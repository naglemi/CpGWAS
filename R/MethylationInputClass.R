#' MethylationInput Class
#'
#' This class encapsulates the input data processing for methylation analysis.
#' It includes methods for loading and processing methylation and SNP data, along with covariates.
#'
#' @slot methylations Matrix of methylation data.
#' @slot genotype_IDs Vector of genotype IDs.
#' @slot pvar_pointer Object representing variant information.
#' @slot pvar_dt Data table of variant information.
#' @slot pgen Object representing genotype information.
#' @slot psam Data table of sample information.
#' @slot cov Data frame of covariates.
#'
#' @export
setClass(
  "MethylationInput",
  slots = c(
    methylations = "matrix",
    methylations_positions = "integer",
    genotype_IDs = "character",
    pvar_pointer = "ANY",
    pvar_dt = "ANY",
    pgen = "ANY",
    psam = "data.frame",
    cov = "matrix"
  )
)

#' Constructor for MethylationInput Class
#'
#' Initializes a MethylationInput object by processing methylation and SNP data.
#'
#' @param BSseq_obj BSseq The BSseq object containing methylation data.
#' @param snp_data_path character Path to the SNP data.
#' @param args list Additional arguments.
#' @return MethylationInput object.
#'
#' @import pgenlibr
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom data.table fread
#' @importFrom pgenlibr NewPvar NewPgen
#' @importFrom SummarizedExperiment colData rowRanges
#' @importFrom GenomicRanges ranges start
#'
#' @examples
#'
#' # Where is SNP data in pgen format to be loaded?
#' pgen_path = system.file("extdata", "chr1_sample_subset.pgen", package = "CpGWAS")
#'
#' data(chr1_methylation_sample_subset, package = "CpGWAS") # Load BSseq_sample
#'
#' methInput <- new("MethylationInput", BSseq_obj = BSobj_sample,
#'                  snp_data_path = pgen_path)
#'
setMethod(
  "initialize",
  "MethylationInput",
  function(.Object, BSseq_obj, snp_data_path, no_cores = detectCores(), start_site = NULL, end_site = NULL) {
    if (is.null(snp_data_path) || is.null(BSseq_obj)) {
      stop("A BSseq object and the path to SNP data are required.")
    }

    if (!inherits(BSseq_obj, "BSseq")) {
      stop("BSseq_obj must be a BSseq object.")
    }
    
    scaffold_name <- tools::file_path_sans_ext(basename(snp_data_path))
    pgen_path <- gsub(snp_data_path, pattern = "pvar", replacement = "pgen")
    pvar_path <- gsub(snp_data_path, pattern = "pgen", replacement = "pvar")
    psam_path <- gsub(pvar_path, pattern = "pvar", replacement = "psam")
    
    if (!file.exists(pgen_path) || !file.exists(pvar_path) || !file.exists(psam_path)) {
      stop("One or more SNP data files not found at the specified paths.")
    }
    
    .Object@pvar_pointer <- pgenlibr::NewPvar(pvar_path)
    .Object@pvar_dt <- fread(pvar_path)[, 1:3]
    .Object@pgen <- pgenlibr::NewPgen(pgen_path, pvar = .Object@pvar_pointer)
    .Object@psam <- fread(psam_path)
    
    #recover()
    
    if (!is.null(start_site) && !is.null(end_site)) {
      
      .Object@methylations <- t(as.matrix(getMeth(BSseq_obj, type = "smooth", what = "perBase")))[, (start_site:end_site)]
      .Object@methylations_positions <- start(ranges(granges(BSseq_obj)))[(start_site:end_site)]     
      
    } else {
      .Object@methylations <- t(as.matrix(getMeth(BSseq_obj, type = "smooth", what = "perBase")))
      .Object@methylations_positions <- start(ranges(granges(BSseq_obj)))
    }
    
    colnames(.Object@methylations) <- paste0("pos_", .Object@methylations_positions)
    

    psam_in_wgbs <- .Object@psam[which(.Object@psam$`#IID` %in% rownames(.Object@methylations))]
    genotype_IDs <- psam_in_wgbs$`#IID`
    genotype_IDs <- intersect(rownames(.Object@methylations), genotype_IDs)
    .Object@genotype_IDs <- genotype_IDs[order(genotype_IDs)]

    .Object@cov <- processCovariates(dataFrame = colData(BSseq_obj),
                                     colsToExclude = c("ID.", "DNum", "brnum", "BrNum", "brnumerical"),
                                     genotype_IDs = genotype_IDs)

    .Object@methylations <- .Object@methylations[which(rownames(.Object@methylations) %in% genotype_IDs), ]
    .Object@methylations <- regress_out_cov_parallel(.Object@methylations, .Object@cov,
                                                     no_cores = no_cores)
    .Object
  }
)





processCovariates <- function(dataFrame,
                              colsToExclude = c("ID.", "DNum", "BrNum",
                                                "brnum", "brnumerical"),
                              genotype_IDs = NULL) {

  if(is.null(genotype_IDs)){
    stop(paste0("`processCovariates must receive a vector of genotype IDs",
                " with order matching methylation trait file."))
  }

  # Sanity check: ensure rownames match brnum
  if (!all(rownames(dataFrame) == dataFrame$brnum)) {
    stop("Row names do not match 'brnum'")
  }

  # Exclude specific columns
  dataFrame <- dataFrame[, !(names(dataFrame) %in% colsToExclude)]

  # Exclude columns with only one factor level and print a warning
  single_factor <- sapply(dataFrame, function(x) length(levels(factor(x))) == 1)

  if (any(single_factor)) {
    # warning("Removing columns with only one factor level: ",
    #         paste(names(dataFrame)[single_factor], collapse = ", "))
    dataFrame <- dataFrame[, !single_factor]
  }

  dataFrame <- dataFrame[which(rownames(dataFrame) %in% genotype_IDs), ]
  dataFrame <- dataFrame[order(rownames(dataFrame)),]

  if(!all(rownames(dataFrame) == genotype_IDs)){
    stop("Mismatch with covariate data and genotype IDs")
  }

  # Get the indices of the factor columns
  factor_columns <- sapply(dataFrame, is.factor)

  # One-hot encode the factor columns using model.matrix
  one_hot_encoded <- model.matrix(~., data = dataFrame[, factor_columns])

  # Combine the one-hot encoded columns with the numeric columns
  final_data <- as.matrix(cbind(one_hot_encoded, dataFrame[, !factor_columns]))

  return(final_data)
}

reorder_and_filter_geno <- function(geno, genotype_IDs) {
  # Filter out rows from geno that are not in genotype_IDs
  idx_geno <- which(row.names(geno) %in% genotype_IDs)
  geno_filtered <- geno[idx_geno, ]

  # Check if both geno and genotype_IDs have the same row names
  if (length(genotype_IDs) != nrow(geno_filtered) || !all(genotype_IDs %in% row.names(geno_filtered))) {
    stop("Row names do not match 100% between geno and genotype_IDs")
  }

  # Find the matching indices for reordering
  match_indices <- match(genotype_IDs, row.names(geno_filtered))

  # Reorder geno matrix to match the order of genotype_IDs
  reordered_geno <- geno_filtered[match_indices, ]

  # Edge case where we only have one SNP in window
  if(!is.matrix(geno_filtered)){
    geno_filtered <- as.matrix(geno_filtered)
  }

  return(reordered_geno)
}

# regress_out_cov <- function(methylations, cov, n_benchmarks = NULL) {
#   print("We just entered regress_out_cov()")
#
#   if(is.null(methylations)){
#     stop("Error: methylation data not found")
#   }
#
#   cat("Dimensions of methylations: ", dim(methylations), "\n")
#
#   # Creating the model formula
#   colnames(cov) <- gsub("\\(Intercept\\)", "Intercept", colnames(cov))
#   cov <- as.data.frame(cov)
#   model_formula <- as.formula(paste("y ~ ", paste(colnames(cov), collapse=" + ")))
#
#   n_tests <- if (is.null(n_benchmarks)) ncol(methylations) else n_benchmarks
#   residuals_matrix <- matrix(NA, nrow = nrow(methylations), ncol = n_tests)
#
#   for(i in 1:n_tests) {
#     y <- methylations[, i]
#     lm_model <- lm(model_formula, data = cbind(y, cov))
#     residuals_matrix[, i] <- residuals(lm_model)
#   }
#
#   cat("Residuals computed for ", n_tests, " tests.\n")
#   return(residuals_matrix)
# }


#' Regress trait (methylation) values on covariates and output residuals
#'
#' @param methylations matrix of methylation rates, n x m for n samples, m sites
#' @param cov_matrix matrix of covariates, n x p for n samples, p covariates
#'
#' @return matrix of residuals, n x m for n samples, m sites
#' @export
#'
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
#'
#' @examples
#' \dontrun{
#' adjusted_methylations <- regress_out_cov_parallel(methylations, cov_matrix)
#' }
regress_out_cov_parallel <- function(methylations, cov_matrix,
                                     no_cores = detectCores() - 1) {
  
  residuals_computation <- function(chunk, cov_matrix, pseudoinv) {
    fitted_values <- cov_matrix %*% (pseudoinv %*% chunk)
    chunk - fitted_values
  }

  if(is.null(methylations)){
    stop("Error: methylation data not found")
  }

  pseudoinv <- solve(t(cov_matrix) %*% cov_matrix) %*% t(cov_matrix)

  chunk_size <- ceiling(ncol(methylations) / no_cores)
  
  if(no_cores > 1 && (ncol(methylations) > no_cores)){
    chunks <- lapply(1:min(no_cores, ncol(methylations)), function(i) {
      start_col <- (i - 1) * chunk_size + 1
      end_col <- min(i * chunk_size, ncol(methylations))
      data = methylations[, start_col:end_col]
    })
    results <- future_lapply(chunks, residuals_computation, cov_matrix, pseudoinv)
    residuals_matrix <- do.call(cbind, results)
  } else {
    residuals_matrix <- residuals_computation(methylations, cov_matrix, pseudoinv)
  }
  
  if(!all(colnames(methylations) == colnames(residuals_matrix))) {
    stop("Error: `regress_out_cov` is broken. Residuals in different order than inputs.")
  }

  return(residuals_matrix)
}
