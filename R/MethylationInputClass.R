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
    methylations = "ANY",
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
  function(.Object, BSseq_obj, snp_data_path, start_site = NULL, end_site = NULL) {
    if (is.null(snp_data_path) || is.null(BSseq_obj)) {
      stop("A BSseq object and the path to SNP data are required.")
    }

    if (!inherits(BSseq_obj, "BSseq")) {
      stop("BSseq_obj must be a BSseq object.")
    }

    methylations_full <- t(as.matrix(getMeth(BSseq_obj, type = "smooth", what = "perBase")))
    colnames(methylations_full) <- paste0("pos_", GenomicRanges::start(GenomicRanges::ranges(SummarizedExperiment::rowRanges(BSseq_obj))))
    
    scaffold_name <- tools::file_path_sans_ext(basename(snp_data_path))
    pgen_path <- gsub(snp_data_path, pattern = "pvar", replacement = "pgen")
    pvar_path <- gsub(snp_data_path, pattern = "pgen", replacement = "pvar")
    psam_path <- gsub(pvar_path, pattern = "pvar", replacement = "psam")

    if (!file.exists(pgen_path) || !file.exists(pvar_path) || !file.exists(psam_path)) {
      stop("One or more SNP data files not found at the specified paths.")
    }

    pvar_pointer <- pgenlibr::NewPvar(pvar_path)
    pvar_dt_full <- fread(pvar_path)
    pgen <- pgenlibr::NewPgen(pgen_path, pvar = pvar_pointer)
    psam <- fread(psam_path)
    psam_in_wgbs <- psam[which(psam$`#IID` %in% rownames(methylations_full))]
    genotype_IDs <- psam_in_wgbs$`#IID`
    genotype_IDs <- intersect(rownames(methylations_full), genotype_IDs)
    genotype_IDs <- genotype_IDs[order(genotype_IDs)]

    cov <- processCovariates(dataFrame = colData(BSseq_obj),
                             colsToExclude = c("ID.", "DNum", "brnum", "BrNum", "brnumerical"),
                             genotype_IDs = genotype_IDs)

    methylations_full <- methylations_full[which(rownames(methylations_full) %in% genotype_IDs), ]
    methylations_full <- regress_out_cov_parallel(methylations_full, cov)
    
    if (!is.null(start_site) && !is.null(end_site)) {
      
      if (start_site > 1) {
        methylations_full[, 1:(start_site-1)] <- 0
      }
      
      if (end_site < ncol(methylations_full)) {
        methylations_full[, (end_site+1):ncol(methylations_full)] <- 0
      }

      .Object@methylations <- Matrix(methylations_full, sparse = TRUE)


      # # Sparse matrix for methylations
      # for (col in start_site:end_site) {
      #   rows_with_data <- which(methylations_full[, col] != 0)
      #   .Object@methylations[rows_with_data, col] <- methylations_full[rows_with_data, col]
      # }
      # 
      # # Assuming pvar_dt_full is a data.frame or data.table
      # for (row in start_site:end_site) {
      #   cols_with_data <- which(!is.na(pvar_dt_full[row, ]))
      #   # Convert the subset of pvar_dt_full to a numeric vector
      #   numeric_values <- as.numeric(pvar_dt_full[row, ..cols_with_data])
      #   # Assign to sparse matrix
      #   .Object@pvar_dt[row, cols_with_data] <- numeric_values
      # }

    } else {
      .Object@methylations <- methylations_full
    }
    
    .Object@pvar_dt <- pvar_dt_full[, 1:3]
    .Object@genotype_IDs <- genotype_IDs
    .Object@pvar_pointer <- pvar_pointer
    .Object@pgen <- pgen
    .Object@psam <- psam
    .Object@cov <- cov
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
#' @param n_benchmarks integer number of sites to compute residuals for, if benchmarking
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
regress_out_cov_parallel <- function(methylations, cov_matrix, n_benchmarks = NULL,
                                     no_cores = detectCores() - 1) {

  if(is.null(methylations)){
    stop("Error: methylation data not found")
  }

  #cat("Dimensions of methylations: ", dim(methylations), "\n")

  #cat("Dimensions of cov_matrix: ", dim(cov_matrix), "\n")

  pseudoinv <- solve(t(cov_matrix) %*% cov_matrix) %*% t(cov_matrix)
  #cat("Dimensions of pseudoinv: ", dim(pseudoinv), "\n")

  n_tests <- if (is.null(n_benchmarks)) ncol(methylations) else n_benchmarks
  chunk_size <- ceiling(n_tests / no_cores)
  chunks <- lapply(1:no_cores, function(i) {
    start_col <- (i - 1) * chunk_size + 1
    end_col <- min(i * chunk_size, n_tests)
    methylations[, start_col:end_col]
  })

  residuals_computation <- function(chunk, cov_matrix, pseudoinv) {
    tryCatch({
      # Compute the fitted values for the entire matrix
      fitted_values <- cov_matrix %*% (pseudoinv %*% chunk)
      # Compute the residuals for the entire matrix
      residuals <- chunk - fitted_values
      return(residuals)
    }, error = function(e) {
      cat("Error occurred: ", e$message, "\n")
      cat("Dimensions of cov_matrix: ", dim(cov_matrix), "\n")
      cat("Dimensions of pseudoinv: ", dim(pseudoinv), "\n")
      cat("Dimensions of chunk: ", dim(chunk), "\n")
      matrix(NA, nrow = nrow(chunk), ncol = ncol(chunk))
    })
  }

  start_time <- Sys.time()
  results <- future_lapply(chunks, residuals_computation, cov_matrix, pseudoinv)
  residuals_matrix <- do.call(cbind, results)

  if (!is.null(n_benchmarks)) {
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    cat("Elapsed time: ", elapsed_time, "\n")
  }

  return(residuals_matrix)
}
