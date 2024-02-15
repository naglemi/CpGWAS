#' MethylationInput Class
#'
#' Encapsulates the processing of input data for methylation analysis in a genomic context. This class facilitates
#' the integration of methylation data, SNP (Single Nucleotide Polymorphism) data in PLINK2 binary format, and
#' covariate information. It supports operations such as data loading, preprocessing, and alignment for subsequent
#' analysis stages.
#'
#' @slot methylations A matrix containing methylation beta values or intensities for genomic positions.
#' @slot methylations_positions An integer vector indicating the genomic positions of methylation sites.
#' @slot genotype_IDs A character vector listing the IDs of genotypes corresponding to the samples in the methylation data.
#' @slot pvar_pointer An external pointer to variant information, compatible with PLINK2 binary format.
#' @slot pvar_dt A data table summarizing variant information, typically including variant IDs, chromosomes, and positions.
#' @slot pgen An external pointer to genotype information, facilitating access to SNP data in PLINK2 binary format.
#' @slot psam A data frame containing sample information extracted from PLINK2 sample (PSAM) file.
#' @slot cov A matrix of covariates to be considered in the analysis, potentially including age, sex, and batch effects.
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
#' Creates a new MethylationInput object by loading and processing methylation and SNP data from specified paths.
#' This method integrates methylation data from a BSseq object and SNP data in PLINK2 binary format, aligning
#' both data types alongside specified covariates for comprehensive analysis.
#'
#' @param BSseq_obj A BSseq object containing methylation data.
#' @param snp_data_path A character string specifying the path to the directory containing SNP data in PLINK2 binary format.
#' @param no_cores Integer, the number of cores to use for parallel operations (defaults to available core count).
#' @param start_site Optional, the start site for subsetting methylation data based on genomic positions.
#' @param end_site Optional, the end site for subsetting methylation data based on genomic positions.
#' @return An object of class MethylationInput, ready for methylation analysis.
#'
#' @import pgenlibr
#' @importFrom tools file_path_sans_ext  # To manipulate file paths
#' @importFrom data.table fread  # To read data tables efficiently
#' @importFrom SummarizedExperiment colData, rowRanges  # For working with experimental data summaries
#' @importFrom GenomicRanges ranges, start  # For handling genomic ranges
#'
#' @examples
#' # Example of initializing a MethylationInput object with methylation and SNP data
#' pgen_path = system.file("extdata", "chr1_sample_subset.pgen", package = "YourPackage")
#' data(chr1_methylation_sample_subset, package = "YourPackage")
#' methInput <- new("MethylationInput", BSseq_obj = BSseq_sample, snp_data_path = pgen_path)
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

setGeneric("sampleMethylationSites", function(object, num_sites) {
  standardGeneric("sampleMethylationSites")
})

setMethod(
  "sampleMethylationSites",
  "MethylationInput",
  function(object, num_sites) {
    if (num_sites <= 0 || num_sites > length(object@methylations_positions)) {
      stop("The number of sites must be greater than 0 and less than or equal to the total number of methylation sites.")
    }
# 
    set.seed(42)
    selected_indices <- sample(x = 1:length(object@methylations_positions),
                               size = num_sites, replace = FALSE)
    
    # Subset methylations and methylations_positions using selected_indices
    object@methylations <- object@methylations[, selected_indices, drop = FALSE]
    object@methylations_positions <- object@methylations_positions[selected_indices]
    
    object
  }
)



#' Reinitialize MethylationInput Object
#'
#' Reloads and updates a MethylationInput object from a saved RDS file, refreshing its SNP data links
#' to reflect new paths or updates in the SNP data stored in PLINK2 binary format.
#'
#' @param rds_path Character string specifying the path to the RDS file containing a previously saved MethylationInput object.
#' @param snp_data_path Character string specifying the new path to the SNP data in PLINK2 binary format.
#' @param no_cores Integer specifying the number of cores to use for parallel operations (defaults to the number of available cores).
#' @return A reinitialized MethylationInput object with updated links to SNP data.
#'
#' @export
#'
#' @examples
#' # Reinitialize a MethylationInput object with a new SNP data path
#' reinitMethInput <- reinitializeMethylationInput("path/to/saved/object.rds", "new/path/to/snp_data", no_cores = 4)
#'
reinitializeMethylationInput <- function(rds_path, snp_data_path, no_cores = detectCores()) {
  if (!file.exists(rds_path)) {
    stop("RDS file does not exist: ", rds_path)
  }
  
  # Load the MethylationInput object from RDS
  loadedObject <- readRDS(rds_path)
  
  # Validate loaded object
  if (!inherits(loadedObject, "MethylationInput")) {
    stop("Loaded object is not a MethylationInput.")
  }
  
  # Re-setup the paths
  pgen_path <- gsub(snp_data_path, pattern = "pvar", replacement = "pgen")
  pvar_path <- gsub(snp_data_path, pattern = "pgen", replacement = "pvar")
  psam_path <- gsub(pvar_path, pattern = "pvar", replacement = "psam")
  
  if (!file.exists(pgen_path) || !file.exists(pvar_path) || !file.exists(psam_path)) {
    stop("One or more SNP data files not found at the specified paths.")
  }
  
  # Reinitialize external pointers
  loadedObject@pvar_pointer <- pgenlibr::NewPvar(pvar_path)
  loadedObject@pvar_dt <- fread(pvar_path)[, 1:3]
  loadedObject@pgen <- pgenlibr::NewPgen(pgen_path, pvar = loadedObject@pvar_pointer)
  loadedObject@psam <- fread(psam_path)
  
  # Reinitialize genotype_IDs based on intersection with methylations
  psam_in_wgbs <- loadedObject@psam[which(loadedObject@psam$`#IID` %in% rownames(loadedObject@methylations))]
  genotype_IDs <- psam_in_wgbs$`#IID`
  genotype_IDs <- intersect(rownames(loadedObject@methylations), genotype_IDs)
  loadedObject@genotype_IDs <- genotype_IDs[order(genotype_IDs)]
  
  # Ensure methylations are filtered and ordered according to the new genotype_IDs, if necessary
  loadedObject@methylations <- loadedObject@methylations[which(rownames(loadedObject@methylations) %in% genotype_IDs), ]
  
  return(loadedObject)
}

    
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

sanityCheckMethylationInput <- function(methInput) {
  # Extract methylations_positions
  positions <- methInput@methylations_positions
  
  # Extract column names from methylations, removing the 'pos_' prefix
  colnames_without_prefix <- gsub("pos_", "", colnames(methInput@methylations))
  
  # Convert column names to numeric for comparison
  colnames_as_numeric <- as.numeric(colnames_without_prefix)
  
  # Check if the positions match the numeric column names
  all_equal <- all(positions == colnames_as_numeric)
  
  if (all_equal) {
    message("Sanity check passed: Positions match column names.")
  } else {
    stop("Sanity check failed: Positions do not match column names.")
  }
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
