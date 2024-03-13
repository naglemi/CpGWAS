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
#' @param start_site Optional, the start site for subsetting methylation data based on indices in the file.
#' @param end_site Optional, the end site for subsetting methylation data based on indices in the file.
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
  function(.Object, BSseq_obj, snp_data_path, no_cores = detectCores(), cov_path = NULL, start_site = NULL, end_site = NULL) {

    if (is.null(snp_data_path) || is.null(BSseq_obj)) {
      stop("A BSseq object and the path to SNP data are required.")
    }
    
    if (!inherits(BSseq_obj, "BSseq")) {
      stop("BSseq_obj must be a BSseq object.")
    }
    
    if (is.null(cov_path)) {
      message("No covariate file provided. Covariates will not be included in the analysis.")
    } else {
      if (!file.exists(cov_path)) {
        stop("Covariate file does not exist: ", cov_path)
      }
      if (file.exists(cov_path)) {
        cov_data <- fread(cov_path)
#        colnames(cov_data) <- gsub("\\(Intercept\\)", "Intercept", colnames(cov_data))
      }
    }
    
    paths <- constructSNPFilePaths(snp_data_path)
    checkFilesExist(c(paths$pgen_path, paths$pvar_path, paths$psam_path))
    
    snp_data <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
    .Object@pvar_pointer <- snp_data$pvar_pointer
    .Object@pvar_dt <- snp_data$pvar_dt
    .Object@pgen <- snp_data$pgen
    .Object@psam <- snp_data$psam
    
    meth_data <- processMethylationData(BSseq_obj, start_site, end_site)
    .Object@methylations <- meth_data$methylations
    .Object@methylations_positions <- meth_data$methylations_positions
    colnames(.Object@methylations) <- paste0("pos_", .Object@methylations_positions)
    
    .Object@cov <- processCovariates(cov_path)

    .Object@genotype_IDs <- selectGenotypeIDs(.Object@psam, .Object@methylations,
                                              .Object@cov)
    
    # Filter and order methylations by genotype IDs
    if(!identical(rownames(.Object@methylations), .Object@genotype_IDs)){
      .Object@methylations <- .Object@methylations[which(rownames(.Object@methylations) %in% .Object@genotype_IDs), ]
    }

    if(!identical(rownames(.Object@cov), .Object@genotype_IDs)){
      .Object@cov <- .Object@cov[which(rownames(.Object@cov) %in% .Object@genotype_IDs), ]
    }
    
    #.Object@cov <- processCovariatesFromBSseq(dataFrame = colData(BSseq_obj),
    #                                 colsToExclude = c("ID.", "DNum", "brnum", "BrNum", "brnumerical"),
    #                                 genotype_IDs = .Object@genotype_IDs)
    
    # Regress out covariates in parallel
    .Object@methylations <- regress_out_cov_parallel(.Object@methylations, .Object@cov, no_cores = no_cores)
    
    return(.Object)
    }
  )

#' @export
setGeneric("sampleMethylationSites", function(object, num_sites, seed) {
  standardGeneric("sampleMethylationSites")
})

setMethod(
  "sampleMethylationSites",
  "MethylationInput",
  function(object, num_sites, seed = NULL) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (num_sites <= 0 || num_sites > length(object@methylations_positions)) {
      recover()
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

#' @export
setGeneric(
  "validatePositionOverlap",
  signature = "object",
  def = function(object, window_size) {
    standardGeneric("validatePositionOverlap")
  }
)

#' Validate Position Overlap Between Methylation Sites and SNPs
#'
#' This method checks if there is an overlap between methylation sites and SNP positions
#' given a specified window size. It ensures that the methylation data being used for analysis
#' has relevant SNP data within the defined proximity.
#' 
#' @param object MethylationInput object containing methylation positions and SNP data.
#' @param window_size Numeric value defining the window size around methylation sites to check for SNP overlap.
#' 
#' @return Invisible. The function stops with an error message if there is no overlap within the specified window size.
#' 
#' @examples
#' # Assuming `methInput` is a MethylationInput object
#' validatePositionOverlap(methInput, 1000)
#' 
#' @details
#' The function calculates the maximum and minimum positions for methylation sites and SNP data,
#' adjusted by the window size, to determine if any SNP falls within the vicinity of the methylation sites.
#' If not, then there's no overlap and the function stops with an error message.
#'
#' @note
#' This method is specific to MethylationInput objects. Ensure that your MethylationInput object
#' is properly initialized with valid methylation positions and SNP data before calling this method.
setMethod(
  "validatePositionOverlap",
  signature(object = "MethylationInput"),
  definition = function(object, window_size) {
    last_possible_position_given_methylation <- max(object@methylations_positions) + window_size
    first_possible_position_given_methylation <- min(object@methylations_positions) - window_size
    
    last_possible_position_given_SNPs <- max(object@pvar_dt$POS)
    first_possible_position_given_SNPs <- min(object@pvar_dt$POS)
    
    if (last_possible_position_given_methylation < first_possible_position_given_SNPs) {
      stop(paste0("There is no overlap between selected methylation sites and given SNP data\n",
                  "with window size of ", window_size, ".\n\n",
                  "The last possible position given methylation sites is ", last_possible_position_given_methylation,
                  ", and the first possible position given SNPs is ", first_possible_position_given_SNPs, "."))
    }
    
    if (first_possible_position_given_methylation > last_possible_position_given_SNPs) {
      stop(paste0("There is no overlap between selected methylation sites and given SNP data\n",
                  "with window size of ", window_size, ".\n\n",
                  "The first possible position given methylation sites is ", first_possible_position_given_methylation,
                  ", and the last possible position given SNPs is ", last_possible_position_given_SNPs, "."))
    }
  }
)

#' @export
processCovariates <- function(cov_path){
  cov <- read.csv(cov_path, stringsAsFactors = TRUE)
  rownames(cov) <- cov$ID
  cov <- cov[, -1]
  
  # Get the indices of the factor columns
  factor_columns <- sapply(cov, is.factor)
  
  # One-hot encode the factor columns using model.matrix
  one_hot_encoded <- model.matrix(~., data = cov[, factor_columns])
  
  # Combine the one-hot encoded columns with the numeric columns
  cov <- as.matrix(cbind(one_hot_encoded, cov[, !factor_columns]))
  
  cov <- as.matrix(cov)
}

#' Reinitialize MethylationInput Object
#'
#' Reloads and updates a MethylationInput object from a saved RDS file, refreshing its SNP data links
#' to reflect new paths or updates in the SNP data stored in PLINK2 binary format.
#'
#' @param rds_path Character string specifying the path to the RDS file containing a previously saved MethylationInput object.
#' @param snp_data_path Character string specifying the new path to the SNP data in PLINK2 binary format.
#' @param no_cores Integer specifying the number of cores to use for parallel operations (defaults to the number of available cores).
#' @param start_site Optional, the start site for subsetting methylation data based on indices in the file.
#' @param end_site Optional, the end site for subsetting methylation data based on indices in the file.
#' @return A reinitialized MethylationInput object with updated links to SNP data.
#'
#' @export
#'
#' @examples
#' # Reinitialize a MethylationInput object with a new SNP data path
#' reinitMethInput <- reinitializeMethylationInput("path/to/saved/object.rds", "new/path/to/snp_data", no_cores = 4)
#'
reinitializeMethylationInput <- function(rds_path, snp_data_path, no_cores = detectCores(), start_site = NULL, end_site = NULL) {
  if (!file.exists(rds_path)) {
    stop("RDS file does not exist: ", rds_path)
  }
  
  # Load the MethylationInput object from RDS
  loadedObject <- readRDS(rds_path)
  
  # Validate loaded object
  if (!inherits(loadedObject, "MethylationInput")) {
    stop("Loaded object is not a MethylationInput.")
  }
  
  # Utilize helper functions for SNP path setup and validation
  paths <- constructSNPFilePaths(snp_data_path)
  checkFilesExist(c(paths$pgen_path, paths$pvar_path, paths$psam_path))
  
  # Load SNP data using helper function
  snp_data <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
  loadedObject@pvar_pointer <- snp_data$pvar_pointer
  loadedObject@pvar_dt <- snp_data$pvar_dt
  loadedObject@pgen <- snp_data$pgen
  loadedObject@psam <- snp_data$psam
  
  if (!is.null(start_site) && !is.null(end_site)) {
    loadedObject@methylations <- loadedObject@methylations[, (start_site:end_site)]
    loadedObject@methylations_positions <- loadedObject@methylations_positions[(start_site:end_site)]
  }
  
  # Utilize helper functions for genotype IDs update
  loadedObject@genotype_IDs <- selectGenotypeIDs(loadedObject@psam, loadedObject@methylations)
  
  # Filter and order methylations by genotype IDs
  if(!identical(rownames(loadedObject@methylations), loadedObject@genotype_IDs)){
    loadedObject@methylations <- loadedObject@methylations[which(rownames(loadedObject@methylations) %in% loadedObject@genotype_IDs), ]
  }
  
  if(!identical(rownames(loadedObject@cov), loadedObject@genotype_IDs)){
    loadedObject@cov <- loadedObject@cov[which(rownames(loadedObject@cov) %in% loadedObject@genotype_IDs), ]
  }
  
  return(loadedObject)
}

# Helper function to select genotype IDs based on methylation data
selectGenotypeIDs <- function(psam_data, methylation_data, cov_data = NULL) {
  psam_in_wgbs <- psam_data[which(psam_data$`#IID` %in% rownames(methylation_data))]
  if (!is.null(cov_data)){
    psam_in_wgbs <- psam_in_wgbs[which(psam_in_wgbs$`#IID` %in% rownames(cov_data))]
  }
  genotype_IDs <- psam_in_wgbs$`#IID`
  intersected_genotype_IDs <- intersect(rownames(methylation_data), genotype_IDs)
  ordered_genotype_IDs <- intersected_genotype_IDs[order(intersected_genotype_IDs)]
  
  return(ordered_genotype_IDs)
}

# Helper function to construct SNP file paths
constructSNPFilePaths <- function(base_path) {
  pgen_path <- gsub(base_path, pattern = "pvar", replacement = "pgen")
  pvar_path <- gsub(base_path, pattern = "pgen", replacement = "pvar")
  psam_path <- gsub(pvar_path, pattern = "pvar", replacement = "psam")
  
  list(pgen_path = pgen_path, pvar_path = pvar_path, psam_path = psam_path)
}

# Helper function to check file existence
checkFilesExist <- function(paths) {
  if (!all(file.exists(paths))) {
    stop("One or more SNP data files not found at the specified paths.")
  }
}

# Helper function to load SNP data
loadSNPData <- function(pvar_path, pgen_path, psam_path) {
  pvar_pointer <- pgenlibr::NewPvar(pvar_path)
  list(
    pvar_pointer = pvar_pointer,
    pvar_dt = fread(pvar_path)[, 1:3],
    pgen = pgenlibr::NewPgen(pgen_path, pvar = pvar_pointer),
    psam = fread(psam_path)
  )
}

# Helper function to process methylation data
processMethylationData <- function(BSseq_obj, start_site = NULL, end_site = NULL) {
  if (!is.null(start_site) && !is.null(end_site)) {
    methylations <- t(as.matrix(getMeth(BSseq_obj, type = "smooth", what = "perBase")))[, (start_site:end_site)]
    methylations_positions <- start(ranges(granges(BSseq_obj)))[(start_site:end_site)]
  } else {
    methylations <- t(as.matrix(getMeth(BSseq_obj, type = "smooth", what = "perBase")))
    methylations_positions <- start(ranges(granges(BSseq_obj)))
  }
  
  list(methylations = methylations, methylations_positions = methylations_positions)
}

    
processCovariatesFromBSseq <- function(dataFrame,
                              colsToExclude = c("ID.", "DNum", "BrNum",
                                                "brnum", "brnumerical"),
                              genotype_IDs = NULL) {

  if(is.null(genotype_IDs)){
    stop(paste0("`processCovariatesFromBSseq must receive a vector of genotype IDs",
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

#' @export
regress_out_cov <- function(methylations, cov, n_benchmarks = NULL) {
  #print("We just entered regress_out_cov()")

  if(is.null(methylations)){
    stop("Error: methylation data not found")
  }

  #cat("Dimensions of methylations: ", dim(methylations), "\n")

  # Creating the model formula
  colnames(cov) <- gsub("\\(Intercept\\)", "Intercept", colnames(cov))
  cov <- as.data.frame(cov)
  model_formula <- as.formula(paste("y ~ ", paste(colnames(cov), collapse=" + ")))

  n_tests <- if (is.null(n_benchmarks)) ncol(methylations) else n_benchmarks
  residuals_matrix <- matrix(NA, nrow = nrow(methylations), ncol = n_tests)

  for(i in 1:n_tests) {
    y <- methylations[, i]
    lm_model <- lm(model_formula, data = cbind(y, cov))
    residuals_matrix[, i] <- residuals(lm_model)
  }
  
  colnames(residuals_matrix) <- colnames(methylations)
  rownames(residuals_matrix) <- rownames(methylations)

  #cat("Residuals computed for ", n_tests, " tests.\n")
  return(residuals_matrix)
}


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
regress_out_cov_parallel <- function(methylations, cov_matrix, no_cores = detectCores() - 1) {

  if(!all(rownames(methylations) == rownames(cov_matrix))) {
    stop("Error: Row names of methylation data and covariates do not match.")
  }
  
  residuals_computation <- function(chunk, cov_matrix, pseudoinv) {
    fitted_values <- cov_matrix %*% (pseudoinv %*% chunk)
    return(chunk - fitted_values)
  }
  
  if(is.null(methylations)){
    stop("Error: methylation data not found")
  }

  pseudoinv <- solve(t(cov_matrix) %*% cov_matrix) %*% t(cov_matrix)
  
  # Calculate the number of columns per chunk
  # Ensure there's no division by zero or incorrect chunk size when no_cores > ncol(methylations)
  no_cores <- min(no_cores, ncol(methylations))
  chunk_size <- max(1, ceiling(ncol(methylations) / no_cores))
  
  chunks <- list()
  for (i in 1:no_cores) {
    start_col <- (i - 1) * chunk_size + 1
    end_col <- min(i * chunk_size, ncol(methylations))
    if (start_col <= ncol(methylations)) { # Ensure start_col is within the column range
      chunks[[i]] <- methylations[, start_col:end_col, drop = FALSE]
    }
  }
  
  if(no_cores > 1 && (ncol(methylations) > no_cores)){
    results <- future_lapply(chunks, residuals_computation, cov_matrix = cov_matrix, pseudoinv = pseudoinv)
    residuals_matrix <- do.call(cbind, results)
  } else {
    residuals_matrix <- residuals_computation(methylations, cov_matrix, pseudoinv)
  }
  
  # Ensure the column names of the residuals_matrix match those of the original methylations matrix
  if(!all(colnames(methylations) == colnames(residuals_matrix))) {
    stop("Error: Residuals in different order than inputs.")
  }
  
  return(residuals_matrix)
}

