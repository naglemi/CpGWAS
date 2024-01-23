#' Build Prediction Model for Methylation Sites
#'
#' Constructs models to predict methylation as a function of single-nucleotide polymorphisms (SNPs)
#' across different window sizes surrounding methylation sites. The function processes a range of
#' methylation sites, fits a model for each, and compiles the results into a MethylationScaff object.
#'
#' @param BSobj A BSseq object containing methylation data.
#' @param methInput A MethylationInput object containing SNP and methylation data.
#' @param window_sizes A vector of integers representing the sizes of windows surrounding methylation sites.
#' @param chunk1 The starting index of methylation sites to process.
#' @param chunk2 The ending index of methylation sites to process.
#' @param n_fold Integer, the number of folds for cross-validation.
#' @param cv_nesting Character, type of cross-validation nesting ("double" or "triple").
#' @param scaffoldIdentifier Character, a unique identifier for the scaffold being processed.
#' @param outdir Character, directory path to save the output files.
#' @param record_runtime Logical, indicating whether to record the runtime of the cross-validation.
#'
#' @return An object of class MethylationScaff containing the fitted models.
#' @export
#'
#' @importFrom methods new
#' @importFrom stats as.formula coef lm model.matrix na.omit predict residuals
#' @importFrom bsseq getMeth granges
#' @importFrom SummarizedExperiment colData rowRanges
#' @importFrom GenomicRanges ranges start
#'
#' @examples
#' \dontrun{
#' build_prediction_model(BSobj, methInput, c(1000, 2000), 1, 10, 5, "double", "scaffold1", tempdir(), TRUE)
#' }
#'
build_prediction_model <- function(BSobj, methInput, window_sizes, chunk1, chunk2, n_fold, cv_nesting, scaffoldIdentifier, outdir, record_runtime) {
  requireNamespace("data.table")

  methBaseModels <- list()

  for(i in chunk1:chunk2) {
    cat("Processing methylation site index: ", i, "\n")
    methylation <- methInput@methylations[, i]  # Extracting methylation for current site
    meth_site_pos <- start(ranges(granges(BSobj)))[i]

    cat("Testing methylation site at position ", meth_site_pos, "\n")

    for(window_size in window_sizes) {

      cat("Testing window size: ", window_size, "\n")

      SNPs <- extract_SNPs(methInput,
                           meth_site_pos = meth_site_pos,
                           window_size = window_size)

      cat("Dimensions of SNP matrix: ", dim(SNPs), "\n")

      if (is.null(SNPs)) {
        cat(paste0("For site at position ", meth_site_pos,
        ", no SNPs were found in the window of size ", window_size, "\n"))
        next
      }

      cv_results <- switch(cv_nesting,
                           "triple" = cv.pred.3nest(SNPs, methylation, n_fold = n_fold, record_runtime = record_runtime),
                           "double" = cv.pred(SNPs, methylation, n_fold = n_fold, record_runtime = record_runtime)
      )

      methBase <- new("MethylationBase",
                      methylationPosition = meth_site_pos,
                      windowSize = window_size,
                      n_SNPs = ncol(SNPs),
                      glmnetModel = cv_results$model,
                      snpWeights = cv_results$snp_weights,
                      cor = cv_results$cor,
                      mse = cv_results$mse,
                      alpha = cv_results$alpha,
                      lambda = cv_results$lambda,
                      runtime = cv_results$runtime)
      methBaseModels[[length(methBaseModels) + 1]] <- methBase
      cat("\n")
    }
  }

  methScaff <- new("MethylationScaff",
                   scaffoldIdentifier = scaffoldIdentifier,
                   models = methBaseModels)

  saveMethylationScaff(methScaff, outputDir = outdir)
  data.table::fwrite(convertToDataFrame(methScaff),
                     file.path(outdir, paste0(scaffoldIdentifier, ".csv")), row.names = FALSE)

  return(methScaff)
}
