#' Build Prediction Model for Methylation Sites
#'
#' Constructs models to predict methylation as a function of single-nucleotide polymorphisms (SNPs)
#' across different window sizes surrounding methylation sites. The function processes a range of
#' methylation sites, fits a model for each, and compiles the results into a MethylationScaff object.
#'
#' @param methInput A MethylationInput object containing SNP and methylation data.
#' @param window_sizes A vector of integers representing the sizes of windows surrounding methylation sites.
#' @param chunk1 The starting index of methylation sites to process.
#' @param chunk2 The ending index of methylation sites to process.
#' @param n_fold Integer, the number of folds for cross-validation.
#' @param scaffoldIdentifier Character, a unique identifier for the scaffold being processed.
#' @param outdir Character, directory path to save the output files.
#' @param maf Numeric, minor allele frequency threshold for filtering SNPs.
#' @param ... Additional arguments passed to `glmnet_tune_alpha` and `cv_eval`.
#' @inheritParams glmnet_tune_alpha
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
#' build_prediction_model(methInput, c(1000, 2000), 1, 10, 5, "double", "scaffold1", tempdir(), TRUE)
#' }
#'
fit_MWAS_models <- function(methInput, window_sizes, chunk1, chunk2,
                            n_fold, scaffoldIdentifier, outdir, alphas,
                            save_evaluation_results_each_fold,
                            save_glmnet_object, cores_per_alpha, num_cores,
                            allow_inefficient_parallelization, maf, na.action,
                            omit_folds_with_na_r = TRUE,
                            cv_eval_mode = "dynamic", verbose = FALSE, ...) {
  requireNamespace("data.table")

  total_iterations <- length(chunk1:chunk2) * length(window_sizes)
  methBaseModels <- vector("list", total_iterations)
  counter <- 1

  for(i in chunk1:chunk2) {

    methylation <- methInput@methylations[, i]  # Extracting methylation for current site
    meth_site_pos <- methInput@methylations_positions[i]
    
    if(verbose == TRUE) {
      cat(paste0("Running methylation site #", i, " at position ", meth_site_pos))
    }

    for(window_size in window_sizes) {

      # 1. Extract SNPs for this window

      SNPs <- extract_SNPs(methInput,
                           meth_site_pos = meth_site_pos,
                           window_size = window_size,
                           verbose = verbose,
                           maf = maf,
                           na.action = na.action)

      if (is.null(SNPs)) {
        if (verbose) {
          message(paste0("For site at position ", meth_site_pos,
                         ", no SNPs meeting criteria found in the window of size ",
                         window_size, "\n\n"))
        }
        next
      }
      
      # Save the original function before making changes
      original_cv_glmnet_raw <- get("cv.glmnet.raw", envir = asNamespace('glmnet'))

      # Step 1: Set the environment of your custom function to the namespace of glmnet
      environment(my_cv_glmnet_raw) <- asNamespace('glmnet')

      # Step 2: Replace the original function in the glmnet namespace with your custom function
      assignInNamespace("cv.glmnet.raw", my_cv_glmnet_raw, ns = "glmnet")

      # 2. Find best parameters model with glmnet.tune.alpha.
      #   This function will also fit final model on full data with best parameters
      tuning_results <- glmnet_tune_alpha(X = SNPs,
                                          y = methylation,
                                          n_fold = n_fold,
                                          verbose = verbose,
                                          lambda_choice = "1se",  # or pass it if it should be variable
                                          alphas = alphas,
                                          cores_per_alpha = cores_per_alpha,
                                          num_cores = num_cores,
                                          allow_inefficient_parallelization = allow_inefficient_parallelization)
      
      # Step 3: Revert to the default version of the function
      assignInNamespace("cv.glmnet.raw", original_cv_glmnet_raw, ns = "glmnet")
      
      if (is.null(tuning_results)) {
          message(paste0("For site at position ", meth_site_pos,
                         "with window size ", window_size,
                         "early termination during k-fold fitting due to low R",
                         " for all alpha values tested\n\n"))
        next
      }

      # 3. Run `cv_eval` (formerly named `cv.pred`) to obtain metrics for model performance

      if(cv_eval_mode == "static"){
        evaluation_results <- cv_eval(X = SNPs,
                                      y = methylation,
                                      n_fold = n_fold,
                                      cv_eval_mode = "static",
                                      best_alpha = tuning_results$para$alpha,
                                      best_lambda = tuning_results$para$lambda,
                                      verbose = verbose,
                                      omit_folds_with_na_r = omit_folds_with_na_r,
                                      ...)
      }
      if(cv_eval_mode == "dynamic"){
        evaluation_results <- cv_eval(X = SNPs,
                                      y = methylation,
                                      n_fold = n_fold,
                                      cv_eval_mode = "dynamic",
                                      verbose = verbose,
                                      alphas = alphas,
                                      cores_per_alpha, num_cores,
                                      allow_inefficient_parallelization,
                                      omit_folds_with_na_r = omit_folds_with_na_r,
                                      ...)
      }

      if(save_evaluation_results_each_fold == TRUE){
        evaluation_results_to_save <- evaluation_results
      } else { # only save mean results across folds
        evaluation_results_to_save <- evaluation_results$eval_cvm
      }

      if(save_glmnet_object == TRUE){
        model_to_save <- tuning_results$model
      } else { # only save model coefficients
        model_to_save <- NULL
      }
      
      #recover()
      
      # Set up predict on full model here
      pred <- predict(tuning_results$model, SNPs)
      
      if(length(levels(factor(pred))) > 1){
        r <- cor(pred, methylation)
      } else {
        r <- NA
      }
      
      mse <- mean((pred - methylation)^2)
      
      full_model_metrics <- c(r = r, mse = mse)
      
      methBase <- new("MethylationBase",
                      methylationPosition = meth_site_pos,
                      windowSize = window_size,
                      n_SNPs = ncol(SNPs),
                      alpha = tuning_results$para$alpha,
                      lambda = tuning_results$para$lambda,
                      glmnetModel = model_to_save,
                      snpWeights = tuning_results$features,
                      intercept = tuning_results$model$a0,
                      evaluation_results = evaluation_results_to_save,
                      cv_eval_mode = cv_eval_mode,
                      full_model_metrics = full_model_metrics)

      methBaseModels[[counter]] <- methBase
      counter <- counter + 1
    }
  }

  # # Calculate size of each slot
  # model_slots <- slotNames(methBase)
  # sizes <- sapply(model_slots, function(slot) object.size(slot(methBase, slot)))

  if (counter >= 2 + length(methBaseModels)) {
    stop("This should not happen. There is a bug; 'counter' should not outpace object 'methBaseModels' by >1.")
  } else if (counter == length(methBaseModels) + 1) {
    # This is okay, proceed without further action
  } else if(!is.null(methBaseModels[[counter]])){
    stop("This should not happen. There is a bug involving the exclusion of windows with no SNPs.")
  }

  #recover()
  
  if(!all(sapply(methBaseModels, is.null))){
    methScaff <- new("MethylationScaff",
                     scaffoldIdentifier = scaffoldIdentifier,
                     models = methBaseModels[1:(counter - 1)])
    
    saveMethylationScaff(methScaff, outputDir = outdir)
    # data.table::fwrite(convertToDataFrame(methScaff),
    #                    file.path(outdir, paste0(scaffoldIdentifier, ".csv")), row.names = FALSE)
  } else {
    # This is not very good coding practice
    #. but it's an acceptable patch since we're relying on checking for RDS files
    #. to see if job completed or not
    message("No models were fitted.")
    fileName <- paste0(scaffoldIdentifier, "-empty.rds")
    filePath <- file.path(outdir, fileName)
    
    filler <- "No models in region"
    
    saveRDS(filler, file = filePath)
    
    message(paste("(Dummy) MethylationScaff object saved to", filePath))
  }

}
