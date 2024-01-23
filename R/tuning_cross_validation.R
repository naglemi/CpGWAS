#' Impute Missing SNP Values with Column Means
#'
#' Handles missing values in SNP data by replacing them with the mean of their respective column.
#' This method is a standard approach in genomic studies to maintain data integrity when dealing with
#' missing data, common in genetic datasets.
#'
#' @param SNPs A matrix or data frame representing SNP data, possibly containing missing values (NAs).
#'             Expected Dimensions: Observations (rows) x SNPs (columns).
#'
#' @return A matrix or data frame with the same dimensions as input, where missing values have been
#'         imputed with column means.
#'
#' @examples
#' \dontrun{
#' imputed_SNPs <- impute_missing_values(SNPs)
#' }
impute_missing_values <- function(SNPs) {
  if (!is.matrix(SNPs) && !is.data.frame(SNPs)) {
    stop("SNPs should be a matrix or data frame.")
  }
  column_means <- colMeans(SNPs, na.rm = TRUE)
  for (column_name in names(SNPs)) {
    SNPs[[column_name]][is.na(SNPs[[column_name]])] <- column_means[column_name]
  }
  return(SNPs)
}


#' Tuning Function for Elastic Net Hyperparameters
#'
#' Performs hyperparameter tuning for elastic net models, specifically tailored for
#' methylome-wide association studies. It iterates over a range of alpha values to find the optimal
#' lambda using cross-validation.
#'
#' @param alpha Numeric value representing the alpha parameter in elastic net, controlling the mix
#'              of L1 and L2 regularization. Expected range: [0, 1].
#' @param X_train Training predictors, a matrix or data frame of SNP data.
#'                Dimensions: Observations (rows) x SNPs (columns).
#' @param yhat_train Numeric vector representing the response variable in the training set.
#'                   Length should match the number of rows in X_train.
#' @param fold_id Numeric vector specifying fold assignments for cross-validation.
#' @param lambda_choice Character string indicating the method for choosing lambda
#'                      ('min' for minimum MSE, '1se' for one-standard-error rule).
#'
#' @return A numeric vector containing the mean squared error, optimal lambda, and alpha.
#'
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @examples
#' \dontrun{
#' optimal_parameters <- tuning_function(alpha, X_train, yhat_train, fold_id, lambda_choice)
#' }
tuning_function <- function(alpha, X_train, yhat_train, fold_id, lambda_choice) {
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Alpha should be a numeric value between 0 and 1.")
  }
  if (!is.matrix(X_train) && !is.data.frame(X_train)) {
    stop("X_train should be a matrix or data frame.")
  }
  if (!is.numeric(yhat_train)) {
    stop("yhat_train should be a numeric vector.")
  }
  if (!is.character(lambda_choice) || !lambda_choice %in% c("min", "1se")) {
    stop("lambda_choice should be 'min' or '1se'.")
  }
  cv <- cv.glmnet(X_train, yhat_train, foldid = fold_id, type.measure = "mse", alpha = alpha)
  lambda_index <- if (lambda_choice == "min") {
    cv$lambda == cv$lambda.min
  } else {
    cv$lambda == cv$lambda.1se
  }
  cv_mean_mse <- cv$cvm[lambda_index]
  lambda <- cv$lambda[lambda_index]
  return(c(cv_mean_mse, lambda, alpha))
}

#' Remove Observations with Missing Response Values
#'
#' Removes any observations in the training set where the response variable (methylation levels)
#' is missing. This function ensures that both predictor (SNPs) and response datasets are aligned,
#' which is crucial in subsequent analyses of methylome-wide association studies.
#'
#' @param yhat_train Numeric vector representing the response variable (methylation levels) in the training set.
#' @param X_train Matrix or data frame of predictors (SNP data) in the training set.
#'                Dimensions: Observations (rows) x SNPs (columns).
#' @param verbose Logical indicating whether to print warnings about removed NAs.
#'
#' @return A list containing two elements: 'yhat_train' and 'X_train', both filtered to exclude
#'         observations with missing response values.
#'
#' @examples
#' \dontrun{
#' filtered_data <- rm_missing_observations(yhat_train, X_train, verbose = TRUE)
#' }
rm_missing_observations <- function(yhat_train, X_train, verbose = FALSE) {
  if (!is.numeric(yhat_train)) {
    stop("yhat_train should be a numeric vector.")
  }
  if (!is.matrix(X_train) && !is.data.frame(X_train)) {
    stop("X_train should be a matrix or data frame.")
  }
  n_na <- sum(is.na(yhat_train))
  if(n_na > 0) {
    if(verbose) {
      warning(paste("There are", n_na, "NA values in yhat_train, which will be removed"))
    }
    idx <- !is.na(yhat_train)
    yhat_train <- yhat_train[idx]
    X_train <- X_train[idx, ]
  }
  return(list(yhat_train = yhat_train, X_train = X_train))
}

#' Extract Non-Zero Coefficients from an Elastic Net Model
#'
#' Post-model fitting in a methylome-wide association study, this function extracts the non-zero
#' coefficients from an elastic net model. These coefficients represent the SNPs that have a
#' significant impact on the response variable (methylation levels), indicating relevant genetic
#' markers.
#'
#' @param model_fit The fitted elastic net model object from 'glmnet'.
#'
#' @return A data frame where each row represents a SNP with a non-zero coefficient.
#'         Columns include SNP names and their corresponding coefficients.
#'
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @examples
#' \dontrun{
#' significant_snps <- extract_non_zero_coefs(fitted_model)
#' }
extract_non_zero_coefs <- function(model_fit) {
  if (!inherits(model_fit, "glmnet")) {
    stop("model_fit must be an object of class 'glmnet'.")
  }
  coef_data <- coef(model_fit)
  valid_coefs <- coef_data[-1]  # Drop intercept and keep non-zero coefficients
  non_zero_indices <- valid_coefs != 0
  data.frame(v = rownames(coef_data)[-1][non_zero_indices], coefs = valid_coefs[non_zero_indices])
}

# Calculate correlation between predicted and observed responses
#
# This function uses a model to predict values, then provides the correlation
# between predicted and observed values, as a measure of model performance.
#
# Args:
#   fit: The fitted elastic net model object.
#   X_test: Subset of SNP matrix assigned for validation.
#   y: Observed response variable (portion in validation dataset).
#
# Returns:
#   The correlation value, or NA if the variance of predictions is zero or NA.

#' Calculate Correlation Between Predicted and Observed Responses
#'
#' Uses a fitted elastic net model to predict values and then computes the correlation
#' between these predicted values and the observed response values. This function is used
#' as a measure of model performance in methylome-wide association studies.
#'
#' @param fit The fitted elastic net model object from 'glmnet'.
#' @param X_test Test set predictors, a subset of SNP matrix assigned for validation.
#'               Dimensions: Observations (rows) x SNPs (columns).
#' @param yhat_test Numeric vector of observed response variable in the validation dataset.
#'
#' @return The correlation value between predicted and observed responses, or NA
#'         if the variance of predictions is zero or NA.
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats cor
#' @importFrom stats var
#'
#' @examples
#' \dontrun{
#' model_performance <- calculate_correlation(fitted_model, X_test, yhat_test)
#' }
calculate_correlation <- function(fit, X_test, yhat_test) {
  if (!inherits(fit, "glmnet")) {
    stop("fit must be an object of class 'glmnet'.")
  }
  if (!is.matrix(X_test) && !is.data.frame(X_test)) {
    stop("X_test should be a matrix or data frame.")
  }
  if (!is.numeric(yhat_test)) {
    stop("yhat_test should be a numeric vector.")
  }
  predY <- predict(fit, X_test)
  if (var(predY) > 0) {
    return(cor(yhat_test, predY))
  } else {
    warning("Variance of predY is 0 or NA, setting correlation to NA")
    return(NA)
  }
}

#' Hyperparameter Tuning for Elastic Net Model
#'
#' This function performs hyperparameter tuning for an elastic net model. It is particularly
#' designed for methylome-wide association studies where the relationship between SNPs (Single
#' Nucleotide Polymorphisms) and methylation levels is being modeled. The function iterates over
#' a range of alpha values to find the optimal lambda using cross-validation.
#'
#' @param X_train A matrix or data frame of SNP data used as predictors in the training set.
#'                Dimensions: Number of Observations (rows) x Number of SNPs (columns).
#' @param yhat_train A numeric vector representing the response variable (methylation levels)
#'                   in the training set. Length should match the number of rows in X_train.
#' @param n_fold Integer, specifying the number of folds for cross-validation.
#' @param verbose Logical, indicating whether to print detailed tuning results.
#' @param lambda_choice Character, specifying the method for choosing lambda
#'                      ('min' for minimum MSE, '1se' for one-standard-error rule).
#' @param alphas Numeric vector, specifying the alpha values to iterate over in tuning.
#'
#' @return A list containing the fitted model, non-zero coefficients, optimal parameters,
#'         and correlation on test data.
#'
#' @importFrom future.apply future_sapply
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @examples
#' \dontrun{
#' fit <- glmnet_tune_alpha(X_train, yhat_train, n_fold = 5, verbose = TRUE,
#'                          lambda_choice = "1se", alphas = seq(0, 1, 0.1))
#' }
glmnet_tune_alpha <- function(X_train, yhat_train, n_fold = 5, verbose = FALSE,
                              lambda_choice = "1se",
                              alphas = seq(0, 1, 0.1)) {
  # Choice of lambda:
  # lambda.min: Minimum CV error; empirically best model.
  #lambda.1se: Largest lambda within 1 SE of minimum; simpler, more generalizable model.
  set.seed(2023)

  if (nrow(X_train) != length(yhat_train)) {
    stop("Number of samples is different")
  }

  results <- rm_missing_observations(yhat_train, X_train, verbose = TRUE)
  yhat_train <- results$yhat_train
  X_train <- results$X_train
  results <- NULL

  fold_id <- sample(rep(1:n_fold, length.out = length(yhat_train)))

  # Pre-allocate data fram
  tuning_results <- data.frame(cv_mean_mse = numeric(length(alphas)),
                     lambda = numeric(length(alphas)),
                     alpha = alphas)

  # Parallel processing
  tuning_results_to_be_parsed <-
    future_sapply(X = alphas, FUN = tuning_function,
                  X_train = X_train, yhat_train = yhat_train,
                  fold_id = fold_id, lambda_choice = lambda_choice,
                  simplify = "array")

  # Assign results to pre-allocated data frame
  tuning_results$cv_mean_mse <- tuning_results_to_be_parsed[1,]
  tuning_results$lambda <- tuning_results_to_be_parsed[2,]

  optimal_model <- tuning_results[which.min(tuning_results$cv_mean_mse),]

  # Debug: Check tuning results
  if (verbose) {
    cat("Tuning results - Lambda:", optimal_model$lambda,
        "Alpha:", optimal_model$alpha, "\n")
  }

  # fit final model on the whole training set
  fit = glmnet(X_train, yhat_train,
               lambda = optimal_model$lambda,
               alpha = optimal_model$alpha)

  fs <- extract_non_zero_coefs(fit)

  return(list(model = fit, features = fs, para = optimal_model))
}

#' Cross-Validation for Elastic Net Model in Methylation Studies
#'
#' Conducts k-fold cross-validation on SNP data against methylation levels using an elastic net model.
#' This function splits the data into k folds, trains the model on k-1 folds, and evaluates it on the
#' remaining fold, iterating through all folds. It selects the best model based on MSE, and outputs
#' the corresponding correlation, MSE, and SNP weight vector for the best model.
#'
#' @param SNPs A matrix or data frame representing SNP data.
#'             Dimensions: Number of Observations (rows) x Number of SNPs (columns).
#' @param y A numeric vector representing methylation levels corresponding to each SNP observation.
#'          Length should match the number of rows in SNPs.
#' @param n_fold Integer, specifying the number of folds for cross-validation.
#' @param omit_folds_with_na_r Logical, indicating whether to omit folds with NA in results.
#' @param record_runtime Logical, indicating whether to record the runtime of the function.
#'
#' @return A list containing the best model, its MSE, correlation, SNP weights (with names),
#'         and total runtime if record_runtime is TRUE.
#'
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' cv_results <- cv.pred(SNPs, y, n_fold = 5, omit_folds_with_na_r = FALSE, record_runtime = TRUE)
#' }

cv.pred <- function(SNPs, y, n_fold = 5, omit_folds_with_na_r = FALSE, record_runtime = TRUE) {
  set.seed(2018)

  if (nrow(SNPs) != length(y)) {
    stop("Number of observations is different")
  }

  # Tune alpha and lambda using the entire dataset and get the best model
  tuned_params <- glmnet_tune_alpha(SNPs, y, n_fold = 5)
  best_model <- tuned_params$model
  best_alpha <- tuned_params$para$alpha
  best_lambda <- tuned_params$para$lambda

  # Set up for cross-validation
  fold_id <- sample(rep(1:n_fold, length.out = length(y)))
  cv <- matrix(NA, nrow = n_fold, ncol = 2)
  colnames(cv) <- c("cor", "mse")

  if (record_runtime) {
    start_time <- Sys.time()
  }

  # Outer k-fold cross-validation
  for (fold in 1:n_fold) {
    testIndices <- which(fold_id == fold, arr.ind = TRUE)
    X_test <- SNPs[testIndices,]
    yhat_test <- y[testIndices]

    predictions <- predict(best_model, X_test)
    cv[fold, 1] <- cor(predictions, yhat_test)
    cv[fold, 2] <- mean((predictions - yhat_test) ^ 2)
  }

  if(omit_folds_with_na_r == TRUE) {
    cv <- na.omit(cv)
  }

  cvm <- apply(cv, 2, mean)
  snp_weights <- coef(best_model)

  if (record_runtime) {
    end_time <- Sys.time()
    total_runtime <- as.numeric(end_time - start_time)
    cvm <- c(cvm, total_runtime = total_runtime)
  }

  return(list(
    cor = cvm[1],
    mse = cvm[2],
    alpha = best_alpha,
    lambda = best_lambda,
    runtime = total_runtime,
    snp_weights = snp_weights,
    model = best_model
  ))
}


cv.pred.3nest <- function(SNPs, y, n_fold = 5, omit_folds_with_na_r = FALSE, record_runtime = TRUE) {
  set.seed(2018)
  if (record_runtime) {
    start_time <- Sys.time()
  }

  if (nrow(SNPs) != length(y)) {
    stop("Number of observations is different")
  }

  #cat("Initial NA check - SNPs:", sum(is.na(SNPs)), "\n")
  if (sum(is.na(SNPs)) > 0) {
    SNPs <- impute_missing_values(SNPs)
  }
  #cat("Post imputation - SNPs:", sum(is.na(SNPs)), "\n")

  fold_id <- sample(rep(1:n_fold, length.out = length(y)))
  cv <- matrix(NA, nrow = n_fold, ncol = 3)
  colnames(cv) <- c("cor", "mse", "alpha")
  models <- list()

  for (fold in 1:n_fold) {
    testIndices <- which(fold_id == fold, arr.ind = TRUE)
    X_train <- SNPs[-testIndices,]
    yhat_train <- y[-testIndices]
    X_test <- SNPs[testIndices,]
    yhat_test <- y[testIndices]

    #print(paste0("Outer fold: ", fold))
    fit <- glmnet_tune_alpha(X_train, yhat_train, n_fold = 5)
    #print(fit$para)
    models[[fold]] <- fit$model
    predictions <- predict(fit$model, X_test)

    cv[fold, 1] <- cor(predictions, yhat_test)
    cv[fold, 2] <- mean((predictions - yhat_test) ^ 2)
    cv[fold, 3] <- fit$para$alpha
  }

  if(omit_folds_with_na_r == TRUE) {
    cv <- na.omit(cv)
  }
  #browser() # Why is alpha NULL?
  # Select the best model based on performance metrics
  best_model_idx <- which.min(cv[, "mse"])
  best_model <- models[[best_model_idx]]
  best_mse <- cv[best_model_idx, "mse"]
  best_cor <- cv[best_model_idx, "cor"]
  best_alpha <- cv[best_model_idx, "alpha"]
  best_lambda <- best_model$lambda

  # Extract and store SNP weights (coefficients) from the best model
  snp_weights <- coef(best_model)

  if (record_runtime) {
    end_time <- Sys.time()
    total_runtime <- as.numeric(end_time - start_time)
  }

  return(list(
    cor = best_cor,
    mse = best_mse,
    alpha = best_alpha,
    lambda = best_lambda,
    runtime = total_runtime,
    snp_weights = snp_weights,
    model = best_model
  ))
}
