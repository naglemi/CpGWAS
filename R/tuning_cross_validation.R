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

#' Remove Observations with Missing Response Values
#'
#' Removes any observations in the training set where the response variable (methylation levels)
#' is missing. This function ensures that both predictor (SNPs) and response datasets are aligned,
#' which is crucial in subsequent analyses of methylome-wide association studies.
#'
#' @param y Numeric vector representing the response variable (methylation levels).
#' @param X Matrix or data frame of predictors (SNP data).
#'          Dimensions: Observations (rows) x SNPs (columns).
#' @param verbose Logical indicating whether to print warnings about removed NAs.
#'
#' @return A list containing two elements: 'y' and 'X', both filtered to exclude
#'         observations with missing response values.
#'
#' @examples
#' \dontrun{
#' filtered_data <- rm_missing_observations(y, X, verbose = TRUE)
#' }
rm_missing_observations <- function(y, X, verbose) {
  if (!is.numeric(y)) {
    stop("y should be a numeric vector.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X should be a matrix or data frame.")
  }
  n_na <- sum(is.na(y))
  if(n_na > 0) {
    if(verbose) {
      warning(paste("There are", n_na, "NA values in y, which will be removed"))
    }
    idx <- !is.na(y)
    y <- y[idx]
    X <- X[idx, ]
  }
  return(list(y = y, X = X))
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
#' features <- extract_non_zero_coefs(fitted_model)
#' }
extract_non_zero_coefs <- function(fitted_model) {
  if (!inherits(fitted_model, "glmnet")) {
    stop("fitted_model must be an object of class 'glmnet'.")
  }
  # extract feature names and effects
  coef <- coef(fitted_model)[-1, ]
  if(sum(coef == 0) >= 1){
    coef <- coef[which(coef != 0)]
  }

  return(coef)
}

#' Hyperparameter Tuning for Elastic Net Model
#'
#' This function performs hyperparameter tuning for an elastic net model using glmnet.
#' It iterates over a range of alpha values to find the optimal lambda using cross-validation.
#' The function adapts its parallelization strategy based on the cores_per_alpha parameter
#' and the number of cores specified.
#'
#' @param X A matrix or data frame of predictors.
#'          Dimensions: Number of Observations (rows) x Number of Predictors (columns).
#' @param y A numeric vector representing the response variable.
#'          Length should match the number of rows in X.
#' @param n_fold Integer, specifying the number of folds for cross-validation.
#' @param verbose Logical, indicating whether to print detailed tuning results.
#' @param lambda_choice Character, specifying the method for choosing lambda
#'                      ('min' for minimum MSE, '1se' for one-standard-error rule).
#' @param alphas Numeric vector, specifying the alpha values to iterate over in tuning.
#' @param cores_per_alpha Can be "all" or 1. "all" uses all cores for parallel processing
#'                        within cv.glmnet (default setting). 1 uses all cores for parallel
#'                        processing of alpha values.
#' @param num_cores Integer, specifying the number of cores to use. Defaults to all
#'                  available cores if not provided.
#' @param allow_inefficient_parallelization Logical, allows inefficient parallelization
#'                                          when cores_per_alpha is 1 and there are
#'                                          more available cores than alpha values.
#'                                          Default is FALSE.
#'
#' @return A list containing the fitted model, non-zero coefficients, optimal parameters,
#'         and correlation on test data.
#'
#' @importFrom future plan availableCores
#' @importFrom future.apply future_lapply
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' \dontrun{
#' # Default example using all cores within cv.glmnet
#' fit_default <- glmnet_tune_alpha(X, y, n_fold = 5, verbose = TRUE,
#'                                  lambda_choice = "1se", alphas = seq(0, 1, 0.1))
#'
#' # Example with parallel processing of alpha values using all available cores
#' fit_parallel <- glmnet_tune_alpha(X, y, n_fold = 5, verbose = TRUE,
#'                                   lambda_choice = "1se", alphas = seq(0, 1, 0.1),
#'                                   cores_per_alpha = 1)
#' }
glmnet_tune_alpha <- function(X, y, n_fold, verbose, lambda_choice, alphas,
                              cores_per_alpha, num_cores,
                              allow_inefficient_parallelization, ...) {

  #set.seed(2023)

  if (cores_per_alpha == "all") {
    plan("sequential")
    registerDoParallel(cores = num_cores)
    if (verbose) {
      message("Parallelization within cv.glmnet using all available cores (",
              num_cores, ").")
    }
  } else if (cores_per_alpha == 1) {
    if(length(alphas) < num_cores && !allow_inefficient_parallelization) {
      stop(paste("Parallelization scheme is inefficient: number of alphas is",
                 "less than available cores. Consider using cores_per_alpha =",
                 "'all' or setting allow_inefficient_parallelization to TRUE."))
    }
    plan("multisession", workers = num_cores)
    if (verbose) {
      message("Parallel processing of alpha values using all available cores (",
              num_cores, ").")
    }
  } else {
    stop("Invalid value for cores_per_alpha. Only 'all' and 1 are allowed.")
  }

  if (nrow(X) != length(y)) {
    stop("Number of samples is different from length of response variable")
  }

  # remove missing data
  idx <- !is.na(y)
  y <- y[idx]
  X <- X[idx,]

  # define fold membership of each sample for k-fold cross-validation
  fold_id <- sample(rep(1:n_fold, length.out = length(y)))

  tuning_results <- future_lapply(alphas, function(alpha) {
    cv <- cv.glmnet(
      X,
      y,
      foldid = fold_id,
      type.measure = "mse",
      parallel = (cores_per_alpha == "all"),
      alpha = alpha
    )

    if(lambda_choice == "1se") {
      lambda_selected <- cv$lambda.1se
    } else if(lambda_choice == "min") {
      lambda_selected <- cv$lambda.min
    } else {
      stop("Invalid lambda_choice: choose either '1se' or 'min'")
    }

    data.frame(
      cvm = cv$cvm[cv$lambda == lambda_selected],
      lambda = lambda_selected,
      alpha = alpha
    )
  })
  
  # loop version so we can use debugger (commented out when not debugging)
  
  # Initialize an empty list to store the results
  # tuning_results <- list()
  # 
  # # Loop through each alpha
  # for(i in seq_along(alphas)) {
  #   alpha <- alphas[i]
  #   cv <- cv.glmnet(
  #     X,
  #     y,
  #     foldid = fold_id,
  #     type.measure = "mse",
  #     parallel = FALSE,
  #     alpha = alpha
  #   )
  #   
  #   if(lambda_choice == "1se") {
  #     lambda_selected <- cv$lambda.1se
  #   } else if(lambda_choice == "min") {
  #     lambda_selected <- cv$lambda.min
  #   } else {
  #     stop("Invalid lambda_choice: choose either '1se' or 'min'")
  #   }
  #   
  #   # Store the result in the list
  #   tuning_results[[i]] <- data.frame(
  #     cvm = cv$cvm[cv$lambda == lambda_selected],
  #     lambda = lambda_selected,
  #     alpha = alpha
  #   )
  # }
  
  # end loop version to be commented out when we're not using debugger

  # Combine the results
  tuning_results <- do.call(rbind, tuning_results)

  cv.opt <- tuning_results[which.min(tuning_results$cvm),]

  # fit final model
  fitted_model <- glmnet(
    X,
    y,
    lambda = cv.opt$lambda,
    alpha = cv.opt$alpha
  )

  features <- extract_non_zero_coefs(fitted_model)

  pred <- predict(fitted_model, X)

  if(length(levels(factor(pred))) > 1){
    r <- cor(pred, y)
  } else {
    r <- NA
  }
  
  # Debug: Check tuning results
  if (verbose) {
    cat("Tuning results - Lambda:", cv.opt$lambda,
        "Alpha:", cv.opt$alpha, "\n")
  }

  return(list(model = fitted_model, features = features,
              para = cv.opt, cor = r))
}

cv_eval <- function(X, y, n_fold, cv_eval_mode, verbose, alphas, cores_per_alpha, num_cores,
                    allow_inefficient_parallelization, omit_folds_with_na_r,
                    best_alpha = NULL, best_lambda = NULL, ...) {
  #set.seed(2018)

  if (nrow(X) != length(y)) {
    stop("Number of observations is different")
  }

  #cat("Initial NA check - X:", sum(is.na(X)), "\n")
  if (sum(is.na(X)) > 0) {
    X <- impute_missing_values(X)
  }
  #cat("Post imputation - X:", sum(is.na(X)), "\n")

  fold_id <- sample(rep(1:n_fold, length.out = length(y)))

  if(cv_eval_mode == "dynamic"){
    evaluation_results <- cv_eval_dynamic(X = X, y = y, n_fold = n_fold, fold_id,
                                          verbose = verbose, alphas = alphas,
                                          cores_per_alpha = cores_per_alpha,
                                          num_cores = num_cores,
                                          omit_folds_with_na_r = omit_folds_with_na_r,
                                          allow_inefficient_parallelization = allow_inefficient_parallelization, ...)
  } else if(cv_eval_mode == "static"){
    evaluation_results <- cv_eval_static(X = X, y = y, n_fold = n_fold, fold_id,
                                         best_alpha = best_alpha, best_lambda = best_lambda,
                                         omit_folds_with_na_r = omit_folds_with_na_r,
                                         verbose = verbose, ...)
  } else {
    stop("Invalid cv_eval_mode: choose either 'dynamic' or 'static'")
  }
}

cv_eval_static <- function(X, y, n_fold, fold_id, cores_per_alpha,
                           num_cores, allow_inefficient_parallelization,
                           omit_folds_with_na_r,
                           best_alpha = NULL, best_lambda = NULL, ...) {

  if(is.null(best_lambda) || is.null(best_alpha)){
    stop(paste("Error: best_lambda and best_alpha must be specified for k-fold",
               "validation to evaluate model in `static` mode. If you want to",
               "evaluate performance without specifying best_lambda and best_alpha,",
               "use cv_eval_mode = `dynamic` instead."))
  }

  cv <- matrix(NA, nrow = n_fold, ncol = 4,
               dimnames = list(NULL, c("cor", "mse", "alpha", "lambda")))

  cv[, 3] <- best_alpha
  cv[, 4] <- best_lambda

  # Let's not doing this part in parallel because it will run in less than
  #  a second on one core in this mode. Here, parameters are known and we only
  #  need to fit one model per fold.

  for (fold in 1:n_fold) {
    testIndices <- which(fold_id == fold, arr.ind = TRUE)
    X_train <- X[-testIndices,]
    y_train <- y[-testIndices]
    X_test <- X[testIndices,]
    y_test <- y[testIndices]

    fit <- glmnet(x = X,
                  y = y,
                  alpha = best_alpha,
                  lambda = best_lambda)

    pred <- predict(fit,
                    X_test)
    
    if(length(levels(factor(pred))) > 1){
      cv[fold, 1] <- cor(pred, y_test)
    } else {
      #recover()
      cv[fold, 1] <- NA
    }

    cv[fold, 2] <- mean((pred - y_test)^2)

    # Make sure we get same result with calculate_correlation(), and make sure
    # fit$para$alpha and lambda match best.
  }

  if(omit_folds_with_na_r == TRUE) {
    cv <- na.omit(cv)
  }

  cvm <- colMeans(cv[, 1:2])

  return(list(eval_cv = cv,
              eval_cvm = cvm))
}


cv_eval_dynamic <- function(X, y, n_fold, fold_id, verbose, alphas, cores_per_alpha,
                            num_cores, allow_inefficient_parallelization,
                            omit_folds_with_na_r, ...) {

  cv <- matrix(NA, nrow = n_fold, ncol = 4,
               dimnames = list(NULL, c("cor", "mse", "alpha", "lambda")))

  #set.seed(2018)

  for (fold in 1:n_fold) {
    testIndices <- which(fold_id == fold, arr.ind = TRUE)
    X_train <- X[-testIndices,]
    y_train <- y[-testIndices]
    X_test <- X[testIndices,]
    y_test <- y[testIndices]

    fit <- glmnet_tune_alpha(X = X_train, y = y_train, n_fold = n_fold,
                             verbose = verbose, alphas = alphas,
                             cores_per_alpha = cores_per_alpha,
                             num_cores = num_cores,
                             allow_inefficient_parallelization = allow_inefficient_parallelization, ...)

    pred <- predict(fit$model,
                    X_test)

    if(length(levels(factor(pred))) > 1){
      cv[fold, 1] <- cor(pred, y_test)
    } else {
      cv[fold, 1] <- NA
    }

    cv[fold, 2] <- mean((pred - y_test)^2)

    cv[fold, 3] <- fit$para$alpha

    cv[fold, 4] <- fit$para$lambda
  }

  if(omit_folds_with_na_r == TRUE) {
    cv <- na.omit(cv)
  }

  cvm <- colMeans(cv[, 1:2])

  return(list(eval_cv = cv,
              eval_cvm = cvm))
}
