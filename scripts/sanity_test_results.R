# Compare outputs from before and after memory optimization to make sure we didn't break anything

before1 <- readRDS("output/libd_chr1-chr1_AA-20240130-121410.rds")
before2 <- readRDS("output/libd_chr1-chr1_AA-20240130-153613.rds")
before3 <- readRDS("output/libd_chr1-chr1_AA-20240130-153840.rds")
after <- readRDS("output/libd_chr1-chr1_AA-20240130-164421.rds")

compareMethylationScaffs <- function(scaff1, scaff2) {
  if (!inherits(scaff1, "MethylationScaff") || !inherits(scaff2, "MethylationScaff")) {
    stop("Both inputs must be MethylationScaff objects.")
  }
  
  # Compare scaffoldIdentifier (expecting differences due to timestamps)
  cat("Comparing scaffoldIdentifier...\n")
  cat("scaff1: ", scaff1@scaffoldIdentifier, "\n")
  cat("scaff2: ", scaff2@scaffoldIdentifier, "\n")
  
  # Compare the length of models
  cat("Comparing number of models...\n")
  if (length(scaff1@models) == length(scaff2@models)) {
    cat("Both have the same number of models: ", length(scaff1@models), "\n")
  } else {
    cat("Different number of models: scaff1 - ", length(scaff1@models), ", scaff2 - ", length(scaff2@models), "\n")
  }
  
  # Compare each model
  for (i in seq(along = scaff1@models)) {
    cat("Comparing model ", i, "\n")
    
    model1 <- scaff1@models[[i]]
    model2 <- scaff2@models[[i]]
    
    # Compare methylationPosition, windowSize, n_SNPs, alpha, lambda
    cat("methylationPosition: ", model1@methylationPosition == model2@methylationPosition, "\n")
    cat("windowSize: ", model1@windowSize == model2@windowSize, "\n")
    cat("n_SNPs: ", model1@n_SNPs == model2@n_SNPs, "\n")
    cat("alpha: ", model1@alpha == model2@alpha, "\n")
    cat("lambda: ", model1@lambda == model2@lambda, "\n")
    
    # Compare snpWeights - expecting slight differences due to stochasticity
    cat("snpWeights comparison...\n")
    diff_weights <- sum(abs(model1@snpWeights - model2@snpWeights))
    cat("Sum of absolute differences in snpWeights: ", diff_weights, "\n")
    
    # Compare intercept
    cat("intercept: ", all.equal(model1@intercept, model2@intercept), "\n")
    
    # Compare evaluation_results
    cat("evaluation_results comparison...\n")
    eval_diffs <- sum(abs(model1@evaluation_results - model2@evaluation_results))
    cat("Sum of absolute differences in evaluation_results: ", eval_diffs, "\n")
  }
}

printLargeObjects <- function(minSizeMB = 1) {
  # Get all object names in the global environment
  all_objects <- ls(envir = .GlobalEnv)
  
  # Loop through each object and check its size
  for (obj_name in all_objects) {
    object_size <- object.size(get(obj_name, envir = .GlobalEnv)) / (1024^2) # Convert to MB
    if (object_size > minSizeMB) {
      cat(obj_name, ": ", round(object_size, 2), "MB\n")
    }
  }
}

printLargeObjects(0.1)

compareMethylationScaffs(before1, before2)
