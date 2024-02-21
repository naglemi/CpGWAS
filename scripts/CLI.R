#!/usr/bin/env Rscript
#profvis({

#options(error = recover)

start_time <- Sys.time()  # Start time capture

library(CpGWAS)
library(optparse)

# Pt. 1: Load libraries and accept, parse user arguments -------------------------------------

# Check if running in RStudio or via command line
if(Sys.getenv("RSTUDIO") != "1") {
  # Define command line options
  option_list <- list(
    make_option(c("--outdir"), type = "character", default = "./output/",
                help = "Output directory, default is './output/'"),
    make_option(c("--chunk1"), type = "integer", default = 1,
                help = "Starting methylation site index for processing, default is 1"),
    make_option(c("--chunk2"), type = "integer", default = NA,
                help = "Ending methylation site index for processing, default is NA (all sites)"),
    make_option(c("-s", "--snp_data_path"), type = "character", default = NULL,
                help = "Path to SNP data (required)"),
    make_option(c("-m", "--methylation_data_path"), type = "character", default = NULL,
                help = "Path to methylation data (required)"),
    make_option(c("-v", "--verbose"), type = "logical", default = FALSE,
                help = "Logical, indicating whether to print detailed tuning results, default is FALSE"),
    make_option(c("-l", "--lambda_choice"), type = "character", default = "1se",
                help = "Method for choosing lambda ('min' for minimum MSE, '1se' for one-standard-error rule), default is '1se'"),
    make_option(c("-a", "--alphas"), type = "character", default = "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1",
                help = "Comma-separated list of alpha values to iterate over in tuning, default is '0,0.1,...,1'"),
    make_option(c("-c", "--cores_per_alpha"), type = "character", default = "all",
                help = "Can be 'all' or '1'. 'all' uses all cores for parallel processing within cv.glmnet (default). '1' uses all cores for parallel processing of alpha values"),
    make_option(c("-n", "--num_cores"), type = "integer", default = future::availableCores(),
                help = "Number of cores to use, defaults to all available cores"),
    make_option(c("--allow_inefficient_parallelization"), type = "logical", default = FALSE,
                help = "Logical, allows inefficient parallelization when cores_per_alpha is 1 and there are more available cores than alpha values. Default is FALSE"),
    make_option(c("-f", "--n_fold"), type = "integer", default = 5,
                help = "Integer, specifying the number of folds for cross-validation, default is 5"),
    make_option(c("-w", "--window_sizes"), type = "character", default = "1000,2000,5000,10000,20000,50000,100000,200000,500000",
                help = "Comma-separated list of window sizes, default is '1000,2000,5000,10000,20000,50000,100000,200000,500000'"),
    make_option(c("-t", "--tag"), type = "character", default = format(Sys.time(), "%Y%m%d-%H%M%S"),
                help = "Tag to append to scaffoldIdentifier, default is current date-time stamp"),
    make_option(c("-e", "--save_evaluation_results_each_fold"), type = "logical", default = FALSE,
                help = "Logical, indicating whether to save MSE, R^2, alpha, lambda for each evaluation fold, default is FALSE"),
    make_option(c("-g", "--save_glmnet_object"), type = "logical", default = FALSE,
                help = "Logical, indicating whether to save the glmnet object, default is FALSE"),
    make_option(c("--cv_eval_mode"), type = "character", default = "dynamic",
                help = paste("Character, indicating whether to use 'dynamic' or 'static' cross-validation for",
                "model evaluation; allow lambda and alpha to vary (dynamic) or keep static after optimization, default is 'dynamic'")),
    make_option(c("--omit_folds_with_na_r"), type = "logical", default = TRUE,
                help = "Logical, indicating whether to omit folds with NA R values, default is TRUE"),
    make_option(c("-r", "--methInput_rds_path"), type = "character", default = NULL,
                help = "Path to an RDS file containing a pre-existing MethylationInput object. If provided, this object will be loaded instead of creating a new one.")
  )

  # Parse options
  args <- parse_args(OptionParser(option_list = option_list))
  args$alphas <- as.numeric(unlist(strsplit(args$alphas, ",")))
  args$window_sizes <- as.numeric(unlist(strsplit(args$window_sizes, ",")))

  if(args$verbose) {
    print("Args from command line: ")
  }

} else {
  args <- list(
    outdir = "./output/",
    chunk1 = 7751,
    chunk2 = 8000,
    snp_data_path = "/Users/mnagle6/data/libd_chr1.pgen",
    methylation_data_path = "/Users/mnagle6/data/chr1_AA.rda",
    verbose = TRUE,
    lambda_choice = "1se",
    alphas = c(0.5),#seq(0, 1, 0.25),
    cores_per_alpha = "all",
    num_cores = "all", #future::availableCores(),
    allow_inefficient_parallelization = FALSE,
    n_fold = 5,
    window_sizes = c(2000,4000,6000,8000),# 5000, 10000, 20000, 50000, 100000, 200000, 500000),
    tag = format(Sys.time(), "%Y%m%d-%H%M%S"),
    save_evaluation_results_each_fold = FALSE,
    save_glmnet_object = FALSE,
    cv_eval_mode = "dynamic",
    omit_folds_with_na_r = TRUE,
    methInput_rds_path = "~/data/chr1_AA_methylation_10k_samples.rds"
  )

  if(args$verbose) {
    print("Args from RStudio: ")
  }
}

#saveRDS(args, file = file.path(args$outdir, paste0(args$tag, "-args.rds")))
#args <- readRDS("output/libd_chr1-chr1_AA-static-1core-20240129-123107-args.rds")
if(args$num_cores == "all"){
  args$num_cores <- future::availableCores()
}

if(args$verbose) {
  print(args)
}

if(!dir.exists(args$outdir)) {
  dir.create(args$outdir)
}

# Check required arguments
if (is.null(args$snp_data_path) || is.null(args$methylation_data_path)) {
  stop("Paths to both SNP and methylation data are required.")
}

load(args$methylation_data_path)

# Pt. 2: Initialize (or load) MethylationInput object -------------------------------

if (!is.null(args$methInput_rds_path) && file.exists(args$methInput_rds_path)) {
  if(args$verbose) {
    message("Loading MethylationInput object from RDS file: ", args$methInput_rds_path)
  }
  methInput <- reinitializeMethylationInput(rds_path = args$methInput_rds_path,
                                            snp_data_path = args$snp_data_path,
                                            start_site = args$chunk1,
                                            end_site = args$chunk2,
                                            no_cores = args$num_cores)
} else {
  if(args$verbose) {
    message("Creating new MethylationInput object")
  }
  methInput <- new("MethylationInput",
                   BSseq_obj = BSobj2,
                   snp_data_path = args$snp_data_path,
                   start_site = args$chunk1,
                   end_site = args$chunk2,
                   no_cores = args$num_cores)
  BSobj2 <- means <- sds <- NULL
  #saveRDS(methInput, "~/data/chr1_AA_methylation_10k_samples.rds")
}

# Pt. 3: Main loop to process SNP data for each methylation site ----------

scaffoldIdentifier <- paste0(tools::file_path_sans_ext(basename(args$snp_data_path)),
                             "-",
                             tools::file_path_sans_ext(basename(args$methylation_data_path)),
                             "-",
                             args$tag)

scaffold_models <- fit_MWAS_models(
  BSobj = BSobj2,
  methInput = methInput,
  window_sizes = args$window_sizes,
  chunk1 = 1,
  chunk2 = length(methInput@methylations_positions),
  n_fold = args$n_fold,
  scaffoldIdentifier = scaffoldIdentifier,
  outdir = args$outdir,
  verbose = args$verbose,
  lambda_choice = args$lambda_choice,
  alphas = args$alphas,
  cores_per_alpha = args$cores_per_alpha,
  num_cores = args$num_cores,
  allow_inefficient_parallelization = args$allow_inefficient_parallelization,
  save_evaluation_results_each_fold = args$save_evaluation_results_each_fold,
  save_glmnet_object = args$save_glmnet_object,
  cv_eval_mode = args$cv_eval_mode,
  omit_folds_with_na_r = args$omit_folds_with_na_r
)
#df <- as.data.frame(scaffold_models)

pgenlibr::ClosePgen(methInput@pgen)
pgenlibr::ClosePvar(methInput@pvar_pointer)

end_time <- Sys.time()  # End time capture

# Calculate total runtime
total_runtime <- end_time - start_time

# Convert 'total_runtime' to numeric seconds
total_runtime_seconds <- as.numeric(total_runtime, units = "secs")

# Calculate hours, minutes, and seconds
hours <- total_runtime_seconds %/% 3600
minutes <- (total_runtime_seconds %% 3600) %/% 60
seconds <- total_runtime_seconds %% 60
# Convert to integers
hours_int <- as.integer(hours)
minutes_int <- as.integer(minutes)
seconds_int <- as.integer(seconds)

# Format the runtime
runtime_formatted <- sprintf("%02d:%02d:%02d", hours_int, minutes_int, seconds_int)

# Report the runtime
cat(sprintf("Processed chunks %d through %d in %d hours, %d minutes, %d seconds.\n",
            args$chunk1, args$chunk2, hours_int, minutes_int, seconds_int))


# Get system information
logical_cores <- future::availableCores()
physical_cores_cmd <- ifelse(Sys.info()["sysname"] == "Darwin",
                             "sysctl -n hw.physicalcpu",
                             "lscpu | grep 'Core(s) per socket:' | awk '{print $4}'")
physical_cores <- as.integer(system(physical_cores_cmd, intern = TRUE))
if (Sys.info()["sysname"] != "Darwin") {
  physical_cores <- physical_cores * as.integer(system(
    "lscpu | grep 'Socket(s):' | awk '{print $2}'",
    intern = TRUE))
}
hyperthreading_enabled <- ifelse(logical_cores > physical_cores, "Yes", "No")
type_CPU_cmd <- ifelse(Sys.info()["sysname"] == "Darwin",
                       "sysctl -n machdep.cpu.brand_string",
                       "lscpu | grep 'Model name:' | awk -F ':' '{print $2}'")
type_CPU <- system(type_CPU_cmd, intern = TRUE)
amount_RAM_cmd <- ifelse(Sys.info()["sysname"] == "Darwin",
                         "sysctl -n hw.memsize",
                         "grep MemTotal /proc/meminfo | awk '{print $2}'")
amount_RAM <- as.numeric(system(amount_RAM_cmd, intern = TRUE)) / 1024^2
if (Sys.info()["sysname"] != "Darwin") {
  amount_RAM <- amount_RAM / 1024  # Adjust for Linux units
}

# Prepare system information for output
system_info <- data.frame(
  time_started = format(start_time, "%Y-%m-%d %H:%M:%S"),
  time_finished = format(end_time, "%Y-%m-%d %H:%M:%S"),
  runtime = runtime_formatted,
  type_CPU = type_CPU,
  amount_RAM = amount_RAM,
  number_cores = logical_cores,
  physical_cores = physical_cores,
  hyperthreading_enabled = hyperthreading_enabled,
  scaffold_ID = scaffoldIdentifier,
  stringsAsFactors = FALSE
)

# Combine arguments and system info
args_df <- as.data.frame(t(unlist(args)), stringsAsFactors = FALSE)
# Transform args_df to a long format
args_long <- stack(args_df)
names(args_long) <- c("Value", "Parameter")

# Transform system_info to a long format
system_info_long <- stack(system_info)
names(system_info_long) <- c("Value", "Parameter")

# Combine args_long and system_info_long
summary_df <- rbind(args_long, system_info_long)

# Switch the order of the columns
summary_df <- summary_df[, c("Parameter", "Value")]

# Define the output file path and save to CSV
summary_file_path <- file.path(args$outdir, paste0(scaffoldIdentifier, "-summary.csv"))
write.csv(summary_df, summary_file_path, row.names = FALSE, quote = TRUE)
#})