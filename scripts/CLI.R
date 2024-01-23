library(CpGWAS)

#rm(list = ls())

record_runtime <- TRUE

# Pt. 1: Load libraries and accept, parse user arguments -------------------------------------

if(Sys.getenv("RSTUDIO") == 0){

  # Define command line options
  option_list <- list(
    make_option(c("-o", "--outdir"), type = "character", default = "./output/",
                help = "Output directory, default is './output/'"),
    make_option(c("-c1", "--chunk1"), type = "integer", default = 1,
                help = "Starting methylation site index for processing"),
    make_option(c("-c2", "--chunk2"), type = "integer", default = NA,
                help = "Ending methylation site index for processing"),
    make_option(c("-s", "--snp_data_path"), type = "character", default = NULL,
                help = "Path to SNP data (required)"),
    make_option(c("-m", "--methylation_data_path"), type = "character", default = NULL,
                help = "Path to methylation data (required)")
  )

  # Parse options
  args <- parse_args(OptionParser(option_list = option_list))
} else {
  args <- list(
    outdir = "./output/",
    chunk1 = 1000000,
    chunk2 = 1000100,
    snp_data_path = "/Users/michaelnagle/code/mwas/gwas/libd_chr1.pgen",
    methylation_data_path = "/Users/michaelnagle/code/mwas/pheno/dlpfc/out/chr1_AA.rda"
  )
}

if(!dir.exists(args$outdir)) {
  dir.create(args$outdir)
}

# Check required arguments
if (is.null(args$snp_data_path) || is.null(args$methylation_data_path)) {
  stop("Paths to both SNP and methylation data are required.")
}

load(args$methylation_data_path)

# Pt. 2: Initialize MethylationInput object -------------------------------

methInput <- new("MethylationInput",
                 BSseq_obj = BSobj2,
                 snp_data_path = args$snp_data_path,
                 args = args)

# Pt. 3: Main loop to process SNP data for each methylation site ----------

start_time <- Sys.time()  # Start time capture

window_sizes <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)

scaffoldIdentifier <- paste0(tools::file_path_sans_ext(basename(args$snp_data_path)),
                             "-",
                             tools::file_path_sans_ext(basename(args$methylation_data_path)))

scaffold_models <- build_prediction_model(
  BSobj = BSobj2,
  methInput = methInput,
  window_sizes = c(1000, 2000),
  chunk1 = 10^6,
  chunk2 = 10^6 + 10,
  n_fold = 5,
  cv_nesting = "double",
  scaffoldIdentifier = scaffoldIdentifier,
  outdir = "output/",
  record_runtime = TRUE
)

df <- as.data.frame(scaffold_models)

pgenlibr::ClosePgen(methInput@pgen)
pgenlibr::ClosePvar(methInput@pvar1)

end_time <- Sys.time()  # End time capture
total_runtime <- end_time - start_time
total_runtime_seconds <- as.numeric(total_runtime, units = "secs")
hours <- total_runtime_seconds %/% 3600
minutes <- (total_runtime_seconds %% 3600) %/% 60
seconds <- total_runtime_seconds %% 60

# Report the runtime
cat(sprintf("Processed chunks %d through %d in %d hours, %d minutes and %d seconds.\n",
            args$chunk1, args$chunk2, as.integer(hours), as.integer(minutes), as.integer(seconds)))
