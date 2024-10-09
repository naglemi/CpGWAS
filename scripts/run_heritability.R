# run_heritability.R

# Load necessary libraries
library(data.table)
library(stringr)
library(bsseq)
library(tools)

# Define a helper function to check file existence and print messages
check_and_report <- function(file_path, description) {
  if (file.exists(file_path)) {
    full_path <- normalizePath(file_path)
    current_wd <- getwd()
    # Since all paths are absolute, relative_path is same as full_path
    relative_path <- full_path
    message(paste0(description, " saved to full path ", full_path,
                   " from current wd ", current_wd,
                   " using relative path ", relative_path))
  } else {
    message(paste0("Failed to save ", description, " to ", file_path))
  }
}

# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  args_list <- list()

  # Iterate over arguments and assign to list
  for (arg in args) {
    if (grepl("^--", arg)) {
      key_value <- strsplit(sub("^--", "", arg), "=")[[1]]
      key <- key_value[1]
      value <- ifelse(length(key_value) > 1, key_value[2], TRUE)
      args_list[[key]] <- value
    }
  }

  return(args_list)
}

# Parse arguments
opt <- parse_args()

# Assign variables from arguments
Chr <- as.integer(opt$Chr)
SNP_data <- opt$SNP_data
methylation_data <- opt$methylation_data
cov_file <- opt$cov_file
modified_methylation_data <- opt$modified_methylation_data
wind <- as.numeric(opt$wind)
chunk_start <- as.integer(opt$chunk_start)
chunk_end <- as.integer(opt$chunk_end)
gwas <- opt$gwas
gcta <- opt$gcta
tag <- opt$tag

# Extract the base name of the input file without extension
input_file_base <- file_path_sans_ext(basename(modified_methylation_data))

# Define the parent directory for all GCTA outputs
parent_outdir <- file.path("/dcs04/lieber/statsgen/mnagle/mwas/CpGWAS/scripts", "all_gcta_outputs")

# Create the parent directory if it doesn't exist
if (!dir.exists(parent_outdir)) {
  dir.create(parent_outdir, recursive = TRUE)
  message(paste0("Parent output directory '", parent_outdir, "' created."))
}

# Modify the output directory to be inside the parent directory
outdir <- file.path(parent_outdir, paste0("gcta_output_", input_file_base))

# Create output directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  message(paste0("Output directory '", outdir, "' created inside '", parent_outdir, "'."))
}

# Set error handling to save dumps in the unique 'outdir'
options(error = function() {
  dump_file <- file.path(outdir, paste0("error_dump_", Sys.getenv("SLURM_JOB_ID"), ".rda"))
  dump.frames(dump_file, to.file = TRUE)
  message(paste0("Error dump saved to ", dump_file))
  q(status = 1)
})

# Verify that all input files exist
input_files <- list(
  "SNP data" = SNP_data,
  "Methylation data" = methylation_data,
  "Covariate file" = cov_file,
  "Modified methylation data" = modified_methylation_data
)

for (desc in names(input_files)) {
  file_path <- input_files[[desc]]
  if (!file.exists(file_path)) {
    stop(paste0(desc, " file not found at path: ", file_path))
  } else {
    check_and_report(file_path, desc)
  }
}

# Load necessary data
p <- readRDS(modified_methylation_data)
check_and_report(modified_methylation_data, "Modified methylation data")

CpG_positions <- p@methylations_positions
p <- p@methylations

ind <- fread(cov_file)$ID  # Adjust if necessary
check_and_report(cov_file, "Covariate file")

# Create 'id' file in 'outdir' based on 'ind'
id_file <- file.path(outdir, "id")
write.table(data.frame(FID = ind, IID = ind), file = id_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
check_and_report(id_file, "ID file")

# Initialize result list and index to avoid copying data when appending
res <- list()
res_index <- 1  # Index to track the position in the result list

# Loop over CpG positions using local indices
for (idx in seq_along(CpG_positions)) {
  i <- idx  # Local index within CpG_positions
  CpG_index <- chunk_start + idx - 1  # Original global index
  chr_num <- Chr

  for (w in seq_along(wind)) {
    # PLINK subset
    CpG_pos <- CpG_positions[idx]  # Use local index
    p1 <- ifelse(CpG_pos - wind[w] > 0, CpG_pos - wind[w], 0)
    p2 <- CpG_pos + wind[w]
    gwas_prefix <- paste0(gwas, "libd_chr", chr_num)

    # Define paths for intermediate files within the unique 'outdir'
    plink_out_prefix <- file.path(outdir, paste0("temp_", i, "_", w))

    command_plink <- paste("plink2 --pfile", shQuote(gwas_prefix), "--silent --keep", shQuote(id_file),
                           "--chr", chr_num,
                           "--from-bp", p1, "--to-bp", p2,
                           "--snps-only 'just-acgt' --make-bed --threads 1 --out", shQuote(plink_out_prefix))
    system(command_plink)

    # Check for PLINK output files
    bim_file <- paste0(plink_out_prefix, ".bim")
    bed_file <- paste0(plink_out_prefix, ".bed")
    fam_file <- paste0(plink_out_prefix, ".fam")

    check_and_report(bim_file, "PLINK BIM file")
    check_and_report(bed_file, "PLINK BED file")
    check_and_report(fam_file, "PLINK FAM file")

    if (!file.exists(bim_file)) {
      message("temp.bim not found. Skipping this window for CpG index ", i, ".")
      next
    }

    # Phenotype file
    pheno <- data.frame(FID = ind, IID = ind, Phenotype = p[, i])
    pheno_path <- file.path(outdir, paste0("pheno_", i, "_", w, ".txt"))
    write.table(pheno, pheno_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    check_and_report(pheno_path, "Phenotype file")

    # GRM
    command_grm <- paste(shQuote(gcta), "--bfile", shQuote(plink_out_prefix),
                         "--make-grm-bin --thread-num 1 --out", shQuote(plink_out_prefix))
    system(command_grm, intern = FALSE)

    # Check for GRM output files
    grm_bin <- paste0(plink_out_prefix, ".grm.bin")
    grm_n <- paste0(plink_out_prefix, ".grm.N.bin")
    grm_id <- paste0(plink_out_prefix, ".grm.id")

    check_and_report(grm_bin, "GRM BIN file")
    check_and_report(grm_n, "GRM N BIN file")
    check_and_report(grm_id, "GRM ID file")

    qcovar_file <- sub("\\.csv$", ".qcovar", cov_file)
    covar_file <- sub("\\.csv$", ".covar", cov_file)

    # Heritability estimation
    command_reml <- paste(shQuote(gcta), "--reml --grm-bin", shQuote(plink_out_prefix),
                          "--pheno", shQuote(pheno_path),
                          "--mpheno 1 --qcovar", shQuote(qcovar_file),
                          "--covar", shQuote(covar_file),
                          "--thread-num 1 --out", shQuote(plink_out_prefix))
    system(command_reml, intern = FALSE)

    # Check for REML output file
    hsq_file <- paste0(plink_out_prefix, ".hsq")
    check_and_report(hsq_file, "REML HSQ file")

    if (!file.exists(hsq_file)) {
      message("temp.hsq not found after REML. Skipping CpG index ", i, ", window ", wind[w], ".")
      next
    }

    # Read and reshape temp.hsq
    temp_before_transpose <- fread(hsq_file, header = TRUE, fill = TRUE)

    # Extract relevant values using single square brackets
    temp2 <- list(
      V_G = temp_before_transpose[1, 2],
      SE_V_G = temp_before_transpose[1, 3],
      V_e = temp_before_transpose[2, 2],
      SE_V_e = temp_before_transpose[2, 3],
      Vp = temp_before_transpose[3, 2],
      SE_Vp = temp_before_transpose[3, 3],
      V_G_Vp = temp_before_transpose[4, 2],
      SE_V_G_Vp = temp_before_transpose[4, 3],
      logL = temp_before_transpose[5, 2],
      logL0 = temp_before_transpose[6, 2],
      LRT = temp_before_transpose[7, 2],
      df = temp_before_transpose[8, 2],
      Pval = temp_before_transpose[9, 2],
      n = temp_before_transpose[10, 2],
      site = paste0("chr", chr_num, "_", CpG_positions[i]),
      wind = wind[w]
    )

    if (is.null(temp2) || temp2$n <= 1) {
      message("Invalid result for CpG index ", i, ", window ", wind[w], ". Skipping.")
      next
    }

    # Append to results using direct assignment to avoid copying the entire list
    res[[res_index]] <- temp2
    res_index <- res_index + 1  # Increment the index

    # Remove temporary files specific to this job
    temp_files_pattern <- paste0("^temp_", i, "_", w)
    temp_files <- list.files(outdir, pattern = temp_files_pattern, full.names = TRUE)
    if (length(temp_files) > 0) {
      file.remove(temp_files)
      message(paste0("Temporary files removed: ", paste(basename(temp_files), collapse = ", ")))
    }

    # Remove phenotype file
    if (file.exists(pheno_path)) {
      file.remove(pheno_path)
      message(paste0("Phenotype file removed: ", basename(pheno_path)))
    }
  }
  cat("\n")
}

if (length(res) > 0) {
  res_dt <- rbindlist(res, fill = TRUE)
} else {
  res_dt <- data.table()
  warning("No valid results to save.")
}

# Define output file path inside 'outdir' and add '.csv' extension
output_file <- file.path(outdir, paste0(file_path_sans_ext(basename(modified_methylation_data)), "_heritability_results.csv"))

# Save results as CSV
fwrite(res_dt, file = output_file)
check_and_report(output_file, "Final heritability results CSV")

message(paste0("Heritability estimation completed. Results saved to ", output_file))
