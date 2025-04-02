#!/usr/bin/env Rscript
# fast_19_plot.R - Optimized version of 19_plot.ipynb and 19_plot_a2.ipynb
# This script efficiently recreates the plots in the 19-OUT_plots directory
# with proper handling of input data files from the OLD directory

# Libraries needed for visualization
library(data.table)
library(qqman)
library(stringr)
library(dplyr)

# Create output directory with proper naming convention
output_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/19_plots_fast"
if (!dir.exists(output_dir)) {
  if (!dir.create(output_dir, recursive = TRUE)) {
    stop("CRITICAL ERROR: Failed to create output directory: ", output_dir,
         "\nThis is a critical error - we cannot proceed without a valid output directory!")
  }
  message(paste("Created output directory:", output_dir))
} else {
  message(paste("Using existing output directory:", output_dir))
}

# Define which traits to process
traits <- c("scz", "bp", "mdd")
message("Starting fast_19_plot.R script...")

# Function to clean GWAS summary statistics column names - FAIL FAST if requirements not met
clean_and_standardize_colnames <- function(summary_stats) {
  # Check if the header is tab-delimited while the rest is space-delimited
  if (grepl("\t", colnames(summary_stats)[1])) {
    real_colnames <- str_split(colnames(summary_stats)[1], "\t")[[1]]
    colnames(summary_stats) <- real_colnames
  }
  
  # Standardize column names - FAIL FAST if required columns are missing
  required_cols <- c("chr", "pos", "p", "BETA")
  missing_cols <- setdiff(required_cols, tolower(colnames(summary_stats)))
  if (length(missing_cols) > 0) {
    stop("CRITICAL ERROR: Missing required columns: ",
         paste(missing_cols, collapse = ", "),
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Standardize column names
  colnames(summary_stats) <- gsub("chr", "CHR", colnames(summary_stats))
  colnames(summary_stats) <- gsub("#CHROM", "CHR", colnames(summary_stats))
  colnames(summary_stats) <- gsub("pos", "BP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("POS", "BP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("MarkerName", "SNP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("ID", "SNP", colnames(summary_stats))
  
  # Convert summary_stats to a keyed data.table for fast lookups
  setkey(summary_stats, SNP)
  
  return(summary_stats)
}

# Function to create Manhattan plot - FAIL FAST if requirements not met
create_manhattan_plot <- function(data_subset, trait, population, region, output_dir) {
  output_file <- file.path(output_dir,
    paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_manhattan.png"))
  
  # Check for required columns - FAIL FAST if missing
  required_cols <- c("CHR", "BP", "P", "SNP")
  missing_cols <- setdiff(required_cols, colnames(data_subset))
  if (length(missing_cols) > 0) {
    stop("CRITICAL ERROR: Missing required columns for Manhattan plot: ",
         paste(missing_cols, collapse = ", "),
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Sample rows if there are too many data points (speed optimization)
  if (nrow(data_subset) > 1000000) {
    set.seed(42)  # For reproducibility
    data_subset <- data_subset[sample(.N, 1000000)]
  }
  
  # Check for valid p-values - FAIL FAST if invalid
  if (any(is.na(data_subset$P)) || any(!is.finite(data_subset$P)) || any(data_subset$P <= 0)) {
    stop("CRITICAL ERROR: Invalid p-values found in data",
         "\nThis is a critical error - we cannot proceed with invalid data!")
  }
  
  # Ensure chromosome values are numeric - FAIL FAST if invalid
  if (is.character(data_subset$CHR)) {
    data_subset$CHR <- as.numeric(gsub("[^0-9]", "", data_subset$CHR))
    if (any(is.na(data_subset$CHR))) {
      stop("CRITICAL ERROR: Invalid chromosome values found",
           "\nThis is a critical error - we cannot proceed with invalid data!")
    }
  }
  
  # Check for valid chromosomes - FAIL FAST if invalid
  invalid_chromosomes <- !data_subset$CHR %in% 1:22
  if (any(invalid_chromosomes)) {
    stop("CRITICAL ERROR: Invalid chromosome numbers found: ",
         paste(unique(data_subset$CHR[invalid_chromosomes]), collapse = ", "),
         "\nThis is a critical error - we cannot proceed with invalid data!")
  }
  
  # Create an improved Manhattan plot with better aesthetics - FAIL FAST if creation fails
  if (!png(filename = output_file, width = 1200, height = 800, res = 150)) {
    stop("CRITICAL ERROR: Failed to create PNG file: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  # Customize colors for better readability
  col_palette <- c("#4286f4", "#f44182")
  
  manhattan(data_subset, chr = "CHR", bp = "BP", p = "P",
          main = paste0("MWAS ", toupper(trait), " - ", toupper(region), " - ", population),
          cex = 0.6, cex.axis = 0.8, col = col_palette,
          suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
          chrlabs = as.character(1:22))
  
  dev.off()
  
  # Verify plot was created - FAIL FAST if not
  if (!file.exists(output_file)) {
    stop("CRITICAL ERROR: Failed to save Manhattan plot: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  message("Created Manhattan plot: ", basename(output_file))
  return(output_file)
}

# Function to create QQ plot - FAIL FAST if requirements not met
create_qq_plot <- function(data_subset, trait, population, region, output_dir, max_points = 5000) {
  output_file <- file.path(output_dir,
    paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_qq.png"))
  
  # Check for required columns - FAIL FAST if missing
  if (!"p" %in% colnames(data_subset)) {
    stop("CRITICAL ERROR: Missing p-value column for QQ plot",
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Extract p-values and check for validity - FAIL FAST if invalid
  p_values <- data_subset$p
  if (any(is.na(p_values)) || any(!is.finite(p_values)) || any(p_values <= 0) || any(p_values > 1)) {
    stop("CRITICAL ERROR: Invalid p-values found in data",
         "\nThis is a critical error - we cannot proceed with invalid data!")
  }
  
  # Sample if we have more values than needed
  if (length(p_values) > max_points) {
    set.seed(42)  # For reproducibility
    p_values <- sample(p_values, max_points)
  }
  
  # Create QQ plot - FAIL FAST if creation fails
  if (!png(filename = output_file, width = 800, height = 800, res = 150)) {
    stop("CRITICAL ERROR: Failed to create PNG file: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  qqman::qq(p_values, main = paste0("QQ Plot: ", toupper(trait), " - ",
                                   toupper(region), " - ", population))
  
  dev.off()
  
  # Verify plot was created - FAIL FAST if not
  if (!file.exists(output_file)) {
    stop("CRITICAL ERROR: Failed to save QQ plot: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  message("Created QQ plot: ", basename(output_file))
  return(output_file)
}

# Function to create histogram of p-values - FAIL FAST if requirements not met
create_histogram_plot <- function(data_subset, trait, population, region, output_dir) {
  output_file <- file.path(output_dir,
    paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_hist.png"))
  
  # Check for required columns - FAIL FAST if missing
  if (!"p" %in% colnames(data_subset)) {
    stop("CRITICAL ERROR: Missing p-value column for histogram",
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Extract p-values and check for validity - FAIL FAST if invalid
  p_values <- data_subset$p
  if (any(is.na(p_values)) || any(!is.finite(p_values)) || any(p_values <= 0) || any(p_values > 1)) {
    stop("CRITICAL ERROR: Invalid p-values found in data",
         "\nThis is a critical error - we cannot proceed with invalid data!")
  }
  
  # Sample if we have too many points
  if (length(p_values) > 100000) {
    set.seed(42)
    p_values <- sample(p_values, 100000)
  }
  
  # Create histogram - FAIL FAST if creation fails
  if (!png(filename = output_file, width = 800, height = 800, res = 150)) {
    stop("CRITICAL ERROR: Failed to create PNG file: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  hist(p_values, breaks = 100, col = "skyblue",
       main = paste0("P-value distribution: ", toupper(trait), " - ",
                    toupper(region), " - ", population),
       xlab = "P-value")
  
  dev.off()
  
  # Verify plot was created - FAIL FAST if not
  if (!file.exists(output_file)) {
    stop("CRITICAL ERROR: Failed to save histogram: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  message("Created histogram: ", basename(output_file))
  return(output_file)
}

# Main processing loop for each trait
for (trait in traits) {
  message("\nProcessing trait: ", trait)
  
  # Input file path - FAIL FAST if file doesn't exist
  input_file <- paste0("/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/OLD_scripts/16a9par-OUT_stage2_MWAS_",
                      trait, ".csv")
  
  if (!file.exists(input_file)) {
    stop("CRITICAL ERROR: Input file not found: ", input_file,
         "\nThis is a critical error - we cannot proceed without input data!")
  }
  
  message("Loading file: ", input_file)
  
  # Load data efficiently - FAIL FAST if required columns are missing
  file_in <- fread(input_file, select = c("chr", "pos", "p", "n", "population", "region"),
                  nThread = 4)
  
  # Check for required columns - FAIL FAST if missing
  required_cols <- c("chr", "pos", "p", "n", "population", "region")
  missing_cols <- setdiff(required_cols, colnames(file_in))
  if (length(missing_cols) > 0) {
    stop("CRITICAL ERROR: Missing required columns in input file: ",
         paste(missing_cols, collapse = ", "),
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Check for empty data - FAIL FAST if empty
  if (nrow(file_in) == 0) {
    stop("CRITICAL ERROR: Input file is empty: ", input_file,
         "\nThis is a critical error - we cannot proceed without data!")
  }
  
  message("Successfully loaded data with ", nrow(file_in), " rows.")
  
  # Get unique combinations of population and region
  pop_region_combinations <- unique(file_in[, .(population, region)])
  message("Found ", nrow(pop_region_combinations), " population-region combinations.")
  
  # Process each combination
  for (i in 1:nrow(pop_region_combinations)) {
    pop <- pop_region_combinations$population[i]
    reg <- pop_region_combinations$region[i]
    message("\nProcessing combination ", i, "/", nrow(pop_region_combinations), 
            ": ", trait, " - ", pop, " - ", reg)
    
    # Filter data for this population and region
    file_in_subset <- file_in[population == pop & region == reg]
    
    # Check for empty subset - FAIL FAST if empty
    if (nrow(file_in_subset) == 0) {
      stop("CRITICAL ERROR: No data found for combination: ", trait, " - ", pop, " - ", reg,
           "\nThis is a critical error - we cannot proceed with empty data!")
    }
    
    # Limit to a reasonable number of rows for processing
    if (nrow(file_in_subset) > 1000000) {
      set.seed(42)
      file_in_subset <- file_in_subset[sample(.N, 1000000)]
    }
    
    # Check for duplicates - FAIL FAST if found
    if (any(duplicated(file_in_subset))) {
      stop("CRITICAL ERROR: Duplicate rows found in data subset",
           "\nThis is a critical error - we cannot proceed with duplicate data!")
    }
    
    message("Data subset has ", nrow(file_in_subset), " rows after sampling.")
    
    # Create plots - each function now has its own error handling
    create_manhattan_plot(file_in_subset, trait, pop, reg, output_dir)
    create_qq_plot(file_in_subset, trait, pop, reg, output_dir)
    create_histogram_plot(file_in_subset, trait, pop, reg, output_dir)
    
    # Free memory
    rm(file_in_subset)
  }
  
  # Free memory
  rm(file_in)
  gc()
}

message("\nFast plot generation complete! All outputs saved to ", output_dir)
message("Successfully processed ", length(traits), " traits.")
