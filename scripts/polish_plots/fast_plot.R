#!/usr/bin/env Rscript
# Fast plot generation script - optimized version
# This script is optimized for speed to quickly recreate plots from existing data files

# Libraries needed for visualization
library(qqman)
library(ggplot2)
library(data.table)

# Define absolute paths for input and output directories
BASE_DIR <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts"
INPUT_DIR <- file.path(BASE_DIR, "plot_outputs/manhattan_qq")
OUTPUT_DIR <- file.path(BASE_DIR, "polish_plots/plot_outputs/manhattan_qq_fast")

# Ensure output directory exists
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Print paths for verification
message("Base directory: ", BASE_DIR)
message("Input directory: ", INPUT_DIR)
message("Output directory: ", OUTPUT_DIR)

# Define variables - now processing all combinations
desired_subpopulations <- c("AA", "EA", "all")
desired_regions <- c("caud", "hippo", "dlpfc")
traits <- c("scz", "bp", "mdd")

# Function to create Manhattan plot (optimized)
create_manhattan_plot <- function(file_path, trait, population, region, output_dir) {
  output_file <- file.path(output_dir, paste0("manhattan_", trait, "_", population, "_", region, ".png"))
  
  # Use fread with select to only load the columns we need
  message("Loading minimal data for Manhattan plot...")
  data <- fread(file_path, select = c("chr", "cg", "p"))
  
  # Add SNP column efficiently
  data[, SNP := paste0(chr, ":", cg)]
  
  # Create Manhattan plot with minimal settings
  png(filename = output_file, width = 1200, height = 800, res = 150)
  manhattan(data, chr = "chr", bp = "cg", p = "p", 
            main = paste("Manhattan Plot:", trait, "-", region, "-", population),
            suggestiveline = -log10(1e-5), 
            genomewideline = -log10(5e-8))
  dev.off()
  
  message("Created Manhattan plot: ", output_file)
  return(TRUE)
}

# Improved QQ plot function with simpler sampling
create_qq_plot <- function(file_path, trait, population, region, output_dir, sample_size = 5000) {
  output_file <- file.path(output_dir, paste0("qqplot_", trait, "_", population, "_", region, ".png"))
  
  # Sample the p-values more efficiently
  message("Sampling p-values for QQ plot...")
  
  # Just read the p column with fread, which is already efficient
  p_values <- fread(file_path, select = "p")$p
  
  # Sample if we have more values than needed
  if (length(p_values) > sample_size) {
    set.seed(42)  # For reproducibility
    p_values <- sample(p_values, sample_size)
  }
  
  # Create QQ plot
  png(filename = output_file, width = 800, height = 800, res = 150)
  qq(p_values, main = paste("QQ Plot:", trait, "-", region, "-", population))
  dev.off()
  
  message("Created QQ plot: ", output_file)
  return(TRUE)
}

# Counter for progress tracking
total_combinations <- length(desired_subpopulations) * length(desired_regions) * length(traits)
current <- 0

# Process combinations
for (trait in traits) {
  for (population in desired_subpopulations) {
    for (region in desired_regions) {
      current <- current + 1
      message("\nProcessing combination ", current, "/", total_combinations, 
              ": ", trait, " - ", population, " - ", region)
      
      # File paths with absolute paths
      filtered_file <- file.path(INPUT_DIR, 
                              paste0("filtered_", trait, "_results_", population, "_", region, ".csv"))
      
      unfiltered_file <- file.path(INPUT_DIR, 
                                paste0("unfiltered_", trait, "_results_", population, "_", region, ".csv"))
      
      # Process Manhattan plot - fail fast if file doesn't exist
      if (!file.exists(filtered_file)) {
        stop("CRITICAL ERROR: Filtered file not found: ", filtered_file)
      }
      
      # No try-catch - let errors propagate and stop execution
      create_manhattan_plot(filtered_file, trait, population, region, OUTPUT_DIR)
      gc()  # Force garbage collection
      
      # Process QQ plot - fail fast if file doesn't exist
      if (!file.exists(unfiltered_file)) {
        stop("CRITICAL ERROR: Unfiltered file not found: ", unfiltered_file)
      }
      
      # No try-catch - let errors propagate and stop execution
      create_qq_plot(unfiltered_file, trait, population, region, OUTPUT_DIR)
      gc()  # Force garbage collection
    }
  }
}

message("\nFast plot generation complete! All outputs saved to ", OUTPUT_DIR)
message("Processed ", current, " combinations out of ", total_combinations, " total.")