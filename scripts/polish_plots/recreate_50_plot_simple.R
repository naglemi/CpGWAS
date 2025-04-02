#!/usr/bin/env Rscript
# recreate_50_plot_simple.R - Simplified script to recreate plots from 50_make_bed_a6.ipynb
# This script skips the data processing steps and focuses only on generating the plots

# Load necessary libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(gtools)
  library(RColorBrewer)
  library(gridExtra)
})

# Initialize timestamp function for logging
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Create output directories with proper naming convention
output_dir <- "./plot_outputs/50_plots_fast"
output_dir_50out <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/50-OUT"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir_50out)) {
  dir.create(output_dir_50out, recursive = TRUE, showWarnings = FALSE)
}
cat(timestamp(), "- Created output directories:", output_dir, "and", output_dir_50out, "\n")

# Define the path to the final merged data file
merged_data_file <- "../final_merged_DNAm_data.csv"

# Check if the file exists
if (!file.exists(merged_data_file)) {
  # Try alternative locations
  alt_locations <- c(
    "./final_merged_DNAm_data.csv",
    "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/final_merged_DNAm_data.csv"
  )
  
  found <- FALSE
  for (loc in alt_locations) {
    if (file.exists(loc)) {
      merged_data_file <- loc
      found <- TRUE
      cat(timestamp(), "- Found merged data file at:", merged_data_file, "\n")
      break
    }
  }
  
  if (!found) {
    stop("Error: Required input file not found: ", merged_data_file)
  }
}

# Load the data in chunks
cat(timestamp(), "- Loading merged data file in chunks (this may take a while)...\n")

# First, count the number of lines in the file
cmd <- paste("wc -l", merged_data_file)
lines_count <- as.numeric(strsplit(system(cmd, intern = TRUE), " ")[[1]][1])
cat(timestamp(), "- File has", lines_count, "lines\n")

# Calculate number of chunks (aim for ~1GB per chunk)
chunk_size <- 1000000  # 1 million rows per chunk
num_chunks <- ceiling(lines_count / chunk_size)
cat(timestamp(), "- Will process file in", num_chunks, "chunks of", chunk_size, "rows each\n")

# Process a small sample first to get feature types
cat(timestamp(), "- Loading a small sample to identify feature types...\n")
sample_data <- fread(merged_data_file, nrows = 100000)
cat(timestamp(), "- Loaded sample with", nrow(sample_data), "rows\n")

# Extract feature type from Feature_INFO
sample_data[, feature_type := str_split_fixed(Feature_INFO, "\\|", 2)[,1]]

# Get list of unique features from sample
unique_features <- unique(sample_data$feature_type)
unique_features <- unique_features[!is.na(unique_features)]
unique_features <- unique_features[!grepl("ttribute", unique_features)]
cat(timestamp(), "- Found", length(unique_features), "unique features\n")

# Create a data structure to store aggregated counts
feature_r2_counts <- list()
feature_dnam_counts <- list()

for (feature in unique_features) {
  feature_r2_counts[[feature]] <- data.table(
    R_squared_bin = character(),
    Overlapping_CpGs = integer(),
    Total_CpGs = integer()
  )
  
  feature_dnam_counts[[feature]] <- data.table(
    DNAm_Level_bin = character(),
    Overlapping_CpGs = integer(),
    Total_CpGs = integer()
  )
}

# Process the file in chunks - process ALL chunks, not just the first 3
for (chunk_idx in 1:num_chunks) {  # Process all chunks to ensure complete data processing
  skip_rows <- (chunk_idx - 1) * chunk_size
  cat(timestamp(), "- Processing chunk", chunk_idx, "of", num_chunks, "(rows", skip_rows + 1, "to", skip_rows + chunk_size, ")...\n")
  
  # Read chunk
  chunk_data <- fread(merged_data_file, skip = skip_rows, nrows = chunk_size)
  cat(timestamp(), "- Loaded chunk with", nrow(chunk_data), "rows\n")
  
  # Extract feature type from Feature_INFO
  chunk_data[, feature_type := str_split_fixed(Feature_INFO, "\\|", 2)[,1]]
  
  # Remove rows with attribute features
  chunk_data <- chunk_data[!grepl("ttribute", feature_type)]
  
  # Ensure feature_type is character
  chunk_data[, feature_type := as.character(feature_type)]
  
  # Create R_squared_bin based on 'all_caud_cor'
  chunk_data[, R_squared_bin := cut(all_caud_cor, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)]
  
  # Create DNAm_Level_bin based on 'all_caud_Mean_DNAm_Level'
  dnam_min <- min(chunk_data$all_caud_Mean_DNAm_Level, na.rm = TRUE)
  dnam_max <- max(chunk_data$all_caud_Mean_DNAm_Level, na.rm = TRUE)
  chunk_data[, DNAm_Level_bin := cut(all_caud_Mean_DNAm_Level, breaks = seq(dnam_min, dnam_max, length.out = 11), include.lowest = TRUE)]
  
  # Process each feature
  for (feature in unique_features) {
    # Count by R² bin
    r2_counts <- chunk_data[feature_type == feature, .N, by = R_squared_bin]
    setnames(r2_counts, "N", "Overlapping_CpGs")
    
    # Count total CpGs per R² bin
    total_r2_counts <- chunk_data[, .N, by = R_squared_bin]
    setnames(total_r2_counts, "N", "Total_CpGs")
    
    # Merge and add to feature_r2_counts
    r2_counts <- merge(r2_counts, total_r2_counts, by = "R_squared_bin", all = TRUE)
    r2_counts[is.na(Overlapping_CpGs), Overlapping_CpGs := 0]
    
    # Append to existing counts
    if (nrow(feature_r2_counts[[feature]]) == 0) {
      feature_r2_counts[[feature]] <- r2_counts
    } else {
      # Merge with existing counts
      existing_counts <- feature_r2_counts[[feature]]
      merged_counts <- merge(existing_counts, r2_counts, by = "R_squared_bin", all = TRUE)
      merged_counts[, Overlapping_CpGs := Overlapping_CpGs.x + Overlapping_CpGs.y]
      merged_counts[, Total_CpGs := Total_CpGs.x + Total_CpGs.y]
      merged_counts[, c("Overlapping_CpGs.x", "Overlapping_CpGs.y", "Total_CpGs.x", "Total_CpGs.y") := NULL]
      feature_r2_counts[[feature]] <- merged_counts
    }
    
    # Count by DNAm Level bin
    dnam_counts <- chunk_data[feature_type == feature, .N, by = DNAm_Level_bin]
    setnames(dnam_counts, "N", "Overlapping_CpGs")
    
    # Count total CpGs per DNAm Level bin
    total_dnam_counts <- chunk_data[, .N, by = DNAm_Level_bin]
    setnames(total_dnam_counts, "N", "Total_CpGs")
    
    # Merge and add to feature_dnam_counts
    dnam_counts <- merge(dnam_counts, total_dnam_counts, by = "DNAm_Level_bin", all = TRUE)
    dnam_counts[is.na(Overlapping_CpGs), Overlapping_CpGs := 0]
    
    # Append to existing counts
    if (nrow(feature_dnam_counts[[feature]]) == 0) {
      feature_dnam_counts[[feature]] <- dnam_counts
    } else {
      # Merge with existing counts
      existing_counts <- feature_dnam_counts[[feature]]
      merged_counts <- merge(existing_counts, dnam_counts, by = "DNAm_Level_bin", all = TRUE)
      merged_counts[, Overlapping_CpGs := Overlapping_CpGs.x + Overlapping_CpGs.y]
      merged_counts[, Total_CpGs := Total_CpGs.x + Total_CpGs.y]
      merged_counts[, c("Overlapping_CpGs.x", "Overlapping_CpGs.y", "Total_CpGs.x", "Total_CpGs.y") := NULL]
      feature_dnam_counts[[feature]] <- merged_counts
    }
  }
  
  # Clean up to free memory
  rm(chunk_data)
  gc()
}

cat(timestamp(), "- Finished processing chunks\n")

# Extract feature type from Feature_INFO
merged_data[, feature_type := str_split_fixed(Feature_INFO, "\\|", 2)[,1]]

# Remove rows with attribute features
merged_data <- merged_data[!grepl("ttribute", feature_type)]

# Ensure feature_type is character
merged_data[, feature_type := as.character(feature_type)]

# Create R_squared_bin based on 'all_caud_cor'
cat(timestamp(), "- Binning R² values into deciles...\n")
merged_data[, R_squared_bin := cut(all_caud_cor, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)]

# Ensure bins are ordered correctly
merged_data$R_squared_bin <- factor(merged_data$R_squared_bin, levels = mixedsort(unique(merged_data$R_squared_bin)))

# Add "Missing" to the levels and assign
levels(merged_data$R_squared_bin) <- c(levels(merged_data$R_squared_bin), "Missing")
merged_data$R_squared_bin[is.na(merged_data$R_squared_bin)] <- 'Missing'

# Create DNAm_Level_bin based on 'all_caud_Mean_DNAm_Level'
cat(timestamp(), "- Binning DNAm levels into deciles...\n")
dnam_min <- min(merged_data$all_caud_Mean_DNAm_Level, na.rm = TRUE)
dnam_max <- max(merged_data$all_caud_Mean_DNAm_Level, na.rm = TRUE)
merged_data[, DNAm_Level_bin := cut(all_caud_Mean_DNAm_Level, breaks = seq(dnam_min, dnam_max, length.out = 11), include.lowest = TRUE)]

# Ensure bins are ordered correctly
merged_data$DNAm_Level_bin <- factor(merged_data$DNAm_Level_bin, levels = mixedsort(unique(merged_data$DNAm_Level_bin)))

# Add "Missing" to the levels and assign
levels(merged_data$DNAm_Level_bin) <- c(levels(merged_data$DNAm_Level_bin), "Missing")
merged_data$DNAm_Level_bin[is.na(merged_data$DNAm_Level_bin)] <- 'Missing'

# Function to sanitize feature names for file names
sanitize_filename <- function(filename) {
  gsub("[/\\?%*:|\"<> ]", "_", filename)
}

# Function to get human-readable feature name
get_human_readable_feature_name <- function(feature) {
  feature <- gsub("_", " ", feature)
  feature <- sub("E[0-9]+_[0-9]+_coreMarks_mnemonics\\|feature:", "", feature)
  return(feature)
}

# Create merged dataset with one row per CpG
merged_data_unique <- merged_data[, .(
  all_caud_cor = unique(all_caud_cor),
  R_squared_bin = unique(R_squared_bin),
  all_caud_Mean_DNAm_Level = unique(all_caud_Mean_DNAm_Level),
  DNAm_Level_bin = unique(DNAm_Level_bin),
  feature_list = list(feature_type)
), by = cg]

# Flatten the feature_list and remove NAs
merged_data_unique[, feature_list := lapply(feature_list, function(x) unique(x[!is.na(x)]))]

# Get list of unique features
unique_features <- unique(unlist(merged_data$feature_type))
unique_features <- unique_features[!is.na(unique_features)]

# For each feature, create an indicator column
for (feature in unique_features) {
  col_name <- sanitize_filename(feature)
  merged_data_unique[, (paste0('overlap_', col_name)) := sapply(feature_list, function(x) as.integer(feature %in% x))]
}

# Initialize list for statistical results
statistical_results <- list()

# Loop through each feature to generate plots and perform statistical analysis
for (feature in unique_features) {
  col_name <- sanitize_filename(feature)
  feature_hr <- get_human_readable_feature_name(feature)
  
  cat(timestamp(), "- Processing feature:", feature_hr, "\n")
  
  # Extract data for this feature
  overlap_col <- paste0('overlap_', col_name)
  feature_data <- merged_data_unique[, .(cg, R_squared_bin, DNAm_Level_bin, Overlap = get(overlap_col))]
  
  # Calculate total CpGs per R² bin
  total_counts_r2 <- feature_data[, .N, by = R_squared_bin]
  setnames(total_counts_r2, "N", "Total_CpGs")
  
  # Calculate number of overlapping CpGs per R² bin
  overlap_counts_r2 <- feature_data[, .(Overlapping_CpGs = sum(Overlap, na.rm = TRUE)), by = R_squared_bin]
  
  # Merge to get total counts for R² bins
  overlap_counts_r2 <- merge(overlap_counts_r2, total_counts_r2, by = "R_squared_bin", all.x = TRUE)
  
  # Calculate proportion for R² bins
  overlap_counts_r2[, Proportion := Overlapping_CpGs / Total_CpGs]
  
  # Generate the plot for R² bins
  p_r2 <- ggplot(overlap_counts_r2, aes(x = R_squared_bin, y = Proportion)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal(base_size = 16) +
    labs(
      title = paste0("Proportion of CpGs Overlapping with ", feature_hr, " Across R² Bins"),
      x = expression(R^2 ~ "Bins"),
      y = "Proportion of CpGs Overlapping"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  # Calculate total CpGs per DNAm Level bin
  total_counts_dnam <- feature_data[, .N, by = DNAm_Level_bin]
  setnames(total_counts_dnam, "N", "Total_CpGs")
  
  # Calculate number of overlapping CpGs per DNAm Level bin
  overlap_counts_dnam <- feature_data[, .(Overlapping_CpGs = sum(Overlap, na.rm = TRUE)), by = DNAm_Level_bin]
  
  # Merge to get total counts for DNAm Level bins
  overlap_counts_dnam <- merge(overlap_counts_dnam, total_counts_dnam, by = "DNAm_Level_bin", all.x = TRUE)
  
  # Calculate proportion for DNAm Level bins
  overlap_counts_dnam[, Proportion := Overlapping_CpGs / Total_CpGs]
  
  # Generate the plot for DNAm Level bins
  p_dnam <- ggplot(overlap_counts_dnam, aes(x = DNAm_Level_bin, y = Proportion)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    theme_minimal(base_size = 16) +
    labs(
      title = paste0("Proportion of CpGs Overlapping with ", feature_hr, " Across DNAm Level Bins"),
      x = "Mean DNAm Level Bins",
      y = "Proportion of CpGs Overlapping"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  # Stack the two plots using gridExtra
  combined_plot <- grid.arrange(p_r2, p_dnam, ncol = 1)
  
  # Save the combined plot with the correct naming convention (50-OUT_overlap_plot_...)
  plot_file <- file.path(output_dir, paste0("overlap_plot_", col_name, "_combined.png"))
  plot_file_50out <- file.path(output_dir_50out, paste0("50-OUT_overlap_plot_", col_name, "_combined.png"))
  ggsave(filename = plot_file, plot = combined_plot, width = 12, height = 16, dpi = 300)
  ggsave(filename = plot_file_50out, plot = combined_plot, width = 12, height = 16, dpi = 300)
  
  # Also save the individual plots with the correct naming convention
  plot_file_standard <- file.path(output_dir_50out, paste0("50-OUT_overlap_plot_", col_name, ".png"))
  ggsave(filename = plot_file_standard, plot = p_r2, width = 12, height = 8, dpi = 300)
  
  cat(timestamp(), "- Combined plot saved to", plot_file, "and", plot_file_50out, "\n")
  cat(timestamp(), "- Standard plot saved to", plot_file_standard, "\n")
  
  # Statistical Analysis for R² bins
  contingency_table_r2 <- table(feature_data$Overlap, feature_data$R_squared_bin)
  
  # Proceed if the contingency table has appropriate dimensions
  if (nrow(contingency_table_r2) == 2 && ncol(contingency_table_r2) > 1) {
    # Suppress warnings from chisq.test
    suppressWarnings({
      # Calculate expected counts
      expected_counts_r2 <- chisq.test(contingency_table_r2)$expected
    })
    if (any(expected_counts_r2 < 5)) {
      # Use Monte Carlo simulation
      chi_squared_test_r2 <- chisq.test(contingency_table_r2, simulate.p.value = TRUE, B = 1e5)
      test_used_r2 <- "Chi-squared Test with Monte Carlo Simulation"
      cat("Note: Monte Carlo simulation used for feature:", feature_hr, "in R² bins due to low expected counts.\n")
    } else {
      chi_squared_test_r2 <- chisq.test(contingency_table_r2)
      test_used_r2 <- "Chi-squared Test"
    }
  } else {
    chi_squared_test_r2 <- NULL
    test_used_r2 <- "Insufficient data for statistical test in R² bins"
    cat("Insufficient data for statistical test for feature:", feature_hr, "in R² bins\n")
  }
  
  # Statistical Analysis for DNAm Level bins
  contingency_table_dnam <- table(feature_data$Overlap, feature_data$DNAm_Level_bin)
  
  # Proceed if the contingency table has appropriate dimensions
  if (nrow(contingency_table_dnam) == 2 && ncol(contingency_table_dnam) > 1) {
    # Suppress warnings from chisq.test
    suppressWarnings({
      # Calculate expected counts
      expected_counts_dnam <- chisq.test(contingency_table_dnam)$expected
    })
    if (any(expected_counts_dnam < 5)) {
      # Use Monte Carlo simulation
      chi_squared_test_dnam <- chisq.test(contingency_table_dnam, simulate.p.value = TRUE, B = 1e5)
      test_used_dnam <- "Chi-squared Test with Monte Carlo Simulation"
      cat("Note: Monte Carlo simulation used for feature:", feature_hr, "in DNAm Level bins due to low expected counts.\n")
    } else {
      chi_squared_test_dnam <- chisq.test(contingency_table_dnam)
      test_used_dnam <- "Chi-squared Test"
    }
  } else {
    chi_squared_test_dnam <- NULL
    test_used_dnam <- "Insufficient data for statistical test in DNAm Level bins"
    cat("Insufficient data for statistical test for feature:", feature_hr, "in DNAm Level bins\n")
  }
  
  # Collect results
  statistical_results[[feature_hr]] <- list(
    Feature = feature_hr,
    R2_Test = test_used_r2,
    R2_Statistic = if (!is.null(chi_squared_test_r2)) round(chi_squared_test_r2$statistic, 2) else NA,
    R2_Degrees_of_Freedom = if (!is.null(chi_squared_test_r2)) chi_squared_test_r2$parameter else NA,
    R2_P_value = if (!is.null(chi_squared_test_r2)) signif(chi_squared_test_r2$p.value, 3) else NA,
    DNAm_Test = test_used_dnam,
    DNAm_Statistic = if (!is.null(chi_squared_test_dnam)) round(chi_squared_test_dnam$statistic, 2) else NA,
    DNAm_Degrees_of_Freedom = if (!is.null(chi_squared_test_dnam)) chi_squared_test_dnam$parameter else NA,
    DNAm_P_value = if (!is.null(chi_squared_test_dnam)) signif(chi_squared_test_dnam$p.value, 3) else NA
  )
}

# Convert statistical_results to a data.table
stats_table <- rbindlist(lapply(statistical_results, function(x) as.data.table(x)), fill = TRUE)

# Save statistical results to CSV files in both directories
stats_file <- file.path(output_dir, "statistical_results.csv")
stats_file_50out <- file.path(output_dir_50out, "50-OUT_statistical_results.csv")
fwrite(stats_table, stats_file)
fwrite(stats_table, stats_file_50out)

cat(timestamp(), "- Statistical results saved to", stats_file, "and", stats_file_50out, "\n")
cat(timestamp(), "- Plot recreation completed successfully!\n")