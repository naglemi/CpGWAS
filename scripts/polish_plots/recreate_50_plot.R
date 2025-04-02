#!/usr/bin/env Rscript
# recreate_50_plot.R - Faithful recreation of plots from 50_make_bed_a6.ipynb
# This script processes the data in chunks to handle the 180GB file
# Modified to use absolute paths and fix logical issues
# Updated to use bash commands for extracting samples

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
output_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/50_plots_fast"
output_dir_50out <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/50-OUT"
temp_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/temp"

if (!dir.exists(output_dir)) {
  if (!dir.create(output_dir, recursive = TRUE)) {
    stop("CRITICAL ERROR: Failed to create output directory: ", output_dir)
  }
  cat(timestamp(), "- Created output directory:", output_dir, "\n")
} else {
  cat(timestamp(), "- Using existing output directory:", output_dir, "\n")
}

if (!dir.exists(output_dir_50out) && !dir.create(output_dir_50out, recursive = TRUE)) {
  cat(timestamp(), "- Warning: Could not create 50-OUT directory:", output_dir_50out, "\n")
}

if (!dir.exists(temp_dir)) {
  if (!dir.create(temp_dir, recursive = TRUE)) {
    stop("CRITICAL ERROR: Failed to create temp directory: ", temp_dir)
  }
}

# Define the path to the final merged data file (absolute path)
merged_data_file <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/final_merged_DNAm_data.csv"

# Check if the file exists
if (!file.exists(merged_data_file)) {
  stop("CRITICAL ERROR: Required input file not found: ", merged_data_file)
}
cat(timestamp(), "- Found merged data file at:", merged_data_file, "\n")

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

# Extract headers to get column names
cat(timestamp(), "- Extracting header from data file...\n")
header_file <- file.path(temp_dir, "header.csv")
system(paste("head -n 1", merged_data_file, ">", header_file))
header <- fread(header_file, header = TRUE, nrows = 0)
column_names <- names(header)
cat(timestamp(), "- Extracted header with", length(column_names), "columns\n")

# Use bash to create a small sample file for feature analysis
cat(timestamp(), "- Creating sample file for feature analysis...\n")
sample_file <- file.path(temp_dir, "sample.csv")
system(paste("head -n 1", merged_data_file, ">", sample_file, 
             "&& awk 'NR > 1 && NR <= 5001' ", merged_data_file, ">>", sample_file))
cat(timestamp(), "- Sample file created\n")

# Process small sample to identify feature types
cat(timestamp(), "- Processing sample to identify feature types...\n")
sample_data <- fread(sample_file)
if (!"Feature_INFO" %in% names(sample_data)) {
  stop("CRITICAL ERROR: Required column 'Feature_INFO' not found in data")
}

# Extract feature type from Feature_INFO
sample_data[, feature_type := str_split_fixed(Feature_INFO, "\\|", 2)[,1]]

# Get list of unique features from sample
unique_features <- unique(sample_data$feature_type)
unique_features <- unique_features[!is.na(unique_features)]
unique_features <- unique_features[!grepl("ttribute", unique_features)]
cat(timestamp(), "- Found", length(unique_features), "unique features in sample\n")

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

# Calculate total number of lines in the file
cat(timestamp(), "- Counting lines in data file (this may take a while)...\n")
lines_count <- as.numeric(system(paste("wc -l <", merged_data_file), intern = TRUE))
cat(timestamp(), "- File has", lines_count, "lines\n")

# Calculate number of chunks (aim for ~1GB per chunk)
chunk_size <- 1000000  # 1 million rows per chunk
num_chunks <- ceiling((lines_count - 1) / chunk_size)  # Subtract 1 for header
cat(timestamp(), "- Will process file in", num_chunks, "chunks of", chunk_size, "rows each\n")

# Process the file in chunks using bash to extract chunks
for (chunk_idx in 1:num_chunks) {
  start_line <- (chunk_idx - 1) * chunk_size + 2  # +2 to skip header and start from data
  end_line <- min(start_line + chunk_size - 1, lines_count)
  
  cat(timestamp(), "- Processing chunk", chunk_idx, "of", num_chunks, 
      "(lines", start_line, "to", end_line, ")...\n")
  
  # Extract chunk to temporary file using bash
  chunk_file <- file.path(temp_dir, paste0("chunk_", chunk_idx, ".csv"))
  system(paste("head -n 1", merged_data_file, ">", chunk_file, 
               "&& awk 'NR >=", start_line, "&& NR <=", end_line, "' ", 
               merged_data_file, ">>", chunk_file))
  
  # Read chunk and process
  tryCatch({
    chunk_data <- fread(chunk_file)
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
    chunk_data[, DNAm_Level_bin := cut(all_caud_Mean_DNAm_Level, 
                                       breaks = seq(dnam_min, dnam_max, length.out = 11), 
                                       include.lowest = TRUE)]
    
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
        merged_counts[, Overlapping_CpGs := ifelse(is.na(Overlapping_CpGs.x), 0, Overlapping_CpGs.x) + 
                                         ifelse(is.na(Overlapping_CpGs.y), 0, Overlapping_CpGs.y)]
        merged_counts[, Total_CpGs := ifelse(is.na(Total_CpGs.x), 0, Total_CpGs.x) + 
                                    ifelse(is.na(Total_CpGs.y), 0, Total_CpGs.y)]
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
        merged_counts[, Overlapping_CpGs := ifelse(is.na(Overlapping_CpGs.x), 0, Overlapping_CpGs.x) + 
                                         ifelse(is.na(Overlapping_CpGs.y), 0, Overlapping_CpGs.y)]
        merged_counts[, Total_CpGs := ifelse(is.na(Total_CpGs.x), 0, Total_CpGs.x) + 
                                    ifelse(is.na(Total_CpGs.y), 0, Total_CpGs.y)]
        merged_counts[, c("Overlapping_CpGs.x", "Overlapping_CpGs.y", "Total_CpGs.x", "Total_CpGs.y") := NULL]
        feature_dnam_counts[[feature]] <- merged_counts
      }
    }
    
    # Clean up to free memory
    rm(chunk_data)
    gc()
    
    # Delete the chunk file to save space
    if (file.exists(chunk_file)) {
      file.remove(chunk_file)
    }
  }, error = function(e) {
    cat(timestamp(), "- WARNING: Error processing chunk", chunk_idx, ":", conditionMessage(e), "\n")
    # Continue with next chunk even if this one fails
  })
}

cat(timestamp(), "- Finished processing chunks\n")

# Initialize list for statistical results
statistical_results <- list()

# Generate plots for each feature using the aggregated counts
for (feature in unique_features) {
  col_name <- sanitize_filename(feature)
  feature_hr <- get_human_readable_feature_name(feature)
  
  cat(timestamp(), "- Generating plots for feature:", feature_hr, "\n")
  
  # Get the aggregated counts for this feature
  r2_counts <- feature_r2_counts[[feature]]
  dnam_counts <- feature_dnam_counts[[feature]]
  
  # Check if we have data for this feature
  if (nrow(r2_counts) == 0 || nrow(dnam_counts) == 0) {
    cat(timestamp(), "- WARNING: No data for feature:", feature_hr, "- skipping\n")
    next
  }
  
  # Calculate proportion for R² bins
  r2_counts[, Proportion := Overlapping_CpGs / Total_CpGs]
  
  # Generate the plot for R² bins
  p_r2 <- ggplot(r2_counts, aes(x = R_squared_bin, y = Proportion)) +
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
  
  # Calculate proportion for DNAm Level bins
  dnam_counts[, Proportion := Overlapping_CpGs / Total_CpGs]
  
  # Generate the plot for DNAm Level bins
  p_dnam <- ggplot(dnam_counts, aes(x = DNAm_Level_bin, y = Proportion)) +
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
  
  # Save the combined plot
  plot_file <- file.path(output_dir, paste0("overlap_plot_", col_name, "_combined.png"))
  ggsave(filename = plot_file, plot = combined_plot, width = 12, height = 16, dpi = 300)
  
  cat(timestamp(), "- Combined plot saved to", plot_file, "\n")
  
  # Statistical Analysis for R² bins
  # First, create a contingency table from our aggregated data
  r2_bins <- unique(r2_counts$R_squared_bin)
  contingency_data <- data.frame(
    bin = rep(r2_bins, each = 2),
    overlap = rep(c(0, 1), length(r2_bins)),
    count = numeric(2 * length(r2_bins))
  )
  
  for (i in 1:nrow(contingency_data)) {
    bin <- contingency_data$bin[i]
    overlap <- contingency_data$overlap[i]
    
    if (overlap == 1) {
      contingency_data$count[i] <- r2_counts[R_squared_bin == bin, Overlapping_CpGs]
    } else {
      contingency_data$count[i] <- r2_counts[R_squared_bin == bin, Total_CpGs - Overlapping_CpGs]
    }
  }
  
  contingency_table_r2 <- xtabs(count ~ overlap + bin, data = contingency_data)
  
  # Perform chi-squared test if appropriate
  if (nrow(contingency_table_r2) == 2 && ncol(contingency_table_r2) > 1) {
    # Calculate expected counts
    expected_counts_r2 <- suppressWarnings(chisq.test(contingency_table_r2)$expected)
    
    # Check if expected counts are sufficient for chi-squared test
    if (all(expected_counts_r2 >= 5)) {
      chi_squared_test_r2 <- suppressWarnings(chisq.test(contingency_table_r2))
      test_used_r2 <- "Chi-squared Test"
      chi_stat_r2 <- chi_squared_test_r2$statistic
      chi_df_r2 <- chi_squared_test_r2$parameter
      chi_p_r2 <- chi_squared_test_r2$p.value
    } else {
      # Use Fisher's exact test for small expected counts
      fisher_test_r2 <- suppressWarnings(fisher.test(contingency_table_r2, simulate.p.value = TRUE))
      test_used_r2 <- "Fisher's Exact Test"
      chi_stat_r2 <- NA
      chi_df_r2 <- NA
      chi_p_r2 <- fisher_test_r2$p.value
    }
  } else {
    test_used_r2 <- "No Test (Insufficient Data)"
    chi_stat_r2 <- NA
    chi_df_r2 <- NA
    chi_p_r2 <- NA
  }
  
  # Repeat for DNAm Level bins
  dnam_bins <- unique(dnam_counts$DNAm_Level_bin)
  contingency_data <- data.frame(
    bin = rep(dnam_bins, each = 2),
    overlap = rep(c(0, 1), length(dnam_bins)),
    count = numeric(2 * length(dnam_bins))
  )
  
  for (i in 1:nrow(contingency_data)) {
    bin <- contingency_data$bin[i]
    overlap <- contingency_data$overlap[i]
    
    if (overlap == 1) {
      contingency_data$count[i] <- dnam_counts[DNAm_Level_bin == bin, Overlapping_CpGs]
    } else {
      contingency_data$count[i] <- dnam_counts[DNAm_Level_bin == bin, Total_CpGs - Overlapping_CpGs]
    }
  }
  
  contingency_table_dnam <- xtabs(count ~ overlap + bin, data = contingency_data)
  
  # Perform chi-squared test if appropriate
  if (nrow(contingency_table_dnam) == 2 && ncol(contingency_table_dnam) > 1) {
    # Calculate expected counts
    expected_counts_dnam <- suppressWarnings(chisq.test(contingency_table_dnam)$expected)
    
    # Check if expected counts are sufficient for chi-squared test
    if (all(expected_counts_dnam >= 5)) {
      chi_squared_test_dnam <- suppressWarnings(chisq.test(contingency_table_dnam))
      test_used_dnam <- "Chi-squared Test"
      chi_stat_dnam <- chi_squared_test_dnam$statistic
      chi_df_dnam <- chi_squared_test_dnam$parameter
      chi_p_dnam <- chi_squared_test_dnam$p.value
    } else {
      # Use Fisher's exact test for small expected counts
      fisher_test_dnam <- suppressWarnings(fisher.test(contingency_table_dnam, simulate.p.value = TRUE))
      test_used_dnam <- "Fisher's Exact Test"
      chi_stat_dnam <- NA
      chi_df_dnam <- NA
      chi_p_dnam <- fisher_test_dnam$p.value
    }
  } else {
    test_used_dnam <- "No Test (Insufficient Data)"
    chi_stat_dnam <- NA
    chi_df_dnam <- NA
    chi_p_dnam <- NA
  }
  
  # Collect results
  statistical_results[[feature_hr]] <- list(
    Feature = feature_hr,
    R2_Test = test_used_r2,
    R2_Statistic = round(chi_stat_r2, 2),
    R2_Degrees_of_Freedom = chi_df_r2,
    R2_P_value = signif(chi_p_r2, 3),
    DNAm_Test = test_used_dnam,
    DNAm_Statistic = round(chi_stat_dnam, 2),
    DNAm_Degrees_of_Freedom = chi_df_dnam,
    DNAm_P_value = signif(chi_p_dnam, 3)
  )
}

# Convert statistical_results to a data.table
stats_table <- rbindlist(lapply(statistical_results, function(x) as.data.table(x)), fill = TRUE)

# Save statistical results to a CSV file
stats_file <- file.path(output_dir, "statistical_results.csv")
fwrite(stats_table, stats_file)

# Clean up temp directory
unlink(file.path(temp_dir, "*"), recursive = TRUE, force = TRUE)

cat(timestamp(), "- Statistical results saved to", stats_file, "\n")
cat(timestamp(), "- Plot recreation completed successfully!\n")