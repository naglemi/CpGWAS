#!/usr/bin/env Rscript
# recreate_41_plot_fixed.R - Fixed version of the script to recreate 41 plots
# This script addresses issues with input file paths and improves file finding logic

# Load necessary libraries
library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)       # For data manipulation
library(scales)      # For better scales in plots
library(tidyr)       # For pivot_wider
library(ggrepel)     # For non-overlapping labels

# Set seed for reproducibility
set.seed(123)

# -----------------------------
# 1. User-Defined Variables
# -----------------------------

# Number of samples per population-region group
n_sample <- Inf  # Set to Inf for full data without sampling

# Indicate whether to include MHC region in plots
include_MHC <- FALSE  # Set to TRUE to include MHC region

# -----------------------------
# 2. Create Output Directory
# -----------------------------

# Use absolute path for output directory as required by the mission
output_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/41_plots_fast"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste("Created output directory:", output_dir))
} else {
  message(paste("Output directory already exists:", output_dir))
}

# -----------------------------
# 3. Define Trait Mapping
# -----------------------------

trait_mapping <- list(
  bp = "Bipolar Disorder",
  mdd = "Major Depressive Disorder",
  scz = "Schizophrenia"
  # Add more mappings as needed
)

# -----------------------------
# 4. Define Folder Path and Files
# -----------------------------

# FIXED: Use a more comprehensive approach to find input files
potential_paths <- c(
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/37-OUT/",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/OLD/37-OUT/",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/OLD_scripts/37-OUT/",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/plot_outputs/37-OUT/"
)

# Find the first valid path
folder_path <- NULL
for (path in potential_paths) {
  if (dir.exists(path)) {
    folder_path <- path
    message(paste("Using folder path:", folder_path))
    break
  }
}

if (is.null(folder_path)) {
  # If none of the predefined paths exist, search for 37-OUT directory
  message("Searching for 37-OUT directory...")
  search_cmd <- "find /expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS -name '37-OUT' -type d | head -1"
  folder_path <- system(search_cmd, intern = TRUE)
  
  if (length(folder_path) == 0 || folder_path == "") {
    stop("CRITICAL ERROR: Could not find 37-OUT directory.")
  }
  
  message(paste("Found 37-OUT directory:", folder_path))
}

# FIXED: Use a more flexible pattern to find input files
# Look for any CSV files with 'bonf' in the name
files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
files <- files[grepl("bonf", files)]

if (length(files) == 0) {
  # If no files found, try a broader search
  message("No files found with initial pattern. Trying broader search...")
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("CRITICAL ERROR: No input files found. Please check the folder path.")
  }
}

message(paste("Found", length(files), "input files"))

# -----------------------------
# 5. Identify Unique Traits and Regions
# -----------------------------

# Extract traits and regions from filenames
filename_info <- data.table(filename = basename(files))

# Add a safety check to print the filenames for debugging
message("Sample filenames:")
message(paste(head(filename_info$filename, 3), collapse = "\n"))

# FIXED: Use a more robust approach to extract information from filenames
filename_info[, trait_code := ifelse(grepl("_bp_|_bp\\.", filename), "bp",
                             ifelse(grepl("_mdd_|_mdd\\.", filename), "mdd",
                             ifelse(grepl("_scz_|_scz\\.", filename), "scz", NA_character_)))]

filename_info[, population := ifelse(grepl("_AA_|_AA\\.", filename), "AA",
                             ifelse(grepl("_EA_|_EA\\.", filename), "EA",
                             ifelse(grepl("_all_|_all\\.", filename), "all", NA_character_)))]

filename_info[, region := ifelse(grepl("_caud|_caud\\.", filename), "caud",
                         ifelse(grepl("_dlpfc|_dlpfc\\.", filename), "dlpfc",
                         ifelse(grepl("_hippo|_hippo\\.", filename), "hippo", NA_character_)))]

# Remove any file entries with missing information
filename_info <- filename_info[!is.na(trait_code) & !is.na(population) & !is.na(region)]

# Print summary of extracted information
message("Extracted information from filenames:")
message(paste("Traits:", paste(unique(filename_info$trait_code), collapse = ", ")))
message(paste("Populations:", paste(unique(filename_info$population), collapse = ", ")))
message(paste("Regions:", paste(unique(filename_info$region), collapse = ", ")))

unique_traits <- unique(filename_info$trait_code)
unique_regions <- unique(filename_info$region)

# -----------------------------
# 6. Initialize Lists to Collect Data
# -----------------------------

# Initialize list to collect all mismatched z-stats data
all_mismatched_z_list <- list()

# Initialize list to collect all data for scatter plots
all_scatter_data_list <- list()

# -----------------------------
# 7. Function to Check and Handle Duplicates
# -----------------------------

check_and_handle_duplicates <- function(dt, step_description) {
  # Identify duplicate rows
  dup_rows <- duplicated(dt)
  num_dups <- sum(dup_rows)
  
  if (num_dups > 0) {
    warning(paste0("Duplicate rows detected at step: ", step_description, 
                   ". Number of duplicates: ", num_dups, "."))
    
    # Extract duplicate details
    dup_details <- dt[dup_rows]
    message("Details of duplicates:")
    print(head(dup_details, n = 10))  # Show first 10 duplicates for brevity
    
    # Remove duplicate rows, keeping the first occurrence
    dt_unique <- unique(dt)
    return(dt_unique)
  } else {
    message(paste0("No duplicate rows found at step: ", step_description, "."))
    return(dt)
  }
}

# -----------------------------
# 8. Process Data Per Trait and Region
# -----------------------------

for (trait_code in unique_traits) {
  # Map trait code to full name
  trait <- trait_mapping[[trait_code]]
  if (is.null(trait)) {
    trait <- trait_code
    message(paste("Trait code", trait_code, "not found in mapping. Using trait code as is."))
  }

  for (region in unique_regions) {
    # FIXED: More flexible file matching
    trait_region_files <- files[grepl(paste0(trait_code, ".*", region), files, ignore.case = TRUE)]

    if (length(trait_region_files) == 0) {
      message(paste("No files found for trait", trait_code, "and region", region))
      next
    }

    # Initialize list to store data for the current trait and region
    data_list <- list()

    # Read and process data files for the current trait and region
    for (file in trait_region_files) {
      # Extract filename without path
      filename <- basename(file)

      # FIXED: More flexible population extraction
      population <- NA_character_
      if (grepl("_AA_|_AA\\.", filename)) {
        population <- "AA"
      } else if (grepl("_EA_|_EA\\.", filename)) {
        population <- "EA"
      } else if (grepl("_all_|_all\\.", filename)) {
        population <- "all"
      }

      if (is.na(population)) {
        message(paste("Could not extract population from filename:", filename))
        next
      }

      # Read the data
      tryCatch({
        dt <- fread(file, select = c("chr", "cg", "z", "p", "n", "cor", "mse"))
        
        # Check for duplicates in the raw data
        dt <- check_and_handle_duplicates(dt, paste("Reading file:", filename))
        
        # Clean 'population' and 'region' columns
        dt[, population := str_trim(toupper(population))]
        dt[, region := str_trim(tolower(region))]
        
        # Add trait column
        dt[, trait := trait]
        
        # Append to the list
        data_list[[length(data_list) + 1]] <- dt
        
        message(paste("Successfully read file:", filename))
      }, error = function(e) {
        message(paste("Error reading file:", filename, "-", e$message))
      })
    }

    if (length(data_list) == 0) {
      message(paste("No data files were read successfully for trait", trait, "and region", region))
      next
    }

    # Combine data for the current trait and region
    combined_dt <- rbindlist(data_list, fill = TRUE)
    message(paste("Combined data for trait", trait, "and region", region))

    # Check for duplicates after combining
    combined_dt <- check_and_handle_duplicates(combined_dt, paste("Combining data for trait:", trait, "and region:", region))

    # -----------------------------
    # 9. Sampling Data
    # -----------------------------

    if (is.infinite(n_sample)) {
      message("No sampling applied. Using full dataset.")
      sampled_dt <- combined_dt
    } else if (is.numeric(n_sample) & length(n_sample) == 1 & n_sample > 0) {
      sampled_dt <- combined_dt[, .SD[sample(.N, min(.N, n_sample))], by = population]
      message(paste("Sampling applied. Total samples after sampling:", nrow(sampled_dt)))

      # Check for duplicates after sampling
      sampled_dt <- check_and_handle_duplicates(sampled_dt, paste("Sampling for trait:", trait, "and region:", region))
    } else {
      stop("n_sample must be a single positive number or Inf.")
    }

    # -----------------------------
    # 10. Exclude MHC Region if Required
    # -----------------------------

    if (!include_MHC) {
      pre_filter_rows <- nrow(sampled_dt)
      sampled_dt <- sampled_dt[!(chr == 6 & cg >= 28000000 & cg <= 34000000)]
      post_filter_rows <- nrow(sampled_dt)
      message(paste("Excluded MHC region from the dataset. Rows before:", pre_filter_rows, "Rows after:", post_filter_rows))

      # Check for duplicates after filtering
      sampled_dt <- check_and_handle_duplicates(sampled_dt, paste("Excluding MHC for trait:", trait, "and region:", region))
    }

    # -----------------------------
    # 11. Early Extraction of Mismatched z-stats
    # -----------------------------

    # Pivot data to wide format using data.table's dcast with specified aggregate function
    # Using mean as the aggregate function
    wide_dt <- dcast(sampled_dt, chr + cg + trait + region ~ population, 
                    value.var = c("z", "p", "n", "cor", "mse"), 
                    fun.aggregate = mean, fill = NA)

    message(paste("Pivoted data to wide format for trait:", trait, "and region:", region))

    # Check for duplicates after pivoting
    wide_dt <- check_and_handle_duplicates(wide_dt, paste("Pivoting data for trait:", trait, "and region:", region))

    # Check if both populations are present
    if (!all(c("z_AA", "z_EA") %in% names(wide_dt))) {
      message(paste("Both populations not present for trait", trait, "and region", region))
    } else {
      # Filter rows where z-stat signs mismatch
      mismatched_z_df <- wide_dt[(z_AA * z_EA < 0), .(chr, cg, trait, region,
                                                      AA_z = z_AA, EA_z = z_EA,
                                                      p = p_AA, n = n_AA, cor = cor_AA, mse = mse_AA)]

      # Add population columns
      mismatched_z_df[, population1 := "AA"]
      mismatched_z_df[, population2 := "EA"]

      if (nrow(mismatched_z_df) > 0) {
        # Append to the all mismatched list
        all_mismatched_z_list[[length(all_mismatched_z_list) + 1]] <- mismatched_z_df

        # Save mismatched z-stats table for the specific trait and region with '-v2-dupfilter' suffix
        filename <- paste0("mismatched_zstats_AA_vs_EA_", 
                           gsub(" ", "_", tolower(trait)),
                           "_", region, 
                           "_fulldata-v2-dupfilter.csv")
        filepath <- file.path(output_dir, filename)
        fwrite(mismatched_z_df, filepath)
        message(paste("Saved mismatched z-stat table for trait", trait, "and region", region, ":", filename))

        # Check for duplicates in mismatched_z_df
        mismatched_z_df <- check_and_handle_duplicates(mismatched_z_df, paste("Mismatched z-stats for trait:", trait, "and region:", region))

      } else {
        message(paste("No mismatched z-stat signs found between AA and EA populations for trait", trait, "and region", region))
      }

      # -----------------------------
      # 12. Collect Data for Scatter Plots
      # -----------------------------

      # Collect data for scatter plots
      scatter_data <- wide_dt[, .(chr, cg, trait, region, AA_z = z_AA, EA_z = z_EA)]
      scatter_data <- scatter_data[!is.na(AA_z) & !is.na(EA_z)]

      # Check for duplicates in scatter_data
      scatter_data <- check_and_handle_duplicates(scatter_data, paste("Collecting scatter data for trait:", trait, "and region:", region))

      if (nrow(scatter_data) > 0) {
        all_scatter_data_list[[length(all_scatter_data_list) + 1]] <- scatter_data
      }
    }

    # -----------------------------
    # 13. Generate Plots for Current Trait and Region
    # -----------------------------

    # Sample counts per group
    sample_counts <- sampled_dt[, .N, by = population]
    subtitle_text_cor <- paste("Sample counts per group:\n",
                               paste(sample_counts$population, "=", sample_counts$N, sep="", collapse = "; "))

    # Full title with MHC information
    full_title_cor <- paste0("Boxplot of Correlation by Population (", trait, " - ", region, ")\n",
                             ifelse(include_MHC, "MHC Region Included", "MHC Region Excluded"))

    # Boxplot for 'cor'
    plt_cor <- ggplot(sampled_dt, aes(x = population, y = cor, fill = population)) +
      geom_boxplot(alpha = 0.7) +
      theme_minimal() +
      labs(title = full_title_cor,
           subtitle = subtitle_text_cor,
           x = "Population",
           y = "Correlation",
           fill = "Population") +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )

    # Save the boxplot with '-v2-dupfilter' suffix
    filename_cor <- paste0("boxplot_cor_", 
                           gsub(" ", "_", tolower(trait)), 
                           "_", 
                           region, 
                           "_fulldata-v2-dupfilter.png")
    filepath_cor <- file.path(output_dir, filename_cor)
    ggsave(filename = filepath_cor, plot = plt_cor, width = 8, height = 6, dpi = 300)
    message(paste("Saved boxplot:", filename_cor))
    
    # Also save without the suffix for compatibility
    filename_cor_compat <- paste0("boxplot_cor_", 
                                 gsub(" ", "_", tolower(trait)), 
                                 "_", 
                                 region, 
                                 "_fulldata.png")
    filepath_cor_compat <- file.path(output_dir, filename_cor_compat)
    ggsave(filename = filepath_cor_compat, plot = plt_cor, width = 8, height = 6, dpi = 300)
    message(paste("Saved boxplot (compatibility version):", filename_cor_compat))

    # -----------------------------
    # 14. Generate Scatter Plots for Current Trait and Region
    # -----------------------------

    if (exists("scatter_data") && nrow(scatter_data) > 0) {
      # Calculate the absolute difference
      scatter_data$diff <- abs(scatter_data$AA_z - scatter_data$EA_z)

      # Select the top 6 points with the highest differences
      top_n <- 6
      top_diff <- scatter_data[order(-diff)][1:min(top_n, nrow(scatter_data)), ]

      # Create scatter plot using ggplot2
      plt_scatter <- ggplot(scatter_data, aes(x = AA_z, y = EA_z)) +
        geom_point(alpha = 0.6, color = "blue") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        geom_text_repel(data = top_diff, aes(label = paste0("chr", chr, ":", cg)), size = 3) +
        theme_minimal() +
        labs(title = paste0("Comparison of z-stats: AA vs EA in ", trait, " (", region, ")"),
             x = "AA z-stat",
             y = "EA z-stat") +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )

      # Save the scatter plot with '-v2-dupfilter' suffix
      filename_scatter <- paste0("scatter_zstats_AA_vs_EA_", 
                                 gsub(" ", "_", tolower(trait)), 
                                 "_", 
                                 region, 
                                 "_fulldata-v2-dupfilter.png")
      filepath_scatter <- file.path(output_dir, filename_scatter)
      ggsave(filename = filepath_scatter, plot = plt_scatter, width = 8, height = 6, dpi = 300)
      message(paste("Saved scatter plot:", filename_scatter))
      
      # Also save without the suffix for compatibility
      filename_scatter_compat <- paste0("scatter_zstats_AA_vs_EA_", 
                                       gsub(" ", "_", tolower(trait)), 
                                       "_", 
                                       region, 
                                       "_fulldata.png")
      filepath_scatter_compat <- file.path(output_dir, filename_scatter_compat)
      ggsave(filename = filepath_scatter_compat, plot = plt_scatter, width = 8, height = 6, dpi = 300)
      message(paste("Saved scatter plot (compatibility version):", filename_scatter_compat))
    } else {
      message(paste("No data available for scatter plot in trait", trait, "and region", region))
    }

    # -----------------------------
    # 15. Clean Up Memory
    # -----------------------------

    # Remove datasets to free up memory
    rm(combined_dt, sampled_dt, wide_dt, scatter_data)
    gc()
  }
}

# -----------------------------
# 16. Save the Combined Big File for EA vs AA
# -----------------------------

if (length(all_mismatched_z_list) > 0) {
  big_table <- rbindlist(all_mismatched_z_list, fill = TRUE)
  
  # Define the filename based on MHC inclusion
  sampling_label <- ifelse(is.infinite(n_sample), "fulldata", paste0("mini", n_sample))
  if (include_MHC) {
    sampling_label <- paste0(sampling_label, "_including_MHC")
  } else {
    sampling_label <- paste0(sampling_label, "_excluding_MHC")
  }
  
  # Save with suffix
  filename_big_table <- paste0("mismatched_zstats_AA_vs_EA_big_table_", 
                               sampling_label, 
                               "-v2-dupfilter.csv")
  filepath_big_table <- file.path(output_dir, filename_big_table)
  fwrite(big_table, filepath_big_table)
  message(paste("Saved big table of mismatched z-stats:", filename_big_table))
  
  # Also save without the suffix for compatibility
  filename_big_table_compat <- paste0("mismatched_zstats_AA_vs_EA_big_table_", 
                                     sampling_label, 
                                     ".csv")
  filepath_big_table_compat <- file.path(output_dir, filename_big_table_compat)
  fwrite(big_table, filepath_big_table_compat)
  message(paste("Saved big table of mismatched z-stats (compatibility version):", filename_big_table_compat))
} else {
  message("No mismatched z-stat signs found between AA and EA populations across all traits and regions.")
}

# -----------------------------
# 17. Generate Overall Scatter Plot
# -----------------------------

if (length(all_scatter_data_list) > 0) {
  overall_scatter_data <- rbindlist(all_scatter_data_list, fill = TRUE)

  # Check for duplicates in overall_scatter_data
  overall_scatter_data <- check_and_handle_duplicates(overall_scatter_data, "Combining all scatter data for overall scatter plot")

  # Calculate the absolute difference
  overall_scatter_data$diff <- abs(overall_scatter_data$AA_z - overall_scatter_data$EA_z)

  # Select the top 10 points with the highest differences
  top_n <- 10
  top_diff_overall <- overall_scatter_data[order(-diff)][1:min(top_n, nrow(overall_scatter_data)), ]

  # Create scatter plot using ggplot2
  plt_overall_scatter <- ggplot(overall_scatter_data, aes(x = AA_z, y = EA_z)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text_repel(data = top_diff_overall, 
                    aes(label = paste0(trait, ":", region, "\nchr", chr, ":", cg)), 
                    size = 3) +
    theme_minimal() +
    labs(title = "Comparison of z-stats: AA vs EA across All Traits and Regions",
         x = "AA z-stat",
         y = "EA z-stat") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

  # Save with suffix
  filename_overall_scatter <- paste0("scatter_zstats_AA_vs_EA_overall_", 
                                     sampling_label, 
                                     "-v2-dupfilter.png")
  filepath_overall_scatter <- file.path(output_dir, filename_overall_scatter)
  ggsave(filename = filepath_overall_scatter, plot = plt_overall_scatter, width = 10, height = 8, dpi = 300)
  message(paste("Saved overall scatter plot:", filename_overall_scatter))
  
  # Also save without the suffix for compatibility
  filename_overall_scatter_compat <- paste0("scatter_zstats_AA_vs_EA_overall_", 
                                           sampling_label, 
                                           ".png")
  filepath_overall_scatter_compat <- file.path(output_dir, filename_overall_scatter_compat)
  ggsave(filename = filepath_overall_scatter_compat, plot = plt_overall_scatter, width = 10, height = 8, dpi = 300)
  message(paste("Saved overall scatter plot (compatibility version):", filename_overall_scatter_compat))
} else {
  message("No data available for overall scatter plot.")
}

message("\nScript execution completed successfully!")
