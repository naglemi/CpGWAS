#!/usr/bin/env Rscript
# verify_plot_completeness.R - Script to verify all plots have been recreated
# This script performs a comprehensive analysis of all original images in the workspace
# and compares them to the recreated plots in polish_plots/plot_outputs

# Load necessary libraries
library(data.table)
library(dplyr)
library(stringr)

# Add timeout functionality
start_time <- Sys.time()
timeout_minutes <- 10

# Function to check if timeout has been reached
check_timeout <- function() {
  current_time <- Sys.time()
  elapsed_minutes <- as.numeric(difftime(current_time, start_time, units = "mins"))
  
  if (elapsed_minutes >= timeout_minutes) {
    cat("\n\n========================================================\n")
    cat("                      TIMEOUT ERROR                      \n")
    cat("========================================================\n\n")
    cat("Script execution exceeded the maximum allowed time of", timeout_minutes, "minutes.\n")
    cat("This script is too slow for interactive use by agents.\n\n")
    cat("Recommendations:\n")
    cat("1. Use the faster verification script: fast_verify_plots.R\n")
    cat("2. Limit the scope of verification to specific directories\n")
    cat("3. Use the verification_results.json file for current status\n\n")
    quit(status = 1)
  }
}

# Print header
cat("========================================================\n")
cat("              PLOT RECREATION VERIFICATION              \n")
cat("========================================================\n\n")
cat("Timeout set to", timeout_minutes, "minutes\n\n")

# Define directories
scripts_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts"
outputs_dir <- file.path(scripts_dir, "plot_outputs")
report_dir <- file.path(scripts_dir, "polish_plots")

# Ensure report directory exists
if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
}

# Function to search for all image files recursively, but with optimizations
find_all_images <- function(base_dir) {
  # List of file extensions to consider image files
  image_extensions <- c("png", "jpg", "jpeg", "pdf", "svg", "tiff", "bmp", "gif")
  
  # Create a pattern to match any of these extensions
  pattern <- paste0("\\.(", paste(image_extensions, collapse = "|"), ")$")
  
  # OPTIMIZATION: Focus only on directories that are likely to contain plots
  # This significantly reduces the search space
  plot_related_dirs <- c(
    "19-OUT_plots", "41-OUT_plots", "43-OUT", "47-OUT", "50-OUT",
    "DNAm_violin_plots", "manhattan_qq", "plot_outputs"
  )
  
  # Find all plot-related directories first
  cat("Finding plot-related directories...\n")
  cmd <- paste0("find ", base_dir, " -type d -name '*plots*' -o -name '*OUT*' | grep -v 'polish_plots/plot_outputs'")
  potential_dirs <- system(cmd, intern = TRUE)
  
  # Add the base directory as well
  potential_dirs <- c(base_dir, potential_dirs)
  
  # Initialize an empty vector to store full paths
  all_paths <- character(0)
  
  # For each potential directory, find image files directly (no recursion)
  cat("Searching for image files in", length(potential_dirs), "directories...\n")
  for (dir_path in potential_dirs) {
    check_timeout()
    
    # Use system command for faster file finding
    cmd <- paste0("find ", dir_path, " -maxdepth 1 -type f -regextype posix-extended -regex '.*\\.(",
                 paste(image_extensions, collapse = "|"), ")$'")
    image_files <- system(cmd, intern = TRUE)
    
    # Add to our collection
    all_paths <- c(all_paths, image_files)
    
    # Progress indicator
    if (length(all_paths) %% 100 == 0) {
      cat("Found", length(all_paths), "images so far...\n")
    }
  }
  
  return(all_paths)
}

# Function to identify plot directories
find_plot_directories <- function() {
  # Check for timeout
  check_timeout()
  
  # Get all directories in the scripts folder
  dirs <- list.dirs(scripts_dir, recursive = TRUE)
  
  # Filter for directories that look like plot directories
  plot_dir_pattern <- "OUT_plots|OUT/|plots"
  plot_dirs <- dirs[grepl(plot_dir_pattern, dirs)]
  
  # Remove polish_plots/plot_outputs from the list
  plot_dirs <- plot_dirs[!grepl("polish_plots/plot_outputs", plot_dirs)]
  
  return(plot_dirs)
}

# Function to find recreated plots
find_recreated_plots <- function() {
  # Check for timeout
  check_timeout()
  
  # Find all plot output directories
  output_dirs <- list.dirs(outputs_dir, recursive = FALSE)
  
  # Initialize an empty vector to store full paths
  recreated_paths <- character(0)
  
  # For each output directory, find all image files
  for (dir in output_dirs) {
    # Check for timeout periodically
    check_timeout()
    
    # List image files in the current directory
    image_extensions <- c("png", "jpg", "jpeg", "pdf", "svg", "tiff", "bmp", "gif")
    pattern <- paste0("\\.(", paste(image_extensions, collapse = "|"), ")$")
    
    files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
    recreated_paths <- c(recreated_paths, files)
  }
  
  return(recreated_paths)
}

# Function to extract base name from a file path (for comparison)
extract_base_name <- function(file_path) {
  # Get the file name without path and extension
  base_name <- basename(file_path)
  base_name <- sub("\\.[^.]+$", "", base_name)  # Remove extension
  
  # Remove common prefixes and suffixes for comparison
  base_name <- sub("^(fast_|19-OUT_|41-OUT_|43-OUT_|47-OUT_|50-OUT_)", "", base_name)
  
  # Some patterns are harder to normalize, so we'll use regex patterns
  # Remove numbered prefixes like "19.1_" or "47.3_"
  base_name <- sub("^\\d+(\\.\\d+)?_", "", base_name)
  
  # Handle special cases and transformations
  if (grepl("manhattan_", base_name)) {
    # Extract trait, population, and region for manhattan plots
    parts <- str_match(base_name, "manhattan_([^_]+)_([^_]+)_(.+)")
    if (!is.na(parts[1,1])) {
      # Standardize order to trait_pop_region for comparison
      trait <- parts[1,2]
      pop <- parts[1,3]
      region <- parts[1,4]
      base_name <- paste(trait, pop, region, sep = "_")
    }
  }
  
  return(tolower(base_name))  # Convert to lowercase for case-insensitive comparison
}

# Function to categorize an image based on directory or filename
categorize_image <- function(file_path) {
  if (grepl("19-OUT_plots", file_path)) {
    return("19_plots")
  } else if (grepl("41-OUT_plots|41A_plots", file_path)) {
    return("41_plots")
  } else if (grepl("43-OUT|43_annotate", file_path)) {
    return("43_plots")
  } else if (grepl("47-OUT|47_merge", file_path)) {
    return("47_plots")
  } else if (grepl("50-OUT|50_make", file_path)) {
    return("50_plots")
  } else if (grepl("DNAm_violin|violin_plot", file_path)) {
    return("dnam_violin_plots")
  } else if (grepl("manhattan|qq_plot", file_path)) {
    return("manhattan_qq")
  } else {
    # Try to categorize by filename pattern
    base_name <- tolower(basename(file_path))
    if (grepl("manhattan|qq_plot", base_name)) {
      return("manhattan_qq")
    } else if (grepl("correlation|scatter", base_name)) {
      return("41_plots")
    } else if (grepl("gwas|manhattan", base_name)) {
      return("43_plots")
    } else if (grepl("feature|genomic|nott", base_name)) {
      return("47_plots")
    } else if (grepl("overlap|bed", base_name)) {
      return("50_plots")
    } else if (grepl("methylation|dnam|violin", base_name)) {
      return("dnam_violin_plots")
    } else {
      return("uncategorized")
    }
  }
}

# OPTIMIZATION: Create a lookup table for faster matching
create_recreated_lookup <- function(recreated_paths) {
  # Extract base names for all recreated plots
  recreated_bases <- sapply(recreated_paths, extract_base_name)
  
  # Create a named list for faster lookup
  lookup <- list()
  for (i in 1:length(recreated_bases)) {
    base <- recreated_bases[i]
    lookup[[base]] <- TRUE
  }
  
  # Also store the full list for fuzzy matching if needed
  attr(lookup, "full_list") <- recreated_bases
  attr(lookup, "full_paths") <- recreated_paths
  
  return(lookup)
}

# Function to check if a plot has been recreated (optimized version)
has_been_recreated <- function(original_path, recreated_lookup) {
  # Check for timeout periodically
  if (sample(1:100, 1) == 1) {  # Random sampling to avoid checking too frequently
    check_timeout()
  }
  
  # Extract base name for comparison
  original_base <- extract_base_name(original_path)
  
  # First try exact match (very fast)
  if (!is.null(recreated_lookup[[original_base]])) {
    return(TRUE)
  }
  
  # If no exact match, try fuzzy matching but only on a sample
  # This is a compromise between accuracy and speed
  recreated_bases <- attr(recreated_lookup, "full_list")
  
  # Only check a random sample of 50 plots for fuzzy matching
  # This is much faster and still catches most matches
  sample_size <- min(50, length(recreated_bases))
  if (sample_size > 0) {
    sample_indices <- sample(1:length(recreated_bases), sample_size)
    sample_bases <- recreated_bases[sample_indices]
    
    for (recreated_base in sample_bases) {
      if (grepl(original_base, recreated_base) || grepl(recreated_base, original_base)) {
        return(TRUE)
      }
    }
  }
  
  return(FALSE)
}

# OPTIMIZATION: Faster function to estimate why a plot might be missing
estimate_missing_reason <- function(file_path) {
  # Check for patterns suggesting temporary/intermediate plots
  if (grepl("temp|intermediate|draft|tmp", file_path, ignore.case = TRUE)) {
    return("Likely a temporary or intermediate plot")
  }
  
  # Check for patterns suggesting old versions
  if (grepl("old|backup|archive|deprecated", file_path, ignore.case = TRUE)) {
    return("Likely an old or deprecated version")
  }
  
  # Check if the plot is in a non-standard location
  if (!any(grepl("OUT", dirname(file_path)))) {
    return("Not in a standard plot output directory")
  }
  
  # OPTIMIZATION: Only check file size if we need to (this is expensive)
  # Use a cache to avoid checking the same file multiple times
  if (!exists("file_size_cache")) {
    file_size_cache <<- new.env(hash = TRUE)
  }
  
  # Check if we've already cached this file's size
  if (is.null(file_size_cache[[file_path]])) {
    # Only check if file exists first (faster than file.info on non-existent files)
    if (file.exists(file_path)) {
      file_info <- file.info(file_path)
      file_size_cache[[file_path]] <- file_info$size
    } else {
      file_size_cache[[file_path]] <- -1  # Mark as non-existent
    }
  }
  
  # Check for unusual file size (very small might be empty/template)
  if (file_size_cache[[file_path]] > 0 && file_size_cache[[file_path]] < 5000) {  # Less than 5KB
    return("Very small file size, might be a template or empty plot")
  }
  
  # Default reason if we can't determine
  return("Unknown - May require additional script for recreation")
}

# Start the verification process
cat("Scanning for all image files in scripts directory...\n")
check_timeout()
all_images <- find_all_images(scripts_dir)
cat("Found", length(all_images), "image files\n\n")

cat("Identifying plot directories...\n")
check_timeout()
plot_dirs <- find_plot_directories()
cat("Found", length(plot_dirs), "plot directories\n\n")

cat("Finding recreated plots...\n")
check_timeout()
recreated_plots <- find_recreated_plots()
cat("Found", length(recreated_plots), "recreated plots\n\n")

# OPTIMIZATION: Create lookup table for faster matching
cat("Creating lookup table for faster matching...\n")
check_timeout()
recreated_lookup <- create_recreated_lookup(recreated_plots)
cat("Lookup table created with", length(recreated_lookup), "entries\n\n")

# Analyze the results
cat("Analyzing plot recreation status...\n")
check_timeout()

# OPTIMIZATION: Process in batches to avoid memory issues
cat("Processing images in batches...\n")
batch_size <- 1000
num_batches <- ceiling(length(all_images) / batch_size)
results_list <- list()

for (i in 1:num_batches) {
  check_timeout()
  start_idx <- (i-1) * batch_size + 1
  end_idx <- min(i * batch_size, length(all_images))
  batch_images <- all_images[start_idx:end_idx]
  
  cat("Processing batch", i, "of", num_batches, "(", length(batch_images), "images )\n")
  
  # Create a data frame for this batch
  batch_results <- data.frame(
    original_path = batch_images,
    category = sapply(batch_images, categorize_image),
    base_name = sapply(batch_images, extract_base_name),
    recreated = sapply(batch_images, function(x) has_been_recreated(x, recreated_lookup)),
    stringsAsFactors = FALSE
  )
  
  results_list[[i]] <- batch_results
}

# Combine all batches
results <- do.call(rbind, results_list)

check_timeout()

# Add additional information for missing plots
results$missing_reason <- ifelse(!results$recreated,
                               sapply(results$original_path, estimate_missing_reason),
                               "N/A")

check_timeout()

# Summarize by category
summary_by_category <- results %>%
  group_by(category) %>%
  summarize(
    total_plots = n(),
    recreated_plots = sum(recreated),
    missing_plots = sum(!recreated),
    recreation_percentage = round(100 * sum(recreated) / n(), 1)
  )

check_timeout()

# Generate detailed report of missing plots
missing_plots <- results[!results$recreated, ]

check_timeout()

# Write results to files
write.csv(results, file.path(report_dir, "plot_verification_full_report.csv"), row.names = FALSE)
write.csv(summary_by_category, file.path(report_dir, "plot_verification_summary.csv"), row.names = FALSE)
write.csv(missing_plots, file.path(report_dir, "plot_verification_missing.csv"), row.names = FALSE)

check_timeout()

# Print summary to console
cat("\nSummary by Category:\n")
print(summary_by_category, row.names = FALSE)

cat("\nOverall Statistics:\n")
cat("Total original plots:", nrow(results), "\n")
cat("Successfully recreated:", sum(results$recreated), "\n")
cat("Not recreated:", sum(!results$recreated), "\n")
cat("Overall recreation percentage:", round(100 * sum(results$recreated) / nrow(results), 1), "%\n\n")

cat("Detailed reports saved to:\n")
cat("- Full report:", file.path(report_dir, "plot_verification_full_report.csv"), "\n")
cat("- Summary by category:", file.path(report_dir, "plot_verification_summary.csv"), "\n")
cat("- Missing plots:", file.path(report_dir, "plot_verification_missing.csv"), "\n\n")

cat("Next steps:\n")
cat("1. Review the missing_plots report to identify any critical plots that still need recreation\n")
cat("2. For plots with 'Unknown' reason, determine if they need dedicated script development\n")
cat("3. Plots identified as temporary, old versions, or templates can likely be ignored\n")
cat("4. Update fast scripts as needed to address any gaps in coverage\n\n")

cat("Verification complete!\n")
