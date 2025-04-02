#!/usr/bin/env Rscript
# verify_plot_recreation.R - Script to verify all original plots have been recreated

# Load necessary libraries
library(data.table)

# Print header
cat("========================================================\n")
cat("             PLOT RECREATION VERIFICATION TOOL          \n")
cat("========================================================\n\n")

# Base directories
scripts_base_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts"
output_base_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/plot_outputs"

# Function to find all image files recursively
find_all_images <- function(base_dir, excluded_dirs = NULL) {
  # Common image extensions
  img_extensions <- "\\.(png|jpg|jpeg|pdf|svg|eps|tiff|bmp)$"
  
  # Build exclusion pattern
  exclusion_pattern <- NULL
  if (!is.null(excluded_dirs) && length(excluded_dirs) > 0) {
    exclusion_pattern <- paste0("(", paste(excluded_dirs, collapse = "|"), ")")
  }
  
  # Command to find all image files
  if (is.null(exclusion_pattern)) {
    cmd <- paste0("find ", base_dir, " -type f | grep -E '", img_extensions, "'")
  } else {
    cmd <- paste0("find ", base_dir, " -type f | grep -v '", exclusion_pattern, "' | grep -E '", img_extensions, "'")
  }
  
  # Execute command and get results
  result <- tryCatch({
    system(cmd, intern = TRUE)
  }, error = function(e) {
    cat(paste("Error executing find command:", e$message, "\n"))
    return(character(0))
  }, warning = function(w) {
    # If grep doesn't find anything, it returns status 1, which causes a warning
    # We'll just return an empty character vector
    return(character(0))
  })
  
  return(result)
}

# Get all output directories
output_dirs <- list.dirs(output_base_dir, recursive = FALSE)
if (length(output_dirs) == 0) {
  cat("No output directories found. Creating the expected directory structure...\n")
  
  # Create the expected directory structure
  dir_names <- c("dnam_summary", "dnam_viz", "general", "manhattan_qq", "parameters", "r_stats")
  for (dir in dir_names) {
    dir.create(file.path(output_base_dir, dir), showWarnings = FALSE, recursive = TRUE)
  }
  
  output_dirs <- list.dirs(output_base_dir, recursive = FALSE)
}

cat(paste("Found", length(output_dirs), "output directories:\n"))
for (dir in output_dirs) {
  cat(paste(" -", basename(dir), "\n"))
}
cat("\n")

# Find all image files in the original scripts directory, excluding the polish_plots directory
cat("Scanning for original image files...\n")
original_images <- find_all_images(scripts_base_dir, c("polish_plots"))

# Check if we got any results
if (length(original_images) == 0) {
  cat("No original image files found. This is unusual and might indicate a problem with the search.\n")
  cat("Please verify the scripts directory contains image files.\n")
  quit(save = "no", status = 1)
}

cat(paste("Found", length(original_images), "original image files.\n\n"))

# Find all recreated images in the output directories
cat("Scanning for recreated image files...\n")
recreated_images <- find_all_images(output_base_dir)

# Check if we got any results
if (length(recreated_images) == 0) {
  cat("No recreated image files found. This suggests the fast scripts haven't been run yet.\n")
  cat("Please run the fast scripts to generate plots before verification.\n")
  
  # Group original images by directory
  original_by_dir <- data.table(
    path = original_images,
    dir = basename(dirname(original_images)),
    filename = basename(original_images)
  )

  # Count by directory
  cat("\nOriginal images by directory:\n")
  dir_counts <- original_by_dir[, .N, by = dir][order(-N)]
  print(dir_counts)
  
  # Generate list of original files
  original_files_list <- file.path(output_base_dir, "original_files.txt")
  writeLines(original_images, original_files_list)
  cat(paste("\nList of original files saved to:", original_files_list, "\n"))
  
  cat("\nOVERALL STATUS: CRITICAL (No plots recreated yet)\n")
  cat("Please run the fast scripts before verification.\n")
  quit(save = "no", status = 2)
} else {
  cat(paste("Found", length(recreated_images), "recreated image files.\n\n"))

  # Group original images by directory
  original_by_dir <- data.table(
    path = original_images,
    dir = basename(dirname(original_images)),
    filename = basename(original_images)
  )

  # Count by directory
  cat("Original images by directory:\n")
  dir_counts <- original_by_dir[, .N, by = dir][order(-N)]
  print(dir_counts)
  cat("\n")

  # Create mapping between original directories and output directories
  dir_mapping <- list(
    # Main plot types to recreate
    "19-OUT_plots" = "19_plots_fast",
    "41-OUT_plots" = "41_plots_fast",
    "43-OUT" = "43_plots_fast",
    "47-OUT" = "47_plots_fast",
    "50-OUT" = "50_plots_fast",
    "DNAm_violin_plots" = "dnam_violin_plots_fast",
    
    # Other possible directories
    "plot_outputs" = "manhattan_qq_fast",
    "manhattan_qq" = "manhattan_qq_fast",
    "manhattan_qq_fast" = "manhattan_qq_fast",
    "OLD" = c("19_plots_fast", "41_plots_fast", "43_plots_fast", "47_plots_fast", "50_plots_fast", "manhattan_qq_fast"),
    
    # Catch-all for plots in main directory
    "." = c("19_plots_fast", "41_plots_fast", "43_plots_fast", "47_plots_fast", "50_plots_fast", "manhattan_qq_fast")
  )

  # Add any missing directories from the original images
  for (dir in unique(original_by_dir$dir)) {
    if (!(dir %in% names(dir_mapping))) {
      # If no mapping exists, map to general
      dir_mapping[[dir]] <- "general"
    }
  }

  # Check recreation status for each directory
  cat("Verification by directory:\n")
  cat("---------------------------\n")

  # Create a dataframe to store verification results
  verification_results <- data.table(
    original_dir = character(),
    original_count = integer(),
    output_dir = character(),
    recreated_count = integer(),
    coverage = numeric(),
    status = character()
  )

  # Keep track of unmapped files
  unmapped_files <- character(0)

  # Check each original directory
  for (dir in unique(original_by_dir$dir)) {
    # Count original images in this directory
    original_count <- original_by_dir[dir == dir, .N]
    original_files <- original_by_dir[dir == dir, filename]
    
    # Get mapping
    mapped_dirs <- if (dir %in% names(dir_mapping)) dir_mapping[[dir]] else NULL
    
    if (is.null(mapped_dirs)) {
      # No mapping found
      cat(sprintf("%-20s: %4d originals, NO MAPPING DEFINED\n", dir, original_count))
      
      verification_results <- rbind(verification_results, data.table(
        original_dir = dir,
        original_count = original_count,
        output_dir = NA_character_,
        recreated_count = 0,
        coverage = 0,
        status = "NO MAPPING"
      ))
      
      # Add files to unmapped list
      unmapped_files <- c(unmapped_files, original_by_dir[dir == dir, path])
    } else {
      # Check each mapped directory
      if (is.character(mapped_dirs) && length(mapped_dirs) == 1) {
        mapped_dirs <- c(mapped_dirs)
      }
      
      for (mapped_dir in mapped_dirs) {
        full_mapped_dir <- file.path(output_base_dir, mapped_dir)
        
        # Count recreated images
        if (dir.exists(full_mapped_dir)) {
          recreated_files <- list.files(full_mapped_dir, pattern = "\\.(png|jpg|jpeg|pdf|svg|eps|tiff|bmp)$", 
                                      recursive = TRUE, full.names = TRUE)
          
          recreated_count <- length(recreated_files)
          coverage <- min(100, round(recreated_count / original_count * 100, 1))
          
          status <- if (coverage >= 90) "GOOD" else if (coverage >= 50) "PARTIAL" else "LOW"
          
          cat(sprintf("%-20s: %4d originals, %4d recreated in %s (%0.1f%%) - %s\n", 
                      dir, original_count, recreated_count, mapped_dir, coverage, status))
          
          verification_results <- rbind(verification_results, data.table(
            original_dir = dir,
            original_count = original_count,
            output_dir = mapped_dir,
            recreated_count = recreated_count,
            coverage = coverage,
            status = status
          ))
        } else {
          cat(sprintf("%-20s: %4d originals, OUTPUT DIRECTORY %s NOT FOUND\n", dir, original_count, mapped_dir))
          
          verification_results <- rbind(verification_results, data.table(
            original_dir = dir,
            original_count = original_count,
            output_dir = mapped_dir,
            recreated_count = 0,
            coverage = 0,
            status = "DIR NOT FOUND"
          ))
        }
      }
    }
  }

  # Print summary statistics
  cat("\nSummary:\n")
  cat("--------\n")

  total_original <- sum(dir_counts$N)
  total_recreated <- length(recreated_images)
  total_coverage <- if (total_recreated > 0) round(total_recreated / total_original * 100, 1) else 0

  cat(sprintf("Total original images: %d\n", total_original))
  cat(sprintf("Total recreated images: %d\n", total_recreated))
  cat(sprintf("Overall coverage: %.1f%%\n", total_coverage))

  # Check which directories have the poorest coverage
  cat("\nDirectories needing attention:\n")
  low_coverage <- verification_results[coverage < 50 & !is.na(output_dir)]
  if (nrow(low_coverage) > 0) {
    print(low_coverage[order(coverage)])
  } else {
    cat("All mapped directories have at least 50% coverage\n")
  }

  # Check which directories have no mapping defined
  cat("\nDirectories with no mapping defined:\n")
  no_mapping <- verification_results[status == "NO MAPPING"]
  if (nrow(no_mapping) > 0) {
    print(no_mapping[order(-original_count)])
  } else {
    cat("All directories have mappings defined\n")
  }

  # Generate comprehensive report file
  report_file <- file.path(output_base_dir, "verification_report.csv")
  fwrite(verification_results, report_file)
  cat(paste("\nDetailed report saved to:", report_file, "\n"))

  # Generate list of unmapped files
  if (length(unmapped_files) > 0) {
    unmapped_files_list <- file.path(output_base_dir, "unmapped_files.txt")
    writeLines(unmapped_files, unmapped_files_list)
    cat(paste("List of unmapped files saved to:", unmapped_files_list, "\n"))
  }

  # Generate list of recreated files
  if (length(recreated_images) > 0) {
    recreated_files_list <- file.path(output_base_dir, "recreated_files.txt")
    writeLines(recreated_images, recreated_files_list)
    cat(paste("List of recreated files saved to:", recreated_files_list, "\n"))
  }

  # Generate list of original files
  original_files_list <- file.path(output_base_dir, "original_files.txt")
  writeLines(original_images, original_files_list)
  cat(paste("List of original files saved to:", original_files_list, "\n"))

  cat("\nVerification complete.\n")

  # Exit with status code
  if (total_coverage >= 90) {
    cat("OVERALL STATUS: GOOD (≥90% coverage)\n")
    quit(save = "no", status = 0)
  } else if (total_coverage >= 75) {
    cat("OVERALL STATUS: ACCEPTABLE (≥75% coverage)\n")
    quit(save = "no", status = 0)
  } else if (total_recreated > 0) {
    cat("OVERALL STATUS: NEEDS IMPROVEMENT (<75% coverage)\n")
    quit(save = "no", status = 1)
  } else {
    cat("OVERALL STATUS: CRITICAL (No plots recreated yet)\n")
    cat("Please run the fast scripts before verification.\n")
    quit(save = "no", status = 2)
  }
}
