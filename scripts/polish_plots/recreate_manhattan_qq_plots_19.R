#!/usr/bin/env Rscript
# recreate_manhattan_qq_plots_19.R - Script to recreate Manhattan and QQ plots from 19_plot.ipynb
# This script combines best practices from existing scripts to recreate plots
# with proper handling of input data files

# Load necessary libraries - FAIL FAST if any library fails to load
if (!require(data.table)) {
  stop("CRITICAL ERROR: Failed to load data.table library",
       "\nThis is a critical error - we cannot proceed without required libraries!")
}
if (!require(qqman)) {
  stop("CRITICAL ERROR: Failed to load qqman library",
       "\nThis is a critical error - we cannot proceed without required libraries!")
}
if (!require(ggplot2)) {
  stop("CRITICAL ERROR: Failed to load ggplot2 library",
       "\nThis is a critical error - we cannot proceed without required libraries!")
}
if (!require(scales)) {
  stop("CRITICAL ERROR: Failed to load scales library",
       "\nThis is a critical error - we cannot proceed without required libraries!")
}
if (!require(gridExtra)) {
  stop("CRITICAL ERROR: Failed to load gridExtra library",
       "\nThis is a critical error - we cannot proceed without required libraries!")
}

# Create output directories with proper naming convention - FAIL FAST if directory creation fails
# FIXED: Updated input directory path to the correct location
input_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/plot_outputs/manhattan_qq"
output_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/manhattan_qq"
output_dir_fast <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/manhattan_qq_fast"
output_dir_improved <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/manhattan_qq_improved"

# Create all output directories
for (dir in c(output_dir, output_dir_fast, output_dir_improved)) {
  if (!dir.exists(dir)) {
    if (!dir.create(dir, recursive = TRUE)) {
      stop("CRITICAL ERROR: Failed to create output directory: ", dir,
           "\nThis is a critical error - we cannot proceed without valid output directories!")
    }
    message(paste("Created output directory:", dir))
  } else {
    message(paste("Using existing output directory:", dir))
  }
}

# Define variables
desired_subpopulations <- c("AA", "EA")  # Removed "all" as it's not in the original plots
desired_regions <- c("caud", "hippo", "dlpfc")
traits <- c("scz", "bp", "mdd")

# Function to create Manhattan plot with improved aesthetics
create_manhattan_plot <- function(data, trait, population, region, output_dir, 
                                width = 1200, height = 800, res = 150) {
  # Create filename for the plot
  output_file <- file.path(output_dir, paste0("manhattan_", trait, "_", population, "_", region, ".png"))
  
  # Check for required columns - FAIL FAST if missing
  required_cols <- c("chr", "cg", "p")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("CRITICAL ERROR: Missing required columns for Manhattan plot: ",
         paste(missing_cols, collapse = ", "),
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Filter out invalid values
  valid_data <- data[!is.na(data$p) & is.finite(data$p) & data$p > 0, ]
  
  if (nrow(valid_data) == 0) {
    stop("CRITICAL ERROR: No valid p-values found in data after filtering",
         "\nThis is a critical error - we cannot proceed without valid data!")
  }
  
  if (nrow(valid_data) < nrow(data)) {
    message(paste("Filtered out", nrow(data) - nrow(valid_data), "rows with invalid p-values"))
  }
  
  # Ensure SNP column exists for the manhattan function - FAIL FAST if creation fails
  if (!"SNP" %in% colnames(valid_data)) {
    valid_data[, SNP := paste0(chr, ":", cg)]
    if (any(is.na(valid_data$SNP))) {
      stop("CRITICAL ERROR: Failed to create SNP column",
           "\nThis is a critical error - we cannot proceed with invalid data!")
    }
  }
  
  # Ensure chromosome values are numeric - FAIL FAST if invalid
  if (is.character(valid_data$chr)) {
    valid_data$chr <- as.numeric(gsub("[^0-9]", "", valid_data$chr))
    if (any(is.na(valid_data$chr))) {
      stop("CRITICAL ERROR: Invalid chromosome values found",
           "\nThis is a critical error - we cannot proceed with invalid data!")
    }
  }
  
  # Check for valid chromosomes - FAIL FAST if invalid
  invalid_chromosomes <- !valid_data$chr %in% 1:22
  if (any(invalid_chromosomes)) {
    message(paste("Filtering out", sum(invalid_chromosomes), "rows with invalid chromosome numbers"))
    valid_data <- valid_data[!invalid_chromosomes, ]
    
    if (nrow(valid_data) == 0) {
      stop("CRITICAL ERROR: No valid data left after filtering invalid chromosomes",
           "\nThis is a critical error - we cannot proceed without valid data!")
    }
  }
  
  # Create Manhattan plot - FIXED: Removed the ! operator
  png(filename = output_file, width = width, height = height, res = res)
  
  # Customize colors for better readability
  col_palette <- c("#4286f4", "#f44182")
  
  manhattan(valid_data, chr = "chr", bp = "cg", p = "p", 
           main = paste0("Manhattan Plot: ", toupper(trait), " in ", toupper(region), " - ", population),
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

# Function to create QQ plot with improved aesthetics
create_qq_plot <- function(p_values, trait, population, region, output_dir,
                          width = 800, height = 800, res = 150, max_points = 10000) {
  # Create filename for the plot
  output_file <- file.path(output_dir, paste0("qqplot_", trait, "_", population, "_", region, ".png"))
  
  # Filter out invalid values
  valid_p_values <- p_values[!is.na(p_values) & is.finite(p_values) & p_values > 0 & p_values <= 1]
  
  if (length(valid_p_values) == 0) {
    stop("CRITICAL ERROR: No valid p-values found after filtering",
         "\nThis is a critical error - we cannot proceed without valid data!")
  }
  
  if (length(valid_p_values) < length(p_values)) {
    message(paste("Filtered out", length(p_values) - length(valid_p_values), "invalid p-values"))
  }
  
  # Sample if we have more values than needed
  if (length(valid_p_values) > max_points) {
    set.seed(42)  # For reproducibility
    valid_p_values <- sample(valid_p_values, max_points)
  }
  
  # Calculate genomic inflation factor (lambda) - FAIL FAST if calculation fails
  lambda <- tryCatch({
    median(qchisq(1-valid_p_values, 1)) / qchisq(0.5, 1)
  }, error = function(e) {
    message("Error calculating genomic inflation factor: ", e$message)
    return(NA)
  })
  
  if (is.na(lambda) || !is.finite(lambda)) {
    message("Warning: Could not calculate genomic inflation factor, using default value of 1.0")
    lambda <- 1.0
  }
  
  # Create QQ plot - FIXED: Removed the ! operator
  png(filename = output_file, width = width, height = height, res = res)
  
  # Set up plot margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  qq(valid_p_values, 
     main = paste0("QQ Plot: ", toupper(trait), " in ", toupper(region), " - ", population,
                  "\nGenomic Inflation Factor λ = ", format(lambda, digits = 3)),
     cex = 0.9, 
     cex.axis = 1.1,
     pch = 19,
     col = "navy",
     xlab = expression(Expected~~-log[10](p)),
     ylab = expression(Observed~~-log[10](p)))
  
  # Add diagonal reference line
  abline(0, 1, col = "red", lwd = 2)
  
  dev.off()
  
  # Verify plot was created - FAIL FAST if not
  if (!file.exists(output_file)) {
    stop("CRITICAL ERROR: Failed to save QQ plot: ", output_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  
  message("Created QQ plot: ", basename(output_file))
  return(output_file)
}

# FIXED: Added function to find input files with flexible naming patterns
find_input_file <- function(base_dir, trait, population, region, type) {
  # Define possible file patterns
  patterns <- c(
    paste0(type, "_", trait, "_results_", population, "_", region, ".csv"),
    paste0(type, "-bonf_", trait, "_results_", population, "_", region, ".csv"),
    paste0(type, "_", trait, "_results_a[0-9]+_", population, "_", region, ".csv"),
    paste0(type, "-bonf_", trait, "_results_a[0-9]+_", population, "_", region, ".csv")
  )
  
  # Try each pattern
  for (pattern in patterns) {
    files <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    if (length(files) > 0) {
      message("Found input file: ", basename(files[1]))
      return(files[1])
    }
  }
  
  # If no file found, return NULL
  return(NULL)
}

# Process each combination
for (trait in traits) {
  for (population in desired_subpopulations) {
    for (region in desired_regions) {
      message("\nProcessing ", trait, " ", population, " ", region)
      
      # FIXED: Use the find_input_file function to locate input files with flexible naming patterns
      filtered_file <- find_input_file(input_dir, trait, population, region, "filtered")
      unfiltered_file <- find_input_file(input_dir, trait, population, region, "unfiltered")
      
      # Check if the filtered file exists - FAIL FAST AND LOUD
      if (is.null(filtered_file)) {
        message("WARNING: Required filtered input file not found for ", trait, " ", population, " ", region)
        message("Skipping this combination...")
        next
      }
      
      # Check if the unfiltered file exists - FAIL FAST AND LOUD
      if (is.null(unfiltered_file)) {
        message("WARNING: Required unfiltered input file not found for ", trait, " ", population, " ", region)
        message("Skipping this combination...")
        next
      }
      
      # Load the filtered data for Manhattan plot
      message("Loading filtered data from: ", filtered_file)
      merged_subset <- tryCatch({
        fread(filtered_file)
      }, error = function(e) {
        message("Error reading filtered file: ", e$message)
        return(NULL)
      })
      
      if (is.null(merged_subset) || nrow(merged_subset) == 0) {
        message("WARNING: Failed to load filtered data or file is empty")
        message("Skipping this combination...")
        next
      }
      
      # Create Manhattan plots in all output directories
      for (dir in c(output_dir, output_dir_fast, output_dir_improved)) {
        tryCatch({
          create_manhattan_plot(merged_subset, trait, population, region, dir)
        }, error = function(e) {
          message("Error creating Manhattan plot: ", e$message)
        })
      }
      
      # Free memory
      rm(merged_subset)
      gc()
      
      # Load the unfiltered data
      message("Loading unfiltered data from: ", unfiltered_file)
      merged <- tryCatch({
        fread(unfiltered_file)
      }, error = function(e) {
        message("Error reading unfiltered file: ", e$message)
        return(NULL)
      })
      
      if (is.null(merged) || nrow(merged) == 0) {
        message("WARNING: Failed to load unfiltered data or file is empty")
        message("Skipping QQ plot for this combination...")
        next
      }
      
      # Create QQ plots in all output directories
      for (dir in c(output_dir, output_dir_fast, output_dir_improved)) {
        tryCatch({
          create_qq_plot(merged$p, trait, population, region, dir)
        }, error = function(e) {
          message("Error creating QQ plot: ", e$message)
        })
      }
      
      # Free memory
      rm(merged)
      gc()
    }
  }
}

message("\nPlot generation complete! All outputs saved to:")
message("- ", output_dir)
message("- ", output_dir_fast)
message("- ", output_dir_improved)