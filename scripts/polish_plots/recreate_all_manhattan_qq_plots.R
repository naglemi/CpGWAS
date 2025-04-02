#!/usr/bin/env Rscript
# recreate_all_manhattan_qq_plots.R - Script to recreate ALL Manhattan and QQ plots
# This script combines both 19_plot.ipynb and 37_plot_and_summarize.ipynb sources
# to provide a comprehensive solution for manhattan and qq plots

# Load required libraries
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(data.table)
library(qqman)
library(stringr)
library(dplyr)

# Record start time
start_time <- Sys.time()
cat("Starting manhattan and qq plot recreation at:", format(start_time), "\n")

# Define base paths using absolute paths (as required by mission)
base_scripts_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts"
output_dir <- file.path(base_scripts_dir, "polish_plots", "plot_outputs", "manhattan_qq_fast")

# Create output directory (and parent directories if needed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", output_dir, "\n")

# Define traits, populations, and regions
traits <- c("scz", "bp", "mdd")
populations <- c("AA", "EA", "all")
regions <- c("caud", "dlpfc", "hippo")

# Function to create Manhattan plot with improved aesthetics
create_manhattan_plot <- function(data, trait, population, region, output_dir, 
                                 name_pattern = "default", width = 1200, height = 800) {
  # Create filename for the plot with appropriate naming pattern
  if (name_pattern == "19_format") {
    output_file <- file.path(output_dir, 
      paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_manhattan.png"))
  } else {
    output_file <- file.path(output_dir, 
      paste0("manhattan_", trait, "_", population, "_", region, ".png"))
  }
  
  # Create a copy of the data to avoid modifying the original
  plot_data <- data.table::copy(data)
  
  # Check for required columns and standardize column names
  if ("chr" %in% colnames(plot_data)) {
    plot_data$CHR <- plot_data$chr
  } else if ("mwas_chr" %in% colnames(plot_data)) {
    plot_data$CHR <- plot_data$mwas_chr
  }
  
  if ("bp" %in% colnames(plot_data)) {
    plot_data$BP <- plot_data$bp
  } else if ("mwas_pos" %in% colnames(plot_data)) {
    plot_data$BP <- plot_data$mwas_pos 
  } else if ("pos" %in% colnames(plot_data)) {
    plot_data$BP <- plot_data$pos
  } else if ("cg" %in% colnames(plot_data)) {
    # Use CpG as position if no position column exists
    plot_data$BP <- 1:nrow(plot_data)
  }
  
  if ("p" %in% colnames(plot_data)) {
    plot_data$P <- plot_data$p
  } else if ("pval" %in% colnames(plot_data)) {
    plot_data$P <- plot_data$pval
  } else if ("p_value" %in% colnames(plot_data)) {
    plot_data$P <- plot_data$p_value
  }
  
  # Check for required columns after standardization
  if (!all(c("CHR", "BP", "P") %in% colnames(plot_data))) {
    cat("Error: Required columns (CHR, BP, P) not found in data\n")
    return(NULL)
  }
  
  # Add SNP column if not present (required by manhattan function)
  if (!"SNP" %in% colnames(plot_data)) {
    plot_data$SNP <- paste0("snp", 1:nrow(plot_data))
  }
  
  # Filter out problematic values
  plot_data <- plot_data[!is.na(plot_data$P) & is.finite(plot_data$P) & plot_data$P > 0, ]
  
  # Ensure chromosome values are numeric
  if (is.character(plot_data$CHR)) {
    plot_data$CHR <- as.numeric(gsub("[^0-9]", "", plot_data$CHR))
  }
  
  # Filter out non-numeric chromosomes
  plot_data <- plot_data[!is.na(plot_data$CHR) & plot_data$CHR %in% 1:22, ]
  
  # Ensure all required columns have data
  if (nrow(plot_data) == 0) {
    cat("Error: No valid data points after filtering\n")
    return(NULL)
  }
  
  # Create Manhattan plot
  png(filename = output_file, width = width, height = height, res = 150)
  
  # Customize colors for better readability
  col_palette <- c("#4286f4", "#f44182")
  
  tryCatch({
    manhattan(plot_data, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
             main = paste0("Manhattan Plot: ", toupper(trait), " in ", toupper(region), " - ", population),
             cex = 0.6, cex.axis = 0.8, col = col_palette,
             suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
             chrlabs = as.character(1:22))
    
    dev.off()
    cat("Created Manhattan plot:", basename(output_file), "\n")
    return(output_file)
  }, error = function(e) {
    dev.off()
    cat("Error creating Manhattan plot:", e$message, "\n")
    return(NULL)
  })
}

# Function to create QQ plot with improved aesthetics
create_qq_plot <- function(data, trait, population, region, output_dir,
                          name_pattern = "default", width = 800, height = 800, max_points = 10000) {
  # Create filename for the plot with appropriate naming pattern
  if (name_pattern == "19_format") {
    output_file <- file.path(output_dir, 
      paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_qq.png"))
  } else {
    output_file <- file.path(output_dir, 
      paste0("qqplot_", trait, "_", population, "_", region, ".png"))
  }
  
  # Extract p-values with correct column name
  if ("p" %in% colnames(data)) {
    p_values <- data$p
  } else if ("pval" %in% colnames(data)) {
    p_values <- data$pval
  } else if ("p_value" %in% colnames(data)) {
    p_values <- data$p_value
  } else if ("P" %in% colnames(data)) {
    p_values <- data$P
  } else {
    cat("Error: No p-value column found. Column names:", paste(colnames(data), collapse=", "), "\n")
    return(NULL)
  }
  
  # Clean p-values
  p_values <- p_values[!is.na(p_values) & is.finite(p_values) & p_values > 0 & p_values <= 1]
  
  # Verify we have sufficient p-values
  if (length(p_values) == 0) {
    cat("Error: No valid p-values found\n")
    return(NULL)
  }
  
  # Sample if we have more values than needed
  if (length(p_values) > max_points) {
    set.seed(42)  # For reproducibility
    p_values <- sample(p_values, max_points)
  }
  
  # Calculate genomic inflation factor (lambda)
  lambda <- median(qchisq(1-p_values, 1)) / qchisq(0.5, 1)
  if (is.na(lambda) || !is.finite(lambda)) {
    lambda <- 1.0
    cat("Warning: Could not calculate lambda, using default value of 1.0\n")
  }
  
  # Create QQ plot
  png(filename = output_file, width = width, height = height, res = 150)
  
  # Set up plot margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Create the QQ plot
  qq(p_values, 
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
  cat("Created QQ plot:", basename(output_file), "\n")
  return(output_file)
}

# Function to find p-value data files using a robust search approach
find_pvalue_files <- function(trait, population, region) {
  # List of potential directories and file patterns to search
  search_locations <- list(
    # 19_plot.ipynb generated files
    list(
      dir = file.path(base_scripts_dir, "19-OUT_plots"),
      patterns = c(
        paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_data.csv"),
        paste0("stage2_MWAS_", trait, "_", population, "_", region, ".csv")
      )
    ),
    # 37_plot_and_summarize.ipynb generated files
    list(
      dir = file.path(base_scripts_dir, "37-OUT"),
      patterns = c(
        paste0("filtered_", trait, "_results_", population, "_", region, ".csv"),
        paste0("unfiltered_", trait, "_results_", population, "_", region, ".csv"),
        paste0("filtered-bonf_", trait, "_results_a6_", population, "_", region, ".csv")
      )
    ),
    # Also check in OLD directories
    list(
      dir = file.path(base_scripts_dir, "OLD", "37-OUT"),
      patterns = c(
        paste0("filtered_", trait, "_results_", population, "_", region, ".csv"),
        paste0("unfiltered_", trait, "_results_", population, "_", region, ".csv"),
        paste0("filtered-bonf_", trait, "_results_a6_", population, "_", region, ".csv")
      )
    ),
    list(
      dir = file.path(base_scripts_dir, "OLD", "19-OUT_plots"),
      patterns = c(
        paste0("16a9par-OUT_stage2_MWAS_", trait, ".csv_", population, "_", region, "_data.csv"),
        paste0("stage2_MWAS_", trait, "_", population, "_", region, ".csv")
      )
    ),
    # Search directly in base scripts dir
    list(
      dir = base_scripts_dir,
      patterns = c(
        paste0("stage2_MWAS_", trait, "_", population, "_", region, ".csv"),
        paste0("filtered_", trait, "_results_", population, "_", region, ".csv"),
        paste0("unfiltered_", trait, "_results_", population, "_", region, ".csv")
      )
    )
  )
  
  # Check each location until we find a file
  for (location in search_locations) {
    dir_path <- location$dir
    if (!dir.exists(dir_path)) {
      cat("Directory not found:", dir_path, "\n")
      next
    }
    
    for (pattern in location$patterns) {
      file_path <- file.path(dir_path, pattern)
      
      if (file.exists(file_path)) {
        cat("Found data file:", file_path, "\n")
        return(list(
          path = file_path,
          pattern_type = if(grepl("16a9par|stage2_MWAS", pattern)) "19_format" else "default"
        ))
      }
    }
  }
  
  # Also search more broadly using find command
  cat("Searching more broadly for data files...\n")
  cmd <- sprintf('find %s -type f -name "*%s*%s*%s*.csv" | grep -v "checkpoint" | sort', 
                base_scripts_dir, trait, population, region)
  potential_files <- system(cmd, intern = TRUE)
  
  if (length(potential_files) > 0) {
    # Use the first file found (most likely to be in the right format)
    cat("Found potential data file:", potential_files[1], "\n")
    return(list(
      path = potential_files[1],
      pattern_type = if(grepl("16a9par|stage2_MWAS", potential_files[1])) "19_format" else "default"
    ))
  }
  
  cat("No data file found for:", trait, population, region, "\n")
  return(NULL)
}

# Process each combination of trait, population, and region
cat("\nStarting to process all combinations...\n")
successful_plots <- list(manhattan = 0, qq = 0)

for (trait in traits) {
  for (pop in populations) {
    for (reg in regions) {
      cat("\n--------------------------------------------\n")
      cat("Processing:", trait, pop, reg, "\n")
      
      # Find data file for this combination
      data_file_info <- find_pvalue_files(trait, pop, reg)
      
      if (is.null(data_file_info)) {
        cat("Skipping: No data file found for", trait, pop, reg, "\n")
        next
      }
      
      # Load the data
      cat("Loading data from:", data_file_info$path, "\n")
      data <- tryCatch({
        fread(data_file_info$path)
      }, error = function(e) {
        cat("Error reading file:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(data) || nrow(data) == 0) {
        cat("Error: Empty or invalid data file\n")
        next
      }
      
      cat("Loaded", nrow(data), "rows from data file\n")
      
      # Create Manhattan plot
      manhattan_result <- create_manhattan_plot(
        data, trait, pop, reg, output_dir, 
        name_pattern = data_file_info$pattern_type
      )
      
      if (!is.null(manhattan_result)) {
        successful_plots$manhattan <- successful_plots$manhattan + 1
      }
      
      # Create QQ plot
      qq_result <- create_qq_plot(
        data, trait, pop, reg, output_dir,
        name_pattern = data_file_info$pattern_type
      )
      
      if (!is.null(qq_result)) {
        successful_plots$qq <- successful_plots$qq + 1
      }
    }
  }
}

# Calculate total possible combinations
total_combinations <- length(traits) * length(populations) * length(regions)
total_possible_plots <- total_combinations * 2  # 2 plots per combination (manhattan + qq)

# Print summary of results
cat("\n\n====== SUMMARY ======\n")
cat("Total Manhattan plots created:", successful_plots$manhattan, 
    "(", round(100 * successful_plots$manhattan / total_combinations, 1), "%)\n")
cat("Total QQ plots created:", successful_plots$qq,
    "(", round(100 * successful_plots$qq / total_combinations, 1), "%)\n")
cat("Total plots created:", successful_plots$manhattan + successful_plots$qq, 
    "(", round(100 * (successful_plots$manhattan + successful_plots$qq) / total_possible_plots, 1), "%)\n")
cat("Plots saved to:", output_dir, "\n")

# Print execution time
end_time <- Sys.time()
execution_time <- end_time - start_time
cat("Execution time:", format(execution_time), "\n")
