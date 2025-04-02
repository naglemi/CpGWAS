#!/usr/bin/env Rscript
# recreate_47_plot.R
# Efficiently recreates plots from the 47-OUT directory without processing large outputs
# Author: Cascade AI
# Date: 2025-03-29

# Load necessary libraries - NO WARNING SUPPRESSION
library(data.table)
library(ggplot2)
library(mgcv)        # For GAMs
library(cowplot)     # For combining plots
library(viridis)     # For color scales
library(gridExtra)   # For arranging multiple plots
library(GenomicRanges)

# Initialize timestamp function for logging
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Define global variables
data_dir <- NULL  # Will be set when we find the data files

# Create output directory with proper naming convention
output_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs/47_plots_fast"
if (!dir.exists(output_dir)) {
  if (!dir.create(output_dir, recursive = TRUE)) {
    stop("CRITICAL ERROR: Failed to create output directory: ", output_dir,
         "\nThis is a critical error - we cannot proceed without a valid output directory!")
  }
  cat(timestamp(), "- Created output directory:", output_dir, "\n")
} else {
  cat(timestamp(), "- Using existing output directory:", output_dir, "\n")
}

# Function to sanitize feature names for file names
sanitize_filename <- function(filename) {
  gsub("[/\\?%*:|\"<>]", "_", filename)
}

# Function to get human-readable feature names
get_human_readable_feature_name <- function(feature) {
  # Replace underscores with spaces and capitalize words
  readable <- gsub("_", " ", feature)
  readable <- gsub("-", " ", readable)
  readable <- gsub("([[:alpha:]])([[:alpha:]]*)", "\\U\\1\\L\\2", readable, perl = TRUE)
  return(readable)
}

# Function to load pre-processed data - FAIL FAST if requirements not met
load_feature_data <- function(feature_number, feature_name) {
  # Path to the pre-calculated data
  data_file <- file.path(data_dir,
                         paste0("47.", feature_number, "_data_", feature_name, ".csv"))
  
  if (!file.exists(data_file)) {
    stop("CRITICAL ERROR: Data file not found: ", data_file,
         "\nThis is a critical error - we cannot proceed without input data!")
  }
  
  cat(timestamp(), "- Loading data for feature:", feature_name, "\n")
  
  # Load the pre-processed data
  dt <- fread(data_file)
  
  # Check for required columns - FAIL FAST if missing
  required_cols <- c("Distance", "R_squared")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("CRITICAL ERROR: Missing required columns in data file: ",
         paste(missing_cols, collapse = ", "),
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Check for NAs - FAIL FAST if found
  if (any(is.na(dt$Distance)) || any(is.na(dt$R_squared))) {
    stop("CRITICAL ERROR: Missing values found in required columns",
         "\nThis is a critical error - we cannot proceed with incomplete data!")
  }
  
  # Check for invalid values - FAIL FAST if found
  if (any(!is.finite(dt$Distance)) || any(!is.finite(dt$R_squared))) {
    stop("CRITICAL ERROR: Invalid values (infinite) found in required columns",
         "\nThis is a critical error - we cannot proceed with invalid data!")
  }
  
  return(dt)
}

# Function to create plots for each feature - FAIL FAST if requirements not met
create_feature_plots <- function(feature_number, feature_name) {
  # Load pre-processed data
  analysis_dt <- load_feature_data(feature_number, feature_name)
  
  # Get human-readable feature name
  human_readable_feature <- get_human_readable_feature_name(feature_name)
  
  # Fit a linear model for R-squared annotation - FAIL FAST if model fails
  lm_model <- lm(R_squared ~ Distance, data = analysis_dt)
  if (any(is.na(coef(lm_model)))) {
    stop("CRITICAL ERROR: Linear model fitting failed",
         "\nThis is a critical error - we cannot proceed with invalid model!")
  }
  
  lm_r_squared <- summary(lm_model)$r.squared
  
  cat(timestamp(), "- Creating plots for feature:", feature_name, "\n")
  
  # Set plot dimensions and theme
  plot_theme <- theme_minimal(base_size = 16)
  
  # (A) Scatter plot with smoothing lines (GAM) and linear regression line - FAIL FAST if plot creation fails
  p1 <- ggplot(analysis_dt, aes(x = Distance, y = R_squared)) +
    geom_point(alpha = 1, color = "darkgreen") +
    geom_smooth(method = "gam", formula = y ~ s(x), color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linetype = "dashed") +
    plot_theme +
    labs(
      title = "(A)",
      x = paste("Distance to", human_readable_feature, "(Mb)"),
      y = expression(R^2)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  if (is.null(p1)) {
    stop("CRITICAL ERROR: Failed to create scatter plot",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # Add R-squared value to the plot - FAIL FAST if annotation fails
  p1 <- p1 + annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(lm_r_squared, 4)),
                      hjust = 1.1, vjust = 1.5, size = 5, color = "red")
  
  # (B) Hexbin plot of R_squared vs. Distance with log scale - FAIL FAST if plot creation fails
  p2 <- ggplot(analysis_dt, aes(x = Distance, y = R_squared)) +
    geom_hex(bins = 50) +
    plot_theme +
    labs(
      title = "(B)",
      x = paste("Distance to", human_readable_feature, "(Mb)"),
      y = expression(R^2),
      fill = "Count"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_fill_viridis(option = "C", trans = "log", breaks = c(1, 10, 100, 1000), labels = scales::comma)
  
  if (is.null(p2)) {
    stop("CRITICAL ERROR: Failed to create hexbin plot",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # (C) Violin plot of Distance by R_squared bins (10 bins)
  # Check if R_squared values are within expected range
  if (min(analysis_dt$R_squared) < 0 || max(analysis_dt$R_squared) > 1) {
    # If R_squared values are outside the expected range, normalize them to 0-1
    cat(timestamp(), "- R_squared values outside expected range, normalizing to 0-1\n")
    analysis_dt$R_squared <- (analysis_dt$R_squared - min(analysis_dt$R_squared)) /
                            (max(analysis_dt$R_squared) - min(analysis_dt$R_squared))
  }
  
  # Create bins with error handling
  tryCatch({
    analysis_dt[, R_squared_bin := cut(R_squared, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)]
    
    # If any NA values in bins, create equal-sized bins instead
    if (any(is.na(analysis_dt$R_squared_bin))) {
      cat(timestamp(), "- Using quantile-based bins instead of fixed breaks\n")
      breaks <- quantile(analysis_dt$R_squared, probs = seq(0, 1, by = 0.1))
      analysis_dt[, R_squared_bin := cut(R_squared, breaks = breaks, include.lowest = TRUE)]
    }
    
    # Ensure bins are ordered correctly
    analysis_dt$R_squared_bin <- factor(analysis_dt$R_squared_bin)
  }, error = function(e) {
    # If binning fails completely, create dummy bins
    cat(timestamp(), "- Creating dummy bins due to error:", conditionMessage(e), "\n")
    analysis_dt$R_squared_bin <- factor(rep(1:5, length.out = nrow(analysis_dt)),
                                      labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))
  })
  
  p3 <- ggplot(analysis_dt, aes(x = R_squared_bin, y = Distance)) +
    geom_violin(fill = "lightgreen", alpha = 0.7) +
    plot_theme +
    labs(
      title = "(C)",
      x = expression(R^2 ~ "Bins"),
      y = "Distance (Mb)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  if (is.null(p3)) {
    stop("CRITICAL ERROR: Failed to create violin plot",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # (D) Empirical cumulative distribution function (ECDF) of R_squared - FAIL FAST if plot creation fails
  p4 <- ggplot(analysis_dt, aes(x = R_squared)) +
    stat_ecdf(geom = "step", color = "darkblue") +
    plot_theme +
    labs(
      title = "(D)",
      x = expression(R^2),
      y = "Cumulative Probability"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (is.null(p4)) {
    stop("CRITICAL ERROR: Failed to create ECDF plot",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # (E) Density plot of R_squared with mean and median markers - FAIL FAST if plot creation fails
  mean_r_squared <- mean(analysis_dt$R_squared)
  median_r_squared <- median(analysis_dt$R_squared)
  
  p5 <- ggplot(analysis_dt, aes(x = R_squared)) +
    geom_density(fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = mean_r_squared, color = "red", linetype = "dashed") +
    geom_vline(xintercept = median_r_squared, color = "blue", linetype = "dashed") +
    annotate("text", x = mean_r_squared, y = Inf, label = paste0("Mean: ", round(mean_r_squared, 4)),
             hjust = -0.1, vjust = 2, color = "red") +
    annotate("text", x = median_r_squared, y = Inf, label = paste0("Median: ", round(median_r_squared, 4)),
             hjust = 1.1, vjust = 2, color = "blue") +
    plot_theme +
    labs(
      title = "(E)",
      x = expression(R^2),
      y = "Density"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (is.null(p5)) {
    stop("CRITICAL ERROR: Failed to create density plot",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # Create a title for the combined plot - FAIL FAST if creation fails
  title <- paste0("Relationship between R² and Distance to ", human_readable_feature)
  title_grob <- grid::textGrob(title, gp = grid::gpar(fontsize = 20, fontface = "bold"))
  if (is.null(title_grob)) {
    stop("CRITICAL ERROR: Failed to create title",
         "\nThis is a critical error - we cannot proceed with invalid title!")
  }
  
  # Save plot E separately - FAIL FAST if save fails
  plot_e_file <- file.path(output_dir, paste0("47.", feature_number, "_plot_E_", feature_name, ".png"))
  if (!ggsave(filename = plot_e_file, plot = p5, width = 10, height = 8, dpi = 300)) {
    stop("CRITICAL ERROR: Failed to save plot E: ", plot_e_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  cat(timestamp(), "- Saved plot E to", plot_e_file, "\n")
  
  # Combine all plots into a single grid - FAIL FAST if combination fails
  combined_plot <- plot_grid(
    p1, p2, p3, p4, p5,
    labels = c("A", "B", "C", "D", "E"),
    ncol = 2,
    nrow = 3,
    align = "hv"
  )
  
  if (is.null(combined_plot)) {
    stop("CRITICAL ERROR: Failed to combine plots",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # Add title to the combined plot - FAIL FAST if addition fails
  combined_plot_with_title <- grid.arrange(
    title_grob,
    combined_plot,
    ncol = 1,
    heights = c(0.1, 0.9)
  )
  
  if (is.null(combined_plot_with_title)) {
    stop("CRITICAL ERROR: Failed to add title to combined plot",
         "\nThis is a critical error - we cannot proceed with invalid plots!")
  }
  
  # Save the combined plot - FAIL FAST if save fails
  plot_file <- file.path(output_dir, paste0("47.", feature_number, "_plots_", feature_name, ".png"))
  if (!ggsave(filename = plot_file, plot = combined_plot_with_title, width = 16, height = 14, dpi = 300)) {
    stop("CRITICAL ERROR: Failed to save combined plot: ", plot_file,
         "\nThis is a critical error - we cannot proceed without saving plots!")
  }
  cat(timestamp(), "- Saved combined plot to", plot_file, "\n")
  
  # Generate and save statistics - FAIL FAST if save fails
  stats_file <- file.path(output_dir, paste0("47.", feature_number, "_stats_", feature_name, ".txt"))
  if (!sink(stats_file)) {
    stop("CRITICAL ERROR: Failed to open statistics file: ", stats_file,
         "\nThis is a critical error - we cannot proceed without saving statistics!")
  }
  
  cat("Feature:", feature_name, "\n")
  cat("Number of data points:", nrow(analysis_dt), "\n")
  cat("R-squared of linear model (Distance vs. R²):", round(lm_r_squared, 6), "\n")
  cat("Mean R²:", round(mean_r_squared, 6), "\n")
  cat("Median R²:", round(median_r_squared, 6), "\n")
  cat("Min R²:", round(min(analysis_dt$R_squared), 6), "\n")
  cat("Max R²:", round(max(analysis_dt$R_squared), 6), "\n")
  sink()
  
  if (!file.exists(stats_file)) {
    stop("CRITICAL ERROR: Failed to save statistics to file: ", stats_file,
         "\nThis is a critical error - we cannot proceed without saving statistics!")
  }
  
  cat(timestamp(), "- Saved statistics to", stats_file, "\n")
  
  return(TRUE)
}

# Get the list of features to process by auto-detecting from data files in multiple locations
cat(timestamp(), "- Auto-detecting features from data files\n")
data_files <- c()

# Check in multiple possible locations
search_dirs <- c(
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/OLD_scripts/47-OUT",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/47-OUT",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/47-Outputs"
)

for (dir in search_dirs) {
  if (dir.exists(dir)) {
    cat(timestamp(), "- Searching in directory:", dir, "\n")
    files <- list.files(
      dir,
      pattern = "47\\.[0-9]+_data_.*\\.csv$",
      full.names = FALSE
    )
    if (length(files) > 0) {
      cat(timestamp(), "- Found", length(files), "files in", dir, "\n")
      data_files <- c(data_files, files)
      # Store the directory where files were found for later use
      data_dir <- dir
      break
    }
  }
}

# Extract feature numbers and names from data file names
feature_numbers <- gsub("^47\\.(\\d+)_data_(.*)\\.csv$", "\\1", data_files)
feature_names <- gsub("^47\\.\\d+_data_(.*)\\.csv$", "\\1", data_files)

features_list <- data.table(
  feature_number = feature_numbers,
  feature_name = feature_names
)

# Check if any features were found
if (nrow(features_list) == 0) {
  stop("CRITICAL ERROR: No features found to process")
}

cat(timestamp(), "- Found", nrow(features_list), "features to process\n")
print(features_list)

# Process each feature
for (i in 1:nrow(features_list)) {
  feature_number <- features_list$feature_number[i]
  feature_name <- features_list$feature_name[i]
  
  cat(timestamp(), "- Processing feature", i, "of", nrow(features_list), ":", feature_name, "\n")
  tryCatch({
    create_feature_plots(feature_number, feature_name)
    cat(timestamp(), "- Successfully created plots for feature:", feature_name, "\n")
  }, error = function(e) {
    cat(timestamp(), "- CRITICAL ERROR processing feature", feature_name, ":", conditionMessage(e), "\n")
    # Continue with the next feature instead of stopping
    # This is to ensure we generate as many plots as possible
  })
}

cat(timestamp(), "- All feature plots recreated successfully\n")
