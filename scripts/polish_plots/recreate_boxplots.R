#!/usr/bin/env Rscript

# recreate_boxplots.R
# Purpose: Recreate various boxplot visualizations for the CpGWAS project
# Author: Code Ninja
# Date: 2024-03-30

# Load required libraries with error handling
required_packages <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "viridis", "readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(sprintf("CRITICAL ERROR: Required package '%s' not found. Please install it first.", pkg))
  }
}

# Function to create timestamp for logging
get_timestamp <- function() {
  return(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

# Function to validate data before plotting
validate_data <- function(data, required_cols) {
  if (is.null(data) || nrow(data) == 0) {
    stop(sprintf("CRITICAL ERROR: Input data is NULL or empty at %s", get_timestamp()))
  }
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("CRITICAL ERROR: Missing required columns: %s at %s", 
                 paste(missing_cols, collapse = ", "), get_timestamp()))
  }
  
  # Check for invalid values
  for (col in required_cols) {
    if (any(is.na(data[[col]])) || any(is.infinite(data[[col]]))) {
      stop(sprintf("CRITICAL ERROR: Column '%s' contains NA or infinite values at %s", 
                   col, get_timestamp()))
    }
  }
  
  return(TRUE)
}

# Function to create output directory with validation
create_output_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    if (!dir.create(dir_path, recursive = TRUE)) {
      stop(sprintf("CRITICAL ERROR: Failed to create output directory: %s at %s", 
                   dir_path, get_timestamp()))
    }
  }
  return(dir_path)
}

# Function to save plot with validation
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  if (is.null(plot)) {
    stop(sprintf("CRITICAL ERROR: Plot object is NULL for %s at %s", 
                 filename, get_timestamp()))
  }
  
  # Ensure filename has proper extension
  if (!grepl("\\.png$", filename)) {
    filename <- paste0(filename, ".png")
  }
  
  # Save plot with error handling
  tryCatch({
    ggsave(filename, plot, width = width, height = height, dpi = dpi)
    if (!file.exists(filename)) {
      stop(sprintf("CRITICAL ERROR: Failed to save plot: %s at %s", 
                   filename, get_timestamp()))
    }
    message(sprintf("Successfully saved plot: %s at %s", filename, get_timestamp()))
  }, error = function(e) {
    stop(sprintf("CRITICAL ERROR: Failed to save plot %s: %s at %s", 
                 filename, e$message, get_timestamp()))
  })
}

# Function to load and validate input data
load_input_data <- function() {
  message("Loading input data...")
  
  # Define required input files
  input_files <- list(
    correlation_data = "data/correlation_data.csv",
    feature_data = "data/feature_data.csv",
    metadata = "data/metadata.csv"
  )
  
  # Check if all required files exist
  for (file_name in input_files) {
    if (!file.exists(file_name)) {
      stop(sprintf("CRITICAL ERROR: Required input file not found: %s at %s", 
                   file_name, get_timestamp()))
    }
  }
  
  # Load data with error handling
  tryCatch({
    correlation_data <- read_csv(input_files$correlation_data)
    feature_data <- read_csv(input_files$feature_data)
    metadata <- read_csv(input_files$metadata)
    
    # Validate loaded data
    validate_data(correlation_data, c("group", "value", "facet"))
    validate_data(feature_data, c("group", "value", "facet", "color_var"))
    validate_data(metadata, c("group", "correlation"))
    
    return(list(
      correlation_data = correlation_data,
      feature_data = feature_data,
      metadata = metadata
    ))
  }, error = function(e) {
    stop(sprintf("CRITICAL ERROR: Failed to load input data: %s at %s", 
                 e$message, get_timestamp()))
  })
}

# Function to create basic faceted boxplot
create_basic_faceted_boxplot <- function(data, output_dir) {
  message("Creating basic faceted boxplot...")
  
  # Validate input data
  required_cols <- c("group", "value", "facet")
  validate_data(data, required_cols)
  
  # Create plot with custom theme
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_boxplot(fill = "lightgray", alpha = 0.7) +
    facet_wrap(~facet, scales = "free_y") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(title = "Basic Faceted Boxplot",
         x = "Group",
         y = "Value")
  
  # Save plot
  save_plot(p, file.path(output_dir, "boxplot_version1_basic_faceted"))
}

# Function to create colored faceted boxplot
create_colored_faceted_boxplot <- function(data, output_dir) {
  message("Creating colored faceted boxplot...")
  
  # Validate input data
  required_cols <- c("group", "value", "facet", "color_var")
  validate_data(data, required_cols)
  
  # Create plot with custom theme and colors
  p <- ggplot(data, aes(x = group, y = value, fill = color_var)) +
    geom_boxplot(alpha = 0.8) +
    facet_wrap(~facet, scales = "free_y") +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    labs(title = "Colored Faceted Boxplot",
         x = "Group",
         y = "Value",
         fill = "Category")
  
  # Save plot
  save_plot(p, file.path(output_dir, "boxplot_version2_colored_faceted"))
}

# Function to create jittered faceted boxplot
create_jittered_faceted_boxplot <- function(data, output_dir) {
  message("Creating jittered faceted boxplot...")
  
  # Validate input data
  required_cols <- c("group", "value", "facet")
  validate_data(data, required_cols)
  
  # Create plot with custom theme and jitter
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_boxplot(fill = "lightgray", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "steelblue") +
    facet_wrap(~facet, scales = "free_y") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(title = "Jittered Faceted Boxplot",
         x = "Group",
         y = "Value")
  
  # Save plot
  save_plot(p, file.path(output_dir, "boxplot_version3_jittered_faceted"))
}

# Function to create average correlation plot for caudate
create_average_correlation_plot <- function(data, output_dir) {
  message("Creating average correlation plot for caudate...")
  
  # Validate input data
  required_cols <- c("group", "correlation")
  validate_data(data, required_cols)
  
  # Create plot with custom theme and styling
  p <- ggplot(data, aes(x = group, y = correlation)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(title = "Average Correlation in Caudate",
         x = "Group",
         y = "Correlation")
  
  # Save plot
  save_plot(p, file.path(output_dir, "average_correlation_plot_caud"))
}

# Main execution
main <- function() {
  message(sprintf("Starting boxplot recreation at %s", get_timestamp()))
  
  # Setup output directory
  output_dir <- create_output_dir("plot_outputs/boxplots")
  
  # Load and validate data
  data <- load_input_data()
  
  # Create all plots
  tryCatch({
    create_basic_faceted_boxplot(data$correlation_data, output_dir)
    create_colored_faceted_boxplot(data$feature_data, output_dir)
    create_jittered_faceted_boxplot(data$correlation_data, output_dir)
    create_average_correlation_plot(data$metadata, output_dir)
    
    message(sprintf("Successfully completed all boxplot recreation at %s", get_timestamp()))
  }, error = function(e) {
    stop(sprintf("CRITICAL ERROR: Failed to complete boxplot recreation: %s at %s", 
                 e$message, get_timestamp()))
  })
}

# Run main function
main() 