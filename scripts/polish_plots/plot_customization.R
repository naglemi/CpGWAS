#!/usr/bin/env Rscript
# Advanced plot customization script for MWAS results
# This script provides enhanced versions of manhattan and QQ plots with more customization options

# Libraries needed for visualization
library(qqman)
library(ggplot2)
library(data.table)
library(scales)  # For better color palettes
library(gridExtra)  # For arranging multiple plots

# Enhanced Manhattan plot function
enhanced_manhattan_plot <- function(data, trait, population, region, 
                                   output_dir = "plot_outputs/manhattan_qq_advanced",
                                   width = 1800, height = 1200, 
                                   highlight_threshold = 1e-7,
                                   color_scheme = "viridis") {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create filename for the plot
  output_file <- file.path(output_dir, paste0("manhattan_advanced_", trait, "_", population, "_", region, ".png"))
  
  # Ensure SNP column exists for the manhattan function
  if (!"SNP" %in% colnames(data)) {
    data[, SNP := paste0(chr, ":", cg)]
  }
  
  # Define color schemes
  if (color_scheme == "viridis") {
    col_palette <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")
  } else if (color_scheme == "brewer") {
    col_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
  } else {
    col_palette <- c("#4286f4", "#f44182")  # Default
  }
  
  # Highlight significant points
  highlight <- data[p < highlight_threshold, SNP]
  
  # Open PNG device with high resolution
  png(filename = output_file, width = width, height = height, res = 220)
  
  # Plot with enhanced parameters
  manhattan(data, chr = "chr", bp = "cg", p = "p", snp = "SNP",
            highlight = highlight,
            main = paste0("Manhattan Plot: ", toupper(trait), " in ", toupper(region), " - ", population),
            cex = 0.8, cex.axis = 1.0, 
            col = col_palette,
            suggestiveline = -log10(1e-5), 
            genomewideline = -log10(5e-8),
            chrlabs = as.character(1:22),
            xlab = "Chromosome", 
            ylab = expression(-log[10](p-value)))
  
  # Add subtitle with details
  title(sub = paste0("Points highlighted at p < ", format(highlight_threshold, scientific = TRUE, digits = 2)), 
        cex.sub = 0.9)
  
  # Close device
  dev.off()
  message("Created enhanced Manhattan plot: ", output_file)
  
  return(output_file)
}

# Enhanced QQ plot function with confidence intervals
enhanced_qq_plot <- function(p_values, trait, population, region, 
                            output_dir = "plot_outputs/manhattan_qq_advanced",
                            width = 1200, height = 1200,
                            conf_interval = TRUE,
                            lambda_include = TRUE) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create filename for the plot
  output_file <- file.path(output_dir, paste0("qqplot_advanced_", trait, "_", population, "_", region, ".png"))
  
  # Calculate genomic inflation factor (lambda)
  lambda <- median(qchisq(1-p_values, 1)) / qchisq(0.5, 1)
  
  # Open PNG device with high resolution
  png(filename = output_file, width = width, height = height, res = 220)
  
  # Set up plot margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Create QQ plot with confidence intervals if requested
  qq(p_values, 
     main = if(lambda_include) {
       paste0("QQ Plot: ", toupper(trait), " in ", toupper(region), " - ", population, 
              "\nGenomic Inflation Factor λ = ", format(lambda, digits = 3))
     } else {
       paste0("QQ Plot: ", toupper(trait), " in ", toupper(region), " - ", population)
     },
     cex = 0.9, 
     cex.axis = 1.1, 
     pch = 19, 
     col = "navy",
     ci = conf_interval,  # Add confidence intervals
     ci.col = "lightblue",
     ci.alpha = 0.5,
     xlab = expression(Expected~~-log[10](p)),
     ylab = expression(Observed~~-log[10](p)))
  
  # Add a diagonal reference line
  abline(0, 1, col = "red", lwd = 2)
  
  # Close device
  dev.off()
  message("Created enhanced QQ plot: ", output_file)
  
  return(output_file)
}

# Create combined visualizations (Manhattan + QQ) function
create_combined_plot <- function(manhattan_file, qq_file, trait, population, region,
                                output_dir = "plot_outputs/manhattan_qq_advanced") {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create filename for the combined plot
  output_file <- file.path(output_dir, paste0("combined_", trait, "_", population, "_", region, ".png"))
  
  # Read the individual plot images
  manhattan_img <- png::readPNG(manhattan_file)
  qq_img <- png::readPNG(qq_file)
  
  # Create a grid of the two images
  png(filename = output_file, width = 2400, height = 1200, res = 220)
  
  grid.arrange(
    rasterGrob(manhattan_img, interpolate = TRUE),
    rasterGrob(qq_img, interpolate = TRUE),
    ncol = 2,
    top = grid::textGrob(
      paste0("MWAS Results: ", toupper(trait), " in ", toupper(region), " - ", population),
      gp = grid::gpar(fontsize = 20, fontface = "bold")
    )
  )
  
  dev.off()
  message("Created combined plot: ", output_file)
  
  return(output_file)
}

# Sample code to demonstrate usage:
# trait <- "scz"
# population <- "AA"
# region <- "caud"
# 
# # Load data
# filtered_file <- file.path("plot_outputs/manhattan_qq", 
#                           paste0("filtered_", trait, "_results_", population, "_", region, ".csv"))
# unfiltered_file <- file.path("plot_outputs/manhattan_qq", 
#                             paste0("unfiltered_", trait, "_results_", population, "_", region, ".csv"))
# 
# if (file.exists(filtered_file) && file.exists(unfiltered_file)) {
#   merged_subset <- fread(filtered_file)
#   merged <- fread(unfiltered_file)
#   
#   # Sample for QQ plot
#   set.seed(42)
#   merged_sample <- merged[sample(1:nrow(merged), min(10000, nrow(merged))), ]
#   
#   # Create enhanced plots
#   manhattan_file <- enhanced_manhattan_plot(merged_subset, trait, population, region)
#   qq_file <- enhanced_qq_plot(merged_sample$p, trait, population, region)
#   
#   # Create combined visualization
#   create_combined_plot(manhattan_file, qq_file, trait, population, region)
# }