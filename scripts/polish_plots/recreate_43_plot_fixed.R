#!/usr/bin/env Rscript
# Recreation script for 43_annotate_MWAS_hits plots
# Based on 43_annotate_MWAS_hits_a11.ipynb
# CRITICAL: This script must recreate plots from data, NOT copy existing plots

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(data.table)
library(qqman)

# Define paths
base_path <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts"
output_dir <- file.path(base_path, "polish_plots/plot_outputs/43_plots_fast")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the specific files we need to recreate the original plots
# These are the exact files that were used to create the original plots in 43-OUT
required_files <- list(
  bp = list(
    path = NULL,
    pattern = "gwas.*bp.*results|GWAS.*bp.*results|bp.*gwas.*results|BP.*gwas.*results"
  ),
  mdd = list(
    path = NULL,
    pattern = "gwas.*mdd.*results|GWAS.*mdd.*results|mdd.*gwas.*results|MDD.*gwas.*results"
  )
)

# Function to find data files in multiple possible locations
find_data_file <- function(pattern) {
  possible_locations <- c(
    file.path(base_path, "OLD_scripts"),
    file.path(base_path, "OLD"),
    file.path(base_path),
    file.path(base_path, "43-OUT"),
    file.path(base_path, "37-OUT")
  )
  
  for (location in possible_locations) {
    if (dir.exists(location)) {
      files <- list.files(path = location, pattern = pattern, full.names = TRUE, recursive = TRUE)
      if (length(files) > 0) {
        return(files)
      }
    }
  }
  
  # Return empty vector if no files found
  return(character(0))
}

# Find the required input files
for (trait in names(required_files)) {
  files <- find_data_file(required_files[[trait]]$pattern)
  if (length(files) > 0) {
    required_files[[trait]]$path <- files[1]
    cat("Found file for", trait, ":", files[1], "\n")
  } else {
    stop(paste0("CRITICAL ERROR: Could not find input file for ", trait, 
                ". This script requires the original data files to recreate the plots."))
  }
}

# Function to recreate Manhattan plots from original data
recreate_manhattan_plot <- function(results_file, output_file, title, plot_type) {
  # Read and process data
  cat("Processing file:", results_file, "\n")
  
  # Try to read the file, handle different possible formats
  results <- tryCatch({
    data <- fread(results_file)
    if (ncol(data) > 0) {
      data
    } else {
      stop("Empty data file")
    }
  }, error = function(e) {
    # Try reading as CSV if fread fails
    data <- read.csv(results_file, stringsAsFactors = FALSE)
    if (ncol(data) > 0) {
      as.data.table(data)
    } else {
      stop("Empty data file")
    }
  })
  
  # Extract chromosome and position information
  # Try different column naming schemes commonly found in these files
  
  # First, identify the chromosome and position columns
  chr_col <- NULL
  pos_col <- NULL
  
  if ("probeID" %in% colnames(results)) {
    cat("  Found probeID column - extracting chromosome and position\n")
    # Extract chromosome and position from CpG IDs
    results$chromosome <- gsub("^chr(.*):.*$", "\\1", results$probeID)
    results$position <- as.numeric(gsub("^chr.*:(.*)$", "\\1", results$probeID))
    chr_col <- "chromosome"
    pos_col <- "position"
  } else if ("chr" %in% colnames(results)) {
    cat("  Found chr column\n")
    chr_col <- "chr"
    
    # Check for position column
    if ("pos" %in% colnames(results)) {
      pos_col <- "pos"
    } else if ("position" %in% colnames(results)) {
      pos_col <- "position"
    } else if ("bp" %in% colnames(results)) {
      pos_col <- "bp"
    }
  } else if ("CHR" %in% colnames(results)) {
    cat("  Found CHR column\n")
    chr_col <- "CHR"
    
    # Check for position column
    if ("POS" %in% colnames(results)) {
      pos_col <- "POS"
    } else if ("BP" %in% colnames(results)) {
      pos_col <- "BP"
    }
  }
  
  # Check if we found the necessary columns
  if (is.null(chr_col) || is.null(pos_col)) {
    stop(paste0("CRITICAL ERROR: Could not identify chromosome and position columns in ", results_file, 
                ". Available columns: ", paste(colnames(results), collapse=", ")))
  }
  
  # Identify p-value column
  p_col <- NULL
  if ("pvalue" %in% colnames(results)) {
    p_col <- "pvalue"
  } else if ("p" %in% colnames(results)) {
    p_col <- "p"
  } else if ("P" %in% colnames(results)) {
    p_col <- "P"
  } else if ("p_value" %in% colnames(results)) {
    p_col <- "p_value"
  } else if ("pval" %in% colnames(results)) {
    p_col <- "pval"
  }
  
  if (is.null(p_col)) {
    stop(paste0("CRITICAL ERROR: Could not identify p-value column in ", results_file, 
                ". Available columns: ", paste(colnames(results), collapse=", ")))
  }
  
  # Ensure chromosome is formatted correctly (numeric where possible)
  results[[chr_col]] <- gsub("X", "23", results[[chr_col]])
  results[[chr_col]] <- gsub("Y", "24", results[[chr_col]])
  results[[chr_col]] <- as.numeric(results[[chr_col]])
  
  # Remove missing or invalid values
  results <- results[!is.na(results[[chr_col]]) & 
                    !is.na(results[[pos_col]]) & 
                    !is.na(results[[p_col]]) & 
                    results[[p_col]] > 0 & 
                    results[[p_col]] <= 1]
  
  # Check if we have enough data points
  if (nrow(results) < 100) {
    stop(paste0("CRITICAL ERROR: Not enough valid data points after filtering: ", nrow(results), 
                ". The original data file may be corrupted or in an unexpected format."))
  }
  
  # Sort by chromosome and position
  results <- results[order(results[[chr_col]], results[[pos_col]]),]
  
  # Add SNP column required by manhattan function
  results$SNP <- paste0(results[[chr_col]], ":", results[[pos_col]])
  
  if (plot_type == "qqman") {
    # Create Manhattan plot using qqman
    png(filename = output_file, width = 12, height = 6, units = "in", res = 300)
    manhattan(results, chr = chr_col, bp = pos_col, p = p_col, snp = "SNP",
              main = title,
              suggestiveline = -log10(1e-5),  # Suggestive significance line
              genomewideline = -log10(5e-8),  # Genome-wide significance line
              col = c("dodgerblue4", "darkorange2"),
              cex = 0.6,
              cex.axis = 0.8)
    dev.off()
    
    cat("Created qqman plot:", output_file, "\n")
  } else if (plot_type == "ggplot2") {
    # Create Manhattan plot with ggplot2
    # Prepare data
    results$logp <- -log10(results[[p_col]])
    
    # Calculate cumulative position for x-axis
    results$pos <- 0
    offset <- 0
    chr_labels <- c()
    chr_pos <- c()
    
    for (chr in sort(unique(results[[chr_col]]))) {
      chr_data <- results[results[[chr_col]] == chr,]
      results$pos[results[[chr_col]] == chr] <- chr_data[[pos_col]] + offset
      chr_labels <- c(chr_labels, chr)
      chr_pos <- c(chr_pos, mean(chr_data[[pos_col]]) + offset)
      offset <- max(results$pos, na.rm=TRUE) + max(results[[pos_col]], na.rm=TRUE)/2
    }
    
    # Create the plot
    png(filename = output_file, width = 12, height = 6, units = "in", res = 300)
    
    # Alternate colors by chromosome
    results$color <- factor(results[[chr_col]] %% 2)
    
    p <- ggplot(results, aes(x = pos, y = logp, color = color)) +
      geom_point(alpha = 0.8, size = 0.8) +
      scale_color_manual(values = c("#4477AA", "#66CCEE")) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = title,
        x = "Chromosome",
        y = "-log10(p-value)"
      ) +
      scale_x_continuous(breaks = chr_pos, labels = chr_labels) +
      geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue") +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red")
    
    print(p)
    dev.off()
    
    cat("Created ggplot2 version:", output_file, "\n")
  }
  
  return(TRUE)
}

# Process each trait and create the original plots
for (trait in names(required_files)) {
  if (!is.null(required_files[[trait]]$path)) {
    # Create qqman plot
    output_file <- file.path(output_dir, paste0("gwas_manhattan_", trait, "_qqman.png"))
    title <- paste("Manhattan Plot -", toupper(trait))
    success <- recreate_manhattan_plot(required_files[[trait]]$path, output_file, title, "qqman")
    
    if (success) {
      cat("Successfully created qqman Manhattan plot for", trait, "\n")
    }
    
    # Create ggplot2 plot
    output_file <- file.path(output_dir, paste0("gwas_manhattan_", trait, "_ggplot2.png"))
    success <- recreate_manhattan_plot(required_files[[trait]]$path, output_file, title, "ggplot2")
    
    if (success) {
      cat("Successfully created ggplot2 Manhattan plot for", trait, "\n")
    }
  }
}

cat("Completed recreation of 43 plots\n")