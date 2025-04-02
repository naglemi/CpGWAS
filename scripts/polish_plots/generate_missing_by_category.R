#!/usr/bin/env Rscript
# Generate missing plots by category

# Load required libraries
library(data.table)

# Define directory paths
verify_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/verification_results"

# Read the missing plots file
missing_plots <- fread(file.path(verify_dir, "missing_plots.csv"))

# Generate a separate file for each category
categories <- unique(missing_plots$category)
for (cat in categories) {
  cat_data <- missing_plots[category == cat]
  if (nrow(cat_data) > 0) {
    # Create a new data frame with original_plot, recreated_plot, category, and status columns
    output_data <- data.table(
      original_plot = cat_data$filepath,
      recreated_plot = "NOT RECREATED",
      category = cat_data$category,
      status = "no_match"
    )
    
    # Write to file
    output_file <- file.path(verify_dir, paste0("missing_", cat, ".csv"))
    fwrite(output_data, output_file)
    cat("Generated:", output_file, "\n")
  }
}

cat("Done!\n")
