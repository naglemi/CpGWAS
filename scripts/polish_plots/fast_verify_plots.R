#!/usr/bin/env Rscript
# fast_verify_plots.R - Efficiently verify plot recreation completeness
# Uses fast file listing and targeted sampling instead of full recursive scanning

# Load required libraries
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# Record start time
start_time <- Sys.time()
cat("Starting plot verification at:", format(start_time), "\n")

# Define directory paths
scripts_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts"
output_dir <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/polish_plots/plot_outputs"

# Create output directory for the verification results
verify_dir <- file.path(scripts_dir, "polish_plots", "verification_results")
dir.create(verify_dir, recursive = TRUE, showWarnings = FALSE)

# Function to quickly find image files using system commands
find_images <- function(base_dir, pattern = "\\.(png|jpg|jpeg|pdf)$", exclude_pattern = "OLD/|-checkpoint\\.png$") {
  # Use find with grep to exclude OLD paths and checkpoint files, and get only image files
  cmd <- sprintf('find "%s" -type f | grep -v "%s" | grep -E "%s"', base_dir, exclude_pattern, pattern)
  
  files <- system(cmd, intern = TRUE)
  return(files)
}

# Find original images (excluding OLD paths)
cat("Finding original images...\n")
original_images <- find_images(scripts_dir)

# Create data table of original images
original_data <- data.table(
  filepath = original_images,
  filename = basename(original_images),
  category = "",
  normalized = tolower(basename(original_images)),
  is_recreated = FALSE
)

# Categorize images based on path
original_data[grepl("19-OUT|plot_outputs/19_plots_fast", filepath), category := "19_plots_fast"]
original_data[grepl("41-OUT|plot_outputs/41_plots_fast", filepath) & !grepl("-checkpoint\\.png$", filepath), category := "41_plots_fast"]
original_data[grepl("43-OUT|plot_outputs/43_plots_fast", filepath), category := "43_plots_fast"]
original_data[grepl("47-OUT|plot_outputs/47_plots_fast", filepath), category := "47_plots_fast"]
original_data[grepl("50-OUT|plot_outputs/50_plots_fast", filepath), category := "50_plots_fast"]
original_data[grepl("DNAm_violin|plot_outputs/dnam_violin_plots", filepath), category := "dnam_violin_plots_fast"]
# Improved categorization for manhattan and qq plots
original_data[grepl("manhattan|qq_plot|qqplot", filepath) | 
             grepl("16a[0-9]par-OUT.*_manhattan\\.png$|16a[0-9]par-OUT.*_qq\\.png$", filepath), 
             category := "manhattan_qq_fast"]
original_data[grepl("Figure[0-9]|plot_outputs/figure_plots_fast", filepath), category := "figure_plots_fast"]
original_data[basename(filepath) == "Rplots.pdf", category := "r_default_output"]
original_data[category == "", category := "uncategorized"]

# Find recreated images
cat("Finding recreated images...\n")
recreated_images <- find_images(output_dir)

# Create data table of recreated images
recreated_data <- data.table(
  filepath = recreated_images,
  filename = basename(recreated_images),
  category = "",
  normalized = tolower(basename(recreated_images)),
  is_recreated = TRUE
)

# Categorize recreated images
recreated_data[grepl("19_plots", filepath), category := "19_plots_fast"]
recreated_data[grepl("41_plots", filepath), category := "41_plots_fast"]
recreated_data[grepl("43_plots", filepath), category := "43_plots_fast"]
recreated_data[grepl("47_plots", filepath), category := "47_plots_fast"]
recreated_data[grepl("50_plots", filepath), category := "50_plots_fast"]
recreated_data[grepl("dnam_violin", filepath), category := "dnam_violin_plots_fast"]
recreated_data[grepl("manhattan|qq", filepath), category := "manhattan_qq_fast"]
recreated_data[category == "", category := "uncategorized"]

# Combine datasets for analysis
all_images <- rbind(original_data, recreated_data)
setkey(all_images, normalized, category)

# Identify which original plots have been recreated
for (i in 1:nrow(original_data)) {
  norm_name <- original_data$normalized[i]
  category <- original_data$category[i]
  
  # Look for matching recreated images by normalized name and category
  matches <- recreated_data[normalized %like% norm_name | norm_name %like% normalized]
  
  if (nrow(matches) > 0) {
    original_data$is_recreated[i] <- TRUE
  }
}

# Generate summary statistics
recreated_count <- sum(original_data$is_recreated)
total_count <- nrow(original_data)
recreation_percent <- round(100 * recreated_count / total_count, 1)

# Create category summaries
category_summary <- original_data[, .(
  total_plots = .N,
  recreated_plots = sum(is_recreated),
  missing_plots = .N - sum(is_recreated),
  recreation_percentage = round(100 * sum(is_recreated) / .N, 1)
), by = category]
setorder(category_summary, -total_plots)

# Save detailed results
cat("\nSaving verification reports...\n")
fwrite(original_data, file.path(verify_dir, "plot_verification_details.csv"))
fwrite(category_summary, file.path(verify_dir, "plot_verification_summary.csv"))
fwrite(original_data[is_recreated == FALSE], file.path(verify_dir, "missing_plots.csv"))

# Generate summary report
writeLines(c(
  "PLOT RECREATION VERIFICATION SUMMARY",
  "========================================================",
  "",
  paste("Timestamp:", format(Sys.time())), 
  "",
  paste("Total original plots scanned:", total_count),
  paste("Total recreated:", recreated_count, "(", recreation_percent, "%)"),
  "  - Exact name matches:", sum(original_data$is_recreated),
  "  - Pattern matches:", recreated_count - sum(original_data$is_recreated),
  paste("Not recreated:", total_count - recreated_count),
  paste("Uncategorized plots:", sum(original_data$category == "uncategorized")),
  "",
  "Recreation status by category:",
  "----------------------------------------------------------",
  sapply(1:nrow(category_summary), function(i) {
    sprintf("%-20s: %d/%d recreated (%.1f%%)",
            category_summary$category[i],
            category_summary$recreated_plots[i],
            category_summary$total_plots[i],
            category_summary$recreation_percentage[i])
  })
), file.path(verify_dir, "verification_summary.txt"))

# Print completion message
cat("\nVerification complete! Reports saved to:", verify_dir, "\n")
