#!/usr/bin/env Rscript
# analyze_uncategorized_plots.R - Script to analyze uncategorized plots and generate reports
# This script reads the plot verification results and generates detailed reports for each category of uncategorized plots

# Load necessary libraries
library(data.table)
library(stringr)

# Print header
cat("========================================================\n")
cat("          UNCATEGORIZED PLOTS ANALYSIS TOOL             \n")
cat("========================================================\n\n")

# Define directories
verify_dir <- "./verification_results"
report_dir <- "./verification_results/uncategorized_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# Read the verification results
cat("Reading verification results...\n")
if (file.exists(file.path(verify_dir, "plot_verification_details.csv"))) {
  all_plots <- fread(file.path(verify_dir, "plot_verification_details.csv"))
} else {
  stop("Verification results not found. Run fast_verify_plots.R first.")
}

# Extract uncategorized plots
uncategorized <- all_plots[category == "uncategorized"]
cat("Found", nrow(uncategorized), "uncategorized plots\n")

# Function to extract directory pattern from a set of file paths
extract_pattern <- function(paths) {
  # Extract directories
  dirs <- dirname(paths)
  
  # Find common patterns
  common_patterns <- table(dirs)
  
  # Sort by frequency
  sorted_patterns <- sort(common_patterns, decreasing = TRUE)
  
  return(sorted_patterns)
}

# Function to generate a report for a subset of plots
generate_report <- function(plots, pattern, output_file) {
  cat("Generating report for pattern:", pattern, "\n")
  
  # Create report content
  report <- c(
    paste("# Uncategorized Plots Report -", pattern),
    "",
    paste("Generated on:", format(Sys.time())),
    "",
    paste("Total plots matching this pattern:", nrow(plots)),
    paste("Recreation status: ", sum(plots$is_recreated), "recreated,", 
          nrow(plots) - sum(plots$is_recreated), "missing"),
    "",
    "## Plot List",
    ""
  )
  
  # Add plot details
  for (i in 1:nrow(plots)) {
    status <- ifelse(plots$is_recreated[i], "✓ RECREATED", "✗ MISSING")
    report <- c(report, 
                paste0("### ", i, ". ", basename(plots$filepath[i])),
                paste0("- **Path**: ", plots$filepath[i]),
                paste0("- **Status**: ", status),
                paste0("- **Missing Reason**: ", plots$missing_reason[i]),
                "")
  }
  
  # Add recommendations
  report <- c(report,
              "## Recommendations",
              "",
              "To properly categorize these plots, update the fast_verify_plots.R script with the following pattern:",
              "```r",
              paste0('original_data[grepl("', pattern, '", filepath), category := "', 
                     gsub(".*/", "", pattern), '_plots"]'),
              "```",
              "",
              "To recreate these plots, create a new script following the pattern of other recreation scripts:",
              "```r",
              paste0("#!/usr/bin/env Rscript"),
              paste0("# recreate_", gsub(".*/", "", pattern), "_plots.R - Script to recreate plots from ", pattern),
              "# Load necessary libraries",
              "library(data.table)",
              "library(ggplot2)",
              "",
              "# Create output directory",
              paste0('output_dir <- "./plot_outputs/', gsub(".*/", "", pattern), '_plots_fast"'),
              "dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)",
              "```")
  
  # Write report to file
  writeLines(report, output_file)
  cat("Report saved to:", output_file, "\n")
}

# Analyze patterns in uncategorized plots
cat("\nAnalyzing patterns in uncategorized plots...\n")
patterns <- extract_pattern(uncategorized$filepath)

# Print top patterns
cat("\nTop directory patterns in uncategorized plots:\n")
for (i in 1:min(10, length(patterns))) {
  cat(paste0(i, ". ", names(patterns)[i], " (", patterns[i], " plots)\n"))
}

# Generate reports for top patterns
cat("\nGenerating reports for top patterns...\n")
for (i in 1:min(5, length(patterns))) {
  pattern <- names(patterns)[i]
  plots_subset <- uncategorized[dirname(filepath) == pattern]
  
  # Skip if too few plots
  if (nrow(plots_subset) < 3) next
  
  # Generate report
  pattern_name <- gsub(".*/", "", pattern)
  output_file <- file.path(report_dir, paste0("uncategorized_", pattern_name, ".md"))
  generate_report(plots_subset, pattern, output_file)
}

# Generate report for 37-OUT qqplot files
qqplot_files <- uncategorized[grepl("qqplot", basename(filepath))]
if (nrow(qqplot_files) > 0) {
  output_file <- file.path(report_dir, "uncategorized_qqplots.md")
  generate_report(qqplot_files, "qqplot files", output_file)
}

# Generate report for figure files
figure_files <- uncategorized[grepl("Figure", basename(filepath))]
if (nrow(figure_files) > 0) {
  output_file <- file.path(report_dir, "uncategorized_figures.md")
  generate_report(figure_files, "Figure files", output_file)
}

# Generate summary report
cat("\nGenerating summary report...\n")
summary_report <- c(
  "# Uncategorized Plots Summary Report",
  "",
  paste("Generated on:", format(Sys.time())),
  "",
  paste("Total uncategorized plots:", nrow(uncategorized)),
  paste("Recreation status: ", sum(uncategorized$is_recreated), "recreated,", 
        nrow(uncategorized) - sum(uncategorized$is_recreated), "missing"),
  "",
  "## Directory Patterns",
  ""
)

# Add pattern details
for (i in 1:min(10, length(patterns))) {
  pattern <- names(patterns)[i]
  plots_count <- patterns[i]
  recreated_count <- sum(uncategorized[dirname(filepath) == pattern]$is_recreated)
  
  summary_report <- c(summary_report,
                      paste0("### ", i, ". ", pattern),
                      paste0("- **Count**: ", plots_count, " plots"),
                      paste0("- **Recreated**: ", recreated_count, " (", 
                             round(100 * recreated_count / plots_count, 1), "%)"),
                      paste0("- **Report**: [Link](./uncategorized_reports/uncategorized_", 
                             gsub(".*/", "", pattern), ".md)"),
                      "")
}

# Write summary report
writeLines(summary_report, file.path(report_dir, "uncategorized_summary.md"))
cat("Summary report saved to:", file.path(report_dir, "uncategorized_summary.md"), "\n")

cat("\nAnalysis complete! Reports saved to", report_dir, "\n")