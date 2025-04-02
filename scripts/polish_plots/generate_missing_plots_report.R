#!/usr/bin/env Rscript
# generate_missing_plots_report.R - Script to generate a detailed report of missing plots
# This script analyzes the plot verification results and generates a detailed report of missing plots

# Load necessary libraries
library(data.table)
library(dplyr)
library(stringr)

# Print header
cat("========================================================\n")
cat("          MISSING PLOTS ANALYSIS TOOL                   \n")
cat("========================================================\n\n")

# Define directories
verify_dir <- "./verification_results"
report_file <- "./missing_plots_report.md"

# Read the verification results
cat("Reading verification results...\n")
if (file.exists(file.path(verify_dir, "plot_verification_details.csv"))) {
  all_plots <- fread(file.path(verify_dir, "plot_verification_details.csv"))
} else {
  stop("Verification results not found. Run fast_verify_plots.R first.")
}

# Extract missing plots
missing_plots <- all_plots[is_recreated == FALSE]
cat("Found", nrow(missing_plots), "missing plots out of", nrow(all_plots), "total plots\n")
cat("Overall recreation rate:", round(100 * (nrow(all_plots) - nrow(missing_plots)) / nrow(all_plots), 1), "%\n\n")

# Function to extract information from plot filename
extract_plot_info <- function(filepath, category) {
  filename <- basename(filepath)
  
  # Extract information based on category
  if (category == "19_plots_fast") {
    # Example: 16a9par-OUT_stage2_MWAS_scz.csv_AA_caud_manhattan.png
    parts <- str_match(filename, "(.*?)_(.*?)_(.*?)_(.*?)\\.(png|pdf)")
    if (!is.na(parts[1])) {
      trait <- parts[3]
      population <- parts[4]
      region <- parts[5]
      plot_type <- ifelse(grepl("manhattan", filename), "Manhattan plot", 
                         ifelse(grepl("qq", filename), "QQ plot", 
                               ifelse(grepl("hist", filename), "Histogram", "Other")))
      return(paste0("Trait: ", trait, ", Population: ", population, ", Region: ", region, ", Type: ", plot_type))
    }
  } else if (category == "41_plots_fast") {
    # Example: scatter_zstats_AA_vs_EA_schizophrenia_caud_fulldata-v2-dupfilter.png
    if (grepl("scatter", filename)) {
      parts <- str_match(filename, "scatter_zstats_(.*?)_vs_(.*?)_(.*?)_(.*?)_.*\\.(png|pdf)")
      if (!is.na(parts[1])) {
        pop1 <- parts[2]
        pop2 <- parts[3]
        trait <- parts[4]
        region <- parts[5]
        return(paste0("Z-statistics scatter plot comparing ", pop1, " vs ", pop2, " for ", trait, " in ", region))
      }
    } else if (grepl("boxplot", filename)) {
      parts <- str_match(filename, "boxplot_cor_(.*?)_(.*?)_.*\\.(png|pdf)")
      if (!is.na(parts[1])) {
        trait <- parts[2]
        region <- parts[3]
        return(paste0("Correlation boxplot for ", trait, " in ", region))
      }
    }
  } else if (category == "43_plots_fast") {
    # Example: manhattan_unfiltered_bp_AA.png
    if (grepl("manhattan", filename)) {
      parts <- str_match(filename, "manhattan_(.*?)_(.*?)_(.*?)\\.(png|pdf)")
      if (!is.na(parts[1])) {
        filter_status <- parts[2]
        trait <- parts[3]
        population <- parts[4]
        return(paste0("Manhattan plot for ", trait, " in ", population, " (", filter_status, ")"))
      }
    } else if (grepl("correlation", filename)) {
      parts <- str_match(filename, "correlation_(.*?)_(.*?)_(.*?)\\.(png|pdf)")
      if (!is.na(parts[1])) {
        filter_status1 <- parts[2]
        filter_status2 <- parts[3]
        trait <- parts[4]
        return(paste0("Correlation plot for ", trait, " (", filter_status1, " vs ", filter_status2, ")"))
      }
    }
  } else if (category == "47_plots_fast") {
    # Example: 47.1_plot_E_astrocyte-promoters-Nott2019.png
    parts <- str_match(filename, "(\\d+\\.\\d+)_plot_E_(.*?)\\.(png|pdf)")
    if (!is.na(parts[1])) {
      plot_num <- parts[2]
      feature <- parts[3]
      return(paste0("Genomic feature plot #", plot_num, " for ", feature))
    }
  } else if (category == "50_plots_fast") {
    # Example: 50.1_plot_E_astrocyte-promoters-Nott2019.png
    parts <- str_match(filename, "(\\d+\\.\\d+)_plot_E_(.*?)\\.(png|pdf)")
    if (!is.na(parts[1])) {
      plot_num <- parts[2]
      feature <- parts[3]
      return(paste0("Overlap plot #", plot_num, " for ", feature))
    }
  } else if (category == "dnam_violin_plots_fast") {
    # Example: violin_plot_dnam_level_by_region.png
    if (grepl("violin", filename)) {
      if (grepl("level", filename)) {
        return("Violin plot of DNA methylation levels by region")
      } else if (grepl("variance", filename)) {
        return("Violin plot of DNA methylation variance by region")
      }
    }
  } else if (category == "manhattan_qq_fast") {
    if (grepl("manhattan", filename)) {
      parts <- str_match(filename, ".*_manhattan_(.*?)_(.*?)\\.(png|pdf)")
      if (!is.na(parts[1])) {
        trait <- parts[2]
        region <- parts[3]
        return(paste0("Manhattan plot for ", trait, " in ", region))
      }
    } else if (grepl("qq", filename)) {
      parts <- str_match(filename, ".*_qq_(.*?)_(.*?)\\.(png|pdf)")
      if (!is.na(parts[1])) {
        trait <- parts[2]
        region <- parts[3]
        return(paste0("QQ plot for ", trait, " in ", region))
      }
    }
  } else if (category == "figure_plots_fast") {
    if (grepl("Figure1", filename)) {
      return("DNA methylation level distribution plot")
    } else if (grepl("Figure2", filename)) {
      return("DNA methylation variance distribution plot")
    } else if (grepl("Figure3", filename)) {
      return("Combined DNA methylation level and variance plot")
    } else if (grepl("Figure4", filename)) {
      return("Biological insight plot (CpG density vs. DNA methylation)")
    }
  } else if (category == "uncategorized") {
    if (grepl("Figure", filename)) {
      if (grepl("DNAm_Level", filename)) {
        return("DNA methylation level distribution plot")
      } else if (grepl("DNAm_Variance", filename)) {
        return("DNA methylation variance distribution plot")
      } else if (grepl("Combined", filename)) {
        return("Combined DNA methylation level and variance plot")
      } else if (grepl("Biological", filename)) {
        return("Biological insight plot")
      }
    } else if (grepl("boxplot", filename)) {
      return("Boxplot visualization")
    } else if (grepl("qqplot", filename)) {
      parts <- str_match(filename, "qqplot_(.*?)_(.*?)_(.*?)\\.(png|pdf)")
      if (!is.na(parts[1])) {
        trait <- parts[2]
        population <- parts[3]
        region <- parts[4]
        return(paste0("QQ plot for ", trait, " in ", population, " ", region))
      }
    }
  }
  
  # Default case
  return("Unknown plot type")
}

# Add information column to missing plots
missing_plots$information <- mapply(extract_plot_info, missing_plots$filepath, missing_plots$category)

# Group by category
missing_by_category <- missing_plots %>%
  group_by(category) %>%
  summarize(
    count = n(),
    percentage = round(100 * n() / nrow(missing_plots), 1)
  ) %>%
  arrange(desc(count))

# Generate report
cat("Generating missing plots report...\n")
report <- c(
  "# Missing Plots Report",
  "",
  paste("Generated on:", format(Sys.time())),
  "",
  paste("Total plots:", nrow(all_plots)),
  paste("Missing plots:", nrow(missing_plots), "(", round(100 * nrow(missing_plots) / nrow(all_plots), 1), "%)"),
  paste("Recreation rate:", round(100 * (nrow(all_plots) - nrow(missing_plots)) / nrow(all_plots), 1), "%"),
  "",
  "## Missing Plots by Category",
  ""
)

# Add category summaries
for (i in 1:nrow(missing_by_category)) {
  cat_name <- missing_by_category$category[i]
  cat_count <- missing_by_category$count[i]
  cat_percentage <- missing_by_category$percentage[i]
  
  report <- c(report,
              paste0("### ", cat_name),
              paste0("- **Count**: ", cat_count, " (", cat_percentage, "% of all missing plots)"),
              paste0("- **Recreation Rate**: ", 
                     round(100 * (1 - cat_count / nrow(all_plots[category == cat_name])), 1), "%"),
              "")
}

# Add detailed listings by category
report <- c(report, "## Detailed Missing Plots by Category", "")

for (cat in unique(missing_plots$category)) {
  cat_plots <- missing_plots[category == cat]
  
  report <- c(report,
              paste0("### ", cat, " (", nrow(cat_plots), " missing plots)"),
              "")
  
  # Add table header
  report <- c(report,
              "| Filename | Information | Missing Reason |",
              "|----------|-------------|----------------|")
  
  # Add rows
  for (i in 1:min(nrow(cat_plots), 50)) {  # Limit to 50 plots per category to keep report manageable
    filename <- basename(cat_plots$filepath[i])
    info <- cat_plots$information[i]
    reason <- cat_plots$missing_reason[i]
    
    report <- c(report,
                paste0("| ", filename, " | ", info, " | ", reason, " |"))
  }
  
  # Add note if more than 50 plots
  if (nrow(cat_plots) > 50) {
    report <- c(report,
                paste0("| ... | ... | ... |"),
                paste0("| _(", nrow(cat_plots) - 50, " more plots not shown)_ | | |"))
  }
  
  report <- c(report, "")
}

# Add recommendations
report <- c(report,
            "## Recommendations for Addressing Missing Plots",
            "",
            "### High Priority Categories",
            "",
            "1. **manhattan_qq_fast** - These plots are critical for visualizing GWAS/MWAS results:",
            "   - Check for missing input data files in the source directories",
            "   - Update the recreate_manhattan_qq_plots.R script to handle edge cases",
            "   - Consider using fallback data for plots with missing inputs",
            "",
            "2. **19_plots_fast** - These plots provide trait-specific visualizations:",
            "   - Verify that all input files are available in the correct locations",
            "   - Update the recreate_19_plot.R script to handle all trait/population/region combinations",
            "",
            "3. **43_plots_fast** - These plots show important correlation analyses:",
            "   - Enhance the recreate_43_plot.R script to handle more plot types",
            "   - Check for missing input data and provide fallback options",
            "",
            "### Medium Priority Categories",
            "",
            "4. **50_plots_fast** - These plots show genomic overlaps:",
            "   - Update the recreate_50_plot.R script to handle all feature types",
            "   - Ensure all required input files are available",
            "",
            "5. **dnam_violin_plots_fast** - These plots visualize DNA methylation distributions:",
            "   - Enhance the recreate_violin_plot.R script to handle all plot variations",
            "   - Check for missing input data files",
            "",
            "### Low Priority Categories",
            "",
            "6. **uncategorized** - These plots may include temporary or obsolete visualizations:",
            "   - Review each plot to determine if it's still needed",
            "   - For important plots, create dedicated recreation scripts",
            "   - Update the categorization logic in fast_verify_plots.R to better classify these plots",
            "",
            "## Next Steps",
            "",
            "1. Focus on high-priority categories first",
            "2. For each category, check if the input data files exist",
            "3. Update the recreation scripts to handle edge cases and provide fallback options",
            "4. Re-run the verification process to track progress",
            "5. Document any plots that cannot be recreated due to missing source data")

# Write report to file
writeLines(report, report_file)
cat("Report saved to:", report_file, "\n")

# Also save a CSV file with all missing plots
fwrite(missing_plots, "missing_plots_details.csv")
cat("Detailed CSV saved to: missing_plots_details.csv\n")

cat("\nAnalysis complete!\n")