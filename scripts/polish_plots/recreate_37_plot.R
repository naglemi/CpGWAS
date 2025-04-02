#!/usr/bin/env Rscript

# Recreation script for 37_plot_and_summarize.ipynb plots
# Date: March 29, 2025
# This script recreates the Manhattan and QQ plots from notebook 37

# =========== SETUP ===========
cat("Starting 37_plot recreation script...\n")
start_time <- Sys.time()

# Load necessary libraries
library(data.table)
library(parallel)
library(qqman)
library(ggplot2)
library(stringr)
library(foreach)
library(doParallel)

# Create output directories
output_dir <- "plot_outputs/37_plots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create a log file
log_file <- file.path(output_dir, "recreation_log.txt")
log_conn <- file(log_file, "w")

log_message <- function(msg) {
  cat(paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n"))
  cat(paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n"), file = log_conn)
}

log_message("Starting recreation of 37_plot_and_summarize plots")

# Define variables
desired_subpopulations <- c("AA", "EA", "all")
desired_regions <- c("caud", "hippo", "dlpfc")
traits <- c("scz", "bp", "mdd")

# Load result table if it exists, otherwise create it
result_table_path <- "37-OUT/37-OUT_result_table_no_dups.csv"

if (file.exists(result_table_path)) {
  log_message(paste("Loading existing result table from", result_table_path))
  result_table <- fread(result_table_path)
} else {
  log_message("Result table not found, creating a new one")
  
  # Create scaffolding code from notebook to find the files
  log_message("Finding data files...")
  
  df <- fread("12-OUT_matched_SNP_meth_cov_outputs.csv")
  df <- df[!grepl("empty", df$path), ]
  df$tag <- paste0(df$region, "_", df$population, "_", df$Chr, "_", df$chunk_start, "_", df$chunk_end)
  df$path_unmodified <- df$path
  
  df$path <- gsub("-caud", "", df$path)
  df$path <- gsub("-hippo", "", df$path)
  df$path <- gsub("-dlpfc", "", df$path)
  
  datetime <- str_split_fixed(df$path, "allcorepera-", 2)[,2]
  datetime <- gsub("\\.rds", "", datetime)
  datetime <- str_split_fixed(datetime, "-", 2)
  
  df$date <- as.numeric(datetime[, 1])
  df$time <- as.numeric(datetime[, 2])
  
  df[, datetime := as.POSIXct(paste0(date, " ", time), format = "%Y%m%d %H%M%S")]
  
  # Order the data by datetime in descending order
  setorder(df, region, Chr, population, chunk_start, chunk_end, -datetime)
  
  # Remove duplicate groups, keeping the first (most recent) entry
  df <- df[, .SD[1], by = .(region, Chr, population, chunk_start, chunk_end)]
  
  # Process traits
  for(i in 1:length(traits)) {
    this_trait <- traits[i]
    log_message(paste("Processing trait:", this_trait))
    
    dir1 <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/output_EXPANSE_dt_caud/"
    dir2 <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/output_EXPANSE_dt_dlpfc/"
    dir3 <- "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/output_EXPANSE_dt_hippo/"
    
    data_files.1 <- list.files(dir1, full.names = TRUE, pattern = paste0("gwas_stat_", this_trait, "_alleleprocessed_a3"))
    data_files.2 <- list.files(dir2, full.names = TRUE, pattern = paste0("gwas_stat_", this_trait, "_alleleprocessed_a3"))
    data_files.3 <- list.files(dir3, full.names = TRUE, pattern = paste0("gwas_stat_", this_trait, "_alleleprocessed_a3"))
    data_files <- c(data_files.1, data_files.2, data_files.3)
    
    data_files <- data_files[grepl("101324", data_files)]
    data_files <- data_files[grepl("a3", data_files)]
    data_files <- data_files[grepl("csv", data_files)]
    
    data_files_alt_name <- gsub("-caud", "", data_files)
    data_files_alt_name <- gsub("-dlpfc", "", data_files_alt_name)
    data_files_alt_name <- gsub("-hippo", "", data_files_alt_name)
    
    result_table <- data.table(str_split_fixed(data_files_alt_name, "[-_]", n = Inf))
    result_table[, c(1:3, 5, 8:11, 14:16, 19:20, 22:25)] <- NULL
    
    colnames(result_table) <- c("region", "Chr", "population", "chunk_start", "chunk_end", "date", "time", "trait")
    result_table$chunk_start <- as.numeric(as.character(result_table$chunk_start))
    result_table$chunk_end <- as.numeric(as.character(result_table$chunk_end))
    
    these_paths <- as.data.frame(data_files)
    colnames(these_paths) <- paste0("result_", this_trait)
    result_table <- cbind(result_table, these_paths)
    
    result_table$Chr <- as.numeric(gsub("chr", "", result_table$Chr))
    result_table$region <- gsub("//libd", "", result_table$region)
    result_table$tag <- paste0(result_table$region, "_", result_table$population, "_", result_table$Chr, "_", result_table$chunk_start, "_", result_table$chunk_end, "_", result_table$trait)
    
    # Convert 'date' and 'time' to a single POSIXct datetime column for accurate comparison
    result_table[, datetime := as.POSIXct(paste0(date, " ", time), format = "%Y%m%d %H%M%S")]
    
    # Order the data by datetime in descending order
    setorder(result_table, region, Chr, population, chunk_start, chunk_end, -datetime)
    
    # Remove duplicate groups, keeping the first (most recent) entry
    result_table <- result_table[, .SD[1], by = .(region, Chr, population, chunk_start, chunk_end)]
    
    df <- df[!grepl("empty", df$path), ]
    result_table$trait <- result_table$date <- result_table$time <- result_table$tag <- result_table$datetime <- NULL
    df <- merge(df, result_table, by = c("region", "Chr", "population", "chunk_start", "chunk_end"))
  }
  
  result_table <- df
  dir.create("37-OUT", showWarnings = FALSE)
  fwrite(result_table, "37-OUT/37-OUT_result_table_no_dups.csv")
}

# Detect the number of cores and set up the parallel backend
log_message("Setting up parallel backend...")
num_cores <- min(32, detectCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

log_message(paste("Using", num_cores, "cores for parallel processing"))

# Create a summary table to track progress
summary_table <- data.frame(
  trait = character(),
  population = character(),
  region = character(),
  total_rows = integer(),
  filtered_rows = integer(),
  bonferroni_rows = integer(),
  status = character(),
  stringsAsFactors = FALSE
)

# Loop over subpopulations, regions, and traits
for (desired_subpopulation in desired_subpopulations) {
  for (desired_region in desired_regions) {
    for (trait in traits) {
      log_message(paste("Processing:", trait, desired_subpopulation, desired_region))
      
      tryCatch({
        sub_df <- result_table[which(result_table$population == desired_subpopulation &
                                       result_table$region == desired_region)]
        
        if (nrow(sub_df) == 0) {
          log_message(paste("ERROR: No rows in sub_df for", trait, desired_subpopulation, desired_region))
          summary_table <- rbind(summary_table, data.frame(
            trait = trait,
            population = desired_subpopulation,
            region = desired_region,
            total_rows = 0,
            filtered_rows = 0,
            bonferroni_rows = 0,
            status = "Failed - No data",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        # Prepare variables for parallel processing
        data_files_subset <- sub_df[[paste0("result_", trait)]]
        result_table_chr <- sub_df$Chr
        
        log_message(paste("  Reading", length(data_files_subset), "data files"))
        
        # Parallel data reading and processing
        data_list <- foreach(i = seq_along(data_files_subset), .packages = 'data.table') %dopar% {
          if (file.exists(data_files_subset[i])) {
            dt <- fread(data_files_subset[i], showProgress = FALSE)
            dt[, chr := result_table_chr[i]]
            return(dt)
          } else {
            return(NULL)
          }
        }
        
        # Remove NULL entries
        data_list <- data_list[!sapply(data_list, is.null)]
        
        if (length(data_list) == 0) {
          log_message(paste("ERROR: No data files found for", trait, desired_subpopulation, desired_region))
          summary_table <- rbind(summary_table, data.frame(
            trait = trait,
            population = desired_subpopulation,
            region = desired_region,
            total_rows = 0,
            filtered_rows = 0,
            bonferroni_rows = 0,
            status = "Failed - Files not found",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        # Combine the data tables into one
        combined_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
        log_message(paste("  Combined data rows:", nrow(combined_data)))
        
        # Read and combine correlation files
        cor_files <- str_split_fixed(data_files_subset, "_gwas", 2)[, 1]
        cor_files <- unique(cor_files)
        cor_files <- paste0(cor_files, "_cv-eval_combined.csv")
        cor_files <- cor_files[file.exists(cor_files)]
        
        if (length(cor_files) == 0) {
          log_message(paste("ERROR: No correlation files found for", trait, desired_subpopulation, desired_region))
          summary_table <- rbind(summary_table, data.frame(
            trait = trait,
            population = desired_subpopulation,
            region = desired_region,
            total_rows = nrow(combined_data),
            filtered_rows = 0,
            bonferroni_rows = 0,
            status = "Failed - No correlation files",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        log_message(paste("  Reading", length(cor_files), "correlation files"))
        
        cor_list <- foreach(file = cor_files, .packages = 'data.table') %dopar% {
          if (file.exists(file)) {
            return(fread(file, showProgress = FALSE))
          } else {
            return(NULL)
          }
        }
        
        # Remove NULL entries
        cor_list <- cor_list[!sapply(cor_list, is.null)]
        
        if (length(cor_list) == 0) {
          log_message(paste("ERROR: Failed to read correlation files for", trait, desired_subpopulation, desired_region))
          summary_table <- rbind(summary_table, data.frame(
            trait = trait,
            population = desired_subpopulation,
            region = desired_region,
            total_rows = nrow(combined_data),
            filtered_rows = 0,
            bonferroni_rows = 0,
            status = "Failed - Error reading correlation files",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        combined_cor <- rbindlist(cor_list, use.names = TRUE, fill = TRUE)
        combined_cor <- unique(combined_cor)
        log_message(paste("  Combined correlation rows:", nrow(combined_cor)))
        
        # Ensure consistent column names for merging
        setnames(combined_data, old = "bp", new = "cg", skip_absent = TRUE)
        setnames(combined_cor, old = "V4", new = "cg", skip_absent = TRUE)
        
        # Set keys for efficient merging
        setkey(combined_data, cg)
        setkey(combined_cor, cg)
        
        # Merge datasets
        log_message("  Merging datasets")
        merged <- merge(combined_data, combined_cor, by = c("chr", "cg"))
        
        if (nrow(merged) == 0) {
          log_message(paste("ERROR: No rows after merging for", trait, desired_subpopulation, desired_region))
          summary_table <- rbind(summary_table, data.frame(
            trait = trait,
            population = desired_subpopulation,
            region = desired_region,
            total_rows = nrow(combined_data),
            filtered_rows = 0,
            bonferroni_rows = 0,
            status = "Failed - Empty merge result",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        # Filter merged data
        log_message("  Filtering data")
        merged <- merged[cor >= 0.1]
        log_message(paste("  Filtered rows:", nrow(merged)))
        
        merged_sample <- merged[sample(1:nrow(merged), min(10000, nrow(merged))), ]
        merged_subset <- merged[p <= 0.0001]
        merged_subset[, SNP := paste0(chr, ":", cg)]
        log_message(paste("  Significant rows (p <= 0.0001):", nrow(merged_subset)))
        
        # Bonferroni correction
        bonferroni_threshold <- 0.05 / nrow(merged)
        bonferroni_pass <- merged[p <= bonferroni_threshold]
        log_message(paste("  Bonferroni-significant rows:", nrow(bonferroni_pass)))
        
        # Save filtered and unfiltered data
        log_message("  Saving result files")
        fwrite(merged, file.path(output_dir, paste0("unfiltered_", trait, "_results_a6_", desired_subpopulation, "_", desired_region, ".csv")))
        fwrite(merged_subset, file.path(output_dir, paste0("filtered_", trait, "_results_a6_", desired_subpopulation, "_", desired_region, ".csv")))
        fwrite(bonferroni_pass, file.path(output_dir, paste0("filtered-bonf_", trait, "_results_a6_", desired_subpopulation, "_", desired_region, ".csv")))
        
        # Save Manhattan plot
        log_message("  Creating Manhattan plot")
        png(filename = file.path(output_dir, paste0("manhattan_", trait, "_", desired_subpopulation, "_", desired_region, ".png")), 
            width = 1200, height = 800, res = 150)
        manhattan(merged_subset, chr = "chr", bp = "cg", p = "p", main = paste("Manhattan plot for", trait, desired_subpopulation, desired_region))
        dev.off()
        
        # Create and save QQ plot
        log_message("  Creating QQ plot")
        png(filename = file.path(output_dir, paste0("qqplot_", trait, "_", desired_subpopulation, "_", desired_region, ".png")), 
            width = 1200, height = 800, res = 150)
        qq(merged_sample$p, main = paste("QQ plot for 10k samples from", trait, desired_subpopulation, desired_region))
        dev.off()
        
        # Update summary table with success
        summary_table <- rbind(summary_table, data.frame(
          trait = trait,
          population = desired_subpopulation,
          region = desired_region,
          total_rows = nrow(merged),
          filtered_rows = nrow(merged_subset),
          bonferroni_rows = nrow(bonferroni_pass),
          status = "Success",
          stringsAsFactors = FALSE
        ))
        
        log_message(paste("Completed:", trait, desired_subpopulation, desired_region))
        gc()
      }, error = function(e) {
        log_message(paste("ERROR in", trait, desired_subpopulation, desired_region, ":", e$message))
        summary_table <- rbind(summary_table, data.frame(
          trait = trait,
          population = desired_subpopulation,
          region = desired_region,
          total_rows = NA,
          filtered_rows = NA,
          bonferroni_rows = NA,
          status = paste("Failed -", e$message),
          stringsAsFactors = FALSE
        ))
      })
    }
  }
}

# Stop the cluster
stopCluster(cl)

# Write summary table
fwrite(summary_table, file.path(output_dir, "plot_recreation_summary.csv"))

# Create summary plots
log_message("Creating summary plots")

# Bar plot of recreation success/failure
if (requireNamespace("ggplot2", quietly = TRUE)) {
  status_summary <- as.data.frame(table(summary_table$status))
  colnames(status_summary) <- c("Status", "Count")
  
  png(filename = file.path(output_dir, "recreation_status_summary.png"), 
      width = 1200, height = 800, res = 150)
  print(
    ggplot(status_summary, aes(x = Status, y = Count, fill = Status)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = "Plot Recreation Status", x = "Status", y = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

# End time and total runtime
end_time <- Sys.time()
total_runtime <- difftime(end_time, start_time, units = "mins")
log_message(paste("Script completed in", round(total_runtime, 2), "minutes"))

# Close log file
close(log_conn)

cat("\nRecreation of 37_plot_and_summarize plots completed.\n")
cat("Check", output_dir, "for results and", log_file, "for detailed log.\n")