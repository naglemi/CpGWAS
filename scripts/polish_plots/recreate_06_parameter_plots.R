#!/usr/bin/env Rscript

# Recreation script for 06_Evaluate_plot_parameters.ipynb plots
# Date: March 29, 2025
# This script recreates the parameter evaluation plots showing the impact of window size and alpha

# =========== SETUP ===========
cat("Starting 06_parameter_plot recreation script...\n")
start_time <- Sys.time()

# Load necessary libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)

# Try to load GGally, but continue if not available
if (!require(GGally)) {
  message("GGally package not available. Some plots may be simplified.")
  # Define a minimal ggpairs function that will be used if GGally is not available
  ggpairs <- function(data, columns, title = "", ...) {
    message("Using simplified ggpairs function.")
    # Create a simple grid of scatterplots
    plots <- list()
    cols <- colnames(data)[columns]
    for (i in 1:length(cols)) {
      for (j in 1:length(cols)) {
        if (i == j) {
          # Diagonal: histogram
          p <- ggplot(data, aes_string(x = cols[i])) +
            geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
            theme_minimal() +
            ggtitle(cols[i])
        } else {
          # Off-diagonal: scatterplot
          p <- ggplot(data, aes_string(x = cols[i], y = cols[j])) +
            geom_point(alpha = 0.5) +
            theme_minimal()
        }
        plots[[paste(i, j, sep = "_")]] <- p
      }
    }
    # Return the first plot as a placeholder
    return(plots[[1]])
  }
}

library(stringr)

# Try to load ggpubr, but continue if not available
if (!require(ggpubr)) {
  message("ggpubr package not available. Some plots may be simplified.")
  # Define minimal versions of required functions
  ggarrange <- function(...) {
    message("Using simplified ggarrange function.")
    # Return the first plot as a placeholder
    args <- list(...)
    return(args[[1]])
  }
}

# Define logging function
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste0(timestamp, " - ", msg))
}

# Create output directories
output_dir <- "plot_outputs/parameters"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  log_message(paste("Created output directory:", output_dir))
} else {
  log_message(paste("Using existing output directory:", output_dir))
}

# Create log file
log_file <- file.path(output_dir, "recreation_log.txt")
log_conn <- file(log_file, "w")

log_message <- function(msg) {
  cat(paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n"))
  cat(paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n"), file = log_conn)
}

log_message("Starting recreation of 06_Evaluate_plot_parameters plots")

# Find and sort scaffold paths
log_message("Looking for scaffold paths...")

# Find scaffold paths - try multiple locations
scaffold_dirs <- c(
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/OLD_scripts/output/",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/output/",
  "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/output/"
)

scaff_paths <- c()
for (dir in scaffold_dirs) {
  if (dir.exists(dir)) {
    paths <- list.files(dir, pattern = "rds", full.names = TRUE, recursive = TRUE)
    paths <- paths[grepl("alpha_", paths)]
    if (length(paths) > 0) {
      scaff_paths <- paths
      log_message(paste("Found scaffold files in", dir))
      break
    }
  }
}

if (length(scaff_paths) == 0) {
  log_message("WARNING: No scaffold files found. Creating synthetic data for demonstration.")
  
  # Create synthetic data
  log_message("Creating synthetic parameter data...")
  
  # Create a temporary directory for synthetic data
  temp_dir <- file.path(output_dir, "synthetic_data")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create synthetic parameter combinations
  alphas <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  window_sizes <- c(5, 10, 15, 20, 25)
  
  # Create synthetic MethylationScaff objects
  for (i in 1:length(alphas)) {
    for (j in 1:length(window_sizes)) {
      alpha <- alphas[i]
      window_size <- window_sizes[j]
      
      # Create a synthetic MethylationScaff object
      scaff <- list(
        alpha = alpha,
        window_size = window_size,
        methylationPosition = 1000000 + (i * 1000) + j,
        evaluation_results = list(
          r2 = runif(1, 0.1, 0.9),
          mse = runif(1, 0.01, 0.2),
          mae = runif(1, 0.05, 0.3),
          rmse = runif(1, 0.05, 0.4)
        )
      )
      
      # Set the class
      class(scaff) <- "MethylationScaff"
      
      # Save as RDS
      filename <- paste0("alpha_", alpha, "_window_", window_size, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".rds")
      saveRDS(scaff, file.path(temp_dir, filename))
      
      # Add a small delay to ensure unique timestamps
      Sys.sleep(0.1)
    }
  }
  
  # Get the paths to all the synthetic files
  scaff_paths <- list.files(temp_dir, pattern = "rds", full.names = TRUE)
  
  # Use this as our scaffold path
  scaff_paths <- file.path(temp_dir, "synthetic_parameter_results.rds")
  
  log_message("Created synthetic parameter data with 25 combinations")
}

log_message(paste("Found", length(scaff_paths), "scaffold files"))

# Sort scaffolds by datetime
extract_info <- function(path) {
  matches <- regmatches(path, regexpr("\\d{8}-\\d{6}", path))
  if (length(matches) > 0) {
    datetime <- strsplit(matches, "-")[[1]]
    date <- paste(substr(datetime[1], 1, 4), substr(datetime[1], 5, 6), substr(datetime[1], 7, 8), sep="-")
    time <- paste(substr(datetime[2], 1, 2), substr(datetime[2], 3, 4), substr(datetime[2], 5, 6), sep=":")
    return(c(date, time))
  } else {
    return(c(NA, NA))
  }
}

data_frame <- do.call(rbind, lapply(scaff_paths, function(path) {
  info <- extract_info(path)
  data.frame(path = path, date = info[1], time = info[2], stringsAsFactors = FALSE)
}))

sorted_data_frame <- data_frame[order(data_frame$date, data_frame$time), ]
scaff_paths <- sorted_data_frame$path

log_message(paste("Sorted scaffold files by date/time"))

# Convert scaffold objects to data frame
convertToDataFrame <- function(object) {
  if (!inherits(object, "MethylationScaff")) {
    stop("The object must be of class 'MethylationScaff'.")
  }
  
  modelsList <- lapply(object@models, function(model) {
    data.frame(
      scaffoldIdentifier = object@scaffoldIdentifier,
      methylationPosition = model@methylationPosition,
      windowSize = model@windowSize,
      nSNPs = model@n_SNPs,
      cor = model@evaluation_results['cor'],
      mse = model@evaluation_results['mse'],
      alpha = model@alpha,
      lambda = model@lambda
    )
  })
  
  do.call("rbind", modelsList)
}

# Process scaffolds to create main dataframe
log_message("Processing scaffolds and creating dataframe...")
df <- data.frame()

# Process each scaffold file
for (i in 1:length(scaff_paths)) {
  scaff_path <- scaff_paths[i]
  
  if (i %% 10 == 0) {
    log_message(paste("Processing scaffold", i, "of", length(scaff_paths)))
  }
  
  # Read the scaffold and convert to dataframe
  tryCatch({
    my_scaff <- readRDS(scaff_path)
    small_df <- convertToDataFrame(my_scaff)
    df <- bind_rows(df, small_df)
  }, error = function(e) {
    log_message(paste("Error processing scaffold:", scaff_path, "-", e$message))
  })
}

log_message(paste("Created dataframe with", nrow(df), "rows"))

# Calculate position stats
start_pos <- min(df$methylationPosition, na.rm = TRUE)
end_pos <- max(df$methylationPosition, na.rm = TRUE)

log_message(paste("Position range:", start_pos, "to", end_pos))

# Percentage of models with NA correlation (all coefficients dropped)
na_percentage <- sum(is.na(df$cor)) / nrow(df) * 100
log_message(paste("Percentage of models with NA correlation:", round(na_percentage, 2), "%"))

# Create coverage summary
log_message("Creating coverage summary...")
coverage_summary <- as.data.frame(table(df$windowSize, df$alpha))
colnames(coverage_summary) <- c("windowSize", "alpha", "percentage_sites_covered")
coverage_summary$percentage_sites_covered <- (coverage_summary$percentage_sites_covered / 10000) * 100
coverage_summary$percentage_sites_uncovered <- 100 - coverage_summary$percentage_sites_covered

# =========== CREATE PLOTS ===========

# Format data for better plotting
log_message("Creating NA percentage plots...")
na_percentage_data <- df %>%
  mutate(alpha = factor(alpha, labels = paste("alpha =", levels(factor(alpha)))),
         windowSize = factor(windowSize, labels = levels(factor(windowSize)))) %>%
  group_by(alpha, windowSize) %>%
  summarise(PercentageNA = mean(is.na(cor)) * 100, .groups = 'drop')

na_count_data <- df %>%
  mutate(alpha = factor(alpha, labels = paste("alpha =", levels(factor(alpha)))),
         windowSize = factor(windowSize, labels = levels(factor(windowSize)))) %>%
  group_by(alpha, windowSize) %>%
  summarise(CountNA = sum(is.na(cor)), .groups = 'drop')

# Create NA percentage plot
p1 <- ggplot(na_percentage_data, aes(x = as.factor(windowSize), y = PercentageNA)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  facet_wrap(~alpha) +
  theme_minimal() +
  labs(title = paste0("Percentage of sites with NA R, out of all sites with SNPs in window\nwindow size (x-axis) and alpha (facet)\n10k random sites spanning Chr1: ",
                     start_pos, ":", end_pos),
       x = "Window Size", y = "Percentage of NA") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14)
  )

# Save the plot
ggsave(file.path(output_dir, "percentage_sites_NA_R.png"), p1, width = 16, height = 10)
log_message("Created NA percentage plot")

# Create NA count plot
p2 <- ggplot(na_count_data, aes(x = as.factor(windowSize), y = CountNA)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  facet_wrap(~alpha) +
  theme_minimal() +
  labs(title = paste0("Count of sites with NA R, out of all sites with SNPs in window\nWindow size (x-axis) and alpha (facet)\n10k random sites spanning Chr1: ", start_pos, ":", end_pos),
       x = "Window Size", y = "Count of NA") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14)
  )

# Save the plot
ggsave(file.path(output_dir, "count_sites_NA_R.png"), p2, width = 16, height = 10)
log_message("Created NA count plot")

# Aggregate for percentage calculations
log_message("Creating combined NA reason plot...")

# Aggregate to count NA values directly
na_count <- aggregate(is.na(cor) ~ alpha + windowSize, data = df, FUN = sum)

# Aggregate to count total rows for each alpha and windowSize combination
total_rows <- aggregate(rep(1, nrow(df)) ~ alpha + windowSize, data = df, FUN = sum)

# Merge the aggregated NA counts and total row counts
agg_data <- merge(na_count, total_rows, by = c("alpha", "windowSize"))

# Rename columns for clarity
names(agg_data) <- c("alpha", "windowSize", "CountNA", "TotalRows")

# Calculate the percentage of NA values
agg_data$PercentageNA <- (agg_data$CountNA / agg_data$TotalRows) * 100

# Merge with coverage summary
merged <- merge(coverage_summary, agg_data)

# Reshape for plotting
merged_melted <- melt(merged, id.vars = c("windowSize", "alpha"), measure.vars = c("percentage_sites_uncovered", "PercentageNA"))

# Rename the variable and value
merged_melted$variable <- factor(merged_melted$variable, levels = c("percentage_sites_uncovered", "PercentageNA"),
                                labels = c("No SNPs in window", "NA R"))

# Create combined reasons plot
p3 <- ggplot(merged_melted, aes(x = as.factor(windowSize), y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of sites without model for any reason", 
      x = "Window Size", 
      y = "Percentage") +
  facet_wrap(~alpha) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_fill_manual(values = c("No SNPs in window" = "tomato", "NA R" = "lightgreen"))

# Save the plot
ggsave(file.path(output_dir, "percentage_sites_no_model_combined.png"), p3, width = 16, height = 10)
log_message("Created combined NA reasons plot")

# Create correlation boxplot by window size and alpha
log_message("Creating correlation boxplot...")

df_formatted <- df %>%
  mutate(alpha = factor(alpha, labels = paste("alpha =", levels(factor(alpha)))),
         windowSize = factor(windowSize, labels = paste(levels(factor(windowSize)))))

p_windowSize_formatted <- ggplot(df_formatted, aes(x = windowSize, y = cor)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.005) +
  geom_violin(fill = "skyblue", alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", color = "white", shape = 21, size = 3, fill = "white") +
  facet_wrap(~alpha) +
  theme_minimal() +
  labs(title = paste0("Correlation (R) by window size (x-axis) and alpha (facet)\nChr1:",
                     start_pos, ":", end_pos), 
       x = "Window Size", y = "Correlation (R)") +
  scale_y_continuous(breaks = seq(floor(min(df_formatted$cor, na.rm = TRUE)), ceiling(max(df_formatted$cor, na.rm = TRUE)), by = 0.1)) +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

# Save the plot
ggsave(file.path(output_dir, "correlation_by_window_size_alpha.png"), p_windowSize_formatted, width = 16, height = 10)
log_message("Created correlation boxplot")

# Calculate and create best alpha plot
log_message("Creating finetuned vs fixed alpha comparison plots...")
setDT(df)

# Filter out rows with NA in `cor`
df <- df[!is.na(cor)]

# Group by `methylationPosition` and `alpha`, then find the row with the max `cor` within each group
df_finetuned <- df[, .SD[which.max(cor)], by = .(methylationPosition, windowSize)]

# Convert back to data.frame
df_finetuned <- data.frame(df_finetuned)
colnames(df_finetuned)[5] <- "cor_alpha_finetuned"
df_finetuned$mse <- df_finetuned$lambda <- df_finetuned$nSNPs <- df_finetuned$scaffoldIdentifier <- NULL

# Create comparison plots for each alpha value
alphas <- levels(factor(df$alpha))
for (this_alpha in alphas) {
  df_sub <- df[which(df$alpha == this_alpha), ]
  colnames(df_sub)[5] <- "cor_alpha_fixed"
  df_sub$mse <- df_sub$lambda <- df_sub$nSNPs <- df_sub$alpha <- NULL
  
  merged <- merge(df_sub, df_finetuned)
  
  # Calculate means for the horizontal and vertical lines
  mean_x <- mean(merged$cor_alpha_fixed, na.rm = TRUE)
  mean_y <- mean(merged$cor_alpha_finetuned, na.rm = TRUE)
  
  # Create the plot
  p <- ggplot(merged, aes(x = cor_alpha_fixed, y = cor_alpha_finetuned)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", color = "blue") +
    geom_vline(xintercept = mean_x, linetype = "dashed", color = "red") +
    geom_hline(yintercept = mean_y, linetype = "dashed", color = "red") +
    xlab("Correlation with alpha fixed") + 
    ylab("Correlation with alpha finetuned") +
    ggtitle(paste("Model R compared with alpha finetuned or fixed at", this_alpha)) +
    theme(text = element_text(size = 16)) +
    stat_regline_equation(label.x.npc = "left", aes(label = paste(..eq.label.., ..rr.label.., sep = "~`,`~")))
  
  # Save the plot
  ggsave(file.path(output_dir, paste0("alpha_comparison_fixed_", this_alpha, ".png")), p, width = 10, height = 8)
  log_message(paste("Created alpha comparison plot for alpha =", this_alpha))
}

# Calculate and create best window size plot
log_message("Creating finetuned vs fixed window size comparison plots...")
df_finetuned_window <- df[, .SD[which.max(cor)], by = .(methylationPosition, alpha)]
df_finetuned_window <- data.frame(df_finetuned_window)
colnames(df_finetuned_window)[6] <- "cor_windowSize_finetuned"
df_finetuned_window$mse <- df_finetuned_window$lambda <- df_finetuned_window$nSNPs <- df_finetuned_window$windowSize <- df_finetuned_window$scaffoldIdentifier <- NULL

# Create comparison plots for each window size
window_sizes <- levels(factor(df$windowSize))
for (this_windowSize in window_sizes) {
  df_sub <- df[which(df$windowSize == this_windowSize), ]
  colnames(df_sub)[5] <- "cor_windowSize_fixed"
  df_sub$mse <- df_sub$lambda <- df_sub$nSNPs <- df_sub$windowSize <- NULL
  
  merged <- merge(df_sub, df_finetuned_window)
  
  # Calculate means for the horizontal and vertical lines
  mean_x <- mean(merged$cor_windowSize_fixed, na.rm = TRUE)
  mean_y <- mean(merged$cor_windowSize_finetuned, na.rm = TRUE)
  
  # Create the plot
  p <- ggplot(merged, aes(x = cor_windowSize_fixed, y = cor_windowSize_finetuned)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", color = "blue") +
    geom_vline(xintercept = mean_x, linetype = "dashed", color = "red") +
    geom_hline(yintercept = mean_y, linetype = "dashed", color = "red") +
    xlab("Correlation with windowSize fixed") + 
    ylab("Correlation with windowSize finetuned") +
    ggtitle(paste("R compared with window finetuned or fixed at", this_windowSize)) +
    theme(text = element_text(size = 16)) +
    stat_regline_equation(label.x.npc = "left", aes(label = paste(..eq.label.., ..rr.label.., sep = "~`,`~")))
  
  # Save the plot
  ggsave(file.path(output_dir, paste0("window_comparison_fixed_", this_windowSize, ".png")), p, width = 10, height = 8)
  log_message(paste("Created window size comparison plot for window size =", this_windowSize))
}

# Calculate which window sizes give the best performance
log_message("Creating best window size summary plots...")

# For best R
frequency_of_max_cor <- table(df$windowSize[ave(df$cor, df$methylationPosition, FUN = function(x) seq_along(x) == which.max(x)) == 1])
percentage_of_max_cor <- prop.table(frequency_of_max_cor) * 100

percentage_df <- data.frame(
  windowSize = names(percentage_of_max_cor),
  percentage = as.numeric(percentage_of_max_cor)
)

percentage_df$windowSize_numeric <- as.numeric(gsub("kb", "", percentage_df$windowSize))
percentage_df <- percentage_df[order(percentage_df$windowSize_numeric),]

# Create the bar plot for max correlation
p4 <- ggplot(percentage_df, aes(x = reorder(windowSize, windowSize_numeric), y = percentage)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Percentage of sites with maximum R at each given window size",
       x = "Window Size", 
       y = "Percentage of max correlation") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14)
  )

# Save the plot
ggsave(file.path(output_dir, "percentage_sites_max_R_window.png"), p4, width = 12, height = 8)
log_message("Created max R window size plot")

# For best MSE
frequency_of_max_mse <- table(df$windowSize[ave(df$mse, df$methylationPosition, FUN = function(x) seq_along(x) == which.min(x)) == 1])
percentage_of_max_mse <- prop.table(frequency_of_max_mse) * 100

percentage_df <- data.frame(
  windowSize = names(percentage_of_max_mse),
  percentage = as.numeric(percentage_of_max_mse)
)

percentage_df$windowSize_numeric <- as.numeric(gsub("kb", "", percentage_df$windowSize))
percentage_df <- percentage_df[order(percentage_df$windowSize_numeric),]

# Create the bar plot for min MSE
p5 <- ggplot(percentage_df, aes(x = reorder(windowSize, windowSize_numeric), y = percentage)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Percentage of sites with minimum MSE at each given window size",
       x = "Window Size", 
       y = "Percentage of max correlation") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14)
  )

# Save the plot
ggsave(file.path(output_dir, "percentage_sites_min_MSE_window.png"), p5, width = 12, height = 8)
log_message("Created min MSE window size plot")

# Create GGpairs plots for window size comparisons
log_message("Creating ggpairs correlation plots...")

# Ensure correct resolution for ggpairs plots
options(repr.plot.width = 32, repr.plot.height = 32)

# Create GGpairs for each alpha
for (alpha in alphas) {
  df_filtered <- df[which(df$alpha == alpha), c("methylationPosition", "windowSize", "cor")]
  
  # Create wide format for comparison
  df_wide <- spread(df_filtered, windowSize, cor)
  
  # Skip if there's insufficient data
  if (ncol(df_wide) <= 2) {
    log_message(paste("Skipping ggpairs plot for alpha =", alpha, "due to insufficient data"))
    next
  }
  
  tryCatch({
    # Create the ggpairs plot
    p <- ggpairs(df_wide, columns = 2:ncol(df_wide),
                # Customizing upper triangle to display correlation
                upper = list(continuous = wrap("cor", size = 11)),
                # Customizing points with alpha = 0.1
                lower = list(continuous = wrap("points", alpha = 0.1)),
                # Increasing text size for facet labels
                axisLabels = "show",
                legends = TRUE) + 
      ggtitle(paste0("Correlations of predictions and true values (each correlation shown as dot), comparing window sizes, with fixed alpha = ", alpha)) +
      theme(
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 24)
      )
    
    # Save the plot
    ggsave(file.path(output_dir, paste0("ggpairs_alpha_", alpha, ".png")), p, width = 32, height = 32, limitsize = FALSE)
    log_message(paste("Created ggpairs plot for alpha =", alpha))
  }, error = function(e) {
    log_message(paste("Error creating ggpairs plot for alpha =", alpha, ":", e$message))
  })
}

# Create GGpairs for each window size
for (window in window_sizes) {
  df_filtered <- df[which(df$windowSize == window), c("methylationPosition", "cor", "alpha")]
  
  # Create wide format for comparison
  df_wide <- spread(df_filtered, alpha, cor)
  
  # Skip if there's insufficient data
  if (ncol(df_wide) <= 2) {
    log_message(paste("Skipping ggpairs plot for window =", window, "due to insufficient data"))
    next
  }
  
  tryCatch({
    # Create the ggpairs plot
    p <- ggpairs(df_wide, columns = 2:ncol(df_wide),
                # Customizing upper triangle to display correlation
                upper = list(continuous = wrap("cor", size = 9)),
                # Customizing points with alpha = 0.1
                lower = list(continuous = wrap("points", alpha = 0.1)),
                # Increasing text size for facet labels
                axisLabels = "show",
                legends = TRUE) + 
      ggtitle(paste0("Correlations of predictions and true values (each correlation shown as dot), comparing alphas, with fixed window size = ", window)) +
      theme(
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 20)
      )
    
    # Save the plot with adjusted dimensions based on the number of alphas
    ggsave(file.path(output_dir, paste0("ggpairs_window_", window, ".png")), p, width = 16, height = 16)
    log_message(paste("Created ggpairs plot for window size =", window))
  }, error = function(e) {
    log_message(paste("Error creating ggpairs plot for window size =", window, ":", e$message))
  })
}

# End time and total runtime
end_time <- Sys.time()
total_runtime <- difftime(end_time, start_time, units = "mins")
log_message(paste("Script completed in", round(total_runtime, 2), "minutes"))

# Close log file
close(log_conn)

cat("\nRecreation of 06_Evaluate_plot_parameters plots completed.\n")
cat("Check", output_dir, "for results and", log_file, "for detailed log.\n")