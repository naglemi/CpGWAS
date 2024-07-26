# Script C: script_C.R
library(CpGWAS)
library(data.table)
library(stringr)
library(optparse)

# Command line options
option_list <- list(
  make_option(c("-g", "--genome_file_index"), type = "integer", default = 1,
              help = "Index of genome file to process"),
  make_option(c("-d", "--data_file"), type = "character", default = "/expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/12-OUT_matched_SNP_meth_cov_outputs.csv",
              help = "Path to data file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load genome files
genome_files <- list.files("/expanse/lustre/projects/jhu152/naglemi/mwas/gwas",
                           pattern = "EUR", full.names = TRUE)
genome_files <- genome_files[grepl("pvar", genome_files)]

genome_files <- data.table(path = genome_files, Chr = NA)

genome_files$Chr <- str_split_fixed(genome_files$path, "chr", 2)[, 2]
genome_files$Chr <- gsub(".pvar", "", genome_files$Chr)

genome_files$Chr <- as.integer(genome_files$Chr)
genome_files <- genome_files[order(genome_files$Chr), ]

df <- fread(opt$data_file)

summary_stats_list <- list.files("/expanse/lustre/projects/jhu152/naglemi/mwas/gwas", pattern = "stat", full.names = TRUE)

# Pre-load all summary stats files into a list and clean/standardize column names
summary_stats_data <- lapply(summary_stats_list, function(path) {
  stats <- suppressWarnings(data.table::fread(path))
  colnames(stats) <- gsub("#CHROM", "CHR", colnames(stats))
  clean_and_standardize_colnames(stats)
})

print("Starting genome file processing")
# Process the specified genome file
g <- opt$genome_file_index
print(paste("Processing genome file index:", g))

paths <- list(
  pvar_path = genome_files[g]$path,
  pgen_path = gsub("pvar", "pgen", genome_files[g]$path),
  psam_path = gsub("pvar", "psam", genome_files[g]$path)
)

my_SNPs <- CpGWAS::loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
setkey(my_SNPs$pvar_dt, `#CHROM`, POS)
df_this_chr <- df[which(df$Chr == genome_files[g]$Chr), ]

summary_stats_data <- lapply(summary_stats_data, function(stats) stats[`CHR` == genome_files[g]$Chr])

print("Loaded SNP data")
print("Files for this Chr:")
print(nrow(df_this_chr))
for(j in 1:nrow(df_this_chr)){
  print(paste0("File number: ", j))
  if (grepl("empty", df_this_chr$path[j])) {
    message(paste0("no model for ", df_this_chr$path[j]))
    next
  }

  all_files_exist <- TRUE
  outnames <- vector("character", length(summary_stats_list))

  for (k in 1:length(summary_stats_list)) {
    outnames[k] <- gsub("\\.rds$", paste0("_", basename(tools::file_path_sans_ext(summary_stats_list[[k]])), "_results.rds"), df_this_chr$path[j])
    if (!file.exists(outnames[k])) {
      all_files_exist <- FALSE
      break  # Exit the loop early if any file does not exist
    }
  }

  if (all_files_exist) {
    print(paste("All output files already exist for", df_this_chr$path[j], "- Skipping"))
    next  # Skip to the next file if all outputs already exist
  }
  #
  my_rds <- tryCatch({
    readRDS(df_this_chr$path[j])
  }, error = function(e) {
    # Print an error message and skip this iteration
    message("ALERT!!! Error reading RDS file: ", e$message)
    return(NULL)  # Return NULL to signal failure
  })
    
  # Check if the readRDS call returned NULL (which indicates an error)
  if (is.null(my_rds)) {
    next  # Skip the rest of this loop iteration
  }

  print(paste("Loaded RDS file:", df_this_chr$path[j]))

  for (k in 1:length(summary_stats_list)) {
    outname <- gsub("\\.rds$", paste0("_", basename(tools::file_path_sans_ext(summary_stats_list[[k]])), "_results.rds"), df_this_chr$path[j])
    if(file.exists(outname)) next
    summary_stats <- summary_stats_data[[k]]

    MWASmodels <- vector("list", length(my_rds@models))
    if (is.null(summary_stats)) {
      summary_stats <- suppressWarnings(fread(summary_stats_list[[k]]))
      summary_stats <- clean_and_standardize_colnames(summary_stats)
    }

    for (i in seq_along(my_rds@models)) {
      this_MethylationBase <- my_rds@models[[i]]
      SNP_split <- stringr::str_split_fixed(names(this_MethylationBase@snpWeights), ":", 4)
      SNP_split[, 1] <- gsub("chr", "", SNP_split[, 1])
      SNP_split_dt <- data.table::as.data.table(SNP_split)
      data.table::setnames(SNP_split_dt, c("chr", "post", "ref", "alt"))
      SNP_split_dt[, `:=`(chr = as.integer(chr), post = as.integer(post))]
      data.table::setkey(SNP_split_dt, chr, post)

      relevant_SNP_indices <- my_SNPs$pvar_dt[SNP_split_dt, on = .(`#CHROM` = chr, POS = post), which = TRUE, nomatch = 0]
      relevant_ids <- my_SNPs$pvar_dt$ID[relevant_SNP_indices]
      summary_stats_sub <- summary_stats[relevant_ids, nomatch = 0]

      if (!identical(summary_stats_sub$BP, SNP_split_dt$post)) {
        summary_stats_sub <- summary_stats_sub[order(summary_stats_sub$BP), ]
        if (!identical(summary_stats_sub$BP, SNP_split_dt$post)) {
          unmatched_positions <- !SNP_split_dt$post %in% summary_stats_sub$BP
          if (any(unmatched_positions)) {
            SNP_split_dt <- SNP_split_dt[!unmatched_positions, ]
            this_MethylationBase@snpWeights <- this_MethylationBase@snpWeights[!unmatched_positions]

            relevant_SNP_indices <- my_SNPs$pvar_dt[SNP_split_dt, on = .(`#CHROM` = chr, POS = post), which = TRUE, nomatch = 0]
            if (!identical(summary_stats_sub$BP, SNP_split_dt$post)) {
              stop("SNP order does not match even after removing unmatched positions. This should not happen. Code is broken.")
            }
          }
        }
      }

      if (!identical(SNP_split_dt$alt, summary_stats_sub$A2) | !identical(SNP_split_dt$ref, summary_stats_sub$A1)) {
        not_matching <- which(SNP_split_dt$alt != summary_stats_sub$A2)
        summary_stats_ref_flipped <- SNP_split_dt$ref[not_matching]
        summary_stats_alt_flipped <- SNP_split_dt$alt[not_matching]
        SNP_split_dt[not_matching, `:=`(ref = summary_stats_alt_flipped, alt = summary_stats_ref_flipped)]
        this_MethylationBase@snpWeights[not_matching] <- this_MethylationBase@snpWeights[not_matching] * -1
      }

      G <- pgenlibr::ReadList(my_SNPs$pgen, variant_subset = relevant_SNP_indices)
      #print(paste("Performing MWAS for model index:", i))
      mwas_out <- mwas(z = summary_stats_sub$BETA, w = this_MethylationBase@snpWeights, G = G)

      MWASmodels[[i]] <- mwas_out
    }

    results <- MWASresults(MWASmodels, paths$pvar_path, paths$pgen_path, paths$psam_path, summary_stats_list[[k]], df_this_chr$path[j])
    saveRDS(results, outname)
    print(paste("Saved results to:", outname))
  }
}
