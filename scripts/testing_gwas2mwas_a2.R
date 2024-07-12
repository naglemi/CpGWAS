library(CpGWAS)
library(data.table)

# rds_path <- "output/libd_chr1-chr1_AA-20240405-131028.rds"
# 
# my_rds <- readRDS(rds_path)
# MWASmodels <- list()
# 
# paths <- list(pvar_path = "/Users/michael.nagle/data/ref_EUR_chr1.pvar",
#               pgen_path = "/Users/michael.nagle/data/ref_EUR_chr1.pgen",
#               psam_path = "/Users/michael.nagle/data/ref_EUR_chr1.psam")
# 
# my_SNPs <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
# 
# summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_bp"
# #summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_scz"
# #summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_mdd"
# #summary_stats <- data.table::fread(summary_stats_path)
# 
# #summary_stats <- clean_and_standardize_colnames(summary_stats)
# 
# # Loop over chromosome genome files (pvar/pgen/psam)
# #  make list of chromosome files
# #. levels factor
# #. select and load first set of files
# #  # subset big file-matching df to those for the chromosome of interest
# #. loop over those, and for each....
# ##Loop over summary stat files
# ### Loop over RDS files containing our MethylationBase objects with SNP->CpG models
# results <- process_MWAS_models(my_rds = my_rds, my_SNPs = my_SNPs, paths = paths,
#                                summary_stats_path = summary_stats_path,
#                                rds_path = rds_path)
# 
# # save results object as RDS with filename same as results@rds_path but with _results appended
# saveRDS(results, 
#         gsub("\\.rds$", 
#              paste0("_", basename(tools::file_path_sans_ext(results@summary_stats_path)), "_results.rds"), 
#              results@rds_path))


#results

#profvis::profvis({

rds_files <- list.files("~/code/CpGWAS/data/stage2_testing/stage2_mwas_testing_samples/",
                        pattern = "rds$", full.names = TRUE)

rds_files <- rds_files[!grepl("results", rds_files)]

summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_scz"
summary_stats <- suppressWarnings(data.table::fread(summary_stats_path))

summary_stats <- clean_and_standardize_colnames(summary_stats)

for(rds_file in rds_files[8:length(rds_files)]){
  print(rds_file)
  my_rds <- readRDS(rds_file)
  MWASmodels <- list()
  
  paths <- list(pvar_path = "/Users/michael.nagle/data/ref_EUR_chr1.pvar",
                pgen_path = "/Users/michael.nagle/data/ref_EUR_chr1.pgen",
                psam_path = "/Users/michael.nagle/data/ref_EUR_chr1.psam")
  
  my_SNPs <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
  setkey(my_SNPs$pvar_dt, `#CHROM`, POS)
  
  # pb <- progress_bar$new(
  #   format = "[:bar] :percent eta: :eta",
  #   total = length(my_rds@models), clear = FALSE, width= 60
  # )
  
  for (i in seq_along(my_rds@models)) {
    this_MethylationBase <- my_rds@models[[i]]
    MWASmodels[[i]] <- process_model(this_MethylationBase, my_SNPs, summary_stats)
    #pb$tick()
  }
  
  # Ensure the lengths of my_rds@models and MWASmodels are the same
  length(my_rds@models) == length(MWASmodels)
  
  results <- MWASresults(MWASmodels, paths$pvar_path, paths$pgen_path, paths$psam_path, summary_stats_path, "output/libd_chr1-chr1_AA-20240405-131028.rds")
  
  #results
  
}


#})
