library(CpGWAS)

my_rds <- readRDS("output/libd_chr1-chr1_AA-20240405-131028.rds")
MWASmodels <- list()

paths <- list(pvar_path = "/Users/michael.nagle/data/ref_EUR_chr1.pvar",
              pgen_path = "/Users/michael.nagle/data/ref_EUR_chr1.pgen",
              psam_path = "/Users/michael.nagle/data/ref_EUR_chr1.psam")

my_SNPs <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)

summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_bp"
summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_scz"
#summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_mdd"
#summary_stats <- data.table::fread(summary_stats_path)

#summary_stats <- clean_and_standardize_colnames(summary_stats)

# Loop over chromosome genome files (pvar/pgen/psam)
#  make list of chromosome files
#. levels factor
#. select and load first set of files
#  # subset big file-matching df to those for the chromosome of interest
#. loop over those, and for each....
##Loop over summary stat files
### Loop over RDS files containing our MethylationBase objects with SNP->CpG models
results <- process_MWAS_models(my_rds = my_rds, my_SNPs = my_SNPs, paths = paths,
                               summary_stats_path = summary_stats_path,
                               "output/libd_chr1-chr1_AA-20240405-131028.rds")

# save results object as RDS with filename same as results@rds_path but with _results appended
saveRDS(results, 
        gsub("\\.rds$", 
             paste0("_", basename(tools::file_path_sans_ext(results@summary_stats_path)), "_results.rds"), 
             results@rds_path))


#results

# profvis::profvis({
#   my_rds <- readRDS("output/libd_chr1-chr1_AA-20240405-131028.rds")
#   MWASmodels <- list()
#   
#   paths <- list(pvar_path = "/Users/michael.nagle/data/ref_EUR_chr1.pvar",
#                 pgen_path = "/Users/michael.nagle/data/ref_EUR_chr1.pgen",
#                 psam_path = "/Users/michael.nagle/data/ref_EUR_chr1.psam")
#   
#   my_SNPs <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
#   
#   summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_scz"
#   summary_stats <- data.table::fread(summary_stats_path)
#   
#   # because for some reason header is tab-delimited while everything else has spaces
#   real_colnames <- stringr::str_split(colnames(summary_stats)[1], "\t")[[1]]
#   colnames(summary_stats) <- real_colnames
#   
#   pb <- progress_bar$new(
#     format = "[:bar] :percent eta: :eta",
#     total = length(my_rds@models), clear = FALSE, width= 60
#   )
#   
#   for (i in seq_along(my_rds@models)) {
#     this_MethylationBase <- my_rds@models[[i]]
#     MWASmodels[[i]] <- process_model(this_MethylationBase, my_SNPs, summary_stats)
#     pb$tick()
#   }
#   
#   # Ensure the lengths of my_rds@models and MWASmodels are the same
#   length(my_rds@models) == length(MWASmodels)
#   
#   results <- MWASresults(MWASmodels, paths$pvar_path, paths$pgen_path, paths$psam_path, summary_stats_path, "output/libd_chr1-chr1_AA-20240405-131028.rds")
#   
#   results
# })
