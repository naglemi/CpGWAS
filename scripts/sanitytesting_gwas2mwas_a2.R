set.seed(42)

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

profvis::profvis({

# Expected result

# cg  z	p
# 73418062	-12.52445	5.487211e-36
# 73418161	-15.61757	5.526944e-55
# 73418186	-15.64082	3.837365e-55
# 73418205	-21.06016	1.845558e-98
# 73418313	-19.28061	7.814446e-83

#rds_files <- list.files("~/code/CpGWAS/data/stage2_testing/stage2_mwas_testing_samples/",
#                        pattern = "rds$", full.names = TRUE)

#rds_files <- rds_files[!grepl("results", rds_files)]

summary_stats_path <- "/Users/michael.nagle/data/gwas_stat_scz"
summary_stats <- suppressWarnings(data.table::fread(summary_stats_path))

summary_stats <- clean_and_standardize_colnames(summary_stats)

# for(rds_file in rds_files[8:length(rds_files)]){

#rds_file <- rds_files[grepl("908982-928981", rds_files)]

#rds_file <- "/Users/michael.nagle/data/libd_chr1-chr1_all-libd_chr1-chr1_all-908982-928981-dynamic-1corestotal-allcorepera-caud-20240510-145818.rds"

rds_file <- "/Users/michael.nagle/data/libd_chr7-chr7_AA-libd_chr7-chr7_AA-40001-60000-dynamic-1corestotal-allcorepera-20240422-091624.rds"

print(rds_file)
my_rds <- readRDS(rds_file)
MWASmodels <- list()

paths <- list(pvar_path = "/Users/michael.nagle/data/ref_EUR_chr7.pvar",
              pgen_path = "/Users/michael.nagle/data/ref_EUR_chr7.pgen",
              psam_path = "/Users/michael.nagle/data/ref_EUR_chr7.psam")

my_SNPs <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
setkey(my_SNPs$pvar_dt, `#CHROM`, POS)

# pb <- progress_bar$new(
#   format = "[:bar] :percent eta: :eta",
#   total = length(my_rds@models), clear = FALSE, width= 60
# )

#desired_sites <- c(73418062, 73418161, 73418186, 73418205, 73418313)

cg_list <- fread("/Users/michael.nagle/data/26-OUT_intermediate_cg_list_chr7.csv")
#head(cg_list$cg)
 
desired_sites <- cg_list$cg

indices_list <- c()

for (i in seq_along(my_rds@models)) {
  this_MethylationPosition <- my_rds@models[[i]]@methylationPosition
  if(this_MethylationPosition >= min(desired_sites) & this_MethylationPosition <= max(desired_sites)){
    indices_list <- c(indices_list, i)
  }
  #pb$tick()
}

my_rds@models <- my_rds@models[indices_list]

for (i in seq_along(my_rds@models)) {
  this_MethylationBase <- my_rds@models[[i]]
  
  if(length(this_MethylationBase@snpWeights) == 0){
    next
  }
  
  MWASmodels[[i]] <- process_model(this_MethylationBase, my_SNPs, summary_stats)
  
  if(is.null(MWASmodels[[i]])){
    next
  }
  
  MWASmodels[[i]]['bp'] <- this_MethylationBase@methylationPosition
  #pb$tick()
}

# Ensure the lengths of my_rds@models and MWASmodels are the same
length(my_rds@models) == length(MWASmodels)


results <- MWASresults(MWASmodels, paths$pvar_path, paths$pgen_path, paths$psam_path, summary_stats_path, "output/libd_chr1-chr1_AA-20240405-131028.rds")

# Filter out non-null models and extract 'mwas_out'
filtered_models <- Filter(Negate(is.null), results@MWASmodels)
df <- do.call(rbind, filtered_models)
df <- as.data.frame(df)

#results
  
# }
})

# Extract p-values from results, skipping NULL models
p_values <- sapply(results@MWASmodels, function(model) {
  if (!is.null(model)) {
    return(model@mwas_out['p'])
  } else {
    return(NA)
  }
})

# Remove NA values
p_values <- na.omit(unlist(p_values))

# Display p-values
p_values

# First run without proper flippage
#> min(p_values)
# [1] 1.299095e-59

# second run unflipped again
#> min(p_values)
#[1] 1.299095e-59

# First run with proper flippage
#> min(p_values)
#[1] 3.386754e-24

# Second run with proper flippage
# > min(p_values)
# [1] 3.386754e-24


