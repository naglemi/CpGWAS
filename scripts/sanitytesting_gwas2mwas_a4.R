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
# #. select and load firwt set of files
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

summary_stats_list <- list.files("~/data/gwas_flipped/", pattern = "stat", full.names = TRUE)
summary_stats_list <- summary_stats_list[grepl("alleleprocessed_a2", summary_stats_list)]

# Pre-load all summary stats files into a list and clean/standardize column names
summary_stats_data <- lapply(summary_stats_list, function(path) {
   stats <- suppressWarnings(data.table::fread(path))
   setkey(stats, SNP)
   stats
   #colnames(stats) <- gsub("#CHROM", "CHR", colnames(stats))
   #clean_and_standardize_colnames(stats)
 })

# saveRDS(summary_stats_data, "~/data/stats_bp_mdd_scz.rds", compress = FALSE)

#Sys.time()
#summary_stats_data <- readRDS("~/data/stats_bp_mdd_scz.rds")
#Sys.time()

#})
# for(rds_file in rds_files[8:length(rds_files)]){

#rds_file <- rds_files[grepl("908982-928981", rds_files)]

#rds_file <- "/Users/michael.nagle/data/libd_chr1-chr1_all-libd_chr1-chr1_all-908982-928981-dynamic-1corestotal-allcorepera-caud-20240510-145818.rds"

#profvis::profvis({

summary_stats_path <- summary_stats_list[[3]]
summary_stats <- summary_stats_data[[3]]

csv_file <- "/Users/michael-non-admin/data/libd_chr7-chr7_AA-libd_chr7-chr7_AA-40001-60000-dynamic-1corestotal-allcorepera-20240422-091624_dt_combined.csv"


my_csv <- fread(csv_file)
setkey(my_csv, cg)
# if(all(colnames(my_csv)[2:3] == c("features", "cg"))){
#   colnames(my_csv)[2:3] <- c("cg", "features")
# }

head(my_csv)

MWASmodels <- list()

paths <- list(pvar_path = "/Users/michael-non-admin/data/gwas_flipped/ref_EUR_allele-match-a2-libd_chr7.pvar",
              pgen_path = "/Users/michael-non-admin/data/gwas_flipped/ref_EUR_allele-match-a2-libd_chr7.pgen",
              psam_path = "/Users/michael-non-admin/data/gwas_flipped/ref_EUR_allele-match-a2-libd_chr7.psam")

chr <- as.numeric(gsub("\\.pvar", "",
            stringr::str_split_fixed(paths$pvar_path, "chr", 2)[, 2]))

summary_stats <- summary_stats[which(summary_stats$CHR == chr), ]

my_SNPs <- loadSNPData(paths$pvar_path, paths$pgen_path, paths$psam_path)
setkey(my_SNPs$pvar_dt, `#CHROM`, POS)

# pb <- progress_bar$new(
#   format = "[:bar] :percent eta: :eta",
#   total = length(my_rds@models), clear = FALSE, width= 60
# )

#desired_sites <- c(73418062, 73418161, 73418186, 73418205, 73418313)

cg_list <- fread("/Users/michael-non-admin/data/26-OUT_intermediate_cg_list_chr7.csv")

#head(cg_list$cg)
 
desired_sites <- cg_list$cg

indices_list <- c()

cgs <- levels(factor(my_csv$cg))

#my_csv <- my_csv[which(my_csv$cg >= min(desired_sites) & cg <= max(desired_sites)), ]

# test for one with plenty of SNPs
# i <- 340
# this_cg <- cgs[i]
# 
# this_cg_data <- my_csv[cg == this_cg]
# 
# this_MethylationBase <- "oopppopo"
# MWASmodels[[i]] <- process_model_csv(this_cg_data, my_SNPs, summary_stats)
# 
# if(is.null(MWASmodels[[i]])){
#   next
# }
# 
# MWASmodels[[i]]['bp'] <- this_cg
#pb$tick()

length(cgs)

Sys.time()

for (i in seq_along(cgs)) {
  this_cg <- cgs[i]
  
  this_cg_data <- my_csv[cg == this_cg]

  MWASmodels[[i]] <- process_model_csv(this_cg_data, my_SNPs, summary_stats)

  if(is.null(MWASmodels[[i]])){
    next
  }

  MWASmodels[[i]]['bp'] <- this_cg
  #pb$tick()
}

Sys.time()

# MWASmodels <- foreach(this_MethylationBase = my_rds@models, .packages = c("data.table", "CpGWAS"), .combine = 'c', .export = c("my_SNPs", "summary_stats")) %dopar% {
#   
#   if(length(this_MethylationBase@snpWeights) == 0) {
#     return(NULL)
#   }
#   
#   model <- process_model(this_MethylationBase, my_SNPs, summary_stats)
#   
#   if(is.null(model)) {
#     return(NULL)
#   }
#   
#   model['bp'] <- this_MethylationBase@methylationPosition
#   return(list(model))  # Return as a single list
# }

# Ensure the lengths of my_rds@models and MWASmodels are the same
length(my_rds@models) == length(MWASmodels)


results <- MWASresults(MWASmodels, paths$pvar_path, paths$pgen_path, paths$psam_path, summary_stats_path, "output/libd_chr1-chr1_AA-20240405-131028.rds")

# Filter out non-null models and extract 'mwas_out'
filtered_models <- Filter(Negate(is.null), results@MWASmodels)
df <- do.call(rbind, filtered_models)
df <- as.data.frame(df)

#results
  
# }
#})

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


