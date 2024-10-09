#!/usr/bin/env Rscript

# example use:
# 

library(optparse)
library(data.table)
library(CpGWAS)
library(bsseq)

option_list <- list(
  make_option(c("--df"), type="character", help="Path to data frame"),
  make_option(c("--df_row"), type="integer", help="Data frame row number"),
  make_option(c("--chunk1"), type="integer", help="Start chunk number"),
  make_option(c("--chunk2"), type="integer", help="End chunk number"),
  make_option(c("--output_file"), type="character", help="Name of output file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

df <- fread(opt$df)
df_row <- opt$df_row
chunk1 <- opt$chunk1
chunk2 <- opt$chunk2
output_file <- opt$output_file

methylation_data <- df$methylation_data[df_row]
cov_file <- df$cov_file[df_row]

load(methylation_data)

p <- t(as.matrix(getMeth(BSobj2, type = "smooth", what = "perBase")))
ind <- gsub("Br0", "Br", BSobj2@colData$brnum)
rownames(p) <- ind

covs <- fread(cov_file)
covs$Sex[covs$Sex == "M"] <- 0
covs$Sex[covs$Sex == "F"] <- 1
covs$Dx[covs$Dx == "Control"] <- 0
covs$Dx[covs$Dx == "SCZ"] <- 1

missing_ids <- ind[!ind %in% covs$ID]
if (length(missing_ids) > 0) {
  cat("Missing covariates for IDs:", paste(missing_ids, collapse = ", "), "\n")
  ind <- ind[ind %in% covs$ID]
}

covs <- covs[match(ind, covs$ID), ]
p <- p[match(ind, rownames(p)), ]
id <- colData(BSobj2)$ID[which(colData(BSobj2)$brnum %in% ind)]

covs <- cbind(0, covs)
colnames(covs)[1] <- "intercept"

write.table(ind, "id", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(covs, "covs", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

chr <- unique(as.character(seqnames(rowRanges(BSobj2))))
if (length(chr) > 1) stop("Should be just one chromosome per BSobj")
chr <- gsub("chr", "", chr)

CpG_positions <- start(ranges(granges(BSobj2)))

res <- data.frame()
wind <- c(10000, 100000, 1000000)
gwas <- "/dcs04/lieber/statsgen/shizhong/michael/mwas/gwas/"
gcta <- "/dcs04/lieber/statsgen/shizhong/software/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"

for (i in chunk1:chunk2) {
  cat(i, "\n")
  chr <- gsub("chr", "", chr)
  for (w in 1:length(wind)) {
    p1 <- ifelse(CpG_positions[i] - wind[w] > 0, CpG_positions[i] - wind[w], 0)
    p2 <- CpG_positions[i] + wind[w]
    gwas_prefix <- paste0(gwas, "libd_chr", chr)
    command <- paste("/dcs04/lieber/statsgen/mnagle/mwas/CpGWAS/scripts/plink2 --pfile ", gwas_prefix, "--silent --keep id",
                 "--chr ",chr,
                 "--from-bp",p1,"--to-bp",p2,
                 "--snps-only 'just-acgt' --make-bed --out temp",
                 sep=" ")
    system(command)
    if (!file.exists("temp.bim")) next
    
    pheno <- cbind(0, ind, p[, i])
    write.table(pheno, "pheno", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    command <- paste(gcta, " --bfile temp --make-grm-bin --out temp", sep = "")
    system(command)
    
    command <- paste(gcta, " --reml --grm-bin temp --pheno pheno --mpheno 1 --qcovar covs --out temp", sep = "")
    system(command)
    
    if (!file.exists("temp.hsq")) next
    
    temp <- read.table("temp.hsq", header = TRUE, fill = TRUE)
    temp$site <- paste0("chr", chr, "_", CpG_positions[i])
    temp$wind <- wind[w]
    res <- rbind(res, temp)
    
    system("rm temp*")
  }
}

write.table(res, output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
