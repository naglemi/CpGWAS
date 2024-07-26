library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
args = commandArgs(trailingOnly=TRUE)

chr <- args[1]

#################
# load BS object
##################
#BSobj = loadHDF5SummarizedExperiment('/dcl02/lieber/WGBS/psychENCODE_szControl/Batch3/bs_objs/batch3_combined/', prefix='CpG')
infile <- paste0("/dcs04/lieber/statsgen/shizhong/wgbs/BSobj/caudate/Caudate_chr",chr,"BSobj_GenotypesDosage.rda")
load(infile)
#BSobj = BSobj[seqnames(BSobj) == chr,] #keep chr

# remove blacklist regions
#load("/dcs04/lieber/statsgen/shizhong/wgbs/BSobj/kira_copy_hippo/blacklist.rda")
#bb = findOverlaps(BSobj, blacklist)
#BSobj = BSobj[-queryHits(bb),]

# remove CT snps at CpG sites
f_snp <- paste0("/dcs04/lieber/statsgen/shizhong/AANRI/DEM2/snps_CT/chr",chr)
snp <- fread(f_snp,header=F,data.table=F)
idx <- is.element(start(BSobj), snp[,1])
BSobj <- BSobj[!idx,]

########## 
# AA only
##########
f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
f_demo <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"
ances <- read.table(f_ances,header=T)
demo  <- read.table(f_demo,header=T)

# update brnum, BrXXX -> Br0XXX
for(i in 1:length(ances$id)){
	if (nchar(ances$id[i]) == 5){
		ances$id[i] <- sub("Br","Br0",ances$id[i])
	}
}
for(i in 1:length(demo$BrNum)){
	if (nchar(demo$BrNum[i]) == 5){
		demo$BrNum[i] <- sub("Br","Br0",demo$BrNum[i])
	}
}

id <- intersect(intersect(ances$id[ances$group == "AA"],colData(BSobj)$brnum),demo$BrNum[demo$Age >= 17])
BSobj2 <- BSobj[,is.element(colData(BSobj)$brnum,id)] 

# exlcude low coverage sites
# cov=getCoverage(BSobj2)
# n <- length(colData(BSobj)$brnum)
# keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
# BSobj2 <- BSobj2[keep,]

# Compute means and sds
M=as.matrix(getMeth(BSobj2,type="smooth"))
sds <- rowSds(M)
means <- rowMeans2(M)
save(sds,means,BSobj2,file=paste0("./out/chr",chr,"_AA.rda"))

########## 
# EA only
##########
id <- intersect(intersect(ances$id[ances$group == "CAUC"],colData(BSobj)$brnum),demo$BrNum[demo$Age >= 17])
BSobj2 <- BSobj[,is.element(colData(BSobj)$brnum,id)] 

# exlcude low coverage sites
# cov=getCoverage(BSobj2)
# n <- length(colData(BSobj)$brnum)
# keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
# BSobj2 <- BSobj2[keep,]

# Compute means and sds
M=as.matrix(getMeth(BSobj2,type="smooth"))
sds <- rowSds(M)
means <- rowMeans2(M)
save(sds,means,BSobj2,file=paste0("./out/chr",chr,"_EA.rda"))

########## 
# EA + AA
##########
id <- intersect(intersect(ances$id[ances$group == "CAUC" | ances$group == "AA"],colData(BSobj)$brnum),demo$BrNum[demo$Age >= 17])
BSobj2 <- BSobj[,is.element(colData(BSobj)$brnum,id)] 

# exlcude low coverage sites
# cov=getCoverage(BSobj2)
# n <- length(colData(BSobj)$brnum)
# keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
# BSobj2 <- BSobj2[keep,]

# Compute means and sds
M=as.matrix(getMeth(BSobj2,type="smooth"))
sds <- rowSds(M)
means <- rowMeans2(M)
save(sds,means,BSobj2,file=paste0("./out/chr",chr,"_all.rda"))
