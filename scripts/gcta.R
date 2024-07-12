args = commandArgs(trailingOnly=TRUE)

outdir <- args[1]
chunk1 <- as.numeric(args[2])
chunk2 <- as.numeric(args[3])

wind <- as.numeric(c("1000","2000","5000","10000","20000","50000","100000","200000","500000"))
indir <- "/dcs04/lieber/statsgen/shizhong/AANRI/VMR2/99/caud/aa/"
gwas <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/topmed/merge_H650_1M_2.5M_5M/AA/all/plink/"
gcta <- "/dcs04/lieber/statsgen/shizhong/software/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"

setwd(paste0("out/",outdir))

# load vmrs
load(paste0(indir,"out/chr1_vmr.rda"))
meth2 <- meth
vmrs2 <- vmrs
for(i in 2:22){
	load(paste0(indir,"out/chr",i,"_vmr.rda"))
	meth2 <- rbind(meth2,meth)
	vmrs2 <- rbind(vmrs2,vmrs)
}
p <- t(meth2)
ind <- rownames(p)

# covariates
f_demo <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"
f_pc <- paste0(indir,"/out/sva.csv")
f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
demo  <- read.table(f_demo,header=T)
pc <- read.csv(f_pc)
ances <- read.table(f_ances,header=T)

# align samples
id <- intersect(intersect(demo$BrNum,ind),pc$ind)
demo <- demo[match(id,demo$BrNum),]
pc <- pc[match(id,pc$ind),]
p <- p[match(id,ind),]
ances <- ances[match(id,ances$id),]

# covariates
#covs <- as.data.frame(cbind(Age=demo$Age,Sex=demo$Sex,demo[,11:20],pc[,3:12]))
covs <- as.data.frame(cbind(Age=demo$Age,Sex=demo$Sex,Afr=ances$Afr,pc[,3:12]))
covs$Sex[covs$Sex=="M"] <- 0
covs$Sex[covs$Sex=="F"] <- 1

# overlap samples with genotype data
fam <- paste0(gwas,"AA_chr1.psam")
fam <- read.table(fam,skip = 1, header = F)
id <- intersect(fam[,1],demo$ID) # genotype ID
write.table(id,"id",col.names=F,row.names=F,quote=F,sep="\t")

# align samples
idx <- match(id,demo$ID)
p <- p[idx,]
covs <- covs[idx,]
covs <- cbind(0,id,covs)
write.table(covs,"covs",col.names=F,row.names=F,quote=F,sep="\t")

# loop over vmr between two chunks
res <- c()
for(i in chunk1:chunk2){
	cat(i,"\n")
	chr <- gsub("chr","",vmrs2[i,1])
	# loop over each window size
	for(w in 1:length(wind)){
	# plink subset
	p1 <- ifelse(vmrs2[i,2] - wind[w] > 0, vmrs2[i,2] - wind[w],0)
	p2 <- vmrs2[i,3] + wind[w]
	gwas_prefix <- paste0(gwas,"AA_chr",chr)
	command <- paste("plink2 --pfile ", gwas_prefix, "--silent --keep id", "--chr ",chr,"--from-bp",p1,"--to-bp",p2,"--snps-only 'just-acgt' --make-bed --out temp",sep=" ")
	system(command)	
	if(!file.exists("temp.bim")){
		next;
	}
	# phenotype file 
	pheno <- cbind(0,id,p[,i])
	write.table(pheno,"pheno",col.names=F,row.names=F,quote=F,sep="\t")
	# grm
	command <- paste(gcta, "--bfile temp --make-grm-bin --out temp", sep=" ")
	system(command)
	# h2 estimation
	command <- paste(gcta, "--reml --grm-bin temp --pheno pheno --mpheno 1 --qcovar covs --out temp", sep=" ")
	system(command)
	# collect results
	if(!file.exists("temp.hsq")){
		next;
	}
	temp <- read.table("temp.hsq",header=T, fill=TRUE)
	vmr <- paste0("chr",chr,"_",vmrs2[i,2],"_",vmrs2[i,3])
	temp$vmr <- vmr
	temp$wind <- wind[w]
	res <- rbind(res,temp)
	# remove temp files
	system("rm temp*")
	}
	
}
write.table(res,"res.txt",col.names=T,row.names=F,quote=F,sep="\t") 