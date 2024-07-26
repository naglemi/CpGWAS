library("glmnet")
library("e1071")
library("doParallel")
source("model.R")

args = commandArgs(trailingOnly=TRUE)

set.seed(2018)
wind <- c(5000,10000)
# output directory
#outd <- "/dcl02/lieber/shan/shizhong/finemapping/GWAS/tags/scz3/mwas/chr22/1/"
outd <- as.character(args[1])

# load data for mwas
# load("./rda/caudate_mwas_data_chr22.rda")
load(as.character(args[2]))
load("BSsample.rda")

# candidate cg
cg <- as.numeric(rownames(p))

# regress out covariates
load("covs_for_meqtl.rda")
mat <- match(BSsample$brnum,colnames(covs)) 
covs <- t(covs[,mat])
p.residual=matrix(0,dim(p)[1],dim(p)[2])
for(i in 1:dim(p)[1]){
	dat <- as.data.frame(cbind(p[i,],covs))
	colnames(dat) <- c("y",paste0("x",1:33))
	model.res <- lm(reformulate(paste0("x",1:33), "y"),dat)
	p.residual[i,] = resid(model.res) 
}

# built predition models
idx.ea <- BSsample$race == "CAUC"
for(k in 1:length(wind)){
models.ea <- c()
models.all <- c()
for(i in 1:length(cg)){
	cat(i,"\n")
	range1 <- ifelse(cg[i] - wind[k] > 0,cg[i] - wind[k],0)
	range2 <- cg[i] + wind[k]
	idx <- map2$POS > range1 & map2$POS < range2
	# go to next cg if no snps within window
	if(sum(idx) <= 1){
		next
	}
	geno <- snp2[idx,]
	rownames(geno) <- map2$POS[idx]
	trainX <- t(geno)
	trainY <- p.residual[i,]
	fit <- tryCatch(
		elastic.net(trainX,trainY),
		error = function(e) {return ("err")})
	if(fit == "err"){
		next
	}			
	fit$cg <- cg[i]	
	models.all <- rbind(models.all,fit)
	# EA only
	trainX <- trainX[idx.ea,]
	if(sum(apply(trainX,2,var)!=0) <= 1){
		next
	}
	trainY <- trainY[idx.ea]	
	fit <- tryCatch(
		elastic.net(trainX,trainY),
		error = function(e) {return ("err")})
	if(fit == "err"){
		next
	}			
	fit$cg <- cg[i]
	models.ea <- rbind(models.ea,fit)
}
models.ea <- models.ea[models.ea[,1] != "(Intercept)",]
models.all <- models.all[models.all[,1] != "(Intercept)",]

# mwas by models of all samples
cg2 <- unique(models.all$cg)
mwas.all <- matrix(0,nrow=length(cg2),ncol=2)
for(i in 1:length(cg2)){
	pos <- models.all[models.all$cg == cg2[i],1]
	gwas <- snp.gwas2$z[is.element(snp.gwas2$pos_hg38, pos)]
	weight <- models.all[models.all$cg == cg2[i],2]
	geno <- snp.1kg.eur2[match(pos,map.1kg.eur2$POS),]
	mwas.all[i,] <- MWAS(gwas, weight, t(geno))
}
rownames(mwas.all) <- cg2
colnames(mwas.all) <- c("z","p")

# mwas by models of EA samples
cg2 <- unique(models.ea$cg)
mwas.ea <- matrix(0,nrow=length(cg2),ncol=2)
for(i in 1:length(cg2)){
	pos <- models.ea[models.ea$cg == cg2[i],1]
	gwas <- snp.gwas2$z[is.element(snp.gwas2$pos_hg38, pos)]
	weight <- models.ea[models.ea$cg == cg2[i],2]
	geno <- snp.1kg.eur2[match(pos,map.1kg.eur2$POS),]
	mwas.ea[i,] <- MWAS(gwas, weight, t(geno))
}
rownames(mwas.ea) <- cg2
colnames(mwas.ea) <- c("z","p")

# output models and mwas results
outf <- paste0(outd,"/models.all.wind.",wind[k])
write.csv(models.all,outf)
outf <- paste0(outd,"/models.ea.wind.",wind[k])
write.csv(models.ea,outf)
outf <- paste0(outd,"/mwas.all.wind.",wind[k])
write.csv(mwas.all,outf)
outf <- paste0(outd,"/mwas.ea.wind.",wind[k])
write.csv(mwas.ea,outf)
}
	