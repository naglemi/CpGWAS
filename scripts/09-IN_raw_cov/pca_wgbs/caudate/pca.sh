gwas=/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/topmed/merge_H650_1M_2.5M_5M/EA_AA/all/plink/all_brnum_all

# keep caudate samples
plink2 --pfile $gwas --keep caudate_brnum --make-pgen --out caudate_gwas

# LD pruning of SNPs
plink2 --pfile caudate_gwas --maf 0.05 --indep-pairwise 100 5 0.1 --out snps

# IBD pruning of samples of close relationships
plink2 --pfile caudate_gwas --make-king-table --king-table-filter 0.05 --out kin --extract snps.prune.in
#awk 'NR>1 && $6>0.05 {print $1}' kin.kin0 | sort | uniq > related_samples.txt

# PCA
#plink2 --pfile caudate_gwas --extract snps.prune.in --remove related_samples.txt --pca 20 --out pca
plink2 --pfile caudate_gwas --extract snps.prune.in --pca 20 --out pca