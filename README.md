# CpGWAS
Tools to support Methylome-wide Association Studies (MWAS), not to be confused with Metabolome-Wide Association Studies (the other MWAS)

## Installation

Bioconductor dependencies need to be installed first.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("GenomicRanges")
BiocManager::install("bsseq")
```

Next, install CpGWAS from the GitHub repository.
```r
If `devtools` isn't already installed:
install.packages("devtools")

devtools::install_github("https://github.com/naglemi/CpGWAS.git")
```
