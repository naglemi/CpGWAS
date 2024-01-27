library(data.table)
library(testthat)
library(bsseq)

pgen_path = system.file("extdata", "chr1_sample_subset.pgen", package = "CpGWAS")
pvar_path = system.file("extdata", "chr1_sample_subset.pvar", package = "CpGWAS")
psam_path = system.file("extdata", "chr1_sample_subset.psam", package = "CpGWAS")

scaffold_name <- "unit_test_scaffold"

# read methylations from RDS

data(chr1_methylation_sample_subset, package = "CpGWAS")
# Processing methylation data from BSseq object
methylations <- t(as.matrix(getMeth(BSobj_sample, type = "smooth", what = "perBase")))
colnames(methylations) <- paste0("pos_",
                                 GenomicRanges::start(
                                   GenomicRanges::ranges(
                                     SummarizedExperiment::rowRanges(BSobj_sample))))

pvar_pointer <- pgenlibr::NewPvar(pvar_path)
pvar_dt <- fread(pvar_path)
pgen <- pgenlibr::NewPgen(pgen_path, pvar = pvar_pointer)
psam <- fread(psam_path)
psam_in_wgbs <- psam[which(psam$`#IID` %in% rownames(methylations))]
genotype_IDs <- psam_in_wgbs$`#IID`

genotype_IDs <- intersect(rownames(methylations), genotype_IDs)
genotype_IDs <- genotype_IDs[order(genotype_IDs)]

cov <- processCovariates(dataFrame = colData(BSobj_sample),
                         colsToExclude = c("ID.", "DNum", "brnum",
                                           "BrNum", "brnumerical"),
                         genotype_IDs = genotype_IDs)

# Begin tests for data loaded without use of `MethylationInput`
test_that("Without using `MethylationInput`, methylation data has 112 samples", {
  expect_equal(length(rownames(methylations)), 112)
})

methylations <- methylations[which(rownames(methylations) %in% genotype_IDs), ]

test_that("Without using `MethylationInput`, after filtering, methylation data has 111 samples", {
  expect_equal(length(rownames(methylations)), 111)
})

test_that("Without using `MethylationInput`, SNP data has 2189 samples as expected", {
  expect_equal(length(psam$`#IID`), 2189)
})

test_that("Without using `MethylationInput`, intersection of methylation, SNP data has 111 samples as expected", {
  genotype_IDs <- intersect(rownames(methylations), genotype_IDs)
  expect_equal(length(genotype_IDs), 111)
})

# Three tests for covariate loading and processing

# 1. Test for Correct Column Names and Types
test_that("cov has correct column names and types", {
  expect_true(is.matrix(cov))
  expected_column_names <- c("(Intercept)", "sexM", "primarydxSchizo", "agedeath", "pmi")
  expect_equal(colnames(cov), expected_column_names)
  expect_true(is.numeric(cov[, 2]))
})

# 2. Test for Correct Number of Rows
test_that("cov has the correct number of rows (samples)", {
  expected_row_count <- 111
  expect_equal(nrow(cov), expected_row_count)
})

# 3. Test for Specific Value Consistency
test_that("Selected specific values in cov are as expected", {
  # Replace the indices and expected values with the correct ones for your data
  expect_equal(cov["Br1003", "sexM"], 0)
  expect_equal(cov["Br1004", "primarydxSchizo"], 0)
  # Add more checks as needed for other specific values
})

methylations <- regress_out_cov_parallel(methylations, cov)

# 1. Test for Output Dimensions
test_that("After regressing out covariates, methylations has correct dimensions", {
  expected_dimensions <- c(111, 11) # Replace with expected dimensions
  expect_equal(dim(methylations), expected_dimensions)
})

# 2. Test for Data Type
test_that("After regressing out covariates, methylations contains numeric values", {
  expect_true(is.numeric(methylations))
})

# 3. Test for Specific Value Range
# Replace -1 and 1 with the range you expect based on your data
test_that("After regressing out covariates, values in methylations are within expected range", {
  expect_true(all(methylations >= -0.15 & methylations <= 0.11))
})

# 4. Test for Specific Known Values
test_that("After regressing out covariates, specific values in methylations are as expected", {
  tolerance <- 1e-6
  expect_true(abs(methylations["Br1003", 1] - 0.0427085754) < tolerance)
  expect_true(abs(methylations["Br1004", 2] - (-0.0200391506)) < tolerance)
})

