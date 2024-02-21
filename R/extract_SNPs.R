#' Extract SNPs within a specified window around a methylation site
#'
#' This function extracts SNP data from a specified genomic window around a methylation site.
#' It returns a matrix where each row represents a sample and each column represents an SNP.
#'
#' @param methInput MethylationInput The MethylationInput object containing SNP and methylation data.
#' @param meth_site_pos numeric The base position of the methylation site.
#' @param window_size numeric Window size for SNPs surrounding methylation site.
#' @return matrix Matrix of SNP data, with row for each sample, column for each SNP (alternative allele count)
#'
#' @examples
#' \dontrun{
#' SNPs <- extract_SNPs(methInput, meth_site_pos = 12345, window_size = 1000)
#' }
#' @export
#'
#'

extract_SNPs <- function(methInput, meth_site_pos, window_size) {
  if (!inherits(methInput, "MethylationInput")) {
    stop("methInput must be a MethylationInput object.")
  }
  if (!is.numeric(meth_site_pos) || length(meth_site_pos) != 1) {
    stop("meth_site_pos must be a single numeric value.")
  }
  if (!is.numeric(window_size) || length(window_size) != 1) {
    stop("window_size must be a single numeric value.")
  }

  lower_bound <- meth_site_pos - window_size
  upper_bound <- meth_site_pos + window_size

  snp_indices <- which(methInput@pvar_dt$POS >= lower_bound & methInput@pvar_dt$POS <= upper_bound)
  snp_IDs <- methInput@pvar_dt$ID[snp_indices]

  if (length(snp_indices) <= 1) {
    return(NULL)  # No SNPs found in the specified window
  }
  SNPs <- pgenlibr::ReadList(methInput@pgen, variant_subset = snp_indices)

  rownames(SNPs) <- methInput@psam$`#IID`
  colnames(SNPs) <- snp_IDs
  SNPs <- reorder_and_filter_geno(geno = SNPs, genotype_IDs = methInput@genotype_IDs)

  # Edge case resulting from when multiple alternative alleles at same position,
  #  which can be filtered out and replaced with NA
  if(length(which(is.na(SNPs), arr.ind = TRUE)) > 0) {
    print("removing NA columns")
    SNPs <- SNPs[, colSums(!is.na(SNPs)) > 0]
  }
  
  # If minor allele count for a given SNP is not above 1, exclude
  #  to avoid edge case that causes error:
  #  "from glmnet C++ code (error code 7777); All used predictors have zero variance"
  if(any(colSums(SNPs, na.rm = TRUE) <= 1)) {
    print("removing SNPs with minor allele count <= 1")
    SNPs <- SNPs[, colSums(SNPs, na.rm = TRUE) > 1, drop = FALSE]
  }

  return(SNPs)
}

are_columns_identical <- function(X) {
  # Use rowSums to check if the sum of absolute differences from the first column is zero
  identical_cols <- all(rowSums(abs(X - X[, 1]) == 0) == ncol(X))
  return(identical_cols)
}
