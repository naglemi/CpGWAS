#' Extract SNPs within a specified window around a methylation site
#'
#' This function extracts SNP data from a specified genomic window around a methylation site.
#' It returns a matrix where each row represents a sample and each column represents an SNP.
#'
#' @param methInput MethylationInput The MethylationInput object containing SNP and methylation data.
#' @param meth_site_pos numeric The base position of the methylation site.
#' @param window_size numeric Window size for SNPs surrounding methylation site.
#' @param verbose logical Whether to print messages.
#' @param maf numeric Minor allele frequency threshold for filtering SNPs.
#'
#' @return matrix Matrix of SNP data, with row for each sample, column for each SNP (alternative allele count)
#'
#' @examples
#' \dontrun{
#' SNPs <- extract_SNPs(methInput, meth_site_pos = 12345, window_size = 1000)
#' }
#' @export
#'
#'

extract_SNPs <- function(methInput, meth_site_pos, window_size, verbose, maf) {
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
    return(NULL)  # No more than 1 SNP found in the specified window
  }
  SNPs <- pgenlibr::ReadList(methInput@pgen, variant_subset = snp_indices)

  rownames(SNPs) <- methInput@psam$`#IID`
  colnames(SNPs) <- snp_IDs
  SNPs <- reorder_and_filter_geno(geno = SNPs, genotype_IDs = methInput@genotype_IDs)
  
  # Filter by maf
  if (maf > 0){
    mafs <- colMeans(SNPs, na.rm = TRUE) / 2
    mafs_below_threshold <- mafs < maf
    if (any(mafs_below_threshold)) {
      if (verbose) {
        message(paste0("removing ", sum(mafs_below_threshold), " SNP(s) with MAF < ",
                       maf, " for position ", meth_site_pos,
                      " with window size ", window_size, ".\n\n"))
      }
      SNPs <- SNPs[, !mafs_below_threshold, drop = FALSE] 
    }
  }

  # If only 1 or fewer SNP remains after maf filtering, skip window (return NULL)
  if (ncol(SNPs) <= 1) {
    if (verbose) {
      message(paste0("For site at position ", meth_site_pos,
                     ", only one SNP remains after MAF filtering in the window of size ",
                     window_size, ".\n\n"))
    }
    return(NULL)
  }
  
  # Edge case resulting from when multiple alternative alleles at same position,
  #  which can be filtered out and replaced with NA
  if (length(which(is.na(SNPs), arr.ind = TRUE)) > 0) {
    if (verbose) {
      message(paste0("removing NA columns for site at position ", meth_site_pos,
                     " with window size ", window_size, ".\n\n"))
    }
    SNPs <- SNPs[, colSums(!is.na(SNPs)) > 0]
  }
  
  # If minor allele count for a given SNP is not above 1, exclude
  #  to avoid edge case that causes error:
  #  "from glmnet C++ code (error code 7777); All used predictors have zero variance"
  if (any(colSums(SNPs, na.rm = TRUE) <= 1)) {
    if (verbose) {
      message(paste0("removing SNP(s) with minor allele count <= 1 for position ", meth_site_pos,
                     " with window size ", window_size, ".\n\n"))
    }
    SNPs <- SNPs[, colSums(SNPs, na.rm = TRUE) > 1, drop = FALSE]
  }
  
  # Edge case where all SNPs are identical
  if (are_columns_identical(SNPs)) {
    if (verbose) {
      message(paste0("For site at position ", meth_site_pos,
                     ", all SNPs are identical in the window of size ",
                     window_size, "\n\n"))
    }
    return(NULL)
  }
  
  # Edge case where all SNPs identical except for in one genotype
  n_rows_identical_columns <- sum(apply(SNPs, 1, function(x) all(x == x[1])))
  n_rows_unique_columns <- nrow(SNPs) - n_rows_identical_columns
  if (n_rows_unique_columns == 1) {
    if (verbose) {
      message(paste0("For site at position ", meth_site_pos,
                     ", case where SNPs differ by only one sample in the window of size ",
                     window_size, "\n\n"))
    }
    return(NULL)
  }
  return(SNPs)
}

are_columns_identical <- function(X) {
  # Use rowSums to check if the sum of absolute differences from the first column is zero
  identical_cols <- all(rowSums(abs(X - X[, 1]) == 0) == ncol(X))
  return(identical_cols)
}

