#' MWAS function
#'
#' Performs methylation-wide association study analysis.
#'
#' @param z Z-scores (log(OR)/SE) for effect of SNPs on external phenotype.
#' @param w Weights for effect of SNPs on methylation.
#' @param G SNP genotype matrix.
#' @return Named vector containing z-score, p-value, and the number of weights.
#' @export
mwas <- function(z, w, G){   
  #recover()
  if(length(w) > 1){
    #recover()
    # z-scores for effect of SNPs on external phenotype
    #. are weighted according to weights for effect of SNPs on methylation
    z <- z %*% w
    # compute correlation matrix of SNP matrix, which captures LD structure
    z.cor <- cor(G)
    # add small value to diagonal to avoid singular matrix
    #  which may otherwise happen if two SNPs in perfect LD
    z.cor <- z.cor + diag(dim(z.cor)[1])*0.1 
    # variance of correlated variables is weighted sum 
    # multiplying w by corr matrix once gives a vector representing
    #. the variance of each individual SNP and the extent to which they are
    #. influenced by other SNPs. Multiplying again by w sums up pairwise contributions
    #. and reflects total variance of weighted sum.
    #. the first w is automatically transposed by R so we don't have to write t(w)
    se <- sqrt(w %*%  z.cor %*%  w)
    z <- z/se
    p <- pnorm(abs(z), lower.tail=F)*2
    return(c(z=z, p=p, n=length(w)))
  } else {
    p <- pnorm(abs(z), lower.tail=F)*2
    return(c(z=z, p=p, n=1))
  }
}

#' MWASmodel class
#' @export
setClass(
  "MWASmodel",
  representation(
    methylationBase = "MethylationBase",
    #summary_stats = "data.table",
    mwas_out = "numeric"
  )
)

#' MWASmodel constructor
#' @param methylationBase MethylationBase object
#' @param summary_stats Data table of summary statistics
#' @param mwas_out Numeric vector of MWAS output
#' @return MWASmodel object
#' @export
MWASmodel <- function(methylationBase,
                      #summary_stats,
                      mwas_out) {
  new("MWASmodel",
      methylationBase = methylationBase,
      #summary_stats = summary_stats,
      mwas_out = mwas_out)
}

#' Process a single MWAS model
#'
#' @param methylationBase MethylationBase object
#' @param my_SNPs SNP data
#' @param summary_stats Data table of summary statistics
#' @return MWASmodel object
#' @export
#' @importFrom stringr str_split_fixed
#' @importFrom data.table as.data.table setnames setkey
#' @importFrom data.table `%chin%`
#' @importFrom pgenlibr ReadList
process_model <- function(methylationBase, my_SNPs, summary_stats) {
  
  #recover()
  SNP_split <- stringr::str_split_fixed(names(methylationBase@snpWeights), ":", 4)

  # Convert SNP_split to data.table and set integer types
  SNP_split_dt <- data.table(chr = as.integer(gsub("chr", "", SNP_split[,1])),
                             post = as.integer(SNP_split[, 2]),
                             ref = SNP_split[, 3],
                             alt = SNP_split[, 4],
                             key = c("chr", "post"))
  
  #setkey(SNP_split_dt, chr, post)
  
  # Use a join with the keys
  #recover()
  relevant_SNP_indices <- my_SNPs$pvar_dt[SNP_split_dt, on = .(`#CHROM` = chr, POS = post), which = TRUE, nomatch = 0]
  # We only want summary stats for the specific SNPs contributing to this
  #. methylation site in our model
  relevant_ids <- my_SNPs$pvar_dt$ID[relevant_SNP_indices]
  
  # Subset summary_stats in constant time using a keyed join
  #recover()
  summary_stats_sub <- summary_stats[relevant_ids, nomatch = 0]
  
  #summary_stats_sub <- summary_stats[SNP %chin% relevant_ids]
  # # z is a vector of the SNP weights from GWAS summary statistics
  # z <- summary_stats_sub$logOR
  # 
  # # w is a vector of the SNP weights from the CpGWAS model
  # w <- methylationBase@snpWeights[relevant_SNP_indices]
  
  # Ensuring the order matches and handling unmatched positions
  
  # Check my_SNPs
  #recover()
  
  if(!identical(summary_stats_sub$BP, SNP_split_dt$post)){
    # Order summary_stats_sub by BP
    summary_stats_sub <- summary_stats_sub[order(summary_stats_sub$BP), ]
    if(!identical(summary_stats_sub$BP, SNP_split_dt$post)){
      # Identify positions in SNP_split not found in summary_stats_sub$BP
      unmatched_positions <- !SNP_split_dt$post %in% summary_stats_sub$BP
      if (any(unmatched_positions)) {
        #recover()
        # Remove rows from SNP_split where positions do not match any in summary_stats_sub$BP
        SNP_split_dt <- SNP_split_dt[!unmatched_positions, ]
        methylationBase@snpWeights <- methylationBase@snpWeights[!unmatched_positions]
        # Assuming you would need to recompute the relevant SNP indices and stats
        relevant_SNP_indices <- my_SNPs$pvar_dt[SNP_split_dt, on = .(`#CHROM` = chr, POS = post), which = TRUE, nomatch = 0]
        relevant_ids <- my_SNPs$pvar_dt$ID[relevant_SNP_indices]
        summary_stats_sub <- summary_stats[relevant_ids, nomatch = 0]
      }
      # Check again after removing unmatched positions
      if(!identical(summary_stats_sub$BP, SNP_split_dt$post)) {
        stop("SNP order does not match even after removing unmatched positions. This should not happen. Code is broken.")
      }
    }
  }
  
  
  # # need to make sure direction is right
  # if(!identical(SNP_split[, 4], summary_stats_sub$A2) |
  #    !identical(SNP_split[, 3], summary_stats_sub$A1)){
  #   recover()
  #   not_matching <- which(SNP_split[, 4] != summary_stats_sub$A2)
  #   # Flip our data to match the summary stats for these
  #   summary_stats_ref_flipped <- SNP_split[, 3][not_matching]
  #   summary_stats_alt_flipped <- SNP_split[, 4][not_matching]
  #   SNP_split[, 3][not_matching] <- summary_stats_alt_flipped
  #   SNP_split[, 4][not_matching] <- summary_stats_ref_flipped
  #   methylationBase@snpWeights[not_matching] <-
  #     methylationBase@snpWeights[not_matching] * -1
  # }
  
  #recover()
  
  #print(SNP_split_dt)
  #print(summary_stats_sub)
  #cat("\n")
  
  # need to make sure direction is right but use SNP_split_dt now
  #if(!identical(SNP_split_dt$alt, summary_stats_sub$A2) |
  #   !identical(SNP_split_dt$ref, summary_stats_sub$A1)){
    #recover()
  #  not_matching <- which(SNP_split_dt$alt != summary_stats_sub$A2)

  these_SNPs_pvar_dt <- my_SNPs$pvar_dt[relevant_SNP_indices]
  
  
  if(!all(these_SNPs_pvar_dt$POS == summary_stats_sub$BP)) {
    stop("ERROR! Mismatch between reference dataset SNPs and summary statistics.")
  }
  
  #recover()
  # need to make sure direction is right but use SNP_split_dt now
  if(!identical(these_SNPs_pvar_dt$ALT, summary_stats_sub$A2) |
     !identical(these_SNPs_pvar_dt$REF, summary_stats_sub$A1)){
    #recover()
    #print("We're flipping alleles")
    #print(SNP_split_dt)
    #print(summary_stats_sub)
    cat("\n")
    not_matching <- which(summary_stats_sub$A2 != these_SNPs_pvar_dt$ALT)
    
    if (length(not_matching) > 0) {
      cat("\nWe're flipping alleles due to non-matching entries:\n")
      
      # Print non-matching rows from summary_stats_sub
      print("Non-matching entries in summary_stats_sub:")
      print(summary_stats_sub[not_matching,])
      
      # Print non-matching rows from these_SNPs_pvar_dt
      print("Non-matching entries in these_SNPs_pvar_dt:")
      print(these_SNPs_pvar_dt[not_matching,])
    } else {
      cat("All alleles match.\n")
    }
    
    # Flip our data to match the summary stats for these
    summary_stats_A1_flipped <- summary_stats_sub$A1[not_matching]
    summary_stats_A2_flipped <- summary_stats_sub$A2[not_matching]
    summary_stats_sub$A1[not_matching] <- summary_stats_A2_flipped
    summary_stats_sub$A2[not_matching] <- summary_stats_A1_flipped
    
    # actually instead flip in summary_stats_sub
    #methylationBase@snpWeights[not_matching] <-
    #  methylationBase@snpWeights[not_matching] * -1
    summary_stats_sub$Z[not_matching] <- summary_stats_sub$Z[not_matching] * -1
  }
  
  
  
  
  # Subset the genotype data
  G <- pgenlibr::ReadList(my_SNPs$pgen,
                          variant_subset = relevant_SNP_indices)
  
  mwas_out <- mwas(z = summary_stats_sub$Z,
                   w = methylationBase@snpWeights,
                   G = G)
  
  # Too bulky
  #MWASmodel(methylationBase,
  #          #summary_stats_sub,
  #          mwas_out)
  
  return(mwas_out)
}


#' MWASresults class
#' @export
setClass(
  "MWASresults",
  representation(
    MWASmodels = "list",
    pvar_path = "character",
    pgen_path = "character",
    psam_path = "character",
    summary_stats_path = "character",
    rds_path = "character"
  )
)

#' MWASresults constructor
#' @param MWASmodels List of MWASmodel objects
#' @param pvar_path Path to pvar file
#' @param pgen_path Path to pgen file
#' @param psam_path Path to psam file
#' @param summary_stats_path Path to summary statistics file
#' @param rds_path Path to RDS file
#' @return MWASresults object
#' @export
MWASresults <- function(MWASmodels, pvar_path, pgen_path, psam_path, summary_stats_path, rds_path) {
  new("MWASresults",
      MWASmodels = MWASmodels,
      pvar_path = pvar_path,
      pgen_path = pgen_path,
      psam_path = psam_path,
      summary_stats_path = summary_stats_path,
      rds_path = rds_path)
}

#' Clean and standardize column names
#'
#' @param summary_stats Data table of summary statistics
#' @return Data table with standardized column names
#' @export
#' @importFrom stringr str_split
#' @importFrom data.table setkey
clean_and_standardize_colnames <- function(summary_stats) {
  # Check if the header is tab-delimited while the rest is space-delimited
  if (grepl("\t", colnames(summary_stats)[1])) {
    real_colnames <- str_split(colnames(summary_stats)[1], "\t")[[1]]
    colnames(summary_stats) <- real_colnames
  }
  
  # Standardize column names
  colnames(summary_stats) <- gsub("chr", "CHR", colnames(summary_stats))
  colnames(summary_stats) <- gsub("#CHROM", "CHR", colnames(summary_stats))
  colnames(summary_stats) <- gsub("pos", "BP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("POS", "BP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("MarkerName", "SNP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("ID", "SNP", colnames(summary_stats))
  colnames(summary_stats) <- gsub("LogOR", "logOR", colnames(summary_stats))
  colnames(summary_stats) <- gsub("StdErrLogOR", "SE", colnames(summary_stats))
  colnames(summary_stats) <- gsub("StdErrlogOR", "SE", colnames(summary_stats))
  colnames(summary_stats) <- gsub("BETA", "logOR", colnames(summary_stats))
  
  # If there's no logOR columns, create one, which will be log of OR column
  # but we only do this if there's already an OR column
  if(!"logOR" %in% colnames(summary_stats)) {
    if("OR" %in% colnames(summary_stats)) {
      summary_stats[, logOR := log(OR)]
    }
  }
  
  #colnames(summary_stats) <- gsub("logOR", "BETA", colnames(summary_stats))
  if(!"Z" %in% colnames(summary_stats)) {
    print("Setting Z")
    summary_stats[, Z := logOR/SE]
  } else {
    print("Z set already")
  }
  # Convert summary_stats to a keyed data.table for fast lookups

  setkey(summary_stats, SNP)
  
  return(summary_stats)
}

#' Process MWAS models
#'
#' @param my_rds An object containing models
#' @param my_SNPs SNP data
#' @param summary_stats Data table of summary statistics
#' @param paths List of paths to data files
#' @param summary_stats_path Path to summary statistics file
#' @return MWASresults object
#' @export
#' @importFrom progress progress_bar
#' @importFrom data.table setkey
process_MWAS_models <- function(my_rds, my_SNPs, paths, summary_stats_path, rds_path, summary_stats = NULL) {
  # pb <- progress_bar$new(
  #   format = "[:bar] :percent eta: :eta",
  #   total = length(my_rds@models), clear = FALSE, width = 60
  # )
  
  MWASmodels <- vector("list", length(my_rds@models))
  
  if(is.null(summary_stats)) {
    summary_stats <- suppressWarnings(fread(summary_stats_path))
    summary_stats <- clean_and_standardize_colnames(summary_stats)
  }

  for (i in seq_along(my_rds@models)) {
    this_MethylationBase <- my_rds@models[[i]]
    MWASmodels[[i]] <- process_model(this_MethylationBase, my_SNPs, summary_stats)
    #pb$tick()
  }
  
  # Ensure the lengths of my_rds@models and MWASmodels are the same
  stopifnot(length(my_rds@models) == length(MWASmodels))
  
  results <- MWASresults(MWASmodels, paths$pvar_path, paths$pgen_path, paths$psam_path, summary_stats_path, rds_path)
  
  return(results)
}