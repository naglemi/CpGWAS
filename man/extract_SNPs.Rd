% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_SNPs.R
\name{extract_SNPs}
\alias{extract_SNPs}
\title{Extract SNPs within a specified window around a methylation site}
\usage{
extract_SNPs(methInput, meth_site_pos, window_size, verbose, maf, na.action)
}
\arguments{
\item{methInput}{MethylationInput The MethylationInput object containing SNP and methylation data.}

\item{meth_site_pos}{numeric The base position of the methylation site.}

\item{window_size}{numeric Window size for SNPs surrounding methylation site.}

\item{verbose}{logical Whether to print messages.}

\item{maf}{numeric Minor allele frequency threshold for filtering SNPs.}

\item{na.action}{character Whether to remove or impute missing data. Options are "remove" or "impute" (to mean).}
}
\value{
matrix Matrix of SNP data, with row for each sample, column for each SNP (alternative allele count)
}
\description{
This function extracts SNP data from a specified genomic window around a methylation site.
It returns a matrix where each row represents a sample and each column represents an SNP.
}
\examples{
\dontrun{
SNPs <- extract_SNPs(methInput, meth_site_pos = 12345, window_size = 1000)
}
}
