% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethylationInputClass.R
\name{reinitializeMethylationInput}
\alias{reinitializeMethylationInput}
\title{Reinitialize MethylationInput Object}
\usage{
reinitializeMethylationInput(
  rds_path,
  snp_data_path,
  no_cores = detectCores(),
  start_site = NULL,
  end_site = NULL
)
}
\arguments{
\item{rds_path}{Character string specifying the path to the RDS file containing a previously saved MethylationInput object.}

\item{snp_data_path}{Character string specifying the new path to the SNP data in PLINK2 binary format.}

\item{no_cores}{Integer specifying the number of cores to use for parallel operations (defaults to the number of available cores).}

\item{start_site}{Optional, the start site for subsetting methylation data based on indices in the file.}

\item{end_site}{Optional, the end site for subsetting methylation data based on indices in the file.}
}
\value{
A reinitialized MethylationInput object with updated links to SNP data.
}
\description{
Reloads and updates a MethylationInput object from a saved RDS file, refreshing its SNP data links
to reflect new paths or updates in the SNP data stored in PLINK2 binary format.
}
\examples{
# Reinitialize a MethylationInput object with a new SNP data path
reinitMethInput <- reinitializeMethylationInput("path/to/saved/object.rds", "new/path/to/snp_data", no_cores = 4)

}
