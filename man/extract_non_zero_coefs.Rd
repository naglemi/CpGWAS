% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tuning_cross_validation.R
\name{extract_non_zero_coefs}
\alias{extract_non_zero_coefs}
\title{Extract Non-Zero Coefficients from an Elastic Net Model}
\usage{
extract_non_zero_coefs(fitted_model)
}
\arguments{
\item{model_fit}{The fitted elastic net model object from 'glmnet'.}
}
\value{
A data frame where each row represents a SNP with a non-zero coefficient.
Columns include SNP names and their corresponding coefficients.
}
\description{
Post-model fitting in a methylome-wide association study, this function extracts the non-zero
coefficients from an elastic net model. These coefficients represent the SNPs that have a
significant impact on the response variable (methylation levels), indicating relevant genetic
markers.
}
\examples{
\dontrun{
features <- extract_non_zero_coefs(fitted_model)
}
}
