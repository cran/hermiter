% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hermite_estimator.R
\name{spearmans}
\alias{spearmans}
\title{Estimates the Spearman's rank correlation coefficient}
\usage{
spearmans(this, clipped = FALSE)
}
\arguments{
\item{this}{A hermite_estimator_bivar object.}

\item{clipped}{A boolean value. Indicates whether to clip Spearman's rank 
correlation estimates to lie between -1 and 1.}
}
\value{
A numeric value.
}
\description{
This method calculates the Spearman's rank correlation coefficient value. It 
is only applicable to the bivariate Hermite estimator i.e. est_type = 
"bivariate".
}
\details{
The object must be updated with observations prior to the use of this method.
}
\examples{
hermite_est <- hermite_estimator(N = 10, standardize = TRUE,
est_type="bivariate")
hermite_est <- update_batch(hermite_est, matrix(rnorm(30*2), nrow=30, 
ncol=2, byrow = TRUE))
spearmans_est <- spearmans(hermite_est)
}
