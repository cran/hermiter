% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hermite_estimator.R
\name{spearmans}
\alias{spearmans}
\title{Estimates the Spearman's rank correlation coefficient}
\usage{
spearmans(h_est_obj, clipped = FALSE)
}
\arguments{
\item{h_est_obj}{A hermite_estimator_bivar object.}

\item{clipped}{A boolean value. Indicates whether to clip Spearman's rank 
correlation estimates to lie between -1 and 1.}
}
\value{
A numeric value.
}
\description{
This method utilizes the estimator (8) in the paper Stephanou, Michael and 
Varughese, Melvin. "Sequential estimation of Spearman rank correlation using 
Hermite series estimators." Journal of Multivariate Analysis (2021) 
<doi:10.1016/j.jmva.2021.104783> to calculate the Spearman rank correlation 
coefficient. It is only applicable to the bivariate Hermite estimator 
i.e. est_type = "bivariate".
}
\details{
The object must be updated with observations prior to the use of this method.
}
\examples{
\dontrun{
hermite_est <- hermite_estimator(N = 10, standardize = TRUE,
est_type="bivariate", observations = matrix(rnorm(30*2), nrow=30, 
ncol=2, byrow = TRUE))
spearmans_est <- spearmans(hermite_est)
}
}
