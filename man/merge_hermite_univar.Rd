% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hermite_estimator_univar.R
\name{merge_hermite_univar}
\alias{merge_hermite_univar}
\title{Merges a list of Hermite estimators}
\usage{
merge_hermite_univar(hermite_estimators)
}
\arguments{
\item{hermite_estimators}{A list of hermite_estimator_univar objects.}
}
\value{
An object of class hermite_estimator_univar.
}
\description{
This method allows a list of Hermite based estimators of class
hermite_estimator_univar to be consistently merged.
}
\details{
Note that the N and standardize arguments must be the same for all estimators
in order to merge them. In addition, note that exponentially weighted
estimators cannot be merged. If the Hermite estimators are not
standardized, the merged estimator will be exactly equivalent to
constructing a single estimator on the data set formed by combining the
data sets used to update the respective hermite_estimator_univar inputs.
If the input Hermite estimators are standardized however, then the
equivalence will be approximate but still accurate in most cases.
}
