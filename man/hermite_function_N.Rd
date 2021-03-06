% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{hermite_function_N}
\alias{hermite_function_N}
\title{Convenience function to output orthonormal Hermite functions
 
 
The method calculates the orthonormal Hermite functions, \eqn{h_k(x)} 
from \eqn{k=0,\dots,N} for the vector of values, x.}
\usage{
hermite_function_N(N, x, normalization_hermite = NULL)
}
\arguments{
\item{N}{An integer number.}

\item{x}{A numeric vector.}

\item{normalization_hermite}{A numeric vector. A vector of normalization 
values necessary in the calculation of the Hermite functions.}
}
\value{
A numeric matrix with N+1 rows and length(x) columns.
}
\description{
Convenience function to output orthonormal Hermite functions
 
 
The method calculates the orthonormal Hermite functions, \eqn{h_k(x)} 
from \eqn{k=0,\dots,N} for the vector of values, x.
}
