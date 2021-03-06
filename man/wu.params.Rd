% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp.R
\name{wu.params}
\alias{wu.params}
\title{Non-central chi-square approximation to mixture of 1-DF chi-square distributions}
\usage{
wu.params(lambda, delta, N = 12)
}
\arguments{
\item{lambda}{mixing coefficients for chi-square mixtures}

\item{delta}{vector of non-centrality parameters for chi-square mixtures}

\item{N}{match (N-1,N)-th moments}
}
\value{
\describe{
  \item{l}{ estimated DF }
  \item{d}{ estimated non-centrality parameter }
  \item{muQ,sigmaQ}{ first two moments based on the chi-square mixture distribution}
  \item{muX,sigmaX}{ first two moments based on the estimated non-central chi-square distribution}
}
}
\description{
Given the mixture of 1-DF chi-square distribution, \eqn{\sum_i\lambda_i\chi_1^2(\delta_i)},
the non-central chi-square distribution approximation is based on matching higher moments. See the ref.
}
\references{
Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.

Wu,B. and Pankow,J.S. (2017) On computing the tail probability of non-negative definite quadratic forms in central normal variables. tech rep.
}
