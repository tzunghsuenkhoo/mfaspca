% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mFASPCA.R
\name{.PACE_S}
\alias{.PACE_S}
\title{Internal Function - Estimation of univariate functional principal components under CBS}
\usage{
.PACE_S(
  X,
  Y,
  W,
  Y.pred = NULL,
  nbasis = 10,
  pve = 0.99,
  npc = NULL,
  makePD = FALSE,
  cov.weight.type = "none"
)
}
\arguments{
\item{X}{A vector of xValues}

\item{Y}{A matrix of observed functions (by row).}

\item{W}{A weight matrix}

\item{Y.pred}{A matrix of functions (by row) to be approximated using the
functional principal components.}

\item{nbasis}{An integer, giving the number of B-spline basis to use.
Defaults to \code{10}.}

\item{pve}{A value between 0 and 1, giving the percentage of variance
explained in the data by the functional principal components. This value is
used to choose the number of principal components. Defaults to \code{0.99}}

\item{npc}{The number of principal components to be estimated. Defaults to
\code{NULL}. If given, this overrides \code{pve}.}

\item{makePD}{Logical, should positive definiteness be enforced for the
covariance estimate? Defaults to \code{FALSE}.}

\item{cov.weight.type}{The type of weighting used for the smooth covariance
estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively,
\code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the pointwise estimates of the covariance function
by the number of observation points.}
}
\value{
\item{fit}{The approximation of \code{Y.pred} (if \code{NULL}, the
approximation of \code{Y}) based on the functional principal components.}
\item{scores}{A matrix containing the estimated scores (observations by
row).} \item{mu}{The estimated mean function.} \item{efunctions}{A matrix
containing the estimated eigenfunctions (by row).} \item{evalues}{The
estimated eigenvalues.} \item{npc}{The number of principal comopnents that
were calculated.} \item{sigma2}{The estimated variance of the measurement
error.}  \item{estVar}{The estimated smooth variance function of the data.}
}
\description{
This function gives the estimation of the FPCA under CBS,
}
\references{
Happ, C. (2021). Multivariate functional principal component
analysis for data observed on different dimensional domains, R package version
Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
Multilevel functional principal component analysis. Annals of Applied
Statistics, 3, 458--488. Yao, F., Mueller, H.-G., and Wang, J.-L. (2005).
Functional data analysis for sparse longitudinal data. Journal of the
American Statistical Association, 100, 577--590.
}
