\encoding{UTF-8}
\name{recSplitConnMat}
\alias{recSplitConnMat}
\title{Recursively subdivides a set of subpoplations}
\usage{
recSplitConnMat(subpops.lst, conn.mat, beta, ...)
}
\arguments{
  \item{subpops.lst}{A list whose elements are vectors of
  indices for each subpopulation.  See
  \code{\link{subpopsVectorToList}}.}

  \item{conn.mat}{A square connectivity matrix.  This
  matrix has typically been normalized and made symmetric
  prior to using this function.}

  \item{beta}{Controls degree of splitting of connectivity
  matrix, with larger values generating more
  subpopulations.}

  \item{\dots}{further arguments to be passed to
  \code{\link{splitConnMat}}}
}
\description{
This funtion recursively splits each subpopulation of a
list of subpopulations until none of the subpopulations can
be split further to improve the minimization.
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\references{
Jacobi, M. N., André, C., Döös, K., and Jonsson, P. R.
2012. Identification of subpopulations from connectivity
matrices. Ecography, 35: 1004-1016.
}
\seealso{
See also \code{\link{optimalSplitConnMat}},
\code{\link{splitConnMat}},
\code{\link{subpopsVectorToList}}
}
