% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractKdiagonalCsrMatrix.R
\name{extractKdiagonalCsrMatrix}
\alias{extractKdiagonalCsrMatrix}
\title{Extract K diagonal data.}
\usage{
extractKdiagonalCsrMatrix(spsCsrMat)
}
\arguments{
\item{spsCsrMat}{A symmetric matrix need to be extract K diagonal data.}
}
\value{
A list, the elements are the values of contact with the same
distance in the original symmetric matrix.
}
\description{
Extract K diagonal data.
}
\examples{
mtx <- matrix(runif(100, min = -1, max = 1), nrow = 10)
mtx_Symmetric <- (mtx + t(mtx)) / 2
contactBydistance <- extractKdiagonalCsrMatrix(mtx_Symmetric)
contactBydistance
}
