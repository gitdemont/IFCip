% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_convolve2d}
\alias{cpp_convolve2d}
\title{Image Filtering by Convolution}
\arguments{
\item{mat, }{a Matrix.}

\item{kernel, }{a Nullable Matrix.}

\item{method}{used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].}

\item{k, }{constant used for padding, a double. Default is \code{0.0}.}
}
\value{
a NumericMatrix.
}
\description{
This function applies 2D convolution filtering on image.
}
\keyword{internal}
