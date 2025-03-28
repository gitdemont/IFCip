% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_rec_erode}
\alias{cpp_rec_erode}
\title{Erosion Reconstruction}
\arguments{
\item{markers, }{a NumericMatrix.}

\item{img, }{a NumericMatrix.}

\item{kernel, }{a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
Default is R_NilValue, resulting in 8-connected pixels neighbors computation.}
}
\value{
a NumericMatrix of reconstructed 'img' from 'markers'.
}
\description{
Performs an erosion reconstruction of an image.
}
\details{
adaptation of 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}
}
\keyword{internal}
