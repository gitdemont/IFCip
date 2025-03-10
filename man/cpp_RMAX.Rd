% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_RMAX}
\alias{cpp_RMAX}
\title{Regional Maxima}
\arguments{
\item{img, }{a NumericMatrix.}

\item{n_lev, }{an int determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.}

\item{kernel, }{a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
Default is R_NilValue, resulting in 8-connected pixels neighbors computation.}

\item{msk_, }{a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
Default is R_NilValue, for using all 'img' elements without masking anything.}
}
\value{
a NumericMatrix of regional maxima of 'img'.
}
\description{
Mask connected component of pixels whose values are higher to their external boundaries neighborhood.
}
\details{
see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
RMAX is defined as img > HMAX(img, h_lev = 1).
}
\keyword{internal}
