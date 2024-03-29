% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_thinning_bst}
\alias{cpp_thinning_bst}
\title{Implementation of Ben Boudaoud-Sider-Tari Thinning}
\arguments{
\item{mat}{a LogicalMatrix, containing mask.}
}
\value{
a LogicalMatrix with the mask thinned.
}
\description{
This function is designed to identify mask skeleton.
}
\details{
adaptation of 'A new thinning algorithm for binary images' from L. Ben Boudaoud, A. Sider, A. Tari.
3rd international conference on control, engineering & information technology, May 2015. \doi{10.1109/CEIT.2015.7233099}.
}
\keyword{internal}
