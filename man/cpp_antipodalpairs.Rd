% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_antipodalpairs}
\alias{cpp_antipodalpairs}
\title{Antipodal Pairs of Convex Hull}
\source{
Adaptation from \url{https://escholarship.mcgill.ca/concern/theses/fx719p46g} in Computational geometry with the rotating calipers authored by Pirzadeh, Hormoz under supervision of Toussaint, Godfried T. at McGill University.
}
\arguments{
\item{pts}{a 2-column matrix defining the locations (x and y coordinates, respectively) of points.
It has to be an object of class `IFCip_convexhull`}
}
\value{
an IntegerVector of antipodal pairs of the convex input polygon. Note that this vector of indices is in C so 1st start at 0; add 1 to use it in R.
}
\description{
Computes antipodal pairs of a convex polygon
}
\keyword{internal}
