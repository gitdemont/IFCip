% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractBasic.R
\name{ExtractBasic}
\alias{ExtractBasic}
\title{Basic Features Extraction}
\usage{
ExtractBasic(
  ...,
  objects,
  offsets,
  removal = "masked",
  batch = 20L,
  display_progress = TRUE,
  parallel = FALSE
)
}
\arguments{
\item{...}{arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
If 'param' is provided 'export'(="matrix"), 'mode'(="raw"), 'size'(="c(0,0)"), 'force_width'(="FALSE") and 'removal' will be overwritten.\cr
If 'offsets' are not provided extra arguments can also be passed with ... to \code{\link{getOffsets}}.\cr
/!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
Remaining arguments with the exception of 'strategy', 'envir' and '...' will be passed to \link[future]{plan}.}

\item{objects}{integers, indices of objects to use.
This argument is not mandatory, if missing, the default, all objects will be used.}

\item{offsets}{object of class `IFC_offset`. 
This argument is not mandatory but it may allow to save time for repeated image export on same file.}

\item{removal}{whether to compute features on "masked" object for each individual channels or on the globally detected object "MC".
Allowed are "masked" or "MC". Default is "masked". Please note that it will overwrite 'param' value if provided.}

\item{batch}{positive integer, number of objects to process at the same time. Default is 20L.}

\item{display_progress}{whether to display a progress bar. Default is TRUE.\cr
When NULL, execution will not be wrapped inside \link[progressr]{with_progress} nor \link[progressr]{withProgressShiny}. This allow user to call the function with \link[progressr]{with_progress} nor \link[progressr]{withProgressShiny} or to use global handler see \link[progressr]{handlers}.\cr
When FALSE, execution will be performed inside \link[progressr]{without_progress}.\cr
When TRUE, execution will be wrapped inside \link[progressr]{with_progress} or \link[progressr]{withProgressShiny}
and \link[progressr]{handlers} will be automatically selected (the last available will be chosen between either):\cr
- \link[progressr]{handler_txtprogressbar},\cr
- a customized version of \link[progressr]{handler_winprogressbar}, (if on windows OS),\cr
- \link[progressr]{handler_shiny} (if shiny is detected).}

\item{parallel}{whether to use parallelization. Default is FALSE.\cr
When NULL, current \pkg{future}'s plan 'strategy' will be used.\cr
When FALSE, \link[future]{plan} will be called with \link[future]{sequential} 'strategy'.
When TRUE, \link[future]{plan} will be called with \link[future.callr]{callr} 'strategy'.}
}
\value{
a 3D array of features values whose dimensions are [object, features, channel] of class `IFCip_features`.
}
\description{
Function to extract basic features (hu moments + intensities) from objects stored within rif and cif files.
}
\details{
arguments of objectExtract() from IFC package will be deduced from \code{\link{ExtractBasic}} input arguments.
}
\examples{
if(!requireNamespace("IFCdata", quietly = TRUE)) {
  ## use a cif file
  file_cif <- system.file("extdata", package = "IFCdata", "example.cif")
  ## features extraction:
  ## the extraction is done for objects 1 to 50 only to allow example to run 
  ## in a reasonable amount of time and without parallelization to fulfill CRAN policies
  feat <- ExtractBasic(fileName = file_cif,
                       objects = 1:50,
                       display_progress = TRUE,
                       parallel = FALSE)
} else {
  message(sprintf('Please run `install.packages("IFCdata", repos = "\%s", type = "source")` \%s',
                  'https://gitdemont.github.io/IFCdata/',
                  'to install extra files required to run this example.'))
}
}
