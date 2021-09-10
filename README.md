<!-- badges: start -->
<!--[![R build status](https://github.com/gitdemont/IFCip/workflows/R-CMD-check/badge.svg)](https://github.com/gitdemont/IFCip/actions) -->
<!-- badges: end -->

# Image Processing Tools for IFC (IFCip)

## INSTALLATION (from **in-dev** github master branch)

IFCip package is under development and can be installed from github

#### On Windows (tested 7 and 10)

- install [Rtools for windows](https://cran.r-project.org/bin/windows/Rtools/)

- ensure that Rtools compiler is present in Windows PATH

In R console, if you installed Rtools directly on C: you should see something like C:\\Rtools\\bin and C:\\Rtools\\mingw_32\\bin 

```R
print(unlist(strsplit(Sys.getenv("PATH"), ";")))
```

Otherwise, try to set it

```R
shell('setx PATH "C:\\Rtools\\bin"')
shell('setx PATH "C:\\Rtools\\mingw_32\\bin"')
```

#### in R

- install devtools
```R
install.packages("devtools")
```

- install R dependencies required for IFCip package
"Rcpp", "IFC"

```R
install.packages(c("Rcpp", "IFC"))
```

- install "remotes", to install IFCip package from github remotes is needed.

```R
install.packages("remotes")
```

- install IFCip

```R
remotes::install_github(repo = "gitdemont/IFCip", ref = "master", dependencies = FALSE)
```

## USAGE

Several examples in IFCip package are dependent on data files that can be found in dedicated IFCdata package.

To install IFCdata package and run examples in IFCip:

```R
install.packages("IFCdata", repos = "https://gitdemont.github.io/IFCdata/", type = "source")
```

## DETAILS

IFCip contains a lots of exported and internal functions to deeply process IFC images. 

IFCip is mostly written using `Rcpp` which allows very fast processing of the images. 

On top of all the tools in IFCip, **ExtractFeatures()** is certainly the most useful one since it can be use to easily compute a bulk of features from IFC files with minimal input from user. In addition, **ExtractFeatures()** has been designed to use parallel backend allowing for even faster computation.

**ExtractFeatures()** derives, extends, or relies on:

- **moments_Hu()** to compute [`Hu moments`](https://doi.org/10.1109/TIT.1962.1057692) allowing to extract:
  - Area
  - circularity according to [Joviša Žunić _et_ al.](https://doi.org/10.1016/j.patcog.2009.06.017)
  - Minor Axis
  - Major Axis
  - Aspect Ratio
  - Angle
  - theta
  - eccentricity
  - pix cx: x centroïd in pixels
  - pix cy: y centroïd in pixels
  - pix min axis: minor axis in pixels
  - pix maj axis: major axis in pixels
  - pix count: area in pixels
  - image invariant moments
  - skewness
  - kurtosis
  - and intensity weighted values
  
  
  
- **moments_Zernike()** to compute `Zernike moments` by projecting image on Zernike complex polynomials [EBImage v3.12.0](https://www.bioconductor.org/packages/2.10/bioc/html/EBImage.html)



- **compute_shape()** to retrieve several `shape features` based notably on `contour tracing labeling` [Change, CJ.  _et_ al.](https://doi.org/10.1016/j.cviu.2003.09.002), `convex hull` computation from Andrew's Monotone Chain Algorithm reported by [M. A. Jayaram, Hasan Fleyeh](http://article.sapub.org/10.5923.j.ajis.20160602.03.html) and bounding box determination thanks to `rotating calipers` as described by [Godfried T. Toussaint](https://escholarship.mcgill.ca/concern/theses/fx719p46g):
  - Perimeter
  - Diameter
  - Circularity
  - convexity
  - roundness
  - Height
  - Width
  - Elongatedness
  - convex perimeter
  - convex cx
  - convex cy
  
  
  
- **compute_haralick()** to extract `Haralick textural features` [see Haralick's original paper](http://haralick.org/journals/TexturalFeatures.pdf) and [Löfstedt T. _et_ al.](https://doi.org/10.1371/journal.pone.0212110) used for the computation:
  - autocorrelation
  - H Contrast
  - H Correlation
  - dissimilarity
  - H Energy
  - H Entropy
  - homogeneity
  - maximum probability
  - inverse difference
  - difference entropy
  - difference variance
  - H Homogeneity
  - sum entropy
  - sum variance
  - cluster prominence
  - cluster shade
  - imc1
  - imc2
  
  
  
- **compute_intensity()** to extract image values of:
  - Bkgd Mean/StdDev
  - Min/Max/Mean Pixel
  - Raw Min/Max/Mean Pixel
  - Intensity

## FURTHER WORDS

IFCip is a new born and under active development to optimize computations, include new features and improve.

Current API may change in the future resulting in function renaming and input or output modifications 

Please take this into account and report bugs or address features implementations to help the development.

## DISCLAMER

- You are using this package **on your own risk!**

- We do not guarantee **privacy nor confidentiality**.

- This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**. In no event shall the copyright holders or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
