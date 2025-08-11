/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2021 Yohann Demont                                              
                                                                                
  It is part of IFCip package, please cite:                                     
  -IFCip: An R Package for Imaging Flow Cytometry Image Processing              
  -YEAR: 2021                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Jean-Pierre Marolleau, Loïc Garçon,        
                      CHU Amiens                                                
                                                                                
                                                                                
  DISCLAIMER:                                                                   
  -You are using this package on your own risk!                                 
  -We do not guarantee privacy nor confidentiality.                             
  -This program is distributed in the hope that it will be useful, but WITHOUT  
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or         
  FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or  
  contributors be liable for any direct, indirect, incidental, special,         
  exemplary, or consequential damages (including, but not limited to,           
  procurement of substitute goods or services; loss of use, data, or profits;   
  or business interruption) however caused and on any theory of liability,      
  whether in contract, strict liability, or tort (including negligence or       
  otherwise) arising in any way out of the use of this software, even if        
  advised of the possibility of such damage.                                    
                                                                                
  You should have received a copy of the GNU General Public License             
  along with IFCip. If not, see <http://www.gnu.org/licenses/>.                 
*/

#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "../inst/include/utils.hpp"
#include "../inst/include/hu.hpp"
#include "../inst/include/otsu.hpp"
#include "../inst/include/distance.hpp"
#include "../inst/include/haralick.hpp"
#include "../inst/include/zernike.hpp"
#include "../inst/include/matrix_logic.hpp"
#include "../inst/include/shift.hpp"
#include "../inst/include/flip.hpp"
#include "../inst/include/filter.hpp"
#include "../inst/include/uw.hpp"
#include "../inst/include/morphology.hpp"
#include "../inst/include/geodesic.hpp"
#include "../inst/include/watershed.hpp"
#include "../inst/include/ctl.hpp"
#include "../inst/include/fill.hpp"
#include "../inst/include/mask.hpp"
#include "../inst/include/thinning.hpp"
#include "../inst/include/similarity.hpp"
#include "../inst/include/kernel.hpp"
#include "../inst/include/convexhull.hpp"
#include "../inst/include/caliper.hpp"
#include "../inst/include/scale.hpp"
using namespace Rcpp;


// FROM utils
//' @title Image Background
//' @name cpp_background
//' @description
//' This function is designed to compute image background on a "raw" image
//' @param img a NumericMatrix, containing image intensity values.
//' @param margin R_len_t number of rows margin used to compute background. Default is 4.
//' @param extra R_len_t number of extra columns used to compute background. Default is 0.
//' @param is_cif a bool whether 'ímg' originates from a cif_file or not. Default is false.
//' @return a NumericVector of background mean and sd
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_background(const Rcpp::NumericMatrix img,
                                   const R_len_t margin = 4,
                                   const R_len_t extra = 0,
                                   const bool is_cif = false) {
  return hpp_background(img, margin, extra, is_cif);
}
// END utils

// FROM caliper
//' @title Antipodal Pairs of Convex Hull
//' @name cpp_antipodalpairs
//' @description
//' Computes antipodal pairs of a convex polygon 
//' @param pts a 2-column matrix defining the locations (x and y coordinates, respectively) of points.
//' It has to be an object of class `IFCip_convexhull`
//' @source Adaptation from \url{https://escholarship.mcgill.ca/concern/theses/fx719p46g} in Computational geometry with the rotating calipers authored by Pirzadeh, Hormoz under supervision of Toussaint, Godfried T. at McGill University.
//' @return an IntegerVector of antipodal pairs of the convex input polygon. Note that this vector of indices is in C so 1st start at 0; add 1 to use it in R.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_antipodalpairs(const Rcpp::NumericMatrix pts) {
  return hpp_antipodalpairs(pts);
}

//' @title Bounding Box of Convex Hull
//' @name cpp_bbox
//' @description
//' Computes features from a Convex Hull 
//' @param pts a 2-column matrix defining the locations (x and y coordinates, respectively) of points.
//' It has to be an object of class `IFCip_convexhull`
//' @param scale a double used to scale the returned values.
//' @return a NumericVector of features from convex hull.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_bbox(const Rcpp::NumericMatrix pts,
                             const double scale = 1.0) {
  return hpp_bbox(pts, scale);
}
// END caliper

// FROM scale
//' @title Image Scaling
//' @name cpp_rescale
//' @description
//' This function is designed to scale a SEXP to [0, n_lev - 1]
//' @param img, a SEXP (logical, raw, integer or numeric) vector or matrix containing image intensity values.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as \code{true}.
//' Default is \code{R_NilValue}, for using all \code{'img'} elements without masking anything.
//' @param value, a double; it is the replacement value that will be used when \code{'msk'} element is interpreted as \code{false}. Default is \code{NA_REAL}.
//' @param n_lev, an int determining the number of levels used for the computation. Default is \code{256}.
//' @param invert, a bool determining whether '\code{'img'} should be scaled from min(\code{'img'}) to max(\code{'img'}) (when \code{'false'}, [min(\code{'img'}),max(\code{'img'})] becoming [0,sign(\code{'n_lev'})*abs(\code{'n_lev'}-1)]) or inverted (when \code{true}, with [max(\code{'img'}),min(\code{'img'})] rescaled to [0,sign(\code{'n_lev'})*\code{'n_lev'}-1]) values. Default is \code{false}.
//' @param bin, a bool determining whether \code{'img'} should be binned or if scaling should be continuous. Default is \code{false}.
//' @param clipmin, a double, minimal value under which \code{'img'} intensity values will be clipped to. Default is \code{NA_REAL}, to use no minimal clipping.
//' @param clipmax, a double, maximal value above which '\code{'img'} intensity values will be clipped to. Default is \code{NA_REAL}, to use no maximal clipping.
//' @param method, an uint8_t determining how scaling should be applied. Default is \code{1}.
//' -when \code{1}, on [min(\code{'img'}),max(\code{'img'})]
//' -when \code{2}, on [is_finite(\code{'clipmin'}) ? \code{'clipmin'} : min(\code{'img'}), is_finite(\code{'clipmax'}) ? \code{'clipmax'} : max(\code{'img'})]
//' -when \code{3}, on [is_finite(\code{'clipmin'}) ? max(\code{'clipmin'}, min(\code{'img'})) : min(\code{'img'}), is_finite(\code{'clipmax'}) ? min(\code{'clipmax'}, max(\code{'img'})) : max(\code{'img'})]
//' @details when \code{'msk'} is provided it has to be of the same dimensions as \code{'img'}, otherwise an error will be thrown.\cr
//' an error will be thrown also if \code{'msk'} contains non-finite value.\cr
//' \code{'img'} range will be determined based on indices of non 0 finite \code{'msk'} values and only the values in \code{'img'} at those indices will be scaled; the others will be filled with \code{'value'}.
//' @return a SEXP of same type as \code{'img'} with class `IFCip_rescale`
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_rescale(const SEXP img,
                 const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                 const double value = NA_REAL,
                 const int n_lev = 256,
                 const bool invert = false,
                 const bool bin = false,
                 const double clipmin = NA_REAL,
                 const double clipmax = NA_REAL,
                 const uint8_t method = 1) {
  return hpp_rescale(img, msk_, value, n_lev, invert, bin, clipmin, clipmax, method);
}

//' @title Image Reverse Scaling
//' @name cpp_scalerev
//' @description
//' This function is designed to revert scaling of a SEXP
//' @param img, a SEXP (logical, raw, integer or numeric) vector or matrix containing image intensity values.
//' @param sca_, a Rcpp::NumericVector of length 7 containing scaling information. Default is \code{R_NilValue} to use attr(img, "scale").
//' @return a SEXP of same type as \code{'img'}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_scalerev(const SEXP img,
                  const Rcpp::Nullable<Rcpp::NumericVector> sca_ = R_NilValue) {
  SEXP out = Rcpp::clone(img);
  hpp_scalerev(out, sca_);
  return out;
}
// END scale

// FROM haralick
//' @title Haralick Co-occurrence Matrix
//' @name cpp_cooc
//' @description
//' This function is designed to compute Haralick co-occurrence matrix
//' @param img a Rcpp::IntegerMatrix of class `IFCip_rescale`, containing image intensity values.
//' @param delta a Rcpp::IntegerVector of column and row shifts. If only one value is provided only row will be shifted.
//' @details See 'Textural Features for Image Classification', Haralick et. al (1979),
//' available at: \url{https://haralick.org/journals/TexturalFeatures.pdf}
//' @return a Rcpp::IntegerMatrix Gray-Level Co-occurrence Matrices.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_cooc(const Rcpp::IntegerMatrix img,
                             const Rcpp::IntegerVector delta) {
  return hpp_cooc(img, delta);
}

//' @title Haralick Features
//' @name cpp_h_features
//' @description
//' This function is designed to compute Haralick's features
//' @param cooc a Rcpp::NumericMatrix of class `IFCip_cooc`, co-occurrence matrix to compute Haralick's features from.
//' @param invariant a bool, whether to compute invariant Haralick's texture features. Default is false.
//' Not yet supported.
//' @details Haralick's invariant texture features are described in Löfstedt T, Brynolfsson P, Asklund T, Nyholm T, Garpebring A (2019) Gray-level invariant Haralick texture features.
//' PLoS ONE 14(2): e0212110. \doi{10.1371/journal.pone.0212110}
//' @return a Rcpp::NumericVector of Haralick's texture features
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_h_features(const Rcpp::IntegerMatrix cooc,
                                   const bool invariant = false) {
  return hpp_h_features(cooc, invariant);
}
// END haralick

// FROM hu
//' @title Hu's Centroid
//' @name cpp_centroid
//' @description
//' This function is designed to compute Hu's image centroids.
//' @param img a NumericMatrix, containing image intensity values.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_centroid(const Rcpp::NumericMatrix img) {
  return hpp_centroid(img);
}

//' @title Hu's Raw Moment
//' @name cpp_rmoment
//' @description
//' This function is designed to compute Hu's image raw moment.
//' @param img a NumericMatrix, containing image intensity values.
//' @param p uint8_t: p order. Default is 0.
//' @param q uint8_t: q order. Default is 0.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_rmoment(const Rcpp::NumericMatrix img, 
                                const uint8_t p = 0, 
                                const uint8_t q = 0) {
  return hpp_rmoment(img, p, q);
}

//' @title Hu's Central Moment
//' @name cpp_cmoment
//' @description
//' This function is designed to compute Hu's image central moment.
//' @param img a NumericMatrix, containing image intensity values.
//' @param cx double, x centroid of the img.
//' @param cy double, y centroid of the img.
//' @param p uint8_t: p order. Default is 0.
//' @param q uint8_t: q order. Default is 0.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_cmoment(const Rcpp::NumericMatrix img, 
                                const double cx = 0.0, 
                                const double cy = 0.0, 
                                const uint8_t p = 0, 
                                const uint8_t q = 0) {
  return hpp_cmoment(img, cx, cy, p, q);
}

//' @title Hu's Partial Features
//' @name cpp_features_hu1
//' @description
//' This function is designed to compute Hu's central moments.
//' @param img a NumericMatrix, containing image intensity values.
//' @param mag a double, magnification scale. Default is 1.0. Use:\cr
//' -1.0 for 20x\cr
//' -4.0 for 40x\cr
//' -9.0 for 60x.
//' @return a NumericVector of Hu's moments values.\cr
//' -Area: img's area\cr
//' -circularity: img's circularity\cr
//' -Minor Axis: img's ellipsis minor axis\cr
//' -Major Axis: img's ellipsis major axis\cr
//' -Aspect Ratio: img's ratio of minor axis over major axis\cr
//' -Angle: img's ellipsis angle with x axis (in radians)\cr
//' -theta: img's ellipsis theta angle (in radians)\cr
//' -Eccentricity: img's ellipsis ecentricity\cr
//' -pix cx: img's pixel x centroïd\cr
//' -pix cy: img's pixel y centroïd\cr
//' -pix min axis: img's ellipsis minor axis in pixels\cr
//' -pix_maj_axis: img's ellipsis major axis in pixels\cr
//' -pix count: img's area in pixels.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_features_hu1(const Rcpp::NumericMatrix img,
                                    const double mag = 1.0) {
  return hpp_features_hu1(img, mag);
}

//' @title Hu's Full Features
//' @name cpp_features_hu2
//' @description
//' This function is designed to compute Hu's central moments + 7 invariant moments
//' @param img a NumericMatrix, containing image intensity values.
//' @param mag a double, magnification scale. Default is 1.0. Use:\cr
//' -1.0 for 20x\cr
//' -4.0 for 40x\cr
//' -9.0 for 60x.
//' @return a NumericVector of Hu's moments values.\cr
//' -Area: img's area\cr
//' -circularity: img's circularity\cr
//' -Minor Axis: img's ellipsis minor axis\cr
//' -Major Axis: img's ellipsis major axis\cr
//' -Aspect Ratio: img's ratio of minor axis over major axis\cr
//' -Angle: img's ellipsis angle with x axis (in radians)\cr
//' -theta: img's ellipsis theta angle (in radians)\cr
//' -Eccentricity: img's ellipsis ecentricity\cr
//' -pix cx: img's pixel x centroïd\cr
//' -pix cy: img's pixel y centroïd\cr
//' -pix min axis: img's ellipsis minor axis in pixels\cr
//' -pix_maj_axis: img's ellipsis major axis in pixels\cr
//' -pix count: img's area in pixels\cr
//' -inv[1-7]: image invariant moments.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_features_hu2(const Rcpp::NumericMatrix img,
                                     const double mag = 1.0) {
  return hpp_features_hu2(img, mag);
}

//' @title Basic Features
//' @name cpp_basic
//' @description
//' This function is designed to compute very basic features based on Hu's moments + intensities.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk a LogicalMatrix, containing msk.
//' @param mag a double, magnification scale. Default is 1.0. Use:\cr
//' -1.0 for 20x\cr
//' -4.0 for 40x\cr
//' -9.0 for 60x.
//' @return a NumericVector of basic features values.\cr
//' -Area: msk's area\cr
//' -circularity: msk's circularity\cr
//' -Minor Axis: msk's ellipsis minor axis\cr
//' -Major Axis: msk's ellipsis major axis\cr
//' -Aspect Ratio: msk's ratio of minor_axis over major_axis\cr
//' -Angle: msk's ellipsis angle with x axis (in radians)\cr
//' -theta: msk's ellipsis theta angle (in radians)\cr
//' -eccentricity: msk's ellipsis ecentricity\cr
//' -pix cx: msk's pixel x centroïd\cr
//' -pix cy: msk's pixel y centroïd\cr
//' -pix min_axis: msk's ellipsis minor axis in pixels\cr
//' -pix maj axis: msk's ellipsis major axis in pixels\cr
//' -pix count: msk's area in pixels\cr
//' -Raw Mean Pixel: img's mean pixel intensity\cr
//' -Raw Min Pixel: img's minimal pixel intensity\cr
//' -Raw Max Pixel: img's maximal pixel intensity.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_basic(const Rcpp::NumericMatrix img,
                              const Rcpp::LogicalMatrix msk,
                              const double mag = 1.0) {
  return hpp_basic(img, msk, mag);
}

//' @title Image Features Extraction
//' @name cpp_features_hu3
//' @description
//' This function is designed to compute image features.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk an IntegerMatrix, containing msk components.
//' @param labels a Nullable IntegerVector corresponding to the desired label(s) to retrieve features about.
//' Default is \code{0} to retrieve features for all components.
//' @param mag a double, magnification scale. Default is \code{1.0}. Use:\cr
//' -1.0 for 20x\cr
//' -4.0 for 40x\cr
//' -9.0 for 60x.
//' @return a NumericMatrix whose rows are component numbers and columns are:\cr
//' -Area, area of the component\cr
//' -circularity, circularity of the component\cr
//' -Minor Axis, minor axis of the component\cr
//' -Major Axis, major axis of the component\cr
//' -Aspect Ratio, aspect ratio of the component\cr
//' -Angle, angle of the component\cr
//' -theta, theta of the component\cr
//' -eccentricity, eccentricity of the component\cr
//' -Minor Axis Intensity, intensity weighted minor axis of the component\cr
//' -Major Axis Intensity, intensity weighted major axis of the component\cr
//' -Aspect Ratio Intensity, intensity weighted aspect ratio of the component\cr
//' -Angle Intensity, intensity weighted angle of the component\cr
//' -theta intensity, intensity weighted theta of the component\cr
//' -eccentricity intensity, intensity weighted eccentricity of the component\cr
//' -pix cx, x centroid of the component in pixels\cr
//' -pix cy, y centroid of the component in pixels\cr
//' -pix min axis, minor axis of the component in pixels\cr
//' -pix maj axis, major axis of the component in pixels\cr
//' -pix count, number of pixels occupied by the component\cr
//' -inv[1-8], component Hu's invariant moments\cr
//' -Raw Mean Pixel, mean pixels intensity of the component\cr
//' -Raw Min Pixel, pixels intensity minimum of the component\cr
//' -Raw Max Pixel, pixels intensity maximum of the component\cr
//' -Std Dev, pixels intensity standard variation of the component\cr
//' -skewness, component's skewness\cr
//' -kurtosis, component's kurtosis\cr
//' -Centroid Y, scaled Y centroid\cr
//' -Centroid X, scaled X centroid\cr
//' -Centroid Y Intensity, intensity weighted scaled Y centroid\cr
//' -Centroid X Intensity. intensity weighted scaled X centroid.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_features_hu3(const Rcpp::NumericMatrix img,
                                     const Rcpp::IntegerMatrix msk,
                                     const Rcpp::Nullable<Rcpp::IntegerVector> labels = Rcpp::IntegerVector::create(0),
                                     const double mag = 1.0) {
  return hpp_features_hu3(img, msk, labels, mag);
}

//' @name cpp_features_hu4
//' @description
//' This function is designed to compute image features.
//' @param msk an IntegerMatrix, containing msk components.
//' @param labels a Nullable IntegerVector corresponding to the desired label(s) to retrieve features about.
//' Default is \code{0} to retrieve features for all components.
//' @param mag a double, magnification scale. Default is \code{1.0}. Use:\cr
//' -1.0 for 20x\cr
//' -4.0 for 40x\cr
//' -9.0 for 60x.
//' @return a NumericMatrix whose rows are component numbers and columns are:\cr
//' -Area, area of the component\cr
//' -circularity, circularity of the component\cr
//' -Minor Axis, minor axis of the component\cr
//' -Major Axis, major axis of the component\cr
//' -Aspect Ratio, aspect ratio of the component\cr
//' -Angle, angle of the component\cr
//' -theta, theta of the component\cr
//' -eccentricity, eccentricity of the component\cr
//' -pix cx, x centroid of the component in pixels\cr
//' -pix cy, y centroid of the component in pixels\cr
//' -pix min axis, minor axis of the component in pixels\cr
//' -pix maj axis, major axis of the component in pixels\cr
//' -pix count, number of pixels occupied by the component\cr
//' -Centroid Y, scaled Y centroid\cr
//' -Centroid X, scaled X centroid.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_features_hu4(const Rcpp::IntegerMatrix msk,
                                     const Rcpp::Nullable<Rcpp::IntegerVector> labels = Rcpp::IntegerVector::create(0),
                                     const double mag = 1.0) {
  return hpp_features_hu4(msk, labels, mag);
}
// END hu

// FROM otsu
//' @title Otsu Multi Thresholding
//' @name cpp_multi_otsu
//' @description
//' This function determines best threshold(s) according to Otsu's method.
//' @param img, a NumericVector.
//' @param msk_, a NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @param n_comp, number of components to separate. Default is 2, should be at least 2.\cr
//' Returned thresholds will be of length n_comp - 1.
//' @param n_lev, an unsigned short determining the number of grey levels used for the computation. Default is 256, should be at least 2.
//' Despite being fast thanks to LUT pre-computation, performance will be highly impacted with large 'n_comp' or 'n_lev' values (typically n_comp = 5 and n_lev = 256).
//' Alternatively, you can try to decrease 'n_lev' when 'n_comp' needs to be large (e.g. n_comp = 8 and n_lev = 32).
//' @details adaptation of 'A Fast Algorithm for Multilevel Thresholding' from L. Ping-Sung, C. Tse-Sheng, and C. Pau-Choo
//' in Jounal of Information Science and Engineering. 2001(17), 713-727.
//' \doi{10.6688/JISE.2001.17.5.1}
//' @return an NumericVector of threshold(s).
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_multi_otsu (const Rcpp::NumericVector img,
                                    const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                                    const uint8_t n_comp = 2,
                                    const unsigned short n_lev = 256) {
  return hpp_multi_otsu(img, msk_, n_comp, n_lev);
}
// END otsu

// FROM distance
//' @title Mask Euclidean Distance
//' @name cpp_distance_eucl
//' @description
//' This function is designed to compute Euclidean distance from background to centroïds' foreground
//' @param msk an IntegerMatrix, containing connected components.
//' @return an NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_distance_eucl(const Rcpp::IntegerMatrix msk) {
  return hpp_distance_eucl(msk);
}

//' @title Mask Normalized Euclidean Distance
//' @name cpp_distance_eucl_norm
//' @description
//' This function is designed to compute normalized Euclidean distance from background to centroïds' foreground
//' @param msk an IntegerMatrix, containing connected components.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_distance_eucl_norm(const Rcpp::IntegerMatrix msk) {
  return hpp_distance_eucl_norm(msk);
}

//' @title Mask Manhattan Distance
//' @name cpp_distance_manh
//' @description
//' This function is designed to compute Manhattan distance from background to centroïds' foreground
//' @param msk an IntegerMatrix, containing connected components.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_distance_manh(const Rcpp::IntegerMatrix msk) {
  return hpp_distance_manh(msk);
}

//' @title Mask Normalized Manhattan Distance
//' @name cpp_distance_manh_norm
//' @description
//' This function is designed to compute normalized Manhattan distance from background to centroïds' foreground
//' @param msk an IntegerMatrix, containing connected components.
//' @return an NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_distance_manh_norm(const Rcpp::IntegerMatrix msk) {
  return hpp_distance_manh_norm(msk);
}

//' @title Manhattan Distance Transform
//' @name cpp_disttrans_manh
//' @description
//' This function computes the Manhattan distance transform of an image by implementing A. Meijster algorithm.
//' @param img, a NumericMatrix.
//' @details adaptation of 'A General Algorithm For Computing Distance Transforms In Linear Time' from W.H. Hesselink, A. Meijster, J.B.T.M. Roerdink.
//' Mathematical Morphology and its Applications to Image and Signal Processing. February 2002, Pages 331-340.\doi{10.1007/0-306-47025-X_36}\cr
//' Values > 0 will be considered as foreground whereas all other values will be concidered as background (i.e. 0).
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_disttrans_manh (const Rcpp::NumericMatrix img) {
  return hpp_disttrans_manh(img);
}

//' @title Euclidean Distance Transform
//' @name cpp_disttrans_eucl
//' @description
//' This function computes the Euclidean distance transform of an image by implementing A. Meijster algorithm.
//' @param img, a NumericMatrix.
//' @details adaptation of 'A General Algorithm For Computing Distance Transforms In Linear Time' from W.H. Hesselink, A. Meijster, J.B.T.M. Roerdink.
//' Mathematical Morphology and its Applications to Image and Signal Processing. February 2002, Pages 331-340.\doi{10.1007/0-306-47025-X_36}\cr
//' Values > 0 will be considered as foreground whereas all other values will be concidered as background (i.e. 0).
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_disttrans_eucl (const Rcpp::NumericMatrix img) {
  return hpp_disttrans_eucl(img);
}

//' @title Euclidean Voronoï
//' @name cpp_voronoi_eucl
//' @description
//' This function computes Voronoï diagram using Euclidean distance of a seed image.
//' @param img, a positive IntegerMatrix where values > 0 represent foreground seeds.
//' @details img will be passed to connected component labelling and centroïd of each identified components will be used to build Voronoï diagram.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_voronoi_eucl (const Rcpp::IntegerMatrix img) {
  return hpp_voronoi_eucl(img);
}

//' @title Manhattan Voronoï
//' @name cpp_voronoi_manh
//' @description
//' This function computes Voronoï diagram using Manhattan distance of a seed image.
//' @param img, a positive IntegerMatrix where values > 0 represent foreground seeds.
//' @details img will be passed to connected component labelling and centroïd of each identified components will be used to build Voronoï diagram.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_voronoi_manh (const Rcpp::IntegerMatrix img) {
  return hpp_voronoi_manh(img);
}
// END distance

// FROM zernike
//' @title Zernike's Features
//' @name cpp_zernike1
//' @description
//' This function is designed to compute Zernike's moments from image.
//' It will compute Zernike's moments but will not return image projection.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk_ a Nullable LogicalMatrix. Default is R_NilValue.
//' @param cx a double. X centroid. Default is 0.0.
//' @param cy a double. Y centroid. Default is 0.0.
//' @param zmax a uint8_t, maximal order of Zernike polynomials to be computed. Default is 15. Values outside [0,99] will be clipped.
//' Be aware that computation of Zernike's moments can be quite long when 'zmax' is high.
//' @param radius a numeric, radius of the circle in pixels around object centers from which the features are calculated. Default is 15.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_zernike1(const Rcpp::NumericMatrix img,
                        const Rcpp::Nullable<Rcpp::LogicalMatrix> msk_ = R_NilValue,
                        const double cx = 0.0,
                        const double cy = 0.0,
                        const uint8_t zmax = 15,
                        const double radius = 15.0) {
  return hpp_zernike1(img, msk_, cx, cy, zmax, radius);
}

//' @title Zernike's Features with Projections
//' @name cpp_zernike2
//' @description
//' This function is designed to compute Zernike's moments from image.
//' It will compute Zernike's moments but also return image projection.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk_ a Nullable LogicalMatrix. Default is R_NilValue.
//' @param cx a double. X centroid. Default is 0.0.
//' @param cy a double. Y centroid. Default is 0.0.
//' @param zmax a uint8_t, maximal order of Zernike polynomials to be computed. Default is 15. Values outside [0,99] will be clipped.
//' Be aware that computation of Zernike's moments can be quite long when 'zmax' is high.
//' @param radius a numeric, radius of the circle in pixels around object centers from which the features are calculated. Default is 15.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_zernike2(const Rcpp::NumericMatrix img,
                        const Rcpp::Nullable<Rcpp::LogicalMatrix> msk_ = R_NilValue,
                        const double cx = 0.0, 
                        const double cy = 0.0, 
                        const uint8_t zmax = 15, 
                        const double radius = 15.0) {
  return hpp_zernike2(img, msk_, cx, cy, zmax, radius);
}
// END zernike


// FROM matrix_logic
//' @title Matrix List AND Logic
//' @name cpp_AND_M
//' @description
//' This function takes a list of matrices and returns the AND operation applied on these matrices.
//' @param list a list of logical matrices.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_AND_M(const Rcpp::List list) {
  return hpp_AND_M(list);
}

//' @title Matrix List OR Logic
//' @name cpp_OR_M
//' @description
//' This function takes a list of matrices and returns the OR operation applied on these matrices.
//' @param list a list of logical matrices.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_OR_M(const Rcpp::List list) {
  return hpp_OR_M(list);
}

//' @title Matrix Neg Logic
//' @name cpp_NEG_M
//' @description
//' This function takes a logical matrix and returns its negation.
//' @param mat LogicalMatrix.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_NEG_M(const Rcpp::LogicalMatrix mat) {
  return hpp_NEG_M(mat);
}

//' @title Matrix Logic Equal
//' @name cpp_k_equal_M
//' @description
//' This function takes an NumericMatrix and checks for members equal to k
//' @param mat a NumericMatrix
//' @param k constant to be checked. Default is 3.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_k_equal_M(const Rcpp::NumericMatrix mat, const double k = 3.0) {
  return hpp_k_equal_M(mat, k);
}

//' @title Matrix Logic Superior and Equal
//' @name cpp_k_sup_equal_M
//' @description
//' This function takes an NumericMatrix and checks for members superior or equal to k
//' @param mat a NumericMatrix
//' @param k constant to be checked. Default is 3.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_k_sup_equal_M(const Rcpp::NumericMatrix mat, const double k = 3.0) {
  return hpp_k_sup_equal_M(mat, k);
}

//' @title Matrix Logic Inferior and Equal
//' @name cpp_k_inf_equal_M
//' @description
//' This function takes an NumericMatrix and checks for members inferior or equal to k
//' @param mat a NumericMatrix
//' @param k constant to be checked. Default is 3.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_k_inf_equal_M(const Rcpp::NumericMatrix mat, const double k = 3.0) {
  return hpp_k_inf_equal_M(mat, k);
}
// END matrix_logic

// FROM shift
//' @title Image Shift
//' @name cpp_shift
//' @description
//' Function that shifts mat according to d_row and d_col parameters
//' @param mat a numeric matrix.
//' @param d_row an integer, giving row shift. Default is 0 for no change.
//' @param d_col an integer, giving col shift. Default is 0 for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions
//' Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @return a shifted matrix with additional rows/columns if d_row or d_col are different from 0.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_shift( const Rcpp::NumericMatrix mat,
                               const int d_row = 0,
                               const int d_col = 0,
                               const bool add_noise = true,
                               const double bg = 0.0,
                               const double sd = 0.0) {
  if(d_row == 0) return hpp_shift_col(mat, d_col, add_noise, bg, sd);
  if(d_col == 0) return hpp_shift_row(mat, d_row, add_noise, bg, sd);
  Rcpp::NumericMatrix M0 = hpp_shift_row(mat, d_row, add_noise, bg, sd);
  Rcpp::NumericMatrix M1 = hpp_shift_col(M0, d_col, add_noise, bg, sd);
  return M1;
}
// END shift

// FROM flip
//' @title Image Flip
//' @name cpp_flip
//' @description
//' Function that flips mat 
//' @param mat NumericMatrix.
//' @param which bool. Default is true. Use true to flip horizontally and false vertically.
//' @return a flipped matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_flip(const Rcpp::NumericMatrix mat, const bool which = true) {
  if(which) return hpp_hflip(mat);
  return hpp_vflip(mat);
}
// END flip

// FROM padding
//' @title Image Padding
//' @name cpp_padding
//' @description
//' This function creates a new matrix with extra rows / cols according to input mat, kernel
//' @param mat, a Matrix.
//' @param extra_rows,extra_cols number of extra rows and/or columns to add. Default is 0.
//' @param method, a uint8_t. Default is 1, allowed are [1-8].\cr
//' -1, extra cols / rows will be filled with 'k', returned 'out' will not be filled.\cr
//' -2, extra cols / rows will be filled with the closest col / row, returned 'out' will not be filled.\cr
//' -3, extra cols / rows will be filled mirroring neighbor cols / rows, returned 'out' will not be filled.\cr
//' -4, extra cols / rows will be filled repeating neighbor cols / rows, returned 'out' will not be filled.\cr
//' -5, extra cols / rows will be filled with 'k', returned 'out' will be filled with mat.\cr
//' -6, extra cols / rows will be filled with the closest col / row, returned 'out' will be filled with mat.\cr
//' -7, extra cols / rows will be filled mirroring neighbor cols / rows, returned 'out' will be filled with mat.\cr
//' -8, extra cols / rows will be filled repeating neighbor cols / rows, returned 'out' will be filled with mat.
//' @param k, a double, constant used when method is 1 or 5. Default is 0.0.
//' @return a Matrix of same type as 'mat', with extra cols / rows
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_padding(const SEXP mat,
                 const R_len_t extra_rows = 0,
                 const R_len_t extra_cols = 0,
                 const uint8_t method = 1,
                 const double k = 0.0) {
  return hpp_padding(mat, extra_rows, extra_cols, method, k);
}
// END padding

// FROM filter
//' @title Image Filtering
//' @name cpp_filter
//' @description
//' This function applies filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{NA_REAL}.
//' @param what, type of filtering, s std::string. Default is \code{""}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_filter(const SEXP mat,
                               const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                               const uint8_t method = 5,
                               const double k = NA_REAL,
                               const std::string what = "") {
  return hpp_filter(mat, kernel, method, k, what);
}
 
//' @title Image Standard Deviation Filtering
//' @name cpp_sd
//' @description
//' This function applies standard deviation filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{NA_REAL}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_sd(const SEXP mat,
                           const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                           const uint8_t method = 5,
                           const double k = NA_REAL) {
  return hpp_filter(mat, kernel, method, k, "sd");
}

//' @title Image Mean Filtering
//' @name cpp_mean
//' @description
//' This function applies mean filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{NA_REAL}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_mean(const SEXP mat,
                             const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                             const uint8_t method = 5,
                             const double k = NA_REAL) {
  return hpp_filter(mat, kernel, method, k, "mean");
}

//' @title Image Median Filtering
//' @name cpp_median
//' @description
//' This function applies median filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{NA_REAL}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_median(const SEXP mat,
                               const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                               const uint8_t method = 5,
                               const double k = NA_REAL) {
  return hpp_filter(mat, kernel, method, k, "median");
}

//' @title Image Mode Filtering
//' @name cpp_mode
//' @description
//' This function applies mode filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{NA_REAL}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_mode(const SEXP mat,
                             const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                             const uint8_t method = 5,
                             const double k = NA_REAL) {
  return hpp_filter(mat, kernel, method, k, "mode");
}

//' @title Image Mid Filtering
//' @name cpp_mid
//' @description
//' This function applies mid filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a NUllable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{NA_REAL}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_mid(const SEXP mat,
                            const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                            const uint8_t method = 5,
                            const double k = NA_REAL) {
  return hpp_filter(mat, kernel, method, k, "mid");
}

//' @title Image Filtering by Convolution
//' @name cpp_convolve2d
//' @description
//' This function applies 2D convolution filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{0.0}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_convolve2d(const SEXP mat,
                                   const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                                   const uint8_t method = 5,
                                   const double k = 0.0) {
  return hpp_filter(mat, kernel, method, k, "convolve");
}

//' @title Image Filtering by Correlation
//' @name cpp_correlate2d
//' @description
//' This function applies 2D correlation filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is \code{5}, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is \code{0.0}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_correlate2d(const SEXP mat,
                                    const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
                                    const uint8_t method = 5,
                                    const double k = 0.0) {
  return hpp_filter(mat, kernel, method, k, "correlate");
}
// END filter

// FROM uw
//' @title Urbach-Wilkinson Algorithm for Image Erosion and Dilation
//' @name cpp_uw
//' @description
//' This function applies erosion or dilatation on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param erode, a bool. whether to do image erosion or dilatation. Default is true to perform erosion.
//' @details see 'Efficient 2-D grayscale morphological transformations with arbitrary flat structuring elements' from  E.R. Urbach, M.H.F. Wilkinson.
//' IEEE Transactions on Image Processing, 17(1):1-8, January 2008.\doi{10.1109/tip.2007.912582}
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_uw(const SEXP mat,
            const Rcpp::Nullable<Rcpp::NumericMatrix> kernel,
            const bool erode = true) {
  return hpp_uw(mat, kernel, erode);
}
// END uw

// FROM morphology
//' @title Image Erosion
//' @name cpp_erode
//' @description
//' This function applies erosion on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @details see 'Efficient 2-D grayscale morphological transformations with arbitrary flat structuring elements' from  E.R. Urbach, M.H.F. Wilkinson.
//' IEEE Transactions on Image Processing, 17(1):1-8, January 2008.\doi{10.1109/tip.2007.912582}
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_erode(const SEXP mat,
               const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_erode(mat, kernel);
}

//' @title Image Dilatation
//' @name cpp_dilate
//' @description
//' This function applies dilatation on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @details see 'Efficient 2-D grayscale morphological transformations with arbitrary flat structuring elements' from  E.R. Urbach, M.H.F. Wilkinson.
//' IEEE Transactions on Image Processing, 17(1):1-8, January 2008.\doi{10.1109/tip.2007.912582}
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_dilate(const SEXP mat,
                const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_dilate(mat, kernel);
}

//' @title Image Opening
//' @name cpp_opening
//' @description
//' This function applies opening on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_opening(const SEXP mat,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_opening(mat, kernel);
}

//' @title Image Closing
//' @name cpp_closing
//' @description
//' This function applies closing on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_closing(const SEXP mat,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_closing(mat, kernel);
}

//' @title Image Morphological Gradient
//' @name cpp_gradient
//' @description
//' This function applies morphological gradient on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_gradient(const SEXP mat,
                  const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_gradient(mat, kernel);
}

//' @title Image White Top Hat
//' @name cpp_tophat_white
//' @description
//' This function applies white top hat on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_tophat_white(SEXP mat,
                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_tophat_white(mat, kernel);
}

//' @title Image Black Top Hat
//' @name cpp_tophat_black
//' @description
//' This function applies black top hat on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_tophat_black(SEXP mat,
                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_tophat_black(mat, kernel);
}

//' @title Image Self Complementary Top Hat
//' @name cpp_tophat_self
//' @description
//' This function applies self complementary on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_tophat_self(SEXP mat,
                     const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_tophat_self(mat, kernel);
}

//' @title Image Contrast Enhancement
//' @name cpp_cont
//' @description
//' This function applies contrast enhancement on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_cont(const SEXP mat,
              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_cont(mat, kernel);
}

//' @title Image Laplacian
//' @name cpp_laplacian
//' @description
//' This function applies Laplacian morphology on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_laplacian(const SEXP mat,
                   const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  return hpp_laplacian(mat, kernel);
}

//' @title Brute Force Image Erosion
//' @name cpp_erosion_old
//' @description
//' This function applies erosion on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time erosion should be iterated. Default is 0.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @details Brute force implementation now replaced by Urbach-Wilkinson algorithm.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_erode_old(const Rcpp::NumericMatrix mat,
                                  const Rcpp::NumericMatrix kernel,
                                  const uint8_t iter = 0,
                                  const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue) {
  return hpp_erode_old(mat, kernel, iter, msk_);
}

//' @title Brute Force Image Dilatation
//' @name cpp_dilate_old
//' @description
//' This function applies dilatation on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilatation should be iterated. Default is 0.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @details Brute force implementation now replaced by Urbach-Wilkinson algorithm.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_dilate_old(const Rcpp::NumericMatrix mat,
                                   const Rcpp::NumericMatrix kernel,
                                   const uint8_t iter = 0,
                                   const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue) {
  return hpp_dilate_old(mat, kernel, iter, msk_);
}

//' @title Contours Dilatation
//' @name cpp_dilate_ctl
//' @description
//' This function applies contours dilatation.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_dilate_ctl(const List ctl,
                                   const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_dilate_ctl: 'ctl' should be of class `IFCip_ctl`");
  }
  return hpp_dilate(as<Rcpp::NumericMatrix>(ctl["matrix"]), kernel);
}

//' @title Contours Erosion
//' @name cpp_erode_ctl
//' @description
//' This function applies contours erosion.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_erode_ctl(const List ctl,
                                  const Rcpp::Nullable<Rcpp::NumericMatrix> kernel) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_erode_ctl: 'ctl' should be of class `IFCip_ctl`");
  }
  return hpp_erode(as<Rcpp::NumericMatrix>(ctl["matrix"]), kernel);
}
// END morphology

// FROM geodesic
//' @title Dilatation Reconstruction
//' @name cpp_rec_dilate
//' @description
//' Performs a dilatation reconstruction of an image.
//' @param markers, a NumericMatrix.
//' @param img, a NumericMatrix.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details adaptation of 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}
//' @return a NumericMatrix of reconstructed 'img' from 'markers'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_rec_dilate (const Rcpp::NumericMatrix markers,
                                    const Rcpp::NumericMatrix img,
                                    const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix out = Rcpp::clone(markers);
  rec_dilate(out, img, kernel);
  return out;
}

//' @title Erosion Reconstruction
//' @name cpp_rec_erode
//' @description
//' Performs an erosion reconstruction of an image.
//' @param markers, a NumericMatrix.
//' @param img, a NumericMatrix.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details adaptation of 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}
//' @return a NumericMatrix of reconstructed 'img' from 'markers'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_rec_erode (const Rcpp::NumericMatrix markers,
                                   const Rcpp::NumericMatrix img,
                                   const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix out = Rcpp::clone(markers);
  rec_erode(out, img, kernel);
  return out;
}

//' @title H-Minima transformation
//' @name cpp_HMIN
//' @description
//' Keep deepest valleys from image.
//' @param img, a NumericMatrix.
//' @param h, a double, specifying the minimal depth. Default is \code{NA_REAL}. When not \code{NA/NaN} it will be used instead of 'h_lev'
//' @param h_lev, an int, specifying the minimal depth normalized to n_lev (being h_lev out of n_lev). Default is \code{1}.
//' @param n_lev, an int determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' HMIN is the erosion reconstruction of (img + h) by kernel.
//' @return a NumericMatrix of H-Minima transformation of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_HMIN (const Rcpp::NumericMatrix img,
                              const double h = NA_REAL,
                              const int h_lev = 1,
                              const int n_lev = 65536,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_HMIN(img, h, h_lev, n_lev, kernel, msk_);
}

//' @title H-Maxima transformation
//' @name cpp_HMAX
//' @description
//' Keep highest peaks in image.
//' @param img, a NumericMatrix.
//' @param h, a double, specifying the minimal height. Default is \code{NA_REAL}. When not \code{NA/NaN} it will be used instead of 'h_lev'
//' @param h_lev, an int, specifying the minimal height normalized to n_lev (being h_lev out of n_lev). Default is \code{1}.
//' @param n_lev, an int determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' HMAX is the dilatation reconstruction of (img - h) by kernel.
//' @return a NumericMatrix of H-Maxima transformation of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_HMAX (const Rcpp::NumericMatrix img,
                              const double h = NA_REAL,
                              const int h_lev = 1,
                              const int n_lev = 65536,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_HMAX(img, h, h_lev, n_lev, kernel, msk_);
}

//' @title Regional Minima
//' @name cpp_RMIN
//' @description
//' Mask connected component of pixels whose values are lower to their external boundaries neighborhood.
//' @param img, a NumericMatrix.
//' @param n_lev, an int determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' RMIN is defined as img < HMIN(img, h_lev = 1).
//' @return a NumericMatrix of regional minima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_RMIN (const Rcpp::NumericMatrix img,
                              const int n_lev = 65536,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_RMIN(img, n_lev, kernel, msk_);
}

//' @title Regional Maxima
//' @name cpp_RMAX
//' @description
//' Mask connected component of pixels whose values are higher to their external boundaries neighborhood.
//' @param img, a NumericMatrix.
//' @param n_lev, an int determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' RMAX is defined as img > HMAX(img, h_lev = 1).
//' @return a NumericMatrix of regional maxima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_RMAX (const Rcpp::NumericMatrix img,
                              const int n_lev = 65536,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_RMAX(img, n_lev, kernel, msk_);
}

//' @title Extended Minima
//' @name cpp_EMIN
//' @description
//' Mask the regional minima of the corresponding h-minima transformation.
//' @param img, a NumericMatrix.
//' @param h, a double, specifying the minimal depth. Default is \code{NA_REAL}. When not \code{NA/NaN} it will be used instead of 'h_lev'
//' @param h_lev, an int, specifying the minimal depth normalized to n_lev (being h_lev out of n_lev). Default is \code{1}.
//' @param n_lev, an unsigned short determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' EMIN is defined as RMIN(HMIN(img, h_lev = 1)).
//' @return a NumericMatrix of extended minima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_EMAX (const Rcpp::NumericMatrix img,
                              const double h = NA_REAL,
                              const int h_lev = 1,
                              const int n_lev = 65536,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_RMAX(hpp_HMAX(img, h, h_lev, n_lev, kernel, msk_), n_lev, kernel, msk_);
}

//' @title Extended Maxima
//' @name cpp_EMAX
//' @description
//' Mask the regional maxima of the corresponding h-maxima transformation.
//' @param img, a NumericMatrix.
//' @param h, a double, specifying the minimal height. Default is \code{NA_REAL}. When not \code{NA/NaN} it will be used instead of 'h_lev'
//' @param h_lev, an int, specifying the minimal height normalized to n_lev (being h_lev out of n_lev). Default is \code{1}.
//' @param n_lev, an unsigned short determining the number levels used for 'img' rescaling. Default is 65536, should be at least 2.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' EMAX is defined as RMAX(HMAX(img, h_lev = 1)).
//' @return a NumericMatrix of extended maxima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_EMIN (const Rcpp::NumericMatrix img,
                              const double h = NA_REAL,
                              const int h_lev = 1,
                              const int n_lev = 65536,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_RMIN(hpp_HMIN(img, h, h_lev, n_lev, kernel, msk_), n_lev, kernel, msk_);
}

//' @title Geodesic White Top Hat
//' @name cpp_geo_tophat_white
//' @description
//' This function applies geodesic white top hat on image.
//' @param img, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' Dilation in closing process is replaced by dilation reconstruction.
//' So, we have out = img - rec_dilate(erode(img)).
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_geo_tophat_white (const Rcpp::NumericMatrix img,
                                          const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_geo_tophat_white(img, kernel);
}

//' @title Geodesic Black Top Hat
//' @name cpp_geo_tophat_black
//' @description
//' This function applies geodesic black top hat on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' Erosion in opening process is replaced by erosion reconstruction.
//' So, we have out = rec_erode(dilate(img)) - img.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_geo_tophat_black (const Rcpp::NumericMatrix img,
                                          const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_geo_tophat_black(img, kernel);
}
// END geodesic

// FROM watershed
//' @title Watershed Transformation SV1
//' @name cpp_watershed_sv1
//' @description
//' This function computes the watershed transformation of an image.
//' @param mat, a NumericMatrix; a distance transform matrix is expected.
//' @param n_lev, an unsigned short determining the number of elevation levels. Default is 256, should be at least 2.
//' @param draw_lines, a bool; whether to draw watershed lines or not. Default is true.
//' @param invert, a bool; whether to fill from basins (lowest values) to peaks (highest values). Default is false.
//' When 'mat' is the result of the distance transformation of an image, peaks (highest values) represent largest distances from background.
//' Thus, they are the ones to be filled first; this can be done with 'invert' set to true.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @details adaptation of 'Determining watersheds in digital pictures via flooding simulations' from P. Soille. and L. Vincent.
//' In Proc. SPIE 1360, Visual Communications and Image Processing '90: Fifth in a Series, (1 September 1990) \doi{10.1117/12.24211}.
//' @source MorphoLib plugin for ImageJ presents a Java implementation of the algorithm in  \url{https://github.com/ijpb/MorphoLibJ/blob/master/src/main/java/inra/ijpb/watershed/WatershedTransform2D.java} authored by Ignacio Arganda-Carreras 
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector cpp_watershed_sv1(const Rcpp::NumericMatrix mat,
                                      const unsigned short n_lev = 256,
                                      const bool draw_lines = true,
                                      const bool invert = false,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                                      const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_watershed_sv1(mat, n_lev, draw_lines, invert, kernel, msk_);
}

//' @title Watershed Transformation SV2
//' @name cpp_watershed_sv2
//' @description
//' This function computes the watershed transformation of an image.
//' @param mat, a NumericMatrix; a distance transform matrix is expected.
//' @param n_lev, an unsigned short determining the number of elevation levels. Default is 256, should be at least 2.
//' @param draw_lines, a bool; whether to draw watershed lines or not. Default is true.
//' @param invert, a bool; whether to fill from basins (lowest values) to peaks (highest values). Default is false.
//' When 'mat' is the result of the distance transformation of an image, peaks (highest values) represent largest distances from background.
//' Thus, they are the ones to be filled first; this can be done with 'invert' set to true.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @details adaptation of 'Watersheds in digital spaces: an efficient algorithm based on immersion simulations' from  L. Vincent and P. Soille.
//' In IEEE Transactions on Pattern Analysis and Machine Intelligence, 13(6):583-598, June 1991.\cr
//' @source The algorithm is reviewed in 'The Watershed Transform: Definitions, Algorithms and Parallelization Strategies'
//' from Roerdink, J. B. T. M. and Meijster, A. (2000) in Fundamenta Informaticae, 41, 187-228 \doi{10.3233/FI-2000-411207}
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector cpp_watershed_sv2(const Rcpp::NumericMatrix mat,
                                      const unsigned short n_lev = 256,
                                      const bool draw_lines = true,
                                      const bool invert = false,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                                      const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  return hpp_watershed_sv2(mat, n_lev, draw_lines, invert, kernel, msk_);
}
// END watershed

// FROM ctl
//' @title Contour Tracing Connected Component Labeling
//' @name cpp_ctl
//' @description
//' This function is designed to identify connected component.
//' @param mat a Matrix, converted to LogicalMatrix, where finite non zero values will be considered as mask.
//' @param global whether to compute the perimeter globally or to evaluate the perimeter of each non 8-connected objects. Default is false.
//' When true pixels of overlapping extra borders of objects are counted only once.
//' @details adaptation of 'A linear-time component-labeling algorithm using contour tracing technique' from F. Chang, C.J. Chen and C.J. Lu.
//' Computer Vision and Image Understanding Volume 93, Issue 2, February 2004, Pages 206-220.\doi{10.1016/j.cviu.2003.09.002}
//' @return a list whose members are:\cr
//' -matrix: an IntegerMatrix with connected component labels.\cr
//' -contours: an IntegerMatrix of identified contours, whose columns are x, y, label, direction and type.\cr
//' -nb_lab: the total number of components identified.
//' -perimeter: the number of pixels outside contours.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_ctl(const Rcpp::LogicalMatrix mat,
                   const bool global = false) {
  return hpp_ctl(mat, global);
}
// END ctl

// FROM fill
//' @title Polygon Drawing
//' @name cpp_polydraw
//' @description
//' This function is designed to trace and fill polygon inside a matrix.
//' @param poly, a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param border, an int used to trace polygon border.
//' @param fill, an int used to fill polygon.
//' @param tol, a double, tolerance between fill color and connected pixels. Use \code{NA}, for filling every pixel inside 'poly'.
//' @param edge, a bool whether to close 'poly' at 'mat_' edges. Default is \code{false}. Closing 'poly' at 'edge' is experimental and may fail.
//' @param mat_, a NumericMatrix to be filled.\cr
//' When 'mat_' is provided 'poly' will be drawn in 'mat_' if its vertices are within 'mat_' dimensions and if there is no \code{NA} directly (4-connectedness) connected to internal poly' border.\cr
//' /!\ Note that filling will not propagate on \code{NA}, \code{NaN} values, unless 'tol' is \code{NA}.
//' @return copy of 'mat_' with 'poly' or a new matrix with 'poly'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_polydraw (const Rcpp::IntegerMatrix poly,
                                  const double border = 1.0,
                                  const double fill = 1.0,
                                  const double tol = 0.0,
                                  const bool edge = false,
                                  const Rcpp::Nullable<Rcpp::NumericMatrix> mat_ = R_NilValue) {
  return hpp_polydraw(poly, border, fill, tol, edge, mat_);
}

//' @title Contours Filling
//' @name cpp_fill
//' @description
//' This function is designed to fill contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param labels a Nullable IntegerVector corresponding to the label(s) of desired set of contour to be filled.
//' Default is \code{0} to fill all sets of contours found.
//' @param i_border a bool, to whether or not draw inside contours if some were identified. Default is \code{true}.
//' @param i_fill a bool, to whether or not fill inside contours if some were identified. Default is \code{true}.
//' @param i_neg_border a bool, to whether or not inside border, if drawn, should be negated. Default is \code{false}.
//' @param o_border a bool, to whether or not draw external contours. Default is \code{true}.
//' @param o_fill a bool, to whether or not fill external contours. Default is \code{true}.
//' @param o_neg_border a bool, to whether or not external border, if drawn, should be negated. Default is \code{false}.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_fill(const List ctl,
                             const Rcpp::Nullable<Rcpp::IntegerVector> labels = Rcpp::IntegerVector::create(0),
                             const bool i_border = true,
                             const bool i_fill = true,
                             const bool i_neg_border = false,
                             const bool o_border = true,
                             const bool o_fill = true,
                             const bool o_neg_border = false ) {
  return hpp_fill(ctl, labels, i_border, i_fill, i_neg_border, o_border, o_fill, o_neg_border);
}

//' @title Contours Default Filling
//' @name cpp_fill_default
//' @description
//' This function is designed to apply default contours filling.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param labels a Nullable IntegerVector corresponding to the label(s) of desired set of contour to be filled.
//' Default is \code{0} to fill all sets of contours found.
//' @param i_neg_border a bool, to whether or not inside border, if drawn, should be negated. Default is \code{false}.
//' @param o_neg_border a bool, to whether or not external border, if drawn, should be negated. Default is \code{false}.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_fill_default (const List ctl,
                                      const Rcpp::Nullable<Rcpp::IntegerVector> labels = Rcpp::IntegerVector::create(0),
                                      const bool i_neg_border = false,
                                      const bool o_neg_border = false) {
  return hpp_fill_default(ctl, labels, i_neg_border, o_neg_border);
}

//' @title Contours Filling Outer Only
//' @name cpp_fill_out
//' @description
//' This function is designed to fill the most external contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`.
//' @param labels a Nullable IntegerVector corresponding to the label(s) of desired set of contour to be filled.
//' Default is \code{0} to fill all sets of contours found.
//' @param o_border a bool, to whether or not draw external contours. Default is \code{true}.
//' @param o_fill a bool, to whether or not fill external contours. Default is \code{true}.
//' @param o_neg_border a bool, to whether or not external border, if drawn, should be negated. Default is \code{false}.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_fill_out(const List ctl,
                                 const Rcpp::Nullable<Rcpp::IntegerVector> labels = Rcpp::IntegerVector::create(0),
                                 const bool o_border = true,
                                 const bool o_fill = true,
                                 const bool o_neg_border = false) {
  return hpp_fill_out(ctl, labels, o_border, o_fill, o_neg_border);
}

//' @title Connected Region Flood Filling
//' @name cpp_floodfill
//' @description
//' Flood fills image region.
//' @param img, a NumericMatrix. The image to be modified.
//' @param markers, a NumericMatrix, It should be a matrix with at least 3 columns being "row", "col", and "value", respectively. It represents coordinates of the seeds to start filling 'img', with the new "value". Eventually, an additional column being "tolerance" can be provided.\cr
//' /!\ Note that "row" and "col" should be provided at C-level meaning 1st start at 0.
//' @return a NumericMatrix, the modified image.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_floodfill (const Rcpp::NumericMatrix img,
                                   const Rcpp::NumericMatrix markers) {
  return hpp_floodfill(img, markers);
}
// END fill

// FROM mask
//' @title Matrix Threshold
//' @name cpp_threshold
//' @description
//' This function takes a image and checks for members superior or equal
//' to max(img) - k * diff(range(img)) / 100 within msk
//' @param img a NumericMatrix
//' @param msk a NumericMatrix
//' @param k constant to be checked. Default is 0.
//' @param removal uint8_t, object removal method. Default is 0 for no removal. Otherwise, if\cr
//' -1, for clipped removal, keep non clipped foreground.\cr
//' -2, height clipped removal.\cr
//' -3, width clipped removal.\cr
//' -4, only keep background:.\cr
//' -5, only keep foreground.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_threshold(const Rcpp::NumericMatrix img,
                                  const Rcpp::NumericMatrix msk,
                                  const double k = 0.0,
                                  uint8_t removal = 0) {
  return hpp_threshold(img, msk, k, removal);
}
// END mask

// FROM thinning
//' @title Implementation of Zheng-Suen Thinning
//' @name cpp_thinning_zs
//' @description
//' This function is designed to identify mask skeleton.
//' @param mat a LogicalMatrix, containing mask.
//' @details adaptation of 'A fast parallel algorithm for thinning digital patterns' from T. Y. Zhang, C. Y. Suen.
//' Communications of the ACM, March 1984 \doi{10.1145/357994.358023}.
//' @return a LogicalMatrix with the mask thinned.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_thinning_zs(const Rcpp::LogicalMatrix mat) {
  return hpp_thinning_zs(mat);
}
//' @title Implementation of Ben Boudaoud-Sider-Tari Thinning
//' @name cpp_thinning_bst
//' @description
//' This function is designed to identify mask skeleton.
//' @param mat a LogicalMatrix, containing mask.
//' @details adaptation of 'A new thinning algorithm for binary images' from L. Ben Boudaoud, A. Sider, A. Tari.
//' 3rd international conference on control, engineering & information technology, May 2015. \doi{10.1109/CEIT.2015.7233099}.
//' @return a LogicalMatrix with the mask thinned.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_thinning_bst(const Rcpp::LogicalMatrix mat) {
  return hpp_thinning_bst(mat);
}
// END thinning

// FROM similarity
//' @title Images Similarity Measurement
//' @name cpp_similarity
//' @description
//' This function is designed to score similarity between two images.
//' @param img1 a NumericMatrix, containing image values.
//' @param img2 a NumericMatrix, containing image values.
//' @param msk a LogicalMatrix, containing mask.
//' @details the similarity is the log transformed Pearson's Correlation Coefficient.
//' It is a measure of the degree to which two images are linearly correlated within a masked region.\cr
//' See "Quantitative measurement of nuclear translocation events using similarity analysis of multispectral cellular images obtained in flow"
//' by T.C. George et al. Journal of Immunological Methods Volume 311, Issues 1–2, 20 April 2006, Pages 117-129 \doi{doi.org/10.1016/j.jim.2006.01.018}
//' @return a double, the similarity.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
double cpp_similarity(const Rcpp::NumericMatrix img1,
                      const Rcpp::NumericMatrix img2,
                      const Rcpp::LogicalMatrix msk) {
  return hpp_similarity(img1, img2, msk);
}

//' @title Images Bright Detail Similarity Measurement
//' @name cpp_bright_similarity
//' @description
//' This function is designed to score similarity between two bright detail images.
//' @param img1 a NumericMatrix, containing bright detail image values.
//' @param img2 a NumericMatrix, containing bright detail image values.
//' @param msk a LogicalMatrix, containing mask.
//' @details the bright detail similarity is the non-mean normalized version of the log transformed Pearson's Correlation Coefficient.
//' It is designed to compare the small bright image detail of two images within a masked region.\cr
//' See "Quantitative analysis of protein co-localization on B cells opsonized with rituximab and complement using the ImageStream multispectral imaging flow cytometer"
//' by P.V. Beum et al. Journal of Immunological Methods Volume 317, Issues 1–2, 20 December 2006, Pages 90-99 \doi{doi.org/10.1016/j.jim.2006.09.012}
//' @return a double, the similarity.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
double cpp_bright_similarity(const Rcpp::NumericMatrix img1,
                             const Rcpp::NumericMatrix img2,
                             const Rcpp::LogicalMatrix msk) {
  return hpp_bright_similarity(img1, img2, msk);
}
// END similarity

// FROM kernel
//' @title Create Gaussian Kernel
//' @name cpp_make_gaussian
//' @description
//' This function is designed to create a gaussian kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @param sigma, a double, deviation of the filter used. Default is \code{-0.3}. If negative, \code{'sigma'} will be determined using \code{-1.0 * sigma * ((size - 1) * 0.5)}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_make_gaussian(const uint8_t size = 3,
                                      const double sigma = -0.3) {
  return hpp_make_gaussian(size, sigma);
}
//' @title Create Laplacian Kernel
//' @name cpp_make_laplacian
//' @description
//' This function is designed to create a laplacian kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @param sigma, a double, deviation of the filter used. Default is \code{-0.3}. If negative, \code{'sigma'} will be determined using \code{-1.0 * sigma * ((size - 1) * 0.5)}.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_make_laplacian(const uint8_t size = 3,
                                       const double sigma = -0.3) {
  return hpp_make_laplacian(size, sigma);
}
//' @title Create a Disc
//' @name cpp_make_disc
//' @description
//' This function is designed to create a disc kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @return a LogicalMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_make_disc(const uint8_t size = 3) {
  return hpp_make_disc(size);
}
//' @title Create a Box
//' @name cpp_make_box
//' @description
//' This function is designed to create a box kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @return a LogicalMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_make_box(const uint8_t size = 3) {
  return hpp_make_box(size);
}
//' @title Create a Plus
//' @name cpp_make_plus
//' @description
//' This function is designed to create a plus kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @return a LogicalMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_make_plus(const uint8_t size = 3) {
  return hpp_make_plus(size);
}
//' @title Create a Cross
//' @name cpp_make_cross
//' @description
//' This function is designed to create a cross kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @return a LogicalMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_make_cross(const uint8_t size = 3) {
  return hpp_make_cross(size);
}
//' @title Create a Diamond
//' @name cpp_make_diamond
//' @description
//' This function is designed to create a diamond kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @return a LogicalMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix cpp_make_diamond(const uint8_t size = 3) {
  return hpp_make_diamond(size);
}
// END kernel

// FROM convexhull
//' @title Convex Hull
//' @name cpp_convexhull
//' @description
//' Computes 2D convex hull of a set of points.
//' @param pts a 2-column matrix defining the locations (x and y coordinates, respectively) of points.
//' x has to be in column 0 and y in column 1 (C index, add 1 for R).
//' @return a vector of row indices of 'pts' that constitutes convex hull vertices.
//' @source Adaptation of Andrew's Monotone Chain Algorithm A. M. Andrew, 'Another Efficient Algorithm
//' for Convex Hulls in Two Dimension', Information Processing Letters, 9, 1979, pp216-219 reported
//' in M. A. Jayaram, Hasan Fleyeh, 'Convex Hulls in Image Processing: A Scoping Review', 
//' American Journal of Intelligent Systems, Vol. 6 No. 2, 2016, pp. 48-58.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_convexhull(const Rcpp::NumericMatrix pts) {
  return hpp_convexhull(pts);
}
// END convexhull
