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
#include "../inst/include/morphology.hpp"
#include "../inst/include/geodesic.hpp"
#include "../inst/include/watershed.hpp"
#include "../inst/include/ctl.hpp"
#include "../inst/include/mask.hpp"
#include "../inst/include/thinning.hpp"
#include "../inst/include/similarity.hpp"
#include "../inst/include/kernel.hpp"
#include "../inst/include/convexhull.hpp"
#include "../inst/include/caliper.hpp"
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
// [[Rcpp::export]]
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
//' @return an IntegerVector of antipodal pairs of the convex input polygon. Note that this vector of indices is in C sort 1st start at 0; add 1 to use it in R.
//' @keywords internal
////' @export
// [[Rcpp::export]]
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
// [[Rcpp::export]]
Rcpp::NumericVector cpp_bbox(const Rcpp::NumericMatrix pts,
                             const double scale = 1.0) {
  return hpp_bbox(pts, scale);
}
// END caliper

// FROM haralick
//' @title Right Shift Matrix
//' @name cpp_R_shift_M
//' @description
//' This function is designed to right shift bits of a matrix to [0, 2^bits - 1]
//' @param mat a Rcpp::IntegerMatrix, containing image intensity values.
//' @param bits uint8_t number of bit to shift matrix values. Default is 4. Allowed are [2,10].
//' If shifted value does not respect [0, 2^bits - 1] an error is thrown.
//' @return an Rcpp::IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_R_shift_M(const Rcpp::IntegerMatrix mat,
                                  const uint8_t bits = 4) {
  return hpp_R_shift_M(mat, bits);
}

//' @title hpp_rescale_M
//' @name cpp_rescale_M
//' @description
//' This function is designed to rescale a matrix to [0, 2^bits - 1]
//' @param mat a Rcpp::NumericMatrix, containing image intensity values.
//' @param bits uint8_t number of bit to shift matrix values. Default is 4. Allowed are [2,10].
//' Rescaled values will normalized to [0, 2^bits - 1]
//' @return an Rcpp::IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_rescale_M(const Rcpp::NumericMatrix mat,
                                  const uint8_t bits = 4) {
  return hpp_rescale_M(mat, bits);
}


//' @title Haralick Co-Occurrence Matrix
//' @name cpp_cooc
//' @description
//' This function is designed to compute Haralick co-occurrence matrix
//' @param img a Rcpp::IntegerMatrix of class `IFCip_rescale`, containing image intensity values.
//' @param msk a LogicalMatrix, containing mask.
//' @param delta uint8_t offset from which co-occurence has to be computed to. Default is 1.
//' @details See 'Textural Features for Image Classification', Haralick et. al (1979),
//' available at: \url{https://haralick.org/journals/TexturalFeatures.pdf}
//' @return a list whose members are normalized Gray-Level Co-occurrence Matrices at angles 0, 45, 90 and 315.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List cpp_cooc(const Rcpp::IntegerMatrix img,
                    const Rcpp::LogicalMatrix msk,
                    const uint8_t delta = 1) {
  return hpp_cooc(img, msk, delta);
}

//' @title Haralick Features
//' @name cpp_h_features
//' @description
//' This function is designed to compute Haralick's features
//' @param cooc a Rcpp::NumericMatrix of class `IFCip_cooc`, normalized co-occurrence matrix to compute Haralick's features from.
//' @param invariant a bool, whether to compute invariant Haralick's texture features. Default is false.
//' Not yet supported.
//' @details Haralick's invariant texture features are described in Löfstedt T, Brynolfsson P, Asklund T, Nyholm T, Garpebring A (2019) Gray-level invariant Haralick texture features.
//' PLoS ONE 14(2): e0212110. \doi{10.1371/journal.pone.0212110}
//' @return a Rcpp::NumericVector of Haralick's texture features
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector cpp_h_features(const Rcpp::NumericMatrix cooc,
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
//' @param cx double, x centroid of the img\cr
//' @param cy double, y centroid of the img\cr
//' @param p uint8_t: p order. Default is 0.
//' @param q uint8_t: q order. Default is 0.
//' @keywords internal
////' @export
// [[Rcpp::export]]
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
//' -5.0 for 40x\cr
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
// [[Rcpp::export]]
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
//' -5.0 for 40x\cr
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
// [[Rcpp::export]]
Rcpp::NumericVector cpp_features_hu2(const Rcpp::NumericMatrix img,
                                     const double mag = 1.0) {
  return hpp_features_hu2(img, mag);
}

//' @title Basic Features
//' @name cpp_basic
//' @description
//' This function is designed to compute very basic features based on Hu's moments + intensities.
//' @param img a NumericMatrix, containing image intensity values.
//' @param mag a double, magnification scale. Default is 1.0. Use:\cr
//' -1.0 for 20x\cr
//' -5.0 for 40x\cr
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
// [[Rcpp::export]]
Rcpp::NumericVector cpp_basic(const Rcpp::NumericMatrix img,
                              const Rcpp::NumericMatrix msk,
                              const double mag = 1.0) {
  return hpp_basic(img, msk, mag);
}


//' @title Image Features Extraction
//' @name cpp_features_hu3
//' @description
//' This function is designed to compute image features.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk an IntegerMatrix, containing msk components.
//' @param components an unsigned integer. Maximal component component number to retrieve features about.
//' Default is 0 to retrieve  features for all components.
//' @param mag a double, magnification scale. Default is 1.0. Use:\cr
//' -1.0 for 20x\cr
//' -5.0 for 40x\cr
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
//' -kurtosis, component's kurtosis.  
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_features_hu3(const Rcpp::NumericMatrix img,
                                     const Rcpp::IntegerMatrix msk,
                                     const unsigned int components = 0,
                                     const double mag = 1.0) {
  return hpp_features_hu3(img, msk, components, mag);
}
// END hu

// FROM otsu
//' @title Otsu Multi Thresholding
//' @name cpp_multi_otsu
//' @description
//' This function determines best threshold(s) according to Otsu's method.
//' @param img, a NumericMatrix.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
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
// [[Rcpp::export]]
Rcpp::NumericVector cpp_multi_otsu (const Rcpp::NumericMatrix img,
                                    const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
//' @param cx a double. X centroid.
//' @param cy a double. Y centroid.
//' @param nmax a uint8_t, maximal order of Zernike polynomials to be computed. Default is 15. Values outside [0,99] will be clipped.
//' Be aware that computation of Zernike's moments can be quite long when 'nmax' is high.
//' @param radius a numeric, radius of the circle in pixels around object centers from which the features are calculated. Default is 15.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List cpp_zernike1(const Rcpp::NumericMatrix img,
                        const double cx,
                        const double cy,
                        const uint8_t nmax = 15,
                        const double radius = 15.0) {
  return hpp_zernike1(img, cx, cy, nmax, radius);
}

//' @title Zernike's Features with Projections
//' @name cpp_zernike2
//' @description
//' This function is designed to compute Zernike's moments from image.
//' It will compute Zernike's moments but also return image projection.
//' @param img a NumericMatrix, containing image intensity values.
//' @param cx a double. X centroid.
//' @param cy a double. Y centroid.
//' @param nmax a uint8_t, maximal order of Zernike polynomials to be computed. Default is 15. Values outside [0,99] will be clipped.
//' Be aware that computation of Zernike's moments can be quite long when 'nmax' is high.
//' @param radius a numeric, radius of the circle in pixels around object centers from which the features are calculated. Default is 15.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List cpp_zernike2(const Rcpp::NumericMatrix img, 
                        const double cx, 
                        const double cy, 
                        const uint8_t nmax = 15, 
                        const double radius = 15.0) {
  return hpp_zernike2(img, cx, cy, nmax, radius);
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_flip(const Rcpp::NumericMatrix mat, const bool which = true) {
  if(which) return hpp_hflip(mat);
  return hpp_hflip(mat);
}
// END flip

// FROM padding
//' @title Image Padding
//' @name cpp_padding
//' @description
//' This function creates a new matrix with extra rows / cols according to input mat, kernel
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method, a uint8_t. Default is 1, allowed are [1-8].\cr
//' -1, extra cols / rows will be filled with 'k', returned 'out' will not be filled.\cr
//' -2, extra cols / rows will be filled with the closest col / row, returned 'out' will not be filled.\cr
//' -3, extra cols / rows will be filled mirroring neighbor cols / rows, returned 'out' will not be filled.\cr
//' -4, extra cols / rows will be filled repeating neighbor cols / rows, returned 'out' will not be filled.\cr
//' -5, extra cols / rows will be filled with 'k', returned 'out' will be filled with mat.\cr
//' -6, extra cols / rows will be filled with the closest col / row, returned 'out' will be filled with mat.\cr
//' -7, extra cols / rows will be filled mirroring neighbor cols / rows, returned 'out' will be filled with mat.\cr
//' -8, extra cols / rows will be filled repeating neighbor cols / rows, returned 'out' will be filled with mat.
//' @param k, a double, constant used when method is 1 or 4. Default is 0.0.
//' @return a List whose elements are:\cr
//' -out, a NumericMatrix, with extra cols / rows\cr
//' -ori_c, a R_len_t with x coordinate of the 1st non extra element,\cr
//' -ori_r, a R_len_t with y coordinate of the 1st non extra element.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List cpp_padding(const Rcpp::NumericMatrix mat,
                       const Rcpp::NumericMatrix kernel,
                       const uint8_t method = 1,
                       const double k = 0.0) {
  return hpp_padding(mat, kernel, method, k);
}
// END padding


// FROM filter
//' @title Image Standard Deviation Filtering
//' @name cpp_sd
//' @description
//' This function applies standard deviation filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_sd(const Rcpp::NumericMatrix mat,
                           const Rcpp::NumericMatrix kernel) {
  return hpp_sd(mat, kernel);
}

//' @title Image Median Filtering
//' @name cpp_median
//' @description
//' This function applies median filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_median(const Rcpp::NumericMatrix mat,
                               const Rcpp::NumericMatrix kernel) {
  return hpp_median(mat, kernel);
}

//' @title Image Mode Filtering
//' @name cpp_mode
//' @description
//' This function applies mode filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_mode(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel) {
  return hpp_mode(mat, kernel);
}

//' @title Image Mid Filtering
//' @name cpp_mid
//' @description
//' This function applies mid filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_mid(const Rcpp::NumericMatrix mat,
                            const Rcpp::NumericMatrix kernel) {
  return hpp_mid(mat, kernel);
}

//' @title Image Filtering by Convolution
//' @name cpp_convolve2d
//' @description
//' This function applies 2D convolution filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_convolve2d(const Rcpp::NumericMatrix mat,
                                   const Rcpp::NumericMatrix kernel) {
  return hpp_convolve2d(mat, kernel);
}

//' @title Image Filtering by Correlation
//' @name cpp_correlate2d
//' @description
//' This function applies 2D correlation filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_correlate2d(const Rcpp::NumericMatrix mat,
                                    const Rcpp::NumericMatrix kernel) {
  return hpp_correlate2d(mat, kernel);
}
// END filter

// FROM morphology
//' @title Image Erosion
//' @name cpp_erode
//' @description
//' This function applies erosion on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_erode(const Rcpp::NumericMatrix mat,
                              const Rcpp::NumericMatrix kernel,
                              const uint8_t iter = 0) {
  return hpp_erode(mat, kernel, iter);
}

//' @title Image Dilatation
//' @name cpp_dilate
//' @description
//' This function applies dilatation on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dilate(const Rcpp::NumericMatrix mat,
                               const Rcpp::NumericMatrix kernel,
                               const uint8_t iter = 0) {
  return hpp_dilate(mat, kernel, iter);
}

//' @title Image Opening
//' @name cpp_opening
//' @description
//' This function applies opening on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_opening(const Rcpp::NumericMatrix mat,
                                const Rcpp::NumericMatrix kernel,
                                const uint8_t iter = 0) {
  return hpp_opening(mat, kernel, iter);
}

//' @title Image Closing
//' @name cpp_closing
//' @description
//' This function applies closing on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_closing(const Rcpp::NumericMatrix mat,
                                const Rcpp::NumericMatrix kernel,
                                const uint8_t iter = 0) {
  return hpp_closing(mat, kernel, iter);
}

//' @title Image Morphological Gradient
//' @name cpp_gradient
//' @description
//' This function applies morphological gradient on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_gradient(const Rcpp::NumericMatrix mat,
                                 const Rcpp::NumericMatrix kernel,
                                 const uint8_t iter = 0) {
  return hpp_gradient(mat, kernel, iter);
}

//' @title Image White Top Hat
//' @name cpp_tophat_white
//' @description
//' This function applies white top hat on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_tophat_white(const Rcpp::NumericMatrix mat,
                                     const Rcpp::NumericMatrix kernel,
                                     const uint8_t iter = 0) {
  return hpp_tophat_white(mat, kernel, iter);
}

//' @title Image Black Top Hat
//' @name cpp_tophat_black
//' @description
//' This function applies black top hat on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_tophat_black(const Rcpp::NumericMatrix mat,
                                     const Rcpp::NumericMatrix kernel,
                                     const uint8_t iter = 0) {
  return hpp_tophat_black(mat, kernel, iter);
}

//' @title Image Self Complementary Top Hat
//' @name cpp_tophat_self
//' @description
//' This function applies self complementary on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_tophat_self(const Rcpp::NumericMatrix mat,
                                    const Rcpp::NumericMatrix kernel,
                                    const uint8_t iter = 0) {
  return hpp_tophat_self(mat, kernel, iter);
}

//' @title Image Contrast Enhancement
//' @name cpp_cont
//' @description
//' This function applies contrast enhancement on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_cont(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel,
                             const uint8_t iter = 0) {
  return hpp_cont(mat, kernel, iter);
}

//' @title Image Laplacian
//' @name cpp_laplacian
//' @description
//' This function applies Laplacian morphology on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate/erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_laplacian(const Rcpp::NumericMatrix mat,
                                  const Rcpp::NumericMatrix kernel,
                                  const uint8_t iter = 0) {
  return hpp_laplacian(mat, kernel, iter);
}
// END morphology

// FROM geodesic
//' @title H-Minima transformation
//' @name cpp_HMIN
//' @description
//' Keep deepest valleys from image.
//' @param img, a NumericMatrix.
//' @param h, a double, specifying the minimal depth. Default is 0.0. it should be positive.
//' @param img_min, a double. Minimal value of image range. Default is 0.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual minimal value of 'img'.
//' @param img_max, a double. Maximal value of image range. Default is 1.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual maximal value of 'img'.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' HMIN is the erosion reconstruction of (img + h) by kernel where img + h is clipped to img_max.
//' @return a NumericMatrix of H-Minima transformation of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_HMIN (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_HMIN(img, h, img_min, img_max, kernel);
}

//' @title H-Maxima transformation
//' @name cpp_HMAX
//' @description
//' Keep highest peaks in image.
//' @param img, a NumericMatrix.
//' @param h, a double, specifying the minimal height. Default is 0.0, it should be positive.
//' @param img_min, a double. Minimal value of image range. Default is 0.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual minimal value of 'img'.
//' @param img_max, a double. Maximal value of image range. Default is 1.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual maximal value of 'img'.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' HMAX is the dilatation reconstruction of (img - h) by kernel where img - h is clipped to img_min.
//' @return a NumericMatrix of H-Maxima transformation of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_HMAX (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_HMAX(img, h, img_min, img_max, kernel);
}

//' @title Regional Minima
//' @name cpp_RMIN
//' @description
//' Mask connected component of pixels whose values are lower to their external boundaries neighborhood.
//' @param img, a NumericMatrix.
//' @param img_min, a double. Minimal value of image range. Default is 0.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual minimal value of 'img'.
//' @param img_max, a double. Maximal value of image range. Default is 1.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual maximal value of 'img'.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' RMIN is defined as HMIN(img, 1)
//' @return a NumericMatrix of regional minima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_RMIN (const Rcpp::NumericMatrix img,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_RMIN(img, img_min, img_max, kernel);
}

//' @title Regional Maxima
//' @name cpp_RMAX
//' @description
//' Mask connected component of pixels whose values are higher to their external boundaries neighborhood.
//' @param img, a NumericMatrix.
//' @param img_min, a double. Minimal value of image range. Default is 0.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual minimal value of 'img'.
//' @param img_max, a double. Maximal value of image range. Default is 1.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual maximal value of 'img'.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' RMAX is defined as HMAX(img, 1)
//' @return a NumericMatrix of regional maxima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_RMAX (const Rcpp::NumericMatrix img,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_RMAX(img, img_min, img_max, kernel);
}

//' @title Extended Minima
//' @name cpp_EMIN
//' @description
//' Mask the regional minima of the corresponding h-minima transformation.
//' @param img, a NumericMatrix.
//' @param img_min, a double. Minimal value of image range. Default is 0.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual minimal value of 'img'.
//' @param img_max, a double. Maximal value of image range. Default is 1.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual maximal value of 'img'.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' EMIN is defined as RMIN(HMIN(img, h))
//' @return a NumericMatrix of extended minima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_EMIN (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_EMIN(img, h, img_min, img_max, kernel);
}

//' @title Extended Maxima
//' @name cpp_EMAX
//' @description
//' Mask the regional maxima of the corresponding h-maxima transformation.
//' @param img, a NumericMatrix.
//' @param img_min, a double. Minimal value of image range. Default is 0.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual minimal value of 'img'.
//' @param img_max, a double. Maximal value of image range. Default is 1.0.
//' It is the value provided by the user, according to what the 'img' range should cover and not the actual maximal value of 'img'.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details see 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}\cr
//' EMAX is defined as RMAX(HMAX(img, h))
//' @return a NumericMatrix of extended maxima of 'img'.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_EMAX (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_EMAX(img, h, img_min, img_max, kernel);
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @param n_lev, an unsigned short determining the number of elevation levels. Default is 256, should be at least 2.
//' @param draw_lines, a bool; whether to draw watershed lines or not. Default is true.
//' @param invert, a bool; whether to fill from basins (lowest values) to peaks (highest values). Default is false.
//' When 'mat' is the result of the distance transformation of an image, peaks (highest values) represent largest distances from background.
//' Thus, they are the ones to be filled first; this is the default behavior with 'invert' set to false.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details adaptation of 'Determining watersheds in digital pictures via flooding simulations' from P. Soille. and L. Vincent.
//' In Proc. SPIE 1360, Visual Communications and Image Processing '90: Fifth in a Series, (1 September 1990) \url{https://doi:10.1117/12.24211}.
//' @source MorphoLib plugin for ImageJ presents a Java implementation of the algorithm in  \url{https://github.com/ijpb/MorphoLibJ/blob/master/src/main/java/inra/ijpb/watershed/WatershedTransform2D.java} authored by Ignacio Arganda-Carreras 
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_watershed_sv1(const Rcpp::NumericMatrix mat,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                                      const unsigned short n_lev = 256,
                                      const bool draw_lines = true,
                                      const bool invert = false,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_watershed_sv1(mat, msk_, n_lev, draw_lines, invert, kernel);
}

//' @title Watershed Transformation SV2
//' @name cpp_watershed_sv2
//' @description
//' This function computes the watershed transformation of an image.
//' @param mat, a NumericMatrix; a distance transform matrix is expected.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @param n_lev, an unsigned short determining the number of elevation levels. Default is 256, should be at least 2.
//' @param draw_lines, a bool; whether to draw watershed lines or not. Default is true.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param invert, a bool; whether to fill from basins (lowest values) to peaks (highest values). Default is false.
//' When 'mat' is the result of the distance transformation of an image, peaks (highest values) represent largest distances from background.
//' Thus, they are the ones to be filled first; this is the default behavior with 'invert' set to false.
//' @details adaptation of 'Watersheds in digital spaces: an efficient algorithm based on immersion simulations' from  L. Vincent and P. Soille.
//' In IEEE Transactions on Pattern Analysis and Machine Intelligence, 13(6):583-598, June 1991.\cr
//' @source The algorithm is reviewed in 'The Watershed Transform: Definitions, Algorithms and Parallelization Strategies'
//' from Roerdink, J. B. T. M. and Meijster, A. (2000) in Fundamenta Informaticae, 41, 187-228 \doi{10.3233/FI-2000-411207}
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_watershed_sv2(const Rcpp::NumericMatrix mat,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                                      const unsigned short n_lev = 256,
                                      const bool draw_lines = true,
                                      const bool invert = false,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  return hpp_watershed_sv2(mat, msk_, n_lev, draw_lines, invert, kernel);
}
// END watershed

// FROM ctl
//' @title Contour Tracing Connected Component Labeling
//' @name cpp_ctl
//' @description
//' This function is designed to identify connected component.
//' @param mat a LogicalMatrix, containing mask.
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
// [[Rcpp::export]]
Rcpp::List cpp_ctl(const Rcpp::LogicalMatrix mat,
                   const bool global = false) {
  return hpp_ctl(mat, global);
}

//' @title Contours Filling
//' @name cpp_fill
//' @description
//' This function is designed to fill contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param label an uint32_t corresponding to the label of desired set of contour to be filled.
//' Default is 0 to fill all set of contours found.
//' @param inner a bool, to whether or not fill hole(s) inside contours if some where identified
//' @param outer a bool, to whether or not fill contours outside hole(s) if some where identified
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_fill(const List ctl,
                             const uint32_t label = 0,
                             const bool inner = true,
                             const bool outer = true) {
  return hpp_fill(ctl, label, inner, outer);
}

//' @title Contours Filling Outer Only
//' @name cpp_fill_out
//' @description
//' This function is designed to fill the most external contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_fill_out(const List ctl) {
  return hpp_fill_out(ctl);
}

//' @title Contours Dilatation
//' @name cpp_dilate_ctl
//' @description
//' This function applies contours dilatation.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dilate_ctl(const List ctl,
                                   const Rcpp::NumericMatrix kernel,
                                   const uint8_t iter = 0) {
  return hpp_dilate_ctl(ctl, kernel, iter);
}

//' @title Contours Erosion
//' @name cpp_erode_ctl
//' @description
//' This function applies contours erosion.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_erode_ctl(const List ctl,
                                  const Rcpp::NumericMatrix kernel,
                                  const uint8_t iter = 0) {
  return hpp_erode_ctl(ctl, kernel, iter);
}
// END ctl

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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
//' It is a measure of the degree to which two images are linearly correlated within a masked region.
//' @return a double, the similarity.
//' @keywords internal
////' @export
// [[Rcpp::export]]
double cpp_similarity(const Rcpp::NumericMatrix img1,
                      const Rcpp::NumericMatrix img2,
                      const Rcpp::LogicalMatrix msk) {
  return hpp_similarity(img1, img2, msk);
}
// END similarity

// FROM kernel
//' @title Create a Disc
//' @name cpp_make_disc
//' @description
//' This function is designed to create a disc kernel.
//' @param size, a uint8_t of the desired kernel size.
//' @return a LogicalMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_convexhull(const Rcpp::NumericMatrix pts) {
  return hpp_convexhull(pts);
}
// END convexhull
