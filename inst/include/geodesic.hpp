/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2022 Yohann Demont                                              

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

#ifndef IFCIP_GEODESIC_HPP
#define IFCIP_GEODESIC_HPP

#include <Rcpp.h>
#include "fifo.hpp"
#include "kernel.hpp"
#include "morphology.hpp"
using namespace Rcpp;

// compute min value in max neighborhood
// M marker image matrix
// n neighbors vector
// u marker value
// v mask value
double nbr_min (const Rcpp::NumericMatrix M,
                const Rcpp::IntegerVector n,
                const double u,
                const double v) {
  double res = u;
  for(uint8_t i = 1; i <= n[0]; i++) if(res < M[n[i]]) res =  M[n[i]];
  return std::min(v, res);
}

// compute max value in min neighborhood
// M marker image matrix
// n neighbors vector
// u marker value
// v mask value
double nbr_max (const Rcpp::NumericMatrix M,
                const Rcpp::IntegerVector n,
                const double u,
                const double v) {
  double res = u;
  for(uint8_t i = 1; i <= n[0]; i++) if(res > M[n[i]]) res =  M[n[i]];
  return std::max(v, res);
}

//' @title Dilatation Reconstruction
//' @name rec_dilate
//' @description
//' Perform a dilatation reconstruction of an image.
//' @param r, a NumericMatrix, containing the dilated values.
//' @param s, a NumericMatrix, specifying the max allowed values.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details adaptation of 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}
//' @return nothing, but r will be modified in-place
//' @keywords internal
void rec_dilate (Rcpp::NumericMatrix r,
                 const Rcpp::NumericMatrix s,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  R_len_t mat_r = r.nrow();
  R_len_t mat_c = r.ncol();
  R_len_t MAX_SIZ = mat_r * mat_c;
  if(MAX_SIZ <= 0) Rcpp::stop("rec_dilate: marker image 'r' should be at least 1px row and 1px col");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 3)) Rcpp::stop("rec_dilate: marker image 'r' is too large");
  if((s.nrow() != mat_r) || (s.ncol() != mat_c)) Rcpp::stop("rec_dilate: marker image 'r' and mask image 's' should have same dimensions");

  Rcpp::IntegerMatrix o = offset_kernel(kernel);               // matrix of offsets in raster order
  Rcpp::IntegerMatrix op = offset_backward(o);                 // before cur_pos, in raster order
  Rcpp::IntegerMatrix on = offset_forward(o);                  // to be scanned after cur_pos is encountered, in raster order
  Rcpp::IntegerVector n(o.ncol() + 1);                         // vector of neighbors idx
  Rcpp::IntegerVector Q = fifo_create(MAX_SIZ + 3, NA_INTEGER);// hierarchical priority queue
  
  // forward scan
  for(R_len_t p = 0; p < r.size(); p++) {
    offset_nbr(p, mat_r, mat_c, op, n);
    r[p] = nbr_min(r, n, r[p], s[p]);
  }
  // backward scan
  for(R_len_t p = r.size() - 1; p >= 0; p--) {
    offset_nbr(p, mat_r, mat_c, on, n);
    r[p] = nbr_min(r, n, r[p], s[p]);
    for(uint8_t i = 1; i <= n[0]; i++) {
      if((r[n[i]] < r[p]) && (r[n[i]] < s[n[i]])) {
        fifo_add(Q, p);
        break;
      }
    }
  }
  // propagate
  while(Q[0] != 0) {
    R_len_t p = fifo_pop(Q);
    offset_nbr(p, mat_r, mat_c, o, n);
    for(uint8_t i = 1; i <= n[0]; i++) {
      if((r[n[i]] < r[p]) && (s[n[i]] != r[n[i]])) {
        r[n[i]] = std::min(r[p], s[n[i]]);
        fifo_add(Q, n[i]);
      }
    }
  }
}

//' @title Erosion Reconstruction
//' @name rec_erode
//' @description
//' Perform an erosion reconstruction of an image.
//' @param r, a NumericMatrix, marker image containing the eroded values.
//' @param s, a NumericMatrix, mask image specifying the min allowed values.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @details adaptation of 'Morphological grayscale reconstruction in image analysis: applications and efficient algorithms' from  L. Vincent.
//' IEEE Transactions on Image Processing, 2(2):176-201, April 1993.\doi{10.1109/83.217222}
//' @return nothing, but r will be modified in-place
//' @keywords internal
void rec_erode (Rcpp::NumericMatrix r,
                const Rcpp::NumericMatrix s,
                const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  R_len_t mat_r = r.nrow();
  R_len_t mat_c = r.ncol();
  R_len_t MAX_SIZ = mat_r * mat_c;
  if(MAX_SIZ <= 0) Rcpp::stop("rec_erode: marker image 'r' should be at least 1px row and 1px col");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 3)) Rcpp::stop("rec_erode: marker image 'r' is too large");
  if((s.nrow() != mat_r) || (s.ncol() != mat_c)) Rcpp::stop("rec_erode: marker image 'r' and mask image 's' should have same dimensions");
  
  Rcpp::IntegerMatrix o = offset_kernel(kernel);               // matrix of offsets in raster order
  Rcpp::IntegerMatrix op = offset_backward(o);                 // before cur_pos, in raster order
  Rcpp::IntegerMatrix on = offset_forward(o);                  // to be scanned after cur_pos is encountered, in raster order
  Rcpp::IntegerVector n(o.ncol() + 1);                         // vector of neighbors idx
  Rcpp::IntegerVector Q = fifo_create(MAX_SIZ + 3, NA_INTEGER);// hierarchical priority queue
  
  // forward scan
  for(R_len_t p = 0; p < r.size(); p++) {
    offset_nbr(p, mat_r, mat_c, op, n);
    r[p] = nbr_max(r, n, r[p], s[p]);
  }
  // backward scan
  for(R_len_t p = r.size() - 1; p >= 0; p--) {
    offset_nbr(p, mat_r, mat_c, on, n);
    r[p] = nbr_max(r, n, r[p], s[p]);
    for(uint8_t i = 1; i <= n[0]; i++) {
      if((r[n[i]] > r[p]) && (r[n[i]] > s[n[i]])) {
        fifo_add(Q, p);
        break;
      }
    }
  }
  // propagate
  while(Q[0] != 0) {
    R_len_t p = fifo_pop(Q);
    offset_nbr(p, mat_r, mat_c, o, n);
    for(uint8_t i = 1; i <= n[0]; i++) {
      if((r[n[i]] > r[p]) && (s[n[i]] != r[n[i]])) {
        r[n[i]] = std::max(r[p], s[n[i]]);
        fifo_add(Q, n[i]);
      }
    }
  }
}

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
Rcpp::NumericMatrix hpp_HMIN (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  if(h < 0) Rcpp::stop("HMIN: 'h' should be positive");
  Rcpp::NumericMatrix out = no_init(img.nrow(), img.ncol());
  for(R_len_t i = 0; i < img.size(); i++) {
    out[i] = std::min(img[i] + h, img_max);
  }
  rec_erode(out, img, kernel);
  return out;
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
Rcpp::NumericMatrix hpp_HMAX (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  if(h < 0) Rcpp::stop("HMAX: 'h' should be positive");
  Rcpp::NumericMatrix out = no_init(img.nrow(), img.ncol());
  for(R_len_t i = 0; i < img.size(); i++) {
    out[i] = std::max(img[i] - h, img_min);
  }
  rec_dilate(out, img, kernel);
  return out;
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
Rcpp::NumericMatrix hpp_RMIN (const Rcpp::NumericMatrix img,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix out = hpp_HMIN(img, (img_max - img_min) / (pow(2.0, 32.0) - 1), img_min, img_max, kernel);
  for(R_len_t i = 0; i < img.size(); i++) out[i] = out[i] - img[i];
  return out;
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
Rcpp::NumericMatrix hpp_RMAX (const Rcpp::NumericMatrix img,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix out = hpp_HMAX(img, (img_max - img_min) / (pow(2.0, 32.0) - 1), img_min, img_max, kernel);
  for(R_len_t i = 0; i < img.size(); i++) out[i] = img[i] - out[i];
  return out;
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
Rcpp::NumericMatrix hpp_EMIN (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix out = hpp_RMIN(hpp_HMIN(img, h, img_min, img_max, kernel), img_min, img_max, kernel);
  return out;
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
Rcpp::NumericMatrix hpp_EMAX (const Rcpp::NumericMatrix img,
                              const double h = 0.0,
                              const double img_min = 0.0,
                              const double img_max = 1.0,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix out = hpp_RMAX(hpp_HMAX(img, h, img_min, img_max, kernel), img_min, img_max, kernel);
  return out;
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
Rcpp::NumericMatrix hpp_geo_tophat_white (const Rcpp::NumericMatrix img,
                                          const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix kk = get_kernel(kernel);
  Rcpp::NumericMatrix out = hpp_erode(img, kk);
  rec_dilate(out, img, kk);
  for(R_len_t i = 0; i < out.size(); i++) out[i] = img[i] - out[i];
  return out;
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
Rcpp::NumericMatrix hpp_geo_tophat_black (const Rcpp::NumericMatrix img,
                                          const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix kk = get_kernel(kernel);
  Rcpp::NumericMatrix out = hpp_dilate(img, kk);
  rec_erode(out, img, kk);
  for(R_len_t i = 0; i < out.size(); i++) out[i] = out[i] - img[i];
  return out;
}

#endif
