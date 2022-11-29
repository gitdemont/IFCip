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

#ifndef IFCIP_HU_HPP
#define IFCIP_HU_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Hu's Centroid
//' @name cpp_centroid
//' @description
//' This function is designed to compute Hu's image centroids.
//' @param img a NumericMatrix, containing image intensity values.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_centroid (const Rcpp::NumericMatrix img) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();
  R_len_t i_row, i_col, i_col_1;
  
  // determine centroid
  double area = 0.0, cx = 0.0, cy = 0.0;
  for(i_col = 0; i_col < mat_c; i_col++) {
    i_col_1 = i_col + 1;
    for(i_row = 0; i_row < mat_r; i_row++) {
      area += img(i_row, i_col);
      cy += i_col_1 * img(i_row, i_col) ;
      cx += (i_row + 1) * img(i_row, i_col);
    }
  }
  return(Rcpp::NumericVector::create(_["pix cx"] = cx / area, _["pix cy"] = cy / area));
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
Rcpp::NumericVector hpp_rmoment(const Rcpp::NumericMatrix img, 
                          const uint8_t p = 0, 
                          const uint8_t q = 0) {
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  R_len_t i_row, i_col, i_col_1;
  double rmoment = 0.0;
  Rcpp::NumericMatrix out(mat_r, mat_c) ;
  out.fill(0.0);
  for(i_col = 0; i_col < mat_c; i_col++) {
    i_col_1 = i_col + 1;
    for(i_row = 0; i_row < mat_r; i_row++) {
      out(i_row, i_col) = std::pow(i_row + 1, p) * std::pow(i_col_1, q) * img(i_row, i_col) + out(i_row, i_col);
      rmoment += out(i_row, i_col);
    }
  }
  Rcpp::NumericVector foo = Rcpp::NumericVector::create(_["rmoment"] = rmoment);
  foo.attr("order") = Rcpp::NumericVector::create(_["p"] = p, _["q"] = q);
  foo.attr("matrix") = out;
  return foo;
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
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_cmoment(const Rcpp::NumericMatrix img, 
                          const double cx = 0.0, 
                          const double cy = 0.0, 
                          const uint8_t p = 0, 
                          const uint8_t q = 0) {
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  R_len_t i_row, i_col, i_col_1;
  double cmoment = 0.0;
  Rcpp::NumericMatrix out(mat_r, mat_c) ;
  out.fill(0.0);
  for(i_col = 0; i_col < mat_c; i_col++) {
    i_col_1 = i_col + 1;
    for(i_row = 0; i_row < mat_r; i_row++) {
      out(i_row, i_col) = std::pow(i_row + 1 - cx, p) * std::pow(i_col_1 - cy, q) * img(i_row, i_col) + out(i_row, i_col);
      cmoment += out(i_row, i_col);
    }
  }
  Rcpp::NumericVector foo = Rcpp::NumericVector::create(_["cmoment"] = cmoment);
  foo.attr("order") = Rcpp::NumericVector::create(_["p"] = p, _["q"] = q);
  foo.attr("matrix") = out;
  return foo;
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
//' -Aspect Ratio: img's ratio of minor_axis over major_axis\cr
//' -Angle: img's ellipsis angle with x axis (in radians)\cr
//' -theta: img's ellipsis theta angle (in radians)\cr
//' -eccentricity: img's ellipsis ecentricity\cr
//' -pix cx: img's pixel x centroïd\cr
//' -pix cy: img's pixel y centroïd\cr
//' -pix min_axis: img's ellipsis minor axis in pixels\cr
//' -pix maj axis: img's ellipsis major axis in pixels\cr
//' -pix count: img's area in pixels.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_features_hu1(const Rcpp::NumericMatrix img,
                              const double mag = 1.0) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();
  
  // define raw moments
  double cx = 0.0, cy = 0.0, U00 = 0.0, U02 = 0.0, U11 = 0.0, U20 = 0.0;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t i_col_1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_row_1 = i_row + 1;
      U00 += img(i_row, i_col);
      cy += i_col_1 * img(i_row, i_col);
      U02 += i_col_1 * i_col_1 * img(i_row, i_col);
      cx += i_row_1 * img(i_row, i_col);
      U11 += i_row_1 * i_col_1 * img(i_row, i_col);
      U20 += i_row_1 * i_row_1 * img(i_row, i_col);
    }
  }
  
  // size invariant
  cx /= U00;
  cy /= U00;
  
  // from https://doi.org/10.1016/j.patcog.2009.06.017
  // A Hu moment invariant as a shape circularity measure
  // Joviša Žunić, Kaoru Hirota, Paul L. Rosin
  double circularity = U00 * U00 / (2 * M_PI * ( U20 + U02 ) );
  
  // covariance
  double Up02 = U02 / U00 - cy * cy;
  double Up11 = U11 / U00 - cx * cy;
  double Up20 = U20 / U00 - cx * cx;
  // if(Up20 == Up02) {
  //   Rcpp::stop("Bad image segmentation");
  // }
  // ellipse
  double d = Up20 - Up02;
  double s = Up20 + Up02;
  double det = std::sqrt(4 * Up11 * Up11 + d * d);
  double pix_maj_axis = std::sqrt( (s + det) / 2) * 4;
  double pix_min_axis = std::sqrt( (s - det) / 2) * 4;
  
  // real size
  double Up02r = Up02 / mag ;
  double Up11r = Up11 / mag ;
  double Up20r = Up20 / mag ;
  
  // real size ellipse
  double area = U00 / mag;
  double dr = Up20r - Up02r;
  double sr = Up20r + Up02r;
  double detr = std::sqrt(4 * Up11r * Up11r + dr * dr);
  double maj_axis = std::sqrt( (sr + detr) / 2) * 4;
  double min_axis = std::sqrt( (sr - detr) / 2) * 4;
  double eccentricity = std::sqrt( 1 - (min_axis * min_axis) / (maj_axis * maj_axis) );
  double theta = std::atan2(2*Up11r, dr) / 2;
  double angle = theta < 0.0 ? -theta:theta;
  angle += M_PI /2;
  return Rcpp::NumericVector::create(_["Area"] = area,
                               _["circularity"] = circularity,
                               _["Minor Axis"] = min_axis,
                               _["Major Axis"] = maj_axis,
                               _["Aspect Ratio"] = min_axis/maj_axis,
                               _["Angle"] = angle,
                               _["theta"] = theta,
                               _["eccentricity"] = eccentricity,
                               _["pix cx"] = cx,
                               _["pix cy"] = cy,
                               _["pix min axis"] = pix_min_axis,
                               _["pix maj axis"] = pix_maj_axis,
                               _["pix count"] = U00);
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
//' -Aspect Ratio: img's ratio of minor_axis over major_axis\cr
//' -Angle: img's ellipsis angle with x axis (in radians)\cr
//' -theta: img's ellipsis theta angle (in radians)\cr
//' -eccentricity: img's ellipsis ecentricity\cr
//' -pix cx: img's pixel x centroïd\cr
//' -pix cy: img's pixel y centroïd\cr
//' -pix min_axis: img's ellipsis minor axis in pixels\cr
//' -pix maj axis: img's ellipsis major axis in pixels\cr
//' -pix count: img's area in pixels\cr
//' -inv[1-7]: image invariant moments.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_features_hu2(const Rcpp::NumericMatrix img,
                               const double mag = 1.0) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();

  // define raw moments
  double m00 = 0.0, m01 = 0.0, m02 = 0.0, m03 = 0.0, m10 = 0.0, m11 = 0.0, m12 = 0.0, m20 = 0.0, m21 = 0.0, m30 = 0.0;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t i_col_1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_row_1 = i_row + 1;
      // for centroid and ellipse
      m00 += img(i_row, i_col) ;
      m01 += i_col_1 * img(i_row, i_col) ;
      m10 += i_row_1 * img(i_row, i_col);
      m11 += i_row_1 * i_col_1 * img(i_row, i_col);
      m02 += i_col_1 * i_col_1 * img(i_row, i_col);
      m20 += i_row_1 * i_row_1 * img(i_row, i_col);
      // for invariants
      m03 += i_col_1 * i_col_1 * i_col_1 * img(i_row, i_col);
      m12 += i_row_1 * i_col_1 * i_col_1 * img(i_row, i_col);
      m21 += i_row_1 * i_row_1 * i_col_1 * img(i_row, i_col);
      m30 += i_row_1 * i_row_1 * i_row_1 * img(i_row, i_col);
    }
  }
  
  // size invariant
  double cx = m10 / m00;
  double cy = m01 / m00;

  // from https://doi.org/10.1016/j.patcog.2009.06.017
  // A Hu moment invariant as a shape circularity measure
  // Joviša Žunić, Kaoru Hirota, Paul L. Rosin
  double circularity = m00 * m00 / (2 * M_PI * ( m20 + m02 ) );
  
  // define covariance
  double Up02 = (m02 / m00 - cy * cy);
  double Up11 = (m11 / m00 - cx * cy);
  double Up20 = (m20 / m00 - cx * cx);
  // if(Up20 == Up02) {
  //   Rcpp::stop("Bad image segmentation");
  // }
  
  // ellipse
  double d = Up20 - Up02;
  double s = Up20 + Up02;
  double det = std::sqrt(4 * Up11 * Up11 + d * d);
  double pix_maj_axis = std::sqrt( (s + det) / 2) * 4;
  double pix_min_axis = std::sqrt( (s - det) / 2) * 4;
  
  // real size
  double Up02r = Up02 / mag ;
  double Up11r = Up11 / mag ;
  double Up20r = Up20 / mag ;
  
  // real size ellipse
  double area = m00 / mag;
  double dr = Up20r - Up02r;
  double sr = Up20r + Up02r;
  double detr = std::sqrt(4 * Up11r * Up11r + dr * dr);
  double maj_axis = std::sqrt( (sr + detr) / 2) * 4;
  double min_axis = std::sqrt( (sr - detr) / 2) * 4;
  double eccentricity = std::sqrt( 1 - (min_axis * min_axis) / (maj_axis * maj_axis) );
  double theta = std::atan2(2 * Up11r, dr) / 2;
  double angle = theta < 0.0 ? -theta:theta;
  angle += M_PI /2;
  
  // define central moments of order 3
  double U00 = m00, U11, U20, U02, U21, U12, U30, U03;
  U11 = m11 - cx * m01;
  U20 = m20 - cx * m10;
  U02 = m02 - cy * m01;
  U21 = m21 - 2 * cx * m11 - cy * m20 + 2 * cx * cx * m01;
  U12 = m12 - 2 * cy * m11 - cx * m02 + 2 * cy * cy * m10;
  U30 = m30 - 3 * cx * m20 + 2 * cx * cx * m10;
  U03 = m03 - 3 * cy * m02 + 2 * cy * cy * m01;
  
  // define invariants
  // scale invariants
  double n02, n03, n11, n12, n20, n21, n30;
  n02 = U02 / U00 / U00; // std::pow(U00, 2.0); // 1 + (0+2)/2 
  n03 = U03 / std::pow(U00, 2.5); // 1 + (0+3)/2 
  n11 = U11 / U00 / U00; // std::pow(U00, 2.0); // 1 + (1+1)/2 
  n12 = U12 / std::pow(U00, 2.5); // 1 + (1+2)/2
  n20 = U20 / U00 / U00; // std::pow(U00, 2.0); // 1 + (2+0)/2
  n21 = U21 / std::pow(U00, 2.5); // 1 + (2+1)/2
  n30 = U30 / U00 / U00; // std::pow(U00, 2.5); // 1 + (3+0)/2 
  
  // rotation invariants
  double I1, I2, I3, I4, I5, I6, I7;
  I1 = n20 + n02;
  I2 = (n20 - n02) * (n20 - n02) + 4 * n11 * n11;
  I3 = (n30 - 3 * n12) * (n30 - 3 * n12) + (3 * n21 - n03) * (3 * n21 - n03);
  I4 = (n30 + n12) * (n30 + n12) + (n21 * n03) * (n21 * n03);
  I5 = (n30 - 3 * n12) * (n30 + n12) * ( (n30 + n12) * (n30 + n12) - 3 * (n21 + n03) * (n21 + n03) )
    + (3 * n21 - n03) * (n21 + n03) * ( 3 * (n30 + n12) * (n30 + n12) - (n21 + n03) * (n21 + n03) );
  I6 = (n20 - n02) * ( (n30 + n12) * (n30 + n12) - (n21 + n03) * (n21 + n03) )
    + 4 * n11 * (n30 + n12) * (n21 + n03);
  I7 = (3 * n21 - n03) * (n30 + n12) *( (n30 + n12) * (n30 + n12) - 3 * (n21 + n03) * (n21 + n03) )
    - (n30 - 3 * n12) * (n21 + n03) * ( 3 * (n30 + n12) * (n30 + n12) - (n21 + n03) * (n21 + n03) );
  
  // return result
  return Rcpp::NumericVector::create(_["Area"] = area,
                               _["circularity"] = circularity,
                               _["Minor Axis"] = min_axis,
                               _["Major Axis"] = maj_axis,
                               _["Aspect Ratio"] = min_axis/maj_axis,
                               _["Angle"] = angle,
                               _["theta"] = theta,
                               _["eccentricity"] = eccentricity,
                               _["pix cx"] = cx,
                               _["pix cy"] = cy,
                               _["pix min axis"] = pix_min_axis,
                               _["pix maj axis"] = pix_maj_axis,
                               _["pix count"] = U00,
                               _["inv1"] = I1,
                               _["inv2"] = I2,
                               _["inv3"] = I3,
                               _["inv4"] = I4,
                               _["inv5"] = I5,
                               _["inv6"] = I6,
                               _["inv7"] = I7);
}

//' @title Basic Features
//' @name cpp_basic
//' @description
//' This function is designed to compute very basic features based on Hu's moments + intensities.
//' @param img a NumericMatrix, containing image intensity values.
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
Rcpp::NumericVector hpp_basic(const Rcpp::NumericMatrix img,
                              const Rcpp::NumericMatrix msk,
                              const double mag = 1.0) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();
  if(msk.ncol() != mat_c || msk.nrow() != mat_r) Rcpp::stop("hpp_basic: please extract 'img' with raw size");
  
  // define raw moments
  double cx = 0.0, cy = 0.0, U00 = 0.0, U02 = 0.0, U11 = 0.0, U20 = 0.0;
  double imin = R_PosInf, imax = R_NegInf, imean = 0.0;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t i_col_1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      if(msk(i_row, i_col)) continue; 
      // for msk moments
      R_len_t i_row_1 = i_row + 1;
      U00 += 1;
      cy += i_col_1;
      U02 += i_col_1 * i_col_1;
      cx += i_row_1;
      U11 += i_row_1 * i_col_1;
      U20 += i_row_1 * i_row_1;
      // for intenties
      imean += img(i_row, i_col);
      if(imin > img(i_row, i_col)) imin = img(i_row, i_col);
      if(imax < img(i_row, i_col)) imax = img(i_row, i_col);
    }
  }
  
  // size invariant
  cx /= U00;
  cy /= U00;
  
  // from https://doi.org/10.1016/j.patcog.2009.06.017
  // A Hu moment invariant as a shape circularity measure
  // Joviša Žunić, Kaoru Hirota, Paul L. Rosin
  double circularity = U00 * U00 / (2 * M_PI * ( U20 + U02 ) );
  
  // covariance
  double Up02 = U02 / U00 - cy * cy;
  double Up11 = U11 / U00 - cx * cy;
  double Up20 = U20 / U00 - cx * cx;
  // if(Up20 == Up02) {
  //   Rcpp::stop("Bad image segmentation");
  // }
  // ellipse
  double d = Up20 - Up02;
  double s = Up20 + Up02;
  double det = std::sqrt(4 * Up11 * Up11 + d * d);
  double pix_maj_axis = std::sqrt( (s + det) / 2) * 4;
  double pix_min_axis = std::sqrt( (s - det) / 2) * 4;
  
  // real size
  double Up02r = Up02 / mag ;
  double Up11r = Up11 / mag ;
  double Up20r = Up20 / mag ;
  
  // real size ellipse
  double area = U00 / mag;
  double dr = Up20r - Up02r;
  double sr = Up20r + Up02r;
  double detr = std::sqrt(4 * Up11r * Up11r + dr * dr);
  double maj_axis = std::sqrt( (sr + detr) / 2) * 4;
  double min_axis = std::sqrt( (sr - detr) / 2) * 4;
  double eccentricity = std::sqrt( 1 - (min_axis * min_axis) / (maj_axis * maj_axis) );
  double theta = std::atan2(2*Up11r, dr) / 2;
  double angle = theta < 0.0 ? -theta:theta;
  angle += M_PI /2;
  return Rcpp::NumericVector::create(_["Area"] = area,
                                     _["circularity"] = circularity,
                                     _["Minor Axis"] = min_axis,
                                     _["Major Axis"] = maj_axis,
                                     _["Aspect Ratio"] = min_axis/maj_axis,
                                     _["Angle"] = angle,
                                     _["theta"] = theta,
                                     _["eccentricity"] = eccentricity,
                                     _["pix cx"] = cx,
                                     _["pix cy"] = cy,
                                     _["pix min axis"] = pix_min_axis,
                                     _["pix maj axis"] = pix_maj_axis,
                                     _["pix count"] = U00,
                                     _["Raw Mean Pixel"] = imean/U00,
                                     _["Raw Min Pixel"] = imin,
                                     _["Raw Max Pixel"] = imax);
}

//' @title Image Features Extraction
//' @name cpp_features_hu3
//' @description
//' This function is designed to compute image features.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk an IntegerMatrix, containing msk components.
//' @param components an unsigned integer. Maximal component number to retrieve features about.
//' Default is 0 to retrieve features for all components.
//' @param mag a double, magnification scale. Default is 1.0. Use:\cr
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
//' -inv[1-7], component Hu's invariant moments\cr
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
Rcpp::NumericMatrix hpp_features_hu3(const Rcpp::NumericMatrix img,
                                     const Rcpp::IntegerMatrix msk,
                                     const unsigned int components = 0,
                                     const double mag = 1.0) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();
  if(mat_r != msk.nrow() || mat_c != msk.ncol()) {
    Rcpp::stop("hpp_features_hu3: 'img' and 'msk' should have same dimensions");
  }
  
  R_len_t nC = components; // == 0 ? 1:components;
  R_len_t alw = (std::ceil(mat_r / 2) + 1) * (std::ceil(mat_c / 2) + 1);
  
  if(nC == 0) {
    for(R_len_t i_msk = 0; i_msk < mat_c * mat_r; i_msk++) {
      if(msk[i_msk] < 0) Rcpp::stop("hpp_features_hu3: invalid negative value for 'msk'");
      if(msk[i_msk] >= alw) Rcpp::stop("hpp_features_hu3: invalid max number of components for 'msk'");
      nC = std::max(nC, msk[i_msk]); 
    }
  } else {
    for(R_len_t i_msk = 0; i_msk < mat_c * mat_r; i_msk++) {
      if(msk[i_msk] < 0) Rcpp::stop("hpp_features_hu3: invalid negative value for 'msk'");
      if(msk[i_msk] >= alw) Rcpp::stop("hpp_features_hu3: invalid max number of components for 'msk'");
    }
  }
  
  // 32 values to return [00-31] 
  // + 16 slots [40-55] for temporary computation
  // +  8 extra slots for future computation ... [32-39]
  
  // initialize matrix
  Rcpp::NumericMatrix foo(nC, 56);
  
  // // create colnames
  Rcpp::StringVector N = Rcpp::StringVector({
    "Area",
    "circularity",
    "Minor Axis","Major Axis","Aspect Ratio",
    "Angle","theta","eccentricity",
    "Minor Axis Intensity","Major Axis Intensity","Aspect Ratio Intensity",
    "Angle Intensity","theta intensity","eccentricity intensity",
    "pix cx","pix cy","pix min axis","pix maj axis",
    "pix count",
    "inv1","inv2","inv3","inv4","inv5","inv6","inv7",
    "Raw Mean Pixel",
    "Raw Min Pixel","Raw Max Pixel",
    "Std Dev","skewness","kurtosis",
    "Centroid Y","Centroid X","Centroid Y Intensity","Centroid X Intensity"});
  
  if(nC == 0) {
    Rcpp::NumericMatrix out(nC, 32);
    Rcpp::colnames(out) = N;
    return out;
  } 
  
  // fill min, max
  for(R_len_t i_comp = 0; i_comp < nC; i_comp++) {
    foo(i_comp, 27) = R_PosInf;
    foo(i_comp, 28) = R_NegInf;
  }
  
  // define raw moments
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    if(Rcpp::is_true(Rcpp::any(Rcpp::is_na(msk(Rcpp::_,i_col))))) Rcpp::stop("hpp_features_hu3: NA - NaN value are not allowed in 'msk'");
    R_len_t i_col_1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_row_1 = i_row + 1;
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp >= 0) && (i_comp < nC)) {
        // for centroid and ellipse
        // M00
        foo(i_comp, 40)++;
        // M01
        foo(i_comp, 41) += i_col_1;
        // M10
        foo(i_comp, 42) += i_row_1;
        // M11
        foo(i_comp, 43) += i_row_1 * i_col_1;
        // M02
        foo(i_comp, 44) += i_col_1 * i_col_1;
        // M20
        foo(i_comp, 45) += i_row_1 * i_row_1;
        
        // for invariants
        // M03
        foo(i_comp, 46) += i_col_1 * i_col_1 * i_col_1 * img(i_row, i_col);    
        // M12
        foo(i_comp, 47) += i_row_1 * i_col_1 * i_col_1 * img(i_row, i_col);
        // M21
        foo(i_comp, 48) += i_row_1 * i_row_1 * i_col_1 * img(i_row, i_col);
        // M30
        foo(i_comp, 49) += i_row_1 * i_row_1 * i_row_1 * img(i_row, i_col);
        // M00
        foo(i_comp, 50) += img(i_row, i_col);
        // M01
        foo(i_comp, 51) += i_col_1 * img(i_row, i_col);
        // M10
        foo(i_comp, 52) += i_row_1 * img(i_row, i_col);
        // M11
        foo(i_comp, 53) += i_row_1 * i_col_1 * img(i_row, i_col);
        // M02
        foo(i_comp, 54) += i_col_1 * i_col_1 * img(i_row, i_col);
        // M20
        foo(i_comp, 55) += i_row_1 * i_row_1 * img(i_row, i_col);
        
        // for intensity
        // Mean
        foo(i_comp, 26) += img(i_row, i_col);
        // Min
        if(foo(i_comp, 27) > img(i_row, i_col)) foo(i_comp, 27) = img(i_row, i_col);
        // Max
        if(foo(i_comp, 28) < img(i_row, i_col)) foo(i_comp, 28) = img(i_row, i_col);
        // pix_count
        foo(i_comp, 18)++;
      }
    }
  }
  
  for(R_len_t i_comp = 0; i_comp < nC; i_comp++) {
    // mean
    foo(i_comp, 26) /= foo(i_comp, 18);
    
    // size invariant
    double cx = foo(i_comp, 42) / foo(i_comp, 40);
    double cy = foo(i_comp, 41) / foo(i_comp, 40);
    double cx_int = foo(i_comp, 52) / foo(i_comp, 50);
    double cy_int = foo(i_comp, 51) / foo(i_comp, 50);
    
    // from https://doi.org/10.1016/j.patcog.2009.06.017
    // A Hu moment invariant as a shape circularity measure
    // Joviša Žunić, Kaoru Hirota, Paul L. Rosin
    double circularity = foo(i_comp, 40) * foo(i_comp, 40) / (2 * M_PI * ( foo(i_comp, 45) + foo(i_comp, 44) ) );
    
    // define covariance
    double Up02 = (foo(i_comp, 44) / foo(i_comp, 40) - cy * cy);
    double Up11 = (foo(i_comp, 43) / foo(i_comp, 40) - cx * cy);
    double Up20 = (foo(i_comp, 45) / foo(i_comp, 40) - cx * cx);
    
    // define covariance intensity weighted
    double Up02i = (foo(i_comp, 54) / foo(i_comp, 50) - cy_int * cy_int);
    double Up11i = (foo(i_comp, 53) / foo(i_comp, 50) - cx_int * cy_int);
    double Up20i = (foo(i_comp, 55) / foo(i_comp, 50) - cx_int * cx_int);
    
    // if(Up20 == Up02) {
    //   circle or square
    //   Rcpp::stop("Bad image segmentation");
    // }
    
    // ellipse
    double d = Up20 - Up02;
    double s = Up20 + Up02;
    double det = std::sqrt(4 * Up11 * Up11 + d * d);
    double pix_maj_axis = std::sqrt( (s + det) / 2) * 4;
    double pix_min_axis = std::sqrt( (s - det) / 2) * 4;
    
    // real size
    double Up02r = Up02 / mag ;
    double Up11r = Up11 / mag ;
    double Up20r = Up20 / mag ;
    double Up02ri = Up02i / mag ;
    double Up11ri = Up11i / mag ;
    double Up20ri = Up20i / mag ;
    double area = foo(i_comp, 40) / mag;
    
    // real size ellipse
    double dr = Up20r - Up02r;
    double sr = Up20r + Up02r;
    double detr = std::sqrt(4 * Up11r * Up11r + dr * dr);
    double maj_axis = std::sqrt( (sr + detr) / 2) * 4;
    double min_axis = std::sqrt( (sr - detr) / 2) * 4;
    double eccentricity = std::sqrt( 1 - (min_axis * min_axis) / (maj_axis * maj_axis) );
    double theta = std::atan2(2 * Up11r, dr) / 2;
    double angle = theta < 0.0 ? -theta:theta;
    angle += M_PI /2;
    
    // intensity weighted real size ellipse
    double dri = Up20ri - Up02ri;
    double sri = Up20ri + Up02ri;
    double detri = std::sqrt(4 * Up11ri * Up11ri + dri * dri);
    double maj_axisi = std::sqrt( (sri + detri) / 2) * 4;
    double min_axisi = std::sqrt( (sri - detri) / 2) * 4;
    double eccentricityi = std::sqrt( 1 - (min_axisi * min_axisi) / (maj_axisi * maj_axisi) );
    double thetai = std::atan2(2 * Up11ri, dri) / 2;
    double anglei = thetai < 0.0 ? -thetai:thetai;
    anglei += M_PI /2;
    
    // define central moments of order 3
    double U11, U20, U02, U21, U12, U30, U03;
    U11 = foo(i_comp, 53) - cx_int * foo(i_comp, 51);
    U20 = foo(i_comp, 55) - cx_int * foo(i_comp, 52);
    U02 = foo(i_comp, 54) - cy_int * foo(i_comp, 51);
    U21 = foo(i_comp, 48) - 2 * cx_int * foo(i_comp, 53) - cy_int * foo(i_comp, 55) + 2 * cx_int * cx_int * foo(i_comp, 51);
    U12 = foo(i_comp, 47) - 2 * cy_int * foo(i_comp, 53) - cx_int * foo(i_comp, 44) + 2 * cy_int * cy_int * foo(i_comp, 52);
    U30 = foo(i_comp, 49) - 3 * cx_int * foo(i_comp, 55) + 2 * cx_int * cx_int * foo(i_comp, 52);
    U03 = foo(i_comp, 46) - 3 * cy_int * foo(i_comp, 54) + 2 * cy_int * cy_int * foo(i_comp, 51);
    
    // define invariants
    // scale invariants
    double n02, n03, n11, n12, n20, n21, n30;
    double U00_2 = foo(i_comp, 40) * foo(i_comp, 40), U00_25 = std::pow(foo(i_comp, 40), 2.5);
    n02 = U02 / U00_2;  // 1 + (0+2)/2 
    n03 = U03 / U00_25; // 1 + (0+3)/2 
    n11 = U11 / U00_2;  // 1 + (1+1)/2 
    n12 = U12 / U00_25; // 1 + (1+2)/2
    n20 = U20 / U00_2;  // 1 + (2+0)/2
    n21 = U21 / U00_25; // 1 + (2+1)/2
    n30 = U30 / U00_25; // 1 + (3+0)/2 
    
    // temp variables
    double foo1, foo2, foo3, foo4, foo5, foo6, foo7;
    foo1 = n30 + n12;
    foo2 = foo1 * foo1;
    foo3 = n21 + n03;
    foo4 = foo3 * foo3;
    foo5 = 3 * n21 - n03;
    foo6 = n30 - 3 * n12;
    foo7 = n20 - n02;
    
    // rotation invariants
    double I1, I2, I3, I4, I5, I6, I7;
    I1 = n20 + n02;
    I2 = foo7 * foo7 + 4 * n11 * n11;
    I3 = foo6 * foo6 + foo5 * foo5;
    I4 = foo2 + foo4;
    I5 = (foo6) * (foo1) * ( foo2 - 3 * foo4 )
      + (foo5) * (foo3) * ( 3 * foo2 - foo4 );
    I6 = foo7 * ( foo2 - foo4 )
      + 4 * n11 * (foo1) * (foo3);
    I7 = (foo5) * (foo1) * ( foo2 - 3 * foo4 )
      - (foo6) * (foo3) * ( 3 * foo2 - foo4 );
    
    foo(i_comp,  0) = area;
    foo(i_comp,  1) = circularity;
    foo(i_comp,  2) = min_axis;
    foo(i_comp,  3) = maj_axis;
    foo(i_comp,  4) = min_axis / maj_axis;
    foo(i_comp,  5) = angle;
    foo(i_comp,  6) = theta;
    foo(i_comp,  7) = eccentricity;
    foo(i_comp,  8) = min_axisi;
    foo(i_comp,  9) = maj_axisi;
    foo(i_comp, 10) = min_axisi / maj_axisi;
    foo(i_comp, 11) = anglei;
    foo(i_comp, 12) = thetai;
    foo(i_comp, 13) = eccentricityi;
    foo(i_comp, 14) = cx;
    foo(i_comp, 15) = cy;
    foo(i_comp, 16) = pix_min_axis;
    foo(i_comp, 17) = pix_maj_axis;
    // foo(i_comp, 18) is pix_count;
    foo(i_comp, 19) = I1 < 0 ? -std::log(-I1):std::log(I1);
    foo(i_comp, 20) = I2 < 0 ? -std::log(-I2):std::log(I2);
    foo(i_comp, 21) = I3 < 0 ? -std::log(-I3):std::log(I3);
    foo(i_comp, 22) = I4 < 0 ? -std::log(-I4):std::log(I4);
    foo(i_comp, 23) = I5 < 0 ? -std::log(-I5):std::log(I5);
    foo(i_comp, 24) = I6 < 0 ? -std::log(-I6):std::log(I6);
    foo(i_comp, 25) = I7 < 0 ? -std::log(-I7):std::log(I7);
    // foo(i_comp, 26) is Mean int;
    // foo(i_comp, 27) is Min int;
    // foo(i_comp, 28) is Max int;
    foo(i_comp, 32) = (cx - 1) / std::sqrt(mag);
    foo(i_comp, 33) = (cy - 1) / std::sqrt(mag);
    foo(i_comp, 34) = (cx_int - 1) / std::sqrt(mag);
    foo(i_comp, 35) = (cy_int - 1) / std::sqrt(mag);
  }
  
  // FIXME: should kurtosis and skewness be intensity weighted ?
  for(R_len_t i = 0; i < img.size(); i++) { 
    R_len_t i_comp = msk[i] - 1;
    if((i_comp >= 0) && (i_comp < nC)) {
      double diff = img[i] - foo(i_comp, 26);
      foo(i_comp, 29) += diff * diff;
      foo(i_comp, 30) += diff * diff * diff;
      foo(i_comp, 31) += diff * diff * diff * diff;
    }
  }
  
  for(R_len_t i_comp = 0; i_comp < nC; i_comp++) {
    foo(i_comp, 29) = std::sqrt(foo(i_comp, 29) / foo(i_comp, 18));
    foo(i_comp, 30) = foo(i_comp, 30) / (foo(i_comp, 29) * foo(i_comp, 29) * foo(i_comp, 29)) / foo(i_comp, 18);
    foo(i_comp, 31) = foo(i_comp, 31) / (foo(i_comp, 29) * foo(i_comp, 29) * foo(i_comp, 29) * foo(i_comp, 29)) / foo(i_comp, 18);
  } 
  
  // return results
  Rcpp::NumericMatrix out = foo(Rcpp::_, Rcpp::Range(0,35));
  Rcpp::colnames(out) = N;
  return out;
}

#endif

