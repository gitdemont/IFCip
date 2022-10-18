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
  
#ifndef IFCIP_DISTANCE_HPP
#define IFCIP_DISTANCE_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Distance Feature Extraction
//' @name hpp_dist_cen
//' @description
//' This function is designed to extract features for distance computation.
//' @param msk an IntegerMatrix, containing connected components.
//' @return a NumericMatrix whose rows are component numbers and columns are:\cr
//' -pix cx, x centroid of the component in pixels\cr
//' -pix cy, y centroid of the component in pixels\cr
//' -pix count, number of pixels occupied by the component.
//' @keywords internal
Rcpp::NumericMatrix hpp_dist_cen(const Rcpp::IntegerMatrix msk) {
  R_len_t mat_r = msk.nrow(), mat_c = msk.ncol();
  R_len_t nC = 0;
  R_len_t alw = (std::ceil(mat_r / 2) + 1) * (std::ceil(mat_c / 2) + 1);
  for(R_len_t i_msk = 0; i_msk < mat_c * mat_r; i_msk++) {
    if(msk[i_msk] < 0) Rcpp::stop("hpp_dist_cen: invalid negative value for 'msk'");
    if(msk[i_msk] >= alw) Rcpp::stop("hpp_dist_cen: invalid max number of components for 'msk'");
    nC = std::max(nC, msk[i_msk]); 
  }
  
  // initialize matrix
  Rcpp::NumericMatrix out(nC, 3);
  
  // create colnames  (not needed)
  // Rcpp::StringVector N = Rcpp::StringVector::create("pix cx", "pix cy", "pix count");
  // Rcpp::colnames(out) = N;
  if(nC == 0) return out;

  // define raw moments
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    if(Rcpp::is_true(Rcpp::any(Rcpp::is_na(msk(Rcpp::_,i_col))))) Rcpp::stop("hpp_dist_cen: NA - NaN value are not allowed in 'msk'");
    R_len_t i_col_1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_row_1 = i_row + 1;
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp >= 0) && (i_comp < nC)) {
        // M01
        out(i_comp, 0) += i_col_1;
        // M10
        out(i_comp, 1) += i_row_1;
        // M00
        out(i_comp, 2)++;
      }
    }
  }
  
  for(R_len_t i_comp = 0; i_comp < nC; i_comp++) {
    out(i_comp, 0) = out(i_comp, 0) / out(i_comp, 2);
    out(i_comp, 1) = out(i_comp, 1) / out(i_comp, 2);
  }
  
  // return results
  return out;
}

//' @title Distance Intensity Extraction
//' @name hpp_dist_int
//' @description
//' This function is designed to extract intensity for distance computation.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk an IntegerMatrix, containing connected components.
//' @param nC an unsigned integer. Maximal component number to retrieve centroids about.
//' Default is 0 to retrieve centroids for all components.
//' @return a NumericMatrix whose rows are component numbers and columns are:\cr
//' -pix cx, x centroid of the component in pixels\cr
//' -pix cy, y centroid of the component in pixels\cr
//' -pix count, number of pixels occupied by the component.
//' @keywords internal
Rcpp::NumericMatrix hpp_dist_int(const Rcpp::NumericMatrix img,
                                 const Rcpp::IntegerMatrix msk,
                                 const R_len_t nC = 0) {
  R_len_t mat_r = msk.nrow(), mat_c = msk.ncol();
  
  // initialize matrix
  Rcpp::NumericMatrix out(nC, 2);
  
  // create colnames (not needed)
  // Rcpp::StringVector N = Rcpp::StringVector::create("Raw Min Pixel","Raw Max Pixel");
  // Rcpp::colnames(out) = N;
  if(nC == 0) return out;
  
  // fill min, max
  for(R_len_t i_comp = 0; i_comp < nC; i_comp++) {
    out(i_comp, 0) = R_PosInf;
    out(i_comp, 1) = R_NegInf;
  }
  
  // define raw moments
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp >= 0) && (i_comp < nC)) {
        // Min
        if(out(i_comp, 0) > img(i_row, i_col)) out(i_comp, 0) = img(i_row, i_col);
        // Max
        if(out(i_comp, 1) < img(i_row, i_col)) out(i_comp, 1) = img(i_row, i_col);
      }
    }
  }
  // return results
  return out;
}

//' @title Mask Euclidean Distance
//' @name cpp_distance_eucl
//' @description
//' This function is designed to compute Euclidean distance from background to centroïds' foreground
//' @param msk an IntegerMatrix, containing connected components.
//' @return an NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_distance_eucl(const Rcpp::IntegerMatrix msk) {
  Rcpp::NumericMatrix cen = hpp_dist_cen(msk);
  R_len_t cen_r = cen.nrow(), mat_r = msk.nrow(), mat_c = msk.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  if(cen_r == 0) return out;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp < cen_r) && (i_comp >= 0)) {
        out(i_row, i_col) = std::sqrt((cen(i_comp, 0) - i_col - 1) * (cen(i_comp, 0) - i_col - 1) +
          (cen(i_comp, 1) - i_row - 1) * (cen(i_comp, 1) - i_row - 1));
      }
    }
  }
  return out;
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
Rcpp::NumericMatrix hpp_distance_eucl_norm(const Rcpp::IntegerMatrix msk) {
  Rcpp::NumericMatrix cen = hpp_dist_cen(msk);
  R_len_t cen_r = cen.nrow(), mat_r = msk.nrow(), mat_c = msk.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  if(cen_r == 0) return out;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp < cen_r) && (i_comp >= 0)) {
        out(i_row, i_col) = std::sqrt((cen(i_comp, 0) - i_col - 1) * (cen(i_comp, 0) - i_col - 1) +
          (cen(i_comp, 1) - i_row - 1) * (cen(i_comp, 1) - i_row - 1));
      }
    }
  }
  Rcpp::NumericMatrix fea = hpp_dist_int(out, msk, cen_r);
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp < cen_r) && (i_comp >= 0)) {
        out(i_row, i_col) = (out(i_row, i_col) - fea(i_comp,0)) / (fea(i_comp,1) - fea(i_comp,0));
      }
    }
  }
  return out;
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
Rcpp::NumericMatrix hpp_distance_manh(const Rcpp::IntegerMatrix msk) {
  Rcpp::NumericMatrix cen = hpp_dist_cen(msk);
  R_len_t cen_r = cen.nrow(), mat_r = msk.nrow(), mat_c = msk.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  if(cen_r == 0) return out;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp < cen_r) && (i_comp >= 0)) {
        out(i_row, i_col) = std::abs(cen(i_comp, 0) - i_col - 1) + std::abs(cen(i_comp, 1) - i_row - 1);
      }
    }
  }
  return out;
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
Rcpp::NumericMatrix hpp_distance_manh_norm(const Rcpp::IntegerMatrix msk) {
  Rcpp::NumericMatrix cen = hpp_dist_cen(msk);
  R_len_t cen_r = cen.nrow(), mat_r = msk.nrow(), mat_c = msk.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  if(cen_r == 0) return out;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp < cen_r) && (i_comp >= 0)) {
        out(i_row, i_col) = std::abs(cen(i_comp, 0) - i_col - 1) + std::abs(cen(i_comp, 1) - i_row - 1);
      }
    }
  }
  Rcpp::NumericMatrix fea = hpp_dist_int(out, msk, cen_r);
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t i_comp = msk(i_row, i_col) - 1;
      if((i_comp < cen_r) && (i_comp >= 0)) {
        out(i_row, i_col) = (out(i_row, i_col) - fea(i_comp,0)) / (fea(i_comp,1) - fea(i_comp,0));
      }
    }
  }
  return out;
}

//' Implementation of Distance Transform Algorithm from A. Meijster et al.
//' https://doi.org/10.1007/0-306-47025-X_36

//' @title Manhattan Distance Transform
//' @name hpp_disttrans_1st
//' @description
//' This function performs the 1st pass of the distance transform of an image by implementing A. Meijster algorithm.
//' @param img, a NumericMatrix. /!\ Elements with values <=0 are considered as background.
//' @details adaptation of 'A General Algorithm For Computing Distance Transforms In Linear Time' from W.H. Hesselink, A. Meijster, J.B.T.M. Roerdink.
//' Mathematical Morphology and its Applications to Image and Signal Processing. February 2002, Pages 331-340.\doi{10.1007/0-306-47025-X_36}\cr
//' Values > 0 will be considered as foreground whereas all other values will be considered as background (i.e. 0).
//' @return a NumericMatrix.
//' @keywords internal
Rcpp::NumericMatrix hpp_disttrans_1st (const Rcpp::NumericMatrix img) {
  R_len_t n = img.nrow();
  R_len_t m = img.ncol();
  R_len_t mat_m = n + m;
  Rcpp::NumericMatrix out = Rcpp::no_init(n, m);
  
  for(R_len_t x = 0; x < m; x++) {
    out(0, x) = img(0, x) <= 0 ? 0 : mat_m;      // zero and negative values are considered as background !
    for(R_len_t y = 1; y < n; y++) {     // forward scan  * scan 1 *
      out(y, x) = img(y, x) <= 0 ? 0 : (1 + out(y - 1, x));
    }
    for(R_len_t y = n - 2; y >= 0; y--) {// backward scan * scan 2 *
      if(out(y + 1, x) < out(y, x)) out(y, x) = 1 + out(y + 1, x);
    }
  }
  return out;
}

R_len_t eucl_fun (R_len_t x, R_len_t i, R_len_t Gi) {
  return (x - i) * (x - i) + Gi * Gi;
}

double eucl_sep (R_len_t i, R_len_t u, R_len_t Gi, R_len_t Gu) {
  return (u * u - i * i + Gu * Gu - Gi * Gi) / (2 * (u - i));
}

R_len_t manh_fun (R_len_t x, R_len_t i, R_len_t Gi) {
  return std::max(std::abs(x - i), Gi);
}

double manh_sep (R_len_t i, R_len_t u, R_len_t Gi, R_len_t Gu) {
  if(Gi <= Gu) {
    return std::max(i + Gu, (i + u) / 2);
  } else {
    return std::min(u - Gi, (i + u) / 2);
  }
}

//' @title Manhattan Distance Transform
//' @name cpp_disttrans_manh
//' @description
//' This function computes the Manhattan distance transform of an image by implementing A. Meijster algorithm.
//' @param img, a NumericMatrix.
//' @details adaptation of 'A General Algorithm For Computing Distance Transforms In Linear Time' from W.H. Hesselink, A. Meijster, J.B.T.M. Roerdink.
//' Mathematical Morphology and its Applications to Image and Signal Processing. February 2002, Pages 331-340.\doi{10.1007/0-306-47025-X_36}\cr
//' Values > 0 will be considered as foreground whereas all other values will be considered as background (i.e. 0).
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_disttrans_manh (const Rcpp::NumericMatrix img) {
  R_len_t n = img.nrow();
  R_len_t m = img.ncol();
  
  Rcpp::IntegerVector s = Rcpp::no_init(m);
  Rcpp::IntegerVector t = Rcpp::no_init(m);
  Rcpp::NumericMatrix ori = hpp_disttrans_1st(img);
  Rcpp::NumericMatrix out = Rcpp::no_init(n, m);

  for(R_len_t y = 0; y < n; y++) {         // forward scan  * scan 3 *
    R_len_t q = 0;
    s[0] = 0;
    t[0] = 0;
    for(R_len_t u = 1; u < m; u++) {
      while((q >= 0) && (manh_fun(t[q], s[q], ori(y, s[q])) > manh_fun(t[q], u, ori(y, u)))) q--;
      if(q < 0) {
        q = 0;
        s[0] = u;
      } else {
        R_len_t w = 1 + manh_sep(s[q], u, ori(y, s[q]), ori(y, u));
        if(w < m) {
          q++;
          s[q] = u;
          t[q] = w;
        }
      }
    }
    for(R_len_t u = m - 1; u >= 0; u--) { // backward scan * scan 4 *
      out(y, u) = manh_fun(u, s[q], ori(y, s[q]));
      if(u == t[q]) q--;
    }
  }
  return out;
}

//' @title Euclidean Distance Transform
//' @name cpp_disttrans_eucl
//' @description
//' This function computes the Euclidean distance transform of an image by implementing A. Meijster algorithm.
//' @param img, a NumericMatrix.
//' @details adaptation of 'A General Algorithm For Computing Distance Transforms In Linear Time' from W.H. Hesselink, A. Meijster, J.B.T.M. Roerdink.
//' Mathematical Morphology and its Applications to Image and Signal Processing. February 2002, Pages 331-340.\doi{10.1007/0-306-47025-X_36}\cr
//' Values > 0 will be considered as foreground whereas all other values will be considered as background (i.e. 0).
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_disttrans_eucl (const Rcpp::NumericMatrix img) {
  R_len_t n = img.nrow();
  R_len_t m = img.ncol();
  
  Rcpp::IntegerVector s = Rcpp::no_init(m);
  Rcpp::IntegerVector t = Rcpp::no_init(m);
  Rcpp::NumericMatrix ori = hpp_disttrans_1st(img);
  Rcpp::NumericMatrix out = Rcpp::no_init(n, m);
  
  for(R_len_t y = 0; y < n; y++) {         // forward scan  * scan 3 *
    R_len_t q = 0;
    s[0] = 0;
    t[0] = 0;
    for(R_len_t u = 1; u < m; u++) {
      while((q >= 0) && (eucl_fun(t[q], s[q], ori(y, s[q])) > eucl_fun(t[q], u, ori(y, u)))) q--;
      if(q < 0) {
        q = 0;
        s[0] = u;
      } else {
        R_len_t w = 1 + eucl_sep(s[q], u, ori(y, s[q]), ori(y, u));
        if(w < m) {
          q++;
          s[q] = u;
          t[q] = w;
        }
      }
    }
    for(R_len_t u = m - 1; u >= 0; u--) { // backward scan * scan 4 *
      out(y, u) = std::sqrt(eucl_fun(u, s[q], ori(y, s[q])));
      if(u == t[q]) q--;
    }
  }
  return out;
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
Rcpp::IntegerMatrix hpp_voronoi_eucl (const Rcpp::IntegerMatrix img) {
  Rcpp::NumericMatrix cen = hpp_dist_cen(img);
  R_len_t cen_r = cen.nrow(), mat_r = img.nrow(), mat_c = img.ncol();
  Rcpp::IntegerMatrix out = Rcpp::no_init(mat_r, mat_c);
  if(cen_r == 0) return out;
  
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t i_col_1 = i_col - 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t lab = 0, dist = i_col * i_col + i_row * i_row;
      for(R_len_t i_cen = 0; i_cen < cen_r; i_cen++) {
        R_len_t d = pow(cen(i_cen, 0) - i_col_1, 2.0) + pow(cen(i_cen, 1) - i_row - 1, 2.0);
        if(d < dist) {
          dist = d;
          lab = i_cen;
        }
      }
      out(i_row, i_col) = lab;
    }
  }
  
  return out;
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
Rcpp::IntegerMatrix hpp_voronoi_manh (const Rcpp::IntegerMatrix img) {
  Rcpp::NumericMatrix cen = hpp_dist_cen(img);
  R_len_t cen_r = cen.nrow(), mat_r = img.nrow(), mat_c = img.ncol();
  Rcpp::IntegerMatrix out = Rcpp::no_init(mat_r, mat_c);
  if(cen_r == 0) return out;
  
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t i_col_1 = i_col - 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t lab = 0, dist = i_col + i_row;
      for(R_len_t i_cen = 0; i_cen < cen_r; i_cen++) {
        R_len_t d = std::abs(cen(i_cen, 0) - i_col_1) + std::abs(cen(i_cen, 1) - i_row - 1);
        if(d < dist) {
          dist = d;
          lab = i_cen;
        }
      }
      out(i_row, i_col) = lab;
    }
  }
  
  return out;
}

#endif
