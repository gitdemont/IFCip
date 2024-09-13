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

#ifndef IFCIP_FILTER_HPP
#define IFCIP_FILTER_HPP

#include <Rcpp.h>
#include "padding.hpp"
using namespace Rcpp;

//' @title Image Standard Deviation Filtering
//' @name cpp_sd
//' @description
//' This function applies standard deviation filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = NA_REAL, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_sd(const Rcpp::NumericMatrix mat,
                           const Rcpp::NumericMatrix kernel,
                           const uint8_t method = 5) {
  // get kernel dimension 
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                    //APPLY OFFSET CORRECTION
  if(kr) kr = !kr;                                                               //APPLY OFFSET CORRECTION
  if(kc) kc = !kc;                                                               //APPLY OFFSET CORRECTION
  // create out padded with NA, NA is important to allow removal of extra values from final result (na_omit)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, method, NA_REAL);//APPLY OFFSET CORRECTION
  // create vector for storage of matrix value at kernel position
  Rcpp::NumericVector K(kernel.size(), NA_REAL);
  
  Rcpp::NumericMatrix foo = clone(out);
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      R_len_t i_ker = -1;
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += 1;
          K[i_ker] = kernel[i_ker] ? foo(f_row, f_col) : NA_REAL;
        }
      }
      Rcpp::NumericVector KK = Rcpp::na_omit(K);
      out(i_row, i_col) = KK.size() == 0 ? NA_REAL : Rcpp::sd(KK);
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Mean Filtering
//' @name cpp_mean
//' @description
//' This function applies mean filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = NA_REAL, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_mean(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel,
                             const uint8_t method = 5) {
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                    //APPLY OFFSET CORRECTION
  if(kr) kr = !kr;                                                               //APPLY OFFSET CORRECTION
  if(kc) kc = !kc;                                                               //APPLY OFFSET CORRECTION
  // create out padded with NA, NA is important to allow removal of extra values from final result (na_omit)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, method, NA_REAL);//APPLY OFFSET CORRECTION
  // create vector for storage of matrix value at kernel position
  Rcpp::NumericVector K(kernel.size(), NA_REAL);
  
  Rcpp::NumericMatrix foo = clone(out);
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      R_len_t i_ker = -1;
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += 1;
          K[i_ker] = kernel[i_ker] ? foo(f_row, f_col) : NA_REAL;
        }
      }
      Rcpp::NumericVector KK = Rcpp::na_omit(K);
      out(i_row, i_col) = KK.size() == 0 ? NA_REAL : Rcpp::mean(KK);
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Median Filtering
//' @name cpp_median
//' @description
//' This function applies median filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = NA_REAL, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_median(const Rcpp::NumericMatrix mat,
                               const Rcpp::NumericMatrix kernel,
                               const uint8_t method = 5) {
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                    //APPLY OFFSET CORRECTION
  if(kr) kr = !kr;                                                               //APPLY OFFSET CORRECTION
  if(kc) kc = !kc;                                                               //APPLY OFFSET CORRECTION
  // create out padded with NA, NA is important to allow removal of extra values from final result (na_omit)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, method, NA_REAL);//APPLY OFFSET CORRECTION
  // create vector for storage of matrix value at kernel position
  Rcpp::NumericVector K(kernel.size(), NA_REAL);
  
  Rcpp::NumericMatrix foo = clone(out);
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      R_len_t i_ker = -1;
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += 1;
          K[i_ker] = kernel[i_ker] ? foo(f_row, f_col) : NA_REAL;
        }
      }
      Rcpp::NumericVector KK = Rcpp::na_omit(K);
      out(i_row, i_col) = KK.size() == 0 ? NA_REAL : Rcpp::median(KK);
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Mode Filtering
//' @name cpp_mode
//' @description
//' This function applies mode filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = NA_REAL, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_mode(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel,
                             const uint8_t method = 5) {
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                    //APPLY OFFSET CORRECTION
  if(kr) kr = !kr;                                                               //APPLY OFFSET CORRECTION
  if(kc) kc = !kc;                                                               //APPLY OFFSET CORRECTION
  // create out padded with NA, NA is important to allow removal of extra values from final result (na_omit)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, method, NA_REAL);//APPLY OFFSET CORRECTION
  // create vector for storage of matrix value at kernel position
  Rcpp::NumericVector K(kernel.size(), NA_REAL);
  
  Rcpp::NumericMatrix foo = clone(out);
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      R_len_t i_ker = -1;
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += 1;
          K[i_ker] = kernel[i_ker] ? foo(f_row, f_col) : NA_REAL;
        }
      }
      Rcpp::NumericVector KK = Rcpp::na_omit(K);
      out(i_row, i_col) = KK.size() == 0 ? NA_REAL : 3 * median(KK) - 2 * mean(KK);
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Mid Filtering
//' @name cpp_mid
//' @description
//' This function applies mid filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = NA_REAL, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_mid(const Rcpp::NumericMatrix mat,
                            const Rcpp::NumericMatrix kernel,
                            const uint8_t method = 5) {
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                    //APPLY OFFSET CORRECTION
  if(kr) kr = !kr;                                                               //APPLY OFFSET CORRECTION
  if(kc) kc = !kc;                                                               //APPLY OFFSET CORRECTION
  // create out padded with NA, NA is important to allow removal of extra values from final result (na_omit)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, method, NA_REAL);//APPLY OFFSET CORRECTION
  // create vectors for storage of matrix value at kernel position
  Rcpp::NumericVector MIN(kernel.size(), NA_REAL);
  Rcpp::NumericVector MAX(kernel.size(), NA_REAL);
  
  Rcpp::NumericMatrix foo = clone(out);
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      R_len_t i_ker = -1;
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += 1;
          if(kernel[i_ker]) {
            MAX[i_ker] = MIN[i_ker] = foo(f_row, f_col);
          } else {
            MAX[i_ker] = MIN[i_ker] = NA_REAL;
          }
        }
      }
      Rcpp::NumericVector KMAX = Rcpp::na_omit(MAX);
      Rcpp::NumericVector KMIN = Rcpp::na_omit(MIN);
      out(i_row, i_col) = (KMAX.size() * KMIN.size()) == 0 ? NA_REAL : (Rcpp::max(KMAX) - Rcpp::min(KMIN))/2;
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Filtering by Convolution
//' @name cpp_convolve2d
//' @description
//' This function applies 2D convolution filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = 0.0, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_convolve2d(const Rcpp::NumericMatrix mat,
                                   const Rcpp::NumericMatrix kernel,
                                   const uint8_t method = 5) {
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                //APPLY OFFSET CORRECTION
  if(kr) kr = !kr;                                                           //APPLY OFFSET CORRECTION
  if(kc) kc = !kc;                                                           //APPLY OFFSET CORRECTION
  // create out padded with 0.0, 0.0 is important to allow removal of extra values from final result (multiply by 0.0)
  Rcpp::NumericMatrix foo = hpp_padding(mat, pad_r + kc, pad_c + kr, method, 0.0);//APPLY OFFSET CORRECTION
  
  Rcpp::NumericMatrix out(foo.nrow(), foo.ncol());
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      R_len_t i_ker = -1;
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += 1;
          out(i_row, i_col) = kernel[i_ker] * foo(f_row, f_col) + out(i_row, i_col);
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Filtering by Correlation
//' @name cpp_correlate2d
//' @description
//' This function applies 2D correlation filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method used for padding, a uint8_t. Default is 5 (with k = 0.0, see hpp_padding), allowed are [1-8].
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_correlate2d(const Rcpp::NumericMatrix mat,
                                    const Rcpp::NumericMatrix kernel, 
                                    const uint8_t method = 5) {
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                 //APPLY OFFSET CORRECTION
  // create out padded with 0.0, 0.0 is important to allow removal of extra values from final result (multiply by 0.0)
  Rcpp::NumericMatrix foo = hpp_padding(mat, pad_r + kc, pad_c + kr, method,  0.0);//APPLY OFFSET CORRECTION
  
  Rcpp::NumericMatrix out(foo.nrow(), foo.ncol());
  for(R_len_t i_col = out.ncol() - pad_c - 1; i_col >= pad_c; i_col--) {
    for(R_len_t i_row = out.nrow() - pad_r - 1; i_row >= pad_r; i_row--) {
      R_len_t i_ker = kernel.size();
      for(R_len_t f_col = i_col + pad_c_1 - 1; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--) {
          i_ker += -1;
          out(i_row, i_col) = kernel[i_ker] * foo(f_row, f_col) + out(i_row, i_col);
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

# endif
