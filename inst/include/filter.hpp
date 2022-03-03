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
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_sd(const Rcpp::NumericMatrix mat,
                           const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // R_len_t i_row, i_col, f_row, f_col, k = 0;
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0);
  NumericMatrix out = pad["out"];
  NumericMatrix foo = clone(out);
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  R_len_t k = 0;
  for(R_len_t i = 0; i < kernel.size(); i++) if(kernel[i]) k++;
  NumericVector K(k);
  
  for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
      k = 0;
      for(R_len_t f_col = -pad_c; f_col <= pad_c; f_col++) {
        for(R_len_t f_row = -pad_r; f_row <= pad_r; f_row++) {
          if(kernel(pad_r + f_row, pad_c + f_col)) {
            K[k++] = foo(i_row + f_row, i_col + f_col);
          }
        }
      }
      out(i_row, i_col) = Rcpp::sd(K);
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
}

//' @title Image Mean Filtering
//' @name cpp_mean
//' @description
//' This function applies mean filtering on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_mean(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // R_len_t i_row, i_col, f_row, f_col, k = 0;
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0);
  NumericMatrix out = pad["out"];
  NumericMatrix foo = clone(out);
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  R_len_t k = 0;
  for(R_len_t i = 0; i < kernel.size(); i++) if(kernel[i]) k++;
  NumericVector K(k);
  
  for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
      k = 0;
      for(R_len_t f_col = -pad_c; f_col <= pad_c; f_col++) {
        for(R_len_t f_row = -pad_r; f_row <= pad_r; f_row++) {
          if(kernel(pad_r + f_row, pad_c + f_col)) {
            K[k++] = foo(i_row + f_row, i_col + f_col);
          }
        }
      }
      out(i_row, i_col) = mean(K);
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_median(const Rcpp::NumericMatrix mat,
                               const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // R_len_t i_row, i_col, f_row, f_col, k = 0;
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0);
  NumericMatrix out = pad["out"];
  NumericMatrix foo = clone(out);
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  R_len_t k = 0;
  for(R_len_t i = 0; i < kernel.size(); i++) if(kernel[i]) k++;
  NumericVector K(k);
  
  for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
      k = 0;
      for(R_len_t f_col = -pad_c; f_col <= pad_c; f_col++) {
        for(R_len_t f_row = -pad_r; f_row <= pad_r; f_row++) {
          if(kernel(pad_r + f_row, pad_c + f_col)) {
            K[k++] = foo(i_row + f_row, i_col + f_col);
          }
        }
      }
      out(i_row, i_col) = median(K);
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_mode(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // R_len_t i_row, i_col, f_row, f_col;
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0);
  NumericMatrix out = pad["out"];
  NumericMatrix foo = clone(out);
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  R_len_t k = 0;
  for(R_len_t i = 0; i < kernel.size(); i++) if(kernel[i]) k++;
  NumericVector K(k);
  
  for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
      k = 0;
      for(R_len_t f_col = -pad_c; f_col <= pad_c; f_col++) {
        for(R_len_t f_row = -pad_r; f_row <= pad_r; f_row++) {
          if(kernel(pad_r + f_row, pad_c + f_col)) {
            K[k++] = foo(i_row + f_row, i_col + f_col);
          }
        }
      }
      out(i_row, i_col) = 3 * median(K) - 2 * mean(K);
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_mid(const Rcpp::NumericMatrix mat,
                            const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // R_len_t i_row, i_col, f_row, f_col;
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0);
  NumericMatrix out = pad["out"];
  NumericMatrix foo = clone(out);
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
      double MIN = R_PosInf, MAX = R_NegInf;
      for(R_len_t f_col = -pad_c; f_col <= pad_c; f_col++) {
        for(R_len_t f_row = -pad_r; f_row <= pad_r; f_row++) {
          if(kernel(pad_r + f_row, pad_c + f_col)) {
            MIN = std::min(foo(i_row + f_row, i_col + f_col), MIN);
            MAX = std::max(foo(i_row + f_row, i_col + f_col), MAX);
          }
        }
      }
      out(i_row, i_col) = (MAX - MIN)/2;
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_convolve2d(const Rcpp::NumericMatrix mat,
                                   const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  
  // create padding
  // 6 add closest row / col + copy mat in the center
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0); 
  NumericMatrix foo = pad["out"];
  NumericMatrix out(foo.nrow(), foo.ncol());
  // NumericMatrix out = foo;
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  // apply filtering
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      for(R_len_t f_col = i_col + pad_c, i_ker = 0; f_col >= i_col - pad_c; f_col--) {
        for(R_len_t f_row = i_row + pad_r; f_row >= i_row - pad_r; f_row--) {
          out(i_row, i_col) = kernel[i_ker++] * foo(f_row, f_col) + out(i_row, i_col);
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_correlate2d(const Rcpp::NumericMatrix mat,
                                    const Rcpp::NumericMatrix kernel) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();

  // create padding
  // 6 add closest row / col + copy mat in the center
  Rcpp::List pad = hpp_padding(mat, kernel, 6, 0.0);
  NumericMatrix foo = pad["out"];
  NumericMatrix out(foo.nrow(), foo.ncol());
  // NumericMatrix out = foo;
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];

  // apply filtering
  for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
    for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
      for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
        for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
          out(i_row, i_col) = kernel[i_ker++] * foo(f_row, f_col) + out(i_row, i_col);
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
}

# endif
