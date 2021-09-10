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

#ifndef IFCIP_MORPHOLOGY_HPP
#define IFCIP_MORPHOLOGY_HPP

#include <Rcpp.h>
#include "padding.hpp"
using namespace Rcpp;

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
Rcpp::NumericMatrix hpp_erode(const Rcpp::NumericMatrix mat,
                              const Rcpp::NumericMatrix kernel,
                              const uint8_t iter = 0) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::List pad = hpp_padding(mat, kernel, 5, R_PosInf);
  NumericMatrix out = pad["out"];
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = R_PosInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) < K) K = foo(f_row, f_col);
            }
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_dilate(const Rcpp::NumericMatrix mat,
                               const Rcpp::NumericMatrix kernel,
                               const uint8_t iter = 0) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::List pad = hpp_padding(mat, kernel, 5, R_NegInf);
  NumericMatrix out = pad["out"];
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = R_NegInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) > K) K = foo(f_row, f_col);
            }
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_opening(const Rcpp::NumericMatrix mat,
                                const Rcpp::NumericMatrix kernel,
                                const uint8_t iter = 0) {
  // return hpp_dilate(hpp_erode(mat, kernel, iter), kernel, iter);
  uint8_t i_iter;
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::List pad = hpp_padding(mat, kernel, 5, R_PosInf);
  NumericMatrix out = pad["out"];
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  // erode
  i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = R_PosInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) < K) K = foo(f_row, f_col);
            }
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  // change border value
  // 1st cols
  for(R_len_t i_col = 0; i_col < pad_c; i_col++) {
    for(R_len_t i_row = 0; i_row < out.nrow(); i_row++) {
      out(i_row, i_col) = R_NegInf;
    }
  }
  // last cols
  for(R_len_t i_col = pad_c + mat_c; i_col < out.ncol(); i_col++) {
    for(R_len_t i_row = 0; i_row < out.nrow(); i_row++) {
      out(i_row, i_col) = R_NegInf;
    }
  }
  
  // 1st rows
  for(R_len_t i_col = 0; i_col < out.ncol(); i_col++) {
    for(R_len_t i_row = 0; i_row < pad_r; i_row++) {
      out(i_row, i_col) = R_NegInf; 
    }
  }
  
  // last rows
  for(R_len_t i_col = 0; i_col < out.ncol(); i_col++) {
    for(R_len_t i_row = pad_r + mat_r; i_row < out.nrow(); i_row++) {
      out(i_row, i_col) = R_NegInf; 
    }
  }
  
  // dilate
  i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = R_NegInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) > K) K = foo(f_row, f_col);
            }
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_closing(const Rcpp::NumericMatrix mat,
                                const Rcpp::NumericMatrix kernel,
                                const uint8_t iter = 0) {
  // return hpp_erode(hpp_dilate(mat, kernel, iter), kernel, iter); 
  uint8_t i_iter;
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::List pad = hpp_padding(mat, kernel, 5, R_NegInf);
  NumericMatrix out = pad["out"];
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  // dilate
  i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = R_NegInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) > K) K = foo(f_row, f_col);
            }
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  
  // change border value
  // 1st cols
  for(R_len_t i_col = 0; i_col < pad_c; i_col++) {
    for(R_len_t i_row = 0; i_row < out.nrow(); i_row++) {
      out(i_row, i_col) = R_PosInf;
    }
  }
  // last cols
  for(R_len_t i_col = pad_c + mat_c; i_col < out.ncol(); i_col++) {
    for(R_len_t i_row = 0; i_row < out.nrow(); i_row ++) {
      out(i_row, i_col) = R_PosInf;
    }
  }
  
  // 1st rows
  for(R_len_t i_col = 0; i_col < out.ncol(); i_col++) {
    for(R_len_t i_row = 0; i_row < pad_r; i_row++) {
      out(i_row, i_col) = R_PosInf; 
    }
  }
  
  // last rows
  for(R_len_t i_col = 0; i_col < out.ncol(); i_col++) {
    for(R_len_t i_row = pad_r + mat_r; i_row < out.nrow(); i_row++) {
      out(i_row, i_col) = R_PosInf; 
    }
  }
  
  // erode
  i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = R_PosInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) < K) K = foo(f_row, f_col);
            }
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
Rcpp::NumericMatrix hpp_gradient(const Rcpp::NumericMatrix mat,
                                 const Rcpp::NumericMatrix kernel,
                                 const uint8_t iter = 0) {
  Rcpp::NumericMatrix dil = hpp_dilate(mat, kernel, iter);
  Rcpp::NumericMatrix ero = hpp_erode(mat, kernel, iter);
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  for(R_len_t i_out = 0; i_out < mat_r * mat_c; i_out++) out[i_out] = dil[i_out] - ero[i_out];
  return out;
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
Rcpp::NumericMatrix hpp_tophat_white(const Rcpp::NumericMatrix mat,
                                     const Rcpp::NumericMatrix kernel,
                                     const uint8_t iter = 0) {
  Rcpp::NumericMatrix ope = hpp_opening(mat, kernel, iter);
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  for(R_len_t i_out = 0; i_out < mat_r * mat_c; i_out++) out[i_out] = mat[i_out] - ope[i_out];
  return out;
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
Rcpp::NumericMatrix hpp_tophat_black(const Rcpp::NumericMatrix mat,
                                     const Rcpp::NumericMatrix kernel,
                                     const uint8_t iter = 0) {
  Rcpp::NumericMatrix clo = hpp_closing(mat, kernel, iter);
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  for(R_len_t i_out = 0; i_out < mat_r * mat_c; i_out++) out[i_out] = clo[i_out] - mat[i_out];
  return out;
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
Rcpp::NumericMatrix hpp_tophat_self(const Rcpp::NumericMatrix mat,
                                    const Rcpp::NumericMatrix kernel,
                                    const uint8_t iter = 0) {
  Rcpp::NumericMatrix ope = hpp_closing(mat, kernel, iter);
  Rcpp::NumericMatrix clo = hpp_opening(mat, kernel, iter);
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  for(R_len_t i_out = 0; i_out < mat_r * mat_c; i_out++) out[i_out] = clo[i_out] - ope[i_out];
  return out;
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
Rcpp::NumericMatrix hpp_cont(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel,
                             const uint8_t iter = 0) {
  Rcpp::NumericMatrix ope = hpp_closing(mat, kernel, iter);
  Rcpp::NumericMatrix clo = hpp_opening(mat, kernel, iter);
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::NumericMatrix out(mat_r, mat_c);
  for(R_len_t i_out = 0; i_out < mat_r * mat_c; i_out++) out[i_out] = 3 * mat[i_out] - clo[i_out] - ope[i_out];
  return out;
}

# endif
