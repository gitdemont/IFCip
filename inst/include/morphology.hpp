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
#include <utils.hpp>
#include "kernel.hpp"
#include "padding.hpp"
#include "uw.hpp"
using namespace Rcpp;

//' @title Image Erosion
//' @name cpp_erode
//' @description
//' This function applies erosion on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @details see 'Efficient 2-D grayscale morphological transformations with arbitrary flat structuring elements' from  E.R. Urbach, M.H.F. Wilkinson.
//' IEEE Transactions on Image Processing, 17(1):1-8, January 2008.\doi{10.1109/tip.2007.912582}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_erode(const Rcpp::NumericMatrix mat,
                              const Rcpp::NumericMatrix kernel) {
  return hpp_uw(mat, kernel, true);
}

//' @title Image Dilatation
//' @name cpp_dilate
//' @description
//' This function applies dilatation on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @details see 'Efficient 2-D grayscale morphological transformations with arbitrary flat structuring elements' from  E.R. Urbach, M.H.F. Wilkinson.
//' IEEE Transactions on Image Processing, 17(1):1-8, January 2008.\doi{10.1109/tip.2007.912582}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_dilate(const Rcpp::NumericMatrix mat,
                               const Rcpp::NumericMatrix kernel) {
  return hpp_uw(mat, kernel, false);
}

//' @title Image Opening
//' @name cpp_opening
//' @description
//' This function applies opening on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_opening(const Rcpp::NumericMatrix mat,
                                const Rcpp::NumericMatrix kernel) {
  return hpp_dilate(hpp_erode(mat, kernel), kernel);
}

//' @title Image Closing
//' @name cpp_closing
//' @description
//' This function applies closing on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_closing(const Rcpp::NumericMatrix mat,
                                const Rcpp::NumericMatrix kernel) {
  return hpp_erode(hpp_dilate(mat, kernel), kernel);
}

//' @title Image Morphological Gradient
//' @name cpp_gradient
//' @description
//' This function applies morphological gradient on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_gradient(const Rcpp::NumericMatrix mat,
                                 const Rcpp::NumericMatrix kernel) {
  Rcpp::NumericMatrix dil = hpp_dilate(mat, kernel);
  Rcpp::NumericMatrix ero = hpp_erode(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) dil[i] -= ero[i];
  return dil;
}

//' @title Image White Top Hat
//' @name cpp_tophat_white
//' @description
//' This function applies white top hat on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_tophat_white(const Rcpp::NumericMatrix mat,
                                     const Rcpp::NumericMatrix kernel) {
  Rcpp::NumericMatrix ope = hpp_opening(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) ope[i] = mat[i] - ope[i];
  return ope;
}

//' @title Image Black Top Hat
//' @name cpp_tophat_black
//' @description
//' This function applies black top hat on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_tophat_black(const Rcpp::NumericMatrix mat,
                                     const Rcpp::NumericMatrix kernel) {
  Rcpp::NumericMatrix clo = hpp_closing(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) clo[i] -= mat[i];
  return clo;
}

//' @title Image Self Complementary Top Hat
//' @name cpp_tophat_self
//' @description
//' This function applies self complementary on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_tophat_self(const Rcpp::NumericMatrix mat,
                                    const Rcpp::NumericMatrix kernel) {
  Rcpp::NumericMatrix clo = hpp_opening(mat, kernel);
  Rcpp::NumericMatrix ope = hpp_closing(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) clo[i] -= ope[i];
  return clo;
}

//' @title Image Contrast Enhancement
//' @name cpp_cont
//' @description
//' This function applies contrast enhancement on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_cont(const Rcpp::NumericMatrix mat,
                             const Rcpp::NumericMatrix kernel) {
  Rcpp::NumericMatrix ope = hpp_closing(mat, kernel);
  Rcpp::NumericMatrix clo = hpp_opening(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) ope[i] = 3 * mat[i] - clo[i] - ope[i];
  return ope;
}

//' @title Image Laplacian
//' @name cpp_laplacian
//' @description
//' This function applies Laplacian morphology on image.
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_laplacian(const Rcpp::NumericMatrix mat,
                                  const Rcpp::NumericMatrix kernel) {
  Rcpp::NumericMatrix dil = hpp_dilate(mat, kernel);
  Rcpp::NumericMatrix ero = hpp_erode(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) dil[i] += ero[i] - 2 * mat[i];
  return dil;
}

//' @title Brute Force Image Erosion
//' @name cpp_erode_old
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
Rcpp::NumericMatrix hpp_erode_old(const Rcpp::NumericMatrix mat,
                                  const Rcpp::NumericMatrix kernel,
                                  const uint8_t iter = 0,
                                  const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                     //APPLY OFFSET CORRECTION
  // create out padded with P_PosInf, P_PosInf is important to allow removal of extra values from final result (we will keep minimum)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, 5, R_PosInf);//APPLY OFFSET CORRECTION
  
  // check msk_ is provided and its dimension
  if(msk_.isNotNull()) {
    Rcpp::NumericMatrix m(msk_.get());
    if(m.size() != 0) {
      if((m.nrow() != mat_r) || (m.ncol() != mat_c)) {
        Rcpp::stop("hpp_erode: when 'msk' is provided 'img' and 'msk' should have same dimensions");
      }
      Rcpp::NumericMatrix msk = hpp_padding(m, pad_r + kc, pad_c + kr, 5, 1.0);
      uint8_t i_iter = 0;
      while(i_iter <= iter) {
        i_iter++;
        Rcpp::NumericMatrix foo = Rcpp::clone(out);
        for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
          for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
            if(msk(i_row, i_col)) {
              double K = R_PosInf;
              for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col < i_col + pad_c_1; f_col++) {
                for(R_len_t f_row = i_row - pad_r; f_row < i_row + pad_r_1; f_row++) {
                  if(kernel[i_ker++] && msk(f_row, f_col) && ((foo(f_row, f_col) < K))) K = foo(f_row, f_col);
                }
              }
              out(i_row, i_col) = K;
            }
          }
        }
      }
      return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
    }
  }
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    Rcpp::NumericMatrix foo = Rcpp::clone(out);
    for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
        double K = R_PosInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col < i_col + pad_c_1; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row < i_row + pad_r_1; f_row++) {
            if(kernel[i_ker++] && ((foo(f_row, f_col) < K))) K = foo(f_row, f_col);
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
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
Rcpp::NumericMatrix hpp_dilate_old(const Rcpp::NumericMatrix mat,
                                   const Rcpp::NumericMatrix kernel,
                                   const uint8_t iter = 0,
                                   const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  // get kernel dimension
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;                                                     //APPLY OFFSET CORRECTION
  // create out padded with R_NegInf, R_NegInf is important to allow removal of extra values from final result (we will keep maximum)
  Rcpp::NumericMatrix out = hpp_padding(mat, pad_r + kc, pad_c + kr, 5, R_NegInf);//APPLY OFFSET CORRECTION
  
  // check msk_ is provided and its dimension
  if(msk_.isNotNull()) {
    Rcpp::NumericMatrix m(msk_.get());
    if(m.size() != 0) {
      if((m.nrow() != mat_r) || (m.ncol() != mat_c)) {
        Rcpp::stop("hpp_dilate: when 'msk' is provided 'img' and 'msk' should have same dimensions");
      }
      Rcpp::NumericMatrix msk = hpp_padding(m, pad_r + kc, pad_c + kr, 5, 1.0);
      uint8_t i_iter = 0;
      while(i_iter <= iter) {
        i_iter++;
        Rcpp::NumericMatrix foo = Rcpp::clone(out);
        for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
          for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
            if(msk(i_row, i_col)) {
              double K = R_NegInf;
              for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col < i_col + pad_c_1; f_col++) {
                for(R_len_t f_row = i_row - pad_r; f_row < i_row + pad_r_1; f_row++) {
                  if(kernel[i_ker++] && msk(f_row, f_col) && (foo(f_row, f_col) > K)) K = foo(f_row, f_col);
                }
              }
              out(i_row, i_col) = K;
            }
          }
        }
      }
      return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
    }
  }
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    Rcpp::NumericMatrix foo = Rcpp::clone(out);
    for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
        double K = R_NegInf;
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col < i_col + pad_c_1; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row < i_row + pad_r_1; f_row++) {
            if(kernel[i_ker++] && (foo(f_row, f_col) > K)) K = foo(f_row, f_col);
          }
        }
        out(i_row, i_col) = K;
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

# endif
