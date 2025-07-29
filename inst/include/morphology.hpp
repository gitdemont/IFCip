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
#include "kernel.hpp"
#include "padding.hpp"
#include "uw.hpp"
using namespace Rcpp;

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> erode_T (Rcpp::Matrix<RTYPE> mat,
                             Rcpp::Matrix<RTYPEk> kernel) {
  return hpp_uw(mat, kernel, true);
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> dilate_T (Rcpp::Matrix<RTYPE> mat,
                              Rcpp::Matrix<RTYPEk> kernel) {
  return hpp_uw(mat, kernel, false);
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> opening_T (Rcpp::Matrix<RTYPE> mat,
                               Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::NumericMatrix rev_kernel = Rcpp::no_init_matrix(kernel.nrow(), kernel.ncol());
  for(R_len_t i = kernel.size() - 1, k = 0; i >= 0; i--) rev_kernel[k++] = kernel[i];
  return dilate_T(erode_T(mat, kernel), rev_kernel);
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> closing_T (Rcpp::Matrix<RTYPE> mat,
                               Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::NumericMatrix rev_kernel = Rcpp::no_init_matrix(kernel.nrow(), kernel.ncol());
  for(R_len_t i = kernel.size() - 1, k = 0; i >= 0; i--) rev_kernel[k++] = kernel[i];
  return erode_T(dilate_T(mat, kernel), rev_kernel);
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> gradient_T (Rcpp::Matrix<RTYPE> mat,
                                Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::Matrix<RTYPE> dil = dilate_T(mat, kernel);
  Rcpp::Matrix<RTYPE> ero = erode_T(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) dil[i] -= ero[i];
  return dil;
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> tophat_white_T (Rcpp::Matrix<RTYPE> mat,
                                    Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::Matrix<RTYPE> ope = opening_T(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) ope[i] = mat[i] - ope[i];
  return ope;
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> tophat_black_T (Rcpp::Matrix<RTYPE> mat,
                                    Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::Matrix<RTYPE> clo = closing_T(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) clo[i] -= mat[i];
  return clo;
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> tophat_self_T (Rcpp::Matrix<RTYPE> mat,
                                   Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::Matrix<RTYPE> clo = opening_T(mat, kernel);
  Rcpp::Matrix<RTYPE> ope = closing_T(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) clo[i] -= ope[i];
  return clo;
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> cont_T (Rcpp::Matrix<RTYPE> mat,
                            Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::Matrix<RTYPE> ope = closing_T(mat, kernel);
  Rcpp::Matrix<RTYPE> clo = opening_T(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) ope[i] = 3 * mat[i] - clo[i] - ope[i];
  return ope;
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> laplacian_T (Rcpp::Matrix<RTYPE> mat,
                                 Rcpp::Matrix<RTYPEk> kernel) {
  Rcpp::Matrix<RTYPE> dil = dilate_T(mat, kernel);
  Rcpp::Matrix<RTYPE> ero = erode_T(mat, kernel);
  for(R_len_t i = 0; i < mat.size(); i++) dil[i] += ero[i] - 2 * mat[i];
  return dil;
}

template <int RTYPE, int RTYPEk>
Rcpp::Matrix<RTYPE> fun_T(Rcpp::Matrix<RTYPE> mat, Rcpp::Matrix<RTYPEk> kernel, int what_n, std::string what) {
  switch(what_n) {
  case 1: return dilate_T(mat, kernel);
  case 2: return erode_T(mat, kernel);
  case 3: return opening_T(mat, kernel);
  case 4: return closing_T(mat, kernel);
  case 5: return tophat_white_T(mat, kernel);
  case 6: return tophat_black_T(mat, kernel);
  case 7: return tophat_self_T(mat, kernel);
  case 8: return gradient_T(mat, kernel);
  case 9: return cont_T(mat, kernel);
  case 10: return laplacian_T(mat, kernel);
  default: Rcpp::stop("hpp_mopho: morpho function['%s'] is not supported", what);
  }
}

// [[Rcpp::export(rng = false)]]
SEXP hpp_morpho(SEXP mat,
                SEXP kernel,
                const std::string what = "") {
  // match morpho function
  int fun = Rcpp::match(Rcpp::CharacterVector::create(what),
                        Rcpp::CharacterVector::create("dilate", "erode", "opening", "closing", "tophat_white", "tophat_black", "tophat_self", "gradient", "contrast", "laplacian"))[0];

  switch( TYPEOF(mat) ) {
  // case NILSXP : return R_NilValue; // prefer throwing error on NILSXP
  case LGLSXP : return fun_T(as<Rcpp::LogicalMatrix>(mat), get_kernel(kernel), fun, what);
  case RAWSXP : return fun_T(as<Rcpp::RawMatrix>(mat), get_kernel(kernel), fun, what);
  case INTSXP : return fun_T(as<Rcpp::IntegerMatrix>(mat), get_kernel(kernel), fun, what);
  case REALSXP :return fun_T(as<Rcpp::NumericMatrix>(mat), get_kernel(kernel), fun, what);
  default : Rcpp::stop("hpp_morpho: not supported SEXP[%i] in 'mat'", TYPEOF(mat));
  }
}

// [[Rcpp::export(rng = false)]]
SEXP hpp_dilate(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "dilate"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_erode(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "erode"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_opening(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "opening"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_closing(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "closing"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_tophat_white(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "tophat_white"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_tophat_black(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "tophat_black"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_tophat_self(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "tophat_self"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_gradient(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "gradient"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_cont(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "contrast"); }
// [[Rcpp::export(rng = false)]]
SEXP hpp_laplacian(SEXP mat, SEXP kernel) { return hpp_morpho(mat, kernel, "laplacian"); }

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
        Rcpp::stop("hpp_erode_old: when 'msk' is provided 'img' and 'msk' should have same dimensions");
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
        Rcpp::stop("hpp_dilate_old: when 'msk' is provided 'img' and 'msk' should have same dimensions");
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
