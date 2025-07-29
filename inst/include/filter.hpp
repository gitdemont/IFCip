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
#include "kernel.hpp"
#include "padding.hpp"
using namespace Rcpp;

double filter_fun(const Rcpp::NumericVector x, int what_n, std::string what) {
  switch(what_n) {
  case 1: return x.size() == 0 ? NA_REAL : Rcpp::sum(x);
  case 2: return x.size() == 0 ? NA_REAL : Rcpp::sum(x);
  case 3: return x.size() == 0 ? NA_REAL : Rcpp::sum(x);
  case 4: return x.size() == 0 ? NA_REAL : Rcpp::mean(x);
  case 5: return x.size() == 0 ? NA_REAL : Rcpp::median(x);
  case 6: return x.size() == 0 ? NA_REAL : Rcpp::sd(x);
  case 7: return x.size() == 0 ? NA_REAL : (Rcpp::max(x) - Rcpp::min(x)) / 2;
  case 8: return x.size() == 0 ? NA_REAL : 3 * Rcpp::median(x) - 2 * Rcpp::mean(x);
  case 9: return x.size() == 0 ? NA_REAL : Rcpp::min(x);
  case 10: return x.size() == 0 ? NA_REAL : Rcpp::max(x);
  default: Rcpp::stop("hpp_filter: filter function['%s']is not supported", what);
  }
}

template <int RTYPE, int RTYPEk>
Rcpp::NumericMatrix filter_T(Rcpp::Matrix<RTYPE> mat,
                             Rcpp::Matrix<RTYPEk> kernel,
                             const uint8_t method = 5,
                             const double k = NA_REAL,
                             const std::string what = "") {
  R_len_t mat_r = mat.nrow();  
  R_len_t mat_c = mat.ncol();
  R_len_t ks = kernel.size();
  if(mat_r == 0 || mat_c == 0 || ks == 0) return Rcpp::as<Rcpp::NumericMatrix>(mat);
  
  // get kernel dimension 
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc; //APPLY OFFSET CORRECTION
  if(what != "correlate") {
    if(kr) kr = !kr;          //APPLY OFFSET CORRECTION
    if(kc) kc = !kc;          //APPLY OFFSET CORRECTION
  }
  
  // match filter function
  int fun = Rcpp::match(Rcpp::CharacterVector::create(what), Rcpp::CharacterVector::create("correlate", "convolve", "convolve0", "mean", "median", "sd", "mid", "mode", "min", "max"))[0];
  
  // initialize out that will be padded using method and k
  // 0.0 is important to allow removal of extra values from final result (multiply by 0.0) for correlate and convolve
  // NA is important to allow removal of extra values from final result (na_omit) for filtering
  // R_NegInf/R_PosInf eventually for max/min 
  Rcpp::NumericMatrix out, foo;
  if(what == "correlate") {
    foo = padding_T(mat, pad_r + kc, pad_c + kr, method, k);
    out = Rcpp::NumericMatrix(foo.nrow(), foo.ncol());
    for(R_len_t i_col = out.ncol() - pad_c - 1; i_col >= pad_c; i_col--) {
      for(R_len_t i_row = out.nrow() - pad_r - 1; i_row >= pad_r; i_row--) {
        for(R_len_t f_col = i_col + pad_c_1 - 1, i_ker = kernel.size() - 1; f_col >= i_col - pad_c; f_col--) {
          for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--, i_ker--) {
            out(i_row, i_col) = kernel[i_ker] * foo(f_row, f_col) + out(i_row, i_col);
          }
        }
      }
    }
  } else {
    if(what == "convolve") {
      foo = padding_T(mat, pad_r + kc, pad_c + kr, method, k);
      out = Rcpp::NumericMatrix(foo.nrow(), foo.ncol());
      for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
        for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
          for(R_len_t f_col = i_col + pad_c_1 - 1, i_ker = 0; f_col >= i_col - pad_c; f_col--) {
            for(R_len_t f_row = i_row + pad_r_1 - 1; f_row >= i_row - pad_r; f_row--, i_ker++) {
              out(i_row, i_col) = kernel[i_ker] * foo(f_row, f_col) + out(i_row, i_col);
            }
          }
        }
      }
    } else { // here we apply filtering
      foo = padding_T(mat, pad_r + kc, pad_c + kr, method, k);
      out = Rcpp::clone(foo);
      // matrix of offsets in reverse raster order to allow to improve speed with when there are FALSE values in kernel
      Rcpp::IntegerMatrix o = offset_kernel(kernel, true, true, true);
      R_len_t ker_c = o.ncol(); 
      Rcpp::NumericVector K(ker_c, NA_REAL);
      
      for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
        for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
          for(R_len_t i_k = 0; i_k < ker_c ; i_k++) {
            R_len_t m_col = o[3 * i_k] + i_col, m_row = o[3 * i_k + 1] + i_row;
            K[i_k] = foo(m_row, m_col);
          }
          // TODO use XPtr instead of fun switch
          out(i_row, i_col) = filter_fun(Rcpp::na_omit(K), fun, what);
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Image Filtering
//' @name cpp_filter
//' @description
//' This function applies filtering on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param method used for padding, a uint8_t. Default is 5, allowed are [1-8].
//' @param k, constant used for padding, a double. Default is NA_REAL.
//' @param what, type of filtering, s std::string. Default is "".
//' @return a NumericMatrix.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_filter (SEXP mat,
                                SEXP kernel,
                                const uint8_t method = 5,
                                const double k = NA_REAL,
                                const std::string what = "") {
  switch( TYPEOF(mat) ) {
  // case NILSXP : return R_NilValue; // prefer throwing error on NILSXP
  case LGLSXP : return filter_T(Rcpp::as<NumericMatrix>(mat), get_kernel(kernel), method, k, what);
  case RAWSXP : return filter_T(Rcpp::as<NumericMatrix>(mat), get_kernel(kernel), method, k, what);
  case INTSXP : return filter_T(Rcpp::as<NumericMatrix>(mat), get_kernel(kernel), method, k, what);
  case REALSXP : return filter_T(Rcpp::as<NumericMatrix>(mat), get_kernel(kernel), method, k, what);
  default : Rcpp::stop("hpp_filter: not supported SEXP[%i] in 'mat'", TYPEOF(mat));
  }
}

# endif
