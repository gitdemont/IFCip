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

#ifndef IFCIP_SHIFT_HPP
#define IFCIP_SHIFT_HPP

#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>  
Rcpp::Matrix<RTYPE> shift_row_T (Rcpp::Matrix<RTYPE> mat, 
                                 const int d_row = 0,
                                 const bool add_noise = true, 
                                 const double bg = 0.0,
                                 const double sd = 0.0,
                                 const bool keep_size = false) {
  if(d_row == 0) return mat;
  R_len_t M_row = mat.nrow();
  R_len_t M_col = mat.ncol();
  R_len_t i_row;
  R_len_t new_row = M_row + std::abs(d_row);
  
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init_matrix(new_row, M_col);
  if(d_row < 0) {
    if(add_noise) {
      for(i_row = 0; i_row < -d_row; i_row++) out(i_row, Rcpp::_) = Rcpp::rnorm(M_col, bg, sd);
    } else {
      for(i_row = 0; i_row < -d_row; i_row++) out(i_row, Rcpp::_) = Rcpp::Vector<RTYPE>(M_col, bg);
    }
    for(; i_row < new_row; i_row++) out(i_row, Rcpp::_) = mat(i_row + d_row, Rcpp::_);
  } else {
    for(i_row = 0; i_row < M_row; i_row++) out(i_row, Rcpp::_) = mat(i_row, Rcpp::_);
    if(add_noise) {
      for(; i_row < new_row; i_row++) out(i_row, Rcpp::_) = Rcpp::rnorm(M_col, bg, sd);
    } else {
      for(; i_row < new_row; i_row++) out(i_row, Rcpp::_) = Rcpp::Vector<RTYPE>(M_col, bg);
    }
  }
  if(keep_size) return out(Range(std::max(0, d_row), std::max(M_row - 1, M_row - 1 + d_row)), Rcpp::_ );
  return out;
}

template <int RTYPE>  
Rcpp::Matrix<RTYPE> shift_col_T (Rcpp::Matrix<RTYPE> mat, 
                                 const int d_col = 0,
                                 const bool add_noise = true, 
                                 const double bg = 0.0,
                                 const double sd = 0.0,
                                 const bool keep_size = false) {
  if(d_col == 0) return mat;
  R_len_t M_row = mat.nrow();
  R_len_t M_col = mat.ncol();
  R_len_t i_col;
  R_len_t new_col = M_col + std::abs(d_col);
  
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init_matrix(M_row, new_col);
  if(d_col < 0) {
    if(add_noise) {
      for(i_col = 0; i_col < -d_col; i_col++) out(Rcpp::_, i_col) = Rcpp::rnorm(M_row, bg, sd);
    } else {
      for(i_col = 0; i_col < -d_col; i_col++) out(Rcpp::_, i_col) = Rcpp::Vector<RTYPE>(M_row, bg);
    }
    for(; i_col < new_col; i_col++) out(Rcpp::_,i_col) = mat(Rcpp::_, i_col + d_col);
  } else {
    for(i_col = 0; i_col < M_col; i_col++) out(Rcpp::_, i_col) = mat(Rcpp::_, i_col);
    if(add_noise) {
      for(; i_col < new_col; i_col++) out(Rcpp::_, i_col) = Rcpp::rnorm(M_row, bg, sd);
    } else {
      for(; i_col < new_col; i_col++) out(Rcpp::_, i_col) = Rcpp::Vector<RTYPE>(M_row, bg);
    }
  }
  if(keep_size) return out(Rcpp::_, Range(std::max(0, d_col), std::max(M_col - 1, M_col - 1 + d_col)));
  return out;
}

//' @title Image Row Shift
//' @name cpp_shift_row
//' @description Function to shift matrix's rows.
//' @param mat a numeric matrix.
//' @param d_row an integer, giving row shift. Default is \code{0} for no change.
//' @param add_noise logical, if \code{true} adds normal noise when at least one new dimension is larger than original mat dimensions
//' Rcpp::rnorm() function is used. Default is \code{true}.
//' @param bg double, mean value of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param sd double, standard deviation of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param keep_size, a bool whether initial \code{'mat'} dimension should be kept. Default is \code{false}.
//' @return a shifted matrix with additional rows if \code{'d_row'} is different from \code{0}.
//' @keywords internal
////' @export
// [[Rcpp::export]]
SEXP hpp_shift_row (SEXP mat,
                    const int d_col = 0,
                    const bool add_noise = true, 
                    const double bg = 0.0,
                    const double sd = 0.0,
                    const bool keep_size = false) {
  switch(TYPEOF(mat)) {
  case NILSXP : return mat;
  case RAWSXP : return shift_row_T<RAWSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  case LGLSXP : return shift_row_T<LGLSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  case INTSXP : return shift_row_T<INTSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  case REALSXP : return shift_row_T<REALSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  default: Rcpp::stop("hpp_shift_row: not supported type in 'mat'");
  }
}

//' @title Image Column Shift
//' @name cpp_shift_col
//' @description Function to shift matrix's columns.
//' @param mat a numeric matrix.
//' @param d_col an integer, giving col shift. Default is \code{0} for no change.
//' @param add_noise logical, if \code{true} adds normal noise when at least one new dimension is larger than original mat dimensions
//' Rcpp::rnorm() function is used. Default is \code{true}.
//' @param bg double, mean value of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param sd double, standard deviation of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param keep_size, a bool whether initial \code{'mat'} dimension should be kept. Default is \code{false}.
//' @return a shifted matrix with additional columns if \code{'d_col'} is different from \code{0}.
//' @keywords internal
////' @export
// [[Rcpp::export]]
SEXP hpp_shift_col (SEXP mat,
                    const int d_col = 0,
                    const bool add_noise = true, 
                    const double bg = 0.0,
                    const double sd = 0.0,
                    const bool keep_size = false) {
  switch(TYPEOF(mat)) {
  case NILSXP : return mat;
  case RAWSXP : return shift_col_T<RAWSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  case LGLSXP : return shift_col_T<LGLSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  case INTSXP : return shift_col_T<INTSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  case REALSXP : return shift_col_T<REALSXP>(mat, d_col, add_noise, bg, sd, keep_size);
  default: Rcpp::stop("hpp_shift_col: not supported type in 'mat'");
  }
}

#endif
