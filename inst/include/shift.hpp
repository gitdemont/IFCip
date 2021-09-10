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

//' @title Image Row Shift
//' @name cpp_shift_row
//' @description Function to shift matrix's rows.
//' @param mat a numeric matrix.
//' @param d_row an integer, giving row shift. Default is 0 for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions
//' Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @return a shifted matrix with additional rows/columns if d_row or d_col are different from 0.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_shift_row( const Rcpp::NumericMatrix mat, 
                                   const int d_row = 0,
                                   const bool add_noise = true, 
                                   const double bg = 0.0,
                                   const double sd = 0.0) {
  if(d_row == 0) return mat;
  R_len_t M_row = mat.nrow();
  R_len_t M_col = mat.ncol();
  R_len_t i_row;
  R_len_t new_row = M_row + std::abs(d_row);
  
  Rcpp::NumericMatrix out(new_row, M_col);
  if(d_row < 0) {
    if(add_noise) {
      for(i_row = 0; i_row < -d_row; i_row++) out(i_row, Rcpp::_) = Rcpp::rnorm(M_col, bg, sd);
    } else {
      for(i_row = 0; i_row < -d_row; i_row++) out(i_row, Rcpp::_) = Rcpp::NumericVector(M_col, bg);
    }
    for(; i_row < new_row; i_row++) out(i_row, Rcpp::_) = mat(i_row + d_row, Rcpp::_);
  } else {
    for(i_row = 0; i_row < M_row; i_row++) out(i_row, Rcpp::_) = mat(i_row, Rcpp::_);
    if(add_noise) {
      for(; i_row < new_row; i_row++) out(i_row, Rcpp::_) = Rcpp::rnorm(M_col, bg, sd);
    } else {
      for(; i_row < new_row; i_row++) out(i_row, Rcpp::_) = Rcpp::NumericVector(M_col, bg);
    }
  }
  return out;
}

//' @title Image Column Shift
//' @name cpp_shift_col
//' @description Function to shift matrix's columns.
//' @param mat a numeric matrix.
//' @param d_col an integer, giving col shift. Default is 0 for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions
//' Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @return a shifted matrix with additional rows/columns if d_row or d_col are different from 0.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_shift_col( const Rcpp::NumericMatrix mat, 
                                   const int d_col = 0,
                                   const bool add_noise = true, 
                                   const double bg = 0.0,
                                   const double sd = 0.0) {
  if(d_col == 0) return mat;
  R_len_t M_row = mat.nrow();
  R_len_t M_col = mat.ncol();
  R_len_t i_col;
  R_len_t new_col = M_col + std::abs(d_col);
  
  Rcpp::NumericMatrix out(M_row, new_col);
  if(d_col < 0) {
    if(add_noise) {
      for(i_col = 0; i_col < -d_col; i_col++) out(Rcpp::_, i_col) = Rcpp::rnorm(M_row, bg, sd);
    } else {
      for(i_col = 0; i_col < -d_col; i_col++) out(Rcpp::_, i_col) = Rcpp::NumericVector(M_row, bg);
    }
    for(; i_col < new_col; i_col++) out(Rcpp::_,i_col) = mat(Rcpp::_, i_col + d_col);
  } else {
    for(i_col = 0; i_col < M_col; i_col++) out(Rcpp::_, i_col) = mat(Rcpp::_, i_col);
    if(add_noise) {
      for(; i_col < new_col; i_col++) out(Rcpp::_, i_col) = Rcpp::rnorm(M_row, bg, sd);
    } else {
      for(; i_col < new_col; i_col++) out(Rcpp::_, i_col) = Rcpp::NumericVector(M_row, bg);
    }
  }
  return out;
}

#endif
