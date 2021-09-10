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

#ifndef IFCIP_MATRIX_LOGIC_HPP
#define IFCIP_MATRIX_LOGIC_HPP

#include <Rcpp.h>
// [[Rcpp::depends(IFC)]]
#include <matrix_logic.hpp>
using namespace Rcpp;

//' @title Matrix Logic Equal
//' @name cpp_k_equal_M
//' @description
//' This function takes an NumericMatrix and checks for members equal to k
//' @param mat a NumericMatrix
//' @param k constant to be checked. Default is 3.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_k_equal_M(const Rcpp::NumericMatrix mat, const double k = 3.0 ) {
  R_len_t mat_r = mat.nrow();
  Rcpp::LogicalMatrix OUT_M(mat_r, mat.ncol());
  for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
    OUT_M(i_row, Rcpp::_) = mat(i_row, Rcpp::_) == k;
  }
  return OUT_M;
}

//' @title Matrix Logic Superior and Equal
//' @name cpp_k_sup_equal_M
//' @description
//' This function takes an NumericMatrix and checks for members superior or equal to k
//' @param mat a NumericMatrix
//' @param k constant to be checked. Default is 3.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_k_sup_equal_M(const Rcpp::NumericMatrix mat, const double k = 3.0) {
  R_len_t mat_r = mat.nrow();
  Rcpp::LogicalMatrix OUT_M(mat_r, mat.ncol());
  for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
    OUT_M(i_row, Rcpp::_) = mat(i_row, Rcpp::_) >= k;
  }
  return OUT_M;
}

//' @title Matrix Logic Inferior and Equal
//' @name cpp_k_inf_equal_M
//' @description
//' This function takes an NumericMatrix and checks for members inferior or equal to k
//' @param mat a NumericMatrix
//' @param k constant to be checked. Default is 3.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_k_inf_equal_M(const Rcpp::NumericMatrix mat, const double k = 3.0) {
  R_len_t mat_r = mat.nrow();
  Rcpp::LogicalMatrix OUT_M(mat_r, mat.ncol());
  for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
    OUT_M(i_row, Rcpp::_) = mat(i_row, Rcpp::_) <= k;
  }
  return OUT_M;
}

#endif
