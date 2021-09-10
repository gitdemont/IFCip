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

#ifndef IFCIP_FLIP_HPP
#define IFCIP_FLIP_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Image Horizontal Flip
//' @name cpp_hflip
//' @description
//' Function that flips mat horizontally
//' @param mat NumericMatrix.
//' @return a flipped matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_hflip(const Rcpp::NumericMatrix mat) {
  // R_len_t img_r = mat.nrow();
  R_len_t img_c = mat.ncol();
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(mat.nrow(), img_c);
  for(R_len_t i_col = img_c -1; i_col >= 0; i_col--) {
    out(Rcpp::_, i_col) = mat(Rcpp::_, img_c - i_col - 1);
  }
  return out;
}

//' @title Image Vertical Flip
//' @name cpp_vflip
//' @description
//' Function that flips mat vertically
//' @param mat NumericMatrix.
//' @return a flipped matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_vflip(const Rcpp::NumericMatrix mat) {
  R_len_t img_r = mat.nrow();
  // R_len_t img_c = mat.ncol();
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(img_r, mat.ncol());
  for(R_len_t i_row = img_r - 1; i_row >= 0 ; i_row--) {
    out(i_row, Rcpp::_) = mat(img_r - i_row - 1, Rcpp::_);
  }
  return out;
}

#endif
