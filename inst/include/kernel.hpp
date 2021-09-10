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

#ifndef IFCIP_KERNEL_HPP
#define IFCIP_KERNEL_HPP

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_make_disc(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double foo = i_col - half;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      out(i_row, i_col) = std::sqrt(foo * foo + pow(i_row - half, 2.0)) <= half;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_make_box(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  out.fill(true);
  return out;
}

// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_make_plus(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  out(half, Rcpp::_) = Rcpp::rep(true, size);
  out(Rcpp::_, half) = Rcpp::rep(true, size);
  return out;
}

// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_make_cross(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      out(i_row, i_col) = i_row == i_col || i_row == (size - 1 - i_col);
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_make_diamond(uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  double half = size >> 1;
  R_len_t i = 0;
  if(size % 2) {
    for(R_len_t i_col = -half; i_col <= half; i_col++) {
      for(R_len_t i_row = -half; i_row <= half; i_row++) {
        if(i >= size * size) Rcpp::stop("Not allowed");
        out[i++] = (abs(i_col) + abs(i_row)) <= half;
      }
    }
    
  } else {
    for(R_len_t i_col = -half; i_col <= half; i_col++) {
      if(i_col == 0) i_col++;
      for(R_len_t i_row = -half; i_row <= half; i_row++) {
        if(i_row == 0) i_row++;
        if(i >= size * size) Rcpp::stop("Not allowed");
        out[i++] = (abs(i_col) + abs(i_row)) <= half;
      }
    }
  }
  return out;
}
#endif
