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

#ifndef IFCIP_PADDING_HPP
#define IFCIP_PADDING_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Image Padding
//' @name cpp_padding
//' @description
//' This function creates a new matrix with extra rows / cols according to input mat, kernel
//' @param mat, a NumericMatrix.
//' @param kernel, a NumericMatrix.
//' @param method, a uint8_t. Default is 1, allowed are [1-8].\cr
//' -1, extra cols / rows will be filled with 'k', returned 'out' will not be filled.\cr
//' -2, extra cols / rows will be filled with the closest col / row, returned 'out' will not be filled.\cr
//' -3, extra cols / rows will be filled mirroring neighbor cols / rows, returned 'out' will not be filled.\cr
//' -4, extra cols / rows will be filled repeating neighbor cols / rows, returned 'out' will not be filled.\cr
//' -5, extra cols / rows will be filled with 'k', returned 'out' will be filled with mat.\cr
//' -6, extra cols / rows will be filled with the closest col / row, returned 'out' will be filled with mat.\cr
//' -7, extra cols / rows will be filled mirroring neighbor cols / rows, returned 'out' will be filled with mat.\cr
//' -8, extra cols / rows will be filled repeating neighbor cols / rows, returned 'out' will be filled with mat.
//' @param k, a double, constant used when method is 1 or 4. Default is 0.0.
//' @return a List whose elements are:\cr
//' -out, a NumericMatrix, with extra cols / rows\cr
//' -ori_c, a R_len_t with x coordinate of the 1st non extra element,\cr
//' -ori_r, a R_len_t with y coordinate of the 1st non extra element.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_padding(const Rcpp::NumericMatrix mat,
                       const Rcpp::NumericMatrix kernel,
                       const uint8_t method = 1,
                       const double k = 0.0) {
  if((method < 1) || (method > 8)) {
    Rcpp::stop("hpp_padding: not allowed 'method' value");
  }
  uint8_t typ = (method > 4) ? (method - 4): method;
  // input dimension
  R_len_t mat_r, mat_c, fil_r, fil_c, ori_r, ori_c, fin_height, fin_width;
  // init var
  mat_r = mat.nrow();
  mat_c = mat.ncol();
  fil_r = kernel.nrow();
  fil_c = kernel.ncol();
  ori_r = (fil_r >> 1);
  ori_c = (fil_c >> 1);
  fin_height = mat_r + ori_r * 2;
  fin_width = mat_c + ori_c * 2;
  
  // create output matrix with expanded dimension
  Rcpp::NumericMatrix out(fin_height, fin_width);
  
  switch(typ) {
  case 1: { // fill with constant
    // 1st cols
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
    for(R_len_t i_row = 0; i_row < fin_height; i_row++) {
      out(i_row, i_col) = k;
    }
  }
    // last cols
    for(R_len_t i_col = ori_c + mat_c; i_col < fin_width; i_col++) {
      for(R_len_t i_row = 0; i_row < fin_height; i_row++) {
        out(i_row, i_col) = k;
      }
    }
    
    // 1st rows
    for(R_len_t i_col = 0; i_col < fin_width; i_col++) {
      for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
        out(i_row, i_col) = k; 
      }
    }
    
    // last rows
    for(R_len_t i_col = 0; i_col < fin_width; i_col++) {
      for(R_len_t i_row = ori_r + mat_r; i_row < fin_height; i_row++) {
        out(i_row, i_col) = k; 
      }
    }
  }
    break;
  case 2: { // fill with the closest row / col
    // add 1st extra rows
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
      out(i_row, ori_c + i_col) = mat(0, i_col);
    }
  }
    // add last extra rows
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
        out(i_row + ori_r + mat_r, ori_c + i_col) = mat(mat_r - 1, i_col);
      }
    }
    // add 1st extra col
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out(i_row + ori_r, i_col) = mat(i_row, 0);
      }
    }
    // add last extra col
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out(i_row + ori_r, i_col + ori_c + mat_c) = mat(i_row, mat_c - 1);
      }
    }
  }
    break;
  case 3: { // fill mirroring neighbor rows / cols
    // add 1st extra rows
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
      out(i_row, ori_c + i_col) = mat(ori_r - i_row - 1, i_col);
    }
  }
    // add last extra rows
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
        out(i_row + ori_r + mat_r, ori_c + i_col) = mat(mat_r - i_row - 1, i_col);
      }
    }
    // add 1st extra col
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out(i_row + ori_r, i_col) = mat(i_row, ori_c - i_col - 1);
      }
    }
    // add last extra col
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out(i_row + ori_r, i_col + ori_c + mat_c) = mat(i_row, mat_c - i_col - 1);
      }
    }
  }
    break;
  case 4: { // fill repeating neighbor rows / cols
    // add 1st extra rows
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
      out(i_row, ori_c + i_col) = mat(i_row, i_col);
    }
  } 
    // add last extra rows
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
        out(i_row + ori_r + mat_r, ori_c + i_col) = mat(mat_r - (ori_r - i_row), i_col);
      }
    }
    // add 1st extra col
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out(i_row + ori_r, i_col) = mat(i_row, i_col);
      }
    }
    // add last extra col
    for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out(i_row + ori_r, i_col + ori_c + mat_c) = mat(i_row, mat_c - (ori_c - i_col));
      }
    }
  }
    break;
  }
  
  // add corners
  double cornerUL, cornerUR, cornerLL, cornerLR;
  if(typ > 1) {
    cornerUL = mat(0, 0);
    cornerUR = mat(0, mat_c - 1);
    cornerLL = mat(mat_r - 1, 0);
    cornerLR = mat(mat_r - 1, mat_c - 1);
  } else {
    cornerUL = k;
    cornerUR = k;
    cornerLL = k;
    cornerLR = k;
  }
  for(R_len_t i_col = 0; i_col < ori_c; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
      out(i_row, i_col) = cornerUL;
    }
    for(R_len_t i_row = ori_r + mat_r; i_row < fin_height; i_row++) {
      out(i_row, i_col) = cornerLL;
    }
  }
  for(R_len_t i_col = mat_c + ori_c; i_col < fin_width; i_col++) {
    for(R_len_t i_row = 0; i_row < ori_r; i_row++) {
      out(i_row, i_col) = cornerUR;
    }
    for(R_len_t i_row = ori_r + mat_r; i_row < fin_height; i_row++) {
      out(i_row, i_col) = cornerLR;
    }
  }
  
  
  if(method > 4) {
    // copy mat into out
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      R_len_t i_out = (i_col + ori_c) * fin_height + ori_r;
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        out[i_out++] = mat(i_row, i_col);
      }
    }
  } 
  return Rcpp::List::create(Rcpp::_["out"] = out,
                            Rcpp::_["ori_c"] = ori_c,
                            Rcpp::_["ori_r"] = ori_r);
}

#endif
