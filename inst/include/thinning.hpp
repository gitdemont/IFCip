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

#ifndef IFCIP_THINNING_HPP
#define IFCIP_THINNING_HPP

#include <Rcpp.h>
using namespace Rcpp;

static int ifcip_thinning_dx[8]={ 0, 1, 1, 1, 0,-1,-1,-1};
static int ifcip_thinning_dy[8]={-1,-1, 0, 1, 1, 1, 0,-1};

//' @title Implementation of Zheng-Suen Thinning
//' @name cpp_thinning_zs
//' @description
//' This function is designed to identify mask skeleton.
//' @param mat a LogicalMatrix, containing mask.
//' @details adaptation of 'A fast parallel algorithm for thinning digital patterns' from T. Y. Zhang, C. Y. Suen.
//' Communications of the ACM, March 1984 \url{https://doi.org/10.1145/357994.358023}.
//' @return a LogicalMatrix with the mask thinned.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_thinning_zs(const Rcpp::LogicalMatrix mat) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  R_len_t i_del = 0, del1 = 0, del2 = 0, i_out = 0, i_row = 0, i_col = 0;
  unsigned short j, N, T;
  
  // create output matrix with extra cols / rows
  Rcpp::LogicalMatrix out(mat_r + 2, mat_c + 2);
  // add padding
  out(0, _) = Rcpp::LogicalMatrix(mat_c + 2, false);
  out(mat_r + 1, _) = Rcpp::LogicalMatrix(mat_c + 2, false);
  out(_, 0) = Rcpp::LogicalMatrix(mat_r + 2, false);
  out(_, mat_c + 1) = Rcpp::LogicalMatrix(mat_r + 2, false);
  for(i_col = 0; i_col < mat_c; i_col++) {
    i_out = (i_col + 1) * (mat_r + 2) + 1;
    for(i_row = 0; i_row < mat_r; i_row++) {
      out[i_out++] = mat(i_row, i_col);
    }
  }
  Rcpp::IntegerVector removal_col((mat_r + 2) * (mat_c + 2));
  Rcpp::IntegerVector removal_row((mat_r + 2) * (mat_c + 2));
  Rcpp::LogicalVector nbr(8);
  
  do {
    // Rcpp::checkUserInterrupt();
    
    // step 1
    del1 = 0;
    for(i_col = 1; i_col <= mat_c; i_col++) {
      // Rcpp::checkUserInterrupt();
      for(i_row = 1; i_row <= mat_r; i_row++) {
        // Rcpp::checkUserInterrupt();
        if(out(i_row, i_col)) {
          N = 0; // number of 8-connected neighbors with non zero value
          T = 0; // number of 0-1 transitions in the ordered sequence
          nbr[0] = out(i_row + ifcip_thinning_dy[0], i_col + ifcip_thinning_dx[0]);
          for(j = 1; j < 8; j++) {
            // Rcpp::checkUserInterrupt();
            nbr[j] = out(i_row + ifcip_thinning_dy[j], i_col + ifcip_thinning_dx[j]);
            if(nbr[j]) {
              N++;
              if(!nbr[j-1]) T++;
            }
          }
          if(nbr[0]) {
            N++;
            if(!nbr[7]) T++;
          } 
          if((N >= 2) && (N <= 6)                  // condition a
               && (T == 1)                         // condition b
               && !(nbr[0] && nbr[2] && nbr[4])    // condition c
               && !(nbr[2] && nbr[4] && nbr[6])) { // condition d
               removal_col[del1] = i_col;
            removal_row[del1++] = i_row;
          }
        }
      }
    }
    
    for(i_del = 0; i_del < del1; i_del++) {
      out(removal_row[i_del], removal_col[i_del]) = 0;
    }
    
    // step 2
    del2 = 0;
    for(i_col = 1; i_col <= mat_c; i_col++) {
      // Rcpp::checkUserInterrupt();
      for(i_row = 1; i_row <= mat_r; i_row++) {
        // Rcpp::checkUserInterrupt();
        if(out(i_row, i_col)) {
          N = 0; // number of 8-connected neighbors with non zero value
          T = 0; // number of 0-1 transitions in the ordered sequence
          nbr[0] = out(i_row + ifcip_thinning_dy[0], i_col + ifcip_thinning_dx[0]);
          for(j = 1; j < 8; j++) {
            // Rcpp::checkUserInterrupt();
            nbr[j] = out(i_row + ifcip_thinning_dy[j], i_col + ifcip_thinning_dx[j]);
            if(nbr[j]) {
              N++;
              if(!nbr[j-1]) T++;
            }
          }
          if(nbr[0]) {
            N++;
            if(!nbr[7]) T++;
          } 
          if((N >= 2) && (N <= 6)                  // condition a
               && (T == 1)                         // condition b
               && !(nbr[0] && nbr[2] && nbr[6])    // condition c
               && !(nbr[0] && nbr[4] && nbr[6])) { // condition d
               removal_col[del2] = i_col;
            removal_row[del2++] = i_row;
          }
        }
      }
    }
    
    for(i_del = 0; i_del < del2; i_del++) {
      out(removal_row[i_del], removal_col[i_del]) = 0;
    }
  } while(del1 + del2);
  
  return out(Rcpp::Range(1, mat_r), Rcpp::Range(1, mat_c));
}

//' @title Implementation of Ben Boudaoud-Sider-Tari Thinning
//' @name cpp_thinning_bst
//' @description
//' This function is designed to identify mask skeleton.
//' @param mat a LogicalMatrix, containing mask.
//' @details adaptation of 'A new thinning algorithm for binary images' from L. Ben Boudaoud, A. Sider, A. Tari.
//' 3rd international conference on control, engineering & information technology, May 2015. \url{https://doi.org/10.1109/CEIT.2015.7233099}.
//' @return a LogicalMatrix with the mask thinned.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_thinning_bst(const Rcpp::LogicalMatrix mat) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  R_len_t i_del = 0, del1 = 0, del2 = 0, i_out = 0, i_row = 0, i_col = 0;
  unsigned short j, N, T;
  
  // create output matrix with extra cols / rows
  Rcpp::LogicalMatrix out(mat_r + 2, mat_c + 2);
  // add padding
  out(0, _) = Rcpp::LogicalMatrix(mat_c + 2, false);
  out(mat_r + 1, _) = Rcpp::LogicalMatrix(mat_c + 2, false);
  out(_, 0) = Rcpp::LogicalMatrix(mat_r + 2, false);
  out(_, mat_c + 1) = Rcpp::LogicalMatrix(mat_r + 2, false);
  for(i_col = 0; i_col < mat_c; i_col++) {
    i_out = (i_col + 1) * (mat_r + 2) + 1;
    for(i_row = 0; i_row < mat_r; i_row++) {
      out[i_out++] = mat(i_row, i_col);
    }
  }
  Rcpp::IntegerVector removal_col((mat_r + 2) * (mat_c + 2));
  Rcpp::IntegerVector removal_row((mat_r + 2) * (mat_c + 2));              
  Rcpp::LogicalVector nbr(8);
  
  do {
    // Rcpp::checkUserInterrupt();
    
    // step 1
    del1 = 0;
    for(i_col = 1; i_col <= mat_c; i_col++) {
      // Rcpp::checkUserInterrupt();
      for(i_row = 1; i_row <= mat_r; i_row++) {
        // Rcpp::checkUserInterrupt();
        if((((i_col + i_row) % 2) == 0)            // condition a
             && out(i_row, i_col)) {
          N = 0; // number of 8-connected neighbors with non zero value
          T = 0; // number of 8-connected components of "1" in the 8-connected neighborhood
          for(j = 0; j < 8; j++) {
            // Rcpp::checkUserInterrupt();
            nbr[j] = out(i_row + ifcip_thinning_dy[j], i_col + ifcip_thinning_dx[j]);
            if(nbr[j]) N++;
          }
          T = (!nbr[0] && (nbr[1] || nbr[2]))
            + (!nbr[2] && (nbr[3] || nbr[4]))
            + (!nbr[4] && (nbr[5] || nbr[6]))
            + (!nbr[6] && (nbr[7] || nbr[0]));
          if((N >= 2) && (N <= 7)                  // condition b
               && (T == 1)                         // condition c
               && !(nbr[0] && nbr[2] && nbr[4])    // condition d
               && !(nbr[2] && nbr[4] && nbr[6])) { // condition e
                removal_col[del1] = i_col;
                removal_row[del1++] = i_row;
          }
        }
      }
    }
    
    for(i_del = 0; i_del < del1; i_del++) {
      out(removal_row[i_del], removal_col[i_del]) = 0;
    }
    
    // step 2
    del2 = 0;
    for(i_col = 1; i_col <= mat_c; i_col++) {
      // Rcpp::checkUserInterrupt();
      for(i_row = 1; i_row <= mat_r; i_row++) {
        // Rcpp::checkUserInterrupt();
        if((((i_col + i_row) % 2) != 0)            // condition a
             && out(i_row, i_col)) {
          N = 0; // number of 8-connected neighbors with non zero value
          T = 0; // number of 8-connected components of "1" in the 8-connected neighborhood
          for(j = 0; j < 8; j++) {
            // Rcpp::checkUserInterrupt();
            nbr[j] = out(i_row + ifcip_thinning_dy[j], i_col + ifcip_thinning_dx[j]);
            if(nbr[j]) N++;
          }
          T = (!nbr[0] && (nbr[1] || nbr[2]))
            + (!nbr[2] && (nbr[3] || nbr[4]))
            + (!nbr[4] && (nbr[5] || nbr[6]))
            + (!nbr[6] && (nbr[7] || nbr[0]));
          if((N >= 2) && (N <= 7)                  // condition b
               && (T == 1)                         // condition c
               && !(nbr[0] && nbr[2] && nbr[6])    // condition d
               && !(nbr[0] && nbr[4] && nbr[6])) { // condition e
                removal_col[del2] = i_col;
                removal_row[del2++] = i_row;
          }
        }
      }
    }
    
    for(i_del = 0; i_del < del2; i_del++) {
      out(removal_row[i_del], removal_col[i_del]) = 0;
    }
  } while(del1 + del2);
  
  return out(Rcpp::Range(1, mat_r), Rcpp::Range(1, mat_c));
}

#endif
