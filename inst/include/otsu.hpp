/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2022 Yohann Demont                                              

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
  
#ifndef IFCIP_OTSU_HPP
#define IFCIP_OTSU_HPP

#include <Rcpp.h>
#include "scale.hpp"
using namespace Rcpp;

//' @title Otsu Multi Thresholding
//' @name cpp_multi_otsu
//' @description
//' This function determines best threshold(s) according to Otsu's method.
//' @param img, a NumericMatrix.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @param n_comp, number of components to separate. Default is 2, should be at least 2.\cr
//' Returned thresholds will be of length n_comp - 1.
//' @param n_lev, an unsigned short determining the number of grey levels used for the computation. Default is 256, should be at least 2.
//' Despite being fast thanks to LUT pre-computation, performance will be highly impacted with large 'n_comp' or 'n_lev' values (typically n_comp = 5 and n_lev = 256).
//' Alternatively, you can try to decrease 'n_lev' when 'n_comp' needs to be large (e.g. n_comp = 8 and n_lev = 32).
//' @details adaptation of 'A Fast Algorithm for Multilevel Thresholding' from L. Ping-Sung, C. Tse-Sheng, and C. Pau-Choo
//' in Jounal of Information Science and Engineering. 2001(17), 713-727.
//' \doi{10.6688/JISE.2001.17.5.1}
//' @return an NumericVector of threshold(s).
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_multi_otsu (const Rcpp::NumericMatrix img,
                                    const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                                    const uint8_t n_comp = 2,
                                    const unsigned short n_lev = 256) {
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  R_len_t MAX_SIZ = mat_r * mat_c;
  if(MAX_SIZ <= 0) Rcpp::stop("'img' should be at least 1px row and 1px col");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 2)) Rcpp::stop("'img' is too large");
  if(n_lev < 2) Rcpp::stop("'n_lev' should be at least >= 2");
  int MAX_LEV = n_lev;
  int MAX_LEV_1 = MAX_LEV - 1;
  int MAX_LEV_2 = MAX_LEV - 2;
  if((n_comp < 2) || (n_comp > MAX_LEV - 2)) Rcpp::stop("'n_comp' value is not allowed");
  
  // // determines image range
  // double img_min = R_PosInf, img_max = R_NegInf;
  // for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
  //   Rcpp::NumericVector col_ran = range(img(Rcpp::_, i_col));
  //   if(col_ran[0] < img_min) {
  //     img_min = col_ran[0];
  //   } else {
  //     if(col_ran[1] > img_max) img_max = col_ran[1];
  //   }
  // }
  // 
  // // creates scaled mat ranging from [0, MAX_LEV - 1]
  // Rcpp::IntegerMatrix sca = Rcpp::no_init(mat_r, mat_c);
  // double MAX_LEV_SCA = MAX_LEV_1 / (img_max - img_min);
  // for(R_len_t i = 0; i < MAX_SIZ; i++) sca[i] = MAX_LEV_SCA * (img[i] - img_min);
  
  Rcpp::NumericMatrix mat = Rcpp::clone(img);
  Rcpp::NumericVector Q = hpp_scale(mat, msk_, -1.0, n_lev, true);
  Rcpp::IntegerMatrix sca = as<Rcpp::IntegerMatrix>(mat);
  
  // fill bins
  // in the article LUT start at 1 and end at L, here there is an extra [0] = 0.0 to handle this
  Rcpp::NumericVector B(MAX_LEV + 1);
  for (int i=0; i < MAX_SIZ; ++i) if(sca[i] >= 0) B[sca[i] + 1]++;
  
  // compute LUT
  // in the article LUT start at 1 and end at L, here there is an extra row0 and col0 with value 0.0 to handle this
  Rcpp::NumericMatrix P(MAX_LEV + 1, MAX_LEV + 1); // P(1,0) = 0.0
  Rcpp::NumericMatrix S(MAX_LEV + 1, MAX_LEV + 1); // S(1,0) = 0.0
  Rcpp::NumericMatrix H = no_init(MAX_LEV + 1, MAX_LEV + 1);      // H(0,_) and H(_,0) are never used
  for (int v = 0; v < MAX_LEV; v++) {
    int vv = v + 1;
    P(1, vv) = P(1, v) + B[vv];                                   // eq 22, 0th-order moment
    S(1, vv) = S(1, v) + vv * B[vv];                              // eq 23, 1st-order moment
  }
  for (int u = 1; u <= MAX_LEV; u++) {
    int uu = u - 1;
    for (int v = 1; v <= MAX_LEV; v++) {
      P(u, v) = P(1, v) - P(1, uu);                               // eq 24
      S(u, v) = S(1, v) - S(1, uu);                               // eq 25
      H(u, v) = P(u, v) == 0 ? 0.0 : S(u, v) * S(u, v) / P(u, v); // eq 29
    }
  }
  
  // initialize threshold idx vector
  uint8_t n_comp_1 = n_comp - 1;
  uint8_t n_comp_2 = n_comp - 2;
  Rcpp::IntegerVector curr_tr(n_comp_1);
  Rcpp::IntegerVector best_tr(n_comp_1);
  curr_tr[0] = 1;                                                                 // In the article indices start at 1
  for(uint8_t i = 1; i < n_comp_1; i++) curr_tr[i] = curr_tr[i - 1] + 1;          // each additional threshold is +1 from previous one
  
  // compute sigma for each combination of threshold(s)
  // and store combination with max sigma
  double smax = 0.0;
  do {
    double s = 0.0;                                                               // eq 28
    s += H(1, curr_tr[0]) + H(curr_tr[n_comp_2] + 1, MAX_LEV_1);                  // H(1,t1) + ... + H(tm-1, L)
    for(uint8_t i = 0; i < n_comp_2; i++) s += H(curr_tr[i] + 1, curr_tr[i + 1]); // all the remaining ...
    if(smax < s) {
      std::copy(curr_tr.begin(), curr_tr.end(), best_tr.begin());
      smax = s;
    }
    // compute all possible combinations of threshold(s)
    for(uint8_t n = n_comp_1; n > 0;) {
      n--;
      if(curr_tr[n] < (MAX_LEV_2 - n_comp_1 + n)) {
        curr_tr[n]++;
        for(uint8_t k = n; k < n_comp_2; k++) curr_tr[k + 1] = curr_tr[k] + 1;
        break;
      }
    }
  } while (curr_tr[0] < MAX_LEV_2 - n_comp_1);
  // one last to compute and test
  double s = 0.0;                                                                 // eq 28
  s += H(1, curr_tr[0]) + H(curr_tr[n_comp_2] + 1, MAX_LEV_1);                    // H(1,t1) + ... + H(tm-1, L)
  for(uint8_t i = 0; i < n_comp_2; i++) s += H(curr_tr[i] + 1, curr_tr[i + 1]);   // all the remaining  ...
  if(smax < s) {
    std::copy(curr_tr.begin(), curr_tr.end(), best_tr.begin());
    smax = s;
  }
  
  // returned threshold(s) need to be rescaled back to initial image value
  Rcpp::NumericVector out = no_init(n_comp_1);
  for(uint8_t i = 0; i < n_comp_1; i++) out[i] = Q[1] - best_tr[n_comp_1 - i - 1] / Q[2];
  return out;
}

#endif
