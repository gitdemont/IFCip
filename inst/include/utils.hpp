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

#ifndef IFCIP_UTILS_HPP
#define IFCIP_UTILS_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Image Background
//' @name cpp_background
//' @description
//' This function is designed to compute image background on a "raw" image
//' @param img a NumericMatrix, containing image intensity values.
//' @param margin R_len_t number of rows margin used to compute background. Default is 4.
//' @param extra R_len_t number of extra columns used to compute background. Default is 0.
//' @param is_cif a bool whether 'ímg' originates from a cif_file or not. Default is false.
//' @return a NumericVector of background mean and sd
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_background(const Rcpp::NumericMatrix img,
                                   const R_len_t margin = 4,
                                   const R_len_t extra = 0,
                                   const bool is_cif = false) {
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  if(mat_r <= 2 * margin) {
    Rcpp::Rcerr << "cpp_background: 'margin' is too high for 'img'" << std::endl;
    Rcpp::stop("cpp_background: 'margin' is too high for 'img'");
  }
  if(mat_c <= 2 * extra) {
    Rcpp::Rcerr << "cpp_background: 'extra' is too high for 'img'" << std::endl;
    Rcpp::stop("cpp_background: 'extra' is too high for 'img'");
  }
  // compute mean
  R_len_t count_b = 0;
  double bkg_mean = 0.0;
  for(R_len_t i_col = extra; i_col < mat_c - extra; i_col++) { 
    for(R_len_t i_row = 0; i_row < margin; i_row++) {
      bkg_mean += img(i_row, i_col);
      count_b++;
    }
    for(R_len_t i_row = mat_r - margin; i_row < mat_r; i_row++) {
      bkg_mean += img(i_row, i_col);
      count_b++;
    }
  }
  bkg_mean /= count_b;
  
  // compute variance
  double bkg_var = 0.0;
  for(R_len_t i_col = extra; i_col < mat_c - extra; i_col++) {
    for(R_len_t i_row = 0; i_row < margin; i_row++) bkg_var += pow(img(i_row, i_col) - bkg_mean, 2);
    for(R_len_t i_row = mat_r - margin; i_row < mat_r; i_row++) bkg_var += pow(img(i_row, i_col) - bkg_mean, 2);
  }
  return NumericVector::create(_["BG_MEAN"] = bkg_mean, _["BG_STD"] = sqrt(bkg_var/(count_b - is_cif))); // n-1 for background
}

#endif
