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

#ifndef IFCIP_SIMILARITY_HPP
#define IFCIP_SIMILARITY_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Images Similarity Measurement
//' @name cpp_similarity
//' @description
//' This function is designed to score similarity between two images.
//' @param img1 a NumericMatrix, containing image values.
//' @param img2 a NumericMatrix, containing image values.
//' @param msk a LogicalMatrix, containing mask.
//' @details the similarity is the log transformed Pearson's Correlation Coefficient.
//' It is a measure of the degree to which two images are linearly correlated within a masked region.
//' @return a double, the similarity.
//' @keywords internal
////' @export
// [[Rcpp::export]]
double hpp_similarity(const Rcpp::NumericMatrix img1,
                      const Rcpp::NumericMatrix img2,
                      const Rcpp::LogicalMatrix msk) {
  R_len_t mat_r = img1.nrow();
  R_len_t mat_c = img1.ncol();
  if(mat_r != img2.nrow() || mat_c != img2.ncol() || mat_r != msk.nrow() || mat_c != msk.ncol()) {
    Rcpp::stop("hpp_similarity: 'img1', 'img2' and 'msk' should have same dimensions");
  }
  // formula in 'Highlighting curcumin-induced crosstalk between autophagy and apoptosis: A biochemical approach coupling impedancemetry, imaging, and flow cytometry.'
  // from F. J. Sala de Oyanguren,  N. E. Rainey, A. Moustapha, A. Saric, F. Sureau, J-E O�Connor, P. X. Petit. https://doi.org/10.1101/827279 Fig.5.
  
  // compute means
  R_len_t i = 0, count = 0;
  double avg1 = 0.0, avg2 = 0.0, dif = 0.0, sq_dif1 = 0.0, sq_dif2 = 0.0, foo = 0.0, bar = 0.0;
  for(i = 0; i < mat_r * mat_c; i++) {
    if(msk[i]) {
      count++;
      avg1 += img1[i];
      avg2 += img2[i];
    }
  }
  avg1 /= count;
  avg2 /= count;
  
  // compute  diff
  for(i = 0; i < mat_r * mat_c; i++) {
    if(msk[i]) {
      foo = img1[i] - avg1;
      bar = img2[i] - avg2;
      dif += foo * bar;
      sq_dif1 += foo * foo;
      sq_dif2 += bar * bar;
    }
  }
  
  // log transform
  double ret = dif / std::sqrt(sq_dif1 * sq_dif2);
  return std::log((1+ret)/(1-ret));
}

#endif
