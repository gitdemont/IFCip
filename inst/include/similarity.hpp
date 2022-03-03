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
//' It is a measure of the degree to which two images are linearly correlated within a masked region.\cr
//' See "Quantitative measurement of nuclear translocation events using similarity analysis of multispectral cellular images obtained in flow"
//' by T.C. George et al. Journal of Immunological Methods Volume 311, Issues 1–2, 20 April 2006, Pages 117-129 \doi{doi.org/10.1016/j.jim.2006.01.018}
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

//' @title Images Bright Detail Similarity Measurement
//' @name cpp_bright_similarity
//' @description
//' This function is designed to score similarity between two bright detail images.
//' @param img1 a NumericMatrix, containing bright detail image values.
//' @param img2 a NumericMatrix, containing bright detail image values.
//' @param msk a LogicalMatrix, containing mask.
//' @details the bright detail similarity is the non-mean normalized version of the log transformed Pearson's Correlation Coefficient.
//' It is designed to compare the small bright image detail of two images within a masked region.\cr
//' See "Quantitative analysis of protein co-localization on B cells opsonized with rituximab and complement using the ImageStream multispectral imaging flow cytometer"
//' by P.V. Beum et al. Journal of Immunological Methods Volume 317, Issues 1–2, 20 December 2006, Pages 90-99 \doi{doi.org/10.1016/j.jim.2006.09.012}
//' @return a double, the similarity.
//' @keywords internal
////' @export
// [[Rcpp::export]]
double hpp_bright_similarity(const Rcpp::NumericMatrix img1,
                             const Rcpp::NumericMatrix img2,
                             const Rcpp::LogicalMatrix msk) {
  R_len_t mat_r = img1.nrow();
  R_len_t mat_c = img1.ncol();
  if(mat_r != img2.nrow() || mat_c != img2.ncol() || mat_r != msk.nrow() || mat_c != msk.ncol()) {
    Rcpp::stop("hpp_bright_similarity: 'img1', 'img2' and 'msk' should have same dimensions");
  }
  
  // compute score
  double dif = 0.0, sq_dif1 = 0.0, sq_dif2 = 0.0;
  for(R_len_t i = 0; i < mat_r * mat_c; i++) {
    if(msk[i]) {
      dif += img1[i] * img2[i];
      sq_dif1 += img1[i] * img1[i];
      sq_dif2 += img2[i] * img2[i];
    }
  }
  
  // log transform
  double ret = dif / std::sqrt(sq_dif1 * sq_dif2);
  return std::log((1+ret)/(1-ret));
}

#endif
