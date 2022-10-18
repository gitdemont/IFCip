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

#ifndef IFCIP_MASK_HPP
#define IFCIP_MASK_HPP

#include <Rcpp.h>
#include "matrix_logic.hpp"
using namespace Rcpp;

//' @title Matrix Threshold
//' @name cpp_threshold
//' @description
//' This function takes a image and checks for members superior or equal
//' to max(img) - k * diff(range(img)) / 100 within msk
//' @param img a NumericMatrix
//' @param msk a NumericMatrix
//' @param k constant to be checked. Default is 0.
//' @param removal uint8_t, object removal method. Default is 0 for no removal. Otherwise, if\cr
//' -1, for clipped removal, keep non clipped foreground.\cr
//' -2, height clipped removal.\cr
//' -3, width clipped removal.\cr
//' -4, only keep background:.\cr
//' -5, only keep foreground.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix hpp_threshold(const Rcpp::NumericMatrix img,
                                  const Rcpp::NumericMatrix msk,
                                  const double k = 0.0,
                                  uint8_t removal = 0) {
  if(!(msk.ncol() == img.ncol()) && (msk.nrow() == img.nrow())) Rcpp::stop("hpp_threshold: 'img' and 'msk' should have same dimensions");
  Rcpp::LogicalMatrix OUT_M(msk.nrow(), msk.ncol());
  OUT_M.fill(false);
  Rcpp::NumericVector ran = NumericVector::create(R_NegInf, R_PosInf);
  // determines range within msk
  switch(removal) {
  case 1: { // only keep non clipped foreground
    for(R_len_t i = 0; i < msk.size() ; i++) {
    if(msk[i] == 1) {
      if(img[i] > ran[0]) ran[0] = img[i];
      if(img[i] < ran[1]) ran[1] = img[i];
      OUT_M[i] = true;
    }
  }
    break;
  }
  case 2: { // only keep non height clipped foreground
    for(R_len_t i = 0; i < msk.size() ; i++) {
    if(msk[i] && (msk[i] != 2)) {
      if(img[i] > ran[0]) ran[0] = img[i];
      if(img[i] < ran[1]) ran[1] = img[i];
      OUT_M[i] = true;
    }
  }
    break;
  }
  case 3: { // only keep non width clipped foreground
    for(R_len_t i = 0; i < msk.size() ; i++) {
    if(msk[i] && (msk[i] != 3)) {
      if(img[i] > ran[0]) ran[0] = img[i];
      if(img[i] < ran[1]) ran[1] = img[i];
      OUT_M[i] = true;
    }
  }
    break;
  }
  case 4: { // only keep background
    for(R_len_t i = 0; i < msk.size() ; i++) {
    if(!msk[i]) {
      if(img[i] > ran[0]) ran[0] = img[i];
      if(img[i] < ran[1]) ran[1] = img[i];
      OUT_M[i] = true;
    }
  }
    break;
  }
  case 5: { // foreground
    for(R_len_t i = 0; i < msk.size() ; i++) {
    if(msk[i]) {
      if(img[i] > ran[0]) ran[0] = img[i];
      if(img[i] < ran[1]) ran[1] = img[i];
      OUT_M[i] = true;
    }
  }
    break;
  }
  default: { // no removal
    for(R_len_t i = 0; i < msk.size() ; i++) {
    if(img[i] > ran[0]) ran[0] = img[i];
    if(img[i] < ran[1]) ran[1] = img[i];
    OUT_M[i] = true;
  }
    break;
  }
  }
  double kk = ran[1] - k * (ran[1] - ran[0]) / 100 ;
  
  // apply threshold
  for(R_len_t i = 0; i < msk.size() ; i++) {
    if(OUT_M[i]) OUT_M[i] = (img[i] >= kk);
  }
  return OUT_M;
}

#endif
