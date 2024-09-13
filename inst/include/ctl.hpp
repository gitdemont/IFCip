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

#ifndef IFCIP_CTL_HPP
#define IFCIP_CTL_HPP

#include <Rcpp.h>
#include "padding.hpp"
using namespace Rcpp;

static int ifcip_ctl_dx [8]={ 1, 1, 0,-1,-1,-1, 0, 1};
static int ifcip_ctl_dy [8]={ 0, 1, 1, 1, 0,-1,-1,-1};
static int ifcip_ctl_bk [8]={ 7, 7, 1, 1, 3, 3, 5, 5};

//' @title Contour Tracing Connected Component Labeling
//' @name cpp_ctl
//' @description
//' This function is designed to identify connected component.
//' @param mat a LogicalMatrix, containing mask.
//' @param global whether to compute the perimeter globally or to evaluate the perimeter of each non 8-connected objects. Default is false.
//' When true pixels of overlapping extra borders of objects are counted only once.
//' @details adaptation of 'A linear-time component-labeling algorithm using contour tracing technique' from F. Chang, CJ. Chen and CJ Lu.
//' Computer Vision and Image Understanding Volume 93, Issue 2, February 2004, Pages 206-220.\doi{10.1016/j.cviu.2003.09.002}
//' @return a list whose members are:\cr
//' -matrix: an IntegerMatrix with connected component labels.\cr
//' -contours: an IntegerMatrix of identified contours, whose columns are x, y, label, direction and type.\cr
//' -nb_lab: the total number of components identified.
//' -perimeter: the number of pixels outside contours.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_ctl(const Rcpp::LogicalMatrix mat,
                   const bool global = false) {
  R_len_t mat_r = mat.nrow(), mat_c = mat.ncol();
  R_len_t i_cont = 0, i_per = -1, new_l = 2, max_count = ((mat_r >> 1) + (mat_r % 2) + 1) * ((mat_c >> 1) + (mat_c % 2) + 1) * 10;
  
  // create output matrix with extra cols / rows
  Rcpp::IntegerMatrix out(mat_r + 2, mat_c + 2);
  
  // create output contours vector
  // we consider we max size is when we have only isolated points
  // meaning that each point is separated by 1 background px
  // when recorded we store (5 + 5) values for each foreground px
  // so maximal final length should be
  // Rcpp::IntegerVector contours(((mat_r >> 1) + (mat_r % 2) + 1) * ((mat_c >> 1) + (mat_c % 2) + 1) * 10);
  Rcpp::IntegerVector contours(max_count);
  
  // create output perimeter
  Rcpp::IntegerVector perimeter;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t i_out = (i_col + 1) * (mat_r + 2) + 1;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      out[i_out++] = mat(i_row, i_col);
    }
  }
  
  for(R_len_t i_col = 1; i_col <= mat_c; i_col++) {
    for(R_len_t i_row = 1; i_row <= mat_r; i_row++) {
      if(out(i_row, i_col) >= 1) {
        if((out(i_row, i_col) == 1) && (out(i_row - 1, i_col) <= 0)) { // unprocess external contour
          // Rcpp::checkUserInterrupt();
          // Initialize external contour recording
          i_per++; 
          perimeter.push_back(0);
          
          // register current position + set 1st visiting position
          R_len_t len = i_cont, nbr_x = 0, nbr_y = 0;
          R_len_t ori_r = i_row, ori_c = i_col;
          unsigned short pos = 7, visit, j;
          
          // for debugging
          // Rcout << "ext ori, i_col: " << i_col + 1 << ", i_row: " << i_row + 1 << std::endl;
          
          // draw external contour
          while(out(i_row, i_col) > 0) {
            // Rcpp::checkUserInterrupt();
            for(j = 0; j < 8; j++) {
              // Rcpp::checkUserInterrupt();
              visit = (j + pos) % 8;
              nbr_x = i_col + ifcip_ctl_dx[visit];
              nbr_y = i_row + ifcip_ctl_dy[visit];
              
              // for debugging
              // Rcout << "go: " << visit << ", i_col: " << nbr_x << ", i_row: " << nbr_y << ", val: "<< out(nbr_y, nbr_x) << std::endl;
              
              if(out(nbr_y, nbr_x) > 0) {
                pos = ifcip_ctl_bk[visit];
                out(i_row, i_col) = new_l;
                break;
              } else {
                if(global) {
                  if(out(nbr_y, nbr_x) == 0) perimeter[i_per]++; 
                } else {
                  perimeter[i_per]++; 
                }
                out(nbr_y, nbr_x) = -1;
              }
            }
            // for debugging
            // Rcout << visit << "," << pos << "," << j << std::endl;
            if(j == 8) { 
              out(i_row, i_col) = new_l;
              if(max_count - i_cont < 5) { // safer
                Rcpp::stop("hpp_ctl: buffer overrun");
              }
              contours[i_cont++] = i_col;
              contours[i_cont++] = i_row;
              contours[i_cont++] = new_l - 1;
              contours[i_cont++] = 8;
              contours[i_cont++] = 1;
              break; // isolated point
            }
            if(//(contours[len + 0]== ori_c) && (contours[len + 1]== ori_r) &&
               (contours[len + 5] == nbr_x) && (contours[len + 6] == nbr_y) &&
               (i_row == ori_r) && (i_col == ori_c)) break; // escape contour
            
            // record current contour point coordinates + its label
            if(max_count - i_cont < 5) {
              Rcpp::stop("hpp_ctl: buffer overrun");
            }
            contours[i_cont++] = i_col;
            contours[i_cont++] = i_row;
            contours[i_cont++] = new_l - 1;
            contours[i_cont++] = visit;
            contours[i_cont++] = 1;
            // go to next neighbour point
            i_col = nbr_x;
            i_row = nbr_y;
          }
          new_l++;
        }
        else if(out(i_row + 1, i_col) == 0) { // unprocessed internal contour
          // Rcpp::checkUserInterrupt();
          // Initialize internal contour recording
          
          // register current position + set 1st visiting position
          R_len_t len = i_cont, nbr_x = 0, nbr_y = 0;
          R_len_t ori_r = i_row, ori_c = i_col;
          unsigned short pos = 3, visit, j;
          
          // for debugging
          // Rcout << "int ori, i_col: " << i_col + 1 << ", i_row: " << i_row + 1 << std::endl;
          
          // use label from the preceding point N in the scan line
          if(out(i_row,i_col) == 1) out(i_row,i_col) = out(i_row,i_col - 1);
          R_len_t exg_l = out(i_row,i_col);
          
          // draw internal contour
          while(out(i_row, i_col) > 0) {
            // Rcpp::checkUserInterrupt();
            for(j = 0; j < 8; j++) {
              // Rcpp::checkUserInterrupt();
              visit = (j + pos) % 8;
              nbr_x = i_col + ifcip_ctl_dx[visit];
              nbr_y = i_row + ifcip_ctl_dy[visit];
              
              // for debugging
              // Rcout << "go: " << visit << ", i_col: " << nbr_x << ", i_row: " << nbr_y << ", val: "<< out(nbr_y, nbr_x) << std::endl;
              
              if(out(nbr_y, nbr_x) > 0) {
                pos = ifcip_ctl_bk[visit];
                out(i_row, i_col) = exg_l;
                break;
              } else {
                out(nbr_y, nbr_x) = -1;
              }
            }
            if(//(contours[len + 0]== ori_c) && (contours[len + 1]== ori_r) && 
               (contours[len + 5] == nbr_x) && (contours[len + 6] == nbr_y) &&
               (i_row == ori_r) && (i_col == ori_c)) break; // escape contour
            
            // record current contour point coordinates + its label
            if(max_count - i_cont < 5) { // safer
              Rcpp::stop("hpp_ctl: buffer overrun");
            }
            contours[i_cont++] = i_col;
            contours[i_cont++] = i_row;
            contours[i_cont++] = exg_l - 1;
            contours[i_cont++] = visit;
            contours[i_cont++] = 2;
            // go to next neighbour point
            i_col = nbr_x;
            i_row = nbr_y;
          }
        }
        else if(out(i_row, i_col) == 1) out(i_row, i_col) = out(i_row, i_col - 1); // use label from the preceding point N in the scan line
      }
    }
  }
  // clip label to [0, label - 1]
  for(R_len_t i_col = 1; i_col <= mat_c; i_col++) {
    for(R_len_t i_row = 1; i_row <= mat_r; i_row++) {
      out(i_row, i_col) = (out(i_row, i_col) > 0) ? (out(i_row, i_col) - 1) : 0;
    }
  }
  
  // return contours to user-friendly output
  Rcpp::IntegerMatrix foo(i_cont / 5, 5);
  i_cont = 0;
  for(R_len_t i_row = 0; i_row < foo.nrow(); i_row++) {
    for(R_len_t i_col = 0; i_col < 5; i_col++) {
      foo(i_row, i_col) = contours[i_cont++];
    }
  }
  // Rcpp::IntegerVector dim = Rcpp::IntegerVector::create(mat_r, mat_c);
  colnames(foo) = Rcpp::CharacterVector::create("x","y","label","direction","type");
  Rcpp::List ret = Rcpp::List::create(_["matrix"] = out(Rcpp::Range(1, mat_r), Rcpp::Range(1, mat_c)),
                                      _["dim"] = Rcpp::IntegerVector::create(_["nrow"] = mat_r, _["ncol"] = mat_c), 
                                        _["contours"] = foo,
                                        _["nb_lab"] = new_l - 2,
                                        _["perimeter"] = perimeter);
  ret.attr("class") = "IFCip_ctl";
  return ret;
}

//' @title Contours Dilatation
//' @name cpp_dilate_ctl
//' @description
//' This function applies contours dilatation.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time dilate should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_dilate_ctl(const List ctl,
                                   const Rcpp::NumericMatrix kernel,
                                   const uint8_t iter = 0) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;
  Rcpp::NumericMatrix out = hpp_padding(as<Rcpp::NumericMatrix>(ctl["matrix"]), pad_r + kc, pad_c + kr, 5, 0.0);
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    Rcpp::NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
        double K = foo(i_row, i_col);
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col < i_col + pad_c_1; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row < i_row + pad_r_1; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) == 0) {
                out(f_row, f_col) = K;
              }
            }
          }
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

//' @title Contours Erosion
//' @name cpp_erode_ctl
//' @description
//' This function applies contours erosion.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param kernel, a NumericMatrix.
//' @param iter, an uint8_t, number of time erode should be iterated. Default is 0.
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_erode_ctl(const List ctl,
                                  const Rcpp::NumericMatrix kernel,
                                  const uint8_t iter = 0) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  R_len_t pad_c = kernel.ncol() >> 1;
  R_len_t pad_r = kernel.nrow() >> 1;
  R_len_t kc = kernel.ncol() % 2;
  R_len_t kr = kernel.nrow() % 2;
  R_len_t pad_c_1 = pad_c + kc;
  R_len_t pad_r_1 = pad_r + kr;
  if(kr == kc) kr = kc = !kc;
  Rcpp::NumericMatrix out = hpp_padding(as<Rcpp::NumericMatrix>(ctl["matrix"]), pad_r + kc, pad_c + kr, 5, 0.0);
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    Rcpp::NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < out.ncol() - pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < out.nrow() - pad_r; i_row++) {
        double K = foo(i_row, i_col);
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col < i_col + pad_c_1; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row < i_row + pad_r_1; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col) != 0) {
                out(f_row, f_col) = K;
              }
            }
          }
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r + kc * 2, out.nrow() - 1 - pad_r), Rcpp::Range(pad_c + kr * 2, out.ncol() - 1 - pad_c));
}

#endif
