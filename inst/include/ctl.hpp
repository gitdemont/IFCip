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
// [[Rcpp::depends(IFC)]]
#include <gate.hpp>
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
//' Computer Vision and Image Understanding Volume 93, Issue 2, February 2004, Pages 206-220.\url{https://doi.org/10.1016/j.cviu.2003.09.002}
//' @return a list whose members are:\cr
//' -matrix: an IntegerMatrix with connected component labels.\cr
//' -contours: an IntegerMatrix of identified contours, whose columns are x, y, label, direction and type.\cr
//' -nb_lab: the total number of components identified.
//' -perimeter: the number of pixels outside contours.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_ctl(const Rcpp::LogicalMatrix mat,
                   const bool global = false) {
  R_len_t mat_r = mat.nrow(), mat_c = mat.ncol();
  R_len_t i_cont = 0, i_per = -1, new_l = 2, max_count = (std::ceil(mat_r / 2) + 1) * (std::ceil(mat_c / 2) + 1) * 10;
  
  // create output matrix with extra cols / rows
  Rcpp::IntegerMatrix out(mat_r + 2, mat_c + 2);
  
  // create output contours vector
  // we consider we max size is when we have only isolated points
  // meaning that each point is separated by 1 background px
  // when recorded we store (5 + 5) values for each foreground px
  // so maximal final length should be
  // Rcpp::IntegerVector contours((ceil(mat_r / 2) + 1) * (ceil(mat_c / 2) + 1) * 10);
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
          if(max_count - i_cont < 5) { // safer
            Rcpp::stop("hpp_ctl: buffer overrun");
          }
          contours[i_cont++] = -2;
          contours[i_cont++] = -2;
          contours[i_cont++] = -2;
          contours[i_cont++] =  8;
          contours[i_cont++] =  1;
          
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
          if(max_count - i_cont < 5) { // safer
            Rcpp::stop("hpp_ctl: buffer overrun");
          }
          contours[i_cont++] = -3;
          contours[i_cont++] = -3;
          contours[i_cont++] = -3;
          contours[i_cont++] =  8;
          contours[i_cont++] =  2;
          
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
  Rcpp::NumericMatrix foo(i_cont / 5, 5);
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

//' @title ray_pnt_in_poly
//' @description
//' This function checks if points lie within a polygonusing an adaptation of the Ray Casting algorithm.
//' Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: 
//' 1.Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers. 
//' 2.Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution. 
//' 3.The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 
//' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
//' @source adaptation from W. Randolph Franklin code \url{https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html}
//' @param pnt a NumericVector with x and y coordinates.
//' @param poly a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param include whether to include vertices of the polygon or not
//' @keywords internal
bool ray_pnt_in_poly(const Rcpp::NumericVector pnt,
                     const Rcpp::NumericMatrix poly,
                     const bool include = true) { 
  bool inside = false;
  for(R_len_t i = 0, j = poly.nrow() - 1; i < poly.nrow(); j = i++) {
    // start
    // this part is the main adaptation. to include the vertice or not
    // thanks to ctl every adjacent pixel of the contours are registered
    // so they are contiguous and we can easily include them or not
    if(include && (pnt[0] == poly(i, 0)) && (pnt[1] == poly(i, 1))) return true;
    // end
    if(((poly(i,1) > pnt[1]) != (poly(j,1) > pnt[1])) &&
       (pnt[0] < (poly(j,0) - poly(i,0)) * (pnt[1] - poly(i,1)) / (poly(j,1) - poly(i,1)) + poly(i,0))) 
      inside = !inside;
  }
  return inside;
}

//' @title Contours Filling
//' @name cpp_fill
//' @description
//' This function is designed to fill contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param label an uint32_t corresponding to the label of desired set of contour to be filled.
//' Default is 0 to fill all set of contours found.
//' @param inner a bool, to whether or not fill hole(s) inside contours if some where identified
//' @param outer a bool, to whether or not fill contours outside hole(s) if some where identified
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_fill(const List ctl,
                             const int label = 0,
                             const bool inner = true,
                             const bool outer = true) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  // retrieve max number of sets of contours identified by ctl
  int label_max = as<uint32_t>(ctl["nb_lab"]);
  // check whether to fill only one set or every sets of contours
  int label_cur = label <= 0 ? label_max : label;
  label_max = label_cur == label ? label : 1;
  
  // set up the type of object to fill (i.e) inside contours, w or w/o hole(s)
  Rcpp::LogicalVector type = Rcpp::LogicalVector::create(false, outer, inner); 
  // if inside contours + hole(s) no need to compute for hole since they are 
  // also inside contours by definition
  if(inner && outer) type[2] = false;
  // retrieve contours coordinates
  Rcpp::IntegerMatrix contours = Rcpp::clone(as<Rcpp::IntegerMatrix>(ctl["contours"]));
  // retrieve label and type from contours
  Rcpp::IntegerVector V1 = Rcpp::no_init_vector(contours.nrow());
  Rcpp::IntegerVector V2 = Rcpp::clone(V1);
  for(R_len_t i = 0; i < V1.size(); i++) {
    V1[i] = contours(i, 2); // 2 is label_cur calling contours(i, "label") is not working
    V2[i] = contours(i, 4); // 4 is type
  }
  // retrieve dimension of original image used to determine contours
  Rcpp::IntegerVector dim = Rcpp::clone(as<Rcpp::IntegerVector>(ctl["dim"]));
  // create all possible combinations of x, y coordinates within original image
  // we will check if theses pts coordinates are within contours or not
  Rcpp::NumericMatrix pnts(dim[0] * dim[1], 2);
  for(R_len_t i_col = 1, i = 0; i_col <= dim[1]; i_col++) {
    for(R_len_t i_row = 1; i_row <= dim[0]; i_row++) {
      pnts(i, 0) = i_col;
      pnts(i++, 1) = i_row;
    }
  }
  // create out 
  Rcpp::IntegerMatrix out(dim[0] + 1, dim[1] + 1);
  
  while(label_cur >= label_max) {
    for(short k = 2; k >= 1; k--) { // try internal contour 1st (k=2), then external
      if(type[k]) {
        // retrieve the contours for label and type (i.e. internal or external contour)
        NumericVector V4;
        double xmin = R_PosInf, xmax = R_NegInf, ymin = R_PosInf, ymax = R_NegInf;
        for(R_len_t i = 0; i < V1.size(); i++) {
          if((V1[i] == label_cur) && (V2[i] == k) && (contours(i, 0)) >= 0 && (contours(i, 1)) >= 0) {
            V4.push_back(contours(i, 0) ); 
            V4.push_back(contours(i, 1) );
            if(contours(i, 0) < xmin) xmin = contours(i, 0);
            if(contours(i, 1) < ymin) ymin = contours(i, 1);
            if(contours(i, 0) > xmax) xmax = contours(i, 0);
            if(contours(i, 1) > ymax) ymax = contours(i, 1);
          }
        }
        
        if(V4.size() > 0) { // check if a contour of current label and type was found
          NumericMatrix gate = Rcpp::no_init_matrix(2, V4.size()/2);
          std::copy(V4.begin(), V4.end(), gate.begin());
          NumericMatrix poly = close_polygon(Rcpp::transpose(gate));
          Rcpp::LogicalVector is_in(pnts.nrow());
          if(k == 1) {
            for(R_len_t i = 0; i < pnts.nrow(); i++) {
              is_in[i] = (out(pnts(i, 1), pnts(i, 0)) == 0) &&
                (pnts(i, 1) >= ymin) && (pnts(i, 1) <= ymax) &&
                (pnts(i, 0) >= xmin) && (pnts(i, 0) <= xmax) &&
                ray_pnt_in_poly(pnts(i, Rcpp::_), poly, true);
            } 
          } else {
            for(R_len_t i = 0; i < pnts.nrow(); i++) {
              is_in[i] = (out(pnts(i, 1), pnts(i, 0)) == 0) &&
                (pnts(i, 1) > ymin) && (pnts(i, 1) < ymax) &&
                (pnts(i, 0) > xmin) && (pnts(i, 0) < xmax) &&
                ray_pnt_in_poly(pnts(i, Rcpp::_), poly, false);
            }
          }
          
          if(outer & !inner) {
            NumericVector V5;
            double in_xmin = R_PosInf, in_xmax = R_NegInf, in_ymin = R_PosInf, in_ymax = R_NegInf;
            for(R_len_t i = 0; i < V1.size(); i++) {
              if((V1[i] == label_cur) && (V2[i] == 2) && (contours(i, 0)) >= 0 && (contours(i, 1)) >= 0) {
                V5.push_back(contours(i, 0) ); 
                V5.push_back(contours(i, 1) ); 
                if(contours(i, 0) < in_xmin) in_xmin = contours(i, 0);
                if(contours(i, 1) < in_ymin) in_ymin = contours(i, 1);
                if(contours(i, 0) > in_xmax) in_xmax = contours(i, 0);
                if(contours(i, 1) > in_ymax) in_ymax = contours(i, 1);
              }
            }
            if(V5.size() > 0) {  // check if a contour of current label and type 2 was found
              NumericMatrix gate_out = Rcpp::no_init_matrix(2, V5.size()/2);
              std::copy(V5.begin(), V5.end(), gate_out.begin());
              NumericMatrix poly_out = close_polygon(Rcpp::transpose(gate_out));
              
              for(R_len_t i = 0; i < pnts.nrow(); i++) { // fill inside the external contour but not inside the internal hole
                if(is_in[i] &&
                   !((out(pnts(i, 1), pnts(i, 0)) < label_cur) &&
                   (pnts(i, 1) > in_ymin) && (pnts(i, 1) < in_ymax) &&
                   (pnts(i, 0) > in_xmin) && (pnts(i, 0) < in_xmax) &&
                   ray_pnt_in_poly(pnts(i, Rcpp::_), poly, false))) out(pnts(i, 1), pnts(i, 0)) = label_cur;
              }
            } else { 
              for(R_len_t i = 0; i < pnts.nrow(); i++) { // fill inside the contour 
                if(is_in[i]) out(pnts(i, 1), pnts(i, 0)) = label_cur;
              }
            }
          } else {
            for(R_len_t i = 0; i < pnts.nrow(); i++) { // fill inside the contour 
              if(is_in[i]) out(pnts(i, 1), pnts(i, 0)) = label_cur;
            }
          }
        }
      }
    }
    label_cur--;
  }
  return out(Rcpp::Range(1, dim[0]), Rcpp::Range(1, dim[1]));
}

//' @title Contours Filling Outer Only
//' @name cpp_fill_out
//' @description
//' This function is designed to fill the most external contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_fill_out(const List ctl) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  // retrieve max number of sets of contours identified by ctl
  int label_max = as<uint32_t>(ctl["nb_lab"]);
  
  // retrieve contours coordinates
  Rcpp::IntegerMatrix contours = Rcpp::clone(as<Rcpp::IntegerMatrix>(ctl["contours"]));
  
  // retrieve dimension of original image used to determine contours
  Rcpp::IntegerVector dim = Rcpp::clone(as<Rcpp::IntegerVector>(ctl["dim"]));
  
  // create VV vector of external contours points length
  Rcpp::IntegerVector VV(label_max);
  for(R_len_t label_cur = 1; label_cur <= label_max; label_cur++) {
    for(R_len_t k = Rcpp::sum(VV); k < contours.nrow(); k++) {
      if((contours(k, 2) == label_cur) && (contours(k, 4) == 1)) VV[label_cur - 1]++;
    }
  }
  
  // create out 
  Rcpp::IntegerMatrix out(dim[0] + 1, dim[1] + 1);
  
  // for each external contours
  for(R_len_t i_V = 0; i_V < VV.size(); i_V++) {
    // extract contours points coordinates
    Rcpp::IntegerVector x = no_init_vector(VV[i_V]);
    Rcpp::IntegerVector y = clone(x);
    for(R_len_t i = 0, k = 0; k < contours.nrow(); k++) {
      if((contours(k, 2) == i_V + 1) && (contours(k, 4) == 1)) {
        x[i] = contours(k, 0);
        y[i++] = contours(k, 1);
      }
    }
    // extract all unique x (columns)
    Rcpp::IntegerVector V1 = Rcpp::unique(x);
    for(R_len_t i_V1 = 0; i_V1 < V1.size(); i_V1++) {
      // for each unique x extract and sort all unique y (rows)
      Rcpp::IntegerVector V2;
      for(R_len_t i_x = 0; i_x < x.size(); i_x++) if(x[i_x] == V1[i_V1]) V2.push_back(y[i_x]);
      Rcpp::IntegerVector V3 = sort_unique(V2);
      // check if 1st contour point is background
      // if not it means it belongs, i.e. it is inside, of another contour
      // so we escape (do nothing)
      if(!out(V3[0], V1[i_V1])) {
        // fill contours
        int color = i_V + 1;
        out(V3[0], V1[i_V1]) = color;
        for(R_len_t i_V3 = 1; i_V3 < V3.size(); i_V3++) {
          for(R_len_t i_row = V3[i_V3 - 1]; i_row <= V3[i_V3]; i_row++) out(i_row, V1[i_V1]) = color;
          // check if current i_row and former i_row are separated by at least one px
          // in such case, we swap from foreground to background and reversely
          if(V3[i_V3] - V3[i_V3 - 1] > 1) color = color ? 0 : i_V + 1;
        }
      }
    }
  }
  return out(Rcpp::Range(1, dim[0]), Rcpp::Range(1, dim[1]));
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
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_dilate_ctl(const List ctl,
                                   const Rcpp::NumericMatrix kernel,
                                   const uint8_t iter = 0) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  Rcpp::IntegerVector dim = Rcpp::clone(as<Rcpp::IntegerVector>(ctl["dim"]));
  R_len_t mat_r = dim[0];
  R_len_t mat_c = dim[1];
  Rcpp::List pad = hpp_padding(as<Rcpp::NumericMatrix>(ctl["matrix"]), kernel, 5, 0);
  Rcpp::NumericMatrix out = pad["out"];
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    Rcpp::NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = foo(i_row, i_col);
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
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
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
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
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_erode_ctl(const List ctl,
                                  const Rcpp::NumericMatrix kernel,
                                  const uint8_t iter = 0) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  Rcpp::IntegerVector dim = Rcpp::clone(as<Rcpp::IntegerVector>(ctl["dim"]));
  R_len_t mat_r = dim[0];
  R_len_t mat_c = dim[1];
  Rcpp::List pad = hpp_padding(as<Rcpp::NumericMatrix>(ctl["matrix"]), kernel, 5, 0);
  Rcpp::NumericMatrix out = pad["out"];
  R_len_t pad_r = pad["ori_r"];
  R_len_t pad_c = pad["ori_c"];
  
  uint8_t i_iter = 0;
  while(i_iter <= iter) {
    i_iter++;
    Rcpp::NumericMatrix foo = clone(out);
    for(R_len_t i_col = pad_c; i_col < mat_c + pad_c; i_col++) {
      for(R_len_t i_row = pad_r; i_row < mat_r + pad_r; i_row++) {
        double K = foo(i_row, i_col);
        for(R_len_t f_col = i_col - pad_c, i_ker = 0; f_col <= i_col + pad_c; f_col++) {
          for(R_len_t f_row = i_row - pad_r; f_row <= i_row + pad_r; f_row++) {
            if(kernel[i_ker++]) {
              if(foo(f_row, f_col)) {
                out(f_row, f_col) = K;
              }
            }
          }
        }
      }
    }
  }
  return out(Rcpp::Range(pad_r , mat_r + pad_r - 1), Rcpp::Range(pad_c , mat_c + pad_c - 1));
}

#endif
