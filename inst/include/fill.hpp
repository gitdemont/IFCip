/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2023 Yohann Demont                                              
                                                                                
  It is part of IFCip package, please cite:                                     
  -IFCip: An R Package for Imaging Flow Cytometry Image Processing              
  -YEAR: 2023                                                                   
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

#ifndef IFCIP_FILL_HPP
#define IFCIP_FILL_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Polygon Closing
//' @name close_polygon
//' @description
//' Pushes 1st row at the end of a matrix.
//' @param poly a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix close_polygon (const Rcpp::IntegerMatrix poly) {
  R_len_t i, n = poly.nrow();
  if(n > 0) {
    Rcpp::IntegerMatrix MM = Rcpp::no_init_matrix(n + 1, poly.ncol());
    for(i = 0; i < n; i++) MM(i, Rcpp::_) = poly(i, Rcpp::_);
    MM(i, Rcpp::_) = poly(0, Rcpp::_);
    return MM;
  }
  return poly;
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
//' @param x,y x and y coordinates of point to test.
//' @param poly a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
bool ray_pnt_in_poly(const R_len_t x,
                     const R_len_t y,
                     const Rcpp::IntegerMatrix poly) { 
  bool inside = false;
  for(R_len_t i = 0, j = poly.nrow() - 1; i < poly.nrow(); j = i++) {
    if(((poly(i,1) > y) != (poly(j,1) > y)) &&
       (x < (poly(j,0) - poly(i,0)) * (y - poly(i,1)) / (poly(j,1) - poly(i,1)) + poly(i,0))) 
      inside = !inside;
  }
  return inside;
}

//' @title Distance to Line Computation
//' @name comp_dist
//' @description
//' Compute distance from point with coords x, y to line passing through points x1, y1 and x2, y2.
//' @param x,x1,x2 x-coordinates.
//' @param y,y1,y2 y-coordinates.
//' @return a double.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double comp_dist(const R_len_t x1, const R_len_t y1,
                 const R_len_t x2, const R_len_t y2,
                 const double x, const double y) {
  R_len_t dx = x2-x1;
  R_len_t dy = y2-y1;
  return std::abs(dx*(y1-y)-dy*(x1-x)) / std::sqrt(dx*dx+dy*dy);
}

//' @title Trace Line
//' @name mk_line
//' @description
//' Make line from point x1, y1 to point x2, y2 in inret.
//' @param x1,x2 x-coordinates.
//' @param y1,y2 y-coordinates.
//' @param inret a LogicalMatrix where line will be drawn.
//' @return nothing /!\ inret is modified in place.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
void mk_line(const R_len_t x1, const R_len_t y1,
             const R_len_t x2, const R_len_t y2,
             Rcpp::LogicalMatrix inret) {
  R_len_t mat_r = inret.nrow() - 1, mat_c = inret.ncol() - 1;
  if((mat_r >= 0) && (mat_c >= 0)) {
    R_len_t y = y1, x = x1;
    double offsetx = 0.5, offsety = 0.5;
    if(x2 < x) offsetx = -offsetx;
    if(y2 < y) offsety = -offsety;
    inret(std::min(std::max(y, 0),mat_r), std::min(std::max(x, 0),mat_c)) = true;
    unsigned short count = 1;
    while(!((std::abs(y - y2) < 0.5) && (std::abs(x - x2) < 0.5))) {
      if((count++ % 10000) == 0) {
        count = 1;
        Rcpp::checkUserInterrupt();
      }
      Rcpp::NumericVector d = Rcpp::NumericVector::create(comp_dist(x1, y1, x2, y2, x + offsetx, y + offsety),
                                                          comp_dist(x1, y1, x2, y2, x + offsetx, y),
                                                          comp_dist(x1, y1, x2, y2, x, y + offsety));
      switch(Rcpp::which_min(d))
      {
      case 0 :
        x = x + 2 * offsetx;
        y = y + 2 * offsety;
        break;
      case 1 :
        x = x + 2 * offsetx;
        break;
      case 2 :
        y = y + 2 * offsety;
        break;
      }
      inret(std::min(std::max(y, 0),mat_r), std::min(std::max(x, 0),mat_c)) = true;
    }
  }
}

//' @title Polygon Mask
//' @name polymask
//' @description
//' Create a logical matrix of polygon contour.
//' @param poly a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param nrow, R_len_t of desired returned matrix rows.
//' @param ncol, R_len_t of desired returned matrix columns.
//' @return a LogicalMatrix
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix polymask(const Rcpp::IntegerMatrix poly,
                             const R_len_t nrow = 0,
                             const R_len_t ncol = 0) {
  Rcpp::LogicalMatrix out(std::max(nrow, 0), std::max(ncol, 0));
  if((nrow <= 0) || (ncol <= 0)) return out;
  if(poly.ncol() < 2) Rcpp::stop("'poly' should have at least 2 columns");
  R_len_t pts_row = poly.nrow();
  for(R_len_t i = 0; i < pts_row - 1; i++) {
    mk_line(poly[i] - 1, poly[i + pts_row] - 1, poly[i + 1] - 1, poly[i + 1 + pts_row] - 1, out);
  }
  if(pts_row > 0) mk_line(poly[pts_row - 1] - 1, poly[pts_row + pts_row - 1] - 1, poly[0] - 1, poly[pts_row] - 1, out); 
  return out;
}

//' @title Polygon Drawing
//' @name polydraw
//' @description
//' This function is designed to trace and fill polygon inside a matrix.
//' @param poly, a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param inret, an IntegerMatrix to be filled.
//' @param border, an int used to trace polygon border.
//' @param fill, an int used to fill polygon.
//' @return nothing /!\ inret is modified in place.
//' @keywords internal
void polydraw(const Rcpp::IntegerMatrix poly,
              Rcpp::IntegerMatrix inret,
              const int border = 1,
              const int fill = 1) {
  Rcpp::LogicalMatrix mask = polymask(poly, inret.nrow(), inret.ncol());// make a polygon mask
  Rcpp::IntegerMatrix pts = close_polygon(poly);
  R_len_t mat_r = inret.nrow(), mat_c = inret.ncol();
  
  // bound to polygon
  R_len_t start_r = mat_r + 1, stop_c = 0, start_c = mat_c + 1, stop_r = 0;
  for(R_len_t i = 0; i < poly.nrow(); i++) {
    if((poly(i, 1) - 2) < start_r) start_r = poly(i, 1) - 2;
    if(poly(i, 1) >  stop_r) stop_r  = poly(i, 1);
    if((poly(i, 0) - 2) < start_c) start_c = poly(i, 0) - 2;
    if(poly(i, 0) >  stop_c) stop_c  = poly(i, 0);
  }
  
  // bound at most to inret limits
  start_c = std::max(0, start_c);
  start_r = std::max(0, start_r);
  stop_c = std::min(mat_c, stop_c);
  stop_r = std::min(mat_r, stop_r);
  
  // trace and fill polygon
  for(R_len_t i_col = start_c; i_col < stop_c; i_col++) {
    bool edge = false;
    R_len_t beg = stop_r;
    for(R_len_t i_row = start_r; i_row < stop_r; i_row++) {
      if(mask(i_row, i_col)) {                                          // an edge is found
        for(R_len_t i = beg; i < i_row; i++) {                          // fill line from last registered starting position to current i_row
          if(!mask(i, i_col)) inret(i, i_col) = fill;                   // avoid coloring polygon contour with fill color
        }
        inret(i_row, i_col) = border;                                   // trace polygon border (i.e. current point)
        edge = true;
      } else {
        if(edge) {                                                      // an edge has been encountered
          edge = false;                                                 // reset edge flag
          beg = stop_r;                                                 // reset edge position
          if(ray_pnt_in_poly(i_col + 1, i_row + 1, pts)) beg = i_row;   // check if point right after border (in the scan line) is in the polygon, if so register current position as starting point
        }
      }
    }
  }
}

//' @title Contours Filling
//' @name cpp_fill
//' @description
//' This function is designed to fill contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param label an int corresponding to the label of desired set of contour to be filled.
//' Default is 0 to fill all set of contours found.
//' @param inner a bool, to whether or not fill hole(s) inside contours if some where identified.
//' @param outer a bool, to whether or not fill contours outside hole(s) if some where identified.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
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
  
  // retrieve contours coordinates
  Rcpp::IntegerMatrix contours = as<Rcpp::IntegerMatrix>(ctl["contours"]);
  
  // retrieve dimension of original image used to determine contours
  Rcpp::IntegerVector dim = as<Rcpp::IntegerVector>(ctl["dim"]);
  
  // create out 
  Rcpp::IntegerMatrix out(dim[0], dim[1]);
  
  if(label_cur <= label_max) {
    // count points in inner/outer contours
    Rcpp::IntegerVector inn_s = Rcpp::IntegerVector(label_max + 1);
    Rcpp::IntegerVector out_s = Rcpp::IntegerVector(label_max + 1);
    for(R_len_t k = 0; k < contours.nrow(); k++) {
      if(contours(k, 2) > 0) {
        if(contours(k, 4) == 1) out_s[contours(k, 2)]++;
        if(contours(k, 4) == 2) inn_s[contours(k, 2)]++;
      }
    }
    
    // fill contours with desired label(s)
    label_max = label_cur == label ? label : 1;
    while(label_cur >= label_max) {
      if(outer) { // extract and fill external contour
        Rcpp::IntegerMatrix poly_out = Rcpp::no_init_matrix(out_s[label_cur], 2);
        for(R_len_t i = 0, k = 0; k < contours.nrow(); k++) {
          if((contours(k, 2) == label_cur) && (contours(k, 4) == 1)) {
            poly_out[i] = contours(k, 0);
            poly_out[poly_out.nrow() + i++] = contours(k, 1);
          }
        }
        polydraw(poly_out, out, label_cur, label_cur);
      }
      if(!inner) {// extract and fill internal contour
        Rcpp::IntegerMatrix poly_in = Rcpp::no_init_matrix(inn_s[label_cur], 2);
        for(R_len_t i = 0, k = 0; k < contours.nrow(); k++) {
          if((contours(k, 2) == label_cur) && (contours(k, 4) == 2)) {
            poly_in[i] = contours(k, 0);
            poly_in[poly_in.nrow() + i++] = contours(k, 1);
          }
        }
        polydraw(poly_in, out, label_cur, 0);
      }
      label_cur--;
    }
  }
  return out;
}

//' @title Contours Filling Outer Only
//' @name cpp_fill_out
//' @description
//' This function is designed to fill the most external contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_fill_out(const List ctl) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) {
    Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  }
  // retrieve max number of sets of contours identified by ctl
  int label_max = as<uint32_t>(ctl["nb_lab"]);
  
  // retrieve contours coordinates
  Rcpp::IntegerMatrix contours = as<Rcpp::IntegerMatrix>(ctl["contours"]);
  
  // retrieve dimension of original image used to determine contours
  Rcpp::IntegerVector dim = as<Rcpp::IntegerVector>(ctl["dim"]);
  
  // create out 
  Rcpp::IntegerMatrix out(dim[0], dim[1]);
  
  // count points in inner/outer contours
  Rcpp::IntegerVector out_s = Rcpp::IntegerVector(label_max + 1);
  for(R_len_t k = 0; k < contours.nrow(); k++) {
    if((contours(k, 2) > 0) && (contours(k, 4) == 1)) out_s[contours(k, 2)]++;
  }
  
  // fill all contours
  for(R_len_t label_cur = 1; label_cur < out_s.size(); label_cur++) {
    // extract external contours points coordinates
    Rcpp::IntegerMatrix poly = Rcpp::no_init_matrix(out_s[label_cur], 2);
    for(R_len_t i = 0, k = 0; k < contours.nrow(); k++) {
      if((contours(k, 2) == label_cur) && (contours(k, 4) == 1)) {
        poly[i] = contours(k, 0);
        poly[poly.nrow() + i++] = contours(k, 1);
      }
    }
    // fill polygon of points coordinates
    polydraw(poly, out, label_cur, label_cur);
  }
  return out;
}

#endif
