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
#include "fifo.hpp"
#include "geometry.hpp"
#include "ctl.hpp"
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
    bool same = true;
    for(R_len_t k = 0; k < poly.ncol(); k++) if(poly(0, k) != poly(n - 1, k)) { same = false; break; }
    if(same) return poly;
    Rcpp::IntegerMatrix MM = Rcpp::no_init_matrix(n + 1, poly.ncol());
    for(i = 0; i < n; i++) MM(i, Rcpp::_) = poly(i, Rcpp::_);
    MM(i, Rcpp::_) = poly(0, Rcpp::_);
    return MM;
  }
  return poly;
}

//' @title ray_pnt_in_poly
//' @description
//' This function checks if points lie within a polygon using an adaptation of the Ray Casting algorithm.
//' @source adaptation from W. Randolph Franklin code \url{https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html}\cr
//' /verb{
//' Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: 
//' 1.Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers. 
//' 2.Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution. 
//' 3.The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 
//' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
//' }
//' @param x,y x and y coordinates of point to test.
//' @param poly a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
bool ray_pnt_in_poly (const R_len_t x,
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

//' scanfill template
template <int RTYPE>
void scanfill_T (Rcpp::Matrix<RTYPE> m,
                 const R_len_t i,
                 const double v,
                 const double tol,
                 Rcpp::IntegerVector &Z,
                 const R_len_t buf = 256) {
  if(i < 0 || i >= m.size()) return;
  double o = m[i];
  bool is_na = traits::is_na<14>(tol);
  double n = is_na ? tol : v;
  if(o == n || (is_na && traits::is_na<RTYPE>(o))) return;
  unsigned short ptr = 1;
  R_len_t mat_r = m.nrow(), mat_c = m.ncol();
  Rcpp::IntegerVector Q = queue_create(buf);
  queue_push(Q, i);
  if(is_na) {
    while(Q[0] != 0) {
      if(((ptr)++ % 10000) == 0) {
        if(Q[0] >= 2147463648) Rcpp::stop("scanfill: endless propagation"); // 2^31 - 2 * 10000
        ptr = 1;
        Rcpp::checkUserInterrupt();
      }
      R_len_t p = queue_pop(Q), i_col = p / mat_r, i_row = p - i_col * mat_r, yy = i_row - 1;
      while((yy >= 0) && (!traits::is_na<RTYPE>(m[yy + mat_r * i_col]))) yy--;
      yy++;
      while((yy < mat_r) && (!traits::is_na<RTYPE>(m[yy + mat_r * i_col]))) {
        m[yy + mat_r * i_col] = n;
        queue_push(Z, yy + mat_r * i_col);
        if((i_col > 0) && (!traits::is_na<RTYPE>(m[yy + mat_r * (i_col - 1)]))) queue_push(Q, yy + mat_r * (i_col - 1));
        if((i_col < mat_c - 1) && (!traits::is_na<RTYPE>(m[yy + mat_r * (i_col + 1)]))) queue_push(Q, yy + mat_r * (i_col + 1));
        yy++;
      }
    }
  } else {
    while(Q[0] != 0) {
      if(((ptr)++ % 10000) == 0) {
        if(Q[0] >= 2147463648) Rcpp::stop("scanfill: endless propagation"); // 2^31 - 2 * 10000
        ptr = 1;
        Rcpp::checkUserInterrupt();
      }
      R_len_t p = queue_pop(Q), i_col = p / mat_r, i_row = p - i_col * mat_r, yy = i_row - 1;
      if(traits::is_na<RTYPE>(m[p]) || (m[p] == n)) continue;
      while((yy >= 0) && (std::abs(m[yy + mat_r * i_col] - o) <= tol)) yy--;
      yy++;
      bool goL = false, goR = false;
      while((yy < mat_r) && (std::abs(m[yy + mat_r * i_col] - o) <= tol)) {
        m[yy + mat_r * i_col] = n;
        // queue_push(Z, yy + mat_r * i_col); not needed except if we want to record order of scan
        if(!goL && (i_col > 0) && (std::abs(m[yy + mat_r * (i_col - 1)] - o) <= tol)) {
          queue_push(Q, yy + mat_r * (i_col - 1));
          goL = true;
        } else if(goL && (i_col > 0) && (std::abs(m[yy + mat_r * (i_col - 1)] - o) > tol)) goL = false;
        if(!goR && (i_col < mat_c - 1) && (std::abs(m[yy + mat_r * (i_col + 1)] - o) <= tol)) {
          queue_push(Q, yy + mat_r * (i_col + 1));
          goR = true;
        } else if(goR && (i_col < mat_c - 1) && (std::abs(m[yy + mat_r * (i_col + 1)] - o) > tol)) goR = false;
        yy++;
      }
    }
  }
}

//' @title Scanline Flood Filling
//' @name scanfill
//' @description
//' Flood fills image region.
//' @param m, a Matrix.
//' @param i, a R_len_t, index within 'm' to start filling.
//' @param v, a double, value to assign to connected pixels. Note that when 'tol' is \code{NA}, 'v' will be forced to \code{NA}.
//' @param tol, a double, tolerance between fill color and connected pixels. Use \code{NA}, for filling every pixel inside 'poly'.
//' Note that when 'tol' is \code{NA}, 'v' is forced to \code{NA}, so it is important that connected points are bounded by \code{NA},
//' otherwise the whole 'm' will be scanned.
//' @param Z, a queue, used to hold filled indices. It is mainly used to allow filling back with 'v' when 'tol' is \code{NA}.
//' @param buf, a R_len_t, initial size for the queue. Default is 256.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @return nothing, but m will be modified in-place.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
void scanfill (SEXP m,
               const R_len_t i,
               const double v,
               const double tol,
               Rcpp::IntegerVector &Z,
               const R_len_t buf = 256) {
  switch( TYPEOF(m) ) {
  case LGLSXP : return scanfill_T(as<Rcpp::LogicalMatrix>(m), i, v, tol, Z, buf);
  case INTSXP : return scanfill_T(as<Rcpp::IntegerMatrix>(m), i, v, tol, Z, buf);
  case REALSXP : return scanfill_T(as<Rcpp::NumericMatrix>(m), i, v, tol, Z, buf);
  default : Rcpp::stop("scanfill: not supported SEXP in 'm'");
  }
}

// convert poly x,y coordinates to img idx and push to queue
void polyidx (Rcpp::IntegerVector &Q,
              const R_len_t i_row,
              const R_len_t i_col,
              const R_len_t nrow) {
  if((i_row >= 1) && ((i_row - 1) < nrow) && (i_col >= 1)) queue_push(Q, (i_col - 1) * nrow + i_row - 1);
}

//' @title Polygon Contours
//' @name polycontour
//' @description
//' Determines contours indices
//' @param poly, a 2 (x and y) columns matrix defining the vertices of the polygon of interest. Eventually, a 4th (direction) column can be provided to speed-ud computation.
//' @param nrow, a R_len_t, number of rows of original image.
//' @param typ, a uint8_t, type of contour to return. 'poly' itself (typ=0), internal (typ=1) or external(typ=2).
//' @param inner, whether 'poly' is an inner or outer contour. Default is false.
//' @param conn, a uint8_t, desired connectedness. Default is 4, unless 8 is provided.
//' @return a Rcpp::IntegerVector of contour indices. 
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector polycontour (const Rcpp::IntegerMatrix poly,
                                 const R_len_t nrow = 0,
                                 const uint8_t typ = 2,
                                 const bool inner = false,
                                 const uint8_t conn = 4) {
  if(poly.ncol() < 2) Rcpp::stop("polycontour: 'poly' should have at least 4 columns");
  if(!(poly.nrow())) return 0;
  
  // close polygon;
  Rcpp::IntegerMatrix p = close_polygon(poly);
  
  // create queue;
  Rcpp::IntegerVector Q = queue_create(6 * p.nrow());
  
  // create vector to hold polygon indices
  for(R_len_t i = 0; i < p.nrow(); i++) polyidx(Q, p(i,1), p(i,0), nrow);
  Rcpp::IntegerVector W = Rcpp::no_init_vector(Q[0]);
  for(R_len_t i = 0; i < W.size(); i++) W[i] = Q[i + 1];
  if(typ == 0) return(W); 
  
  // empty queue
  Q[0] = 0;
  
  // define offsets
  Rcpp::IntegerVector dx = Rcpp::IntegerVector::create( 1, 1, 0,-1,-1,-1, 0, 1);
  Rcpp::IntegerVector dy = Rcpp::IntegerVector::create( 0, 1, 1, 1, 0,-1,-1,-1);
  
  if(poly.ncol() >= 4) {
    Rcpp::LogicalVector dd = Rcpp::LogicalVector::create(true,false,true,false,true,false,true,false);
    if(conn == 8) dd.fill(true);
    if(inner) {
      if(typ == 2) {
        short j = poly(poly.nrow() - 1, 3) % 8;
        j += 4;
        if(j > 7) j -= 8;
        for(R_len_t i = 0, ii = 1; i < p.nrow(); i++, ii++) {
          if(ii >= p.nrow()) ii = 0;
          for(short n = 0; n <= 7; n++) {
            if(++j > 7) j = 0;
            R_len_t x = p(i, 0) + dx[j];
            R_len_t y = p(i, 1) + dy[j];
            if((p(ii, 0) == x) && (p(ii, 1) == y)) { 
              j += 4;
              if(j > 7) j -= 8;
              break;
            }
            if(dd[j]) polyidx(Q, y, x, nrow);
          }
        } 
      } else {
        dx = Rcpp::IntegerVector::create(-1,-1,-1, 0,+1,+1,+1, 0);
        dy = Rcpp::IntegerVector::create(+1, 0,-1,-1,-1, 0,+1,+1);
        if(conn != 8) dd = !dd;
        for(R_len_t ii = p.nrow() - 1, i = ii - 1; ii >= 0 ; i--, ii--) {
          if(i < 0) i = p.nrow() - 1;
          short j = poly(ii, 3) % 8;
          j += 4;
          if(++j > 7) j -= 8;
          for(short n = 0; n <= 7; n++) {
            if(++j > 7) j = 0;
            R_len_t x = p(ii, 0) + dx[j];
            R_len_t y = p(ii, 1) + dy[j];
            if((p(i, 0) == x) && (p(i, 1) == y)) break;
            if(dd[j]) polyidx(Q, y, x, nrow);
          }
        }
      }
    } else {
      if(typ == 2) {
        short j = poly(poly.nrow() - 1, 3) % 8;
        j += 4;
        if(j > 7) j -= 8;
        for(R_len_t i = 0, ii = 1; i < p.nrow(); i++, ii++) {
          if(ii >= p.nrow()) ii = 0;
          for(short n = 0; n <= 7; n++) {
            if(--j < 0) j = 7;
            R_len_t x = p(i, 0) + dx[j];
            R_len_t y = p(i, 1) + dy[j];
            if((p(ii, 0) == x) && (p(ii, 1) == y)) {
              --j;
              j += 4;
              if(j > 7) j -= 8;
              break;
            }
            if(dd[j]) polyidx(Q, y, x, nrow);
          }
        }
        
        Rcpp::IntegerVector V = Rcpp::no_init_vector(Q[0]);
        for(R_len_t i = 0; i < V.size(); i++) V[i] = Q[i + 1];
        return Rcpp::setdiff(Rcpp::setdiff(V,W), polycontour(poly, nrow, 1, inner, 8));
      } else {
        short j = poly(poly.nrow() - 1, 3) % 8;
        j += 4;
        if(j > 7) j -= 8;
        for(R_len_t i = 0, ii = 1; i < p.nrow(); i++, ii++) {
          if(ii >= p.nrow()) ii = 0;
          for(short n = 0; n <= 7; n++) {
            if(++j > 7) j = 0;
            R_len_t x = p(i, 0) + dx[j];
            R_len_t y = p(i, 1) + dy[j];
            if((p(ii, 0) == x) && (p(ii, 1) == y)) {
              j += 4;
              if(j > 7) j -= 8;
              break;
            }
            if(dd[j]) polyidx(Q, y, x, nrow);
          }
        }
        
        Rcpp::IntegerVector V = Rcpp::no_init_vector(Q[0]);
        for(R_len_t i = 0; i < V.size(); i++) V[i] = Q[i + 1];
        return Rcpp::unique(V);
      }
    }
  }
  if(conn != 8) {
    dx = Rcpp::IntegerVector::create(+1, 0,-1, 0);
    dy = Rcpp::IntegerVector::create( 0,+1, 0,-1); 
  }
  
  for(R_len_t i = 0; i < p.nrow(); i++) {
    for(R_len_t j = 0; j < dx.size(); j++) { // test points around current polygon vertex with connectedness of 'conn'
      R_len_t x = p(i, 0) + dx[j], y = p(i, 1) + dy[j];
      if(ray_pnt_in_poly(x, y, p) ^ (typ == 1)) polyidx(Q, y, x, nrow);
    }
  }
  
  // return values in queue
  Rcpp::IntegerVector V = Rcpp::no_init(Q[0]);
  for(R_len_t i = 0; i < V.size(); i++) V[i] = Q[i + 1];
  return Rcpp::setdiff(V, W);
}

//' @title Polygon Seeds
//' @name polyseeds
//' @description
//' Determines seeds of polygon contours
//' @param poly, a 2 (x and y) columns matrix defining the vertices of the polygon of interest. Eventually, a 4th (direction) column can be provided to speed-ud computation.
//' @param nrow, a R_len_t, number of rows of original image.
//' @param inner, whether 'poly' is a inner or outer contour. Default is false.
//' @return a Rcpp::IntegerVector of seeds indices. 
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector polyseeds (const Rcpp::IntegerMatrix poly,
                               const R_len_t nrow = 0,
                               const bool inner = false) {
  if(poly.ncol() < 4) return polycontour(poly, nrow, 2, inner, 4);
  
  // close polygon
  Rcpp::IntegerMatrix p = close_polygon(poly);
  
  // init x,y coordinates and dx,dy offsets to be tested by ray_pnt_in_poly
  R_len_t x = 0, y = 0;
  Rcpp::IntegerVector dx1 = Rcpp::IntegerVector::create(0,-1, 0,-1, 0,+1, 0,+1, 0);
  Rcpp::IntegerVector dy2 = Rcpp::IntegerVector::create(0,+1, 0,-1, 0,-1, 0,+1, 0);
  Rcpp::LogicalVector dd  = Rcpp::LogicalVector::create(false,true,false,true,false,true,false,true);
  Rcpp::IntegerVector dy1(9), dx2 = dy1;
  if(inner == 1) { dx1 = -dx1; dy2 = -dy2; }
  
  Rcpp::IntegerVector Q = queue_create(2 * p.nrow());
  for(R_len_t i = 1; i < p.nrow(); i++) {
    int j = p(i, 3) % 8;
    if(dd[j]) { // only test when direction is diag
      x = p(i, 0) + dx1[j]; y = p(i, 1) + dy1[j];
      if(!((p(i - 1, 0) == x) && (p(i - 1, 1) == y))) if(ray_pnt_in_poly(x, y, p)) {
        polyidx(Q, y, x, nrow);
        if(inner) continue;
      }
      x = p(i, 0) + dx2[j]; y = p(i, 1) + dy2[j];
      if(!((p(i - 1, 0) == x) && (p(i - 1, 1) == y))) if(ray_pnt_in_poly(x, y, p)) polyidx(Q, y, x, nrow);
    }
  }
  // in case poly is a rectangle we add seed from diagonal offsets from 1st point
  if(p.nrow()) { // if(p.nrow()) nrow() is at least 2 thanks to close_polygon so p(0,) and p(1,) exist
    for(int j = 1; j <= 7; j += 2) {
      x = p(0, 0) + dx1[j]; y = p(0, 1) + dy2[j];
      if(!((p(1, 0) == x) && (p(1, 1) == y))) if(ray_pnt_in_poly(x, y, p)) {
        polyidx(Q, y, x, nrow);
        break;
      }
    }
  }
  
  // return values in queue
  Rcpp::IntegerVector V = Rcpp::no_init(Q[0]);
  for(R_len_t i = 0; i < V.size(); i++) V[i] = Q[i + 1];
  return V;
}

//' polycheckNA_T template
//' used to check that there is no NA close to polygon border
template <int RTYPE>
bool polycheckNA_T (const Rcpp::IntegerMatrix poly,
                    Rcpp::Matrix<RTYPE> img,
                    const uint8_t typ,
                    const bool inner,
                    const short conn) {
  if(poly.ncol() < 2) Rcpp::stop("polycheckNA: 'poly' should have at least 2 columns");
  Rcpp::IntegerVector nbr = polycontour(poly, img.nrow(), typ, inner, conn);
  for(R_len_t i = 0; i < nbr.size(); i++) if((nbr[i] >= 0) && (nbr[i] < img.size())) if(traits::is_na<RTYPE>(img[nbr[i]])) return true;
  return false;
}
 
//' polyvoid template
template <int RTYPE>
void polyvoid_T (const Rcpp::IntegerMatrix poly,
                 Rcpp::Matrix<RTYPE> img,
                 const double border,
                 const double fill,
                 const double tol,
                 const bool inner,
                 const bool checkNA,
                 const bool polydraw) {
  if(poly.ncol() < 2) Rcpp::stop("polyvoid: 'poly' should have at least 2 columns");
  if(poly.nrow()) {
    // check that 4-connected points outside 'poly' are not NA to avoid propagation outside polygon
    if(checkNA) if(polycheckNA_T(poly, img, 1, inner, 4)) Rcpp::stop("polycheckNA: [NA] found outside 'poly'");
    
    // a temporary polygon border is drawn with NA
    R_len_t mat_r = img.nrow();
    R_len_t mat_c = img.ncol();
    for(R_len_t i_row = 0; i_row < poly.nrow(); i_row++) {
      if(poly(i_row, 1) >= 1 && poly(i_row, 1) <= mat_r &&
         poly(i_row, 0) >= 1 && poly(i_row, 0) <= mat_c) {
        img(poly(i_row, 1) - 1, poly(i_row, 0) - 1) = NA_REAL;
      } else {
        Rcpp::stop("polyvoid: 'poly' is outside of 'img'");
      }
    }
    
    // find internal seeds
    Rcpp::IntegerVector seeds = polyseeds(poly, img.nrow(), inner);
    
    // fill polygon
    // when 'tol' is NA polygon is temporarily filled with NA so Z is used to record modified pix
    // afterwards Z is used to fill 'img' with original 'fill' input value within 'tol'
    Rcpp::IntegerVector Z = queue_create();
    if(polydraw && !traits::is_na<14>(tol)) {
      Rcpp::Matrix<RTYPE> img2 = Rcpp::clone(img);
      for(R_len_t i = 0; i < seeds.size(); i++) scanfill(img2, seeds[i], fill, NA_REAL, Z);
      for(R_len_t i = 1; i <= Z[0]; i++) if(std::abs(img[Z[i]] - fill) <= tol) img[Z[i]] = fill;
    } else {
      for(R_len_t i = 0; i < seeds.size(); i++) scanfill(img, seeds[i], fill, tol, Z);
      if(traits::is_na<14>(tol)) for(R_len_t i = 1; i <= Z[0]; i++) img[Z[i]] = fill;
    }
    
    // finally, overwrite temp NA polygon border with original 'border' input value
    for(R_len_t i_row = 0; i_row < poly.nrow(); i_row++) img(poly(i_row, 1) - 1, poly(i_row, 0) - 1) = border; 
  }
}

//' @name polyvoid
//' @description
//' This function is designed to trace and fill polygon inside a matrix.
//' @param poly, a 2-columns matrix defining the locations (x and y and direction) of vertices of the polygon of interest. Eventually, a 4th (direction) column can be provided to speed-ud computation.
//' @param img, a Matrix to be filled.
//' @param border, a double used to trace polygon 'border'.
//' @param fill, a double used to fill polygon.
//' @param tol, a double, tolerance between fill color and connected pixels. Use NA, for filling every pixel inside 'poly'.
//' @param inner, whether 'poly' is a inner or outer contour.
//' @param checkNA, whether or not to check that 'poly' is not surrounded by NA.
//' @param polydraw, whether polyvoid is called from polydraw. Special case to fill inside polygon using 'tol'
//' @return nothing /!\ img is modified in place.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
void polyvoid (const Rcpp::IntegerMatrix poly,
               SEXP img,
               const double border,
               const double fill,
               const double tol,
               const bool inner,
               const bool checkNA,
               const bool polydraw) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return polyvoid_T(poly, as<Rcpp::LogicalMatrix>(img), border, fill, tol, inner, checkNA, polydraw);
  case INTSXP : return polyvoid_T(poly, as<Rcpp::IntegerMatrix>(img), border, fill, tol, inner, checkNA, polydraw);
  case REALSXP : return polyvoid_T(poly, as<Rcpp::NumericMatrix>(img), border, fill, tol, inner, checkNA, polydraw);
  default : Rcpp::stop("polyvoid: not supported SEXP in 'img'");
  }
}

//' @title Connected Region Flood Filling
//' @name cpp_floodfill
//' @description
//' Flood fills image region.
//' @param img, a NumericMatrix. The image to be modified.
//' @param markers, a NumericMatrix, It should be a matrix with at least 3 columns being "row", "col", and "value", respectively. It represents coordinates of the seeds to start filling 'img', with the new "value". Eventually, an additional column being "tolerance" can be provided.\cr
//' /!\ Note that "row" and "col" should be provided at C-level meaning 1st start at 0.
//' @return a NumericMatrix, the modified image.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_floodfill (const Rcpp::NumericMatrix img,
                                   const Rcpp::NumericMatrix markers) {
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  R_len_t MAX_SIZ = markers.nrow();
  R_len_t m_col = markers.ncol();
  if(m_col < 3) Rcpp::stop("hpp_floodfill: 'markers' should have least 3 columns");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 2)) Rcpp::stop("floodfill: 'markers' is too large");
  
  // loop over seeds
  Rcpp::NumericMatrix out = Rcpp::clone(img);
  Rcpp::IntegerVector Q = queue_create();
  for(R_len_t i = 0; i < markers.nrow(); i++) {
    R_len_t i_row = markers(i, 0);
    R_len_t i_col = markers(i, 1);
    if(i_row < 0 || i_row >= mat_r || i_col < 0 || i_col >= mat_c) continue;
    if(m_col >= 4 && traits::is_na<14>(markers(i, 3))) Rcpp::stop("tolerance should not be [NA] for markers[%i,4]", (i + 1));
    scanfill(out, i_col * mat_r + i_row, markers(i, 2), m_col >= 4 ? markers(i, 3) : 0.0, Q);
  }
  return out;
}

//' @title Contours Filling
//' @name cpp_fill
//' @description
//' This function is designed to fill contours.
//' @param ctl a List, containing contour tracing labeling, object of class `IFCip_ctl`
//' @param label an int corresponding to the label of desired set of contour to be filled.
//' Default is \code{0} to fill all sets of contours found.
//' @param i_border a bool, to whether or not draw inside contours if some were identified.
//' @param i_fill a bool, to whether or not fill inside contours if some were identified.
//' @param o_border a bool, to whether or draw external contours.
//' @param o_fill a bool, to whether or not fill external contours.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_fill (const List ctl,
                              const int label = 0,
                              const bool i_border = true,
                              const bool i_fill = true,
                              const bool o_border = true,
                              const bool o_fill = true) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) Rcpp::stop("hpp_fill: 'ctl' should be of class `IFCip_ctl`");
  
  // retrieve max number of sets of contours identified by ctl
  int label_max = ctl["nb_lab"];
  
  // check whether to fill only one set or every sets of contours
  int label_cur = label <= 0 ? label_max : label;
  
  // retrieve contours coordinates
  Rcpp::IntegerMatrix contours = as<Rcpp::IntegerMatrix>(ctl["contours"]);
  
  // retrieve dimension of original image used to determine contours
  Rcpp::IntegerVector dim = as<Rcpp::IntegerVector>(ctl["dim"]);
  
  // create out 
  Rcpp::IntegerVector Q = queue_create();
  R_len_t p = 0, j = 0;
  bool o = o_fill || o_border;
  bool i = i_fill || i_border || o;
  Rcpp::IntegerMatrix out(dim[0], dim[1]);
  
  if((i || o) && (label_cur <= label_max)) {
    // fill contours with desired label(s)
    label_max = label_cur == label ? label : 1;
    while(label_cur >= label_max) {
      if(o) {
        for(R_len_t k = contours.nrow() - 1; k >= 0; k--) {
          if((contours(k, 2) == label_cur) && (contours(k, 4) == 1)) {
            queue_push(Q, k);
            if(contours(k, 3) > 7) {
              Rcpp::IntegerMatrix poly = Rcpp::no_init_matrix(Q[0], 4);
              j = 0;
              while(Q[0] != 0) {
                p = queue_pop(Q);
                poly[j] = contours(p, 0);
                poly[poly.nrow() + j] = contours(p, 1);
                poly[2 * poly.nrow() + j] = contours(p, 2);
                poly[3 * poly.nrow() + j++] = contours(p, 3);
              }
              polyvoid(poly, out, o_border ? label_cur : 0, o_fill ? label_cur : 0, 0.0, false, false, false);
            }
          }
        }
      }
      if(i) {
        for(R_len_t k = contours.nrow() - 1; k >= 0; k--) {
          if((contours(k, 2) == label_cur) && (contours(k, 4) == 2)) {
            queue_push(Q, k);
            if(contours(k, 3) > 7) {
              Rcpp::IntegerMatrix poly = Rcpp::no_init_matrix(Q[0], 4);
              j = 0;
              while(Q[0] != 0) {
                p = queue_pop(Q);
                poly[j] = contours(p, 0);
                poly[poly.nrow() + j] = contours(p, 1);
                poly[2 * poly.nrow() + j] = contours(p, 2);
                poly[3 * poly.nrow() + j++] = contours(p, 3);
              }
              polyvoid(poly, out, i_border ? label_cur : 0, i_fill ? label_cur : 0, 0.0, true, false, false);
            }
          }
        }
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
//' @param o_border a bool, to whether or draw external contours.
//' @param o_fill a bool, to whether or not fill external contours.
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_fill_out (const List ctl,
                                  const bool o_border = true,
                                  const bool o_fill = true) {
  if(!Rf_inherits(ctl, "IFCip_ctl")) Rcpp::stop("hpp_fill_out: 'ctl' should be of class `IFCip_ctl`");
  
  // retrieve contours coordinates
  Rcpp::IntegerMatrix contours = as<Rcpp::IntegerMatrix>(ctl["contours"]);
  
  // retrieve dimension of original image used to determine contours
  Rcpp::IntegerVector dim = as<Rcpp::IntegerVector>(ctl["dim"]);
  
  // create out 
  Rcpp::IntegerMatrix out(dim[0], dim[1]);
  if(!(o_fill || o_border)) return out;
  
  // fill every external contours with desired label(s)
  Rcpp::IntegerVector Q = queue_create();
  int label_cur = ctl["nb_lab"];
  while(label_cur >= 1) {
    for(R_len_t k = contours.nrow() - 1; k >= 0; k--) {
      if((contours(k, 2) == label_cur) && (contours(k, 4) == 1)) {
        queue_push(Q, k);
        if(contours(k, 3) > 7) {
          Rcpp::IntegerMatrix poly = Rcpp::no_init_matrix(Q[0], 4);
          R_len_t j = 0;
          while(Q[0] != 0) {
            R_len_t p = queue_pop(Q);
            poly[j] = contours(p, 0);
            poly[poly.nrow() + j] = contours(p, 1);
            poly[2 * poly.nrow() + j] = contours(p, 2);
            poly[3 * poly.nrow() + j++] = contours(p, 3);
          }
          polyvoid(poly, out, o_border ? label_cur : 0, o_fill ? label_cur : 0, 0.0, false, false, false);
        }
      }
    }
    label_cur--;
  }
  return out;
}

//// Helpers for polydraw
//' Compute distance from point with coords x, y to line passing through points x1, y1 and x2, y2.
double comp_dist (const R_len_t x1, const R_len_t y1,
                  const R_len_t x2, const R_len_t y2,
                  const double x, const double y) {
  R_len_t dx = x2-x1;
  R_len_t dy = y2-y1;
  return std::abs(dx*(y1-y)-dy*(x1-x)) / std::sqrt(dx*dx+dy*dy);
}

// check whether 2 values are almost equal +/- eps
bool check_same (const double x, const double y, const double eps = 0.000000001) {
  return std::abs(x - y) <= eps;
}

// check whether a point with coords x,y is inside a bounding box of coords [xmin,ymin], [xmax,ymax]
bool check_in (const double x, const double y,
               const double xmax, const double ymax,
               const double ymin = 0.0,
               const double xmin = 0.0) {
  return x >= xmin && y >= ymin && x <= xmax && y <= ymax;
}

//' @title Crossing Points Computation
//' @name compute_cross
//' @description
//' Compute crossing point(s) from point x1, y1 to point x2, y2 at edges of inret.
//' @param x1,x2 x-coordinates.
//' @param y1,y2 y-coordinates.
//' @param inret a LogicalMatrix where line will be drawn.
//' @return an IntegerVector of crossing point(s) edge + position on the edge;
//' edge is 0 for crossing point at y=0, 
//' edge 1 at x=0,
//' edge 2 at y=nrow(inret), and
//' edge 3 at x=ncol(inret)
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector compute_cross(const R_len_t x1, const R_len_t y1,
                                  const R_len_t x2, const R_len_t y2,
                                  Rcpp::LogicalMatrix inret) {
  Rcpp::IntegerVector out = fifo_create(4); // /!\ be carefull to not resize
  R_len_t mat_r = inret.nrow() - 1, mat_c = inret.ncol() - 1;
  if((mat_r < 0) || (mat_c < 0)) return out;
  if(!traits::is_finite<13>(x1) || !traits::is_finite<13>(y1) ||
     !traits::is_finite<13>(x2) || !traits::is_finite<13>(y2)) return(out);
     
  // variables
  short cross = 0;
  Rcpp::LogicalVector d = Rcpp::no_init_vector(4);
  Rcpp::NumericVector pt1 = Rcpp::NumericVector::create(x1, y1);
  Rcpp::NumericVector pt2 = Rcpp::NumericVector::create(x2, y2);
  bool pt1_in = check_in(pt1[0], pt1[1], mat_c, mat_r);
  bool pt2_in = check_in(pt2[0], pt2[1], mat_c, mat_r);
  short K1 = -1, K2 = -1, K3 = -1, K4 = -1;                         // edge of crossing point 0 (y=0), 1 (x=0), 2 (y=mat_r), 3 (x=mat_c)
  
  // compute coordinates of points crossing edges
  Rcpp::NumericMatrix M = Rcpp::no_init_matrix(4,2);
  M[0] = pt_projection_x(pt1, pt2, 0.0);  M[4] =0.0;
  M[1] = 0.0;                             M[5] = pt_projection_y(pt1, pt2, 0.0);
  M[2] = pt_projection_x(pt1, pt2, mat_r);M[6] = mat_r;
  M[3] = mat_c;                           M[7] = pt_projection_y(pt1, pt2, mat_c);
  
  // register crossing points
  for(R_len_t k = 0; k < 4; k++) {
    d[k] = traits::is_finite<14>(pt_shortest(M(k, _), pt1, pt2)) && // check crossing point is inside [pt1,pt2] segment
      check_in(M(k, 0), M(k, 1), mat_c, mat_r);                     // check crossing point is inside image
    if(d[k]) {
      cross++;
      if(K1 < 0) { K1 = k; } else { if(K2 < 0) { K2 = k; } else { if(K3 < 0) { K3 = k; } else { K4 = k; } }  }
    }
  }
  
  // keep only 0, 1 or 2 crossing points
  if(cross == 3) {            // 3 crossing points, so 2 of them are identical (same corner), only 2 different crossing points are kept
    if(check_same(M(K1, 0), M(K2, 0)) && check_same(M(K1, 1), M(K2, 1))) K2 = K3;
    K3 = -1; cross = 2;
  }
  if(cross == 4) {            // 4 crossing points, so there are 2 corners, only 2 different (they can still be identical) crossing points are kept
    if(check_same(M(K1, 0), M(K2, 0)) && check_same(M(K1, 1), M(K2, 1))) {
      K2 = check_same(M(K1, 0), M(K3, 0)) && check_same(M(K1, 1), M(K3, 1)) ? K4 : K3;
    }
    K3 = -1; K4 = -1; cross = 2;
  }
  
  // store crossing point in returned fifo queue
  if(cross == 1) {            // only 1 crossing point
    if(!(pt1_in && pt2_in)) { // if both pt1 and pt2 are in it means that crossing point is either pt1 or pt2 so it is not returned
      fifo_add(out, K1); fifo_add(out, M(K1, K1 % 2));
    }
  }
  if(cross == 2) {            // 2 crossing points
    double d11 = pt_distance_fast(pt1, M(K1, _));
    double d12 = pt_distance_fast(pt1, M(K2, _));
    bool same11 = check_same(pt1[0], M(K1, 0)) && check_same(pt1[1], M(K1, 1));
    bool same21 = check_same(pt2[0], M(K1, 0)) && check_same(pt2[1], M(K1, 1));
    bool same12 = check_same(pt1[0], M(K2, 0)) && check_same(pt1[1], M(K2, 1));
    bool same22 = check_same(pt2[0], M(K2, 0)) && check_same(pt2[1], M(K2, 1));
    
    if(pt_distance_fast(M(K1,_), M(K2, _)) > 2) {                      // check that crossing points K1 and K2 are not +/-1 px identical
      if(d12 > d11) {                                                  // sort returned crossing point so that first-in copied to out is the closest one to pt1 (unless it is pt1)
        if(!same11) {fifo_add(out, K1); fifo_add(out, M(K1, K1 % 2));} // add K1 (closest to pt1), if it is not pt1
        if(!same22) {fifo_add(out, K2); fifo_add(out, M(K2, K2 % 2));} // add K2 (farthest from pt1), if it is not pt2
      } else {
        if(!same12) {fifo_add(out, K2); fifo_add(out, M(K2, K2 % 2));} // add K2 (closest to pt1), if it is not pt1
        if(!same21) {fifo_add(out, K1); fifo_add(out, M(K1, K1 % 2));} // add K1 (farthest from pt1), if it is not pt2
      }
    } else {
      if(pt1_in || pt2_in) if(!same11 && !same21) {fifo_add(out, K1); fifo_add(out, M(K1, K1 % 2));}
    }
  }
  return out;
}

//' @title Trace Edge
//' @name trace_edge
//' @description
//' Trace segment from point x1, y1 to point x2, y2  at edge of inret.
//' @param x1,x2 x-coordinates.
//' @param y1,y2 y-coordinates.
//' @param inret a LogicalMatrix where line will be drawn.
//' @return nothing /!\ inret is modified in place.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
void trace_edge (const R_len_t x1, const R_len_t y1,
                 const R_len_t x2, const R_len_t y2,
                 Rcpp::LogicalMatrix inret) {
  R_len_t mat_r = inret.nrow() - 1, mat_c = inret.ncol() - 1;
  if((mat_r < 0) || (mat_c < 0)) return ;
  if(!traits::is_finite<13>(x1) || !traits::is_finite<13>(y1) ||
     !traits::is_finite<13>(x2) || !traits::is_finite<13>(y2)) return ;
     
  if(x1 == x2 && y1 != y2) {
    if(x1 >= 0 && x1 <= mat_c) {
      R_len_t beg = y1, end = y2;
      if(y2 < y1) { beg = y2; end = y1; }
      beg = std::max(beg, 0);
      end = std::min(end, mat_r);
      for(R_len_t y = beg ; y <= end; y++) inret(y, x1) = !inret(y, x1); // trick to avoid tracing edge twice
    }
  }
  if(y1 == y2 && x1 != x2) {
    if(y1 >= 0 && y1 <= mat_r) {
      R_len_t beg = x1, end = x2;
      if(x2 < x1) { beg = x2; end = x1; }
      beg = std::max(beg, 0);
      end = std::min(end, mat_c);
      for(R_len_t x = beg ; x <= end; x++) inret(y1, x) = !inret(y1, x);
    }
  }
}

//' @title Trace Line
//' @name trace_seg
//' @description
//' Trace segment from point x1, y1 to point x2, y2 in inret.
//' @param x1,x2 x-coordinates.
//' @param y1,y2 y-coordinates.
//' @param inret a LogicalMatrix where segment will be drawn.
//' @param k an IntegerVector of crossing points values.
//' @return Nothing but /!\ inret is modified in place.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
void trace_seg (const R_len_t x1, const R_len_t y1,
                const R_len_t x2, const R_len_t y2,
                Rcpp::LogicalMatrix inret,
                Rcpp::IntegerVector K) {
  R_len_t mat_r = inret.nrow() - 1, mat_c = inret.ncol() - 1,
    X1 = x1, Y1 = y1, X2 = x2, Y2 = y2;
  if((mat_r < 0) || (mat_c < 0)) return ;
  if(!traits::is_finite<13>(x1) || !traits::is_finite<13>(y1) ||
     !traits::is_finite<13>(x2) || !traits::is_finite<13>(y2)) return ;
     
  // variables
  short cross = K[0] / 2;
  Rcpp::NumericMatrix M(4, 4);
  Rcpp::IntegerVector edg = Rcpp::IntegerVector::create(0, 0, mat_r, mat_c);
  short K1 = -1, K2 = -1;
  while(K[0] != 0) { // get K1/K2 and positions of crossing points from K
    short k = fifo_pop(K);
    M(k, k % 2) = fifo_pop(K);
    M(k, !(k % 2)) = edg[k];
    if(K1 < 0) { K1 = k; } else { K2 = k; }
  }
  
  Rcpp::NumericVector pt1 = Rcpp::NumericVector::create(x1, y1);
  Rcpp::NumericVector pt2 = Rcpp::NumericVector::create(x2, y2);
  bool pt1_in = check_in(pt1[0], pt1[1], mat_c, mat_r);
  bool pt2_in = check_in(pt2[0], pt2[1], mat_c, mat_r);
  
  // restrict pt1 and pt2 coordinates to image range (e.g. there is no need to start tracing at x1 = -1e9)
  if(!pt1_in && !pt2_in) {        // pt1 and pt2 are outside of image
    if(cross == 2) {              // there should be 2 crossing points 
      pt1 = M(K1, _);
      pt2 = M(K2, _);
    } else {                      // otherwise it means that there is no segment to draw
      return ;
    }
  } else {
    if(cross == 1) {              // one crossing point was found
      if(!(pt1_in && pt2_in)) {   // if both pt1 and pt2 are in it means that crossing point is one of these, so no need to replace
        if(pt1_in) pt2 = M(K1, _);// pt1 is in, so pt2 is out and pt2 is replaced by crossing point
        if(pt2_in) pt1 = M(K1, _);// pt2 is in, so pt1 is out and pt1 is replaced by crossing point
      }
    }
    if(cross == 2) {              // 2 crossing points found
      if(!(pt1_in && pt2_in)) {   // if both pt1 and pt2 are out we replace them
        if(!pt1_in) pt1 = pt_distance_fast(pt2, M(K1, _)) > pt_distance_fast(pt2, M(K2, _)) ? M(K1, _) : M(K2, _);
        if(!pt2_in) pt2 = pt_distance_fast(pt1, M(K1, _)) > pt_distance_fast(pt1, M(K2, _)) ? M(K1, _) : M(K2, _);
      } else {                    // one can can be exactly crossing point and not the other, so we need to replace but only if it is out
        bool same11 = check_same(pt1[0], M(K1, 0)) && check_same(pt1[1], M(K1, 1));
        bool same21 = check_same(pt2[0], M(K1, 0)) && check_same(pt2[1], M(K1, 1));
        bool same12 = check_same(pt1[0], M(K2, 0)) && check_same(pt1[1], M(K2, 1));
        bool same22 = check_same(pt2[0], M(K2, 0)) && check_same(pt2[1], M(K2, 1));
        if(!pt2_in) if(same11 || same12) pt2 = same11 ? M(K2, _) : M(K1, _);
        if(!pt1_in) if(same21 || same22) pt1 = same21 ? M(K2, _) : M(K1, _);
      }
    }
  }
  
  // sort new values for pt1 and pt2 so that pt1 will be the closest one to (x1, y1)
  if(pt_distance_fast(Rcpp::NumericVector::create(x1, y1), pt2) < pt_distance_fast(Rcpp::NumericVector::create(x1, y1), pt1)) {
    Rcpp::NumericVector tmp = clone(pt1);
    pt1 = Rcpp::clone(pt2);
    pt2 = Rcpp::clone(tmp);
  }
  
  // start tracing line
  X1 = pt1[0]; Y1 = pt1[1]; X2 = pt2[0]; Y2 = pt2[1];
  R_len_t y = Y1, x = X1;
  double offsetx = 0.5, offsety = 0.5;
  if(X2 < x) offsetx = -offsetx;
  if(Y2 < y) offsety = -offsety;
  if(x >= 0 && y >= 0 && x <= mat_c && y <= mat_r) inret(y, x) = true;
  if(X1 == X2 && Y1 == Y2) return ;
  unsigned short count = 1;
  while(!((std::abs(y - Y2) < 0.5) && (std::abs(x - X2) < 0.5))) {
    if((count++ % 10000) == 0) {
      count = 1;
      Rcpp::checkUserInterrupt();
    }
    Rcpp::NumericVector dist = Rcpp::NumericVector::create(comp_dist(X1, Y1, X2, Y2, x + offsetx, y + offsety),
                                                           comp_dist(X1, Y1, X2, Y2, x + offsetx, y),
                                                           comp_dist(X1, Y1, X2, Y2, x, y + offsety));
    switch(Rcpp::which_min(dist))
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
    if(x >= 0 && y >= 0 && x <= mat_c && y <= mat_r) inret(y, x) = true;
  }
  return ;
}

// dequeue points to draw edge segments
// Q: queue with at least 4 values, edge1, position1, edge2, position2 of points pt1 and pt2
// M: image
// invert: used for inverting pt1 and pt2
// [[Rcpp::export(rng = false)]]
void dq_edge (IntegerVector &Q,
              LogicalMatrix &inret,
              short invert = 0) {
  if((Q.size() == 0) || (Q[0] < 4)) return ; // safer
  Rcpp::IntegerVector edg = Rcpp::IntegerVector::create(0, 0, inret.nrow() - 1, inret.ncol() - 1);
  short k1, k2;
  R_len_t idx1, idx2;
  if(invert) {
    k2 = fifo_pop(Q); idx2 = fifo_pop(Q);
    k1 = fifo_pop(Q); idx1 = fifo_pop(Q);
  } else {
    k1 = fifo_pop(Q); idx1 = fifo_pop(Q);
    k2 = fifo_pop(Q); idx2 = fifo_pop(Q); 
  }
  if(k1 == k2) { // pt1 and pt2 are at same edge
    if(k1 % 2) {
      trace_edge(edg[k1], idx1, edg[k1], idx2, inret);
    } else {
      trace_edge(idx1, edg[k1], idx2, edg[k1], inret);
    }
  } else {       // pt1 and pt2 are 2 different edges so we need to draw portion of edge in between them
    Rcpp::NumericVector pt1 = Rcpp::no_init_vector(2);
    Rcpp::NumericVector pt1o= Rcpp::no_init_vector(2);
    Rcpp::NumericVector pt2 = Rcpp::no_init_vector(2);
    Rcpp::NumericVector pt2o= Rcpp::no_init_vector(2);
    if(k2 % 2) {
      pt2[0] = edg[k2]; pt2[1] = idx2;
      pt2o[0] = edg[k2]; pt2o[1] = edg[k1];
      pt2o[1] = pt2o[1] == 0 ? 1 : pt2o[1] - 1; // pt2o is +/- 1 at corner to prevent drawing corner twice
    } else {
      pt2[0] = idx2; pt2[1] = edg[k2];
      pt2o[1] = edg[k2]; pt2o[0] = edg[k1];
      pt2o[0] = pt2o[0] == 0 ? 1 : pt2o[0] - 1; // pt2o is +/- 1 at corner to prevent drawing corner twice 
    }
    if((std::abs(k1 - k2) == 1) || ((k1 == (short) 3) && (k2 == (short) 0)) || ((k1 == (short) 0) && (k2 == (short) 3))) { // both edges are connected
      if(k1 % 2) {
        pt1[0] = edg[k1]; pt1[1] = idx1;
        pt1o[0] = edg[k1]; pt1o[1] = edg[k2];
      } else {
        pt1[0] = idx1; pt1[1] = edg[k1];
        pt1o[1] = edg[k1]; pt1o[0] = edg[k2];
      }
      trace_edge(pt1[0], pt1[1], pt1o[0], pt1o[1], inret);
      trace_edge(pt2[0], pt2[1], pt2o[0], pt2o[1], inret);
    } else {  // edges are parallel
      // TODO the problem is which edge to go between k1 vs k2. This question is handled by 'invert'
      if(k1 % 2) {
        pt1[0] = edg[k1]; pt1[1] = idx1;
        pt1o[0] = edg[k1]; pt1o[1] = edg[k1];
      } else {
        pt1[0] = idx1; pt1[1] = edg[k1];
        pt1o[1] = edg[k1]; pt1o[0] = edg[k1];
      }
      trace_edge(pt1[0], pt1[1], pt1o[0], pt1o[1], inret);
      trace_edge(pt2[0], pt2[1], pt2o[0], pt2o[1], inret);
      trace_edge(pt1o[0], pt2o[1], pt1o[1], pt2o[0], inret);
    }
  }
}

// dequeue crossing points to draw edge segments
// [[Rcpp::export(rng = false)]]
void polyedge (const IntegerMatrix M,
               const bool invert,
               LogicalMatrix &inret) {
  Rcpp::IntegerVector Q = fifo_create(8, NA_INTEGER);      // hold crossing points to trace the edge (load and unload every 2 points)
  Rcpp::IntegerVector H = fifo_create(2, NA_INTEGER);      // hold the very 1st crossing point of a series of segments with 2 crossing points
  Rcpp::IntegerVector B = fifo_create(6, NA_INTEGER);      // hold series of 2 crossing points (load and unload every 2 points but start unloading at 2).
  bool one = false;
  bool swap = false;
  for(R_len_t i = 0; i < M.nrow(); i++) {
    Rcpp::IntegerVector X = M(i, _);
    if(X[0] == 4) {                                        // 2 crossing points
      swap = !swap;
      if(H[0] == 0) {                                      // series start
        fifo_add(H, fifo_pop(X)); fifo_add(H, fifo_pop(X));// hold very 1st crossing point out of 2 in the series
        fifo_add(B, fifo_pop(X)); fifo_add(B, fifo_pop(X));// start banking X in B
      } else {                                             // series follow (unload X and dequeue)
        while(X[0] != 0) {fifo_add(B, fifo_pop(X));}       // unload X into B
        while(B[0] >= 6) dq_edge(B, inret, true);          // dequeue B
      }
    } else {                                               // no or 1 crossing point
      if(X[0] == 2) one = true;
      while(B[0] != 0) fifo_add(Q, fifo_pop(B));           // load remaining point of 2 crossing point series (if any) in Q
      while(H[0] != 0) fifo_add(Q, fifo_pop(H));           // load H (if any) into Q
      while(X[0] != 0) fifo_add(Q, fifo_pop(X));           // copy remaining crossing point to Q
    }
    while(Q[0] >= 4) dq_edge(Q, inret, true);              // dequeue Q
  }
  while(B[0] != 0) fifo_add(Q, fifo_pop(B));               // load remaining point of 2 crossing point series (if any) in Q
  while(H[0] != 0) fifo_add(Q, fifo_pop(H));               // load H (if any) into Q
  while(Q[0] >= 4) dq_edge(Q, inret, swap & one);          // dequeue Q
}

//' @title Polygon Mask
//' @name polymask
//' @description
//' Create a logical matrix of polygon contour.
//' @param poly a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param nrow, R_len_t of desired returned matrix rows.
//' @param ncol, R_len_t of desired returned matrix columns.
//' @param edge, a bool whether to close 'poly' at edges
//' @return a LogicalMatrix
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix polymask (const Rcpp::IntegerMatrix poly,
                              const R_len_t nrow = 0,
                              const R_len_t ncol = 0,
                              const bool edge = false) {
  Rcpp::LogicalMatrix out(std::max(nrow, 0), std::max(ncol, 0));
  if((nrow <= 0) || (ncol <= 0)) return out;
  if(poly.ncol() < 2) Rcpp::stop("'poly' should have at least 2 columns");
  
  // close polygon
  Rcpp::IntegerMatrix p = close_polygon(poly);
  R_len_t pts_row = p.nrow();
  
  // matrix to store crossing points happening when tracing lines between poly vertices
  Rcpp::IntegerMatrix M(pts_row - 1, 6); // /!\ trace_seg will return a fifo_queue of size = 6 (with at most 4 values stored + ele 0 describing number of value(s) stored + last ele describing index of 1st value stored if any)
  Rcpp::IntegerMatrix MM(pts_row - 1, 6);
  
  // compute crossing points
  for(R_len_t i = 0; i < pts_row - 1; i++) {
    M(i, _) = compute_cross(p[i] - 1, p[i + pts_row] - 1, p[i + 1] - 1, p[i + 1 + pts_row] - 1, out);
  }
  
  // draw edge segments between crossing point
  if(edge) polyedge(M, true, out);
  
  //draw segments between polygon vertices
  for(R_len_t i = 0; i < pts_row - 1; i++) {
    Rcpp::IntegerVector K = M(i, _);
    trace_seg(p[i] - 1, p[i + pts_row] - 1, p[i + 1] - 1, p[i + 1 + pts_row] - 1, out, K);
  }
  
  return out;
}

//' @title Polygon Drawing
//' @name cpp_polydraw
//' @description
//' This function is designed to trace and fill polygon inside a matrix.
//' @param poly, a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param border, a double used to trace polygon border.
//' @param fill, a double used to fill polygon.
//' @param tol, a double, tolerance between fill color and connected pixels. Use \code{NA}, for filling every pixel inside 'poly'.
//' @param edge, a bool whether to close 'poly' at 'mat_' edges. Default is \code{false}. Closing 'poly' at 'edge' is experimental and may fail.
//' @param mat_, a NumericMatrix to be filled.\cr
//' When 'mat_' is provided 'poly' will be drawn in 'mat_' if its vertices are within 'mat_' dimensions and if there is no \code{NA} directly (4-connectedness) connected to internal part of 'poly' border.\cr
//' /!\ Note that filling will not propagate on \code{NA}, \code{NaN} values, unless 'tol' is \code{NA}.
//' @return copy of 'mat_' with 'poly' or a new matrix with 'poly'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_polydraw (const Rcpp::IntegerMatrix poly,
                                  const double border = 1.0,
                                  const double fill = 1.0,
                                  const double tol = 0.0,
                                  const bool edge = false,
                                  const Rcpp::Nullable<Rcpp::NumericMatrix> mat_ = R_NilValue) {
  if(poly.ncol() < 0) Rcpp::stop("hpp_polydraw: 'poly' should have at least 2 columns");
  R_len_t pnrow = 0, pncol = 0;
  for(R_len_t i = 0; i < poly.nrow(); i++) {
    if(poly(i,0) > pncol) pncol = poly(i,0);
    if(poly(i,1) > pnrow) pnrow = poly(i,1);
  }
  Rcpp::NumericMatrix out(pnrow, pncol);
  if(mat_.isNotNull()) {
    Rcpp::NumericMatrix mat(mat_.get());  
    out = Rcpp::clone(mat);
  }
  Rcpp::List ctl = hpp_ctl(polymask(poly, out.nrow(), out.ncol(), edge));
  Rcpp::IntegerMatrix full_poly = ctl["contours"];
  Rcpp::IntegerVector Q = queue_create();
  R_len_t j = 0, p = 0, label_cur = 0, label_max = ctl["nb_lab"];
  while((++label_cur) <= label_max) {
    for(R_len_t k = full_poly.nrow() - 1; k >= 0; k--) {
      if((full_poly(k, 2) == label_cur) && (full_poly(k, 4) == 1)) {
        queue_push(Q, k);
        if(full_poly(k, 3) > 7) {
          Rcpp::IntegerMatrix poly_ = Rcpp::no_init_matrix(Q[0], 4);
          j = 0;
          while(Q[0] != 0) {
            p = queue_pop(Q);
            poly_[j] = full_poly(p, 0);
            poly_[poly_.nrow() + j] = full_poly(p, 1);
            poly_[2 * poly_.nrow() + j] = full_poly(p, 2);
            poly_[3 * poly_.nrow() + j++] = full_poly(p, 3);
          }
          polyvoid(poly_, out, border, fill, tol, false, true, true);
        }
      }
    }
  }
  return out;
}

#endif
