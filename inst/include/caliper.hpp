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

#ifndef IFCIP_CALIPER_HPP
#define IFCIP_CALIPER_HPP

#include <Rcpp.h>
#include "fifo.hpp"
#include "geometry.hpp"
using namespace Rcpp;

// helper
R_len_t get_next (R_len_t idx,
                  const Rcpp::NumericMatrix pts) {
  idx++;
  if(idx >= pts.nrow()) idx = 0;
  return idx;
}

// helper
double get_area (const R_len_t idx1,
                 const R_len_t idx2,
                 const R_len_t idx3,
                 const Rcpp::NumericMatrix pts) {
  return area(pts(idx1, Rcpp::_), pts(idx2, Rcpp::_), pts(idx3, Rcpp::_));
}

//' @title Antipodal Pairs of Convex Hull
//' @name cpp_antipodalpairs
//' @description
//' Computes antipodal pairs of a convex polygon 
//' @param pts a 2-column matrix defining the locations (x and y coordinates, respectively) of points.
//' It has to be an object of class `IFCip_convexhull`
//' @source Adaptation from \url{https://escholarship.mcgill.ca/concern/theses/fx719p46g} in Computational geometry with the rotating calipers authored by Pirzadeh, Hormoz under supervision of Toussaint, Godfried T. at McGill University.
//' @return an IntegerVector of antipodal pairs of the convex input polygon. Note that this vector of indices is in C so 1st start at 0; add 1 to use it in R.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_antipodalpairs (const Rcpp::NumericMatrix pts) {
  if(!Rf_inherits(pts, "IFCip_convexhull")) {
    Rcpp::stop("hpp_antipodalpairs: 'pts' should be of class `IFCip_convexhull`");
  }
  if(pts.nrow() < 3) {
    Rcpp::stop("hpp_antipodalpairs: 'pts' should have at least 3 rows");
  }
  R_len_t q = 0, p0 = 0, p = pts.nrow() - 1;
  Rcpp::IntegerVector out_p = Rcpp::no_init_vector(2 * pts.nrow() + 1);
  Rcpp::IntegerVector out_q = Rcpp::no_init_vector(2 * pts.nrow() + 1);
  out_p[0] = 0;
  out_q[0] = 0;
  
  while(get_area(p, get_next(p, pts), get_next(q, pts), pts) > get_area(p, get_next(p, pts), q, pts)) q = get_next(q, pts);
  R_len_t q0 = q;
  
  while (q != p0) {
    p = get_next(p, pts);
    queue_push(out_p, p);
    queue_push(out_q, q);;
    while(get_area(p, get_next(p, pts), get_next(q, pts), pts) > get_area(p, get_next(p, pts), q, pts)) {
      q = get_next(q, pts);
      if((p != q0) || (q != p0)) {
        queue_push(out_p, p);;
        queue_push(out_q, q);;
      } else {
        break;
      }
    }
    // parallel case
    if(get_area(p, get_next(p, pts), get_next(q, pts), pts) == get_area(p, get_next(p, pts), q, pts)) {
      if((p != q0) || (q != (pts.nrow() - 1))) {
        queue_push(out_p, p);;
        queue_push(out_q, get_next(q, pts));
      } else {
        queue_push(out_p, get_next(p, pts));
        queue_push(out_q, q);;
        break;
      }
    }
  }
  
  Rcpp::IntegerMatrix out(out_p[0], 2);
  for(R_len_t i = 0; i < out_p[0]; i++) {
    out(i, 0) = out_p[i + 1];
    out(i, 1) = out_q[i + 1];
  }
  out.attr("class") = "IFCip_antipodal";
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_antipodalpairs_width (const Rcpp::IntegerMatrix pairs,
                                              const Rcpp::NumericMatrix pts,
                                              const double scale = 1.0) {
  if(!Rf_inherits(pts, "IFCip_convexhull")) {
    Rcpp::stop("hpp_antipodalpairs: 'pts' should be of class `IFCip_convexhull`");
  }
  if(!Rf_inherits(pairs, "IFCip_antipodal")) {
    Rcpp::stop("hpp_antipodalpairs: 'pts' should be of class `IFCip_antipodal`");
  }
  Rcpp::NumericMatrix out(pts.nrow(), 4);
  for(R_len_t vertex = 0; vertex < pts.nrow(); vertex++) {
    Rcpp::IntegerVector V1 = Rcpp::no_init_vector(2 * pairs.nrow() + 1);
    V1[0] = 0;
    for(R_len_t anti = 0; anti < pairs.nrow(); anti++) {
      if(pairs(anti, 0) == vertex) queue_push(V1, pairs(anti, 1));
      if(pairs(anti, 1) == vertex) queue_push(V1, pairs(anti, 0));
    }
    if(V1[0] < 2) {
      out(vertex, 0) = vertex;
      out(vertex, 1) = V1[1];
      out(vertex, 2) = NA_REAL;
      out(vertex, 3) = R_PosInf;
    } else {
      Rcpp::NumericVector V2(V1[0] - 1);
      for(R_len_t i = 0; i < V2.size(); i++) {
        V2[i] = pt_shortest(pts(vertex, Rcpp::_), pts(V1[i + 1], Rcpp::_), pts(V1[i + 2], Rcpp::_), scale);
      }
      R_len_t m = which_min(V2);
      out(vertex, 0) = vertex;
      out(vertex, 1) = V1[m + 1];
      out(vertex, 2) = V1[m + 2];
      out(vertex, 3) = V2[m];
    }
  }
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("pt", "S1", "S2", "d");
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_antipodalpairs_feat (const Rcpp::IntegerMatrix pairs,
                                             const Rcpp::NumericMatrix pts,
                                             const double scale = 1.0) {
  if(!Rf_inherits(pts, "IFCip_convexhull")) {
    Rcpp::stop("hpp_antipodalpairs: 'pts' should be of class `IFCip_convexhull`");
  }
  if(!Rf_inherits(pairs, "IFCip_antipodal")) {
    Rcpp::stop("hpp_antipodalpairs: 'pts' should be of class `IFCip_antipodal`");
  }
  Rcpp::NumericVector out =  Rcpp::NumericVector::create(R_NegInf, R_PosInf, 0.0, 0.0, 0.0, 0.0);
  R_len_t vertex = 0;
  for(; vertex < pts.nrow() - 1; vertex++) {
    out[4] += pts(vertex, 0);
    out[5] += pts(vertex, 1);
    out[3] += pt_distance_accu(pts(vertex, Rcpp::_), pts(vertex + 1, Rcpp::_), scale);
  }
  out[4] += pts(vertex, 0);
  out[5] += pts(vertex, 1);
  out[3] += pt_distance_accu(pts(vertex, Rcpp::_), pts(0, Rcpp::_), scale);
  
  out[4] /= pts.nrow();
  out[5] /= pts.nrow();
  
  Rcpp::NumericVector dis(pairs.nrow());
  for(R_len_t anti = 0; anti < pairs.nrow(); anti++) {
    if((pairs(anti, 0) >= pts.nrow()) || (pairs(anti, 1) >= pts.nrow())) {
      Rcpp::stop("hpp_antipodalpairs_feat: 'pairs' and 'pts' are incompatible");
    }
    dis[anti] = pt_distance_accu(pts(pairs(anti, 0), Rcpp::_), pts(pairs(anti, 1), Rcpp::_), scale);
    if(out[0] < dis[anti]) out[0] = dis[anti];
    if(out[1] > dis[anti]) out[1] = dis[anti];
  }
  out[2] = out[0] / out[1];
  out.attr("names") = Rcpp::CharacterVector::create("Height", "Width", "Elongatedness", "convex perimeter", "convex cx", "convex cy");
  return out;
}

//' @title Bounding Box of Convex Hull
//' @name cpp_bbox
//' @description
//' Computes features from a Convex Hull 
//' @param pts a 2-column matrix defining the locations (x and y coordinates, respectively) of points.
//' It has to be an object of class `IFCip_convexhull`
//' @param scale a double used to scale the returned values.
//' @return a NumericVector of features from convex hull.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_bbox (const Rcpp::NumericMatrix pts,
                              const double scale = 1.0) {
  return hpp_antipodalpairs_feat(hpp_antipodalpairs(pts), pts, scale);
}

#endif
