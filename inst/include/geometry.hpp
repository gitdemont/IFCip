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

#ifndef IFCIP_GEOMETRY_HPP
#define IFCIP_GEOMETRY_HPP

#include <Rcpp.h>
using namespace Rcpp;

// ptx is supposed to be an Rcpp::NumericVector where
// -element 0 is x
// -element 1 is y

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_translate(const Rcpp::NumericVector pt1,
                                 const double x = 0.0,
                                 const double y = 0.0) {
  return Rcpp::NumericVector::create(_["x"] = pt1[0] + x,
                                      _["y"] = pt1[1] + y);
}

// [[Rcpp::export(rng = false)]]
double to_radians(const double angle = 0.0) {
  return angle * M_PI / 180;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_rotate_ori(const Rcpp::NumericVector pt1,
                                  const double theta = 0.0) {
  double s = std::sin(theta);
  double c = std::cos(theta);
  return Rcpp::NumericVector::create(_["x"] = pt1[0] * c - pt1[1] * s,
                                      _["y"] = pt1[1] * c + pt1[0] * s);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_rotate(const Rcpp::NumericVector pt1,
                              const Rcpp::NumericVector pt2 = Rcpp::NumericVector::create(0,0),
                              const double theta = 0.0) {
  return pt_translate(pt_rotate_ori(pt_translate(pt1, -pt2[0], -pt2[1]), theta), pt2[0], pt2[1]);
}

// helper to determine pt_slope of a line given 2 points 
// [[Rcpp::export(rng = false)]]
double pt_slope(const Rcpp::NumericVector pt1,
                const Rcpp::NumericVector pt2) {
  return (pt1[1] - pt2[1]) / (pt1[0] - pt2[0]);
}

// helper to determine pt_angle of a line given 2 points
// [[Rcpp::export(rng = false)]]
double pt_angle(const Rcpp::NumericVector pt1,
                const Rcpp::NumericVector pt2) {
  return 180 * std::atan(pt_slope(pt1, pt2[1])) / M_PI;
}

// helper to determine x coordinate of a line passing from
// points pt1 and pt2 and intersecting with line y
// [[Rcpp::export(rng = false)]]
double pt_projection_x(const Rcpp::NumericVector pt1,
                       const Rcpp::NumericVector pt2,
                       const double y = 0.0) { // cross with y
  if(pt1[0] == pt2[0]) return pt2[0];
  double m = pt_slope(pt1, pt2);
  double p = pt2[1] - m * pt2[0];
  return (y - p) / m;
}

// helper to determine y coordinate of a line passing from
// points pt1 and pt2 and intersecting with line x
// [[Rcpp::export(rng = false)]]
double pt_projection_y(const Rcpp::NumericVector pt1,
                       const Rcpp::NumericVector pt2,
                       const double x = 0.0) { // cross with x
  if(pt1[1] == pt2[1]) return pt2[1];
  double m = pt_slope(pt1, pt2);
  double p = pt2[1] - m * pt2[0];
  return m * x + p;
}

// helper to determine point central symmetry of pt1 against pt2.
// It returns pt' x and y coordinates of pt1 reflected against pt2.
// As a result pt2 will be in the middle of the segment formed by pt' and pt1.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_inversion(const Rcpp::NumericVector pt1,
                                 const Rcpp::NumericVector pt2) {
  double xo = pt2[0], yo = pt2[1], xx = pt1[0], yy = pt1[1];
  if(yo == yy) return Rcpp::NumericVector::create(_["x"] = 2 * xo - xx, _["y"] = yy);
  if(xo == xx) return Rcpp::NumericVector::create(_["x"] = xx, _["y"] = 2 * yo - yy);
  double m = pt_slope(pt1, pt2);
  double a = -1 / m;
  double c = yo + xo / m;
  double d = 2 * (xx + a * (yy - c) ) / (1 + a * a);
  return Rcpp::NumericVector::create(_["x"] = d - xx,
                                     _["y"] = d * a - yy + 2 * c);
}

// pt_inversion is faster than pt_inversion2
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_inversion2(const Rcpp::NumericVector pt1,
                                  const Rcpp::NumericVector pt2) {
  return pt_translate(pt_rotate_ori(pt_translate(pt1, -pt2[0], -pt2[1]), M_PI), pt2[0], pt2[1]);
}

// helper to determine reflection of point pt1 by the line defined by pt2 and pt3.
// It returns pt' x and y coordinates of pt1 pt_reflection accross line pt2 - pt3.
// As a result pt'- pt1 and pt2 - pt3 segments will be perpendicular.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_reflection(const Rcpp::NumericVector pt1,  // pt of interest
                                  const Rcpp::NumericVector pt2,  // 1st coord of the segment
                                  const Rcpp::NumericVector pt3) {// 2nd coord of the segment
  // line ax + by + c; b = -1
  if(pt2[0] == pt3[0]) {
    return Rcpp::NumericVector::create(_["x"] = 2 * pt2[0] - pt1[0],
                                       _["y"] = pt1[1]);
  }
  if(pt2[1] == pt3[1]) {
    return Rcpp::NumericVector::create(_["x"] = pt1[0],
                                       _["y"] = 2 * pt2[1] - pt1[1]);
  }
  double a = pt_slope(pt2, pt3);
  double c = pt_projection_y(pt2, pt3, 0);
  double d = 2 * (a * pt1[0] - pt1[1] + c) / (a * a + 1);
  return Rcpp::NumericVector::create(_["x"] = pt1[0] - d * a,
                                     _["y"] = pt1[1] + d);
}

// helper to determine orthogonal projection of point pt1 on the line defined by pt2 and pt3.
// It returns pt' x and y coordinates of pt1 projected on the line pt2 - pt3.
// As a result pt'- pt1 and pt2 - pt3 segments will be perpendicular.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector pt_projection(const Rcpp::NumericVector pt1,  // pt of interest
                                  const Rcpp::NumericVector pt2,  // 1st coord of the segment
                                  const Rcpp::NumericVector pt3) {// 2nd coord of the segment
  // pt' is at the middle of pt1 and its pt_projection
  if(pt2[0] == pt3[0]) {
    return Rcpp::NumericVector::create(_["x"] = pt2[0],
                                       _["y"] = pt1[1]);
  }
  if(pt2[1] == pt3[1]) {
    return Rcpp::NumericVector::create(_["x"] = pt1[0],
                                       _["y"] = pt2[1]);
  }
  double a = pt_slope(pt2, pt3);
  double c = pt_projection_y(pt2, pt3, 0);
  double d = (a * pt1[0] - pt1[1] + c) / (a * a + 1);
  return Rcpp::NumericVector::create(_["x"] = pt1[0] - d * a,
                                     _["y"] = pt1[1] + d);
}

// [[Rcpp::export(rng = false)]]
double pt_distance_fast(const Rcpp::NumericVector pt1,
                        const Rcpp::NumericVector pt2,
                        const double scale = 1.0) {
  return std::pow(pt2[0] - pt1[0], 2.0) + std::pow(pt2[1] - pt1[1], 2.0);
}

// [[Rcpp::export(rng = false)]]
double pt_distance_accu(const Rcpp::NumericVector pt1,
                        const Rcpp::NumericVector pt2,
                        const double scale = 1.0) {
  return scale * std::sqrt(std::pow(pt2[0] - pt1[0], 2.0) + std::pow(pt2[1] - pt1[1], 2.0));
}

// Helper to determine shortest distance of a point pt1 to a segment formed by points pt2 and pt3.
// It returns R_PosInf if it projects outside of the segment
// [[Rcpp::export(rng = false)]]
double pt_shortest(const Rcpp::NumericVector pt1,  // pt of interest
                   const Rcpp::NumericVector pt2,  // 1st coord of the segment
                   const Rcpp::NumericVector pt3,  // 2nd coord of the segment
                   const double scale = 1.0,
                   const double eps = 0.000000001) {
  if((pt1[0] == pt2[0] && pt1[1] == pt2[1]) || (pt1[0] == pt3[0] && pt1[1] == pt3[1])) return 0;
  double dx = pt3[0] - pt2[0], dy = pt3[1] - pt2[1];
  double dot = (pt1[0] - pt2[0]) * dx + (pt1[1] - pt2[1]) * dy;
  if((dot < 0) && (dot > (-1.0 * eps))) dot = 0;
  double dist = dx * dx + dy * dy;
  double n = -1.0;
  n = dot / dist;
  if(n < 0.0) return R_PosInf;
  if((n > 1.0) && (n < (1 + eps))) n = 1.0;
  if(n > 1.0) return R_PosInf;
  return scale * std::sqrt(pow(pt1[0] - (pt2[0] + n * dx), 2.0) + pow(pt1[1] - (pt2[1] + n * dy), 2.0));
}

// pt_shortest is faster than pt_shortest2
// /!\ Note that if pt_shortest returns R_PosInf when pt1 projects outside of the segment formed by pt2 - pt3
// pt_shortest2 returns the distance of pt1 to the line, no matter if projection is outside of this segment
// [[Rcpp::export(rng = false)]]
double pt_shortest2(const Rcpp::NumericVector pt1,  // pt of interest
                    const Rcpp::NumericVector pt2,  // 1st coord of the segment
                    const Rcpp::NumericVector pt3,  // 2nd coord of the segment
                    const double scale = 1.0) {
  return pt_distance_accu(pt1, pt_projection(pt1, pt2, pt3), scale);
}

// area returns the area of a parallelogram defined by pt1, pt2 and pt3
// [[Rcpp::export(rng = false)]]
double area(const Rcpp::NumericVector pt1,
            const Rcpp::NumericVector pt2,
            const Rcpp::NumericVector pt3) {
  return std::abs(pt2[0] * (pt3[1] - pt1[1]) + pt3[0] * (pt1[1] - pt2[1]) + pt1[0] * (pt2[1] - pt3[1]));
}

// // [[Rcpp::export(rng = false)]]
// double area2(const R_len_t idx1,
//              const R_len_t idx2,
//              const R_len_t idx3,
//              const Rcpp::NumericMatrix pts) {
//   Rcpp::NumericVector pt1 = pts(idx1, Rcpp::_), pt2 = pts(idx2, Rcpp::_), pt3 = pts(idx3, Rcpp::_);
//   return std::abs(pt1[0]*pt2[1] - pt1[1]*pt2[0]
//                 + pt1[1]*pt3[0] - pt1[0]*pt3[1]
//                 + pt2[0]*pt3[1] - pt3[0]*pt2[1]);
// }

// // [[Rcpp::export(rng = false)]]
// double area3(const R_len_t idx1,
//              const R_len_t idx2,
//              const R_len_t idx3,
//              const Rcpp::NumericMatrix pts) {
//   Rcpp::NumericVector pt1 = pts(idx1, Rcpp::_), pt2 = pts(idx2, Rcpp::_), pt3 = pts(idx3, Rcpp::_);
//   double dx1 = pt2[0] - pt1[0], dy1 = pt2[1] - pt1[1],
//          dx2 = pt3[0] - pt1[0], dy2 = pt3[1] - pt1[1];
//   return std::abs(dx1*dy2 - dy1*dx2);
// }

#endif
