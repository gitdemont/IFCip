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

#ifndef IFCIP_CONVEX_HPP
#define IFCIP_CONVEX_HPP

#include <Rcpp.h>
using namespace Rcpp;

// helper to sort one vector against 2 others
Rcpp::IntegerVector sort_1vs2(const Rcpp::IntegerVector V,
                              const Rcpp::NumericVector vec1,
                              const Rcpp::NumericVector vec2) {
  Rcpp::IntegerVector idx = seq_along(V) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){
    if ( vec1[i] == vec1[j] ) {
      return vec2[i] < vec2[j];
    }
    return vec1[i] < vec1[j];
  });
  return V[idx];
}

// helper to compute cross product. it returns false if 3 points make a counter clockwise turn
bool cross (const Rcpp::NumericVector pt1, 
                   const Rcpp::NumericVector pt2, 
                   const Rcpp::NumericVector pt3) {
  return ((pt2[0] - pt1[0]) * (pt3[1] - pt1[1]) <= (pt2[1] - pt1[1]) * (pt3[0] - pt1[0]));
}

//' @title Convex Hull
//' @name cpp_convexhull
//' @description
//' Computes 2D convex hull of a set of points.
//' @param pts a 2-columns matrix defining the locations (x and y coordinates, respectively) of points.
//' x has to be in column 0 and y in column 1 (C index, add 1 for R).
//' @return a 2-column matrix subset of points within 'pts' that constitutes convex hull vertices.
//' @source Adaptation of Andrew's Monotone Chain Algorithm A. M. Andrew, 'Another Efficient Algorithm
//' for Convex Hulls in Two Dimension', Information Processing Letters, 9, 1979, pp216-219 reported
//' in M. A. Jayaram, Hasan Fleyeh, 'Convex Hulls in Image Processing: A Scoping Review', 
//' American Journal of Intelligent Systems, Vol. 6 No. 2, 2016, pp. 48-58.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_convexhull(const Rcpp::NumericMatrix pts) {
  Rcpp::NumericVector x = pts(Rcpp::_, 0);
  Rcpp::NumericVector y = pts(Rcpp::_, 1);
  Rcpp::IntegerVector idx = Rcpp::seq(0, pts.nrow() - 1);
  Rcpp::IntegerVector oo =  sort_1vs2(idx, x, y);
  
  // reorder input according to sorted x-y
  Rcpp::NumericMatrix M(pts.nrow(), pts.ncol());
  for(R_len_t i = 0; i < M.nrow(); i++) {
    M(i, Rcpp::_) = pts(oo[i], Rcpp::_);
  }
  
  // compute convex hull upper and lower
  Rcpp::IntegerVector L, U;
  for(R_len_t i = 0; i < M.nrow(); i++) {
    while(L.size() >= 2 && cross(M(L[L.size() - 2], Rcpp::_), M(L[L.size()- 1], Rcpp::_), M(i, Rcpp::_))) {
      L.erase(L.size() - 1);
    }
    L.push_back(i);
  }
  for(R_len_t i = M.nrow() - 1; i >= 0; i--) {
    while(U.size() >= 2 && cross(M(U[U.size() - 2], Rcpp::_), M(U[U.size() - 1], Rcpp::_), M(i, Rcpp::_))) {
      U.erase(U.size() - 1);
    }
    U.push_back(i);
  }
  
  // return result
  Rcpp::IntegerVector ret = Rcpp::no_init_vector(L.size() + U.size() - 2);
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(ret.size(), 2);
  R_len_t i = 0;
  for(; i < L.size() - 1; i++) {
    ret[i] = oo[L[i]] + 1;
    out(i, 0) = x[oo[L[i]]];
    out(i, 1) = y[oo[L[i]]];
  }
  for(R_len_t j = 0; j <  U.size() - 1; j++) {
    ret[i + j] = oo[U[j]] + 1;
    out(i + j, 0) = x[oo[U[j]]];
    out(i + j, 1) = y[oo[U[j]]];
  }
  colnames(out) = CharacterVector::create("x", "y");
  out.attr("subset") = ret;
  out.attr("class") = "IFCip_convexhull";
  return out;
}

#endif
