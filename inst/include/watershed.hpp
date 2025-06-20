/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2022 Yohann Demont                                              

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
  
#ifndef IFCIP_WATERSHED_HPP
#define IFCIP_WATERSHED_HPP

#include <Rcpp.h>
#include "fifo.hpp"
#include "kernel.hpp"
#include "scale.hpp"
using namespace Rcpp;

//' @title Watershed Transformation SV1
//' @name cpp_watershed_sv1
//' @description
//' This function computes the watershed transformation of an image.
//' @param mat, a NumericMatrix; a distance transform matrix is expected.
//' @param n_lev, an unsigned short determining the number of elevation levels. Default is 256, should be at least 2.
//' @param draw_lines, a bool; whether to draw watershed lines or not. Default is true.
//' @param invert, a bool; whether to fill from basins (lowest values) to peaks (highest values). Default is false.
//' When 'mat' is the result of the distance transformation of an image, peaks (highest values) represent largest distances from background.
//' Thus, they are the ones to be filled first; this can be done with 'invert' set to true.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @details adaptation of 'Determining watersheds in digital pictures via flooding simulations' from P. Soille. and L. Vincent.
//' In Proc. SPIE 1360, Visual Communications and Image Processing '90: Fifth in a Series, (1 September 1990) \doi{10.1117/12.24211}.
//' @source MorphoLib plugin for ImageJ presents a Java implementation of the algorithm in  \url{https://github.com/ijpb/MorphoLibJ/blob/master/src/main/java/inra/ijpb/watershed/WatershedTransform2D.java} authored by Ignacio Arganda-Carreras 
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector hpp_watershed_sv1(const Rcpp::NumericMatrix mat,
                                      const unsigned short n_lev = 256,
                                      const bool draw_lines = true,
                                      const bool invert = false,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                                      const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  R_len_t MAX_SIZ = mat_r * mat_c;
  if(MAX_SIZ <= 0) Rcpp::stop("'mat' should be at least 1px row and 1px col");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 3)) Rcpp::stop("'mat' is too large");
  if(n_lev < 2) Rcpp::stop("'n_lev' should be at least >= 2");
  int MAX_LEV = n_lev;
  
  // create scaled image
  Rcpp::NumericMatrix img = Rcpp::clone(mat);
  int hop = -1;
  hpp_scale(img, msk_, hop, n_lev, invert, true);
  Rcpp::IntegerMatrix sca = as<Rcpp::IntegerMatrix>(img);
  
  // idx for sorting
  Rcpp::IntegerVector idx = seq_along(mat) - 1;
  
  // initialization
  int wshed = 0,
    init = -1,
    mask = -2,
    inqueue = -3;
  Rcpp::IntegerMatrix out = Rcpp::no_init(mat_r, mat_c);       // returned segmented image
  out.fill(init);
  R_len_t cur_lab = 0;
  
  // working var
  Rcpp::IntegerMatrix o = offset_kernel(kernel);               // matrix of offsets in raster order
  Rcpp::IntegerVector nbr(o.ncol() + 1);                       // vector of neighbors idx
  Rcpp::IntegerVector Q = fifo_create(MAX_SIZ + 3, NA_INTEGER);// hierarchical priority queue
  unsigned short count = 1;
  
  bool flag = true;
  int h, k_start = 0, k = 0, k_stop = mat.size();
  h = MAX_LEV - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return sca[i] > sca[j];});
  // start flooding
  while(h != 0) {
    h += hop;
    k_stop = mat.size();
    for(k = k_start; k < k_stop; k++) {
      if(sca[idx[k]] < h) {
        k_stop = k;
        break;
      }
    }
    for(k = k_start; k < k_stop; k++) {
      int p = idx[k];
      out[p] = mask;
      offset_nbr(p, mat_r, mat_c, o, nbr, &count);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        if((out[nbr[i]] > 0) || (out[p] == wshed)) {
          out[p] = inqueue;
          fifo_add(Q, p);
          break;
        }
      }
    }
    while(Q[0] != 0) {
      int p = fifo_pop(Q);
      offset_nbr(p, mat_r, mat_c, o, nbr, &count);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        if(out[nbr[i]] > 0) {
          if((out[p] == inqueue) || (flag && (out[p] == wshed))) {
            out[p] = out[nbr[i]];
          } else {
            if((out[p] > 0) && (out[p] != out[nbr[i]])) {
              out[p] = wshed;
              flag = false;
            }
          }
        } else {
          if(out[nbr[i]] == wshed) {
            if(out[p] == inqueue) {
              out[p] = wshed;
              flag = true;
            }
          } else {
            if(out[nbr[i]] == mask) {
              out[nbr[i]] = inqueue;
              fifo_add(Q, nbr[i]);
            }
          }
        }
      }
    }
    // search for new minima
    k = k_start;
    for(k = k_start; k < k_stop; k++) {
      int p = idx[k];
      if(out[p] == mask) {
        cur_lab++;
        out[p] = cur_lab;
        fifo_add(Q, p);
        while(Q[0] != 0) {
          offset_nbr(fifo_pop(Q), mat_r, mat_c, o, nbr, &count);
          for(R_len_t i = 1; i <= nbr[0]; i++) {
            if(out[nbr[i]] == mask) {
              fifo_add(Q, nbr[i]);
              out[nbr[i]] = cur_lab;
            }
          }
        }
      }
    }
    k_start = k_stop;
  }
  
  // remove watershed lines
  if(!draw_lines) {
    for(k = 0; k < MAX_SIZ; k++) {
      if(out[k] == wshed) {
        offset_nbr(k, mat_r, mat_c, o, nbr, &count);
        for(R_len_t i = 1; i <= nbr[0]; i++) {
          if(out[nbr[i]] > wshed) {
            fifo_add(Q, k);  
            break;
          }
        }
      }
    }
    while(Q[0] != 0) {
      int p = fifo_pop(Q);
      double d, dmax = R_NegInf;
      offset_nbr(p, mat_r, mat_c, o, nbr, &count);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        if(out[nbr[i]] > wshed) {
          d = std::abs(sca[nbr[i]] - sca[p]);
          if(d > dmax) {
            out[p] = out[nbr[i]];
            dmax = d;
          }
        }
      }
      if(out[p] == wshed) fifo_add(Q, p);
    }
  }
  for(k = k_start; k < MAX_SIZ; k++) {
    out[idx[k]] = 0;
  }
  return out;
}

//' @title Watershed Transformation SV2
//' @name cpp_watershed_sv2
//' @description
//' This function computes the watershed transformation of an image.
//' @param mat, a NumericMatrix; a distance transform matrix is expected.
//' @param n_lev, an unsigned short determining the number of elevation levels. Default is 256, should be at least 2.
//' @param draw_lines, a bool; whether to draw watershed lines or not. Default is true.
//' @param invert, a bool; whether to fill from basins (lowest values) to peaks (highest values). Default is false.
//' When 'mat' is the result of the distance transformation of an image, peaks (highest values) represent largest distances from background.
//' Thus, they are the ones to be filled first; this can be done with 'invert' set to true.
//' @param kernel, a NumericMatrix; the structuring shape determining neighborhood. All non-zero elements will be considered as neighbors (except center).\cr
//' Default is R_NilValue, resulting in 8-connected pixels neighbors computation.
//' @param msk_, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'mat' elements without masking anything.
//' @details adaptation of 'Watersheds in digital spaces: an efficient algorithm based on immersion simulations' from  L. Vincent and P. Soille.
//' In IEEE Transactions on Pattern Analysis and Machine Intelligence, 13(6):583-598, June 1991.\cr
//' @source The algorithm is reviewed in 'The Watershed Transform: Definitions, Algorithms and Parallelization Strategies'
//' from Roerdink, J. B. T. M. and Meijster, A. (2000) in Fundamenta Informaticae, 41, 187-228 \doi{10.3233/FI-2000-411207}
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_watershed_sv2(const Rcpp::NumericMatrix mat,
                                      const unsigned short n_lev = 256,
                                      const bool draw_lines = true,
                                      const bool invert = false,
                                      const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                                      const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  R_len_t MAX_SIZ = mat_r * mat_c;
  if(MAX_SIZ <= 0) Rcpp::stop("'mat' should be at least 1px row and 1px col");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 3)) Rcpp::stop("'mat' is too large");
  if(n_lev < 2) Rcpp::stop("'n_lev' should be at least >= 2");
  int MAX_LEV = n_lev;
  Rcpp::NumericMatrix kk = get_kernel(kernel);
  
  // create scaled image
  Rcpp::NumericMatrix img = Rcpp::clone(mat);
  int hop = -1;
  hpp_scale(img, msk_, hop, n_lev, invert, true);
  Rcpp::IntegerMatrix sca = as<Rcpp::IntegerMatrix>(img);
  
  // idx for sorting
  Rcpp::IntegerVector idx = seq_along(mat) - 1;
  
  // initialization
  int mask = -2,     // px is initially found at threshold level
    wshed = 0,       // px belongs to a watershed line
    init = -1;       // initial value of out
  
  Rcpp::IntegerMatrix out = Rcpp::no_init(mat_r, mat_c);       // returned segmented image
  out.fill(init);
  R_len_t cur_lab = 0;
  Rcpp::IntegerMatrix dst(mat_r, mat_c);                       // matrix to store distances
  
  // working var
  Rcpp::IntegerMatrix o = offset_kernel(kernel);               // matrix of offsets in raster order
  Rcpp::IntegerVector nbr(o.ncol() + 1);                       // vector of neighbors idx
  Rcpp::IntegerVector Q = fifo_create(MAX_SIZ + 3, NA_INTEGER);// hierarchical priority queue
  unsigned short count = 1;
  
  int h, k_start = 0, k = 0, k_stop = mat.size();
  h = MAX_LEV - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return sca[i] > sca[j];});
  
  // start flooding
  while(h != 0) {
    h += hop;
    k_stop = mat.size();
    for(k = k_start; k < k_stop; k++) {
      if(sca[idx[k]] < h) {
        k_stop = k;
        break;
      }
    }
    for(k = k_start; k < k_stop; k++) {
      int p = idx[k];
      out[p] = mask;
      offset_nbr(p, mat_r, mat_c, o, nbr, &count);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        if((out[nbr[i]] > 0) || (out[nbr[i]] == wshed)) {
          dst[p]  = 1;
          fifo_add(Q, p);
          break;
        }
      }
    }
    int cur_dist = 1; 
    fifo_add(Q, NA_INTEGER);
    while(true) {
      int p = fifo_pop(Q);
      if(p == NA_INTEGER) {
        if(Q[0] == 0) {
          break;
        } else {
          fifo_add(Q, NA_INTEGER);
          cur_dist++;
          p = fifo_pop(Q);
        }
      }
      offset_nbr(p, mat_r, mat_c, o, nbr, &count);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        int pp = nbr[i];
        if((dst[pp] < cur_dist) && ((out[pp] > 0) || (out[pp] == wshed))) {
          if((out[pp] > 0)) {
            if((out[p] == mask) || (out[p] == wshed)) {
              out[p] = out[pp];
            } else {
              if(out[p] != out[pp]) out[p] = wshed;
            }
          } else {
            if(out[p] == mask) out[p] = wshed;
          }
        } else {
          if((out[pp] == mask) && (dst[pp] == 0)) {
            dst[pp] = cur_dist + 1;
            fifo_add(Q, pp);
          }
        }
      }
    }
    // search for new minima
    k = k_start;
    for(k = k_start; k < k_stop; k++) {
      int p = idx[k];
      dst[p] = 0;
      if(out[p] == mask) {
        cur_lab++;
        fifo_add(Q, p); out[p] = cur_lab;
        while(Q[0] != 0) {
          offset_nbr(fifo_pop(Q), mat_r, mat_c, o, nbr, &count);
          for(R_len_t i = 1; i <= nbr[0]; i++) {
            if(out[nbr[i]] == mask) {
              fifo_add(Q, nbr[i]);
              out[nbr[i]] = cur_lab;
            }
          }
        }
      }
    }
    k_start = k_stop;
  }
  
  // remove watershed lines
  if(!draw_lines) {
    for(k = 0; k < MAX_SIZ; k++) {
      if(out[k] == wshed) {
        offset_nbr(k, mat_r, mat_c, o, nbr, &count);
        for(R_len_t i = 1; i <= nbr[0]; i++) {
          if(out[nbr[i]] > wshed) {
            fifo_add(Q, k);  
            break;
          }
        }
      }
    }
    while(Q[0] != 0) {
      int p = fifo_pop(Q);
      double d, dmax = R_NegInf;
      offset_nbr(p, mat_r, mat_c, o, nbr, &count);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        if(out[nbr[i]] > wshed) {
          d = std::abs(sca[nbr[i]] - sca[p]);
          if(d > dmax) {
            out[p] = out[nbr[i]];
            dmax = d;
          }
        }
      }
      if(out[p] == wshed) fifo_add(Q, p);
    }
  }
  for(k = k_start; k < MAX_SIZ; k++) {
    out[idx[k]] = 0;
  }
  return out;
}

#endif
