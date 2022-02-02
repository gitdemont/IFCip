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
using namespace Rcpp;

static int ifcip_ws_dx [8]={ 1, 0,-1, 0,-1,-1, 1, 1};
static int ifcip_ws_dy [8]={ 0, 1, 0,-1, 1,-1,-1, 1};

// not a real pop since pop/push are quite slow but rather a pre-allocated vector
// used circularly where:
// x[0]             is current element(s) count
// x[x.size() - 1]  is current starting position
int fifo_pop (Rcpp::IntegerVector x) {
  int xx = x.size() - 1;
  int n = x[0];            // count
  int s = x[xx];           // start index
  int out = x[s];
  x[s] = NA_INTEGER;       // should we do it ?
  if(n > 0) {
    x[0]--;
    x[xx] = (s < (xx - 1)) ? (s + 1) : 1;
  } else { // should never happen
    Rcpp::stop("fifo_pop: can't pop vector");
  }
  return out;
}

// not a real push since pop/push are quite slow but rather a pre-allocated vector
// used circularly where:
// x[0]             is current element(s) count
// x[x.size() - 1]  is current starting position
void fifo_add (Rcpp::IntegerVector x,
               const int value) {
  int xx = x.size() - 1;
  int n = x[0];            // count
  int s = x[xx];           // start index
  if(n < (xx - 1)) {
    x[0]++;
    x[((s + n) < xx) ? (s + n) : (s + n - xx + 1)] = value;
  } else {  // should never happen
    Rcpp::stop("fifo_add: can't add value");
  }
}

// neighbor position computation
// Each time a neighbor is found, it is added to nbr vector starting at nbr[1]
// while nbr[0] store the total count of neighbors found.
// Only possible positions are stored i.e. inside nrow / ncol ranges.
// Depending on nbr.size(), 4 or 8 neighbors are retrieved:
// -4 neighbors = east, south, west, north
// -8 neighbors = 4 neigbors + diagonals
void get_neighbor(const R_len_t cur_idx,
                  const R_len_t nrow, 
                  const R_len_t ncol,
                  Rcpp::IntegerVector nbr) {
  R_len_t n = 0;
  if((cur_idx >= 0) && (cur_idx < nrow * ncol)) {
    for(R_len_t i = 0, i_col = cur_idx / nrow; i < (nbr.size() - 1); i++) {
      R_len_t x = ifcip_ws_dx[i] + i_col;
      R_len_t y = ifcip_ws_dy[i] + cur_idx - i_col * nrow;
      if((x >= 0) && (x < ncol) &&
         (y >= 0) && (y < nrow)) {
        nbr[++n] = y + x * nrow;
      }
    }
  } else { // should never happen
    Rcpp::stop("get_neighbor: cur_idx is out of matrix dimensions");
  }
  nbr[0] = n;
}

//' @title Watershed Transformation
//' @name cpp_watershed
//' @description
//' This function computes the watershed transformation of an image.
//' @param mat, a NumericMatrix; a distance transform matrix is expected.
//' @param connectivity, an uint8_t either 4 or 8 describing pixel neighborhood. Default is 8.
//' @param levels, an unsigned short determining the number elevation levels. Default is 256, should be at least 2.
//' @details adaptation of 'Watersheds in digital spaces: an efficient algorithm based on immersion simulations' from  L. Vincent and P. Soille.
//' In IEEE Transactions on Pattern Analysis and Machine Intelligence, 13(6):583-598, June 1991.\cr
//' The algorithm is reviewed in 'The Watershed Transform: Definitions, Algorithms and Parallelization Strategies'
//' from Roerdink, J. B. T. M., & Meijster, A. (2000) in Fundamenta Informaticae, 41, 187-228 \url{https://doi.org/10.3233/FI-2000-411207}
//' @return an IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_watershed(const Rcpp::NumericMatrix mat,
                                  const uint8_t connectivity = 8,
                                  const unsigned short levels = 256) {
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  R_len_t MAX_SIZ = mat_r * mat_c;
  if(MAX_SIZ <= 0) Rcpp::stop("'mat' should be at least 1px row and 1px col");
  if(MAX_SIZ >= (std::pow(2.0,31.0) - 2)) Rcpp::stop("'mat' is too large");
  if(!((connectivity==4)||(connectivity==8))) Rcpp::stop("'connectivity' should be either 4 or 8");
  if(levels < 2) Rcpp::stop("'levels' should be at least >= 2");
  int MAX_LEV = levels;
  
  // determines image range
  double mat_min = R_PosInf, mat_max = R_NegInf;
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    Rcpp::NumericVector col_ran = range(mat(Rcpp::_, i_col));
    if(col_ran[0] < mat_min) {
      mat_min = col_ran[0];
    } else {
      if(col_ran[1] > mat_max) mat_max = col_ran[1];
    }
  }
  
  // creates scaled mat ranging from [1 - MAX_LEV, 0]
  Rcpp::IntegerMatrix sca = Rcpp::no_init(mat_r, mat_c);
  double MAX_LEV_SCA = (1 - MAX_LEV) / mat_max;
  for(R_len_t i = 0; i < MAX_SIZ; i++) sca[i] = MAX_LEV_SCA * (mat[i] - mat_min);
  
  // sort values
  Rcpp::IntegerVector idx = seq_along(mat) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return sca[i] < sca[j];});
  
  // initialization
  R_len_t mask = -2, // px is initially found at threshold level
    wshed = 0,       // px belongs to a watershed line
    init = -1;       // initial value of out
  
  Rcpp::IntegerMatrix out = Rcpp::no_init(mat_r, mat_c);       // returned segmented image
  out.fill(init);
  R_len_t cur_lab = 0;
  Rcpp::IntegerMatrix dst(mat_r, mat_c);                       // matrix to store distances
  
  // working var
  Rcpp::IntegerVector nbr(connectivity + 1);                   // vector of neighbors idx
  Rcpp::IntegerVector Q(MAX_SIZ + 2, NA_INTEGER);              // hierarchical priority queue
  Q[0] = 0;
  Q[MAX_SIZ + 1] = 1;
  int k_start = 0, k = 0, k_stop = MAX_SIZ;
  
  // start flooding
  for(int h = 2-MAX_LEV; h <= 0; h++) {
    k_stop = MAX_SIZ;
    for(k = k_start; k < k_stop; k++) {
      if(sca[idx[k]] >= h) {
        k_stop = k;
        break;
      }
    }
    for(k = k_start; k < k_stop; k++) {
      int p = idx[k];
      out[p] = mask;
      get_neighbor(p, mat_r, mat_c, nbr);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        if((out[nbr[i]] > 0) || (out[nbr[i]] == wshed)) {
          dst[p]  = 1;
          fifo_add(Q, p);
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
      get_neighbor(p, mat_r, mat_c, nbr);
      for(R_len_t i = 1; i <= nbr[0]; i++) {
        int pp = nbr[i];
        if((dst[pp] < cur_dist) && ((out[pp] > 0) || (out[pp] == wshed))) {
          if(out[pp] > 0) {
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
          get_neighbor(fifo_pop(Q), mat_r, mat_c, nbr);
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
  // TODO-FIXME:
  // should we replace 0 (wshed values) by something else to keep track of watershed lines
  // before replacing remaining non visited -1 values ?
  for(k = k_start; k < mat.size(); k++) {
    out[idx[k]] = 0;
  }
  return out;
}

#endif
