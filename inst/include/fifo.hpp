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

#ifndef IFCIP_FIFO_HPP
#define IFCIP_FIFO_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Clear FIFO Queue
//' @name fifo_clear
//' @description
//' This function clears a fifo queue.
//' @param Q, an IntegerVector; the fifo queue.
//' @details Q[0] will be set to 0 and Q[Q.size() - 1] to 1.
//' @return it returns nothing but Q will be modified in-place.
void fifo_clear (Rcpp::IntegerVector Q) {
  Q[0] = 0;
  Q[Q.size() - 1] = 1;
}

//' @title Create FIFO Queue
//' @name fifo_create
//' @description
//' This function creates a fifo queue.
//' @param size, a R_len_t; the desired size. Default is 1. Should be at least 1.
//' @param value, an int; the value used to fill the queue. Default is NA_INTEGER.
//' @details Q[0] will be set to 0 and Q[Q.size() - 1] to 1.
//' @return an IntegerVector.
Rcpp::IntegerVector fifo_create (const R_len_t size = 1,
                                 const int value = NA_INTEGER) {
  if(size < 1) Rcpp::stop("fifo queue: size should be at least 1, queue can't be created");
  Rcpp::IntegerVector out(size + 2, value);
  fifo_clear(out);
  return out;
}

//' @title Pop FIFO Queue
//' @name fifo_pop
//' @description
//' This function pops a fifo queue.
//' @param Q, an IntegerVector; the fifo queue.
//' @details not a real pop since pop/push are quite slow but rather a pre-allocated vector
//' used circularly where:
//' Q[0]             is current element(s) count,
//' Q[Q.size() - 1]  is current starting position.
//' Q[0] will be decreased by 1 and Q[Q.size() - 1] will be shift by -1
//' @return an int, the 1st element of the queue.
int fifo_pop (Rcpp::IntegerVector Q) {
  R_len_t xx = Q.size() - 1;
  R_len_t n = Q[0];            // count
  R_len_t s = Q[xx];           // start index
  int out = Q[s];
  Q[s] = NA_INTEGER;       // should we do it ?
  if(n > 0) {
    Q[0]--;
    Q[xx] = (s < (xx - 1)) ? (s + 1) : 1;
  } else { // should never happen
    Rcpp::stop("fifo_pop: can't pop vector");
  }
  return out;
}

//' @title Add to FIFO Queue
//' @name fifo_add
//' @description
//' This function push value to a fifo queue.
//' @param Q, an IntegerVector; the fifo queue.
//' @param value, an int; the value used to fill the queue. Default is NA_INTEGER. 
//' @details not a real pop since pop/push are quite slow but rather a pre-allocated vector
//' used circularly where:
//' Q[0]             is current element(s) count,
//' Q[Q.size() - 1]  is current starting position.
//' Q[0] will be increased by 1 and Q[Q.size() - 1] will be shift by +1
//' @return it returns nothing but Q will be modified in-place.
void fifo_add (Rcpp::IntegerVector Q,
               const int value) {
  R_len_t xx = Q.size() - 1;
  R_len_t n = Q[0];            // count
  R_len_t s = Q[xx];           // start index
  if(n < (xx - 1)) {
    Q[0]++;
    Q[((s + n) < xx) ? (s + n) : (s + n - xx + 1)] = value;
  } else {  // should never happen
    Rcpp::stop("fifo_add: can't add value");
  }
}

#endif
