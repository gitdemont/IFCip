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
// [[Rcpp::export(rng = false)]]
void fifo_clear (Rcpp::IntegerVector &Q) {
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
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector fifo_create (const R_len_t size = 1,
                                 const int value = NA_INTEGER) {
  if(size < 1) Rcpp::stop("fifo_create: size should be at least 1, fifo queue can't be created");
  Rcpp::IntegerVector out(size + 2, value);
  fifo_clear(out);
  return out;
}

//' @title Resize FIFO Queue
//' @name fifo_resize
//' @description
//' This function expand a fifo queue by doubling its capacity.
//' @param Q, an IntegerVector; the fifo queue.
//' @param value, an int; the value used to fill the queue. Default is NA_INTEGER.
//' @details Q capacity will be doubled and all values copied and current starting position, Q[Q.size()], reset to 1.
//' @return it returns nothing but Q will be modified in-place.
// [[Rcpp::export(rng = false)]]
void fifo_resize (Rcpp::IntegerVector &Q,
                  const int value = NA_INTEGER) {
  R_len_t SIZE = Q.size() > 3 ? 4 + 2 * (Q.size() - 3) : 4;
  if(SIZE > (std::pow(2.0,31.0) - 3)) Rcpp::stop("fifo_resize: can't resize queue");
  Rcpp::IntegerVector QQ = Rcpp::no_init_vector(SIZE);
  QQ[0] = Q[0];
  R_len_t j = 1;
  for(R_len_t i = Q[Q.size() - 1]; i < Q.size() - 1; i++) QQ[j++] = Q[i];
  for(R_len_t i = 1; i < Q[Q.size() - 1]; i++) QQ[j++] = Q[i];
  for(R_len_t i = Q.size(); i < QQ.size(); i++) QQ[i] = value;
  QQ[QQ.size() - 1] = 1;
  Q = QQ;
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
// [[Rcpp::export(rng = false)]]
int fifo_pop (Rcpp::IntegerVector &Q) {
  R_len_t xx = Q.size() - 1;
  R_len_t n = Q[0];            // count
  R_len_t s = Q[xx];           // start index
  int out = Q[s];
  Q[s] = NA_INTEGER; // should we do it ?
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
//' @param value, an int; the value used to fill the queue. 
//' @details not a real pop since pop/push are quite slow but rather a pre-allocated vector
//' used circularly where:
//' Q[0]             is current element(s) count,
//' Q[Q.size() - 1]  is current starting position.
//' Q[0] will be increased by 1 and Q[Q.size() - 1] will be shift by +1
//' @return it returns nothing but Q will be modified in-place.
// [[Rcpp::export(rng = false)]]
void fifo_add (Rcpp::IntegerVector &Q,
               const int value) {
  R_len_t xx = Q.size() - 1;
  R_len_t n = Q[0];            // count
  R_len_t s = Q[xx];           // start index
  if((n < (xx - 1)) && (n >= 0)) {
    Q[0]++;
    Q[((s + n) < xx) ? (s + n) : (s + n - xx + 1)] = value;
  } else { // dynamically expand queue
    if(n < 0) Rcpp::stop("fifo_add: can't add value");
    fifo_resize(Q);
    fifo_add(Q, value);
  }
}

//' @title Create Queue
//' @name queue_create
//' @description
//' This function creates a vector queue.
//' @param size, a R_len_t; the desired size. Default is 1. Should be at least 1.
//' @details Q[0] will be set to 0.
//' @return an IntegerVector.
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector queue_create (const R_len_t size = 1) {
  if(size < 1) Rcpp::stop("queue_create: size should be at least 1, vector queue can't be created");
  Rcpp::IntegerVector out = Rcpp::no_init_vector(size + 1);
  out[0] = 0;
  return out;
}

//' @title Resize Queue
//' @name queue_resize
//' @description
//' This function expand a queue by doubling its capacity.
//' @param Q, an IntegerVector; the queue.
//' @param value, an int; the value used to fill the queue. Default is NA_INTEGER.
//' @details Q capacity will be doubled and all values copied and current starting position, Q[Q.size()], reset to 1.
//' @return it returns nothing but Q will be modified in-place.
// [[Rcpp::export(rng = false)]]
void queue_resize (Rcpp::IntegerVector &Q) {
   R_len_t SIZE = Q.size() > 2 ? 3 + 2 * (Q.size() - 2) : 3;
   if(SIZE > (std::pow(2.0,31.0) - 2)) Rcpp::stop("queue_resize: can't resize queue");
   Rcpp::IntegerVector QQ = Rcpp::no_init_vector(SIZE);
   std::copy(Q.begin(), Q.end(), QQ.begin());
   Q = QQ;
 }

//' @title Pop Queue
//' @name queue_pop
//' @description
//' This function pops a vector queue.
//' @param Q, an IntegerVector; the vector queue.
//' @details not a real pop since pop/push are quite slow but rather a pre-allocated vector where:\cr
//' Q[0] is used to store current element(s) count.\cr
//' On pop, Q[0] will be decreased by 1 while last entered element will be returned.
//' @return an int, the most recently pushed element of the vector queue.
// [[Rcpp::export(rng = false)]]
int queue_pop (Rcpp::IntegerVector &Q) {
  R_len_t n = Q[0];
  if(n > 0) {
    Q[0]--;
  } else { // should never happen
    Rcpp::stop("queue_pop: can't pop vector");
  }
  return Q[n];
}

//' @title Push to Queue
//' @name queue_push
//' @description
//' This function pushes a value to a vector queue.
//' @param Q, an IntegerVector; the vector queue.
//' @param value, an int; the value used to fill the queue.
//' @details not a real pop since pop/push are quite slow but rather a pre-allocated vector where:\cr
//' Q[0] is used to store current element(s) count.\cr
//' On push, Q[0] will be increased by 1 while value will be positioned in Q at new Q[0].
//' @return it returns nothing but Q will be modified in-place.
// [[Rcpp::export(rng = false)]]
void queue_push (Rcpp::IntegerVector &Q, const int value) {
  if((Q[0] < (Q.size() - 1)) && (Q[0] >= 0)) {
    Q[0]++;
    Q[Q[0]] = value;
  } else { // dynamically expand queue
    if(Q[0] < 0) Rcpp::stop("queue_push: can't push value");
    queue_resize(Q);
    queue_push(Q, value);
  }
}

#endif
