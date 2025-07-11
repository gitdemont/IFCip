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

#ifndef IFCIP_KERNEL_HPP
#define IFCIP_KERNEL_HPP

#include <Rcpp.h>
using namespace Rcpp;

// retrieve kernel matrix or set it to default [[1,1,1],[1,1,1],[1,1,1]] if NULL or of size 0
Rcpp::NumericMatrix get_kernel(const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue) {
  Rcpp::NumericMatrix kk;
  if(kernel.isNotNull()) {
    Rcpp::NumericMatrix k(kernel.get());
    kk = k;
  }
  if(kk.size() == 0) {
    kk = Rcpp::NumericMatrix(3,3);
    kk.fill(1.0);
  }
  return kk;
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix rev_mat(const Rcpp::IntegerMatrix mat,
                            bool reverse = true) {
  if(!reverse) return mat;
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  for(R_len_t i = 0; i < mat.size(); i++) out[i] = mat[mat.size() - 1 - i];
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix swap_mat(const Rcpp::IntegerMatrix mat,
                             const bool swap = true) {
  if(!swap) return mat;
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  for(R_len_t i = 0, j = 0; i < mat.ncol(); i++) out(Rcpp::_, j++) = mat(Rcpp::_, i);
  return out;
}

// compute matrix of x, y, offsets relative to kernel center
// nozero: according to non-zero kernel elements
// center: including center
// reverse: should matrix be reversed, this is useful to adjust center when kernel has even dimension depending on reading direction
// NULL kernel will result in 8 connected neighbors offsets
// return an IntegerMatrix whose rows are [x, y, index], respectively and columns are elements
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix offset_kernel(const Rcpp::Nullable<Rcpp::NumericMatrix> kernel = R_NilValue,
                                  bool nozero = true,
                                  bool center = false,
                                  bool reverse = false) {
  Rcpp::LogicalMatrix k = Rcpp::as<Rcpp::LogicalMatrix>(get_kernel(kernel));
  Rcpp::IntegerMatrix kk = rev_mat(Rcpp::as<Rcpp::IntegerMatrix>(k), reverse);
  R_len_t rr = kk.nrow() >> 1;
  R_len_t rc = kk.ncol() >> 1;
  R_len_t kc = kk.ncol() % 2;
  R_len_t kr = kk.nrow() % 2;
  R_len_t rc_1 = rc + kc;
  R_len_t rr_1 = rr + kr;
  if(nozero) {
    R_len_t n = 0;
    for(R_len_t c = -rc, j = 0; c < rc_1; c++) {
      for(R_len_t r = -rr; r < rr_1; r++, j++) {
        if(kk[j] == 1) {
          if((r == 0) && (c == 0) && !center) continue;
          n++;
        }
      }
    }
    Rcpp::IntegerMatrix out(3, n);
    for(R_len_t c = -rc, i = 0, j = 0; c < rc_1; c++) {
      for(R_len_t r = -rr; r < rr_1; r++, j++) {
        if(kk[j] == 1) {
          if((r == 0) && (c == 0) && !center) continue;
          out[i++] = c;
          out[i++] = r;
          out[i++] = j;
        }
      }
    }
    return swap_mat(out, reverse);
  }
  R_len_t s = kk.size() - !center;
  Rcpp::IntegerMatrix out(3, std::max(0, s));
  for(R_len_t c = -rc, i = 0, j = 0; c < rc_1; c++) {
    for(R_len_t r = -rr; r < rr_1; r++, j++) {
      if((r == 0) && (c == 0) && !center) continue;
      out[i++] = c;
      out[i++] = r;
      out[i++] = j;
    }
  }
  return swap_mat(out, reverse);
}

// determine kernel [x, y, index] that are scanned before (kernel) center in raster order
// return an IntegerMatrix whose rows are [x, y, index], respectively and columns are elements
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix offset_backward(const Rcpp::IntegerMatrix offset) {
  if(offset.nrow() < 2) Rcpp::stop("offset_backward: 'offset' should be at least 2 rows");
  Rcpp::LogicalVector L = Rcpp::no_init_vector(offset.ncol());
  for(R_len_t i = 0; i < offset.ncol(); i++) L[i] = !((offset(0, i) > 0) || ((offset(0, i) == 0) && (offset(1, i) > 0)));
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(offset.nrow(), sum(L));
  for(R_len_t i = 0, j = 0; i < offset.ncol(); i++) if(L[i]) out(Rcpp::_, j++) = offset(Rcpp::_, i);
  return out;
}

// determine kernel [x, y, index] that are scanned after (kernel) center in raster order
// return an IntegerMatrix whose rows are [x, y, index], respectively and columns are elements
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix offset_forward(const Rcpp::IntegerMatrix offset) {
  if(offset.nrow() < 2) Rcpp::stop("offset_forward: 'offset' should be at least 2 rows");
  Rcpp::LogicalVector L = Rcpp::no_init_vector(offset.ncol());
  for(R_len_t i = offset.ncol() - 1; i >= 0; i--) L[i] = ((offset(0, i) > 0) || ((offset(0, i) == 0) && (offset(1, i) >= 0)));
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(offset.nrow(), sum(L));
  for(R_len_t i = offset.ncol() - 1, j = 0; i >= 0; i--) if(L[i]) out(Rcpp::_, j++) = offset(Rcpp::_, i);
  return out;
}

// neighbor position computation according to offset
// Each time a neighbor is found, it is added to nbr vector starting at nbr[1]
// while nbr[0] stores the total count of neighbors found.
// Only possible positions are stored i.e. inside nrow / ncol ranges.
void offset_nbr(const R_len_t idx,
                const R_len_t nrow, 
                const R_len_t ncol,
                const Rcpp::IntegerMatrix offset,
                Rcpp::IntegerVector nbr,
                unsigned short *ptr) {
  if(((*ptr)++ % 10000) == 0) {
    *ptr = 1;
    Rcpp::checkUserInterrupt();
  }
  R_len_t n = 0;
  if((idx >= 0) && (idx < nrow * ncol)) {
    for(R_len_t i = 0, i_col = idx / nrow; i < offset.ncol(); i++) {
      R_len_t x = offset(0, i) + i_col;
      R_len_t y = offset(1, i) + idx - i_col * nrow;
      if((x >= 0) && (x < ncol) &&
         (y >= 0) && (y < nrow)) {
        nbr[++n] = y + x * nrow;
      }
    }
  } else { // should never happen
    Rcpp::stop("offset_nbr: 'idx' is out of matrix dimensions");
  }
  nbr[0] = n;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_make_sobel() {
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(3,3);
  out[0]=-1.0; out[3]=0.0; out[6]=1.0;
  out[1]=-2.0; out[4]=0.0; out[7]=2.0;
  out[2]=-1.0; out[5]=0.0; out[8]=1.0;
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_make_scharr() {
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(3,3);
  out[0]= -3.0; out[3]=0.0; out[6]= 3.0;
  out[1]=-10.0; out[4]=0.0; out[7]=10.0;
  out[2]= -3.0; out[5]=0.0; out[8]= 3.0;
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_make_gaussian(const uint8_t size,
                                      const double sigma = -0.3) {
  double s = sigma < 0 ? -1.0 * sigma * ((size - 1) * 0.5) : sigma;
  R_len_t h = size >> 1;
  double c = size % 2 ? (2.0 * h) : (2.0 * h - 1);
  c *= 0.5;
  double k = s * s;
  k *= 2.0;
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(size, size);
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double x = i_col - c;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      double y = i_row - c;  
      out(i_row, i_col) = std::exp(-1.0 * (x * x + y * y) / k);
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_make_laplacian(const uint8_t size,
                                       const double sigma = -0.3) {
  double s = sigma < 0 ? -1.0 * sigma * ((size - 1) * 0.5) : sigma;
  R_len_t h = size >> 1;
  double c = size % 2 ? (2.0 * h) : (2.0 * h - 1);
  c *= 0.5;
  double k = s * s;
  double kk = k * k;
  k *= 2.0;
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(size, size);
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double x = i_col - c;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      double y = i_row - c; 
      double foo = -1.0 * (x * x + y * y) / k;
      out(i_row, i_col) = -(1.0 + foo) * std::exp(foo) / kk;
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix hpp_make_disc(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out = Rcpp::no_init_matrix(size, size);
  if(size == 0) return out;
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double foo = i_col - half;
    foo = foo < 0 ? foo + 0.3 : foo - 0.3;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      double bar = i_row - half;
      bar = bar < 0 ? bar + 0.3 : bar - 0.3;
      out(i_row, i_col) = std::sqrt(foo * foo + bar * bar) <= half;
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix hpp_make_box(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out = Rcpp::no_init_matrix(size, size);
  out.fill(true);
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix hpp_make_plus(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  if(half > 0.0) {
    out(half, Rcpp::_) = Rcpp::rep(true, size);
    out(Rcpp::_, half) = Rcpp::rep(true, size);
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix hpp_make_cross(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out = Rcpp::no_init_matrix(size, size);
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      out(i_row, i_col) = i_row == i_col || i_row == (size - 1 - i_col);
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalMatrix hpp_make_diamond(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size >= 1) {
    double half = size >> 1;
    R_len_t i = 0;
    if(size % 2) {
      for(R_len_t i_col = -half; i_col <= half; i_col++) {
        for(R_len_t i_row = -half; i_row <= half; i_row++) {
          if(i >= size * size) Rcpp::stop("Not allowed");
          out[i++] = (abs(i_col) + abs(i_row)) <= half;
        }
      }
    } else {
      for(R_len_t i_col = -half; i_col <= half; i_col++) {
        if(i_col == 0) i_col++;
        for(R_len_t i_row = -half; i_row <= half; i_row++) {
          if(i_row == 0) i_row++;
          if(i >= size * size) Rcpp::stop("Not allowed");
          out[i++] = (abs(i_col) + abs(i_row)) <= half;
        }
      }
    }
  }
  return out;
}

#endif
