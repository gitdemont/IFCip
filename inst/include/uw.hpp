/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2023 Yohann Demont                                              
                                                                                
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
  FITNESS FOR arr PARTICULAR PURPOSE. In no event shall the copyright holders or
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
  
#ifndef IFCIP_UW_HPP
#define IFCIP_UW_HPP

#include <Rcpp.h>
#include "kernel.hpp"
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix kernel_coords(const Rcpp::NumericMatrix kernel,
                                  const bool erode = true) {
  double kcx = (kernel.ncol() - 1) / 2;
  double kcy = (kernel.nrow() - 1) / 2;
  
  R_len_t n = 0;
  for(R_len_t j = 0; j < kernel.size(); j++) if(kernel[j] && R_finite(kernel[j])) n++;
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(2, n);
  
  if(erode) {
    R_len_t kc_x = std::floor(kcx);
    R_len_t kc_y = std::floor(kcy);
    for(R_len_t c = 0, i = 0, j = 0; c < kernel.ncol(); c++) {
      for(R_len_t r = 0; r < kernel.nrow(); r++, j++) {
        if(kernel[j] && R_finite(kernel[j])) {
          out[i++] = c - kc_x;
          out[i++] = r - kc_y;
        }
      }
    } 
  } else {
    R_len_t kc_x = std::ceil(kcx);
    R_len_t kc_y = std::ceil(kcy);
    for(R_len_t c = 0, i = 0, j = kernel.size() - 1; c < kernel.ncol(); c++) {
      for(R_len_t r = 0; r < kernel.nrow(); r++, j--) {
        if(kernel[j] && R_finite(kernel[j])) {
          out[i++] = c - kc_x;
          out[i++] = r - kc_y;
        }
      }
    } 
  }
  return out;
}

//' @title Chordset Computation
//' @name chordset
//' @description
//' Computes chordset from kernel.
//' @param kern, a NumericMatrix.
//' @return a 4 columns IntegetMatrix whose columns are:\cr
//' - x starting position,\cr
//' - y starting position,\cr
//' - y ending position,\cr
//' - and length\cr
//' of the chords constituting the kernel in rows.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix chordset(const Rcpp::NumericMatrix kernel,
                             const bool erode = true) {
  // compute coordinates of non zero pixels from kernel center
  Rcpp::IntegerMatrix o = kernel_coords(kernel);
  if(o.ncol() == 0) { // no non zero pixel
    Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(0, 4);
    return out;
  }
  // compute chords
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(o.ncol() + 1, 4);
  R_len_t count = -1, k = 0;
  while(k < o.ncol()) {
    count++;
    if(count >= out.nrow()) Rcpp::stop("chordset: buffer overrun");
    R_len_t n = 0;
    out(count, 0) = o(0, k);     // x start
    out(count, 1) = o(1, k);     // y start
    while(k < o.ncol() && (out(count, 0) == o(0, k)) && (out(count, 1) == o(1, k) - n)) {
      n++;
      k++;
    }
    out(count, 2) = o(1, k - 1); // y stop
    out(count, 3) = n;           // length
  }
  return out(Range(0, count), Rcpp::_);
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector chordset_R(const IntegerMatrix C) {
  // sort chords length
  Rcpp::IntegerVector U = sort_unique(C(Rcpp::_, 3));
  
  // prepare returned vector
  Rcpp::IntegerVector out = Rcpp::no_init_vector(U[U.size() - 1]);
  
  // computes R so that for every chords of length A in chordset C
  // it exists a chord of length B that respects A/2 <= B < A
  R_len_t i = 0, count = 0, n = 1;
  while(i < U.size()) {
    while(U[i] > n) {
      if(count >= out.size()) Rcpp::stop("chordset_R: buff over while computing R array");
      out[count++] = n;
      n *= 2;
    }
    n = U[i++];
    if(count >= out.size()) Rcpp::stop("chordset_R: buff over while computing R array");
    out[count++] = n;
    n *= 2;
  }
  return out[Range(0, count - 1)];
}

template <int RTYPE, typename TYPE>
void LUT_erode_T (Rcpp::Vector<RTYPE> T,
                  const Rcpp::IntegerVector O,
                  const Rcpp::IntegerVector R,
                  const Rcpp::Matrix<RTYPE> mat,
                  const R_len_t coln,
                  const R_len_t xoff,
                  const R_len_t xmin,
                  const R_len_t ymax,
                  const R_len_t ymin,
                  const TYPE limit) {
  // Rcpp::checkUserInterrupt();
  R_len_t mat_r = mat.nrow(),
    d0 = mat_r + ymax - ymin + 1,
    d1 = R.size(),
    x = coln + xoff,
    rr = O[xoff - xmin];
  
  // fill for R[0]
  for(typename Rcpp::Vector<RTYPE>::iterator T_it = T.begin() + rr; T_it != T.begin() + rr + d0 * d1; ++T_it) *T_it = limit;
  if(x >= 0 && x < mat.ncol()) for(R_len_t y = 0, yy = rr - ymin; y < mat_r; y++, yy++) T[yy] = mat[y + x * mat_r];
  for(R_len_t i = 1; i < d1; i++) {
    R_len_t d = R[i] - R[i - 1];
    for(R_len_t yy = i * d0 + rr; yy < i * d0 + rr + d0 - d; yy++) {
      T[yy] = std::min(T[yy - d0], T[yy - d0 + d]);
    }
  }
}

template <int RTYPE, typename TYPE>
void LUT_dilate_T (Rcpp::Vector<RTYPE> T,
                   const Rcpp::IntegerVector O,
                   const Rcpp::IntegerVector R,
                   const Rcpp::Matrix<RTYPE> mat,
                   const R_len_t coln,
                   const R_len_t xoff,
                   const R_len_t xmin,
                   const R_len_t ymax,
                   const R_len_t ymin,
                   const TYPE limit) {
  // Rcpp::checkUserInterrupt();
  R_len_t mat_r = mat.nrow(),
    d0 = mat_r + ymax - ymin + 1,
    d1 = R.size(),
    x = coln + xoff,
    rr = O[xoff - xmin];
  
  // fill for R[0]
  for(typename Rcpp::Vector<RTYPE>::iterator T_it = T.begin() + rr; T_it != T.begin() + rr + d0 * d1; ++T_it) *T_it = limit;
  if(x >= 0 && x < mat.ncol()) for(R_len_t y = 0, yy = rr - ymin; y < mat_r; y++, yy++) T[yy] = mat[y + x * mat_r];
  
  // fill for R[1-R.size()[
  for(R_len_t i = 1; i < d1; i++) {
    R_len_t d = R[i] - R[i - 1];
    for(R_len_t yy = i * d0 + rr; yy < i * d0 + rr + d0 - d; yy++) {
      T[yy] = std::max(T[yy - d0], T[yy - d0 + d]);
    }
  }
}

template <int RTYPE>
void erode_column_T(Rcpp::Matrix<RTYPE> out,
                    Rcpp::Vector<RTYPE> T,
                    Rcpp::IntegerVector idx,
                    const R_len_t beg,
                    const R_len_t coln) {
  // Rcpp::checkUserInterrupt();
  for(Rcpp::IntegerVector::iterator idx_it = idx.begin(); idx_it != idx.end(); ++idx_it) {
    for(typename Rcpp::Vector<RTYPE>::iterator T_it = T.begin() + (*idx_it + beg) % T.size(), out_it = out.begin() + coln * out.nrow();
        out_it != out.begin() + (coln + 1) * out.nrow();
        ++T_it, ++out_it) {
      *out_it = std::min(*out_it, *T_it);
    }
  }
}

template <int RTYPE>
void dilate_column_T(Rcpp::Matrix<RTYPE> out,
                     Rcpp::Vector<RTYPE> T,
                     Rcpp::IntegerVector idx,
                     const R_len_t beg,
                     const R_len_t coln) {
  // Rcpp::checkUserInterrupt();
  for(Rcpp::IntegerVector::iterator idx_it = idx.begin(); idx_it != idx.end(); ++idx_it) {
    for(typename Rcpp::Vector<RTYPE>::iterator T_it = T.begin() + (*idx_it + beg) % T.size(), out_it = out.begin() + coln * out.nrow();
        out_it != out.begin() + (coln + 1) * out.nrow();
        ++T_it, ++out_it) {
      *out_it = std::max(*out_it, *T_it);
    }
  }
}

template <int RTYPE, int RTYPEk, typename TYPE>
Rcpp::Matrix<RTYPE> uw_T(Rcpp::Matrix<RTYPE> mat,
                         Rcpp::Matrix<RTYPEk> kernel,
                         const TYPE limit) {
  Rcpp::IntegerMatrix C = chordset(kernel, limit > 0);
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  std::copy(mat.begin(), mat.end(), out.begin());
  if(C.nrow() == 0) return out;
  // compute R
  Rcpp::IntegerVector R = chordset_R(C);
  
  // some vars
  R_len_t ymin = Rcpp::min(C(Rcpp::_, 1)),
    ymax = Rcpp::max(C(Rcpp::_, 2)),
    xmin = Rcpp::min(C(Rcpp::_, 0)),
    xmax = Rcpp::max(C(Rcpp::_, 0)),
    d0 = mat.nrow() + ymax - ymin + 1,
    d1 = R.size(),
    d2 = xmax - xmin + 1,
    dd = d0 * d1;
  
  // declare LUT
  Rcpp::Vector<RTYPE> T = Rcpp::no_init_vector(d0 * d1 * d2);
  Rcpp::IntegerVector O = dd * (seq_len(d2) - 1);
  
  // vector of idx offsets
  Rcpp::IntegerVector idx = Rcpp::no_init_vector(C.nrow());
  for(R_len_t i = 0; i < C.nrow(); i++) {
    R_len_t R_id = 0;
    while(C(i, 3) != R[R_id]) if(++R_id >= R.size()) Rcpp::stop("uw_T: no matching chord length found");
    idx[i] = C(i, 1) - ymin + R_id * d0 + (C(i, 0) - xmin) * dd;
  }
  
  // erode / dilate
  if(limit > 0) {
    // init LUT for column 0
    for(R_len_t xoff = xmin; xoff <= xmax; xoff++) LUT_erode_T(T, O, R, mat, 0, xoff, xmin, ymax, ymin, limit);
    // erode column 0
    erode_column_T(out, T, idx, 0, 0);
    for(R_len_t coln = 1; coln < out.ncol(); coln++) {
      // rotate offsets so that T(;;r) = T(;;r-1) for [xmin - xmax[
      std::rotate(O.begin(), O.begin() + 1, O.end());
      // compute LUT for xmax
      LUT_erode_T(T, O, R, mat, coln, xmax, xmin, ymax, ymin, limit);
      // erode column
      erode_column_T(out, T, idx, O[0], coln); 
    }
  } else {
    // init LUT for column 0
    for(R_len_t xoff = xmin; xoff <= xmax; xoff++) LUT_dilate_T(T, O, R, mat, 0, xoff, xmin, ymax, ymin, limit);
    // dilate column 0
    dilate_column_T(out, T, idx, 0, 0);
    for(R_len_t coln = 1; coln < out.ncol(); coln++) {
      // rotate offsets so that T(;;r) = T(;;r-1) for [xmin - xmax[
      std::rotate(O.begin(), O.begin() + 1, O.end());
      // compute LUT for xmax
      LUT_dilate_T(T, O, R, mat, coln, xmax, xmin, ymax, ymin, limit);
      // dilate column
      dilate_column_T(out, T, idx, O[0], coln); 
    }
  }
  return out;
}

//' @title Urbach-Wilkinson Algorithm for Image Erosion and Dilation
//' @name cpp_uw
//' @description
//' This function applies erosion or dilatation on image.
//' @param mat, a Matrix.
//' @param kernel, a Nullable Matrix.
//' @param erode, a bool. whether to do image erosion or dilatation. Default is true to perform erosion.
//' @details see 'Efficient 2-D grayscale morphological transformations with arbitrary flat structuring elements' from  E.R. Urbach, M.H.F. Wilkinson.
//' IEEE Transactions on Image Processing, 17(1):1-8, January 2008.\doi{10.1109/tip.2007.912582}
//' @return a Matrix of same type as 'mat'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP hpp_uw (SEXP mat,
             SEXP kernel,
             const bool erode = true) {
  switch( TYPEOF(mat) ) {
  // case NILSXP : return R_NilValue; // prefer throwing error on NILSXP
  case LGLSXP : {
    int limit = erode ? 1 : 0;
    return uw_T(as<Rcpp::LogicalMatrix>(mat), get_kernel(kernel), limit);
  }
  case RAWSXP : {
    int limit = erode ? 255 : 0;
    return uw_T(as<Rcpp::RawMatrix>(mat), get_kernel(kernel), limit);
  }
  case INTSXP : {
    int limit = erode ? 2147483647 : -2147483647;
    return uw_T(as<Rcpp::IntegerMatrix>(mat), get_kernel(kernel), limit);
  }
  case REALSXP : {
    double limit = erode ? R_PosInf : R_NegInf;
    return uw_T(as<Rcpp::NumericMatrix>(mat), get_kernel(kernel), limit);
  }
  default : Rcpp::stop("hpp_uw: not supported SEXP[%i] in 'mat'", TYPEOF(mat));
  }
}

#endif
