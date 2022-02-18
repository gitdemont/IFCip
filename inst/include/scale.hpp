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

#ifndef IFCIP_SCALE_HPP
#define IFCIP_SCALE_HPP

#include <Rcpp.h>
using namespace Rcpp;

Rcpp::NumericVector hpp_n_scale(Rcpp::NumericMatrix img,
                                const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                                const double value = NA_REAL,
                                const double n_lev = 256,
                                const bool positive = true) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();
  double MAX_LEV_SCA = positive ? (n_lev - 1) : (1 - n_lev);
  double mat_min = R_PosInf, mat_max = R_NegInf;
  if(msk_.isNotNull()) {
    Rcpp::NumericMatrix msk(msk_.get());
    if(mat_r != msk.nrow() || mat_c != msk.ncol()) {
      Rcpp::stop("hpp_scale: when 'msk' is provided 'img' and 'msk' should have same dimensions");
    }
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        if((msk(i_row, i_col) != NA_REAL) &&
           (msk(i_row, i_col) != R_NaN) &&
           (msk(i_row, i_col) != R_PosInf) &&
           (msk(i_row, i_col) != R_NegInf)) {
          if(msk(i_row, i_col)) {
            if((img(i_row, i_col) != NA_REAL) &&
               (img(i_row, i_col) != R_NaN) &&
               (img(i_row, i_col) != R_PosInf) &&
               (img(i_row, i_col) != R_NegInf)) {
              if(img(i_row, i_col) < mat_min) {
                mat_min = img(i_row, i_col);
              } else {
                if(img(i_row, i_col) > mat_max) mat_max = img(i_row, i_col);
              }
            } else {
              Rcpp::stop("hpp_scale: non-finite values found in 'img'");
            }
          }
        } else {
          Rcpp::stop("hpp_scale: non-finite values found in 'msk'");
        }
      }
    }
    MAX_LEV_SCA /= (mat_max - mat_min);
    if(positive) {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? MAX_LEV_SCA * (mat_max -img[i]) : value; 
    } else {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? MAX_LEV_SCA * (img[i] - mat_min) : value;
    }
    return Rcpp::NumericVector::create(mat_min, mat_max, MAX_LEV_SCA, positive, value);
  } else {
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        if((img(i_row, i_col) != NA_REAL) &&
           (img(i_row, i_col) != R_NaN) &&
           (img(i_row, i_col) != R_PosInf) &&
           (img(i_row, i_col) != R_NegInf)) {
          if(img(i_row, i_col) < mat_min) {
            mat_min = img(i_row, i_col);
          } else {
            if(img(i_row, i_col) > mat_max) mat_max = img(i_row, i_col);
          }
        } else {
          Rcpp::stop("hpp_scale: non-finite values found in 'img'");
        }
      }
    }
    MAX_LEV_SCA /= (mat_max - mat_min);
    if(positive) {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = MAX_LEV_SCA * (mat_max - img[i]);
    } else {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = MAX_LEV_SCA * (img[i] - mat_min);
    }
    return Rcpp::NumericVector::create(mat_min, mat_max, MAX_LEV_SCA, positive, value);
  }
}

Rcpp::NumericVector hpp_i_scale(Rcpp::IntegerMatrix img,
                                const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                                const int value = NA_INTEGER,
                                const int n_lev = 256,
                                const bool positive = true) {
  R_len_t mat_r = img.nrow(), mat_c = img.ncol();
  double MAX_LEV_SCA = positive ? (n_lev - 1) : (1 - n_lev);
  double mat_min = R_PosInf, mat_max = R_NegInf;
  if(msk_.isNotNull()) {
    Rcpp::NumericMatrix msk(msk_.get());
    if(mat_r != msk.nrow() || mat_c != msk.ncol()) {
      Rcpp::stop("hpp_scale: when 'msk' is provided 'img' and 'msk' should have same dimensions");
    }
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        if(msk(i_row, i_col) != NA_INTEGER) {
          if(msk(i_row, i_col)) {
            if(img(i_row, i_col) != NA_INTEGER) {
              if(img(i_row, i_col) < mat_min) {
                mat_min = img(i_row, i_col);
              } else {
                if(img(i_row, i_col) > mat_max) mat_max = img(i_row, i_col);
              }
            } else {
              Rcpp::stop("hpp_scale: non-finite values found in 'img'");
            }
          }
        } else {
          Rcpp::stop("hpp_scale: non-finite values found in 'msk'");
        }
      }
    }
    MAX_LEV_SCA /= (mat_max - mat_min);
    if(positive) {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? MAX_LEV_SCA * (mat_max -img[i]) : value; 
    } else {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? MAX_LEV_SCA * (img[i] - mat_min) : value;
    }
    return Rcpp::NumericVector::create(mat_min, mat_max, MAX_LEV_SCA, positive, value);
  } else {
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
        if(img(i_row, i_col) != NA_INTEGER) {
          if(img(i_row, i_col) < mat_min) {
            mat_min = img(i_row, i_col);
          } else {
            if(img(i_row, i_col) > mat_max) mat_max = img(i_row, i_col);
          }
        } else {
          Rcpp::stop("hpp_scale: non-finite values found in 'img'");
        }
      }
    }
    MAX_LEV_SCA /= (mat_max - mat_min);
    if(positive) {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = MAX_LEV_SCA * (mat_max - img[i]);
    } else {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = MAX_LEV_SCA * (img[i] - mat_min);
    }
    return Rcpp::NumericVector::create(mat_min, mat_max, MAX_LEV_SCA, positive, value);
  }
}

void hpp_n_scalerev(Rcpp::NumericMatrix img,
                    const Rcpp::NumericVector sca) {
  if(sca.size() != 5) Rcpp::stop("hpp_rescale: bad 'sca' argument");
  if(sca[3]) {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = sca[1] - img[i] / sca[2];
  } else {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = img[i] / sca[2] + sca[0];
  }
}

void hpp_i_scalerev(Rcpp::IntegerMatrix img,
                    const Rcpp::NumericVector sca) {
  if(sca.size() != 5) Rcpp::stop("hpp_rescale: bad 'sca' argument");
  if(sca[3]) {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = sca[1] - img[i] / sca[2];
  } else {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = img[i] / sca[2] + sca[0];
  }
}

//' @title Image Scaling
//' @name cpp_scale
//' @description
//' This function scales a image.
//' @param img, either non-null IntegerMatrix or NumericMatrix, with finite values. Non-finite values will trigger an error.
//' @param msk, a NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @param value, a double; it is the replacement value that will be used when 'msk' element is false. Default is NA_REAL..
//' @param n_lev, a double determining the number of grey levels used for the computation. Default is 256.
//' @param positive, a bool determining whether 'img' should be scaled with positive (when true) or negative(when false) values. Default is true.
//' @details when 'msk' is provided it has to be of the same dimensions as 'img', otherwise an error will be thrown.\cr
//' an error will be thrown also if 'msk' contains non-finite value.\cr
//' 'img' range will be determined based on indexes of non 0 'msk' values and only the values in 'img' at those indexes will be scaled; the other will be filled with 'value'.
//' @return a NumericVector containing min(img), max(img), scaling factor and 'positive' flag as masked according to 'msk' when provided.\cr
//' /!\ in addition 'img' will be modified (scaled) in-place.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_scale(SEXP img,
                              const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                              const double value = NA_REAL,
                              const double n_lev = 256,
                              const bool positive = true) {
  switch( TYPEOF(img) ) {
  case INTSXP: {
    return hpp_i_scale(img, msk_, value, n_lev, positive);
  }
    break;
  case REALSXP: {
    return hpp_n_scale(img, msk_, value, n_lev, positive);
  } break;
  default: Rcpp::stop("hpp_scale: not supported SEXP type in 'img'");
  }
}

//' @title Image Scaling Reversion
//' @name cpp_scalerev
//' @description
//' This function reverts image scaling done by hpp_scale.
//' @param img, an IntegerMatrix or NumericMatrix.
//' @param sca, a NumericVector as returned by hpp_scale.
//' @details all values in 'img' that are different from sca[4] will be reverted back to original values.
//' @return Nothing, but 'img' will be modified (scale reversion) in-place.
//' @keywords internal
////' @export
// [[Rcpp::export]]
void hpp_scalerev(SEXP img,
                  const Rcpp::NumericVector sca) {
  switch( TYPEOF(img) ) {
  case INTSXP: {
    hpp_i_scalerev(img, sca);
  }
    break;
  case REALSXP: {
    hpp_n_scalerev(img, sca);
  } break;
  default: Rcpp::stop("hpp_rescale: not supported SEXP type in 'img'");
  }
}

#endif
