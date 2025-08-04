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
#include "utils.hpp"
using namespace Rcpp;

template <int RTYPE>
Rcpp::NumericVector range_T(const Rcpp::Vector<RTYPE>& img,
                            const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  Rcpp::NumericVector out = Rcpp::NumericVector::create(R_PosInf, R_NegInf);
  if(msk_.isNotNull()) {
    Rcpp::NumericVector msk = getmask_T(img, msk_);
    if(RTYPE == REALSXP) {
      for(R_len_t i = 0; i < img.size(); i++) {
        if(R_finite(msk[i])) {
          if(msk[i]) {
            if(R_finite(img[i])) {
              if(img[i] < out[0]) out[0] = img[i];
              if(img[i] > out[1]) out[1] = img[i];
            } else {
              Rcpp::stop("hpp_range: non-finite values found in 'img'");
            }
          }
        } else {
          Rcpp::stop("hpp_range: non-finite values found in 'msk'");
        }
      }
    } else {
      for(R_len_t i = 0; i < img.size(); i++) {
        if(R_finite(msk[i])) {
          if(msk[i]) {
            if(img[i] != NA_INTEGER) {
              if(img[i] < out[0]) out[0] = img[i];
              if(img[i] > out[1]) out[1] = img[i];
            } else {
              Rcpp::stop("hpp_range: non-finite values found in 'img'");
            }
          }
        } else {
          Rcpp::stop("hpp_range: non-finite values found in 'msk'");
        }
      }
    }
  } else {
    if(RTYPE == REALSXP) {
      for(R_len_t i = 0; i < img.size(); i++) {
        if(R_finite(img[i])) {
          if(img[i] < out[0]) out[0] = img[i];
          if(img[i] > out[1]) out[1] = img[i];
        } else {
          Rcpp::stop("hpp_range: non-finite values found in 'img'");
        }
      }
    } else {
      for(R_len_t i = 0; i < img.size(); i++) {
        if(img[i] != NA_INTEGER) {
          if(img[i] < out[0]) out[0] = img[i];
          if(img[i] > out[1]) out[1] = img[i];
        } else {
          Rcpp::stop("hpp_range: non-finite values found in 'img'");
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector hpp_range(const SEXP img,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return range_T(as<Rcpp::LogicalVector>(img), msk_);
  case INTSXP : return range_T(as<Rcpp::IntegerVector>(img), msk_);
  case REALSXP : return range_T(as<Rcpp::NumericVector>(img), msk_);
  case RAWSXP : return range_T(as<Rcpp::RawVector>(img), msk_);
  default : Rcpp::stop("hpp_range: not supported SEXP in 'img'");
  }
}

template <int RTYPE>
void scalerev_T(Rcpp::Vector<RTYPE> img,
                const Rcpp::Nullable<Rcpp::NumericVector> sca_ = R_NilValue) {
  Rcpp::NumericVector sca;
  if(sca_.isNotNull()) {
    sca = sca_.get();
  } else {
    Rcpp::Nullable<Rcpp::NumericVector> foo = img.attr("scale");
    if(foo.isNotNull()) { 
      sca = foo.get();
    } else {
      Rcpp::stop("hpp_scalerev: bad 'sca' attribute");
    }
  }
  if(sca.size() != 5) Rcpp::stop("hpp_scalerev: bad 'sca' argument");
  if(sca[3]) {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = sca[1] - img[i] / sca[2];
  } else {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = img[i] / sca[2] + sca[0];
  }
  if(img.hasAttribute("scale")) img.attr("scale") = R_NilValue;
  if(img.hasAttribute("msk")) img.attr("msk") = R_NilValue;
  if(img.hasAttribute("levels")) img.attr("levels") = R_NilValue;
}

// [[Rcpp::export(rng = false)]]
void hpp_scalerev(SEXP img,
                  const Rcpp::Nullable<Rcpp::NumericVector> sca_ = R_NilValue) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return scalerev_T(as<Rcpp::LogicalVector>(img), sca_);
  case INTSXP : return scalerev_T(as<Rcpp::IntegerVector>(img), sca_);
  case REALSXP : return scalerev_T(as<Rcpp::NumericVector>(img), sca_);
  case RAWSXP : return scalerev_T(as<Rcpp::RawVector>(img), sca_);
  default : Rcpp::stop("hpp_scalerev: not supported SEXP in 'img'");
  }
}

template <int RTYPE>
Rcpp::NumericVector scale_T(Rcpp::Vector<RTYPE> img,
                            const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                            const double value = NA_REAL,
                            const int n_lev = 256,
                            const bool invert = false,
                            const bool bin = true) {
  if(std::abs(n_lev) < 2) Rcpp::stop("hpp_scale: abs(n_lev) should be at least 2.0");
  Rcpp::NumericVector out = Rcpp::NumericVector::create(_["min"]=R_PosInf, _["max"]=R_NegInf, _["factor"]=1.0, _["invert"]=invert, _["value"]=value);
  double MAX_LEV_SCA = bin + std::abs(n_lev) - 1;
  double MAX_LEV_BIN = std::abs(n_lev) - 1;
  if(n_lev < 0) {
    MAX_LEV_SCA *= -1.0;
    MAX_LEV_BIN *= -1.0;
  }
  Rcpp::NumericVector ran = range_T(img, msk_);
  out[0] = ran[0];
  out[1] = ran[1];
  out[2] = bin ? MAX_LEV_BIN : MAX_LEV_SCA;
  out[2] /= (out[1] - out[0]);
  double K = MAX_LEV_SCA / (out[1] - out[0]);
  double MX = out[1], MN = out[0];
  if(invert) {
    K *= -1.0;
    MX = out[0];
    MN = out[1];
  }
  if(msk_.isNotNull()) {
    Rcpp::NumericVector msk = getmask_T(img, msk_);
    if(bin) {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? (img[i] == MX ? MAX_LEV_BIN : std::min(MAX_LEV_BIN, std::trunc(K * (img[i] - MN)))) : value;
    } else {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? K * (img[i] - MN) : value;
    }
  } else {
    if(bin) {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = img[i] == MX ? MAX_LEV_BIN : std::min(MAX_LEV_BIN, std::trunc(K * (img[i] - MN)));
    } else {
      for(R_len_t i = 0; i < img.size(); i++) img[i] = K * (img[i] - MN);
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_scale(SEXP img,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                              const double value = NA_REAL,
                              const int n_lev = 256,
                              const bool invert = false,
                              const bool bin = true) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return scale_T(as<Rcpp::LogicalVector>(img), msk_, value, n_lev, invert, bin);
  case INTSXP : return scale_T(as<Rcpp::IntegerVector>(img), msk_, value, n_lev, invert, bin);
  case REALSXP : return scale_T(as<Rcpp::NumericVector>(img), msk_, value, n_lev, invert, bin);
  case RAWSXP : return scale_T(as<Rcpp::RawVector>(img), msk_, value, n_lev, invert, bin);
  default : Rcpp::stop("hpp_scale: not supported SEXP in 'img'");
  }
}

// Template to allow safe-use of hpp scale in R i.e. without in place modification of img
template <int RTYPE>
Rcpp::Vector<RTYPE> rescale_T(const Rcpp::Vector<RTYPE>& img,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                              const double value = NA_REAL,
                              const int n_lev = 256,
                              const bool invert = false,
                              const bool bin = false) {
  Rcpp::Vector<RTYPE> out = Rcpp::clone(img);
  out.attr("class") = "IFCip_rescale";
  out.attr("msk") = msk_;
  out.attr("scale") = hpp_scale(out, msk_, value, n_lev, invert, bin); 
  out.attr("levels") = n_lev;
  return out;
}

//' @title Image Scaling
//' @name cpp_rescale
//' @description
//' This function is designed to scale a SEXP to [0, n_lev - 1]
//' @param img, a SEXP (logical, raw, integer or numeric) vector or matrix containing image intensity values.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @param value, a double; it is the replacement value that will be used when 'msk' element is false. Default is NA_REAL.
//' @param n_lev, an int determining the number of levels used for the computation. Default is 256.
//' @param invert, a bool determining whether 'img' should be scaled from min to max (when false, [min(img),max(img)] becoming [0,n_lev-1]) or inverted (when true, with [max(img),min(img)] rescaled to [0,n_lev-1]) values. Default is false.
//' @param bin, a bool determining whether 'img' should be binned or if scaling should be continuous. Default is true to return discrete values.
//' @details when 'msk' is provided it has to be of the same dimensions as 'img', otherwise an error will be thrown.\cr
//' an error will be thrown also if 'msk' contains non-finite value.\cr
//' 'img' range will be determined based on indices of non 0 'msk' values and only the values in 'img' at those indices will be scaled; the others will be filled with 'value'.
//' @return a SEXP of same type as 'img' with class `IFCip_rescale`
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP hpp_rescale(const SEXP img,
                 const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                 const double value = NA_REAL,
                 const int n_lev = 256,
                 const bool invert = false,
                 const bool bin = false) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return rescale_T(as<Rcpp::RawVector>(img), msk_, value, n_lev, invert, bin);
  case INTSXP: return rescale_T(as<Rcpp::IntegerVector>(img), msk_, value, n_lev, invert, bin);
  case REALSXP: return rescale_T(as<Rcpp::NumericVector>(img), msk_, value, n_lev, invert, bin);
  case RAWSXP : return rescale_T(as<Rcpp::RawVector>(img), msk_, value, n_lev, invert, bin);
  default: Rcpp::stop("hpp_rescale: not supported SEXP type in 'img'");
  }
}

#endif