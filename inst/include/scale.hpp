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
bool is_scaled_T(const Rcpp::Vector<RTYPE>& img) {
  if(Rf_inherits(img, "IFCip_rescale")) {
    if(img.hasAttribute("scale")) {
      Rcpp::NumericVector sca = img.attr("scale"); 
      if(sca.size() == 7) return true;
      Rcpp::stop("hpp_is_scaled: bad 'scale' attribute size");
    } else {
      Rcpp::stop("hpp_is_scaled: missing 'scale' attribute");
    }
  }
  return false;
}

// [[Rcpp::export(rng = false)]]
bool hpp_is_scaled(SEXP img) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return is_scaled_T(as<Rcpp::LogicalVector>(img));
  case INTSXP : return is_scaled_T(as<Rcpp::IntegerVector>(img));
  case REALSXP : return is_scaled_T(as<Rcpp::NumericVector>(img));
  case RAWSXP : return is_scaled_T(as<Rcpp::RawVector>(img));
  default : Rcpp::stop("hpp_is_scaled: not supported SEXP in 'img'");
  }
}

template <int RTYPE>
void scalerev_T(Rcpp::Vector<RTYPE> img,
                const Rcpp::Nullable<Rcpp::NumericVector> sca_ = R_NilValue) {
  Rcpp::NumericVector sca;
  if(sca_.isNotNull()) {
    sca = sca_.get();
    if(sca.size() != 7) Rcpp::stop("hpp_scalerev: bad 'sca_' argument");
  } else {
    if(!img.hasAttribute("scale")) Rcpp::stop("hpp_scalerev: missing 'scale' attribute");
    sca = img.attr("scale");
    if(sca.size() != 7) Rcpp::stop("hpp_scalerev: bad 'scale' attribute size");
  }
  if(sca[3]) {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = sca[1] - img[i] / sca[2];
  } else {
    for(R_len_t i = 0; i < img.size(); i++) if(img[i] != sca[4]) img[i] = img[i] / sca[2] + sca[0];
  }
  if(img.hasAttribute("class")) img.attr("class") = Rcpp::setdiff(Rcpp::as<Rcpp::CharacterVector>(img.attr("class")), Rcpp::CharacterVector::create("IFCip_rescale"));
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

template <typename T, int RTYPE>
Rcpp::NumericVector scale_T(Rcpp::Vector<RTYPE> img,
                            const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                            const double value = NA_REAL,
                            const int n_lev = 256,
                            const bool invert = false,
                            const bool bin = true,
                            const double clipmin = NA_REAL,
                            const double clipmax = NA_REAL,
                            const uint8_t method = 1) {
  bool as_bin = bin || RTYPE != REALSXP;
  T val =
    RTYPE == INTSXP ? value > 2147483647.0 || value < -2147483647.0 ? NA_INTEGER : value :
      RTYPE == LGLSXP ? Rcpp::traits::is_na<REALSXP>(value) ? NA_LOGICAL : value != 0 :
        RTYPE == RAWSXP ? value > 255.0 || value < 0.0 ? 0.0 : value : value;
  double val_out = RTYPE == REALSXP ? value : Rcpp::traits::is_na<RTYPE>(val) ? NA_REAL : val;
  if(is_scaled_T(img)) {
    Rcpp::NumericVector sca = img.attr("scale");
    if(sca[3] == invert && sca[5] == as_bin && sca[6] == n_lev) {
      return img.attr("scale"); 
    } else {
      hpp_scalerev(img, sca);
    }
    if(img.hasAttribute("class")) img.attr("class") = Rcpp::setdiff(Rcpp::as<Rcpp::CharacterVector>(img.attr("class")), Rcpp::CharacterVector::create("IFCip_rescale"));
  } else {
  }
  if(std::abs(n_lev) < 2) Rcpp::stop("hpp_scale: abs(n_lev) should be at least 2.0");
  Rcpp::NumericVector out = Rcpp::NumericVector::create(_["min"]=R_NegInf, _["max"]=R_PosInf, _["factor"]=1.0, _["invert"]=invert, _["value"]=val_out, _["bin"]=as_bin, _["levels"]=n_lev);
  T MAX_LEV_SCA = as_bin + std::abs(n_lev) - 1;
  T MAX_LEV_BIN = std::abs(n_lev) - 1;
  Rcpp::NumericVector ran = range_T(img, msk_);
  double MN = std::max(ran[0],R_finite(clipmin) && method != 1 ? clipmin : R_NegInf);
  double MX = std::min(ran[1],R_finite(clipmax) && method != 1 ? clipmax : R_PosInf);
  switch(method) {
  case 1 :{
    // don't modify ran
  }
    break;
  case 2 :{
    if(R_finite(clipmin)) ran[0] = clipmin;
    if(R_finite(clipmax)) ran[1] = clipmax;
  }
    break;
  case 3 :{
    ran[0] = MN;
    ran[1] = MX;
  }
    break;
  default: Rcpp::stop("hpp_scale: 'method'[%i] is not supported", method);
  }
  if(n_lev < 0) {
    MAX_LEV_SCA *= -1.0;
    MAX_LEV_BIN *= -1.0;
  }
  out[0] = ran[0];
  out[1] = ran[1];
  out[2] = as_bin ? MAX_LEV_BIN : MAX_LEV_SCA;
  out[2] /= (out[1] - out[0]);
  double K = MAX_LEV_SCA / (out[1] - out[0]), Q = out[0];
  T V = MAX_LEV_BIN, W = invert ? MAX_LEV_BIN : 0, Z = invert ? 0 : MAX_LEV_BIN;
  if(invert) { K *= -1.0; Q = out[1]; }
  if(n_lev < 0) {
    if(msk_.isNotNull()) {
      Rcpp::NumericVector msk = getmask_T(img, msk_);
      if(as_bin) {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? img[i] >= MX ? Z : img[i] <= MN ? W : std::max((T) std::trunc(K * (img[i] - Q)), V) : val; 
      } else {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? img[i] >= MX ? Z : img[i] <= MN ? W : std::max((T) (K * (img[i] - Q)), V) : val;
      }
    } else {
      if(as_bin) {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = img[i] >= MX ? Z : img[i] <= MN ? W : std::max((T) std::trunc(K * (img[i] - Q)), V); 
      } else {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = img[i] >= MX ? Z : img[i] <= MN ? W : std::max((T) (K * (img[i] - Q)), V);
      }
    }
  } else {
    if(msk_.isNotNull()) {
      Rcpp::NumericVector msk = getmask_T(img, msk_);
      if(as_bin) {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? img[i] >= MX ? Z : img[i] <= MN ? W : std::min((T) std::trunc(K * (img[i] - Q)), V) : val;
      } else {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = msk[i] ? img[i] >= MX ? Z : img[i] <= MN ? W : std::min((T) (K * (img[i] - Q)), V) : val;
      }
    } else {
      if(as_bin) {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = img[i] >= MX ? Z : img[i] <= MN ? W : std::min((T) std::trunc(K * (img[i] - Q)), V);
      } else {
        for(R_len_t i = 0; i < img.size(); i++) img[i] = img[i] >= MX ? Z : img[i] <= MN ? W : std::min((T) (K * (img[i] - Q)), V);
      }
    }
  }
  if(img.hasAttribute("class")) {
    img.attr("class") = Rcpp::union_(Rcpp::as<Rcpp::CharacterVector>(img.attr("class")), Rcpp::CharacterVector::create("IFCip_rescale"));
  } else {
    img.attr("class") = Rcpp::CharacterVector::create("IFCip_rescale");
  }
  img.attr("msk") = msk_;
  img.attr("scale") = out; 
  img.attr("levels") = n_lev;
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_scale(SEXP img,
                              const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                              const double value = NA_REAL,
                              const int n_lev = 256,
                              const bool invert = false,
                              const bool bin = true,
                              const double clipmin = NA_REAL,
                              const double clipmax = NA_REAL,
                              const uint8_t method = 1) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return scale_T<int>(as<Rcpp::LogicalVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method);
  case INTSXP : return scale_T<int>(as<Rcpp::IntegerVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method);
  case REALSXP : return scale_T<double>(as<Rcpp::NumericVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method);
  case RAWSXP : return scale_T<uint8_t>(as<Rcpp::RawVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method);
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
                              const bool bin = false,
                              const double clipmin = NA_REAL,
                              const double clipmax = NA_REAL,
                              const uint8_t method = 1) {
  Rcpp::Vector<RTYPE> out = Rcpp::clone(img);
  hpp_scale(out, msk_, value, n_lev, invert, bin, clipmin, clipmax, method); 
  return out;
}

//' @title Image Scaling
//' @name cpp_rescale
//' @description
//' This function is designed to scale a SEXP to [0, n_lev - 1]
//' @param img, a SEXP (logical, raw, integer or numeric) vector or matrix containing image intensity values.
//' @param msk_, a Rcpp::NumericVector with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as \code{true}.
//' Default is \code{R_NilValue}, for using all \code{'img'} elements without masking anything.
//' @param value, a double; it is the replacement value that will be used when \code{'msk'} element is interpreted as \code{false}. Default is \code{NA_REAL}. /!\ Note that it will be silently cast to \code{'img'} type.
//' @param n_lev, an int determining the number of levels used for the computation. Default is \code{256}.
//' @param invert, a bool determining whether \code{'img'} should be scaled from min(\code{'img'}) to max(\code{'img'}) (when \code{false}, [min(\code{'img'}),max(\code{'img'})] becoming [0,sign(\code{'n_lev'})*abs(\code{'n_lev'}-1)]) or inverted (when \code{true}, with [max(\code{'img'}),min(\code{'img'})] rescaled to [0,sign(\code{'n_lev'})*\code{'n_lev'}-1]) values. Default is \code{false}.
//' @param bin, a bool determining whether \code{'img'} should be binned or if scaling should be continuous. Default is \code{false}. It will be forced to \code{true} if \code{'img'} is not of numeric type.
//' @param clipmin, a double, minimal value under which \code{'img'} intensity values will be clipped to. Default is \code{NA_REAL}, to use no minimal clipping.
//' @param clipmax, a double, maximal value above which '\code{'img'} intensity values will be clipped to. Default is \code{NA_REAL}, to use no maximal clipping.
//' @param method, an uint8_t determining how scaling should be applied. Default is \code{1}.
//' -when \code{1}, on [min(\code{'img'}),max(\code{'img'})]
//' -when \code{2}, on [is_finite(\code{'clipmin'}) ? \code{'clipmin'} : min(\code{'img'}), is_finite(\code{'clipmax'}) ? \code{'clipmax'} : max(\code{'img'})]
//' -when \code{3}, on [is_finite(\code{'clipmin'}) ? max(\code{'clipmin'}, min(\code{'img'})) : min(\code{'img'}), is_finite(\code{'clipmax'}) ? min(\code{'clipmax'}, max(\code{'img'})) : max(\code{'img'})]
//' @details when \code{'msk'} is provided it has to be of the same dimensions as \code{'img'}, otherwise an error will be thrown.\cr
//' an error will be thrown also if \code{'msk'} contains non-finite value.\cr
//' \code{'img'} range will be determined based on indices of non 0 finite \code{'msk'} values and only the values in \code{'img'} at those indices will be scaled; the others will be filled with \code{'value'}.
//' @return a SEXP of same type as \code{'img'} with class `IFCip_rescale`
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP hpp_rescale(const SEXP img,
                 const Rcpp::Nullable<Rcpp::NumericVector> msk_ = R_NilValue,
                 const double value = NA_REAL,
                 const int n_lev = 256,
                 const bool invert = false,
                 const bool bin = false,
                 const double clipmin = NA_REAL,
                 const double clipmax = NA_REAL,
                 const uint8_t method = 1) {
  switch( TYPEOF(img) ) {
  case LGLSXP : return rescale_T(as<Rcpp::LogicalVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method); 
  case INTSXP: return rescale_T(as<Rcpp::IntegerVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method); 
  case REALSXP: return rescale_T(as<Rcpp::NumericVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method); 
  case RAWSXP : return rescale_T(as<Rcpp::RawVector>(img), msk_, value, n_lev, invert, bin, clipmin, clipmax, method); 
  default: Rcpp::stop("hpp_rescale: not supported SEXP type in 'img'");
  }
}

#endif