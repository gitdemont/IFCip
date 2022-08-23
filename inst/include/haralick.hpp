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

#ifndef IFCIP_HARALICK_HPP
#define IFCIP_HARALICK_HPP

#include <Rcpp.h>
// [[Rcpp::depends(IFC)]]
#include <utils.hpp>
using namespace Rcpp;

//' @title Right Shift Matrix
//' @name cpp_R_shift_M
//' @description
//' This function is designed to right shift bits of a matrix to [0, 2^bits - 1]
//' @param mat a Rcpp::IntegerMatrix, containing image intensity values.
//' @param bits uint8_t number of bit to shift matrix values. Default is 4. Allowed are [2,10].
//' If shifted value does not respect [0, 2^bits - 1] an error is thrown.
//' @return an Rcpp::IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_R_shift_M(const Rcpp::IntegerMatrix mat,
                                  const uint8_t bits = 4) {
  if(bits > 12 || bits < 2) {
    Rcpp::stop("hpp_R_shift_M: 'bits' should be a [2,10] integer");
  }
  int min = 0, max = std::pow(2.0, bits) - 1;
  Rcpp::IntegerMatrix out(mat.nrow(), mat.ncol());
  for(R_len_t i = 0; i < mat.size(); i++) {
    out[i] = mat[i] >> (bits - 1);
    if(out[i] < min || out[i] > max) {
      Rcpp::stop("hpp_R_shift_M: shifted value does not respect allowed interval");
    }
  }
  out.attr("class") = "IFCip_rescale";
  out.attr("bits") = bits;
  return out;
}

//' @title hpp_rescale_M
//' @name cpp_rescale_M
//' @description
//' This function is designed to rescale a matrix to [0, 2^bits - 1]
//' @param img a Rcpp::IntegerMatrix, containing image intensity values.
//' @param msk_, a Rcpp::NumericMatrix with finite values. Non-finite values will trigger an error. All non 0 values will be interpreted as true.
//' Default is R_NilValue, for using all 'img' elements without masking anything.
//' @param bits uint8_t number of bit to shift matrix values. Default is 4. Allowed are [2,10].
//' @return a Rcpp::IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_rescale_M(const Rcpp::IntegerMatrix img,
                                  const Rcpp::Nullable<Rcpp::NumericMatrix> msk_ = R_NilValue,
                                  const uint8_t bits = 4) {
  if(bits > 12 || bits < 2) {
    Rcpp::stop("hpp_rescale_M: 'bits' should be a [2,10] integer");
  }
  Rcpp::IntegerMatrix out = Rcpp::clone(img);
  Rcpp::NumericVector sca = hpp_scale(out, msk_, NA_REAL, std::pow(2.0, bits), true);
  out.attr("class") = "IFCip_rescale";
  out.attr("scale") = sca;
  out.attr("bits") = bits;
  return out;
}

//' @title Haralick Co-occurrence Matrix
//' @name cpp_cooc
//' @description
//' This function is designed to compute Haralick co-occurrence matrix
//' @param img a Rcpp::IntegerMatrix of class `IFCip_rescale`, containing image intensity values.
//' @param msk a LogicalMatrix, containing mask.
//' @param delta uint8_t offset from which co-occurence has to be computed to. Default is 1.
//' @details See 'Textural Features for Image Classification', Haralick et. al (1979),
//' available at: \url{https://haralick.org/journals/TexturalFeatures.pdf}
//' @return a list whose members are normalized Gray-Level Co-occurrence Matrices at angles 0, 45, 90 and 315.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_cooc(const Rcpp::IntegerMatrix img,
                    const Rcpp::LogicalMatrix msk,
                    const uint8_t delta = 1) {
  if(!Rf_inherits(img, "IFCip_rescale")) {
   Rcpp::stop("hpp_cooc: 'img' should be of class `IFCip_rescale`");
  }

  uint16_t depth = std::pow(2.0, as<int>(img.attr("bits")));
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  if(mat_r != msk.nrow() || mat_c != msk.ncol()) {
    Rcpp::stop("hpp_cooc: 'img' and 'msk' should have same dimensions");
  }
  if((delta >= mat_c) || (delta >= mat_r)) {
    Rcpp::warning("hpp_cooc: 'delta'[%i] should be smaller than 'img' size [%i,%i]", delta, mat_r, mat_c);
  }
  
  // Non normalized Gray-Level Co-occurence Matrix, GLCM
  // DO NOT USE IntegerMatrix 
  Rcpp::NumericMatrix G1(depth, depth);
  Rcpp::NumericMatrix G2 = Rcpp::clone(G1);
  Rcpp::NumericMatrix G3 = Rcpp::clone(G1);
  Rcpp::NumericMatrix G4 = Rcpp::clone(G1);
  
  // Initialize normalization constant
  R_len_t K1 = 0, K2 = 0, K3 = 0, K4 = 0;
  
  // Compute co-occurence
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t t_col = i_col + delta;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      if(msk(i_row, i_col)) {
        R_len_t t_row = i_row + delta;
        if(t_row < mat_r) {
          // 270 identical to 90
          if(msk(t_row, i_col)) {
            G1(img(i_row, i_col), img(t_row, i_col))++;
            G1(img(t_row, i_col), img(i_row, i_col))++;
            K1 += 2;
          }
          // 315 identical to 135
          if(t_col < mat_c &&
             msk(t_row, t_col)) {
            G2(img(i_row, i_col), img(t_row, t_col))++;
            G2(img(t_row, t_col), img(i_row, i_col))++;
            K2 += 2;
          }
        }
        // 0 identical to 180
        if(t_col < mat_c &&
           msk(i_row, t_col)) {
          G3(img(i_row, i_col), img(i_row, t_col))++;
          G3(img(i_row, t_col), img(i_row, i_col))++;
          K3 += 2;
        }
        // 45 identical to 225
        t_row = i_row - delta;
        if(t_row >= 0 &&
           t_col < mat_c &&
           msk(t_row, t_col)) {
          G4(img(i_row, i_col), img(t_row, t_col))++;
          G4(img(t_row, t_col), img(i_row, i_col))++;
          K4 += 2;
        }
      }
    }
  }
  
  // Normalization
  Rcpp::NumericMatrix P1 = Rcpp::no_init_matrix(depth, depth);
  P1.attr("class") = "IFCip_cooc";
  P1.attr("delta") = delta;
  Rcpp::NumericMatrix P2 = Rcpp::clone(P1);
  Rcpp::NumericMatrix P3 = Rcpp::clone(P1);
  Rcpp::NumericMatrix P4 = Rcpp::clone(P1);
  R_len_t dd = P1.size() - 1;
  if(K1 > 0){
    for(R_len_t i = 0; i <= dd; i++) P1[i] = G1[dd - i] / K1;
  } else {
    P1.fill(0.0);
  }
  if(K2 > 0) {
    for(R_len_t i = 0; i <= dd; i++) P2[i] = G2[dd - i] / K2;
  } else {
    P2.fill(0.0);
  }
  if(K3 > 0) {
    for(R_len_t i = 0; i <= dd; i++) P3[i] = G3[dd - i] / K3;
  } else {
    P3.fill(0.0);
  }
  if(K4 > 0) {
    for(R_len_t i = 0; i <= dd; i++) P4[i] = G4[dd - i] / K4;
  } else {
    P4.fill(0.0);
  }
  return Rcpp::List::create(_["P+x0y"] = P3,  // 000, 180
                            _["P+x-y"] = P4,  // 045, 225
                            _["P0x+y"] = P1,  // 090, 270
                            _["P+x+y"] = P2); // 135, 315
}

//' @title Haralick Features
//' @name cpp_h_features
//' @description
//' This function is designed to compute Haralick's features
//' @param cooc a Rcpp::NumericMatrix of class `IFCip_cooc`, normalized co-occurrence matrix to compute Haralick's features from.
//' @param invariant a bool, whether to compute invariant Haralick's texture features. Default is false.
//' Not yet supported.
//' @details Haralick's invariant texture features are described in Löfstedt T, Brynolfsson P, Asklund T, Nyholm T, Garpebring A (2019) Gray-level invariant Haralick texture features.
//' PLoS ONE 14(2): e0212110. \doi{10.1371/journal.pone.0212110}
//' @return a Rcpp::NumericVector of Haralick's texture features
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_h_features(const Rcpp::NumericMatrix cooc,
                                   const bool invariant = false) {
  if(!Rf_inherits(cooc, "IFCip_cooc")) {
    Rcpp::stop("hpp_h_features: 'cooc' should be of class `IFCip_cooc`");
  }
  R_len_t N = cooc.ncol();
  double mu_x, mu_y, sig_x, sig_y, mu_xpy, mu_xmy,
         HX, HY, HXY, HXY1, HXY2;
  double a_cor, con, cor, d_ent, d_var, dis, nrj, ent, hom,
         imc1, imc2, i_dif, p_max, s_ent, s_var;
  double mu, c_pro, c_sha, s_sqr;
  double eps = 0.0000001; // to prevent log2(0)
  
  // unknown lambda dependent + computationally instable Q
  // lambda, Q, mcc
  
  Rcpp::NumericMatrix p(N + 1, N + 1);
  // for invariant
  double d, d_xpy, d_ij;
  if(invariant) {
    // TODO check to add invariant
    // d = 1/N, d_xpy = 1 / (2 * N - 1), d_ij = 1 / (N * N);
    d = 1, d_xpy = 1, d_ij = 1;
  } else {
    d = 1, d_xpy = 1, d_ij = 1;
  }
  // copy and weight cooc to p, in addition add +1 index offset
  for(R_len_t i_col = 0; i_col < N; i_col++) {
    for(R_len_t i_row = 0; i_row < N; i_row++) {
      p(i_row + 1, i_col + 1) = cooc(i_row, i_col) * d_ij;
    }
  }
  
  // compute variables
  Rcpp::NumericVector p_x(N + 1, 0.0);
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    p_x[i_row1] += sum(p(i_row1, Rcpp::_));
  }
  Rcpp::NumericVector p_y(N + 1, 0.0);
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    p_y[i_col1] += sum(p(Rcpp::_, i_col1));
  }
  
  Rcpp::NumericVector p_xpy(2 * N + 2, 0.0);
  Rcpp::NumericVector p_xmy(N, 0.0);
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
      p_xpy[i_row1 + i_col1] += p(i_row1, i_col1) * d;
      p_xmy[std::abs(i_row1 - i_col1)] += p(i_row1, i_col1) * d;
    }
  }
  
  mu_x = 1.0;
  HX   = 0.0;
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    mu_x += (i_row1 - 1) * d * p_x[i_row1];
    HX += p_x[i_row1] ? p_x[i_row1] * std::log2(p_x[i_row1]) : p_x[i_row1] * std::log2(p_x[i_row1] + eps);
  }
  HX *= -1.0;
  
  mu_y = 1.0;
  HY   = 0.0;
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    mu_y += (i_col1 - 1) * d * p_y[i_col1];
    HY += p_y[i_col1] ? p_y[i_col1] * std::log2(p_y[i_col1]) : p_y[i_col1] * std::log2(p_y[i_col1] + eps);
  }
  HY *= -1.0;
  
  sig_x = 0.0;
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    sig_x += std::pow(i_row1 * d - mu_x, 2.0) * p_x[i_row1];
  }
  sig_x = std::sqrt(sig_x);
  
  sig_y = 0.0;
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    sig_y += std::pow(i_col1 * d - mu_y, 2.0) * p_y[i_col1];
  }
  sig_y = std::sqrt(sig_y);
  
  mu_xpy = 0.0;
  mu_xmy = 0.0;
  d_ent  = 0.0;
  s_ent  = 0.0;
  d_var  = 0.0;
  s_var  = 0.0;
  if(invariant) {
    for(double k = 2; k <= 2 * N; k++) {
      mu_xpy += 2 * (k - 1) * d_xpy * p_xpy[k];
      s_ent += p_xpy[k] ? p_xpy[k] * std::log2(p_xpy[k]) : p_xpy[k] * std::log2(p_xpy[k] + eps);
    }
    for(double k = 0; k < N; k++) {
      mu_xmy += (k + 1) * d * p_xmy[k];
      d_ent += p_xmy[k] ? p_xmy[k] * std::log2(p_xmy[k]) : p_xmy[k] * std::log2(p_xmy[k] + eps);
    }
    for(double k = 0; k < N; k++) d_var += std::pow((k + 1) * d - mu_xmy, 2.0) * p_xmy[k];
    for(double k = 2; k <= 2 * N; k++) s_var += std::pow(2 * (k - 1) * d_xpy - mu_xpy, 2.0) * p_xpy[k];
  } else {
    for(double k = 2; k <= 2 * N; k++) {
      mu_xpy += k * p_xpy[k];
      s_ent += p_xpy[k] ? p_xpy[k] * std::log2(p_xpy[k]) : p_xpy[k] * std::log2(p_xpy[k] + eps);
    }
    for(double k = 0; k < N; k++) {
      mu_xmy += k * p_xmy[k];
      d_ent += p_xmy[k] ? p_xmy[k] * std::log2(p_xmy[k]) : p_xmy[k] * std::log2(p_xmy[k] + eps);
    }
    for(double k = 0; k < N; k++) d_var += std::pow(k - mu_xmy, 2.0) * p_xmy(k);
    for(double k = 2; k <= 2 * N; k++) s_var += std::pow(k - mu_xpy, 2.0) * p_xpy(k);
  }
  s_ent *= -1.0;
  d_ent *= -1.0;
  
  HXY  = 0.0;
  HXY1 = 0.0;
  HXY2 = 0.0;
  // Rcpp::NumericMatrix Q(N, N);
  // Q.fill(0.0);
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
      HXY += p(i_row1, i_col1) ? p(i_row1, i_col1) * std::log2(p(i_row1, i_col1)) : p(i_row1, i_col1) * std::log2(p(i_row1, i_col1) + eps);
      if(p_x[i_row1] * p_y[i_col1] == 0.0) {
        HXY1 += p(i_row1, i_col1) * std::log2(p_x[i_row1] * p_y[i_col1] + eps);
        HXY2 += p_x[i_row1] * p_y[i_col1] * std::log2(p_x[i_row1] * p_y[i_col1] + eps);
      } else {
        HXY1 += p(i_row1, i_col1) * std::log2(p_x[i_row1] * p_y[i_col1]);
        HXY2 += p_x[i_row1] * p_y[i_col1] * std::log2(p_x[i_row1] * p_y[i_col1]);
      }
      // for(k = 1; k <= N; k++) {
      //   Q(i_row, i_col) += p(i_row, k - 1) * p(i_col, k - 1) / (p_x(i_row) * p_y(k));
      // }
    }
  }
  HXY  *= -1.0;
  HXY1 *= -1.0;
  HXY2 *= -1.0;
  
  imc1 = (HXY - HXY1) / std::max(HX, HY);
  imc2 = std::sqrt(1 - std::exp(-2*(HXY2 - HXY)));
  // mcc = std::sqrt(lambda(Q(i_col1, i_row1)))
  
  c_pro = 0.0;
  c_sha = 0.0;
  a_cor = 0.0;
  con   = 0.0;
  cor   = 0.0;
  dis   = 0.0;
  nrj   = 0.0;
  ent   = 0.0;
  hom   = 0.0;
  i_dif = 0.0;
  s_sqr = 0.0;
  mu = (mu_x + mu_y) / 2;
  p_max = p(0, 0);
  
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
      // AUTOCORRELATION
      a_cor += (i_row1 * d * i_col1 * d) * p(i_row1, i_col1); 
      // CLUSTER PROMINENCE
      c_pro += std::pow(i_row1 * d + i_col1 * d - 2 * mu, 3.0) * p(i_row1, i_col1);
      // CLUSTER SHAPE
      c_sha += std::pow(i_row1 * d + i_col1 * d - 2 * mu, 4.0) * p(i_row1, i_col1);
      // H Contrast
      con += std::pow(i_row1 * d - i_col1 * d, 2.0) * p(i_row1, i_col1);
      // CORRELATION
      cor += (i_row1 * d - mu_x) * (i_col1 * d - mu_y) * p(i_row1, i_col1) / sig_x / sig_y;
      // CONTRAST
      dis += std::abs(i_row1 * d - i_col1 * d) * p(i_row1, i_col1);
      // ANGULAR SECOND MOMENT
      nrj += std::pow(p(i_row1, i_col1), 2.0); 
      // ENTROPY
      ent += p(i_row1, i_col1) ? p(i_row1, i_col1) * std::log2(p(i_row1, i_col1)) : p(i_row1, i_col1) * std::log2(p(i_row1, i_col1) + eps);
      // HOMOGENEITY (it is INVERSE DIFFERENCE MOMENT in Haralick paper)
      hom += p(i_row1, i_col1) / (1 + std::pow(i_row1 * d - i_col1 * d, 2.0));
      // INVERSE DIFFERENT MOMENT
      i_dif += p(i_row1, i_col1) / (1 + std::abs(i_row1 * d - i_col1 * d));
      // MAXIMUM PROBABILITY
      p_max = std::max(p(i_row1, i_col1), p_max);
      // SUM OF SQUARE
      s_sqr += (i_row1 - mu) * (i_row1 - mu) * p(i_row1, i_col1);
    }
  }
  
  // to handle potential division by 0
  if(std::max(HX, HY) == 0.0) imc1 = 0.0; 
  if((sig_x == 0.0) || (sig_y == 0.0)) cor = 0.0;
  
  return Rcpp::NumericVector::create(_["autocorrelation"] = a_cor,
                                     _["H Contrast"] = con / N / N, // IDEAS H_Contrast con is divided by probability since IDEAS claims contrast is [0,1]
                                     _["H Correlation"] = cor, // IDEAS H_Correlation;
                                     _["dissimilarity"] = dis,
                                     _["H Energy"] = nrj, // IDEAS H_Energy
                                     _["H Entropy"] = -1.0 * ent, // IDEAS H_Entropy
                                     _["homogeneity"] = hom,
                                     // _["H Homogeneity"] = hom2, // IDEAS H_Homogeneity
                                     _["maximum probability"] = p_max,
                                     _["inverse difference"] = i_dif,
                                     _["difference entropy"] = d_ent,
                                     _["difference variance"] = d_var,
                                     _["H Homogeneity"] = mu_xpy, //IDEAS H_Homogeneity = sum average
                                     _["sum entropy"] = s_ent,
                                     _["sum variance"] = s_var,
                                     _["H Variance"] = s_sqr,  // IDEAS H_Variance
                                     _["cluster prominence"] = c_pro,
                                     _["cluster shade"] = c_sha,
                                     _["imc1"] = imc1,
                                     _["imc2"] = imc2);
}

#endif
