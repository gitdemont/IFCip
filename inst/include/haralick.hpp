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
#include "utils.hpp"
using namespace Rcpp;

//' @title Haralick Co-occurrence Matrix
//' @name cpp_cooc
//' @description
//' This function is designed to compute Haralick co-occurrence matrix
//' @param img a Rcpp::IntegerMatrix of class `IFCip_rescale`, containing image intensity values.
//' @param delta a Rcpp::IntegerVector of column and row shifts. If only one value is provided only row will be shifted.
//' @details See 'Textural Features for Image Classification', Haralick et. al (1979),
//' available at: \url{https://haralick.org/journals/TexturalFeatures.pdf}
//' @return a Rcpp::IntegerMatrix Gray-Level Co-occurrence Matrices.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_cooc(const Rcpp::IntegerMatrix img,
                             const Rcpp::IntegerVector delta) {
  if(!Rf_inherits(img, "IFCip_rescale")) {
    Rcpp::stop("hpp_cooc: 'img' should be of class `IFCip_rescale`");
  }
  if(!img.hasAttribute("levels")) {
    Rcpp::stop("hpp_cooc: 'img' should have `levels` attribute");
  }
  
  uint16_t depth = as<int>(img.attr("levels"));
  R_len_t mat_r = img.nrow();
  R_len_t mat_c = img.ncol();
  
  R_len_t dy = delta[0];
  R_len_t dx = delta.size() > 1 ? delta[1] : 0;
  
  if(std::abs(dx) >= mat_c) {
    Rcpp::warning("hpp_cooc: 'dx'[%i] should be smaller than 'img' ncol [%i,%i]", dx, mat_r, mat_c);
  }
  if(std::abs(dy) >= mat_r) {
    Rcpp::warning("hpp_cooc: 'dy'[%i] should be smaller than 'img' nrow [%i,%i]", dy, mat_r, mat_c);
  }
  
  // Ensure no NA will be encountered
  Rcpp::LogicalMatrix msk = get_mask(img.attr("msk"), mat_r, mat_c);
  Rcpp::LogicalMatrix M = Rcpp::no_init(mat_r, mat_c);
  for(R_len_t i = 0; i < img.size(); i++) {
    if(img[i] == NA_INTEGER ||
       msk[i] == NA_LOGICAL) {
      M[i] = false;
    } else {
      M[i] = msk[i];
    }
  }
  
  // Non normalized Gray-Level Co-occurence Matrix, GLCM
  Rcpp::IntegerMatrix out(depth, depth);
  out.attr("class") = "IFCip_cooc";
  out.attr("delta") = delta;
  
  // Initialize count
  R_len_t K = 0;
  
  // Compute co-occurence
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t t_col = i_col + dx;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      if(M(i_row, i_col)) {
        R_len_t t_row = i_row + dy;
        if(t_row >= 0 && t_row < mat_r && 
           t_col >= 0 && t_col < mat_c) {
          if(M(t_row, t_col)) {
            if(img(i_row, i_col) >= 0 &&
               img(t_row, t_col) >= 0 &&
               img(i_row, i_col) < depth &&
               img(t_row, t_col) < depth) {
              out(img(i_row, i_col), img(t_row, t_col))++;
              out(img(t_row, t_col), img(i_row, i_col))++;
              K++;
            } else {
              Rcpp::stop("hpp_cooc: values are out of range");
            }
          }
        }
      }
    }
  }
  
  out.attr("count") = K;
  return out;
}

//' @title Haralick Features
//' @name cpp_h_features
//' @description
//' This function is designed to compute Haralick's features
//' @param cooc a Rcpp::NumericMatrix of class `IFCip_cooc`, co-occurrence matrix to compute Haralick's features from.
//' @param invariant a bool, whether to compute invariant Haralick's texture features. Default is false.
//' Not yet supported.
//' @details Haralick's invariant texture features are described in Löfstedt T, Brynolfsson P, Asklund T, Nyholm T, Garpebring A (2019) Gray-level invariant Haralick texture features.
//' PLoS ONE 14(2): e0212110. \doi{10.1371/journal.pone.0212110}
//' @return a Rcpp::NumericVector of Haralick's texture features
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_h_features(const Rcpp::IntegerMatrix cooc,
                                   const bool invariant = false) {
  if(!Rf_inherits(cooc, "IFCip_cooc")) {
    Rcpp::stop("hpp_h_features: 'cooc' should be of class `IFCip_cooc`");
  }
  R_len_t N = cooc.ncol();
  if(N < 4) {
    Rcpp::stop("hpp_h_features: 'cooc' is too small");
  }
  if(N != cooc.nrow() || (N % 2)) {
    Rcpp::stop("hpp_h_features: 'cooc' should be a square matrix");
  }
  double S = cooc.attr("count"); S *= 2;
  double mu_x, mu_y, sig_x, sig_y, sig_xy, mu_xpy, mu_xmy,
         HX, HY, HXY, HXY1, HXY2;
  double con, cor, d_ent, d_var, dis, nrj, ent, hom,
         imc1, imc2, i_dif, p_max, s_ent, s_var;
  double mu, c_pro, c_sha, s_sqr;
  
  // unknown lambda dependent + computationally instable Q
  // lambda, Q, mcc
  
  Rcpp::NumericMatrix p(N + 1, N + 1);
  // for invariant
  double d, d_xpy, d_ij;
  if(invariant) {
    // TODO check to add invariant
    // d = 1/N, d_xpy = 1 / (2 * N - 1), d_ij = 1 / (N * N) / S;
    d = 1, d_xpy = 1, d_ij = 1 / S;
  } else {
    d = 1, d_xpy = 1, d_ij = 1 / S;
  }
  // copy and weight cooc to p, in addition add +1 index offset
  for(R_len_t i = 0, i_col1 = 1; i_col1 <= N; i_col1++) {
    for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++, i++) {
      p(i_row1, i_col1) = cooc[i] * d_ij;
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
  
  mu_x = 0.0;
  HX   = 0.0;
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    mu_x += (i_row1) * d * p_x[i_row1];
    if(p_x[i_row1]) HX -= p_x[i_row1] * std::log2(p_x[i_row1]);
  }
  
  mu_y = 0.0;
  HY   = 0.0;
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    mu_y += (i_col1) * d * p_y[i_col1];
    if(p_y[i_col1]) HY -= p_y[i_col1] * std::log2(p_y[i_col1]);
  }
  
  sig_x = 0.0;
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    sig_x += std::pow((i_row1) * d - mu_x, 2.0) * p_x[i_row1];
  }
  sig_x = std::sqrt(sig_x);
  
  sig_y = 0.0;
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    sig_y += std::pow((i_col1) * d - mu_y, 2.0) * p_y[i_col1];
  }
  sig_y = std::sqrt(sig_y);
  sig_xy = sig_x * sig_y;
  
  mu_xpy = 0.0;
  mu_xmy = 0.0;
  d_ent  = 0.0;
  s_ent  = 0.0;
  d_var  = 0.0;
  s_var  = 0.0;
  if(invariant) {
    for(double k = 2; k <= 2 * N; k++) {
      // 11
      mu_xpy += 2 * (k - 1) * d_xpy * p_xpy[k];
      // 12
      if(p_xpy[k]) s_ent -= p_xpy[k] * std::log2(p_xpy[k]);
    }
    for(double k = 0; k < N; k++) {
      // 20
      mu_xmy += (k + 1) * d * p_xmy[k];
      // 09
      if(p_xmy[k]) d_ent -= p_xmy[k] * std::log2(p_xmy[k]);
    }
    // 10
    for(double k = 0; k < N; k++) d_var += std::pow((k + 1) * d - mu_xmy, 2.0) * p_xmy[k];
    // 13
    for(double k = 2; k <= 2 * N; k++) s_var += std::pow(2 * (k - 1) * d_xpy - mu_xpy, 2.0) * p_xpy[k];
  } else {
    for(double k = 2; k <= 2 * N; k++) {
      // 11
      mu_xpy += k * p_xpy[k];
      // 12
      if(p_xpy[k]) s_ent -= p_xpy[k] * std::log2(p_xpy[k]);
    }
    for(double k = 0; k < N; k++) {
      // 20
      mu_xmy += k * p_xmy[k];
      // 09
      if(p_xmy[k]) d_ent -= p_xmy[k] * std::log2(p_xmy[k]);
    }
    // 10
    for(double k = 0; k < N; k++) d_var += std::pow(k - mu_xmy, 2.0) * p_xmy(k);
    // 13
    for(double k = 2; k <= 2 * N; k++) s_var += std::pow(k - mu_xpy, 2.0) * p_xpy(k);
  }
  
  HXY  = 0.0;
  HXY1 = 0.0;
  HXY2 = 0.0;
  // Rcpp::NumericMatrix Q(N, N);
  // Q.fill(0.0);
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
      if(p(i_row1, i_col1)) HXY += p(i_row1, i_col1) * std::log2(p(i_row1, i_col1));
      if(p_x[i_row1] * p_y[i_col1] != 0.0) {
        HXY1 += p(i_row1, i_col1) * std::log2(p_x[i_row1] * p_y[i_col1]);
        HXY2 += p_x[i_row1] * p_y[i_col1] * std::log2(p_x[i_row1] * p_y[i_col1]);
      }
      // for(k = 1; k <= N; k++) {
      //   Q(i_row, i_col) += p(i_row, k - 1) * p(i_col, k - 1) / (p_x(i_row) * p_y(k));
      // }
    }
  }
  
  imc1 = (HXY - HXY1) / std::max(HX, HY);
  imc2 = std::sqrt(1 - std::exp(2*(HXY2 - HXY))); //should be std::sqrt(1 - std::exp(-2*(HXY2 - HXY)))
  // mcc = std::sqrt(lambda(Q(i_col1, i_row1)))
  
  c_pro = 0.0;
  c_sha = 0.0;
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
    double i_cold = (i_col1) * d;
    for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
      double i_rowd = (i_row1) * d;
      double pij = p(i_row1, i_col1);
      double m = i_rowd - i_cold;
      double a = std::abs(m);
      double b = std::pow(m, 2.0);
      // 01: H Contrast
      con += b * pij;
      // 02: CORRELATION
      cor += (i_rowd - mu_x) * (i_cold - mu_y) * pij / sig_xy;
      // 03: CONTRAST
      dis += a * pij;
      // 04: ANGULAR SECOND MOMENT
      nrj += std::pow(pij, 2.0); 
      // 05: ENTROPY
      if(pij) ent -= pij * std::log2(pij);
      // 06: HOMOGENEITY (it is INVERSE DIFFERENCE MOMENT in Haralick paper)
      hom += pij / (1 + b);
      // 07: MAXIMUM PROBABILITY
      p_max = std::max(pij, p_max);
      // 08: INVERSE DIFFERENT MOMENT
      i_dif += pij / (1 + a);
      
      double v = i_rowd + i_cold - 2 * mu;
      double w = std::pow(v, 3.0) * pij;
      // 14: SUM OF SQUARE
      s_sqr += (i_rowd - mu) * (i_rowd - mu) * pij;
      // 15: CLUSTER PROMINENCE
      c_pro += w;
      // 16: CLUSTER SHAPE
      c_sha += v * w;
    }
  }
  
  // to handle potential division by 0
  if(std::max(HX, HY) == 0.0) imc1 = 0.0;
  if(imc2 == R_NaN) imc2 = 0.0;
  if((sig_x == 0.0) || (sig_y == 0.0)) cor = 0.0;
  
  return Rcpp::NumericVector::create(_["H Contrast"] = con / N / N, // IDEAS H_Contrast con is divided by probability since IDEAS claims contrast is [0,1]
                                     _["H Correlation"] = cor, // IDEAS H_Correlation;
                                     _["dissimilarity"] = dis,
                                     _["H Energy"] = nrj, // IDEAS H_Energy
                                     _["H Entropy"] = ent, // IDEAS H_Entropy
                                     _["homogeneity"] = hom,
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
