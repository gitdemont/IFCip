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
//' @param mat a Rcpp::NumericMatrix, containing image intensity values.
//' @param bits uint8_t number of bit to shift matrix values. Default is 4. Allowed are [2,10].
//' Rescaled values will normalized to [0, 2^bits - 1]
//' @return an Rcpp::IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_rescale_M(const Rcpp::NumericMatrix mat,
                                  const uint8_t bits = 4) {
  if(bits > 12 || bits < 2) {
    Rcpp::stop("hpp_rescale_M: 'bits' should be a [2,10] integer");
  }
  uint16_t d = std::pow(2.0, bits);
  uint16_t dd = d - 1;
  Rcpp::IntegerMatrix out = no_init_matrix(mat.nrow(), mat.ncol());
  Rcpp::NumericVector ran = hpp_check_range(hpp_check_range(mat));
  double diff = ran[1] - ran[0];
  for(R_len_t i = 0; i < mat.size(); i++) {
    // normalize to [0-2^bits]
    // besides setting out as IntegerMatrix will automatically factorize out[i]
    // to an interval equal to:
    // -0 for [0, 1[,
    // -1 for [1, 2[, 
    // -and so on up to [2^bits-1, 2^bits[,
    // and a last one for exact [2^bits]
    out[i] = d * (mat[i] - ran[0]) / diff;
    // we need to pass 2^bits to the former interval [2^bits-1, 2^bits[
    if(out[i] == d) out[i] = dd;
  }
  out.attr("class") = "IFCip_rescale";
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
  if(delta >= mat_c) {
    Rcpp::stop("hpp_cooc: 'delta' should be smaller than 'img' width");
  }
  if(delta >= mat_r) {
    Rcpp::stop("hpp_cooc: 'delta' should be smaller than 'img' height");
  }
  
  // Non normalized Gray-Level Co-occurence Matrix, GLCM
  // DO NOT USE IntegerMatrix 
  Rcpp::NumericMatrix G1 = Rcpp::no_init_matrix(depth, depth);
  G1.fill(0);
  Rcpp::NumericMatrix G2 = Rcpp::clone(G1);
  Rcpp::NumericMatrix G3 = Rcpp::clone(G1);
  Rcpp::NumericMatrix G4 = Rcpp::clone(G1);
  
  // Initialize normalization constant
  R_len_t K1 = 0, K2 = 0, K3 = 0, K4 = 0;
  
  // Compute co-occurence
  for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
    R_len_t t_col = i_col + delta;
    for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
      R_len_t t_row = i_row + delta;
      // 270 identical to 90
      if(t_row < mat_r) {
        if(msk(i_row, i_col) && msk(t_row, i_col)) {
          G1(img(i_row, i_col), img(t_row, i_col))++;
          G1(img(t_row, i_col), img(i_row, i_col))++;
          K1 += 2;
        }
      }
      // 315 identical to 135
      if(t_row < mat_r && t_col < mat_c) {
        if(msk(i_row, i_col) && msk(t_row, t_col)) {
          G2(img(i_row, i_col), img(t_row, t_col))++;
          G2(img(t_row, t_col), img(i_row, i_col))++;
          K2 += 2;
        }
      }
      // 0 identical to 180
      if(t_col < mat_c) {
        if(msk(i_row, i_col) && msk(i_row, t_col)) {
          G3(img(i_row, i_col), img(i_row, t_col))++;
          G3(img(i_row, t_col), img(i_row, i_col))++;
          K3 += 2;
        }
      }
      // 45 identical to 225
      t_row = i_row - delta;
      if(t_row >= 0 && t_col < mat_c) {
        if(msk(i_row, i_col) && msk(t_row, t_col)) {
          G4(img(i_row, i_col), img(t_row, t_col))++;
          G4(img(t_row, t_col), img(i_row, i_col))++;
          K4 += 2;
        }
      }
    }
  }
  
  // Normalization
  Rcpp::NumericMatrix P1 = Rcpp::no_init_matrix(depth, depth);
  Rcpp::NumericMatrix P2 = Rcpp::no_init_matrix(depth, depth);
  Rcpp::NumericMatrix P3 = Rcpp::no_init_matrix(depth, depth);
  Rcpp::NumericMatrix P4 = Rcpp::no_init_matrix(depth, depth);
  // Rcpp::NumericMatrix P5 = Rcpp::no_init_matrix(depth, depth);
  // Rcpp::NumericMatrix P6 = Rcpp::no_init_matrix(depth, depth);
  
  for(R_len_t i_col = 0; i_col < depth; i_col++) {
    for(R_len_t i_row = 0; i_row < depth; i_row++) {
      /* 
       It could have been faster to determine mean but we need the 4 cooc matrices
       to compute std afterwards
       Rcpp::NumericVector V = Rcpp::NumericVector::create(G1(i_row, i_col)/ K1, 
                                                     G2(i_row, i_col)/ K2, 
                                                     G3(i_row, i_col)/ K3, 
                                                     G4(i_row, i_col)/ K4);
       P1(i_row, i_col) = V[0];
       P2(i_row, i_col) = V[1];
       P3(i_row, i_col) = V[2];
       P4(i_row, i_col) = V[3];
       P5(i_row, i_col) = Rcpp::mean(V);
       P6(i_row, i_col) = Rcpp::sd(V);
       Rcpp::List out = Rcpp::List::create(_["P+x0y"] = P3, // 000, 180
                                           _["P+x-y"] = P4, // 045, 225
                                           _["P0x+y"] = P1, // 090, 270
                                           _["P+x+y"] = P2, // 135, 315
                                           _["mean"] = P5,  // all directions
                                           _["sd"] = P6);
      */
      P1(i_row, i_col) = G1(i_row, i_col)/ K1;
      P2(i_row, i_col) = G2(i_row, i_col)/ K2;
      P3(i_row, i_col) = G3(i_row, i_col)/ K3;
      P4(i_row, i_col) = G4(i_row, i_col)/ K4;
    }
  }
  P1.attr("class") = "IFCip_cooc";
  P1.attr("delta") = delta;
  P2.attr("class") = "IFCip_cooc";
  P2.attr("delta") = delta;
  P3.attr("class") = "IFCip_cooc";
  P3.attr("delta") = delta;
  P4.attr("class") = "IFCip_cooc";
  P4.attr("delta") = delta;
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
//' PLoS ONE 14(2): e0212110. \url{https://doi.org/10.1371/journal.pone.0212110}
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
  double a_cor, con, cor, d_ent, d_var, dis, nrj, ent, hom, // hom2, 
         imc1, imc2, i_dif, p_max, s_ent, s_var;
  double mu, c_pro, c_sha, s_sqr;
  
  // unknown lambda dependent + computationally instable Q
  // lambda, Q, mcc
  
  Rcpp::NumericMatrix p(N, N);
  // for invariant
  double d, d_xpy, d_ij;
  if(invariant) {
    // TODO check to add invariant
    // d = 1/N, d_xpy = 1 / (2 * N - 1), d_ij = 1 / (N * N);
    d = 1, d_xpy = 1, d_ij = 1;
    for(R_len_t i_col = 0; i_col < N; i_col++) {
      for(R_len_t i_row = 0; i_row < N; i_row++) {
        p(i_row, i_col) = cooc(i_row, i_col) * d_ij;
      }
    }
  } else {
    d = 1, d_xpy = 1, d_ij = 1;
    p = Rcpp::clone(cooc);
  }
  
  // compute variables
  Rcpp::NumericVector p_x(N + 1, 0.0);
  for(R_len_t i_col = 0; i_col < N; i_col++) {
    p_x[i_col + 1] += sum(p(Rcpp::_, i_col));
  }
  
  Rcpp::NumericVector p_y(N + 1, 0.0);
  for(R_len_t i_row = 0; i_row < N; i_row++) {
    p_y[i_row + 1] += sum(p(i_row, Rcpp::_));
  }
  Rcpp::NumericVector p_xpy(2 * N + 2, 0.0);
  Rcpp::NumericVector p_xmy(N, 0.0);
  for(R_len_t i_col = 0; i_col < N; i_col++) {
    for(R_len_t i_row = 0; i_row < N; i_row++) {
      p_xpy[i_row + 1 + i_col + 1] += p(i_row, i_col) * d;
      // p_xmy[std::abs(i_row1 - i_col1)] += p(i_row, i_col);
      // but i_row1 - icol1 = i_row + 1 - (i_col + 1) = i_row - i_col;
      p_xmy[std::abs(i_row - i_col)] += p(i_row, i_col) * d; 
    }
  }
  
  mu_x = 0.0;
  HX = 0.0;
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    mu_x += i_row1 * d * p_x[i_row1];
    if(p_x[i_row1]) HX -= p_x[i_row1] * std::log(p_x[i_row1]);
  }
  
  mu_y = 0.0;
  HY = 0.0;
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    mu_y += i_col1 * d * p_y[i_col1];
    if(p_y[i_col1]) HY -= p_y[i_col1] * std::log(p_y[i_col1]);
  }
  
  sig_x = 0.0;
  for(R_len_t i_row1 = 1; i_row1 <= N; i_row1++) {
    sig_x += std::pow(i_row1 * d - mu_x, 2.0) * p_x[i_row1]; // IDEAS H_Variance
  }
  
  sig_y = 0.0;
  for(R_len_t i_col1 = 1; i_col1 <= N; i_col1++) {
    sig_y += std::pow(i_col1 * d - mu_y, 2.0) * p_y[i_col1];
  }
  
  mu_xpy = 0.0;
  s_ent = 0.0;
  mu_xmy = 0.0;
  d_ent = 0.0;
  d_var = 0.0;
  s_var = 0.0;
  if(invariant) {
    for(double k = 2; k <= 2 * N; k++) {
      mu_xpy += 2 * (k - 1) * d_xpy * p_xpy[k];
      if(p_xpy[k]) s_ent -= p_xpy[k] * std::log(p_xpy[k]);
    }
    for(double k = 0; k < N; k++) {
      mu_xmy += (k + 1) * d * p_xmy[k];
      if(p_xmy[k]) d_ent -= p_xmy[k] * std::log(p_xmy[k]);
    }
    for(double k = 0; k < N; k++) d_var += std::pow((k + 1) * d - mu_xmy, 2.0) * p_xmy[k];
    for(double k = 2; k <= 2 * N; k++) s_var += std::pow(2 * (k - 1) * d_xpy - mu_xpy, 2.0) * p_xpy[k];
  } else {
    for(double k = 2; k <= 2 * N; k++) {
      mu_xpy += k * p_xpy[k];
      if(p_xpy[k]) s_ent -= p_xpy[k] * std::log(p_xpy[k]);
    }
    for(double k = 0; k < N; k++) {
      mu_xmy += k * p_xmy[k];
      if(p_xmy[k]) d_ent -= p_xmy[k] * std::log(p_xmy[k]);
    }
    for(double k = 0; k < N; k++) d_var += std::pow(k - mu_xmy, 2.0) * p_xmy(k);
    for(double k = 2; k <= 2 * N; k++) s_var += std::pow(k - mu_xpy, 2.0) * p_xpy(k);
  }
  
  HXY = 0.0;
  HXY1 = 0.0;
  HXY2 = 0.0;
  // Rcpp::NumericMatrix Q(N, N);
  // Q.fill(0.0);
  for(R_len_t i_col = 0; i_col < N; i_col++) {
    R_len_t i_col1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < N; i_row++) {
      R_len_t i_row1 = i_row + 1;
      if(p(i_row, i_col)) HXY -= p(i_row, i_col) * std::log(p(i_row, i_col));
      if(p_x[i_row1] * p_y[i_col1] != 0) {
        HXY1 -= p(i_row, i_col) * std::log(p_x[i_row1] * p_y[i_col1]);
        HXY2 -= p_x[i_row1] * p_y[i_col1] * std::log(p_x[i_row1] * p_y[i_col1]);
      }
      // for(k = 1; k <= N; k++) {
      //   Q(i_row, i_col) += p(i_row, k - 1) * p(i_col, k - 1) / (p_x(i_row) * p_y(k));
      // }
    }
  }
  imc1 = (HXY - HXY1) / std::max(HX, HY);
  imc2 = std::sqrt(1 - std::exp(-2*(HXY2 - HXY)));
  // mcc = std::sqrt(lambda(Q(i_col1, i_row1)))
  
  c_pro = 0.0;
  c_sha = 0.0;
  a_cor = 0.0;
  con = 0.0;
  cor = 0.0;
  dis = 0.0;
  nrj = 0.0;
  ent = 0.0;
  hom = 0.0;
  // hom2 = 0.0;
  i_dif = 0.0;
  s_sqr = 0.0;
  mu = (mu_x + mu_y) / 2;
  p_max = p(0, 0);
  
  for(R_len_t i_col = 0; i_col < N; i_col++) {
    R_len_t i_col1 = i_col + 1;
    for(R_len_t i_row = 0; i_row < N; i_row++) {
      R_len_t i_row1 = i_row + 1;
      // AUTOCORRELATION
      a_cor += (i_row1 * d * i_col1 * d) * p(i_row, i_col); 
      // CLUSTER PROMINENCE
      c_pro += std::pow(i_row1 * d + i_col1 * d - mu, 3.0) * p(i_row, i_col);
      // CLUSTER SHAPE
      c_sha += std::pow(i_row1 * d + i_col1 * d - mu, 4.0) * p(i_row, i_col);
      con += std::pow(i_row1 * d - i_col1 * d, 2.0) * p(i_row, i_col);
      // CORRELATION
      cor += (i_row1 * d - mu_x) * (i_col1 * d - mu_y) * p(i_row, i_col) / std::sqrt(sig_x * sig_y);
      // CONTRAST
      dis += std::abs(i_row1 * d - i_col1 * d) * p(i_row, i_col);
      // ANGULAR SECOND MOMENT
      nrj += std::pow(p(i_row, i_col), 2.0); 
      // ENTROPY
      if(p(i_row, i_col)) ent -= p(i_row, i_col) * std::log(p(i_row, i_col));
      // HOMOGENEITY (it is INVERSE DIFFERENCE MOMENT in Haralick paper)
      hom += p(i_row, i_col) / (1 + std::pow(i_row1 * d - i_col1 * d, 2.0));
      // H_Homogeneity from IDEAS
      // hom2 += (i_row1 * d + i_col1 * d) * p(i_row, i_col); 
      // INVERSE DIFFERENT MOMENT
      i_dif += p(i_row, i_col) / (1 + std::abs(i_row1 * d - i_col1 * d));
      // MAXIMUM PROBABILITY
      p_max = std::max(p(i_row, i_col), p_max);
      // SUM OF SQUARE
      s_sqr += (i_row1 - mu) * (i_row1 - mu) * p(i_row, i_col);
    }
  }
  
  return Rcpp::NumericVector::create(_["autocorrelation"] = a_cor,
                                     _["H Contrast"] = con / N / N, // IDEAS H_Contrast con is divided by probability since IDEAS claims contrast is [0,1]
                                     _["H Correlation"] = cor, // IDEAS H_Correlation;
                                     _["dissimilarity"] = dis,
                                     _["H Energy"] = nrj, // IDEAS H_Energy
                                     _["H Entropy"] = ent, // IDEAS H_Entropy
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
