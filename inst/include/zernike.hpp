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

#ifndef IFCIP_ZERNIKE_HPP
#define IFCIP_ZERNIKE_HPP

#include <Rcpp.h>
using namespace Rcpp;

// generates Z index from n, l
int nl_index(int n, int l) {
  int i, sum = 0;
  if ( n < 2 ) return n;
  for ( i = 1; i <= n / 2; i++ ) sum += i;
  if ( n % 2 == 0 ) return 2 * sum + l / 2 + 0;
  return 2 * sum + l / 2 + n / 2 + 1;
}

int p_nl_index(int n, int l) {
  return ((n * (n + 1)) / 2 + l);
}

double factorial(int x) {
  return std::floor( Rf_gammafn( (x + 1) ) );
}

R_len_t ZIND(uint8_t n, uint8_t l, uint8_t N1) {
  return l + n * N1;
}

R_len_t VIND(uint8_t n, uint8_t l, uint8_t m, uint8_t N1) {
  return l + n * N1 + m * N1 * N1;
}

//' @title Zernike's Features
//' @name cpp_zernike1
//' @description
//' This function is designed to compute Zernike's moments from image.
//' It will compute Zernike's moments but will not return image projection.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk_ a Nullable LogicalMatrix. Default is R_NilValue.
//' @param cx a double. X centroid. Default is 0.0.
//' @param cy a double. Y centroid. Default is 0.0.
//' @param zmax a uint8_t, maximal order of Zernike polynomials to be computed. Default is 15. Values outside [0,99] will be clipped.
//' Be aware that computation of Zernike's moments can be quite long when 'zmax' is high.
//' @param radius a numeric, radius of the circle in pixels around object centers from which the features are calculated. Default is 15.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_zernike1(const Rcpp::NumericMatrix img,
                        const Rcpp::Nullable<Rcpp::LogicalMatrix> msk_ = R_NilValue,
                        const double cx = 0.0,
                        const double cy = 0.0,
                        const uint8_t zmax = 15,
                        const double radius = 15.0) {
  R_len_t nx, ny, x, y;
  double d2, newx, newy, vnl, theta;
  char label[8] = "zn00m00";
  uint8_t Nmax = zmax > 99 ? 99:zmax;
  uint8_t l, m, n, i, msize, N1 = Nmax + 1;
  
  nx = img.nrow();
  ny = img.ncol();
  msize = nl_index(Nmax, Nmax) + 1; // inlcude n = 0
  Rcpp::NumericVector v(N1 * N1 * (Nmax/2 + 1));
  Rcpp::IntegerVector z(N1 * N1);
  
  Rcpp::CharacterVector nm(msize);
  for ( n = 0; n <= Nmax; n++ ) {
    for ( l = 0; l <= n; l++ ) {
      if ( (n - l) % 2 != 0 ) continue;
      i = nl_index(n, l);
      /* z(n,l) cached */
      z[ ZIND(n,l,N1) ] = i;
      label[2] = (char)(n / 10 + 48);
      label[3] = (char)(n - (n / 10) * 10 + 48);
      label[5] = (char)(l / 10 + 48);
      label[6] = (char)(l - (l / 10) * 10 + 48);
      nm(i) = label;
      for (m = 0; m <= (n - l) / 2; m++ ) {
        /* v(n,l,m) cached */
        if (n == l) {
          v[ VIND(n,l,m,N1) ] = (n + 1) / M_PI;
        } else {
          // n <- n, m <- l, s <- m,   
          v[ VIND(n,l,m,N1) ] = (n + 1) * ((m%2==0)?(1.0):(-1.0)) * factorial(n - m) / 
            ( M_PI * factorial(m) * factorial(0.5*(n - 2*m + l)) 
                * factorial(0.5*(n - 2*m - l)) );
        }
      }
    }
  }
  
  Rcpp::LogicalMatrix msk = get_mask(msk_, nx, ny);
  Rcpp::NumericVector zmoment(msize, 0.0);
  Rcpp::NumericVector mRe(msize, 0.0);
  Rcpp::NumericVector mIm(msize, 0.0);
  
  for ( x = 0; x < nx; x++ ) {
    newx = (x - cx) / radius;
    for ( y = 0; y < ny; y++ ) {
      newy = (y - cy) / radius;
      d2 = newx * newx + newy * newy;
      if ( !msk(x, y) || (d2 > 1.0) ) continue; // only use pixels within normalized unit circle
      theta = std::atan2(newy, newx);
      
      for ( n = 0; n <= Nmax; n++ ) {
        for ( l = 0; l <= n; l++ ) {
          if ( (n - l) % 2 != 0 ) continue; // only want values when (n - l) is even
          vnl = 0.0;
          for (m = 0; m  <= (n - l) / 2; m++ ) vnl += v[ VIND(n,l,m,N1) ] * pow(d2, 0.5 * n - m);
          vnl *= img(x, y);
          // no need to do e^i*l*theta if l=0 e^j*n*theta
          mRe[z[ZIND(n,l,N1)]] += vnl * std::cos(l * theta);
          mIm[z[ZIND(n,l,N1)]] += vnl * std::sin(l * theta);
        }
      }
    }
  }
  // loop through mRe, mIm, get abs and mult by (n+1)/pi, save in mRe, which is returned! 
  for ( i = 0; i < msize; i++ ) {
    zmoment[i] = std::sqrt(mRe[i]*mRe[i] + mIm[i]*mIm[i]);
  }
  zmoment.names() = nm;
  return List::create(_["zmoment"] = zmoment);   
}

//' @title Zernike's Features with Projections
//' @name cpp_zernike2
//' @description
//' This function is designed to compute Zernike's moments from image.
//' It will compute Zernike's moments but also return image projection.
//' @param img a NumericMatrix, containing image intensity values.
//' @param msk_ a Nullable LogicalMatrix. Default is R_NilValue.
//' @param cx a double. X centroid. Default is 0.0.
//' @param cy a double. Y centroid. Default is 0.0.
//' @param zmax a uint8_t, maximal order of Zernike polynomials to be computed. Default is 15. Values outside [0,99] will be clipped.
//' Be aware that computation of Zernike's moments can be quite long when 'zmax' is high.
//' @param radius a numeric, radius of the circle in pixels around object centers from which the features are calculated. Default is 15.
//' @source Adaptation from \url{https://github.com/aoles/EBImage} in v3.12.0, authored by Andrzej Oles, Gregoire Pau, Mike Smith, Oleg Sklyar, Wolfgang Huber, with contributions from Joseph Barry and Philip A. Marais \email{andrzej.oles@embl.de}.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_zernike2(const Rcpp::NumericMatrix img,
                        const Rcpp::Nullable<Rcpp::LogicalMatrix> msk_ = R_NilValue,
                        const double cx = 0.0, 
                        const double cy = 0.0, 
                        const uint8_t zmax = 15, 
                        const double radius = 15.0) {
  R_len_t nx, ny, x, y;
  double d2, newx, newy, vnl, theta;
  char label[8] = "zn00m00";
  uint8_t Nmax = zmax > 99 ? 99:zmax;
  uint8_t l, m, n, i, msize, N1 = Nmax + 1;
  
  nx = img.nrow();
  ny = img.ncol();
  msize = nl_index(Nmax, Nmax) + 1;
  Rcpp::NumericVector v(N1 * N1 * (Nmax/2 + 1));
  Rcpp::IntegerVector z(N1 * N1);
  
  Rcpp::CharacterVector nm(msize);
  for ( n = 0; n <= Nmax; n++ ) {
    for ( l = 0; l <= n; l++ ) {
      if ( (n - l) % 2 != 0 ) continue;
      i = nl_index(n, l);
      z[ ZIND(n,l,N1) ] = i;
      label[2] = (char)(n / 10 + 48);
      label[3] = (char)(n - (n / 10) * 10 + 48);
      label[5] = (char)(l / 10 + 48);
      label[6] = (char)(l - (l / 10) * 10 + 48);
      nm(i) = label;
      for (m = 0; m <= (n - l) / 2; m++ ) {
        if (n == l) {
          v[ VIND(n,l,m,N1) ] = (n + 1) / M_PI;
        } else {
          v[ VIND(n,l,m,N1) ] = (n + 1) * ((m%2==0)?(1.0):(-1.0)) * factorial(n - m) /
            ( M_PI * factorial(m) * factorial((n - 2 * m + l) / 2)  * factorial((n - 2 * m - l) / 2) );
        }
      }
    }
  }
  
  Rcpp::LogicalMatrix msk = get_mask(msk_, nx, ny);
  Rcpp::NumericVector zmoment(msize, 0.0);
  Rcpp::NumericVector mRe(msize, 0.0);
  Rcpp::NumericVector mIm(msize, 0.0);
  Rcpp::NumericVector nZm(msize * nx * ny, 0.0);
  Rcpp::NumericVector nZ_m(msize * nx * ny, 0.0);
  R_len_t idx;
  
  for ( x = 0; x < nx; x++ ) {
    newx = (x - cx) / radius;
    for ( y = 0; y < ny; y++ ) {
      newy = (y - cy) / radius;
      d2 = newx * newx + newy * newy;
      if ( !msk(x, y) || (d2 > 1.0) ) continue;
      theta = std::atan2(newy, newx);
      
      for ( n = 0; n <= Nmax; n++ ) {
        for ( l = 0; l <= n; l++ ) {
          if ( (n - l) % 2 != 0 ) continue;
          vnl = 0.0;
          for (m = 0; m  <= (n - l) / 2; m++ ) vnl += v[ VIND(n,l,m,N1) ] * pow(d2, 0.5 * n - m);
          vnl *= img(x,y);
          idx = z[ZIND(n,l,N1)] * nx * ny + y * nx + x;
          nZm[idx] = vnl * std::cos(l * theta);
          nZ_m[idx] = vnl * std::sin(l * theta);
          mRe[z[ZIND(n,l,N1)]] += vnl * std::cos(l * theta);
          mIm[z[ZIND(n,l,N1)]] += vnl * std::sin(l * theta);
        }
      }
    }
  }
  for ( i = 0; i < msize; i++ ) {
    zmoment[i] = std::sqrt(mRe[i]*mRe[i] + mIm[i]*mIm[i]);
  }
  zmoment.names() = nm;
  mRe.names() = nm;
  mIm.names() = nm;
  return Rcpp::List::create(_["zmoment"] = zmoment, 
                            _["Re"] = mRe,
                            _["Im"] = mIm,
                            _["even"] = nZm,
                            _["odd"] = nZ_m);
}

#endif
