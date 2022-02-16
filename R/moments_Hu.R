################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
#                                                                              #
# It is part of IFCip package, please cite:                                    #
# -IFCip: An R Package for Imaging Flow Cytometry Image Processing             #
# -YEAR: 2021                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Jean-Pierre Marolleau, Loïc Garçon,       #
#                     CHU Amiens                                               #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFCip. If not, see <http://www.gnu.org/licenses/>.                #
################################################################################

#' @title Hu's Moments
#' @description This function is designed to compute Hu's central + invariant moments.
#' @param img an object of class `IFC_img` or `IFC_msk`.
#' @param mag magnification. Default is 40. Allowed are 20, 40 and 60.
#' @param full a logical. Controls the quantity of information that will be returned. Default is TRUE.
#' if full == TRUE, the seven first invariant moments will also be returned.
#' @details returned moments are based on Visual pattern recognition by moment invariants authored
#' by Ming-Kuei Hu \doi{10.1109/TIT.1962.1057692}
#' circularity is computed according to A Hu moment invariant as a shape circularity
#' measure authored by Joviša Žunić, Kaoru Hirota, Paul L. Rosin \doi{10.1016/j.patcog.2009.06.017}
#' @return a named vector of Hu's moments values.\cr
#' -Area: img's area\cr
#' -circularity: img's circularity\cr
#' -Minor Axis: img's ellipsis minor axis\cr
#' -Major Axis: img's ellipsis major axis\cr
#' -Aspect Ratio: img's ratio of minor_axis over major_axis\cr
#' -Angle: img's ellipsis angle with x axis (in radians)\cr
#' -theta: img's ellipsis theta angle (in radians)\cr
#' -eccentricity: img's ellipsis ecentricity\cr
#' -pix cx: img's pixel x centroïd\cr
#' -pix cy: img's pixel y centroïd\cr
#' -pix min axis: img's ellipsis minor axis in pixels\cr
#' -pix maj axis: img's ellipsis major axis in pixels\cr
#' -pix count: img's area in pixels\cr
#' -inv[1-7]: image invariant moments, depending on 'full' parameter.
#' @export
moments_Hu <- function(img, mag = 40, full = TRUE) {
  if(!any(inherits(img, what = c("IFC_img", "IFC_msk")))) stop("'img' should be of class `IFC_img` or `IFC_msk`")
  mag = as.character(as.integer(mag)); assert(mag, len = 1, alw = c("20", "40", "60"))
  mag = switch(mag, "20" = 1, "40" = 4, "60" = 9)
  full = na.omit(as.logical(full)); assert(full, len = 1, alw = c(TRUE, FALSE))
  if(full) {
    foo = cpp_features_hu2(img = img, mag = mag)
  } else {
    foo = cpp_features_hu1(img = img, mag = mag)
  }
  class(foo) = "IFCip_hu"
  return(foo)
}
