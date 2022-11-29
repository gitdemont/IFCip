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

#' @title Zernike's Moments
#' @description Computation of Zernike moment features from image.
#' @param img an object of class `IFC_img` or `IFC_msk`.
#' @param msk an object of class `IFC_msk`, the mask matrix. If missing, the default mask information will be extracted from 'img'.
#' Otherwise, only elements superior to zero will be considered.
#' @param centroid a length 2 numeric vector of img centroids.
#' When missing, the default, it will be computed from img.
#' @param zmax an integer. Indicates the maximal order of Zernike polynomials
#' to be computed. Default value is 15. Values outside [0,99] will be clipped.
#' Be aware that computation of Zernike's moments can be quite long when 'zmax' is high.
#' @param radius a numeric. Defines the radius of the circle in pixels around
#' component center from which the features are calculated. Default is 30.
#' @param full a logical. Controls the quantity of information that will be returned. Default is FALSE.
#' @return a named list zernike moments values.\cr
#' if full == TRUE, then Re and Im parts are also returned + zernike matrix projection of the image.
#' @details Zernike features are computed by projecting image on Zernike complex polynomials.
#' @source adapted from File src/features_zernike.c in EBImage R package in its version 3.12.0 available at \url{https://www.bioconductor.org/packages/2.10/bioc/html/EBImage.html}
#' @export
moments_Zernike <- function(img, msk, centroid, zmax = 15, radius = 30, full = FALSE) {
  if(!any(inherits(img, what = c("IFC_img", "IFC_msk")))) stop("'img' be should of class `IFC_img` or `IFC_msk`")
  if(missing(msk)) {
    msk = attr(img, "mask")
    if(attr(msk, "removal") == "raw") {
      msk = (msk == 1)
    } else {
      msk = !msk
    }
  } else {
    if(!any(inherits(msk, what = "IFC_msk"))) stop("When provided 'msk' should be of class `IFC_msk`")
    msk = (msk > 0)
  }
  if(missing(centroid)) {
    centroid = cpp_centroid(msk)
  } else {
    centroid = na.omit(as.numeric(centroid));
    if(length(centroid) < 2) stop("when provided 'centroid' should be a numeric vector of at least length 2")
  }
  zmax = na.omit(as.integer(zmax)); zmax = zmax[(zmax>=0) && (zmax<=99)]; assert(zmax, len = 1, alw = 0:99)
  radius = na.omit(as.integer(radius)); assert(radius, len = 1, typ = "integer")
  full = na.omit(as.logical(full)); assert(full, len = 1, alw = c(TRUE, FALSE))
  if(full) {
    foo = cpp_zernike2(img = img, msk_ = msk, cx = centroid[1] - 1, cy = centroid[2] - 1, zmax = zmax, radius = radius)
    L = length(foo$zmoment)
    even = array(foo[["even"]], dim = c(dim(img), L))
    odd = array(foo[["odd"]], dim = c(dim(img), L))
    even = lapply(1:L, FUN=function(i) even[,,i] )
    names(even) = gsub("m", "+", names(foo$zmoment))
    odd = lapply(1:L, FUN=function(i) odd[,,i] )
    names(odd) = gsub("m", "-", names(foo$zmoment))
    foo[["even"]] = even
    foo[["odd"]] = odd
  }  else{
    foo = cpp_zernike1(img = img, msk_ = msk, cx = centroid[1] - 1, cy = centroid[2] - 1, zmax = zmax, radius = radius)
  }
  class(foo) = "IFCip_zernike"
  return(foo)
}
