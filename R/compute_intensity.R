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

#' @title Intensity Computation
#' @description
#' Computes intensity of image
#' @param img an object of class `IFC_img`, the image matrix.
#' @param msk an object of class `IFC_msk`, the mask matrix. If missing, the default mask information will be extracted from 'img'.
#' Otherwise, only elements superior to zero will be considered.
#' @return intensity wihtin 'msk'.
#' @export
compute_intensity <- function(img, msk) {
  if(!any(inherits(img, what = "IFC_img"))) stop("'img' should be of class `IFC_img`")
  if(missing(msk)) msk = attr(img, "mask")
  if(!any(inherits(msk, what = "IFC_msk"))) stop("'msk' should be of class `IFC_msk`")
  if(length(attr(msk, "removal")) != 0) {
    if(identical(attr(msk, "removal"), "raw")) {
      msk = (msk == 1)
    } else {
      msk = !msk
    }
  } else {
    msk = msk != 0
  }
  observed = sum(img * (msk == 1))
  expected = sum(msk == 1) * attr(img, "BG_MEAN")
  return(structure(observed - expected, names = "Intensity"))
} 
