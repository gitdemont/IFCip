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

#' @title Threshold Mask
#' @description
#' The threshold mask identifies intensity areas from an image that are superior or equal to max(img) - threshold * diff(range(img)) / 100 within msk.
#' @param img an object of class `IFC_img`, the image matrix.
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param threshold constant to be checked. Default is 0.
#' @param removal uint8_t, object removal method. Default is 0 for no removal. Otherwise, if\cr
#' -1, for clipped removal, keep non clipped foreground.\cr
#' -2, height clipped removal.\cr
#' -3, width clipped removal.\cr
#' -4, only keep background:.\cr
#' -5, only keep foreground.
#' @return a logical matrix of area(s) above threshold found.
mask_threshold <- function(msk, img, threshold = 0, removal = 0) {
  assert(msk, cla = "IFC_msk")
  assert(img, cla = "IFC_img")
  threshold = na.omit(as.integer(threshold)); threshol=threshold[(threshold>=0) & (threshold<=100)]
  assert(threshold, len=1, typ="integer")
  
  # TODO ask amnis what they used to detemine img range (within msk or not)
  foo = cpp_threshold(img = img, msk = msk, k = threshold, removal = removal)
  
  ATT = list("IFC_msk", "threshold", threshold, c(nrow(msk), ncol(msk)))
  names(ATT) = c("class", "type", "area", "dim")
  attributes(foo) <- ATT
  return(foo)
}
