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

#' @title Haralick Texture
#' @description This function is designed to compute Haralick's Features.
#' @param img an object of class `IFC_img`, the image matrix.
#' @param msk an object of class `IFC_msk`, the mask matrix. If missing, the default mask information will be extracted from 'img'.
#' Otherwise, only elements superior to zero will be considered.
#' @param granularity a integer vector. Controls the grain of the texture. Default is 3. Allowed are [1-20].
#' For very fine textures, this value is small (1-3 pixels), while for very coarse textures, it is large (>10).
#' @param bits gray levels depth. Default is 4 (i.e. 2^4 = 16 gray levels). Allowed are [2-10].
#' Before applying compensation `IFC_img` have 2^12 = 4096 gray levels.
#' @details See Haralick's original paper \url{https://haralick.org/journals/TexturalFeatures.pdf} and Löfstedt T. et al. url{https://doi.org/10.1371/journal.pone.0212110} used for the computation:
#' @return array of Haralick's texture Features with dimensions corresponding to:\cr
#' -1st: the type of value returned,\cr
#' -2nd: the Haralick feature,\cr
#' -3rd: the granularity.
#' @export
compute_haralick = function(img, msk, granularity = 3, bits = 4) {
  if(!any(inherits(img, what = "IFC_img"))) stop("'img' should be of class `IFC_img`")
  if(!attr(img, "mode") == "raw") stop("'img' should have been extracted using 'mode' = \"raw\"")
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
  assert(img, cla = "IFC_img")
  assert(bits, len = 1, alw = 2:10)
  granularity = as.integer(granularity); granularity = granularity[granularity>0]; granularity = na.omit(granularity[is.finite(granularity)])
  assert(granularity, alw = c(1:20))
  
  # clip img to 0-4095
  # due to compensation raw img may exceed this range
  # TODO ask Amnis what they choose: normalization / clip...
  # otherwise use cpp_R_shift_M which generate error when img is outside [0,4095]
  rescaled = cpp_rescale_M(img, bits = bits)

  ans = lapply(granularity, FUN = function(i) {
    apply(sapply(cpp_cooc(img = rescaled, msk = msk, delta = i), cpp_h_features), 1, FUN = function(x) c(mean(x), sd(x)))
  })
  N = dimnames(ans[[1]])[[2]]
  return(array(unlist(ans), dim = c(2, length(N), length(granularity)), dimnames = list(value = c("Mean","Std"), features = N,  granularity = granularity)))
}
