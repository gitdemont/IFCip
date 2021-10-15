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

################################################################################
#              function described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Identify Mask
#' @description
#' The identify mask identifies objects within image.
#' @param img the image matrix.
#' @param threshold quantile above which object(s) should be kept.
#' @param size the size of the kernel used to identify object(s). Default is 5. Should be higher than 3.
#' @return an integer matrix of object(s) found.
#' @keywords internal
mask_identify <- function(img, threshold = 0.95, size = 5) {
  # TODO find better method to identify objects
  threshold = na.omit(threshold[threshold >= 0 & threshold <= 1]);assert(threshold, len=1)
  size = as.integer(size); size = na.omit(size[size > 3]); assert(size, len=1)
  
  # blur image
  sigma = sqrt(size) / 2
  kg = make_kernel(size, type = "gaussian", sigma = sigma)
  img_blur = IFCip:::cpp_convolve2d(img, kernel = kg/sum(kg))
  
  # apply laplacian
  kl = make_kernel(size = size + 2, type = "laplacian", sigma = sigma)
  LOG = IFCip:::cpp_convolve2d(img_blur, kernel = kl)
  
  # compute background
  # back = cpp_background(img)
  # int_tresh = back["BG_MEAN"] + 0.5 * back["BG_STD"]
  
  # apply threshold
  threshold = sort(c(threshold, 1 - threshold))
  thresh = quantile(LOG, threshold)
  # a = ((LOG > thresh[2]) & (img < int_tresh)) |
      # ((LOG < thresh[1]) & (img > int_tresh))
  a = ((LOG > thresh[2]) | (LOG < thresh[1]))
  
  # dilate and fill object(s) found
  b = IFCip:::cpp_dilate(a, kernel = make_kernel(size = size + 2, type = "disc"))
  ctl = IFCip:::cpp_ctl(b)
  aa = IFCip:::cpp_fill(ctl, inner = TRUE, outer = TRUE)

  # erode and identify object(s)
  bb = IFCip:::cpp_erode(aa, kernel = make_kernel(size = size, type = "box"), iter = 1)
  ctl = IFCip:::cpp_ctl(bb)
  foo = IFCip:::cpp_fill(ctl, inner = TRUE, outer = TRUE)
  
  ATT = list("IFC_msk", "identify", threshold, size, c(nrow(foo), ncol(foo)))
  names(ATT) = c("class", "type", "threshold", "size", "dim")
  attributes(foo) <- ATT
  return(foo)
}
