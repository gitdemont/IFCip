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
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

# TODO find better method to identify objects

#' @title Identify Mask Version 1
#' @description
#' The identify mask identifies objects within image using algorithm 1.
#' @param img the image matrix.
#' @param threshold quantile above which object(s) should be kept.
#' @param size the size of the kernel used to identify object(s). Default is 5. Should be higher than 3.
#' @return an integer matrix of object(s) found.
#' @keywords internal
mask_identify1 <- function(img, threshold = 0.95, size = 5) {
  # blur image
  sigma = sqrt(size) / 2
  kg = make_kernel(size, type = "gaussian", sigma = sigma)
  img_blur = cpp_convolve2d(img, kernel = kg/sum(kg))

  # apply laplacian
  kl = make_kernel(size = size + 2, type = "laplacian", sigma = sigma)
  LOG = cpp_convolve2d(img_blur, kernel = kl)

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
  b = cpp_dilate(a, kernel = make_kernel(size = size + 2, type = "disc"))
  ctl = cpp_ctl(b)
  aa = cpp_fill(ctl, inner = TRUE, outer = TRUE)

  # erode and identify object(s)
  bb = cpp_erode(cpp_erode(aa, kernel = make_kernel(size = size, type = "box")), kernel = make_kernel(size = size, type = "box"))
  ctl = cpp_ctl(bb)
  foo = cpp_fill(ctl, inner = TRUE, outer = TRUE)
  
  ATT = list("IFC_msk", "identify", threshold, size, ctl$perimeter, c(nrow(foo), ncol(foo)))
  names(ATT) = c("class", "type", "threshold", "size", "perimeter", "dim")
  attributes(foo) <- ATT
  return(foo)
}

#' @title Identify Mask Version 2
#' @description
#' The identify mask identifies objects within image using algorithm 2.
#' @param img the image matrix.
#' @param threshold standard deviation above which object(s) should be kept. Default is 0.
#' @param size the size of the kernel used to identify object(s). Default is 5. Should be higher than 3.
#' @return an integer matrix of object(s) found.
#' @keywords internal
mask_identify2 <- function(img, threshold = 0, size = 5) {
  k = make_kernel(size, "box")
  ctl = cpp_ctl(cpp_sd(img, k) > threshold)
  foo = cpp_fill_out(ctl)
  
  ATT = list("IFC_msk", "identify", threshold, size, ctl$perimeter, c(nrow(foo), ncol(foo)))
  names(ATT) = c("class", "type", "threshold", "size", "perimeter","dim")
  attributes(foo) <- ATT
  return(foo)
}

#' @title Identify Mask
#' @description
#' The identify mask identifies objects within image.
#' @param img the image matrix.
#' @param threshold threshold above which object(s) should be kept. Default is 0.
#' @param size the size of the kernel used to identify object(s). Default is 5. Should be higher than 3.
#' @param version version to be used. Default is 2. Allowed are 1 and 2.
#' @details With version 1 object identification is done by percentage thresholding on Laplacian of Gaussian after Gaussian denoising.\cr
#' Version 2 will perform thresholding on local standard deviation.
#' @return an integer matrix of object(s) found.
#' @keywords internal
mask_identify <- function(img, threshold = 0, size = 5, version = 2) {
  version = as.integer(version); assert(version, len = 1, alw = c(1,2))
  size = as.integer(size); size = na.omit(size[size > 3]); assert(size, len=1)
  if(version == 1) {
    threshold = na.omit(threshold[threshold >= 0 & threshold <= 1]);assert(threshold, len=1) 
    return(mask_identify1(img = img, threshold = threshold, size = size))
  }
  threshold = na.omit(threshold[threshold > 0]);assert(threshold, len=1)
  return(mask_identify2(img = img, threshold = threshold, size = size))
}
