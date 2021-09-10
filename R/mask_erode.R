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

#' @title Erode Mask
#' @description
#' The erode mask removes the selected number of pixels from all edges of the starting mask.
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param n a positive integer corresponding to number of pixels to erode. Default is 1.
#' @return a logical matrix corresponding to the eroded area.
mask_erode <- function(msk, n = 1) {
  assert(msk, cla = "IFC_msk")
  n = as.integer(n); n = n[n>=0]; n = na.omit(n[is.finite(n)])
  assert(n, len=1)
  # TODO check / ask the kernel type used, the number of iteration
  foo = cpp_erode(mat = msk, kernel = make_kernel(size = n, type = "box"), iter = 0)
  ATT = list("IFC_msk", "erode", n, c(nrow(msk), ncol(msk)))
  names(ATT) = c("class", "type", "n", "dim")
  attributes(foo) <- ATT
  return(foo)
}
