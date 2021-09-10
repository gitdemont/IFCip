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

#' @title Perimeter Computation
#' @description
#' Computes perimeter of a mask
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param mag magnification. Default is 40. Allowed are 20, 40 and 60.
#' @param global whether to compute the perimeter globally or to evaluate the perimeter of each non 8-connected objects. Default is TRUE.
#' @details when global is TRUE the value returned is the total perimeter of all objects found in the mask.
#' If pixels from 2 different objects overlap, they are counted only once.
#' However, if FALSE, value returned is the number of perimeter of each single objects.
#' In case out border of 2 different objects overlap they are counted for each objects.
#' @return a data.frame of 'msk' perimeter.
#' @export
compute_perimeter <- function(msk, mag = 40, global = TRUE) {
  if(!any(inherits(msk, what = "IFC_msk"))) stop("When provided 'msk' should be of class `IFC_msk`")
  global = as.logical(global); assert(global, len = 1, alw = c(TRUE, FALSE))
  mag = as.character(as.integer(mag)); assert(mag, len = 1, alw = c("20", "40", "60"))
  mag = switch(mag, "20" = 1, "40" = 2, "60" = 3)
  if(global) return(data.frame(perimeter = sum(cpp_ctl(msk, global = TRUE)$perimeter)/mag))
  return(data.frame(perimeter = cpp_ctl(msk, global = TRUE)$perimeter)/mag)
}
