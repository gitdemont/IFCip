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

#' @title Clipped Mask
#' @description
#' This function is aimed to detect all components that are bound to edges.
#' Except if it is principal component
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param borderLR number of pixels considered as Left/Right edge.
#' @param borderTB number of pixels considered as Top/Bottom edge.
#' @return a logical matrix of clipped components found.
#' @export
mask_clipped <- function(msk, borderLR = 1, borderTB = 6) {
  borderLR = na.omit(as.integer(borderLR)); assert(borderLR, len=1, typ="integer")
  borderTB = na.omit(as.integer(borderTB)); assert(borderTB, len=1, typ="integer")
  mat_r = nrow(msk);
  mat_c = ncol(msk);
  if(borderLR >= mat_c) stop("mask_clipped: 'borderLR' can't be superior than image width");
  if(borderTB >= mat_r) stop("mask_clipped: 'borderTB' can't be superior than image height");
  if(borderLR < 1) stop("mask_clipped: 'borderLR' can't be less than 1px");
  if(borderTB < 1) stop("mask_clipped: 'borderTB' can't be less than 1px");

  ALL_M = suppressWarnings(mask_component(msk));

  is_on_edge = sapply(ALL_M, FUN=function(CUR_M) {
    is_on_top = any(CUR_M[1:borderTB,])
    is_on_bottom = any(CUR_M[(mat_r - borderTB):mat_r, ])
    is_on_left = any(CUR_M[,1:borderLR])
    is_on_right = any(CUR_M[,(mat_c - borderLR):mat_c])
    return(any(is_on_top, is_on_bottom, is_on_left, is_on_right))
  })
  
  NONE_mask = matrix(FALSE, nrow = mat_r, ncol = mat_c);

  Idx = which(is_on_edge);
  if(length(Idx) == 0) return(NONE_mask);

  CLIPPED_mask = cpp_OR_M(ALL_M[Idx]);
  N_CLIPPED_mask = cpp_NEG_M(CLIPPED_mask);
  val = sum(cpp_AND_M(list(msk, N_CLIPPED_mask)));

  if(val == 0) return(NONE_mask);
  return(CLIPPED_mask);
}
