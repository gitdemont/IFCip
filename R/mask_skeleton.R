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

#' @title Skeleton Mask
#' @description
#' The skeleton mask provides the barebone structure of the object from the starting mask
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param option whether "Thin" or "Thick". Default is "Thin". Not yet implemented "Thick".
#' @details See 'A new thinning algorithm for binary images' from L. Ben Boudaoud, A. Sider, A. Tari.'
#' 3rd international conference on control, engineering & information technology, May 2015. \doi{10.1109/CEIT.2015.7233099}.
#' @return a logical matrix corresponding to the skeleton area.
mask_skeleton <- function(msk, option = "Thin") {
  assert(msk, cla = "IFC_msk")
  # Thin:0, Thick:1
  assert(option, len = 1, alw = c("Thick", "Thin"))
  if(option == "Thin") {
    foo = cpp_thinning_bst(mat = msk)
  } else {
    stop("Skeleton with 'option' = \"Thick\" is not yet implemented")
  }
  ATT = list("IFC_msk", "skeleton", option, c(nrow(msk), ncol(msk)))
  names(ATT) = c("class", "type", "option", "dim")
  attributes(foo) <- ATT
  return(foo)
}
