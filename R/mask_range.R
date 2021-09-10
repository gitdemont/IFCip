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

#' @title Range Mask
#' @description
#' The range mask identifies region of desired area and aspect ration.
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param area a length 2 numeric vector of area range defining the mask. Default is c(0,5000).
#' @param aspect_ratio a length 2 numeric vector of aspect ratio defining the mask. Default is c(0,1).
#' @return a logical matrix of area(s) respecting range parameters found.
mask_range <- function(msk, area = c(0,5000), aspect_ratio = c(0,1)) {
  assert(msk, cla = "IFC_msk")
  assert(area, len=2)
  assert(aspect_ratio, len=2)
  components = cpp_ctl(msk)
  nC = setdiff(components[["nb_lab"]],0)
  if(length(nC) == 0) {
    foo = matrix(FALSE, ncol = ncol(msk), nrow = nrow(msk))
  } else {
    components = lapply(1:nC, FUN = function(lab) { components[["matrix"]] == lab })
    bar = sapply(components, FUN=function(i_comp) { 
      moments = moments_Hu(i_comp, full = FALSE)
      aspect = moments["Minor Axis"] / moments["Major Axis"]
      return((moments["Area"] >= area[1]) & (moments["Area"] <= area[2]) & (aspect >= aspect_ratio[1]) & (aspect<= aspect_ratio[2]))
    })
    if(any(bar)) {
      foo = cpp_OR_M(components[bar])
    } else {
      foo = matrix(FALSE, ncol = ncol(msk), nrow = nrow(msk))
    }
  }
  ATT = list("IFC_msk", "range", area, aspect_ratio, c(nrow(msk), ncol(msk)))
  names(ATT) = c("class", "type", "area", "aspect_ratio", "dim")
  attributes(foo) <- ATT
  return(foo)
}
