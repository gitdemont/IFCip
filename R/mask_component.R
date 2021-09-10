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

#' @title Component Mask
#' @description
#' Apply component masking on an image. 
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param img an object of class `IFC_img`, the image matrix.\cr
#' Can be missing when feature is "Area", "Aspect Ratio", "Circularity", "Thickness Max" or "Thickness Min".
#' @param order Either "Ascending" or "Descending". Default is "Descending".
#' @param feature Ranking feature. Default is "Area".\cr
#' Allowed are: "Area", "Aspect Ratio", "Bright Detail Intensity R3", "Bright Detail Intensity R7", "Bright Detail Similarity R3", "Circularity",
#' "Contrast", "Gradient RMS", "Intensity", "Max Pixel", "Mean Pixel", "Median Pixel", "Min Pixel", "Saturation Count", "Thickness Max", "Thickness Min".\cr
#' Only 'feature' = "Area" is possible for the moment
#' @return list of mask(s) sorted according to 'order' and 'feature' selected.\cr
#' If no component is found (i.e. 'msk' is filled with 0), 'msk' is returned.
mask_component <- function(msk, img, order="Descending", feature="Area") {
  assert(msk, cla = "IFC_msk")
  assert(feature, len=1, alw=c("Area", "Aspect Ratio", "Bright Detail Intensity R3", "Bright Detail Intensity R7", "Bright Detail Similarity R3", "Circularity",
                               "Contrast", "Gradient RMS", "Intensity", "Max Pixel", "Mean Pixel", "Median Pixel", "Min Pixel", "Saturation Count", "Thickness Max", "Thickness Min"))
  if(feature%in%c("Bright Detail Intensity R3", "Bright Detail Intensity R7", "Bright Detail Similarity R3", "Circularity",
                  "Contrast", "Gradient RMS", "Intensity", "Max Pixel", "Mean Pixel", "Median Pixel", "Min Pixel", "Saturation Count")) {
    assert(img, cla = "IFC_img")
  }
  assert(order, len=1, alw=c("Descending", "Ascending"))
  components = cpp_ctl(msk)
  nC = setdiff(components[["nb_lab"]],0)
  foo = lapply(1:nC, FUN = function(lab) components[["matrix"]] == lab)
  if(length(foo)>0) {
    # TODO add cpp_features_hu3 to be able to sort by other than "Area"
    if(feature=="Area") {
      val = unlist(lapply(foo, sum))
    } else {
      stop("only 'feature'=='Area' is available for the moment")
    }
    O = order(val)
  } else {
    warning("no component found")
    O = 1
    foo = list(msk)
    val = 0
  }
  if(order=="Descending") O = rev(O)
  foo = foo[O]
  ATT = list("IFC_msk", "component", order, val[O])
  names(ATT) = c("class", "type", "order", feature)
  attributes(foo) <- ATT
  return(foo)
}
