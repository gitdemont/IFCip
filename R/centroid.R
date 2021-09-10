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

#' @title Hu's Centroid
# This function is designed to compute Hu's image centroids.
#' @param img an object of class `IFC_img` or `IFC_msk`.
#' @details See Visual pattern recognition by moment invariants authored
#' by Ming-Kuei Hu \url{https://doi.org/10.1109/TIT.1962.1057692}
#' @return a named vector. where pix cx is x image's centroid and pix cy y image's centroid.
centroid <- function(img) {
  if(!any(inherits(img, what = c("IFC_img", "IFC_msk")))) stop("'img' be should of class `IFC_img` or `IFC_msk`")
  return(cpp_centroid(img))
}
