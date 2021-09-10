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

#' @title Images Similarity
#' @description
#' This function is designed to score similarity between two images.
#' @param img1 a NumericMatrix, containing image values.
#' @param img2 a NumericMatrix, containing image values.
#' @param msk a LogicalMatrix, containing mask. If missing, the default mask information will be extracted from 'img1'.
#' Otherwise, only elements superior to zero will be considered.
#' @details the similarity is the log transformed Pearson's Correlation Coefficient.
#' It is a measure of the degree to which two images are linearly correlated within a masked region.
#' @return similarity score within 'msk'.
#' @export
compute_similarity <- function(img1, img2, msk) {
  if(!any(inherits(img1, what = "IFC_img"))) stop("'img1' should be of class `IFC_img`")
  if(!any(inherits(img2, what = "IFC_img"))) stop("'img2' should be of class `IFC_img`")
  if(missing(msk)) {
    msk = attr(img1, "mask")
    if(attr(msk, "removal") == "raw") {
      msk = (msk == 1)
    } else {
      msk = !msk
    }
  } else {
    if(!any(inherits(msk, what = "IFC_msk"))) stop("When provided 'msk' should be of class `IFC_msk`")
    msk = (msk > 0)
  }
  return(structure(cpp_similarity(img1 = img1, img2 = img2, msk = msk), names = "Similarity"))
} 
