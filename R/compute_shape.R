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

#' @title Shape Computation
#' @description
#' Computes shape of a mask
#' @param msk an object of class `IFC_msk`, the mask matrix.
#' @param mag magnification. Default is 40. Allowed are 20, 40 and 60.
#' @details Based on contour tracing labeling by Chang F. et al. \doi{10.1016/j.cviu.2003.09.002},
#' convex hull computation from Andrew's Monotone Chain Algorithm \doi{10.1016/0020-0190(79)90072-3}
#' and bounding box determination thanks to rotating calipers as described by Pirzadeh, Hormoz under supervision of Godfried T. Toussaint \url{https://escholarship.mcgill.ca/concern/theses/fx719p46g}
#' @return shape features from 'msk'.
#' @export
compute_shape <- function(msk, mag = 40) {
  if(!any(inherits(msk, what = "IFC_msk"))) stop("When provided 'msk' should be of class `IFC_msk`")
  out_names = c("Perimeter", "Diameter", "Circularity", "convexity", "roundness", 
                "Height", "Width", "Elongatedness", "convex perimeter", "convex cx", "convex cy") # from bbox
  ctl = cpp_ctl(msk, global = TRUE)
  # no object masked
  if(ctl$nb_lab == 0) return(structure(rep(NA, length(out_names)), names = out_names))
  
  mag = as.character(as.integer(mag)); assert(mag, len = 1, alw = c("20", "40", "60"))
  k = switch(mag, "20" = 1, "40" = 0.5, "60" = 0.3)
  hu = cpp_features_hu2(msk, switch(mag, "20" = 1, "40" = 4, "60" = 9))
  contours = ctl$contours
  contours = by(contours[, c(1,2,4,5), drop = FALSE], contours[, 3, drop = FALSE], FUN =function(d) by(d[,c(1,2,3), drop = FALSE], d[,4, drop = FALSE], FUN = function(dd) dd))
  contours = contours[as.integer(names(contours)) > 0] 
  contours = contours[[1]] # we only keep 1 object detected
  if(inherits(contours, what = "by")) contours = contours[[1]] # we remove internal contour if any

  perimeter = k * sum(ctl$perimeter)
  diameter = 2 * sqrt(hu["Area"] / pi)
  # center = apply(contours[,1:2], 2, mean)
  center = hu[c("pix cy", "pix cx")]
  distance = k * apply(contours[,1:2], 1, FUN =function(coord)  sqrt((coord[1] - center[1])^2 + (coord[2] - center[2])^2))
  radius = mean(distance)
  circularity = radius / sd(distance)
  bbox = cpp_bbox(cpp_convexhull(as.matrix(contours)), k)
  convexity = bbox["convex perimeter"] / perimeter
  roundness = 4 * pi * hu["Area"] / bbox["convex perimeter"]^2
  
  return(structure(c(perimeter, diameter, circularity, convexity, roundness, bbox), names = out_names))
}
