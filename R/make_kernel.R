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

#' @title Create a Kernel
#' @description
#' Create a kernel of the desired 'type'. 
#' @param size the size of the desired kernel. Default is \code{3}.\cr
#' It does not apply to 'type' "sobelx", "sobely", "scharr", "n4" nor "n8" which are of 'size' \code{3}.
#' @param type a string indicating the type of the kernel Default is "box". Allowed are: 
#' "box", "cross", "plus", "disc", "diamond", "mean", "gaussian", "laplacian", "sobelx", "sobely", "scharrx",  "scharry", "n4", "n8".
#' @param sigma standard deviation of the "gaussian" or "laplacian" \code{'type'}. Default is \code{-0.3}. If negative, \code{'sigma'} will be determined using \code{-1.0 * sigma * ((size - 1) * 0.5)}.
#' @details "box", "cross", "plus", "disc", "diamond" will result in a logical matrix whereas the other will produce numeric a matrix.
#' @source 'gaussian' is from EBImage::makeBrush R package and 'Laplacian' adapted from 'gaussian'
#' @return a kernel matrix.
make_kernel <- function(size = 3, type = "box", sigma = -0.3) {
  size = na.omit(as.integer(size)); assert(size, len = 1)
  type = tolower(type);assert(type, len = 1, alw = c("box", "cross", "plus", "disc", "diamond", 
                                                     "mean", "gaussian", "laplacian",
                                                     "sobelx", "sobely", "scharrx", "scharry",
                                                     "n4","n8"))
  return(switch(type, 
                "box" = cpp_make_box(size),
                "cross" =  cpp_make_cross(size),
                "plus" = cpp_make_plus(size),
                "disc" = cpp_make_disc(size),
                "diamond" = cpp_make_diamond(size),
                "mean" = matrix(1/(size*size), ncol = size, nrow = size),
                "gaussian" = cpp_make_gaussian(size, sigma),
                "laplacian" = cpp_make_laplacian(size, sigma),
                "sobelx" = matrix(c( -1,  0,  1,
                                     -2,  0,  2,
                                     -1,  0,  1), ncol = 3, byrow = TRUE),
                
                "sobely" = matrix(c( -1, -2, -1,
                                     0,  0,  0,
                                     1,  2,  1), ncol = 3, byrow = TRUE),
                
                "scharrx" = matrix(c( -3,  0,  3,
                                      -10,  0, 10,
                                      -3,  0,  3), ncol = 3, byrow = TRUE),
                
                "scharry" = matrix(c(  -3, -10, -3,
                                       0,   0,  0,
                                       3,  10,  3), ncol = 3, byrow = TRUE),
                
                "n4" = matrix(c(  0, -1,  0,
                                  -1,  4, -1,
                                  0, -1,  0), ncol = 3, byrow = TRUE),
                
                "n8" = matrix(c( -1, -1, -1,
                                 -1,  8, -1,
                                 -1, -1, -1), ncol = 3, byrow = TRUE)
  ))
}
