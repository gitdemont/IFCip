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

#' @title Convert Features
#' @name as_IFC_features
#' @description
#' Function to convert an object of class `IFCip_features` to a list of `IFC_features` and `IFC_features_def`
#' @param x an object of class `IFCip_features` 
#' @return a list containing 2 members named:
#' -features, data.frame, object of class `IFC_features` 
#' -features_def, a list, object of class `IFC_features_def`
#' @export
as_IFC_features <- function(x) {
  assert(x, cla = "IFCip_features")
  ans = lapply(1:(dim(x)[3]), FUN = function(chan) {
    if(attr(x, "removal")[chan] == "MC") {
      mask = "MC"
    } else {
      mask = sprintf("M%02i", as.integer(attr(x, "channel_id")[chan]))
    }
    image = attr(x, "channel_names")[chan]
    foo = lapply(dimnames(x)[[2]], FUN = function(feat) {
      name = character()
      # check if 1st letter of feature's name is major case
      # in such case it can also be computed by IDEAS so we can create the def
      if(grepl("^[A-Z]", feat, ignore.case = FALSE)) {
        # special for H haralick feature we need to reorder to retrieve granularity
        if(grepl("^H ", feat, ignore.case = FALSE)) {
          name = sprintf(gsub("^(.*) (\\d.*)$", "\\1|%s|%s|\\2", feat), mask, image)
        } else {
          # feature is based on image
          if(grepl("intensity|pixel|modulation|contrast|gradient|spot|std|bkgd|raw|saturation|bright|internalization|similarity|xcorr", feat, ignore.case = TRUE)) {
            # feature is image only
            if(grepl("^Bkgd", feat, ignore.case = FALSE)) {
              name = sprintf("%s|%s", feat, image)
            } else {
              name = sprintf("%s|%s|%s", feat, mask, image)
            }
          } else {
            name = sprintf("%s|%s", feat, mask)
          }
        }
      }
      if(length(name) == 0) {
        buildFeature(name = sprintf("%s_%s_%s", feat, mask, image),
                     val = x[,feat,chan])
      } else {
        buildFeature(name = gsub("|", "_", name, fixed = T), type = "single", def = paste0(name, ifelse(grepl("H ", feat, ignore.case = FALSE), "|Granularity:|1|20", "")),
                     val = x[,feat,chan])
      }
    })
  })
  if(length(ans) > 1) {
    ans = do.call(what = "c", args = ans)
  }
  # remove duplicated, which came from using MC mask
  ans = ans[!duplicated(sapply(ans, simplify = TRUE, USE.NAMES = FALSE, FUN = function(i) i$name))]
  return(list(features = structure(data.frame("Object Number" = as.numeric(dimnames(x)[[1]]), 
                                              do.call(what = cbind, args = lapply(ans, FUN = function(i) i$val)),
                                              stringsAsFactors = FALSE, check.names = FALSE),
                                   names = c("Object Number", sapply(ans, simplify = TRUE, USE.NAMES = FALSE, FUN = function(i) i$name)),
                                   class = c("data.frame", "IFC_features")),
              features_def = structure(ans, class = c("list", "IFC_features_def"),
                                       names = sapply(ans, simplify = TRUE, USE.NAMES = FALSE, FUN = function(i) i$name))))
             
}
