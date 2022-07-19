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

#' @name assert
#' @keywords internal
assert <- getFromNamespace("assert", "IFC")

#' @name whoami
#' @keywords internal
whoami <- getFromNamespace("whoami", "IFC")

#' @name cpp_getTAGS
#' @keywords internal
cpp_getTAGS <- getFromNamespace("cpp_getTAGS", "IFC")

#' @name num_to_string
#' @keywords internal
num_to_string <- getFromNamespace("num_to_string", "IFC")

#' @name ifcip_handler_winprogressbar
#' @source modified from \link[progressr]{handler_winprogressbar} to allow to pass label argument
#' @keywords internal
ifcip_handler_winprogressbar <- function (intrusiveness = getOption("progressr.intrusiveness.gui", 1), target = "gui", ...) {
  backend_args <- list(...)
  winProgressBar <- utils::winProgressBar
  setWinProgressBar <- utils::setWinProgressBar
  reporter <- local({
    pb <- NULL
    make_pb <- function(...) {
      if (!is.null(pb)) return(pb)
      args <- c(backend_args, list(...))
      alw_names <- methods::formalArgs(winProgressBar)
      args <- args[names(args) %in% alw_names]
      args <- args[!duplicated(names(args))]
      #### trick
      # Empty title or label are replaced to avoid error while creating the progress bar
      # In addition, if the progress bar has been created with default label="" value label,
      # it won't be possible to modify with setProgressBar afterwards,
      # so as a trick label value is replaced with " " when NULL or equal to ""
      if(length(args$title) == 0) args$title = " "
      if((length(args$label) == 0) || (args$label == "")) args$label = " "
      pb <<- do.call(winProgressBar, args = args)
      pb
    }
    list(reset = function(...) { pb <<- NULL },
         initiate = function(config, state, progression, ...) {
           if (!state$enabled || config$times == 1L) return()
           make_pb(max = config$max_steps, label = state$message, ...) },
         update = function(config, state, progression, ...) {
           if (!state$enabled || progression$amount == 0 || config$times <= 2L) return()
           make_pb(max = config$max_steps, label = state$message, ...)
           setWinProgressBar(pb, value = state$step, label = paste0(state$message, "")) },
         finish = function(config, state, progression, ...) {
           if (is.null(pb)) return()
           if (!state$enabled) return()
           if (config$clear) {
             close(pb)
             pb <<- NULL
           } else {
             setWinProgressBar(pb, value = state$step, label = paste0(state$message, ""))
           }
         })
  })
  progressr::make_progression_handler("winprogressbar", reporter, intrusiveness = intrusiveness, target = target, ...)
}