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

# #' @name assert
# #' @keywords internal
assert <- getFromNamespace("assert", "IFC")

# #' @name whoami
# #' @keywords internal
whoami <- getFromNamespace("whoami", "IFC")

# #' @name cpp_getTAGS
# #' @keywords internal
cpp_getTAGS <- getFromNamespace("cpp_getTAGS", "IFC")

# #' @name num_to_string
# #' @keywords internal
num_to_string <- getFromNamespace("num_to_string", "IFC")

# #' @name ifcip_handler_winprogressbar
# #' @source modified from \link[progressr]{handler_winprogressbar} to allow to pass label argument
# #' @keywords internal
ifcip_handler_winprogressbar <- function(intrusiveness = getOption("progressr.intrusiveness.gui", 1), target = "gui", inputs = list(title = "sticky_message", label = "non_sticky_message"), ...) {
  backend_args <- list(...)
  not_fake <- TRUE
  if (not_fake) {
    if (.Platform$OS.type != "windows") {
      stop("handler_winprogressbar requires MS Windows: ",
           sQuote(.Platform$OS.type))
    }
    ## Import functions
    winProgressBar <- utils::winProgressBar
    setWinProgressBar <- utils::setWinProgressBar
  } else {
    winProgressBar <- function(...) rawConnection(raw(0L))
    setWinProgressBar <- function(...) NULL
  }
  alw_ini <- methods::formalArgs(winProgressBar)
  alw_upd <- methods::formalArgs(setWinProgressBar)
  
  reporter <- local({
    pb <- NULL
    
    stopifnot(
      is.list(inputs),
      !is.null(names(inputs)),
      all(names(inputs) %in% c("title", "label")),
      all(vapply(inputs, FUN = function(x) {
        if (is.null(x)) return(TRUE)
        if (!is.character(x)) return(FALSE)
        x %in% c("message", "non_sticky_message", "sticky_message")
      }, FUN.VALUE = FALSE))
    )
    
    ## Expand 'message' => c("non_sticky_message", "sticky_message")
    for (name in names(inputs)) {
      input <- inputs[[name]]
      if ("message" %in% input) {
        input <- setdiff(input, "message")
        input <- c(input, "non_sticky_message", "sticky_message")
      }
      inputs[[name]] <- unique(input)
    }
    
    ## Update winProgressBar
    set_pb <- function(state, progression) {
      args <- list()
      for (target in c("title", "label")) {
        if (inherits(progression, "sticky")) {
          if ("sticky_message" %in% inputs[[target]])
            args[[target]] <- progression$message
        } else {
          if ("non_sticky_message" %in% inputs[[target]])
            args[[target]] <- progression$message
        }
      }
      for (target in c("title", "label")) if (length(args[[target]]) == 0) args[[target]] <- ifelse(length(pb$args[[target]]) == 0, " ", pb$args[[target]])
      args <- c(list(pb = pb$bar, value = state$step), args)
      if(not_fake) args <- args[names(args) %in% alw_upd]
      pb$args <<- args[setdiff(names(args), "pb")]
      do.call(what = setWinProgressBar, args = args)
    }
    
    list(
      reset = function(...) {
        pb <<- NULL
      },
      
      initiate = function(config, state, progression, ...) {
        if (!state$enabled || config$times == 1L) return()
        ## NOTE: 'pb' may be re-used for winProgressBar:s
        if (config$clear) stopifnot(is.null(pb))
        args <- c(backend_args, list(max = config$max_steps, initial = state$step), list(...))
        if(not_fake) args <- args[names(args) %in% alw_ini]
        # Empty title or label are replaced to avoid error while creating the progress bar
        # In addition, if the progress bar has been created with default label="" value label,
        # it won't be possible to modify with setProgressBar afterwards,
        # so as a trick label value is replaced with " " when NULL or equal to ""
        if (length(args$title) == 0) args$title = " "
        if (length(args$label) == 0 || args$label == "") args$label = " "
        pb <<- c(list(bar = do.call(winProgressBar, args = args)), list(args = args))
        pb
      },
      
      update = function(config, state, progression, ...) {
        if (!state$enabled || config$times <= 2L) return()
        set_pb(state, progression)
      },
      
      finish = function(config, state, progression, ...) {
        ## Already finished?
        if (is.null(pb)) return()
        if (!state$enabled) return()
        if (config$clear) {
          close(pb$bar)
          pb <<- NULL
        } else {
          set_pb(state, progression)
        }
      }
    )
  })
  
  progressr::make_progression_handler("winprogressbar", reporter, intrusiveness = intrusiveness, target = target, ...)
}
