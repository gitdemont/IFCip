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

#' @title Features Extraction
#' @name ExtractFeatures
#' @description
#' Function to extract features from objects stored within rif and cif files.
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
#' If 'param' is provided 'export'(="matrix"), 'mode'(="raw"), 'size'(="c(0,0)"), 'force_width'(="FALSE") and 'removal' will be overwritten.\cr
#' If 'offsets' are not provided extra arguments can also be passed with ... to \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @param objects integers, indices of objects to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param removal whether to compute features on "masked" object fo each individual channels or on the globally detected object "MC".
#' Allowed are "masked" or "MC". Default is "masked". Please note that it will overwrite 'param' value if provided.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param nmax maximal order of Zernike polynomials to be computed. Default value is 9. Values outside [0,99] will be clipped.
#' Be aware that computation of Zernike's Moments can be quite long when 'nmax' is high.
#' @param granularity an integer vector. Controls the grain of the Haralick texture. Default is 3. Allowed are [1-20].
#' For very fine textures, this value is small (1-3 pixels), while for very coarse textures, it is large (>10).
#' @param batch number of objects to process at the same time. Default is 20.
#' @param parallel whether to use parallelization. Default is FALSE. Note that parallelization requires a parallel backend to be registered (see example).
#' In addition when parallelization is possible display_progress will be turned to FALSE.
#' @examples
#' msg_dat = character()
#' msg_par = character()
#' if(!requireNamespace("IFCdata", quietly = TRUE)) msg_dat = "IFCdata"
#' if(!requireNamespace("parallel", quietly = TRUE)) msg_par = c(msg_par, "parallel")
#' if(!requireNamespace("doParallel", quietly = TRUE)) msg_par = c(msg_par, "doParallel")
#' if(!requireNamespace("foreach", quietly = TRUE)) msg_par = c(msg_par, "foreach")
#' if(length(msg_dat) == 0) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", package = "IFCdata", "example.cif")
#'   ## features extraction
#'   ## the extraction is run for only objects 1 to 50 to allow example to run 
#'   ## in a reasonable amount of time to fulfill CRAN policies
#'   ## in current usage 'objects' argument may be missing to allow features extraction for all objects
#'   time_seq <- system.time({
#'     feat_seq <- ExtractFeatures(fileName = file_cif, 
#'                                 objects = 1:50,
#'                                 display_progress = TRUE,
#'                                 parallel = FALSE)
#'   })
#'   if(length(msg_par) == 0) {
#'     ## same extraction with a parallel backend,
#'     ## use of parallelization can clearly speed up the process,
#'     ## notably when Zernike features are extracted on a large amount of objects
#'     ## here is a small example that requires 'parallel', 'doParallel', and 'foreach' packages
#'     ## not installed along with 'IFCip' package
#'     library(parallel)
#'     library(doParallel)
#'     library(foreach)
#'     no_cores <- max(1, parallel::detectCores() - 1)
#'     ## the following is for R CMD cran check which allows to use at most 2 cores
#'     if("TRUE" %in% Sys.getenv("_R_CHECK_LIMIT_CORES_", "")) no_cores = min(2, no_cores)
#'     cl <- parallel::makePSOCKcluster(no_cores)
#'     doParallel::registerDoParallel(cl)
#'     time_par <- system.time({
#'        feat_par <- ExtractFeatures(fileName = file_cif, 
#'                                    objects = 1:50, 
#'                                    display_progress = TRUE, 
#'                                    parallel = TRUE)
#'     })
#'     parallel::stopCluster(cl)
#'     foreach::registerDoSEQ()
#'   } else {
#'     message(sprintf('Please run `install.packages(%s)` %s',
#'                     paste0("c(", paste0('"',msg_par,'"', collapse = ", "), ")"),
#'                    'to use parallelization.'))
#'   }
#' } else {
#'   message(paste0(sprintf('Please run `install.packages("IFCdata",repos="%s",type="source")` %s',
#'                  'https://gitdemont.github.io/IFCdata/',
#'                  'to install extra files required to run this example.'),
#'                  ifelse(length(msg_par) == 0, "" ,
#'                         sprintf('Please run `install.packages(%s)` %s',
#'                                 paste0("c(", paste0('"',msg_par,'"', collapse = ", "),")"),
#'                                'to use example of features extraction with parallelization.'))))
#' }
#' @details arguments of objectExtract() from IFC package will be deduced from \code{\link{ExtractFeatures}} input arguments.
#' @return a 3D array of features values whose dimensions are [object, features, channel] of class `IFCip_features`.
#' @export
ExtractFeatures <- function(...,
                            objects,
                            offsets,
                            removal = "masked",
                            display_progress = TRUE,
                            nmax = 9,
                            granularity = 3,
                            batch = 20,
                            parallel = FALSE)  {
  dots=list(...)
  
  # check input
  input = whoami(entries = as.list(match.call()))
  if(!any(sapply(input, FUN = function(i) length(i) != 0))) {
    stop("can't determine what to extract with provided parameters.\n try to input at least one of: 'fileName', 'info', 'param' or 'offsets'")
  }
  
  # reattribute needed param
  offsets = input[["offsets"]]
  param = input[["param"]]
  if(length(offsets) == 0) {
    fileName = input[["fileName"]]
  } else {
    fileName = attr(offsets, "fileName_image")
  }
  # check mandatory param
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  
  # process extra parameters
  if(length(dots[["verbose"]]) == 0) { 
    verbose = FALSE
  } else {
    verbose = dots[["verbose"]]
  }
  if(length(dots[["verbosity"]]) == 0) { 
    verbosity = 1
  } else {
    verbosity = dots[["verbosity"]]
  }
  if(length(dots[["fast"]]) == 0) { 
    fast = TRUE
  } else {
    fast = dots[["fast"]]
  }
  fast = as.logical(fast); assert(fast, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1, 2))
  batch = na.omit(as.integer(batch)); assert(verbosity, len = 1, typ = "integer")
  assert(removal, len=1, alw = c("masked", "MC"))
  param_extra = names(dots) %in% c("ifd","param","mode","export","size","force_width","removal","bypass","verbose")
  dots = dots[!param_extra] # remove not allowed param
  param_param = names(dots) %in% c("write_to","mode","base64_id","base64_att","overwrite",
                                   "composite","selection","random_seed","size","force_width",
                                   "removal","add_noise","full_range","force_range","spatial_correction")
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # compute object param
  # 1: prefer using 'param' if found,
  # 2: otherwise use 'info' if found,
  # 3: finally look at fileName
  if(length(param) == 0) {
    if(length(input$info) == 0) { 
      param = do.call(what = "objectParam",
                      args = c(list(fileName = fileName,
                                    mode = "raw",
                                    size = c(0,0),
                                    force_width = FALSE,
                                    removal = removal,
                                    export = "matrix"), dots_param))
    } else {
      param = do.call(what = "objectParam",
                      args = c(list(info = input$info,
                                    mode = "raw",
                                    size = c(0,0),
                                    force_width = FALSE,
                                    removal = removal,
                                    export = "matrix"), dots_param))
    }
  } else {
    param = input$param
    param$mode = "raw"
    param$export = "matrix"
    param$size = c(0,0)
    param$removal = rep(removal, length(param$chan_to_keep))
    param$channels$removal = rep(ifelse(removal == "masked", 3, 4), length(param$channels$removal))
    param$extract_msk = ifelse(removal == "masked", 3, 4)
  }
  fileName = param$fileName
  title_progress = basename(fileName)
  
  # check input offsets if any
  compute_offsets = TRUE
  if(length(offsets) != 0) {
    if(!("IFC_offset" %in% class(offsets))) {
      warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != checksumIFC(param$fileName_image)) {
        warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = param$fileName_image, fast = fast, display_progress = display_progress, verbose = verbose))
  }
  
  compute_mask <- FALSE
  if(param$XIF_test != 1) {
    compute_mask <- TRUE
  } else {
    ifd = getIFD(fileName = param$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = 0, image_type = "msk"), display_progress = FALSE)
    msk = objectExtract(ifd = ifd, param = param,  verbose = FALSE, bypass = TRUE)
    if(all(unname(unlist(msk)) == 0))  compute_mask <- TRUE
  }
  if(compute_mask) {
    param$removal = rep("none", length(param$chan_to_keep))
    param$channels$removal = rep(0, length(param$channels$removal))
    param$extract_msk = 0
    message("ExtractFeatures: can't find masks within file. They will be computed.")
  }
  
  # check objects to extract
  nobj = as.integer(attr(x = offsets, which = "obj_count"))
  
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  
  magnification = as.character(param$magnification)
  # for hu
  mag = switch(magnification, "20" = 1, "40" = 4, "60" = 9)
  # for shape
  k = switch(magnification, "20" = 1, "40" = 0.5, "60" = 0.3)
  
  ##### PAY ATTENTION TO MODIFY THIS PART IF NAMES ARE CHANGED IF DEDICATED FUNCTIONS
  # pre compute names in case object is not masked
  names_zernike = unlist(lapply(0:(nmax+1), FUN = function(a) { 
    foo = as.logical(sapply(0:a, FUN = function(b) {
      return (b %% 2)
    }))
    if(a %% 2) foo = !foo
    sprintf("zn%02im%02i", a-1, which(foo)-1)
  }))
  no_zernike = structure(rep(NA, length(names_zernike)), names = names_zernike)
  names_shape = c("Perimeter", "Diameter", "Circularity", "convexity", "roundness", 
                  "Height", "Width", "Elongatedness", "convex perimeter",
                  "convex cx", "convex cy")
  no_shape = structure(rep(NA, length(names_shape)), names = names_shape)
  names_hu = c("Area", #0.000
               "circularity", 
               "Minor Axis","Major Axis","Aspect Ratio",
               "Angle","theta","eccentricity",
               "Minor Axis Intensity","Major Axis Intensity","Aspect Ratio Intensity",
               "Angle Intensity","theta intensity","eccentricity intensity",
               "pix cx","pix cy","pix min axis","pix maj axis",
               "pix count", #0.000
               "inv1","inv2","inv3","inv4","inv5","inv6","inv7",
               "Raw Mean Pixel",
               "Raw Min Pixel","Raw Max Pixel", #+Inf,-Inf
               "Std Dev","skewness","kurtosis")
  no_hu = structure(c(0.000, rep(NaN, 17), 0.000, rep(NaN, 8), +Inf, -Inf, NaN, NaN, NaN), names = names_hu)
  #####
  show_pb = display_progress
  if(display_progress && 
     parallel &&
     requireNamespace("foreach", quietly = TRUE) &&
     (foreach::getDoParWorkers() > 1)) {
    show_pb = FALSE
    message("A parallel backend has been detected. 'display_progress' has been turned to FALSE")
  }
  
  if(parallel && requireNamespace("foreach", quietly = TRUE) && (foreach::getDoParWorkers() > 1)) {
    `%op%` <- foreach::`%dopar%`
  } else {
    `%op%` <- foreach::`%do%`
  }
  
  # extract objects
  sel = subsetOffsets(offsets = offsets, objects = objects, image_type = "img")
  sel = split(sel, ceiling(seq_along(sel)/batch))
  L=length(sel)
  if(L == 0) {
    warning("ExtractFeatures: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
    return(NULL)
  }
  tryCatch({
    if(show_pb) {
      pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
      ans = lapply(1:L, FUN=function(ifcip_iter) {
        setPB(pb, value = ifcip_iter, title = title_progress, label = "computing features from images")
        img = do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[ifcip_iter]],
                                                                         FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                         offset = off,
                                                                                                         trunc_bytes = 1, 
                                                                                                         force_trunc = TRUE, 
                                                                                                         verbose = verbose)),
                                                            param = param,
                                                            verbose = verbose,
                                                            bypass = TRUE),
                                                       dots))
        bar = lapply(img, FUN=function(i_img) {
          foo = lapply(i_img, FUN=function(i_chan) {
            if(compute_mask) {
              back = cpp_background(i_chan)
              bg_mean = back["BG_MEAN"]
              bg_sd = back["BG_STD"]
              msk = mask_identify(i_chan, 2 * bg_sd)
              msk_i = which.max(attr(msk, "perimeter"))
              if(length(msk_i) != 0) {
                msk = cpp_k_equal_M(msk, msk_i)
              } else {
                msk = msk
              }
            } else {
              msk = !attr(i_chan, "mask")
            }
            class(msk) = "IFC_msk"
            hu = cpp_features_hu3(img = i_chan, msk = msk, components = 1, mag = mag)
            if((nrow(hu) == 0) || !is.finite(hu[1,1]) || (hu[1,1] == 0)) {
              hu = no_hu
              shape = no_shape
            } else {
              hu = hu[1,]
              ctl = cpp_ctl(msk, global = TRUE)
              contours = ctl$contours
              contours = by(contours[, c(1,2,4,5)], contours[, 3], FUN =function(d) by(d[,c(1,2,3)], d[,4], FUN = function(dd) dd))
              contours = contours[as.integer(names(contours)) > 0] 
              contours = contours[[1]]
              if(inherits(contours, what = "by")) contours = contours[[1]]
              
              perimeter = k * sum(ctl$perimeter)
              # if(length(perimeter) == 0) perimeter = NA
              
              diameter = 2 * sqrt(hu["Area"] / pi)
              
              center = apply(contours[,1:2], 2, mean)
              center = hu[c("pix cy", "pix cx")]
              distance = k * apply(contours[,1:2], 1, FUN =function(coord)  sqrt((coord[1] - center[1])^2 + (coord[2] - center[2])^2))
              radius = mean(distance)
              circularity = radius / sd(distance)
              
              bbox = try(cpp_bbox(cpp_convexhull(as.matrix(contours)), k), silent = TRUE)
              if(inherits(x = bbox, what = "try-error")) {
                shape = structure(c(perimeter, diameter, circularity, rep(NA, 8)), names = names_shape)
              } else {
                convexity = bbox["convex perimeter"] / perimeter
                roundness = 4 * pi * hu["Area"] / bbox["convex perimeter"]^2
                shape = structure(c(perimeter, diameter, circularity, convexity, roundness, bbox), names = names_shape) 
              }
            }
            avg_intensity = hu["Raw Mean Pixel"] - attr(i_chan, "BG_MEAN")
            min_intensity = hu["Raw Min Pixel"] - attr(i_chan, "BG_MEAN")
            max_intensity = hu["Raw Max Pixel"] - attr(i_chan, "BG_MEAN")
            intensities = structure(c(attr(i_chan, "BG_MEAN"), attr(i_chan, "BG_STD"),
                                      min_intensity, max_intensity, avg_intensity, avg_intensity * hu["pix count"]), 
                                    names = c("Bkgd Mean", "Bkgd StdDev", "Min Pixel", "Max Pixel", "Mean Pixel", "Intensity"))
            # modulation TODO ask Amnis
            # max_intensity -  min_intensity / max_intensity + min_intensity is not working
            # modulation = (attr(img, "BG_MEAN") - (hu["Raw Max Pixel"] - hu["Raw Min Pixel"])) / ((hu["Raw Max Pixel"] + hu["Raw Min Pixel"])) 
            ze = try(moments_Zernike(img = i_chan, centroid = c(hu["pix cx"], hu["pix cy"]), radius = max(2, hu["pix maj axis"]/2+1), nmax = nmax, full = FALSE)$zmoment, silent = TRUE)
            if(inherits(x = ze, what = "try-error")) {
              ze = no_zernike
            }
            # return(c(hu, intensity, ze))
            har = compute_haralick(img = i_chan, msk = msk, granularity = granularity, bits = 4)
            return(c(hu, shape, intensities, ze, structure(unlist(har), names = paste0(apply(expand.grid(dimnames(har))[,c(2,1,3)], 1, paste0, collapse = " ")))))
          })
          attr(foo, "object_id") <- attr(i_img, "object_id")
          attr(foo, "offset_id") <- attr(i_img, "offset_id")
          attr(foo, "channel_id") <- attr(i_img, "channel_id")
          attr(foo, "removal") <- attr(i_img, "removal")
          return(foo)
        })
      })
    } else {
      if(show_pb) pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
      if(show_pb) setPB(pb, value = 0, title = title_progress, label = "progress bar will not update with parallel work but it is computing features from images")
      ans = list(foreach::foreach(ifcip_iter = 1:L, .combine = "c", .verbose = FALSE, .packages = c("IFC","IFCip"),
                         .export = c("cpp_features_hu3", "compute_haralick","assert","cpp_h_features","cpp_cooc","cpp_ctl",
                                     "cpp_background","cpp_closing","cpp_sd","cpp_fill","cpp_k_equal_M","mask_identify","mask_identify2","make_kernel",
                                     "cpp_getTAGS")) %op% { 
                           img = do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[ifcip_iter]],
                                                                                            FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                                            offset = off,
                                                                                                                            trunc_bytes = 1, 
                                                                                                                            force_trunc = TRUE, 
                                                                                                                            verbose = verbose)),
                                                                               param = param,
                                                                               verbose = verbose,
                                                                               bypass = TRUE),
                                                                          dots))
        bar = lapply(img, FUN=function(i_img) {
          foo = lapply(i_img, FUN=function(i_chan) {
            if(compute_mask) {
              back = cpp_background(i_chan)
              bg_mean = back["BG_MEAN"]
              bg_sd = back["BG_STD"]
              msk = mask_identify(i_chan, 2 * bg_sd)
              msk_i = which.max(attr(msk, "perimeter"))
              if(length(msk_i) != 0) {
                msk = cpp_k_equal_M(msk, msk_i)
              } else {
                msk = msk
              }
            } else {
              msk = !attr(i_chan, "mask")
            }
            class(msk) = "IFC_msk"
            hu = cpp_features_hu3(img = i_chan, msk = msk, components = 1, mag = mag)
            if((nrow(hu) == 0) || !is.finite(hu[1,1]) || (hu[1,1] == 0)) {
              hu = no_hu
              shape = no_shape
            } else {
              hu = hu[1,]
              ctl = cpp_ctl(msk, global = TRUE)
              contours = ctl$contours
              contours = by(contours[, c(1,2,4,5)], contours[, 3], FUN =function(d) by(d[,c(1,2,3)], d[,4], FUN = function(dd) dd))
              contours = contours[as.integer(names(contours)) > 0] 
              contours = contours[[1]]
              if(inherits(contours, what = "by")) contours = contours[[1]]
              
              perimeter = k * sum(ctl$perimeter)
              # if(length(perimeter) == 0) perimeter = 0
              
              diameter = 2 * sqrt(hu["Area"] / pi)
              
              center = apply(contours[,1:2], 2, mean)
              center = hu[c("pix cy", "pix cx")]
              distance = k * apply(contours[,1:2], 1, FUN =function(coord)  sqrt((coord[1] - center[1])^2 + (coord[2] - center[2])^2))
              radius = mean(distance)
              circularity = radius / sd(distance)
              
              bbox = try(cpp_bbox(cpp_convexhull(as.matrix(contours)), k), silent = TRUE)
              if(inherits(x = bbox, what = "try-error")) {
                shape = structure(c(perimeter, diameter, circularity, rep(NA, 8)), names = names_shape)
              } else {
                convexity = bbox["convex perimeter"] / perimeter
                roundness = 4 * pi * hu["Area"] / bbox["convex perimeter"]^2
                shape = structure(c(perimeter, diameter, circularity, convexity, roundness, bbox), names = names_shape) 
              }
            }
            
            avg_intensity = hu["Raw Mean Pixel"] - attr(i_chan, "BG_MEAN")
            min_intensity = hu["Raw Min Pixel"] - attr(i_chan, "BG_MEAN")
            max_intensity = hu["Raw Max Pixel"] - attr(i_chan, "BG_MEAN")
            intensities = structure(c(attr(i_chan, "BG_MEAN"), attr(i_chan, "BG_STD"),
                                      min_intensity, max_intensity, avg_intensity, avg_intensity * hu["pix count"]), 
                                    names = c("Bkgd Mean", "Bkgd StdDev", "Min Pixel", "Max Pixel", "Mean Pixel", "Intensity"))
            # modulation TODO ask Amnis
            # max_intensity -  min_intensity / max_intensity + min_intensity is not working
            # modulation = (attr(img, "BG_MEAN") - (hu["Raw Max Pixel"] - hu["Raw Min Pixel"])) / ((hu["Raw Max Pixel"] + hu["Raw Min Pixel"])) 
            ze = try(moments_Zernike(img = i_chan, centroid = c(hu["pix cx"], hu["pix cy"]), radius = max(2, hu["pix maj axis"]/2+1), nmax = nmax, full = FALSE)$zmoment, silent = TRUE)
            if(inherits(x = ze, what = "try-error")) {
              ze = no_zernike
            }
            har = compute_haralick(img = i_chan, msk = msk, granularity = granularity, bits = 4)
            return(c(hu, shape, intensities, ze, structure(unlist(har), names = paste0(apply(expand.grid(dimnames(har))[,c(2,1,3)], 1, paste0, collapse = " ")))))
          })
          attr(foo, "object_id") <- attr(i_img, "object_id")
          attr(foo, "offset_id") <- attr(i_img, "offset_id")
          attr(foo, "channel_id") <- attr(i_img, "channel_id")
          attr(foo, "removal") <- attr(i_img, "removal")
          return(foo)
        })
      })
    }
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  }, finally = {
    if(show_pb) endPB(pb)
  })
  channel_id = attr(ans[[1]][[1]], "channel_id")
  channel_removal = attr(ans[[1]][[1]], "removal")
  if(L > 1) {
    ans = do.call(what = "c", args = ans)
  } else {
    ans = ans[[1]]
  }
  channel_names = names(ans[[1]])
  ids = sapply(ans, attr, which = "object_id")
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  ret = aperm(array(unlist(ans), dim = c(length(ans[[1]][[1]]), 
                                         length(ans[[1]]),
                                         length(objects))),
              perm = c(3, 1, 2))
  dimnames(ret) = list("object" = ids,
                       "feature" = names(ans[[1]][[1]]),
                       "channel" = channel_id)
  attr(ret, "offset_id") <- sapply(ans, attr, which = "offset_id")
  attr(ret, "channel_id") <- channel_id
  attr(ret, "channel_names") <- channel_names
  attr(ret, "removal") <- channel_removal
  class(ret) = "IFCip_features"
  return(ret)
}
