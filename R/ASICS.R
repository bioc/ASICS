#' Automatic Statistical Identification in Complex Spectra
#'
#' Quantification of 1D 1H NMR spectra with ASICS method using a library of
#' pure metabolite spectra. The method is presented in Tardivel et al. (2017).
#'
#' @param spectra_obj An object of class \linkS4class{Spectra} obtained with the
#' function \link{createSpectra}.
#' @param exclusion.areas Definition domain of spectra that has to be excluded
#' for the quantification (ppm). By default, the water region is excluded
#' (4.5-5.1 ppm).
#' @param max.shift Maximum chemical shift allowed (in ppm). Default to 0.02.
#' @param pure.library An object of class \linkS4class{PureLibrary} containing
#' the reference spectra (pure metabolite spectra). If \code{NULL}, the library
#' included in the package (that contains 191 reference spectra) is used.
#' @param noise.thres Threshold for signal noise. Default to 0.02.
#' @param threshold.noise DEPRECATED, use \code{noise.thres} instead.
#' @param joint.align Logical. If \code{TRUE}, information from all spectra is
#' taken into account to align individual library.
#' @param combine DEPRECATED, use \code{joint.align} instead.
#' @param add.noise,mult.noise additive and multiplicative noises. To set these
#' noises, you can compute the standard deviation in a noisy area for
#' \code{add.noise} or the standard deviation in a peak area for
#' \code{mult.noise} when several spectra of the same sample are available.
#' By default, \code{add.noise = 0.15} and \code{mult.noise = 0.172}
#' @param quantif.method either \code{"FWER"} to perform an independent
#' quantification (the method available in ASICS since the beginning),
#' \code{"Lasso"} to perform a joint quantification (all the spectra together)
#' or \code{"both"} to perform a joint quantification after the FWER selection
#' of the independent quantification. More details can be founded in the
#' user's guide.
#' @param clean.thres if \code{quantif.method == "both"} the percentage of
#' spectra in which the metabolite needs to be identified by the FWER selection.
#' Default to 1, \emph{i.e.} metabolite is quantified if it was identified
#' in at least 1\% of the spectra.
#' @param ref.spectrum index of the reference spectrum used for the alignment.
#' Default to \code{NULL}, \emph{i.e.} the reference spectrum is automatically
#' detected.
#' @param seed Random seed to control randomness in the algorithm (used in the
#' estimation of the significativity of a given metabolite concentration).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param verbose A Boolean value to allow print out process information.
#'
#' @note Since version 2.3.1 small changes were applied in order to improve the
#' speed of metabolite selection algorithm, which can slightly impact outputs
#' of the method.
#'
#' @return An object of type \linkS4class{ASICSResults} containing the
#' quantification results.
#'
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SerialParam
#' @importFrom stats reshape
#' @export
#'
#' @seealso \linkS4class{ASICSResults} \code{\link{pure_library}}
#' \code{\link{createSpectra}}
#'
#' @references Tardivel P., Canlet C., Lefort G., Tremblay-Franco M., Debrauwer
#' L., Concordet D., Servien R. (2017). ASICS: an automatic method for
#' identification and quantification of metabolites in complex 1D 1H NMR
#' spectra. \emph{Metabolomics}, \strong{13}(10): 109.
#' \url{https://doi.org/10.1007/s11306-017-1244-5}
#'
#' @examples
#' # Import data and create object
#' current_path <- system.file("extdata", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path,
#'                      name.file = "spectra_example.txt", type.import = "txt")
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Estimation of relative quantifications
#' to_exclude <- matrix(c(4.5, 10), ncol = 2)
#' resASICS <- ASICS(spectra_obj, exclusion.areas = to_exclude,
#'                   joint.align = FALSE, quantif.method = "FWER")
ASICS <- function(spectra_obj,
                  exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                  max.shift = 0.02, pure.library = NULL,
                  noise.thres = 0.02, joint.align = TRUE,
                  threshold.noise = NULL, combine = NULL,
                  add.noise = 0.15, mult.noise = 0.172,
                  quantif.method = c("FWER", "Lasso", "both"),
                  clean.thres = 1, ref.spectrum = NULL,
                  seed = 1234, ncores = 1, verbose = TRUE) {

  quantif.method <- match.arg(quantif.method)
  
  if(!is.null(exclusion.areas) &&
     (!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2)){
    stop("'exclusion.areas' must be a matrix with 2 columns.")
  }

  if(max.shift < 0){
    stop("'max.shift' must be non negative.")
  }

  if(!is.null(threshold.noise)){
    message("'threshold.noise' is depreceated, use 'noise.thres' instead")
    noise.thres <- threshold.noise
  }

  if(!is.null(combine)){
    message("'combine' is depreceated, use 'joint.align' instead")
    joint.align <- combine
  }

  if(noise.thres < 0){
    stop("'noise.thres' must be non negative.")
  }

  if(!is(pure.library, "PureLibrary") & !is.null(pure.library)){
    stop(paste("'pure.library' must be either NULL or an object of class",
               "'PureLibrary'."))
  }

  res_estimation <- .ASICSInternal(spectra_obj, exclusion.areas, max.shift,
                                   pure.library, noise.thres,  seed, ncores,
                                   joint.align, verbose, add.noise, mult.noise,
                                   quantif.method, clean.thres, ref.spectrum)

  return(res_estimation)
}



#' @importFrom methods new
.ASICSInternal <- function(spectra_obj_raw,
                           exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                           max.shift = 0.02, pure.library = NULL,
                           noise.thres = 0.02, seed = 1234, ncores = 1,
                           joint.align = TRUE, verbose = TRUE,
                           add.noise = 0.15, mult.noise = 0.172,
                           quantif.method = c("FWER", "Lasso", "both"),
                           clean.thres = 1, ref.spectrum = NULL) {

    # seed and parallel environment
    set.seed(seed)

    # default library or not
    if(is.null(pure.library)){
      pure.library <- ASICS::pure_library
    }

    # spectra object as a list where each element are 1 spectrum
    spectra_list <- lapply(seq_along(spectra_obj_raw),
                           function(x) spectra_obj_raw[x])


    #-----------------------------------------------------------------------------
    #### Remove areas from spectrum and library ####
    if (verbose) cat("Remove areas from spectrum and library \n")
    spectra_obj <- bplapply(spectra_list, .removeAreas, exclusion.areas,
                            pure.library,
                            BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                                 verbose))

    # number of points on library grid corresponding to maximum shift
    if (length(spectra_list) == 1 | !joint.align) {
      nb_points_shift <-
        floor(max.shift / (spectra_obj[[1]][["cleaned_library"]]@ppm.grid[2] -
                             spectra_obj[[1]][["cleaned_library"]]@ppm.grid[1]))
    } else {
      max.shift <- seq_len(5) * max.shift / 5
      nb_points_shift <-
        floor(max.shift / (spectra_obj[[1]][["cleaned_library"]]@ppm.grid[2] -
                             spectra_obj[[1]][["cleaned_library"]]@ppm.grid[1]))
    }

    # compute weights
    if (verbose) cat("Compute weights \n")
    spectra_obj <-
      bplapply(spectra_obj,
               function(x){x[["mixture_weights"]] <-
                 as.numeric(1 / (abs(x[["cleaned_spectrum"]]@spectra) *
                                   mult.noise ^ 2 + add.noise ^ 2)); return(x)},
               BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                    verbose))

    #-----------------------------------------------------------------------------
    #### Cleaning step: remove metabolites that cannot belong to the mixture ####
    pure_lib_raw <- spectra_obj[[1]]$cleaned_library
    if (quantif.method %in% c("FWER", "both")) {
      if (verbose) cat("Remove metabolites that cannot belong to the mixture \n")
      spectra_obj <- bplapply(spectra_obj, .cleanLibrary, noise.thres,
                              nb_points_shift[length(nb_points_shift)],
                              BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                                   verbose))
    }

    if (quantif.method == "Lasso") {
      if (verbose) cat("Remove metabolites that cannot belong to the mixture \n")
      spectra_obj <- bplapply(spectra_obj, .cleanLibrary, noise.thres,
                              nb_points_shift[length(nb_points_shift)],
                              BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                                   verbose))

      metab_list <- unlist(bplapply(spectra_obj,
                                    function(x) x$cleaned_library@sample.name,
                                    BPPARAM = .createEnv(ncores,
                                                         length(spectra_obj_raw),
                                                         verbose)))
      metab_list_th <-
        names(table(metab_list)[table(metab_list) >= clean.thres])

      spectra_obj <- bplapply(spectra_obj,
            function(x) {
              x$cleaned_library <-
                pure_lib_raw[pure_lib_raw@sample.name %in% metab_list_th];
              return(x)},
            BPPARAM = .createEnv(ncores, length(spectra_obj_raw), verbose))
    }

      #-----------------------------------------------------------------------------
      #### Find the best shift between each pure spectra and mixture ####
      #and sort metabolites by regression residuals
      if (quantif.method %in% c("FWER", "Lasso", "both")) {
        if (verbose) cat("Translate library \n")
        if (length(spectra_list) == 1 | !joint.align) {
          spectra_obj <- bplapply(spectra_obj, .translateLibrary,
                                    nb_points_shift[length(nb_points_shift)],
                                    max.shift[length(max.shift)],
                                    BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                                         verbose))
        } else {
          # spectra binning
          spectra_to_bin <- data.frame(as.matrix(getSpectra(spectra_obj_raw)))
          rownames(spectra_to_bin) <- spectra_obj_raw@ppm.grid
          norm.param <- c(list(spectra = spectra_to_bin,
                               exclusion.areas = exclusion.areas,
                               ncores = ncores,
                               bin = 0.001,
                               verbose = FALSE,
                               type.norm = spectra_obj_raw@norm.method),
                          spectra_obj_raw@norm.params)
          spec_bin <- do.call("binning", norm.param)
          spec_bin <- spec_bin[rowSums(spec_bin) != 0, ]

          spectra_obj <- .translateLibrary_combineVersion(spectra_obj, max.shift,
                                                          nb_points_shift, spec_bin,
                                                          pure.library, ncores,
                                                          length(spectra_obj_raw),
                                                          verbose)
        }
      }

      #-----------------------------------------------------------------------------
      #### Localized deformations of pure spectra ####
      if (quantif.method %in% c("FWER", "both")) {
        if (verbose) cat("Deform library peaks \n")
        spectra_obj <- bplapply(spectra_obj, .deformLibrary,
                                BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                                     verbose))

      }

    #-----------------------------------------------------------------------------
    #### Threshold and concentration optimisation for each metabolites
    # with FWER method ####
    if (quantif.method %in% c("FWER", "both")) {
      if (verbose) cat("Compute quantifications \n")
      spectra_obj <- bplapply(spectra_obj, .concentrationOpti_FWER,
                            BPPARAM = .createEnv(ncores, length(spectra_obj_raw),
                                                 verbose))
    }

    #-----------------------------------------------------------------------------
    #### Cleaning step (on FWER results): remove metabolites that cannot ####
    #### belong to the mixture ####
    if (quantif.method == "both") {
      if (verbose) cat("Remove metabolites that cannot belong to the mixture \n")

      metab_list <- unlist(bplapply(spectra_obj,
                                    function(x) x$cleaned_library@sample.name,
                                    BPPARAM = .createEnv(ncores,
                                                         length(spectra_obj_raw),
                                                         verbose)))
      metab_list_th <-
        names(table(metab_list)[table(metab_list) >= clean.thres])
    }

  #-----------------------------------------------------------------------------
  #### (method: Lasso + both) Mean translation between each pure
  #spectra and mixture and sort metabolites by regression residuals
  if (quantif.method %in% c("Lasso", "both")) {
    pure_lib <- pure_lib_raw[pure_lib_raw@sample.name %in% metab_list_th]

    if (verbose) cat("Translate library of the new cleaned spectra \n")
    shift <- bplapply(as.list(seq_along(spectra_obj)),
                      function(x) {
        res <- data.frame(metab = spectra_obj[[x]]$cleaned_library@sample.name,
                          shift = spectra_obj[[x]]$shift);
        colnames(res)[2] <- paste0("shift", x);
        return(res)},
                      BPPARAM = .createEnv(ncores,
                                           length(spectra_obj_raw),
                                           verbose))

    shifts <- join_all(shift, by = "metab", type = "full")
    rownames(shifts) <- shifts$metab
    shifts$metab <- NULL

    shift_by_metab <- apply(shifts, 1, median, na.rm = TRUE)
    shift_by_metab <-
      floor(shift_by_metab /
              (spectra_obj[[1]][["cleaned_library"]]@ppm.grid[2] -
                 spectra_obj[[1]][["cleaned_library"]]@ppm.grid[1]))

    shift_by_metab <-
      shift_by_metab[names(shift_by_metab) %in% pure_lib@sample.name]

    pure_lib@spectra <-
      Matrix(do.call("cbind",
                     lapply(seq_along(pure_lib),
                            function(idx)
                              .doShift(pure_lib@spectra[, idx],
                                       shift_by_metab[pure_lib@sample.name[idx]],
                                       0, 0))))
  }

  #### Localized deformations of pure spectra ####
  if (quantif.method %in% c("Lasso", "both")) {
    if (verbose & is.null(ref.spectrum)) cat("Find the reference spectrum \n")

    if (is.null(ref.spectrum)) {
      spectra <-
        as.matrix(do.call("cbind",
                        lapply(spectra_obj,
                               function(x) return(x[["cleaned_spectrum"]]@spectra))))
      ref.spectrum <- .findReference(spectra, ncores, verbose)
    }

    #max_shift for reference
    max_shift_all <- bplapply(as.list(seq_along(spectra_obj)),
                              function(x) {res <- data.frame(metab = spectra_obj[[x]]$cleaned_library@sample.name,
                                                             max_shift = spectra_obj[[x]]$max_shift); colnames(res)[2] <- paste0("shift", x);
                                                             return(res)},
                              BPPARAM = .createEnv(ncores,
                                                   length(spectra_obj_raw),
                                                   verbose))

    max_shift_all <- join_all(max_shift_all, by = "metab", type = "full")
    rownames(max_shift_all) <- max_shift_all$metab
    max_shift_all$metab <- NULL

    max_shift_by_metab <- apply(max_shift_all, 1, max, na.rm = TRUE)
    max_shift_by_metab <-
      max_shift_by_metab[names(max_shift_by_metab) %in% pure_lib@sample.name]


    #nb_point_shift for reference
    nb_point_all <- bplapply(as.list(seq_along(spectra_obj)),
                              function(x) {res <- data.frame(metab = spectra_obj[[x]]$cleaned_library@sample.name,
                                                             max_nb_points_shift = spectra_obj[[x]]$max_nb_points_shift); colnames(res)[2] <- paste0("shift", x);
                                                             return(res)},
                              BPPARAM = .createEnv(ncores,
                                                   length(spectra_obj_raw),
                                                   verbose))

    nb_point_all <- join_all(nb_point_all, by = "metab", type = "full")
    rownames(nb_point_all) <- nb_point_all$metab
    nb_point_all$metab <- NULL

    nb_point_by_metab <- apply(nb_point_all, 1, max, na.rm = TRUE)
    nb_point_by_metab <-
      nb_point_by_metab[names(nb_point_by_metab) %in% pure_lib@sample.name]

    spectra_obj[[ref.spectrum]]$cleaned_library <- pure_lib
    spectra_obj[[ref.spectrum]]$shift <- shift_by_metab[pure_lib@sample.name]
    spectra_obj[[ref.spectrum]]$max_shift <- max_shift_by_metab[pure_lib@sample.name]
    spectra_obj[[ref.spectrum]]$max_nb_points_shift <- nb_point_by_metab[pure_lib@sample.name]

    if (verbose) cat("Deform library peaks \n")
    spectra_obj[[ref.spectrum]] <- .deformLibrary(spectra_obj[[ref.spectrum]])
  }


  #### Threshold and concentration optimisation for each metabolites
  # with Lasso method ####
  if (quantif.method %in% c("Lasso", "both")) {
    if (verbose) cat("Compute quantifications \n")
    spectra_obj <- .concentrationOpti_MRLasso(spectra_obj, ref.spectrum,
                                              ncores, verbose)
  }

  #-----------------------------------------------------------------------------
  #### Results ####
  if (verbose) cat("Format global results... \n")
  sample_name <-
    unlist(vapply(spectra_obj,
                  function(x) return(x[["cleaned_spectrum"]]@sample.name),
                  "character"))
  spectra <-
    do.call("cbind",
            lapply(spectra_obj,
                   function(x) return(x[["cleaned_spectrum"]]@spectra)))
  rec_spectra <-
    do.call("cbind", lapply(spectra_obj,
                            function(x) return(x[["est_mixture"]])))

  rel_conc <- lapply(spectra_obj, function(x) {
    x[["relative_concentration"]]$row_names <-
      x[["cleaned_library"]]@sample.name ;
    return(x[["relative_concentration"]])})
  metab_conc <- join_all(rel_conc, by = "row_names", type = "full")
  rownames(metab_conc) <- metab_conc$row_names
  metab_conc$row_names <- NULL
  metab_conc[is.na(metab_conc)] <- 0
  pure_lib_format <- do.call("rbind",
                             lapply(spectra_obj,
                                    function(x) return(x[["format_library"]])))
  # Object to return
  res_object <- new(Class = "ASICSResults",
                    sample.name = sample_name,
                    ppm.grid = spectra_obj[[1]][["cleaned_spectrum"]]@ppm.grid,
                    spectra = spectra,
                    reconstructed.spectra = rec_spectra,
                    quantification = metab_conc,
                    deformed.library = pure_lib_format)

  return(res_object)
}
