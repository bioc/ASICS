#' Simulate a set of spectra
#'
#' Simulate a set of spectra based on the default library with shifts
#'
#' @param n.spectra Number of spectra to simulate.
#' @param max.shift Maximum shift allowed for artificial deformation of pure
#' spectra (default to 0.02).
#' @param metab.percent Percentage of present metabolites in complex spectra
#' (default to 0.5).
#' @param metab.different Number of metabolites that are different between each 
#' complex spectra (default to 4).
#' @param add.noise,mult.noise additive and multiplicative noises. By default,
#' \code{add.noise = 0.15} and \code{mult.noise = 0.172}
#'
#' @return A list with a data frame of simulated spectra in columns and
#' a data frame of simulated quantifications.
#'
#' @export
#' @importFrom plyr ldply adply
#' @importFrom stats rbinom runif rnorm
#'
#' @examples
#' spectra <- simulate_spectra(n.spectra = 10)
simulate_spectra <- function(n.spectra, max.shift = 0.02, metab.percent = 0.5,
                             metab.different = 4, add.noise = 0.07,
                             mult.noise = 0.09) {
  # number of metabolite in the library
  p <- length(ASICS::pure_library)

  # present metabolites in the complex mixture
  present <- rbinom(p, 1, metab.percent)
  # for each sample : add or remove some metabolites
  present_all <- t(ldply(as.list(1:n.spectra),
      function(x) {
        to_add_or_remove <- sample(1:length(ASICS::pure_library), metab.different)
        present[to_add_or_remove] <-
          ifelse(present[to_add_or_remove] == 0, 1, 0)
        return(present)} ))

  ##### Simulate metabolite concentrations

  # vector of log normal means and standard deviations
  log_mu <- rlnorm(p, -8, 2)
  log_sd <- 0.3 * log_mu

  # matrix of coefficients
  params_coeff <- data.frame(idx = 1:p, mu = log_mu, sd = log_sd,
                             present = present_all)
  coeff <- adply(params_coeff, 1, .simu_coeff,
                 n.spectra, metab.different,
                 .expand = FALSE)[, 2:(n.spectra + 1)]

  ##### Simulate shifted library and create simulated spectra

  # maximum shift in points
  max_points <- max.shift / diff(ASICS::pure_library@ppm.grid)[1] #max shift in points

  spectra <- matrix(NA, nrow = length(ASICS::pure_library@ppm.grid), ncol = n.spectra)

  all_shift <- matrix(NA, nrow = p, ncol = n.spectra)
  all_shift_lib <- list()
  for (i in seq_len(n.spectra)) {
    res <- llply(as.list(1:p), function(x) .shift_spec(x, max_points, coeff[, i]))

    shift_lib <- t(ldply(as.list(1:p), function(x) res[[x]]$spectrum))
    shift_lib_raw <- t(ldply(as.list(1:p), function(x) res[[x]]$raw_spectrum))
    all_shift_lib[[i]] <- Matrix(shift_lib_raw)
    all_shift[, i] <- laply(as.list(1:p), function(x) res[[x]]$shift)

    spectra[, i] <- apply(shift_lib, 1, sum)
    spectra[, i] <- spectra[, i]  / .AUC(ASICS::pure_library@ppm.grid, spectra[, i])
    spectra[, i] <- spectra[, i] + spectra[, i] *
      rnorm(length(spectra[, i]), 0, mult.noise^2) +
      rnorm(length(spectra[, i]), 0, add.noise^2)
  }

  spectra <- as.data.frame(spectra)
  rownames(spectra) <- ASICS::pure_library@ppm.grid

  return(list(spectra = spectra, sim_quantif = coeff, all_shift = all_shift,
              all_shift_lib = all_shift_lib))
}

## Simulate coefficients
#' @importFrom stats rlnorm
#' @keywords internal
.simu_coeff <- function(param, n.spectra, metab.different) {
  present <- param[, 4:(n.spectra + 3)]

  coeff_sim <- rnorm(n.spectra, param$mu, param$sd)
  # remove values that are too high
  coeff_sim[coeff_sim > 1] <- 0
  # remove values that are too low (metabolites non present)
  coeff_sim[coeff_sim < 0] <- 0

  return(ifelse(present, coeff_sim, 0))
}

## Create library shifts
#' @keywords internal
#' @importFrom stats rnbinom
.shift_spec <- function(idx, max_points, coeff_ind) {
  spectrum <- ASICS::pure_library@spectra[, idx]
  coeff_spec <- coeff_ind[idx]

  pos_neg <- rbinom(1, 1, 1/2)

  global_shift <- rnbinom(1, 2, 1/4)
  global_shift <- ifelse(global_shift > max_points, max_points, global_shift)
  global_shift <- ifelse(pos_neg == 0, - global_shift, global_shift)

  spectrum <- .doShift(spectrum, global_shift)

  # # detect peak bounds for local shift
  # # expanded connected components
  signal_peak_lib <- which(spectrum > 1)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] -
    floor(max_points / 10)
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] +
    floor(max_points / 10)
  peaks_extremities <-
    cbind(vapply(min_extremities, max, 1,
                 FUN.VALUE = numeric(1)),
          vapply(max_extremities, min, length(ASICS::pure_library@ppm.grid),
                 FUN.VALUE = numeric(1)))

  # remove overlapping
  long_signal <- unique(unlist(apply(peaks_extremities, 1,
                                     function(x) x[[1]]:x[[2]])))

  min_extremities <- long_signal[!((long_signal - 1) %in% long_signal)]
  max_extremities <- long_signal[!((long_signal + 1) %in% long_signal)]
  peaks_extremities <- cbind(min_extremities, max_extremities)

  # remove small peak (less than 3 points)
  peaks_extremities <-
    matrix(peaks_extremities[peaks_extremities[, 2] -
                               peaks_extremities[, 1] >= 3, ], ncol = 2)

  if (nrow(peaks_extremities) > 0) {
    # deform on each connected component
    for(peak in seq_len(nrow(peaks_extremities))){
      # simulate a parameter for deformation
      a <- rnorm(1, mean = 0, sd = 0.09)
      a[a < -1] <- -1
      a[a > 1] <- 1

      # area of peak to deforme
      peak_area <- peaks_extremities[peak, 1]:peaks_extremities[peak, 2]
      # peak to deform
      to_deform <- spectrum[peak_area]
      # grid to deform
      grid_to_deform <- ASICS::pure_library@ppm.grid[peak_area]

      spectrum[peak_area] <- .deforme(grid_to_deform, to_deform, a)
    }
  }

  spectrum_coeff <- spectrum * coeff_spec * ASICS::pure_library@nb.protons[idx]

  return(list(spectrum = spectrum_coeff, raw_spectrum = spectrum,
              shift = global_shift))
}
