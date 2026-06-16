#' Generate FDR Decoy Spectra Intensities
#'
#' @description
#' Prepares mass spectrometry data for False Discovery Rate (FDR) estimation by manipulating
#' peak intensities. When using the default `"supervised"` approach, it replaces the original
#' centroided peak intensities with simulated values that mimic an upward trend over cycles
#' defined by expected isotopic patterns.
#'
#' @param full_mgf \code{data.frame}. The data imported using the \code{\link{import_spectra}}
#' function.
#' @param isotope_peaks_included \code{integer}. The number of peaks expected per isotopic cluster
#' phase (used as the cycle length/bucket count for the simulated intensity trends).
#' @param fdr_approach \code{character}. The methodology used for decoy generation. Currently supports
#' \code{"supervised"} to simulate bounded, cyclic trend intensities. Default is \code{"supervised"}.
#'
#' @return A \code{data.frame} filtered for valid centroid entries, containing either the original
#' or the newly simulated decoy \code{centroided_intensity} values grouped by spectra ID and slice.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming `processed_mgf` is the output from import_spectra
#' decoy_mgf <- generate_fdr_spectra(
#'   full_mgf = processed_mgf,
#'   isotope_peaks_included = 3,
#'   fdr_approach = "supervised"
#' )
#' }
generate_fdr_spectra <- function(full_mgf, isotope_peaks_included, fdr_approach = "supervised") {
  full_mgf <- full_mgf |>
    dplyr::filter(!is.na(.data$centroided_intensity))

  if (fdr_approach == "supervised") {
    full_mgf <- full_mgf |>
      dplyr::mutate(
                    .by = c("id", "slice"),
                    centroided_intensity =
                      supervised_intensity(
                                           .data$centroided_intensity,
                                           isotope_peaks_included))
  }

  full_mgf
}

supervised_intensity <- function(intensities, isotope_peaks_included) {
  n <- length(intensities)

  min_i  <- min(intensities, na.rm = TRUE)
  max_i  <- max(intensities, na.rm = TRUE)
  range_i <- max_i - min_i

  if (range_i <= 0) {
    return(intensities)
  }

  new_i <- numeric(n)
  num_phases <- isotope_peaks_included
  phase_width <- range_i / num_phases

  for (i in seq_len(n)) {
    # Determine which phase we are in (0 to num_phases - 1)
    phase_idx <- (i - 1) %% num_phases

    # Calculate the lower bound:
    # For the first phase of a cycle, use min_i.
    # Otherwise, stay above the previous value to maintain the upward trend.
    lower_bound <- if (phase_idx == 0) min_i else new_i[i - 1]

    # Calculate the upper bound:
    # (phase_idx + 1) * phase_width gives the top of the current bucket
    upper_bound <- min_i + (phase_idx + 1) * phase_width

    # Ensure numerical stability (bounds can't cross)
    if (lower_bound > upper_bound) lower_bound <- upper_bound - 1e-10

    new_i[i] <- stats::runif(1, lower_bound, upper_bound)
  }

  # for (i in seq_len(n)) {
  #   phase <- (i - 1) %% 3
  #
  #   if (phase == 0) {
  #     new_i[i] <- runif(1, min_i, min_i + range_i / 3)
  #
  #   } else if (phase == 1) {
  #     new_i[i] <- runif(1, new_i[i - 1], min_i + 2 * range_i / 3)
  #
  #   } else {
  #     new_i[i] <- runif(1, new_i[i - 1], max_i)
  #   }
  # }

  new_i
}
