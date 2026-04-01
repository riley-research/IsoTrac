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

    new_i[i] <- runif(1, lower_bound, upper_bound)
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
