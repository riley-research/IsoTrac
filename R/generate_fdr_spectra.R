generate_fdr_spectra <- function(full_mgf, fdr_approach = "supervised") {
  full_mgf <- full_mgf |>
    dplyr::filter(!is.na(.data$centroided_intensity))

  if (fdr_approach == "supervised") {
    full_mgf <- full_mgf |>
      dplyr::mutate(
                    .by = c("id", "slice"),
                    centroided_intensity =
                      supervised_intensity(
                                           .data$centroided_intensity))
  }

  full_mgf
}

supervised_intensity <- function(intensities) {
  n <- length(intensities)

  min_i  <- min(intensities, na.rm = TRUE)
  max_i  <- max(intensities, na.rm = TRUE)
  range_i <- max_i - min_i

  if (range_i <= 0) {
    return(intensities)
  }

  new_i <- numeric(n)

  for (i in seq_len(n)) {
    phase <- (i - 1) %% 3

    if (phase == 0) {
      new_i[i] <- runif(1, min_i, min_i + range_i / 3)

    } else if (phase == 1) {
      new_i[i] <- runif(1, new_i[i - 1], min_i + 2 * range_i / 3)

    } else {
      new_i[i] <- runif(1, new_i[i - 1], max_i)
    }
  }

  new_i
}
