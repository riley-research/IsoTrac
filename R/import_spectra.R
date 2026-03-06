import_spectra <- function(mgf_path, intensity_threshold, charge_range,
                           separation_threshold = 20, min_peaks_per_slice = 6,
                           min_charge_percentage_matched = 0.2,
                           scan_range = c(0, 0)) {
  full_mgf <- utils::read.delim(mgf_path, sep = ",")

  if (!identical(scan_range, c(0, 0))) {
    full_mgf <- scan_range_filter(full_mgf, scan_range)
  }

  full_mgf <- intensity_filter(full_mgf,
                               intensity_threshold = intensity_threshold)

  full_mgf <- min_peak_filter(full_mgf)

  full_mgf <- full_mgf |>
    dplyr::mutate(.by = "id",
                  slice = cumsum(c(TRUE,
                                   diff(.data$mz) >= separation_threshold)))

  full_mgf <- min_peaks_per_slice_filter(full_mgf, min_peaks_per_slice)

  full_mgf <- centroid_data(full_mgf)

  full_mgf <- calculate_differences(full_mgf, charge_range,
                                    min_charge_percentage_matched)

  full_mgf <- full_mgf |>
    dplyr::mutate(
      estimated_intact_mass = estimate_intact_mass(dplyr::pick("mz", "charge")),
      .by = c("id", "slice")
    )


  full_mgf
}

scan_range_filter <- function(spectrum, scan_range) {
  to_keep <- spectrum |>
    dplyr::mutate(spectrum_nr = as.numeric(gsub("Spectrum", "", .data$id))) |>
    dplyr::filter(.data$spectrum_nr >= scan_range[1] &
                    .data$spectrum_nr <= scan_range[2]) |>
    dplyr::pull(.data$id)

  spectrum <- spectrum |>
    dplyr::filter(.data$id %in% to_keep)

  spectrum
}

intensity_filter <- function(spectrum, intensity_threshold) {
  keep <- spectrum |>
    dplyr::filter(.data$intensity >= intensity_threshold) |>
    dplyr::pull(.data$peakgroup)

  spectrum <- spectrum |>
    dplyr::filter(.data$peakgroup %in% keep)

  spectrum
}

min_peak_filter <- function(spectrum) {
  num_to_keep <- 2

  keep_these <- spectrum |>
    dplyr::summarise(.by = c("id", "peakgroup"),
                     count = dplyr::n()) |>
    dplyr::filter(.data$count >= num_to_keep) |>
    dplyr::distinct(.data$id, .data$peakgroup)

  filtered_df <- spectrum |>
    dplyr::semi_join(keep_these, by = c("id", "peakgroup"))

  filtered_df
}

min_peaks_per_slice_filter <- function(spectrum, min_peaks_per_slice) {
  num_to_keep <- min_peaks_per_slice

  keep_these <- spectrum |>
    dplyr::distinct(.data$id, .data$peakgroup, .data$slice) |>
    dplyr::summarise(.by = c("id", "slice"),
                     count = dplyr::n()) |>
    dplyr::filter(.data$count >= num_to_keep) |>
    dplyr::distinct(.data$id, .data$slice)

  filtered_df <- spectrum |>
    dplyr::semi_join(keep_these, by = c("id", "slice"))

  filtered_df
}

centroid_data <- function(spectrum) {
  # Ensure column exists
  if (!"centroided_intensity" %in% names(spectrum)) {
    spectrum$centroided_intensity <- NA_real_
  }

  # Group by id + peakgroup and process each group
  results <- spectrum |>
    dplyr::group_by(.data$id, .data$peakgroup) |>
    dplyr::group_map(~ {
      g <- .
      spline_interp <- stats::splinefun(g$mz, g$intensity, method = "fmm")

      opt_result <- optimize(
        function(x) -spline_interp(x),
        interval = range(g$mz)
      )

      # Take the row with max intensity to copy other metadata
      line_to_add <- g |> dplyr::slice(1)
      line_to_add$centroided_intensity <- -opt_result$objective
      line_to_add$intensity <- line_to_add$centroided_intensity
      line_to_add$mz <- opt_result$minimum
      line_to_add
    }, .keep = TRUE)

  # Combine all groups
  results_df <- dplyr::bind_rows(results)

  # Add centroided points to original spectrum
  spectrum <- dplyr::bind_rows(spectrum, results_df) |>
    dplyr::mutate(spectrum_nr = as.numeric(gsub("Spectrum", "", .data$id))) |>
    dplyr::arrange(.data$spectrum_nr, .data$mz) |>
    dplyr::select(-"spectrum_nr")

  spectrum
}

calculate_differences <- function(spectrum, charge_range,
                                  min_charge_percentage_matched) {
  w <- data.frame(charge = charge_range,
                  delta = 1.00784 / charge_range)

  tempdf <- spectrum |>
    dplyr::filter(!is.na(.data$centroided_intensity)) |>
    dplyr::distinct(.data$id, .data$slice, .data$mz) |>
    dplyr::summarise(
      .by = c("id", "slice"),
      mz = list(as.double(.data$mz))
    ) |>
    dplyr::mutate(
      difVal = purrr::map(.data$mz, diff)
    )

  what_col_names <- as.character(w$charge)
  for (i in seq_len(nrow(w))) {
    tempdf <- tempdf |>
      dplyr::mutate(
        !!rlang::sym(what_col_names[i]) :=
          purrr::map_int(
            .data$difVal,
            ~ sum(abs(.x - w$delta[i]) < 0.01)
          )
      )
  }

  tempdf <- tempdf |>
    dplyr::mutate(
      numberOfPeaks = lengths(.data$mz),
      highestValue = do.call(
        pmax,
        c(
          dplyr::pick(dplyr::all_of(what_col_names)),
          na.rm = TRUE
        )
      ),
      charge = purrr::pmap_chr(
        dplyr::pick(dplyr::all_of(what_col_names)),
        function(...) {
          vals <- c(...)
          names(vals)[which.max(vals)]
        }
      ),
      charge = as.integer(.data$charge),
      percentageMatched = .data$highestValue / .data$numberOfPeaks
    ) |>
    dplyr::distinct(.data$id, .data$slice,
                    .data$charge, .data$percentageMatched)

  spectrum <- spectrum |>
    dplyr::left_join(tempdf, dplyr::join_by("id", "slice"))

  spectrum <- spectrum |>
    dplyr::filter(.data$percentageMatched >= min_charge_percentage_matched)

  spectrum
}

estimate_intact_mass <- function(sliced_df) {
  averagemz <- (min(sliced_df$mz) + max(sliced_df$mz)) / 2

  averagemz * sliced_df$charge[1]
}
