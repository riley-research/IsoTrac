#' Extract Isotope Fit Patterns from Mass Spectrometry Data
#'
#' Determines intact masses from the dataset imported using the import_spectra function,
#' by iteratively matching molecular isotope patterns. Processing can be
#' executed sequentially or in parallel using a \code{future} backend. Matched and unmatched
#' spectra profiles are written out as binned CSV tables and visual PDF spectrum plots.
#'
#' @param full_mgf \code{data.frame}. The data imported using the \code{\link{import_spectra}}
#' function.
#' @param aa_seq \code{character}. Primary amino acid sequence of the target protein backbone
#'   used for matching theoretical distribution envelopes.
#' @param protein_mass \code{numeric}. Average mass of the bare protein backbone.
#' @param threshold \code{numeric}. Minimum similarity fit score required to accept and
#'   subtract an isotopic envelope from the spectrum.
#' @param output_folder \code{character}. Path to the parent directory where the timestamped
#'   results folder and diagnostic files will be created.
#' @param score \code{character}. Vector choice containing the similarity metric to optimize.
#'   Must be one of \code{"pearson"}, \code{"cosine"}, or \code{"NRMSE"}.
#' @param number_of_cores \code{integer}. Number of parallel workers to allocate for the
#'   \code{future} multisession execution. Default is \code{1}.
#' @param isotope_peaks_included \code{integer}. Total number of isotopic peaks trailing the
#'   most intense centroided peak tracked within the calculation frame. Default is \code{3}.
#' @param cutoff \code{numeric}. General filtering score threshold limit constraint. Default is \code{0.05}.
#' @param ppm_tolerance \code{numeric}. Maximum allowable mass error assignment threshold in
#'   parts-per-million (\code{PPM}) between experimental and template centroids. Default is \code{20}.
#' @param cosine_score_correction \code{numeric}. Intensity multiplier offset constant applied
#'   during manual vector translation checks. Default is \code{1}.
#' @param minimal_isotope_sequence \code{numeric vector}. Array of required continuous isotope
#'   peaks (e.g., \code{c(-2, -1, 0, 1, 2)}).
#' @param nrmse_cutoff \code{numeric}. Normalized Root Mean Square Error limit cutoff filter.
#'   Default is \code{0}.
#'
#' @return \code{list}. A structured list containing two data frames:
#'   \itemize{
#'     \item \code{founddf}: Data frame containing peaks matching the template threshold bounds.
#'     \item \code{nonefounddf}: Data frame containing spectrum features failing template criteria.
#'   }
#'
#' @seealso \code{\link[future:plan]{future::plan}}, \code{\link[ggplot2:ggsave]{ggplot2::ggsave}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' extraction_results <- extract_isotopes(
#'   full_mgf = experimental_data,
#'   aa_seq = "PEPTIDEK",
#'   protein_mass = 910.45,
#'   threshold = 0.85,
#'   output_folder = "./results/run_1",
#'   score = "cosine",
#'   number_of_cores = 4,
#'   ppm_tolerance = 15
#' )
#' }
extract_isotopes <- function(full_mgf, aa_seq, protein_mass,
                             threshold, output_folder,
                             score = c("pearson", "cosine", "NRMSE"), number_of_cores = 1,
                             isotope_peaks_included = 3, cutoff = 0.05,
                             ppm_tolerance = 20, cosine_score_correction = 1,
                             minimal_isotope_sequence = c(-2, -1, 0, 1, 2),
                             nrmse_cutoff = 0) {

  path_to_export <- create_output_folder(output_folder)

  full_mgf <- full_mgf |>
    dplyr::filter(!is.na(.data$centroided_intensity))

  all_frames <- split(full_mgf, list(full_mgf$id, full_mgf$slice), drop = TRUE)

  if (number_of_cores != 1) {
    the_results <- runner_parallel(
      all_frames,
      aa_seq,
      protein_mass,
      path_to_export,
      number_of_cores,
      isotope_peaks_included,
      ppm_tolerance,
      minimal_isotope_sequence,
      score,
      cosine_score_correction,
      threshold,
      nrmse_cutoff
    )
  }else {
    the_results <- runner_sequential(
      all_frames,
      aa_seq,
      protein_mass,
      path_to_export,
      number_of_cores,
      isotope_peaks_included,
      ppm_tolerance,
      minimal_isotope_sequence,
      score,
      cosine_score_correction,
      threshold,
      nrmse_cutoff
    )
  }

  mergepdfs(paste0(path_to_export, "/matched_spectra.pdf"),
            paste0(path_to_export, "/separated_matched_spectra/"))
  mergepdfs(paste0(path_to_export, "/unmatched_spectra.pdf"),
            paste0(path_to_export, "/separated_unmatched_spectra/"))

  founddf <- do.call(rbind, lapply(the_results, `[[`, "founddf"))
  nonefounddf <- do.call(rbind, lapply(the_results, `[[`, "nonefounddf"))

  utils::write.csv(founddf, paste0(path_to_export, "/matched_peaks.csv"),
            row.names = FALSE)

  utils::write.csv(nonefounddf, paste0(path_to_export, "/unmatched_peaks.csv"),
            row.names = FALSE)

  return(list(founddf = founddf, nonefounddf = nonefounddf))

}

create_output_folder <- function(output_folder) {
  folder_name <- strsplit(as.character(Sys.time()), "\\.")[[1]][1]
  folder_name <- gsub(" ", "_", folder_name)
  folder_name <- gsub(":", "-", folder_name)
  path_to_export <- paste0(output_folder, "/", folder_name)

  dir.create(path_to_export, showWarnings = FALSE)
  dir.create(paste0(path_to_export, "/separated_matched_spectra/"))
  dir.create(paste0(path_to_export, "/separated_unmatched_spectra/"))

  path_to_export
}

runner_parallel <- function(all_frames,
                            aa_seq,
                            protein_mass,
                            path_to_export,
                            number_of_cores,
                            isotope_peaks_included,
                            ppm_tolerance,
                            minimal_isotope_sequence,
                            score,
                            cosine_score_correction,
                            threshold,
                            nrmse_cutoff) {

  # 1. Setup future plan
  future::plan(future::multisession, workers = number_of_cores)
  # Ensure cleanup of workers when function exits
  on.exit(future::plan(future::sequential))

  # Fetch isotopes reference once for export
  isotopes_ref <- get("isotopes", envir = asNamespace("enviPat"))

  # 2. Parallel execution using future_lapply
  results <- future.apply::future_lapply(all_frames, function(frame) {

    spectrum <- frame
    spectrum_id <- spectrum$id[1]

    local_found_df <- data.frame()
    local_none_found_df <- data.frame()

    for (i in unique(spectrum$slice)) {
      slice_counter <- 1
      inaction <- TRUE
      isotopeidentified <- FALSE

      temp_spectrum <- spectrum |>
        dplyr::filter(.data$slice == i)

      spectrum_slice <- i

      isotope_template <- get_isotope_template(
        temp_spectrum,
        aa_seq,
        protein_mass,
        isotope_peaks_included,
        isotopes = isotopes_ref
      )

      isotope_template_average_mass <- get_isotope_template(temp_spectrum,
                                                            aa_seq,
                                                            protein_mass,
                                                            1000000,
                                                            isotopes = isotopes_ref)

      while (inaction) {
        score_matrix <- vapply(temp_spectrum$mz, function(x) {
          isotopefit_opt_cpp(
            spectrum = temp_spectrum,
            isotope_template = isotope_template,
            poi = x,
            ppm_tolerance = ppm_tolerance,
            cosine_score_correction = cosine_score_correction,
            minimal_isotope_sequence = minimal_isotope_sequence,
            score = score
          )
        }, numeric(3))

        temp_spectrum[, c("pearson", "cosine", "NRMSE")] <- t(score_matrix)
        if (score == "pearson") {
          temp_spectrum$cor <- temp_spectrum$pearson
        }else if (score == "cosine") {
          temp_spectrum$cor <- temp_spectrum$cosine
        } else if (score == "NRMSE") {
          temp_spectrum$cor <- temp_spectrum$NRMSE
        }

        best_nrmse <- temp_spectrum$NRMSE[which.max(temp_spectrum$cor)]

        #If the maximum score is higher than the set threshold
        #Replot and substract the isotope fit
        #Return the spectrum, p1, lineToAdd
        if ((score %in% c("pearson", "cosine") & max(temp_spectrum$cor) >= threshold &
             best_nrmse <= nrmse_cutoff) |
            (score == "NRMSE" & min(temp_spectrum$cor) <= threshold)) {
          if (score == "NRMSE") {
            poi <- temp_spectrum$mz[temp_spectrum$cor == min(temp_spectrum$cor, na.rm = TRUE)]
          } else {
            poi <- temp_spectrum$mz[temp_spectrum$cor == max(temp_spectrum$cor, na.rm = TRUE)]
          }

          temp <- substract_isotope_fit_opt(
            spectrum = temp_spectrum,
            isotope_template = isotope_template,
            poi = poi,
            ppm_tolerance = ppm_tolerance,
            slice = i,
            path_to_export = path_to_export,
            id = spectrum_id,
            score = score,
            cosine_score_correction =
              cosine_score_correction,
            slice_counter = slice_counter,
            isotope_template_am = isotope_template_average_mass
          )

          slice_counter <- slice_counter + 1
          temp_spectrum <- temp[[1]]

          local_found_df <- rbind(local_found_df, temp[[2]])

          isotopeidentified <- TRUE
          if (nrow(temp_spectrum) == 0) inaction <- FALSE

        } else {
          inaction = FALSE

          if (!isotopeidentified & score == "NRMSE") {
            if (nrow(local_none_found_df) == 0) {
              local_none_found_df <- temp_spectrum[which(temp_spectrum$cor == min(temp_spectrum$cor, na.rm = TRUE))[1], ]
            } else {
              local_none_found_df <- rbind(
                local_none_found_df,
                temp_spectrum[which(temp_spectrum$cor == min(temp_spectrum$cor, na.rm = TRUE))[1], ]
              )
            }
          } else {
            if (nrow(local_none_found_df) == 0) {
              local_none_found_df <- temp_spectrum[which(temp_spectrum$cor == max(temp_spectrum$cor, na.rm = TRUE))[1], ]
            } else {
              local_none_found_df <- rbind(
                local_none_found_df,
                temp_spectrum[which(temp_spectrum$cor == max(temp_spectrum$cor, na.rm = TRUE))[1], ]
              )
            }
          }

          p1 <-  ggplot2::ggplot(temp_spectrum,
                                 ggplot2::aes(x = .data$mz, y = .data$intensity)) +
            ggplot2::geom_bar(
              data = temp_spectrum,
              ggplot2::aes(
                x = .data$mz,
                y = .data$centroided_intensity
              ),
              stat = "identity",
              width = 0.05,
              fill = "#444444"
            ) +
            ggplot2::labs(
              title = paste0(temp_spectrum$id[1], "; Slice ", i, " ; ", max(temp_spectrum$cor)),
              x = "m/z",
              y = "Intensity"
            ) +
            ggplot2::theme_classic() +
            ggplot2::theme(panel.grid = ggplot2::element_blank()) +
            ggplot2::scale_y_continuous(expand = c(0, 0),
                               limits =
                                 c(0,
                                   max(temp_spectrum$intensity) * 1.05))

          plot_path <- file.path(path_to_export, "separated_unmatched_spectra")
          if(!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

          ggplot2::ggsave(
            filename = file.path(plot_path, paste0(spectrum_id, "_", spectrum_slice, ".pdf")),
            plot = p1, width = 6, height = 6, device = "pdf"
          )
        }
      }
    }

    # Return local results to be combined by future_lapply
    return(list(founddf = local_found_df, nonefounddf = local_none_found_df))

  },
  future.seed = TRUE,
  future.globals = c(
    "aa_seq",
    "protein_mass",
    "isotopes_ref",
    "ppm_tolerance",
    "cosine_score_correction",
    "minimal_isotope_sequence",
    "score",
    "threshold",
    "nrmse_cutoff",
    "path_to_export",
    "substract_isotope_fit_opt",
    "get_isotope_template",
    "isotopefit_opt_cpp"
  ),
  future.packages = c("dplyr", "ggplot2", "tidyr", "enviPat", "IsotopeExtractor", "Rcpp"))

  return(results)
}

runner_sequential <- function(all_frames,
                           aa_seq,
                           protein_mass,
                           path_to_export,
                           number_of_cores,
                           isotope_peaks_included,
                           ppm_tolerance,
                           minimal_isotope_sequence,
                           score,
                           cosine_score_correction,
                           threshold,
                           nrmse_cutoff) {
  message("Use  the number_of_cores argument to use multiple cores.
          It can take a while when using only one core.")

  isotopes <- get("isotopes", envir = asNamespace("enviPat"))

  frame <- NULL

  results <- {
    foreach::foreach(
      frame = all_frames,
      .packages = c("dplyr", "ggplot2", "tidyr", "enviPat"),
      .export = c(
        "get_isotope_template",
        "isotopefit_opt_cpp",
        "substract_isotope_fit_opt"
      )
    ) %do% {
      spectrum <- frame
      spectrum_id <- spectrum$id[1]

      # Local storage for each thread
      local_found_df <- data.frame()
      local_none_found_df <- data.frame()

      for (i in unique(spectrum$slice)) {
        slice_counter = 1
        inaction <- TRUE
        isotopeidentified <- FALSE
        temp_spectrum <- spectrum |>
          dplyr::filter(.data$slice == i)

        spectrum_slice <- i
        isotope_template <- get_isotope_template(temp_spectrum,
                                                 aa_seq,
                                                 protein_mass,
                                                 isotope_peaks_included,
                                                 isotopes = isotopes)

        isotope_template_average_mass <- get_isotope_template(temp_spectrum,
                                                 aa_seq,
                                                 protein_mass,
                                                 1000000,
                                                 isotopes = isotopes)

        while (inaction) {
          score_matrix <- vapply(temp_spectrum$mz, function(x) {
            isotopefit_opt_cpp(
              spectrum = temp_spectrum,
              isotope_template = isotope_template,
              poi = x,
              ppm_tolerance = ppm_tolerance,
              cosine_score_correction = cosine_score_correction,
              minimal_isotope_sequence = minimal_isotope_sequence,
              score = score
            )
          }, numeric(3))

          temp_spectrum[, c("pearson", "cosine", "NRMSE")] <- t(score_matrix)
          if (score == "pearson") {
            temp_spectrum$cor <- temp_spectrum$pearson
          }else if (score == "cosine") {
            temp_spectrum$cor <- temp_spectrum$cosine
          } else if (score == "NRMSE") {
            temp_spectrum$cor <- temp_spectrum$NRMSE
          }

          best_nrmse <- temp_spectrum$NRMSE[which.max(temp_spectrum$cor)]

          #If the maximum score is higher than the set threshold
          #Replot and substract the isotope fit
          #Return the spectrum, p1, lineToAdd
          if ((score %in% c("pearson", "cosine") & max(temp_spectrum$cor) >= threshold &
               best_nrmse <= nrmse_cutoff) |
              (score == "NRMSE" & min(temp_spectrum$cor) <= threshold)) {
            if (score == "NRMSE") {
              poi <- temp_spectrum$mz[temp_spectrum$cor == min(temp_spectrum$cor, na.rm = TRUE)]
            } else {
              poi <- temp_spectrum$mz[temp_spectrum$cor == max(temp_spectrum$cor, na.rm = TRUE)]
            }

            temp <- substract_isotope_fit_opt(
              spectrum = temp_spectrum,
              isotope_template = isotope_template,
              poi = poi,
              ppm_tolerance = ppm_tolerance,
              slice = i,
              path_to_export = path_to_export,
              id = spectrum_id,
              score = score,
              cosine_score_correction =
                cosine_score_correction,
              slice_counter = slice_counter,
              isotope_template_am = isotope_template_average_mass
            )

            slice_counter = slice_counter + 1

            temp_spectrum <- temp[[1]]

            if (nrow(local_found_df) == 0) {
              local_found_df <- temp[[2]]
            } else {
              local_found_df <- rbind(local_found_df, temp[[2]])
            }

            isotopeidentified <- TRUE
            if (nrow(temp_spectrum) == 0) {
              inaction = FALSE
            }

          }else {
            inaction = FALSE

            if (!isotopeidentified & score == "NRMSE") {
              if (nrow(local_none_found_df) == 0) {
                local_none_found_df <- temp_spectrum[which(temp_spectrum$cor == min(temp_spectrum$cor, na.rm = TRUE))[1], ]
              } else {
                local_none_found_df <- rbind(
                  local_none_found_df,
                  temp_spectrum[which(temp_spectrum$cor == min(temp_spectrum$cor, na.rm = TRUE))[1], ]
                )
              }
            } else {
              if (nrow(local_none_found_df) == 0) {
                local_none_found_df <- temp_spectrum[which(temp_spectrum$cor == max(temp_spectrum$cor, na.rm = TRUE))[1], ]
              } else {
                local_none_found_df <- rbind(
                  local_none_found_df,
                  temp_spectrum[which(temp_spectrum$cor == max(temp_spectrum$cor, na.rm = TRUE))[1], ]
                )
              }
            }

            p1 <-  ggplot2::ggplot(temp_spectrum,
                                   ggplot2::aes(x = .data$mz, y = .data$intensity)) +
              ggplot2::geom_bar(
                data = temp_spectrum,
                ggplot2::aes(
                  x = .data$mz,
                  y = .data$centroided_intensity
                ),
                stat = "identity",
                width = 0.05,
                fill = "#444444"
              ) +
              ggplot2::labs(
                title = paste0(temp_spectrum$id[1], "; Slice ", i, " ; ", max(temp_spectrum$cor)),
                x = "m/z",
                y = "Intensity"
              ) +
              ggplot2::theme_classic() +
              ggplot2::theme(panel.grid = ggplot2::element_blank()) +
              ggplot2::scale_y_continuous(expand = c(0, 0),
                                 limits =
                                   c(0,
                                     max(temp_spectrum$intensity) * 1.05))

            invisible(suppressMessages(ggplot2::ggsave(
              paste0(
                path_to_export,
                "/separated_unmatched_spectra/",
                spectrum_id,
                "_",
                spectrum_slice,
                ".pdf"
              ),
              plot = p1,
              width = 6,
              height = 6
            )))
          }
        }
      }

      # Return results from each thread
      list(founddf = local_found_df, nonefounddf = local_none_found_df)
    }
  }

  results
}

substract_isotope_fit_opt <- function(spectrum,
                                          isotope_template,
                                          poi,
                                          slice,
                                          ppm_tolerance,
                                          path_to_export,
                                          id,
                                          score,
                                          #cosine,
                                          cosine_score_correction,
                                          slice_counter,
                                          isotope_template_am) {
  # --- Step 1: Align template and Scale (Matches C++ Step 1) ---
  # Find 100% peak for scaling (difint)
  #idx_100 <- which.min(abs(isotope_template$percent - 100))
  idx_100 <- which.max(isotope_template$intensity)
  intensity_at_100pct <- isotope_template$intensity[idx_100]
  mz_at_max_temp <- isotope_template$mz[idx_100]

  # Find max spectrum intensity exactly at POI
  # C++ uses: std::abs(spec_mz[i] - poi) < 1e-6
  max_spectrum_intensity <- spectrum$centroided_intensity[which.min(abs(spectrum$mz - poi))]

  difint <- max_spectrum_intensity / intensity_at_100pct
  difmz <- poi - mz_at_max_temp

  # Update templates
  isotope_template$mz <- isotope_template$mz + difmz
  isotope_template$intensity <- isotope_template$intensity * difint
  isotope_template_am$mz <- isotope_template_am$mz + difmz

  # --- Step 2: Nearest m/z search (Matches C++ std::lower_bound logic) ---
  # Vectorized approach to find the absolute closest peak
  spec_mz <- spectrum$mz
  temp_mz <- isotope_template$mz

  # findInterval finds index i such that spec_mz[i] <= temp_mz < spec_mz[i+1]
  idx <- findInterval(temp_mz, spec_mz)

  # We need to check both idx and idx + 1 to find the absolute closest
  idx_plus_1 <- pmin(idx + 1, length(spec_mz))
  idx_safe <- pmax(idx, 1)

  dist_idx <- abs(spec_mz[idx_safe] - temp_mz)
  dist_idx_plus_1 <- abs(spec_mz[idx_plus_1] - temp_mz)

  closest <- ifelse(dist_idx_plus_1 < dist_idx, idx_plus_1, idx_safe)

  isotope_template$spectrum_mz <- spec_mz[closest]
  isotope_template$centroided_intensity <- spectrum$centroided_intensity[closest]

  # PPM calculation using experimental mz as denominator (Matches C++)
  isotope_template$ppmError <- (isotope_template$mz - isotope_template$spectrum_mz) /
    isotope_template$spectrum_mz * 1e6

  # --- Step 3: Filter by ppm_tolerance ---
  isotope_template <- isotope_template[abs(isotope_template$ppmError) <= ppm_tolerance, ]

  # --- Step 4: Validate and Contiguous Trim (Matches C++ Step 4) ---
  existing_seq <- isotope_template$isotope_seq
  # if (!all(minimal_isotope_sequence %in% existing_seq)) return(0)

  # C++ logic: Find 0-isotope and expand outwards until a gap is found
  zero_idx <- which(isotope_template$isotope_seq == 0)
  if (length(zero_idx) == 0) return(0)

  # Replicate the contiguous check logic
  # Find sequence gaps relative to 0
  seqs <- sort(existing_seq)
  # Forward from 0
  fwd <- seqs[seqs >= 0]
  if(length(fwd) > 1) {
    break_fwd <- which(diff(fwd) > 1)
    if(length(break_fwd) > 0) fwd <- fwd[1:break_fwd[1]]
  }
  # Backward from 0
  bwd <- sort(seqs[seqs <= 0], decreasing = TRUE)
  if(length(bwd) > 1) {
    break_bwd <- which(diff(bwd) < -1)
    if(length(break_bwd) > 0) bwd <- bwd[1:break_bwd[1]]
  }

  keep_seq <- c(fwd, bwd)
  isotope_template <- isotope_template[isotope_template$isotope_seq %in% keep_seq, ]

  # --- Step 5: Compute fit (Matches C++ Step 5) ---
  fit <- 0
  final_temp <- isotope_template$intensity
  final_spec <- isotope_template$centroided_intensity

  if (length(final_temp) >= 5) {
    if (score == "pearson") {
      fit <- stats::cor(final_temp, final_spec)
    } else if (score == "cosine") {
      # Use manual manual calculation to match C++ exactly
      min_val <- min(c(final_temp, final_spec))
      offset <- min_val * cosine_score_correction

      a <- final_temp - offset
      b <- final_spec - offset

      # Explicit dot product / magnitude calculation
      num <- sum(a * b)
      den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
      fit <- if(den > 0) num / den else 0
    }else if (score == "NRMSE") {
      min_val <- min(c(final_temp, final_spec))
      offset <- min_val * cosine_score_correction

      a <- final_temp - offset
      b <- final_spec - offset

      diff_sq <- (a - b)^2
      sum_sq_diff <- sum(diff_sq)

      rmse_val <- sqrt(sum_sq_diff / length(final_temp))

      peak_range <- max(a) - min(b)

      if (peak_range > 1e-12) {
        fit <- (rmse_val / peak_range) * 100.0
      } else {
        fit <- 999.0
      }
    }
  }

  p1 <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = spectrum,
      ggplot2::aes(x = .data$mz, y = .data$centroided_intensity),
      stat = "identity",
      fill = "#808080",
      alpha = 1,
      width = 0.05,
      position = "identity"
    ) +
    ggplot2::geom_bar(
      data = isotope_template,
      ggplot2::aes(x = .data$mz, y = .data$intensity),
      stat = "identity",
      fill = "#840032",
      alpha = 0.5,
      width = 0.05,
      position = "identity"
    ) +
    ggplot2::labs(title =
                    paste0(
                      id, "; Slice ", slice, "-", slice_counter, " ; ", fit),
                  x = "m/z",
                  y = "Intensity") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                limits = c(0, max(spectrum$intensity) *
                                             1.05))

  invisible(suppressMessages(ggplot2::ggsave(
    paste0(
      path_to_export,
      "/separated_matched_spectra/",
      id,
      "_",
      slice,
      "_",
      slice_counter,
      ".pdf"
    ),
    plot = p1,
    width = 5,
    height = 5
  )))

  spectrum <- isotope_template |>
    dplyr::select("spectrum_mz", "intensity") |>
    dplyr::rename("mz" = "spectrum_mz", "intensity_isotope" = "intensity") |>
    dplyr::right_join(spectrum, dplyr::join_by("mz")) |>
    dplyr::arrange(.data$mz)

  spectrum$intensity_isotope[is.na(spectrum$intensity_isotope)] <- 0

  line_to_add <- spectrum[which(spectrum$cor == max(spectrum$cor, na.rm = TRUE)), ]
  line_to_add$averageMass <-
    calculate_average_mass(isotope_template_am, spectrum$charge[1])
  #####

  spectrum$centroided_intensity <-
    spectrum$centroided_intensity - spectrum$intensity_isotope

  spectrum <- spectrum |>
    dplyr::select(-c("intensity_isotope", "cor")) |>
    dplyr::filter(.data$centroided_intensity > 0)

  list(spectrum, line_to_add)
}

calculate_average_mass <- function(isotope_template, charge) {
  isotope_template$deconvMz <- sapply(isotope_template$mz, function(x) {
    (x * charge) - (charge * 1.00784)
  })

  isotope_template$mt <-
    isotope_template$deconvMz * (isotope_template$percent * 0.01)

  average_mass <-
    sum(isotope_template$mt) / sum(isotope_template$percent * 0.01)

  isotope_template <- isotope_template |>
    dplyr::select(-c("mt"))

  average_mass
}

mergepdfs <- function(outputfile, inputfolder, chunk_size = 10){
  # List and sort PDF files
  pdf_files <- sort(list.files(inputfolder, pattern = "\\.pdf$", full.names = TRUE))

  # Step 1: Merge in small batches
  temp_files <- c()
  for(i in seq(1, length(pdf_files), by = chunk_size)) {
    chunk <- pdf_files[i:min(i + chunk_size - 1, length(pdf_files))]
    temp_out <- tempfile(fileext = ".pdf")
    qpdf::pdf_combine(input = chunk, output = temp_out)
    temp_files <- c(temp_files, temp_out)
  }

  # Step 2: Merge the temporary chunk PDFs into final output
  while(length(temp_files) > 1){
    new_temp <- c()
    for(i in seq(1, length(temp_files), by = chunk_size)){
      chunk <- temp_files[i:min(i + chunk_size - 1, length(temp_files))]
      temp_out <- tempfile(fileext = ".pdf")
      qpdf::pdf_combine(input = chunk, output = temp_out)
      new_temp <- c(new_temp, temp_out)
    }
    temp_files <- new_temp
  }

  # Rename final single PDF
  file.rename(temp_files, outputfile)

  message("Merged PDF created at: ", outputfile)
}
