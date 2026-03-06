extract_isotopes <- function(full_mgf, aa_seq, protein_mass,
                             threshold, output_folder,
                             score = "cosine", number_of_cores = 1,
                             isotope_peaks_included = 3, cutoff = 0.05,
                             ppm_tolerance = 20, cosine_score_correction = 1,
                             minimal_isotope_sequence = c(-2, -1, 0, 1, 2)) {

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
      threshold
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
      threshold
    )
  }

  mergepdfs(paste0(path_to_export, "/matched_spectra.pdf"),
            paste0(path_to_export, "/separated_matched_spectra/"))
  mergepdfs(paste0(path_to_export, "/unmatched_spectra.pdf"),
            paste0(path_to_export, "/separated_unmatched_spectra/"))

  founddf <- do.call(rbind, lapply(the_results, `[[`, "founddf"))
  nonefounddf <- do.call(rbind, lapply(the_results, `[[`, "nonefounddf"))

  write.csv(founddf, paste0(path_to_export, "/matched_peaks.csv"),
            row.names = FALSE)

  write.csv(nonefounddf, paste0(path_to_export, "/unmatched_peaks.csv"),
            row.names = FALSE)
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
                            threshold) {
  cl <- parallel::makeCluster(number_of_cores)
  doParallel::registerDoParallel(cl)

  isotopes <- get("isotopes", envir = asNamespace("enviPat"))

  results <- {
    foreach::foreach(
      frame = all_frames,
      .packages = c("dplyr", "ggplot2", "tidyr", "enviPat"),
      .export = c(
        "get_isotope_template",
        "isotopefit_opt",
        "substract_isotope_fit_opt"
      )
    ) %dopar% {
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
        while (inaction) {
          #Calculate correlation values
          temp_spectrum$cor <- sapply(temp_spectrum$mz, function(x) {
            isotopefit_opt(
              spectrum = temp_spectrum,
              isotope_template = isotope_template,
              poi = x,
              ppm_tolerance = ppm_tolerance,
              cosine_score_correction = cosine_score_correction,
              minimal_isotope_sequence = minimal_isotope_sequence,
              score = score
            )
          })

          #If the maximum score is higher than the set threshold
          #Replot and substract the isotope fit
          #Return the spectrum, p1, lineToAdd
          if (max(temp_spectrum$cor) >= threshold) {
            poi <- subset(temp_spectrum, cor == max(temp_spectrum$cor))$mz

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
              slice_counter = slice_counter
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

            if (!isotopeidentified) {
              if (nrow(local_none_found_df) == 0) {
                local_none_found_df <-
                  subset(
                    temp_spectrum,
                    cor == max(temp_spectrum$cor)
                  )[1, ]
              }else {
                local_none_found_df <-
                  rbind(
                    local_none_found_df,
                    subset(temp_spectrum, cor == max(temp_spectrum$cor))[1, ]
                  )
              }
            }

            p1 <-  ggplot(temp_spectrum,
                          aes(x = .data$mz, y = .data$intensity)) +
              geom_bar(
                data = temp_spectrum,
                aes(
                  x = .data$mz,
                  y = .data$centroided_intensity
                ),
                stat = "identity",
                width = 0.05,
                fill = "#444444"
              ) +
              labs(
                title = paste0("Slice ", i, " ; ", max(temp_spectrum$cor)),
                x = "m/z",
                y = "Intensity"
              ) +
              theme_classic() +
              theme(panel.grid = element_blank()) +
              scale_y_continuous(expand = c(0, 0),
                                 limits =
                                   c(0,
                                     max(temp_spectrum$intensity) * 1.05))

            invisible(suppressMessages(ggsave(
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
                           threshold) {
  message("Use  the number_of_cores argument to use multiple cores.
          It can take a while when using only one core.")

  isotopes <- get("isotopes", envir = asNamespace("enviPat"))

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

        while (inaction) {
          #Calculate correlation values
          temp_spectrum$cor <- sapply(temp_spectrum$mz, function(x) {
            isotopefit_opt_cpp(
              spectrum = temp_spectrum,
              isotope_template = isotope_template,
              poi = x,
              ppm_tolerance = ppm_tolerance,
              cosine_score_correction = cosine_score_correction,
              minimal_isotope_sequence = minimal_isotope_sequence,
              score = score
            )
          })

          #If the maximum score is higher than the set threshold
          #Replot and substract the isotope fit
          #Return the spectrum, p1, lineToAdd
          if (max(temp_spectrum$cor) >= threshold) {
            poi <- subset(temp_spectrum, cor == max(temp_spectrum$cor))$mz

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
              slice_counter = slice_counter
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

            if (!isotopeidentified) {
              if (nrow(local_none_found_df) == 0) {
                local_none_found_df <-
                  subset(
                    temp_spectrum,
                    cor == max(temp_spectrum$cor)
                  )[1, ]
              }else {
                local_none_found_df <-
                  rbind(
                    local_none_found_df,
                    subset(temp_spectrum, cor == max(temp_spectrum$cor))[1, ]
                  )
              }
            }

            p1 <-  ggplot(temp_spectrum,
                          aes(x = .data$mz, y = .data$intensity)) +
              geom_bar(
                data = temp_spectrum,
                aes(
                  x = .data$mz,
                  y = .data$centroided_intensity
                ),
                stat = "identity",
                width = 0.05,
                fill = "#444444"
              ) +
              labs(
                title = paste0("Slice ", i, " ; ", max(temp_spectrum$cor)),
                x = "m/z",
                y = "Intensity"
              ) +
              theme_classic() +
              theme(panel.grid = element_blank()) +
              scale_y_continuous(expand = c(0, 0),
                                 limits =
                                   c(0,
                                     max(temp_spectrum$intensity) * 1.05))

            invisible(suppressMessages(ggsave(
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
                                      cosine,
                                      cosine_score_correction,
                                      slice_counter) {
  # --- Step 1: Align template to poi ---
  max_intensity_row <- isotope_template |>
    dplyr::filter(.data$intensity == max(.data$intensity)) |>
    dplyr::slice(1)
  difmz <- poi - max_intensity_row$mz
  isotope_template$mz <- isotope_template$mz + difmz

  # Scale intensities
  max_spectrum_intensity <-
    max(spectrum$centroided_intensity[spectrum$mz == poi], na.rm = TRUE)
  difint <- max_spectrum_intensity /
    isotope_template$intensity[isotope_template$percent == 100][1]
  isotope_template$intensity <- isotope_template$intensity * difint

  # --- Step 2: Vectorized nearest m/z search ---
  # This replaces sapply(...) with fast vectorized approximation
  spectrum_mz <- spectrum$mz
  template_mz <- isotope_template$mz

  # For each template mz, find closest spectrum mz
  idx <- tryCatch({
    findInterval(template_mz, spectrum_mz)
  }, error = function(e) {
    message(
      "Error in findInterval with template_mz length: ",
      length(template_mz),
      ", spectrum_mz length: ",
      length(spectrum_mz),
      " -> ",
      e$message
    )
    NA
  })
  idx[idx == 0] <- 1
  idx[idx > length(spectrum_mz)] <- length(spectrum_mz)

  # Choose the closer of idx and idx+1
  closest <- ifelse(
    idx < length(spectrum_mz) &
      abs(spectrum_mz[idx + 1] - template_mz) <
        abs(spectrum_mz[idx] - template_mz),
    idx + 1,
    idx
  )

  isotope_template$spectrum_mz <- spectrum_mz[closest]
  isotope_template$centroided_intensity <-
    spectrum$centroided_intensity[closest]
  isotope_template$ppmError <- (isotope_template$mz -
                                  isotope_template$spectrum_mz) /
    isotope_template$spectrum_mz * 1e6

  # --- Step 3: Filter by ppm_tolerance ---
  isotope_template <-
    isotope_template[abs(isotope_template$ppmError) <= ppm_tolerance, ]

  # --- Step 4: Trim missing isotopeSeq without loops ---
  existing_seq <- isotope_template$isotope_seq

  # Trim front
  trim_front <- min(existing_seq):max(existing_seq)
  first_missing <- trim_front[!trim_front %in% existing_seq][1]
  if (!is.na(first_missing)) {
    isotope_template <-
      isotope_template[isotope_template$isotope_seq < first_missing, ]
  }
  # Trim back
  trim_back <- max(existing_seq):min(existing_seq)
  last_missing <- trim_back[!trim_back %in% existing_seq][1]

  if (!is.na(last_missing)) {
    isotope_template <-
      isotope_template[isotope_template$isotope_seq > last_missing, ]
  }

  # --- Step 5: Compute fit ---
  fit <- 0
  if (length(isotope_template$intensity) >= 5) {
    if (score == "pearson") {
      fit <- cor(isotope_template$intensity,
                 isotope_template$centroided_intensity)
    } else if (score == "cosine") {
      template_intensity <- isotope_template$intensity
      isotope_intensity <- isotope_template$centroided_intensity
      offset <- min(template_intensity, isotope_intensity, na.rm = TRUE) *
        cosine_score_correction
      template_intensity <- template_intensity - offset
      isotope_intensity <- isotope_intensity - offset
      fit <- lsa::cosine(template_intensity, isotope_intensity)[[1]]
    }
  } else {
    stop("score type unkown")
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

  invisible(suppressMessages(ggsave(
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
    dplyr::right_join(spectrum, join_by("mz")) |>
    dplyr::arrange(.data$mz)

  spectrum$intensity_isotope[is.na(spectrum$intensity_isotope)] <- 0

  line_to_add <- subset(spectrum, cor == max(spectrum$cor))
  line_to_add$averageMass <-
    IsotopeExtractor:::calculate_average_mass(isotope_template, spectrum$charge[1])

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
