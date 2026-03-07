get_fdr_threshold <- function(full_mgf, full_fdr_mgf, aa_seq, protein_mass,
                              score = "cosine", number_of_cores = 1,
                              isotope_peaks_included = 3, cutoff = 0.05,
                              ppm_tolerance = 20, cosine_score_correction = 1,
                              minimal_isotope_sequence = c(-2, -1, 0, 1, 2),
                              show_plot = TRUE) {
  full_mgf <- full_mgf |>
    dplyr::filter(!is.na(.data$centroided_intensity))

  if (number_of_cores == 1) {
    message("Calculating the FDR is slow. Please consider setting the
            number_of_cores argument to parallelize the calculating.")
    message("Calculating target scores.")
    targets <- return_max_cor_scores_cpp(full_mgf, aa_seq, protein_mass,
                                     isotope_peaks_included, ppm_tolerance,
                                     cosine_score_correction,
                                     minimal_isotope_sequence, score)

    # write.csv(targets, "Z:/Tim/26-02-23_eqQ/EpCAM/EpCAM_DIA-PTCR_2800-3350_20260217172110/Analysis/classic.csv",
    #           row.names = FALSE)
    #
    # write.csv(targets2, "Z:/Tim/26-02-23_eqQ/EpCAM/EpCAM_DIA-PTCR_2800-3350_20260217172110/Analysis/cpp.csv",
    #           row.names = FALSE)

    message("Calculating decoy scores scores.")
    decoys <- return_max_cor_scores_cpp(full_fdr_mgf, aa_seq, protein_mass,
                                    isotope_peaks_included, ppm_tolerance,
                                    cosine_score_correction,
                                    minimal_isotope_sequence, score)
  }else {
    message("Calculating target scores.")
    targets <- return_max_cor_scores_parallel_cpp(
      full_mgf,
      aa_seq,
      protein_mass,
      isotope_peaks_included,
      ppm_tolerance,
      cosine_score_correction,
      minimal_isotope_sequence,
      score,
      number_of_cores
    )

    message("Calculating decoy scores scores.")
    decoys <- return_max_cor_scores_parallel_cpp(
      full_fdr_mgf,
      aa_seq,
      protein_mass,
      isotope_peaks_included,
      ppm_tolerance,
      cosine_score_correction,
      minimal_isotope_sequence,
      score,
      number_of_cores
    )
  }

  cutoff_value <- compute_qvalues(targets, decoys, cutoff)


  if (show_plot) {
    plot_fdr_results(targets, decoys, cutoff_value)
  }

  cutoff_value
}

return_max_cor_scores_parallel_cpp <- function(spectrum,
                                               aa_seq,
                                               protein_mass,
                                               isotope_peaks_included,
                                               ppm_tolerance,
                                               cosine_score_correction,
                                               minimal_isotope_sequence,
                                               score,
                                               number_of_cores) {

  # 1. Setup future plan
  # multisession is the equivalent of a socket cluster
  future::plan(future::multisession, workers = number_of_cores)

  # Ensure cleanup of workers when function exits
  on.exit(future::plan(future::sequential))

  # 2. Pre-process data
  spectrum_list <- split(spectrum, list(spectrum$id, spectrum$slice), drop = TRUE)

  # Efficiency Fix: Calculate template ONCE if it's constant for all spectra
  # If it depends on temp_spectrum, move this back inside the loop
  isotopes <- get("isotopes", envir = asNamespace("enviPat"))

  # 3. Parallel execution using future_lapply
  # future_lapply handles export of globals automatically
  max_vec <- future.apply::future_lapply(spectrum_list, function(temp_spectrum) {

    if (nrow(temp_spectrum) == 0) return(0.0)

    isotope_template <- get_isotope_template(
      temp_spectrum, # Use the actual slice
      aa_seq,
      protein_mass,
      isotope_peaks_included,
      isotopes = isotopes
    )

    cor_vec <- vapply(temp_spectrum$mz, function(x) {
      IsotopeExtractor:::isotopefit_opt_cpp(
        spectrum = temp_spectrum,
        isotope_template = isotope_template,
        poi = x,
        ppm_tolerance = ppm_tolerance,
        cosine_score_correction = cosine_score_correction,
        minimal_isotope_sequence = minimal_isotope_sequence,
        score = score
      )
    }, numeric(1))

    res <- max(cor_vec, na.rm = TRUE)
    return(if(is.infinite(res)) 0.0 else res)

  }, future.seed = TRUE, future.packages = c("IsotopeExtractor", "enviPat", "Rcpp"))

  # Return as a combined vector
  return(unlist(max_vec))
}

return_max_cor_scores_parallel_cpp_old <- function(spectrum,
                                           aa_seq,
                                           protein_mass,
                                           isotope_peaks_included,
                                           ppm_tolerance,
                                           cosine_score_correction,
                                           minimal_isotope_sequence,
                                           score,
                                           number_of_cores) {

  cl <- parallel::makeCluster(number_of_cores, outfile = "Z:/Tim/26-02-23_eqQ/EpCAM/EpCAM_DIA-PTCR_2800-3350_20260217172110/Analysis/cpp.log")
  doParallel::registerDoParallel(cl)

  spectrum_list <- split(spectrum, list(spectrum$id, spectrum$slice),
                         drop = TRUE)

  isotopes <- get("isotopes", envir = asNamespace("enviPat"))

  max_vec <- foreach::foreach(
    temp_spectrum = spectrum_list,
    .packages = c("dplyr", "ggplot2", "tidyr", "enviPat", "IsotopeExtractor"),
    .export = c(
      "get_isotope_template",
      "isotopefit_opt_cpp",
      "amino_acid_composition"
    ),
    .combine = c
  ) %dopar% {
    isotope_template <- get_isotope_template(temp_spectrum,
                                             aa_seq,
                                             protein_mass,
                                             isotope_peaks_included,
                                             isotopes = isotopes)

    cor_vec <- vapply(temp_spectrum$mz, function(x) {
      isotopefit_opt_cpp(
        spectrum = temp_spectrum,
        isotope_template = isotope_template,
        poi = x,
        ppm_tolerance,
        cosine_score_correction,
        minimal_isotope_sequence,
        score
      )
    }, numeric(1))

    max(cor_vec, na.rm = TRUE)

  }

  parallel::stopCluster(cl)

  max_vec
}

return_max_cor_scores_parallel <- function(spectrum,
                                           aa_seq,
                                           protein_mass,
                                           isotope_peaks_included,
                                           ppm_tolerance,
                                           cosine_score_correction,
                                           minimal_isotope_sequence,
                                           score,
                                           number_of_cores) {

  cl <- parallel::makeCluster(number_of_cores)
  doParallel::registerDoParallel(cl)

  spectrum_list <- split(spectrum, list(spectrum$id, spectrum$slice),
                         drop = TRUE)

  isotopes <- get("isotopes", envir = asNamespace("enviPat"))

  max_vec <- foreach::foreach(
    temp_spectrum = spectrum_list,
    .packages = c("dplyr", "ggplot2", "tidyr", "enviPat", "IsotopeExtractor"),
    .export = c(
      "get_isotope_template",
      "isotopefit_opt",
      "amino_acid_composition"
    ),
    .combine = c
  ) %dopar% {
    isotope_template <- get_isotope_template(temp_spectrum,
                                             aa_seq,
                                             protein_mass,
                                             isotope_peaks_included,
                                             isotopes = isotopes)

    cor_vec <- vapply(temp_spectrum$mz, function(x) {
      isotopefit_opt(
        spectrum = temp_spectrum,
        isotope_template = isotope_template,
        poi = x,
        ppm_tolerance,
        cosine_score_correction,
        minimal_isotope_sequence,
        score
      )
    }, numeric(1))

    max(cor_vec, na.rm = TRUE)

  }

  parallel::stopCluster(cl)

  max_vec
}

return_max_cor_scores <- function(spectrum, aa_seq, protein_mass,
                                  isotope_peaks_included, ppm_tolerance,
                                  cosine_score_correction,
                                  minimal_isotope_sequence, score) {

  spectrum_list <- split(spectrum, list(spectrum$id, spectrum$slice),
                         drop = TRUE)

  max_vec <- vapply(spectrum_list, function(temp_spectrum) {

    isotope_template <- get_isotope_template(temp_spectrum,
                                             aa_seq, protein_mass,
                                             isotope_peaks_included)
    cor_vec <- vapply(temp_spectrum$mz, function(x) {
      isotopefit_opt(
        spectrum = temp_spectrum,
        isotope_template = isotope_template,
        poi = x,
        ppm_tolerance,
        cosine_score_correction,
        minimal_isotope_sequence,
        score
      )
    }, numeric(1))

    max(cor_vec, na.rm = TRUE)

  }, numeric(1))

  max_vec
}

return_max_cor_scores_cpp <- function(spectrum, aa_seq, protein_mass,
                                  isotope_peaks_included, ppm_tolerance,
                                  cosine_score_correction,
                                  minimal_isotope_sequence, score) {

  spectrum_list <- split(spectrum, list(spectrum$id, spectrum$slice),
                         drop = TRUE)

  max_vec <- vapply(spectrum_list, function(temp_spectrum) {

    isotope_template <- get_isotope_template(temp_spectrum,
                                             aa_seq, protein_mass,
                                             isotope_peaks_included)

    cor_vec <- vapply(temp_spectrum$mz, function(x) {
      isotopefit_opt_cpp(
        spectrum = temp_spectrum,
        isotope_template = isotope_template,
        poi = x,
        ppm_tolerance,
        cosine_score_correction,
        minimal_isotope_sequence,
        score
      )
    }, numeric(1))

    max(cor_vec, na.rm = TRUE)

  }, numeric(1))

  max_vec
}

get_isotope_template <- function(spectrum, aa_seq, protein_mass,
                                 isotope_peaks_included, isotopes = FALSE) {
  aa_comp <- IsotopeExtractor:::amino_acid_composition(sequence = aa_seq)

  if (spectrum$estimated_intact_mass[1] > protein_mass) {
    glycan_comp <- add_glycan(spectrum$estimated_intact_mass[1] - protein_mass)

    total_comp <- get_clean_total_comp(aa_comp +
                                         glycan_comp +
                                         c(C = 0, H = spectrum$charge[1],
                                           N = 0, O = 0, S = 0, P = 0))
  }else {
    total_comp <- get_clean_total_comp(aa_comp +
                                         c(C = 0, H = spectrum$charge[1],
                                           N = 0, O = 0, S = 0, P = 0))
  }

  if (identical(isotopes, FALSE)) {
    isotopes <- get("isotopes", envir = asNamespace("enviPat"))
  }

  pattern <- enviPat::isopattern(isotopes,
                                 total_comp,
                                 charge = spectrum$charge[1],
                                 threshold = 0.01,
                                 emass = 0.00054858,
                                 verbose = FALSE,
                                 plotit = FALSE,
                                 algo = 1)

  pattern <- as.data.frame(pattern[[1]][, 1:2])
  names(pattern) <- c("mz", "intensity")

  pattern <- pattern |>
    dplyr::arrange(dplyr::desc(.data$intensity))

  if(FALSE){
  cleandf <- data.frame()
  mz_seen <- c()
  for (i in seq_len(nrow(pattern))) {
    mzoi <- pattern$mz[i]

    if (mzoi %in% mz_seen) {
      next
    }

    tempdf <- pattern |>
      dplyr::filter(abs(.data$mz - mzoi) <= 0.01)
    mz_seen <- append(mz_seen, tempdf$mz)

    if (nrow(cleandf) == 0) {
      cleandf <- data.frame(mz = weighted.mean(tempdf$mz, tempdf$intensity),
                            intensity = sum(tempdf$intensity, na.rm = TRUE))
    }else {
      cleandf <- rbind(cleandf, data.frame(
        mz = weighted.mean(tempdf$mz, tempdf$intensity),
        intensity = sum(tempdf$intensity, na.rm = TRUE)
      ))
    }
  }
  }

  if(TRUE){
    results_list <- list()
    counter <- 1

    # Create a logical mask of 'unprocessed' peaks
    # This replaces the 'mz_seen' vector and is much faster
    unprocessed_mask <- rep(TRUE, nrow(pattern))

    while (any(unprocessed_mask)) {

      # 1. Pick the FIRST TRUE index (This matches the 'i' in your for loop)
      first_idx <- which(unprocessed_mask)[1]
      mzoi <- pattern$mz[first_idx]

      # 2. Identify all peaks in the CURRENT pattern within tolerance
      # We check against the WHOLE pattern to match your filter(abs(...) <= 0.01)
      in_window <- abs(pattern$mz - mzoi) <= 0.01

      # 3. Create tempdf from the window
      tempdf <- pattern[in_window, ]

      # 4. Store result
      results_list[[counter]] <- data.frame(
        mz = weighted.mean(tempdf$mz, tempdf$intensity, na.rm = TRUE),
        intensity = sum(tempdf$intensity, na.rm = TRUE)
      )

      # 5. Flag these peaks as 'processed' so they aren't picked as anchors again
      # This replicates 'mz_seen <- append(mz_seen, tempdf$mz)'
      unprocessed_mask[in_window] <- FALSE

      counter <- counter + 1
    }

    # Final result
    cleandf <- do.call(rbind, results_list)
  }

  # if(!identical(cleandf, cleandf1)){
  #   View(cleandf)
  #   View(cleandf1)
  #   stop()
  # }

  cleandf <- cleandf |>
    dplyr::arrange(.data$mz)

  cleandf <- cleandf[complete.cases(cleandf), ]

  cleandf$percent <- as.numeric(cleandf$intensity /
                                  max(cleandf$intensity) * 100)

  isotope_template <- clean_isotope_duplicates(cleandf) #Remove any duplicates

  for (j in seq_len(nrow(isotope_template))) {
    isotope_template$isotope_seq[j] <- j -
      which(isotope_template$percent == 100)[1]
  }

  if (isotope_peaks_included != 0) {
    isotope_template <- isotope_template |>
      dplyr::filter(abs(.data$isotope_seq) <= isotope_peaks_included)
  }

  isotope_template
}

amino_acid_composition <- function(sequence, iaa = TRUE) {

  # Initialize elemental counts for the sequence
  elements <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)

  # Define the elemental composition for each amino acid
  amino_acids <- list(
    A = c(C = 3, H = 5, N = 1, O = 1, S = 0, P = 0),
    R = c(C = 6, H = 12, N = 4, O = 1, S = 0, P = 0),
    N = c(C = 4, H = 6, N = 2, O = 2, S = 0, P = 0),
    D = c(C = 4, H = 5, N = 1, O = 3, S = 0, P = 0),
    E = c(C = 5, H = 7, N = 1, O = 3, S = 0, P = 0),
    Q = c(C = 5, H = 8, N = 2, O = 2, S = 0, P = 0),
    G = c(C = 2, H = 3, N = 1, O = 1, S = 0, P = 0),
    H = c(C = 6, H = 7, N = 3, O = 1, S = 0, P = 0),
    I = c(C = 6, H = 11, N = 1, O = 1, S = 0, P = 0),
    L = c(C = 6, H = 11, N = 1, O = 1, S = 0, P = 0),
    K = c(C = 6, H = 12, N = 2, O = 1, S = 0, P = 0),
    M = c(C = 5, H = 9, N = 1, O = 1, S = 1, P = 0),
    F = c(C = 9, H = 9, N = 1, O = 1, S = 0, P = 0),
    P = c(C = 5, H = 7, N = 1, O = 1, S = 0, P = 0),
    S = c(C = 3, H = 5, N = 1, O = 2, S = 0, P = 0),
    T = c(C = 4, H = 7, N = 1, O = 2, S = 0, P = 0),
    W = c(C = 11, H = 10, N = 2, O = 1, S = 0, P = 0),
    Y = c(C = 9, H = 9, N = 1, O = 2, S = 0, P = 0),
    V = c(C = 5, H = 9, N = 1, O = 1, S = 0, P = 0),
    C_iaa = c(C = 5, H = 8, N = 1, O = 2, S = 1, P = 0),
    C_no_iaa = c(C = 3, H = 5, N = 1, O = 1, S = 1, P = 0)
  )

  # Convert sequence into a vector of characters
  seq_vector <- strsplit(sequence, split = "")[[1]]

  # Loop through each amino acid in the sequence
  for (aa in seq_vector) {
    # Determine if it's C and choose composition based on IAA
    if (aa == "C") {
      elements <- elements + if (iaa)
        amino_acids$C_iaa
      else
        amino_acids$C_no_iaa
    } else if (aa %in% names(amino_acids)) {
      # Add the elemental composition for other amino acids
      elements <- elements + amino_acids[[aa]]
    }
  }

  elements
}

add_glycan <- function(toadd) {
  elements <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)

  #This mass is the sum of Hex and HexNac, normalized to a 1 Da unit.
  a <- list(A = c(
    C = 0.038321192,
    H = 0.062956245,
    N = 0.002737228,
    O = 0.027372280,
    S = 0,
    P = 0
  ))

  seq_vector <- rep("A", toadd)

  for (aa in seq_vector) {
    if (aa %in% names(a)) {
      elements <- elements + a[[aa]]
    }
  }

  elements
}

get_clean_total_comp <- function(total_comp) {
  total_comp <- total_comp[total_comp != 0]
  total_comp <- round(total_comp)
  total_comp <- paste0(names(total_comp), total_comp, collapse = "")
  total_comp
}

clean_isotope_duplicates <- function(iso_tempate) {
  iso_tempate |>
    dplyr::mutate(isotope_seq = NA) |>
    dplyr::arrange(.data$mz) |>
    dplyr::mutate(group = cumsum(c(TRUE, diff(.data$mz) > 0.01))) |>
    dplyr::group_by(.data$group) |>
    dplyr::slice_max(.data$intensity, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-"group")
}

#This is the R version producing identical results as the cpp code
isotopefit_opt <- function(spectrum, isotope_template, poi,
                           ppm_tolerance, cosine_score_correction,
                           minimal_isotope_sequence, score, plot = FALSE) {

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

  # Update template
  isotope_template$mz <- isotope_template$mz + difmz
  isotope_template$intensity <- isotope_template$intensity * difint

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
  if (!all(minimal_isotope_sequence %in% existing_seq)) return(0)

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
      fit <- cor(final_temp, final_spec)
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
    }
  }

  return(fit)
}

isotopefit_opt_old <- function(spectrum, isotope_template, poi,
                           ppm_tolerance, cosine_score_correction,
                           minimal_isotope_sequence, score, plot = FALSE) {
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
  idx <- tryCatch(
    {
      findInterval(template_mz, spectrum_mz)
    },
    error = function(e) {
      message("Error in findInterval with template_mz length: ",
              length(template_mz),
              ", spectrum_mz length: ", length(spectrum_mz),
              " -> ", e$message)
      NA
    }
  )
  idx[idx == 0] <- 1
  idx[idx > length(spectrum_mz)] <- length(spectrum_mz)

  # Choose the closer of idx and idx+1
  closest <- ifelse(
    idx < length(spectrum_mz) & abs(spectrum_mz[idx + 1] - template_mz) <
      abs(spectrum_mz[idx] - template_mz),
    idx + 1, idx
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

  #Return 0 if no minimal isotope sequence
  if (!all(minimal_isotope_sequence %in% existing_seq)) {
    return(0)
  }

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
      plot_title <- paste0("Pearson correlation coefficient: ", fit)
    }else if (score == "cosine") {
      template_intensity <- isotope_template$intensity
      isotope_intensity <- isotope_template$centroided_intensity
      offset <- min(template_intensity, isotope_intensity, na.rm = TRUE) *
        cosine_score_correction
      template_intensity <- template_intensity - offset
      isotope_intensity <- isotope_intensity - offset
      fit <- lsa::cosine(template_intensity, isotope_intensity)[[1]]
      plot_title <- paste0("Cosine score: ", fit)
    }
  }else {
    plot_title <- ""
  }

  # --- Step 6: Optional plot ---
  if (plot) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_bar(
        data = spectrum,
        ggplot2::aes(x = .data$mz, y = .data$centroided_intensity),
        stat = "identity",
        fill = "red",
        alpha = 0.7,
        width = 0.05
      ) +
      ggplot2::geom_bar(
        data = isotope_template,
        ggplot2::aes(x = .data$mz, y = .data$intensity),
        stat = "identity",
        fill = "blue",
        alpha = 0.7,
        width = 0.05
      ) +
      ggplot2::labs(title = plot_title, x = "m/z", y = "Intensity") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(
          0,
          max(spectrum$centroided_intensity) * 1.05
        )
      )

    suppressWarnings(print(p))
  }

  fit
}

compute_qvalues <- function(targets, decoys, cutoff) {
  correction <- 1
  targets <- as.numeric(as.vector(unlist(targets)))
  decoys <- as.numeric(as.vector(unlist(decoys)))

  scores <- c(targets, decoys)
  labels <- c(rep(1, length(targets)),
              rep(0, length(decoys)))

  df <- data.frame(score = scores, target = labels)
  df <- df[order(df$score, decreasing = TRUE), ]

  df$cum_target <- cumsum(df$target)
  df$cum_decoy  <- cumsum(1 - df$target)

  df$fdr <- (correction * df$cum_decoy) / df$cum_target
  df$qvalue <- rev(cummin(rev(df$fdr)))

  return_val <- df |>
    dplyr::filter(.data$qvalue <= cutoff) |>
    dplyr::pull(.data$score) |>
    min(na.rm = TRUE)

  return_val
}

plot_fdr_results <- function(target_scores, decoy_scores, cutoff_value) {
  target_df <- data.frame(score = as.vector(unlist(target_scores)),
                          type = "target") |>
    dplyr::filter(.data$score != 0) |>
    dplyr::mutate()

  target_df$bin <- binner(target_df$score, 0.006)
  target_df <- target_df |>
    dplyr::count(.data$bin, name = "count")

  decoy_df <- data.frame(score = as.vector(unlist(decoy_scores)),
                         type = "decoy") |>
    dplyr::filter(.data$score != 0)

  decoy_df$bin <- binner(decoy_df$score, 0.006)
  decoy_df <- decoy_df |>
    dplyr::count(.data$bin, name = "count")

  print(
    ggplot2::ggplot() +
      ggplot2::geom_vline(
        xintercept = cutoff_value,
        color = "grey50",
        linetype = "dashed"
      ) +
      ggplot2::geom_line(
        data = decoy_df,
        ggplot2::aes(x = .data$bin, y = .data$count),
        color = "black"
      ) +
      ggplot2::geom_line(
        data = target_df,
        ggplot2::aes(x = .data$bin, y = .data$count),
        color = "red"
      ) +
      ggplot2::labs(x = "Score", y = "Number of matches") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  )
}

binner <- function(x, binwidth) {
  bins <- as.character(cut(x, breaks = seq(floor(min(x)), ceiling(max(x)) + 1,
                                           by = binwidth),
                           include.lowest = TRUE))
  bins <- gsub("\\(|\\]", "", bins)

  minbin <- sapply(bins, function(x) as.numeric(strsplit(x, ",")[[1]][1]))
  maxbin <- sapply(bins, function(x) as.numeric(strsplit(x, ",")[[1]][2]))

  meanbin <- (minbin + maxbin) / 2

  meanbin
}
