plot_spectra <- function(full_mgf, output_folder) {
  outputfile <- paste(normalizePath(output_folder,
                                    winslash = "/", mustWork = FALSE),
                      "all_spectra_plot.pdf", sep = "/")

  pdf(outputfile, width = 6, height = 6)

  for (i in unique(full_mgf$id)){
    tempdf <- full_mgf |>
      dplyr::filter(.data$id == i)

    rect <- data.frame(minx = as.numeric(),
                       maxx = as.numeric(),
                       miny = as.numeric(),
                       maxy = as.numeric(),
                       charge = as.character())

    for (i in unique(tempdf$slice)) {
      tempdf_sub <- tempdf |>
        dplyr::filter(.data$slice == i)

      minx <- min(tempdf_sub$mz)
      maxx <- max(tempdf_sub$mz)
      miny <- min(tempdf_sub$intensity)
      maxy <- max(tempdf_sub$intensity)
      charge <- tempdf_sub$charge[1]

      rect <- rbind(rect, data.frame(minx = minx,
                                     maxx = maxx,
                                     miny = miny,
                                     maxy = maxy,
                                     charge = charge))
    }

    max_gap <- 1

    gap_rows <- tempdf |>
      dplyr::mutate(next_mz = dplyr::lead(.data$mz),
                    next_intensity = dplyr::lead(.data$intensity)) |>
      dplyr::filter(!is.na(.data$next_mz) &
                      (.data$next_mz - .data$mz) > max_gap) |>
      # for each gap, create two rows: 0 after current and 0 before next
      dplyr::reframe(
        mz = c(.data$mz + 1e-6, .data$next_mz - 1e-6),
        intensity = 0,
        id = .data$id,
        .by = c("id", "mz", "next_mz")
      )

    # Combine with original
    tempdf_filled <- dplyr::bind_rows(tempdf, gap_rows) |>
      dplyr::arrange(.data$mz)

    p <- ggplot2::ggplot(data = tempdf_filled,
                         ggplot2::aes(x = .data$mz, y = .data$intensity)) +
      ggplot2::geom_line() +
      ggplot2::geom_text(data = rect, ggplot2::aes(x = (minx + maxx) / 2,
                                                   y = maxy + (0.05 * maxy),
                                                   label = charge),
                         color = "red") +
      ggplot2::labs(title = tempdf$id[1]) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_x_continuous(n.breaks = 10)

    print(p)
  }

  dev.off()
}
