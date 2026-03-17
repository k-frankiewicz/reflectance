safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }
  
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) tibble::tibble()
  )
}

make_empty_typed_tibble <- function(chr_cols = character(), num_cols = character()) {
  out <- tibble::tibble()
  
  for (nm in chr_cols) {
    out[[nm]] <- character()
  }
  
  for (nm in num_cols) {
    out[[nm]] <- numeric()
  }
  
  out
}

ensure_typed_cols <- function(df, chr_cols = character(), num_cols = character()) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    df <- tibble::tibble()
  }
  
  for (nm in chr_cols) {
    if (!nm %in% names(df)) {
      df[[nm]] <- rep(NA_character_, nrow(df))
    }
    df[[nm]] <- as.character(df[[nm]])
  }
  
  for (nm in num_cols) {
    if (!nm %in% names(df)) {
      df[[nm]] <- rep(NA_real_, nrow(df))
    }
    df[[nm]] <- suppressWarnings(readr::parse_double(as.character(df[[nm]])))
  }
  
  tibble::as_tibble(df)
}

read_typed_csv <- function(path, chr_cols = character(), num_cols = character()) {
  df <- safe_read_csv(path)
  
  if (nrow(df) == 0 && ncol(df) == 0) {
    df <- make_empty_typed_tibble(chr_cols = chr_cols, num_cols = num_cols)
  } else {
    df <- ensure_typed_cols(df, chr_cols = chr_cols, num_cols = num_cols)
  }
  
  df
}

read_comparison_with_block <- function(filename, block_name, input_tables_dir) {
  comparison_chr_cols <- c(
    "comparison_type",
    "comparison_id",
    "line",
    "individual",
    "drying",
    "ageing",
    "sample_group_sample",
    "sample_group_ref",
    "timepoint_sample",
    "timepoint_ref",
    "comparison_block"
  )
  
  comparison_num_cols <- c(
    "rmse",
    "sam",
    "iauc",
    "n_replicates_sample",
    "n_replicates_ref",
    "delta_mfdre",
    "delta_rep",
    "delta_res700_740",
    "delta_mean_blue_450_500",
    "delta_mean_green_500_570",
    "delta_mean_red_650_680",
    "delta_mean760_900",
    "delta_datt",
    "delta_pri",
    "delta_sipi",
    "delta_psri"
  )
  
  path <- file.path(input_tables_dir, filename)
  
  df <- read_typed_csv(
    path,
    chr_cols = comparison_chr_cols,
    num_cols = comparison_num_cols
  )
  
  df$comparison_block <- block_name
  df
}

read_delta_with_block <- function(filename, block_name, input_tables_dir) {
  delta_chr_cols <- c(
    "comparison_type",
    "comparison_id",
    "line",
    "individual",
    "drying",
    "ageing",
    "comparison_block"
  )
  
  delta_num_cols <- c(
    "wavelength",
    "reflectance_ref",
    "reflectance_sample",
    "delta_reflectance"
  )
  
  path <- file.path(input_tables_dir, filename)
  
  df <- read_typed_csv(
    path,
    chr_cols = delta_chr_cols,
    num_cols = delta_num_cols
  )
  
  df$comparison_block <- block_name
  df
}

read_sample_spectra_table <- function(path) {
  chr_cols <- c(
    "sample_group",
    "timepoint",
    "line",
    "individual",
    "drying",
    "ageing"
  )
  
  num_cols <- c(
    "wavelength",
    "reflectance",
    "mean_reflectance",
    "sd_reflectance",
    "cv_reflectance",
    "n_replicates_used"
  )
  
  read_typed_csv(path, chr_cols = chr_cols, num_cols = num_cols)
}

read_emmeans_primary_table <- function(path) {
  chr_cols <- c(
    "comparison_block",
    "response_family",
    "response",
    "term",
    "level"
  )
  
  num_cols <- c(
    "emmean",
    "SE",
    "df",
    "lower_cl",
    "upper_cl",
    "statistic",
    "p_value",
    "p_value_adjusted"
  )
  
  read_typed_csv(path, chr_cols = chr_cols, num_cols = num_cols)
}

save_plot_if_has_data <- function(plot_obj, data, filename, width = 11, height = 7, dpi = 300) {
  if (inherits(plot_obj, "ggplot") && nrow(data) > 0) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi
    )
  }
}

make_treatment_label <- function(drying, ageing = NA_character_, drying_labels, ageing_labels) {
  drying <- as.character(drying)
  ageing <- as.character(ageing)
  
  out <- unname(drying_labels[drying])
  
  has_ageing <- !is.na(ageing) & ageing != ""
  out[has_ageing] <- paste0(
    unname(drying_labels[drying[has_ageing]]),
    " / ",
    unname(ageing_labels[ageing[has_ageing]])
  )
  
  out
}

extract_legend_grob <- function(plot_obj) {
  gt <- ggplot2::ggplotGrob(plot_obj)
  
  grob_names <- vapply(
    gt$grobs,
    function(x) {
      if (!is.null(x$name)) x$name else ""
    },
    character(1)
  )
  
  guide_idx <- which(grob_names == "guide-box")
  
  if (length(guide_idx) == 0) {
    return(NULL)
  }
  
  gt$grobs[[guide_idx[1]]]
}

save_arranged_plots_with_shared_legend <- function(
    plots,
    filename,
    title = NULL,
    legend_grob = NULL,
    widths = NULL,
    width = 12,
    height = 8,
    dpi = 300
) {
  plots <- Filter(function(x) inherits(x, "ggplot"), plots)
  
  if (length(plots) == 0) {
    return(invisible(FALSE))
  }
  
  if (is.null(widths)) {
    widths <- rep(1, length(plots))
  }
  
  widths <- rep(widths, length.out = length(plots))
  
  has_title <- !is.null(title) && nzchar(title)
  has_legend <- !is.null(legend_grob)
  
  row_heights <- list()
  
  if (has_title) {
    row_heights <- c(row_heights, list(grid::unit(0.9, "lines")))
  }
  
  row_heights <- c(row_heights, list(grid::unit(1, "null")))
  
  if (has_legend) {
    row_heights <- c(row_heights, list(grid::unit(1.4, "lines")))
  }
  
  layout <- grid::grid.layout(
    nrow = length(row_heights),
    ncol = length(plots),
    widths = grid::unit(widths, "null"),
    heights = do.call(grid::unit.c, row_heights)
  )
  
  grDevices::png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = layout))
  
  current_row <- 1
  
  if (has_title) {
    grid::grid.text(
      label = title,
      x = 0.01,
      hjust = 0,
      gp = grid::gpar(fontface = "bold", cex = 1.2),
      vp = grid::viewport(
        layout.pos.row = current_row,
        layout.pos.col = seq_along(plots)
      )
    )
    current_row <- current_row + 1
  }
  
  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      vp = grid::viewport(
        layout.pos.row = current_row,
        layout.pos.col = i
      )
    )
  }
  
  if (has_legend) {
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = current_row + 1,
        layout.pos.col = seq_along(plots)
      )
    )
    grid::grid.draw(legend_grob)
    grid::popViewport()
  }
  
  grid::popViewport()
  
  invisible(TRUE)
}