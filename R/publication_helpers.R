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
    "se",
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
    legend_height <- grid::grobHeight(legend_grob) +
      grid::unit(0.5, "lines")
    
    row_heights <- c(
      row_heights,
      list(legend_height)
    )
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

# ------------------------------------------------------------------------------
# Shared figure design (used by scripts/07_publication_outputs.R)
# ------------------------------------------------------------------------------

# Colour-vision-safe palette for drying methods. Chosen by an exhaustive search over
# Okabe-Ito/Tol colours for the largest minimum pairwise distance (CIELAB) under normal,
# deuteranopic and protanopic vision, ordered light -> dark so that the three methods
# also differ in lightness (readable in grayscale).
pub_colors_drying <- c(P = "#E69F00", C = "#117733", L = "#332288")
pub_color_fresh   <- "#666666"
pub_color_dried   <- "#111111"   # dried reference in spectra plots (drying method is given by the column)

# Ageing regimes get their own colours (temperature = red-orange, humidity = blue,
# both = purple, i.e. the mixture); these differ from the drying-method colours and are
# always combined with point shapes where points are drawn.
pub_colors_ageing <- c(T = "#D55E00", H = "#0072B2", B = "#CC79A7")

pub_labels_drying <- c(P = "Air-dried", C = "Oven-dried", L = "Lyophilized")
pub_labels_drying_short <- c(P = "Air", C = "Oven", L = "Lyoph.")
pub_labels_ageing <- c(T = "Temperature", H = "Humidity", B = "Both")
# full wording for figures (labels on the figure instead of abbreviations plus legend)
pub_labels_ageing_full <- c(
  T = "Ageing with temperature",
  H = "Ageing with humidity",
  B = "Ageing with temperature and humidity (both)"
)
pub_labels_ageing_with <- c(T = "With temperature", H = "With humidity", B = "With temperature and humidity")
pub_labels_ageing_short <- c(T = "Temp.", H = "Humid.", B = "Both")

# Ageing regimes: point shapes (filled, so they take fill + outline)
pub_shapes_ageing <- c(T = 24, H = 22, B = 21)

theme_pub <- function(base_size = 7) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.ticks = ggplot2::element_line(colour = "grey40", linewidth = 0.3),
      strip.text = ggplot2::element_text(size = base_size, colour = "black"),
      strip.background = ggplot2::element_rect(fill = "grey93", colour = NA),
      panel.border = ggplot2::element_rect(colour = "grey40", linewidth = 0.4, fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey91", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.key.size = grid::unit(3.5, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(2, 3, 2, 2, "mm")
    )
}

# Text size given in points (geom_text/annotate take mm)
pt_to_mm <- function(pt) pt / ggplot2::.pt

# Vector PDF (cairo) + high-resolution PNG; sizes in mm (figures are designed at final size)
save_publication_figure <- function(plot, path_base, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  ggplot2::ggsave(
    paste0(path_base, ".pdf"), plot,
    width = w, height = h, device = grDevices::cairo_pdf, bg = "white"
  )

  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(
    paste0(path_base, ".png"), plot,
    width = w, height = h, dpi = dpi, device = png_device, bg = "white"
  )

  invisible(path_base)
}

format_p_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "P < 0.001",
    p < 0.01 ~ sprintf("P = %.3f", p),
    p >= 0.995 ~ "P = 1",
    TRUE ~ sprintf("P = %.2f", p)
  )
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

# Compact letter display from pairwise contrasts ("A - B" with adjusted p-values).
# `levels_ordered` should be sorted by decreasing estimated marginal mean so that "a"
# marks the largest mean. Returns NA for all levels when no contrast is significant.
compact_letters <- function(levels_ordered, contrasts_df, alpha = 0.05) {
  n <- length(levels_ordered)
  sig <- matrix(FALSE, n, n, dimnames = list(levels_ordered, levels_ordered))

  for (i in seq_len(nrow(contrasts_df))) {
    parts <- trimws(strsplit(contrasts_df$contrast[i], " - ", fixed = TRUE)[[1]])

    if (length(parts) == 2 && all(parts %in% levels_ordered)) {
      is_sig <- isTRUE(contrasts_df$p_value_adjusted[i] < alpha)
      sig[parts[1], parts[2]] <- is_sig
      sig[parts[2], parts[1]] <- is_sig
    }
  }

  if (!any(sig)) {
    return(stats::setNames(rep(NA_character_, n), levels_ordered))
  }

  subsets <- unlist(
    lapply(seq_len(n), function(k) utils::combn(levels_ordered, k, simplify = FALSE)),
    recursive = FALSE
  )

  no_sig_pair <- function(s) {
    if (length(s) == 1) return(TRUE)
    m <- sig[s, s, drop = FALSE]
    !any(m[upper.tri(m)])
  }

  ok <- Filter(no_sig_pair, subsets)
  maximal <- Filter(
    function(s) !any(vapply(ok, function(o) length(o) > length(s) && all(s %in% o), logical(1))),
    ok
  )

  first_pos <- vapply(maximal, function(s) min(match(s, levels_ordered)), numeric(1))
  maximal <- maximal[order(first_pos, -lengths(maximal))]

  out <- stats::setNames(character(n), levels_ordered)

  for (k in seq_along(maximal)) {
    out[maximal[[k]]] <- paste0(out[maximal[[k]]], letters[k])
  }

  out
}
