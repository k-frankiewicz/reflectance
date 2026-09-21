# Figures use non-ASCII symbols (minus sign, multiplication sign, Delta, en dash).
# Under a non-UTF-8 locale (e.g. Rscript started with LANG unset) they would be drawn
# as dots, so switch the character type to UTF-8 when necessary.
if (!l10n_info()[["UTF-8"]]) {
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

source(file.path("R", "publication_helpers.R"))

# ------------------------------------------------------------------------------
# Paths and global settings
# ------------------------------------------------------------------------------

input_tables_dir <- file.path("output", "tables")
output_figures_dir <- file.path("output", "publication", "figures")
output_tables_dir  <- file.path("output", "publication", "tables")

dir.create(output_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_tables_dir,  recursive = TRUE, showWarnings = FALSE)

# Figures are designed at final print size (mm). 174 mm = full page width (two columns);
# use 84 mm for a single column. Text is >= 6 pt at this size.
fig_width_mm <- 174
base_size <- 7

# Figure numbers as in the manuscript (Fig. 2 is a photograph that is not produced by this
# workflow). Publication files are named after these numbers: "Fig. 1.pdf", "Fig. 1.png", ...
# and the matching source-data tables "Fig. 1 source data.csv", ...
figure_numbers <- c(
  experimental_design = 1,
  spectral_trajectories = 3,
  primary_metrics = 4,
  index_heatmap = 5
)
fig_path <- function(key) file.path(output_figures_dir, paste0("Fig. ", figure_numbers[[key]]))
source_path <- function(key, suffix = "") {
  file.path(output_tables_dir, paste0("Fig. ", figure_numbers[[key]], " source data", suffix, ".csv"))
}

drying_levels <- c("P", "C", "L")
ageing_levels <- c("T", "H", "B")
comparison_levels <- c("drying_vs_fresh", "ageing_vs_dried", "total_vs_fresh")
primary_responses <- c("rmse", "sam", "iauc")

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

sample_spectra <- read_sample_spectra_table(
  file.path(input_tables_dir, "analysis_sample_spectra.csv")
)

sample_indices <- safe_read_csv(file.path(input_tables_dir, "analysis_sample_indices.csv"))

comparison_all <- bind_rows(
  read_comparison_with_block("comparison_drying_vs_fresh.csv", "drying_vs_fresh", input_tables_dir),
  read_comparison_with_block("comparison_ageing_vs_dried.csv", "ageing_vs_dried", input_tables_dir),
  read_comparison_with_block("comparison_total_vs_fresh.csv",  "total_vs_fresh", input_tables_dir)
)

delta_all <- bind_rows(
  read_delta_with_block("delta_spectra_drying_vs_fresh.csv", "drying_vs_fresh", input_tables_dir),
  read_delta_with_block("delta_spectra_ageing_vs_dried.csv", "ageing_vs_dried", input_tables_dir),
  read_delta_with_block("delta_spectra_total_vs_fresh.csv",  "total_vs_fresh", input_tables_dir)
)

emmeans_primary <- read_emmeans_primary_table(
  file.path(input_tables_dir, "stats_emmeans_primary.csv")
)

model_overview <- safe_read_csv(file.path(input_tables_dir, "stats_model_overview.csv"))
tests_primary  <- safe_read_csv(file.path(input_tables_dir, "stats_tests_primary.csv"))
contrasts_primary <- safe_read_csv(file.path(input_tables_dir, "stats_contrasts_primary.csv"))
tests_secondary <- safe_read_csv(file.path(input_tables_dir, "stats_tests_secondary_indices.csv"))
aged_design_check <- safe_read_csv(file.path(input_tables_dir, "analysis_aged_design_check.csv"))
file_qc <- safe_read_csv(file.path(input_tables_dir, "qc_file_level_results.csv"))

# ------------------------------------------------------------------------------
# Metadata on available outputs
# ------------------------------------------------------------------------------

available_blocks <- comparison_all %>%
  distinct(comparison_block) %>%
  mutate(comparison_block = factor(comparison_block, levels = comparison_levels)) %>%
  arrange(comparison_block)

write_csv(
  available_blocks,
  file.path(output_tables_dir, "07_available_blocks.csv")
)

if (nrow(model_overview) > 0) {
  write_csv(
    model_overview,
    file.path(output_tables_dir, "07_model_overview_copy.csv")
  )
}

if (nrow(tests_primary) > 0) {
  write_csv(
    tests_primary,
    file.path(output_tables_dir, "07_primary_tests_copy.csv")
  )
}

# ==============================================================================
# FIGURE: experimental design (schematic)
#   A  workflow: fresh leaf -> three drying protocols -> Latin-square allocation to three
#      accelerated-ageing regimes, a spectral measurement after every stage, and the three
#      comparisons analysed (drying, ageing and total effect)
#   B  allocation of dried pieces to ageing regimes (samples per cell; black dots = one
#      example individual)
#   C  programmed 8-h cycles of the three ageing regimes (repeated for 28 days)
#   Protocol values (temperatures, humidities, durations) are taken from the Methods;
#   sample numbers and the number of replicate spectra come from the analysis tables.
# ==============================================================================

if (!requireNamespace("cowplot", quietly = TRUE)) {
  message("Package 'cowplot' is needed for the experimental-design figure: skipped.")
} else {
  n_of <- function(tp) sum(sample_indices$timepoint == tp, na.rm = TRUE)

  n_individuals <- sample_indices %>%
    filter(timepoint == "fresh") %>%
    distinct(line, individual) %>%
    nrow()
  n_lines <- n_distinct(sample_indices$line)
  rep_range <- range(sample_indices$n_replicates_used, na.rm = TRUE)
  rep_label <- paste0(rep_range[1], "–", rep_range[2], " per sample")

  design_counts <- tibble(
    stage = c("fresh", "dried", "aged"),
    n_samples = c(n_of("fresh"), n_of("dried"), n_of("aged")),
    replicate_spectra_min = rep_range[1],
    replicate_spectra_max = rep_range[2]
  )

  # blend a colour with white (amount = share of white)
  tint <- function(cols, amount = 0.86) {
    vapply(cols, function(cl) {
      v <- grDevices::col2rgb(cl)[, 1] / 255
      v <- v + (1 - v) * amount
      grDevices::rgb(v[1], v[2], v[3])
    }, character(1), USE.NAMES = FALSE)
  }

  # rounded rectangle as polygon vertices
  rr_poly <- function(xmin, xmax, ymin, ymax, r = 1.5, id = 1, n = 8) {
    th <- seq(0, pi / 2, length.out = n)
    corner <- function(cx, cy, a0) tibble(x = cx + r * cos(a0 + th), y = cy + r * sin(a0 + th))
    bind_rows(
      corner(xmax - r, ymax - r, 0),
      corner(xmin + r, ymax - r, pi / 2),
      corner(xmin + r, ymin + r, pi),
      corner(xmax - r, ymin + r, 3 * pi / 2)
    ) %>% mutate(id = id)
  }

  # curly brace spanning y0..y1 at x = xa, tip pointing in direction dir (+1 right, -1 left)
  curly <- function(xa, y0, y1, dir = 1, r = 1.9, n = 10) {
    ym <- (y0 + y1) / 2
    arc <- function(cx, cy, a0, a1) {
      t <- seq(a0, a1, length.out = n)
      tibble(x = cx + r * cos(t), y = cy + r * sin(t))
    }
    pts <- bind_rows(
      arc(xa, y1 - r, pi / 2, 0),
      arc(xa + 2 * r, ym + r, pi, 3 * pi / 2),
      arc(xa + 2 * r, ym - r, pi / 2, pi),
      arc(xa, y0 + r, 0, -pi / 2)
    )
    pts$x <- xa + dir * (pts$x - xa)
    pts
  }

  # mini reflectance spectrum icon (schematic leaf spectrum) fitted into a box
  mini_spectrum <- function(x0, y0, w, h, id) {
    t <- seq(0, 1, length.out = 60)
    v <- 0.10 + 0.14 * exp(-((t - 0.30) / 0.10)^2) + 0.62 * plogis((t - 0.58) * 28)
    tibble(x = x0 + t * w, y = y0 + v / max(v) * h, id = id)
  }

  # ---- geometry --------------------------------------------------------------------
  box_h <- 13.5
  box_gap <- 1.5
  y_top <- 75
  box_ymin <- y_top - box_h - (0:2) * (box_h + box_gap)   # 61.5, 46.5, 31.5
  box_ymax <- box_ymin + box_h
  stack_ymin <- min(box_ymin)
  stack_ymax <- max(box_ymax)
  stack_mid <- (stack_ymin + stack_ymax) / 2

  x_fresh <- c(0.8, 33.2)
  x_dry <- c(44, 88)
  x_age <- c(128, 173)
  drying_cols <- unname(pub_colors_drying[drying_levels])
  ageing_cols <- unname(pub_colors_ageing[ageing_levels])
  txt_pt <- base_size - 0.5

  boxes <- bind_rows(
    rr_poly(x_fresh[1], x_fresh[2], stack_ymin, stack_ymax, id = 1),
    bind_rows(lapply(1:3, function(i) rr_poly(x_dry[1], x_dry[2], box_ymin[i], box_ymax[i], id = 1 + i))),
    bind_rows(lapply(1:3, function(i) rr_poly(x_age[1], x_age[2], box_ymin[i], box_ymax[i], id = 4 + i)))
  )
  box_style <- tibble(
    id = 1:7,
    fill = c("#F2F2F2", tint(drying_cols), tint(ageing_cols)),
    edge = c(pub_color_fresh, drying_cols, ageing_cols)
  )
  boxes <- left_join(boxes, box_style, by = "id")

  # leaf icon with cut lines (fresh leaf cut into three pieces)
  leaf_x0 <- 4.5; leaf_len <- 26; leaf_y <- 62.6; leaf_w <- 4.6
  lx <- seq(0, leaf_len, length.out = 70)
  lh <- leaf_w * sin(pi * (lx / leaf_len)^0.75)^0.7
  leaf <- tibble(x = c(leaf_x0 + lx, rev(leaf_x0 + lx)), y = c(leaf_y + lh, rev(leaf_y - lh)))
  leaf_cuts <- tibble(
    x = leaf_x0 + leaf_len * c(1 / 3, 2 / 3), y = leaf_y - leaf_w - 1, yend = leaf_y + leaf_w + 1
  )

  # ---- texts -----------------------------------------------------------------------
  title_y <- box_ymin + 10.1
  proto_y <- box_ymin + 4.0
  flow_text <- bind_rows(
    tibble(x = 17, y = 70.6, label = "Fresh leaf", face = "bold", hjust = 0.5),
    tibble(x = 17, y = 55.6, label = "cut into three pieces", face = "italic", hjust = 0.5),
    tibble(x = 17, y = 43.4, hjust = 0.5, face = "plain",
           label = paste0("third fully expanded\nleaf of ", n_individuals, " individuals\nfrom ", n_lines, " maize inbred lines")),
    tibble(x = 66, y = 79, label = "Drying (three protocols)", face = "bold", hjust = 0.5),
    tibble(x = 150.5, y = 79, label = "Accelerated ageing (28 days)", face = "bold", hjust = 0.5),
    tibble(x = x_dry[1] + 4, y = c(title_y[1], proto_y[1]), hjust = 0, face = c("bold", "plain"),
           label = c("Air-dried", "room temperature, blotting paper")),
    tibble(x = x_dry[1] + 4, y = c(title_y[2], proto_y[2]), hjust = 0, face = c("bold", "plain"),
           label = c("Oven-dried", "38 °C for 72 h")),
    tibble(x = x_dry[1] + 4, y = c(title_y[3], proto_y[3]), hjust = 0, face = c("bold", "plain"),
           label = c("Lyophilized", "freeze-dried to constant mass")),
    tibble(x = x_age[1] + 4, y = c(title_y[1], proto_y[1]), hjust = 0, face = c("bold", "plain"),
           label = c(pub_labels_ageing_full[["T"]], "15 and 40 °C at 50% RH")),
    tibble(x = x_age[1] + 4, y = c(title_y[2], proto_y[2]), hjust = 0, face = c("bold", "plain"),
           label = c(pub_labels_ageing_full[["H"]], "25 and 75% RH at 20 °C")),
    tibble(x = x_age[1] + 4, y = c(box_ymin[3] + 9.5, proto_y[3] - 0.7), hjust = 0, face = c("bold", "plain"),
           label = c("Ageing with temperature\nand humidity (both)", "15 °C/75% RH and 40 °C/25% RH")),
    tibble(x = 108, y = 61, label = "Latin-square\nallocation (B)", face = "bold", hjust = 0.5),
    tibble(x = 108, y = 45, label = "every drying method\nreaches every regime", face = "plain", hjust = 0.5),
    # sample numbers
    tibble(x = 17, y = 28.6, label = paste0(design_counts$n_samples[1], " leaves"), face = "plain", hjust = 0.5),
    tibble(x = 66, y = 28.6, label = paste0(design_counts$n_samples[2], " samples (", design_counts$n_samples[2] / 3, " per drying method)"),
           face = "plain", hjust = 0.5),
    tibble(x = 150.5, y = 28.6, label = paste0(design_counts$n_samples[3], " samples (", design_counts$n_samples[3] / 3, " per ageing regime)"),
           face = "plain", hjust = 0.5),
    # comparison brackets
    tibble(x = 42, y = c(10.6, 7.6), label = c("Drying effect", "dried − fresh"), face = c("bold", "plain"), hjust = 0.5),
    tibble(x = 108.5, y = c(10.6, 7.6), label = c("Ageing effect", "aged − dried"), face = c("bold", "plain"), hjust = 0.5),
    tibble(x = 83.5, y = c(1.6, -1.4), label = c("Total effect", "aged − fresh"), face = c("bold", "plain"), hjust = 0.5)
  )

  # ---- spectra measurement pills (one after every stage) ------------------------------
  pill_x <- list(c(0.8, 33.2), c(49.8, 82.2), c(134.3, 166.7))
  pill_y <- c(19, 26.5)
  pills <- bind_rows(lapply(seq_along(pill_x), function(i) rr_poly(pill_x[[i]][1], pill_x[[i]][2], pill_y[1], pill_y[2], r = 2.2, id = i)))
  icons <- bind_rows(lapply(seq_along(pill_x), function(i) mini_spectrum(pill_x[[i]][1] + 2, pill_y[1] + 1.6, 6.5, 4.3, i)))
  pill_text <- bind_rows(lapply(seq_along(pill_x), function(i) tibble(
    x = pill_x[[i]][1] + 9.8, y = c(pill_y[2] - 2.4, pill_y[1] + 2.4),
    label = c("reflectance spectra", rep_label), face = c("bold", "plain"), hjust = 0
  )))

  # ---- arrows, braces, brackets --------------------------------------------------------
  box_mid <- box_ymin + box_h / 2
  fan <- tibble(x = x_fresh[2] + 0.6, xend = x_dry[1] - 0.6, y = stack_mid, yend = box_mid)
  main_arrow <- tibble(x = 94.8, xend = 121.2, y = stack_mid, yend = stack_mid)
  brace_l <- curly(90, stack_ymin, stack_ymax, dir = 1) %>% mutate(g = "dry")
  brace_r <- curly(126, stack_ymin, stack_ymax, dir = -1) %>% mutate(g = "age")
  braces <- bind_rows(brace_l, brace_r)

  bracket_y <- c(13.5, 4.6)
  brackets <- tribble(
    ~x, ~xend, ~y, ~yend,
    17, 65, bracket_y[1], bracket_y[1],
    17, 17, bracket_y[1], 16.6,
    65, 65, bracket_y[1], 16.6,
    67, 150.5, bracket_y[1], bracket_y[1],
    67, 67, bracket_y[1], 16.6,
    150.5, 150.5, bracket_y[1], 16.6,
    17, 150.5, bracket_y[2], bracket_y[2],
    17, 17, bracket_y[2], bracket_y[1],
    150.5, 150.5, bracket_y[2], bracket_y[1]
  )

  p_flow <- ggplot() +
    geom_polygon(data = boxes, aes(x = x, y = y, group = id), fill = boxes$fill, colour = boxes$edge, linewidth = 0.6) +
    geom_polygon(data = leaf, aes(x = x, y = y), fill = "#A8CC86", colour = "#4F7A3A", linewidth = 0.4) +
    geom_segment(aes(x = leaf_x0, xend = leaf_x0 + leaf_len, y = leaf_y, yend = leaf_y), colour = "#4F7A3A", linewidth = 0.3) +
    geom_segment(data = leaf_cuts, aes(x = x, xend = x, y = y, yend = yend), colour = "grey20", linewidth = 0.4, linetype = "dashed") +
    geom_polygon(data = pills, aes(x = x, y = y, group = id), fill = "#F6F6F6", colour = "grey65", linewidth = 0.35) +
    geom_path(data = icons, aes(x = x, y = y, group = id), colour = "grey25", linewidth = 0.5) +
    geom_path(data = braces, aes(x = x, y = y, group = g), colour = "grey35", linewidth = 0.55) +
    geom_segment(data = fan, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = grid::unit(1.7, "mm"), type = "closed"), colour = "grey35", linewidth = 0.45) +
    geom_segment(data = main_arrow, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = grid::unit(2.1, "mm"), type = "closed"), colour = "grey35", linewidth = 0.7) +
    geom_segment(data = brackets, aes(x = x, xend = xend, y = y, yend = yend), colour = "grey30", linewidth = 0.45) +
    geom_text(data = bind_rows(flow_text, pill_text),
              aes(x = x, y = y, label = label, fontface = face, hjust = hjust),
              size = pt_to_mm(txt_pt), lineheight = 0.95, vjust = 0.5, colour = "grey10") +
    coord_fixed(xlim = c(0, 174), ylim = c(-3, 81), expand = FALSE) +
    theme_void()

  # ---- panel B: Latin-square allocation ---------------------------------------------------
  alloc <- aged_design_check %>%
    mutate(
      x = match(ageing, ageing_levels),
      y = 4 - match(drying, drying_levels),
      example = paste(drying, ageing) %in% c("P T", "C H", "L B")
    )

  alloc_head_rows <- tibble(
    x = -0.28, y = 3:1, label = unname(pub_labels_drying[drying_levels]),
    fill = drying_cols, tcol = c("black", "white", "white")
  )
  alloc_head_cols <- tibble(
    x = 1:3, y = 3.95, label = c("Temperature", "Humidity", "Temperature\n+ humidity"),
    fill = ageing_cols, tcol = c("white", "white", "black")
  )

  p_alloc <- ggplot() +
    geom_tile(data = alloc, aes(x = x, y = y), width = 0.94, height = 0.94,
              fill = "white", colour = "grey60", linewidth = 0.35) +
    geom_point(data = filter(alloc, example), aes(x = x - 0.33, y = y + 0.31),
               shape = 16, size = 2.2, colour = "black") +
    geom_text(data = alloc, aes(x = x, y = y, label = n_samples),
              size = pt_to_mm(base_size + 2), fontface = "bold") +
    geom_tile(data = alloc_head_rows, aes(x = x, y = y), width = 1.24, height = 0.94,
              fill = alloc_head_rows$fill, colour = NA) +
    geom_text(data = alloc_head_rows, aes(x = x, y = y, label = label),
              colour = alloc_head_rows$tcol, size = pt_to_mm(base_size - 0.5)) +
    geom_tile(data = alloc_head_cols, aes(x = x, y = y), width = 0.99, height = 0.6,
              fill = alloc_head_cols$fill, colour = NA) +
    geom_text(data = alloc_head_cols, aes(x = x, y = y, label = label),
              colour = alloc_head_cols$tcol, size = pt_to_mm(base_size - 1)) +
    annotate("text", x = 2, y = 4.55, label = "Ageing with:", size = pt_to_mm(base_size - 0.5), fontface = "italic") +
    annotate("text", x = -0.28, y = 4.0, label = "Drying:", size = pt_to_mm(base_size - 0.5), fontface = "italic") +
    annotate("text", x = 1.4, y = 5.1, label = "Allocation of dried samples", fontface = "bold", size = pt_to_mm(base_size)) +
    annotate("text", x = 1.4, y = 0.05,
             label = "numbers: samples per cell\n\u25cf the three pieces of one example individual:\none per drying method, one per ageing regime",
             size = pt_to_mm(base_size - 1), lineheight = 0.95) +
    coord_cartesian(xlim = c(-0.95, 3.5), ylim = c(-0.45, 5.4), expand = FALSE) +
    theme_void()

  # ---- panel C: programmed cycles ---------------------------------------------------------------
  cycle_defs <- tribble(
    ~regime, ~variable,                     ~v1, ~v2, ~varied,
    "T", "Temperature (°C)",            15,  40,  TRUE,
    "T", "Relative humidity (%)",            50,  50,  FALSE,
    "H", "Temperature (°C)",            20,  20,  FALSE,
    "H", "Relative humidity (%)",            25,  75,  TRUE,
    "B", "Temperature (°C)",            15,  40,  TRUE,
    "B", "Relative humidity (%)",            75,  25,  TRUE
  )

  regime_strip <- c(
    T = "Ageing with\ntemperature",
    H = "Ageing with\nhumidity",
    B = "Ageing with temperature\nand humidity (both)"
  )

  # two 4-h blocks per cycle; 0.5-h ramp then 3.5-h hold; the cycle starts from the end state
  cycle_paths <- cycle_defs %>%
    rowwise() %>%
    reframe(
      regime = regime, variable = variable, varied = varied,
      t = c(0, 0.5, 4, 4.5, 8),
      value = c(v2, v1, v1, v2, v2)
    ) %>%
    mutate(
      regime_label = factor(regime, levels = ageing_levels, labels = unname(regime_strip[ageing_levels])),
      variable = factor(variable, levels = c("Temperature (°C)", "Relative humidity (%)"),
                        labels = c("Temperature\n(°C)", "Relative\nhumidity (%)")),
      line_colour = if_else(varied, unname(pub_colors_ageing[regime]), "grey55")
    )

  p_cycles <- ggplot(cycle_paths, aes(x = t, y = value, group = interaction(regime, variable))) +
    geom_line(aes(colour = line_colour, linewidth = varied)) +
    scale_colour_identity() +
    scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.5), guide = "none") +
    facet_grid(variable ~ regime_label, scales = "free_y", switch = "y") +
    scale_x_continuous(breaks = c(0, 4, 8), expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(x = "Time within the 8-h cycle (h)", y = NULL,
         title = "Programmed cycles (repeated for 28 days)") +
    theme_pub(base_size) +
    theme(
      strip.placement = "outside",
      plot.title = element_text(face = "bold", size = base_size, hjust = 0.5),
      plot.background = element_blank(),
      panel.spacing = grid::unit(2, "mm")
    )

  fig1 <- cowplot::plot_grid(
    p_flow,
    cowplot::plot_grid(
      p_alloc, p_cycles, nrow = 1, rel_widths = c(1.4, 2),
      labels = c("B", "C"), label_size = base_size + 3, label_fontface = "bold"
    ),
    ncol = 1, rel_heights = c(84, 66),
    labels = c("A", ""), label_size = base_size + 3, label_fontface = "bold"
  )

  save_publication_figure(
    fig1,
    fig_path("experimental_design"),
    width_mm = fig_width_mm, height_mm = 150
  )

  write_csv(design_counts, source_path("experimental_design"))
}

# ==============================================================================
# FIGURE: spectral trajectories across specimen history
#   columns = drying method; rows = mean spectra and the three difference spectra
#   colour = ageing regime (temperature / humidity / both), black = dried, grey = fresh;
#   the drying method is given by the column (and the coloured bar on top of it)
#   ribbons = 95% CI of the mean difference across individuals
#   NOTE: spectra are averaged into 1-nm bins for display only (raw data are sampled
#   at ~0.2 nm); statistics are not based on binned data.
# ==============================================================================

wl_bin_nm <- 1
ci_level <- 0.95
bin_wl <- function(x) round(x / wl_bin_nm) * wl_bin_nm

f2_row_labels <- c(
  spectra = "A  Mean reflectance (%)",
  drying  = "B  Dried − fresh (Δ %)",
  ageing  = "C  Aged − dried (Δ %)",
  total   = "D  Aged − fresh (Δ %)"
)

fig2_data <- tibble()

if (nrow(sample_spectra) > 0) {
  spectra_binned <- sample_spectra %>%
    filter(!is.na(reflectance)) %>%
    mutate(wl = bin_wl(wavelength)) %>%
    group_by(timepoint, drying, ageing, sample_group, wl) %>%
    summarise(r = mean(reflectance), .groups = "drop") %>%
    group_by(timepoint, drying, ageing, wl) %>%
    summarise(y = mean(r), n = n(), .groups = "drop")

  fresh_series <- spectra_binned %>%
    filter(timepoint == "fresh") %>%
    select(wl, y) %>%
    tidyr::crossing(col = drying_levels) %>%
    mutate(row = "spectra", series = "fresh")

  dried_series <- spectra_binned %>%
    filter(timepoint == "dried", !is.na(drying)) %>%
    transmute(row = "spectra", col = drying, series = "dried", wl, y)

  aged_series <- spectra_binned %>%
    filter(timepoint == "aged", !is.na(drying), !is.na(ageing)) %>%
    transmute(row = "spectra", col = drying, series = paste0("aged_", ageing), wl, y)

  fig2_data <- bind_rows(fresh_series, dried_series, aged_series)
}

if (nrow(delta_all) > 0) {
  delta_binned <- delta_all %>%
    filter(!is.na(delta_reflectance), !is.na(drying)) %>%
    mutate(wl = bin_wl(wavelength)) %>%
    group_by(comparison_block, comparison_id, drying, ageing, wl) %>%
    summarise(d = mean(delta_reflectance), .groups = "drop") %>%
    group_by(comparison_block, drying, ageing, wl) %>%
    summarise(n = n(), y = mean(d), sd = sd(d), .groups = "drop") %>%
    mutate(
      half = qt(1 - (1 - ci_level) / 2, df = pmax(n - 1, 1)) * sd / sqrt(n),
      lo = y - half,
      hi = y + half
    )

  delta_series <- delta_binned %>%
    transmute(
      row = recode(comparison_block,
                   drying_vs_fresh = "drying",
                   ageing_vs_dried = "ageing",
                   total_vs_fresh = "total"),
      col = drying,
      series = if_else(comparison_block == "drying_vs_fresh", "dried", paste0("aged_", ageing)),
      wl, y, lo, hi, n
    )

  fig2_data <- bind_rows(fig2_data, delta_series)
}

if (nrow(fig2_data) > 0) {
  fig2_data <- fig2_data %>%
    mutate(
      colour_key = as.character(series),
      row = factor(row, levels = names(f2_row_labels), labels = unname(f2_row_labels)),
      col = factor(col, levels = drying_levels, labels = unname(pub_labels_drying)),
      series = factor(series, levels = c("fresh", "dried", "aged_T", "aged_H", "aged_B"))
    )

  # wavelength ranges used for band means / red-edge indices (labelled in the first panel)
  bands <- tibble(
    label = c("B", "G", "R", "RE", "NIR"),
    xmin = c(450, 500, 650, 680, 760),
    xmax = c(500, 570, 680, 750, 900),
    shade = c(0.07, 0.03, 0.07, 0.03, 0.07)
  ) %>%
    mutate(xmid = (xmin + xmax) / 2)

  band_labels <- bands %>%
    mutate(
      row = factor(f2_row_labels["spectra"], levels = levels(fig2_data$row)),
      col = factor(pub_labels_drying["P"], levels = levels(fig2_data$col))
    )

  zero_lines <- tibble(row = factor(f2_row_labels[c("drying", "ageing", "total")],
                                    levels = levels(fig2_data$row)))

  # coloured cap on top of each column = drying method (same colours as in the other figures)
  drying_caps <- tibble(
    row = factor(f2_row_labels["spectra"], levels = levels(fig2_data$row)),
    col = factor(unname(pub_labels_drying), levels = levels(fig2_data$col)),
    cap = unname(pub_colors_drying[drying_levels]),
    x0 = 400,
    x1 = 950
  )

  series_colours <- c(
    fresh = pub_color_fresh,
    dried = pub_color_dried,
    aged_T = unname(pub_colors_ageing["T"]),
    aged_H = unname(pub_colors_ageing["H"]),
    aged_B = unname(pub_colors_ageing["B"])
  )

  p_fig2 <- ggplot(fig2_data, aes(x = wl, y = y, group = series)) +
    geom_rect(
      data = bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = shade),
      inherit.aes = FALSE, fill = "black"
    ) +
    scale_alpha_identity() +
    geom_hline(
      data = zero_lines, aes(yintercept = 0),
      colour = "grey35", linewidth = 0.25, linetype = 2
    ) +
    geom_ribbon(
      data = filter(fig2_data, !is.na(lo)),
      aes(ymin = lo, ymax = hi, fill = colour_key),
      alpha = 0.13, colour = NA
    ) +
    geom_line(aes(colour = colour_key, linewidth = series)) +
    geom_segment(
      data = drying_caps,
      aes(x = x0, xend = x1, y = Inf, yend = Inf),
      inherit.aes = FALSE, colour = drying_caps$cap, linewidth = 2.4, lineend = "butt"
    ) +
    geom_text(
      data = band_labels,
      aes(x = xmid, y = Inf, label = label),
      inherit.aes = FALSE, vjust = 1.9,
      size = pt_to_mm(base_size - 1.5), colour = "grey25"
    ) +
    facet_grid(row ~ col, scales = "free_y") +
    scale_colour_manual(
      values = series_colours,
      breaks = names(series_colours),
      labels = c("Fresh", "Dried", unname(pub_labels_ageing_full[ageing_levels])),
      name = NULL
    ) +
    scale_fill_manual(values = series_colours, guide = "none") +
    scale_linewidth_manual(
      values = c(fresh = 0.9, dried = 0.65, aged_T = 0.55, aged_H = 0.55, aged_B = 0.55),
      guide = "none"
    ) +
    scale_x_continuous(limits = c(400, 950), breaks = seq(400, 900, 100), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = "Wavelength (nm)", y = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1))) +
    theme_pub(base_size) +
    theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(8, "mm"),
      strip.text.y = element_text(hjust = 0),
      panel.spacing = grid::unit(2, "mm")
    )

  save_publication_figure(
    p_fig2,
    fig_path("spectral_trajectories"),
    width_mm = fig_width_mm, height_mm = 172
  )

  write_csv(
    fig2_data %>% mutate(across(where(is.factor), as.character)),
    source_path("spectral_trajectories")
  )
}

# ==============================================================================
# FIGURE: global spectral divergence metrics
#   raw observations + estimated marginal means (95% CI), omnibus P (Holm-adjusted)
#   and compact letters from Holm-adjusted pairwise contrasts (shown only when at
#   least one contrast is significant).
# ==============================================================================

fig3_panels <- tribble(
  ~comparison_block,   ~term,    ~panel,
  "drying_vs_fresh",   "drying", "Dried vs fresh\nby drying method",
  "ageing_vs_dried",   "drying", "Aged vs dried\nby drying method",
  "ageing_vs_dried",   "ageing", "Aged vs dried\nby ageing regime",
  "total_vs_fresh",    "drying", "Aged vs fresh\nby drying method",
  "total_vs_fresh",    "ageing", "Aged vs fresh\nby ageing regime"
)

fig3_response_labels <- c(
  rmse = "A  RMSE\n(% reflectance)",
  sam  = "B  SAM\n(rad)",
  iauc = "C  iAUC\n(% × nm)"
)

fig3_level_order <- c(drying_levels, ageing_levels)

fig3_raw <- tibble()
fig3_emm <- tibble()
fig3_omnibus <- tibble()
fig3_letters <- tibble()

if (nrow(comparison_all) > 0 && nrow(emmeans_primary) > 0) {
  raw_long <- comparison_all %>%
    select(comparison_block, line, drying, ageing, all_of(primary_responses)) %>%
    pivot_longer(all_of(primary_responses), names_to = "response", values_to = "value") %>%
    filter(!is.na(value))

  fig3_raw <- bind_rows(
    raw_long %>% filter(!is.na(drying)) %>% mutate(term = "drying", level = drying),
    raw_long %>% filter(!is.na(ageing)) %>% mutate(term = "ageing", level = ageing)
  ) %>%
    inner_join(fig3_panels, by = c("comparison_block", "term")) %>%
    mutate(
      level = factor(level, levels = fig3_level_order),
      panel = factor(panel, levels = fig3_panels$panel),
      response = factor(response, levels = names(fig3_response_labels), labels = unname(fig3_response_labels))
    )

  emm_sel <- emmeans_primary %>%
    filter(response %in% primary_responses, term %in% c("drying", "ageing"))

  fig3_emm <- emm_sel %>%
    inner_join(fig3_panels, by = c("comparison_block", "term")) %>%
    mutate(
      level = factor(level, levels = fig3_level_order),
      panel = factor(panel, levels = fig3_panels$panel),
      response = factor(response, levels = names(fig3_response_labels), labels = unname(fig3_response_labels))
    )

  letters_raw <- emm_sel %>%
    group_by(comparison_block, response, term) %>%
    group_modify(function(.x, .y) {
      ord <- .x %>% arrange(desc(emmean))

      cx <- contrasts_primary %>%
        filter(
          comparison_block == .y$comparison_block,
          response == .y$response,
          term == .y$term
        )

      tibble(
        level = ord$level,
        letters = unname(compact_letters(ord$level, cx))
      )
    }) %>%
    ungroup()

  # letters are shown only where the Holm-adjusted omnibus test is significant, so that
  # pairwise differences are not over-interpreted when the omnibus effect is not supported
  omnibus_significant <- tests_primary %>%
    filter(response %in% primary_responses, term %in% c("drying", "ageing")) %>%
    transmute(comparison_block, response, term, omnibus_sig = p_value_adjusted < 0.05)

  fig3_letters <- letters_raw %>%
    filter(!is.na(letters)) %>%
    inner_join(omnibus_significant, by = c("comparison_block", "response", "term")) %>%
    filter(omnibus_sig %in% TRUE) %>%
    inner_join(fig3_panels, by = c("comparison_block", "term")) %>%
    mutate(
      level = factor(level, levels = fig3_level_order),
      panel = factor(panel, levels = fig3_panels$panel),
      response = factor(response, levels = names(fig3_response_labels), labels = unname(fig3_response_labels))
    )

  fig3_omnibus <- tests_primary %>%
    filter(response %in% primary_responses, term %in% c("drying", "ageing")) %>%
    transmute(
      comparison_block, response, term,
      p_adjusted = p_value_adjusted,
      label = format_p_label(p_value_adjusted),
      face = if_else(p_value_adjusted < 0.05, "bold", "plain")
    ) %>%
    inner_join(fig3_panels, by = c("comparison_block", "term")) %>%
    mutate(
      panel = factor(panel, levels = fig3_panels$panel),
      response = factor(response, levels = names(fig3_response_labels), labels = unname(fig3_response_labels))
    )
}

if (nrow(fig3_emm) > 0) {
  fill_values <- c(pub_colors_drying, pub_colors_ageing)
  shape_values <- c(setNames(rep(21, length(drying_levels)), drying_levels), pub_shapes_ageing)
  x_labels <- c(pub_labels_drying, pub_labels_ageing_with)
  legend_breaks <- c(drying_levels, ageing_levels)
  legend_labels <- c(unname(pub_labels_drying[drying_levels]), unname(pub_labels_ageing_full[ageing_levels]))

  p_fig3 <- ggplot() +
    geom_point(
      data = fig3_raw,
      aes(x = level, y = value, fill = level, shape = level),
      position = position_jitter(width = 0.14, height = 0, seed = 20240921),
      size = 1.05, alpha = 0.7, colour = "white", stroke = 0.15
    ) +
    geom_linerange(
      data = fig3_emm,
      aes(x = level, ymin = lower_cl, ymax = upper_cl),
      linewidth = 0.5, colour = "black"
    ) +
    geom_point(
      data = fig3_emm,
      aes(x = level, y = emmean, fill = level),
      shape = 23, size = 2.4, colour = "black", stroke = 0.5, show.legend = FALSE
    ) +
    geom_text(
      data = fig3_omnibus,
      aes(x = 2, y = Inf, label = label, fontface = face),
      vjust = 1.5, size = pt_to_mm(base_size - 1)
    ) +
    geom_text(
      data = fig3_letters,
      aes(x = level, y = Inf, label = letters),
      vjust = 3.1, size = pt_to_mm(base_size), fontface = "bold"
    ) +
    facet_grid(response ~ panel, scales = "free", switch = "y") +
    scale_fill_manual(values = fill_values, breaks = legend_breaks, labels = legend_labels, name = NULL) +
    scale_shape_manual(values = shape_values, breaks = legend_breaks, labels = legend_labels, name = NULL) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.30))) +
    labs(x = NULL, y = NULL) +
    guides(
      fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.6, alpha = 1)),
      shape = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.6, alpha = 1))
    ) +
    theme_pub(base_size) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90),
      axis.text.x = element_text(angle = 40, hjust = 1),
      panel.spacing.x = grid::unit(c(4, 1.5, 4, 1.5), "mm"),
      panel.spacing.y = grid::unit(2, "mm"),
      legend.position = "bottom",
      legend.key.width = grid::unit(5, "mm")
    )

  save_publication_figure(
    p_fig3,
    fig_path("primary_metrics"),
    width_mm = fig_width_mm, height_mm = 144
  )

  write_csv(
    fig3_raw %>% mutate(across(where(is.factor), as.character)),
    source_path("primary_metrics", " - raw points")
  )

  write_csv(
    fig3_emm %>%
      mutate(across(where(is.factor), as.character)) %>%
      left_join(
        fig3_letters %>%
          transmute(comparison_block, response = as.character(response), term, level = as.character(level), letters),
        by = c("comparison_block", "response", "term", "level")
      ),
    source_path("primary_metrics", " - estimated marginal means")
  )
}

# ==============================================================================
# FIGURE: mean changes in derived spectral indices
#   colour and numbers = mean change standardised by the SD of the index among fresh leaves
#   (so that indices with different units share one scale; raw mean changes are in
#   "Fig. 5 source data.csv");
#   "Tests" columns = Benjamini-Hochberg-adjusted omnibus P (*** < 0.001, ** < 0.01, * < 0.05).
# ==============================================================================

index_defs <- tribble(
  ~metric,                     ~label,               ~group,
  "delta_mfdre",               "MFDRE",              "Red edge",
  "delta_rep",                 "REP",                "Red edge",
  "delta_res700_740",          "RES700–740",    "Red edge",
  "delta_mean_blue_450_500",   "Mean blue",          "Visible bands",
  "delta_mean_green_500_570",  "Mean green",         "Visible bands",
  "delta_mean_red_650_680",    "Mean red",           "Visible bands",
  "delta_datt",                "Datt",               "Pigment indices",
  "delta_pri",                 "PRI",                "Pigment indices",
  "delta_sipi",                "SIPI",               "Pigment indices",
  "delta_psri",                "PSRI",               "Pigment indices",
  "delta_mean760_900",         "Mean760–900",   "NIR"
)

fig4_z_limit <- 8

# compact tile label for the standardised change (one decimal; integer from 9.5 upwards),
# so that labels stay within a ~4.4 mm tile. Raw mean changes are in the source-data table.
fmt_z <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    a <- abs(v)
    f <- if (a >= 9.5) formatC(a, format = "f", digits = 0) else formatC(a, format = "f", digits = 1)
    paste0(if (v < 0 && f != "0.0") "−" else "", f)
  }, character(1))
}

# full x labels (rotated 90 degrees): drying methods for the dried-vs-fresh block, ageing regimes elsewhere
fig4_x_labels <- c(
  P = unname(pub_labels_drying[["P"]]), C = unname(pub_labels_drying[["C"]]), L = unname(pub_labels_drying[["L"]]),
  T = unname(pub_labels_ageing_with[["T"]]), H = unname(pub_labels_ageing_with[["H"]]), B = unname(pub_labels_ageing_with[["B"]])
)

fig4_tiles <- tibble()
fig4_tests <- tibble()

selected_metrics <- index_defs$metric[index_defs$metric %in% names(comparison_all)]

if (nrow(comparison_all) > 0 && nrow(sample_indices) > 0 && length(selected_metrics) > 0) {
  sd_fresh <- sample_indices %>%
    filter(timepoint == "fresh") %>%
    summarise(across(all_of(sub("^delta_", "", selected_metrics)), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "index", values_to = "sd_fresh") %>%
    mutate(metric = paste0("delta_", index))

  # panel layout: one "effects" panel per block (dried vs fresh) or per drying method
  # (ageing blocks), each followed by a "Tests" panel
  panel_layout <- bind_rows(
    tibble(comparison_block = "drying_vs_fresh", kind = "effects", drying = NA_character_),
    tibble(comparison_block = "drying_vs_fresh", kind = "tests",   drying = NA_character_),
    tibble(comparison_block = "ageing_vs_dried", kind = "effects", drying = drying_levels),
    tibble(comparison_block = "ageing_vs_dried", kind = "tests",   drying = NA_character_),
    tibble(comparison_block = "total_vs_fresh",  kind = "effects", drying = drying_levels),
    tibble(comparison_block = "total_vs_fresh",  kind = "tests",   drying = NA_character_)
  ) %>%
    mutate(panel_index = row_number()) %>%
    group_by(comparison_block) %>%
    mutate(first_in_block = panel_index == min(panel_index)) %>%
    ungroup() %>%
    mutate(
      block_label = c(
        drying_vs_fresh = "Dried vs fresh",
        ageing_vs_dried = "Aged vs dried",
        total_vs_fresh  = "Aged vs fresh"
      )[comparison_block],
      # Nested strips: the block title is a plotmath label (bold) shown on the first panel of
      # each block; the other panels get unique invisible titles (phantom) so that facets stay
      # distinct but the title is printed only once (see labeller in facet_grid below).
      block_title = if_else(
        first_in_block,
        paste0("bold('", block_label, "')"),
        paste0("phantom(", strrep("0", panel_index), ")")
      ),
      panel_label = case_when(
        kind == "tests" ~ "Tests",
        comparison_block == "drying_vs_fresh" ~ "Drying",
        TRUE ~ unname(pub_labels_drying[drying])
      )
    )

  block_title_levels <- panel_layout$block_title

  means_long <- comparison_all %>%
    select(comparison_block, drying, ageing, all_of(selected_metrics)) %>%
    pivot_longer(all_of(selected_metrics), names_to = "metric", values_to = "delta") %>%
    group_by(comparison_block, drying, ageing, metric) %>%
    summarise(mean_delta = mean(delta, na.rm = TRUE), n = sum(!is.na(delta)), .groups = "drop") %>%
    left_join(sd_fresh %>% select(metric, sd_fresh), by = "metric") %>%
    mutate(z = mean_delta / sd_fresh)

  fig4_tiles <- bind_rows(
    # dried vs fresh: x = drying method
    means_long %>%
      filter(comparison_block == "drying_vs_fresh") %>%
      mutate(kind = "effects", panel_drying = NA_character_, x_code = drying),
    # ageing blocks: panel = drying method, x = ageing regime
    means_long %>%
      filter(comparison_block != "drying_vs_fresh", !is.na(ageing)) %>%
      mutate(kind = "effects", panel_drying = drying, x_code = ageing)
  ) %>%
    inner_join(
      panel_layout %>% filter(kind == "effects") %>%
        select(comparison_block, kind, panel_drying = drying, block_title, panel_label),
      by = c("comparison_block", "kind", "panel_drying")
    ) %>%
    inner_join(index_defs, by = "metric") %>%
    mutate(
      x = factor(
        x_code,
        levels = names(fig4_x_labels),
        labels = unname(fig4_x_labels)
      ),
      metric_label = factor(label, levels = rev(index_defs$label)),
      group = factor(group, levels = unique(index_defs$group)),
      block_title = factor(block_title, levels = block_title_levels),
      raw_label = fmt_z(z),
      text_colour = if_else(abs(z) > 4.5, "white", "black")
    )

  term_labels <- c(drying = "Drying", ageing = "Ageing", `drying:ageing` = "Interaction")

  fig4_tests <- tests_secondary %>%
    filter(response %in% selected_metrics, term %in% names(term_labels)) %>%
    inner_join(
      panel_layout %>% filter(kind == "tests") %>% select(comparison_block, block_title, panel_label),
      by = "comparison_block"
    ) %>%
    inner_join(index_defs, by = c("response" = "metric")) %>%
    mutate(
      x = factor(unname(term_labels[term]), levels = unname(term_labels)),
      metric_label = factor(label, levels = rev(index_defs$label)),
      group = factor(group, levels = unique(index_defs$group)),
      block_title = factor(block_title, levels = block_title_levels),
      stars = p_to_stars(p_value_adjusted),
      is_sig = stars != ""
    )
}

if (nrow(fig4_tiles) > 0) {
  n_panels <- nrow(panel_layout)
  # wider gaps between blocks than between panels of the same block
  gaps <- vapply(seq_len(n_panels - 1), function(i) {
    if (panel_layout$comparison_block[i] != panel_layout$comparison_block[i + 1]) 5
    else if (panel_layout$kind[i + 1] == "tests") 2.5
    else 1
  }, numeric(1))

  p_fig4 <- ggplot() +
    geom_tile(
      data = fig4_tiles,
      aes(x = x, y = metric_label, fill = z),
      colour = "white", linewidth = 0.3
    ) +
    geom_text(
      data = fig4_tiles,
      aes(x = x, y = metric_label, label = raw_label, colour = text_colour),
      size = pt_to_mm(base_size - 1.5)
    ) +
    scale_colour_identity() +
    geom_tile(
      data = fig4_tests,
      aes(x = x, y = metric_label),
      fill = "white", colour = "grey85", linewidth = 0.3
    ) +
    geom_text(
      data = filter(fig4_tests, is_sig),
      aes(x = x, y = metric_label, label = stars),
      size = pt_to_mm(base_size + 1), vjust = 0.75, fontface = "bold"
    ) +
    geom_text(
      data = filter(fig4_tests, !is_sig),
      aes(x = x, y = metric_label, label = "ns"),
      size = pt_to_mm(base_size - 2), colour = "grey55"
    ) +
    facet_grid(
      group ~ block_title + panel_label,
      scales = "free", space = "free",
      labeller = labeller(block_title = label_parsed, .default = label_value)
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
      limits = c(-fig4_z_limit, fig4_z_limit), oob = scales::squish,
      breaks = seq(-fig4_z_limit, fig4_z_limit, length.out = 5),
      na.value = "grey90",
      name = "Mean Δ index (in SD of fresh leaves)"
    ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colourbar(
      title.position = "top", barwidth = grid::unit(55, "mm"), barheight = grid::unit(2.5, "mm")
    )) +
    theme_pub(base_size) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size - 1.5),
      strip.clip = "off",
      strip.text.x = element_text(hjust = 0),
      strip.text.y = element_text(hjust = 0.5),
      strip.background = element_blank(),
      panel.spacing.x = grid::unit(gaps, "mm"),
      panel.spacing.y = grid::unit(1.5, "mm"),
      legend.position = "bottom"
    )

  save_publication_figure(
    p_fig4,
    fig_path("index_heatmap"),
    width_mm = fig_width_mm, height_mm = 150
  )

  write_csv(
    fig4_tiles %>%
      transmute(
        comparison_block, drying, ageing, index = metric, group = as.character(group),
        mean_delta, sd_fresh, standardised_delta = z, n
      ),
    source_path("index_heatmap")
  )

  write_csv(
    fig4_tests %>%
      transmute(
        comparison_block, response, term, F = statistic, df1, df2,
        p_value, p_value_adjusted, stars
      ),
    source_path("index_heatmap", " - omnibus tests")
  )
}

message("07_publication_outputs.R finished successfully.")
message("Figures saved to: ", normalizePath(output_figures_dir, winslash = "/", mustWork = FALSE))
message("Tables saved to: ", normalizePath(output_tables_dir, winslash = "/", mustWork = FALSE))
