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
# Paths
# ------------------------------------------------------------------------------

input_tables_dir <- file.path("output", "tables")
output_figures_dir <- file.path("output", "publication", "figures")
output_tables_dir  <- file.path("output", "publication", "tables")

dir.create(output_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_tables_dir,  recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Visual settings
# ------------------------------------------------------------------------------

drying_levels <- c("P", "C", "L")
ageing_levels <- c("T", "H", "B")
comparison_levels <- c("drying_vs_fresh", "ageing_vs_dried", "total_vs_fresh")

drying_labels <- c(
  P = "Air-dried",
  C = "Oven-dried",
  L = "Lyophilized"
)

ageing_labels <- c(
  T = "Temperature",
  H = "Humidity",
  B = "Both"
)

comparison_labels <- c(
  drying_vs_fresh = "Dried vs fresh",
  ageing_vs_dried = "Aged vs dried",
  total_vs_fresh  = "Aged vs fresh"
)

response_labels <- c(
  rmse = "RMSE",
  sam  = "SAM",
  iauc = "iAUC"
)

index_labels <- c(
  delta_mfdre       = "MFDRE",
  delta_rep         = "REP",
  delta_res700_740  = "RES700-740",
  delta_mean760_900 = "Mean760-900",
  delta_datt        = "Datt",
  delta_pri         = "PRI",
  delta_sipi        = "SIPI",
  delta_psri        = "PSRI"
)

drying_colors <- c(
  P = "#C49A00",
  C = "#D55E00",
  L = "#0072B2"
)

reference_colors <- c(
  fresh = "#222222",
  ageing = "#555555"
)

ageing_linetypes <- c(
  T = "dashed",
  H = "dotted",
  B = "dotdash"
)

heatmap_low  <- "#2166AC"
heatmap_mid  <- "#F7F7F7"
heatmap_high <- "#B2182B"

treatment_label_order <- c(
  drying_labels["P"], drying_labels["C"], drying_labels["L"],
  paste0(drying_labels["P"], " / ", ageing_labels["T"]),
  paste0(drying_labels["P"], " / ", ageing_labels["H"]),
  paste0(drying_labels["P"], " / ", ageing_labels["B"]),
  paste0(drying_labels["C"], " / ", ageing_labels["T"]),
  paste0(drying_labels["C"], " / ", ageing_labels["H"]),
  paste0(drying_labels["C"], " / ", ageing_labels["B"]),
  paste0(drying_labels["L"], " / ", ageing_labels["T"]),
  paste0(drying_labels["L"], " / ", ageing_labels["H"]),
  paste0(drying_labels["L"], " / ", ageing_labels["B"])
)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

sample_spectra <- read_sample_spectra_table(
  file.path(input_tables_dir, "analysis_sample_spectra.csv")
)

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
# FIGURE 2. Spectral trajectories across specimen history
# ==============================================================================

fig2_spectra_data <- tibble()
fig2_delta_data <- tibble()

if (nrow(sample_spectra) > 0) {
  fresh_mean <- sample_spectra %>%
    filter(timepoint == "fresh") %>%
    group_by(wavelength) %>%
    summarise(y_value = mean(reflectance, na.rm = TRUE), .groups = "drop")
  
  dried_mean <- sample_spectra %>%
    filter(timepoint == "dried", !is.na(drying)) %>%
    group_by(drying, wavelength) %>%
    summarise(y_value = mean(reflectance, na.rm = TRUE), .groups = "drop")
  
  aged_mean <- sample_spectra %>%
    filter(timepoint == "aged", !is.na(drying), !is.na(ageing)) %>%
    group_by(drying, ageing, wavelength) %>%
    summarise(y_value = mean(reflectance, na.rm = TRUE), .groups = "drop")
  
  if (nrow(fresh_mean) > 0 && nrow(dried_mean) > 0) {
    fig2_spectra_data <- bind_rows(
      fig2_spectra_data,
      fresh_mean %>%
        mutate(
          panel = "Fresh vs dried spectra",
          panel_order = 1,
          series_id = "fresh",
          color_key = "fresh",
          linetype_key = "solid",
          linewidth_key = "reference",
          plot_group = "spectra",
          y_measure = "reflectance"
        ),
      dried_mean %>%
        mutate(
          panel = "Fresh vs dried spectra",
          panel_order = 1,
          series_id = paste0("dried__", drying),
          color_key = as.character(drying),
          linetype_key = "solid",
          linewidth_key = "treatment",
          plot_group = "spectra",
          y_measure = "reflectance"
        )
    )
  }
  
  if (nrow(fresh_mean) > 0 && nrow(aged_mean) > 0) {
    fig2_spectra_data <- bind_rows(
      fig2_spectra_data,
      fresh_mean %>%
        mutate(
          panel = "Fresh vs aged spectra",
          panel_order = 2,
          series_id = "fresh",
          color_key = "fresh",
          linetype_key = "solid",
          linewidth_key = "reference",
          plot_group = "spectra",
          y_measure = "reflectance"
        ),
      aged_mean %>%
        mutate(
          panel = "Fresh vs aged spectra",
          panel_order = 2,
          series_id = paste0("aged__", drying, "__", ageing),
          color_key = as.character(drying),
          linetype_key = as.character(ageing),
          linewidth_key = "treatment",
          plot_group = "spectra",
          y_measure = "reflectance"
        )
    )
  }
}

if (nrow(delta_all) > 0) {
  drying_delta_mean <- delta_all %>%
    filter(comparison_block == "drying_vs_fresh", !is.na(drying)) %>%
    group_by(drying, wavelength) %>%
    summarise(y_value = mean(delta_reflectance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      panel = "Dried - fresh",
      panel_order = 1,
      series_id = paste0("delta_drying__", drying),
      color_key = as.character(drying),
      linetype_key = "solid",
      linewidth_key = "treatment",
      plot_group = "delta",
      y_measure = "delta_reflectance"
    )
  
  ageing_delta_mean <- delta_all %>%
    filter(comparison_block == "ageing_vs_dried", !is.na(drying), !is.na(ageing)) %>%
    group_by(drying, ageing, wavelength) %>%
    summarise(y_value = mean(delta_reflectance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      panel = "Aged - dried",
      panel_order = 2,
      series_id = paste0("delta_ageing__", drying, "__", ageing),
      color_key = as.character(drying),
      linetype_key = as.character(ageing),
      linewidth_key = "treatment",
      plot_group = "delta",
      y_measure = "delta_reflectance"
    )
  
  total_delta_mean <- delta_all %>%
    filter(comparison_block == "total_vs_fresh", !is.na(drying), !is.na(ageing)) %>%
    group_by(drying, ageing, wavelength) %>%
    summarise(y_value = mean(delta_reflectance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      panel = "Aged - fresh",
      panel_order = 3,
      series_id = paste0("delta_total__", drying, "__", ageing),
      color_key = as.character(drying),
      linetype_key = as.character(ageing),
      linewidth_key = "treatment",
      plot_group = "delta",
      y_measure = "delta_reflectance"
    )
  
  fig2_delta_data <- bind_rows(
    drying_delta_mean,
    ageing_delta_mean,
    total_delta_mean
  )
}

fig2_source_data <- bind_rows(fig2_spectra_data, fig2_delta_data)

if (nrow(fig2_source_data) > 0) {
  fig2_source_data <- fig2_source_data %>%
    mutate(
      color_key = factor(color_key, levels = c("fresh", drying_levels)),
      linetype_key = as.character(linetype_key),
      linewidth_key = factor(linewidth_key, levels = c("reference", "treatment"))
    )
  
  write_csv(
    fig2_source_data,
    file.path(output_tables_dir, "07_Figure2_source_data.csv")
  )
}

fig2_plots <- list()

if (nrow(fig2_spectra_data) > 0) {
  fig2_spectra_data <- fig2_spectra_data %>%
    mutate(
      panel = factor(panel, levels = unique(panel[order(panel_order)])),
      color_key = factor(color_key, levels = c("fresh", drying_levels)),
      linetype_key = as.character(linetype_key),
      linewidth_key = factor(linewidth_key, levels = c("reference", "treatment"))
    )
  
  p_fig2_spectra <- ggplot(
    fig2_spectra_data,
    aes(
      x = wavelength,
      y = y_value,
      group = series_id,
      color = color_key,
      linetype = linetype_key,
      linewidth = linewidth_key
    )
  ) +
    geom_line() +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_color_manual(
      values = c(
        fresh = reference_colors["fresh"],
        drying_colors
      )
    ) +
    scale_linetype_manual(
      values = c(
        solid = "solid",
        ageing_linetypes
      )
    ) +
    scale_linewidth_manual(
      values = c(reference = 0.95, treatment = 0.60),
      guide = "none"
    ) +
    labs(
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 5.5)
    )
  
  fig2_plots <- c(fig2_plots, list(p_fig2_spectra))
}

if (nrow(fig2_delta_data) > 0) {
  fig2_delta_data <- fig2_delta_data %>%
    mutate(
      panel = factor(panel, levels = unique(panel[order(panel_order)])),
      color_key = factor(color_key, levels = c("fresh", drying_levels)),
      linetype_key = as.character(linetype_key),
      linewidth_key = factor(linewidth_key, levels = c("reference", "treatment"))
    )
  
  delta_hline_data <- fig2_delta_data %>%
    distinct(panel)
  
  p_fig2_delta <- ggplot(
    fig2_delta_data,
    aes(
      x = wavelength,
      y = y_value,
      group = series_id,
      color = color_key,
      linetype = linetype_key,
      linewidth = linewidth_key
    )
  ) +
    geom_hline(
      data = delta_hline_data,
      aes(yintercept = 0),
      inherit.aes = FALSE,
      color = "grey70",
      linetype = 2,
      linewidth = 0.35
    ) +
    geom_line() +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_color_manual(
      values = c(
        fresh = reference_colors["fresh"],
        drying_colors
      )
    ) +
    scale_linetype_manual(
      values = c(
        solid = "solid",
        ageing_linetypes
      )
    ) +
    scale_linewidth_manual(
      values = c(reference = 0.95, treatment = 0.60),
      guide = "none"
    ) +
    labs(
      x = "Wavelength (nm)",
      y = expression(Delta * " reflectance")
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 8)
    )
  
  fig2_plots <- c(fig2_plots, list(p_fig2_delta))
}

if (length(fig2_plots) > 0) {
  fig2_legend_plot <- ggplot(
    fig2_source_data,
    aes(
      x = wavelength,
      y = y_value,
      group = series_id,
      color = color_key,
      linetype = linetype_key,
      linewidth = linewidth_key
    )
  ) +
    geom_line() +
    scale_color_manual(
      values = c(
        fresh = reference_colors["fresh"],
        drying_colors
      ),
      breaks = c("fresh", drying_levels),
      labels = c("Fresh", unname(drying_labels))
    ) +
    scale_linetype_manual(
      values = c(
        solid = "solid",
        ageing_linetypes
      ),
      breaks = ageing_levels,
      labels = unname(ageing_labels)
    ) +
    scale_linewidth_manual(
      values = c(reference = 0.95, treatment = 0.60),
      guide = "none"
    ) +
    labs(
      color = "Reference / drying",
      linetype = "Ageing"
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  fig2_legend_grob <- extract_legend_grob(fig2_legend_plot)
  
  save_arranged_plots_with_shared_legend(
    plots = fig2_plots,
    filename = file.path(output_figures_dir, "07_Figure2_spectral_trajectories.png"),
    title = "Figure 2. Spectral trajectories across specimen history",
    legend_grob = fig2_legend_grob,
    widths = rep(1, length(fig2_plots)),
    width = 12,
    height = if (length(fig2_plots) == 2) 8 else 7,
    dpi = 300
  )
}

# ==============================================================================
# FIGURE 3. Global spectral divergence metrics
# ==============================================================================

fig3_raw <- tibble()
fig3_emm <- tibble()

if (nrow(comparison_all) > 0) {
  raw_primary_base <- comparison_all %>%
    select(comparison_block, drying, ageing, rmse, sam, iauc) %>%
    pivot_longer(
      cols = c(rmse, sam, iauc),
      names_to = "response",
      values_to = "value"
    )
  
  raw_primary_drying <- raw_primary_base %>%
    filter(!is.na(drying)) %>%
    mutate(
      term = "drying",
      level = as.character(drying),
      plot_color = as.character(drying)
    )
  
  raw_primary_ageing <- raw_primary_base %>%
    filter(!is.na(ageing)) %>%
    mutate(
      term = "ageing",
      level = as.character(ageing),
      plot_color = "ageing"
    )
  
  fig3_raw <- bind_rows(raw_primary_drying, raw_primary_ageing) %>%
    mutate(
      comparison_block = factor(comparison_block, levels = comparison_levels, labels = unname(comparison_labels)),
      response = factor(response, levels = names(response_labels), labels = unname(response_labels)),
      term = factor(term, levels = c("drying", "ageing"), labels = c("Drying effect", "Ageing effect")),
      level = as.character(level),
      plot_color = as.character(plot_color)
    )
}

if (nrow(emmeans_primary) > 0) {
  fig3_emm <- emmeans_primary %>%
    filter(response %in% c("rmse", "sam", "iauc")) %>%
    mutate(
      level = as.character(level),
      plot_color = if_else(term == "drying", level, "ageing", missing = "ageing"),
      comparison_block = factor(comparison_block, levels = comparison_levels, labels = unname(comparison_labels)),
      response = factor(response, levels = names(response_labels), labels = unname(response_labels)),
      term = factor(term, levels = c("drying", "ageing"), labels = c("Drying effect", "Ageing effect"))
    )
}

if (nrow(fig3_emm) > 0) {
  x_label_map <- c(drying_labels, ageing_labels)
  
  p_fig3 <- ggplot() +
    geom_point(
      data = fig3_raw,
      aes(x = level, y = value, color = plot_color),
      position = position_jitter(width = 0.10, height = 0),
      alpha = 0.45,
      size = 1.2
    ) +
    geom_pointrange(
      data = fig3_emm,
      aes(x = level, y = emmean, ymin = lower_cl, ymax = upper_cl, color = plot_color),
      linewidth = 0.40,
      fatten = 1.4
    ) +
    facet_grid(response ~ comparison_block + term, scales = "free_y", space = "free_x") +
    scale_color_manual(
      values = c(
        drying_colors,
        ageing = reference_colors["ageing"]
      ),
      guide = "none"
    ) +
    scale_x_discrete(labels = x_label_map) +
    labs(
      title = "Figure 3. Global spectral divergence metrics",
      x = NULL,
      y = "Estimated marginal mean (with raw observations)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot_if_has_data(
    plot_obj = p_fig3,
    data = fig3_emm,
    filename = file.path(output_figures_dir, "07_Figure3_primary_metrics.png"),
    width = 12,
    height = 8
  )
  
  write_csv(
    fig3_raw,
    file.path(output_tables_dir, "07_Figure3_raw_points.csv")
  )
  
  write_csv(
    fig3_emm,
    file.path(output_tables_dir, "07_Figure3_emmeans.csv")
  )
}

# ==============================================================================
# FIGURE 4. Mean changes in derived spectral indices
# ==============================================================================

selected_indices <- names(index_labels)

fig4_data <- tibble()

if (nrow(comparison_all) > 0) {
  fig4_data <- comparison_all %>%
    select(comparison_block, drying, ageing, any_of(selected_indices)) %>%
    mutate(
      treatment = if_else(
        comparison_block == "drying_vs_fresh",
        unname(drying_labels[drying]),
        make_treatment_label(
          drying = drying,
          ageing = ageing,
          drying_labels = drying_labels,
          ageing_labels = ageing_labels
        )
      )
    ) %>%
    pivot_longer(
      cols = any_of(selected_indices),
      names_to = "metric",
      values_to = "delta_value"
    ) %>%
    group_by(comparison_block, treatment, metric) %>%
    summarise(
      mean_delta = if (all(is.na(delta_value))) NA_real_ else mean(delta_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      comparison_block = factor(comparison_block, levels = comparison_levels, labels = unname(comparison_labels)),
      treatment = factor(treatment, levels = treatment_label_order),
      metric = factor(metric, levels = rev(selected_indices), labels = rev(unname(index_labels)))
    )
}

if (nrow(fig4_data) > 0) {
  heat_limit <- quantile(abs(fig4_data$mean_delta), 0.7, na.rm = TRUE)
  p_fig4 <- ggplot(
    fig4_data,
    aes(x = treatment, y = metric, fill = mean_delta)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    facet_wrap(~ comparison_block, ncol = 1, scales = "free_x") +
    scale_fill_gradient2(
      low = heatmap_low,
      mid = heatmap_mid,
      high = heatmap_high,
      midpoint = 0,
      limits = c(-heat_limit, heat_limit),
      oob = scales::squish,
      na.value = "grey90"
    ) +
    labs(
      title = "Figure 4. Mean changes in derived spectral indices",
      x = NULL,
      y = NULL,
      fill = "Mean \u0394 index"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot_if_has_data(
    plot_obj = p_fig4,
    data = fig4_data,
    filename = file.path(output_figures_dir, "07_Figure4_index_heatmap.png"),
    width = 11,
    height = 8
  )
  
  write_csv(
    fig4_data,
    file.path(output_tables_dir, "07_Figure4_source_data.csv")
  )
}

message("07_publication_outputs.R finished successfully.")
message("Figures saved to: ", normalizePath(output_figures_dir, winslash = "/", mustWork = FALSE))
message("Tables saved to: ", normalizePath(output_tables_dir, winslash = "/", mustWork = FALSE))