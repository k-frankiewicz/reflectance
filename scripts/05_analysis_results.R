suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

source(file.path("R", "replicate_qc_helpers.R"))
source(file.path("R", "analysis_helpers.R"))

manifest_path <- file.path("data", "metadata", "sample_manifest.csv")
processed_dir <- file.path("data", "processed")
file_qc_path <- file.path("output", "tables", "qc_file_level_results.csv")
replicate_qc_path <- file.path("output", "tables", "replicate_qc_results.csv")
output_tables_dir <- file.path("output", "tables")
output_figures_dir <- file.path("output", "figures")

analysis_min_wl <- 400
analysis_max_wl <- 950

# By default we keep file-level QC pass/warn and exclude only replicate-level outlier candidates
exclude_replicate_qc_status <- c("outlier_candidate")

index_names <- spectral_index_names()
comparison_metric_names <- c("rmse", "sam", "iauc")
summary_metrics <- c(comparison_metric_names, paste0("delta_", index_names))

heatmap_metrics <- paste0(
  "delta_",
  c("mfdre", "rep", "res700_740", "mean760_900", "datt", "pri", "sipi", "psri")
)

heatmap_labels <- c(
  delta_mfdre = "MFDRE",
  delta_rep = "REP",
  delta_res700_740 = "RES700-740",
  delta_mean760_900 = "Mean760-900",
  delta_datt = "Datt",
  delta_pri = "PRI",
  delta_sipi = "SIPI",
  delta_psri = "PSRI"
)

dir.create(output_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(manifest_path, file_qc_path, replicate_qc_path)) {
  if (!file.exists(path)) {
    stop(
      "Required input not found: ",
      normalizePath(path, winslash = "/", mustWork = FALSE)
    )
  }
}

manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
  mutate(file = basename(file))

file_qc <- read_csv(file_qc_path, show_col_types = FALSE) %>%
  mutate(file = basename(file))

replicate_qc <- read_csv(replicate_qc_path, show_col_types = FALSE) %>%
  mutate(file = basename(file))

eligible_files <- file_qc %>%
  filter(qc_status %in% c("pass", "warn")) %>%
  pull(file) %>%
  unique()

processed_files <- file.path(processed_dir, eligible_files)
processed_files <- processed_files[file.exists(processed_files)]

if (length(processed_files) == 0) {
  stop("No eligible processed files found for analysis.")
}

# ------------------------------------------------------------------------------
# Load spectra retained for analysis
# ------------------------------------------------------------------------------

spectra <- map_dfr(
  processed_files,
  read_one_processed_spectrum,
  analysis_min_wl = analysis_min_wl,
  analysis_max_wl = analysis_max_wl
) %>%
  left_join(manifest, by = "file") %>%
  left_join(
    file_qc %>%
      select(file, file_qc_status = qc_status, file_qc_reasons = qc_reasons),
    by = "file"
  ) %>%
  left_join(
    replicate_qc %>%
      select(
        file,
        replicate_qc_status,
        pearson_r,
        rmse_replicate_qc = rmse,
        n_replicates,
        assessed_for_replicate_qc
      ),
    by = "file"
  ) %>%
  mutate(
    sample_group = build_sample_group(
      timepoint = timepoint,
      line = line,
      individual = individual,
      drying = drying,
      ageing = ageing
    ),
    analysis_keep = file_qc_status %in% c("pass", "warn") &
      (is.na(replicate_qc_status) | !replicate_qc_status %in% exclude_replicate_qc_status)
  )

if (any(is.na(spectra$sample_group))) {
  stop("Could not build sample_group for some files. Check manifest fields.")
}

analysis_inclusion_summary <- spectra %>%
  distinct(file, timepoint, file_qc_status, replicate_qc_status, analysis_keep) %>%
  count(timepoint, file_qc_status, replicate_qc_status, analysis_keep, name = "n_files")

write_csv(
  analysis_inclusion_summary,
  file.path(output_tables_dir, "analysis_inclusion_summary.csv")
)

retained_spectra <- spectra %>%
  filter(analysis_keep)

if (nrow(retained_spectra) == 0) {
  stop("All spectra were excluded by the current analysis filters.")
}

write_csv(
  retained_spectra,
  file.path(output_tables_dir, "analysis_retained_spectra.csv")
)

# ------------------------------------------------------------------------------
# Build representative sample spectra
# ------------------------------------------------------------------------------

sample_replicate_counts <- retained_spectra %>%
  distinct(sample_group, timepoint, line, individual, drying, ageing, file) %>%
  count(sample_group, timepoint, line, individual, drying, ageing, name = "n_replicates_used")

representative_spectra <- retained_spectra %>%
  group_by(sample_group, timepoint, line, individual, drying, ageing, wavelength) %>%
  summarise(
    # replicate-level statistics must be computed BEFORE `reflectance` is
    # overwritten with the median (summarise() evaluates arguments sequentially)
    mean_reflectance = mean(reflectance, na.rm = TRUE),
    sd_reflectance = if (dplyr::n() >= 2) stats::sd(reflectance, na.rm = TRUE) else NA_real_,
    cv_reflectance = if (!is.na(mean_reflectance) && !is.na(sd_reflectance) && mean_reflectance != 0) {
      sd_reflectance / abs(mean_reflectance)
    } else {
      NA_real_
    },
    reflectance = median(reflectance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  relocate(reflectance, .before = mean_reflectance) %>%
  left_join(
    sample_replicate_counts,
    by = c("sample_group", "timepoint", "line", "individual", "drying", "ageing")
  )

write_csv(
  representative_spectra,
  file.path(output_tables_dir, "analysis_sample_spectra.csv")
)

# ------------------------------------------------------------------------------
# Compute indices: replicate-level and representative-sample-level
# ------------------------------------------------------------------------------

replicate_indices <- retained_spectra %>%
  group_by(file, sample_group, timepoint, line, individual, drying, ageing) %>%
  group_modify(~ compute_spectral_indices(.x %>% select(wavelength, reflectance))) %>%
  ungroup()

write_csv(
  replicate_indices,
  file.path(output_tables_dir, "analysis_replicate_indices.csv")
)

replicate_index_summary <- summarise_numeric_by_group(
  df = replicate_indices,
  group_cols = c("sample_group", "timepoint", "line", "individual", "drying", "ageing"),
  value_cols = index_names
)

write_csv(
  replicate_index_summary,
  file.path(output_tables_dir, "analysis_replicate_index_summary.csv")
)

sample_indices <- representative_spectra %>%
  distinct(sample_group, timepoint, line, individual, drying, ageing, n_replicates_used) %>%
  left_join(
    representative_spectra %>%
      group_by(sample_group) %>%
      group_modify(~ compute_spectral_indices(.x %>% select(wavelength, reflectance))) %>%
      ungroup(),
    by = "sample_group"
  )

write_csv(
  sample_indices,
  file.path(output_tables_dir, "analysis_sample_indices.csv")
)

message(
  "Timepoints available after QC filtering: ",
  paste(sort(unique(sample_indices$timepoint)), collapse = ", ")
)

# ------------------------------------------------------------------------------
# Check completeness of aged experimental design
# ------------------------------------------------------------------------------

expected_drying_levels <- c("P", "C", "L")
expected_ageing_levels <- c("T", "H", "B")

aged_design_counts <- sample_indices %>%
  filter(timepoint == "aged") %>%
  count(drying, ageing, name = "n_samples")

aged_design_check <- tidyr::expand_grid(
  drying = expected_drying_levels,
  ageing = expected_ageing_levels
) %>%
  left_join(
    aged_design_counts,
    by = c("drying", "ageing")
  ) %>%
  mutate(
    n_samples = dplyr::coalesce(n_samples, 0L)
  )

write_csv(
  aged_design_check,
  file.path(output_tables_dir, "analysis_aged_design_check.csv")
)

missing_design_cells <- aged_design_check %>%
  filter(n_samples == 0)

if (nrow(missing_design_cells) > 0) {
  stop(
    "Incomplete aged experimental design after QC. Missing drying × ageing cells: ",
    paste(
      paste0(
        missing_design_cells$drying,
        "/",
        missing_design_cells$ageing
      ),
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------------------------
# Convenience lookup tables for pairwise comparisons
# ------------------------------------------------------------------------------

sample_spectra_ref <- representative_spectra %>%
  select(sample_group, wavelength, reflectance) %>%
  rename(sample_group_ref = sample_group, reflectance_ref = reflectance)

sample_spectra_sample <- representative_spectra %>%
  select(sample_group, wavelength, reflectance) %>%
  rename(sample_group_sample = sample_group, reflectance_sample = reflectance)

empty_metrics <- empty_comparison_metrics(index_names)
empty_delta <- empty_delta_spectra()

# ==============================================================================
# BLOCK A. DRYING EFFECT: dried vs fresh
# ==============================================================================

fresh_idx <- sample_indices %>%
  filter(timepoint == "fresh") %>%
  transmute(
    line,
    individual,
    sample_group_ref = sample_group,
    timepoint_ref = timepoint,
    n_replicates_ref = n_replicates_used,
    across(all_of(index_names), ~ .x, .names = "{.col}_ref")
  )

dried_idx <- sample_indices %>%
  filter(timepoint == "dried") %>%
  transmute(
    line,
    individual,
    drying,
    sample_group_sample = sample_group,
    timepoint_sample = timepoint,
    n_replicates_sample = n_replicates_used,
    across(all_of(index_names), ~ .x, .names = "{.col}_sample")
  )

if (nrow(fresh_idx) > 0 && nrow(dried_idx) > 0) {
  drying_pair_indices <- dried_idx %>%
    inner_join(fresh_idx, by = c("line", "individual")) %>%
    mutate(
      comparison_type = "drying_vs_fresh",
      comparison_id = paste(comparison_type, line, individual, drying, sep = "__"),
      ageing = NA_character_
    ) %>%
    make_index_delta_table(index_names = index_names)
  
  drying_delta_spectra <- drying_pair_indices %>%
    select(
      comparison_type, comparison_id, line, individual, drying, ageing,
      sample_group_sample, sample_group_ref
    ) %>%
    distinct() %>%
    left_join(sample_spectra_sample, by = "sample_group_sample") %>%
    left_join(sample_spectra_ref, by = c("sample_group_ref", "wavelength")) %>%
    mutate(delta_reflectance = reflectance_sample - reflectance_ref) %>%
    select(
      comparison_type, comparison_id, line, individual, drying, ageing,
      wavelength, reflectance_ref, reflectance_sample, delta_reflectance
    )
  
  drying_metrics <- drying_delta_spectra %>%
    group_by(comparison_type, comparison_id, line, individual, drying, ageing) %>%
    summarise(
      rmse = spectral_rmse(reflectance_sample, reflectance_ref),
      sam = spectral_sam(reflectance_sample, reflectance_ref),
      iauc = spectral_iauc(wavelength, delta_reflectance),
      .groups = "drop"
    ) %>%
    left_join(
      drying_pair_indices %>%
        select(
          comparison_type, comparison_id, line, individual, drying, ageing,
          sample_group_sample, sample_group_ref,
          timepoint_sample, timepoint_ref,
          n_replicates_sample, n_replicates_ref,
          starts_with("delta_")
        ),
      by = c("comparison_type", "comparison_id", "line", "individual", "drying", "ageing")
    )
} else {
  drying_delta_spectra <- empty_delta
  drying_metrics <- empty_metrics
}

write_csv(
  drying_delta_spectra,
  file.path(output_tables_dir, "delta_spectra_drying_vs_fresh.csv")
)

write_csv(
  drying_metrics,
  file.path(output_tables_dir, "comparison_drying_vs_fresh.csv")
)

write_csv(
  summarise_numeric_by_group(
    drying_metrics,
    group_cols = c("drying"),
    value_cols = summary_metrics
  ),
  file.path(output_tables_dir, "summary_drying_vs_fresh_by_drying.csv")
)

# ==============================================================================
# BLOCK B. AGEING EFFECT: aged vs dried
# ==============================================================================

aged_idx <- sample_indices %>%
  filter(timepoint == "aged") %>%
  transmute(
    line,
    individual,
    drying,
    ageing,
    sample_group_sample = sample_group,
    timepoint_sample = timepoint,
    n_replicates_sample = n_replicates_used,
    across(all_of(index_names), ~ .x, .names = "{.col}_sample")
  )

dried_idx_ref <- sample_indices %>%
  filter(timepoint == "dried") %>%
  transmute(
    line,
    individual,
    drying,
    sample_group_ref = sample_group,
    timepoint_ref = timepoint,
    n_replicates_ref = n_replicates_used,
    across(all_of(index_names), ~ .x, .names = "{.col}_ref")
  )

if (nrow(aged_idx) > 0 && nrow(dried_idx_ref) > 0) {
  ageing_pair_indices <- aged_idx %>%
    inner_join(dried_idx_ref, by = c("line", "individual", "drying")) %>%
    mutate(
      comparison_type = "ageing_vs_dried",
      comparison_id = paste(comparison_type, line, individual, drying, ageing, sep = "__")
    ) %>%
    make_index_delta_table(index_names = index_names)
  
  ageing_delta_spectra <- ageing_pair_indices %>%
    select(
      comparison_type, comparison_id, line, individual, drying, ageing,
      sample_group_sample, sample_group_ref
    ) %>%
    distinct() %>%
    left_join(sample_spectra_sample, by = "sample_group_sample") %>%
    left_join(sample_spectra_ref, by = c("sample_group_ref", "wavelength")) %>%
    mutate(delta_reflectance = reflectance_sample - reflectance_ref) %>%
    select(
      comparison_type, comparison_id, line, individual, drying, ageing,
      wavelength, reflectance_ref, reflectance_sample, delta_reflectance
    )
  
  ageing_metrics <- ageing_delta_spectra %>%
    group_by(comparison_type, comparison_id, line, individual, drying, ageing) %>%
    summarise(
      rmse = spectral_rmse(reflectance_sample, reflectance_ref),
      sam = spectral_sam(reflectance_sample, reflectance_ref),
      iauc = spectral_iauc(wavelength, delta_reflectance),
      .groups = "drop"
    ) %>%
    left_join(
      ageing_pair_indices %>%
        select(
          comparison_type, comparison_id, line, individual, drying, ageing,
          sample_group_sample, sample_group_ref,
          timepoint_sample, timepoint_ref,
          n_replicates_sample, n_replicates_ref,
          starts_with("delta_")
        ),
      by = c("comparison_type", "comparison_id", "line", "individual", "drying", "ageing")
    )
} else {
  ageing_delta_spectra <- empty_delta
  ageing_metrics <- empty_metrics
}

write_csv(
  ageing_delta_spectra,
  file.path(output_tables_dir, "delta_spectra_ageing_vs_dried.csv")
)

write_csv(
  ageing_metrics,
  file.path(output_tables_dir, "comparison_ageing_vs_dried.csv")
)

write_csv(
  summarise_numeric_by_group(
    ageing_metrics,
    group_cols = c("drying", "ageing"),
    value_cols = summary_metrics
  ),
  file.path(output_tables_dir, "summary_ageing_vs_dried_by_drying_ageing.csv")
)

# ==============================================================================
# BLOCK C. TOTAL EFFECT: aged vs fresh
# ==============================================================================

fresh_idx_ref <- sample_indices %>%
  filter(timepoint == "fresh") %>%
  transmute(
    line,
    individual,
    sample_group_ref = sample_group,
    timepoint_ref = timepoint,
    n_replicates_ref = n_replicates_used,
    across(all_of(index_names), ~ .x, .names = "{.col}_ref")
  )

aged_idx_total <- sample_indices %>%
  filter(timepoint == "aged") %>%
  transmute(
    line,
    individual,
    drying,
    ageing,
    sample_group_sample = sample_group,
    timepoint_sample = timepoint,
    n_replicates_sample = n_replicates_used,
    across(all_of(index_names), ~ .x, .names = "{.col}_sample")
  )

if (nrow(fresh_idx_ref) > 0 && nrow(aged_idx_total) > 0) {
  total_pair_indices <- aged_idx_total %>%
    inner_join(fresh_idx_ref, by = c("line", "individual")) %>%
    mutate(
      comparison_type = "total_vs_fresh",
      comparison_id = paste(comparison_type, line, individual, drying, ageing, sep = "__")
    ) %>%
    make_index_delta_table(index_names = index_names)
  
  total_delta_spectra <- total_pair_indices %>%
    select(
      comparison_type, comparison_id, line, individual, drying, ageing,
      sample_group_sample, sample_group_ref
    ) %>%
    distinct() %>%
    left_join(sample_spectra_sample, by = "sample_group_sample") %>%
    left_join(sample_spectra_ref, by = c("sample_group_ref", "wavelength")) %>%
    mutate(delta_reflectance = reflectance_sample - reflectance_ref) %>%
    select(
      comparison_type, comparison_id, line, individual, drying, ageing,
      wavelength, reflectance_ref, reflectance_sample, delta_reflectance
    )
  
  total_metrics <- total_delta_spectra %>%
    group_by(comparison_type, comparison_id, line, individual, drying, ageing) %>%
    summarise(
      rmse = spectral_rmse(reflectance_sample, reflectance_ref),
      sam = spectral_sam(reflectance_sample, reflectance_ref),
      iauc = spectral_iauc(wavelength, delta_reflectance),
      .groups = "drop"
    ) %>%
    left_join(
      total_pair_indices %>%
        select(
          comparison_type, comparison_id, line, individual, drying, ageing,
          sample_group_sample, sample_group_ref,
          timepoint_sample, timepoint_ref,
          n_replicates_sample, n_replicates_ref,
          starts_with("delta_")
        ),
      by = c("comparison_type", "comparison_id", "line", "individual", "drying", "ageing")
    )
} else {
  total_delta_spectra <- empty_delta
  total_metrics <- empty_metrics
}

write_csv(
  total_delta_spectra,
  file.path(output_tables_dir, "delta_spectra_total_vs_fresh.csv")
)

write_csv(
  total_metrics,
  file.path(output_tables_dir, "comparison_total_vs_fresh.csv")
)

write_csv(
  summarise_numeric_by_group(
    total_metrics,
    group_cols = c("drying", "ageing"),
    value_cols = summary_metrics
  ),
  file.path(output_tables_dir, "summary_total_vs_fresh_by_drying_ageing.csv")
)

# ------------------------------------------------------------------------------
# Combined outputs
# ------------------------------------------------------------------------------

all_comparison_metrics <- bind_rows(
  drying_metrics,
  ageing_metrics,
  total_metrics
)

write_csv(
  all_comparison_metrics,
  file.path(output_tables_dir, "comparison_all_blocks.csv")
)

write_csv(
  summarise_numeric_by_group(
    all_comparison_metrics,
    group_cols = c("comparison_type"),
    value_cols = summary_metrics
  ),
  file.path(output_tables_dir, "summary_all_blocks_by_comparison_type.csv")
)

# ==============================================================================
# Exploratory plots
# ==============================================================================

# ------------------------------------------------------------------------------
# Drying block plots
# ------------------------------------------------------------------------------

if (nrow(drying_delta_spectra) > 0) {
  drying_spectra_plot_data <- drying_delta_spectra %>%
    select(drying, wavelength, reflectance_ref, reflectance_sample) %>%
    pivot_longer(
      cols = c(reflectance_ref, reflectance_sample),
      names_to = "state_key",
      values_to = "reflectance"
    ) %>%
    mutate(
      state = dplyr::recode(
        state_key,
        reflectance_ref = "fresh",
        reflectance_sample = "dried"
      )
    ) %>%
    group_by(drying, state, wavelength) %>%
    summarise(mean_reflectance = mean(reflectance, na.rm = TRUE), .groups = "drop")
  
  p_drying_spectra <- ggplot(
    drying_spectra_plot_data,
    aes(x = wavelength, y = mean_reflectance, color = state)
  ) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ drying, ncol = 1) +
    labs(
      title = "Representative spectra: dried vs fresh",
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "State"
    ) +
    theme_bw()
  
  save_plot_if_has_data(
    p_drying_spectra,
    file.path(output_figures_dir, "05_drying_vs_fresh_spectra.png"),
    width = 8,
    height = 9
  )
  
  drying_delta_plot_data <- drying_delta_spectra %>%
    group_by(drying, wavelength) %>%
    summarise(mean_delta = mean(delta_reflectance, na.rm = TRUE), .groups = "drop")
  
  p_drying_delta <- ggplot(
    drying_delta_plot_data,
    aes(x = wavelength, y = mean_delta, color = drying)
  ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.6) +
    labs(
      title = "Difference spectra: dried - fresh",
      x = "Wavelength (nm)",
      y = expression(Delta * "Reflectance"),
      color = "Drying"
    ) +
    theme_bw()
  
  save_plot_if_has_data(
    p_drying_delta,
    file.path(output_figures_dir, "05_drying_vs_fresh_delta_spectra.png"),
    width = 9,
    height = 6
  )
  
  drying_metric_plot_data <- drying_metrics %>%
    pivot_longer(
      cols = all_of(comparison_metric_names),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  if (nrow(drying_metric_plot_data) > 0) {
    p_drying_metrics <- ggplot(
      drying_metric_plot_data,
      aes(x = drying, y = value)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.8) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title = "Global spectral change metrics: dried vs fresh",
        x = "Drying method",
        y = "Value"
      ) +
      theme_bw()
    
    save_plot_if_has_data(
      p_drying_metrics,
      file.path(output_figures_dir, "05_drying_vs_fresh_metrics.png"),
      width = 10,
      height = 6
    )
  }
  
  drying_heatmap_data <- drying_metrics %>%
    group_by(drying) %>%
    summarise(
      across(
        any_of(heatmap_metrics),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = any_of(heatmap_metrics),
      names_to = "metric",
      values_to = "mean_delta"
    ) %>%
    mutate(
      metric_label = dplyr::recode(metric, !!!heatmap_labels),
      metric_label = factor(metric_label, levels = rev(unname(heatmap_labels)))
    )
  
  if (nrow(drying_heatmap_data) > 0) {
    p_drying_heatmap <- ggplot(
      drying_heatmap_data,
      aes(x = drying, y = metric_label, fill = mean_delta)
    ) +
      geom_tile() +
      labs(
        title = "Mean index changes: dried vs fresh",
        x = "Drying method",
        y = "Index",
        fill = "Mean delta"
      ) +
      theme_bw()
    
    save_plot_if_has_data(
      p_drying_heatmap,
      file.path(output_figures_dir, "05_drying_vs_fresh_index_heatmap.png"),
      width = 7,
      height = 5
    )
  }
}

# ------------------------------------------------------------------------------
# Ageing block plots
# ------------------------------------------------------------------------------

if (nrow(ageing_delta_spectra) > 0) {
  ageing_spectra_plot_data <- ageing_delta_spectra %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    select(treatment, wavelength, reflectance_ref, reflectance_sample) %>%
    pivot_longer(
      cols = c(reflectance_ref, reflectance_sample),
      names_to = "state_key",
      values_to = "reflectance"
    ) %>%
    mutate(
      state = dplyr::recode(
        state_key,
        reflectance_ref = "dried",
        reflectance_sample = "aged"
      )
    ) %>%
    group_by(treatment, state, wavelength) %>%
    summarise(mean_reflectance = mean(reflectance, na.rm = TRUE), .groups = "drop")
  
  p_ageing_spectra <- ggplot(
    ageing_spectra_plot_data,
    aes(x = wavelength, y = mean_reflectance, color = state)
  ) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ treatment) +
    labs(
      title = "Representative spectra: aged vs dried",
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "State"
    ) +
    theme_bw()
  
  save_plot_if_has_data(
    p_ageing_spectra,
    file.path(output_figures_dir, "05_ageing_vs_dried_spectra.png"),
    width = 10,
    height = 7
  )
  
  ageing_delta_plot_data <- ageing_delta_spectra %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    group_by(treatment, wavelength) %>%
    summarise(mean_delta = mean(delta_reflectance, na.rm = TRUE), .groups = "drop")
  
  p_ageing_delta <- ggplot(
    ageing_delta_plot_data,
    aes(x = wavelength, y = mean_delta, color = treatment)
  ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.6) +
    labs(
      title = "Difference spectra: aged - dried",
      x = "Wavelength (nm)",
      y = expression(Delta * "Reflectance"),
      color = "Treatment"
    ) +
    theme_bw()
  
  save_plot_if_has_data(
    p_ageing_delta,
    file.path(output_figures_dir, "05_ageing_vs_dried_delta_spectra.png"),
    width = 9,
    height = 6
  )
  
  ageing_metric_plot_data <- ageing_metrics %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    pivot_longer(
      cols = all_of(comparison_metric_names),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  if (nrow(ageing_metric_plot_data) > 0) {
    p_ageing_metrics <- ggplot(
      ageing_metric_plot_data,
      aes(x = treatment, y = value)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.8) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title = "Global spectral change metrics: aged vs dried",
        x = "Drying_ageing treatment",
        y = "Value"
      ) +
      theme_bw()
    
    save_plot_if_has_data(
      p_ageing_metrics,
      file.path(output_figures_dir, "05_ageing_vs_dried_metrics.png"),
      width = 11,
      height = 6
    )
  }
  
  ageing_heatmap_data <- ageing_metrics %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    group_by(treatment) %>%
    summarise(
      across(
        any_of(heatmap_metrics),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = any_of(heatmap_metrics),
      names_to = "metric",
      values_to = "mean_delta"
    ) %>%
    mutate(
      metric_label = dplyr::recode(metric, !!!heatmap_labels),
      metric_label = factor(metric_label, levels = rev(unname(heatmap_labels)))
    )
  
  if (nrow(ageing_heatmap_data) > 0) {
    p_ageing_heatmap <- ggplot(
      ageing_heatmap_data,
      aes(x = treatment, y = metric_label, fill = mean_delta)
    ) +
      geom_tile() +
      labs(
        title = "Mean index changes: aged vs dried",
        x = "Drying_ageing treatment",
        y = "Index",
        fill = "Mean delta"
      ) +
      theme_bw()
    
    save_plot_if_has_data(
      p_ageing_heatmap,
      file.path(output_figures_dir, "05_ageing_vs_dried_index_heatmap.png"),
      width = 8,
      height = 5
    )
  }
}

# ------------------------------------------------------------------------------
# Total block plots
# ------------------------------------------------------------------------------

if (nrow(total_delta_spectra) > 0) {
  total_spectra_plot_data <- total_delta_spectra %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    select(treatment, wavelength, reflectance_ref, reflectance_sample) %>%
    pivot_longer(
      cols = c(reflectance_ref, reflectance_sample),
      names_to = "state_key",
      values_to = "reflectance"
    ) %>%
    mutate(
      state = dplyr::recode(
        state_key,
        reflectance_ref = "fresh",
        reflectance_sample = "aged"
      )
    ) %>%
    group_by(treatment, state, wavelength) %>%
    summarise(mean_reflectance = mean(reflectance, na.rm = TRUE), .groups = "drop")
  
  p_total_spectra <- ggplot(
    total_spectra_plot_data,
    aes(x = wavelength, y = mean_reflectance, color = state)
  ) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ treatment) +
    labs(
      title = "Representative spectra: aged vs fresh",
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "State"
    ) +
    theme_bw()
  
  save_plot_if_has_data(
    p_total_spectra,
    file.path(output_figures_dir, "05_total_vs_fresh_spectra.png"),
    width = 10,
    height = 7
  )
  
  total_delta_plot_data <- total_delta_spectra %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    group_by(treatment, wavelength) %>%
    summarise(mean_delta = mean(delta_reflectance, na.rm = TRUE), .groups = "drop")
  
  p_total_delta <- ggplot(
    total_delta_plot_data,
    aes(x = wavelength, y = mean_delta, color = treatment)
  ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.6) +
    labs(
      title = "Difference spectra: aged - fresh",
      x = "Wavelength (nm)",
      y = expression(Delta * "Reflectance"),
      color = "Treatment"
    ) +
    theme_bw()
  
  save_plot_if_has_data(
    p_total_delta,
    file.path(output_figures_dir, "05_total_vs_fresh_delta_spectra.png"),
    width = 9,
    height = 6
  )
  
  total_metric_plot_data <- total_metrics %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    pivot_longer(
      cols = all_of(comparison_metric_names),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  if (nrow(total_metric_plot_data) > 0) {
    p_total_metrics <- ggplot(
      total_metric_plot_data,
      aes(x = treatment, y = value)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.8) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(
        title = "Global spectral change metrics: aged vs fresh",
        x = "Drying_ageing treatment",
        y = "Value"
      ) +
      theme_bw()
    
    save_plot_if_has_data(
      p_total_metrics,
      file.path(output_figures_dir, "05_total_vs_fresh_metrics.png"),
      width = 11,
      height = 6
    )
  }
  
  total_heatmap_data <- total_metrics %>%
    mutate(treatment = paste(drying, ageing, sep = "_")) %>%
    group_by(treatment) %>%
    summarise(
      across(
        any_of(heatmap_metrics),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = any_of(heatmap_metrics),
      names_to = "metric",
      values_to = "mean_delta"
    ) %>%
    mutate(
      metric_label = dplyr::recode(metric, !!!heatmap_labels),
      metric_label = factor(metric_label, levels = rev(unname(heatmap_labels)))
    )
  
  if (nrow(total_heatmap_data) > 0) {
    p_total_heatmap <- ggplot(
      total_heatmap_data,
      aes(x = treatment, y = metric_label, fill = mean_delta)
    ) +
      geom_tile() +
      labs(
        title = "Mean index changes: aged vs fresh",
        x = "Drying_ageing treatment",
        y = "Index",
        fill = "Mean delta"
      ) +
      theme_bw()
    
    save_plot_if_has_data(
      p_total_heatmap,
      file.path(output_figures_dir, "05_total_vs_fresh_index_heatmap.png"),
      width = 8,
      height = 5
    )
  }
}

message("05_analysis_results.R finished successfully.")