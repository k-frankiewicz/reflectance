suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

source(file.path("R", "replicate_qc_helpers.R"))

manifest_path <- file.path("data", "metadata", "sample_manifest.csv")
processed_dir <- file.path("data", "processed")
file_qc_path <- file.path("output", "tables", "qc_file_level_results.csv")
output_dir <- file.path("output", "tables")

analysis_min_wl <- 400
analysis_max_wl <- 950
min_n_for_group_threshold <- 20

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manifest_path)) {
  stop(
    "Manifest not found: ",
    normalizePath(manifest_path, winslash = "/", mustWork = FALSE)
  )
}

if (!file.exists(file_qc_path)) {
  stop(
    "File-level QC results not found: ",
    normalizePath(file_qc_path, winslash = "/", mustWork = FALSE)
  )
}

manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
  mutate(file = basename(file))

file_qc <- read_csv(file_qc_path, show_col_types = FALSE) %>%
  mutate(file = basename(file))

eligible_files <- file_qc %>%
  filter(qc_status %in% c("pass", "warn")) %>%
  pull(file) %>%
  unique()

processed_files <- file.path(processed_dir, eligible_files)
processed_files <- processed_files[file.exists(processed_files)]

if (length(processed_files) == 0) {
  stop(
    "No eligible processed files found after filtering by script 03 QC."
  )
}

spectra <- map_dfr(
  processed_files,
  read_one_processed_spectrum,
  analysis_min_wl = analysis_min_wl,
  analysis_max_wl = analysis_max_wl
) %>%
  left_join(manifest, by = "file") %>%
  left_join(
    file_qc %>%
      select(file, qc_status_03 = qc_status, qc_reasons_03 = qc_reasons),
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
    threshold_group = case_when(
      timepoint == "fresh" ~ "fresh",
      timepoint == "dried" ~ paste("dried", drying, sep = "__"),
      timepoint == "aged" ~ paste("aged", drying, ageing, sep = "__"),
      TRUE ~ NA_character_
    )
  )

if (any(is.na(spectra$sample_group))) {
  stop("Could not build sample_group for some files. Check manifest fields.")
}

replicate_counts <- spectra %>%
  distinct(sample_group, file) %>%
  count(sample_group, name = "n_replicates")

median_spectra <- spectra %>%
  group_by(sample_group, wavelength) %>%
  summarise(
    median_reflectance = median(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

replicate_qc <- spectra %>%
  left_join(
    median_spectra,
    by = c("sample_group", "wavelength")
  ) %>%
  group_by(
    sample_group,
    threshold_group,
    timepoint,
    line,
    individual,
    drying,
    ageing,
    file,
    replicate,
    qc_status_03,
    qc_reasons_03
  ) %>%
  summarise(
    n_wavelengths = n(),
    pearson_r = safe_cor(reflectance, median_reflectance),
    rmse = rmse_vec(reflectance, median_reflectance),
    .groups = "drop"
  ) %>%
  left_join(replicate_counts, by = "sample_group") %>%
  mutate(
    assessed_for_replicate_qc = n_replicates >= 2,
    pearson_r = if_else(assessed_for_replicate_qc, pearson_r, NA_real_),
    rmse = if_else(assessed_for_replicate_qc, rmse, NA_real_)
  )

thresholds_detailed_r <- replicate_qc %>%
  filter(assessed_for_replicate_qc, !is.na(timepoint), !is.na(threshold_group)) %>%
  group_by(timepoint, threshold_group) %>%
  summarise(
    metric = "pearson_r",
    threshold_level = "detailed",
    threshold_direction = "low",
    n_metric_values = sum(!is.na(pearson_r)),
    threshold_value = if (sum(!is.na(pearson_r)) >= min_n_for_group_threshold) {
      compute_outlier_threshold(pearson_r, direction = "low")
    } else {
      NA_real_
    },
    .groups = "drop"
  )

thresholds_detailed_rmse <- replicate_qc %>%
  filter(assessed_for_replicate_qc, !is.na(timepoint), !is.na(threshold_group)) %>%
  group_by(timepoint, threshold_group) %>%
  summarise(
    metric = "rmse",
    threshold_level = "detailed",
    threshold_direction = "high",
    n_metric_values = sum(!is.na(rmse)),
    threshold_value = if (sum(!is.na(rmse)) >= min_n_for_group_threshold) {
      compute_outlier_threshold(rmse, direction = "high")
    } else {
      NA_real_
    },
    .groups = "drop"
  )

thresholds_timepoint_r <- replicate_qc %>%
  filter(assessed_for_replicate_qc, !is.na(timepoint)) %>%
  group_by(timepoint) %>%
  summarise(
    metric = "pearson_r",
    threshold_level = "timepoint",
    threshold_direction = "low",
    n_metric_values = sum(!is.na(pearson_r)),
    threshold_value = if (sum(!is.na(pearson_r)) >= min_n_for_group_threshold) {
      compute_outlier_threshold(pearson_r, direction = "low")
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(threshold_group = NA_character_) %>%
  select(
    timepoint,
    threshold_group,
    metric,
    threshold_level,
    threshold_direction,
    n_metric_values,
    threshold_value
  )

thresholds_timepoint_rmse <- replicate_qc %>%
  filter(assessed_for_replicate_qc, !is.na(timepoint)) %>%
  group_by(timepoint) %>%
  summarise(
    metric = "rmse",
    threshold_level = "timepoint",
    threshold_direction = "high",
    n_metric_values = sum(!is.na(rmse)),
    threshold_value = if (sum(!is.na(rmse)) >= min_n_for_group_threshold) {
      compute_outlier_threshold(rmse, direction = "high")
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(threshold_group = NA_character_) %>%
  select(
    timepoint,
    threshold_group,
    metric,
    threshold_level,
    threshold_direction,
    n_metric_values,
    threshold_value
  )

thresholds <- bind_rows(
  thresholds_detailed_r,
  thresholds_detailed_rmse,
  thresholds_timepoint_r,
  thresholds_timepoint_rmse
)

r_thresholds_detailed <- thresholds %>%
  filter(metric == "pearson_r", threshold_level == "detailed") %>%
  select(
    timepoint,
    threshold_group,
    r_threshold_detailed = threshold_value,
    r_n_detailed = n_metric_values
  )

rmse_thresholds_detailed <- thresholds %>%
  filter(metric == "rmse", threshold_level == "detailed") %>%
  select(
    timepoint,
    threshold_group,
    rmse_threshold_detailed = threshold_value,
    rmse_n_detailed = n_metric_values
  )

r_thresholds_timepoint <- thresholds %>%
  filter(metric == "pearson_r", threshold_level == "timepoint") %>%
  select(
    timepoint,
    r_threshold_timepoint = threshold_value,
    r_n_timepoint = n_metric_values
  )

rmse_thresholds_timepoint <- thresholds %>%
  filter(metric == "rmse", threshold_level == "timepoint") %>%
  select(
    timepoint,
    rmse_threshold_timepoint = threshold_value,
    rmse_n_timepoint = n_metric_values
  )

replicate_qc <- replicate_qc %>%
  left_join(
    r_thresholds_detailed,
    by = c("timepoint", "threshold_group")
  ) %>%
  left_join(
    rmse_thresholds_detailed,
    by = c("timepoint", "threshold_group")
  ) %>%
  left_join(
    r_thresholds_timepoint,
    by = "timepoint"
  ) %>%
  left_join(
    rmse_thresholds_timepoint,
    by = "timepoint"
  ) %>%
  mutate(
    r_threshold = coalesce(r_threshold_detailed, r_threshold_timepoint),
    rmse_threshold = coalesce(rmse_threshold_detailed, rmse_threshold_timepoint),
    r_threshold_source = case_when(
      !is.na(r_threshold_detailed) ~ "detailed",
      !is.na(r_threshold_timepoint) ~ "timepoint",
      TRUE ~ NA_character_
    ),
    rmse_threshold_source = case_when(
      !is.na(rmse_threshold_detailed) ~ "detailed",
      !is.na(rmse_threshold_timepoint) ~ "timepoint",
      TRUE ~ NA_character_
    ),
    low_r_flag = case_when(
      !assessed_for_replicate_qc ~ NA,
      is.na(r_threshold) ~ NA,
      is.na(pearson_r) ~ NA,
      TRUE ~ pearson_r < r_threshold
    ),
    high_rmse_flag = case_when(
      !assessed_for_replicate_qc ~ NA,
      is.na(rmse_threshold) ~ NA,
      is.na(rmse) ~ NA,
      TRUE ~ rmse > rmse_threshold
    ),
    outlier_candidate = case_when(
      !assessed_for_replicate_qc ~ NA,
      is.na(low_r_flag) | is.na(high_rmse_flag) ~ NA,
      TRUE ~ low_r_flag & high_rmse_flag
    ),
    replicate_qc_status = case_when(
      !assessed_for_replicate_qc ~ "not_assessed",
      is.na(outlier_candidate) ~ "thresholds_unavailable",
      outlier_candidate ~ "outlier_candidate",
      TRUE ~ "ok"
    )
  )

replicate_qc_summary <- replicate_qc %>%
  count(replicate_qc_status, name = "n_files")

replicate_qc_summary_by_timepoint <- replicate_qc %>%
  count(timepoint, replicate_qc_status, name = "n_files") %>%
  arrange(timepoint, replicate_qc_status)

sample_qc_summary <- replicate_qc %>%
  group_by(
    sample_group,
    timepoint,
    line,
    individual,
    drying,
    ageing
  ) %>%
  summarise(
    n_replicates = first(n_replicates),
    n_outlier_candidates = sum(outlier_candidate %in% TRUE, na.rm = TRUE),
    mean_pearson_r = mean(pearson_r, na.rm = TRUE),
    median_pearson_r = median(pearson_r, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    median_rmse = median(rmse, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  replicate_qc,
  file.path(output_dir, "replicate_qc_results.csv")
)

write_csv(
  replicate_qc_summary,
  file.path(output_dir, "replicate_qc_summary.csv")
)

write_csv(
  replicate_qc_summary_by_timepoint,
  file.path(output_dir, "replicate_qc_summary_by_timepoint.csv")
)

write_csv(
  sample_qc_summary,
  file.path(output_dir, "sample_qc_summary.csv")
)

write_csv(
  thresholds,
  file.path(output_dir, "replicate_qc_thresholds.csv")
)

message("Replicate QC summary:")
print(replicate_qc_summary)

message("Replicate QC thresholds:")
print(thresholds)

if (interactive()) {
  View(replicate_qc_summary)
  View(replicate_qc)
  View(sample_qc_summary)
  View(thresholds)
}