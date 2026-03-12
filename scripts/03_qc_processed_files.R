suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

source(file.path("R", "qc_helpers.R"))

manifest_path <- file.path("data", "metadata", "sample_manifest.csv")
processed_dir <- file.path("data", "processed")
output_dir <- file.path("output", "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

processed_files <- list.files(
  processed_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (!file.exists(manifest_path)) {
  stop(
    "Manifest not found: ",
    normalizePath(manifest_path, winslash = "/", mustWork = FALSE)
  )
}

if (length(processed_files) == 0) {
  stop(
    "No processed .csv files found: ",
    normalizePath(processed_dir, winslash = "/", mustWork = FALSE)
  )
}

manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
  mutate(file = basename(file))

qc_results <- map_dfr(processed_files, assess_file)

roughness_reference <- qc_results %>%
  filter(qc_status != "fail", !is.na(roughness)) %>%
  pull(roughness)

if (length(roughness_reference) >= 4) {
  q1 <- as.numeric(quantile(roughness_reference, 0.25, na.rm = TRUE))
  q3 <- as.numeric(quantile(roughness_reference, 0.75, na.rm = TRUE))
  iqr_val <- q3 - q1
  roughness_threshold <- q3 + 3 * iqr_val
  
  qc_results <- qc_results %>%
    mutate(
      warn_reasons = if_else(
        !is.na(roughness) & roughness > roughness_threshold,
        if_else(
          is.na(warn_reasons) | warn_reasons == "",
          "high_roughness",
          paste(warn_reasons, "high_roughness", sep = "; ")
        ),
        warn_reasons
      )
    )
}

qc_results <- qc_results %>%
  mutate(
    qc_status = case_when(
      !is.na(fail_reasons) & fail_reasons != "" ~ "fail",
      !is.na(warn_reasons) & warn_reasons != "" ~ "warn",
      TRUE ~ "pass"
    ),
    qc_reasons = map2_chr(fail_reasons, warn_reasons, combine_two_reasons)
  ) %>%
  left_join(manifest, by = "file") %>%
  relocate(
    file, line, individual, replicate, timepoint, drying, ageing, notes,
    qc_status, qc_reasons
  )

qc_summary <- qc_results %>%
  count(qc_status, name = "n_files") %>%
  mutate(qc_status = factor(qc_status, levels = c("pass", "warn", "fail"))) %>%
  arrange(qc_status)

qc_reasons_summary <- qc_results %>%
  filter(!is.na(qc_reasons), qc_reasons != "") %>%
  separate_rows(qc_reasons, sep = ";\\s*") %>%
  count(qc_reasons, sort = TRUE, name = "n_files")

qc_summary_by_timepoint <- qc_results %>%
  count(timepoint, qc_status, name = "n_files") %>%
  arrange(timepoint, qc_status)

write_csv(
  qc_results,
  file.path(output_dir, "qc_file_level_results.csv")
)

write_csv(
  qc_summary,
  file.path(output_dir, "qc_file_level_summary.csv")
)

write_csv(
  qc_reasons_summary,
  file.path(output_dir, "qc_file_level_reasons_summary.csv")
)

write_csv(
  qc_summary_by_timepoint,
  file.path(output_dir, "qc_file_level_summary_by_timepoint.csv")
)

message("QC summary:")
print(qc_summary)

message("QC reasons:")
print(qc_reasons_summary)

if (interactive()) {
  View(qc_summary)
  View(qc_results)
  View(qc_reasons_summary)
  View(qc_summary_by_timepoint)
}