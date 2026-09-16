suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

required_pkgs <- c("lme4", "lmerTest", "emmeans")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running script 06."
  )
}

source(file.path("R", "stats_helpers.R"))

input_dir <- file.path("output", "tables")
output_tables_dir <- file.path("output", "tables")
output_figures_dir <- file.path("output", "figures")

dir.create(output_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures_dir, recursive = TRUE, showWarnings = FALSE)

comparison_files <- tibble(
  comparison_block = c(
    "drying_vs_fresh",
    "ageing_vs_dried",
    "total_vs_fresh"
  ),
  path = c(
    file.path(input_dir, "comparison_drying_vs_fresh.csv"),
    file.path(input_dir, "comparison_ageing_vs_dried.csv"),
    file.path(input_dir, "comparison_total_vs_fresh.csv")
  )
)

primary_responses <- c("rmse", "sam", "iauc")

index_base_names <- c(
  "mfdre",
  "rep",
  "res700_740",
  "mean_blue_450_500",
  "mean_green_500_570",
  "mean_red_650_680",
  "mean760_900",
  "datt",
  "pri",
  "sipi",
  "psri"
)

secondary_responses <- paste0("delta_", index_base_names)

primary_test_adjust_method <- "holm"
primary_contrast_adjust_method <- "holm"
secondary_test_adjust_method <- "BH"
secondary_contrast_adjust_method <- "BH"

comparison_data <- comparison_files %>%
  mutate(
    raw_data = map(path, safe_read_comparison),
    data = map2(raw_data, comparison_block, prepare_block_data),
    block_summary = pmap(
      list(data, comparison_block, path),
      function(data, comparison_block, path) {
        build_block_summary(
          df = data,
          comparison_block = comparison_block,
          path = path,
          primary_responses = primary_responses,
          secondary_responses = secondary_responses
        )
      }
    )
  )

block_summary <- bind_rows(comparison_data$block_summary)

all_model_overview <- empty_model_overview()
all_diagnostics <- empty_diagnostics()

all_tests_primary <- empty_tests()
all_emmeans_primary <- empty_emmeans()
all_contrasts_primary <- empty_contrasts()

all_tests_secondary <- empty_tests()
all_emmeans_secondary <- empty_emmeans()
all_contrasts_secondary <- empty_contrasts()

for (i in seq_len(nrow(comparison_data))) {
  block_name <- comparison_data$comparison_block[[i]]
  df <- comparison_data$data[[i]]
  
  require_individual_random <- block_name %in% c(
    "drying_vs_fresh",
    "ageing_vs_dried",
    "total_vs_fresh"
  )
  
  candidate_terms <- switch(
    block_name,
    drying_vs_fresh = c("drying"),
    ageing_vs_dried = c("drying", "ageing"),
    total_vs_fresh = c("drying", "ageing"),
    character()
  )
  
  primary_available <- intersect(primary_responses, names(df))
  secondary_available <- intersect(secondary_responses, names(df))
  
  primary_out <- run_models_for_family(
    df = df,
    comparison_block = block_name,
    response_family = "primary",
    responses = primary_available,
    candidate_terms = candidate_terms,
    require_individual_random = require_individual_random
  )
  
  secondary_out <- run_models_for_family(
    df = df,
    comparison_block = block_name,
    response_family = "secondary",
    responses = secondary_available,
    candidate_terms = candidate_terms,
    require_individual_random = require_individual_random
  )
  
  all_model_overview <- bind_rows(
    all_model_overview,
    primary_out$model_overview,
    secondary_out$model_overview
  )
  
  all_diagnostics <- bind_rows(
    all_diagnostics,
    primary_out$diagnostics,
    secondary_out$diagnostics
  )
  
  all_tests_primary <- bind_rows(all_tests_primary, primary_out$tests)
  all_emmeans_primary <- bind_rows(all_emmeans_primary, primary_out$emmeans)
  all_contrasts_primary <- bind_rows(all_contrasts_primary, primary_out$contrasts)
  
  all_tests_secondary <- bind_rows(all_tests_secondary, secondary_out$tests)
  all_emmeans_secondary <- bind_rows(all_emmeans_secondary, secondary_out$emmeans)
  all_contrasts_secondary <- bind_rows(all_contrasts_secondary, secondary_out$contrasts)
}

all_tests_primary <- adjust_pvalues(
  all_tests_primary,
  group_cols = c("comparison_block"),
  method = primary_test_adjust_method
)

all_contrasts_primary <- adjust_pvalues(
  all_contrasts_primary,
  group_cols = c("comparison_block", "response", "term"),
  method = primary_contrast_adjust_method
)

all_tests_secondary <- adjust_pvalues(
  all_tests_secondary,
  group_cols = c("comparison_block"),
  method = secondary_test_adjust_method
)

all_contrasts_secondary <- adjust_pvalues(
  all_contrasts_secondary,
  group_cols = c("comparison_block", "response", "term"),
  method = secondary_contrast_adjust_method
)

safe_write_csv(
  block_summary,
  file.path(output_tables_dir, "stats_block_summary.csv")
)

safe_write_csv(
  all_model_overview,
  file.path(output_tables_dir, "stats_model_overview.csv")
)

safe_write_csv(
  all_diagnostics,
  file.path(output_tables_dir, "stats_diagnostics.csv")
)

safe_write_csv(
  all_tests_primary,
  file.path(output_tables_dir, "stats_tests_primary.csv")
)

safe_write_csv(
  all_emmeans_primary,
  file.path(output_tables_dir, "stats_emmeans_primary.csv")
)

safe_write_csv(
  all_contrasts_primary,
  file.path(output_tables_dir, "stats_contrasts_primary.csv")
)

safe_write_csv(
  all_tests_secondary,
  file.path(output_tables_dir, "stats_tests_secondary_indices.csv")
)

safe_write_csv(
  all_emmeans_secondary,
  file.path(output_tables_dir, "stats_emmeans_secondary_indices.csv")
)

safe_write_csv(
  all_contrasts_secondary,
  file.path(output_tables_dir, "stats_contrasts_secondary_indices.csv")
)

save_primary_figures(
  emmeans_df = all_emmeans_primary,
  contrasts_df = all_contrasts_primary,
  output_dir = output_figures_dir
)

message("06_inferential_statistics.R finished successfully.")

message("Blocks summary:")
print(block_summary)

message("Model overview:")
print(
  all_model_overview %>%
    count(comparison_block, response_family, status, name = "n_models")
)

message("Primary tests:")
print(
  all_tests_primary %>%
    select(comparison_block, response, term, statistic, p_value, p_value_adjusted)
)

message("Secondary tests:")
print(
  all_tests_secondary %>%
    select(comparison_block, response, term, statistic, p_value, p_value_adjusted)
)