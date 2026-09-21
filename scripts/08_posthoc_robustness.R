## Post-hoc robustness analyses for the results produced by scripts 05 and 06.
##
## These analyses are sensitivity checks, not part of the primary analysis. They
## answer two questions that the primary models cannot answer on their own:
##
##   1. How much of an observed spectral distance is measurement noise?
##      Every distance metric used here is non-negative, so measurement noise
##      inflates it even when nothing has changed. A split-half resampling of
##      the retained replicates gives an empirical null: the distance between
##      two halves of the same sample, measured in the same session under the
##      same calibration, with no treatment effect between them. The primary
##      models are then refitted on noise-corrected responses.
##
##   2. Does any single inbred line drive the primary results? Each primary
##      model is refitted eight times, dropping one line at a time.
##
## Run after scripts 05 and 06, from the project root.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(rlang)
})

required_pkgs <- c("lme4", "lmerTest", "emmeans")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running script 08."
  )
}

source(file.path("R", "posthoc_helpers.R"))

input_dir <- file.path("output", "tables")
output_tables_dir <- file.path("output", "tables")

dir.create(output_tables_dir, recursive = TRUE, showWarnings = FALSE)

analysis_min_wl <- 400
analysis_max_wl <- 950

## Number of random split-half draws averaged per sample group. Averaging over
## draws removes most of the sampling variability of a single split.
n_splits <- 20L

## Fixed with set.seed so the split-half null is reproducible.
set.seed(20260921)

primary_responses <- c("rmse", "sam", "iauc")

## Metrics for which the quadrature-based noise correction is applied. See
## noise_corrected() in R/posthoc_helpers.R for why iAUC is excluded.
corrected_responses <- c("rmse", "sam")

## Terms that form the multiplicity-correction family. This matches the family
## used in script 06, where `line` enters as a blocking factor and is excluded
## from the exported omnibus tests; keeping the two families identical is what
## makes the raw columns below reproduce the numbers reported in the manuscript.
treatment_terms <- c("drying", "ageing", "drying:ageing")

## The split-half null is the slow step (it reads every retained spectrum and
## resamples it n_splits times). It is cached, and recomputed only when the
## cache is missing, was built with a different number of draws, or is older
## than the retained-spectra table. Set to TRUE to force a recomputation.
force_recompute_null <- FALSE

comparison_blocks <- tibble(
  comparison_block = c("drying_vs_fresh", "ageing_vs_dried", "total_vs_fresh"),
  fixed_formula = c("drying", "drying * ageing", "drying * ageing"),
  path = file.path(
    input_dir,
    c(
      "comparison_drying_vs_fresh.csv",
      "comparison_ageing_vs_dried.csv",
      "comparison_total_vs_fresh.csv"
    )
  )
)

retained_path <- file.path(input_dir, "analysis_retained_spectra.csv")

if (!file.exists(retained_path)) {
  stop("Missing ", retained_path, ". Run script 05 first.")
}

missing_blocks <- comparison_blocks$path[!file.exists(comparison_blocks$path)]
if (length(missing_blocks) > 0) {
  stop(
    "Missing comparison tables: ", paste(missing_blocks, collapse = ", "),
    ". Run script 05 first."
  )
}

# ---------------------------------------------------------------------------
# 1. Replicate-level spread, as reported by the replicate QC step
# ---------------------------------------------------------------------------

replicate_qc_path <- file.path(input_dir, "replicate_qc_results.csv")

if (file.exists(replicate_qc_path)) {
  replicate_noise_summary <- read_csv(replicate_qc_path, show_col_types = FALSE) %>%
    filter(.data$replicate_qc_status == "ok", !is.na(.data$rmse)) %>%
    group_by(.data$timepoint, .data$drying) %>%
    summarise(
      n_replicates = n(),
      rmse_median = stats::median(.data$rmse),
      rmse_q75 = stats::quantile(.data$rmse, 0.75),
      rmse_q95 = stats::quantile(.data$rmse, 0.95),
      rmse_max = max(.data$rmse),
      .groups = "drop"
    )

  write_csv(
    replicate_noise_summary,
    file.path(output_tables_dir, "posthoc_replicate_spread.csv")
  )

  message("Replicate-level spread against the sample median:")
  print(replicate_noise_summary)
}

# ---------------------------------------------------------------------------
# 2. Split-half null distribution of the three divergence metrics
# ---------------------------------------------------------------------------

null_cache_path <- file.path(output_tables_dir, "posthoc_splithalf_null.csv")

cache_is_usable <- function(path) {
  if (force_recompute_null || !file.exists(path)) {
    return(FALSE)
  }
  if (file.mtime(path) < file.mtime(retained_path)) {
    message("Split-half cache is older than the retained spectra; recomputing.")
    return(FALSE)
  }
  cached <- try(read_csv(path, show_col_types = FALSE, n_max = 1), silent = TRUE)
  if (inherits(cached, "try-error") || !"n_splits" %in% names(cached)) {
    return(FALSE)
  }
  if (cached$n_splits[1] != n_splits) {
    message("Split-half cache was built with a different number of draws; recomputing.")
    return(FALSE)
  }
  TRUE
}

if (cache_is_usable(null_cache_path)) {
  message("Reusing the cached split-half null from ", null_cache_path, ".")
  null_by_sample <- read_csv(null_cache_path, show_col_types = FALSE)
} else {
  message("Reading retained spectra for the split-half null ...")

  retained_spectra <- read_csv(retained_path, show_col_types = FALSE) %>%
    filter(
      .data$analysis_keep,
      .data$wavelength >= analysis_min_wl,
      .data$wavelength <= analysis_max_wl
    ) %>%
    select(
      "sample_group", "timepoint", "line", "individual", "drying", "ageing",
      "replicate", "wavelength", "reflectance"
    )

  message("Computing split-half null over ", n_splits, " draws per sample group ...")

  null_by_sample <- split_half_null(retained_spectra, n_splits = n_splits)

  ## Replicates actually used when the representative spectrum was built.
  replicates_used <- retained_spectra %>%
    distinct(.data$sample_group, .data$replicate) %>%
    count(.data$sample_group, name = "n_used")

  null_by_sample <- null_by_sample %>%
    left_join(replicates_used, by = "sample_group") %>%
    mutate(
      across(
        all_of(primary_responses),
        ~ noise_contribution(.x, .data$n_half_a, .data$n_half_b, .data$n_used),
        .names = "contribution_{.col}"
      )
    )

  write_csv(null_by_sample, null_cache_path)
}

## Reported on the scale of a real comparison, assuming two endpoints with the
## same replicate structure, so that the numbers are directly comparable with
## the observed distances in the comparison tables.
noise_floor_summary <- null_by_sample %>%
  group_by(.data$timepoint, .data$drying) %>%
  summarise(
    n_samples = n(),
    across(
      all_of(paste0("contribution_", primary_responses)),
      ~ stats::median(sqrt(2 * .x)),
      .names = "null_{.col}"
    ),
    .groups = "drop"
  ) %>%
  rename_with(~ sub("^null_contribution_", "null_", .x))

write_csv(
  noise_floor_summary,
  file.path(output_tables_dir, "posthoc_noise_floor.csv")
)

message("Measurement-noise floor, on the scale of a pairwise comparison:")
print(noise_floor_summary)

# ---------------------------------------------------------------------------
# 3. Noise-corrected primary models
# ---------------------------------------------------------------------------

## Identifiers of the two endpoints of each comparison, so that the per-sample
## noise contributions can be attached to the comparison rows.
endpoint_contributions <- null_by_sample %>%
  select("sample_group", starts_with("contribution_"))

attach_null <- function(comparisons) {
  comparisons %>%
    left_join(
      endpoint_contributions %>%
        rename_with(~ paste0(.x, "_sample"), starts_with("contribution_")),
      by = c("sample_group_sample" = "sample_group")
    ) %>%
    left_join(
      endpoint_contributions %>%
        rename_with(~ paste0(.x, "_ref"), starts_with("contribution_")),
      by = c("sample_group_ref" = "sample_group")
    ) %>%
    mutate(
      !!!set_names(
        map(
          primary_responses,
          function(metric) {
            expr(
              sqrt(
                .data[[!!paste0("contribution_", metric, "_sample")]] +
                  .data[[!!paste0("contribution_", metric, "_ref")]]
              )
            )
          }
        ),
        paste0("null_", primary_responses)
      )
    ) %>%
    mutate(
      !!!set_names(
        map(
          corrected_responses,
          function(metric) {
            expr(
              noise_corrected(
                .data[[!!metric]],
                .data[[!!paste0("null_", metric)]]
              )
            )
          }
        ),
        paste0(corrected_responses, "_corrected")
      )
    )
}

all_corrected_tests <- list()
all_corrected_emmeans <- list()
all_line_sensitivity <- list()
all_noise_summary <- list()

for (i in seq_len(nrow(comparison_blocks))) {
  block <- comparison_blocks$comparison_block[i]
  fixed_formula <- comparison_blocks$fixed_formula[i]

  message("Block ", block, " ...")

  comparisons <- read_csv(comparison_blocks$path[i], show_col_types = FALSE) %>%
    mutate(individual_id = paste(.data$line, .data$individual, sep = "__")) %>%
    attach_null()

  all_noise_summary[[block]] <- comparisons %>%
    group_by(.data$drying) %>%
    summarise(
      comparison_block = block,
      n = n(),
      across(
        all_of(c(primary_responses, paste0("null_", primary_responses),
                 paste0(corrected_responses, "_corrected"))),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  ## Omnibus tests on the raw and on the noise-corrected responses, so that the
  ## two can be compared directly. Holm correction is applied across all
  ## responses and terms within the block, matching script 06.
  responses_to_test <- c(primary_responses, paste0(corrected_responses, "_corrected"))

  tests <- map_dfr(responses_to_test, function(response) {
    fit_primary_model(comparisons, response, fixed_formula)
  })

  tests <- tests %>% filter(.data$term %in% treatment_terms)

  if (nrow(tests) > 0) {
    tests <- tests %>%
      select(-"dropped_line") %>%
      mutate(
        comparison_block = block,
        response_kind = ifelse(grepl("_corrected$", .data$response), "noise_corrected", "raw"),
        base_response = sub("_corrected$", "", .data$response)
      ) %>%
      group_by(.data$response_kind) %>%
      mutate(
        p_value_adjusted = stats::p.adjust(.data$p_value, method = "holm"),
        p_adjust_method = "holm"
      ) %>%
      ungroup()

    all_corrected_tests[[block]] <- tests

    emm <- map_dfr(responses_to_test, function(response) {
      n_line <- n_distinct(comparisons$line)
      rhs <- if (n_line >= 2) paste(fixed_formula, "+ line") else fixed_formula
      form <- stats::as.formula(paste(response, "~", rhs, "+ (1 | individual_id)"))
      fit <- try(lmerTest::lmer(form, data = comparisons, REML = TRUE), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(tibble())
      }
      terms_to_report <- if (grepl("\\*", fixed_formula)) c("drying", "ageing") else "drying"
      map_dfr(terms_to_report, function(term) {
        as.data.frame(emmeans::emmeans(fit, stats::as.formula(paste("~", term)))) %>%
          as_tibble() %>%
          transmute(
            comparison_block = block,
            response = response,
            response_kind = ifelse(grepl("_corrected$", response), "noise_corrected", "raw"),
            term = term,
            level = as.character(.data[[term]]),
            emmean = .data$emmean,
            se = .data$SE,
            lower_cl = .data$lower.CL,
            upper_cl = .data$upper.CL
          )
      })
    })

    all_corrected_emmeans[[block]] <- emm
  }

  # -------------------------------------------------------------------------
  # 4. Leave-one-line-out sensitivity of the primary models
  # -------------------------------------------------------------------------

  sensitivity <- leave_one_line_out(comparisons, primary_responses, fixed_formula) %>%
    filter(.data$term %in% treatment_terms)

  if (nrow(sensitivity) > 0) {
    all_line_sensitivity[[block]] <- sensitivity %>%
      mutate(
        comparison_block = block,
        dropped_line = ifelse(is.na(.data$dropped_line), "(none)", .data$dropped_line)
      )
  }
}

# ---------------------------------------------------------------------------
# 5. Distance between the drying methods themselves
# ---------------------------------------------------------------------------

sample_spectra_path <- file.path(input_dir, "analysis_sample_spectra.csv")

if (file.exists(sample_spectra_path)) {
  message("Distances between drying methods ...")

  between_drying <- read_csv(sample_spectra_path, show_col_types = FALSE) %>%
    between_drying_distances(analysis_min_wl, analysis_max_wl)

  if (nrow(between_drying) > 0) {
    between_drying <- between_drying %>%
      mutate(individual_id = paste(.data$line, .data$individual, sep = "__"))

    write_csv(
      between_drying,
      file.path(output_tables_dir, "posthoc_between_drying_distances.csv")
    )

    between_drying_summary <- between_drying %>%
      group_by(.data$pair) %>%
      summarise(
        n = n(),
        across(all_of(primary_responses), list(mean = mean, sd = stats::sd)),
        .groups = "drop"
      )

    ## Paired across individuals, so the same random intercept as elsewhere.
    between_drying_tests <- map_dfr(primary_responses, function(response) {
      form <- stats::as.formula(
        paste(response, "~ pair + line + (1 | individual_id)")
      )
      fit <- try(lmerTest::lmer(form, data = between_drying, REML = TRUE), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(tibble())
      }
      emm <- emmeans::emmeans(fit, ~ pair)
      contrasts <- as.data.frame(
        stats::confint(emmeans::contrast(emm, "pairwise", adjust = "holm"))
      )
      tests <- as.data.frame(emmeans::contrast(emm, "pairwise", adjust = "holm"))
      tibble(
        response = response,
        contrast = tests$contrast,
        estimate = tests$estimate,
        se = tests$SE,
        lower_cl = contrasts$lower.CL,
        upper_cl = contrasts$upper.CL,
        p_value_adjusted = tests$p.value,
        p_adjust_method = "holm"
      )
    })

    write_csv(
      between_drying_summary,
      file.path(output_tables_dir, "posthoc_between_drying_summary.csv")
    )
    write_csv(
      between_drying_tests,
      file.path(output_tables_dir, "posthoc_between_drying_tests.csv")
    )

    message("Distance between drying methods (same individual, dried stage):")
    print(between_drying_summary)
    print(between_drying_tests)
  }
}

# ---------------------------------------------------------------------------
# 6. Thermal acceleration of the ageing regimes
# ---------------------------------------------------------------------------
#
# Accelerated ageing standards for paper convert a chamber exposure into an
# equivalent number of years at storage temperature through the Arrhenius
# equation. This section applies that conversion to the three regimes, so that
# the size of the exposure can be stated in the terms collections care uses.
#
# Two limits of the conversion are worth stating with the numbers themselves.
# First, the activation energy is not estimated here: a proper Arrhenius
# extrapolation measures a degradation rate at three or more temperatures and
# fits the slope, which this design does not allow, so a published range for
# cellulose is used instead and the answer is reported across that range.
# Second, the conversion covers only temperature-driven chemistry. It does not
# describe damage caused by repeated swelling and shrinking of the tissue as
# moisture enters and leaves it, which conservation practice treats as a
# fatigue process specified by the amplitude and the number of fluctuations
# rather than as an equivalent age. The humidity regime below is the
# informative case: it ran at the reference temperature, so its thermal
# acceleration is exactly 1 by construction.
#
# Protocol values are typed in from the Methods, as in script 07.

reference_storage_c <- 20

## Published range of activation energies for cellulose degradation.
activation_energies_kj <- c(100, 113, 130)

ramp_h <- 0.5
hold_h <- 3.5
cycle_h <- 8
treatment_days <- 28

ageing_regimes <- list(
  temperature = list(
    label = "Temperature ageing (RH fixed at 50%)",
    blocks = list(
      list(target_c = 15, ramp_h = ramp_h, hold_h = hold_h),
      list(target_c = 40, ramp_h = ramp_h, hold_h = hold_h)
    ),
    rh_amplitude_pp = 0,
    temperature_amplitude_c = 25
  ),
  humidity = list(
    label = "Humidity ageing (temperature fixed at 20 C)",
    blocks = list(
      list(target_c = 20, ramp_h = ramp_h, hold_h = hold_h),
      list(target_c = 20, ramp_h = ramp_h, hold_h = hold_h)
    ),
    rh_amplitude_pp = 50,
    temperature_amplitude_c = 0
  ),
  combined = list(
    label = "Combined ageing (crossed cycle)",
    blocks = list(
      list(target_c = 15, ramp_h = ramp_h, hold_h = hold_h),
      list(target_c = 40, ramp_h = ramp_h, hold_h = hold_h)
    ),
    rh_amplitude_pp = 50,
    temperature_amplitude_c = 25
  )
)

n_cycles <- treatment_days * 24 / cycle_h

arrhenius_summary <- map_dfr(names(ageing_regimes), function(regime) {
  spec <- ageing_regimes[[regime]]
  map_dfr(activation_energies_kj, function(ea_kj) {
    factor <- arrhenius_factor(spec$blocks, ea_kj * 1000, reference_storage_c)
    tibble(
      regime = regime,
      regime_label = spec$label,
      reference_temperature_c = reference_storage_c,
      activation_energy_kj_mol = ea_kj,
      mean_cycle_temperature_c = mean_cycle_temperature(spec$blocks),
      temperature_amplitude_c = spec$temperature_amplitude_c,
      rh_amplitude_pp = spec$rh_amplitude_pp,
      n_cycles = n_cycles,
      treatment_days = treatment_days,
      thermal_acceleration_factor = factor,
      equivalent_days_at_reference = treatment_days * factor,
      equivalent_years_at_reference = treatment_days * factor / 365.25
    )
  })
})

write_csv(
  arrhenius_summary,
  file.path(output_tables_dir, "posthoc_thermal_acceleration.csv")
)

message("Thermal (Arrhenius) acceleration of the ageing regimes:")
print(
  arrhenius_summary %>%
    select(
      "regime", "activation_energy_kj_mol", "thermal_acceleration_factor",
      "equivalent_days_at_reference", "equivalent_years_at_reference"
    )
)

noise_summary <- bind_rows(all_noise_summary)
corrected_tests <- bind_rows(all_corrected_tests)
corrected_emmeans <- bind_rows(all_corrected_emmeans)
line_sensitivity <- bind_rows(all_line_sensitivity)

write_csv(noise_summary, file.path(output_tables_dir, "posthoc_noise_summary_by_drying.csv"))
write_csv(corrected_tests, file.path(output_tables_dir, "posthoc_noise_corrected_tests.csv"))
write_csv(corrected_emmeans, file.path(output_tables_dir, "posthoc_noise_corrected_emmeans.csv"))
write_csv(line_sensitivity, file.path(output_tables_dir, "posthoc_leave_one_line_out.csv"))

message("Observed distances against the noise floor, by drying method:")
print(noise_summary)

message("Omnibus tests, raw versus noise-corrected responses:")
print(
  corrected_tests %>%
    select("comparison_block", "response", "term", "statistic", "p_value", "p_value_adjusted")
)

message("Leave-one-line-out sensitivity of the primary models:")
print(
  line_sensitivity %>%
    select("comparison_block", "response", "term", "dropped_line", "n_obs", "statistic", "p_value")
)
