## Helper functions for post-hoc robustness analyses (script 08).
##
## Two families of helpers are provided:
##   * measurement-noise helpers, which estimate how much of an observed
##     spectral distance is attributable to measurement noise rather than to
##     treatment, using a split-half resampling of the retained replicates;
##   * sensitivity helpers, which refit the primary models while dropping one
##     inbred line at a time.

spectral_angle <- function(a, b) {
  denom <- sqrt(sum(a^2) * sum(b^2))
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  acos(pmin(1, pmax(-1, sum(a * b) / denom)))
}

integrated_abs_area <- function(wavelength, delta) {
  n <- length(wavelength)
  if (n < 2) {
    return(NA_real_)
  }
  d <- abs(delta)
  sum((d[-1] + d[-n]) / 2 * diff(wavelength))
}

## One split-half draw for a single sample group: the retained replicates are
## split at random into two disjoint halves, each half is aggregated with the
## same median-per-wavelength rule used for the representative spectra, and the
## two half-spectra are compared with the same three metrics as the real
## pairwise comparisons. Both halves come from the same sample, the same
## session and the same calibration, so the resulting distance is measurement
## noise plus within-leaf heterogeneity, with no treatment effect.
split_half_once <- function(df) {
  reps <- unique(df$replicate)
  n_rep <- length(reps)
  if (n_rep < 4) {
    return(NULL)
  }

  shuffled <- sample(reps)
  size_a <- floor(n_rep / 2)
  half_a <- shuffled[seq_len(size_a)]

  wide <- df %>%
    dplyr::mutate(half = ifelse(.data$replicate %in% half_a, "a", "b")) %>%
    dplyr::group_by(.data$half, .data$wavelength) %>%
    dplyr::summarise(reflectance = stats::median(.data$reflectance), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = "half", values_from = "reflectance")

  wide <- wide[stats::complete.cases(wide), ]
  if (nrow(wide) < 2) {
    return(NULL)
  }

  tibble::tibble(
    rmse = sqrt(mean((wide$a - wide$b)^2)),
    sam = spectral_angle(wide$a, wide$b),
    iauc = integrated_abs_area(wide$wavelength, wide$a - wide$b),
    n_replicates = n_rep,
    n_half_a = size_a,
    n_half_b = n_rep - size_a
  )
}

split_half_null <- function(spectra, n_splits = 20L) {
  spectra %>%
    dplyr::group_by(
      .data$sample_group, .data$timepoint, .data$line, .data$individual,
      .data$drying, .data$ageing
    ) %>%
    dplyr::group_modify(~ {
      draws <- purrr::map_dfr(seq_len(n_splits), function(i) split_half_once(.x))
      if (nrow(draws) == 0) {
        return(tibble::tibble())
      }
      dplyr::summarise(
        draws,
        dplyr::across(c("rmse", "sam", "iauc"), mean),
        n_replicates = dplyr::first(.data$n_replicates),
        n_half_a = dplyr::first(.data$n_half_a),
        n_half_b = dplyr::first(.data$n_half_b),
        n_splits = dplyr::n()
      )
    }) %>%
    dplyr::ungroup()
}

## Rescale a split-half distance to the geometry of a real pairwise comparison.
##
## A split-half distance compares a median of n_half_a replicates with a median
## of n_half_b replicates, whereas a real comparison uses the full replicate set
## on each side. All three metrics scale with the standard deviation of the
## aggregated spectrum, so the per-sample noise contribution is
##   metric_split^2 / ((1 / n_half_a + 1 / n_half_b) * n_replicates_used)
## and the null for a comparison is the sum of the contributions of its two
## endpoints.
noise_contribution <- function(metric_split, n_half_a, n_half_b, n_used) {
  metric_split^2 / ((1 / n_half_a + 1 / n_half_b) * n_used)
}

## Remove the measurement-noise component from an observed distance. RMSE and
## SAM combine signal and noise approximately in quadrature, so the corrected
## value is sqrt(observed^2 - null^2), floored at zero. This is a moment-based
## approximation and is intended as a sensitivity check, not as a replacement
## for the primary analysis. It is deliberately not applied to iAUC, whose
## absolute-area construction does not combine in quadrature; for iAUC only the
## null is reported.
noise_corrected <- function(observed, null) {
  sqrt(pmax(0, observed^2 - null^2))
}

## Spectral distance between the drying methods themselves, within an
## individual. The comparison blocks of script 05 measure each dried sample
## against its own fresh or dried reference, and therefore never say how far
## the drying methods sit from one another. That distance is what a trait model
## has to cross when it is trained on material prepared one way and applied to
## material prepared another way.
between_drying_distances <- function(sample_spectra, min_wl, max_wl) {
  wide <- sample_spectra %>%
    dplyr::filter(
      .data$timepoint == "dried",
      .data$wavelength >= min_wl,
      .data$wavelength <= max_wl
    ) %>%
    dplyr::select("line", "individual", "drying", "wavelength", "reflectance") %>%
    tidyr::pivot_wider(names_from = "drying", values_from = "reflectance")

  methods <- intersect(c("P", "C", "L"), names(wide))
  if (length(methods) < 2) {
    return(tibble::tibble())
  }

  wide <- wide[stats::complete.cases(wide), ]

  pairs <- utils::combn(methods, 2, simplify = FALSE)

  wide %>%
    dplyr::group_by(.data$line, .data$individual) %>%
    dplyr::group_modify(~ {
      purrr::map_dfr(pairs, function(p) {
        a <- .x[[p[1]]]
        b <- .x[[p[2]]]
        tibble::tibble(
          drying_a = p[1],
          drying_b = p[2],
          pair = paste(p[1], p[2], sep = "_vs_"),
          rmse = sqrt(mean((a - b)^2)),
          sam = spectral_angle(a, b),
          iauc = integrated_abs_area(.x$wavelength, a - b)
        )
      })
    }) %>%
    dplyr::ungroup()
}

## Thermal acceleration of an ageing regime, in the sense used by accelerated
## ageing standards for paper (ASTM D6819, ISO 5630): the rate of a
## temperature-driven degradation reaction is assumed to follow the Arrhenius
## equation, and the exposure is expressed as the equivalent time at a
## reference storage temperature.
##
## The factor is integrated numerically over one programmed cycle, ramps
## included, because the Arrhenius factor is convex in temperature and the rate
## at the mean temperature of a ramp is not the mean rate over that ramp.
##
## `blocks` is a list of blocks, each with `target_c`, `ramp_h` and `hold_h`.
## The regime is assumed to be in steady-state cycling, so the ramp of the
## first block starts from the target of the last one.
arrhenius_factor <- function(blocks, activation_energy_j, reference_c, steps = 200L) {
  gas_constant <- 8.314462618
  ref_k <- reference_c + 273.15

  rate_ratio <- function(temp_c) {
    exp(-activation_energy_j / gas_constant * (1 / (temp_c + 273.15) - 1 / ref_k))
  }

  previous_target <- blocks[[length(blocks)]]$target_c

  weighted <- vapply(blocks, function(b) {
    ## Linear ramp from the previous hold temperature to this block's target.
    ramp_temps <- seq(previous_target, b$target_c, length.out = steps)
    ramp_rate <- mean(rate_ratio(ramp_temps))
    previous_target <<- b$target_c
    ramp_rate * b$ramp_h + rate_ratio(b$target_c) * b$hold_h
  }, numeric(1))

  total_h <- sum(vapply(blocks, function(b) b$ramp_h + b$hold_h, numeric(1)))
  sum(weighted) / total_h
}

mean_cycle_temperature <- function(blocks, steps = 200L) {
  previous_target <- blocks[[length(blocks)]]$target_c
  weighted <- vapply(blocks, function(b) {
    ramp_mean <- mean(seq(previous_target, b$target_c, length.out = steps))
    previous_target <<- b$target_c
    ramp_mean * b$ramp_h + b$target_c * b$hold_h
  }, numeric(1))
  total_h <- sum(vapply(blocks, function(b) b$ramp_h + b$hold_h, numeric(1)))
  sum(weighted) / total_h
}

## Directional consistency of a signed change across inbred lines.
##
## The primary metrics are non-negative distances, so only the signed index
## changes (`delta_*`) carry a direction. For each of them this reports the
## dominant direction, the share of comparisons that follow it overall and
## within the least consistent line, and how many lines agree on it. A change
## that keeps one direction in every genetic background is a bias that a model
## can absorb; one whose sign depends on the line is not.
direction_consistency <- function(comparisons, responses, block) {
  purrr::map_dfr(responses, function(response) {
    values <- comparisons[[response]]
    if (is.null(values) || all(is.na(values))) {
      return(tibble::tibble())
    }

    overall_mean <- mean(values, na.rm = TRUE)
    dominant <- sign(overall_mean)
    if (dominant == 0) {
      return(tibble::tibble())
    }

    follows <- sign(values) == dominant

    per_line <- tapply(
      seq_along(values), comparisons$line,
      function(idx) {
        v <- values[idx]
        v <- v[!is.na(v)]
        c(share = mean(sign(v) == dominant), line_mean = mean(v))
      }
    )
    shares <- vapply(per_line, function(z) z[["share"]], numeric(1))
    line_means <- vapply(per_line, function(z) z[["line_mean"]], numeric(1))

    tibble::tibble(
      comparison_block = block,
      response = response,
      overall_mean = overall_mean,
      dominant_direction = if (dominant > 0) "increase" else "decrease",
      share_following_overall = mean(follows, na.rm = TRUE),
      n_lines = length(shares),
      share_following_min = min(shares),
      share_following_max = max(shares),
      n_lines_agreeing = sum(sign(line_means) == dominant)
    )
  })
}

## Refit one primary model and return the omnibus tests, optionally after
## dropping a single inbred line.
fit_primary_model <- function(data, response, fixed_formula, drop_line = NA_character_) {
  if (!is.na(drop_line)) {
    data <- data[data$line != drop_line, , drop = FALSE]
  }
  if (nrow(data) == 0 || all(is.na(data[[response]]))) {
    return(tibble::tibble())
  }

  n_line <- dplyr::n_distinct(data$line)
  rhs <- if (n_line >= 2) paste(fixed_formula, "+ line") else fixed_formula
  form <- stats::as.formula(
    paste(response, "~", rhs, "+ (1 | individual_id)")
  )

  fit <- try(
    lmerTest::lmer(form, data = data, REML = TRUE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(tibble::tibble())
  }

  ## joint_tests(), not anova(), so that these omnibus tests are produced by
  ## exactly the same machinery as in script 06 and the uncorrected results
  ## below reproduce the values reported in the manuscript. The two differ
  ## slightly for the interaction term under this unbalanced design.
  tests <- try(emmeans::joint_tests(fit), silent = TRUE)
  if (inherits(tests, "try-error")) {
    return(tibble::tibble())
  }

  tests <- as.data.frame(tests)
  names(tests) <- make.names(names(tests))

  tibble::tibble(
    dropped_line = drop_line,
    response = response,
    term = as.character(tests[["model.term"]]),
    n_obs = nrow(data),
    n_line = n_line,
    df1 = as.numeric(tests[["df1"]]),
    df2 = as.numeric(tests[["df2"]]),
    statistic = as.numeric(tests[["F.ratio"]]),
    p_value = as.numeric(tests[["p.value"]]),
    singular_fit = lme4::isSingular(fit)
  )
}

leave_one_line_out <- function(data, responses, fixed_formula) {
  lines <- sort(unique(data$line))
  purrr::map_dfr(responses, function(response) {
    purrr::map_dfr(c(NA_character_, lines), function(drop_line) {
      fit_primary_model(data, response, fixed_formula, drop_line)
    })
  })
}
