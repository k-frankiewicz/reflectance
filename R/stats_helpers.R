empty_block_summary <- function() {
  tibble(
    comparison_block = character(),
    input_file = character(),
    file_exists = logical(),
    n_rows = integer(),
    n_lines = integer(),
    n_individuals = integer(),
    n_drying_levels = integer(),
    n_ageing_levels = integer(),
    primary_responses_available = character(),
    secondary_responses_available = character()
  )
}

empty_model_overview <- function() {
  tibble(
    comparison_block = character(),
    response_family = character(),
    response = character(),
    status = character(),
    message = character(),
    n_obs = integer(),
    n_line = integer(),
    n_individual = integer(),
    active_fixed_terms = character(),
    engine = character(),
    formula = character(),
    singular_fit = logical()
  )
}

empty_diagnostics <- function() {
  tibble(
    comparison_block = character(),
    response_family = character(),
    response = character(),
    engine = character(),
    formula = character(),
    singular_fit = logical(),
    n_obs = integer(),
    n_line = integer(),
    n_individual = integer(),
    sigma = double(),
    aic = double(),
    bic = double(),
    logLik = double(),
    convergence_message = character()
  )
}

empty_tests <- function() {
  tibble(
    comparison_block = character(),
    response_family = character(),
    response = character(),
    engine = character(),
    formula = character(),
    term = character(),
    df1 = double(),
    df2 = double(),
    statistic = double(),
    p_value = double(),
    p_value_adjusted = double(),
    p_adjust_method = character()
  )
}

empty_emmeans <- function() {
  tibble(
    comparison_block = character(),
    response_family = character(),
    response = character(),
    engine = character(),
    formula = character(),
    term = character(),
    level = character(),
    emmean = double(),
    se = double(),
    df = double(),
    lower_cl = double(),
    upper_cl = double()
  )
}

empty_contrasts <- function() {
  tibble(
    comparison_block = character(),
    response_family = character(),
    response = character(),
    engine = character(),
    formula = character(),
    term = character(),
    contrast = character(),
    estimate = double(),
    se = double(),
    df = double(),
    statistic = double(),
    lower_cl = double(),
    upper_cl = double(),
    p_value = double(),
    p_value_adjusted = double(),
    p_adjust_method = character()
  )
}

safe_read_comparison <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  
  readr::read_csv(path, show_col_types = FALSE)
}

count_non_missing_levels <- function(x) {
  dplyr::n_distinct(x[!is.na(x)])
}

safe_char <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }
  
  paste(stats::na.omit(as.character(x)), collapse = " | ")
}

prepare_block_data <- function(df, comparison_block) {
  if (nrow(df) == 0) {
    return(tibble())
  }
  
  df %>%
    dplyr::mutate(
      comparison_block = comparison_block,
      line = if ("line" %in% names(.)) as.character(line) else NA_character_,
      individual = if ("individual" %in% names(.)) as.character(individual) else NA_character_,
      drying = if ("drying" %in% names(.)) as.character(drying) else NA_character_,
      ageing = if ("ageing" %in% names(.)) as.character(ageing) else NA_character_,
      individual_id = interaction(line, individual, drop = TRUE, sep = "__")
    )
}

build_block_summary <- function(
    df,
    comparison_block,
    path,
    primary_responses,
    secondary_responses
) {
  primary_available <- intersect(primary_responses, names(df))
  secondary_available <- intersect(secondary_responses, names(df))
  
  tibble(
    comparison_block = comparison_block,
    input_file = path,
    file_exists = file.exists(path),
    n_rows = nrow(df),
    n_lines = if ("line" %in% names(df)) count_non_missing_levels(df$line) else 0L,
    n_individuals = if (all(c("line", "individual") %in% names(df))) {
      count_non_missing_levels(interaction(df$line, df$individual, drop = TRUE, sep = "__"))
    } else {
      0L
    },
    n_drying_levels = if ("drying" %in% names(df)) count_non_missing_levels(df$drying) else 0L,
    n_ageing_levels = if ("ageing" %in% names(df)) count_non_missing_levels(df$ageing) else 0L,
    primary_responses_available = if (length(primary_available) == 0) NA_character_ else paste(primary_available, collapse = "; "),
    secondary_responses_available = if (length(secondary_available) == 0) NA_character_ else paste(secondary_available, collapse = "; ")
  )
}

active_fixed_terms <- function(data, candidate_terms) {
  candidate_terms[
    vapply(
      candidate_terms,
      function(term) {
        term %in% names(data) && count_non_missing_levels(data[[term]]) >= 2
      },
      logical(1)
    )
  ]
}

candidate_formula_strings <- function(
    response,
    fixed_terms,
    data,
    require_individual_random = FALSE
) {
  fixed_part <- paste(fixed_terms, collapse = " + ")
  
  n_line <- if ("line" %in% names(data)) {
    count_non_missing_levels(data$line)
  } else {
    0L
  }
  
  n_individual <- if ("individual_id" %in% names(data)) {
    count_non_missing_levels(data$individual_id)
  } else {
    0L
  }
  
  formulas <- character()
  
  if (require_individual_random) {
    if (n_line >= 2 && n_individual >= 2) {
      formulas <- c(
        formulas,
        paste(response, "~", fixed_part, "+ (1 | line) + (1 | individual_id)")
      )
    }
    
    if (n_individual >= 2) {
      formulas <- c(
        formulas,
        paste(response, "~", fixed_part, "+ (1 | individual_id)")
      )
    }
    
    return(unique(formulas))
  }
  
  if (n_line >= 2 && n_individual >= 2) {
    formulas <- c(
      formulas,
      paste(response, "~", fixed_part, "+ (1 | line) + (1 | individual_id)")
    )
  }
  
  if (n_individual >= 2) {
    formulas <- c(
      formulas,
      paste(response, "~", fixed_part, "+ (1 | individual_id)")
    )
  }
  
  if (n_line >= 2) {
    formulas <- c(
      formulas,
      paste(response, "~", fixed_part, "+ (1 | line)")
    )
  }
  
  formulas <- c(formulas, paste(response, "~", fixed_part))
  
  unique(formulas)
}

safe_fit_formula <- function(formula_string, data) {
  form <- stats::as.formula(formula_string)
  is_mixed <- grepl("\\|", formula_string)
  
  if (is_mixed) {
    fit <- tryCatch(
      lmerTest::lmer(form, data = data, REML = TRUE),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) {
      return(list(
        success = FALSE,
        fit = NULL,
        engine = "lmer",
        error_message = conditionMessage(fit),
        singular_fit = NA
      ))
    }
    
    singular <- tryCatch(
      lme4::isSingular(fit, tol = 1e-5),
      error = function(e) NA
    )
    
    return(list(
      success = TRUE,
      fit = fit,
      engine = "lmer",
      error_message = NA_character_,
      singular_fit = singular
    ))
  }
  
  fit <- tryCatch(
    stats::lm(form, data = data),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      success = FALSE,
      fit = NULL,
      engine = "lm",
      error_message = conditionMessage(fit),
      singular_fit = FALSE
    ))
  }
  
  list(
    success = TRUE,
    fit = fit,
    engine = "lm",
    error_message = NA_character_,
    singular_fit = FALSE
  )
}

fit_best_model <- function(
    data,
    response,
    candidate_terms,
    require_individual_random = FALSE
) {
  needed <- unique(c(response, candidate_terms, "line", "individual_id"))
  needed <- needed[needed %in% names(data)]
  
  dat <- data %>%
    dplyr::select(dplyr::any_of(needed)) %>%
    dplyr::filter(!is.na(.data[[response]]))
  
  if ("line" %in% names(dat)) {
    dat <- dat %>% dplyr::filter(!is.na(line), line != "")
  }
  
  if ("individual_id" %in% names(dat)) {
    dat <- dat %>% dplyr::filter(!is.na(individual_id), individual_id != "")
  }
  
  if (nrow(dat) < 3) {
    return(list(
      status = "skipped",
      message = "Too few observations after filtering.",
      fit = NULL,
      data = dat,
      engine = NA_character_,
      formula = NA_character_,
      singular_fit = NA,
      active_terms = character()
    ))
  }
  
  active_terms <- active_fixed_terms(dat, candidate_terms)
  
  if (length(active_terms) == 0) {
    return(list(
      status = "skipped",
      message = "No fixed effect with at least two observed levels.",
      fit = NULL,
      data = dat,
      engine = NA_character_,
      formula = NA_character_,
      singular_fit = NA,
      active_terms = character()
    ))
  }
  
  dat_model <- dat %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(response, active_terms)), ~ !is.na(.x)))
  
  if (nrow(dat_model) < 3) {
    return(list(
      status = "skipped",
      message = "Too few complete cases for model fitting.",
      fit = NULL,
      data = dat_model,
      engine = NA_character_,
      formula = NA_character_,
      singular_fit = NA,
      active_terms = active_terms
    ))
  }
  
  if (require_individual_random) {
    n_individual <- if ("individual_id" %in% names(dat_model)) {
      count_non_missing_levels(dat_model$individual_id)
    } else {
      0L
    }
    
    if (n_individual < 2) {
      return(list(
        status = "skipped",
        message = "Too few individual_id levels for models requiring random intercept for individual_id.",
        fit = NULL,
        data = dat_model,
        engine = NA_character_,
        formula = NA_character_,
        singular_fit = NA,
        active_terms = active_terms
      ))
    }
  }
  
  formula_strings <- candidate_formula_strings(
    response = response,
    fixed_terms = active_terms,
    data = dat_model,
    require_individual_random = require_individual_random
  )
  
  first_success <- NULL
  first_non_singular <- NULL
  
  for (formula_string in formula_strings) {
    one_fit <- safe_fit_formula(formula_string, dat_model)
    
    if (!isTRUE(one_fit$success)) {
      next
    }
    
    candidate <- list(
      status = "fitted",
      message = NA_character_,
      fit = one_fit$fit,
      data = dat_model,
      engine = one_fit$engine,
      formula = formula_string,
      singular_fit = isTRUE(one_fit$singular_fit),
      active_terms = active_terms
    )
    
    if (is.null(first_success)) {
      first_success <- candidate
    }
    
    if (!isTRUE(one_fit$singular_fit)) {
      first_non_singular <- candidate
      break
    }
  }
  
  if (!is.null(first_non_singular)) {
    return(first_non_singular)
  }
  
  if (!is.null(first_success)) {
    return(first_success)
  }
  
  list(
    status = "failed",
    message = "All candidate models failed.",
    fit = NULL,
    data = dat_model,
    engine = NA_character_,
    formula = NA_character_,
    singular_fit = NA,
    active_terms = active_terms
  )
}

extract_model_overview_row <- function(fit_result, comparison_block, response_family, response) {
  dat <- fit_result$data
  
  tibble(
    comparison_block = comparison_block,
    response_family = response_family,
    response = response,
    status = fit_result$status,
    message = fit_result$message,
    n_obs = if (!is.null(dat)) nrow(dat) else NA_integer_,
    n_line = if (!is.null(dat) && "line" %in% names(dat)) count_non_missing_levels(dat$line) else NA_integer_,
    n_individual = if (!is.null(dat) && "individual_id" %in% names(dat)) count_non_missing_levels(dat$individual_id) else NA_integer_,
    active_fixed_terms = if (length(fit_result$active_terms) == 0) NA_character_ else paste(fit_result$active_terms, collapse = "; "),
    engine = fit_result$engine,
    formula = fit_result$formula,
    singular_fit = fit_result$singular_fit
  )
}

extract_diagnostics_row <- function(fit_result, comparison_block, response_family, response) {
  fit <- fit_result$fit
  dat <- fit_result$data
  
  if (is.null(fit)) {
    return(empty_diagnostics())
  }
  
  sigma_val <- tryCatch(stats::sigma(fit), error = function(e) NA_real_)
  aic_val <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
  bic_val <- tryCatch(stats::BIC(fit), error = function(e) NA_real_)
  logLik_val <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
  
  conv_msg <- NA_character_
  if (inherits(fit, "merMod")) {
    conv_raw <- tryCatch(
      fit@optinfo$conv$lme4$messages,
      error = function(e) NULL
    )
    conv_msg <- safe_char(conv_raw)
  }
  
  tibble(
    comparison_block = comparison_block,
    response_family = response_family,
    response = response,
    engine = fit_result$engine,
    formula = fit_result$formula,
    singular_fit = fit_result$singular_fit,
    n_obs = nrow(dat),
    n_line = if ("line" %in% names(dat)) count_non_missing_levels(dat$line) else NA_integer_,
    n_individual = if ("individual_id" %in% names(dat)) count_non_missing_levels(dat$individual_id) else NA_integer_,
    sigma = sigma_val,
    aic = aic_val,
    bic = bic_val,
    logLik = logLik_val,
    convergence_message = conv_msg
  )
}

extract_joint_tests <- function(fit_result, comparison_block, response_family, response) {
  fit <- fit_result$fit
  
  if (is.null(fit)) {
    return(empty_tests())
  }
  
  jt <- tryCatch(
    emmeans::joint_tests(fit),
    error = function(e) NULL
  )
  
  if (is.null(jt)) {
    return(empty_tests())
  }
  
  jt_df <- as.data.frame(jt)
  if (nrow(jt_df) == 0) {
    return(empty_tests())
  }
  
  names(jt_df) <- make.names(names(jt_df))
  
  term_col <- if ("model.term" %in% names(jt_df)) {
    "model.term"
  } else if ("term" %in% names(jt_df)) {
    "term"
  } else {
    NULL
  }
  
  stat_col <- if ("F.ratio" %in% names(jt_df)) {
    "F.ratio"
  } else if ("Chisq" %in% names(jt_df)) {
    "Chisq"
  } else if ("t.ratio" %in% names(jt_df)) {
    "t.ratio"
  } else {
    NULL
  }
  
  p_col <- if ("p.value" %in% names(jt_df)) {
    "p.value"
  } else if ("P.value" %in% names(jt_df)) {
    "P.value"
  } else {
    NULL
  }
  
  if (is.null(term_col) || is.null(stat_col) || is.null(p_col)) {
    return(empty_tests())
  }
  
  out <- jt_df %>%
    dplyr::mutate(
      term = .data[[term_col]],
      statistic = .data[[stat_col]],
      p_value = .data[[p_col]],
      df1 = if ("df1" %in% names(.)) df1 else NA_real_,
      df2 = if ("df2" %in% names(.)) df2 else NA_real_
    ) %>%
    dplyr::filter(
      !is.na(term),
      term != "(confounded)",
      term != "1",
      term %in% fit_result$active_terms
    ) %>%
    dplyr::transmute(
      comparison_block = comparison_block,
      response_family = response_family,
      response = response,
      engine = fit_result$engine,
      formula = fit_result$formula,
      term = as.character(term),
      df1 = as.numeric(df1),
      df2 = as.numeric(df2),
      statistic = as.numeric(statistic),
      p_value = as.numeric(p_value),
      p_value_adjusted = NA_real_,
      p_adjust_method = NA_character_
    )
  
  if (nrow(out) == 0) {
    return(empty_tests())
  }
  
  out
}

extract_emmeans_one_term <- function(fit_result, comparison_block, response_family, response, term) {
  fit <- fit_result$fit
  
  emm_obj <- tryCatch(
    emmeans::emmeans(fit, specs = stats::as.formula(paste("~", term))),
    error = function(e) NULL
  )
  
  if (is.null(emm_obj)) {
    return(list(
      emmeans = empty_emmeans(),
      contrasts = empty_contrasts()
    ))
  }
  
  emm_df <- as.data.frame(emm_obj)
  
  if (nrow(emm_df) == 0 || !term %in% names(emm_df)) {
    emm_out <- empty_emmeans()
  } else {
    emm_out <- emm_df %>%
      dplyr::transmute(
        comparison_block = comparison_block,
        response_family = response_family,
        response = response,
        engine = fit_result$engine,
        formula = fit_result$formula,
        term = term,
        level = as.character(.data[[term]]),
        emmean = if ("emmean" %in% names(.)) as.numeric(emmean) else NA_real_,
        se = if ("SE" %in% names(.)) as.numeric(SE) else NA_real_,
        df = if ("df" %in% names(.)) as.numeric(df) else NA_real_,
        lower_cl = if ("lower.CL" %in% names(.)) as.numeric(lower.CL) else NA_real_,
        upper_cl = if ("upper.CL" %in% names(.)) as.numeric(upper.CL) else NA_real_
      )
  }
  
  contrast_obj <- tryCatch(
    emmeans::pairs(emm_obj, adjust = "none"),
    error = function(e) NULL
  )
  
  if (is.null(contrast_obj)) {
    contrast_out <- empty_contrasts()
  } else {
    contrast_df <- tryCatch(
      as.data.frame(summary(contrast_obj, infer = c(TRUE, TRUE))),
      error = function(e) NULL
    )
    
    if (is.null(contrast_df) || nrow(contrast_df) == 0) {
      contrast_out <- empty_contrasts()
    } else {
      names(contrast_df) <- make.names(names(contrast_df))
      
      stat_col <- if ("t.ratio" %in% names(contrast_df)) {
        "t.ratio"
      } else if ("z.ratio" %in% names(contrast_df)) {
        "z.ratio"
      } else {
        NULL
      }
      
      contrast_out <- contrast_df %>%
        dplyr::transmute(
          comparison_block = comparison_block,
          response_family = response_family,
          response = response,
          engine = fit_result$engine,
          formula = fit_result$formula,
          term = term,
          contrast = as.character(contrast),
          estimate = if ("estimate" %in% names(.)) as.numeric(estimate) else NA_real_,
          se = if ("SE" %in% names(.)) as.numeric(SE) else NA_real_,
          df = if ("df" %in% names(.)) as.numeric(df) else NA_real_,
          statistic = if (!is.null(stat_col)) as.numeric(.data[[stat_col]]) else NA_real_,
          lower_cl = if ("lower.CL" %in% names(.)) as.numeric(lower.CL) else NA_real_,
          upper_cl = if ("upper.CL" %in% names(.)) as.numeric(upper.CL) else NA_real_,
          p_value = if ("p.value" %in% names(.)) as.numeric(p.value) else NA_real_,
          p_value_adjusted = NA_real_,
          p_adjust_method = NA_character_
        )
    }
  }
  
  list(
    emmeans = emm_out,
    contrasts = contrast_out
  )
}

adjust_pvalues <- function(df, group_cols, method) {
  if (nrow(df) == 0 || !"p_value" %in% names(df)) {
    return(df)
  }
  
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      p_value_adjusted = stats::p.adjust(p_value, method = method),
      p_adjust_method = method
    ) %>%
    dplyr::ungroup()
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

save_primary_figures <- function(emmeans_df, contrasts_df, output_dir) {
  if (nrow(emmeans_df) > 0) {
    p_emm <- ggplot2::ggplot(
      emmeans_df,
      ggplot2::aes(x = level, y = emmean, ymin = lower_cl, ymax = upper_cl)
    ) +
      ggplot2::geom_pointrange() +
      ggplot2::facet_grid(response ~ comparison_block + term, scales = "free_y") +
      ggplot2::labs(
        title = "Estimated marginal means for primary spectral metrics",
        x = "Level",
        y = "Estimated marginal mean"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    ggplot2::ggsave(
      filename = file.path(output_dir, "06_primary_emmeans.png"),
      plot = p_emm,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
  
  if (nrow(contrasts_df) > 0) {
    p_contr <- ggplot2::ggplot(
      contrasts_df,
      ggplot2::aes(x = contrast, y = estimate, ymin = lower_cl, ymax = upper_cl)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_pointrange() +
      ggplot2::facet_grid(response ~ comparison_block + term, scales = "free_y") +
      ggplot2::labs(
        title = "Pairwise contrasts for primary spectral metrics",
        x = "Contrast",
        y = "Estimated difference"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    ggplot2::ggsave(
      filename = file.path(output_dir, "06_primary_contrasts.png"),
      plot = p_contr,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
}

run_models_for_family <- function(
    df,
    comparison_block,
    response_family,
    responses,
    candidate_terms,
    require_individual_random = FALSE
) {
  model_overview <- empty_model_overview()
  diagnostics <- empty_diagnostics()
  tests <- empty_tests()
  emmeans_out <- empty_emmeans()
  contrasts_out <- empty_contrasts()
  
  if (nrow(df) == 0 || length(responses) == 0) {
    return(list(
      model_overview = model_overview,
      diagnostics = diagnostics,
      tests = tests,
      emmeans = emmeans_out,
      contrasts = contrasts_out
    ))
  }
  
  for (response in responses) {
    fit_result <- fit_best_model(
      data = df,
      response = response,
      candidate_terms = candidate_terms,
      require_individual_random = require_individual_random
    )
    
    model_overview <- dplyr::bind_rows(
      model_overview,
      extract_model_overview_row(
        fit_result = fit_result,
        comparison_block = comparison_block,
        response_family = response_family,
        response = response
      )
    )
    
    if (!identical(fit_result$status, "fitted") || is.null(fit_result$fit)) {
      next
    }
    
    diagnostics <- dplyr::bind_rows(
      diagnostics,
      extract_diagnostics_row(
        fit_result = fit_result,
        comparison_block = comparison_block,
        response_family = response_family,
        response = response
      )
    )
    
    tests <- dplyr::bind_rows(
      tests,
      extract_joint_tests(
        fit_result = fit_result,
        comparison_block = comparison_block,
        response_family = response_family,
        response = response
      )
    )
    
    for (term in fit_result$active_terms) {
      term_out <- extract_emmeans_one_term(
        fit_result = fit_result,
        comparison_block = comparison_block,
        response_family = response_family,
        response = response,
        term = term
      )
      
      emmeans_out <- dplyr::bind_rows(emmeans_out, term_out$emmeans)
      contrasts_out <- dplyr::bind_rows(contrasts_out, term_out$contrasts)
    }
  }
  
  list(
    model_overview = model_overview,
    diagnostics = diagnostics,
    tests = tests,
    emmeans = emmeans_out,
    contrasts = contrasts_out
  )
}