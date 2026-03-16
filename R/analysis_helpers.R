spectral_index_names <- function() {
  c(
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
}

safe_divide <- function(num, den) {
  if (length(num) == 0 || length(den) == 0) {
    return(NA_real_)
  }
  
  if (is.na(num) || is.na(den) || den == 0) {
    return(NA_real_)
  }
  
  num / den
}

nearest_reflectance <- function(wavelength, reflectance, target_wavelength) {
  ok <- !is.na(wavelength) & !is.na(reflectance)
  
  if (sum(ok) == 0) {
    return(NA_real_)
  }
  
  wavelength <- wavelength[ok]
  reflectance <- reflectance[ok]
  
  idx <- which.min(abs(wavelength - target_wavelength))
  
  if (length(idx) == 0) {
    return(NA_real_)
  }
  
  reflectance[idx[1]]
}

band_mean <- function(wavelength, reflectance, wl_min, wl_max) {
  ok <- !is.na(wavelength) & !is.na(reflectance) &
    wavelength >= wl_min & wavelength <= wl_max
  
  if (sum(ok) == 0) {
    return(NA_real_)
  }
  
  mean(reflectance[ok], na.rm = TRUE)
}

compute_first_derivative <- function(wavelength, reflectance) {
  ok <- !is.na(wavelength) & !is.na(reflectance)
  wavelength <- wavelength[ok]
  reflectance <- reflectance[ok]
  
  if (length(wavelength) < 2) {
    return(tibble::tibble(fd_wavelength = numeric(), fd = numeric()))
  }
  
  ord <- order(wavelength)
  wavelength <- wavelength[ord]
  reflectance <- reflectance[ord]
  
  tibble::tibble(
    fd_wavelength = (wavelength[-1] + wavelength[-length(wavelength)]) / 2,
    fd = diff(reflectance) / diff(wavelength)
  )
}

compute_spectral_indices <- function(spectrum_df) {
  stopifnot(all(c("wavelength", "reflectance") %in% names(spectrum_df)))
  
  dat <- spectrum_df %>%
    dplyr::transmute(
      wavelength = as.numeric(wavelength),
      reflectance = as.numeric(reflectance)
    ) %>%
    dplyr::filter(!is.na(wavelength), !is.na(reflectance)) %>%
    dplyr::arrange(wavelength)
  
  if (nrow(dat) == 0) {
    return(tibble::tibble(
      mfdre = NA_real_,
      rep = NA_real_,
      res700_740 = NA_real_,
      mean_blue_450_500 = NA_real_,
      mean_green_500_570 = NA_real_,
      mean_red_650_680 = NA_real_,
      mean760_900 = NA_real_,
      datt = NA_real_,
      pri = NA_real_,
      sipi = NA_real_,
      psri = NA_real_
    ))
  }
  
  r445 <- nearest_reflectance(dat$wavelength, dat$reflectance, 445)
  r500 <- nearest_reflectance(dat$wavelength, dat$reflectance, 500)
  r531 <- nearest_reflectance(dat$wavelength, dat$reflectance, 531)
  r570 <- nearest_reflectance(dat$wavelength, dat$reflectance, 570)
  r680 <- nearest_reflectance(dat$wavelength, dat$reflectance, 680)
  r700 <- nearest_reflectance(dat$wavelength, dat$reflectance, 700)
  r710 <- nearest_reflectance(dat$wavelength, dat$reflectance, 710)
  r740 <- nearest_reflectance(dat$wavelength, dat$reflectance, 740)
  r750 <- nearest_reflectance(dat$wavelength, dat$reflectance, 750)
  r800 <- nearest_reflectance(dat$wavelength, dat$reflectance, 800)
  r850 <- nearest_reflectance(dat$wavelength, dat$reflectance, 850)
  
  fd_df <- compute_first_derivative(dat$wavelength, dat$reflectance) %>%
    dplyr::filter(fd_wavelength >= 680, fd_wavelength <= 750)
  
  if (nrow(fd_df) > 0 && any(!is.na(fd_df$fd))) {
    max_idx <- which.max(fd_df$fd)
    mfdre <- fd_df$fd[max_idx]
    rep <- fd_df$fd_wavelength[max_idx]
  } else {
    mfdre <- NA_real_
    rep <- NA_real_
  }
  
  tibble::tibble(
    mfdre = mfdre,
    rep = rep,
    res700_740 = safe_divide(r740 - r700, 40),
    mean_blue_450_500 = band_mean(dat$wavelength, dat$reflectance, 450, 500),
    mean_green_500_570 = band_mean(dat$wavelength, dat$reflectance, 500, 570),
    mean_red_650_680 = band_mean(dat$wavelength, dat$reflectance, 650, 680),
    mean760_900 = band_mean(dat$wavelength, dat$reflectance, 760, 900),
    datt = safe_divide(r850 - r710, r850 - r680),
    pri = safe_divide(r531 - r570, r531 + r570),
    sipi = safe_divide(r800 - r445, r800 - r680),
    psri = safe_divide(r680 - r500, r750)
  )
}

spectral_rmse <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  
  if (sum(ok) == 0) {
    return(NA_real_)
  }
  
  sqrt(mean((x[ok] - y[ok])^2))
}

spectral_sam <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  
  if (sum(ok) == 0) {
    return(NA_real_)
  }
  
  x <- x[ok]
  y <- y[ok]
  
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  
  if (denom == 0) {
    return(NA_real_)
  }
  
  cosine_val <- sum(x * y) / denom
  cosine_val <- max(min(cosine_val, 1), -1)
  
  acos(cosine_val)
}

spectral_iauc <- function(wavelength, delta_reflectance) {
  ok <- !is.na(wavelength) & !is.na(delta_reflectance)
  
  if (sum(ok) < 2) {
    return(NA_real_)
  }
  
  wavelength <- wavelength[ok]
  delta_reflectance <- abs(delta_reflectance[ok])
  
  ord <- order(wavelength)
  wavelength <- wavelength[ord]
  delta_reflectance <- delta_reflectance[ord]
  
  sum((delta_reflectance[-1] + delta_reflectance[-length(delta_reflectance)]) / 2 * diff(wavelength))
}

make_index_delta_table <- function(df, index_names, sample_suffix = "_sample", ref_suffix = "_ref") {
  out <- df
  
  for (nm in index_names) {
    sample_col <- paste0(nm, sample_suffix)
    ref_col <- paste0(nm, ref_suffix)
    delta_col <- paste0("delta_", nm)
    
    if (all(c(sample_col, ref_col) %in% names(out))) {
      out[[delta_col]] <- out[[sample_col]] - out[[ref_col]]
    } else {
      out[[delta_col]] <- NA_real_
    }
  }
  
  out
}

summarise_numeric_by_group <- function(df, group_cols, value_cols) {
  value_cols <- intersect(value_cols, names(df))
  
  if (length(value_cols) == 0) {
    return(tibble::tibble())
  }
  
  df %>%
    dplyr::select(dplyr::any_of(c(group_cols, value_cols))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(value_cols),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), metric) %>%
    dplyr::summarise(
      n = sum(!is.na(value)),
      mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
      sd = if (sum(!is.na(value)) >= 2) stats::sd(value, na.rm = TRUE) else NA_real_,
      cv = if (!is.na(mean) && !is.na(sd) && mean != 0) sd / abs(mean) else NA_real_,
      median = if (all(is.na(value))) NA_real_ else stats::median(value, na.rm = TRUE),
      iqr = if (sum(!is.na(value)) >= 1) stats::IQR(value, na.rm = TRUE) else NA_real_,
      min = if (all(is.na(value))) NA_real_ else min(value, na.rm = TRUE),
      max = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
      .groups = "drop"
    )
}

empty_comparison_metrics <- function(index_names = spectral_index_names()) {
  base <- tibble::tibble(
    comparison_type = character(),
    comparison_id = character(),
    line = character(),
    individual = character(),
    drying = character(),
    ageing = character(),
    sample_group_sample = character(),
    sample_group_ref = character(),
    timepoint_sample = character(),
    timepoint_ref = character(),
    n_replicates_sample = integer(),
    n_replicates_ref = integer(),
    rmse = double(),
    sam = double(),
    iauc = double()
  )
  
  delta_cols <- stats::setNames(
    rep(list(numeric()), length(index_names)),
    paste0("delta_", index_names)
  )
  
  tibble::add_column(base, !!!delta_cols)
}

empty_delta_spectra <- function() {
  tibble::tibble(
    comparison_type = character(),
    comparison_id = character(),
    line = character(),
    individual = character(),
    drying = character(),
    ageing = character(),
    wavelength = double(),
    reflectance_ref = double(),
    reflectance_sample = double(),
    delta_reflectance = double()
  )
}

save_plot_if_has_data <- function(plot_obj, filename, width = 9, height = 6, dpi = 300) {
  if (inherits(plot_obj, "ggplot")) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi
    )
  }
}