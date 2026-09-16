collapse_unique <- function(x) {
  x <- unique(na.omit(as.character(x)))
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = " | ")
}

add_reason <- function(existing, new_reason) {
  if (is.na(existing) || existing == "") {
    return(new_reason)
  }
  
  paste(existing, new_reason, sep = "; ")
}

combine_two_reasons <- function(a, b) {
  x <- c(a, b)
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = "; ")
}

check_expected_value <- function(x, expected) {
  vals <- unique(x[!is.na(x)])
  
  if (length(vals) == 0) {
    return(FALSE)
  }
  
  length(vals) == 1 && vals[[1]] %in% expected
}

empty_qc_output <- function(file_name) {
  tibble(
    file = file_name,
    n_rows = NA_integer_,
    n_columns = NA_integer_,
    has_spectrum = NA,
    has_calibration = NA,
    n_spectrum_rows = NA_integer_,
    n_calibration_rows = NA_integer_,
    n_points_used = NA_integer_,
    wl_min_used = NA_real_,
    wl_max_used = NA_real_,
    min_reflectance = NA_real_,
    max_reflectance = NA_real_,
    pct_clipped = NA_real_,
    min_cal_span = NA_real_,
    pct_cal_span_lt_2000 = NA_real_,
    roughness = NA_real_,
    mode_found = NA_character_,
    integration_time_found = NA_character_,
    boxcar_width_found = NA_character_,
    scans_to_average_found = NA_character_,
    fail_reasons = NA_character_,
    warn_reasons = NA_character_,
    qc_status = "pass"
  )
}

assess_file <- function(
    path,
    analysis_min_wl = 400,
    analysis_max_wl = 950,
    expected_mode = "Reflectance",
    expected_integration_time = 400,
    expected_boxcar_width = 2,
    expected_scans_to_average = 1,
    clip_threshold_low = 0.1,
    clip_threshold_high = 99.9,
    warn_pct_clipped = 0.5,
    low_cal_threshold = 2000,
    warn_pct_low_cal = 1
) {
  file_name <- basename(path)
  out <- empty_qc_output(file_name)
  
  df <- tryCatch(
    read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(df)) {
    out$fail_reasons <- "read_error"
    out$qc_status <- "fail"
    return(out)
  }
  
  out$n_rows <- nrow(df)
  out$n_columns <- ncol(df)
  
  required_cols <- c(
    "section", "wavelength", "mode",
    "integration_time", "boxcar_width", "scans_to_average",
    "raw_spectrometer_data", "calibrated_and_averaged_data",
    "light_calibration_value", "darkness_calibration_value"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    out$fail_reasons <- paste0(
      "missing_columns: ",
      paste(missing_cols, collapse = ", ")
    )
    out$qc_status <- "fail"
    return(out)
  }
  
  out$has_spectrum <- any(df$section == "spectrum", na.rm = TRUE)
  out$has_calibration <- any(df$section == "calibration", na.rm = TRUE)
  out$n_spectrum_rows <- sum(df$section == "spectrum", na.rm = TRUE)
  out$n_calibration_rows <- sum(df$section == "calibration", na.rm = TRUE)
  
  out$mode_found <- collapse_unique(df$mode)
  out$integration_time_found <- collapse_unique(df$integration_time)
  out$boxcar_width_found <- collapse_unique(df$boxcar_width)
  out$scans_to_average_found <- collapse_unique(df$scans_to_average)
  
  if (!out$has_spectrum) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "missing_spectrum_section"
    )
  }
  
  if (!out$has_calibration) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "missing_calibration_section"
    )
  }
  
  if (!check_expected_value(as.character(df$mode), expected_mode)) {
    out$fail_reasons <- add_reason(out$fail_reasons, "unexpected_mode")
  }
  
  if (!check_expected_value(
    parse_number(as.character(df$integration_time)),
    expected_integration_time
  )) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "unexpected_integration_time"
    )
  }
  
  if (!check_expected_value(
    parse_number(as.character(df$boxcar_width)),
    expected_boxcar_width
  )) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "unexpected_boxcar_width"
    )
  }
  
  if (!check_expected_value(
    parse_number(as.character(df$scans_to_average)),
    expected_scans_to_average
  )) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "unexpected_scans_to_average"
    )
  }
  
  spec <- df %>%
    filter(section == "spectrum") %>%
    transmute(
      wavelength = as.numeric(wavelength),
      reflectance = as.numeric(calibrated_and_averaged_data),
      raw_signal = as.numeric(raw_spectrometer_data)
    )
  
  cal <- df %>%
    filter(section == "calibration") %>%
    transmute(
      wavelength = as.numeric(wavelength),
      light_calibration_value = as.numeric(light_calibration_value),
      darkness_calibration_value = as.numeric(darkness_calibration_value),
      cal_span = light_calibration_value - darkness_calibration_value
    )
  
  spec_use <- spec %>%
    filter(!is.na(wavelength)) %>%
    filter(wavelength >= analysis_min_wl, wavelength <= analysis_max_wl)
  
  cal_use <- cal %>%
    filter(!is.na(wavelength)) %>%
    filter(wavelength >= analysis_min_wl, wavelength <= analysis_max_wl)
  
  out$n_points_used <- nrow(spec_use)
  
  if (nrow(spec_use) == 0) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "no_spectrum_points_in_range"
    )
  } else {
    out$wl_min_used <- min(spec_use$wavelength, na.rm = TRUE)
    out$wl_max_used <- max(spec_use$wavelength, na.rm = TRUE)
    
    if (any(is.na(spec_use$reflectance))) {
      out$fail_reasons <- add_reason(
        out$fail_reasons,
        "missing_reflectance_values"
      )
    }
    
    if (anyDuplicated(spec_use$wavelength) > 0) {
      out$fail_reasons <- add_reason(
        out$fail_reasons,
        "duplicate_wavelengths"
      )
    }
    
    if (any(diff(spec_use$wavelength) <= 0, na.rm = TRUE)) {
      out$fail_reasons <- add_reason(
        out$fail_reasons,
        "non_monotonic_wavelengths"
      )
    }
    
    if (!all(is.na(spec_use$reflectance))) {
      out$min_reflectance <- min(spec_use$reflectance, na.rm = TRUE)
      out$max_reflectance <- max(spec_use$reflectance, na.rm = TRUE)
      
      clipped <- spec_use$reflectance <= clip_threshold_low |
        spec_use$reflectance >= clip_threshold_high
      
      out$pct_clipped <- mean(clipped, na.rm = TRUE) * 100
      
      if (!is.na(out$pct_clipped) && out$pct_clipped > warn_pct_clipped) {
        out$warn_reasons <- add_reason(
          out$warn_reasons,
          "high_clipping"
        )
      }
      
      if (nrow(spec_use) >= 3) {
        out$roughness <- median(abs(diff(spec_use$reflectance)), na.rm = TRUE)
      }
    }
  }
  
  if (nrow(cal_use) == 0) {
    out$fail_reasons <- add_reason(
      out$fail_reasons,
      "no_calibration_points_in_range"
    )
  } else {
    if (any(is.na(cal_use$cal_span))) {
      out$fail_reasons <- add_reason(
        out$fail_reasons,
        "missing_calibration_values"
      )
    } else {
      out$min_cal_span <- min(cal_use$cal_span, na.rm = TRUE)
      out$pct_cal_span_lt_2000 <- mean(
        cal_use$cal_span < low_cal_threshold,
        na.rm = TRUE
      ) * 100
      
      if (!is.na(out$min_cal_span) && out$min_cal_span <= 0) {
        out$fail_reasons <- add_reason(
          out$fail_reasons,
          "non_positive_cal_span"
        )
      }
      
      if (!is.na(out$pct_cal_span_lt_2000) &&
          out$pct_cal_span_lt_2000 > warn_pct_low_cal) {
        out$warn_reasons <- add_reason(
          out$warn_reasons,
          "low_cal_span"
        )
      }
    }
  }
  
  if (!is.na(out$fail_reasons) && out$fail_reasons != "") {
    out$qc_status <- "fail"
  } else if (!is.na(out$warn_reasons) && out$warn_reasons != "") {
    out$qc_status <- "warn"
  } else {
    out$qc_status <- "pass"
  }
  
  out
}