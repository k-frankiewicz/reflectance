build_sample_group <- function(timepoint, line, individual, drying, ageing) {
  case_when(
    timepoint == "fresh" ~ paste(line, individual, sep = "__"),
    timepoint == "dried" ~ paste(line, individual, drying, sep = "__"),
    timepoint == "aged" ~ paste(line, individual, drying, ageing, sep = "__"),
    TRUE ~ NA_character_
  )
}

read_one_processed_spectrum <- function(
    file,
    analysis_min_wl = 400,
    analysis_max_wl = 950
) {
  df <- read_csv(file, show_col_types = FALSE, progress = FALSE)
  
  required_cols <- c(
    "section",
    "wavelength",
    "calibrated_and_averaged_data"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in file ",
      basename(file),
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    filter(section == "spectrum") %>%
    transmute(
      file = basename(file),
      wavelength = as.numeric(wavelength),
      reflectance = as.numeric(calibrated_and_averaged_data)
    ) %>%
    filter(!is.na(wavelength), !is.na(reflectance)) %>%
    filter(wavelength >= analysis_min_wl, wavelength <= analysis_max_wl)
}

safe_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  
  cor(x[ok], y[ok], method = "pearson")
}

rmse_vec <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  
  if (sum(ok) == 0) {
    return(NA_real_)
  }
  
  sqrt(mean((x[ok] - y[ok])^2))
}

compute_outlier_threshold <- function(x, direction = c("high", "low")) {
  direction <- match.arg(direction)
  x <- x[!is.na(x)]
  
  if (length(x) < 4) {
    return(NA_real_)
  }
  
  q1 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
  q3 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
  iqr_val <- q3 - q1
  
  if (direction == "high") {
    q3 + 3 * iqr_val
  } else {
    q1 - 3 * iqr_val
  }
}