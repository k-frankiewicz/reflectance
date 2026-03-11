suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

raw_files <- list.files(raw_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(raw_files) == 0) {
  stop(
    "No .csv files found: ",
    normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  )
}

files <- raw_files[
  vapply(raw_files, function(file) {
    out_file <- file.path(processed_dir, basename(file))
    
    if (!file.exists(out_file)) {
      return(TRUE)
    }
    
    file.info(file)$mtime > file.info(out_file)$mtime
  }, logical(1))
]

if (length(files) == 0) {
  message("No new/updated files to be processed.")
}

empty_output <- function() {
  tibble(
    source_file = character(),
    specimen_id = character(),
    layer_title = character(),
    color = character(),
    mode = character(),
    integration_time = character(),
    boxcar_width = character(),
    scans_to_average = character(),
    section = character(),
    record_id = integer(),
    wavelength = double(),
    raw_spectrometer_data = double(),
    calibrated_and_averaged_data = double(),
    peak_wavelength = double(),
    light_calibration_value = double(),
    darkness_calibration_value = double()
  )
}

read_lines_clean <- function(file) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  if (length(x) > 0) {
    x[1] <- sub("^\ufeff", "", x[1])
  }
  x
}

parse_simple_row <- function(line, expected_n) {
  vals <- trimws(strsplit(line, ",", fixed = TRUE)[[1]])
  if (length(vals) < expected_n) {
    vals <- c(vals, rep(NA_character_, expected_n - length(vals)))
  }
  vals[seq_len(expected_n)]
}

get_block_lines <- function(lines, header_pattern, stop_patterns = character()) {
  header_idx <- which(str_detect(lines, header_pattern))[1]
  
  if (is.na(header_idx) || header_idx >= length(lines)) {
    return(character(0))
  }
  
  after <- lines[(header_idx + 1):length(lines)]
  stop_idx <- which(trimws(after) == "")
  
  if (length(stop_patterns) > 0) {
    for (pat in stop_patterns) {
      stop_idx <- c(stop_idx, which(str_detect(after, pat)))
    }
  }
  
  stop_idx <- sort(unique(stop_idx))
  
  if (length(stop_idx) == 0) {
    block <- after
  } else {
    block <- after[seq_len(min(stop_idx) - 1)]
  }
  
  block[trimws(block) != ""]
}

split_tokens <- function(line) {
  trimws(strsplit(line, ",", fixed = TRUE)[[1]])
}

join_decimal <- function(a, b) {
  paste0(a, ".", b)
}

to_num <- function(x) {
  as.numeric(str_remove(x, "%$"))
}

extract_metadata <- function(lines, file) {
  specimen_fallback <- tools::file_path_sans_ext(basename(file))
  
  layer_header_idx <- which(
    str_detect(lines, "^Layer Title\\s*,\\s*Color\\s*,\\s*Specimen ID\\s*$")
  )[1]
  
  mode_header_idx <- which(
    str_detect(lines, "^Mode\\s*,\\s*Integration Time\\s*,\\s*Boxcar Width\\s*,\\s*Scans to Average\\s*$")
  )[1]
  
  layer_vals <- c(NA_character_, NA_character_, specimen_fallback)
  mode_vals  <- c(NA_character_, NA_character_, NA_character_, NA_character_)
  
  if (!is.na(layer_header_idx) && layer_header_idx < length(lines)) {
    layer_vals <- parse_simple_row(lines[layer_header_idx + 1], 3)
  }
  
  if (!is.na(mode_header_idx) && mode_header_idx < length(lines)) {
    mode_vals <- parse_simple_row(lines[mode_header_idx + 1], 4)
  }
  
  list(
    layer_title = layer_vals[1],
    color = layer_vals[2],
    specimen_id = ifelse(
      is.na(layer_vals[3]) || layer_vals[3] == "",
      specimen_fallback,
      layer_vals[3]
    ),
    mode = mode_vals[1],
    integration_time = mode_vals[2],
    boxcar_width = mode_vals[3],
    scans_to_average = mode_vals[4]
  )
}

parse_spectrum_line <- function(line) {
  tok <- split_tokens(line)
  
  if (length(tok) == 3) {
    return(tibble(
      wavelength = to_num(tok[1]),
      raw_spectrometer_data = to_num(tok[2]),
      calibrated_and_averaged_data = to_num(tok[3])
    ))
  }
  
  if (length(tok) == 5) {
    return(tibble(
      wavelength = to_num(join_decimal(tok[1], tok[2])),
      raw_spectrometer_data = to_num(join_decimal(tok[3], tok[4])),
      calibrated_and_averaged_data = to_num(tok[5])
    ))
  }
  
  if (length(tok) == 6) {
    return(tibble(
      wavelength = to_num(join_decimal(tok[1], tok[2])),
      raw_spectrometer_data = to_num(join_decimal(tok[3], tok[4])),
      calibrated_and_averaged_data = to_num(join_decimal(tok[5], tok[6]))
    ))
  }
  
  stop("Unexpected number of tokens in spectrum section: ", length(tok), " | ", line)
}

parse_peak_line <- function(line) {
  tok <- split_tokens(line)
  
  if (length(tok) == 1) {
    return(tibble(peak_wavelength = to_num(tok[1])))
  }
  
  if (length(tok) == 2) {
    return(tibble(peak_wavelength = to_num(join_decimal(tok[1], tok[2]))))
  }
  
  stop("Unexpected number of tokens in peaks section: ", length(tok), " | ", line)
}

parse_calibration_line <- function(line) {
  tok <- split_tokens(line)
  
  if (length(tok) == 3) {
    return(tibble(
      wavelength = to_num(tok[1]),
      light_calibration_value = to_num(tok[2]),
      darkness_calibration_value = to_num(tok[3])
    ))
  }
  
  if (length(tok) == 6) {
    return(tibble(
      wavelength = to_num(join_decimal(tok[1], tok[2])),
      light_calibration_value = to_num(join_decimal(tok[3], tok[4])),
      darkness_calibration_value = to_num(join_decimal(tok[5], tok[6]))
    ))
  }
  
  stop("Unexpected number of tokens in calibration section: ", length(tok), " | ", line)
}

read_ci710_file <- function(file) {
  lines <- read_lines_clean(file)
  meta <- extract_metadata(lines, file)
  
  spectrum_lines <- get_block_lines(
    lines,
    "^Wavelength\\s*,\\s*Raw Spectrometer Data\\s*,\\s*Calibrated and Averaged Data\\s*$",
    stop_patterns = c(
      "^Peak Wavelengths\\s*$",
      "^Wavelength\\s*,\\s*Light Calibration Value\\s*,\\s*Darkness Calibration Value\\s*$"
    )
  )
  
  peaks_lines <- get_block_lines(
    lines,
    "^Peak Wavelengths\\s*$",
    stop_patterns = c(
      "^Wavelength\\s*,\\s*Light Calibration Value\\s*,\\s*Darkness Calibration Value\\s*$"
    )
  )
  
  calibration_lines <- get_block_lines(
    lines,
    "^Wavelength\\s*,\\s*Light Calibration Value\\s*,\\s*Darkness Calibration Value\\s*$"
  )
  
  out <- empty_output()
  
  if (length(spectrum_lines) > 0) {
    spectrum_df <- map_dfr(seq_along(spectrum_lines), function(i) {
      vals <- parse_spectrum_line(spectrum_lines[i])
      
      tibble(
        source_file = basename(file),
        specimen_id = meta$specimen_id,
        layer_title = meta$layer_title,
        color = meta$color,
        mode = meta$mode,
        integration_time = meta$integration_time,
        boxcar_width = meta$boxcar_width,
        scans_to_average = meta$scans_to_average,
        section = "spectrum",
        record_id = i,
        wavelength = vals$wavelength,
        raw_spectrometer_data = vals$raw_spectrometer_data,
        calibrated_and_averaged_data = vals$calibrated_and_averaged_data,
        peak_wavelength = NA_real_,
        light_calibration_value = NA_real_,
        darkness_calibration_value = NA_real_
      )
    })
    
    out <- bind_rows(out, spectrum_df)
  }
  
  if (length(peaks_lines) > 0) {
    peaks_df <- map_dfr(seq_along(peaks_lines), function(i) {
      vals <- parse_peak_line(peaks_lines[i])
      
      tibble(
        source_file = basename(file),
        specimen_id = meta$specimen_id,
        layer_title = meta$layer_title,
        color = meta$color,
        mode = meta$mode,
        integration_time = meta$integration_time,
        boxcar_width = meta$boxcar_width,
        scans_to_average = meta$scans_to_average,
        section = "peaks",
        record_id = i,
        wavelength = NA_real_,
        raw_spectrometer_data = NA_real_,
        calibrated_and_averaged_data = NA_real_,
        peak_wavelength = vals$peak_wavelength,
        light_calibration_value = NA_real_,
        darkness_calibration_value = NA_real_
      )
    })
    
    out <- bind_rows(out, peaks_df)
  }
  
  if (length(calibration_lines) > 0) {
    calibration_df <- map_dfr(seq_along(calibration_lines), function(i) {
      vals <- parse_calibration_line(calibration_lines[i])
      
      tibble(
        source_file = basename(file),
        specimen_id = meta$specimen_id,
        layer_title = meta$layer_title,
        color = meta$color,
        mode = meta$mode,
        integration_time = meta$integration_time,
        boxcar_width = meta$boxcar_width,
        scans_to_average = meta$scans_to_average,
        section = "calibration",
        record_id = i,
        wavelength = vals$wavelength,
        raw_spectrometer_data = NA_real_,
        calibrated_and_averaged_data = NA_real_,
        peak_wavelength = NA_real_,
        light_calibration_value = vals$light_calibration_value,
        darkness_calibration_value = vals$darkness_calibration_value
      )
    })
    
    out <- bind_rows(out, calibration_df)
  }
  
  if (nrow(out) == 0) {
    stop("No sections found in file: ", basename(file))
  }
  
  out
}

walk(files, function(file) {
  message("Processing: ", basename(file))
  
  one_file <- read_ci710_file(file)
  
  out_file <- file.path(processed_dir, basename(file))
  write_csv(one_file, out_file, na = "")
  
  message(
    "Saved: ",
    normalizePath(out_file, winslash = "/", mustWork = FALSE)
  )
  
  print(one_file %>% count(section))
  print(one_file %>% filter(section == "spectrum") %>% slice_head(n = 15))
  print(one_file %>% filter(section == "peaks"))
  print(one_file %>% filter(section == "calibration") %>% slice_head(n = 15))
})