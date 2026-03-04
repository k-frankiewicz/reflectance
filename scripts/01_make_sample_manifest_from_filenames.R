library(tidyverse)

# Paths
raw_dir  <- "data/raw"
out_path <- "data/metadata/sample_manifest.csv"

# List all CSV files in data/raw
files <- list.files(raw_dir, pattern = "\\.csv$", full.names = TRUE)

# Parse metadata from filename
parse_file <- function(path) {
  fname <- basename(path)
  stem  <- sub("\\.csv$", "", fname)
  parts <- strsplit(stem, ".", fixed = TRUE)[[1]]
  
  line <- parts[1]
  individual <- parts[2]
  
  if (length(parts) == 3) {
    # LINE.INDIVIDUAL.REP 
    drying <- NA_character_
    ageing <- NA_character_
    replicate <- as.integer(parts[3])
    timepoint <- "fresh"
    
  } else if (length(parts) == 4) {
    # LINE.INDIVIDUAL.DRYING.REP
    drying <- parts[3]
    ageing <- NA_character_
    replicate <- as.integer(parts[4])
    timepoint <- "dried"
    
  } else if (length(parts) == 5) {
    # LINE.INDIVIDUAL.DRYING.AGEING.REP
    drying <- parts[3]
    ageing <- parts[4]
    replicate <- as.integer(parts[5])
    timepoint <- "aged"
    
  } else {
    stop("Unexpected filename format: ", fname)
  }
  
  # Validate codes
  if (!is.na(drying) && !drying %in% c("P", "C", "L")) {
    stop("Unknown drying code (allowed: P/C/L) in file: ", fname)
  }
  if (!is.na(ageing) && !ageing %in% c("T", "H", "B")) {
    stop("Unknown ageing code (allowed: T/H/B) in file: ", fname)
  }
  if (is.na(replicate)) {
    stop("Invalid replicate (should be an integer) in file: ", fname)
  }
  
  tibble(
    file = fname,
    line = line,
    individual = individual,
    replicate = replicate,
    timepoint = timepoint,
    drying = drying,
    ageing = ageing,
    notes = NA_character_
  )
}

# Build a fresh manifest from current files
manifest_new <- purrr::map_dfr(files, parse_file) %>%
  mutate(timepoint = factor(timepoint, levels = c("fresh", "dried", "aged"), ordered = TRUE)) %>%
  arrange(timepoint, line, individual, replicate, file) %>%
  mutate(timepoint = as.character(timepoint))

# Ensure output directory exists
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

# Update mode: keep any existing manual columns (e.g. notes) and carry them over by 'file'
if (file.exists(out_path)) {
  manifest_old <- readr::read_csv(out_path, show_col_types = FALSE)
  
  # Extra columns that exist in the old manifest but not in the new one
  extra_cols <- setdiff(names(manifest_old), names(manifest_new))
  
  # Keep notes + any extra columns from the old file (if they exist)
  carry_cols <- c("file", intersect(names(manifest_old), c("notes", extra_cols)))
  manifest_old_carry <- manifest_old %>% select(any_of(carry_cols))
  
  # Prefer parsed columns from the new manifest, but carry over manual fields
  manifest_updated <- manifest_new %>%
    select(-notes) %>%
    left_join(manifest_old_carry, by = "file")
  
  # If notes column did not exist before, create it
  if (!"notes" %in% names(manifest_updated)) {
    manifest_updated <- manifest_updated %>% mutate(notes = NA_character_)
  }
} else {
  manifest_updated <- manifest_new
}

# Final ordering (fresh -> dried -> aged, then line -> individual -> replicate)
manifest_updated <- manifest_updated %>%
  mutate(timepoint = factor(timepoint, levels = c("fresh", "dried", "aged"), ordered = TRUE)) %>%
  arrange(timepoint, line, individual, replicate, file) %>%
  mutate(timepoint = as.character(timepoint))

# Save
readr::write_csv(manifest_updated, out_path)

cat("Manifest updated:", nrow(manifest_updated), "rows ->", out_path, "\n")