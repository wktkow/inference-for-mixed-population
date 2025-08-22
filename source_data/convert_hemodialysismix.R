#!/usr/bin/env Rscript

suppressPackageStartupMessages({
})

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

ensure_package("haven")
suppressPackageStartupMessages(library(haven))

input_path <- "hemodialysismix.sas7bdat"
output_rds_path <- "hemodialysismix.rds"
output_csv_path <- "hemodialysismix.csv"
report_path <- "hemodialysismix_validation.txt"

if (!file.exists(input_path)) {
  stop(sprintf("Input file not found: %s", input_path))
}

read_warnings <- character(0)
dataset <- NULL

start_time <- Sys.time()
dataset <- withCallingHandlers(
  read_sas(input_path),
  warning = function(w) {
    read_warnings <<- c(read_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
read_ms <- as.integer(round(as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000))

validate_dataset <- function(df) {
  issues <- character(0)
  notes <- character(0)

  if (!is.data.frame(df)) {
    issues <- c(issues, "Object is not a data.frame/tibble.")
  }

  nr <- tryCatch(nrow(df), error = function(e) NA_integer_)
  nc <- tryCatch(ncol(df), error = function(e) NA_integer_)
  if (is.na(nr) || is.na(nc)) {
    issues <- c(issues, "Failed to compute dataset dimensions.")
  } else {
    if (nr <= 0) issues <- c(issues, sprintf("Row count is %d (expected > 0).", nr))
    if (nc <= 0) issues <- c(issues, sprintf("Column count is %d (expected > 0).", nc))
  }

  col_names <- names(df)
  if (anyDuplicated(col_names) > 0) {
    dupes <- unique(col_names[duplicated(col_names)])
    issues <- c(issues, sprintf("Duplicate column names: %s", paste(dupes, collapse = ", ")))
  }

  all_missing_cols <- names(df)[vapply(df, function(x) all(is.na(x)), logical(1))]
  if (length(all_missing_cols) > 0) {
    notes <- c(notes, sprintf("Columns with all values missing (%d): %s",
                              length(all_missing_cols), paste(all_missing_cols, collapse = ", ")))
  }

  dup_rows_count <- tryCatch(sum(duplicated(df)), error = function(e) NA_integer_)
  if (!is.na(dup_rows_count) && dup_rows_count > 0) {
    notes <- c(notes, sprintf("Found %d duplicate rows.", dup_rows_count))
  }

  classes <- vapply(df, function(x) paste(class(x), collapse = "/"), character(1))
  unusual <- names(classes)[!grepl("integer|numeric|double|character|logical|Date|POSIXct|haven_labelled", classes)]
  if (length(unusual) > 0) {
    notes <- c(notes, sprintf("Unusual column classes: %s", paste(sprintf("%s[%s]", unusual, classes[unusual]), collapse = ", ")))
  }

  list(
    ok = length(issues) == 0,
    issues = issues,
    notes = notes,
    nrow = nr,
    ncol = nc,
    classes = classes
  )
}

val <- validate_dataset(dataset)

save_warnings <- character(0)
save_ok <- FALSE
save_time <- Sys.time()
withCallingHandlers({
  saveRDS(dataset, file = output_rds_path, compress = "xz")
  save_ok <- TRUE
}, warning = function(w) {
  save_warnings <<- c(save_warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})
save_ms <- as.integer(round(as.numeric(difftime(Sys.time(), save_time, units = "secs")) * 1000))

csv_warnings <- character(0)
csv_ok <- FALSE
csv_time <- Sys.time()
to_csv_value <- function(x) {
  if (inherits(x, "haven_labelled")) {
    return(as.character(as_factor(x, levels = "labels")))
  }
  if (inherits(x, "POSIXct")) {
    return(format(x, "%Y-%m-%d %H:%M:%S%z"))
  }
  if (inherits(x, "Date")) {
    return(format(x, "%Y-%m-%d"))
  }
  x
}
dataset_csv <- as.data.frame(lapply(dataset, to_csv_value), stringsAsFactors = FALSE)
withCallingHandlers({
  write.csv(dataset_csv, file = output_csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  csv_ok <- TRUE
}, warning = function(w) {
  csv_warnings <<- c(csv_warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})
csv_ms <- as.integer(round(as.numeric(difftime(Sys.time(), csv_time, units = "secs")) * 1000))

con <- file(report_path, open = "wt", encoding = "UTF-8")
on.exit(close(con), add = TRUE)

write_line <- function(...) {
  cat(paste0(..., "\n"), file = con, append = TRUE)
}

write_line("hemodialysismix conversion and validation report")
write_line(sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
write_line(sprintf("Input: %s", normalizePath(input_path)))
write_line(sprintf("Output (RDS): %s", normalizePath(output_rds_path, mustWork = FALSE)))
write_line(sprintf("Output (CSV): %s", normalizePath(output_csv_path, mustWork = FALSE)))
in_size <- tryCatch(file.info(input_path)$size, error = function(e) NA_real_)
out_size <- tryCatch(file.info(output_rds_path)$size, error = function(e) NA_real_)
out_csv_size <- tryCatch(file.info(output_csv_path)$size, error = function(e) NA_real_)
write_line(sprintf("Input size (bytes): %s", ifelse(is.na(in_size), "NA", format(in_size, big.mark = ","))))
write_line(sprintf("Output size (bytes): %s", ifelse(is.na(out_size), "NA", format(out_size, big.mark = ","))))
write_line(sprintf("Output size (CSV bytes): %s", ifelse(is.na(out_csv_size), "NA", format(out_csv_size, big.mark = ","))))
write_line("")

write_line("Read:")
write_line(sprintf("- Duration: %d ms", read_ms))
if (length(read_warnings) > 0) {
  write_line(sprintf("- Warnings (%d):", length(read_warnings)))
  for (w in read_warnings) write_line(sprintf("  * %s", w))
} else {
  write_line("- Warnings: none")
}
write_line("")

write_line("Validation:")
write_line(sprintf("- Rows: %s", ifelse(is.na(val$nrow), "NA", format(val$nrow, big.mark = ","))))
write_line(sprintf("- Columns: %s", ifelse(is.na(val$ncol), "NA", format(val$ncol, big.mark = ","))))
write_line(sprintf("- Issues: %d", length(val$issues)))
for (i in val$issues) write_line(sprintf("  ! %s", i))
if (length(val$notes) > 0) {
  write_line(sprintf("- Notes: %d", length(val$notes)))
  for (n in val$notes) write_line(sprintf("  - %s", n))
}
write_line("")

write_line("Save RDS:")
write_line(sprintf("- Success: %s", save_ok))
write_line(sprintf("- Duration: %d ms", save_ms))
if (length(save_warnings) > 0) {
  write_line(sprintf("- Warnings (%d):", length(save_warnings)))
  for (w in save_warnings) write_line(sprintf("  * %s", w))
} else {
  write_line("- Warnings: none")
}

write_line("")
write_line("Save CSV:")
write_line(sprintf("- Success: %s", csv_ok))
write_line(sprintf("- Duration: %d ms", csv_ms))
if (length(csv_warnings) > 0) {
  write_line(sprintf("- Warnings (%d):", length(csv_warnings)))
  for (w in csv_warnings) write_line(sprintf("  * %s", w))
} else {
  write_line("- Warnings: none")
}

if (!val$ok) {
  message("Validation found issues. See report: ", report_path)
} else {
  message("Conversion complete. RDS written to: ", output_rds_path)
  message("CSV written to: ", output_csv_path)
  message("Validation report: ", report_path)
}


