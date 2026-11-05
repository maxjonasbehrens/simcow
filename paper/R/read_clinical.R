# Clinical application (SENIOR cohort). The data are not public and must never enter this
# repository: scripts take the data path as a command-line argument and write every output to
# a directory outside the repository.

# Patients aged 75 or older; the outcome (radiotherapy dose) enters the models squared;
# numeric columns with fewer than five distinct values become factors.
read_clinical <- function(path) {
  d <- readr::read_csv2(path, show_col_types = FALSE) |>
    dplyr::filter(Age >= 75) |>
    dplyr::mutate(DosisEQDW = DosisEQDW^2)
  nm <- names(d)
  d[] <- lapply(nm, function(col) { v <- d[[col]]; if (is.numeric(v) && length(unique(v)) < 5) as.factor(v) else v })
  names(d) <- nm
  d
}

# Stop if the output directory lies inside the repository; otherwise create it
outdir_outside_repo <- function(dir, repo) {
  d <- if (dir.exists(dir)) normalizePath(dir) else file.path(normalizePath(dirname(dir), mustWork = TRUE), basename(dir))
  r <- normalizePath(repo)
  inside <- tolower(d) == tolower(r) || startsWith(tolower(d), paste0(tolower(r), .Platform$file.sep))
  if (inside) stop("The output directory must lie outside the repository (patient-level output): ", dir)
  dir.create(d, showWarnings = FALSE)
  d
}
