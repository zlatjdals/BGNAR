#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))

local_lib <- file.path(script_dir, ".Rlib")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(local_lib, file.path(dirname(script_dir), ".Rlib"), .libPaths())))

required <- c("igraph", "GNAR", "BVAR")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
install_missing <- tolower(Sys.getenv("INSTALL_MISSING", "0")) %in%
  c("1", "true", "yes")
if (length(missing) && install_missing) {
  install.packages(
    missing, lib = local_lib,
    repos = Sys.getenv("CRAN_REPO", "https://cloud.r-project.org"),
    dependencies = TRUE
  )
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
}
if (length(missing)) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    ". Install them centrally or rerun with INSTALL_MISSING=1."
  )
}

cat("R version: ", R.version.string, "\n", sep = "")
cat("Library paths:\n", paste0("  ", .libPaths(), collapse = "\n"), "\n", sep = "")
for (pkg in required) {
  cat(pkg, " ", as.character(utils::packageVersion(pkg)), "\n", sep = "")
}
cat("Environment check passed.\n")
