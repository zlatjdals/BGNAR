#!/usr/bin/env Rscript

############################################################
## 10_run_wind_one.R
## One rolling-origin wind-speed analysis job.
############################################################

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))

.libPaths(unique(c(
  file.path(script_dir, ".Rlib"), file.path(dirname(script_dir), ".Rlib"),
  .libPaths()
)))
source(file.path(script_dir, "00_utils.R"))
source(file.path(script_dir, "01_bgnar_model.R"))
source(file.path(script_dir, "02_compare_methods.R"))
source(file.path(script_dir, "03_cluster_functions.R"))
source(file.path(script_dir, "09_wind_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 10_run_wind_one.R <origin_id>")
origin_id <- suppressWarnings(as.integer(args[[1L]]))
if (!is.finite(origin_id) || origin_id < 1L) stop("origin_id must be positive.")

cfg <- read_wind_config()
required <- unique(c(
  "igraph", "GNAR",
  if ("bvar" %in% cfg$methods) "BVAR" else character()
))
missing_packages <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

design <- load_wind_design(cfg)
cfg$N <- cfg$wind_n_stations
dat <- wind_window_data(design, origin_id, cfg)

result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")), mustWork = FALSE
)
checkpoint_dir <- file.path(result_root, "wind", "checkpoints")
full_fit_dir <- file.path(result_root, "wind", "full_fits")
spatial_mean_dir <- file.path(result_root, "wind", "spatial_mean")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_mean_dir, recursive = TRUE, showWarnings = FALSE)
if (cfg$save_full_fits) dir.create(full_fit_dir, recursive = TRUE, showWarnings = FALSE)

job_id <- sprintf("wind_origin_%03d", origin_id)
checkpoint <- file.path(checkpoint_dir, paste0(job_id, ".rds"))
if (file.exists(checkpoint) && !cfg$overwrite) {
  message("Checkpoint exists; skipping: ", checkpoint)
  quit(save = "no", status = 0L)
}

started <- Sys.time()
message(
  "Starting ", job_id, " at observation ", dat$origin_index,
  " with methods: ", paste(cfg$methods, collapse = ", ")
)
compact_methods <- list()
failures <- list()
spatial_mean_rows <- list()
for (method in cfg$methods) {
  message("  fitting ", method)
  seed <- wind_method_seed(cfg, origin_id, method)
  tick <- proc.time()[[3L]]
  fit_result <- tryCatch(
    fit_requested_method(method, dat, cfg$wind_window, cfg, seed),
    error = function(e) e
  )
  elapsed <- proc.time()[[3L]] - tick
  if (inherits(fit_result, "error") || is.null(fit_result)) {
    msg <- if (inherits(fit_result, "error")) conditionMessage(fit_result) else "NULL fit"
    failures[[length(failures) + 1L]] <- data.frame(
      method = method, error = msg, elapsed_seconds = elapsed
    )
    message("  FAILED ", method, ": ", msg)
    next
  }
  compact_methods[[method]] <- compact_wind_method_result(
    fit_result, dat, cfg, method, elapsed
  )
  spatial_summary <- compact_methods[[method]]$spatial_mean_predictive
  if (nrow(spatial_summary)) {
    if (nrow(spatial_summary) != cfg$horizon) {
      stop("Spatial-mean predictive summary has the wrong horizon length.")
    }
    spatial_mean_rows[[method]] <- transform(
      spatial_summary,
      job_id = job_id,
      origin_id = origin_id,
      origin_index = dat$origin_index,
      test_index = dat$test_index,
      method = method
    )
  }
  if (cfg$save_full_fits) {
    full_path <- file.path(full_fit_dir, paste0(job_id, "_", method, ".rds"))
    tmp_full <- tempfile(pattern = paste0(job_id, "_"), tmpdir = full_fit_dir)
    saveRDS(fit_result, tmp_full, compress = "xz")
    if (!file.rename(tmp_full, full_path)) stop("Could not finalize ", full_path)
  }
  rm(fit_result)
  invisible(gc(FALSE))
}

failure_df <- if (length(failures)) do.call(rbind, failures) else data.frame(
  method = character(), error = character(), elapsed_seconds = numeric()
)
result <- list(
  metadata = data.frame(
    job_id = job_id, origin_id = origin_id,
    origin_index = dat$origin_index,
    window_start = dat$window_start, window_end = dat$window_end,
    window_length = cfg$wind_window,
    N = cfg$N, horizon = cfg$horizon,
    gnar_mode = cfg$gnar_mode,
    started = format(started, "%Y-%m-%d %H:%M:%S %z"),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  ),
  config = list(
    methods = cfg$methods, gnar_mode = cfg$gnar_mode,
    preprocessing = "raw data; training-window nodewise centering for centered GNAR only",
    bgnar_p = cfg$p_fit, bgnar_s = cfg$s_fit,
    alpha_mode = "hierarchical local alpha",
    reduce_soc = cfg$reduce_soc,
    fixed_tau2 = cfg$fixed_tau2, fixed_tau_scale = cfg$fixed_tau_scale,
    adaptive_tau_a = cfg$adaptive_a, adaptive_tau_b = cfg$adaptive_b,
    spatial_mean_interval = paste(
      "For each posterior predictive draw, average across stations;",
      "then take the 0.025 and 0.975 quantiles across drawwise means."
    )
  ),
  methods = compact_methods,
  failures = failure_df,
  package_versions = vapply(
    required, function(pkg) as.character(utils::packageVersion(pkg)), character(1L)
  )
)

tmp_checkpoint <- tempfile(pattern = paste0(job_id, "_"), tmpdir = checkpoint_dir)
saveRDS(result, tmp_checkpoint, compress = "xz")
if (!file.rename(tmp_checkpoint, checkpoint)) stop("Could not finalize ", checkpoint)

## Write one origin-specific CSV.  Unique files avoid concurrent-write races
## when the 15 rolling origins are run in parallel.  The same information is
## also retained inside each compact checkpoint.
if (length(spatial_mean_rows)) {
  spatial_df <- do.call(rbind, spatial_mean_rows)
  keep <- c(
    "job_id", "origin_id", "origin_index", "test_index", "method",
    "horizon", "mean", "q025", "q500", "q975"
  )
  spatial_path <- file.path(
    spatial_mean_dir, paste0(job_id, "_spatial_mean.csv")
  )
  tmp_spatial <- tempfile(pattern = paste0(job_id, "_"),
                          tmpdir = spatial_mean_dir)
  utils::write.csv(spatial_df[keep], tmp_spatial, row.names = FALSE)
  if (!file.rename(tmp_spatial, spatial_path)) {
    stop("Could not finalize ", spatial_path)
  }
}
message("Completed ", job_id, ": ", checkpoint)
if (nrow(failure_df)) quit(save = "no", status = 2L)
