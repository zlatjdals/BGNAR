#!/usr/bin/env Rscript

############################################################
## 11_collect_wind_results.R
## Aggregate and validate independent wind-origin checkpoints.
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
source(file.path(script_dir, "03_cluster_functions.R"))
source(file.path(script_dir, "09_wind_functions.R"))

cfg <- read_wind_config()
design <- load_wind_design(cfg)
result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")), mustWork = FALSE
)
checkpoint_dir <- file.path(result_root, "wind", "checkpoints")
summary_dir <- file.path(result_root, "wind", "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
files <- sort(list.files(checkpoint_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("No wind checkpoints in ", checkpoint_dir)

add_meta <- function(df, meta) {
  if (is.null(df) || !nrow(df)) return(data.frame())
  keys <- c("job_id", "origin_id", "origin_index", "window_start", "window_end",
            "window_length", "N", "horizon", "gnar_mode")
  for (key in rev(keys)) {
    value <- meta[[key]]
    if (is.null(value) || !length(value)) value <- NA
    df <- cbind(setNames(data.frame(value), key), df)
  }
  df
}
bind_nonempty <- function(x) {
  x <- x[vapply(x, function(z) !is.null(z) && nrow(z) > 0L, logical(1L))]
  if (!length(x)) data.frame() else do.call(rbind, x)
}

manifest <- fit_metrics <- overall <- horizon <- point <- coefficients <- hyperparameters <-
  orders <- failures <- list()
for (i in seq_along(files)) {
  z <- readRDS(files[[i]])
  meta <- z$metadata[1L, , drop = FALSE]
  meta$checkpoint <- files[[i]]
  meta$n_methods_ok <- length(z$methods)
  meta$n_methods_failed <- nrow(z$failures)
  manifest[[i]] <- meta
  if (nrow(z$failures)) failures[[length(failures) + 1L]] <- add_meta(z$failures, meta)
  for (method in names(z$methods)) {
    one <- z$methods[[method]]
    if (!is.null(one$fit) && nrow(as.data.frame(one$fit)) > 0L) {
      ff <- as.data.frame(one$fit)
      ff$method <- method
      fit_metrics[[length(fit_metrics) + 1L]] <- add_meta(ff, meta)
    }
    ov <- as.data.frame(as.list(one$forecast$overall), check.names = FALSE)
    ov$method <- method
    ov$elapsed_seconds <- one$elapsed_seconds
    overall[[length(overall) + 1L]] <- add_meta(ov, meta)
    hh <- one$forecast$horizon
    hh$method <- method
    horizon[[length(horizon) + 1L]] <- add_meta(hh, meta)
    pp <- one$forecast$point
    pp$method <- method
    pp$test_index <- meta$origin_index + pp$horizon
    pp$station <- design$stations$station[pp$node]
    point[[length(point) + 1L]] <- add_meta(pp, meta)
    if (nrow(one$coefficients)) {
      coefficients[[length(coefficients) + 1L]] <- add_meta(one$coefficients, meta)
    }
    if (nrow(one$hyperparameters)) {
      hyperparameters[[length(hyperparameters) + 1L]] <-
        add_meta(one$hyperparameters, meta)
    }
    orders[[length(orders) + 1L]] <- add_meta(one$selected_order, meta)
  }
}

manifest <- bind_nonempty(manifest)
fit_metrics <- bind_nonempty(fit_metrics)
overall <- bind_nonempty(overall)
horizon <- bind_nonempty(horizon)
point <- bind_nonempty(point)
coefficients <- bind_nonempty(coefficients)
hyperparameters <- bind_nonempty(hyperparameters)
orders <- bind_nonempty(orders)
failures <- bind_nonempty(failures)

write.csv(manifest, file.path(summary_dir, "wind_run_manifest.csv"), row.names = FALSE)
write.csv(fit_metrics, file.path(summary_dir, "wind_fit_by_origin.csv"), row.names = FALSE)
write.csv(overall, file.path(summary_dir, "wind_forecast_by_origin.csv"), row.names = FALSE)
write.csv(horizon, file.path(summary_dir, "wind_forecast_by_horizon.csv"), row.names = FALSE)
write.csv(point, file.path(summary_dir, "wind_point_forecasts.csv"), row.names = FALSE)
write.csv(orders, file.path(summary_dir, "wind_order_selection.csv"), row.names = FALSE)
write.csv(failures, file.path(summary_dir, "wind_failures.csv"), row.names = FALSE)
if (nrow(coefficients)) {
  write.csv(coefficients, file.path(summary_dir, "wind_coefficients.csv"), row.names = FALSE)
}
if (nrow(hyperparameters)) {
  write.csv(hyperparameters, file.path(summary_dir, "wind_hyperparameters.csv"),
            row.names = FALSE)
}

method_summary <- do.call(rbind, lapply(split(overall, overall$method), function(z) {
  finite_cov <- is.finite(z$coverage)
  data.frame(
    method = z$method[[1L]], origins = nrow(z),
    rmse = sqrt(mean(z$rmse^2)),
    mean_origin_rmse = mean(z$rmse),
    sd_origin_rmse = stats::sd(z$rmse),
    mae = mean(z$mae),
    coverage = if (any(finite_cov)) mean(z$coverage[finite_cov]) else NA_real_,
    avg_width = if (any(is.finite(z$avg_width))) {
      mean(z$avg_width[is.finite(z$avg_width)])
    } else NA_real_
  )
}))
method_summary <- method_summary[order(method_summary$rmse), ]
write.csv(method_summary, file.path(summary_dir, "wind_forecast_summary.csv"),
          row.names = FALSE)

if (nrow(fit_metrics)) {
  fit_summary <- do.call(rbind, lapply(split(fit_metrics, fit_metrics$method), function(z) {
    data.frame(
      method = z$method[[1L]], origins = nrow(z),
      fit_rmse = sqrt(mean(z$fit_rmse^2, na.rm = TRUE)),
      mean_origin_fit_rmse = mean(z$fit_rmse, na.rm = TRUE),
      sd_origin_fit_rmse = stats::sd(z$fit_rmse, na.rm = TRUE),
      fit_mae = mean(z$fit_mae, na.rm = TRUE),
      fit_start = unique(z$fit_start)[1L]
    )
  }))
  fit_summary <- fit_summary[order(fit_summary$fit_rmse), ]
} else {
  fit_summary <- data.frame()
}
write.csv(fit_summary, file.path(summary_dir, "wind_fit_summary.csv"), row.names = FALSE)

if (nrow(fit_summary)) {
  wind_performance <- merge(
    fit_summary, method_summary, by = c("method", "origins"), all = TRUE
  )
  names(wind_performance)[names(wind_performance) == "rmse"] <- "forecast_rmse"
  names(wind_performance)[names(wind_performance) == "mean_origin_rmse"] <-
    "mean_origin_forecast_rmse"
  names(wind_performance)[names(wind_performance) == "sd_origin_rmse"] <-
    "sd_origin_forecast_rmse"
  names(wind_performance)[names(wind_performance) == "mae"] <- "forecast_mae"
  names(wind_performance)[names(wind_performance) == "coverage"] <-
    "forecast_coverage"
  names(wind_performance)[names(wind_performance) == "avg_width"] <-
    "forecast_interval_length"
  wind_performance <- wind_performance[order(wind_performance$forecast_rmse), ]
} else {
  wind_performance <- data.frame()
}
write.csv(
  wind_performance,
  file.path(summary_dir, "wind_performance_summary.csv"),
  row.names = FALSE
)

winners <- do.call(rbind, lapply(split(overall, overall$origin_id), function(z) {
  z <- z[is.finite(z$rmse), , drop = FALSE]
  best <- min(z$rmse)
  z <- z[abs(z$rmse - best) < 1e-12, , drop = FALSE]
  data.frame(origin_id = z$origin_id, origin_index = z$origin_index,
             method = z$method, rmse = z$rmse, fractional_win = 1 / nrow(z))
}))
write.csv(winners, file.path(summary_dir, "wind_winner_by_origin.csv"), row.names = FALSE)
winner_counts <- stats::aggregate(fractional_win ~ method, winners, sum)
names(winner_counts)[2L] <- "wins"
write.csv(winner_counts, file.path(summary_dir, "wind_winner_counts.csv"), row.names = FALSE)

write.csv(design$stations, file.path(summary_dir, "wind_station_info.csv"), row.names = FALSE)
write.csv(design$adjacency, file.path(summary_dir, "wind_network_adjacency.csv"),
          row.names = FALSE)

expected_ids <- seq_len(min(cfg$wind_n_origins, length(design$origins)))
present_ids <- sort(unique(manifest$origin_id))
missing_ids <- setdiff(expected_ids, present_ids)
requested_methods <- cfg$methods
method_issues <- do.call(rbind, lapply(expected_ids, function(id) {
  found <- unique(overall$method[overall$origin_id == id])
  miss <- setdiff(requested_methods, found)
  if (!length(miss)) return(NULL)
  data.frame(origin_id = id, issue = "missing successful method",
             detail = paste(miss, collapse = ","))
}))
fit_issues <- do.call(rbind, lapply(expected_ids, function(id) {
  rows <- fit_metrics[fit_metrics$origin_id == id, , drop = FALSE]
  missing_fit <- setdiff(requested_methods, rows$method[is.finite(rows$fit_rmse)])
  if (!length(missing_fit)) return(NULL)
  data.frame(
    origin_id = id, issue = "missing/non-finite fit RMSE",
    detail = paste(missing_fit, collapse = ",")
  )
}))
issues <- rbind(
  if (length(missing_ids)) data.frame(
    origin_id = missing_ids, issue = "missing checkpoint", detail = ""
  ) else NULL,
  method_issues,
  fit_issues,
  if (nrow(failures)) data.frame(
    origin_id = failures$origin_id,
    issue = paste0("method failure: ", failures$method), detail = failures$error
  ) else NULL
)
if (is.null(issues)) issues <- data.frame(
  origin_id = integer(), issue = character(), detail = character()
)
write.csv(issues, file.path(summary_dir, "wind_validation_issues.csv"), row.names = FALSE)
writeLines(c(
  paste0("expected_origins=", length(expected_ids)),
  paste0("present_origins=", length(present_ids)),
  paste0("methods=", paste(requested_methods, collapse = ",")),
  paste0("validation_issues=", nrow(issues))
), file.path(summary_dir, "wind_validation_report.txt"))

message("Collected ", length(present_ids), " wind origins into ", summary_dir)
if (env_flag("STRICT_VALIDATION", TRUE) && nrow(issues)) {
  quit(save = "no", status = 2L)
}
