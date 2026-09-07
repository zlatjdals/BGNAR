#!/usr/bin/env Rscript

############################################################
## 05_collect_results.R
## Merge independent checkpoints; no concurrent file writes.
############################################################

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))

result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")),
  mustWork = FALSE
)
checkpoint_dir <- file.path(result_root, "checkpoints")
summary_dir <- file.path(result_root, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
files <- sort(list.files(checkpoint_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("No checkpoints in ", checkpoint_dir)

add_meta <- function(df, meta) {
  if (is.null(df) || !nrow(df)) return(data.frame())
  keys <- c("job_id", "scenario", "structure", "train_T", "replication",
            "N", "horizon", "intercept_mode", "gnar_mode")
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

manifest <- list()
fit_metrics <- list()
overall <- list()
horizon <- list()
point <- list()
coefficients <- list()
hyperparameters <- list()
orders <- list()
failures <- list()

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

    if (isTRUE(tolower(Sys.getenv("WRITE_POINT_FORECASTS", "0")) %in%
               c("1", "true", "yes"))) {
      pp <- one$forecast$point
      pp$method <- method
      point[[length(point) + 1L]] <- add_meta(pp, meta)
    }

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

write.csv(manifest, file.path(summary_dir, "run_manifest.csv"), row.names = FALSE)
write.csv(fit_metrics, file.path(summary_dir, "fit_by_run.csv"), row.names = FALSE)
write.csv(overall, file.path(summary_dir, "forecast_by_run.csv"), row.names = FALSE)
write.csv(horizon, file.path(summary_dir, "forecast_by_horizon.csv"), row.names = FALSE)
write.csv(orders, file.path(summary_dir, "order_selection_by_run.csv"), row.names = FALSE)
write.csv(failures, file.path(summary_dir, "failures.csv"), row.names = FALSE)
if (nrow(coefficients)) {
  write.csv(coefficients, file.path(summary_dir, "coefficient_by_run.csv"), row.names = FALSE)
}
if (nrow(hyperparameters)) {
  write.csv(hyperparameters, file.path(summary_dir, "hyperparameter_by_run.csv"),
            row.names = FALSE)
}
if (nrow(point)) {
  con <- gzfile(file.path(summary_dir, "point_forecasts.csv.gz"), open = "wt")
  utils::write.csv(point, con, row.names = FALSE)
  close(con)
}

summarize_forecasts <- function(df, group_cols) {
  key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(df, key), function(z) {
    first <- z[1L, group_cols, drop = FALSE]
    finite_cov <- is.finite(z$coverage)
    cbind(
      first,
      data.frame(
        runs = nrow(z),
        rmse_mean = mean(z$rmse, na.rm = TRUE),
        rmse_sd = stats::sd(z$rmse, na.rm = TRUE),
        rmse_median = stats::median(z$rmse, na.rm = TRUE),
        mae_mean = mean(z$mae, na.rm = TRUE),
        coverage_mean = if (any(finite_cov)) mean(z$coverage[finite_cov]) else NA_real_,
        width_mean = if (any(is.finite(z$avg_width))) {
          mean(z$avg_width[is.finite(z$avg_width)])
        } else NA_real_,
        elapsed_mean_seconds = mean(z$elapsed_seconds, na.rm = TRUE)
      )
    )
  })
  do.call(rbind, rows)
}

forecast_summary <- summarize_forecasts(
  overall, c("scenario", "structure", "train_T", "method")
)
write.csv(forecast_summary, file.path(summary_dir, "forecast_summary.csv"),
          row.names = FALSE)

if (nrow(fit_metrics)) {
  fit_key <- interaction(
    fit_metrics[c("scenario", "structure", "train_T", "method")],
    drop = TRUE, lex.order = TRUE
  )
  fit_summary <- do.call(rbind, lapply(split(fit_metrics, fit_key), function(z) {
    first <- z[1L, c("scenario", "structure", "train_T", "method"), drop = FALSE]
    cbind(
      first,
      data.frame(
        runs = nrow(z),
        fit_rmse_mean = mean(z$fit_rmse, na.rm = TRUE),
        fit_rmse_sd = stats::sd(z$fit_rmse, na.rm = TRUE),
        fit_mae_mean = mean(z$fit_mae, na.rm = TRUE),
        fit_start = unique(z$fit_start)[1L]
      )
    )
  }))
} else {
  fit_summary <- data.frame()
}
write.csv(fit_summary, file.path(summary_dir, "fit_summary.csv"), row.names = FALSE)

## Paired winner counts: each job is judged only among successfully fitted methods.
job_split <- split(overall, overall$job_id)
winners <- do.call(rbind, lapply(job_split, function(z) {
  z <- z[is.finite(z$rmse), , drop = FALSE]
  if (!nrow(z)) return(NULL)
  best <- min(z$rmse)
  z <- z[abs(z$rmse - best) < 1e-12, , drop = FALSE]
  data.frame(
    job_id = z$job_id, scenario = z$scenario, structure = z$structure,
    train_T = z$train_T, replication = z$replication,
    method = z$method, rmse = z$rmse,
    fractional_win = 1 / nrow(z)
  )
}))
write.csv(winners, file.path(summary_dir, "winner_by_run.csv"), row.names = FALSE)
winner_counts <- stats::aggregate(
  fractional_win ~ scenario + structure + train_T + method,
  data = winners, FUN = sum
)
names(winner_counts)[names(winner_counts) == "fractional_win"] <- "wins"
write.csv(winner_counts, file.path(summary_dir, "winner_counts.csv"), row.names = FALSE)

if (nrow(coefficients)) {
  coefficients$lag_group <- ifelse(is.na(coefficients$lag), 0L, coefficients$lag)
  coefficients$order_group <- ifelse(
    is.na(coefficients$neighbor_order), 0L, coefficients$neighbor_order
  )
  coef_key <- interaction(
    coefficients[c("scenario", "structure", "train_T", "method", "parameter",
                   "lag_group", "order_group")],
    drop = TRUE, lex.order = TRUE
  )
  coefficient_summary <- do.call(rbind, lapply(split(coefficients, coef_key), function(z) {
    first <- z[1L, c("scenario", "structure", "train_T", "method", "parameter",
                     "lag_group", "order_group"), drop = FALSE]
    cbind(
      first,
      data.frame(
        estimates = nrow(z),
        rmse = sqrt(mean(z$error^2, na.rm = TRUE)),
        bias = mean(z$error, na.rm = TRUE),
        coverage = mean(z$covered, na.rm = TRUE),
        mean_interval_width = mean(z$q975 - z$q025, na.rm = TRUE)
      )
    )
  }))
  names(coefficient_summary)[names(coefficient_summary) == "lag_group"] <- "lag"
  names(coefficient_summary)[names(coefficient_summary) == "order_group"] <-
    "neighbor_order"
  write.csv(coefficient_summary,
            file.path(summary_dir, "coefficient_summary.csv"), row.names = FALSE)

  ## Parameter recovery is first computed within each replication and then
  ## summarized across replications.  Active and truly-zero coefficients are
  ## kept separate so a large envelope cannot obtain an artificially small
  ## aggregate RMSE merely by containing many zeros.
  scored <- coefficients[coefficients$parameter %in% c("alpha", "beta"), , drop = FALSE]
  scored$active <- abs(scored$truth) > 1e-12
  score_frames <- list(
    transform(scored, score_set = "all"),
    transform(scored[scored$active, , drop = FALSE], score_set = "active"),
    transform(scored[!scored$active, , drop = FALSE], score_set = "zero")
  )
  score_frames <- score_frames[vapply(score_frames, nrow, integer(1L)) > 0L]
  scored <- do.call(rbind, score_frames)
  run_cols <- c(
    "job_id", "scenario", "structure", "train_T", "replication",
    "method", "parameter", "score_set"
  )
  parameter_key <- interaction(scored[run_cols], drop = TRUE, lex.order = TRUE)
  parameter_by_run <- do.call(rbind, lapply(split(scored, parameter_key), function(z) {
    first <- z[1L, run_cols, drop = FALSE]
    widths <- z$q975 - z$q025
    finite_coverage <- !is.na(z$covered)
    cbind(
      first,
      data.frame(
        n_parameters = nrow(z),
        rmse = sqrt(mean(z$error^2, na.rm = TRUE)),
        bias = mean(z$error, na.rm = TRUE),
        coverage = if (any(finite_coverage)) {
          mean(z$covered[finite_coverage])
        } else NA_real_,
        interval_length = if (any(is.finite(widths))) {
          mean(widths[is.finite(widths)])
        } else NA_real_
      )
    )
  }))
  write.csv(
    parameter_by_run, file.path(summary_dir, "parameter_by_run.csv"),
    row.names = FALSE
  )

  summary_cols <- c(
    "scenario", "structure", "train_T", "method", "parameter", "score_set"
  )
  parameter_summary_key <- interaction(
    parameter_by_run[summary_cols], drop = TRUE, lex.order = TRUE
  )
  parameter_summary <- do.call(rbind, lapply(
    split(parameter_by_run, parameter_summary_key),
    function(z) {
      first <- z[1L, summary_cols, drop = FALSE]
      finite_coverage <- is.finite(z$coverage)
      finite_width <- is.finite(z$interval_length)
      cbind(
        first,
        data.frame(
          runs = nrow(z),
          rmse_mean = mean(z$rmse, na.rm = TRUE),
          rmse_sd = stats::sd(z$rmse, na.rm = TRUE),
          bias_mean = mean(z$bias, na.rm = TRUE),
          coverage_mean = if (any(finite_coverage)) {
            mean(z$coverage[finite_coverage])
          } else NA_real_,
          interval_length_mean = if (any(finite_width)) {
            mean(z$interval_length[finite_width])
          } else NA_real_
        )
      )
    }
  ))
  write.csv(
    parameter_summary, file.path(summary_dir, "parameter_summary.csv"),
    row.names = FALSE
  )
}

if (nrow(hyperparameters)) {
  hp_key <- interaction(
    hyperparameters[c("scenario", "structure", "train_T", "method", "parameter")],
    drop = TRUE, lex.order = TRUE
  )
  hp_summary <- do.call(rbind, lapply(split(hyperparameters, hp_key), function(z) {
    first <- z[1L, c("scenario", "structure", "train_T", "method", "parameter"),
               drop = FALSE]
    cbind(first, data.frame(runs = nrow(z), posterior_mean = mean(z$mean),
                            posterior_median = mean(z$q500)))
  }))
  write.csv(hp_summary, file.path(summary_dir, "hyperparameter_summary.csv"),
            row.names = FALSE)
}

message("Collected ", nrow(manifest), " checkpoints into ", summary_dir)
message("Successful method fits: ", nrow(overall), "; failures: ", nrow(failures))
