#!/usr/bin/env Rscript

############################################################
## 14_make_simulation_tables.R
## Three graph-specific performance tables and parameter tables.
############################################################

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))
result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")), mustWork = FALSE
)
summary_dir <- file.path(result_root, "summary")
report_dir <- file.path(summary_dir, "report_tables")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

required <- file.path(summary_dir, c(
  "fit_by_run.csv", "forecast_by_run.csv", "parameter_summary.csv"
))
if (any(!file.exists(required))) {
  stop(
    "Missing metric-enabled summaries: ",
    paste(basename(required[!file.exists(required)]), collapse = ", "),
    ". Rerun simulation checkpoints, then run 05_collect_results.R."
  )
}

fit <- read.csv(required[[1L]], stringsAsFactors = FALSE)
forecast <- read.csv(required[[2L]], stringsAsFactors = FALSE)
parameter <- read.csv(required[[3L]], stringsAsFactors = FALSE)

keys <- c("job_id", "scenario", "structure", "train_T", "replication", "method")
fit_keep <- fit[c(keys, "fit_rmse", "fit_mae")]
forecast_keep <- forecast[c(keys, "rmse", "mae", "coverage", "avg_width")]
names(forecast_keep)[names(forecast_keep) == "rmse"] <- "forecast_rmse"
names(forecast_keep)[names(forecast_keep) == "mae"] <- "forecast_mae"
perf <- merge(fit_keep, forecast_keep, by = keys, all = FALSE)

method_label <- function(x) {
  out <- x
  out[out == "bgnar_nondummy"] <- "BGNAR (no dummy)"
  out[grepl("^bgnar", out) & out != "BGNAR (no dummy)"] <- "BGNAR"
  out[out == "gnar"] <- "GNAR"
  out[out == "bvar"] <- "BVAR"
  out[out == "rw"] <- "RW"
  out
}

safe_sd <- function(x) if (sum(is.finite(x)) > 1L) stats::sd(x, na.rm = TRUE) else NA_real_
safe_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
group_cols <- c("scenario", "structure", "train_T", "method")
perf_key <- interaction(perf[group_cols], drop = TRUE, lex.order = TRUE)
perf_summary <- do.call(rbind, lapply(split(perf, perf_key), function(z) {
  first <- z[1L, group_cols, drop = FALSE]
  bayesian <- grepl("^bgnar", first$method) || first$method == "bvar"
  cbind(
    first,
    data.frame(
      runs = nrow(z),
      fit_rmse_mean = safe_mean(z$fit_rmse),
      fit_rmse_sd = safe_sd(z$fit_rmse),
      forecast_rmse_mean = safe_mean(z$forecast_rmse),
      forecast_rmse_sd = safe_sd(z$forecast_rmse),
      forecast_interval_length = if (bayesian) safe_mean(z$avg_width) else NA_real_,
      forecast_coverage = if (bayesian) safe_mean(z$coverage) else NA_real_
    )
  )
}))
perf_summary$method_label <- method_label(perf_summary$method)
write.csv(
  perf_summary, file.path(report_dir, "simulation_performance_all.csv"),
  row.names = FALSE
)

fmt <- function(x, digits = 3L) {
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "--")
}
fmt_mean_sd <- function(m, s) {
  ifelse(is.finite(m), paste0(fmt(m), " (", fmt(s), ")"), "--")
}
tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x
}
bold_min <- function(value, label, group) {
  out <- label
  for (g in unique(group)) {
    ii <- which(group == g & is.finite(value))
    if (!length(ii)) next
    best <- min(value[ii])
    hit <- ii[abs(value[ii] - best) < 1e-12]
    out[hit] <- paste0("\\textbf{", out[hit], "}")
  }
  out
}

write_performance_tex <- function(z, network, path) {
  z <- z[order(z$scenario, z$train_T,
               match(z$method_label, c("BGNAR", "BGNAR (no dummy)", "GNAR", "BVAR", "RW"))), ]
  group <- interaction(z$scenario, z$train_T, drop = TRUE)
  fit_text <- fmt_mean_sd(z$fit_rmse_mean, z$fit_rmse_sd)
  forecast_text <- fmt_mean_sd(z$forecast_rmse_mean, z$forecast_rmse_sd)
  fit_text <- bold_min(z$fit_rmse_mean, fit_text, group)
  forecast_text <- bold_min(z$forecast_rmse_mean, forecast_text, group)
  lines <- c(
    "\\begin{longtable}{lllrrrr}",
    paste0("\\caption{", network,
           " network: common-start fit RMSE and five-step forecast performance.}",
           "\\label{tab:performance-", tolower(network), "}\\\\"),
    "\\toprule",
    "DGP & $T$ & Method & Fit RMSE & Forecast RMSE & Interval length & Coverage \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    "DGP & $T$ & Method & Fit RMSE & Forecast RMSE & Interval length & Coverage \\\\",
    "\\midrule",
    "\\endhead"
  )
  body <- vapply(seq_len(nrow(z)), function(i) {
    paste(
      tex_escape(z$scenario[i]), z$train_T[i], tex_escape(z$method_label[i]),
      fit_text[i], forecast_text[i],
      fmt(z$forecast_interval_length[i]), fmt(z$forecast_coverage[i]),
      sep = " & "
    )
  }, character(1L))
  body <- paste0(body, " \\\\")
  writeLines(c(lines, body, "\\bottomrule", "\\end{longtable}"), path)
}

for (network in c("ER", "SBM", "SWN")) {
  z <- perf_summary[perf_summary$structure == network, , drop = FALSE]
  write.csv(
    z, file.path(report_dir, paste0("simulation_performance_", network, ".csv")),
    row.names = FALSE
  )
  write_performance_tex(
    z, network,
    file.path(report_dir, paste0("simulation_performance_", network, ".tex"))
  )
}

## Parameter recovery: BGNAR and GNAR only.  CSV files retain all/active/zero
## score sets; the compact LaTeX table reports active coefficients.
parameter <- parameter[
  grepl("^bgnar", parameter$method) | parameter$method == "gnar",
  , drop = FALSE
]
parameter$method_label <- method_label(parameter$method)
write.csv(
  parameter, file.path(report_dir, "simulation_parameters_all.csv"),
  row.names = FALSE
)

write_parameter_tex <- function(z, network, path) {
  z <- z[z$score_set == "active", , drop = FALSE]
  z <- z[order(z$scenario, z$train_T, z$parameter,
               match(z$method_label, c("BGNAR", "BGNAR (no dummy)", "GNAR"))), ]
  group <- interaction(z$scenario, z$train_T, z$parameter, drop = TRUE)
  rmse_text <- bold_min(
    z$rmse_mean, fmt_mean_sd(z$rmse_mean, z$rmse_sd), group
  )
  lines <- c(
    "\\begin{longtable}{lllrrrr}",
    paste0("\\caption{", network,
           " network: active-parameter recovery.}",
           "\\label{tab:parameters-", tolower(network), "}\\\\"),
    "\\toprule",
    "DGP & $T$ & Parameter & Method & RMSE & Interval length & Coverage \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    "DGP & $T$ & Parameter & Method & RMSE & Interval length & Coverage \\\\",
    "\\midrule",
    "\\endhead"
  )
  body <- vapply(seq_len(nrow(z)), function(i) {
    paste(
      tex_escape(z$scenario[i]), z$train_T[i], tex_escape(z$parameter[i]),
      tex_escape(z$method_label[i]), rmse_text[i],
      fmt(z$interval_length_mean[i]), fmt(z$coverage_mean[i]),
      sep = " & "
    )
  }, character(1L))
  body <- paste0(body, " \\\\")
  writeLines(c(lines, body, "\\bottomrule", "\\end{longtable}"), path)
}

for (network in c("ER", "SBM", "SWN")) {
  z <- parameter[parameter$structure == network, , drop = FALSE]
  write.csv(
    z, file.path(report_dir, paste0("simulation_parameters_", network, ".csv")),
    row.names = FALSE
  )
  write_parameter_tex(
    z, network,
    file.path(report_dir, paste0("simulation_parameters_", network, ".tex"))
  )
}

message("Simulation tables written to ", report_dir)
