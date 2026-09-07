#!/usr/bin/env Rscript

############################################################
## 15_make_wind_outputs.R
## Wind fit/forecast table and one two-panel rolling-origin figure.
############################################################

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))
result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")), mustWork = FALSE
)
summary_dir <- file.path(result_root, "wind", "summary")
report_dir <- file.path(summary_dir, "report_outputs")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

fit_path <- file.path(summary_dir, "wind_fit_by_origin.csv")
forecast_path <- file.path(summary_dir, "wind_forecast_by_origin.csv")
performance_path <- file.path(summary_dir, "wind_performance_summary.csv")
required <- c(fit_path, forecast_path, performance_path)
if (any(!file.exists(required))) {
  stop(
    "Missing metric-enabled wind summaries: ",
    paste(basename(required[!file.exists(required)]), collapse = ", "),
    ". Rerun wind checkpoints, then run 11_collect_wind_results.R."
  )
}

fit <- read.csv(fit_path, stringsAsFactors = FALSE)
forecast <- read.csv(forecast_path, stringsAsFactors = FALSE)
performance <- read.csv(performance_path, stringsAsFactors = FALSE)

method_label <- function(x) {
  out <- x
  out[out == "bgnar_nondummy"] <- "BGNAR (no dummy)"
  out[grepl("^bgnar", out) & out != "BGNAR (no dummy)"] <- "BGNAR"
  out[out == "gnar"] <- "GNAR"
  out[out == "bvar"] <- "BVAR"
  out[out == "rw"] <- "RW"
  out
}
fit$method_label <- method_label(fit$method)
forecast$method_label <- method_label(forecast$method)
performance$method_label <- method_label(performance$method)

write.csv(
  performance,
  file.path(report_dir, "wind_fit_forecast_table.csv"),
  row.names = FALSE
)

fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "--")
rows <- performance[order(match(performance$method_label,
                                c("BGNAR", "BGNAR (no dummy)", "GNAR", "BVAR", "RW"))), ]
tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Wind-data common-start fit RMSE and five-step forecast RMSE.}",
  "\\label{tab:wind-fit-forecast}",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Method & Fit RMSE & Forecast RMSE \\\\",
  "\\midrule"
)
body <- vapply(seq_len(nrow(rows)), function(i) {
  paste(
    rows$method_label[i], fmt(rows$fit_rmse[i]),
    fmt(rows$forecast_rmse[i]), sep = " & "
  )
}, character(1L))
body <- paste0(body, " \\\\")
writeLines(
  c(tex, body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(report_dir, "wind_fit_forecast_table.tex")
)

plot_df <- rbind(
  data.frame(
    origin_index = fit$origin_index, method = fit$method_label,
    metric = "Common-start fit RMSE", value = fit$fit_rmse
  ),
  data.frame(
    origin_index = forecast$origin_index, method = forecast$method_label,
    metric = "Five-step forecast RMSE", value = forecast$rmse
  )
)
method_levels <- c("BGNAR", "BGNAR (no dummy)", "GNAR", "BVAR", "RW")
method_levels <- method_levels[method_levels %in% unique(plot_df$method)]
plot_df$method <- factor(plot_df$method, levels = method_levels)
plot_df$metric <- factor(
  plot_df$metric,
  levels = c("Common-start fit RMSE", "Five-step forecast RMSE")
)

figure_path <- file.path(report_dir, "wind_fit_forecast_rmse_by_origin.pdf")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = origin_index, y = value, color = method, group = method)
  ) +
    ggplot2::geom_line(linewidth = 0.65, alpha = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~metric, ncol = 1L, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(
      BGNAR = "#0072B2", `BGNAR (no dummy)` = "#56B4E9",
      GNAR = "#D55E00", BVAR = "#009E73", RW = "#7F7F7F"
    )) +
    ggplot2::labs(x = "Forecast origin", y = "RMSE", color = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(figure_path, p, width = 8.0, height = 6.0)
} else {
  grDevices::pdf(figure_path, width = 8.0, height = 6.0)
  old_par <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4, 2, 1))
  on.exit({graphics::par(old_par); grDevices::dev.off()}, add = TRUE)
  colors <- c(
    BGNAR = "#0072B2", `BGNAR (no dummy)` = "#56B4E9",
    GNAR = "#D55E00", BVAR = "#009E73", RW = "#7F7F7F"
  )
  for (metric_name in levels(plot_df$metric)) {
    z <- plot_df[plot_df$metric == metric_name, , drop = FALSE]
    graphics::plot(
      range(z$origin_index), range(z$value, finite = TRUE), type = "n",
      xlab = "Forecast origin", ylab = "RMSE", main = metric_name
    )
    for (method_name in levels(plot_df$method)) {
      zz <- z[z$method == method_name, , drop = FALSE]
      graphics::lines(zz$origin_index, zz$value, type = "o", pch = 16,
                      col = colors[[method_name]])
    }
    graphics::legend("topright", levels(plot_df$method),
                     col = colors[levels(plot_df$method)],
                     lty = 1, pch = 16, bty = "n", horiz = TRUE)
  }
}

message("Wind table and figure written to ", report_dir)
