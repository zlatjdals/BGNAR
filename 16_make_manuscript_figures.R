#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))
input_root <- if (length(args) >= 1L) args[[1L]] else
  file.path(script_dir, "results")
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(input_root, "figures")
wind_input_root <- if (length(args) >= 3L) args[[3L]] else
  file.path(input_root, "wind")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_label <- function(x) {
  out <- x
  out[out == "bgnar_soc_dio_adaptive"] <- "BGNAR"
  out[out == "bvar"] <- "BVAR"
  out[out == "gnar"] <- "GNAR"
  out
}

method_order <- c("BGNAR", "GNAR", "BVAR")
method_colors <- c(BGNAR = "#0072B2", GNAR = "#D55E00", BVAR = "#009E73")
method_pch <- c(BGNAR = 16, GNAR = 17, BVAR = 15)

draw_method_lines <- function(z, x_name, y_name, ylim = NULL,
                              xlab = "Training length", ylab = "RMSE",
                              add_legend = FALSE) {
  x_values <- sort(unique(z[[x_name]]))
  if (is.null(ylim)) ylim <- range(z[[y_name]], finite = TRUE)
  graphics::plot(
    range(x_values), ylim, type = "n", xaxt = "n",
    xlab = xlab, ylab = ylab, bty = "l"
  )
  graphics::axis(1, at = x_values)
  for (method in method_order) {
    zz <- z[z$Method == method, , drop = FALSE]
    if (!nrow(zz)) next
    zz <- zz[order(zz[[x_name]]), , drop = FALSE]
    graphics::lines(
      zz[[x_name]], zz[[y_name]], type = "o",
      col = method_colors[[method]], pch = method_pch[[method]],
      lwd = 1.5, cex = 0.75
    )
  }
  graphics::grid(col = "grey90", lty = 1)
  if (add_legend) {
    shown <- method_order[method_order %in% z$Method]
    graphics::legend(
      "topright", shown, col = method_colors[shown], pch = method_pch[shown],
      lty = 1, lwd = 1.5, bty = "n", cex = 0.75
    )
  }
}

sim_dir <- file.path(input_root, "summary")
forecast <- read.csv(file.path(sim_dir, "forecast_by_run.csv"),
                     stringsAsFactors = FALSE)
parameter <- read.csv(file.path(sim_dir, "parameter_by_run.csv"),
                      stringsAsFactors = FALSE)
coefficient <- read.csv(file.path(sim_dir, "coefficient_by_run.csv"),
                        stringsAsFactors = FALSE)

forecast$Method <- method_label(forecast$method)
fplot <- aggregate(rmse ~ scenario + structure + train_T + Method,
                   forecast, mean)

grDevices::pdf(file.path(output_dir, "simulation_forecast_rmse.pdf"),
               width = 9.2, height = 6.7, family = "Helvetica")
old_par <- graphics::par(
  mfrow = c(3, 4), mar = c(2.8, 3.1, 2.2, 0.8), oma = c(2.5, 2.8, 0.5, 0.3),
  mgp = c(1.8, 0.55, 0), tcl = -0.25, cex.axis = 0.72, cex.lab = 0.78
)
for (structure in c("ER", "SBM", "SWN")) {
  for (scenario in c("M11", "M12", "M13", "M2")) {
    z <- fplot[fplot$structure == structure & fplot$scenario == scenario, ]
    draw_method_lines(z, "train_T", "rmse", xlab = "", ylab = "",
                      add_legend = structure == "ER" && scenario == "M11")
    graphics::title(main = paste(structure, scenario, sep = " - "),
                    cex.main = 0.85)
  }
}
graphics::mtext("Training length", side = 1, outer = TRUE, line = 1.0)
graphics::mtext("Five-step forecast RMSE", side = 2, outer = TRUE, line = 1.2)
graphics::par(old_par)
grDevices::dev.off()

parameter <- subset(parameter,
                    parameter %in% c("alpha", "beta") &
                      score_set %in% c("active", "zero", "all") &
                      method %in% c("bgnar_soc_dio_adaptive", "gnar"))
parameter$Method <- method_label(parameter$method)
pplot <- aggregate(rmse ~ structure + train_T + Method + parameter + score_set,
                   parameter, mean)

for (parameter_name in c("alpha", "beta")) {
  grDevices::pdf(
    file.path(output_dir, paste0("simulation_parameter_rmse_", parameter_name, ".pdf")),
    width = 8.7, height = 6.3, family = "Helvetica"
  )
  old_par <- graphics::par(
    mfrow = c(3, 3), mar = c(2.8, 3.2, 2.2, 0.8),
    oma = c(2.5, 2.8, 0.5, 0.3), mgp = c(1.8, 0.55, 0),
    tcl = -0.25, cex.axis = 0.72, cex.lab = 0.78
  )
  for (structure in c("ER", "SBM", "SWN")) {
    for (score_name in c("active", "zero", "all")) {
      z <- pplot[
        pplot$structure == structure & pplot$parameter == parameter_name &
          pplot$score_set == score_name,
      ]
      draw_method_lines(z, "train_T", "rmse", xlab = "", ylab = "",
                        add_legend = structure == "ER" && score_name == "active")
      label <- switch(score_name, active = "active", zero = "inactive",
                      all = "full envelope")
      graphics::title(
        main = bquote(.(structure) ~ .(parameter_name) ~ "(" * .(label) * ")"),
        cex.main = 0.78
      )
    }
  }
  graphics::mtext("Training length", side = 1, outer = TRUE, line = 1.0)
  graphics::mtext("Parameter RMSE", side = 2, outer = TRUE, line = 1.2)
  graphics::par(old_par)
  grDevices::dev.off()
}

one <- subset(
  coefficient,
  job_id == "M2_SWN_T040_r001" & parameter %in% c("alpha", "beta") &
    method %in% c("bgnar_soc_dio_adaptive", "gnar")
)
grDevices::pdf(file.path(output_dir, "single_replication_parameter_recovery.pdf"),
               width = 12.0, height = 5.9, family = "Helvetica")
old_par <- graphics::par(
  mfrow = c(2, 5), mar = c(3.5, 3.0, 2.1, 0.5), oma = c(0.4, 2.2, 0.3, 0.2),
  mgp = c(1.8, 0.55, 0), tcl = -0.25, cex.axis = 0.62, cex.lab = 0.70
)
for (parameter_name in c("alpha", "beta")) {
  for (lag in 1:5) {
    zb <- one[one$parameter == parameter_name & one$lag == lag &
                one$method == "bgnar_soc_dio_adaptive", ]
    zg <- one[one$parameter == parameter_name & one$lag == lag &
                one$method == "gnar", ]
    if (parameter_name == "alpha") {
      zb <- zb[order(zb$node), ]
      zg <- zg[order(zg$node), ]
      x <- zb$node
      xlab <- "Node"
    } else {
      zb <- zb[order(zb$neighbor_order), ]
      zg <- zg[order(zg$neighbor_order), ]
      x <- zb$neighbor_order
      xlab <- "Network order"
    }
    ylim <- range(c(zb$q025, zb$q975, zg$mean, zb$truth), finite = TRUE)
    graphics::plot(
      range(x), ylim, type = "n", xlab = xlab, ylab = "", bty = "l",
      xaxt = "n"
    )
    graphics::axis(1, at = x)
    if (parameter_name == "alpha") {
      graphics::abline(h = unique(zb$truth), col = "grey45", lty = 2, lwd = 1)
    } else {
      graphics::lines(x, zb$truth, col = "grey45", lty = 2, lwd = 1)
    }
    graphics::segments(x, zb$q025, x, zb$q975, col = method_colors[["BGNAR"]])
    graphics::points(x - 0.06, zb$mean, pch = 16,
                     col = method_colors[["BGNAR"]], cex = 0.72)
    graphics::points(x + 0.06, zg$mean, pch = 17,
                     col = method_colors[["GNAR"]], cex = 0.72)
    graphics::grid(col = "grey92")
    graphics::title(main = paste(parameter_name, "lag", lag), cex.main = 0.85)
    if (parameter_name == "alpha" && lag == 1L) {
      graphics::legend(
        "bottomright", c("BGNAR mean and 95% CI", "GNAR estimate", "truth"),
        col = c(method_colors[["BGNAR"]], method_colors[["GNAR"]], "grey45"),
        pch = c(16, 17, NA), lty = c(1, NA, 2), bty = "n", cex = 0.63
      )
    }
  }
}
graphics::mtext("Coefficient estimate", side = 2, outer = TRUE, line = 1.0)
graphics::par(old_par)
grDevices::dev.off()

wind_dir <- if (dir.exists(file.path(wind_input_root, "summary"))) {
  file.path(wind_input_root, "summary")
} else {
  file.path(wind_input_root, "wind", "summary")
}
w_origin <- read.csv(file.path(wind_dir, "wind_forecast_by_origin.csv"),
                     stringsAsFactors = FALSE)
w_point <- read.csv(file.path(wind_dir, "wind_point_forecasts.csv"),
                    stringsAsFactors = FALSE)
w_station <- read.csv(file.path(wind_dir, "wind_station_info.csv"),
                      stringsAsFactors = FALSE)
w_coef <- read.csv(file.path(wind_dir, "wind_coefficients.csv"),
                   stringsAsFactors = FALSE)
w_adj <- as.matrix(read.csv(file.path(wind_dir, "wind_network_adjacency.csv"),
                            check.names = FALSE))
storage.mode(w_adj) <- "numeric"
w_origin$Method <- method_label(w_origin$method)

local_library_paths <- c(
  file.path(script_dir, ".Rlib"),
  file.path(dirname(script_dir), ".Rlib")
)
.libPaths(unique(c(local_library_paths[dir.exists(local_library_paths)], .libPaths())))
wind_env <- new.env(parent = emptyenv())
utils::data("vswind", package = "GNAR", envir = wind_env)
wind_full <- as.matrix(wind_env$vswindts)
wind_names <- make.unique(trimws(wind_env$vswindnames))
station_columns <- match(w_station$station, wind_names)
if (anyNA(station_columns)) {
  stop("Could not match all summarized stations to GNAR::vswind.")
}
actual_index <- seq.int(min(w_point$test_index), max(w_point$test_index))

grDevices::pdf(file.path(output_dir, "wind_analysis.pdf"),
               width = 9.2, height = 3.8, family = "Helvetica")
old_par <- graphics::par(
  mfrow = c(1, 2), mar = c(3.8, 4.1, 2.5, 0.8),
  mgp = c(2.3, 0.72, 0), tcl = -0.25,
  cex.axis = 0.84, cex.lab = 0.94
)

graphics::plot(
  w_station$x, w_station$y, type = "n", axes = FALSE, xlab = "", ylab = "",
  asp = 1, main = "(a) Wind-station network", cex.main = 1.03
)
edge_idx <- which(upper.tri(w_adj) & w_adj > 0, arr.ind = TRUE)
graphics::segments(
  w_station$x[edge_idx[, 1L]], w_station$y[edge_idx[, 1L]],
  w_station$x[edge_idx[, 2L]], w_station$y[edge_idx[, 2L]],
  col = "grey75", lwd = 0.7
)
node_cex <- 0.8 + 0.25 * w_station$degree
graphics::points(w_station$x, w_station$y, pch = 21, bg = "#56B4E9",
                 col = "#005A8D", cex = node_cex)
capel_index <- match("CAPEL", w_station$station)
if (!is.na(capel_index)) {
  graphics::points(
    w_station$x[capel_index], w_station$y[capel_index],
    pch = 21, bg = "#E69F00", col = "#8C510A",
    lwd = 1.6, cex = node_cex[capel_index] + 0.45
  )
  graphics::text(
    w_station$x[capel_index], w_station$y[capel_index], "CAPEL",
    pos = 4, offset = 0.65, cex = 0.78, font = 2, col = "#6B3A00"
  )
}
graphics::box(col = "grey70")

yr <- range(w_origin$rmse, finite = TRUE)
graphics::plot(
  range(w_origin$origin_index), yr, type = "n",
  xlab = "Forecast origin", ylab = "Five-step RMSE", bty = "l",
  main = "(b) Origin-specific forecast error", cex.main = 1.03
)
method_offset <- c(BGNAR = -0.9, GNAR = 0, BVAR = 0.9)
for (method in c(setdiff(method_order, "BGNAR"), "BGNAR")) {
  z <- w_origin[w_origin$Method == method, ]
  if (!nrow(z)) next
  z <- z[order(z$origin_index), ]
  graphics::points(
    z$origin_index + method_offset[[method]], z$rmse,
    col = method_colors[[method]], pch = method_pch[[method]],
    cex = if (method == "BGNAR") 1.05 else 0.92
  )
}
graphics::grid(col = "grey90")
graphics::legend(
  "topright", method_order, col = method_colors[method_order],
  pch = method_pch[method_order], lty = 0, bty = "n", cex = 0.78
)
graphics::par(old_par)
grDevices::dev.off()

w_point$Method <- method_label(w_point$method)
draw_station_forecast <- function(selected_station, output_file) {
  if (!selected_station %in% w_station$station) {
    stop(selected_station, " is not present in the selected wind-station network.")
  }
  w_station_point <- w_point[w_point$station == selected_station, ]
  w_station_bgnar <- w_station_point[
    w_station_point$Method == "BGNAR" &
      is.finite(w_station_point$lo) & is.finite(w_station_point$hi),
  ]
  selected_station_column <- station_columns[
    match(selected_station, w_station$station)
  ]
  w_station_actual <- data.frame(
    test_index = actual_index,
    truth = wind_full[actual_index, selected_station_column]
  )

  grDevices::pdf(file.path(output_dir, output_file),
                 width = 9.2, height = 3.9, family = "Helvetica")
  old_par <- graphics::par(
    mar = c(3.8, 4.0, 2.2, 1.0), mgp = c(2.15, 0.65, 0),
    tcl = -0.25, cex.axis = 0.78, cex.lab = 0.84
  )
  graphics::plot(
    range(w_station_actual$test_index),
    range(c(w_station_actual$truth, w_station_point$prediction,
            w_station_bgnar$lo, w_station_bgnar$hi), finite = TRUE),
    type = "n", xlab = "Observation index",
    ylab = "Wind-speed series (package scale)",
    bty = "l", main = ""
  )
  graphics::grid(col = "grey90")
  for (origin in sort(unique(w_station_bgnar$origin_id))) {
    zi <- w_station_bgnar[w_station_bgnar$origin_id == origin, ]
    zi <- zi[order(zi$test_index), ]
    graphics::polygon(
      c(zi$test_index, rev(zi$test_index)), c(zi$lo, rev(zi$hi)),
      col = grDevices::adjustcolor(method_colors[["BGNAR"]], alpha.f = 0.16),
      border = NA
    )
  }
  for (method in c(setdiff(method_order, "BGNAR"), "BGNAR")) {
    z <- w_station_point[w_station_point$Method == method, ]
    if (!nrow(z)) next
    for (origin in sort(unique(z$origin_id))) {
      zz <- z[z$origin_id == origin, ]
      zz <- zz[order(zz$test_index), ]
      graphics::lines(
        zz$test_index, zz$prediction, col = method_colors[[method]],
        lwd = if (method == "BGNAR") 2.3 else 1.45
      )
    }
  }
  graphics::lines(w_station_actual$test_index, w_station_actual$truth,
                  col = "black", lwd = 2.3)
  plot_usr <- graphics::par("usr")
  graphics::legend(
    x = mean(plot_usr[1:2]),
    y = plot_usr[4] + 0.10 * diff(plot_usr[3:4]),
    c("Observed", method_order, "BGNAR 95% interval"),
    col = c("black", method_colors[method_order],
            grDevices::adjustcolor(method_colors[["BGNAR"]], alpha.f = 0.35)),
    lty = c(rep(1, 1 + length(method_order)), NA),
    lwd = c(2.5, 2.3, 1.45, 1.45, NA),
    pch = c(rep(NA, 1 + length(method_order)), 15),
    pt.cex = 1.5, bty = "n", cex = 0.72,
    horiz = TRUE, xjust = 0.5, yjust = 0.5, xpd = NA
  )
  graphics::par(old_par)
  grDevices::dev.off()
}

draw_station_forecast("CAPEL", "wind_station_forecast.pdf")

# Three-panel BGNAR forecast diagnostics across all stations.
w_diag <- w_point[w_point$Method == "BGNAR", ]
station_order <- w_station$station
w_diag$station_order <- match(w_diag$station, station_order)
w_diag <- w_diag[order(
  w_diag$station_order, w_diag$origin_id, w_diag$test_index
), ]
origin_lookup <- unique(w_diag[c("origin_id", "origin_index")])
origin_lookup <- origin_lookup[order(origin_lookup$origin_id), ]
n_station <- length(station_order)
n_origin <- nrow(origin_lookup)
n_horizon <- length(unique(w_diag$test_index[w_diag$origin_id == origin_lookup$origin_id[1L]]))
n_forecast <- n_origin * n_horizon
if (nrow(w_diag) != n_station * n_forecast) {
  stop("The BGNAR wind diagnostic grid is incomplete.")
}
matrix_by_station <- function(x) {
  matrix(x, nrow = n_station, ncol = n_forecast, byrow = TRUE)
}
observed_matrix <- matrix_by_station(w_diag$truth)
prediction_matrix <- matrix_by_station(w_diag$prediction)
coverage_matrix <- matrix_by_station(
  as.integer(w_diag$truth >= w_diag$lo & w_diag$truth <= w_diag$hi)
)
common_level_range <- range(observed_matrix, prediction_matrix, finite = TRUE)
grDevices::pdf(file.path(output_dir, "wind_forecast_diagnostics.pdf"),
               width = 12.6, height = 9.2, family = "Helvetica")
graphics::layout(
  matrix(c(1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L),
         nrow = 4L, ncol = 2L, byrow = TRUE),
  heights = c(1, 0.12, 1, 0.10)
)
old_par <- graphics::par(
  mar = c(4.8, 6.6, 4.1, 1.2),
  oma = c(0.4, 0.4, 0.6, 0.4), mgp = c(3.6, 0.80, 0),
  tcl = -0.25, cex.axis = 1.20, cex.lab = 1.35
)
station_ticks <- unique(round(seq(1, n_station, length.out = 10L)))
station_tick_at <- n_station - station_ticks + 1L
origin_ticks <- unique(c(1L, seq(1L, n_origin, by = 3L), n_origin))
origin_tick_at <- (origin_ticks - 1L) * n_horizon + (n_horizon + 1) / 2
origin_tick_labels <- origin_lookup$origin_index[origin_ticks]
block_boundaries <- seq(n_horizon + 0.5, n_forecast - 0.5, by = n_horizon)

draw_diagnostic_heatmap <- function(z, main, subtitle, colors, zlim,
                                    breaks = NULL) {
  graphics::par(
    mar = c(4.8, 6.6, 4.1, 1.2),
    mgp = c(3.6, 0.80, 0), tcl = -0.25,
    cex.axis = 1.20, cex.lab = 1.35
  )
  z_plot <- z[n_station:1L, , drop = FALSE]
  image_args <- list(
    x = seq_len(n_forecast), y = seq_len(n_station), z = t(z_plot),
    col = colors, xlim = c(0.5, n_forecast + 0.5),
    ylim = c(0.5, n_station + 0.5), axes = FALSE,
    xlab = "Forecast origin (five horizons per block)",
    ylab = "Station", useRaster = TRUE
  )
  if (is.null(breaks)) {
    image_args$zlim <- zlim
  } else {
    image_args$breaks <- breaks
  }
  do.call(graphics::image, image_args)
  graphics::axis(1, at = origin_tick_at, labels = origin_tick_labels)
  graphics::axis(
    2, at = station_tick_at, labels = station_order[station_ticks],
    las = 1, cex.axis = 1.10
  )
  graphics::abline(v = block_boundaries, col = "white", lwd = 0.75)
  graphics::box(col = "grey35")
  graphics::title(main = main, cex.main = 1.45, line = 1.45)
  if (nzchar(subtitle)) {
    graphics::mtext(subtitle, side = 3, line = 0.20, cex = 0.90,
                    col = "grey35")
  }
}

draw_continuous_legend <- function(colors, limits, label, digits = 2L) {
  graphics::par(mar = c(1.8, 6.6, 0.15, 1.2))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  bar_at <- c(0.17, 0.50, 0.83)
  graphics::rasterImage(
    as.raster(matrix(colors, nrow = 1L)),
    bar_at[1L], 0.40, bar_at[3L], 0.78,
    interpolate = TRUE
  )
  tick_values <- c(limits[1L], mean(limits), limits[2L])
  graphics::axis(
    1, at = bar_at, labels = formatC(tick_values, format = "f", digits = digits),
    pos = 0.40, tcl = -0.22, cex.axis = 1.00
  )
  if (nzchar(label)) {
    graphics::mtext(label, side = 1, line = 1.05, cex = 0.90)
  }
}

draw_coverage_legend <- function() {
  graphics::par(mar = c(1.2, 6.6, 0.15, 1.2))
  graphics::plot.new()
  graphics::legend(
    "center", c("Covered", "Missed"),
    fill = c("#56B4E9", "#D73027"), border = "grey45",
    horiz = TRUE, bty = "n", cex = 1.08
  )
}

level_colors <- grDevices::hcl.colors(101, "viridis")
draw_diagnostic_heatmap(
  observed_matrix, "(a) Observed wind-speed series", "",
  level_colors, common_level_range
)
draw_diagnostic_heatmap(
  prediction_matrix, "(b) BGNAR posterior predictive median", "",
  level_colors, common_level_range
)
draw_continuous_legend(
  level_colors, common_level_range, "", digits = 2L
)
draw_diagnostic_heatmap(
  coverage_matrix, "(c) BGNAR 95% predictive coverage",
  "",
  c("#D73027", "#56B4E9"), c(0, 1), breaks = c(-0.5, 0.5, 1.5)
)
draw_coverage_legend()
graphics::par(old_par)
graphics::layout(1)
grDevices::dev.off()

final_origin_id <- max(w_coef$origin_id)
alpha_last <- w_coef[
  w_coef$parameter == "alpha" & w_coef$origin_id == final_origin_id,
]
beta_all <- w_coef[w_coef$parameter == "beta", ]
beta_last <- beta_all[beta_all$origin_id == final_origin_id, ]
beta_last <- beta_last[order(beta_last$lag, beta_last$neighbor_order), ]
beta_last$coef_label <- sprintf("(%d,%d)", beta_last$lag,
                                beta_last$neighbor_order)
coef_levels <- unique(beta_last$coef_label)
beta_all$coef_label <- factor(
  sprintf("(%d,%d)", beta_all$lag, beta_all$neighbor_order),
  levels = coef_levels
)
origin_values <- sort(unique(beta_all$origin_index))
beta_all$origin_factor <- factor(beta_all$origin_index, levels = origin_values)
beta_matrix <- as.matrix(stats::xtabs(
  mean ~ origin_factor + coef_label, data = beta_all
))

lag_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
grDevices::pdf(file.path(output_dir, "wind_parameter_estimates.pdf"),
               width = 11.2, height = 4.8, family = "Helvetica")
old_par <- graphics::par(
  mfrow = c(1, 3), mar = c(5.4, 5.0, 3.2, 1.0),
  mgp = c(2.9, 0.90, 0), tcl = -0.30,
  cex.axis = 1.18, cex.lab = 1.25
)

graphics::boxplot(
  mean ~ lag, data = alpha_last, outline = FALSE,
  col = grDevices::adjustcolor(lag_colors, alpha.f = 0.28),
  border = lag_colors, xlab = "Temporal lag", ylab = "Posterior mean",
  main = expression(paste("(a) Posterior means of ", alpha[i * "," * j], " by station")),
  cex.main = 1.35
)
for (lag in sort(unique(alpha_last$lag))) {
  za <- alpha_last[alpha_last$lag == lag, ]
  za <- za[order(za$node), ]
  xj <- lag + seq(-0.17, 0.17, length.out = nrow(za))
  graphics::points(xj, za$mean, pch = 16, cex = 0.36,
                   col = grDevices::adjustcolor(lag_colors[[lag]], alpha.f = 0.60))
}
graphics::abline(h = 0, col = "grey55", lty = 2)
graphics::grid(col = "grey92")

x_beta <- seq_len(nrow(beta_last))
beta_ylim <- range(c(beta_last$q025, beta_last$q975), finite = TRUE)
graphics::plot(
  x_beta, beta_last$mean, type = "n", xaxt = "n", ylim = beta_ylim,
  xlab = "Coefficient (j,r)", ylab = "Estimate",
  main = expression(paste("(b) Posterior estimates of ", beta[j * "," * r], " at the final origin")),
  cex.main = 1.35, bty = "l"
)
graphics::abline(h = 0, col = "grey55", lty = 2)
graphics::segments(
  x_beta, beta_last$q025, x_beta, beta_last$q975,
  col = lag_colors[beta_last$lag], lwd = 2.3
)
graphics::points(x_beta, beta_last$mean, pch = 16,
                 col = lag_colors[beta_last$lag], cex = 0.95)
graphics::axis(1, at = x_beta, labels = beta_last$coef_label,
               las = 2, cex.axis = 0.96)
graphics::grid(col = "grey92")

zmax <- max(abs(beta_matrix), na.rm = TRUE)
heat_colors <- grDevices::colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)
graphics::image(
  x = origin_values, y = seq_along(coef_levels), z = beta_matrix,
  col = heat_colors, zlim = c(-zmax, zmax),
  xlab = "Forecast origin", ylab = "Coefficient (j,r)",
  yaxt = "n", bty = "l",
  main = expression(paste("(c) Posterior means of ", beta[j * "," * r], " by origin")),
  cex.main = 1.35
)
graphics::axis(2, at = seq_along(coef_levels), labels = coef_levels,
               las = 2, cex.axis = 0.96)
graphics::legend(
  "topright", c("Positive", "Near zero", "Negative"),
  fill = c("#B2182B", "#F7F7F7", "#2166AC"),
  border = "grey70", bty = "n", cex = 1.15
)

graphics::par(old_par)
grDevices::dev.off()

message("Manuscript figures written to ", output_dir)
