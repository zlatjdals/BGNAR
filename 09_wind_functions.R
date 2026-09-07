############################################################
## 09_wind_functions.R
## Wind-data network, rolling origins and compact summaries.
############################################################

read_wind_config <- function() {
  cfg <- read_cluster_config()
  wind_methods <- csv_tokens(Sys.getenv("WIND_METHODS", ""))
  if (length(wind_methods)) {
    if (!all(wind_methods %in% all_method_names())) {
      stop("Unknown WIND_METHODS entry.")
    }
    cfg$methods <- unique(wind_methods)
  }
  cfg$wind_n_stations <- env_int("WIND_N_STATIONS", 50L)
  cfg$wind_window <- env_int("WIND_WINDOW", 240L)
  cfg$wind_first_origin <- env_int("WIND_FIRST_ORIGIN", 541L)
  cfg$wind_origin_spacing <- env_int("WIND_ORIGIN_SPACING", 12L)
  cfg$wind_tail_buffer <- env_int("WIND_TAIL_BUFFER", 12L)
  cfg$wind_n_origins <- env_int("WIND_N_ORIGINS", 15L)
  cfg
}

load_wind_design <- function(cfg) {
  wind_env <- new.env(parent = emptyenv())
  data("vswind", package = "GNAR", envir = wind_env)
  X_full <- as.matrix(wind_env$vswindts)
  station_full <- make.unique(trimws(wind_env$vswindnames))
  n_full <- ncol(X_full)
  A_full <- matrix(0L, n_full, n_full)
  for (i in seq_len(n_full)) {
    nei <- as.integer(wind_env$vswindnet$edges[[i]])
    if (length(nei)) A_full[i, nei] <- 1L
  }
  A_full <- 1L * ((A_full + t(A_full)) > 0)
  diag(A_full) <- 0L
  g_full <- make_graph_from_adjacency(A_full, station_full)
  root <- match("CAPEL", station_full)
  if (is.na(root)) stop("CAPEL station not found in GNAR::vswind.")
  bfs <- as.integer(igraph::bfs(g_full, root = root, unreachable = FALSE)$order)
  bfs <- bfs[!is.na(bfs)]
  if (length(bfs) < cfg$wind_n_stations) {
    stop("Requested more connected wind stations than are available.")
  }
  selected <- bfs[seq_len(cfg$wind_n_stations)]
  station_names <- station_full[selected]
  X <- X_full[, selected, drop = FALSE]
  colnames(X) <- station_names
  A <- A_full[selected, selected, drop = FALSE]
  g <- make_graph_from_adjacency(A, station_names)
  W_list <- make_neighbor_matrices(g, max_r = max(cfg$s_fit))
  coords <- as.data.frame(wind_env$vswindcoords[selected, , drop = FALSE])
  names(coords) <- c("x", "y")
  station_info <- data.frame(
    station = station_names, coords, degree = as.numeric(igraph::degree(g))
  )
  last_origin <- nrow(X) - cfg$wind_tail_buffer
  origins <- seq.int(
    cfg$wind_first_origin, last_origin, by = cfg$wind_origin_spacing
  )
  origins <- head(origins, cfg$wind_n_origins)
  if (!length(origins) || any(origins + cfg$horizon > nrow(X))) {
    stop("Wind origin configuration exceeds the available data.")
  }
  list(
    X = X, g = g, W_list = W_list, adjacency = A,
    stations = station_info, origins = origins,
    source_rows = nrow(X_full)
  )
}

wind_window_data <- function(design, origin_id, cfg) {
  if (origin_id < 1L || origin_id > length(design$origins)) {
    stop("origin_id must be between 1 and ", length(design$origins), ".")
  }
  origin <- design$origins[[origin_id]]
  start <- origin - cfg$wind_window + 1L
  if (start < 1L) stop("WIND_WINDOW is too long for the first origin.")
  rows <- start:(origin + cfg$horizon)
  list(
    X = design$X[rows, , drop = FALSE],
    g = design$g, W_list = design$W_list,
    origin_id = origin_id, origin_index = origin,
    window_start = start, window_end = origin,
    test_index = origin + seq_len(cfg$horizon),
    station_names = design$stations$station
  )
}

wind_method_seed <- function(cfg, origin_id, method) {
  cfg$base_seed + 80000000L + 100L * as.integer(origin_id) +
    match(method, all_method_names())
}

compact_wind_coefficients <- function(fit_result, method, station_names) {
  ps <- summarize_bgnar_parameters(fit_result$fit)
  mu <- transform(
    ps$mu, method = method, parameter = "intercept",
    station = station_names[node], lag = NA_integer_,
    neighbor_order = NA_integer_
  )
  alpha <- transform(
    ps$alpha, method = method, parameter = "alpha",
    station = station_names[node], neighbor_order = NA_integer_
  )
  beta <- transform(
    ps$beta, method = method, parameter = "beta",
    node = NA_integer_, station = NA_character_
  )
  keep <- c("method", "parameter", "node", "station", "lag",
            "neighbor_order", "mean", "q025", "q500", "q975")
  rbind(mu[keep], alpha[keep], beta[keep])
}

## Posterior predictive distribution of the spatial mean.  The aggregation
## must be performed within each joint posterior predictive draw before taking
## quantiles.  Averaging nodewise marginal interval endpoints does not produce
## an interval for the spatial mean and fails to reflect uncertainty reduction
## from spatial aggregation.
compact_spatial_mean_predictive <- function(fit_result, method) {
  if (!method %in% all_bgnar_variants()) return(data.frame())
  pred_y <- fit_result$pred$pred_y
  if (is.null(pred_y) || length(dim(pred_y)) != 3L) {
    stop("BGNAR posterior predictive draws are unavailable.")
  }
  spatial_draws <- apply(pred_y, c(1L, 3L), mean)
  if (is.null(dim(spatial_draws))) {
    spatial_draws <- matrix(spatial_draws, nrow = dim(pred_y)[1L])
  }
  data.frame(
    horizon = seq_len(nrow(spatial_draws)),
    mean = rowMeans(spatial_draws),
    q025 = apply(spatial_draws, 1L, stats::quantile, probs = 0.025,
                 names = FALSE),
    q500 = apply(spatial_draws, 1L, stats::quantile, probs = 0.500,
                 names = FALSE),
    q975 = apply(spatial_draws, 1L, stats::quantile, probs = 0.975,
                 names = FALSE)
  )
}

compact_wind_method_result <- function(fit_result, dat, cfg, method, elapsed) {
  list(
    method = method,
    elapsed_seconds = as.numeric(elapsed),
    forecast = compact_forecast(
      fit_result, dat, cfg$wind_window, method, cfg
    ),
    coefficients = if (method %in% all_bgnar_variants()) {
      compact_wind_coefficients(fit_result, method, dat$station_names)
    } else data.frame(),
    hyperparameters = if (method %in% all_bgnar_variants()) {
      compact_bgnar_hyperparameters(fit_result, method)
    } else data.frame(),
    spatial_mean_predictive = compact_spatial_mean_predictive(
      fit_result, method
    ),
    selected_order = compact_selected_order(fit_result, method, cfg)
  )
}
