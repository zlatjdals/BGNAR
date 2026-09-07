############################################################
## 03_cluster_functions.R
## Reproducible cluster simulation functions for BGNAR.
############################################################

`%||%` <- function(x, y) if (is.null(x)) y else x

csv_tokens <- function(x) {
  z <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  z[nzchar(z)]
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!length(value) || is.na(value)) stop("Invalid integer environment variable: ", name)
  value
}

env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (!length(value) || !is.finite(value)) stop("Invalid numeric environment variable: ", name)
  value
}

env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, if (default) "1" else "0"))
  if (!value %in% c("0", "1", "false", "true", "no", "yes")) {
    stop("Invalid logical environment variable: ", name)
  }
  value %in% c("1", "true", "yes")
}

simulation_scenarios <- function() {
  list(
    M11 = list(
      alpha = c(0.7), beta = list(0.2),
      s = c(1L), local_sd = 0,
      label = "(1,[1]) case1"
    ),
    M12 = list(
      alpha = c(0.45), beta = list(0.45),
      s = c(1L), local_sd = 0,
      label = "(1,[1]) case2"
    ),
    M13 = list(
      alpha = c(0.2), beta = list(0.7),
      s = c(1L), local_sd = 0,
      label = "(1,[1]) case3"
    ),
    M2 = list(
      alpha = c(0.3, 0.2),
      beta = list(c(0.2, 0.1), c(0.1, 0)),
      s = c(2L, 2L), local_sd = 0,
      label = "(2,[2,2])"
    )
  )
}

all_bgnar_variants <- function() {
  c(
    "bgnar_nondummy",
    "bgnar_soc_fixed", "bgnar_soc_adaptive",
    "bgnar_dio_fixed", "bgnar_dio_adaptive",
    "bgnar_soc_dio_fixed", "bgnar_soc_dio_adaptive"
  )
}

all_method_names <- function() c(all_bgnar_variants(), "gnar", "bvar", "rw")

read_cluster_config <- function() {
  p_fit <- env_int("BGNAR_P", 5L)
  s_fit <- as.integer(csv_tokens(Sys.getenv("BGNAR_S", "3,3,3,3,3")))
  if (length(s_fit) != p_fit || anyNA(s_fit) || any(s_fit < 0L)) {
    stop("BGNAR_S must contain one nonnegative order for each BGNAR lag.")
  }
  methods <- csv_tokens(Sys.getenv(
    "METHODS", "bgnar_soc_dio_adaptive,gnar,bvar"
  ))
  if (!length(methods) || !all(methods %in% all_method_names())) {
    stop("Unknown METHODS entry. Valid names: ", paste(all_method_names(), collapse = ", "))
  }
  gnar_mode <- tolower(Sys.getenv("GNAR_MODE", "centered"))
  if (!gnar_mode %in% c("centered", "raw")) {
    stop("GNAR_MODE must be 'centered' or 'raw'.")
  }
  intercept_mode <- tolower(Sys.getenv("INTERCEPT_MODE", "common"))
  if (!intercept_mode %in% c("none", "common", "nodewise")) {
    stop("INTERCEPT_MODE must be none, common, or nodewise.")
  }
  fixed_tau_scale <- tolower(Sys.getenv("FIXED_TAU_SCALE", "absolute"))
  if (!fixed_tau_scale %in% c("absolute", "pooled")) {
    stop("FIXED_TAU_SCALE must be absolute or pooled.")
  }
  list(
    N = env_int("N_NODES", 20L),
    horizon = env_int("FORECAST_HORIZON", 5L),
    p_fit = p_fit,
    s_fit = s_fit,
    dgp_burn = env_int("DGP_BURN", 300L),
    innovation_sd = env_num("INNOVATION_SD", sqrt(0.1)),
    intercept_mode = intercept_mode,
    intercept_value = env_num("INTERCEPT_VALUE", 0.1),
    intercept_spread = env_num("INTERCEPT_SPREAD", 0),
    methods = unique(methods),
    gnar_mode = gnar_mode,
    gnar_p_max = env_int("GNAR_P_MAX", 5L),
    gnar_s_max = env_int("GNAR_S_MAX", 3L),
    ## Despite the legacy name, BVAR_P_MAX is the fixed fitted BVAR lag order.
    bvar_p_max = env_int("BVAR_P_MAX", 5L),
    mcmc_iter = env_int("MCMC_ITER", 2000L),
    mcmc_burn = env_int("MCMC_BURN", 1000L),
    mcmc_thin = env_int("MCMC_THIN", 2L),
    bvar_draws = env_int("BVAR_DRAWS", 1000L),
    bvar_burn = env_int("BVAR_BURN", 500L),
    reduce_soc = env_flag("REDUCE_SOC", TRUE),
    fixed_tau2 = env_num("FIXED_TAU2", 1),
    fixed_tau_scale = fixed_tau_scale,
    adaptive_a = env_num("ADAPTIVE_TAU_A", 2.5),
    adaptive_b = env_num("ADAPTIVE_TAU_B", 1),
    save_full_fits = env_flag("SAVE_FULL_FITS", FALSE),
    overwrite = env_flag("OVERWRITE", FALSE),
    base_seed = env_int("BASE_SEED", 260804L)
  )
}

make_structure_graph_cluster <- function(structure, N, seed) {
  structure <- toupper(structure)
  if (!structure %in% c("ER", "SBM", "SWN")) stop("Unknown structure: ", structure)
  if (structure == "SBM" && N %% 4L != 0L) {
    stop("The four-block SBM requires N_NODES divisible by four.")
  }
  for (attempt in seq_len(500L)) {
    set.seed(seed + attempt - 1L)
    g <- switch(
      structure,
      ER = igraph::sample_gnp(N, p = 0.2, directed = FALSE, loops = FALSE),
      SBM = {
        K <- 4L
        pref_mat <- matrix(0.075, nrow = K, ncol = K)
        diag(pref_mat) <- 0.70
        
        igraph::sample_sbm(
          N,
          pref.matrix = pref_mat,
          block.sizes = rep(N / 4L, 4L),
          directed = FALSE, loops = FALSE
        )
      },
      SWN = igraph::sample_smallworld(
        dim = 1L,
        size = N,
        nei = 2,
        p = 0.1,
        loops = FALSE,
        multiple = FALSE
      )
    )
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    if (igraph::is_connected(g) && igraph::diameter(g) >= 3L) {
      igraph::V(g)$name <- paste0("Node", seq_len(N))
      return(g)
    }
  }
  stop("Could not generate a connected ", structure, " graph.")
}

make_intercept_vector <- function(N, mode = "common", value = 0.05, spread = 0.02) {
  switch(
    mode,
    none = rep(0, N),
    common = rep(value, N),
    nodewise = value + seq(-spread, spread, length.out = N),
    stop("Unknown intercept mode: ", mode)
  )
}

pad_dgp_truth <- function(alpha_true, beta_true, s_true, N, p_fit, s_fit) {
  alpha_pad <- matrix(0, N, p_fit)
  alpha_pad[, seq_len(ncol(alpha_true))] <- alpha_true
  beta_pad <- matrix(0, p_fit, max(s_fit))
  for (j in seq_along(beta_true)) {
    if (s_true[j] > 0L) {
      beta_pad[j, seq_len(s_true[j])] <- beta_true[[j]][seq_len(s_true[j])]
    }
  }
  list(alpha = alpha_pad, beta = beta_pad)
}

simulate_cluster_dgp <- function(scenario, structure, train_end, replication, cfg) {
  specs <- simulation_scenarios()
  if (!scenario %in% names(specs)) stop("Unknown DGP: ", scenario)
  spec <- specs[[scenario]]
  structure <- toupper(structure)
  scenario_id <- match(scenario, names(specs))
  structure_id <- match(structure, c("ER", "SBM", "SWN"))
  graph_seed <- cfg$base_seed + 10000L * structure_id + as.integer(replication)
  data_seed <- cfg$base_seed + 1000000L * scenario_id +
    10000L * structure_id + as.integer(replication)
  
  g <- make_structure_graph_cluster(structure, cfg$N, graph_seed)
  W_list <- make_neighbor_matrices(g, max_r = max(cfg$s_fit))
  alpha_true <- matrix(rep(spec$alpha, each = cfg$N), nrow = cfg$N)
  
  if (spec$local_sd > 0) {
    accepted <- FALSE
    for (attempt in seq_len(2000L)) {
      set.seed(data_seed + 70000L + attempt)
      u <- stats::rnorm(cfg$N, 0, spec$local_sd)
      u <- u - mean(u)
      trial <- alpha_true
      trial[, 1L] <- trial[, 1L] + u
      trial[, 2L] <- trial[, 2L] - u
      if (gnar_spectral_radius(trial, spec$beta, W_list, spec$s) < 1) {
        alpha_true <- trial
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("Could not draw a stationary local-alpha DGP.")
  }
  
  rho <- gnar_spectral_radius(alpha_true, spec$beta, W_list, spec$s)
  if (!is.finite(rho) || rho >= 1) {
    stop("Nonstationary fixed DGP for ", scenario, "/", structure, ": rho=", rho)
  }
  
  intercept <- make_intercept_vector(
    cfg$N, cfg$intercept_mode, cfg$intercept_value, cfg$intercept_spread
  )
  p_true <- ncol(alpha_true)
  total_keep <- train_end + cfg$horizon
  total_sim <- cfg$dgp_burn + total_keep
  set.seed(data_seed)
  X_full <- matrix(0, nrow = p_true + total_sim, ncol = cfg$N)
  X_full[seq_len(p_true), ] <- matrix(
    stats::rnorm(p_true * cfg$N, 0, cfg$innovation_sd),
    nrow = p_true, ncol = cfg$N
  )
  for (tt in (p_true + 1L):nrow(X_full)) {
    mean_tt <- intercept
    for (j in seq_len(p_true)) {
      xlag <- X_full[tt - j, ]
      mean_tt <- mean_tt + alpha_true[, j] * xlag
      if (spec$s[j] > 0L) {
        for (r in seq_len(spec$s[j])) {
          mean_tt <- mean_tt + spec$beta[[j]][r] *
            as.vector(W_list[[r]] %*% xlag)
        }
      }
    }
    X_full[tt, ] <- mean_tt + stats::rnorm(cfg$N, 0, cfg$innovation_sd)
  }
  first <- p_true + cfg$dgp_burn + 1L
  X <- X_full[first:(first + total_keep - 1L), , drop = FALSE]
  colnames(X) <- paste0("Node", seq_len(cfg$N))
  truth <- pad_dgp_truth(
    alpha_true, spec$beta, spec$s, cfg$N, cfg$p_fit, cfg$s_fit
  )
  A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  list(
    X = X, g = g, W_list = W_list, adjacency = A,
    alpha_true = truth$alpha, beta_true = truth$beta,
    intercept_true = intercept, rho = rho,
    graph_seed = graph_seed, data_seed = data_seed,
    graph_density = igraph::edge_density(g),
    mean_degree = mean(igraph::degree(g)),
    dgp_label = spec$label
  )
}

## Public selector for the requested BGNAR dummy/tau combinations.
## tau_mode is ignored when dummy="none".
make_bgnar_spec <- function(dummy = c("none", "soc", "dio", "soc_dio"),
                            tau_mode = c("fixed", "adaptive"),
                            X_train = NULL, p = 3L, W_list = NULL,
                            s_vec = c(3L, 3L, 3L), alpha_mode = "local",
                            reduce_soc = TRUE, fixed_tau2 = 1,
                            fixed_tau_scale = c("absolute", "pooled"),
                            adaptive_a = 2.5, adaptive_b = 1) {
  dummy <- match.arg(dummy)
  tau_mode <- match.arg(tau_mode)
  fixed_tau_scale <- match.arg(fixed_tau_scale)
  use_soc <- dummy %in% c("soc", "soc_dio")
  use_dio <- dummy %in% c("dio", "soc_dio")
  if (dummy == "none") {
    return(list(
      use_dummy = FALSE, use_soc = FALSE, use_dio = FALSE,
      learn_dummy_tightness = FALSE, tau_soc = 1, tau_dio = 1
    ))
  }
  if (is.null(X_train) || is.null(W_list)) {
    stop("X_train and W_list are required for fixed/adaptive tau calibration.")
  }
  raw_scale <- stats::sd(as.vector(X_train))
  if (!is.finite(raw_scale) || raw_scale <= 0) raw_scale <- 1
  base <- list(
    use_dummy = TRUE, use_soc = use_soc, use_dio = use_dio,
    reduce_soc_dummy = isTRUE(reduce_soc),
    dummy_error_covariance = "identity"
  )
  if (tau_mode == "fixed") {
    tau2 <- fixed_tau2 * if (fixed_tau_scale == "pooled") raw_scale^2 else 1
    return(c(base, list(
      learn_dummy_tightness = FALSE,
      tau_soc = sqrt(tau2), tau_dio = sqrt(tau2)
    )))
  }
  
  dummy_dim <- make_dummy_observations(
    X_train, p, W_list, s_vec, alpha_mode = alpha_mode,
    reduce_soc = isTRUE(reduce_soc)
  )
  q <- length(dummy_dim$x_soc)
  N <- ncol(X_train)
  a_raw <- adaptive_a * raw_scale^2
  shape_soc <- 1 + q / 2 * (1 / adaptive_b - 1)
  shape_dio <- 1 + N / 2 * (1 / adaptive_b - 1)
  scale_soc <- a_raw * q / (2 * adaptive_b)
  scale_dio <- a_raw * N / (2 * adaptive_b)
  if (any(!is.finite(c(shape_soc, shape_dio, scale_soc, scale_dio))) ||
      min(shape_soc, shape_dio, scale_soc, scale_dio) <= 0) {
    stop("Invalid dimension/scale-calibrated adaptive tau prior.")
  }
  c(base, list(
    learn_dummy_tightness = TRUE,
    tau_soc = sqrt(a_raw), tau_dio = sqrt(a_raw),
    a_soc = shape_soc, b_soc = scale_soc,
    a_dio = shape_dio, b_dio = scale_dio
  ))
}

parse_bgnar_variant <- function(method) {
  switch(
    method,
    bgnar_nondummy = list(dummy = "none", tau_mode = "fixed"),
    bgnar_soc_fixed = list(dummy = "soc", tau_mode = "fixed"),
    bgnar_soc_adaptive = list(dummy = "soc", tau_mode = "adaptive"),
    bgnar_dio_fixed = list(dummy = "dio", tau_mode = "fixed"),
    bgnar_dio_adaptive = list(dummy = "dio", tau_mode = "adaptive"),
    bgnar_soc_dio_fixed = list(dummy = "soc_dio", tau_mode = "fixed"),
    bgnar_soc_dio_adaptive = list(dummy = "soc_dio", tau_mode = "adaptive"),
    stop("Unknown BGNAR variant: ", method)
  )
}

fit_bgnar_variant <- function(dat, train_end, method, cfg, seed) {
  variant <- parse_bgnar_variant(method)
  X_train <- dat$X[seq_len(train_end), , drop = FALSE]
  raw_scale <- stats::sd(as.vector(X_train))
  if (!is.finite(raw_scale) || raw_scale <= 0) raw_scale <- 1
  dummy_args <- make_bgnar_spec(
    dummy = variant$dummy, tau_mode = variant$tau_mode,
    X_train = X_train, p = cfg$p_fit, W_list = dat$W_list,
    s_vec = cfg$s_fit, alpha_mode = "local",
    reduce_soc = cfg$reduce_soc, fixed_tau2 = cfg$fixed_tau2,
    fixed_tau_scale = cfg$fixed_tau_scale,
    adaptive_a = cfg$adaptive_a, adaptive_b = cfg$adaptive_b
  )
  base_args <- list(
    X_scaled = dat$X, train_end = train_end,
    p = cfg$p_fit, W_list = dat$W_list, s_vec = cfg$s_fit,
    horizon = cfg$horizon,
    n_iter = cfg$mcmc_iter, burn = cfg$mcmc_burn, thin = cfg$mcmc_thin,
    b = 0, b_mode = "fixed", alpha_mode = "local",
    c_kappa = 3, d_kappa = 2,
    A_grid = c(1, 1.5, 2, 2.5, 3),
    D_grid = c(1, 1.5, 2, 2.5, 3),
    a_alpha = 3, b_alpha = 0.16,
    a_beta = 3, b_beta = 0.16,
    c_sigma = 2, d_sigma = raw_scale^2,
    c_mu = 1e7 * raw_scale^2,
    simulate = TRUE, verbose = FALSE
  )
  set.seed(seed)
  do.call(fit_predict_bgnar, c(base_args, dummy_args))
}

## Public GNAR selector. Centered GNAR is the production default.
fit_gnar_by_mode <- function(dat, train_end, cfg, mode = cfg$gnar_mode) {
  mode <- match.arg(mode, c("centered", "raw"))
  if (mode == "centered") {
    fit_predict_gnar_bic_centered(
      X_scaled = dat$X, train_end = train_end, g = dat$g,
      horizon = cfg$horizon, p_max = cfg$gnar_p_max,
      s_max = cfg$gnar_s_max, globalalpha_grid = c(TRUE, FALSE)
    )
  } else {
    fit_predict_gnar_bic(
      X_scaled = dat$X, train_end = train_end, g = dat$g,
      horizon = cfg$horizon, p_max = cfg$gnar_p_max,
      s_max = cfg$gnar_s_max, globalalpha_grid = c(TRUE, FALSE)
    )
  }
}

fit_requested_method <- function(method, dat, train_end, cfg, seed) {
  if (method %in% all_bgnar_variants()) {
    return(fit_bgnar_variant(dat, train_end, method, cfg, seed))
  }
  if (method == "gnar") return(fit_gnar_by_mode(dat, train_end, cfg))
  if (method == "bvar") {
    set.seed(seed)
    return(fit_predict_bvar(
      X_scaled = dat$X, train_end = train_end, horizon = cfg$horizon,
      p = cfg$bvar_p_max, draws = cfg$bvar_draws,
      burn = cfg$bvar_burn, prior_mean = 1, use_dummy = FALSE
    ))
  }
  if (method == "rw") {
    return(fit_predict_random_walk(dat$X, train_end, cfg$horizon))
  }
  stop("Unknown method: ", method)
}

forecast_component <- function(fit_result) {
  if (!is.null(fit_result$forecast$pred_med)) return(fit_result$forecast)
  if (!is.null(fit_result$pred$pred_med)) return(fit_result$pred)
  if (!is.null(fit_result$summary$pred_med)) return(fit_result$summary)
  stop("Could not locate a summarized forecast matrix in the fitted object.")
}

compact_forecast <- function(fit_result, dat, train_end, method, cfg) {
  fc <- forecast_component(fit_result)
  pred <- as.matrix(fc$pred_med)
  truth <- dat$X[(train_end + 1L):(train_end + cfg$horizon), , drop = FALSE]
  if (!all(dim(pred) == dim(truth))) stop("Forecast dimension mismatch for ", method)
  as_interval <- function(x) {
    if (is.null(x)) matrix(NA_real_, nrow(truth), ncol(truth)) else as.matrix(x)
  }
  lo <- as_interval(fc$pred_lo)
  hi <- as_interval(fc$pred_hi)
  overall <- evaluate_forecast(truth, pred, fc$pred_lo, fc$pred_hi)
  horizon_metrics <- do.call(rbind, lapply(seq_len(cfg$horizon), function(h) {
    z <- evaluate_forecast(
      truth[h, , drop = FALSE], pred[h, , drop = FALSE],
      if (all(is.finite(lo[h, ]))) lo[h, , drop = FALSE] else NULL,
      if (all(is.finite(hi[h, ]))) hi[h, , drop = FALSE] else NULL
    )
    data.frame(horizon = h, rmse = z$rmse, mae = z$mae,
               coverage = z$coverage, avg_width = z$avg_width)
  }))
  point <- data.frame(
    horizon = rep(seq_len(cfg$horizon), each = cfg$N),
    node = rep(seq_len(cfg$N), times = cfg$horizon),
    truth = as.vector(t(truth)), prediction = as.vector(t(pred)),
    lo = as.vector(t(lo)), hi = as.vector(t(hi))
  )
  point$error <- point$prediction - point$truth
  list(overall = overall, horizon = horizon_metrics, point = point)
}

## Compare in-sample one-step fitted values on the same time range for every
## method.  With the production settings this scores t = 4, ..., T because the
## largest candidate/fitted temporal order is three.  Centered GNAR fitted
## values are returned to the original scale before scoring.
compact_fit_metrics <- function(fit_result, dat, train_end, method, cfg) {
  common_p <- max(cfg$p_fit, cfg$gnar_p_max, cfg$bvar_p_max)
  if (train_end <= common_p) {
    return(data.frame(
      fit_rmse = NA_real_, fit_mae = NA_real_,
      fit_start = common_p + 1L, fit_n_time = 0L
    ))
  }
  
  truth <- dat$X[(common_p + 1L):train_end, , drop = FALSE]
  
  if (method == "rw") {
    pred <- dat$X[common_p:(train_end - 1L), , drop = FALSE]
  } else {
    fs <- fit_result$fitted_summary
    if (is.null(fs) || is.null(fs$fit_med)) {
      return(data.frame(
        fit_rmse = NA_real_, fit_mae = NA_real_,
        fit_start = common_p + 1L, fit_n_time = nrow(truth)
      ))
    }
    p_method <- if (method %in% all_bgnar_variants()) {
      cfg$p_fit
    } else {
      as.integer(fit_result$selected_order$p[[1L]])
    }
    pred_all <- as.matrix(fs$fit_med)
    first_row <- common_p - p_method + 1L
    if (first_row < 1L || first_row > nrow(pred_all)) {
      stop("Could not align fitted values for ", method, ".")
    }
    pred <- pred_all[first_row:nrow(pred_all), , drop = FALSE]
    if (method == "gnar" && cfg$gnar_mode == "centered") {
      pred <- sweep(pred, 2L, fit_result$training_node_means, FUN = "+")
    }
  }
  
  if (!all(dim(pred) == dim(truth))) {
    stop(
      "Common-start fitted-value dimension mismatch for ", method,
      ": truth=", paste(dim(truth), collapse = "x"),
      ", fitted=", paste(dim(pred), collapse = "x"), "."
    )
  }
  metric <- evaluate_forecast(truth, pred)
  data.frame(
    fit_rmse = metric$rmse, fit_mae = metric$mae,
    fit_start = common_p + 1L, fit_n_time = nrow(truth)
  )
}

compact_bgnar_coefficients <- function(fit_result, dat, method) {
  ps <- summarize_bgnar_parameters(fit_result$fit)
  mu <- transform(
    ps$mu, method = method, parameter = "intercept", lag = NA_integer_,
    neighbor_order = NA_integer_, truth = dat$intercept_true[node]
  )
  alpha <- transform(
    ps$alpha, method = method, parameter = "alpha",
    neighbor_order = NA_integer_, truth = dat$alpha_true[cbind(node, lag)]
  )
  beta <- transform(
    ps$beta, method = method, parameter = "beta", node = NA_integer_,
    truth = dat$beta_true[cbind(lag, neighbor_order)]
  )
  keep <- c("method", "parameter", "node", "lag", "neighbor_order",
            "truth", "mean", "q025", "q500", "q975")
  out <- rbind(mu[keep], alpha[keep], beta[keep])
  out$covered <- out$truth >= out$q025 & out$truth <= out$q975
  out$error <- out$mean - out$truth
  out
}

## Convert a BIC-selected GNAR coefficient vector to the same fixed envelope
## used for BGNAR.  Omitted coefficients are scored as zero so order-selection
## errors are included in parameter recovery.  A selected global alpha is
## repeated over nodes; a selected local alpha retains its node index.
compact_gnar_coefficients <- function(fit_result, dat, method, cfg) {
  coef_hat <- stats::coef(fit_result$fit)
  coef_hat <- as.numeric(coef_hat)
  names(coef_hat) <- names(stats::coef(fit_result$fit))
  
  alpha <- expand.grid(
    node = seq_len(cfg$N), lag = seq_len(cfg$p_fit),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  alpha$parameter <- "alpha"
  alpha$neighbor_order <- NA_integer_
  alpha$truth <- dat$alpha_true[cbind(alpha$node, alpha$lag)]
  alpha$mean <- 0
  
  beta_rows <- lapply(seq_len(cfg$p_fit), function(j) {
    if (cfg$s_fit[j] < 1L) return(NULL)
    data.frame(
      node = NA_integer_, lag = j,
      neighbor_order = seq_len(cfg$s_fit[j]),
      parameter = "beta", stringsAsFactors = FALSE
    )
  })
  beta <- do.call(rbind, beta_rows[!vapply(beta_rows, is.null, logical(1L))])
  beta$truth <- dat$beta_true[cbind(beta$lag, beta$neighbor_order)]
  beta$mean <- 0
  
  for (nm in names(coef_hat)) {
    value <- unname(coef_hat[[nm]])
    global_match <- regmatches(nm, regexec("^dmatalpha([0-9]+)$", nm))[[1L]]
    if (length(global_match)) {
      j <- as.integer(global_match[[2L]])
      if (j <= cfg$p_fit) alpha$mean[alpha$lag == j] <- value
      next
    }
    local_match <- regmatches(
      nm, regexec("^dmatalpha([0-9]+)node([0-9]+)$", nm)
    )[[1L]]
    if (length(local_match)) {
      j <- as.integer(local_match[[2L]])
      i <- as.integer(local_match[[3L]])
      if (j <= cfg$p_fit && i <= cfg$N) {
        alpha$mean[alpha$lag == j & alpha$node == i] <- value
      }
      next
    }
    beta_match <- regmatches(
      nm, regexec("^dmatbeta([0-9]+)\\.([0-9]+)$", nm)
    )[[1L]]
    if (length(beta_match)) {
      j <- as.integer(beta_match[[2L]])
      r <- as.integer(beta_match[[3L]])
      hit <- beta$lag == j & beta$neighbor_order == r
      if (any(hit)) beta$mean[hit] <- value
    }
  }
  
  make_output <- function(z) {
    z$method <- method
    z$q025 <- z$q500 <- z$q975 <- NA_real_
    z$covered <- NA
    z$error <- z$mean - z$truth
    z[c("method", "parameter", "node", "lag", "neighbor_order",
        "truth", "mean", "q025", "q500", "q975", "covered", "error")]
  }
  rbind(make_output(alpha), make_output(beta))
}

summarize_draw_vector <- function(x, name) {
  if (is.null(x) || !length(x) || all(!is.finite(x))) return(NULL)
  x <- x[is.finite(x)]
  data.frame(
    parameter = name, mean = mean(x), sd = stats::sd(x),
    q025 = stats::quantile(x, 0.025), q500 = stats::median(x),
    q975 = stats::quantile(x, 0.975)
  )
}

compact_bgnar_hyperparameters <- function(fit_result, method) {
  fit <- fit_result$fit
  rows <- list(
    summarize_draw_vector(fit$tau_soc2_draws, "tau_soc2"),
    summarize_draw_vector(fit$tau_dio2_draws, "tau_dio2"),
    summarize_draw_vector(fit$kappa_alpha_draws, "kappa_alpha"),
    summarize_draw_vector(fit$lambda_alpha_draws, "lambda_alpha2"),
    summarize_draw_vector(fit$lambda_beta_draws, "lambda_beta2")
  )
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(data.frame())
  transform(do.call(rbind, rows), method = method)
}

compact_selected_order <- function(fit_result, method, cfg) {
  if (method %in% all_bgnar_variants()) {
    return(data.frame(
      method = method, p = cfg$p_fit,
      s_vec = paste(cfg$s_fit, collapse = ";"),
      globalalpha = FALSE, criterion = "fixed large envelope"
    ))
  }
  if (method == "rw") {
    return(data.frame(method = method, p = NA_integer_, s_vec = NA_character_,
                      globalalpha = NA, criterion = "random walk"))
  }
  z <- fit_result$selected_order
  data.frame(
    method = method,
    p = as.integer(z$p[[1L]]),
    s_vec = if ("s_vec" %in% names(z)) as.character(z$s_vec[[1L]]) else NA_character_,
    globalalpha = if ("globalalpha" %in% names(z)) as.logical(z$globalalpha[[1L]]) else NA,
    criterion = as.character(z$criterion[[1L]] %||% "BIC")
  )
}

compact_method_result <- function(fit_result, dat, train_end, method, cfg, elapsed) {
  list(
    method = method,
    elapsed_seconds = as.numeric(elapsed),
    fit = compact_fit_metrics(fit_result, dat, train_end, method, cfg),
    forecast = compact_forecast(fit_result, dat, train_end, method, cfg),
    coefficients = if (method %in% all_bgnar_variants()) {
      compact_bgnar_coefficients(fit_result, dat, method)
    } else if (method == "gnar") {
      compact_gnar_coefficients(fit_result, dat, method, cfg)
    } else data.frame(),
    hyperparameters = if (method %in% all_bgnar_variants()) {
      compact_bgnar_hyperparameters(fit_result, method)
    } else data.frame(),
    selected_order = compact_selected_order(fit_result, method, cfg)
  )
}

method_seed <- function(cfg, scenario, structure, train_end, replication, method) {
  cfg$base_seed + 50000000L + 1000000L * match(scenario, names(simulation_scenarios())) +
    100000L * match(toupper(structure), c("ER", "SBM", "SWN")) +
    1000L * as.integer(train_end) + 10L * as.integer(replication) +
    match(method, all_method_names())
}
