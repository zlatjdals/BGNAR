############################################################
## 00_utils.R
## Utilities for Bayesian GNAR simulation and evaluation
############################################################

## Keep package loading local to the functions that need it.  This makes the
## numerical code usable in batch jobs even when optional plotting packages
## (notably ggraph) are unavailable.

check_packages <- function(pkgs, install_missing = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    msg <- paste0("Missing packages: ", paste(missing, collapse = ", "), 
                  ". Install with install.packages(c(", 
                  paste(sprintf('"%s"', missing), collapse = ", "), ")).")
    if (install_missing) {
      install.packages(missing)
    } else {
      warning(msg, call. = FALSE)
    }
  }
  invisible(missing)
}

rinvgamma <- function(n, shape, scale) {
  1 / stats::rgamma(n, shape = shape, rate = scale)
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

sample_categorical_logprob <- function(logw) {
  prob <- exp(logw - log_sum_exp(logw))
  sample(seq_along(prob), size = 1, prob = prob)
}

gamma_shape_rate_from_mode_sd <- function(mode, sd) {
  if (mode <= 0 || sd <= 0) stop("mode and sd must be positive.")
  scale <- (sqrt(mode^2 + 4 * sd^2) - mode) / 2
  shape <- 1 + mode / scale
  c(shape = shape, rate = 1 / scale)
}

gamma_grid_prob <- function(grid, mode, sd) {
  pars <- gamma_shape_rate_from_mode_sd(mode, sd)
  w <- stats::dgamma(grid, shape = pars[["shape"]], rate = pars[["rate"]])
  if (!all(is.finite(w)) || sum(w) <= 0) stop("Invalid discretized Gamma weights.")
  w / sum(w)
}

row_normalize <- function(A) {
  A <- as.matrix(A)
  rs <- rowSums(A)
  W <- A
  for (i in seq_len(nrow(A))) {
    if (is.finite(rs[i]) && rs[i] > 0) W[i, ] <- A[i, ] / rs[i] else W[i, ] <- 0
  }
  W
}

make_connected_er_graph <- function(N, prob = 0.25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  repeat {
    g <- igraph::sample_gnp(N, prob, directed = FALSE, loops = FALSE)
    if (igraph::is_connected(g)) break
  }
  igraph::V(g)$name <- as.character(seq_len(N))
  g
}

make_neighbor_matrices <- function(g, max_r = 2) {
  N <- igraph::vcount(g)
  dist_mat <- igraph::distances(g)
  W_list <- vector("list", max_r)
  for (r in seq_len(max_r)) {
    A_r <- matrix(0, N, N)
    A_r[dist_mat == r] <- 1
    diag(A_r) <- 0
    W_list[[r]] <- row_normalize(A_r)
  }
  W_list
}

plot_network <- function(g, title = "Simulated network") {
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_network() requires the ggraph and ggplot2 packages.")
  }
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.5) +
    ggraph::geom_node_point(size = 4) +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 4) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(title)
}

make_graph_from_adjacency <- function(A, node_names = NULL) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The igraph package is required.")
  }
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  A <- 1L * ((A + t(A)) > 0)
  diag(A) <- 0L
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  if (is.null(node_names)) node_names <- as.character(seq_len(nrow(A)))
  igraph::V(g)$name <- node_names
  g
}

standardize_ts <- function(X_train, X_all = NULL) {
  X_train <- as.matrix(X_train)
  center <- colMeans(X_train)
  scalev <- apply(X_train, 2, sd)
  scalev[scalev == 0 | !is.finite(scalev)] <- 1
  if (is.null(X_all)) X_all <- X_train
  X_scaled <- sweep(sweep(as.matrix(X_all), 2, center, "-"), 2, scalev, "/")
  list(X = X_scaled, center = center, scale = scalev)
}

standardize_network_ts <- function(X_train, X_all = NULL) {
  ## Use one common affine transformation for all nodes.  This is preferable
  ## when nodes share physical units because it preserves common GNAR beta
  ## coefficients; node-wise scaling generally does not.
  X_train <- as.matrix(X_train)
  center <- mean(X_train)
  scalev <- stats::sd(as.vector(X_train))
  if (!is.finite(scalev) || scalev == 0) scalev <- 1
  if (is.null(X_all)) X_all <- X_train
  list(X = (as.matrix(X_all) - center) / scalev,
       center = center, scale = scalev)
}

center_scale_network_ts <- function(X_train, X_all = NULL) {
  ## Match the uploaded d analysis: remove a separate training-sample mean
  ## from each node, then divide every node by one common scale.  The common
  ## scale preserves the interpretation of shared network coefficients.
  X_train <- as.matrix(X_train)
  center <- colMeans(X_train)
  centered_train <- sweep(X_train, 2, center, "-")
  common_scale <- stats::sd(as.vector(centered_train))
  if (!is.finite(common_scale) || common_scale == 0) common_scale <- 1
  if (is.null(X_all)) X_all <- X_train
  list(
    X = sweep(sweep(as.matrix(X_all), 2, center, "-"), 2, common_scale, "/"),
    center = center,
    scale = rep(common_scale, ncol(X_train))
  )
}

scale_network_levels <- function(X_train, X_all = NULL) {
  ## Rescale all nodes by one common training-sample standard deviation while
  ## preserving the level.  SOC/DIO dummy observations are level restrictions,
  ## so subtracting the training mean would change their substantive target.
  X_train <- as.matrix(X_train)
  scalev <- stats::sd(as.vector(X_train))
  if (!is.finite(scalev) || scalev == 0) scalev <- 1
  if (is.null(X_all)) X_all <- X_train
  list(X = as.matrix(X_all) / scalev, center = 0, scale = scalev)
}

compute_metrics <- function(pred, truth, lo = NULL, hi = NULL) {
  pred <- as.matrix(pred)
  truth <- as.matrix(truth)
  err <- pred - truth
  out <- list(
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    mae = mean(abs(err), na.rm = TRUE)
  )
  if (!is.null(lo) && !is.null(hi)) {
    out$coverage <- mean(truth >= lo & truth <= hi, na.rm = TRUE)
    out$avg_width <- mean(hi - lo, na.rm = TRUE)
  }
  out
}

plot_prediction_node <- function(pred_med, pred_lo, pred_hi, truth, node = 1,
                                 title_prefix = "Posterior predictive interval") {
  df <- data.frame(
    h = seq_len(nrow(pred_med)),
    truth = truth[seq_len(nrow(pred_med)), node],
    median = pred_med[, node],
    lo = pred_lo[, node],
    hi = pred_hi[, node]
  )
  ggplot2::ggplot(df, ggplot2::aes(x = h)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = truth), linewidth = 0.9) +
    ggplot2::geom_line(ggplot2::aes(y = median), linetype = "dashed", linewidth = 0.9) +
    ggplot2::labs(
      title = paste0(title_prefix, ": node ", node),
      x = "Forecast horizon",
      y = "Standardized value"
    ) +
    ggplot2::theme_minimal()
}

extract_forecast_mats_from_result <- function(res_m, horizon, N) {
  ## Return list(pred_med, pred_lo, pred_hi)
  ## Handles BGNAR, BVAR, GNAR-like objects with different storage structures.
  
  pred_med <- NULL
  pred_lo <- NULL
  pred_hi <- NULL
  
  ## ------------------------------------------------------------
  ## 1. Try summary slot
  ## ------------------------------------------------------------
  if (!is.null(res_m$summary)) {
    if (!is.null(res_m$summary$pred_med)) pred_med <- res_m$summary$pred_med
    if (!is.null(res_m$summary$pred_lo))  pred_lo  <- res_m$summary$pred_lo
    if (!is.null(res_m$summary$pred_hi))  pred_hi  <- res_m$summary$pred_hi
  }
  
  ## ------------------------------------------------------------
  ## 2. Try forecast slot
  ## ------------------------------------------------------------
  if (is.null(pred_med) && !is.null(res_m$forecast)) {
    if (!is.null(res_m$forecast$pred_med)) pred_med <- res_m$forecast$pred_med
    if (!is.null(res_m$forecast$pred_lo))  pred_lo  <- res_m$forecast$pred_lo
    if (!is.null(res_m$forecast$pred_hi))  pred_hi  <- res_m$forecast$pred_hi
  } else if (!is.null(res_m$forecast)) {
    ## Even if pred_med was found in summary, intervals may be in forecast.
    if (is.null(pred_lo) && !is.null(res_m$forecast$pred_lo)) pred_lo <- res_m$forecast$pred_lo
    if (is.null(pred_hi) && !is.null(res_m$forecast$pred_hi)) pred_hi <- res_m$forecast$pred_hi
  }
  
  ## ------------------------------------------------------------
  ## 3. Try pred slot directly, useful for GNAR ts forecast
  ## ------------------------------------------------------------
  if (is.null(pred_med) && !is.null(res_m$pred)) {
    pred_med <- res_m$pred
  }
  
  ## ------------------------------------------------------------
  ## Convert ts/data.frame/vector-like objects to matrix
  ## ------------------------------------------------------------
  to_mat <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    x
  }
  
  pred_med <- to_mat(pred_med)
  pred_lo  <- to_mat(pred_lo)
  pred_hi  <- to_mat(pred_hi)
  
  ## ------------------------------------------------------------
  ## Check dimensions.
  ## Expected: horizon x N.
  ## If transposed, fix it.
  ## ------------------------------------------------------------
  fix_dim <- function(mat, name) {
    if (is.null(mat)) return(NULL)
    
    if (all(dim(mat) == c(horizon, N))) {
      return(mat)
    }
    
    if (all(dim(mat) == c(N, horizon))) {
      return(t(mat))
    }
    
    warning(
      name, " has unexpected dimension: ",
      paste(dim(mat), collapse = " x "),
      ". Expected ", horizon, " x ", N,
      " or ", N, " x ", horizon, ". Returning NULL."
    )
    NULL
  }
  
  pred_med <- fix_dim(pred_med, "pred_med")
  pred_lo  <- fix_dim(pred_lo, "pred_lo")
  pred_hi  <- fix_dim(pred_hi, "pred_hi")
  
  list(
    pred_med = pred_med,
    pred_lo = pred_lo,
    pred_hi = pred_hi
  )
}


plot_prediction_node_compare_from_cmp <- function(cmp, X_test, node = 1,
                                                  interval_methods = "BGNAR",
                                                  method_order = NULL,
                                                  title = NULL,
                                                  y_label = "Standardized value",
                                                  method_labels = NULL) {
  if (is.null(cmp$results)) {
    stop("cmp$results is NULL.")
  }
  
  horizon <- nrow(X_test)
  N <- ncol(X_test)
  h <- seq_len(horizon)
  
  method_names <- names(cmp$results)
  
  if (is.null(method_order)) {
    method_order <- method_names
  } else {
    method_order <- intersect(method_order, method_names)
  }
  
  ## ------------------------------------------------------------
  ## Actual observed values
  ## ------------------------------------------------------------
  actual_df <- data.frame(
    h = h,
    value = X_test[, node],
    method = "Actual"
  )
  
  ## ------------------------------------------------------------
  ## Prediction lines and ribbons
  ## ------------------------------------------------------------
  pred_df_list <- list()
  ribbon_df_list <- list()
  
  for (method_name in method_order) {
    res_m <- cmp$results[[method_name]]
    
    mats <- extract_forecast_mats_from_result(
      res_m = res_m,
      horizon = horizon,
      N = N
    )
    
    if (is.null(mats$pred_med)) {
      warning("Skipping ", method_name, ": pred_med not found.")
      next
    }
    
    pred_df_list[[method_name]] <- data.frame(
      h = h,
      value = mats$pred_med[, node],
      method = method_name
    )
    
    if (method_name %in% interval_methods &&
        !is.null(mats$pred_lo) &&
        !is.null(mats$pred_hi)) {
      ribbon_df_list[[method_name]] <- data.frame(
        h = h,
        lo = mats$pred_lo[, node],
        hi = mats$pred_hi[, node],
        method = method_name
      )
    }
  }
  
  if (length(pred_df_list) == 0) {
    stop("No valid predictions found in cmp$results.")
  }
  
  if (is.null(method_labels)) {
    method_labels <- setNames(method_order, method_order)
  }
  display_order <- unname(method_labels[method_order])
  pred_df <- do.call(rbind, pred_df_list)
  pred_df$method <- factor(unname(method_labels[as.character(pred_df$method)]),
                           levels = display_order)
  
  ribbon_df <- if (length(ribbon_df_list) > 0) {
    out <- do.call(rbind, ribbon_df_list)
    out$method <- factor(unname(method_labels[as.character(out$method)]),
                         levels = display_order)
    out
  } else {
    NULL
  }
  
  if (is.null(title)) {
    title <- paste0("Forecast comparison: node ", node)
  }
  
  ## ------------------------------------------------------------
  ## Plot
  ## ------------------------------------------------------------
  p <- ggplot2::ggplot()
  
  if (!is.null(ribbon_df)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = ribbon_df,
        ggplot2::aes(x = h, ymin = lo, ymax = hi, fill = method),
        alpha = 0.15
      )
  }
  
  p <- p +
    ggplot2::geom_line(
      data = actual_df,
      ggplot2::aes(x = h, y = value),
      color = "black",
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = pred_df,
      ggplot2::aes(x = h, y = value, color = method, linetype = method),
      linewidth = 0.85
    ) +
    ggplot2::labs(
      title = title,
      x = "Forecast horizon",
      y = y_label,
      color = "Method",
      linetype = "Method",
      fill = "Prediction interval"
    ) +
    ggplot2::theme_minimal()
  
  p
}

simulate_gnar <- function(Tn, N, p, W_list, s_vec,
                          mu_true, alpha_true_mat, beta_true_list,
                          sigma_true, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(0, Tn, N)
  X[seq_len(p), ] <- matrix(rnorm(p * N, 0, 0.5), p, N)
  for (t in (p + 1):Tn) {
    mean_t <- mu_true
    for (j in seq_len(p)) {
      xlag <- X[t - j, ]
      mean_t <- mean_t + alpha_true_mat[, j] * xlag
      for (r in seq_len(s_vec[j])) {
        mean_t <- mean_t + beta_true_list[[j]][r] * as.vector(W_list[[r]] %*% xlag)
      }
    }
    X[t, ] <- mean_t + rnorm(N, 0, sigma_true)
  }
  X
}

evaluate_forecast <- function(truth, pred_med, pred_lo = NULL, pred_hi = NULL) {
  truth <- as.matrix(truth)
  pred_med <- as.matrix(pred_med)
  
  if (!all(dim(truth) == dim(pred_med))) {
    stop(
      "Dimension mismatch in evaluate_forecast(): truth is ",
      paste(dim(truth), collapse = " x "),
      ", pred_med is ",
      paste(dim(pred_med), collapse = " x ")
    )
  }
  
  err <- pred_med - truth
  
  out <- list(
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    mae = mean(abs(err), na.rm = TRUE)
  )
  
  if (!is.null(pred_lo) && !is.null(pred_hi)) {
    pred_lo <- as.matrix(pred_lo)
    pred_hi <- as.matrix(pred_hi)
    
    if (all(dim(truth) == dim(pred_lo)) && all(dim(truth) == dim(pred_hi))) {
      out$coverage <- mean(truth >= pred_lo & truth <= pred_hi, na.rm = TRUE)
      out$avg_width <- mean(pred_hi - pred_lo, na.rm = TRUE)
    } else {
      out$coverage <- NA_real_
      out$avg_width <- NA_real_
    }
  } else {
    out$coverage <- NA_real_
    out$avg_width <- NA_real_
  }
  
  out
}


check_gnar_stationarity <- function(alpha_true_mat, beta_true_list) {
  N <- nrow(alpha_true_mat)
  p <- ncol(alpha_true_mat)
  

  if (length(beta_true_list) != p) {
    stop("length(beta_true_list) must equal p.")
  }
  
  sum_beta <- sum( abs(unlist( beta_true_list )) )
  
  companion <- NULL
  for(i in seq_len(N)) {
    sum_abs_alpha_i <- sum(abs(alpha_true_mat[i, ]))
    
    companion <- c(companion, sum_abs_alpha_i + sum_beta)
  }
  rho = max(companion)
  
  list(
    rho = rho,
    stationary = rho < 1
  )
}

gnar_spectral_radius <- function(alpha_true_mat, beta_true_list, W_list, s_vec) {
  ## Exact companion-matrix spectral radius for a time-invariant GNAR DGP.
  Phi <- make_phi_true(alpha_true_mat, beta_true_list, W_list, s_vec)
  N <- nrow(Phi[[1]])
  p <- length(Phi)
  if (p == 1L) {
    companion <- Phi[[1]]
  } else {
    companion <- rbind(
      do.call(cbind, Phi),
      cbind(diag(N * (p - 1L)), matrix(0, N * (p - 1L), N))
    )
  }
  max(Mod(eigen(companion, only.values = TRUE)$values))
}

make_phi_true <- function(alpha_true_mat, beta_true_list, W_list, s_vec) {
  alpha_true_mat <- as.matrix(alpha_true_mat)
  N <- nrow(alpha_true_mat)
  p <- ncol(alpha_true_mat)

  if (length(s_vec) != p) {
    stop("make_phi_true(): length(s_vec) must equal ncol(alpha_true_mat).")
  }
  if (length(beta_true_list) != p) {
    stop("make_phi_true(): length(beta_true_list) must equal ncol(alpha_true_mat).")
  }
  if (any(s_vec < 0) || any(s_vec != as.integer(s_vec))) {
    stop("make_phi_true(): s_vec must contain nonnegative integers.")
  }
  max_order <- if (length(s_vec)) max(s_vec) else 0L
  if (max_order > length(W_list)) {
    stop(
      "make_phi_true(): s_vec requires ", max_order,
      " neighborhood matrices, but W_list contains ", length(W_list), "."
    )
  }
  if (max_order > 0L) {
    bad_w <- which(!vapply(
      W_list[seq_len(max_order)],
      function(W) is.matrix(W) && identical(dim(W), c(N, N)),
      logical(1)
    ))
    if (length(bad_w)) {
      dims <- vapply(
        W_list[bad_w],
        function(W) paste(dim(as.matrix(W)), collapse = " x "),
        character(1)
      )
      stop(
        "make_phi_true(): W_list dimension mismatch at order ",
        paste(bad_w, collapse = ", "), "; expected ", N, " x ", N,
        ", found ", paste(dims, collapse = ", "), "."
      )
    }
  }
  too_short <- which(lengths(beta_true_list) < s_vec)
  if (length(too_short)) {
    stop(
      "make_phi_true(): beta_true_list is shorter than s_vec at temporal lag ",
      paste(too_short, collapse = ", "), "."
    )
  }
  
  Phi_true <- vector("list", p)
  
  for (j in seq_len(p)) {
    Phi_j <- diag(alpha_true_mat[, j], nrow = N, ncol = N)
    
    if (s_vec[j] > 0) {
      for (r in seq_len(s_vec[j])) {
        Phi_j <- Phi_j + beta_true_list[[j]][r] * W_list[[r]]
      }
    }
    
    Phi_true[[j]] <- Phi_j
  }
  
  Phi_true
}


make_phi_posterior_mean <- function(fit) {
  N <- fit$N
  p <- fit$p
  s_vec <- fit$s_vec
  W_list <- fit$W_list
  idx <- fit$idx
  G_mean <- colMeans(fit$gamma_draws)
  
  Phi_hat <- vector("list", p)
  
  for (j in seq_len(p)) {
    alpha_j <- G_mean[idx[[j]]$alpha]
    beta_j <- G_mean[idx[[j]]$beta]
    
    Phi_j <- diag(alpha_j, nrow = N, ncol = N)
    
    if (s_vec[j] > 0) {
      for (r in seq_len(s_vec[j])) {
        Phi_j <- Phi_j + beta_j[r] * W_list[[r]]
      }
    }
    
    Phi_hat[[j]] <- Phi_j
  }
  
  Phi_hat
}


compare_phi <- function(Phi_true, Phi_hat) {
  p <- length(Phi_true)
  
  data.frame(
    lag = seq_len(p),
    frob_error = sapply(seq_len(p), function(j) {
      sqrt(sum((Phi_hat[[j]] - Phi_true[[j]])^2))
    }),
    true_frob = sapply(seq_len(p), function(j) {
      sqrt(sum(Phi_true[[j]]^2))
    }),
    relative_frob_error = sapply(seq_len(p), function(j) {
      numerator <- sqrt(sum((Phi_hat[[j]] - Phi_true[[j]])^2))
      denominator <- sqrt(sum(Phi_true[[j]]^2))
      if (denominator > 0) numerator / denominator else NA_real_
    }),
    spectral_error = sapply(seq_len(p), function(j) {
      max(svd(Phi_hat[[j]] - Phi_true[[j]])$d)
    })
  )
}

extract_generic_fitted <- function(fit, X_train, p) {
  ## Generic extractor for package models (BVAR/GNAR).
  ## Returns NULL matrices if fitted() is not available or dimensions cannot be parsed.
  X_train <- as.matrix(X_train)
  T_fit <- nrow(X_train) - p
  N <- ncol(X_train)

  fit_obj <- tryCatch(stats::fitted(fit), error = function(e) NULL)
  if (is.null(fit_obj)) {
    return(list(fit_med = NULL, fit_lo = NULL, fit_hi = NULL))
  }

  fit_med <- tryCatch(
    as_prediction_matrix(fit_obj, horizon = T_fit, N = N),
    error = function(e) NULL
  )

  list(fit_med = fit_med, fit_lo = NULL, fit_hi = NULL)
}

extract_fitted_mats_from_result <- function(res_m, T_fit, N) {
  fit_med <- fit_lo <- fit_hi <- NULL

  if (!is.null(res_m$fitted_summary)) {
    if (!is.null(res_m$fitted_summary$fit_med)) fit_med <- res_m$fitted_summary$fit_med
    if (!is.null(res_m$fitted_summary$fit_lo))  fit_lo  <- res_m$fitted_summary$fit_lo
    if (!is.null(res_m$fitted_summary$fit_hi))  fit_hi  <- res_m$fitted_summary$fit_hi
  }

  to_mat <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    x
  }

  fix_dim <- function(mat, name) {
    mat <- to_mat(mat)
    if (is.null(mat)) return(NULL)
    if (all(dim(mat) == c(T_fit, N))) return(mat)
    if (all(dim(mat) == c(N, T_fit))) return(t(mat))
    warning(name, " has unexpected dimension: ", paste(dim(mat), collapse = " x "),
            ". Expected ", T_fit, " x ", N, " or ", N, " x ", T_fit,
            ". Returning NULL.")
    NULL
  }

  list(
    fit_med = fix_dim(fit_med, "fit_med"),
    fit_lo = fix_dim(fit_lo, "fit_lo"),
    fit_hi = fix_dim(fit_hi, "fit_hi")
  )
}

plot_fitted_node_compare_from_cmp <- function(cmp, X_scaled, train_end, p, node = 1,
                                              interval_methods = "BGNAR",
                                              method_order = NULL,
                                              title = NULL) {
  if (is.null(cmp$results)) stop("cmp$results is NULL.")

  X_train <- as.matrix(X_scaled[seq_len(train_end), , drop = FALSE])
  X_fit <- X_train[(p + 1):train_end, , drop = FALSE]
  T_fit <- nrow(X_fit)
  N <- ncol(X_fit)
  tt <- seq_len(T_fit)

  if (is.null(method_order)) {
    method_order <- names(cmp$results)
  } else {
    method_order <- intersect(method_order, names(cmp$results))
  }

  actual_df <- data.frame(
    t = tt,
    value = X_fit[, node],
    method = "Actual"
  )

  fit_df_list <- list()
  ribbon_df_list <- list()

  for (method_name in method_order) {
    res_m <- cmp$results[[method_name]]
    mats <- extract_fitted_mats_from_result(res_m = res_m, T_fit = T_fit, N = N)

    if (is.null(mats$fit_med)) {
      warning("Skipping ", method_name, ": fit_med not found.")
      next
    }

    fit_df_list[[method_name]] <- data.frame(
      t = tt,
      value = mats$fit_med[, node],
      method = method_name
    )

    if (method_name %in% interval_methods &&
        !is.null(mats$fit_lo) &&
        !is.null(mats$fit_hi)) {
      ribbon_df_list[[method_name]] <- data.frame(
        t = tt,
        lo = mats$fit_lo[, node],
        hi = mats$fit_hi[, node],
        method = method_name
      )
    }
  }

  if (length(fit_df_list) == 0) stop("No valid fitted values found in cmp$results.")

  fit_df <- do.call(rbind, fit_df_list)
  fit_df$method <- factor(fit_df$method, levels = method_order)

  ribbon_df <- if (length(ribbon_df_list) > 0) {
    out <- do.call(rbind, ribbon_df_list)
    out$method <- factor(out$method, levels = method_order)
    out
  } else {
    NULL
  }

  if (is.null(title)) title <- paste0("Fitted value comparison: node ", node)

  p_plot <- ggplot2::ggplot()

  if (!is.null(ribbon_df)) {
    p_plot <- p_plot +
      ggplot2::geom_ribbon(
        data = ribbon_df,
        ggplot2::aes(x = t, ymin = lo, ymax = hi, fill = method),
        alpha = 0.15
      )
  }

  p_plot +
    ggplot2::geom_line(
      data = actual_df,
      ggplot2::aes(x = t, y = value),
      color = "black",
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = fit_df,
      ggplot2::aes(x = t, y = value, color = method, linetype = method),
      linewidth = 0.85
    ) +
    ggplot2::labs(
      title = title,
      x = "Training time index after lag p",
      y = paste0("Node ", node, " value"),
      color = "Method",
      linetype = "Method",
      fill = "Fitted interval"
    ) +
    ggplot2::theme_minimal()
}
