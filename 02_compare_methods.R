############################################################
## 02_compare_methods.R
## Wrappers for BGNAR, BVAR, and GNAR package comparisons
############################################################

############################################################
## BGNAR wrapper
############################################################

fit_predict_bgnar <- function(X_scaled, train_end, p, W_list, s_vec,
                              horizon = nrow(X_scaled) - train_end,
                              n_iter = 3000, burn = 1000, thin = 5,
                              use_dummy = TRUE,
                              b = 0,
                              forecast_type = c("recursive", "one_step"),
                              simulate = TRUE,
                              point_source = c("pred_mean", "pred_y"),
                              ...) {
  forecast_type <- match.arg(forecast_type)
  point_source <- match.arg(point_source)
  
  fit <- fit_bgnar(
    X_train = X_scaled[seq_len(train_end), , drop = FALSE],
    p = p,
    W_list = W_list,
    s_vec = s_vec,
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    use_dummy = use_dummy,
    b = b,
    ...
  )
  
  if (forecast_type == "recursive") {
    pred <- predict_bgnar_recursive(
      fit = fit,
      X_history = X_scaled[seq_len(train_end), , drop = FALSE],
      horizon = horizon,
      simulate = simulate
    )
  } else {
    pred <- predict_bgnar_one_step(
      fit = fit,
      X_all = X_scaled,
      train_end = train_end,
      horizon = horizon,
      simulate = simulate
    )
  }
  
  X_train <- X_scaled[seq_len(train_end), , drop = FALSE]
  X_test <- X_scaled[(train_end + 1):(train_end + horizon), , drop = FALSE]
  
  sum_pred <- summarize_bgnar_prediction(
    pred,
    X_test,
    point_source = point_source
  )
  
  fitted <- fitted_bgnar(
    fit = fit,
    X_train = X_train,
    simulate = simulate
  )
  
  sum_fit <- summarize_bgnar_fitted(
    fitted,
    X_train = X_train,
    p = p,
    point_source = "fit_mean"
  )
  
  list(
    method = paste0("BGNAR_", forecast_type),
    fit = fit,
    pred = pred,
    fitted = fitted,
    fitted_summary = sum_fit,
    summary = sum_pred,
    performance = sum_pred
  )
}

############################################################
## BVAR wrapper
############################################################

extract_bvar_forecast <- function(pred_obj, horizon, N) {
  sm <- summary(pred_obj)
  
  ## BVAR summary stores quantiles as:
  ## sm$quants: quantile x horizon x variable
  ## dim = 3 x horizon x N
  if (!is.null(sm$quants)) {
    qarr <- sm$quants
  } else if (!is.null(sm$fcast) && is.list(sm$fcast) && !is.null(sm$fcast$quants)) {
    qarr <- sm$fcast$quants
  } else {
    stop("Could not find BVAR forecast quantiles in summary(pred_obj).")
  }
  
  qarr <- as.array(qarr)
  dd <- dim(qarr)
  
  if (length(dd) != 3) {
    stop("BVAR quantile object is not a 3D array. dim = ", paste(dd, collapse = " x "))
  }
  
  ## Expected: quantile x horizon x variable
  if (!(dd[2] == horizon && dd[3] == N)) {
    stop(
      "Unexpected BVAR quantile dimensions. dim = ",
      paste(dd, collapse = " x "),
      ", expected 3 x ", horizon, " x ", N
    )
  }
  
  qnames <- dimnames(qarr)[[1]]
  
  if (!is.null(qnames)) {
    q_lower <- grep("2.5|0.025|lower|low", qnames, ignore.case = TRUE)
    q_mid   <- grep("50|0.5|median|mean", qnames, ignore.case = TRUE)
    q_upper <- grep("97.5|0.975|upper|high", qnames, ignore.case = TRUE)
    
    if (length(q_lower) == 0) q_lower <- 1
    if (length(q_mid) == 0) q_mid <- 2
    if (length(q_upper) == 0) q_upper <- 3
  } else {
    q_lower <- 1
    q_mid <- 2
    q_upper <- 3
  }
  
  q_lower <- q_lower[1]
  q_mid <- q_mid[1]
  q_upper <- q_upper[1]
  
  ## Convert to horizon x N
  pred_lo_raw  <- as.matrix(qarr[q_lower, , ])
  pred_med     <- as.matrix(qarr[q_mid, , ])
  pred_hi_raw  <- as.matrix(qarr[q_upper, , ])
  
  ## Safety: ensure lower <= upper
  pred_lo <- pmin(pred_lo_raw, pred_hi_raw)
  pred_hi <- pmax(pred_lo_raw, pred_hi_raw)
  
  if (!all(dim(pred_med) == c(horizon, N))) {
    stop(
      "Parsed BVAR forecast has wrong dimension. pred_med dim = ",
      paste(dim(pred_med), collapse = " x "),
      ", expected ",
      horizon, " x ", N
    )
  }
  
  list(
    pred_med = pred_med,
    pred_lo = pred_lo,
    pred_hi = pred_hi
  )
}

fit_predict_bvar <- function(X_scaled, train_end, p, horizon,
                             draws = 3000L, burn = 1000L,
                             prior_mean = 1, use_dummy = TRUE) {
  if (!requireNamespace("BVAR", quietly = TRUE)) {
    warning("Package BVAR is not installed. Skipping BVAR.")
    return(NULL)
  }
  
  X_train <- as.matrix(X_scaled[seq_len(train_end), , drop = FALSE])
  X_test <- as.matrix(X_scaled[(train_end + 1):(train_end + horizon), , drop = FALSE])
  
  make_priors <- function(psi = BVAR::bv_psi()) {
    mn <- BVAR::bv_mn(
      lambda = BVAR::bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
      alpha = BVAR::bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3),
      psi = psi,
      var = 1e7,
      b = prior_mean
    )
    if (use_dummy) {
      BVAR::bv_priors(
        hyper = "auto",
        mn = mn,
        soc = BVAR::bv_soc(mode = 1, sd = 1, min = 0.0001, max = 50),
        sur = BVAR::bv_sur(mode = 1, sd = 1, min = 0.0001, max = 50)
      )
    } else {
      BVAR::bv_priors(hyper = "auto", mn = mn)
    }
  }
  
  run_bvar <- function(priors) {
    BVAR::bvar(
      data = X_train,
      lags = p,
      n_draw = as.integer(draws),
      n_burn = as.integer(burn),
      priors = priors,
      verbose = FALSE
    )
  }
  
  ## BVAR's automatic psi routine fits a univariate AR(p) to every column.
  ## Very persistent short samples can make that initialization fail even
  ## though the Bayesian VAR itself is estimable.  Retry only those cases with
  ## the same pre-declared AR(1)-residual rule for every affected series.
  fit <- tryCatch(
    run_bvar(make_priors()),
    error = function(e) {
      warning(
        "BVAR automatic psi failed; retrying with robust AR(1) residual scales: ",
        conditionMessage(e)
      )
      psi_mode <- vapply(seq_len(ncol(X_train)), function(j) {
        x <- X_train[, j]
        ar1 <- stats::lm(x[-1L] ~ x[-length(x)])
        value <- stats::var(stats::residuals(ar1))
        if (!is.finite(value) || value <= 0) value <- stats::var(diff(x))
        max(value, 1e-6, na.rm = TRUE)
      }, numeric(1L))
      run_bvar(make_priors(BVAR::bv_psi(mode = psi_mode)))
    }
  )
  
  pred <- predict(
    fit,
    horizon = as.integer(horizon),
    conf_bands = 0.025
  )
  
  fc <- extract_bvar_forecast(
    pred_obj = pred,
    horizon = horizon,
    N = ncol(X_train)
  )
  
  perf <- evaluate_forecast(
    truth = X_test,
    pred_med = fc$pred_med,
    pred_lo = fc$pred_lo,
    pred_hi = fc$pred_hi
  )
  
  fitted <- extract_generic_fitted(
    fit = fit,
    X_train = X_train,
    p = p
  )
  fit_perf <- if (!is.null(fitted$fit_med)) {
    evaluate_forecast(
      truth = X_train[(p + 1):nrow(X_train), , drop = FALSE],
      pred_med = fitted$fit_med,
      pred_lo = fitted$fit_lo,
      pred_hi = fitted$fit_hi
    )
  } else {
    list(rmse = NA_real_, mae = NA_real_, coverage = NA_real_, avg_width = NA_real_)
  }
  
  list(
    method = "BVAR",
    fit = fit,
    pred = pred,
    forecast = fc,
    fitted_summary = c(fitted, fit_perf),
    summary = perf,
    performance = perf,
    selected_order = data.frame(
      p = p,
      criterion = "fixed from BVAR_P_MAX"
    )
  )
}

############################################################
## Random-walk benchmark
############################################################

fit_predict_random_walk <- function(X_scaled, train_end, horizon,
                                    interval_level = 0.95) {
  X_train <- as.matrix(X_scaled[seq_len(train_end), , drop = FALSE])
  X_test <- as.matrix(X_scaled[(train_end + 1):(train_end + horizon), , drop = FALSE])
  pred_med <- matrix(
    X_train[nrow(X_train), ],
    nrow = horizon,
    ncol = ncol(X_train),
    byrow = TRUE
  )
  delta <- apply(X_train, 2L, diff)
  if (is.null(dim(delta))) delta <- matrix(delta, ncol = ncol(X_train))
  innovation_sd <- apply(delta, 2L, stats::sd, na.rm = TRUE)
  pooled_sd <- stats::sd(as.vector(delta), na.rm = TRUE)
  if (!is.finite(pooled_sd) || pooled_sd <= 0) pooled_sd <- 1e-6
  innovation_sd[!is.finite(innovation_sd) | innovation_sd <= 0] <- pooled_sd
  z <- stats::qnorm((1 + interval_level) / 2)
  step_sd <- outer(sqrt(seq_len(horizon)), innovation_sd)
  pred_lo <- pred_med - z * step_sd
  pred_hi <- pred_med + z * step_sd
  perf <- evaluate_forecast(X_test, pred_med, pred_lo, pred_hi)
  list(
    method = "Random_walk",
    forecast = list(pred_med = pred_med, pred_lo = pred_lo, pred_hi = pred_hi),
    fitted_summary = list(
      fit_med = NULL, fit_lo = NULL, fit_hi = NULL,
      rmse = NA_real_, mae = NA_real_, coverage = NA_real_, avg_width = NA_real_
    ),
    summary = perf,
    performance = perf
  )
}

############################################################
## GNAR wrapper
############################################################

igraph_to_gnar_net <- function(g) {
  if (!requireNamespace("GNAR", quietly = TRUE)) {
    stop("Package GNAR is not installed.")
  }
  
  ## GNAR has both matrixtoGNAR() and igraphtoGNAR().
  ## igraphtoGNAR() is more direct, but matrixtoGNAR() is kept as fallback.
  net <- tryCatch(
    GNAR::igraphtoGNAR(g),
    error = function(e) {
      A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
      GNAR::matrixtoGNAR(A)
    }
  )
  
  net
}

as_prediction_matrix <- function(pred_obj, horizon, N) {
  obj <- pred_obj
  
  if (is.data.frame(obj)) obj <- as.matrix(obj)
  if (is.ts(obj)) obj <- as.matrix(obj)
  
  if (is.vector(obj) && length(obj) == horizon * N) {
    obj <- matrix(obj, nrow = horizon, ncol = N, byrow = TRUE)
  }
  
  if (is.matrix(obj)) {
    if (nrow(obj) == horizon && ncol(obj) == N) {
      return(obj)
    }
    if (nrow(obj) == N && ncol(obj) == horizon) {
      return(t(obj))
    }
  }
  
  if (is.array(obj)) {
    dd <- dim(obj)
    h_dim <- which(dd == horizon)[1]
    n_dim <- which(dd == N)[1]
    if (!is.na(h_dim) && !is.na(n_dim)) {
      other_dim <- setdiff(seq_along(dd), c(h_dim, n_dim))
      if (length(other_dim) == 0) {
        mat <- aperm(obj, c(h_dim, n_dim))
        return(as.matrix(mat))
      } else {
        ## Use first slice of remaining dimension
        perm <- c(h_dim, n_dim, other_dim)
        arr <- aperm(obj, perm)
        mat <- arr[, , 1]
        return(as.matrix(mat))
      }
    }
  }
  
  stop("Cannot convert prediction object to horizon x N matrix.")
}

fit_predict_gnar <- function(X_scaled, train_end, g, p, s_vec, horizon,
                             globalalpha = FALSE) {
  if (!requireNamespace("GNAR", quietly = TRUE)) {
    warning("Package GNAR is not installed. Skipping GNAR.")
    return(NULL)
  }
  
  X_train <- as.matrix(X_scaled[seq_len(train_end), , drop = FALSE])
  X_test <- as.matrix(X_scaled[(train_end + 1):(train_end + horizon), , drop = FALSE])
  
  net <- igraph_to_gnar_net(g)
  
  fit <- GNAR::GNARfit(
    vts = X_train,
    net = net,
    alphaOrder = p,
    betaOrder = s_vec,
    globalalpha = globalalpha,
    ErrorIfNoNei = FALSE
  )
  
  pred_obj <- predict(
    fit,
    n.ahead = horizon
  )
  
  pred_med <- as_prediction_matrix(
    pred_obj = pred_obj,
    horizon = horizon,
    N = ncol(X_train)
  )
  
  perf <- evaluate_forecast(
    truth = X_test,
    pred_med = pred_med,
    pred_lo = NULL,
    pred_hi = NULL
  )
  
  fitted <- extract_generic_fitted(
    fit = fit,
    X_train = X_train,
    p = p
  )
  fit_perf <- if (!is.null(fitted$fit_med)) {
    evaluate_forecast(
      truth = X_train[(p + 1):nrow(X_train), , drop = FALSE],
      pred_med = fitted$fit_med,
      pred_lo = fitted$fit_lo,
      pred_hi = fitted$fit_hi
    )
  } else {
    list(rmse = NA_real_, mae = NA_real_, coverage = NA_real_, avg_width = NA_real_)
  }
  
  list(
    method = "GNAR",
    fit = fit,
    pred = pred_obj,
    forecast = list(
      pred_med = pred_med,
      pred_lo = NULL,
      pred_hi = NULL
    ),
    fitted_summary = c(fitted, fit_perf),
    summary = perf,
    performance = perf
  )
}

make_gnar_order_candidates <- function(p_max = 3L, s_max = 3L) {
  out <- list()
  k <- 1L
  for (p in seq_len(p_max)) {
    grid <- expand.grid(rep(list(0:s_max), p))
    for (i in seq_len(nrow(grid))) {
      s <- as.integer(grid[i, ])
      ## A weak hierarchy avoids implausibly complex distant-lag structures.
      if (length(s) > 1L && any(diff(s) > 0)) next
      out[[k]] <- list(p = p, s_vec = s)
      k <- k + 1L
    }
  }
  out
}

fit_predict_gnar_bic <- function(X_scaled, train_end, g, horizon,
                                 p_max = 3L, s_max = 3L,
                                 globalalpha_grid = c(TRUE, FALSE)) {
  if (!requireNamespace("GNAR", quietly = TRUE)) {
    warning("Package GNAR is not installed. Skipping GNAR.")
    return(NULL)
  }
  X_train <- as.matrix(X_scaled[seq_len(train_end), , drop = FALSE])
  X_test <- as.matrix(X_scaled[(train_end + 1):(train_end + horizon), , drop = FALSE])
  net <- igraph_to_gnar_net(g)
  candidates <- make_gnar_order_candidates(p_max, s_max)
  fits <- list()
  score_rows <- list()
  counter <- 1L
  
  for (cand in candidates) {
    for (ga in globalalpha_grid) {
      ## Compare BIC values on a common effective sample.  Without this
      ## alignment, larger p discards more initial rows and its likelihood/BIC
      ## is not directly comparable with lower-order candidates.
      score_start <- p_max - cand$p + 1L
      X_score <- X_train[score_start:nrow(X_train), , drop = FALSE]
      fit_score <- tryCatch(
        GNAR::GNARfit(
          vts = X_score, net = net,
          alphaOrder = cand$p, betaOrder = cand$s_vec,
          globalalpha = ga, ErrorIfNoNei = FALSE
        ),
        error = function(e) NULL
      )
      score <- if (is.null(fit_score)) Inf else tryCatch(
        as.numeric(stats::BIC(fit_score)), error = function(e) Inf
      )
      fits[[counter]] <- fit_score
      score_rows[[counter]] <- data.frame(
        candidate = counter, p = cand$p,
        s_vec = paste(cand$s_vec, collapse = ";"),
        globalalpha = ga, score_start = score_start,
        effective_time = nrow(X_score) - cand$p, bic = score
      )
      counter <- counter + 1L
    }
  }
  scores <- do.call(rbind, score_rows)
  if (!any(is.finite(scores$bic))) stop("All GNAR BIC candidates failed.")
  best_id <- scores$candidate[which.min(scores$bic)]
  best <- scores[scores$candidate == best_id, , drop = FALSE]
  p_selected <- as.integer(best$p[[1L]])
  s_selected <- as.integer(strsplit(
    best$s_vec[[1L]], ";", fixed = TRUE
  )[[1L]])
  global_selected <- as.logical(best$globalalpha[[1L]])
  ## Refit the selected specification to all training observations before
  ## forecasting; the truncated fit above is used only for comparable BIC.
  fit <- GNAR::GNARfit(
    vts = X_train, net = net,
    alphaOrder = p_selected, betaOrder = s_selected,
    globalalpha = global_selected, ErrorIfNoNei = FALSE
  )
  pred_obj <- predict(fit, n.ahead = horizon)
  pred_med <- as_prediction_matrix(pred_obj, horizon, ncol(X_train))
  perf <- evaluate_forecast(X_test, pred_med)
  fitted <- extract_generic_fitted(fit, X_train, p_selected)
  fit_perf <- if (!is.null(fitted$fit_med)) {
    evaluate_forecast(
      X_train[(p_selected + 1):nrow(X_train), , drop = FALSE],
      fitted$fit_med
    )
  } else {
    list(rmse = NA_real_, mae = NA_real_, coverage = NA_real_, avg_width = NA_real_)
  }
  list(
    method = "GNAR_BIC", fit = fit, pred = pred_obj,
    forecast = list(pred_med = pred_med, pred_lo = NULL, pred_hi = NULL),
    fitted_summary = c(fitted, fit_perf),
    summary = perf, performance = perf,
    selected_order = best,
    selection_table = scores
  )
}

## GNAR has no intercept in the package implementation.  For comparisons that
## explicitly request centering, remove each node's training-window mean before
## BIC selection and fitting, then add that same training mean back to every
## point forecast.  Held-out observations never enter the centering constants.
fit_predict_gnar_bic_centered <- function(
    X_scaled, train_end, g, horizon, p_max = 3L, s_max = 3L,
    globalalpha_grid = c(TRUE, FALSE)) {
  X_scaled <- as.matrix(X_scaled)
  train_mean <- colMeans(X_scaled[seq_len(train_end), , drop = FALSE])
  X_centered <- sweep(X_scaled, 2L, train_mean, FUN = "-")
  out <- fit_predict_gnar_bic(
    X_scaled = X_centered, train_end = train_end, g = g,
    horizon = horizon, p_max = p_max, s_max = s_max,
    globalalpha_grid = globalalpha_grid
  )
  out$forecast$pred_med <- sweep(
    out$forecast$pred_med, 2L, train_mean, FUN = "+"
  )
  out$forecast$pred_lo <- NULL
  out$forecast$pred_hi <- NULL
  truth <- X_scaled[(train_end + 1L):(train_end + horizon), , drop = FALSE]
  out$summary <- evaluate_forecast(truth, out$forecast$pred_med)
  out$performance <- out$summary
  out$method <- "GNAR_BIC_nodewise_centered"
  out$training_node_means <- train_mean
  out
}
