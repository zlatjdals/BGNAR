############################################################
## 01_bgnar_model.R
## Bayesian GNAR with Minnesota-type prior + dummy observations
############################################################

build_H_one <- function(x_lags, W_list, s_vec,
                        alpha_mode = c("independent", "global", "local")) {
  alpha_mode <- match.arg(alpha_mode)
  p <- length(x_lags)
  N <- length(x_lags[[1]])
  H <- diag(N)
  for (j in seq_len(p)) {
    xlag <- x_lags[[j]]
    G_alpha <- if (alpha_mode == "global") {
      matrix(xlag, nrow = N, ncol = 1L)
    } else {
      diag(xlag, nrow = N, ncol = N)
    }
    G_beta <- NULL
    if (s_vec[j] > 0) {
      for (r in seq_len(s_vec[j])) {
        G_beta <- cbind(G_beta, W_list[[r]] %*% xlag)
      }
    }
    H <- cbind(H, cbind(G_alpha, G_beta))
  }
  H
}

make_observed_regression <- function(X, p, W_list, s_vec,
                                     alpha_mode = "independent") {
  Tn <- nrow(X)
  H_list <- list()
  x_list <- list()
  idx <- 1
  for (t in (p + 1):Tn) {
    x_lags <- lapply(seq_len(p), function(j) X[t - j, ])
    H_list[[idx]] <- build_H_one(x_lags, W_list, s_vec, alpha_mode)
    x_list[[idx]] <- X[t, ]
    idx <- idx + 1
  }
  list(H = do.call(rbind, H_list), x = as.vector(do.call(c, x_list)))
}

make_dummy_observations <- function(X_train, p, W_list, s_vec,
                                    tau_soc = NULL, tau_dio = NULL,
                                    alpha_mode = "independent",
                                    reduce_soc = FALSE) {
  ## tau_soc and tau_dio are kept only for backward compatibility.
  ## In the updated model, dummy observations are NOT divided by tau.
  ## Instead, tau_soc^2 and tau_dio^2 enter the pseudo-error variances.
  
  N <- ncol(X_train)
  xbar0 <- colMeans(X_train[seq_len(p), , drop = FALSE])
  
  ## ------------------------------------------------------------
  ## Sum-of-coefficients dummy restriction
  ## X_soc = diag(xbar0)
  ## Y_soc = [0; X_soc; ...; X_soc]
  ## ------------------------------------------------------------
  X_soc <- diag(xbar0, nrow = N, ncol = N)
  
  H_soc_list <- list()
  x_soc_list <- list()
  
  for (k in seq_len(N)) {
    y0 <- 0
    lag_vec <- X_soc[, k]
    x_lags <- replicate(p, lag_vec, simplify = FALSE)
    
    H_y <- build_H_one(x_lags, W_list, s_vec, alpha_mode)
    
    ## Replace the intercept block by y0 * I_N.
    H_y[, seq_len(N)] <- y0 * diag(N)
    
    H_soc_list[[k]] <- H_y
    x_soc_list[[k]] <- X_soc[, k]
  }
  
  H_soc <- do.call(rbind, H_soc_list)
  x_soc <- as.vector(do.call(c, x_soc_list))
  soc_original_dim <- length(x_soc)

  ## Compress the N^2-row SOC block to the nonzero singular directions of
  ## [x_soc, H_soc].  For every gamma this preserves the original BVAR-style
  ## penalty ||x_soc - H_soc gamma||^2, including its xbar0 dependence.
  if (isTRUE(reduce_soc)) {
    M_soc <- cbind(x_soc, H_soc)
    sv_soc <- svd(M_soc, nu = 0L, nv = min(dim(M_soc)))
    tol_soc <- max(dim(M_soc)) * .Machine$double.eps * sv_soc$d[1L]
    rank_soc <- sum(sv_soc$d > tol_soc)
    B_soc <- sweep(
      t(sv_soc$v[, seq_len(rank_soc), drop = FALSE]),
      1L, sv_soc$d[seq_len(rank_soc)], "*"
    )
    x_soc <- as.vector(B_soc[, 1L])
    H_soc <- B_soc[, -1L, drop = FALSE]
  }
  soc_reduced_dim <- length(x_soc)
  
  ## ------------------------------------------------------------
  ## Dummy-initial-observation restriction
  ## X_dio = xbar0
  ## Y_dio = [1; xbar0; ...; xbar0]
  ## ------------------------------------------------------------
  X_dio <- matrix(xbar0, ncol = 1)
  
  y0 <- 1
  lag_vec <- xbar0
  x_lags <- replicate(p, lag_vec, simplify = FALSE)
  
  H_dio <- build_H_one(x_lags, W_list, s_vec, alpha_mode)
  H_dio[, seq_len(N)] <- y0 * diag(N)
  
  x_dio <- as.vector(X_dio)
  
  H_d <- rbind(H_soc, H_dio)
  x_d <- c(x_soc, x_dio)
  
  list(
    H = H_d,
    x = x_d,
    H_soc = H_soc,
    x_soc = x_soc,
    H_dio = H_dio,
    x_dio = x_dio,
    M_soc = N,
    M_dio = 1,
    M_d = N + 1,
    soc_original_dim = soc_original_dim,
    soc_reduced_dim = soc_reduced_dim
  )
}

make_theta_index <- function(N, p, s_vec,
                             alpha_mode = c("independent", "global", "local")) {
  alpha_mode <- match.arg(alpha_mode)
  idx <- list()
  start <- N + 1
  n_alpha <- if (alpha_mode == "global") 1L else N
  
  for (j in seq_len(p)) {
    alpha_idx <- start:(start + n_alpha - 1L)
    
    if (s_vec[j] > 0) {
      beta_idx <- (start + n_alpha):(start + n_alpha + s_vec[j] - 1)
    } else {
      beta_idx <- integer(0)
    }
    
    idx[[j]] <- list(alpha = alpha_idx, beta = beta_idx)
    start <- start + n_alpha + s_vec[j]
  }
  
  idx
}

make_prior_mean_var <- function(N, p, s_vec,
                                lambda_alpha2, lambda_beta2, a, d,
                                b = 0, center_all_lags = FALSE,
                                c_mu = 1e7,
                                alpha_mode = "independent",
                                kappa_alpha = 1,
                                b_lag = rep(0, p)) {
  alpha_mode <- match.arg(alpha_mode, c("independent", "global", "local"))
  n_alpha <- if (alpha_mode == "global") 1L else N
  K_theta <- sum(n_alpha + s_vec)
  K_gamma <- N + K_theta
  m <- numeric(K_gamma)
  V_diag <- numeric(K_gamma)
  m[seq_len(N)] <- 0
  V_diag[seq_len(N)] <- c_mu
  idx <- make_theta_index(N, p, s_vec, alpha_mode)
  for (j in seq_len(p)) {
    a_idx <- idx[[j]]$alpha
    b_idx <- idx[[j]]$beta
    m[a_idx] <- if (alpha_mode == "local") {
      rep(b_lag[[j]], length(a_idx))
    } else if (alpha_mode == "global") {
      0
    } else if (center_all_lags || j == 1) {
      rep(b[[min(j, length(b))]], length.out = N)
    } else {
      0
    }
    V_diag[a_idx] <- lambda_alpha2 * j^(-a) *
      if (alpha_mode == "local") kappa_alpha else 1
    r_seq <- seq_len(s_vec[j])
    m[b_idx] <- 0
    V_diag[b_idx] <- lambda_beta2 * j^(-a) * r_seq^(-d)
  }
  list(m = m, V_diag = V_diag, idx = idx)
}

## One common hierarchical prior is used for every simulation and real-data
## application. Hyperparameters are selected once on independent simulation
## calibration seeds and pre-holdout wind origins.
fit_bgnar <- function(X_train, p, W_list, s_vec,
                      n_iter = 3000, burn = 1000, thin = 2,
                      tau_soc = 1, tau_dio = 1, use_dummy = TRUE,
                      use_soc = use_dummy, use_dio = use_dummy,
                      dummy_error_covariance = c("identity", "innovation"),
                      reduce_soc_dummy = FALSE,
                      learn_dummy_tightness = TRUE,
                      dummy_precision_selection = FALSE,
                      dummy_precision_shape = 2,
                      dummy_precision_rate = 1,
                      dummy_tightness_mixture = FALSE,
                      tau_active = 0.5, tau_inactive = 2,
                      pi_soc_active = 0.5, pi_dio_active = 0.5,
                      initial_soc_active = 1L, initial_dio_active = 1L,
                      learn_dummy_mixture_probability = FALSE,
                      pi_soc_alpha = 1, pi_soc_beta = 1,
                      pi_dio_alpha = 1, pi_dio_beta = 1,
                      soc_variance_multiplier = 1,
                      A_grid = c(1, 1.5, 2, 2.5, 3),
                      D_grid = c(1, 1.5, 2, 2.5, 3),
                      A_prob = NULL, D_prob = NULL,
                      b = 0, c_mu = 1e7,
                      b_mode = c("fixed", "lag", "common", "node"),
                      b_prior_mean = 0, b_prior_var = 1,
                      alpha_mode = c("global", "local", "independent"),
                      c_kappa = 3, d_kappa = 2,
                      c_sigma = 2, d_sigma = 1,
                      a_alpha = 3, b_alpha = 0.16,
                      a_beta = 3, b_beta = 0.16,
                      ## SOC contributes N^2 scalar pseudo residuals whereas
                      ## DIO contributes N.  With a_soc=a_dio=3, the default
                      ## scales make E(tau^2 | Q=0)=1 for both restrictions.
                      a_soc = 3, b_soc = NULL,
                      a_dio = 3, b_dio = NULL,
                      verbose = TRUE) {
  X_train <- as.matrix(X_train)
  Tn <- nrow(X_train)
  N <- ncol(X_train)
  Tstar <- Tn - p
  b_mode <- match.arg(b_mode)
  alpha_mode <- match.arg(alpha_mode)
  dummy_error_covariance <- match.arg(dummy_error_covariance)
  reduce_soc_dummy <- isTRUE(reduce_soc_dummy)
  if (reduce_soc_dummy && dummy_error_covariance != "identity") {
    stop("Reduced SOC dummy requires dummy_error_covariance='identity'.")
  }
  if (alpha_mode != "independent" && b_mode != "fixed") {
    stop("b_mode must be 'fixed' for global- and local-alpha BGNAR.")
  }
  if (c_kappa <= 0 || d_kappa <= 0) {
    stop("c_kappa and d_kappa must be positive.")
  }
  if (!is.finite(b_prior_var) || b_prior_var <= 0) {
    stop("b_prior_var must be positive.")
  }
  use_soc <- isTRUE(use_soc)
  use_dio <- isTRUE(use_dio)
  use_dummy <- use_soc || use_dio
  dummy_precision_selection <- isTRUE(dummy_precision_selection)
  dummy_tightness_mixture <- isTRUE(dummy_tightness_mixture)
  learn_dummy_mixture_probability <-
    isTRUE(learn_dummy_mixture_probability)
  if (dummy_precision_selection &&
      (learn_dummy_tightness || dummy_tightness_mixture)) {
    stop(paste0(
      "Set learn_dummy_tightness=FALSE and dummy_tightness_mixture=FALSE ",
      "when dummy_precision_selection=TRUE."
    ))
  }
  if (dummy_precision_selection &&
      (!is.finite(dummy_precision_shape) || dummy_precision_shape <= 0 ||
       !is.finite(dummy_precision_rate) || dummy_precision_rate <= 0)) {
    stop("The Gamma slab shape and rate must be positive.")
  }
  if (dummy_tightness_mixture && learn_dummy_tightness) {
    stop(paste0(
      "Set learn_dummy_tightness=FALSE when ",
      "dummy_tightness_mixture=TRUE."
    ))
  }
  if (dummy_tightness_mixture) {
    if (!initial_soc_active %in% 0:1 || !initial_dio_active %in% 0:1) {
      stop("Initial dummy states must be zero or one.")
    }
    if (!is.finite(tau_active) || tau_active <= 0 ||
        !is.finite(tau_inactive) || tau_inactive <= tau_active) {
      stop("Require 0 < tau_active < tau_inactive.")
    }
    if (!is.finite(soc_variance_multiplier) ||
        soc_variance_multiplier <= 0) {
      stop("soc_variance_multiplier must be positive.")
    }
    if (!is.finite(pi_soc_active) || pi_soc_active <= 0 ||
        pi_soc_active >= 1 || !is.finite(pi_dio_active) ||
        pi_dio_active <= 0 || pi_dio_active >= 1) {
      stop("The active-state probabilities must lie strictly between 0 and 1.")
    }
    if (learn_dummy_mixture_probability &&
        (pi_soc_alpha <= 0 || pi_soc_beta <= 0 ||
         pi_dio_alpha <= 0 || pi_dio_beta <= 0)) {
      stop("Beta prior parameters for mixture probabilities must be positive.")
    }
  }

  if (length(s_vec) != p) {
    stop("length(s_vec) must equal p.")
  }
  if (is.null(A_prob)) A_prob <- rep(1 / length(A_grid), length(A_grid))
  if (is.null(D_prob)) D_prob <- rep(1 / length(D_grid), length(D_grid))
  if (length(A_prob) != length(A_grid) || any(A_prob < 0) || sum(A_prob) <= 0) {
    stop("A_prob must be nonnegative and match A_grid.")
  }
  if (length(D_prob) != length(D_grid) || any(D_prob < 0) || sum(D_prob) <= 0) {
    stop("D_prob must be nonnegative and match D_grid.")
  }
  A_prob <- A_prob / sum(A_prob)
  D_prob <- D_prob / sum(D_prob)
  
  ## ------------------------------------------------------------
  ## Observed regression
  ## ------------------------------------------------------------
  obs <- make_observed_regression(X_train, p, W_list, s_vec, alpha_mode)
  H_o <- obs$H
  x_o <- obs$x
  
  K_gamma <- ncol(H_o)
  
  ## ------------------------------------------------------------
  ## Dummy regressions, unscaled
  ## ------------------------------------------------------------
  if (use_dummy) {
    dum <- make_dummy_observations(
      X_train, p, W_list, s_vec, alpha_mode = alpha_mode,
      reduce_soc = reduce_soc_dummy
    )
    H_soc <- dum$H_soc
    x_soc <- dum$x_soc
    H_dio <- dum$H_dio
    x_dio <- dum$x_dio
  } else {
    H_soc <- NULL
    x_soc <- NULL
    H_dio <- NULL
    x_dio <- NULL
  }
  d_soc <- if (use_soc) length(x_soc) else N^2
  d_dio <- if (use_dio) length(x_dio) else N
  if (is.null(b_soc)) b_soc <- d_soc / 2 + 2
  if (is.null(b_dio)) b_dio <- d_dio / 2 + 2
  
  idx <- make_theta_index(N, p, s_vec, alpha_mode)
  
  ## ------------------------------------------------------------
  ## Initial values
  ## ------------------------------------------------------------
  gamma <- rep(0, K_gamma)
  sigma2 <- rep(1, N)
  
  lambda_alpha2 <- 0.04
  lambda_beta2 <- 0.01
  kappa_alpha <- 1
  
  a_cur <- 2
  d_cur <- 2
  ## The submitted zero-centered Minnesota specification uses fixed b=0.
  ## Learned-center modes are retained only as optional diagnostic variants;
  ## they are not used by the simulation or empirical application.
  b_cur <- if (alpha_mode == "local") {
    rep(0, p)
  } else {
    switch(
      b_mode,
      lag = c(b, rep(0, p - 1L)),
      node = rep(b, N),
      b
    )
  }
  
  ## tau_soc and tau_dio are initial standard deviations.  Their squares are
  ## subsequently updated from conjugate inverse-gamma full conditionals.
  z_soc <- if (use_soc &&
               (dummy_tightness_mixture || dummy_precision_selection)) {
    as.integer(initial_soc_active)
  } else NA_integer_
  z_dio <- if (use_dio &&
               (dummy_tightness_mixture || dummy_precision_selection)) {
    as.integer(initial_dio_active)
  } else NA_integer_
  omega_soc <- if (use_soc && dummy_precision_selection) {
    dummy_precision_shape / dummy_precision_rate
  } else NA_real_
  omega_dio <- if (use_dio && dummy_precision_selection) {
    dummy_precision_shape / dummy_precision_rate
  } else NA_real_
  tau_soc2 <- if (use_soc && dummy_precision_selection) {
    1 / omega_soc
  } else if (use_soc && dummy_tightness_mixture) {
    soc_variance_multiplier *
      c(tau_inactive^2, tau_active^2)[z_soc + 1L]
  } else {
    tau_soc^2
  }
  tau_dio2 <- if (use_dio && dummy_precision_selection) {
    1 / omega_dio
  } else if (use_dio && dummy_tightness_mixture) {
    c(tau_inactive^2, tau_active^2)[z_dio + 1L]
  } else {
    tau_dio^2
  }
  
  keep_idx <- seq(burn + 1, n_iter, by = thin)
  n_keep <- length(keep_idx)
  
  gamma_draws <- matrix(NA_real_, n_keep, K_gamma)
  sigma_draws <- matrix(NA_real_, n_keep, N)
  
  lambda_alpha_draws <- numeric(n_keep)
  lambda_beta_draws <- numeric(n_keep)
  kappa_alpha_draws <- rep(NA_real_, n_keep)
  
  a_draws <- numeric(n_keep)
  d_draws <- numeric(n_keep)
  
  tau_soc2_draws <- numeric(n_keep)
  tau_dio2_draws <- numeric(n_keep)
  z_soc_draws <- rep(NA_integer_, n_keep)
  z_dio_draws <- rep(NA_integer_, n_keep)
  pi_soc_draws <- rep(NA_real_, n_keep)
  pi_dio_draws <- rep(NA_real_, n_keep)
  omega_soc_draws <- rep(NA_real_, n_keep)
  omega_dio_draws <- rep(NA_real_, n_keep)
  n_b <- if (alpha_mode == "local") {
    p
  } else {
    switch(b_mode, lag = p, node = N, 1L)
  }
  b_draws <- matrix(NA_real_, n_keep, n_b)
  colnames(b_draws) <- switch(
    if (alpha_mode == "local") "lag" else b_mode,
    lag = paste0("b_lag", seq_len(p)),
    node = paste0("b_node", seq_len(N)),
    "b"
  )
  
  keep_counter <- 1
  
  ## Row-node identifiers.
  ## Observed data are stacked by time block, each block has N rows.
  row_node_obs <- rep(seq_len(N), times = Tstar)
  
  ## SOC dummy has N blocks, each block has N rows.
  row_node_soc <- if (reduce_soc_dummy) {
    rep(NA_integer_, d_soc)
  } else {
    rep(seq_len(N), times = N)
  }
  
  ## DIO dummy has one N-vector.
  row_node_dio <- seq_len(N)
  
  for (iter in seq_len(n_iter)) {
    
    ## ------------------------------------------------------------
    ## Build prior for gamma
    ## ------------------------------------------------------------
    prior <- make_prior_mean_var(
      N = N,
      p = p,
      s_vec = s_vec,
      lambda_alpha2 = lambda_alpha2,
      lambda_beta2 = lambda_beta2,
      a = a_cur,
      d = d_cur,
      b = b_cur,
      center_all_lags = b_mode == "lag",
      c_mu = c_mu,
      alpha_mode = alpha_mode,
      kappa_alpha = kappa_alpha,
      b_lag = if (alpha_mode == "local") b_cur else rep(0, p)
    )
    
    m0 <- prior$m
    Q0_diag <- 1 / prior$V_diag
    
    ## ------------------------------------------------------------
    ## Step 1. Sample gamma
    ## Weighted regression representation:
    ## observed: variance sigma_i^2
    ## SOC/DIO: variance tau^2 under the main identity-covariance
    ## specification.  The legacy innovation-scaled option is retained only
    ## for reproducibility of earlier results.
    ## ------------------------------------------------------------
    
    # inv_sigma_obs <- 1 / sigma2[row_node_obs]
    # 
    # H_weighted <- H_o * sqrt(inv_sigma_obs)
    # x_weighted <- x_o * sqrt(inv_sigma_obs)
    # 
    # if (use_dummy) {
    #   inv_sigma_soc <- 1 / (tau_soc2 * sigma2[row_node_soc])
    #   inv_sigma_dio <- 1 / (tau_dio2 * sigma2[row_node_dio])
    #   
    #   H_soc_weighted <- H_soc * sqrt(inv_sigma_soc)
    #   x_soc_weighted <- x_soc * sqrt(inv_sigma_soc)
    #   
    #   H_dio_weighted <- H_dio * sqrt(inv_sigma_dio)
    #   x_dio_weighted <- x_dio * sqrt(inv_sigma_dio)
    #   
    #   H_weighted <- rbind(H_weighted, H_soc_weighted, H_dio_weighted)
    #   x_weighted <- c(x_weighted, x_soc_weighted, x_dio_weighted)
    # }
    # 
    # Q_post <- crossprod(H_weighted) + diag(Q0_diag, K_gamma)
    # h_post <- crossprod(H_weighted, x_weighted) + Q0_diag * m0
    
    w_obs <- 1 / sigma2[row_node_obs]
    
    Q_lik <- crossprod(H_o, sweep(H_o, 1, w_obs, "*"))
    h_lik <- as.vector(crossprod(H_o, x_o * w_obs))
    
    soc_is_active <- use_soc &&
      (!dummy_precision_selection || identical(z_soc, 1L))
    dio_is_active <- use_dio &&
      (!dummy_precision_selection || identical(z_dio, 1L))

    if (soc_is_active) {
      w_soc <- if (dummy_error_covariance == "identity") {
        rep(1 / tau_soc2, length(x_soc))
      } else {
        1 / (tau_soc2 * sigma2[row_node_soc])
      }
      Q_lik <- Q_lik + crossprod(H_soc, sweep(H_soc, 1, w_soc, "*"))
      h_lik <- h_lik + as.vector(crossprod(H_soc, x_soc * w_soc))
    }
    if (dio_is_active) {
      w_dio <- if (dummy_error_covariance == "identity") {
        rep(1 / tau_dio2, length(x_dio))
      } else {
        1 / (tau_dio2 * sigma2[row_node_dio])
      }
      Q_lik <- Q_lik + crossprod(H_dio, sweep(H_dio, 1, w_dio, "*"))
      h_lik <- h_lik + as.vector(crossprod(H_dio, x_dio * w_dio))
    }
    
    Q_post <- Q_lik + diag(Q0_diag, K_gamma)
    h_post <- h_lik + Q0_diag * m0
    
    Q_post <- 0.5 * (Q_post + t(Q_post))
    
    # V_post <- solve(Q_post)
    # m_post <- as.vector(V_post %*% h_post)
    # 
    # gamma <- as.vector(MASS::mvrnorm(1, mu = m_post, Sigma = V_post))
    
    Rchol <- chol(Q_post)
    
    m_post <- backsolve(
      Rchol,
      forwardsolve(t(Rchol), h_post)
    )
    
    z <- rnorm(K_gamma)
    
    gamma <- as.vector(
      m_post + backsolve(Rchol, z)
    )
    
    ## ------------------------------------------------------------
    ## Residuals
    ## ------------------------------------------------------------
    resid_o <- x_o - as.vector(H_o %*% gamma)
    E_o <- matrix(resid_o, nrow = N, byrow = FALSE)
    
    if (use_soc) {
      resid_soc <- x_soc - as.vector(H_soc %*% gamma)
      E_soc <- if (dummy_error_covariance == "innovation") {
        matrix(resid_soc, nrow = N, byrow = FALSE)
      } else NULL
    }
    if (use_dio) {
      resid_dio <- x_dio - as.vector(H_dio %*% gamma)
    }
    
    ## ------------------------------------------------------------
    ## Step 2. Sample sigma_i^2
    ## ------------------------------------------------------------
    for (i in seq_len(N)) {
      ss_i <- sum(E_o[i, ]^2)
      
      shape_i <- c_sigma + Tstar / 2
      
      if (soc_is_active && dummy_error_covariance == "innovation") {
        ss_i <- ss_i + sum(E_soc[i, ]^2) / tau_soc2
        shape_i <- shape_i + N / 2
      }
      if (dio_is_active && dummy_error_covariance == "innovation") {
        ss_i <- ss_i + resid_dio[i]^2 / tau_dio2
        shape_i <- shape_i + 1 / 2
      }
      
      sigma2[i] <- rinvgamma(
        1,
        shape = shape_i,
        scale = d_sigma + 0.5 * ss_i
      )
    }
    
    ## ------------------------------------------------------------
    ## Extract alpha and beta
    ## ------------------------------------------------------------
    alpha_list <- vector("list", p)
    beta_list <- vector("list", p)
    
    for (j in seq_len(p)) {
      alpha_list[[j]] <- gamma[idx[[j]]$alpha]
      beta_list[[j]] <- gamma[idx[[j]]$beta]
    }

    ## ------------------------------------------------------------
    ## Step 3. Optional learned-center update.  In the main specification
    ## b_mode="fixed" and b=0, so this block performs no update and every
    ## own-lag coefficient is shrunk directly toward zero.
    ## ------------------------------------------------------------
    if (alpha_mode == "local") {
      ## alpha_ij = b_j + u_ij, with
      ## b_j ~ N(0, lambda_alpha2 * j^(-a)) and
      ## u_ij ~ N(0, kappa_alpha * lambda_alpha2 * j^(-a)).
      for (j in seq_len(p)) {
        b_post_var <- lambda_alpha2 * j^(-a_cur) *
          kappa_alpha / (kappa_alpha + N)
        b_post_mean <- sum(alpha_list[[j]]) / (kappa_alpha + N)
        b_cur[[j]] <- stats::rnorm(1L, b_post_mean, sqrt(b_post_var))
      }
    } else if (b_mode == "lag") {
      for (j in seq_len(p)) {
        alpha_var_j <- lambda_alpha2 * j^(-a_cur)
        b_post_var <- 1 / (1 / b_prior_var + N / alpha_var_j)
        b_post_mean <- b_post_var * (
          b_prior_mean / b_prior_var + sum(alpha_list[[j]]) / alpha_var_j
        )
        b_cur[[j]] <- stats::rnorm(1L, b_post_mean, sqrt(b_post_var))
      }
    } else if (b_mode == "common") {
      b_post_var <- 1 / (1 / b_prior_var + N / lambda_alpha2)
      b_post_mean <- b_post_var * (
        b_prior_mean / b_prior_var + sum(alpha_list[[1]]) / lambda_alpha2
      )
      b_cur <- stats::rnorm(1L, b_post_mean, sqrt(b_post_var))
    } else if (b_mode == "node") {
      b_post_var <- 1 / (1 / b_prior_var + 1 / lambda_alpha2)
      b_post_mean <- b_post_var * (
        b_prior_mean / b_prior_var + alpha_list[[1]] / lambda_alpha2
      )
      b_cur <- stats::rnorm(N, b_post_mean, sqrt(b_post_var))
    }
    alpha_center <- function(j) {
      if (alpha_mode == "local") {
        rep(b_cur[[j]], N)
      } else if (alpha_mode == "global") {
        0
      } else if (b_mode == "lag") {
        rep(b_cur[[j]], N)
      } else if (j == 1L) {
        rep(b_cur, length.out = N)
      } else {
        rep(0, N)
      }
    }
    
    ## ------------------------------------------------------------
    ## Step 4. Sample lambda_alpha^2
    ## ------------------------------------------------------------
    R_alpha <- 0
    if (alpha_mode == "local") {
      for (j in seq_len(p)) {
        u_j <- alpha_list[[j]] - b_cur[[j]]
        R_alpha <- R_alpha + j^a_cur * (
          b_cur[[j]]^2 + sum(u_j^2) / kappa_alpha
        )
      }
      lambda_alpha2 <- rinvgamma(
        1,
        shape = a_alpha + ((N + 1) * p) / 2,
        scale = b_alpha + 0.5 * R_alpha
      )
    } else {
      n_alpha <- if (alpha_mode == "global") 1L else N
      for (j in seq_len(p)) {
        m_alpha_j <- alpha_center(j)
        R_alpha <- R_alpha +
          j^a_cur * sum((alpha_list[[j]] - m_alpha_j)^2)
      }
      lambda_alpha2 <- rinvgamma(
        1,
        shape = a_alpha + (n_alpha * p) / 2,
        scale = b_alpha + 0.5 * R_alpha
      )
    }

    ## ------------------------------------------------------------
    ## Step 4b. Sample kappa_alpha for the local-alpha hierarchy
    ## ------------------------------------------------------------
    if (alpha_mode == "local") {
      R_kappa <- 0
      for (j in seq_len(p)) {
        u_j <- alpha_list[[j]] - b_cur[[j]]
        R_kappa <- R_kappa + j^a_cur * sum(u_j^2)
      }
      kappa_alpha <- rinvgamma(
        1,
        shape = c_kappa + (N * p) / 2,
        scale = d_kappa + R_kappa / (2 * lambda_alpha2)
      )
    }
    
    ## ------------------------------------------------------------
    ## Step 5. Sample lambda_beta^2
    ## ------------------------------------------------------------
    R_beta <- 0
    
    for (j in seq_len(p)) {
      if (s_vec[j] > 0) {
        for (r in seq_len(s_vec[j])) {
          R_beta <- R_beta + j^a_cur * r^d_cur * beta_list[[j]][r]^2
        }
      }
    }
    
    lambda_beta2 <- rinvgamma(
      1,
      shape = a_beta + sum(s_vec) / 2,
      scale = b_beta + 0.5 * R_beta
    )
    
    ## ------------------------------------------------------------
    ## Step 6. Sample a from categorical grid
    ## ------------------------------------------------------------
    logw_a <- numeric(length(A_grid))
    
    for (aa in seq_along(A_grid)) {
      a_val <- A_grid[aa]
      
      R_alpha_a <- 0
      if (alpha_mode == "local") {
        for (j in seq_len(p)) {
          u_j <- alpha_list[[j]] - b_cur[[j]]
          R_alpha_a <- R_alpha_a + j^a_val * (
            b_cur[[j]]^2 + sum(u_j^2) / kappa_alpha
          )
        }
      } else {
        for (j in seq_len(p)) {
          m_alpha_j <- alpha_center(j)
          R_alpha_a <- R_alpha_a +
            j^a_val * sum((alpha_list[[j]] - m_alpha_j)^2)
        }
      }
      
      R_beta_a <- 0
      for (j in seq_len(p)) {
        if (s_vec[j] > 0) {
          for (r in seq_len(s_vec[j])) {
            R_beta_a <- R_beta_a +
              j^a_val * r^d_cur * beta_list[[j]][r]^2
          }
        }
      }
      
      det_a <- 0.5 * sum(vapply(
        seq_len(p),
        function(j) {
          n_alpha_prior <- switch(
            alpha_mode, global = 1L, local = N + 1L, independent = N
          )
          (n_alpha_prior + s_vec[j]) * a_val * log(j)
        },
        numeric(1)
      ))
      
      logw_a[aa] <-
        log(A_prob[aa]) +
        det_a -
        R_alpha_a / (2 * lambda_alpha2) -
        R_beta_a / (2 * lambda_beta2)
    }
    
    a_cur <- A_grid[sample_categorical_logprob(logw_a)]
    
    ## ------------------------------------------------------------
    ## Step 7. Sample d from categorical grid
    ## ------------------------------------------------------------
    logw_d <- numeric(length(D_grid))
    
    for (dd in seq_along(D_grid)) {
      d_val <- D_grid[dd]
      
      R_beta_d <- 0
      det_d <- 0
      
      for (j in seq_len(p)) {
        if (s_vec[j] > 0) {
          for (r in seq_len(s_vec[j])) {
            R_beta_d <- R_beta_d +
              j^a_cur * r^d_val * beta_list[[j]][r]^2
          }
          
          det_d <- det_d +
            0.5 * d_val * sum(log(seq_len(s_vec[j])))
        }
      }
      
      logw_d[dd] <-
        log(D_prob[dd]) +
        det_d -
        R_beta_d / (2 * lambda_beta2)
    }
    
    d_cur <- D_grid[sample_categorical_logprob(logw_d)]
    
    ## ------------------------------------------------------------
    ## Step 8. Sample tau_soc^2 and tau_dio^2.  Large learned values make an
    ## incompatible pseudo restriction diffuse instead of forcing it.
    ## ------------------------------------------------------------
    if (use_soc && dummy_precision_selection) {
      ## Spike-and-slab selection on the dummy precision:
      ##   delta ~ Bernoulli(pi), omega_tilde ~ Gamma(shape, rate),
      ##   omega = delta * omega_tilde.
      ## Integrating omega_tilde out of the active Gaussian pseudo-error
      ## density gives a stable collapsed Bernoulli update and avoids the
      ## nearly absorbing states of an uncollapsed indicator update.
      resid_soc <- x_soc - as.vector(H_soc %*% gamma)
      E_soc <- if (dummy_error_covariance == "innovation") {
        matrix(resid_soc, nrow = N, byrow = FALSE)
      } else NULL
      quad_soc <- if (dummy_error_covariance == "identity") {
        sum(resid_soc^2)
      } else {
        sum(vapply(seq_len(N), function(i) {
          sum(E_soc[i, ]^2) / sigma2[i]
        }, numeric(1)))
      }
      dim_soc <- d_soc
      log_marginal_soc <-
        lgamma(dummy_precision_shape + dim_soc / 2) -
        lgamma(dummy_precision_shape) +
        dummy_precision_shape * log(dummy_precision_rate) -
        (dummy_precision_shape + dim_soc / 2) *
          log(dummy_precision_rate + quad_soc / 2) -
        (dim_soc / 2) * log(2 * pi)
      z_soc <- sample_categorical_logprob(c(
        log1p(-pi_soc_active),
        log(pi_soc_active) + log_marginal_soc
      )) - 1L
      omega_soc <- if (z_soc == 1L) {
        rgamma(
          1,
          shape = dummy_precision_shape + dim_soc / 2,
          rate = dummy_precision_rate + quad_soc / 2
        )
      } else {
        rgamma(1, shape = dummy_precision_shape,
               rate = dummy_precision_rate)
      }
      tau_soc2 <- 1 / omega_soc
    } else if (use_soc && dummy_tightness_mixture) {
      ## z_soc=1 selects the active (tight) state.  SOC uses the effective
      ## variance multiplier * tau_z^2.  The reduced representation uses
      ## d_soc effective pseudo residuals instead of the original N^2 rows.
      resid_soc <- x_soc - as.vector(H_soc %*% gamma)
      E_soc <- if (dummy_error_covariance == "innovation") {
        matrix(resid_soc, nrow = N, byrow = FALSE)
      } else NULL
      quad_soc <- if (dummy_error_covariance == "identity") {
        sum(resid_soc^2)
      } else {
        sum(vapply(seq_len(N), function(i) {
          sum(E_soc[i, ]^2) / sigma2[i]
        }, numeric(1)))
      }
      var_soc <- soc_variance_multiplier *
        c(tau_inactive^2, tau_active^2)
      logw_soc <- c(
        log1p(-pi_soc_active), log(pi_soc_active)
      ) - (d_soc / 2) * log(var_soc) - quad_soc / (2 * var_soc)
      z_soc <- sample_categorical_logprob(logw_soc) - 1L
      tau_soc2 <- var_soc[z_soc + 1L]
      if (learn_dummy_mixture_probability) {
        pi_soc_active <- rbeta(
          1, pi_soc_alpha + z_soc, pi_soc_beta + 1L - z_soc
        )
      }
    } else if (use_soc && learn_dummy_tightness) {
      ## Recompute dummy residuals after possible sigma update.
      resid_soc <- x_soc - as.vector(H_soc %*% gamma)
      E_soc <- if (dummy_error_covariance == "innovation") {
        matrix(resid_soc, nrow = N, byrow = FALSE)
      } else NULL
      
      quad_soc <- if (dummy_error_covariance == "identity") {
        sum(resid_soc^2)
      } else {
        sum(vapply(seq_len(N), function(i) {
          sum(E_soc[i, ]^2) / sigma2[i]
        }, numeric(1)))
      }
      
      tau_soc2 <- rinvgamma(
        1,
        shape = a_soc + d_soc / 2,
        scale = b_soc + 0.5 * quad_soc
      )
      
    }
    if (use_dio && dummy_precision_selection) {
      resid_dio <- x_dio - as.vector(H_dio %*% gamma)
      quad_dio <- if (dummy_error_covariance == "identity") {
        sum(resid_dio^2)
      } else {
        sum(resid_dio^2 / sigma2)
      }
      dim_dio <- N
      log_marginal_dio <-
        lgamma(dummy_precision_shape + dim_dio / 2) -
        lgamma(dummy_precision_shape) +
        dummy_precision_shape * log(dummy_precision_rate) -
        (dummy_precision_shape + dim_dio / 2) *
          log(dummy_precision_rate + quad_dio / 2) -
        (dim_dio / 2) * log(2 * pi)
      z_dio <- sample_categorical_logprob(c(
        log1p(-pi_dio_active),
        log(pi_dio_active) + log_marginal_dio
      )) - 1L
      omega_dio <- if (z_dio == 1L) {
        rgamma(
          1,
          shape = dummy_precision_shape + dim_dio / 2,
          rate = dummy_precision_rate + quad_dio / 2
        )
      } else {
        rgamma(1, shape = dummy_precision_shape,
               rate = dummy_precision_rate)
      }
      tau_dio2 <- 1 / omega_dio
    } else if (use_dio && dummy_tightness_mixture) {
      resid_dio <- x_dio - as.vector(H_dio %*% gamma)
      quad_dio <- if (dummy_error_covariance == "identity") {
        sum(resid_dio^2)
      } else {
        sum(resid_dio^2 / sigma2)
      }
      var_dio <- c(tau_inactive^2, tau_active^2)
      logw_dio <- c(
        log1p(-pi_dio_active), log(pi_dio_active)
      ) - (N / 2) * log(var_dio) - quad_dio / (2 * var_dio)
      z_dio <- sample_categorical_logprob(logw_dio) - 1L
      tau_dio2 <- var_dio[z_dio + 1L]
      if (learn_dummy_mixture_probability) {
        pi_dio_active <- rbeta(
          1, pi_dio_alpha + z_dio, pi_dio_beta + 1L - z_dio
        )
      }
    } else if (use_dio && learn_dummy_tightness) {
      resid_dio <- x_dio - as.vector(H_dio %*% gamma)
      quad_dio <- if (dummy_error_covariance == "identity") {
        sum(resid_dio^2)
      } else {
        sum(resid_dio^2 / sigma2)
      }
      
      tau_dio2 <- rinvgamma(
        1,
        shape = a_dio + N / 2,
        scale = b_dio + 0.5 * quad_dio
      )
    }
    
    ## ------------------------------------------------------------
    ## Save draws
    ## ------------------------------------------------------------
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      gamma_draws[keep_counter, ] <- gamma
      sigma_draws[keep_counter, ] <- sigma2
      
      lambda_alpha_draws[keep_counter] <- lambda_alpha2
      lambda_beta_draws[keep_counter] <- lambda_beta2
      kappa_alpha_draws[keep_counter] <- if (alpha_mode == "local") {
        kappa_alpha
      } else {
        NA_real_
      }
      
      a_draws[keep_counter] <- a_cur
      d_draws[keep_counter] <- d_cur
      
      tau_soc2_draws[keep_counter] <- if (use_soc) tau_soc2 else NA_real_
      tau_dio2_draws[keep_counter] <- if (use_dio) tau_dio2 else NA_real_
      z_soc_draws[keep_counter] <- z_soc
      z_dio_draws[keep_counter] <- z_dio
      pi_soc_draws[keep_counter] <- if (
        use_soc && dummy_tightness_mixture
      ) pi_soc_active else NA_real_
      pi_dio_draws[keep_counter] <- if (
        use_dio && dummy_tightness_mixture
      ) pi_dio_active else NA_real_
      omega_soc_draws[keep_counter] <- omega_soc
      omega_dio_draws[keep_counter] <- omega_dio
      b_draws[keep_counter, ] <- b_cur
      
      keep_counter <- keep_counter + 1
    }
    
    if (verbose && iter %% 500 == 0) {
      cat(
        "BGNAR iter =", iter,
        "| a =", a_cur,
        "| d =", d_cur,
        "| lambda_alpha =", round(sqrt(lambda_alpha2), 3),
        "| lambda_beta =", round(sqrt(lambda_beta2), 3)
      )
      if (alpha_mode == "local") {
        cat("| kappa_alpha =", round(kappa_alpha, 3))
      }
      
      if (use_soc) {
        cat("| tau_soc =", round(sqrt(tau_soc2), 3))
      }
      if (use_dio) {
        cat(
          "| tau_dio =", round(sqrt(tau_dio2), 3)
        )
      }
      
      cat("\n")
    }
  }
  
  structure(
    list(
      gamma_draws = gamma_draws,
      sigma_draws = sigma_draws,
      
      lambda_alpha_draws = lambda_alpha_draws,
      lambda_beta_draws = lambda_beta_draws,
      kappa_alpha_draws = kappa_alpha_draws,
      
      a_draws = a_draws,
      d_draws = d_draws,
      
      tau_soc2_draws = tau_soc2_draws,
      tau_dio2_draws = tau_dio2_draws,
      z_soc_draws = z_soc_draws,
      z_dio_draws = z_dio_draws,
      pi_soc_draws = pi_soc_draws,
      pi_dio_draws = pi_dio_draws,
      omega_soc_draws = omega_soc_draws,
      omega_dio_draws = omega_dio_draws,
      b_draws = b_draws,
      
      idx = idx,
      W_list = W_list,
      s_vec = s_vec,
      p = p,
      N = N,
      alpha_mode = alpha_mode,
      A_grid = A_grid,
      D_grid = D_grid,
      A_prob = A_prob,
      D_prob = D_prob,
      
    use_dummy = use_dummy,
    use_soc = use_soc,
    use_dio = use_dio,
    dummy_error_covariance = dummy_error_covariance,
    reduce_soc_dummy = reduce_soc_dummy,
    dummy_dimensions = c(SOC = d_soc, DIO = d_dio),
    learn_dummy_tightness = learn_dummy_tightness,
    dummy_precision_selection = dummy_precision_selection,
    dummy_tightness_mixture = dummy_tightness_mixture,
      
      tau_soc = sqrt(tau_soc2),
      tau_dio = sqrt(tau_dio2),
      
      a_soc = a_soc,
      b_soc = b_soc,
      a_dio = a_dio,
      b_dio = b_dio,
      prior_specification = paste0(
        "bvar_inspired_", alpha_mode, "_alpha_",
        if (dummy_precision_selection) {
          "bernoulli_gamma_dummy_precision"
        } else if (dummy_tightness_mixture) {
          "bernoulli_dummy_tightness"
        } else if (learn_dummy_tightness) {
          "adaptive_dummy_tightness"
        } else {
          "fixed_dummy_tightness"
        }
      ),
      hyperparameters = list(
        b = b, b_mode = b_mode,
        b_prior_mean = b_prior_mean, b_prior_var = b_prior_var,
        alpha_mode = alpha_mode,
        c_kappa = c_kappa, d_kappa = d_kappa,
        c_mu = c_mu,
        c_sigma = c_sigma, d_sigma = d_sigma,
        a_alpha = a_alpha, b_alpha = b_alpha,
        a_beta = a_beta, b_beta = b_beta,
        tau_soc = tau_soc, tau_dio = tau_dio,
        dummy_error_covariance = dummy_error_covariance,
        reduce_soc_dummy = reduce_soc_dummy,
        learn_dummy_tightness = learn_dummy_tightness,
        dummy_precision_selection = dummy_precision_selection,
        dummy_precision_shape = dummy_precision_shape,
        dummy_precision_rate = dummy_precision_rate,
        dummy_tightness_mixture = dummy_tightness_mixture,
        tau_active = tau_active, tau_inactive = tau_inactive,
        pi_soc_active = pi_soc_active, pi_dio_active = pi_dio_active,
        initial_soc_active = initial_soc_active,
        initial_dio_active = initial_dio_active,
        learn_dummy_mixture_probability = learn_dummy_mixture_probability,
        pi_soc_alpha = pi_soc_alpha, pi_soc_beta = pi_soc_beta,
        pi_dio_alpha = pi_dio_alpha, pi_dio_beta = pi_dio_beta,
        soc_variance_multiplier = soc_variance_multiplier,
        a_soc = a_soc, b_soc = b_soc,
        a_dio = a_dio, b_dio = b_dio,
        A_grid = A_grid, A_prob = A_prob,
        D_grid = D_grid, D_prob = D_prob
      )
    ),
    class = "bgnar_fit"
  )
}

predict_bgnar_one_step <- function(fit, X_all, train_end, horizon,
                                   simulate = TRUE) {
  N <- fit$N
  p <- fit$p
  W_list <- fit$W_list
  s_vec <- fit$s_vec
  
  G <- fit$gamma_draws
  S <- fit$sigma_draws
  n_draw <- nrow(G)
  
  pred_mean <- array(NA_real_, c(horizon, N, n_draw))
  pred_y <- array(NA_real_, c(horizon, N, n_draw))
  
  for (h in seq_len(horizon)) {
    t <- train_end + h
    
    x_lags <- lapply(seq_len(p), function(j) {
      X_all[t - j, ]
    })
    
    H_t <- build_H_one(x_lags, W_list, s_vec, fit$alpha_mode)
    
    for (m in seq_len(n_draw)) {
      mu_m <- as.vector(H_t %*% G[m, ])
      pred_mean[h, , m] <- mu_m
      
      if (simulate) {
        pred_y[h, , m] <- stats::rnorm(N, mean = mu_m, sd = sqrt(S[m, ]))
      } else {
        pred_y[h, , m] <- mu_m
      }
    }
  }
  
  list(pred_mean = pred_mean, pred_y = pred_y)
}

predict_bgnar_recursive <- function(fit, X_history, horizon,
                                    simulate = TRUE) {
  N <- fit$N
  p <- fit$p
  W_list <- fit$W_list
  s_vec <- fit$s_vec
  
  G <- fit$gamma_draws
  S <- fit$sigma_draws
  n_draw <- nrow(G)
  
  X_history <- as.matrix(X_history)
  
  pred_mean <- array(NA_real_, c(horizon, N, n_draw))
  pred_y <- array(NA_real_, c(horizon, N, n_draw))
  
  for (m in seq_len(n_draw)) {
    ## Keep two paths.  The deterministic path gives the exact conditional
    ## mean for a parameter draw and is used for squared-error point forecasts.
    ## The stochastic path propagates future innovations and is retained for
    ## posterior predictive intervals.
    X_mean_path <- X_history
    X_y_path <- X_history
    
    for (h in seq_len(horizon)) {
      current_T <- nrow(X_mean_path)
      
      mean_lags <- lapply(seq_len(p), function(j) {
        X_mean_path[current_T - j + 1, ]
      })
      H_mean <- build_H_one(mean_lags, W_list, s_vec, fit$alpha_mode)
      mu_mean <- as.vector(H_mean %*% G[m, ])
      
      pred_mean[h, , m] <- mu_mean
      X_mean_path <- rbind(X_mean_path, mu_mean)
      
      if (simulate) {
        y_lags <- lapply(seq_len(p), function(j) {
          X_y_path[current_T - j + 1, ]
        })
        H_y <- build_H_one(y_lags, W_list, s_vec, fit$alpha_mode)
        mu_y <- as.vector(H_y %*% G[m, ])
        y_t <- stats::rnorm(N, mean = mu_y, sd = sqrt(S[m, ]))
      } else {
        y_t <- mu_mean
      }
      
      pred_y[h, , m] <- y_t
      X_y_path <- rbind(X_y_path, y_t)
    }
  }
  
  list(pred_mean = pred_mean, pred_y = pred_y)
}



fitted_bgnar <- function(fit, X_train, simulate = TRUE) {
  ## In-sample one-step fitted values for t = p+1, ..., T_train.
  ## Uses the same posterior draws of gamma as prediction functions.
  X_train <- as.matrix(X_train)
  N <- fit$N
  p <- fit$p
  W_list <- fit$W_list
  s_vec <- fit$s_vec
  G <- fit$gamma_draws
  S <- fit$sigma_draws
  n_draw <- nrow(G)
  T_fit <- nrow(X_train) - p

  fit_mean <- array(NA_real_, c(T_fit, N, n_draw))
  fit_y <- array(NA_real_, c(T_fit, N, n_draw))

  for (tt in seq_len(T_fit)) {
    t <- p + tt
    x_lags <- lapply(seq_len(p), function(j) X_train[t - j, ])
    H_t <- build_H_one(x_lags, W_list, s_vec, fit$alpha_mode)

    for (m in seq_len(n_draw)) {
      mu_m <- as.vector(H_t %*% G[m, ])
      fit_mean[tt, , m] <- mu_m

      if (simulate) {
        fit_y[tt, , m] <- stats::rnorm(N, mean = mu_m, sd = sqrt(S[m, ]))
      } else {
        fit_y[tt, , m] <- mu_m
      }
    }
  }

  list(fit_mean = fit_mean, fit_y = fit_y)
}

summarize_bgnar_fitted <- function(fitted_obj, X_train, p,
                                   point_source = c("fit_mean", "fit_y")) {
  point_source <- match.arg(point_source)

  T_fit <- dim(fitted_obj$fit_mean)[1]
  N <- dim(fitted_obj$fit_mean)[2]
  X_fit <- as.matrix(X_train)[(p + 1):nrow(X_train), , drop = FALSE]

  med <- lo <- hi <- matrix(NA_real_, T_fit, N)

  for (tt in seq_len(T_fit)) {
    for (i in seq_len(N)) {
      vals_y <- fitted_obj$fit_y[tt, i, ]
      vals_point <- if (point_source == "fit_mean") {
        fitted_obj$fit_mean[tt, i, ]
      } else {
        fitted_obj$fit_y[tt, i, ]
      }

      med[tt, i] <- median(vals_point)
      lo[tt, i] <- stats::quantile(vals_y, 0.025)
      hi[tt, i] <- stats::quantile(vals_y, 0.975)
    }
  }

  metrics <- compute_metrics(
    med,
    X_fit[seq_len(T_fit), , drop = FALSE],
    lo,
    hi
  )

  c(list(fit_med = med, fit_lo = lo, fit_hi = hi), metrics)
}

summarize_bgnar_prediction <- function(pred_obj, X_test,
                                       point_source = c("pred_y", "pred_mean")) {
  point_source <- match.arg(point_source)
  
  horizon <- dim(pred_obj$pred_y)[1]
  N <- dim(pred_obj$pred_y)[2]
  
  med <- lo <- hi <- matrix(NA_real_, horizon, N)
  
  for (h in seq_len(horizon)) {
    for (i in seq_len(N)) {
      vals_y <- pred_obj$pred_y[h, i, ]
      vals_point <- if (point_source == "pred_mean") {
        pred_obj$pred_mean[h, i, ]
      } else {
        pred_obj$pred_y[h, i, ]
      }
      
      ## The posterior mean minimizes expected squared forecast loss.  Keep
      ## predictive quantiles from stochastic paths for interval evaluation.
      med[h, i] <- if (point_source == "pred_mean") {
        mean(vals_point)
      } else {
        median(vals_point)
      }
      lo[h, i] <- stats::quantile(vals_y, 0.025)
      hi[h, i] <- stats::quantile(vals_y, 0.975)
    }
  }
  
  metrics <- compute_metrics(
    med,
    X_test[seq_len(horizon), , drop = FALSE],
    lo,
    hi
  )
  
  c(list(pred_med = med, pred_lo = lo, pred_hi = hi), metrics)
}

summarize_bgnar_parameters <- function(fit) {
  N <- fit$N
  p <- fit$p
  s_vec <- fit$s_vec
  idx <- fit$idx
  G <- fit$gamma_draws
  
  M_draw <- G[, seq_len(N), drop = FALSE]
  mu_summary <- data.frame(
    node = seq_len(N),
    mean = colMeans(M_draw),
    q025 = apply(M_draw, 2, stats::quantile, 0.025),
    q500 = apply(M_draw, 2, stats::quantile, 0.5),
    q975 = apply(M_draw, 2, stats::quantile, 0.975)
  )

  alpha_summary <- data.frame()
  beta_summary <- data.frame()
  
  for (j in seq_len(p)) {
    A_draw <- G[, idx[[j]]$alpha, drop = FALSE]
    if (identical(fit$alpha_mode, "global")) {
      A_draw <- matrix(
        rep(A_draw[, 1L], N), nrow = nrow(A_draw), ncol = N
      )
    }
    
    alpha_summary <- rbind(
      alpha_summary,
      data.frame(
        lag = j,
        node = seq_len(N),
        mean = colMeans(A_draw),
        q025 = apply(A_draw, 2, stats::quantile, 0.025),
        q500 = apply(A_draw, 2, stats::quantile, 0.5),
        q975 = apply(A_draw, 2, stats::quantile, 0.975)
      )
    )
    
    if (s_vec[j] > 0L) {
      B_draw <- G[, idx[[j]]$beta, drop = FALSE]
      beta_summary <- rbind(
        beta_summary,
        data.frame(
          lag = j,
          neighbor_order = seq_len(s_vec[j]),
          mean = colMeans(B_draw),
          q025 = apply(B_draw, 2, stats::quantile, 0.025),
          q500 = apply(B_draw, 2, stats::quantile, 0.5),
          q975 = apply(B_draw, 2, stats::quantile, 0.975)
        )
      )
    }
  }
  
  tau_soc_summary <- if (!is.null(fit$tau_soc2_draws) &&
                         !all(is.na(fit$tau_soc2_draws))) {
    summary(sqrt(fit$tau_soc2_draws))
  } else {
    NA
  }
  
  tau_dio_summary <- if (!is.null(fit$tau_dio2_draws) &&
                         !all(is.na(fit$tau_dio2_draws))) {
    summary(sqrt(fit$tau_dio2_draws))
  } else {
    NA
  }
  
  list(
    mu = mu_summary,
    alpha = alpha_summary,
    beta = beta_summary,
    
    a_table = prop.table(table(fit$a_draws)),
    d_table = prop.table(table(fit$d_draws)),
    
    lambda_alpha = summary(sqrt(fit$lambda_alpha_draws)),
    lambda_beta = summary(sqrt(fit$lambda_beta_draws)),
    kappa_alpha = if (!is.null(fit$kappa_alpha_draws) &&
                      !all(is.na(fit$kappa_alpha_draws))) {
      summary(fit$kappa_alpha_draws)
    } else {
      NA
    },
    
    tau_soc = tau_soc_summary,
    tau_dio = tau_dio_summary
  )
}

plot_bgnar_beta <- function(beta_summary, beta_true_list = NULL) {
  
  plot_df <- beta_summary
  
  ## ------------------------------------------------------------
  ## Add true beta values if provided
  ## ------------------------------------------------------------
  if (!is.null(beta_true_list)) {
    beta_true_df <- data.frame(
      lag = rep(seq_along(beta_true_list), times = lengths(beta_true_list)),
      neighbor_order = unlist(lapply(beta_true_list, seq_along)),
      true = unlist(beta_true_list)
    )
    
    plot_df <- merge(
      plot_df,
      beta_true_df,
      by = c("lag", "neighbor_order"),
      all.x = TRUE
    )
  } else {
    plot_df$true <- NA_real_
  }
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = factor(neighbor_order),
      y = mean,
      ymin = q025,
      ymax = q975
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_pointrange() +
    ggplot2::geom_point(
      ggplot2::aes(y = true),
      shape = 4,
      size = 3,
      stroke = 1.1,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~ lag, scales = "free_y") +
    ggplot2::labs(
      title = "BGNAR posterior beta[j,r]",
      subtitle = "Dots: posterior mean and 95% CI; crosses: true values",
      x = "Neighborhood order r",
      y = "Posterior mean and 95% CI"
    ) +
    ggplot2::theme_minimal()
}


plot_bgnar_alpha <- function(alpha_summary, alpha_true_mat = NULL) {
  
  plot_df <- alpha_summary
  
  ## ------------------------------------------------------------
  ## Add true alpha values if provided
  ## alpha_true_mat: N x p matrix
  ## ------------------------------------------------------------
  if (!is.null(alpha_true_mat)) {
    N <- nrow(alpha_true_mat)
    p <- ncol(alpha_true_mat)
    
    alpha_true_df <- data.frame(
      lag = rep(seq_len(p), each = N),
      node = rep(seq_len(N), times = p),
      true = as.vector(alpha_true_mat)
    )
    
    plot_df <- merge(
      plot_df,
      alpha_true_df,
      by = c("lag", "node"),
      all.x = TRUE
    )
  } else {
    plot_df$true <- NA_real_
  }
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = factor(node),
      y = mean,
      ymin = q025,
      ymax = q975
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_pointrange() +
    ggplot2::geom_point(
      ggplot2::aes(y = true),
      shape = 4,
      size = 2.7,
      stroke = 1.0,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~ lag, scales = "free_y") +
    ggplot2::labs(
      title = "BGNAR posterior alpha[i,j]",
      subtitle = "Dots: posterior mean and 95% CI; crosses: true values",
      x = "Node",
      y = "Posterior mean and 95% CI"
    ) +
    ggplot2::theme_minimal()
}
