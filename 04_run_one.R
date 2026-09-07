#!/usr/bin/env Rscript

############################################################
## 04_run_one.R
## One deterministic job: DGP x graph x T x replication.
############################################################

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))

.libPaths(unique(c(
  file.path(script_dir, ".Rlib"),
  file.path(dirname(script_dir), ".Rlib"),
  .libPaths()
)))

source(file.path(script_dir, "00_utils.R"))
source(file.path(script_dir, "01_bgnar_model.R"))
source(file.path(script_dir, "02_compare_methods.R"))
source(file.path(script_dir, "03_cluster_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript 04_run_one.R <DGP> <ER|SBM|SWN> <train_T> <replication>"
  )
}

scenario <- toupper(args[[1L]])
structure <- toupper(args[[2L]])
train_end <- suppressWarnings(as.integer(args[[3L]]))
replication <- suppressWarnings(as.integer(args[[4L]]))
if (!scenario %in% names(simulation_scenarios())) stop("Unknown DGP: ", scenario)
if (!structure %in% c("ER", "SBM", "SWN")) stop("Unknown graph: ", structure)
if (!is.finite(train_end) || train_end < 6L) stop("train_T must be at least 6.")
if (!is.finite(replication) || replication < 1L) stop("replication must be positive.")

cfg <- read_cluster_config()
if (cfg$mcmc_burn >= cfg$mcmc_iter) stop("MCMC_BURN must be smaller than MCMC_ITER.")
if (cfg$bvar_burn >= cfg$bvar_draws) stop("BVAR_BURN must be smaller than BVAR_DRAWS.")

required <- c("igraph")
if ("gnar" %in% cfg$methods) required <- c(required, "GNAR")
if ("bvar" %in% cfg$methods) required <- c(required, "BVAR")
missing_packages <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")),
  mustWork = FALSE
)
checkpoint_dir <- file.path(result_root, "checkpoints")
full_fit_dir <- file.path(result_root, "full_fits")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
if (cfg$save_full_fits) dir.create(full_fit_dir, recursive = TRUE, showWarnings = FALSE)

job_id <- sprintf("%s_%s_T%03d_r%03d", scenario, structure, train_end, replication)
checkpoint <- file.path(checkpoint_dir, paste0(job_id, ".rds"))
if (file.exists(checkpoint) && !cfg$overwrite) {
  message("Checkpoint exists; skipping: ", checkpoint)
  quit(save = "no", status = 0L)
}

started <- Sys.time()
message("Starting ", job_id, " with methods: ", paste(cfg$methods, collapse = ", "))
dat <- simulate_cluster_dgp(scenario, structure, train_end, replication, cfg)

compact_methods <- list()
failures <- list()
for (method in cfg$methods) {
  message("  fitting ", method)
  seed <- method_seed(cfg, scenario, structure, train_end, replication, method)
  tick <- proc.time()[[3L]]
  fit_result <- tryCatch(
    fit_requested_method(method, dat, train_end, cfg, seed),
    error = function(e) e
  )
  elapsed <- proc.time()[[3L]] - tick
  if (inherits(fit_result, "error") || is.null(fit_result)) {
    msg <- if (inherits(fit_result, "error")) conditionMessage(fit_result) else "NULL fit"
    failures[[length(failures) + 1L]] <- data.frame(
      method = method, error = msg, elapsed_seconds = elapsed
    )
    message("  FAILED ", method, ": ", msg)
    next
  }
  compact_methods[[method]] <- compact_method_result(
    fit_result, dat, train_end, method, cfg, elapsed
  )
  if (cfg$save_full_fits) {
    full_path <- file.path(full_fit_dir, paste0(job_id, "_", method, ".rds"))
    tmp_full <- tempfile(pattern = paste0(job_id, "_"), tmpdir = full_fit_dir)
    saveRDS(fit_result, tmp_full, compress = "xz")
    if (!file.rename(tmp_full, full_path)) stop("Could not finalize ", full_path)
  }
  rm(fit_result)
  invisible(gc(FALSE))
}

failure_df <- if (length(failures)) do.call(rbind, failures) else data.frame(
  method = character(), error = character(), elapsed_seconds = numeric()
)

result <- list(
  metadata = data.frame(
    job_id = job_id, scenario = scenario, structure = structure,
    train_T = train_end, replication = replication,
    N = cfg$N, horizon = cfg$horizon,
    dgp_label = dat$dgp_label, spectral_radius = dat$rho,
    graph_density = dat$graph_density, mean_degree = dat$mean_degree,
    graph_seed = dat$graph_seed, data_seed = dat$data_seed,
    intercept_mode = cfg$intercept_mode,
    intercept_value = cfg$intercept_value,
    gnar_mode = cfg$gnar_mode,
    innovation_sd = cfg$innovation_sd,
    started = format(started, "%Y-%m-%d %H:%M:%S %z"),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  ),
  config = list(
    methods = cfg$methods, gnar_mode = cfg$gnar_mode,
    bgnar_p = cfg$p_fit, bgnar_s = cfg$s_fit,
    mcmc_iter = cfg$mcmc_iter, mcmc_burn = cfg$mcmc_burn,
    reduce_soc = cfg$reduce_soc,
    fixed_tau2 = cfg$fixed_tau2, fixed_tau_scale = cfg$fixed_tau_scale,
    adaptive_tau_a = cfg$adaptive_a, adaptive_tau_b = cfg$adaptive_b,
    alpha_mode = "hierarchical local alpha",
    alpha_hyperprior = "kappa_alpha ~ IG(3,2)"
  ),
  truth = list(
    intercept = dat$intercept_true,
    alpha = dat$alpha_true, beta = dat$beta_true,
    adjacency = dat$adjacency
  ),
  methods = compact_methods,
  failures = failure_df,
  package_versions = vapply(
    unique(required), function(pkg) as.character(utils::packageVersion(pkg)),
    character(1L)
  )
)

tmp_checkpoint <- tempfile(pattern = paste0(job_id, "_"), tmpdir = checkpoint_dir)
saveRDS(result, tmp_checkpoint, compress = "xz")
if (!file.rename(tmp_checkpoint, checkpoint)) {
  stop("Could not atomically finalize checkpoint: ", checkpoint)
}
message("Completed ", job_id, ": ", checkpoint)

if (nrow(failure_df)) quit(save = "no", status = 2L)
