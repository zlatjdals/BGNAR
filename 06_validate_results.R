#!/usr/bin/env Rscript

############################################################
## 06_validate_results.R
## Check grid completeness and per-method failures.
############################################################

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else normalizePath(Sys.getenv("BGNAR_SCRIPT_DIR", "."))
source(file.path(script_dir, "03_cluster_functions.R"))

result_root <- normalizePath(
  Sys.getenv("RESULT_ROOT", file.path(script_dir, "results")), mustWork = FALSE
)
checkpoint_dir <- file.path(result_root, "checkpoints")
summary_dir <- file.path(result_root, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

dgps <- toupper(csv_tokens(Sys.getenv("DGP_LIST", "M11,M12,M13,M2")))
structures <- toupper(csv_tokens(Sys.getenv("GRAPH_LIST", "ER,SBM,SWN")))
train_grid <- as.integer(csv_tokens(Sys.getenv("TRAIN_GRID", "20,40,80,160")))
replications <- env_int("REPLICATIONS", 100L)
methods <- csv_tokens(Sys.getenv(
  "METHODS", "bgnar_soc_dio_adaptive,gnar,bvar"
))
strict <- env_flag("STRICT_VALIDATION", TRUE)

expected <- expand.grid(
  scenario = dgps, structure = structures, train_T = train_grid,
  replication = seq_len(replications), stringsAsFactors = FALSE
)
expected$job_id <- sprintf(
  "%s_%s_T%03d_r%03d", expected$scenario, expected$structure,
  expected$train_T, expected$replication
)
expected$checkpoint <- file.path(checkpoint_dir, paste0(expected$job_id, ".rds"))
expected$exists <- file.exists(expected$checkpoint)

issues <- list()
if (any(!expected$exists)) {
  issues[[length(issues) + 1L]] <- data.frame(
    job_id = expected$job_id[!expected$exists],
    issue = "missing checkpoint", detail = expected$checkpoint[!expected$exists]
  )
}

for (path in expected$checkpoint[expected$exists]) {
  z <- tryCatch(readRDS(path), error = function(e) e)
  id <- sub("\\.rds$", "", basename(path))
  if (inherits(z, "error")) {
    issues[[length(issues) + 1L]] <- data.frame(
      job_id = id, issue = "unreadable checkpoint", detail = conditionMessage(z)
    )
    next
  }
  missing_methods <- setdiff(methods, names(z$methods))
  if (length(missing_methods)) {
    issues[[length(issues) + 1L]] <- data.frame(
      job_id = id, issue = "missing successful method",
      detail = paste(missing_methods, collapse = ",")
    )
  }
  if (!is.null(z$failures) && nrow(z$failures)) {
    issues[[length(issues) + 1L]] <- data.frame(
      job_id = id, issue = paste0("method failure: ", z$failures$method),
      detail = z$failures$error
    )
  }
  for (method in intersect(methods, names(z$methods))) {
    rmse <- z$methods[[method]]$forecast$overall$rmse
    if (!length(rmse) || !is.finite(rmse)) {
      issues[[length(issues) + 1L]] <- data.frame(
        job_id = id, issue = paste0("non-finite RMSE: ", method), detail = ""
      )
    }
  }
}

issue_df <- if (length(issues)) do.call(rbind, issues) else data.frame(
  job_id = character(), issue = character(), detail = character()
)
write.csv(expected, file.path(summary_dir, "expected_grid.csv"), row.names = FALSE)
write.csv(issue_df, file.path(summary_dir, "validation_issues.csv"), row.names = FALSE)

report <- c(
  paste0("expected_checkpoints=", nrow(expected)),
  paste0("present_checkpoints=", sum(expected$exists)),
  paste0("requested_methods=", paste(methods, collapse = ",")),
  paste0("validation_issues=", nrow(issue_df)),
  paste0("strict_validation=", strict)
)
writeLines(report, file.path(summary_dir, "validation_report.txt"))
message(paste(report, collapse = "\n"))
if (strict && nrow(issue_df)) quit(save = "no", status = 2L)
