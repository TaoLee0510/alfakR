args <- commandArgs(trailingOnly = TRUE)
command_args <- commandArgs()
script_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE)
} else {
  NA_character_
}

default_repo_dir <- if (is.na(script_path)) {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE)
}

usage <- paste(
  "Usage:",
  "  Rscript benchmark/scr/search_weighted_params.R <alfak_input_rds> [n_cores] [repo_dir]",
  "",
  "Outputs are written under:",
  "  <repo_dir>/benchmark/results/calibration/<sample_name>_weighted_param_search/",
  sep = "\n"
)

if (!length(args) || args[[1L]] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(save = "no", status = if (length(args)) 0L else 1L)
}

input_rds <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (!file.exists(input_rds)) {
  stop("Input Rds not found: ", input_rds)
}

requested_n_cores <- if (length(args) >= 2L) suppressWarnings(as.integer(args[[2L]])) else 1L
if (!is.finite(requested_n_cores) || is.na(requested_n_cores) || requested_n_cores < 1L) {
  stop("`n_cores` must be a positive integer.")
}

repo_dir <- if (length(args) >= 3L) {
  normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
} else {
  default_repo_dir
}
if (!dir.exists(repo_dir)) {
  stop("Repo dir not found: ", repo_dir)
}

r_libs_user <- Sys.getenv("R_LIBS_USER", unset = "")

if (nzchar(r_libs_user)) {
  dir.create(r_libs_user, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(r_libs_user, .libPaths())))
}

available_n_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(available_n_cores) || is.na(available_n_cores)) {
  available_n_cores <- parallel::detectCores()
}
if (!is.finite(available_n_cores) || is.na(available_n_cores)) {
  available_n_cores <- requested_n_cores
}
effective_n_cores <- as.integer(max(1L, min(requested_n_cores, available_n_cores)))

sample_id <- sub("\\.Rds$", "", basename(input_rds), ignore.case = TRUE)
sample_key <- gsub("[^A-Za-z0-9._-]+", "_", sample_id)

calibration_root <- file.path(repo_dir, "benchmark", "results", "calibration")
output_root <- file.path(calibration_root, sprintf("%s_weighted_param_search", sample_key))

fixed_cfg <- list(
  patient_id = sample_id,
  minobs = 20L,
  nboot_stage1 = 20L,
  nboot_stage2 = 45L,
  nboot_stage3 = 45L,
  n0 = 1e5,
  nb = 1e7,
  pm = 5e-05,
  correct_efflux = TRUE,
  nn_prior = "empirical_censored_weighted",
  nn_prior_fit_subset = "hybrid",
  nn_prior_grid_n = 161L,
  benchmark_seed = 31415926L,
  stage3_seeds = c(31415926L, 31415927L, 31415928L)
)

coarse_birth_levels <- c(0.02, 0.10, 0.50)
coarse_cap_levels <- c(NA_real_, 0.05, 0.10, 0.15)
coarse_quantile_levels <- c(0.10, 0.30, 0.50)
coarse_scale_levels <- 0.50

fine_birth_levels <- c(0.02, 0.05, 0.10, 0.15, 0.25, 0.50)
fine_cap_levels <- c(NA_real_, 0.05, 0.10, 0.15, 0.20)
fine_quantile_levels <- c(0.10, 0.20, 0.30, 0.40, 0.50)
fine_scale_levels <- c(0.25, 0.50)

diploid_state <- paste(rep(2, 22), collapse = ".")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "runs"), recursive = TRUE, showWarnings = FALSE)

sample_output_path <- function(stem, ext) {
  file.path(output_root, sprintf("%s_%s.%s", sample_key, stem, ext))
}

load_alfakr <- function(repo_dir) {
  if (requireNamespace("alfakR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(alfakR))
    return("installed")
  }
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_dir, quiet = TRUE)
    return("source")
  }
  stop("Could not load `alfakR`. Install the package or make `pkgload` available.")
}

format_param_value <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) {
    return("adaptive")
  }
  txt <- formatC(as.numeric(x), digits = 6, format = "fg", flag = "#")
  txt <- gsub("\\+", "", txt, fixed = TRUE)
  txt <- gsub("-", "m", txt, fixed = TRUE)
  txt <- gsub("\\.", "p", txt)
  txt
}

combo_key <- function(cfg) {
  paste(
    sprintf("bw_%s", format_param_value(cfg$nn_prior_zero_birth_fallback_weight)),
    sprintf("cap_%s", format_param_value(cfg$nn_prior_zero_weight_cap_ratio)),
    sprintf("q_%s", format_param_value(cfg$nn_prior_zero_exposure_quantile)),
    sprintf("scale_%s", format_param_value(cfg$nn_prior_zero_weight_scale)),
    sep = "__"
  )
}

make_cfg <- function(birth, cap, quantile, scale = 0.50) {
  list(
    nn_prior_zero_birth_fallback_weight = as.numeric(birth),
    nn_prior_zero_weight_cap_ratio = if (is.na(cap)) NA_real_ else as.numeric(cap),
    nn_prior_zero_exposure_quantile = as.numeric(quantile),
    nn_prior_zero_weight_scale = as.numeric(scale)
  )
}

expand_cfg_grid <- function(birth_levels, cap_levels, quantile_levels, scale_levels) {
  grid <- expand.grid(
    nn_prior_zero_birth_fallback_weight = birth_levels,
    nn_prior_zero_weight_cap_ratio = cap_levels,
    nn_prior_zero_exposure_quantile = quantile_levels,
    nn_prior_zero_weight_scale = scale_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  split(grid, seq_len(nrow(grid)))
}

cfg_from_row <- function(row) {
  make_cfg(
    birth = row$nn_prior_zero_birth_fallback_weight,
    cap = row$nn_prior_zero_weight_cap_ratio,
    quantile = row$nn_prior_zero_exposure_quantile,
    scale = row$nn_prior_zero_weight_scale
  )
}

same_level <- function(levels, value) {
  if (is.na(value)) {
    which(is.na(levels))
  } else {
    which(!is.na(levels) & abs(levels - value) < 1e-12)
  }
}

neighbor_levels <- function(levels, value) {
  idx <- same_level(levels, value)
  if (!length(idx)) {
    return(value)
  }
  keep <- unique(c(max(1L, idx - 1L), idx, min(length(levels), idx + 1L)))
  levels[keep]
}

dedupe_cfgs <- function(cfgs) {
  if (!length(cfgs)) {
    return(cfgs)
  }
  keys <- vapply(cfgs, combo_key, character(1))
  cfgs[!duplicated(keys)]
}

stage2_cfgs_from_top <- function(top_rows) {
  out <- list()
  idx <- 1L
  for (i in seq_len(nrow(top_rows))) {
    cfg <- cfg_from_row(top_rows[i, , drop = FALSE])
    out[[idx]] <- cfg
    idx <- idx + 1L

    for (birth in neighbor_levels(fine_birth_levels, cfg$nn_prior_zero_birth_fallback_weight)) {
      out[[idx]] <- make_cfg(
        birth = birth,
        cap = cfg$nn_prior_zero_weight_cap_ratio,
        quantile = cfg$nn_prior_zero_exposure_quantile,
        scale = cfg$nn_prior_zero_weight_scale
      )
      idx <- idx + 1L
    }

    for (cap in neighbor_levels(fine_cap_levels, cfg$nn_prior_zero_weight_cap_ratio)) {
      out[[idx]] <- make_cfg(
        birth = cfg$nn_prior_zero_birth_fallback_weight,
        cap = cap,
        quantile = cfg$nn_prior_zero_exposure_quantile,
        scale = cfg$nn_prior_zero_weight_scale
      )
      idx <- idx + 1L
    }

    for (quantile in neighbor_levels(fine_quantile_levels, cfg$nn_prior_zero_exposure_quantile)) {
      out[[idx]] <- make_cfg(
        birth = cfg$nn_prior_zero_birth_fallback_weight,
        cap = cfg$nn_prior_zero_weight_cap_ratio,
        quantile = quantile,
        scale = cfg$nn_prior_zero_weight_scale
      )
      idx <- idx + 1L
    }

    for (scale in fine_scale_levels) {
      out[[idx]] <- make_cfg(
        birth = cfg$nn_prior_zero_birth_fallback_weight,
        cap = cfg$nn_prior_zero_weight_cap_ratio,
        quantile = cfg$nn_prior_zero_exposure_quantile,
        scale = scale
      )
      idx <- idx + 1L
    }
  }
  dedupe_cfgs(out)
}

prepare_input <- function(input_rds) {
  yi <- readRDS(input_rds)
  yi$x <- as.data.frame(yi$x)
  if (diploid_state %in% rownames(yi$x)) {
    yi$x <- yi$x[rownames(yi$x) != diploid_state, , drop = FALSE]
  }
  yi
}

safe_mean <- function(x) {
  if (!length(x) || all(!is.finite(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (!length(x) || all(!is.finite(x))) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

summarise_diag <- function(diag_df) {
  if (is.null(diag_df) || !nrow(diag_df)) {
    return(list(
      n_replicates = 0L,
      mean_observed_children = NA_real_,
      mean_zero_children_retained = NA_real_,
      mean_birth_fallback_children = NA_real_,
      mean_fallback_share = NA_real_,
      mean_zero_weight_final = NA_real_,
      mean_zero_to_observed_weight_ratio = NA_real_,
      cap_applied_rate = NA_real_,
      mean_prior_sigma_hat = NA_real_,
      median_prior_sigma_hat = NA_real_,
      mean_lower_boundary_rate = NA_real_,
      mean_upper_boundary_rate = NA_real_,
      mean_no_prior_fallback_rate = NA_real_
    ))
  }

  zero_to_observed <- diag_df$sum_zero_weight_final / pmax(diag_df$sum_observed_weight, 1e-12)
  fallback_share <- diag_df$n_zero_children_with_birth_fallback / pmax(diag_df$n_zero_children_retained, 1)

  list(
    n_replicates = nrow(diag_df),
    mean_observed_children = safe_mean(diag_df$n_observed_children),
    mean_zero_children_retained = safe_mean(diag_df$n_zero_children_retained),
    mean_birth_fallback_children = safe_mean(diag_df$n_zero_children_with_birth_fallback),
    mean_fallback_share = safe_mean(fallback_share),
    mean_zero_weight_final = safe_mean(diag_df$sum_zero_weight_final),
    mean_zero_to_observed_weight_ratio = safe_mean(zero_to_observed),
    cap_applied_rate = safe_mean(as.numeric(diag_df$zero_weight_cap_applied)),
    mean_prior_sigma_hat = safe_mean(diag_df$prior_sigma_hat),
    median_prior_sigma_hat = safe_median(diag_df$prior_sigma_hat),
    mean_lower_boundary_rate = safe_mean(diag_df$map_delta_lower_boundary_rate),
    mean_upper_boundary_rate = safe_mean(diag_df$map_delta_upper_boundary_rate),
    mean_no_prior_fallback_rate = safe_mean(as.numeric(diag_df$used_no_prior_fallback_for_this_replicate))
  )
}

parallel_map <- function(tasks, fun) {
  if (!length(tasks)) {
    return(list())
  }
  if (effective_n_cores <= 1L || .Platform$OS.type == "windows") {
    return(lapply(tasks, fun))
  }
  parallel::mclapply(
    tasks,
    fun,
    mc.cores = effective_n_cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

run_one <- function(cfg, stage, seed, nboot, yi) {
  key <- combo_key(cfg)
  run_dir <- file.path(
    output_root,
    "runs",
    stage,
    sprintf("%s__seed_%d__nboot_%d", key, as.integer(seed), as.integer(nboot))
  )
  summary_rds <- file.path(run_dir, "run_summary.rds")
  warning_log <- file.path(run_dir, "warnings.log")

  if (file.exists(summary_rds)) {
    return(readRDS(summary_rds))
  }

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  started_at <- Sys.time()
  warning_messages <- character()
  status <- "ok"
  error_message <- NA_character_
  xval <- NA_real_
  diag_summary <- summarise_diag(NULL)

  tryCatch(
    withCallingHandlers({
      set.seed(seed)
      alfakR::alfak(
        yi = yi,
        outdir = run_dir,
        passage_times = NULL,
        minobs = fixed_cfg$minobs,
        nboot = nboot,
        n0 = fixed_cfg$n0,
        nb = fixed_cfg$nb,
        pm = fixed_cfg$pm,
        correct_efflux = fixed_cfg$correct_efflux,
        nn_prior = fixed_cfg$nn_prior,
        nn_prior_grid_n = fixed_cfg$nn_prior_grid_n,
        nn_prior_fit_subset = fixed_cfg$nn_prior_fit_subset,
        nn_prior_zero_birth_fallback_weight = cfg$nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_weight_cap_ratio = if (is.na(cfg$nn_prior_zero_weight_cap_ratio)) NULL else cfg$nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_exposure_quantile = cfg$nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = cfg$nn_prior_zero_weight_scale
      )
    }, warning = function(w) {
      warning_messages <<- unique(c(warning_messages, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  elapsed_sec <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

  xval_path <- file.path(run_dir, "xval.Rds")
  bootstrap_path <- file.path(run_dir, "bootstrap_res.Rds")
  if (status == "ok" && file.exists(xval_path)) {
    xval <- readRDS(xval_path)
  }
  if (status == "ok" && file.exists(bootstrap_path)) {
    boot_res <- readRDS(bootstrap_path)
    diag_summary <- summarise_diag(boot_res$nn_prior_diagnostics)
  }

  if (length(warning_messages)) {
    writeLines(warning_messages, warning_log)
  }

  row <- c(
    list(
      stage = stage,
      combo_key = key,
      seed = as.integer(seed),
      nboot = as.integer(nboot),
      status = status,
      error_message = error_message,
      elapsed_sec = elapsed_sec,
      warning_count = length(warning_messages),
      warning_messages = paste(warning_messages, collapse = " || "),
      outdir = run_dir,
      xval = xval
    ),
    cfg,
    diag_summary
  )

  saveRDS(row, summary_rds)
  row
}

run_stage <- function(cfgs, stage, seeds, nboot, yi) {
  tasks <- list()
  idx <- 1L
  for (cfg in cfgs) {
    for (seed in seeds) {
      tasks[[idx]] <- list(cfg = cfg, seed = as.integer(seed))
      idx <- idx + 1L
    }
  }

  message(
    sprintf(
      "Running %s for sample %s: %d task(s) on %d core(s).",
      stage,
      sample_id,
      length(tasks),
      effective_n_cores
    )
  )

  parallel_map(tasks, function(task) {
    run_one(
      cfg = task$cfg,
      stage = stage,
      seed = task$seed,
      nboot = nboot,
      yi = yi
    )
  })
}

rows_to_df <- function(rows) {
  if (!length(rows)) {
    return(data.frame())
  }
  as.data.frame(do.call(rbind, lapply(rows, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  })), stringsAsFactors = FALSE)
}

write_table <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
}

aggregate_stage <- function(df, group_cols) {
  if (!nrow(df)) {
    return(df)
  }

  split_keys <- interaction(df[, group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, split_keys)
  out <- lapply(pieces, function(piece) {
    ok_piece <- piece[piece$status == "ok" & is.finite(piece$xval), , drop = FALSE]
    data.frame(
      combo_key = piece$combo_key[[1]],
      nn_prior_zero_birth_fallback_weight = piece$nn_prior_zero_birth_fallback_weight[[1]],
      nn_prior_zero_weight_cap_ratio = piece$nn_prior_zero_weight_cap_ratio[[1]],
      nn_prior_zero_exposure_quantile = piece$nn_prior_zero_exposure_quantile[[1]],
      nn_prior_zero_weight_scale = piece$nn_prior_zero_weight_scale[[1]],
      runs = nrow(piece),
      ok_runs = nrow(ok_piece),
      mean_xval = safe_mean(ok_piece$xval),
      sd_xval = if (nrow(ok_piece) >= 2L) stats::sd(ok_piece$xval) else NA_real_,
      mean_elapsed_sec = safe_mean(piece$elapsed_sec),
      mean_warning_count = safe_mean(piece$warning_count),
      mean_fallback_share = safe_mean(piece$mean_fallback_share),
      mean_zero_to_observed_weight_ratio = safe_mean(piece$mean_zero_to_observed_weight_ratio),
      mean_prior_sigma_hat = safe_mean(piece$mean_prior_sigma_hat),
      mean_lower_boundary_rate = safe_mean(piece$mean_lower_boundary_rate),
      cap_applied_rate = safe_mean(piece$cap_applied_rate),
      stringsAsFactors = FALSE
    )
  })

  out_df <- do.call(rbind, out)
  out_df[order(-out_df$mean_xval, out_df$mean_elapsed_sec), , drop = FALSE]
}

write_best_summary <- function(best_row, stage3_summary, path) {
  lines <- c(
    sprintf("Sample: %s", sample_id),
    sprintf("Best combo key: %s", best_row$combo_key),
    sprintf("Mean xval across stage3 seeds: %.6f", best_row$mean_xval),
    sprintf("SD xval across stage3 seeds: %s", if (is.finite(best_row$sd_xval)) sprintf("%.6f", best_row$sd_xval) else "NA"),
    sprintf("birth_fallback_weight: %.6f", best_row$nn_prior_zero_birth_fallback_weight),
    sprintf("zero_weight_cap_ratio: %s", if (is.na(best_row$nn_prior_zero_weight_cap_ratio)) "adaptive" else sprintf("%.6f", best_row$nn_prior_zero_weight_cap_ratio)),
    sprintf("zero_exposure_quantile: %.6f", best_row$nn_prior_zero_exposure_quantile),
    sprintf("zero_weight_scale: %.6f", best_row$nn_prior_zero_weight_scale),
    sprintf("mean_fallback_share: %.6f", best_row$mean_fallback_share),
    sprintf("mean_zero_to_observed_weight_ratio: %.6f", best_row$mean_zero_to_observed_weight_ratio),
    sprintf("mean_prior_sigma_hat: %.6f", best_row$mean_prior_sigma_hat),
    sprintf("mean_lower_boundary_rate: %.6f", best_row$mean_lower_boundary_rate),
    "",
    "Top stage3 rows:",
    paste(capture.output(print(utils::head(stage3_summary, 10L), row.names = FALSE)), collapse = "\n")
  )
  writeLines(lines, path)
}

message("Input Rds: ", input_rds)
message("Sample: ", sample_id)
message("Output root: ", output_root)
message("Requested cores: ", requested_n_cores, "; effective cores: ", effective_n_cores)

loader_mode <- load_alfakr(repo_dir)
yi <- prepare_input(input_rds)

metadata <- list(
  loader_mode = loader_mode,
  input_rds = input_rds,
  sample_id = sample_id,
  sample_key = sample_key,
  output_root = output_root,
  calibration_root = calibration_root,
  repo_dir = repo_dir,
  requested_n_cores = requested_n_cores,
  effective_n_cores = effective_n_cores,
  fixed_cfg = fixed_cfg,
  search_space = list(
    coarse_birth_levels = coarse_birth_levels,
    coarse_cap_levels = coarse_cap_levels,
    coarse_quantile_levels = coarse_quantile_levels,
    coarse_scale_levels = coarse_scale_levels,
    fine_birth_levels = fine_birth_levels,
    fine_cap_levels = fine_cap_levels,
    fine_quantile_levels = fine_quantile_levels,
    fine_scale_levels = fine_scale_levels
  )
)
saveRDS(metadata, sample_output_path("search_metadata", "rds"))

stage1_cfgs <- lapply(
  expand_cfg_grid(
    birth_levels = coarse_birth_levels,
    cap_levels = coarse_cap_levels,
    quantile_levels = coarse_quantile_levels,
    scale_levels = coarse_scale_levels
  ),
  cfg_from_row
)

stage1_rows <- run_stage(
  cfgs = stage1_cfgs,
  stage = "stage1",
  seeds = fixed_cfg$benchmark_seed,
  nboot = fixed_cfg$nboot_stage1,
  yi = yi
)
stage1_df <- rows_to_df(stage1_rows)
stage1_path <- sample_output_path("stage1_runs", "tsv")
write_table(stage1_df, stage1_path)
stage1_summary <- aggregate_stage(stage1_df, group_cols = "combo_key")
stage1_summary_path <- sample_output_path("stage1_summary", "tsv")
write_table(stage1_summary, stage1_summary_path)

top_stage1 <- utils::head(stage1_summary[stage1_summary$ok_runs > 0L, , drop = FALSE], 5L)
if (!nrow(top_stage1)) {
  stop("Stage 1 produced no successful runs.")
}

stage2_cfgs <- stage2_cfgs_from_top(top_stage1)
stage2_rows <- run_stage(
  cfgs = stage2_cfgs,
  stage = "stage2",
  seeds = fixed_cfg$benchmark_seed,
  nboot = fixed_cfg$nboot_stage2,
  yi = yi
)
stage2_df <- rows_to_df(stage2_rows)
stage2_path <- sample_output_path("stage2_runs", "tsv")
write_table(stage2_df, stage2_path)
stage2_summary <- aggregate_stage(stage2_df, group_cols = "combo_key")
stage2_summary_path <- sample_output_path("stage2_summary", "tsv")
write_table(stage2_summary, stage2_summary_path)

top_stage2 <- utils::head(stage2_summary[stage2_summary$ok_runs > 0L, , drop = FALSE], 3L)
if (!nrow(top_stage2)) {
  stop("Stage 2 produced no successful runs.")
}

stage3_cfgs <- lapply(seq_len(nrow(top_stage2)), function(i) cfg_from_row(top_stage2[i, , drop = FALSE]))
stage3_rows <- run_stage(
  cfgs = stage3_cfgs,
  stage = "stage3",
  seeds = fixed_cfg$stage3_seeds,
  nboot = fixed_cfg$nboot_stage3,
  yi = yi
)
stage3_df <- rows_to_df(stage3_rows)
stage3_path <- sample_output_path("stage3_runs", "tsv")
write_table(stage3_df, stage3_path)
stage3_summary <- aggregate_stage(stage3_df, group_cols = "combo_key")
stage3_summary_path <- sample_output_path("stage3_summary", "tsv")
write_table(stage3_summary, stage3_summary_path)

if (!nrow(stage3_summary)) {
  stop("Stage 3 produced no successful runs.")
}

best_row <- stage3_summary[which.max(stage3_summary$mean_xval), , drop = FALSE]
write_table(best_row, sample_output_path("best_params", "tsv"))
saveRDS(best_row, sample_output_path("best_params", "rds"))
write_best_summary(best_row, stage3_summary, sample_output_path("best_params", "txt"))

message("Search complete.")
message("Best combo: ", best_row$combo_key[[1]])
message("Mean stage3 xval: ", sprintf("%.6f", best_row$mean_xval[[1]]))
