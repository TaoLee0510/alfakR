has_complete_alfak_outputs <- function(outdir) {
  required_paths <- file.path(
    outdir,
    c(
      "landscape.Rds",
      "bootstrap_res.Rds",
      "landscape_posterior_samples.Rds",
      "xval.Rds"
    )
  )
  all(file.exists(required_paths))
}

cohort_transition_output_setting <- function(outdir, field) {
  diag_path <- file.path(outdir, "nn_prior_diagnostics.Rds")
  diag_obj <- safe_read_rds(diag_path)
  if (!is.list(diag_obj) || is.null(diag_obj$replicate) || !is.data.frame(diag_obj$replicate)) {
    return(NA_character_)
  }

  value_vec <- diag_obj$replicate[[field]]
  if (is.null(value_vec)) {
    return(NA_character_)
  }
  value_vec <- unique(as.character(value_vec[!is.na(value_vec) & nzchar(value_vec)]))
  if (length(value_vec) != 1L) {
    return(NA_character_)
  }
  value_vec[[1L]]
}

cohort_transition_output_version <- function(outdir) {
  cohort_transition_output_setting(outdir, "cohort_transition_version")
}

cohort_transition_output_contextual_apply_to <- function(outdir) {
  cohort_transition_output_setting(outdir, "cohort_contextual_apply_to")
}

cohort_transition_output_keep_baseline_when_sparse <- function(outdir) {
  cohort_transition_output_setting(outdir, "cohort_context_keep_baseline_when_sparse")
}

cohort_transition_output_lambda_sparse_unknown <- function(outdir) {
  cohort_transition_output_setting(outdir, "cohort_context_lambda_sparse_unknown")
}

cohort_transition_output_matches <- function(outdir,
                                             nn_prior,
                                             cohort_transition_version = NA_character_,
                                             cohort_contextual_apply_to = NA_character_,
                                             cohort_context_keep_baseline_when_sparse = NA,
                                             cohort_context_lambda_sparse_unknown = NA_real_) {
  if (!identical(as.character(nn_prior), "cohort_transition")) {
    return(TRUE)
  }
  requested_version <- as.character(cohort_transition_version)
  if (length(requested_version) && !is.na(requested_version[[1L]]) && nzchar(requested_version[[1L]]) &&
      !identical(cohort_transition_output_version(outdir), requested_version[[1L]])) {
    return(FALSE)
  }

  requested_apply_to <- as.character(cohort_contextual_apply_to)
  if (length(requested_apply_to) && !is.na(requested_apply_to[[1L]]) && nzchar(requested_apply_to[[1L]]) &&
      !identical(cohort_transition_output_contextual_apply_to(outdir), requested_apply_to[[1L]])) {
    return(FALSE)
  }

  requested_keep_sparse_missing <- is.null(cohort_context_keep_baseline_when_sparse) ||
    !length(cohort_context_keep_baseline_when_sparse) ||
    (length(cohort_context_keep_baseline_when_sparse) == 1L && is.na(cohort_context_keep_baseline_when_sparse))
  if (!requested_keep_sparse_missing) {
    if (!same_optional_logical(
      cohort_transition_output_keep_baseline_when_sparse(outdir),
      cohort_context_keep_baseline_when_sparse
    )) {
      return(FALSE)
    }
  }

  requested_sparse_lambda_missing <- is.null(cohort_context_lambda_sparse_unknown) ||
    !length(cohort_context_lambda_sparse_unknown) ||
    (length(cohort_context_lambda_sparse_unknown) == 1L && is.na(cohort_context_lambda_sparse_unknown))
  if (!requested_sparse_lambda_missing) {
    if (!same_optional_numeric(
      cohort_transition_output_lambda_sparse_unknown(outdir),
      cohort_context_lambda_sparse_unknown
    )) {
      return(FALSE)
    }
  }

  TRUE
}

same_optional_logical <- function(lhs, rhs) {
  lhs_missing <- is.null(lhs) || !length(lhs) || (length(lhs) == 1L && is.na(lhs))
  rhs_missing <- is.null(rhs) || !length(rhs) || (length(rhs) == 1L && is.na(rhs))
  if (lhs_missing && rhs_missing) {
    return(TRUE)
  }
  if (lhs_missing || rhs_missing) {
    return(FALSE)
  }
  lhs_log <- suppressWarnings(as.logical(lhs[[1]]))
  rhs_log <- suppressWarnings(as.logical(rhs[[1]]))
  is.logical(lhs_log) && is.logical(rhs_log) && !is.na(lhs_log) && !is.na(rhs_log) && identical(lhs_log, rhs_log)
}

same_optional_numeric <- function(lhs, rhs, tol = 1e-12) {
  lhs_missing <- is.null(lhs) || !length(lhs) || (length(lhs) == 1L && is.na(lhs))
  rhs_missing <- is.null(rhs) || !length(rhs) || (length(rhs) == 1L && is.na(rhs))
  if (lhs_missing && rhs_missing) {
    return(TRUE)
  }
  if (lhs_missing || rhs_missing) {
    return(FALSE)
  }
  lhs_num <- suppressWarnings(as.numeric(lhs[[1]]))
  rhs_num <- suppressWarnings(as.numeric(rhs[[1]]))
  is.finite(lhs_num) && is.finite(rhs_num) && abs(lhs_num - rhs_num) <= tol
}

weighted_prior_cache_matches <- function(cached,
                                         nn_prior,
                                         nn_prior_grid_n,
                                         nn_prior_fit_subset,
                                         nn_prior_zero_exposure_quantile,
                                         nn_prior_zero_weight_scale,
                                         nn_prior_zero_weight_cap_ratio,
                                         nn_prior_zero_birth_fallback_weight,
                                         nn_prior_zero_birth_child_floor,
                                         nn_prior_zero_birth_child_shape,
                                         nn_prior_zero_birth_replicate_floor,
                                         nn_prior_zero_birth_replicate_shape,
                                         nn_prior_two_step_support,
                                         nn_prior_two_step_support_min,
                                         nn_prior_two_step_cap_floor) {
  if (identical(nn_prior, "empirical_censored") &&
      !same_optional_numeric(cached$nn_prior_grid_n, nn_prior_grid_n)) {
    return(FALSE)
  }

  if (!nn_prior %in% c("empirical_censored_weighted", "empirical_two_shell")) {
    return(TRUE)
  }

  same_optional_numeric(cached$nn_prior_grid_n, nn_prior_grid_n) &&
    identical(as.character(cached$nn_prior_fit_subset), nn_prior_fit_subset) &&
    same_optional_numeric(cached$nn_prior_zero_exposure_quantile, nn_prior_zero_exposure_quantile) &&
    same_optional_numeric(cached$nn_prior_zero_weight_scale, nn_prior_zero_weight_scale) &&
    same_optional_numeric(cached$nn_prior_zero_weight_cap_ratio, nn_prior_zero_weight_cap_ratio) &&
    same_optional_numeric(cached$nn_prior_zero_birth_fallback_weight, nn_prior_zero_birth_fallback_weight) &&
    same_optional_numeric(cached$nn_prior_zero_birth_child_floor, nn_prior_zero_birth_child_floor) &&
    same_optional_numeric(cached$nn_prior_zero_birth_child_shape, nn_prior_zero_birth_child_shape) &&
    same_optional_numeric(cached$nn_prior_zero_birth_replicate_floor, nn_prior_zero_birth_replicate_floor) &&
    same_optional_numeric(cached$nn_prior_zero_birth_replicate_shape, nn_prior_zero_birth_replicate_shape) &&
    identical(as.character(cached$nn_prior_two_step_support), nn_prior_two_step_support) &&
    same_optional_numeric(cached$nn_prior_two_step_support_min, nn_prior_two_step_support_min) &&
    same_optional_numeric(cached$nn_prior_two_step_cap_floor, nn_prior_two_step_cap_floor)
}

cohort_transition_cache_matches <- function(cached,
                                            outdir,
                                            nn_prior,
                                            cohort_transition_version = NA_character_,
                                            cohort_contextual_apply_to = NA_character_,
                                            cohort_context_keep_baseline_when_sparse = NA,
                                            cohort_context_lambda_sparse_unknown = NA_real_) {
  if (!identical(as.character(nn_prior), "cohort_transition")) {
    return(TRUE)
  }
  requested_version <- as.character(cohort_transition_version)
  if (length(requested_version) && !is.na(requested_version[[1L]]) && nzchar(requested_version[[1L]])) {
    cached_version <- as.character(cached$cohort_transition_version)
    if (length(cached_version) && !is.na(cached_version[[1L]]) && nzchar(cached_version[[1L]])) {
      if (!identical(cached_version[[1L]], requested_version[[1L]])) {
        return(FALSE)
      }
    } else if (!identical(cohort_transition_output_version(outdir), requested_version[[1L]])) {
      return(FALSE)
    }
  }

  requested_apply_to <- as.character(cohort_contextual_apply_to)
  if (length(requested_apply_to) && !is.na(requested_apply_to[[1L]]) && nzchar(requested_apply_to[[1L]])) {
    cached_apply_to <- as.character(cached$cohort_contextual_apply_to)
    if (length(cached_apply_to) && !is.na(cached_apply_to[[1L]]) && nzchar(cached_apply_to[[1L]])) {
      if (!identical(cached_apply_to[[1L]], requested_apply_to[[1L]])) {
        return(FALSE)
      }
    } else if (!identical(cohort_transition_output_contextual_apply_to(outdir), requested_apply_to[[1L]])) {
      return(FALSE)
    }
  }

  requested_keep_sparse_missing <- is.null(cohort_context_keep_baseline_when_sparse) ||
    !length(cohort_context_keep_baseline_when_sparse) ||
    (length(cohort_context_keep_baseline_when_sparse) == 1L && is.na(cohort_context_keep_baseline_when_sparse))
  if (!requested_keep_sparse_missing) {
    cached_keep_sparse <- if ("cohort_context_keep_baseline_when_sparse" %in% names(cached)) {
      cached$cohort_context_keep_baseline_when_sparse
    } else {
      cohort_transition_output_keep_baseline_when_sparse(outdir)
    }
    if (!same_optional_logical(cached_keep_sparse, cohort_context_keep_baseline_when_sparse)) {
      return(FALSE)
    }
  }

  requested_sparse_lambda_missing <- is.null(cohort_context_lambda_sparse_unknown) ||
    !length(cohort_context_lambda_sparse_unknown) ||
    (length(cohort_context_lambda_sparse_unknown) == 1L && is.na(cohort_context_lambda_sparse_unknown))
  if (!requested_sparse_lambda_missing) {
    cached_sparse_lambda <- if ("cohort_context_lambda_sparse_unknown" %in% names(cached)) {
      cached$cohort_context_lambda_sparse_unknown
    } else {
      cohort_transition_output_lambda_sparse_unknown(outdir)
    }
    if (!same_optional_numeric(cached_sparse_lambda, cohort_context_lambda_sparse_unknown)) {
      return(FALSE)
    }
  }

  TRUE
}

build_existing_fit_row <- function(patient_id,
                                   outdir,
                                   pm,
                                   minobs,
                                   benchmark_seed,
                                   parameter_label,
                                   nn_prior,
                                   nn_prior_grid_n,
                                   nn_prior_fit_subset,
                                   nn_prior_zero_exposure_quantile,
                                   nn_prior_zero_weight_scale,
                                   nn_prior_zero_weight_cap_ratio,
                                   nn_prior_zero_birth_fallback_weight,
                                   nn_prior_zero_birth_child_floor,
                                   nn_prior_zero_birth_child_shape,
                                   nn_prior_zero_birth_replicate_floor,
                                   nn_prior_zero_birth_replicate_shape,
                                   nn_prior_two_step_support,
                                   nn_prior_two_step_support_min,
                                   nn_prior_two_step_cap_floor,
                                   cohort_transition_version = NA_character_,
                                   cohort_contextual_apply_to = NA_character_,
                                   cohort_context_keep_baseline_when_sparse = NA,
                                   cohort_context_lambda_sparse_unknown = NA_real_,
                                   warning_log_path,
                                   landscape_path,
                                   bootstrap_path,
                                   posterior_path,
                                   xval_path) {
  warning_lines <- if (file.exists(warning_log_path)) readLines(warning_log_path, warn = FALSE) else character()
  xval_obj <- safe_read_rds(xval_path)

  c(
    list(
      patient_id = patient_id,
      outdir = outdir,
      pm = pm,
      pm_label = pm_to_label(pm),
      minobs = minobs,
      benchmark_seed = benchmark_seed,
      parameter_label = parameter_label,
      nn_prior = nn_prior,
      nn_prior_grid_n = nn_prior_grid_n,
      nn_prior_fit_subset = nn_prior_fit_subset,
      nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
      nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
      nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
      nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
      nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
      nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
      nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
      nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
      nn_prior_two_step_support = nn_prior_two_step_support,
      nn_prior_two_step_support_min = nn_prior_two_step_support_min,
      nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
      cohort_transition_version = if (identical(nn_prior, "cohort_transition")) as.character(cohort_transition_version) else NA_character_,
      cohort_contextual_apply_to = if (identical(nn_prior, "cohort_transition")) as.character(cohort_contextual_apply_to) else NA_character_,
      cohort_context_keep_baseline_when_sparse = if (identical(nn_prior, "cohort_transition")) as.logical(cohort_context_keep_baseline_when_sparse) else NA,
      cohort_context_lambda_sparse_unknown = if (identical(nn_prior, "cohort_transition")) as.numeric(cohort_context_lambda_sparse_unknown) else NA_real_,
      status = "ok",
      cached = TRUE,
      error_message = NA_character_,
      elapsed_sec = NA_real_,
      warning_count = length(warning_lines),
      lambda_endpoint_warning_count = sum(grepl("^Grid searches over lambda", warning_lines)),
      warning_messages = if (length(warning_lines)) paste(warning_lines, collapse = " || ") else NA_character_,
      landscape_path = if (file.exists(landscape_path)) landscape_path else NA_character_,
      bootstrap_path = if (file.exists(bootstrap_path)) bootstrap_path else NA_character_,
      posterior_path = if (file.exists(posterior_path)) posterior_path else NA_character_,
      xval_path = if (file.exists(xval_path)) xval_path else NA_character_
    ),
    extract_xval_metrics(xval_obj)
  )
}

extract_xval_metrics <- function(xv) {
  metrics <- list(
    xval_r2 = NA_real_,
    xval_cor = NA_real_,
    xval_rmse = NA_real_,
    xval_mae = NA_real_,
    n_xval = NA_integer_
  )

  if (is.null(xv)) {
    return(metrics)
  }

  if (is.atomic(xv) && is.numeric(xv) && is.null(dim(xv)) && length(xv) >= 1L) {
    metrics$xval_r2 <- suppressWarnings(as.numeric(xv[[1]]))
    return(metrics)
  }

  if (!is.list(xv)) {
    return(metrics)
  }

  if (!is.null(xv$R2R)) {
    metrics$xval_r2 <- suppressWarnings(as.numeric(xv$R2R))
  }

  if (is.null(xv$tmp)) {
    return(metrics)
  }

  tmp_df <- as.data.frame(xv$tmp)
  if (ncol(tmp_df) < 2) {
    return(metrics)
  }

  colnames(tmp_df)[1:2] <- c("f_est", "f_xv")
  tmp_df <- tmp_df[, c("f_est", "f_xv"), drop = FALSE]
  tmp_df$f_est <- suppressWarnings(as.numeric(tmp_df$f_est))
  tmp_df$f_xv <- suppressWarnings(as.numeric(tmp_df$f_xv))
  tmp_df <- tmp_df[is.finite(tmp_df$f_est) & is.finite(tmp_df$f_xv), , drop = FALSE]

  metrics$n_xval <- nrow(tmp_df)
  if (!nrow(tmp_df)) {
    return(metrics)
  }

  diff_vec <- tmp_df$f_est - tmp_df$f_xv
  metrics$xval_rmse <- sqrt(mean(diff_vec^2))
  metrics$xval_mae <- mean(abs(diff_vec))
  if (nrow(tmp_df) >= 2) {
    metrics$xval_cor <- suppressWarnings(stats::cor(tmp_df$f_est, tmp_df$f_xv, use = "complete.obs"))
  }
  metrics
}

refresh_cached_fit_row <- function(cached,
                                   outdir,
                                   patient_id,
                                   pm,
                                   minobs,
                                   benchmark_seed,
                                   parameter_label,
                                   nn_prior,
                                   nn_prior_grid_n,
                                   nn_prior_fit_subset,
                                   nn_prior_zero_exposure_quantile,
                                   nn_prior_zero_weight_scale,
                                   nn_prior_zero_weight_cap_ratio,
                                   nn_prior_zero_birth_fallback_weight,
                                   nn_prior_zero_birth_child_floor,
                                   nn_prior_zero_birth_child_shape,
                                   nn_prior_zero_birth_replicate_floor,
                                   nn_prior_zero_birth_replicate_shape,
                                   nn_prior_two_step_support,
                                   nn_prior_two_step_support_min,
                                   nn_prior_two_step_cap_floor,
                                   cohort_transition_version = NA_character_,
                                   cohort_contextual_apply_to = NA_character_,
                                   cohort_context_keep_baseline_when_sparse = NA,
                                   cohort_context_lambda_sparse_unknown = NA_real_,
                                   warning_log_path,
                                   landscape_path,
                                   bootstrap_path,
                                   posterior_path,
                                   xval_path) {
  if (!is.list(cached)) {
    cached <- list()
  }

  warning_lines <- if (file.exists(warning_log_path)) readLines(warning_log_path, warn = FALSE) else character()
  cached$patient_id <- patient_id
  cached$outdir <- outdir
  cached$pm <- pm
  cached$pm_label <- pm_to_label(pm)
  cached$minobs <- minobs
  cached$benchmark_seed <- benchmark_seed
  cached$parameter_label <- parameter_label
  cached$nn_prior <- nn_prior
  cached$nn_prior_grid_n <- nn_prior_grid_n
  cached$nn_prior_fit_subset <- nn_prior_fit_subset
  cached$nn_prior_zero_exposure_quantile <- nn_prior_zero_exposure_quantile
  cached$nn_prior_zero_weight_scale <- nn_prior_zero_weight_scale
  cached$nn_prior_zero_weight_cap_ratio <- nn_prior_zero_weight_cap_ratio
  cached$nn_prior_zero_birth_fallback_weight <- nn_prior_zero_birth_fallback_weight
  cached$nn_prior_zero_birth_child_floor <- nn_prior_zero_birth_child_floor
  cached$nn_prior_zero_birth_child_shape <- nn_prior_zero_birth_child_shape
  cached$nn_prior_zero_birth_replicate_floor <- nn_prior_zero_birth_replicate_floor
  cached$nn_prior_zero_birth_replicate_shape <- nn_prior_zero_birth_replicate_shape
  cached$nn_prior_two_step_support <- nn_prior_two_step_support
  cached$nn_prior_two_step_support_min <- nn_prior_two_step_support_min
  cached$nn_prior_two_step_cap_floor <- nn_prior_two_step_cap_floor
  cached$cohort_transition_version <- if (identical(nn_prior, "cohort_transition")) as.character(cohort_transition_version) else NA_character_
  cached$cohort_contextual_apply_to <- if (identical(nn_prior, "cohort_transition")) as.character(cohort_contextual_apply_to) else NA_character_
  cached$cohort_context_keep_baseline_when_sparse <- if (identical(nn_prior, "cohort_transition")) as.logical(cohort_context_keep_baseline_when_sparse) else NA
  cached$cohort_context_lambda_sparse_unknown <- if (identical(nn_prior, "cohort_transition")) as.numeric(cohort_context_lambda_sparse_unknown) else NA_real_
  cached$cached <- TRUE
  cached$warning_count <- length(warning_lines)
  cached$lambda_endpoint_warning_count <- sum(grepl("^Grid searches over lambda", warning_lines))
  cached$warning_messages <- if (length(warning_lines)) paste(warning_lines, collapse = " || ") else NA_character_
  cached$landscape_path <- if (file.exists(landscape_path)) landscape_path else NA_character_
  cached$bootstrap_path <- if (file.exists(bootstrap_path)) bootstrap_path else NA_character_
  cached$posterior_path <- if (file.exists(posterior_path)) posterior_path else NA_character_
  cached$xval_path <- if (file.exists(xval_path)) xval_path else NA_character_
  cached
}

fit_path_component_label <- function(value, key) {
  if (is.null(value) || (length(value) == 1L && is.na(value))) {
    return(NA_character_)
  }

  switch(
    key,
    parameter_label = as.character(value),
    pm = paste0("pm_", pm_to_label(as.numeric(value))),
    minobs = paste0("MINOBS_", as.integer(value)),
    patient_id = as.character(value),
    stop("Unsupported fit path component key: ", key)
  )
}

build_fit_outdir <- function(root_dir, patient_id, path_components) {
  component_labels <- purrr::imap_chr(path_components, fit_path_component_label)
  component_labels <- component_labels[is.finite(nchar(component_labels))]
  do.call(file.path, c(list(root_dir), as.list(component_labels), list(patient_id)))
}

infer_parameter_label <- function(nn_prior) {
  paste0("nn_prior_", nn_prior)
}

build_parameter_spec_tbl <- function(parameter_labels) {
  tibble::tibble(
    parameter_label = parameter_labels,
    nn_prior = sub("^nn_prior_", "", parameter_labels)
  ) %>%
    dplyr::mutate(parameter_label = factor(parameter_label, levels = parameter_labels)) %>%
    dplyr::arrange(parameter_label) %>%
    dplyr::mutate(parameter_label = as.character(parameter_label))
}

reconcile_fit_results_tbl <- function(fit_results_tbl, task_tbl = NULL) {
  if (is.null(fit_results_tbl) || !nrow(fit_results_tbl)) {
    return(tibble::tibble())
  }

  fit_results_tbl <- tibble::as_tibble(fit_results_tbl)
  if (!"parameter_label" %in% names(fit_results_tbl)) {
    if ("nn_prior" %in% names(fit_results_tbl)) {
      fit_results_tbl$parameter_label <- infer_parameter_label(fit_results_tbl$nn_prior)
    } else {
      fit_results_tbl$parameter_label <- NA_character_
    }
  }

  if (!is.null(task_tbl) && nrow(task_tbl)) {
    task_idx_tbl <- task_tbl %>%
      dplyr::transmute(
        patient_id,
        minobs,
        pm,
        parameter_label,
        nn_prior,
        expected_outdir = outdir
      )

    fit_results_tbl <- fit_results_tbl %>%
      dplyr::rename(recorded_outdir = outdir) %>%
      dplyr::left_join(task_idx_tbl, by = c("patient_id", "minobs", "pm", "parameter_label", "nn_prior")) %>%
      dplyr::mutate(outdir = dplyr::coalesce(expected_outdir, recorded_outdir)) %>%
      dplyr::select(-recorded_outdir, -expected_outdir)
  }

  fit_results_tbl %>%
    dplyr::mutate(
      landscape_path = ifelse(!is.na(outdir) & nzchar(outdir), file.path(outdir, "landscape.Rds"), NA_character_),
      bootstrap_path = ifelse(!is.na(outdir) & nzchar(outdir), file.path(outdir, "bootstrap_res.Rds"), NA_character_),
      posterior_path = ifelse(!is.na(outdir) & nzchar(outdir), file.path(outdir, "landscape_posterior_samples.Rds"), NA_character_),
      xval_path = ifelse(!is.na(outdir) & nzchar(outdir), file.path(outdir, "xval.Rds"), NA_character_)
    ) %>%
    dplyr::select(-dplyr::any_of("fit_mode")) %>%
    dplyr::relocate(parameter_label, .after = benchmark_seed)
}

run_alfak_fit <- function(patient_id,
                          input_rds,
                          outdir,
                          minobs,
                          pm,
                          nboot,
                          n0,
                          nb,
                          benchmark_seed,
                          parameter_label,
                          diploid_state,
                          correct_efflux = TRUE,
                          nn_prior = "none",
                          nn_prior_grid_n = 81L,
                          nn_prior_fit_subset = "hybrid",
                          nn_prior_zero_exposure_quantile = 0.10,
                          nn_prior_zero_weight_scale = 0.50,
                          nn_prior_zero_weight_cap_ratio = NA_real_,
                          nn_prior_zero_birth_fallback_weight = NA_real_,
                          nn_prior_zero_birth_child_floor = 0.25,
                          nn_prior_zero_birth_child_shape = 1,
                          nn_prior_zero_birth_replicate_floor = 0.50,
                          nn_prior_zero_birth_replicate_shape = 1,
                          nn_prior_two_step_support = "none",
                          nn_prior_two_step_support_min = 0.15,
                          nn_prior_two_step_cap_floor = 0.30,
                          cohort_transition_version = NA_character_,
                          cohort_contextual_apply_to = NA_character_,
                          cohort_context_keep_baseline_when_sparse = NA,
                          cohort_context_lambda_sparse_unknown = NA_real_,
                          force_refit = FALSE) {
  landscape_path <- file.path(outdir, "landscape.Rds")
  bootstrap_path <- file.path(outdir, "bootstrap_res.Rds")
  posterior_path <- file.path(outdir, "landscape_posterior_samples.Rds")
  xval_path <- file.path(outdir, "xval.Rds")
  summary_path <- file.path(outdir, "fit_task_result.rds")
  warning_log_path <- file.path(outdir, "fit_warnings.log")
  pm_label <- pm_to_label(pm)
  task_tag <- paste0(
    patient_id,
    " | minobs=", minobs,
    " | pm=", pm_label,
    " | parameter_label=", parameter_label,
    " | nn_prior=", nn_prior,
    " | grid=", nn_prior_grid_n,
    if (nn_prior %in% c("empirical_censored_weighted", "empirical_two_shell")) paste0(
      " | fit_subset=", nn_prior_fit_subset,
      " | zero_q=", signif(nn_prior_zero_exposure_quantile, 4),
      " | zero_scale=", signif(nn_prior_zero_weight_scale, 4),
      " | zero_cap=", if (is.na(nn_prior_zero_weight_cap_ratio)) "adaptive" else signif(nn_prior_zero_weight_cap_ratio, 4),
      if (!is.na(nn_prior_zero_birth_fallback_weight)) paste0(
        " | zero_birth_fallback_alias=", signif(nn_prior_zero_birth_fallback_weight, 4)
      ) else "",
      " | zero_birth_child_floor=", signif(nn_prior_zero_birth_child_floor, 4),
      " | zero_birth_child_shape=", signif(nn_prior_zero_birth_child_shape, 4),
      " | zero_birth_replicate_floor=", signif(nn_prior_zero_birth_replicate_floor, 4),
      " | zero_birth_replicate_shape=", signif(nn_prior_zero_birth_replicate_shape, 4),
      " | two_step_support=", nn_prior_two_step_support,
      " | two_step_support_min=", signif(nn_prior_two_step_support_min, 4),
      " | two_step_cap_floor=", signif(nn_prior_two_step_cap_floor, 4)
    ) else ""
  )

  if (!force_refit && file.exists(summary_path)) {
    cached <- tryCatch(readRDS(summary_path), error = function(e) NULL)
    if (is.list(cached) &&
        identical(cached$status, "ok") &&
        has_complete_alfak_outputs(outdir) &&
        weighted_prior_cache_matches(
          cached = cached,
          nn_prior = nn_prior,
          nn_prior_grid_n = nn_prior_grid_n,
          nn_prior_fit_subset = nn_prior_fit_subset,
          nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
          nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
          nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
          nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
          nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
          nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
          nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
          nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
          nn_prior_two_step_support = nn_prior_two_step_support,
          nn_prior_two_step_support_min = nn_prior_two_step_support_min,
          nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor
        ) &&
          cohort_transition_cache_matches(
            cached = cached,
            outdir = outdir,
            nn_prior = nn_prior,
            cohort_transition_version = cohort_transition_version,
            cohort_contextual_apply_to = cohort_contextual_apply_to,
            cohort_context_keep_baseline_when_sparse = cohort_context_keep_baseline_when_sparse,
            cohort_context_lambda_sparse_unknown = cohort_context_lambda_sparse_unknown
          )) {
      cached <- refresh_cached_fit_row(
        cached = cached,
        outdir = outdir,
        patient_id = patient_id,
        pm = pm,
        minobs = minobs,
        benchmark_seed = benchmark_seed,
        parameter_label = parameter_label,
        nn_prior = nn_prior,
        nn_prior_grid_n = nn_prior_grid_n,
        nn_prior_fit_subset = nn_prior_fit_subset,
        nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
        nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
        nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
        nn_prior_two_step_support = nn_prior_two_step_support,
        nn_prior_two_step_support_min = nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
        cohort_transition_version = cohort_transition_version,
        cohort_contextual_apply_to = cohort_contextual_apply_to,
        cohort_context_keep_baseline_when_sparse = cohort_context_keep_baseline_when_sparse,
        cohort_context_lambda_sparse_unknown = cohort_context_lambda_sparse_unknown,
        warning_log_path = warning_log_path,
        landscape_path = landscape_path,
        bootstrap_path = bootstrap_path,
        posterior_path = posterior_path,
        xval_path = xval_path
      )
      saveRDS(cached, summary_path)
      alfak_log("ALFA-K cached: ", task_tag)
      return(tibble::as_tibble(cached))
    }
  }

  if (!force_refit &&
      has_complete_alfak_outputs(outdir) &&
      cohort_transition_output_matches(
        outdir,
        nn_prior,
        cohort_transition_version,
        cohort_contextual_apply_to,
        cohort_context_keep_baseline_when_sparse,
        cohort_context_lambda_sparse_unknown
      )) {
    cached <- build_existing_fit_row(
      patient_id = patient_id,
      outdir = outdir,
      pm = pm,
      minobs = minobs,
      benchmark_seed = benchmark_seed,
      parameter_label = parameter_label,
      nn_prior = nn_prior,
      nn_prior_grid_n = nn_prior_grid_n,
      nn_prior_fit_subset = nn_prior_fit_subset,
      nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
      nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
      nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
      nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
      nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
      nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
      nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
      nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
      nn_prior_two_step_support = nn_prior_two_step_support,
      nn_prior_two_step_support_min = nn_prior_two_step_support_min,
      nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
      cohort_transition_version = cohort_transition_version,
      cohort_contextual_apply_to = cohort_contextual_apply_to,
      cohort_context_keep_baseline_when_sparse = cohort_context_keep_baseline_when_sparse,
      cohort_context_lambda_sparse_unknown = cohort_context_lambda_sparse_unknown,
      warning_log_path = warning_log_path,
      landscape_path = landscape_path,
      bootstrap_path = bootstrap_path,
      posterior_path = posterior_path,
      xval_path = xval_path
    )
    saveRDS(cached, summary_path)
    alfak_log("ALFA-K adopted existing outputs: ", task_tag)
    return(tibble::as_tibble(cached))
  }

  alfak_log("ALFA-K start: ", task_tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (identical(nn_prior, "cohort_transition")) {
    stop(
      "Benchmark does not generate `nn_prior = \"cohort_transition\"` fits directly. ",
      "Place precomputed cohort-transition outputs in `", outdir, "` and rerun with `force_refit = FALSE`."
    )
  }

  yi <- readRDS(input_rds)
  yi$x <- as.data.frame(yi$x)
  if (diploid_state %in% rownames(yi$x)) {
    yi$x <- yi$x[rownames(yi$x) != diploid_state, , drop = FALSE]
  }
  if (!nrow(yi$x)) {
    stop("All rows were filtered out for input: ", input_rds)
  }
  if (max(rowSums(yi$x), na.rm = TRUE) < minobs) {
    stop("No frequent karyotypes reach minobs=", minobs, " for input: ", input_rds)
  }

  warning_messages <- character()
  started_at <- Sys.time()

  fit_row <- tryCatch({
    set.seed(benchmark_seed)
    withCallingHandlers({
      alfakR::alfak(
        yi = yi,
        outdir = outdir,
        passage_times = NULL,
        minobs = minobs,
        nboot = nboot,
        n0 = n0,
        nb = nb,
        pm = pm,
        correct_efflux = correct_efflux,
        nn_prior = nn_prior,
        nn_prior_grid_n = nn_prior_grid_n,
        nn_prior_fit_subset = nn_prior_fit_subset,
        nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = if (is.na(nn_prior_zero_weight_cap_ratio)) NULL else nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_birth_fallback_weight = if (is.na(nn_prior_zero_birth_fallback_weight)) NULL else nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
        nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
        nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
        nn_prior_two_step_support = nn_prior_two_step_support,
        nn_prior_two_step_support_min = nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor
      )
    }, warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

    xv <- readRDS(xval_path)
    elapsed_sec <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    warning_messages <- unique(warning_messages)
    if (length(warning_messages)) {
      writeLines(warning_messages, warning_log_path)
    }

    c(
      list(
        patient_id = patient_id,
        outdir = outdir,
        pm = pm,
        pm_label = pm_label,
        minobs = minobs,
        benchmark_seed = benchmark_seed,
        parameter_label = parameter_label,
        nn_prior = nn_prior,
        nn_prior_grid_n = nn_prior_grid_n,
        nn_prior_fit_subset = nn_prior_fit_subset,
        nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
        nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
        nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
        nn_prior_two_step_support = nn_prior_two_step_support,
        nn_prior_two_step_support_min = nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
        cohort_transition_version = if (identical(nn_prior, "cohort_transition")) as.character(cohort_transition_version) else NA_character_,
        cohort_contextual_apply_to = if (identical(nn_prior, "cohort_transition")) as.character(cohort_contextual_apply_to) else NA_character_,
        cohort_context_keep_baseline_when_sparse = if (identical(nn_prior, "cohort_transition")) as.logical(cohort_context_keep_baseline_when_sparse) else NA,
        cohort_context_lambda_sparse_unknown = if (identical(nn_prior, "cohort_transition")) as.numeric(cohort_context_lambda_sparse_unknown) else NA_real_,
        status = "ok",
        cached = FALSE,
        error_message = NA_character_,
        elapsed_sec = elapsed_sec,
        warning_count = length(warning_messages),
        lambda_endpoint_warning_count = sum(grepl("^Grid searches over lambda", warning_messages)),
        warning_messages = if (length(warning_messages)) paste(warning_messages, collapse = " || ") else NA_character_,
        landscape_path = landscape_path,
        bootstrap_path = if (file.exists(bootstrap_path)) bootstrap_path else NA_character_,
        posterior_path = if (file.exists(posterior_path)) posterior_path else NA_character_,
        xval_path = xval_path
      ),
      extract_xval_metrics(xv)
    )
  }, error = function(e) {
    elapsed_sec <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    warning_messages <- unique(warning_messages)
    if (length(warning_messages)) {
      writeLines(warning_messages, warning_log_path)
    }

    alfak_log("ALFA-K error: ", task_tag, " | ", conditionMessage(e))
    list(
      patient_id = patient_id,
      outdir = outdir,
      pm = pm,
      pm_label = pm_label,
      minobs = minobs,
      benchmark_seed = benchmark_seed,
      parameter_label = parameter_label,
      nn_prior = nn_prior,
      nn_prior_grid_n = nn_prior_grid_n,
      nn_prior_fit_subset = nn_prior_fit_subset,
      nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
      nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
      nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
      nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
      nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
      nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
      nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
      nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
      nn_prior_two_step_support = nn_prior_two_step_support,
      nn_prior_two_step_support_min = nn_prior_two_step_support_min,
      nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
      cohort_transition_version = if (identical(nn_prior, "cohort_transition")) as.character(cohort_transition_version) else NA_character_,
      cohort_contextual_apply_to = if (identical(nn_prior, "cohort_transition")) as.character(cohort_contextual_apply_to) else NA_character_,
      cohort_context_keep_baseline_when_sparse = if (identical(nn_prior, "cohort_transition")) as.logical(cohort_context_keep_baseline_when_sparse) else NA,
      cohort_context_lambda_sparse_unknown = if (identical(nn_prior, "cohort_transition")) as.numeric(cohort_context_lambda_sparse_unknown) else NA_real_,
      status = "error",
      cached = FALSE,
      error_message = conditionMessage(e),
      elapsed_sec = elapsed_sec,
      warning_count = length(warning_messages),
      lambda_endpoint_warning_count = sum(grepl("^Grid searches over lambda", warning_messages)),
      warning_messages = if (length(warning_messages)) paste(warning_messages, collapse = " || ") else NA_character_,
      xval_r2 = NA_real_,
      xval_cor = NA_real_,
      xval_rmse = NA_real_,
      xval_mae = NA_real_,
      n_xval = NA_integer_,
      landscape_path = if (file.exists(landscape_path)) landscape_path else NA_character_,
      bootstrap_path = if (file.exists(bootstrap_path)) bootstrap_path else NA_character_,
      posterior_path = if (file.exists(posterior_path)) posterior_path else NA_character_,
      xval_path = if (file.exists(xval_path)) xval_path else NA_character_
    )
  })

  saveRDS(fit_row, summary_path)
  if (identical(fit_row$status, "ok")) {
    alfak_log("ALFA-K done: ", task_tag, " | xval_r2=", signif(fit_row$xval_r2, 4))
  }
  tibble::as_tibble(fit_row)
}

task_outdir_parameter <- function(root_dir, patient_id, minobs, pm, parameter_label) {
  build_fit_outdir(
    root_dir = root_dir,
    patient_id = patient_id,
    path_components = list(
      parameter_label = parameter_label,
      pm = pm,
      minobs = minobs
    )
  )
}

build_parameter_tasks <- function(input_index_tbl,
                                  fit_root,
                                  minobs_values,
                                  pm_values,
                                  parameter_spec_tbl,
                                  nn_prior_grid_n,
                                  nn_prior_fit_subset,
                                  nn_prior_zero_exposure_quantile,
                                  nn_prior_zero_weight_scale,
                                  nn_prior_zero_weight_cap_ratio,
                                  nn_prior_zero_birth_fallback_weight,
                                  nn_prior_zero_birth_child_floor,
                                  nn_prior_zero_birth_child_shape,
                                  nn_prior_zero_birth_replicate_floor,
                                  nn_prior_zero_birth_replicate_shape,
                                  nn_prior_two_step_support,
                                  nn_prior_two_step_support_min,
                                  nn_prior_two_step_cap_floor,
                                  cohort_transition_version,
                                  cohort_contextual_apply_to,
                                  cohort_context_keep_baseline_when_sparse,
                                  cohort_context_lambda_sparse_unknown,
                                  nboot,
                                  n0,
                                  nb,
                                  benchmark_seed,
                                  correct_efflux,
                                  force_refit) {
  if (is.null(parameter_spec_tbl) || !nrow(parameter_spec_tbl)) {
    return(tibble::tibble())
  }

  tidyr::crossing(
    input_index_tbl %>% dplyr::select(patient_id, input_rds, input_row_count_minobs5, input_row_count),
    minobs = minobs_values,
    pm = pm_values
  ) %>%
    tidyr::crossing(parameter_spec_tbl) %>%
    dplyr::mutate(
      nn_prior_grid_n = nn_prior_grid_n,
      nn_prior_fit_subset = nn_prior_fit_subset,
      nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
      nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
      nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
      nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
      nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
      nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
      nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
      nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
      nn_prior_two_step_support = nn_prior_two_step_support,
      nn_prior_two_step_support_min = nn_prior_two_step_support_min,
      nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
      cohort_transition_version = dplyr::if_else(
        nn_prior == "cohort_transition",
        as.character(cohort_transition_version),
        NA_character_
      ),
      cohort_contextual_apply_to = dplyr::if_else(
        nn_prior == "cohort_transition",
        as.character(cohort_contextual_apply_to),
        NA_character_
      ),
      cohort_context_keep_baseline_when_sparse = dplyr::if_else(
        nn_prior == "cohort_transition",
        as.logical(cohort_context_keep_baseline_when_sparse),
        NA
      ),
      cohort_context_lambda_sparse_unknown = dplyr::if_else(
        nn_prior == "cohort_transition",
        as.numeric(cohort_context_lambda_sparse_unknown),
        NA_real_
      ),
      outdir = purrr::pmap_chr(
        list(patient_id, minobs, pm, parameter_label),
        ~ task_outdir_parameter(fit_root, ..1, ..2, ..3, ..4)
      ),
      nboot = nboot,
      n0 = n0,
      nb = nb,
      benchmark_seed = benchmark_seed,
      correct_efflux = correct_efflux,
      force_refit = force_refit
    ) %>%
    dplyr::arrange(
      dplyr::desc(input_row_count_minobs5),
      dplyr::desc(input_row_count),
      factor(patient_id, levels = sort_pid_levels(patient_id)),
      minobs,
      pm,
      factor(parameter_label, levels = parameter_spec_tbl$parameter_label)
    )
}

run_task_table_parallel <- function(task_tbl, n_cores, diploid_state) {
  if (!nrow(task_tbl)) {
    return(tibble::tibble())
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    res <- parallel::mclapply(
      seq_len(nrow(task_tbl)),
      function(i) {
        rr <- task_tbl[i, , drop = FALSE]
        run_alfak_fit(
          patient_id = rr$patient_id,
          input_rds = rr$input_rds,
          outdir = rr$outdir,
          minobs = rr$minobs,
          pm = rr$pm,
          nboot = rr$nboot,
          n0 = rr$n0,
          nb = rr$nb,
          benchmark_seed = rr$benchmark_seed,
          parameter_label = rr$parameter_label,
          diploid_state = diploid_state,
          correct_efflux = rr$correct_efflux,
          nn_prior = rr$nn_prior,
          nn_prior_grid_n = rr$nn_prior_grid_n,
          nn_prior_fit_subset = rr$nn_prior_fit_subset,
          nn_prior_zero_exposure_quantile = rr$nn_prior_zero_exposure_quantile,
          nn_prior_zero_weight_scale = rr$nn_prior_zero_weight_scale,
          nn_prior_zero_weight_cap_ratio = rr$nn_prior_zero_weight_cap_ratio,
          nn_prior_zero_birth_fallback_weight = rr$nn_prior_zero_birth_fallback_weight,
          nn_prior_zero_birth_child_floor = rr$nn_prior_zero_birth_child_floor,
          nn_prior_zero_birth_child_shape = rr$nn_prior_zero_birth_child_shape,
          nn_prior_zero_birth_replicate_floor = rr$nn_prior_zero_birth_replicate_floor,
          nn_prior_zero_birth_replicate_shape = rr$nn_prior_zero_birth_replicate_shape,
          nn_prior_two_step_support = rr$nn_prior_two_step_support,
          nn_prior_two_step_support_min = rr$nn_prior_two_step_support_min,
          nn_prior_two_step_cap_floor = rr$nn_prior_two_step_cap_floor,
          cohort_transition_version = if ("cohort_transition_version" %in% names(rr)) rr$cohort_transition_version else NA_character_,
          cohort_contextual_apply_to = if ("cohort_contextual_apply_to" %in% names(rr)) rr$cohort_contextual_apply_to else NA_character_,
          cohort_context_keep_baseline_when_sparse = if ("cohort_context_keep_baseline_when_sparse" %in% names(rr)) rr$cohort_context_keep_baseline_when_sparse else NA,
          cohort_context_lambda_sparse_unknown = if ("cohort_context_lambda_sparse_unknown" %in% names(rr)) rr$cohort_context_lambda_sparse_unknown else NA_real_,
          force_refit = rr$force_refit
        )
      },
      mc.cores = n_cores,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    res <- lapply(seq_len(nrow(task_tbl)), function(i) {
      rr <- task_tbl[i, , drop = FALSE]
      run_alfak_fit(
        patient_id = rr$patient_id,
        input_rds = rr$input_rds,
        outdir = rr$outdir,
        minobs = rr$minobs,
        pm = rr$pm,
        nboot = rr$nboot,
        n0 = rr$n0,
        nb = rr$nb,
        benchmark_seed = rr$benchmark_seed,
        parameter_label = rr$parameter_label,
        diploid_state = diploid_state,
        correct_efflux = rr$correct_efflux,
        nn_prior = rr$nn_prior,
        nn_prior_grid_n = rr$nn_prior_grid_n,
        nn_prior_fit_subset = rr$nn_prior_fit_subset,
        nn_prior_zero_exposure_quantile = rr$nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = rr$nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = rr$nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_birth_fallback_weight = rr$nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = rr$nn_prior_zero_birth_child_floor,
        nn_prior_zero_birth_child_shape = rr$nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = rr$nn_prior_zero_birth_replicate_floor,
        nn_prior_zero_birth_replicate_shape = rr$nn_prior_zero_birth_replicate_shape,
        nn_prior_two_step_support = rr$nn_prior_two_step_support,
        nn_prior_two_step_support_min = rr$nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = rr$nn_prior_two_step_cap_floor,
        cohort_transition_version = if ("cohort_transition_version" %in% names(rr)) rr$cohort_transition_version else NA_character_,
        cohort_contextual_apply_to = if ("cohort_contextual_apply_to" %in% names(rr)) rr$cohort_contextual_apply_to else NA_character_,
        cohort_context_keep_baseline_when_sparse = if ("cohort_context_keep_baseline_when_sparse" %in% names(rr)) rr$cohort_context_keep_baseline_when_sparse else NA,
        cohort_context_lambda_sparse_unknown = if ("cohort_context_lambda_sparse_unknown" %in% names(rr)) rr$cohort_context_lambda_sparse_unknown else NA_real_,
        force_refit = rr$force_refit
      )
    })
  }

  dplyr::bind_rows(res)
}

build_task_error_rows <- function(task_tbl, error_message) {
  if (is.null(task_tbl) || !nrow(task_tbl)) {
    return(tibble::tibble())
  }

  error_message <- as.character(error_message)
  if (length(error_message) == 1L) {
    error_message <- rep(error_message, nrow(task_tbl))
  }
  if (length(error_message) != nrow(task_tbl)) {
    stop("`error_message` must have length 1 or match `task_tbl` rows.", call. = FALSE)
  }

  dplyr::bind_rows(lapply(seq_len(nrow(task_tbl)), function(i) {
    rr <- task_tbl[i, , drop = FALSE]
    outdir <- as.character(rr$outdir)
    tibble::tibble(
      patient_id = as.character(rr$patient_id),
      outdir = outdir,
      pm = as.numeric(rr$pm),
      pm_label = pm_to_label(as.numeric(rr$pm)),
      minobs = as.integer(rr$minobs),
      benchmark_seed = as.integer(rr$benchmark_seed),
      parameter_label = as.character(rr$parameter_label),
      nn_prior = as.character(rr$nn_prior),
      nn_prior_grid_n = as.integer(rr$nn_prior_grid_n),
      nn_prior_fit_subset = as.character(rr$nn_prior_fit_subset),
      nn_prior_zero_exposure_quantile = as.numeric(rr$nn_prior_zero_exposure_quantile),
      nn_prior_zero_weight_scale = as.numeric(rr$nn_prior_zero_weight_scale),
      nn_prior_zero_weight_cap_ratio = as.numeric(rr$nn_prior_zero_weight_cap_ratio),
      nn_prior_zero_birth_fallback_weight = as.numeric(rr$nn_prior_zero_birth_fallback_weight),
      nn_prior_zero_birth_child_floor = as.numeric(rr$nn_prior_zero_birth_child_floor),
      nn_prior_zero_birth_child_shape = as.numeric(rr$nn_prior_zero_birth_child_shape),
      nn_prior_zero_birth_replicate_floor = as.numeric(rr$nn_prior_zero_birth_replicate_floor),
      nn_prior_zero_birth_replicate_shape = as.numeric(rr$nn_prior_zero_birth_replicate_shape),
      nn_prior_two_step_support = as.character(rr$nn_prior_two_step_support),
      nn_prior_two_step_support_min = as.numeric(rr$nn_prior_two_step_support_min),
      nn_prior_two_step_cap_floor = as.numeric(rr$nn_prior_two_step_cap_floor),
      cohort_transition_version = if ("cohort_transition_version" %in% names(rr)) as.character(rr$cohort_transition_version) else NA_character_,
      cohort_contextual_apply_to = if ("cohort_contextual_apply_to" %in% names(rr)) as.character(rr$cohort_contextual_apply_to) else NA_character_,
      cohort_context_keep_baseline_when_sparse = if ("cohort_context_keep_baseline_when_sparse" %in% names(rr)) as.logical(rr$cohort_context_keep_baseline_when_sparse) else NA,
      cohort_context_lambda_sparse_unknown = if ("cohort_context_lambda_sparse_unknown" %in% names(rr)) as.numeric(rr$cohort_context_lambda_sparse_unknown) else NA_real_,
      status = "error",
      cached = FALSE,
      error_message = error_message[[i]],
      elapsed_sec = NA_real_,
      warning_count = 0L,
      lambda_endpoint_warning_count = 0L,
      warning_messages = NA_character_,
      xval_r2 = NA_real_,
      xval_cor = NA_real_,
      xval_rmse = NA_real_,
      xval_mae = NA_real_,
      n_xval = NA_integer_,
      landscape_path = if (file.exists(file.path(outdir, "landscape.Rds"))) file.path(outdir, "landscape.Rds") else NA_character_,
      bootstrap_path = if (file.exists(file.path(outdir, "bootstrap_res.Rds"))) file.path(outdir, "bootstrap_res.Rds") else NA_character_,
      posterior_path = if (file.exists(file.path(outdir, "landscape_posterior_samples.Rds"))) file.path(outdir, "landscape_posterior_samples.Rds") else NA_character_,
      xval_path = if (file.exists(file.path(outdir, "xval.Rds"))) file.path(outdir, "xval.Rds") else NA_character_
    )
  }))
}

read_benchmark_patient_input <- function(input_rds, diploid_state) {
  yi <- readRDS(input_rds)
  yi$x <- as.data.frame(yi$x)
  if (diploid_state %in% rownames(yi$x)) {
    yi$x <- yi$x[rownames(yi$x) != diploid_state, , drop = FALSE]
  }
  if (!nrow(yi$x)) {
    stop("All rows were filtered out for input: ", input_rds, call. = FALSE)
  }
  yi
}

build_cohort_patient_list <- function(task_tbl, diploid_state) {
  patients <- lapply(task_tbl$input_rds, read_benchmark_patient_input, diploid_state = diploid_state)
  names(patients) <- as.character(task_tbl$patient_id)
  patients
}

lookup_cohort_refit_errors <- function(cohort_outdir, task_tbl) {
  status_tbl <- safe_read_rds(file.path(cohort_outdir, "cohort_transition_refit_status.Rds"))
  if (is.null(status_tbl) || !is.data.frame(status_tbl) || !"patient_id" %in% names(status_tbl)) {
    return(rep("cohort_transition did not produce complete outputs.", nrow(task_tbl)))
  }

  status_tbl <- tibble::as_tibble(status_tbl) %>%
    dplyr::mutate(patient_id = as.character(patient_id)) %>%
    dplyr::select(patient_id, dplyr::any_of("error_message"))
  if (!"error_message" %in% names(status_tbl)) {
    status_tbl$error_message <- NA_character_
  }

  task_tbl %>%
    dplyr::mutate(.row_id = dplyr::row_number()) %>%
    dplyr::left_join(status_tbl, by = "patient_id") %>%
    dplyr::arrange(.row_id) %>%
    dplyr::transmute(
      error_message = dplyr::if_else(
        !is.na(error_message) & nzchar(error_message),
        as.character(error_message),
        "cohort_transition did not produce complete outputs."
      )
    ) %>%
    dplyr::pull(error_message)
}

run_cohort_transition_task_group <- function(task_tbl,
                                             base_results_tbl,
                                             fit_root,
                                             nboot,
                                             n0,
                                             nb,
                                             correct_efflux,
                                             n_cores = 1L,
                                             diploid_state,
                                             force_refit) {
  if (!nrow(task_tbl)) {
    return(tibble::tibble())
  }

  minobs_value <- as.integer(task_tbl$minobs[[1]])
  pm_value <- as.numeric(task_tbl$pm[[1]])
  cohort_transition_version <- if ("cohort_transition_version" %in% names(task_tbl)) {
    as.character(task_tbl$cohort_transition_version[[1]])
  } else {
    "contextual"
  }
  if (!length(cohort_transition_version) || is.na(cohort_transition_version) || !nzchar(cohort_transition_version)) {
    cohort_transition_version <- "contextual"
  }
  cohort_contextual_apply_to <- if ("cohort_contextual_apply_to" %in% names(task_tbl)) {
    as.character(task_tbl$cohort_contextual_apply_to[[1]])
  } else {
    "all"
  }
  if (!length(cohort_contextual_apply_to) || is.na(cohort_contextual_apply_to) || !nzchar(cohort_contextual_apply_to)) {
    cohort_contextual_apply_to <- "all"
  }
  cohort_context_keep_baseline_when_sparse <- if ("cohort_context_keep_baseline_when_sparse" %in% names(task_tbl)) {
    as.logical(task_tbl$cohort_context_keep_baseline_when_sparse[[1]])
  } else {
    TRUE
  }
  if (!length(cohort_context_keep_baseline_when_sparse) || is.na(cohort_context_keep_baseline_when_sparse[[1L]])) {
    cohort_context_keep_baseline_when_sparse <- TRUE
  } else {
    cohort_context_keep_baseline_when_sparse <- cohort_context_keep_baseline_when_sparse[[1L]]
  }
  cohort_context_lambda_sparse_unknown <- if ("cohort_context_lambda_sparse_unknown" %in% names(task_tbl)) {
    suppressWarnings(as.numeric(task_tbl$cohort_context_lambda_sparse_unknown[[1]]))
  } else {
    0
  }
  if (!length(cohort_context_lambda_sparse_unknown) ||
      is.na(cohort_context_lambda_sparse_unknown[[1L]]) ||
      !is.finite(cohort_context_lambda_sparse_unknown[[1L]]) ||
      cohort_context_lambda_sparse_unknown[[1L]] < 0) {
    cohort_context_lambda_sparse_unknown <- 0
  } else {
    cohort_context_lambda_sparse_unknown <- cohort_context_lambda_sparse_unknown[[1L]]
  }

  if (!force_refit) {
    adopted_tbl <- adopt_existing_task_outputs(task_tbl)
    if (nrow(adopted_tbl) == nrow(task_tbl)) {
      return(adopted_tbl)
    }
  }

  base_result_ok_tbl <- if (
    !is.null(base_results_tbl) &&
      nrow(base_results_tbl) &&
      all(c("parameter_label", "status", "patient_id", "minobs", "pm") %in% names(base_results_tbl))
  ) {
    base_results_tbl %>%
      dplyr::filter(
        parameter_label == "nn_prior_empirical_two_shell",
        status == "ok",
        minobs == minobs_value,
        abs(pm - pm_value) <= max(1e-12, abs(pm_value) * 1e-8)
      ) %>%
      dplyr::select(patient_id, minobs, pm)
  } else {
    tibble::tibble(patient_id = character(), minobs = integer(), pm = numeric())
  }
  base_disk_ok_tbl <- task_tbl %>%
    dplyr::transmute(
      patient_id,
      minobs,
      pm,
      base_outdir = purrr::pmap_chr(
        list(patient_id, minobs, pm),
        ~ task_outdir_parameter(fit_root, ..1, ..2, ..3, "nn_prior_empirical_two_shell")
      )
    ) %>%
    dplyr::filter(vapply(base_outdir, has_complete_alfak_outputs, logical(1))) %>%
    dplyr::select(patient_id, minobs, pm)
  base_ok_tbl <- dplyr::bind_rows(base_result_ok_tbl, base_disk_ok_tbl) %>%
    dplyr::distinct(patient_id, minobs, pm)

  eligible_task_tbl <- task_tbl %>%
    dplyr::semi_join(base_ok_tbl, by = c("patient_id", "minobs", "pm"))
  missing_base_task_tbl <- task_tbl %>%
    dplyr::anti_join(base_ok_tbl, by = c("patient_id", "minobs", "pm"))

  missing_base_rows <- build_task_error_rows(
    missing_base_task_tbl,
    "Missing successful nn_prior_empirical_two_shell prerequisite for cohort_transition."
  )
  if (!nrow(eligible_task_tbl)) {
    return(missing_base_rows)
  }
  cohort_refit_cores <- suppressWarnings(as.integer(n_cores))
  if (!length(cohort_refit_cores) || is.na(cohort_refit_cores[[1L]]) || cohort_refit_cores[[1L]] < 1L) {
    cohort_refit_cores <- 1L
  } else {
    cohort_refit_cores <- min(cohort_refit_cores[[1L]], max(1L, nrow(eligible_task_tbl)))
  }

  cohort_outdir <- dirname(as.character(eligible_task_tbl$outdir[[1]]))
  two_shell_root <- file.path(fit_root, "nn_prior_empirical_two_shell")
  cohort_started_at <- Sys.time()
  cohort_error <- NULL

  patients <- tryCatch(
    build_cohort_patient_list(eligible_task_tbl, diploid_state = diploid_state),
    error = function(e) {
      cohort_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(patients)) {
    return(dplyr::bind_rows(
      missing_base_rows,
      build_task_error_rows(eligible_task_tbl, cohort_error)
    ))
  }

  alfak_log(
    "ALFA-K cohort_transition start: minobs=", minobs_value,
    " | pm=", pm_to_label(pm_value),
    " | version=", cohort_transition_version,
    " | contextual_apply_to=", cohort_contextual_apply_to,
    " | keep_sparse_baseline=", cohort_context_keep_baseline_when_sparse,
    " | sparse_unknown_lambda=", signif(cohort_context_lambda_sparse_unknown, 4),
    " | patients=", nrow(eligible_task_tbl),
    " | refit_cores=", cohort_refit_cores
  )
  cohort_result <- tryCatch(
    {
      set.seed(as.integer(eligible_task_tbl$benchmark_seed[[1]]))
      alfakR::alfak_cohort_transition(
        patients = patients,
        patient_ids = as.character(eligible_task_tbl$patient_id),
        outdir = cohort_outdir,
        two_shell_root = two_shell_root,
        two_shell_pm = pm_value,
        two_shell_minobs = minobs_value,
        two_shell_pm_tag = paste0("pm_", pm_to_label(pm_value)),
        two_shell_minobs_tag = paste0("MINOBS_", minobs_value),
        reuse_two_shell = TRUE,
        rerun_missing_two_shell = FALSE,
        rerun_corrupt_two_shell = FALSE,
        two_shell_integrity_check = "strict",
        base_nn_prior = "empirical_two_shell",
        cohort_transition_version = cohort_transition_version,
        cohort_contextual_apply_to = cohort_contextual_apply_to,
        cohort_context_keep_baseline_when_sparse = cohort_context_keep_baseline_when_sparse,
        cohort_context_lambda_sparse_unknown = cohort_context_lambda_sparse_unknown,
        minobs = minobs_value,
        nboot = nboot,
        n0 = n0,
        nb = nb,
        pm = pm_value,
        passage_times = NULL,
        allow_noninteger_counts = FALSE,
        correct_efflux = correct_efflux,
        cohort_refit_cores = cohort_refit_cores,
        cohort_refit_seed = as.integer(eligible_task_tbl$benchmark_seed[[1]]),
        nn_prior_grid_n = as.integer(eligible_task_tbl$nn_prior_grid_n[[1]]),
        nn_prior_fit_subset = as.character(eligible_task_tbl$nn_prior_fit_subset[[1]]),
        nn_prior_zero_exposure_quantile = as.numeric(eligible_task_tbl$nn_prior_zero_exposure_quantile[[1]]),
        nn_prior_zero_weight_scale = as.numeric(eligible_task_tbl$nn_prior_zero_weight_scale[[1]]),
        nn_prior_zero_weight_cap_ratio = if (is.na(eligible_task_tbl$nn_prior_zero_weight_cap_ratio[[1]])) NULL else as.numeric(eligible_task_tbl$nn_prior_zero_weight_cap_ratio[[1]]),
        nn_prior_zero_birth_fallback_weight = if (is.na(eligible_task_tbl$nn_prior_zero_birth_fallback_weight[[1]])) NULL else as.numeric(eligible_task_tbl$nn_prior_zero_birth_fallback_weight[[1]]),
        nn_prior_zero_birth_child_floor = as.numeric(eligible_task_tbl$nn_prior_zero_birth_child_floor[[1]]),
        nn_prior_zero_birth_child_shape = as.numeric(eligible_task_tbl$nn_prior_zero_birth_child_shape[[1]]),
        nn_prior_zero_birth_replicate_floor = as.numeric(eligible_task_tbl$nn_prior_zero_birth_replicate_floor[[1]]),
        nn_prior_zero_birth_replicate_shape = as.numeric(eligible_task_tbl$nn_prior_zero_birth_replicate_shape[[1]]),
        nn_prior_two_step_support = as.character(eligible_task_tbl$nn_prior_two_step_support[[1]]),
        nn_prior_two_step_support_min = as.numeric(eligible_task_tbl$nn_prior_two_step_support_min[[1]]),
        nn_prior_two_step_cap_floor = as.numeric(eligible_task_tbl$nn_prior_two_step_cap_floor[[1]])
      )
    },
    error = function(e) {
      cohort_error <<- conditionMessage(e)
      NULL
    }
  )
  elapsed_sec <- as.numeric(difftime(Sys.time(), cohort_started_at, units = "secs"))

  if (is.null(cohort_result)) {
    alfak_log(
      "ALFA-K cohort_transition error: minobs=", minobs_value,
      " | pm=", pm_to_label(pm_value),
      " | version=", cohort_transition_version,
      " | contextual_apply_to=", cohort_contextual_apply_to,
      " | keep_sparse_baseline=", cohort_context_keep_baseline_when_sparse,
      " | sparse_unknown_lambda=", signif(cohort_context_lambda_sparse_unknown, 4),
      " | ", cohort_error
    )
    return(dplyr::bind_rows(
      missing_base_rows,
      build_task_error_rows(eligible_task_tbl, cohort_error)
    ))
  }

  adopted_tbl <- adopt_existing_task_outputs(eligible_task_tbl)
  if (nrow(adopted_tbl)) {
    adopted_tbl$elapsed_sec <- elapsed_sec / max(1L, nrow(adopted_tbl))
  }
  adopted_key_tbl <- if (
    nrow(adopted_tbl) &&
      all(c("patient_id", "minobs", "pm", "parameter_label", "nn_prior") %in% names(adopted_tbl))
  ) {
    adopted_tbl %>%
      dplyr::select(patient_id, minobs, pm, parameter_label, nn_prior)
  } else {
    tibble::tibble(
      patient_id = character(),
      minobs = integer(),
      pm = numeric(),
      parameter_label = character(),
      nn_prior = character()
    )
  }

  incomplete_task_tbl <- eligible_task_tbl %>%
    dplyr::anti_join(
      adopted_key_tbl,
      by = c("patient_id", "minobs", "pm", "parameter_label", "nn_prior")
    )
  incomplete_rows <- build_task_error_rows(
    incomplete_task_tbl,
    lookup_cohort_refit_errors(cohort_outdir, incomplete_task_tbl)
  )

  alfak_log(
    "ALFA-K cohort_transition done: minobs=", minobs_value,
    " | pm=", pm_to_label(pm_value),
    " | version=", cohort_transition_version,
    " | contextual_apply_to=", cohort_contextual_apply_to,
    " | keep_sparse_baseline=", cohort_context_keep_baseline_when_sparse,
    " | sparse_unknown_lambda=", signif(cohort_context_lambda_sparse_unknown, 4),
    " | ok=", nrow(adopted_tbl),
    " | error=", nrow(incomplete_rows) + nrow(missing_base_rows)
  )

  dplyr::bind_rows(missing_base_rows, adopted_tbl, incomplete_rows)
}

run_cohort_transition_tasks <- function(task_tbl,
                                        base_results_tbl,
                                        fit_root,
                                        nboot,
                                        n0,
                                        nb,
                                        correct_efflux,
                                        n_cores = 1L,
                                        cohort_transition_version = "contextual",
                                        cohort_contextual_apply_to = "all",
                                        cohort_context_keep_baseline_when_sparse = TRUE,
                                        cohort_context_lambda_sparse_unknown = 0,
                                        diploid_state,
                                        force_refit) {
  cohort_transition_version_requested <- as.character(cohort_transition_version)
  if (!length(cohort_transition_version_requested) ||
      is.na(cohort_transition_version_requested[[1L]]) ||
      !nzchar(cohort_transition_version_requested[[1L]])) {
    cohort_transition_version_requested <- "contextual"
  } else {
    cohort_transition_version_requested <- cohort_transition_version_requested[[1L]]
  }
  cohort_contextual_apply_to_requested <- as.character(cohort_contextual_apply_to)
  if (!length(cohort_contextual_apply_to_requested) ||
      is.na(cohort_contextual_apply_to_requested[[1L]]) ||
      !nzchar(cohort_contextual_apply_to_requested[[1L]])) {
    cohort_contextual_apply_to_requested <- "all"
  } else {
    cohort_contextual_apply_to_requested <- cohort_contextual_apply_to_requested[[1L]]
  }
  cohort_context_keep_baseline_when_sparse_requested <- as.logical(cohort_context_keep_baseline_when_sparse)
  if (!length(cohort_context_keep_baseline_when_sparse_requested) ||
      is.na(cohort_context_keep_baseline_when_sparse_requested[[1L]])) {
    cohort_context_keep_baseline_when_sparse_requested <- TRUE
  } else {
    cohort_context_keep_baseline_when_sparse_requested <- cohort_context_keep_baseline_when_sparse_requested[[1L]]
  }
  cohort_context_lambda_sparse_unknown_requested <- suppressWarnings(as.numeric(cohort_context_lambda_sparse_unknown))
  if (!length(cohort_context_lambda_sparse_unknown_requested) ||
      is.na(cohort_context_lambda_sparse_unknown_requested[[1L]]) ||
      !is.finite(cohort_context_lambda_sparse_unknown_requested[[1L]]) ||
      cohort_context_lambda_sparse_unknown_requested[[1L]] < 0) {
    cohort_context_lambda_sparse_unknown_requested <- 0
  } else {
    cohort_context_lambda_sparse_unknown_requested <- cohort_context_lambda_sparse_unknown_requested[[1L]]
  }
  cohort_task_tbl <- task_tbl %>%
    dplyr::filter(nn_prior == "cohort_transition")
  if (!"cohort_transition_version" %in% names(cohort_task_tbl)) {
    cohort_task_tbl$cohort_transition_version <- cohort_transition_version_requested
  }
  if (!"cohort_contextual_apply_to" %in% names(cohort_task_tbl)) {
    cohort_task_tbl$cohort_contextual_apply_to <- cohort_contextual_apply_to_requested
  }
  if (!"cohort_context_keep_baseline_when_sparse" %in% names(cohort_task_tbl)) {
    cohort_task_tbl$cohort_context_keep_baseline_when_sparse <- cohort_context_keep_baseline_when_sparse_requested
  }
  if (!"cohort_context_lambda_sparse_unknown" %in% names(cohort_task_tbl)) {
    cohort_task_tbl$cohort_context_lambda_sparse_unknown <- cohort_context_lambda_sparse_unknown_requested
  }
  cohort_task_tbl <- cohort_task_tbl %>%
    dplyr::mutate(cohort_transition_version = dplyr::if_else(
      is.na(.data$cohort_transition_version) | !nzchar(.data$cohort_transition_version),
      cohort_transition_version_requested,
      as.character(.data$cohort_transition_version)
    ), cohort_contextual_apply_to = dplyr::if_else(
      is.na(.data$cohort_contextual_apply_to) | !nzchar(.data$cohort_contextual_apply_to),
      cohort_contextual_apply_to_requested,
      as.character(.data$cohort_contextual_apply_to)
    ), cohort_context_keep_baseline_when_sparse = dplyr::if_else(
      is.na(.data$cohort_context_keep_baseline_when_sparse),
      cohort_context_keep_baseline_when_sparse_requested,
      as.logical(.data$cohort_context_keep_baseline_when_sparse)
    ), cohort_context_lambda_sparse_unknown = dplyr::if_else(
      is.na(.data$cohort_context_lambda_sparse_unknown),
      cohort_context_lambda_sparse_unknown_requested,
      as.numeric(.data$cohort_context_lambda_sparse_unknown)
    )) %>%
    dplyr::arrange(minobs, pm, factor(patient_id, levels = sort_pid_levels(patient_id)))
  if (!nrow(cohort_task_tbl)) {
    return(tibble::tibble())
  }

  cohort_task_tbl %>%
    dplyr::group_split(pm, minobs, .keep = TRUE) %>%
    lapply(
      run_cohort_transition_task_group,
      base_results_tbl = base_results_tbl,
      fit_root = fit_root,
      nboot = nboot,
      n0 = n0,
      nb = nb,
      correct_efflux = correct_efflux,
      n_cores = n_cores,
      diploid_state = diploid_state,
      force_refit = force_refit
    ) %>%
    dplyr::bind_rows()
}

adopt_existing_task_outputs <- function(task_tbl, fit_results_tbl = NULL) {
  if (is.null(task_tbl) || !nrow(task_tbl)) {
    return(if (is.null(fit_results_tbl)) tibble::tibble() else tibble::as_tibble(fit_results_tbl))
  }

  key_cols <- c("patient_id", "minobs", "pm", "parameter_label", "nn_prior")
  fit_results_tbl <- if (is.null(fit_results_tbl)) tibble::tibble() else tibble::as_tibble(fit_results_tbl)
  incompatible_task_key_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(task_tbl)), function(i) {
    rr <- task_tbl[i, , drop = FALSE]
    if (!identical(as.character(rr$nn_prior), "cohort_transition")) {
      return(NULL)
    }
    if (cohort_transition_output_matches(
      outdir = rr$outdir,
      nn_prior = as.character(rr$nn_prior),
      cohort_transition_version = if ("cohort_transition_version" %in% names(rr)) rr$cohort_transition_version else NA_character_,
      cohort_contextual_apply_to = if ("cohort_contextual_apply_to" %in% names(rr)) rr$cohort_contextual_apply_to else NA_character_,
      cohort_context_keep_baseline_when_sparse = if ("cohort_context_keep_baseline_when_sparse" %in% names(rr)) rr$cohort_context_keep_baseline_when_sparse else NA,
      cohort_context_lambda_sparse_unknown = if ("cohort_context_lambda_sparse_unknown" %in% names(rr)) rr$cohort_context_lambda_sparse_unknown else NA_real_
    )) {
      return(NULL)
    }
    rr[, key_cols, drop = FALSE]
  }))
  if (nrow(incompatible_task_key_tbl) && nrow(fit_results_tbl)) {
    fit_results_tbl <- fit_results_tbl %>%
      dplyr::anti_join(incompatible_task_key_tbl, by = key_cols)
  }

  adopted_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(task_tbl)), function(i) {
    rr <- task_tbl[i, , drop = FALSE]
    if (!has_complete_alfak_outputs(rr$outdir) ||
        !cohort_transition_output_matches(
          outdir = rr$outdir,
          nn_prior = as.character(rr$nn_prior),
          cohort_transition_version = if ("cohort_transition_version" %in% names(rr)) rr$cohort_transition_version else NA_character_,
          cohort_contextual_apply_to = if ("cohort_contextual_apply_to" %in% names(rr)) rr$cohort_contextual_apply_to else NA_character_,
          cohort_context_keep_baseline_when_sparse = if ("cohort_context_keep_baseline_when_sparse" %in% names(rr)) rr$cohort_context_keep_baseline_when_sparse else NA,
          cohort_context_lambda_sparse_unknown = if ("cohort_context_lambda_sparse_unknown" %in% names(rr)) rr$cohort_context_lambda_sparse_unknown else NA_real_
        )) {
      return(NULL)
    }

    tibble::as_tibble(build_existing_fit_row(
      patient_id = rr$patient_id,
      outdir = rr$outdir,
      pm = rr$pm,
      minobs = rr$minobs,
      benchmark_seed = rr$benchmark_seed,
      parameter_label = rr$parameter_label,
      nn_prior = rr$nn_prior,
      nn_prior_grid_n = rr$nn_prior_grid_n,
      nn_prior_fit_subset = rr$nn_prior_fit_subset,
      nn_prior_zero_exposure_quantile = rr$nn_prior_zero_exposure_quantile,
      nn_prior_zero_weight_scale = rr$nn_prior_zero_weight_scale,
      nn_prior_zero_weight_cap_ratio = rr$nn_prior_zero_weight_cap_ratio,
      nn_prior_zero_birth_fallback_weight = rr$nn_prior_zero_birth_fallback_weight,
      nn_prior_zero_birth_child_floor = rr$nn_prior_zero_birth_child_floor,
      nn_prior_zero_birth_child_shape = rr$nn_prior_zero_birth_child_shape,
      nn_prior_zero_birth_replicate_floor = rr$nn_prior_zero_birth_replicate_floor,
      nn_prior_zero_birth_replicate_shape = rr$nn_prior_zero_birth_replicate_shape,
      nn_prior_two_step_support = rr$nn_prior_two_step_support,
      nn_prior_two_step_support_min = rr$nn_prior_two_step_support_min,
      nn_prior_two_step_cap_floor = rr$nn_prior_two_step_cap_floor,
      cohort_transition_version = if ("cohort_transition_version" %in% names(rr)) rr$cohort_transition_version else NA_character_,
      cohort_contextual_apply_to = if ("cohort_contextual_apply_to" %in% names(rr)) rr$cohort_contextual_apply_to else NA_character_,
      cohort_context_keep_baseline_when_sparse = if ("cohort_context_keep_baseline_when_sparse" %in% names(rr)) rr$cohort_context_keep_baseline_when_sparse else NA,
      cohort_context_lambda_sparse_unknown = if ("cohort_context_lambda_sparse_unknown" %in% names(rr)) rr$cohort_context_lambda_sparse_unknown else NA_real_,
      warning_log_path = file.path(rr$outdir, "fit_warnings.log"),
      landscape_path = file.path(rr$outdir, "landscape.Rds"),
      bootstrap_path = file.path(rr$outdir, "bootstrap_res.Rds"),
      posterior_path = file.path(rr$outdir, "landscape_posterior_samples.Rds"),
      xval_path = file.path(rr$outdir, "xval.Rds")
    ))
  }))

  if (!nrow(adopted_tbl)) {
    return(fit_results_tbl)
  }

  if (!nrow(fit_results_tbl)) {
    return(adopted_tbl)
  }

  fit_results_tbl %>%
    dplyr::anti_join(adopted_tbl %>% dplyr::select(dplyr::all_of(key_cols)), by = key_cols) %>%
    dplyr::bind_rows(adopted_tbl)
}

summarize_fit_results <- function(fit_results_tbl, group_cols) {
  if (!nrow(fit_results_tbl)) {
    return(tibble::tibble())
  }

  fit_results_tbl %>%
    dplyr::group_by(dplyr::across(all_of(group_cols))) %>%
    dplyr::summarise(
      n_tasks = dplyr::n(),
      n_ok = sum(status == "ok", na.rm = TRUE),
      n_error = sum(status == "error", na.rm = TRUE),
      success_rate = n_ok / n_tasks,
      mean_xval_r2 = mean(xval_r2, na.rm = TRUE),
      median_xval_r2 = median(xval_r2, na.rm = TRUE),
      mean_elapsed_sec = mean(elapsed_sec, na.rm = TRUE),
      median_elapsed_sec = median(elapsed_sec, na.rm = TRUE),
      total_lambda_endpoint_warnings = sum(lambda_endpoint_warning_count, na.rm = TRUE),
      .groups = "drop"
    )
}
