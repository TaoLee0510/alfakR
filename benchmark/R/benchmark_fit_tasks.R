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

  if (!identical(nn_prior, "empirical_censored_weighted")) {
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
    if (identical(nn_prior, "empirical_censored_weighted")) paste0(
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

  alfak_log("ALFA-K start: ", task_tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
        force_refit = rr$force_refit
      )
    })
  }

  dplyr::bind_rows(res)
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
