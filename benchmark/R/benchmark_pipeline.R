benchmark_load_libraries <- function() {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(knitr)
    library(parallel)
    library(pkgload)
    library(purrr)
    library(readxl)
    library(tibble)
    library(tidyr)
  })
}

normalize_focus_pids <- function(params) {
  focus_raw <- NULL
  if (!is.null(params$focus_pids)) {
    focus_raw <- params$focus_pids
  } else if (!is.null(params$focus_pid)) {
    focus_raw <- params$focus_pid
  }

  focus_pids <- unique(trimws(as.character(unlist(focus_raw))))
  focus_pids <- focus_pids[nzchar(focus_pids)]
  if (!length(focus_pids)) {
    focus_pids <- "P2"
  }
  sort_pid_levels(focus_pids)
}

benchmark_param_scalar <- function(params, name, default = NA) {
  value <- params[[name]]
  if (is.null(value) || !length(value)) {
    return(default)
  }
  value[[1L]]
}

focus_table_stem <- function(ctx, focus_pid, stem) {
  file.path(ctx$tables_dir, paste0("focus_", focus_pid, "_", stem))
}

build_parameter_figure_paths <- function(ctx) {
  list(
    runtime = file.path(ctx$figures_dir, "runtime_by_parameter_label.png"),
    xval = file.path(ctx$figures_dir, "xval_by_parameter_label.png"),
    global_diff = file.path(ctx$figures_dir, "global_landscape_and_beneficial_diff.png"),
    global_state_diff = file.path(ctx$figures_dir, "global_landscape_state_diff.png"),
    retained_landscape_distribution = file.path(ctx$figures_dir, "retained_landscape_distribution_comparison.png"),
    retained_xval_scatter = file.path(ctx$figures_dir, "retained_xval_validation_scatter.png")
  )
}

build_input_figure_paths <- function(ctx) {
  list(
    fq_nn_counts = file.path(ctx$figures_dir, "input_fq_nn_counts_by_minobs.png"),
    fq_nn_up_prop = file.path(ctx$figures_dir, "input_fq_nn_up_proportion_by_minobs.png"),
    fq_group_nn_up_prop = file.path(ctx$figures_dir, "input_fq_group_nn_up_proportion_by_minobs.png")
  )
}

build_input_group_prop_heatmap_data <- function(input_fq_group_nn_summary_tbl, patient_ids, minobs_levels) {
  patient_ids <- sort_pid_levels(patient_ids)
  feature_tbl <- tidyr::expand_grid(
    fq_prop_direction = c("up", "down"),
    minobs = minobs_levels
  ) %>%
    dplyr::mutate(feature = paste(fq_prop_direction, minobs, sep = "__"))

  mat_tbl <- tidyr::expand_grid(
    patient_id = patient_ids,
    feature = feature_tbl$feature
  ) %>%
    dplyr::left_join(
      input_fq_group_nn_summary_tbl %>%
        dplyr::mutate(
          patient_id = as.character(patient_id),
          fq_prop_direction = as.character(fq_prop_direction),
          minobs = as.integer(minobs),
          feature = paste(fq_prop_direction, minobs, sep = "__")
        ) %>%
        dplyr::select(patient_id, feature, prop_group_nn_count_up),
      by = c("patient_id", "feature")
    ) %>%
    tidyr::pivot_wider(names_from = feature, values_from = prop_group_nn_count_up)

  panel_mat <- as.matrix(mat_tbl[, setdiff(names(mat_tbl), "patient_id"), drop = FALSE])
  rownames(panel_mat) <- mat_tbl$patient_id
  storage.mode(panel_mat) <- "numeric"

  panel_mat <- panel_mat[, feature_tbl$feature, drop = FALSE]
  colnames(panel_mat) <- as.character(feature_tbl$minobs)

  list(
    panel_mat = panel_mat,
    column_split = factor(feature_tbl$fq_prop_direction, levels = c("up", "down"), labels = c("fq proportion up", "fq proportion down"))
  )
}

build_focus_figure_paths <- function(ctx, focus_pid) {
  list(
    parity = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_landscape_parity_by_parameter_label.png")),
    beneficial_heatmap = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_beneficial_proportion_heatmap.png")),
    state_fitness = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_state_fitness_decomposition.png")),
    edge_delta = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_edge_delta_decomposition.png")),
    edge_beneficial = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_edge_beneficial_decomposition.png"))
  )
}

build_benchmark_context <- function(params, repo_dir = resolve_repo_dir()) {
  benchmark_load_libraries()
  pkgload::load_all(repo_dir, quiet = TRUE)

  benchmark_dir <- file.path(repo_dir, "benchmark")
  data_dir <- file.path(benchmark_dir, "data")
  meta_path <- file.path(data_dir, "meta_data.xlsx")
  results_dir <- file.path(benchmark_dir, "results")
  input_dir <- file.path(results_dir, "inputs")
  fit_dir <- file.path(results_dir, "fits")
  tables_dir <- file.path(results_dir, "tables")
  figures_dir <- file.path(results_dir, "figures")
  cache_dir <- file.path(results_dir, "cache")

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  patient_subset_use <- unique(trimws(as.character(unlist(params$patient_subset))))
  patient_subset_use <- patient_subset_use[nzchar(patient_subset_use)]
  if (!length(patient_subset_use)) {
    patient_subset_use <- NULL
  }

  minobs_values_use <- sort(unique(as.integer(unlist(params$minobs_values))))
  minobs_values_use <- minobs_values_use[is.finite(minobs_values_use) & minobs_values_use > 0L]
  if (!length(minobs_values_use)) {
    stop("No valid minobs values supplied.")
  }

  pm_values_use <- sort(unique(as.numeric(unlist(params$pm_values))))
  pm_values_use <- pm_values_use[is.finite(pm_values_use) & pm_values_use > 0]
  if (!length(pm_values_use)) {
    stop("No valid pm values supplied.")
  }

  allowed_nn_prior_values <- eval(formals(alfakR::alfak)$nn_prior)
  allowed_parameter_labels <- paste0("nn_prior_", allowed_nn_prior_values)
  parameter_labels_use <- unique(as.character(unlist(params$parameter_labels)))
  parameter_labels_use <- parameter_labels_use[nzchar(parameter_labels_use)]
  parameter_labels_use <- ifelse(
    parameter_labels_use %in% allowed_nn_prior_values,
    paste0("nn_prior_", parameter_labels_use),
    parameter_labels_use
  )
  if (!length(parameter_labels_use) || !all(parameter_labels_use %in% allowed_parameter_labels)) {
    stop(
      "parameter_labels must be chosen from: ",
      paste(allowed_parameter_labels, collapse = ", ")
    )
  }
  if ("nn_prior_cohort_transition" %in% parameter_labels_use) {
    parameter_labels_use <- c(
      parameter_labels_use[parameter_labels_use != "nn_prior_cohort_transition"],
      "nn_prior_cohort_transition"
    )
  }

  selected_grid_n_use <- suppressWarnings(as.integer(params$nn_prior_grid_n))
  if (!is.finite(selected_grid_n_use) || selected_grid_n_use < 3L) {
    stop("nn_prior_grid_n must be a finite integer >= 3.")
  }

  n_cores_use <- suppressWarnings(as.integer(params$n_cores))
  if (!is.finite(n_cores_use) || n_cores_use < 1L) {
    detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (!is.finite(detected_cores) || detected_cores < 1L) {
      detected_cores <- 1L
    }
    n_cores_use <- max(1L, detected_cores - 1L)
  }

  nboot_use <- as.integer(params$nboot)
  n0_use <- as.numeric(params$n0)
  nb_use <- as.numeric(params$nb)
  benchmark_seed_use <- suppressWarnings(as.integer(params$benchmark_seed))
  if (!is.finite(benchmark_seed_use)) {
    stop("benchmark_seed must be a finite integer.")
  }
  correct_efflux_use <- isTRUE(params$correct_efflux)
  nn_prior_fit_subset_use <- match.arg(as.character(params$nn_prior_fit_subset), c("hybrid", "all"))
  nn_prior_zero_exposure_quantile_use <- suppressWarnings(as.numeric(params$nn_prior_zero_exposure_quantile))
  if (!is.finite(nn_prior_zero_exposure_quantile_use) ||
      nn_prior_zero_exposure_quantile_use < 0 ||
      nn_prior_zero_exposure_quantile_use > 1) {
    stop("nn_prior_zero_exposure_quantile must be a finite number in [0, 1].")
  }
  nn_prior_zero_weight_scale_use <- suppressWarnings(as.numeric(params$nn_prior_zero_weight_scale))
  if (!is.finite(nn_prior_zero_weight_scale_use) ||
      nn_prior_zero_weight_scale_use < 0 ||
      nn_prior_zero_weight_scale_use > 1) {
    stop("nn_prior_zero_weight_scale must be a finite number in [0, 1].")
  }
  nn_prior_zero_weight_cap_ratio_raw <- params$nn_prior_zero_weight_cap_ratio
  if (is.null(nn_prior_zero_weight_cap_ratio_raw) ||
      (length(nn_prior_zero_weight_cap_ratio_raw) == 1L && is.na(nn_prior_zero_weight_cap_ratio_raw))) {
    nn_prior_zero_weight_cap_ratio_use <- NA_real_
  } else {
    nn_prior_zero_weight_cap_ratio_use <- suppressWarnings(as.numeric(nn_prior_zero_weight_cap_ratio_raw))
    if (!is.finite(nn_prior_zero_weight_cap_ratio_use) || nn_prior_zero_weight_cap_ratio_use < 0) {
      stop("nn_prior_zero_weight_cap_ratio must be NULL/NA or a finite non-negative number.")
    }
  }
  nn_prior_zero_birth_fallback_weight_raw <- params$nn_prior_zero_birth_fallback_weight
  if (is.null(nn_prior_zero_birth_fallback_weight_raw) ||
      (length(nn_prior_zero_birth_fallback_weight_raw) == 1L && is.na(nn_prior_zero_birth_fallback_weight_raw))) {
    nn_prior_zero_birth_fallback_weight_use <- NA_real_
  } else {
    nn_prior_zero_birth_fallback_weight_use <- suppressWarnings(as.numeric(nn_prior_zero_birth_fallback_weight_raw))
    if (!is.finite(nn_prior_zero_birth_fallback_weight_use) ||
        nn_prior_zero_birth_fallback_weight_use < 0 ||
        nn_prior_zero_birth_fallback_weight_use > 1) {
      stop("nn_prior_zero_birth_fallback_weight must be NULL/NA or a finite number in [0, 1].")
    }
  }
  nn_prior_zero_birth_child_floor_use <- suppressWarnings(as.numeric(params$nn_prior_zero_birth_child_floor))
  if (!is.finite(nn_prior_zero_birth_child_floor_use) ||
      nn_prior_zero_birth_child_floor_use < 0 ||
      nn_prior_zero_birth_child_floor_use > 1) {
    stop("nn_prior_zero_birth_child_floor must be a finite number in [0, 1].")
  }
  nn_prior_zero_birth_child_shape_use <- suppressWarnings(as.numeric(params$nn_prior_zero_birth_child_shape))
  if (!is.finite(nn_prior_zero_birth_child_shape_use) ||
      nn_prior_zero_birth_child_shape_use < 0) {
    stop("nn_prior_zero_birth_child_shape must be a finite non-negative number.")
  }
  nn_prior_zero_birth_replicate_floor_use <- suppressWarnings(as.numeric(params$nn_prior_zero_birth_replicate_floor))
  if (!is.finite(nn_prior_zero_birth_replicate_floor_use) ||
      nn_prior_zero_birth_replicate_floor_use < 0 ||
      nn_prior_zero_birth_replicate_floor_use > 1) {
    stop("nn_prior_zero_birth_replicate_floor must be a finite number in [0, 1].")
  }
  nn_prior_zero_birth_replicate_shape_use <- suppressWarnings(as.numeric(params$nn_prior_zero_birth_replicate_shape))
  if (!is.finite(nn_prior_zero_birth_replicate_shape_use) ||
      nn_prior_zero_birth_replicate_shape_use < 0) {
    stop("nn_prior_zero_birth_replicate_shape must be a finite non-negative number.")
  }
  nn_prior_two_step_support_use <- match.arg(
    as.character(params$nn_prior_two_step_support),
    c("none", "rescue")
  )
  nn_prior_two_step_support_min_use <- suppressWarnings(as.numeric(params$nn_prior_two_step_support_min))
  if (!is.finite(nn_prior_two_step_support_min_use) ||
      nn_prior_two_step_support_min_use < 0 ||
      nn_prior_two_step_support_min_use > 1) {
    stop("nn_prior_two_step_support_min must be a finite number in [0, 1].")
  }
  nn_prior_two_step_cap_floor_use <- suppressWarnings(as.numeric(params$nn_prior_two_step_cap_floor))
  if (!is.finite(nn_prior_two_step_cap_floor_use) ||
      nn_prior_two_step_cap_floor_use < 0 ||
      nn_prior_two_step_cap_floor_use > 1) {
    stop("nn_prior_two_step_cap_floor must be a finite number in [0, 1].")
  }
  cohort_transition_version_use <- if (is.null(params$cohort_transition_version)) {
    "contextual"
  } else {
    match.arg(as.character(params$cohort_transition_version), c("contextual", "v2", "v1"))
  }
  cohort_contextual_apply_to_use <- if (is.null(params$cohort_contextual_apply_to)) {
    "all"
  } else {
    match.arg(as.character(params$cohort_contextual_apply_to), c("zero_only", "low_information", "all"))
  }
  cohort_context_keep_baseline_when_sparse_use <- if (is.null(params$cohort_context_keep_baseline_when_sparse)) {
    TRUE
  } else {
    isTRUE(params$cohort_context_keep_baseline_when_sparse)
  }
  cohort_context_lambda_sparse_unknown_use <- if (is.null(params$cohort_context_lambda_sparse_unknown)) {
    0
  } else {
    suppressWarnings(as.numeric(params$cohort_context_lambda_sparse_unknown))
  }
  if (!is.finite(cohort_context_lambda_sparse_unknown_use) || cohort_context_lambda_sparse_unknown_use < 0) {
    stop("cohort_context_lambda_sparse_unknown must be a non-negative finite number.")
  }

  force_refit_use <- isTRUE(params$force_refit)
  rebuild_inputs_use <- isTRUE(params$rebuild_inputs)
  run_benchmark_use <- isTRUE(params$run_benchmark)
  include_posterior_use <- isTRUE(params$include_posterior_comparison)
  render_figures_use <- isTRUE(params$render_figures)
  top_shift_n_use <- max(1L, as.integer(params$top_shift_n))
  run_nn_identifiability_use <- if (is.null(params$run_nn_identifiability)) TRUE else isTRUE(params$run_nn_identifiability)
  run_nn_stability_use <- if (is.null(params$run_nn_stability)) TRUE else isTRUE(params$run_nn_stability)
  run_nn_holdout_use <- isTRUE(params$run_nn_holdout)
  nn_holdout_repeats_use <- max(1L, suppressWarnings(as.integer(params$nn_holdout_repeats)))
  if (!is.finite(nn_holdout_repeats_use)) {
    nn_holdout_repeats_use <- 5L
  }
  nn_holdout_fraction_use <- suppressWarnings(as.numeric(params$nn_holdout_fraction))
  if (!is.finite(nn_holdout_fraction_use) || nn_holdout_fraction_use <= 0 || nn_holdout_fraction_use > 1) {
    nn_holdout_fraction_use <- 0.25
  }
  nn_holdout_min_count_use <- max(1L, suppressWarnings(as.integer(params$nn_holdout_min_count)))
  if (!is.finite(nn_holdout_min_count_use)) {
    nn_holdout_min_count_use <- 2L
  }
  nn_holdout_seed_use <- suppressWarnings(as.integer(params$nn_holdout_seed))
  if (!is.finite(nn_holdout_seed_use)) {
    nn_holdout_seed_use <- 271828L
  }
  run_nn_simulation_use <- isTRUE(params$run_nn_simulation)
  nn_simulation_n_use <- max(1L, suppressWarnings(as.integer(params$nn_simulation_n)))
  if (!is.finite(nn_simulation_n_use)) {
    nn_simulation_n_use <- 50L
  }
  nn_simulation_seed_use <- suppressWarnings(as.integer(params$nn_simulation_seed))
  if (!is.finite(nn_simulation_seed_use)) {
    nn_simulation_seed_use <- 161803L
  }
  nn_simulation_scenarios_use <- unique(as.character(unlist(params$nn_simulation_scenarios)))
  nn_simulation_scenarios_use <- nn_simulation_scenarios_use[nzchar(nn_simulation_scenarios_use)]
  if (!length(nn_simulation_scenarios_use)) {
    nn_simulation_scenarios_use <- c("sparse_zero_heavy", "moderate_observed", "two_step_supported")
  }
  run_nn_grf_simulation_use <- isTRUE(params$run_nn_grf_simulation)
  nn_grf_simulation_n_use <- max(1L, suppressWarnings(as.integer(benchmark_param_scalar(params, "nn_grf_simulation_n"))))
  if (!is.finite(nn_grf_simulation_n_use)) {
    nn_grf_simulation_n_use <- 12L
  }
  nn_grf_seed_use <- suppressWarnings(as.integer(benchmark_param_scalar(params, "nn_grf_seed")))
  if (!is.finite(nn_grf_seed_use)) {
    nn_grf_seed_use <- 424242L
  }
  nn_grf_lambdas_use <- sort(unique(suppressWarnings(as.numeric(unlist(params$nn_grf_lambdas)))))
  nn_grf_lambdas_use <- nn_grf_lambdas_use[is.finite(nn_grf_lambdas_use) & nn_grf_lambdas_use > 0]
  if (!length(nn_grf_lambdas_use)) {
    nn_grf_lambdas_use <- c(0.2, 0.4, 0.8, 1.6)
  }
  nn_grf_training_windows_use <- sort(unique(suppressWarnings(as.integer(unlist(params$nn_grf_training_windows)))))
  nn_grf_training_windows_use <- nn_grf_training_windows_use[is.finite(nn_grf_training_windows_use) & nn_grf_training_windows_use >= 2L]
  if (!length(nn_grf_training_windows_use)) {
    nn_grf_training_windows_use <- c(2L, 4L, 8L)
  }
  nn_grf_nboot_use <- suppressWarnings(as.integer(benchmark_param_scalar(params, "nn_grf_nboot")))
  if (!is.finite(nn_grf_nboot_use) || nn_grf_nboot_use < 1L) {
    nn_grf_nboot_use <- 15L
  }
  nn_grf_k_dim_use <- suppressWarnings(as.integer(benchmark_param_scalar(params, "nn_grf_k_dim")))
  if (!is.finite(nn_grf_k_dim_use) || nn_grf_k_dim_use < 2L) {
    nn_grf_k_dim_use <- 22L
  }
  nn_grf_n_centroids_use <- suppressWarnings(as.integer(benchmark_param_scalar(params, "nn_grf_n_centroids")))
  if (!is.finite(nn_grf_n_centroids_use) || nn_grf_n_centroids_use < 1L) {
    nn_grf_n_centroids_use <- 64L
  }
  nn_grf_time_max_use <- suppressWarnings(as.numeric(benchmark_param_scalar(params, "nn_grf_time_max")))
  if (!is.finite(nn_grf_time_max_use) || nn_grf_time_max_use <= 0) {
    nn_grf_time_max_use <- 140
  }
  nn_grf_passage_interval_use <- suppressWarnings(as.numeric(benchmark_param_scalar(params, "nn_grf_passage_interval")))
  if (!is.finite(nn_grf_passage_interval_use) || nn_grf_passage_interval_use <= 0) {
    nn_grf_passage_interval_use <- 20
  }
  nn_grf_sample_depth_use <- suppressWarnings(as.integer(benchmark_param_scalar(params, "nn_grf_sample_depth")))
  if (!is.finite(nn_grf_sample_depth_use) || nn_grf_sample_depth_use < 1L) {
    nn_grf_sample_depth_use <- 2000L
  }
  nn_grf_abm_pop_size_use <- suppressWarnings(as.numeric(benchmark_param_scalar(params, "nn_grf_abm_pop_size")))
  if (!is.finite(nn_grf_abm_pop_size_use) || nn_grf_abm_pop_size_use < 1) {
    nn_grf_abm_pop_size_use <- 50000
  }
  nn_grf_abm_delta_t_use <- suppressWarnings(as.numeric(benchmark_param_scalar(params, "nn_grf_abm_delta_t")))
  if (!is.finite(nn_grf_abm_delta_t_use) || nn_grf_abm_delta_t_use <= 0) {
    nn_grf_abm_delta_t_use <- 1
  }
  nn_grf_abm_max_pop_use <- suppressWarnings(as.numeric(benchmark_param_scalar(params, "nn_grf_abm_max_pop")))
  if (!is.finite(nn_grf_abm_max_pop_use)) {
    nn_grf_abm_max_pop_use <- 2e6
  }
  nn_grf_abm_culling_survival_use <- suppressWarnings(as.numeric(benchmark_param_scalar(params, "nn_grf_abm_culling_survival")))
  if (!is.finite(nn_grf_abm_culling_survival_use) ||
      nn_grf_abm_culling_survival_use < 0 ||
      nn_grf_abm_culling_survival_use > 1) {
    nn_grf_abm_culling_survival_use <- 0.01
  }

  focus_pids_use <- normalize_focus_pids(params)
  focus_minobs_use <- suppressWarnings(as.integer(params$focus_minobs))
  if (!is.finite(focus_minobs_use) || !focus_minobs_use %in% minobs_values_use) {
    focus_minobs_use <- max(minobs_values_use)
  }
  focus_pm_use <- suppressWarnings(as.numeric(params$focus_pm))
  if (!is.finite(focus_pm_use)) {
    focus_pm_use <- pm_values_use[1]
  } else {
    focus_pm_use <- pm_values_use[which.min(abs(pm_values_use - focus_pm_use))][1]
  }

  diploid_state <- paste(rep(2, 22), collapse = ".")
  stage_levels <- c("Primary", "Recurrent")
  beneficial_move_levels <- as.vector(rbind(paste0(seq_len(22), "+"), paste0(seq_len(22), "-")))
  parameter_levels_use <- parameter_labels_use

  set.seed(benchmark_seed_use)
  options(mc.cores = n_cores_use)
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 12))

  list(
    repo_dir = repo_dir,
    benchmark_dir = benchmark_dir,
    data_dir = data_dir,
    meta_path = meta_path,
    results_dir = results_dir,
    input_dir = input_dir,
    fit_dir = fit_dir,
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    cache_dir = cache_dir,
    patient_subset_use = patient_subset_use,
    minobs_values_use = minobs_values_use,
    pm_values_use = pm_values_use,
    parameter_labels_use = parameter_labels_use,
    parameter_levels_use = parameter_levels_use,
    selected_grid_n_use = selected_grid_n_use,
    n_cores_use = n_cores_use,
    nboot_use = nboot_use,
    n0_use = n0_use,
    nb_use = nb_use,
    benchmark_seed_use = benchmark_seed_use,
    correct_efflux_use = correct_efflux_use,
    nn_prior_fit_subset_use = nn_prior_fit_subset_use,
    nn_prior_zero_exposure_quantile_use = nn_prior_zero_exposure_quantile_use,
    nn_prior_zero_weight_scale_use = nn_prior_zero_weight_scale_use,
    nn_prior_zero_weight_cap_ratio_use = nn_prior_zero_weight_cap_ratio_use,
    nn_prior_zero_birth_fallback_weight_use = nn_prior_zero_birth_fallback_weight_use,
    nn_prior_zero_birth_child_floor_use = nn_prior_zero_birth_child_floor_use,
    nn_prior_zero_birth_child_shape_use = nn_prior_zero_birth_child_shape_use,
    nn_prior_zero_birth_replicate_floor_use = nn_prior_zero_birth_replicate_floor_use,
    nn_prior_zero_birth_replicate_shape_use = nn_prior_zero_birth_replicate_shape_use,
    nn_prior_two_step_support_use = nn_prior_two_step_support_use,
    nn_prior_two_step_support_min_use = nn_prior_two_step_support_min_use,
    nn_prior_two_step_cap_floor_use = nn_prior_two_step_cap_floor_use,
    cohort_transition_version_use = cohort_transition_version_use,
    cohort_contextual_apply_to_use = cohort_contextual_apply_to_use,
    cohort_context_keep_baseline_when_sparse_use = cohort_context_keep_baseline_when_sparse_use,
    cohort_context_lambda_sparse_unknown_use = cohort_context_lambda_sparse_unknown_use,
    force_refit_use = force_refit_use,
    rebuild_inputs_use = rebuild_inputs_use,
    run_benchmark_use = run_benchmark_use,
    include_posterior_use = include_posterior_use,
    render_figures_use = render_figures_use,
    top_shift_n_use = top_shift_n_use,
    run_nn_identifiability_use = run_nn_identifiability_use,
    run_nn_stability_use = run_nn_stability_use,
    run_nn_holdout_use = run_nn_holdout_use,
    nn_holdout_repeats_use = nn_holdout_repeats_use,
    nn_holdout_fraction_use = nn_holdout_fraction_use,
    nn_holdout_min_count_use = nn_holdout_min_count_use,
    nn_holdout_seed_use = nn_holdout_seed_use,
    run_nn_simulation_use = run_nn_simulation_use,
    nn_simulation_n_use = nn_simulation_n_use,
    nn_simulation_seed_use = nn_simulation_seed_use,
    nn_simulation_scenarios_use = nn_simulation_scenarios_use,
    run_nn_grf_simulation_use = run_nn_grf_simulation_use,
    nn_grf_simulation_n_use = nn_grf_simulation_n_use,
    nn_grf_seed_use = nn_grf_seed_use,
    nn_grf_lambdas_use = nn_grf_lambdas_use,
    nn_grf_training_windows_use = nn_grf_training_windows_use,
    nn_grf_nboot_use = nn_grf_nboot_use,
    nn_grf_k_dim_use = nn_grf_k_dim_use,
    nn_grf_n_centroids_use = nn_grf_n_centroids_use,
    nn_grf_time_max_use = nn_grf_time_max_use,
    nn_grf_passage_interval_use = nn_grf_passage_interval_use,
    nn_grf_sample_depth_use = nn_grf_sample_depth_use,
    nn_grf_abm_pop_size_use = nn_grf_abm_pop_size_use,
    nn_grf_abm_delta_t_use = nn_grf_abm_delta_t_use,
    nn_grf_abm_max_pop_use = nn_grf_abm_max_pop_use,
    nn_grf_abm_culling_survival_use = nn_grf_abm_culling_survival_use,
    focus_pids_use = focus_pids_use,
    focus_minobs_use = focus_minobs_use,
    focus_pm_use = focus_pm_use,
    diploid_state = diploid_state,
    stage_levels = stage_levels,
    beneficial_move_levels = beneficial_move_levels
  )
}

save_parameter_figures <- function(ctx,
                                   parameter_results_all_tbl,
                                   parameter_fit_summary_tbl,
                                   parameter_global_landscape_overview_tbl,
                                   parameter_global_landscape_state_tbl) {
  figure_paths <- build_parameter_figure_paths(ctx)

  if (ctx$render_figures_use && nrow(parameter_fit_summary_tbl)) {
    p_parameter_runtime <- ggplot2::ggplot(
      parameter_results_all_tbl %>% dplyr::filter(status == "ok") %>% dplyr::mutate(parameter_label = factor(parameter_label, levels = ctx$parameter_levels_use)),
      ggplot2::aes(x = parameter_label, y = elapsed_sec)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "#E8D8E8", color = "#633A63") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1.6, color = "#633A63") +
      ggplot2::labs(
        title = paste0("Runtime by parameter label (grid = ", ctx$selected_grid_n_use, ")"),
        x = "Parameter label",
        y = "Elapsed seconds"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
    ggplot2::ggsave(figure_paths$runtime, p_parameter_runtime, width = 9, height = 5, dpi = 150)

    p_parameter_xval <- ggplot2::ggplot(
      parameter_results_all_tbl %>% dplyr::filter(status == "ok") %>% dplyr::mutate(parameter_label = factor(parameter_label, levels = ctx$parameter_levels_use)),
      ggplot2::aes(x = parameter_label, y = xval_r2)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "#D8EBCF", color = "#31572C") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1.6, color = "#31572C") +
      ggplot2::labs(
        title = paste0("xval by parameter label (grid = ", ctx$selected_grid_n_use, ")"),
        x = "Parameter label",
        y = "xval_r2"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
    ggplot2::ggsave(figure_paths$xval, p_parameter_xval, width = 9, height = 5, dpi = 150)
  }

  if (ctx$render_figures_use && nrow(parameter_global_landscape_overview_tbl)) {
    p_parameter_global_diff <- parameter_global_landscape_overview_tbl %>%
      dplyr::mutate(comparison = paste(lhs_label, rhs_label, sep = " vs ")) %>%
      ggplot2::ggplot(ggplot2::aes(x = comparison, y = mean_abs_diff, fill = comparison)) +
      ggplot2::geom_col(width = 0.7, alpha = 0.9, show.legend = FALSE) +
      ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
      ggplot2::labs(
        title = "Benchmark-wide mean absolute differences",
        x = NULL,
        y = "Mean absolute difference"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
    ggplot2::ggsave(figure_paths$global_diff, p_parameter_global_diff, width = 9, height = 10, dpi = 150)
  }

  if (ctx$render_figures_use && nrow(parameter_global_landscape_state_tbl)) {
    p_parameter_state_diff <- parameter_global_landscape_state_tbl %>%
      dplyr::mutate(comparison = paste(lhs_label, rhs_label, sep = " vs ")) %>%
      ggplot2::ggplot(ggplot2::aes(x = state_group, y = mean_abs_diff, fill = comparison)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
      ggplot2::labs(
        title = "Benchmark-wide landscape mean differences by state group",
        x = "State group",
        y = "Mean absolute difference",
        fill = "Comparison"
      ) +
      ggplot2::theme_bw(base_size = 12)
    ggplot2::ggsave(figure_paths$global_state_diff, p_parameter_state_diff, width = 9, height = 5, dpi = 150)
  }

  figure_paths
}

make_all_sample_umap_plot <- function(plot_tbl,
                                      value_col,
                                      title,
                                      subtitle = NULL,
                                      scale_mode = c("raw", "reference"),
                                      raw_limits = NULL,
                                      reference_limits = NULL) {
  scale_mode <- match.arg(scale_mode)
  value_vec <- suppressWarnings(as.numeric(plot_tbl[[value_col]]))
  plot_tbl <- plot_tbl[order(value_vec, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

  if (scale_mode == "reference") {
    color_center <- 1
    color_limits <- if (!is.null(reference_limits) && length(reference_limits) == 2L && all(is.finite(reference_limits))) {
      reference_limits
    } else {
      c(min(value_vec, na.rm = TRUE), max(value_vec, na.rm = TRUE))
    }
    color_limits <- ensure_center_in_limits(color_limits, center = color_center)
    if (!all(is.finite(color_limits)) || diff(color_limits) == 0) {
      color_limits <- c(color_center - 5e-07, color_center + 5e-07)
    }
    color_knots <- c(
      seq(color_limits[1], color_center, length.out = 6),
      seq(color_center, color_limits[2], length.out = 6)[-1]
    )
    color_breaks <- c(color_limits[1], color_center, color_limits[2])
    legend_name <- "Relative fitness\n(vs global diploid)"
  } else {
    color_limits <- if (!is.null(raw_limits) && length(raw_limits) == 2L && all(is.finite(raw_limits))) {
      raw_limits
    } else {
      c(min(value_vec, na.rm = TRUE), max(value_vec, na.rm = TRUE))
    }
    if (!all(is.finite(color_limits)) || diff(color_limits) == 0) {
      color_limits[2] <- color_limits[1] + 1e-06
    }
    color_knots <- seq(color_limits[1], color_limits[2], length.out = 11)
    color_breaks <- c(color_limits[1], mean(color_limits), color_limits[2])
    legend_name <- "Raw fitness"
  }

  color_values <- scales::rescale(color_knots, from = color_limits)
  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = UMAP1, y = UMAP2, color = .data[[value_col]])) +
    ggplot2::geom_point(alpha = 0.85, size = 0.9) +
    ggplot2::scale_color_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
      values = color_values,
      limits = color_limits,
      breaks = color_breaks,
      labels = signif(color_breaks, 3),
      oob = scales::squish,
      name = legend_name
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

compute_validation_scatter_limits <- function(plot_tbl) {
  rng <- range(c(plot_tbl$estimate, plot_tbl$validation), finite = TRUE)
  if (!all(is.finite(rng))) {
    return(NULL)
  }

  if (rng[1] == rng[2]) {
    pad_base <- if (rng[1] == 0) 1 else abs(rng[1]) * 0.1
    lims_raw <- c(rng[1] - pad_base, rng[2] + pad_base)
  } else {
    lims_raw <- rng
  }

  span <- diff(lims_raw)
  pad <- span * 0.05
  c(lims_raw[1] - pad, lims_raw[2] + pad)
}

fit_validation_scatter_trend <- function(plot_tbl, weight_threshold = 0.5) {
  df_plot <- plot_tbl %>%
    dplyr::mutate(outlier_flag = FALSE)
  slope_est <- NA_real_
  intercept_est <- NA_real_
  fit_ribbon_tbl <- NULL

  if (nrow(df_plot) >= 2L && length(unique(df_plot$estimate)) >= 2L) {
    if (requireNamespace("MASS", quietly = TRUE)) {
      hub_fit <- tryCatch(
        MASS::rlm(validation ~ estimate, data = df_plot, psi = MASS::psi.huber),
        error = function(e) NULL
      )
      if (!is.null(hub_fit) && length(hub_fit$coefficients) >= 2L) {
        intercept_est <- as.numeric(hub_fit$coefficients[[1]])
        slope_est <- as.numeric(hub_fit$coefficients[[2]])
        if (length(hub_fit$w) == nrow(df_plot)) {
          df_plot$outlier_flag <- hub_fit$w < weight_threshold
        }
        x_seq <- seq(min(df_plot$estimate, na.rm = TRUE), max(df_plot$estimate, na.rm = TRUE), length.out = 100)
        X <- cbind(1, x_seq)
        hub_summary <- tryCatch(MASS::summary.rlm(hub_fit), error = function(e) NULL)
        cov_unscaled <- if (!is.null(hub_summary) && !is.null(hub_summary$cov.unscaled)) {
          hub_summary$cov.unscaled
        } else {
          NULL
        }
        scale_val <- if (!is.null(hub_summary) && !is.null(hub_summary$sigma)) {
          hub_summary$sigma
        } else {
          hub_fit$s
        }
        if (!is.null(cov_unscaled) && is.finite(scale_val)) {
          cov_beta <- cov_unscaled * (scale_val^2)
          se_fit <- sqrt(rowSums((X %*% cov_beta) * X))
          y_hat <- as.numeric(X %*% hub_fit$coefficients[1:2])
          fit_ribbon_tbl <- tibble::tibble(
            x = x_seq,
            y = y_hat,
            ylower = y_hat - 1.96 * se_fit,
            yupper = y_hat + 1.96 * se_fit
          )
        }
      }
    }

    if (is.null(fit_ribbon_tbl)) {
      lm_fit <- tryCatch(stats::lm(validation ~ estimate, data = df_plot), error = function(e) NULL)
      if (!is.null(lm_fit) && length(stats::coef(lm_fit)) >= 2L) {
        intercept_est <- as.numeric(stats::coef(lm_fit)[[1]])
        slope_est <- as.numeric(stats::coef(lm_fit)[[2]])
        x_seq <- seq(min(df_plot$estimate, na.rm = TRUE), max(df_plot$estimate, na.rm = TRUE), length.out = 100)
        pred <- tryCatch(
          stats::predict(lm_fit, newdata = data.frame(estimate = x_seq), interval = "confidence"),
          error = function(e) NULL
        )
        if (!is.null(pred)) {
          fit_ribbon_tbl <- tibble::tibble(
            x = x_seq,
            y = pred[, "fit"],
            ylower = pred[, "lwr"],
            yupper = pred[, "upr"]
          )
        }
      }
    }
  }

  list(
    df_plot = df_plot,
    fit_ribbon_tbl = fit_ribbon_tbl,
    slope_est = slope_est,
    intercept_est = intercept_est
  )
}

make_validation_scatter_panel <- function(plot_tbl, title_text) {
  lims <- compute_validation_scatter_limits(plot_tbl)
  fit_res <- fit_validation_scatter_trend(plot_tbl)
  df_plot <- fit_res$df_plot
  fit_ribbon_tbl <- fit_res$fit_ribbon_tbl

  annotation_lines <- c(
    if (is.finite(fit_res$slope_est) && is.finite(fit_res$intercept_est)) {
      sprintf("Fit slope = %.3f\nintercept = %.3f", fit_res$slope_est, fit_res$intercept_est)
    } else {
      "Fit slope = NA\nintercept = NA"
    },
    paste0("Total points: ", nrow(df_plot))
  )
  if (any(df_plot$outlier_flag %in% TRUE, na.rm = TRUE)) {
    annotation_lines <- c(annotation_lines, paste0("Outliers: ", sum(df_plot$outlier_flag %in% TRUE, na.rm = TRUE)))
  }

  has_outliers <- any(df_plot$outlier_flag %in% TRUE, na.rm = TRUE)
  annotation_x <- if (!is.null(lims)) lims[1] + diff(lims) * 0.02 else min(df_plot$estimate, na.rm = TRUE)
  annotation_y <- if (!is.null(lims)) lims[2] - diff(lims) * 0.02 else max(df_plot$validation, na.rm = TRUE)

  plot_obj <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = estimate, y = validation, color = outlier_flag, shape = outlier_flag)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::annotate(
      "text",
      x = annotation_x,
      y = annotation_y,
      hjust = 0,
      vjust = 1,
      label = paste(annotation_lines, collapse = "\n"),
      size = 4.4,
      color = "black",
      lineheight = 1.1
    ) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "black", `TRUE` = "red"),
      labels = c(`FALSE` = "Inlier", `TRUE` = "Outlier"),
      guide = if (has_outliers) "legend" else "none"
    ) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17),
      labels = c(`FALSE` = "Inlier", `TRUE` = "Outlier"),
      guide = if (has_outliers) "legend" else "none"
    ) +
    ggplot2::labs(title = title_text, x = "Estimate", y = "Validation") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))

  if (!is.null(fit_ribbon_tbl)) {
    plot_obj <- plot_obj +
      ggplot2::geom_ribbon(
        data = fit_ribbon_tbl,
        ggplot2::aes(x = x, ymin = ylower, ymax = yupper),
        inherit.aes = FALSE,
        fill = "#b2d7ff",
        alpha = 0.35
      ) +
      ggplot2::geom_line(
        data = fit_ribbon_tbl,
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        color = "#2c7fb8",
        linewidth = 0.7
      )
  }

  if (!is.null(lims)) {
    plot_obj <- plot_obj +
      ggplot2::scale_x_continuous(limits = lims) +
      ggplot2::scale_y_continuous(limits = lims)
  }

  plot_obj <- plot_obj + ggplot2::coord_equal()

  if (has_outliers) {
    plot_obj <- plot_obj +
      ggplot2::theme(
        legend.position = c(0.97, 0.03),
        legend.justification = c(1, 0),
        legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.7), color = NA),
        legend.text = ggplot2::element_text(size = 9),
        legend.key.size = grid::unit(0.35, "cm"),
        legend.title = ggplot2::element_blank()
      )
  } else {
    plot_obj <- plot_obj + ggplot2::theme(legend.position = "none")
  }

  plot_obj
}

save_retained_summary_figures <- function(ctx,
                                          retained_fit_tbl,
                                          retained_landscape_long_tbl,
                                          retained_xval_scatter_tbl,
                                          retained_umap_res) {
  figure_paths <- build_parameter_figure_paths(ctx)
  landscape_artifact_tbl <- tibble::tibble(
    parameter_label = ctx$parameter_levels_use,
    title = ctx$parameter_levels_use,
    png_path = file.path(ctx$figures_dir, paste0("all_sample_landscape_distribution_", ctx$parameter_levels_use, ".png")),
    n_patients = NA_integer_,
    n_karyotypes = NA_integer_
  )
  xval_artifact_tbl <- tibble::tibble(
    parameter_label = ctx$parameter_levels_use,
    title = ctx$parameter_levels_use,
    png_path = file.path(ctx$figures_dir, paste0("all_sample_xval_", ctx$parameter_levels_use, ".png")),
    pooled_correlation = NA_real_,
    pooled_r2 = NA_real_,
    n_points = NA_integer_,
    n_patients = NA_integer_
  )

  umap_artifact_tbl <- dplyr::bind_rows(
    tibble::tibble(
      panel_label = paste0(ctx$parameter_levels_use, "__raw"),
      parameter_label = ctx$parameter_levels_use,
      panel_group = "raw",
      title = paste0(ctx$parameter_levels_use, " raw_fitness"),
      png_path = file.path(ctx$figures_dir, paste0("all_sample_umap_", ctx$parameter_levels_use, ".png"))
    ),
    tibble::tibble(
      panel_label = paste0(ctx$parameter_levels_use, "__reference_global_diploid"),
      parameter_label = ctx$parameter_levels_use,
      panel_group = "reference_global_diploid",
      title = paste0(ctx$parameter_levels_use, " reference_global_diploid"),
      png_path = file.path(ctx$figures_dir, paste0("all_sample_umap_reference_global_diploid_", ctx$parameter_levels_use, ".png"))
    )
  )

  if (!ctx$render_figures_use) {
    return(list(
      figure_paths = figure_paths,
      landscape_artifact_tbl = landscape_artifact_tbl,
      xval_artifact_tbl = xval_artifact_tbl,
      umap_artifact_tbl = umap_artifact_tbl
    ))
  }

  if (nrow(retained_landscape_long_tbl)) {
    raw_limits <- range(retained_landscape_long_tbl$mean, na.rm = TRUE)
    sample_stat_tbl <- summarize_selected_landscape_sample_tbl(retained_landscape_long_tbl)

    for (parameter_label_name in ctx$parameter_levels_use) {
      raw_hist_tbl <- retained_landscape_long_tbl %>%
        dplyr::mutate(parameter_label = as.character(parameter_label)) %>%
        dplyr::filter(parameter_label == parameter_label_name)
      if (!nrow(raw_hist_tbl)) {
        next
      }
      sample_stat_sub_tbl <- sample_stat_tbl %>%
        dplyr::mutate(parameter_label = as.character(parameter_label)) %>%
        dplyr::filter(parameter_label == parameter_label_name)

      p_hist <- ggplot2::ggplot(raw_hist_tbl, ggplot2::aes(x = mean)) +
        ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)), bins = 60, fill = "#9ECAE1", color = "#2C7FB8", alpha = 0.45) +
        ggplot2::geom_density(linewidth = 0.8, color = "#08519C", adjust = 1) +
        ggplot2::coord_cartesian(xlim = raw_limits) +
        ggplot2::labs(
          title = parameter_label_name,
          subtitle = "Pooled landscape mean values across that method's positive-xval-selected samples",
          x = "Landscape mean fitness",
          y = "Density"
        ) +
        ggplot2::theme_bw(base_size = 11)

      p_mean <- ggplot2::ggplot(sample_stat_sub_tbl, ggplot2::aes(x = 1, y = landscape_mean_mean)) +
        ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, fill = "#9ECAE1", alpha = 0.8) +
        ggplot2::geom_jitter(width = 0.08, height = 0, alpha = 0.7, size = 1.8, color = "#08519C") +
        ggplot2::labs(
          title = "Per-sample mean",
          x = NULL,
          y = "Mean of landscape mean"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

      p_sd <- ggplot2::ggplot(sample_stat_sub_tbl, ggplot2::aes(x = 1, y = landscape_mean_sd)) +
        ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, fill = "#FDD0A2", alpha = 0.8) +
        ggplot2::geom_jitter(width = 0.08, height = 0, alpha = 0.7, size = 1.8, color = "#A63603") +
        ggplot2::labs(
          title = "Per-sample SD",
          x = NULL,
          y = "SD of landscape mean"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

      p_landscape_distribution <- patchwork::wrap_plots(
        p_hist,
        patchwork::wrap_plots(p_mean, p_sd, nrow = 1),
        ncol = 1,
        heights = c(2, 1)
      )

      png_path <- landscape_artifact_tbl %>%
        dplyr::filter(parameter_label == parameter_label_name) %>%
        dplyr::pull(png_path)
      ggplot2::ggsave(png_path[[1]], p_landscape_distribution, width = 8.5, height = 8.5, dpi = 150)

      landscape_artifact_tbl <- landscape_artifact_tbl %>%
        dplyr::mutate(
          n_patients = dplyr::if_else(parameter_label == parameter_label_name, as.integer(length(unique(as.character(raw_hist_tbl$patient_id)))), n_patients),
          n_karyotypes = dplyr::if_else(parameter_label == parameter_label_name, as.integer(nrow(raw_hist_tbl)), n_karyotypes)
        )
    }
  }

  if (nrow(retained_xval_scatter_tbl)) {
    for (parameter_label_name in ctx$parameter_levels_use) {
      plot_tbl <- retained_xval_scatter_tbl %>%
        dplyr::mutate(
          parameter_label = as.character(parameter_label),
          patient_id = factor(patient_id, levels = sort_pid_levels(patient_id))
        ) %>%
        dplyr::filter(parameter_label == parameter_label_name)
      if (!nrow(plot_tbl)) {
        next
      }

      pooled_correlation_value <- if (nrow(plot_tbl) >= 2L) {
        suppressWarnings(stats::cor(plot_tbl$estimate, plot_tbl$validation, use = "complete.obs"))
      } else {
        NA_real_
      }
      pooled_r2_value <- if (nrow(plot_tbl) >= 2L) {
        suppressWarnings(alfakR:::R2R(plot_tbl$estimate, plot_tbl$validation))
      } else {
        NA_real_
      }
      patient_ids <- sort_pid_levels(unique(as.character(plot_tbl$patient_id)))
      patient_plot_list <- lapply(patient_ids, function(patient_id_name) {
        patient_plot_tbl <- plot_tbl %>%
          dplyr::mutate(patient_id = as.character(patient_id)) %>%
          dplyr::filter(patient_id == patient_id_name)
        patient_title <- sprintf(
          "%s (MINOBS=%s, Signed R²=%.3f)",
          patient_id_name,
          unique(patient_plot_tbl$minobs)[1],
          unique(patient_plot_tbl$xval_r2)[1]
        )
        make_validation_scatter_panel(patient_plot_tbl, patient_title)
      })
      subtitle_text <- sprintf(
        "fq-only cross-validation points; pooled corr = %.3f, pooled R² = %.3f, n = %d",
        pooled_correlation_value,
        pooled_r2_value,
        nrow(plot_tbl)
      )
      p_xval_scatter <- patchwork::wrap_plots(patient_plot_list, ncol = 2) +
        patchwork::plot_annotation(
          title = parameter_label_name,
          subtitle = subtitle_text,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14),
            plot.subtitle = ggplot2::element_text(size = 11)
          )
        )

      png_path <- xval_artifact_tbl %>%
        dplyr::filter(parameter_label == parameter_label_name) %>%
        dplyr::pull(png_path)
      ggplot2::ggsave(
        png_path[[1]],
        p_xval_scatter,
        width = 12.5,
        height = max(6.5, 4.3 * ceiling(length(patient_plot_list) / 2)),
        dpi = 150
      )

      xval_artifact_tbl <- xval_artifact_tbl %>%
        dplyr::mutate(
          pooled_correlation = dplyr::if_else(parameter_label == parameter_label_name, pooled_correlation_value, .data$pooled_correlation),
          pooled_r2 = dplyr::if_else(parameter_label == parameter_label_name, pooled_r2_value, .data$pooled_r2),
          n_points = dplyr::if_else(parameter_label == parameter_label_name, as.integer(nrow(plot_tbl)), n_points),
          n_patients = dplyr::if_else(parameter_label == parameter_label_name, as.integer(length(unique(plot_tbl$patient_id))), n_patients)
        )
    }
  }

  if (!is.null(retained_umap_res$landscape_with_umap) && nrow(retained_umap_res$landscape_with_umap)) {
    raw_limits <- range(retained_umap_res$landscape_with_umap$mean, na.rm = TRUE)
    reference_limits <- NULL
    if (!is.null(retained_umap_res$reference_umap_tbl) && nrow(retained_umap_res$reference_umap_tbl)) {
      reference_limits <- range(retained_umap_res$reference_umap_tbl$reference_value, na.rm = TRUE)
    }
    for (i in seq_along(ctx$parameter_levels_use)) {
      parameter_label_name <- ctx$parameter_levels_use[[i]]
      plot_tbl <- retained_umap_res$landscape_with_umap %>%
        dplyr::filter(parameter_label == parameter_label_name)
      if (!nrow(plot_tbl)) {
        next
      }
      raw_png_path <- umap_artifact_tbl %>%
        dplyr::filter(parameter_label == parameter_label_name, panel_group == "raw") %>%
        dplyr::pull(png_path)
      plot_obj <- make_all_sample_umap_plot(
        plot_tbl = plot_tbl,
        value_col = "mean",
        title = paste(parameter_label_name, "raw_fitness"),
        subtitle = "All per-method positive-xval-selected samples, raw fitness",
        scale_mode = "raw",
        raw_limits = raw_limits
      )
      ggplot2::ggsave(raw_png_path[[1]], plot_obj, width = 5.5, height = 5, dpi = 150)
    }

    if (!is.null(retained_umap_res$reference_umap_tbl) && nrow(retained_umap_res$reference_umap_tbl)) {
      for (parameter_label_name in ctx$parameter_levels_use) {
        ref_plot_tbl <- retained_umap_res$reference_umap_tbl %>%
          dplyr::filter(parameter_label == parameter_label_name)
        if (!nrow(ref_plot_tbl)) {
          next
        }
        ref_png_path <- umap_artifact_tbl %>%
          dplyr::filter(parameter_label == parameter_label_name, panel_group == "reference_global_diploid") %>%
          dplyr::pull(png_path)
        ref_plot <- make_all_sample_umap_plot(
          plot_tbl = ref_plot_tbl,
          value_col = "reference_value",
          title = paste(parameter_label_name, "reference_global_diploid"),
          subtitle = "All per-method positive-xval-selected samples, scaled by that method's global diploid mean",
          scale_mode = "reference",
          reference_limits = reference_limits
        )
        ggplot2::ggsave(ref_png_path[[1]], ref_plot, width = 5.5, height = 5, dpi = 150)
      }
    }
  }

  list(
    figure_paths = figure_paths,
    landscape_artifact_tbl = landscape_artifact_tbl,
    xval_artifact_tbl = xval_artifact_tbl,
    umap_artifact_tbl = umap_artifact_tbl
  )
}

save_input_overview_figures <- function(ctx,
                                        input_fq_nn_summary_tbl,
                                        input_fq_group_nn_summary_tbl) {
  figure_paths <- build_input_figure_paths(ctx)

  if (!ctx$render_figures_use || is.null(input_fq_nn_summary_tbl) || !nrow(input_fq_nn_summary_tbl)) {
    return(figure_paths)
  }

  patient_levels <- rev(sort_pid_levels(unique(as.character(input_fq_nn_summary_tbl$patient_id))))
  plot_height <- max(5.5, 0.28 * length(patient_levels) + 3)
  minobs_levels <- sort(unique(as.integer(input_fq_nn_summary_tbl$minobs)))

  input_count_long_tbl <- input_fq_nn_summary_tbl %>%
    dplyr::select(patient_id, minobs, n_fq, n_nn, n_nn_observed) %>%
    tidyr::pivot_longer(cols = c(n_fq, n_nn, n_nn_observed), names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      patient_id = factor(patient_id, levels = patient_levels),
      minobs = factor(minobs, levels = minobs_levels),
      metric = factor(
        metric,
        levels = c("n_fq", "n_nn", "n_nn_observed"),
        labels = c("Number of fq", "Number of one-step neighbours", "Observed one-step neighbours")
      )
    )

  p_input_counts <- ggplot2::ggplot(input_count_long_tbl, ggplot2::aes(x = minobs, y = patient_id, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = value), size = 2.8) +
    ggplot2::facet_wrap(~ metric, ncol = 1) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "grey90") +
    ggplot2::labs(
      title = "Input fq and one-step-neighbour counts by sample and minobs",
      subtitle = "fq uses the benchmark minobs rule after diploid removal",
      x = "minobs",
      y = NULL,
      fill = "Count"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(figure_paths$fq_nn_counts, p_input_counts, width = 8.5, height = plot_height, dpi = 150)

  input_prop_long_tbl <- input_fq_nn_summary_tbl %>%
    dplyr::select(patient_id, minobs, prop_fq_count_up, prop_nn_count_up, prop_nn_observed) %>%
    tidyr::pivot_longer(cols = c(prop_fq_count_up, prop_nn_count_up, prop_nn_observed), names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      patient_id = factor(patient_id, levels = patient_levels),
      minobs = factor(minobs, levels = minobs_levels),
      metric = factor(
        metric,
        levels = c("prop_fq_count_up", "prop_nn_count_up", "prop_nn_observed"),
        labels = c("fq count-up proportion", "Neighbour count-up proportion", "Observed-neighbour proportion")
      ),
      value_label = ifelse(is.finite(value), sprintf("%.2f", value), "")
    )

  p_input_props <- ggplot2::ggplot(input_prop_long_tbl, ggplot2::aes(x = minobs, y = patient_id, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = value_label), size = 2.7) +
    ggplot2::facet_wrap(~ metric, ncol = 1) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#CB181D", limits = c(0, 1), na.value = "grey90") +
    ggplot2::labs(
      title = "Input count-up proportions by sample and minobs",
      subtitle = "A state counts as up when timepoint 2 count is greater than timepoint 1 count",
      x = "minobs",
      y = NULL,
      fill = "Proportion"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(figure_paths$fq_nn_up_prop, p_input_props, width = 8.5, height = plot_height, dpi = 150)

  if (!is.null(input_fq_group_nn_summary_tbl) && nrow(input_fq_group_nn_summary_tbl)) {
    input_group_heatmap_data <- build_input_group_prop_heatmap_data(
      input_fq_group_nn_summary_tbl = input_fq_group_nn_summary_tbl,
      patient_ids = unique(as.character(input_fq_group_nn_summary_tbl$patient_id)),
      minobs_levels = minobs_levels
    )
    panel_mat <- input_group_heatmap_data$panel_mat
    cluster_mat <- panel_mat
    if (any(!is.finite(cluster_mat))) {
      col_fill <- apply(cluster_mat, 2, function(x) {
        mm <- mean(x, na.rm = TRUE)
        if (is.finite(mm)) mm else 0
      })
      for (j in seq_len(ncol(cluster_mat))) {
        miss <- !is.finite(cluster_mat[, j])
        if (any(miss)) {
          cluster_mat[miss, j] <- col_fill[j]
        }
      }
    }
    row_cluster <- if (nrow(cluster_mat) > 1) {
      stats::as.dendrogram(stats::hclust(stats::dist(cluster_mat)))
    } else {
      FALSE
    }
    col_fun_group <- circlize::colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#9ECAE1", "#238B45"))

    ht_group <- ComplexHeatmap::Heatmap(
      panel_mat,
      name = "Neighbour\ncount-up\nproportion",
      col = col_fun_group,
      cluster_rows = row_cluster,
      cluster_columns = FALSE,
      column_split = input_group_heatmap_data$column_split,
      na_col = "grey90",
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8),
      rect_gp = grid::gpar(col = "white", lwd = 1),
      row_dend_side = "left",
      column_title = "Neighbour count-up proportion split by fq proportion direction",
      column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        val <- panel_mat[i, j]
        if (is.finite(val)) {
          grid::grid.text(sprintf("%.2f", val), x, y, gp = grid::gpar(fontsize = 7))
        }
      }
    )

    grDevices::png(figure_paths$fq_group_nn_up_prop, width = 9 * 150, height = plot_height * 150, res = 150)
    ComplexHeatmap::draw(ht_group, heatmap_legend_side = "right")
    grDevices::dev.off()
  }

  figure_paths
}

build_focus_outputs <- function(focus_pid, ctx, parameter_selected_fit_tbl, parameter_pair_results = NULL) {
  focus_fit_index_tbl <- parameter_selected_fit_tbl %>%
    dplyr::filter(
      patient_id == focus_pid,
      parameter_label %in% ctx$parameter_levels_use
    ) %>%
    dplyr::arrange(match(parameter_label, ctx$parameter_levels_use))
  save_table_bundle(focus_fit_index_tbl, focus_table_stem(ctx, focus_pid, "fit_index"))

  focus_parameter_bundles <- build_focus_parameter_bundles_from_selected(
    selected_fit_tbl = focus_fit_index_tbl,
    beneficial_move_levels = ctx$beneficial_move_levels,
    parameter_levels = ctx$parameter_levels_use
  )
  focus_parameter_levels_available <- names(focus_parameter_bundles)

  focus_landscape_long_tbl <- landscape_long_from_bundles(focus_parameter_bundles)
  focus_landscape_summary_tbl <- summarize_landscape_by_parameter(focus_landscape_long_tbl)
  focus_landscape_variation_tbl <- summarize_focus_landscape_variation(
    focus_landscape_long_tbl,
    parameter_levels = focus_parameter_levels_available,
    top_n = ctx$top_shift_n_use
  )
  focus_parity_tbl <- build_focus_parity_tbl(
    focus_parameter_bundles,
    parameter_levels = focus_parameter_levels_available,
    value_col = "mean"
  )

  focus_pairwise_component_tbl <- tibble::tibble()
  save_table_bundle(focus_pairwise_component_tbl, focus_table_stem(ctx, focus_pid, "pairwise_component_summary"))

  focus_pairwise_landscape_group_tbl <- tibble::tibble()
  save_table_bundle(focus_pairwise_landscape_group_tbl, focus_table_stem(ctx, focus_pid, "pairwise_landscape_group_summary"))

  focus_top_shift_tbl <- tibble::tibble()
  save_table_bundle(focus_top_shift_tbl, focus_table_stem(ctx, focus_pid, "pairwise_top_landscape_shifts"))

  focus_beneficial_profiles <- lapply(focus_parameter_levels_available, function(parameter_name) {
    focus_parameter_bundles[[parameter_name]]$beneficial
  })
  names(focus_beneficial_profiles) <- focus_parameter_levels_available

  focus_beneficial_long_tbl <- beneficial_long_from_profiles(
    beneficial_profiles = focus_beneficial_profiles,
    beneficial_move_levels = ctx$beneficial_move_levels,
    parameter_levels = focus_parameter_levels_available
  )
  focus_beneficial_summary_tbl <- summarize_beneficial_by_parameter(focus_beneficial_long_tbl)
  focus_beneficial_shift_tbl <- build_focus_beneficial_shift_tbl(
    focus_beneficial_long_tbl,
    parameter_levels = focus_parameter_levels_available,
    top_n = ctx$top_shift_n_use
  )
  focus_beneficial_proportion_matrix_tbl <- if (nrow(focus_beneficial_long_tbl)) {
    focus_beneficial_long_tbl %>%
      dplyr::mutate(move = as.character(move)) %>%
      dplyr::select(parameter_label, move, proportion) %>%
      tidyr::pivot_wider(names_from = move, values_from = proportion)
  } else {
    tibble::tibble()
  }
  focus_beneficial_valid_n_matrix_tbl <- if (nrow(focus_beneficial_long_tbl)) {
    focus_beneficial_long_tbl %>%
      dplyr::mutate(move = as.character(move)) %>%
      dplyr::select(parameter_label, move, valid_n) %>%
      tidyr::pivot_wider(names_from = move, values_from = valid_n)
  } else {
    tibble::tibble()
  }

  focus_input_rds <- file.path(ctx$input_dir, paste0(focus_pid, ".Rds"))
  focus_decomposition <- build_focus_observed_latent_decomposition(
    bundle_list = focus_parameter_bundles,
    input_rds = focus_input_rds,
    diploid_state = ctx$diploid_state,
    beneficial_move_levels = ctx$beneficial_move_levels,
    parameter_levels = focus_parameter_levels_available
  )
  focus_state_fitness_summary_tbl <- focus_decomposition$state_summary
  focus_edge_decomposition_summary_tbl <- focus_decomposition$edge_summary
  focus_nn_prior_diag_summary_tbl <- focus_decomposition$nn_prior_diag_summary

  focus_umap_parameter_levels <- intersect(
    ctx$parameter_levels_use,
    focus_parameter_levels_available
  )
  focus_umap_artifact_tbl <- dplyr::bind_rows(lapply(focus_umap_parameter_levels, function(parameter_label_name) {
    landscape_df <- focus_parameter_bundles[[parameter_label_name]]$landscape
    if (is.null(landscape_df) || !nrow(landscape_df)) {
      return(tibble::tibble(
        parameter_label = parameter_label_name,
        scale_mode = "relative",
        n_karyotypes = 0L,
        png_path = NA_character_
      ))
    }

    plot_obj <- make_focus_umap_plot(
      landscape_df = landscape_df,
      patient_id = focus_pid,
      parameter_label = parameter_label_name,
      benchmark_seed = ctx$benchmark_seed_use,
      diploid_state = ctx$diploid_state,
      scale_mode = "relative"
    )
    png_path <- file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_", parameter_label_name, "_relative_umap.png"))
    if (ctx$render_figures_use) {
      ggplot2::ggsave(png_path, plot_obj, width = 5.5, height = 5, dpi = 150)
    }

    tibble::tibble(
      parameter_label = parameter_label_name,
      scale_mode = "relative",
      n_karyotypes = nrow(landscape_df),
      png_path = png_path
    )
  }))

  save_table_bundle(focus_landscape_summary_tbl, focus_table_stem(ctx, focus_pid, "landscape_summary"))
  save_table_bundle(focus_landscape_variation_tbl, focus_table_stem(ctx, focus_pid, "top_variable_karyotypes"))
  save_table_bundle(focus_beneficial_summary_tbl, focus_table_stem(ctx, focus_pid, "beneficial_summary"))
  save_table_bundle(focus_beneficial_shift_tbl, focus_table_stem(ctx, focus_pid, "beneficial_top_shifts"))
  save_table_bundle(focus_beneficial_proportion_matrix_tbl, focus_table_stem(ctx, focus_pid, "beneficial_proportion_matrix"))
  save_table_bundle(focus_beneficial_valid_n_matrix_tbl, focus_table_stem(ctx, focus_pid, "beneficial_valid_n_matrix"))
  save_table_bundle(focus_state_fitness_summary_tbl, focus_table_stem(ctx, focus_pid, "state_fitness_summary"))
  save_table_bundle(focus_edge_decomposition_summary_tbl, focus_table_stem(ctx, focus_pid, "edge_decomposition_summary"))
  save_table_bundle(focus_nn_prior_diag_summary_tbl, focus_table_stem(ctx, focus_pid, "nn_prior_diagnostics_summary"))
  save_table_bundle(focus_umap_artifact_tbl, focus_table_stem(ctx, focus_pid, "umap_artifacts"))

  focus_figure_paths <- build_focus_figure_paths(ctx, focus_pid)
  focus_density_artifact_tbl <- tibble::tibble(
    parameter_label = focus_parameter_levels_available,
    png_path = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_landscape_density_", focus_parameter_levels_available, ".png"))
  )

  if (ctx$render_figures_use && nrow(focus_landscape_long_tbl)) {
    for (parameter_label_name in focus_parameter_levels_available) {
      focus_landscape_density_long_tbl <- focus_landscape_long_tbl %>%
        dplyr::mutate(parameter_label = as.character(parameter_label)) %>%
        dplyr::filter(parameter_label == parameter_label_name) %>%
        tidyr::pivot_longer(cols = c(mean, median), names_to = "metric", values_to = "value")
      if (!nrow(focus_landscape_density_long_tbl)) {
        next
      }

      p_focus_landscape_density <- ggplot2::ggplot(
        focus_landscape_density_long_tbl,
        ggplot2::aes(x = value)
      ) +
        ggplot2::geom_density(fill = "#9ECAE1", color = "#08519C", alpha = 0.35, adjust = 1) +
        ggplot2::facet_wrap(~ metric, scales = "free", ncol = 2) +
        ggplot2::labs(
          title = paste0(focus_pid, " ", parameter_label_name, " landscape distributions"),
          subtitle = "Retained best positive-xval fit for this method",
          x = "Landscape value",
          y = "Density"
        ) +
        ggplot2::theme_bw(base_size = 12)

      png_path <- focus_density_artifact_tbl %>%
        dplyr::filter(parameter_label == parameter_label_name) %>%
        dplyr::pull(png_path)
      ggplot2::ggsave(png_path[[1]], p_focus_landscape_density, width = 8, height = 4.5, dpi = 150)
    }
  }
  save_table_bundle(focus_density_artifact_tbl, focus_table_stem(ctx, focus_pid, "density_artifacts"))

  if (ctx$render_figures_use && nrow(focus_parity_tbl)) {
    p_focus_parity <- ggplot2::ggplot(focus_parity_tbl, ggplot2::aes(x = lhs_value, y = rhs_value, color = state_group)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
      ggplot2::geom_point(alpha = 0.28, size = 0.9) +
      ggplot2::facet_wrap(~ comparison, scales = "free", ncol = 2) +
      ggplot2::scale_color_manual(values = c(fq = "#C53030", nn = "#2B6CB0", other = "#6B7280")) +
      ggplot2::labs(
        title = paste0(focus_pid, " pairwise landscape parity plots"),
        subtitle = "Each facet shows matched landscape means on shared karyotypes",
        x = "Left-hand landscape mean",
        y = "Right-hand landscape mean",
        color = "State group"
      ) +
      ggplot2::theme_bw(base_size = 12)
    ggplot2::ggsave(focus_figure_paths$parity, p_focus_parity, width = 10, height = 8, dpi = 150)
  }

  if (ctx$render_figures_use && nrow(focus_beneficial_long_tbl)) {
    p_focus_beneficial <- focus_beneficial_long_tbl %>%
      dplyr::mutate(parameter_label = factor(parameter_label, levels = focus_parameter_levels_available)) %>%
      ggplot2::ggplot(ggplot2::aes(x = move, y = parameter_label, fill = proportion)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.15) +
      ggplot2::scale_fill_gradient2(
        low = "#2B6CB0",
        mid = "#F7F7F7",
        high = "#C53030",
        midpoint = 0.5,
        limits = c(0, 1),
        na.value = "grey85"
      ) +
      ggplot2::labs(
        title = paste0(focus_pid, " beneficial-karyotype proportion by parameter label"),
        subtitle = "Retained best positive-xval fit per parameter label",
        x = "Chromosome move",
        y = NULL,
        fill = "Beneficial\nproportion"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = ggplot2::element_blank()
      )
    ggplot2::ggsave(focus_figure_paths$beneficial_heatmap, p_focus_beneficial, width = 14, height = 3.5, dpi = 150)
  }

  if (ctx$render_figures_use && nrow(focus_decomposition$state_long)) {
    p_focus_state_fitness <- focus_decomposition$state_long %>%
      dplyr::mutate(
        parameter_label = factor(parameter_label, levels = focus_parameter_levels_available),
        state_class = factor(state_class, levels = focus_state_class_levels())
      ) %>%
      ggplot2::ggplot(ggplot2::aes(x = state_class, y = landscape_mean, fill = state_class)) +
      ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6, color = NA) +
      ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.9) +
      ggplot2::facet_wrap(~ parameter_label, scales = "free_y", ncol = 2) +
      ggplot2::scale_fill_manual(values = c(fq = "#C53030", observed_nn = "#2B6CB0", latent_nn = "#9CA3AF")) +
      ggplot2::labs(
        title = paste0(focus_pid, " fq and one-step-neighbour fitness decomposition"),
        subtitle = "Landscape mean split into fq, observed one-step neighbours, and latent one-step neighbours",
        x = NULL,
        y = "Landscape mean fitness",
        fill = "State class"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
    ggplot2::ggsave(focus_figure_paths$state_fitness, p_focus_state_fitness, width = 10, height = 6, dpi = 150)
  }

  if (ctx$render_figures_use && nrow(focus_decomposition$edge_long)) {
    p_focus_edge_delta <- focus_decomposition$edge_long %>%
      dplyr::mutate(
        parameter_label = factor(parameter_label, levels = focus_parameter_levels_available),
        child_state_class = factor(as.character(child_state_class), levels = focus_edge_class_levels()[-1])
      ) %>%
      ggplot2::ggplot(ggplot2::aes(x = child_state_class, y = delta, fill = child_state_class)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
      ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6, color = NA) +
      ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.35, alpha = 0.9) +
      ggplot2::facet_wrap(~ parameter_label, scales = "free_y", ncol = 2) +
      ggplot2::scale_fill_manual(values = c(fq_target = "#DD6B20", observed_nn = "#2B6CB0", latent_nn = "#6B7280")) +
      ggplot2::labs(
        title = paste0(focus_pid, " edge-level fitness shifts from fq parents"),
        subtitle = "delta = child fitness - fq parent fitness for all valid one-step moves present in the landscape",
        x = "Child state class",
        y = "Fitness delta",
        fill = "Child class"
      ) +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(focus_figure_paths$edge_delta, p_focus_edge_delta, width = 10, height = 6, dpi = 150)
  }

  if (ctx$render_figures_use && nrow(focus_edge_decomposition_summary_tbl)) {
    p_focus_edge_beneficial <- focus_edge_decomposition_summary_tbl %>%
      dplyr::mutate(
        parameter_label = factor(parameter_label, levels = focus_parameter_levels_available),
        edge_class = factor(edge_class, levels = focus_edge_class_levels()),
        edge_label = paste0(edge_class, "\n(n=", n_edges, ")")
      ) %>%
      ggplot2::ggplot(ggplot2::aes(x = edge_label, y = beneficial_proportion, fill = edge_class)) +
      ggplot2::geom_col(width = 0.72, alpha = 0.9) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(is.finite(beneficial_proportion), sprintf("%.2f", beneficial_proportion), "")),
        vjust = -0.4,
        size = 3
      ) +
      ggplot2::facet_wrap(~ parameter_label, ncol = 2) +
      ggplot2::scale_fill_manual(values = c(all_valid_moves = "#4A5568", fq_target = "#DD6B20", observed_nn = "#2B6CB0", latent_nn = "#6B7280")) +
      ggplot2::scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, by = 0.2)) +
      ggplot2::labs(
        title = paste0(focus_pid, " beneficial proportion split by edge source"),
        subtitle = "Comparing all valid moves against fq-target, observed-neighbour, and latent-neighbour subsets",
        x = NULL,
        y = "Beneficial proportion",
        fill = "Edge class"
      ) +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(focus_figure_paths$edge_beneficial, p_focus_edge_beneficial, width = 10, height = 6, dpi = 150)
  }

  list(
    focus_pid = focus_pid,
    focus_fit_index_tbl = focus_fit_index_tbl,
    focus_parameter_bundles = focus_parameter_bundles,
    focus_parameter_levels_available = focus_parameter_levels_available,
    focus_landscape_long_tbl = focus_landscape_long_tbl,
    focus_landscape_summary_tbl = focus_landscape_summary_tbl,
    focus_landscape_variation_tbl = focus_landscape_variation_tbl,
    focus_parity_tbl = focus_parity_tbl,
    focus_pairwise_component_tbl = focus_pairwise_component_tbl,
    focus_pairwise_landscape_group_tbl = focus_pairwise_landscape_group_tbl,
    focus_top_shift_tbl = focus_top_shift_tbl,
    focus_beneficial_long_tbl = focus_beneficial_long_tbl,
    focus_beneficial_summary_tbl = focus_beneficial_summary_tbl,
    focus_beneficial_shift_tbl = focus_beneficial_shift_tbl,
    focus_beneficial_proportion_matrix_tbl = focus_beneficial_proportion_matrix_tbl,
    focus_beneficial_valid_n_matrix_tbl = focus_beneficial_valid_n_matrix_tbl,
    focus_state_fitness_summary_tbl = focus_state_fitness_summary_tbl,
    focus_edge_decomposition_summary_tbl = focus_edge_decomposition_summary_tbl,
    focus_nn_prior_diag_summary_tbl = focus_nn_prior_diag_summary_tbl,
    focus_density_artifact_tbl = focus_density_artifact_tbl,
    focus_umap_artifact_tbl = focus_umap_artifact_tbl,
    focus_figure_paths = focus_figure_paths
  )
}

build_benchmark_artifact_index <- function(ctx, focus_results) {
  base_artifacts <- tibble::tibble(
    artifact = c(
      "benchmark_input_index",
      "input_fq_nn_summary",
      "input_fq_group_nn_summary",
      "parameter_spec",
      "parameter_tasks",
      "parameter_fit_results",
      "parameter_all_results",
      "parameter_fit_summary",
      "parameter_selected_fit_index",
      "parameter_beneficial_all_patient_artifacts",
      "parameter_pair_component_summary",
      "parameter_pair_overview",
      "parameter_pair_landscape_group_summary",
      "parameter_pair_sign_summary",
      "parameter_pair_top_landscape_shifts",
      "parameter_global_landscape_overview",
      "parameter_global_landscape_state_summary",
      "nn_identifiability_summary",
      "nn_identifiability_by_replicate",
      "nn_fitness_stability_by_case",
      "nn_fitness_stability_by_child",
      "nn_fitness_weighted_empirical_comparison",
      "nn_diagnostic_figure_artifacts",
      "nn_holdout_summary",
      "nn_holdout_predictions",
      "nn_simulation_summary",
      "nn_simulation_by_child",
      "nn_grf_simulation_summary",
      "nn_grf_simulation_by_lambda",
      "nn_grf_simulation_by_child",
      "nn_grf_simulation_fit_results"
    ),
    path = file.path(
      ctx$tables_dir,
      c(
        "benchmark_input_index.tsv",
        "input_fq_nn_summary.tsv",
        "input_fq_group_nn_summary.tsv",
        "parameter_spec.tsv",
        "parameter_tasks.tsv",
        "parameter_fit_results.tsv",
        "parameter_all_results.tsv",
        "parameter_fit_summary.tsv",
        "parameter_selected_fit_index.tsv",
        "parameter_beneficial_all_patient_artifacts.tsv",
        "parameter_pair_component_summary.tsv",
        "parameter_pair_overview.tsv",
        "parameter_pair_landscape_group_summary.tsv",
        "parameter_pair_sign_summary.tsv",
        "parameter_pair_top_landscape_shifts.tsv",
        "parameter_global_landscape_overview.tsv",
        "parameter_global_landscape_state_summary.tsv",
        "nn_identifiability_summary.tsv",
        "nn_identifiability_by_replicate.tsv",
        "nn_fitness_stability_by_case.tsv",
        "nn_fitness_stability_by_child.tsv",
        "nn_fitness_weighted_empirical_comparison.tsv",
        "nn_diagnostic_figure_artifacts.tsv",
        "nn_holdout_summary.tsv",
        "nn_holdout_predictions.tsv",
        "nn_simulation_summary.tsv",
        "nn_simulation_by_child.tsv",
        "nn_grf_simulation_summary.tsv",
        "nn_grf_simulation_by_lambda.tsv",
        "nn_grf_simulation_by_child.tsv",
        "nn_grf_simulation_fit_results.tsv"
      )
    )
  )

  focus_artifacts <- dplyr::bind_rows(lapply(names(focus_results), function(focus_pid) {
    tibble::tibble(
      artifact = c(
        paste0("focus_", focus_pid, "_fit_index"),
        paste0("focus_", focus_pid, "_landscape_summary"),
        paste0("focus_", focus_pid, "_pairwise_component_summary"),
        paste0("focus_", focus_pid, "_pairwise_landscape_group_summary"),
        paste0("focus_", focus_pid, "_pairwise_top_landscape_shifts"),
        paste0("focus_", focus_pid, "_top_variable_karyotypes"),
        paste0("focus_", focus_pid, "_beneficial_summary"),
        paste0("focus_", focus_pid, "_beneficial_top_shifts"),
        paste0("focus_", focus_pid, "_beneficial_proportion_matrix"),
        paste0("focus_", focus_pid, "_beneficial_valid_n_matrix"),
        paste0("focus_", focus_pid, "_state_fitness_summary"),
        paste0("focus_", focus_pid, "_edge_decomposition_summary"),
        paste0("focus_", focus_pid, "_nn_prior_diagnostics_summary"),
        paste0("focus_", focus_pid, "_umap_artifacts")
      ),
      path = c(
        paste0(focus_table_stem(ctx, focus_pid, "fit_index"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "landscape_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "pairwise_component_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "pairwise_landscape_group_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "pairwise_top_landscape_shifts"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "top_variable_karyotypes"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "beneficial_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "beneficial_top_shifts"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "beneficial_proportion_matrix"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "beneficial_valid_n_matrix"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "state_fitness_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "edge_decomposition_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "nn_prior_diagnostics_summary"), ".tsv"),
        paste0(focus_table_stem(ctx, focus_pid, "umap_artifacts"), ".tsv")
      )
    )
  }))

  dplyr::bind_rows(base_artifacts, focus_artifacts)
}

run_benchmark_pipeline <- function(ctx) {
  meta_tbl <- readxl::read_xlsx(ctx$meta_path)
  input_index_tbl <- build_benchmark_inputs(
    meta_tbl = meta_tbl,
    base_dir = ctx$data_dir,
    input_dir = ctx$input_dir,
    tables_dir = ctx$tables_dir,
    stage_levels = ctx$stage_levels,
    diploid_state = ctx$diploid_state,
    rebuild_inputs = ctx$rebuild_inputs_use,
    patient_subset = ctx$patient_subset_use
  )
  save_table_bundle(input_index_tbl, file.path(ctx$tables_dir, "benchmark_input_index"))

  input_overview <- summarize_input_fq_nn_overview(
    input_index_tbl = input_index_tbl,
    minobs_values = ctx$minobs_values_use,
    diploid_state = ctx$diploid_state
  )
  input_fq_nn_summary_tbl <- input_overview$input_fq_nn_summary_tbl
  input_fq_group_nn_summary_tbl <- input_overview$input_fq_group_nn_summary_tbl
  save_table_bundle(input_fq_nn_summary_tbl, file.path(ctx$tables_dir, "input_fq_nn_summary"))
  save_table_bundle(input_fq_group_nn_summary_tbl, file.path(ctx$tables_dir, "input_fq_group_nn_summary"))
  input_figure_paths <- save_input_overview_figures(
    ctx = ctx,
    input_fq_nn_summary_tbl = input_fq_nn_summary_tbl,
    input_fq_group_nn_summary_tbl = input_fq_group_nn_summary_tbl
  )

  parameter_spec_tbl <- build_parameter_spec_tbl(
    parameter_labels = ctx$parameter_levels_use
  )
  save_table_bundle(parameter_spec_tbl, file.path(ctx$tables_dir, "parameter_spec"))

  parameter_tasks_tbl <- build_parameter_tasks(
    input_index_tbl = input_index_tbl,
    fit_root = ctx$fit_dir,
    minobs_values = ctx$minobs_values_use,
    pm_values = ctx$pm_values_use,
    parameter_spec_tbl = parameter_spec_tbl,
    nn_prior_grid_n = ctx$selected_grid_n_use,
    nn_prior_fit_subset = ctx$nn_prior_fit_subset_use,
    nn_prior_zero_exposure_quantile = ctx$nn_prior_zero_exposure_quantile_use,
    nn_prior_zero_weight_scale = ctx$nn_prior_zero_weight_scale_use,
    nn_prior_zero_weight_cap_ratio = ctx$nn_prior_zero_weight_cap_ratio_use,
    nn_prior_zero_birth_fallback_weight = ctx$nn_prior_zero_birth_fallback_weight_use,
    nn_prior_zero_birth_child_floor = ctx$nn_prior_zero_birth_child_floor_use,
    nn_prior_zero_birth_child_shape = ctx$nn_prior_zero_birth_child_shape_use,
    nn_prior_zero_birth_replicate_floor = ctx$nn_prior_zero_birth_replicate_floor_use,
    nn_prior_zero_birth_replicate_shape = ctx$nn_prior_zero_birth_replicate_shape_use,
    nn_prior_two_step_support = ctx$nn_prior_two_step_support_use,
    nn_prior_two_step_support_min = ctx$nn_prior_two_step_support_min_use,
    nn_prior_two_step_cap_floor = ctx$nn_prior_two_step_cap_floor_use,
    cohort_transition_version = ctx$cohort_transition_version_use,
    cohort_contextual_apply_to = ctx$cohort_contextual_apply_to_use,
    cohort_context_keep_baseline_when_sparse = ctx$cohort_context_keep_baseline_when_sparse_use,
    cohort_context_lambda_sparse_unknown = ctx$cohort_context_lambda_sparse_unknown_use,
    nboot = ctx$nboot_use,
    n0 = ctx$n0_use,
    nb = ctx$nb_use,
    benchmark_seed = ctx$benchmark_seed_use,
    correct_efflux = ctx$correct_efflux_use,
    force_refit = ctx$force_refit_use
  )
  save_table_bundle(parameter_tasks_tbl, file.path(ctx$tables_dir, "parameter_tasks"))

  parameter_results_path <- file.path(ctx$tables_dir, "parameter_fit_results")
  if (ctx$run_benchmark_use && nrow(parameter_tasks_tbl)) {
    regular_task_tbl <- parameter_tasks_tbl %>%
      dplyr::filter(nn_prior != "cohort_transition")
    cohort_task_tbl <- parameter_tasks_tbl %>%
      dplyr::filter(nn_prior == "cohort_transition")

    regular_results_tbl <- run_task_table_parallel(
      task_tbl = regular_task_tbl,
      n_cores = ctx$n_cores_use,
      diploid_state = ctx$diploid_state
    )
    cohort_results_tbl <- run_cohort_transition_tasks(
      task_tbl = cohort_task_tbl,
      base_results_tbl = regular_results_tbl,
      fit_root = ctx$fit_dir,
      nboot = ctx$nboot_use,
      n0 = ctx$n0_use,
      nb = ctx$nb_use,
      correct_efflux = ctx$correct_efflux_use,
      n_cores = ctx$n_cores_use,
      cohort_transition_version = ctx$cohort_transition_version_use,
      cohort_contextual_apply_to = ctx$cohort_contextual_apply_to_use,
      cohort_context_keep_baseline_when_sparse = ctx$cohort_context_keep_baseline_when_sparse_use,
      cohort_context_lambda_sparse_unknown = ctx$cohort_context_lambda_sparse_unknown_use,
      diploid_state = ctx$diploid_state,
      force_refit = ctx$force_refit_use
    )
    parameter_results_all_tbl <- dplyr::bind_rows(regular_results_tbl, cohort_results_tbl)
  } else {
    parameter_results_all_tbl <- load_saved_table(parameter_results_path)
    if (is.null(parameter_results_all_tbl)) {
      parameter_results_all_tbl <- tibble::tibble()
    }
    parameter_results_all_tbl <- adopt_existing_task_outputs(
      task_tbl = parameter_tasks_tbl,
      fit_results_tbl = parameter_results_all_tbl
    )
  }

  parameter_results_all_tbl <- reconcile_fit_results_tbl(
    fit_results_tbl = parameter_results_all_tbl,
    task_tbl = parameter_tasks_tbl
  ) %>%
    dplyr::filter(parameter_label %in% ctx$parameter_levels_use)
  save_table_bundle(parameter_results_all_tbl, parameter_results_path)

  if (nrow(parameter_results_all_tbl)) {
    parameter_results_all_tbl <- parameter_results_all_tbl %>%
      dplyr::arrange(factor(patient_id, levels = sort_pid_levels(patient_id)), minobs, pm, factor(parameter_label, levels = ctx$parameter_levels_use))
  }
  save_table_bundle(parameter_results_all_tbl, file.path(ctx$tables_dir, "parameter_all_results"))

  parameter_fit_summary_tbl <- summarize_fit_results(parameter_results_all_tbl, group_cols = c("parameter_label"))
  save_table_bundle(parameter_fit_summary_tbl, file.path(ctx$tables_dir, "parameter_fit_summary"))

  parameter_positive_xval_fit_tbl <- select_positive_xval_best_parameter_fit_tbl(
    parameter_results_all_tbl,
    parameter_levels = ctx$parameter_levels_use
  )
  save_table_bundle(parameter_positive_xval_fit_tbl, file.path(ctx$tables_dir, "parameter_positive_xval_fit_index"))

  retained_patient_ids <- select_common_positive_xval_patients(
    parameter_positive_xval_fit_tbl,
    parameter_levels = ctx$parameter_levels_use
  )
  retained_patient_tbl <- tibble::tibble(patient_id = retained_patient_ids)
  save_table_bundle(retained_patient_tbl, file.path(ctx$tables_dir, "retained_patient_index"))

  parameter_selected_fit_tbl <- parameter_positive_xval_fit_tbl %>%
    dplyr::filter(patient_id %in% retained_patient_ids)
  save_table_bundle(parameter_selected_fit_tbl, file.path(ctx$tables_dir, "parameter_selected_fit_index"))

  all_sample_landscape_long_tbl <- build_selected_landscape_long_tbl(
    parameter_positive_xval_fit_tbl,
    beneficial_move_levels = ctx$beneficial_move_levels,
    parameter_levels = ctx$parameter_levels_use
  )
  retained_landscape_sample_summary_tbl <- summarize_selected_landscape_sample_tbl(all_sample_landscape_long_tbl)
  save_table_bundle(retained_landscape_sample_summary_tbl, file.path(ctx$tables_dir, "retained_landscape_sample_summary"))

  retained_xval_scatter_tbl <- build_selected_xval_scatter_tbl(
    parameter_positive_xval_fit_tbl,
    beneficial_move_levels = ctx$beneficial_move_levels,
    benchmark_seed = ctx$benchmark_seed_use,
    parameter_levels = ctx$parameter_levels_use
  )
  save_table_bundle(retained_xval_scatter_tbl, file.path(ctx$tables_dir, "retained_xval_scatter_points"))

  retained_umap_res <- build_selected_umap_long_tbl(
    landscape_long = all_sample_landscape_long_tbl,
    benchmark_seed = ctx$benchmark_seed_use,
    diploid_state = ctx$diploid_state
  )
  all_sample_umap_artifact <- save_retained_summary_figures(
    ctx = ctx,
    retained_fit_tbl = parameter_positive_xval_fit_tbl,
    retained_landscape_long_tbl = all_sample_landscape_long_tbl,
    retained_xval_scatter_tbl = retained_xval_scatter_tbl,
    retained_umap_res = retained_umap_res
  )
  all_sample_landscape_artifact_tbl <- all_sample_umap_artifact$landscape_artifact_tbl
  save_table_bundle(all_sample_landscape_artifact_tbl, file.path(ctx$tables_dir, "all_sample_landscape_artifacts"))
  all_sample_xval_artifact_tbl <- all_sample_umap_artifact$xval_artifact_tbl
  save_table_bundle(all_sample_xval_artifact_tbl, file.path(ctx$tables_dir, "all_sample_xval_artifacts"))
  all_sample_umap_artifact_tbl <- all_sample_umap_artifact$umap_artifact_tbl
  save_table_bundle(all_sample_umap_artifact_tbl, file.path(ctx$tables_dir, "all_sample_umap_artifacts"))

  parameter_selected_landscapes <- load_selected_landscapes_by_parameter(
    parameter_positive_xval_fit_tbl,
    beneficial_move_levels = ctx$beneficial_move_levels,
    parameter_levels = ctx$parameter_levels_use
  )
  parameter_beneficial_artifact_tbl <- build_parameter_beneficial_artifacts(
    selected_landscapes_by_parameter = parameter_selected_landscapes,
    beneficial_move_levels = ctx$beneficial_move_levels,
    parameter_levels = ctx$parameter_levels_use,
    tables_dir = ctx$tables_dir,
    figures_dir = ctx$figures_dir
  )
  save_table_bundle(parameter_beneficial_artifact_tbl, file.path(ctx$tables_dir, "parameter_beneficial_all_patient_artifacts"))

  parameter_pair_results <- build_pairwise_comparisons(
    results_tbl = parameter_results_all_tbl,
    pair_values = ctx$parameter_levels_use,
    setting_col = "parameter_label",
    comparison_set = "parameter",
    beneficial_move_levels = ctx$beneficial_move_levels,
    top_n = ctx$top_shift_n_use,
    include_posterior = ctx$include_posterior_use
  )
  save_table_bundle(parameter_pair_results$component, file.path(ctx$tables_dir, "parameter_pair_component_summary"))
  save_table_bundle(parameter_pair_results$sign, file.path(ctx$tables_dir, "parameter_pair_sign_summary"))
  save_table_bundle(parameter_pair_results$landscape_group, file.path(ctx$tables_dir, "parameter_pair_landscape_group_summary"))
  save_table_bundle(parameter_pair_results$top_shift, file.path(ctx$tables_dir, "parameter_pair_top_landscape_shifts"))

  parameter_pair_overview_tbl <- if (nrow(parameter_pair_results$component)) {
    parameter_pair_results$component %>%
      dplyr::filter(metric %in% c("landscape_mean", "bootstrap_nn_fitness", "xval_r2")) %>%
      dplyr::group_by(lhs_label, rhs_label, metric) %>%
      dplyr::summarise(
        n_pairs = dplyr::n(),
        mean_abs_diff = mean(mean_abs_diff, na.rm = TRUE),
        max_abs_diff = max(max_abs_diff, na.rm = TRUE),
        mean_correlation = mean(correlation, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble::tibble()
  }
  save_table_bundle(parameter_pair_overview_tbl, file.path(ctx$tables_dir, "parameter_pair_overview"))

  parameter_top_shift_preview_tbl <- if (nrow(parameter_pair_results$top_shift)) {
    parameter_pair_results$top_shift %>%
      dplyr::group_by(lhs_label, rhs_label) %>%
      dplyr::slice_head(n = min(10L, ctx$top_shift_n_use)) %>%
      dplyr::ungroup()
  } else {
    tibble::tibble()
  }

  parameter_global_landscape_overview_tbl <- if (nrow(parameter_pair_results$component)) {
    parameter_pair_results$component %>%
      dplyr::filter(metric %in% c("landscape_mean", "landscape_median", "landscape_sd", "beneficial_proportion", "beneficial_valid_n")) %>%
      dplyr::group_by(lhs_label, rhs_label, metric) %>%
      dplyr::summarise(
        n_pairs = dplyr::n(),
        mean_diff = mean(mean_diff, na.rm = TRUE),
        mean_abs_diff = mean(mean_abs_diff, na.rm = TRUE),
        max_abs_diff = max(max_abs_diff, na.rm = TRUE),
        mean_correlation = mean(correlation, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble::tibble()
  }
  save_table_bundle(parameter_global_landscape_overview_tbl, file.path(ctx$tables_dir, "parameter_global_landscape_overview"))

  parameter_global_landscape_state_tbl <- if (nrow(parameter_pair_results$landscape_group)) {
    parameter_pair_results$landscape_group %>%
      dplyr::filter(metric == "landscape_mean") %>%
      dplyr::group_by(lhs_label, rhs_label, state_group) %>%
      dplyr::summarise(
        n_pairs = dplyr::n(),
        mean_diff = mean(mean_diff, na.rm = TRUE),
        mean_abs_diff = mean(mean_abs_diff, na.rm = TRUE),
        max_abs_diff = max(max_abs_diff, na.rm = TRUE),
        mean_correlation = mean(correlation, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble::tibble()
  }
  save_table_bundle(parameter_global_landscape_state_tbl, file.path(ctx$tables_dir, "parameter_global_landscape_state_summary"))

  parameter_figure_paths <- save_parameter_figures(
    ctx = ctx,
    parameter_results_all_tbl = parameter_results_all_tbl,
    parameter_fit_summary_tbl = parameter_fit_summary_tbl,
    parameter_global_landscape_overview_tbl = parameter_global_landscape_overview_tbl,
    parameter_global_landscape_state_tbl = parameter_global_landscape_state_tbl
  )

  nn_diagnostic_results <- build_benchmark_nn_diagnostics(
    ctx = ctx,
    parameter_results_all_tbl = parameter_results_all_tbl,
    input_index_tbl = input_index_tbl,
    parameter_spec_tbl = parameter_spec_tbl
  )

  focus_results <- setNames(lapply(retained_patient_ids, function(focus_pid) {
    build_focus_outputs(
      focus_pid = focus_pid,
      ctx = ctx,
      parameter_selected_fit_tbl = parameter_selected_fit_tbl
    )
  }), retained_patient_ids)

  artifact_index_tbl <- build_benchmark_artifact_index(ctx, focus_results)
  save_table_bundle(artifact_index_tbl, file.path(ctx$tables_dir, "artifact_index"))

  list(
    ctx = ctx,
    input_index_tbl = input_index_tbl,
    input_fq_nn_summary_tbl = input_fq_nn_summary_tbl,
    input_fq_group_nn_summary_tbl = input_fq_group_nn_summary_tbl,
    input_figure_paths = input_figure_paths,
    parameter_spec_tbl = parameter_spec_tbl,
    parameter_tasks_tbl = parameter_tasks_tbl,
    parameter_results_all_tbl = parameter_results_all_tbl,
    parameter_fit_summary_tbl = parameter_fit_summary_tbl,
    retained_patient_tbl = retained_patient_tbl,
    parameter_positive_xval_fit_tbl = parameter_positive_xval_fit_tbl,
    parameter_selected_fit_tbl = parameter_selected_fit_tbl,
    retained_landscape_long_tbl = all_sample_landscape_long_tbl,
    retained_landscape_sample_summary_tbl = retained_landscape_sample_summary_tbl,
    retained_xval_scatter_tbl = retained_xval_scatter_tbl,
    all_sample_landscape_artifact_tbl = all_sample_landscape_artifact_tbl,
    all_sample_xval_artifact_tbl = all_sample_xval_artifact_tbl,
    all_sample_umap_artifact_tbl = all_sample_umap_artifact_tbl,
    parameter_beneficial_artifact_tbl = parameter_beneficial_artifact_tbl,
    parameter_pair_results = parameter_pair_results,
    parameter_pair_overview_tbl = parameter_pair_overview_tbl,
    parameter_top_shift_preview_tbl = parameter_top_shift_preview_tbl,
    parameter_global_landscape_overview_tbl = parameter_global_landscape_overview_tbl,
    parameter_global_landscape_state_tbl = parameter_global_landscape_state_tbl,
    parameter_figure_paths = parameter_figure_paths,
    nn_diagnostic_results = nn_diagnostic_results,
    nn_identifiability_summary_tbl = nn_diagnostic_results$identifiability_summary_tbl,
    nn_identifiability_replicate_tbl = nn_diagnostic_results$identifiability_replicate_tbl,
    nn_stability_case_tbl = nn_diagnostic_results$stability_case_tbl,
    nn_stability_child_tbl = nn_diagnostic_results$stability_child_tbl,
    nn_stability_weighted_empirical_tbl = nn_diagnostic_results$stability_weighted_empirical_tbl,
    nn_diagnostic_figure_artifact_tbl = nn_diagnostic_results$figure_artifact_tbl,
    nn_holdout_summary_tbl = nn_diagnostic_results$holdout_summary_tbl,
    nn_holdout_prediction_tbl = nn_diagnostic_results$holdout_prediction_tbl,
    nn_simulation_summary_tbl = nn_diagnostic_results$simulation_summary_tbl,
    nn_simulation_child_tbl = nn_diagnostic_results$simulation_child_tbl,
    nn_grf_simulation_summary_tbl = nn_diagnostic_results$grf_simulation_summary_tbl,
    nn_grf_simulation_by_lambda_tbl = nn_diagnostic_results$grf_simulation_by_lambda_tbl,
    nn_grf_simulation_child_tbl = nn_diagnostic_results$grf_simulation_child_tbl,
    nn_grf_simulation_fit_tbl = nn_diagnostic_results$grf_simulation_fit_tbl,
    focus_results = focus_results,
    artifact_index_tbl = artifact_index_tbl
  )
}
