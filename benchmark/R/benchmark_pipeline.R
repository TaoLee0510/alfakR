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

focus_table_stem <- function(ctx, focus_pid, stem) {
  file.path(ctx$tables_dir, paste0("focus_", focus_pid, "_", stem))
}

build_parameter_figure_paths <- function(ctx) {
  list(
    runtime = file.path(ctx$figures_dir, "runtime_by_parameter_label.png"),
    xval = file.path(ctx$figures_dir, "xval_by_parameter_label.png"),
    global_diff = file.path(ctx$figures_dir, "global_landscape_and_beneficial_diff.png"),
    global_state_diff = file.path(ctx$figures_dir, "global_landscape_state_diff.png")
  )
}

build_focus_figure_paths <- function(ctx, focus_pid) {
  list(
    density = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_landscape_density_by_parameter_label.png")),
    parity = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_landscape_parity_by_parameter_label.png")),
    beneficial_heatmap = file.path(ctx$figures_dir, paste0("focus_", focus_pid, "_beneficial_proportion_heatmap.png"))
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
  nn_prior_zero_birth_fallback_weight_use <- suppressWarnings(as.numeric(params$nn_prior_zero_birth_fallback_weight))
  if (!is.finite(nn_prior_zero_birth_fallback_weight_use) ||
      nn_prior_zero_birth_fallback_weight_use < 0 ||
      nn_prior_zero_birth_fallback_weight_use > 1) {
    stop("nn_prior_zero_birth_fallback_weight must be a finite number in [0, 1].")
  }

  force_refit_use <- isTRUE(params$force_refit)
  rebuild_inputs_use <- isTRUE(params$rebuild_inputs)
  run_benchmark_use <- isTRUE(params$run_benchmark)
  include_posterior_use <- isTRUE(params$include_posterior_comparison)
  render_figures_use <- isTRUE(params$render_figures)
  top_shift_n_use <- max(1L, as.integer(params$top_shift_n))

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
    force_refit_use = force_refit_use,
    rebuild_inputs_use = rebuild_inputs_use,
    run_benchmark_use = run_benchmark_use,
    include_posterior_use = include_posterior_use,
    render_figures_use = render_figures_use,
    top_shift_n_use = top_shift_n_use,
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

build_focus_outputs <- function(focus_pid, ctx, parameter_results_all_tbl, parameter_pair_results) {
  focus_fit_index_tbl <- parameter_results_all_tbl %>%
    dplyr::filter(
      patient_id == focus_pid,
      minobs == ctx$focus_minobs_use,
      dplyr::near(pm, ctx$focus_pm_use),
      parameter_label %in% ctx$parameter_levels_use
    ) %>%
    dplyr::distinct(parameter_label, .keep_all = TRUE) %>%
    dplyr::arrange(match(parameter_label, ctx$parameter_levels_use))
  save_table_bundle(focus_fit_index_tbl, focus_table_stem(ctx, focus_pid, "fit_index"))

  focus_parameter_bundles <- build_focus_parameter_bundles(
    results_tbl = parameter_results_all_tbl,
    focus_pid = focus_pid,
    focus_minobs = ctx$focus_minobs_use,
    focus_pm = ctx$focus_pm_use,
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

  focus_pairwise_component_tbl <- if (nrow(parameter_pair_results$component)) {
    parameter_pair_results$component %>%
      dplyr::filter(
        patient_id == focus_pid,
        minobs == ctx$focus_minobs_use,
        dplyr::near(pm, ctx$focus_pm_use),
        metric %in% c("landscape_mean", "landscape_median", "landscape_sd", "xval_r2", "beneficial_proportion", "beneficial_valid_n")
      ) %>%
      dplyr::arrange(match(lhs_label, focus_parameter_levels_available), match(rhs_label, focus_parameter_levels_available), metric)
  } else {
    tibble::tibble()
  }
  save_table_bundle(focus_pairwise_component_tbl, focus_table_stem(ctx, focus_pid, "pairwise_component_summary"))

  focus_pairwise_landscape_group_tbl <- if (nrow(parameter_pair_results$landscape_group)) {
    parameter_pair_results$landscape_group %>%
      dplyr::filter(
        patient_id == focus_pid,
        minobs == ctx$focus_minobs_use,
        dplyr::near(pm, ctx$focus_pm_use),
        metric == "landscape_mean"
      ) %>%
      dplyr::arrange(match(lhs_label, focus_parameter_levels_available), match(rhs_label, focus_parameter_levels_available), state_group)
  } else {
    tibble::tibble()
  }
  save_table_bundle(focus_pairwise_landscape_group_tbl, focus_table_stem(ctx, focus_pid, "pairwise_landscape_group_summary"))

  focus_top_shift_tbl <- if (nrow(parameter_pair_results$top_shift)) {
    parameter_pair_results$top_shift %>%
      dplyr::filter(
        patient_id == focus_pid,
        minobs == ctx$focus_minobs_use,
        dplyr::near(pm, ctx$focus_pm_use)
      ) %>%
      dplyr::group_by(lhs_label, rhs_label) %>%
      dplyr::slice_head(n = min(10L, ctx$top_shift_n_use)) %>%
      dplyr::ungroup()
  } else {
    tibble::tibble()
  }
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
  save_table_bundle(focus_umap_artifact_tbl, focus_table_stem(ctx, focus_pid, "umap_artifacts"))

  focus_figure_paths <- build_focus_figure_paths(ctx, focus_pid)

  if (ctx$render_figures_use && nrow(focus_landscape_long_tbl)) {
    focus_landscape_density_long_tbl <- focus_landscape_long_tbl %>%
      tidyr::pivot_longer(cols = c(mean, median), names_to = "metric", values_to = "value") %>%
      dplyr::mutate(parameter_label = factor(parameter_label, levels = focus_parameter_levels_available))

    p_focus_landscape_density <- ggplot2::ggplot(
      focus_landscape_density_long_tbl,
      ggplot2::aes(x = value, color = parameter_label, fill = parameter_label)
    ) +
      ggplot2::geom_density(alpha = 0.18, adjust = 1) +
      ggplot2::facet_wrap(~ metric, scales = "free", ncol = 2) +
      ggplot2::labs(
        title = paste0(focus_pid, " landscape distributions by parameter label"),
        subtitle = paste0("minobs = ", ctx$focus_minobs_use, ", pm = ", pm_to_label(ctx$focus_pm_use)),
        x = "Landscape value",
        y = "Density",
        color = "Parameter label",
        fill = "Parameter label"
      ) +
      ggplot2::theme_bw(base_size = 12)
    ggplot2::ggsave(focus_figure_paths$density, p_focus_landscape_density, width = 10, height = 5, dpi = 150)
  }

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
        subtitle = paste0("minobs = ", ctx$focus_minobs_use, ", pm = ", pm_to_label(ctx$focus_pm_use)),
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
    focus_umap_artifact_tbl = focus_umap_artifact_tbl,
    focus_figure_paths = focus_figure_paths
  )
}

build_benchmark_artifact_index <- function(ctx, focus_results) {
  base_artifacts <- tibble::tibble(
    artifact = c(
      "benchmark_input_index",
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
      "parameter_global_landscape_state_summary"
    ),
    path = file.path(
      ctx$tables_dir,
      c(
        "benchmark_input_index.tsv",
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
        "parameter_global_landscape_state_summary.tsv"
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
    parameter_results_all_tbl <- run_task_table_parallel(
      task_tbl = parameter_tasks_tbl,
      n_cores = ctx$n_cores_use,
      diploid_state = ctx$diploid_state
    )
  } else {
    parameter_results_all_tbl <- load_saved_table(parameter_results_path)
    if (is.null(parameter_results_all_tbl)) {
      parameter_results_all_tbl <- tibble::tibble()
    }
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

  parameter_selected_fit_tbl <- select_best_parameter_fit_tbl(
    parameter_results_all_tbl,
    parameter_levels = ctx$parameter_levels_use
  )
  save_table_bundle(parameter_selected_fit_tbl, file.path(ctx$tables_dir, "parameter_selected_fit_index"))

  parameter_selected_landscapes <- load_selected_landscapes_by_parameter(
    parameter_selected_fit_tbl,
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

  focus_results <- setNames(lapply(ctx$focus_pids_use, function(focus_pid) {
    build_focus_outputs(
      focus_pid = focus_pid,
      ctx = ctx,
      parameter_results_all_tbl = parameter_results_all_tbl,
      parameter_pair_results = parameter_pair_results
    )
  }), ctx$focus_pids_use)

  artifact_index_tbl <- build_benchmark_artifact_index(ctx, focus_results)
  save_table_bundle(artifact_index_tbl, file.path(ctx$tables_dir, "artifact_index"))

  list(
    ctx = ctx,
    input_index_tbl = input_index_tbl,
    parameter_spec_tbl = parameter_spec_tbl,
    parameter_tasks_tbl = parameter_tasks_tbl,
    parameter_results_all_tbl = parameter_results_all_tbl,
    parameter_fit_summary_tbl = parameter_fit_summary_tbl,
    parameter_selected_fit_tbl = parameter_selected_fit_tbl,
    parameter_beneficial_artifact_tbl = parameter_beneficial_artifact_tbl,
    parameter_pair_results = parameter_pair_results,
    parameter_pair_overview_tbl = parameter_pair_overview_tbl,
    parameter_top_shift_preview_tbl = parameter_top_shift_preview_tbl,
    parameter_global_landscape_overview_tbl = parameter_global_landscape_overview_tbl,
    parameter_global_landscape_state_tbl = parameter_global_landscape_state_tbl,
    parameter_figure_paths = parameter_figure_paths,
    focus_results = focus_results,
    artifact_index_tbl = artifact_index_tbl
  )
}
