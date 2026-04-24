numeric_or_na <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- NA_real_
  out
}

safe_median <- function(x) {
  x <- numeric_or_na(x)
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  x <- numeric_or_na(x)
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_fraction <- function(x) {
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x %in% TRUE, na.rm = TRUE)
}

extract_nn_prior_diagnostics_tbl <- function(fit_tbl) {
  if (is.null(fit_tbl) || !nrow(fit_tbl)) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(lapply(seq_len(nrow(fit_tbl)), function(i) {
    rr <- fit_tbl[i, , drop = FALSE]
    if (!identical(as.character(rr$status), "ok")) {
      return(tibble::tibble())
    }

    boot_path <- if ("bootstrap_path" %in% names(rr) && !is.na(rr$bootstrap_path) && nzchar(rr$bootstrap_path)) {
      as.character(rr$bootstrap_path)
    } else {
      file.path(as.character(rr$outdir), "bootstrap_res.Rds")
    }
    boot_obj <- safe_read_rds(boot_path)
    diag_tbl <- boot_obj$nn_prior_diagnostics
    if (is.null(diag_tbl) || !nrow(diag_tbl)) {
      return(tibble::tibble())
    }

    diag_tbl <- tibble::as_tibble(diag_tbl)
    diag_tbl$bootstrap_iter <- seq_len(nrow(diag_tbl))
    diag_tbl %>%
      dplyr::mutate(
        patient_id = as.character(rr$patient_id),
        parameter_label = as.character(rr$parameter_label),
        nn_prior = as.character(rr$nn_prior),
        minobs = as.integer(rr$minobs),
        pm = as.numeric(rr$pm),
        outdir = as.character(rr$outdir),
        .before = 1
      )
  }))
}

classify_nn_identifiability <- function(observed_anchor_rate,
                                        median_n_observed_children,
                                        prior_source_sample_pooled_frac,
                                        median_zero_effective_mass,
                                        median_n_zero_children_with_two_step_support) {
  dplyr::case_when(
    is.finite(observed_anchor_rate) &
      observed_anchor_rate >= 0.8 &
      is.finite(median_n_observed_children) &
      median_n_observed_children >= 3 ~ "learnable",
    (is.finite(observed_anchor_rate) & observed_anchor_rate > 0) |
      (is.finite(prior_source_sample_pooled_frac) & prior_source_sample_pooled_frac > 0) |
      (is.finite(median_zero_effective_mass) & median_zero_effective_mass > 0.5) |
      (is.finite(median_n_zero_children_with_two_step_support) & median_n_zero_children_with_two_step_support > 0) ~ "weakly_learnable",
    TRUE ~ "unidentifiable"
  )
}

build_nn_identifiability_tables <- function(fit_tbl) {
  replicate_tbl <- extract_nn_prior_diagnostics_tbl(fit_tbl)
  if (!nrow(replicate_tbl)) {
    return(list(
      replicate_tbl = tibble::tibble(),
      summary_tbl = tibble::tibble()
    ))
  }

  replicate_tbl <- replicate_tbl %>%
    dplyr::mutate(
      n_observed_children = numeric_or_na(n_observed_children),
      n_zero_children_total = numeric_or_na(n_zero_children_total),
      n_zero_children_retained = numeric_or_na(n_zero_children_retained),
      sum_observed_weight = numeric_or_na(sum_observed_weight),
      sum_zero_weight_final = numeric_or_na(sum_zero_weight_final),
      zero_weight_ratio = dplyr::if_else(
        is.finite(sum_observed_weight) & sum_observed_weight > 0,
        sum_zero_weight_final / sum_observed_weight,
        NA_real_
      ),
      zero_effective_mass_used = numeric_or_na(zero_effective_mass_used),
      n_zero_children_with_two_step_support = numeric_or_na(n_zero_children_with_two_step_support),
      map_delta_lower_boundary_rate = numeric_or_na(map_delta_lower_boundary_rate),
      map_delta_upper_boundary_rate = numeric_or_na(map_delta_upper_boundary_rate),
      nn_prior_source_used = as.character(nn_prior_source_used),
      has_observed_anchor = n_observed_children > 0
    )

  summary_tbl <- replicate_tbl %>%
    dplyr::group_by(parameter_label, nn_prior, patient_id, minobs, pm) %>%
    dplyr::summarise(
      n_bootstrap = dplyr::n(),
      observed_anchor_rate = safe_fraction(has_observed_anchor),
      median_n_observed_children = safe_median(n_observed_children),
      median_n_zero_children_total = safe_median(n_zero_children_total),
      median_n_zero_children_retained = safe_median(n_zero_children_retained),
      median_zero_weight_ratio = safe_median(zero_weight_ratio),
      median_zero_effective_mass = safe_median(zero_effective_mass_used),
      prior_source_observed_replicate_frac = safe_fraction(nn_prior_source_used == "observed_replicate"),
      prior_source_sample_pooled_frac = safe_fraction(nn_prior_source_used == "sample_pooled"),
      prior_source_none_frac = safe_fraction(nn_prior_source_used == "none"),
      median_n_zero_children_with_two_step_support = safe_median(n_zero_children_with_two_step_support),
      median_map_delta_lower_boundary_rate = safe_median(map_delta_lower_boundary_rate),
      median_map_delta_upper_boundary_rate = safe_median(map_delta_upper_boundary_rate),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      nn_identifiability_class = classify_nn_identifiability(
        observed_anchor_rate = observed_anchor_rate,
        median_n_observed_children = median_n_observed_children,
        prior_source_sample_pooled_frac = prior_source_sample_pooled_frac,
        median_zero_effective_mass = median_zero_effective_mass,
        median_n_zero_children_with_two_step_support = median_n_zero_children_with_two_step_support
      ),
      patient_id = factor(patient_id, levels = sort_pid_levels(patient_id))
    ) %>%
    dplyr::arrange(parameter_label, patient_id, minobs, pm)

  list(
    replicate_tbl = replicate_tbl,
    summary_tbl = summary_tbl
  )
}

summarize_bootstrap_matrix_by_child <- function(mat) {
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    finite_vals <- vals[is.finite(vals)]
    tibble::tibble(
      k = colnames(mat)[[j]],
      n_bootstrap = length(vals),
      finite_fraction = if (length(vals)) length(finite_vals) / length(vals) else NA_real_,
      bootstrap_mean = if (length(finite_vals)) mean(finite_vals) else NA_real_,
      bootstrap_sd = if (length(finite_vals) >= 2L) stats::sd(finite_vals) else NA_real_,
      bootstrap_iqr = if (length(finite_vals)) stats::IQR(finite_vals) else NA_real_,
      bootstrap_q025 = if (length(finite_vals)) as.numeric(stats::quantile(finite_vals, 0.025, names = FALSE, na.rm = TRUE)) else NA_real_,
      bootstrap_q975 = if (length(finite_vals)) as.numeric(stats::quantile(finite_vals, 0.975, names = FALSE, na.rm = TRUE)) else NA_real_
    )
  }))
}

build_nn_fitness_stability_tables <- function(fit_tbl) {
  if (is.null(fit_tbl) || !nrow(fit_tbl)) {
    empty <- tibble::tibble()
    return(list(child_tbl = empty, case_tbl = empty, weighted_empirical_tbl = empty))
  }

  child_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(fit_tbl)), function(i) {
    rr <- fit_tbl[i, , drop = FALSE]
    if (!identical(as.character(rr$status), "ok")) {
      return(tibble::tibble())
    }
    boot_obj <- safe_read_rds(file.path(as.character(rr$outdir), "bootstrap_res.Rds"))
    child_summary <- summarize_bootstrap_matrix_by_child(boot_obj$nn_fitness)
    if (!nrow(child_summary)) {
      return(tibble::tibble())
    }
    child_summary %>%
      dplyr::mutate(
        patient_id = as.character(rr$patient_id),
        parameter_label = as.character(rr$parameter_label),
        nn_prior = as.character(rr$nn_prior),
        minobs = as.integer(rr$minobs),
        pm = as.numeric(rr$pm),
        .before = 1
      )
  }))

  if (!nrow(child_tbl)) {
    empty <- tibble::tibble()
    return(list(child_tbl = empty, case_tbl = empty, weighted_empirical_tbl = empty))
  }

  case_tbl <- child_tbl %>%
    dplyr::group_by(parameter_label, nn_prior, patient_id, minobs, pm) %>%
    dplyr::summarise(
      n_nn_children = dplyr::n(),
      median_finite_fraction = safe_median(finite_fraction),
      median_nn_bootstrap_sd = safe_median(bootstrap_sd),
      median_nn_bootstrap_iqr = safe_median(bootstrap_iqr),
      mean_nn_bootstrap_sd = safe_mean(bootstrap_sd),
      mean_nn_bootstrap_iqr = safe_mean(bootstrap_iqr),
      .groups = "drop"
    ) %>%
    dplyr::mutate(patient_id = factor(patient_id, levels = sort_pid_levels(patient_id))) %>%
    dplyr::arrange(parameter_label, patient_id, minobs, pm)

  weighted_empirical_tbl <- child_tbl %>%
    dplyr::filter(parameter_label %in% c(
      "nn_prior_empirical_censored_weighted",
      "nn_prior_empirical_two_shell",
      "nn_prior_empirical"
    )) %>%
    dplyr::select(patient_id, minobs, pm, parameter_label, k, bootstrap_mean, bootstrap_sd, bootstrap_iqr) %>%
    tidyr::pivot_wider(
      names_from = parameter_label,
      values_from = c(bootstrap_mean, bootstrap_sd, bootstrap_iqr),
      names_sep = "__"
    )

  pair_specs <- tibble::tribble(
    ~lhs_label, ~rhs_label, ~comparison_label,
    "nn_prior_empirical_censored_weighted", "nn_prior_empirical", "weighted vs empirical",
    "nn_prior_empirical_two_shell", "nn_prior_empirical", "empirical_two_shell vs empirical",
    "nn_prior_empirical_censored_weighted", "nn_prior_empirical_two_shell", "weighted vs empirical_two_shell"
  )

  if (nrow(weighted_empirical_tbl)) {
    pair_rows <- lapply(seq_len(nrow(pair_specs)), function(i) {
      lhs_label <- pair_specs$lhs_label[[i]]
      rhs_label <- pair_specs$rhs_label[[i]]
      required_cols <- c(
        paste0("bootstrap_mean__", lhs_label),
        paste0("bootstrap_mean__", rhs_label),
        paste0("bootstrap_sd__", lhs_label),
        paste0("bootstrap_sd__", rhs_label),
        paste0("bootstrap_iqr__", lhs_label),
        paste0("bootstrap_iqr__", rhs_label)
      )
      if (!all(required_cols %in% names(weighted_empirical_tbl))) {
        return(NULL)
      }

      weighted_empirical_tbl %>%
        dplyr::transmute(
          patient_id,
          minobs,
          pm,
          k,
          lhs_parameter_label = lhs_label,
          rhs_parameter_label = rhs_label,
          comparison_label = pair_specs$comparison_label[[i]],
          lhs_bootstrap_mean = .data[[paste0("bootstrap_mean__", lhs_label)]],
          rhs_bootstrap_mean = .data[[paste0("bootstrap_mean__", rhs_label)]],
          lhs_bootstrap_sd = .data[[paste0("bootstrap_sd__", lhs_label)]],
          rhs_bootstrap_sd = .data[[paste0("bootstrap_sd__", rhs_label)]],
          lhs_bootstrap_iqr = .data[[paste0("bootstrap_iqr__", lhs_label)]],
          rhs_bootstrap_iqr = .data[[paste0("bootstrap_iqr__", rhs_label)]],
          mean_diff = lhs_bootstrap_mean - rhs_bootstrap_mean,
          sd_ratio = lhs_bootstrap_sd / rhs_bootstrap_sd,
          iqr_ratio = lhs_bootstrap_iqr / rhs_bootstrap_iqr
        )
    })
    weighted_empirical_tbl <- dplyr::bind_rows(pair_rows)
  }

  if (nrow(weighted_empirical_tbl)) {
    weighted_empirical_tbl <- weighted_empirical_tbl %>%
      dplyr::mutate(
        patient_id = factor(as.character(patient_id), levels = sort_pid_levels(patient_id)),
        comparison_label = factor(
          comparison_label,
          levels = c(
            "weighted vs empirical",
            "empirical_two_shell vs empirical",
            "weighted vs empirical_two_shell"
          )
        )
      ) %>%
      dplyr::arrange(comparison_label, patient_id, minobs, pm, k)
  }

  list(
    child_tbl = child_tbl,
    case_tbl = case_tbl,
    weighted_empirical_tbl = weighted_empirical_tbl
  )
}

save_nn_diagnostic_figures <- function(ctx,
                                       nn_identifiability_summary_tbl,
                                       nn_stability_case_tbl,
                                       nn_stability_weighted_empirical_tbl) {
  artifact_tbl <- tibble::tibble(
    artifact = c(
      "nn_identifiability_class",
      "nn_prior_source_fraction",
      "nn_stability_by_method",
      "weighted_empirical_nn_sd_ratio"
    ),
    title = c(
      "NN identifiability class",
      "Weighted prior source fractions",
      "NN bootstrap stability by method",
      "Pairwise NN SD ratio comparisons"
    ),
    png_path = file.path(
      ctx$figures_dir,
      c(
        "nn_identifiability_class.png",
        "nn_prior_source_fraction.png",
        "nn_stability_by_method.png",
        "weighted_empirical_nn_sd_ratio.png"
      )
    )
  )

  if (!ctx$render_figures_use) {
    return(artifact_tbl)
  }

  if (nrow(nn_identifiability_summary_tbl)) {
    plot_tbl <- nn_identifiability_summary_tbl %>%
      dplyr::mutate(
        patient_id = factor(as.character(patient_id), levels = sort_pid_levels(patient_id)),
        minobs = factor(minobs),
        nn_identifiability_class = factor(nn_identifiability_class, levels = c("learnable", "weakly_learnable", "unidentifiable"))
      )
    p_ident <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes(x = minobs, y = patient_id, fill = nn_identifiability_class)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.35) +
      ggplot2::facet_wrap(~ parameter_label, ncol = 2) +
      ggplot2::scale_fill_manual(values = c(learnable = "#2B8CBE", weakly_learnable = "#A6BDDB", unidentifiable = "#FDD49E"), drop = FALSE) +
      ggplot2::labs(x = "MINOBS", y = "Patient", fill = "Class") +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(artifact_tbl$png_path[artifact_tbl$artifact == "nn_identifiability_class"], p_ident, width = 10, height = 6.5, dpi = 150)

    source_tbl <- plot_tbl %>%
      dplyr::filter(parameter_label == "nn_prior_empirical_censored_weighted") %>%
      dplyr::select(patient_id, minobs, prior_source_observed_replicate_frac, prior_source_sample_pooled_frac, prior_source_none_frac) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("prior_source_"),
        names_to = "source",
        values_to = "fraction"
      ) %>%
      dplyr::mutate(
        source = sub("^prior_source_", "", source),
        source = sub("_frac$", "", source),
        panel = paste(patient_id, minobs, sep = " / MINOBS ")
      )
    if (nrow(source_tbl)) {
      source_tbl$panel <- factor(source_tbl$panel, levels = unique(source_tbl$panel))
      p_source <- ggplot2::ggplot(source_tbl, ggplot2::aes(x = panel, y = fraction, fill = source)) +
        ggplot2::geom_col(width = 0.8) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
        ggplot2::scale_fill_manual(values = c(observed_replicate = "#1B9E77", sample_pooled = "#7570B3", none = "#D95F02"), drop = FALSE) +
        ggplot2::labs(x = NULL, y = "Bootstrap fraction", fill = "Prior source") +
        ggplot2::theme_bw(base_size = 11)
      ggplot2::ggsave(artifact_tbl$png_path[artifact_tbl$artifact == "nn_prior_source_fraction"], p_source, width = 9, height = 8, dpi = 150)
    }
  }

  if (nrow(nn_stability_case_tbl)) {
    p_stab <- ggplot2::ggplot(
      nn_stability_case_tbl,
      ggplot2::aes(x = parameter_label, y = median_nn_bootstrap_sd, color = parameter_label)
    ) +
      ggplot2::geom_boxplot(outlier.alpha = 0.25) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Median NN bootstrap SD") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(artifact_tbl$png_path[artifact_tbl$artifact == "nn_stability_by_method"], p_stab, width = 8, height = 5.5, dpi = 150)
  }

  if (nrow(nn_stability_weighted_empirical_tbl) &&
      "sd_ratio" %in% names(nn_stability_weighted_empirical_tbl)) {
    ratio_tbl <- nn_stability_weighted_empirical_tbl %>%
      dplyr::filter(is.finite(sd_ratio), sd_ratio > 0)
    if (nrow(ratio_tbl)) {
      p_ratio <- ggplot2::ggplot(ratio_tbl, ggplot2::aes(x = sd_ratio)) +
        ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
        ggplot2::geom_histogram(bins = 40, fill = "#74A9CF", color = "white") +
        ggplot2::scale_x_log10() +
        ggplot2::facet_wrap(~ comparison_label, ncol = 1, scales = "free_y") +
        ggplot2::labs(x = "Pairwise NN bootstrap SD ratio", y = "NN children") +
        ggplot2::theme_bw(base_size = 11)
      ggplot2::ggsave(artifact_tbl$png_path[artifact_tbl$artifact == "weighted_empirical_nn_sd_ratio"], p_ratio, width = 8.5, height = 8.5, dpi = 150)
    }
  }

  artifact_tbl
}

identify_observed_nn_candidates <- function(input_rds, minobs, diploid_state) {
  yi <- safe_read_rds(input_rds)
  if (is.null(yi) || is.null(yi$x)) {
    return(character(0))
  }
  x <- as.data.frame(yi$x)
  if (diploid_state %in% rownames(x)) {
    x <- x[rownames(x) != diploid_state, , drop = FALSE]
  }
  if (!nrow(x)) {
    return(character(0))
  }

  fq <- alfakR:::get_frequent_karyotypes(x, minobs)
  if (!length(fq)) {
    return(character(0))
  }
  nn_mat <- alfakR:::gen_all_neighbours(fq)
  if (!nrow(nn_mat)) {
    return(character(0))
  }
  nn_ids <- unique(apply(nn_mat, 1, paste, collapse = "."))
  observed_ids <- rownames(x)[rowSums(x, na.rm = TRUE) > 0]
  setdiff(intersect(nn_ids, observed_ids), fq)
}

run_nn_holdout_diagnostics <- function(ctx, input_index_tbl, parameter_spec_tbl) {
  if (!isTRUE(ctx$run_nn_holdout_use)) {
    return(list(summary_tbl = tibble::tibble(), prediction_tbl = tibble::tibble()))
  }

  holdout_root <- file.path(ctx$results_dir, "fits_nn_holdout")
  dir.create(holdout_root, recursive = TRUE, showWarnings = FALSE)
  set.seed(ctx$nn_holdout_seed_use)

  prediction_rows <- list()
  row_idx <- 0L
  for (input_i in seq_len(nrow(input_index_tbl))) {
    input_rr <- input_index_tbl[input_i, , drop = FALSE]
    for (minobs_value in ctx$minobs_values_use) {
      candidates <- identify_observed_nn_candidates(
        input_rds = input_rr$input_rds,
        minobs = minobs_value,
        diploid_state = ctx$diploid_state
      )
      if (length(candidates) < ctx$nn_holdout_min_count_use) {
        next
      }
      n_holdout <- max(1L, floor(length(candidates) * ctx$nn_holdout_fraction_use))
      for (repeat_idx in seq_len(ctx$nn_holdout_repeats_use)) {
        heldout_ids <- sort(sample(candidates, size = min(n_holdout, length(candidates)), replace = FALSE))
        yi_train <- readRDS(input_rr$input_rds)
        yi_train$x <- as.data.frame(yi_train$x)
        heldout_present <- intersect(heldout_ids, rownames(yi_train$x))
        if (!length(heldout_present)) {
          next
        }
        yi_train$x[heldout_present, ] <- 0
        train_rds <- file.path(
          ctx$cache_dir,
          paste0("nn_holdout_", input_rr$patient_id, "_MINOBS_", minobs_value, "_rep", repeat_idx, ".rds")
        )
        saveRDS(yi_train, train_rds)

        for (pm_value in ctx$pm_values_use) {
          for (param_i in seq_len(nrow(parameter_spec_tbl))) {
            param_rr <- parameter_spec_tbl[param_i, , drop = FALSE]
            full_outdir <- task_outdir_parameter(
              root_dir = ctx$fit_dir,
              patient_id = as.character(input_rr$patient_id),
              minobs = minobs_value,
              pm = pm_value,
              parameter_label = as.character(param_rr$parameter_label)
            )
            full_boot <- safe_read_rds(file.path(full_outdir, "bootstrap_res.Rds"))
            target_mat <- full_boot$nn_fitness
            outdir <- build_fit_outdir(
              root_dir = holdout_root,
              patient_id = as.character(input_rr$patient_id),
              path_components = list(
                parameter_label = as.character(param_rr$parameter_label),
                pm = pm_value,
                minobs = minobs_value,
                patient_id = paste0("repeat_", repeat_idx)
              )
            )
            fit_res <- tryCatch(
              run_alfak_fit(
                patient_id = as.character(input_rr$patient_id),
                input_rds = train_rds,
                outdir = outdir,
                minobs = minobs_value,
                pm = pm_value,
                nboot = ctx$nboot_use,
                n0 = ctx$n0_use,
                nb = ctx$nb_use,
                benchmark_seed = ctx$benchmark_seed_use + repeat_idx,
                parameter_label = as.character(param_rr$parameter_label),
                diploid_state = ctx$diploid_state,
                correct_efflux = ctx$correct_efflux_use,
                nn_prior = as.character(param_rr$nn_prior),
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
                force_refit = ctx$force_refit_use
              ),
              error = function(e) tibble::tibble(status = "error", error_message = conditionMessage(e))
            )
            boot_pred <- safe_read_rds(file.path(outdir, "bootstrap_res.Rds"))
            pred_mat <- boot_pred$nn_fitness
            for (kid in heldout_ids) {
              row_idx <- row_idx + 1L
              pred_vals <- if (!is.null(pred_mat) && kid %in% colnames(pred_mat)) {
                as.numeric(pred_mat[, kid])
              } else {
                numeric(0)
              }
              target_vals <- if (!is.null(target_mat) && kid %in% colnames(target_mat)) {
                as.numeric(target_mat[, kid])
              } else {
                numeric(0)
              }
              pred_mean <- if (length(pred_vals) && any(is.finite(pred_vals))) mean(pred_vals, na.rm = TRUE) else NA_real_
              target_mean <- if (length(target_vals) && any(is.finite(target_vals))) mean(target_vals, na.rm = TRUE) else NA_real_
              prediction_rows[[row_idx]] <- tibble::tibble(
                patient_id = as.character(input_rr$patient_id),
                minobs = minobs_value,
                pm = pm_value,
                holdout_repeat = repeat_idx,
                parameter_label = as.character(param_rr$parameter_label),
                nn_prior = as.character(param_rr$nn_prior),
                k = kid,
                status = if ("status" %in% names(fit_res)) as.character(fit_res$status[[1]]) else NA_character_,
                target_mean = target_mean,
                predicted_mean = pred_mean,
                prediction_error = pred_mean - target_mean,
                predicted_sd = if (length(pred_vals) >= 2L) stats::sd(pred_vals, na.rm = TRUE) else NA_real_,
                finite_fraction = if (length(pred_vals)) mean(is.finite(pred_vals)) else 0
              )
            }
          }
        }
      }
    }
  }

  prediction_tbl <- dplyr::bind_rows(prediction_rows)
  summary_tbl <- if (nrow(prediction_tbl)) {
    prediction_tbl %>%
      dplyr::group_by(parameter_label, nn_prior, patient_id, minobs, pm) %>%
      dplyr::summarise(
        n_predictions = dplyr::n(),
        n_finite_predictions = sum(is.finite(predicted_mean)),
        n_targeted_predictions = sum(is.finite(prediction_error)),
        rmse = if (any(is.finite(prediction_error))) sqrt(mean(prediction_error^2, na.rm = TRUE)) else NA_real_,
        mae = if (any(is.finite(prediction_error))) mean(abs(prediction_error), na.rm = TRUE) else NA_real_,
        signed_bias = if (any(is.finite(prediction_error))) mean(prediction_error, na.rm = TRUE) else NA_real_,
        correlation = if (sum(is.finite(predicted_mean) & is.finite(target_mean)) >= 2L) {
          suppressWarnings(stats::cor(predicted_mean, target_mean, use = "complete.obs"))
        } else {
          NA_real_
        },
        median_predicted_sd = safe_median(predicted_sd),
        median_finite_fraction = safe_median(finite_fraction),
        .groups = "drop"
      )
  } else {
    tibble::tibble()
  }

  list(summary_tbl = summary_tbl, prediction_tbl = prediction_tbl)
}

simulate_nn_prior_counts <- function(seed, scenario, depth = c(1200, 1600, 2200)) {
  set.seed(seed)
  k_id <- function(vals) paste(vals, collapse = ".")
  base_a <- rep(2L, 22)
  base_a[1] <- 3L
  base_b <- rep(2L, 22)
  base_b[2] <- 3L
  nn_a <- base_a
  nn_a[3] <- 3L
  nn_b <- base_b
  nn_b[4] <- 1L
  nn_c <- base_a
  nn_c[5] <- 3L
  ids <- c(k_id(base_a), k_id(base_b), k_id(nn_a), k_id(nn_b), k_id(nn_c))

  fitness <- switch(
    scenario,
    sparse_zero_heavy = c(0.04, 0.025, -0.015, -0.025, -0.035),
    moderate_observed = c(0.035, 0.02, 0.005, -0.005, -0.015),
    two_step_supported = c(0.04, 0.015, -0.005, -0.02, 0.0),
    c(0.035, 0.02, -0.01, -0.02, -0.03)
  )
  x0 <- switch(
    scenario,
    sparse_zero_heavy = c(0.62, 0.35, 0.01, 0.01, 0.01),
    moderate_observed = c(0.55, 0.35, 0.04, 0.03, 0.03),
    two_step_supported = c(0.55, 0.34, 0.02, 0.02, 0.07),
    c(0.6, 0.35, 0.02, 0.02, 0.01)
  )
  timepoints <- c(0, 1, 2)
  prob_mat <- sapply(timepoints, function(tt) {
    lv <- log(x0) + fitness * tt
    p <- exp(lv - max(lv))
    p / sum(p)
  })
  counts <- sapply(seq_along(depth), function(j) as.integer(stats::rmultinom(1, size = depth[j], prob = prob_mat[, j])))
  rownames(counts) <- ids
  colnames(counts) <- as.character(timepoints)

  list(
    yi = list(x = counts, dt = 1),
    truth = tibble::tibble(k = ids, true_fitness = fitness, is_true_nn = seq_along(ids) > 2L),
    scenario = scenario
  )
}

run_nn_simulation_diagnostics <- function(ctx, parameter_spec_tbl) {
  if (!isTRUE(ctx$run_nn_simulation_use)) {
    return(list(summary_tbl = tibble::tibble(), child_tbl = tibble::tibble()))
  }

  sim_root <- file.path(ctx$results_dir, "fits_nn_simulation")
  dir.create(sim_root, recursive = TRUE, showWarnings = FALSE)
  scenario_values <- ctx$nn_simulation_scenarios_use
  child_rows <- list()
  row_idx <- 0L

  for (sim_idx in seq_len(ctx$nn_simulation_n_use)) {
    scenario <- scenario_values[((sim_idx - 1L) %% length(scenario_values)) + 1L]
    sim <- simulate_nn_prior_counts(ctx$nn_simulation_seed_use + sim_idx, scenario)
    input_rds <- file.path(ctx$cache_dir, paste0("nn_simulation_", scenario, "_", sim_idx, ".rds"))
    saveRDS(sim$yi, input_rds)

    for (param_i in seq_len(nrow(parameter_spec_tbl))) {
      param_rr <- parameter_spec_tbl[param_i, , drop = FALSE]
      outdir <- build_fit_outdir(
        root_dir = sim_root,
        patient_id = paste0("sim_", sim_idx),
        path_components = list(
          parameter_label = as.character(param_rr$parameter_label),
          pm = ctx$pm_values_use[[1]],
          minobs = min(ctx$minobs_values_use)
        )
      )
      fit_res <- tryCatch(
        run_alfak_fit(
          patient_id = paste0("sim_", sim_idx),
          input_rds = input_rds,
          outdir = outdir,
          minobs = min(ctx$minobs_values_use),
          pm = ctx$pm_values_use[[1]],
          nboot = min(ctx$nboot_use, 15L),
          n0 = ctx$n0_use,
          nb = ctx$nb_use,
          benchmark_seed = ctx$nn_simulation_seed_use + sim_idx,
          parameter_label = as.character(param_rr$parameter_label),
          diploid_state = ctx$diploid_state,
          correct_efflux = ctx$correct_efflux_use,
          nn_prior = as.character(param_rr$nn_prior),
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
          force_refit = ctx$force_refit_use
        ),
        error = function(e) tibble::tibble(status = "error", error_message = conditionMessage(e))
      )
      boot_obj <- safe_read_rds(file.path(outdir, "bootstrap_res.Rds"))
      nn_mat <- boot_obj$nn_fitness
      if (is.null(nn_mat) || !ncol(nn_mat)) {
        next
      }
        for (kid in intersect(colnames(nn_mat), sim$truth$k[sim$truth$is_true_nn])) {
          pred <- as.numeric(nn_mat[, kid])
          target <- sim$truth$true_fitness[match(kid, sim$truth$k)]
          estimate <- if (any(is.finite(pred))) mean(pred, na.rm = TRUE) else NA_real_
          row_idx <- row_idx + 1L
          child_rows[[row_idx]] <- tibble::tibble(
            simulation_id = sim_idx,
            scenario = scenario,
            parameter_label = as.character(param_rr$parameter_label),
            nn_prior = as.character(param_rr$nn_prior),
            k = kid,
            true_fitness = target,
            estimated_fitness = estimate,
            estimation_error = estimate - target,
            status = if ("status" %in% names(fit_res)) as.character(fit_res$status[[1]]) else NA_character_
          )
        }
    }
  }

  child_tbl <- dplyr::bind_rows(child_rows)
  summary_tbl <- if (nrow(child_tbl)) {
    child_tbl %>%
      dplyr::filter(is.finite(estimation_error)) %>%
      dplyr::group_by(scenario, parameter_label, nn_prior) %>%
      dplyr::summarise(
        n_children = dplyr::n(),
        rmse = sqrt(mean(estimation_error^2)),
        mae = mean(abs(estimation_error)),
        signed_bias = mean(estimation_error),
        false_high_rate = mean(estimated_fitness > 0 & true_fitness <= 0, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble::tibble()
  }

  list(summary_tbl = summary_tbl, child_tbl = child_tbl)
}

build_benchmark_nn_diagnostics <- function(ctx,
                                           parameter_results_all_tbl,
                                           input_index_tbl,
                                           parameter_spec_tbl) {
  ident <- if (isTRUE(ctx$run_nn_identifiability_use)) {
    build_nn_identifiability_tables(parameter_results_all_tbl)
  } else {
    list(replicate_tbl = tibble::tibble(), summary_tbl = tibble::tibble())
  }
  save_table_bundle(ident$replicate_tbl, file.path(ctx$tables_dir, "nn_identifiability_by_replicate"))
  save_table_bundle(ident$summary_tbl, file.path(ctx$tables_dir, "nn_identifiability_summary"))

  stability <- if (isTRUE(ctx$run_nn_stability_use)) {
    build_nn_fitness_stability_tables(parameter_results_all_tbl)
  } else {
    list(child_tbl = tibble::tibble(), case_tbl = tibble::tibble(), weighted_empirical_tbl = tibble::tibble())
  }
  save_table_bundle(stability$child_tbl, file.path(ctx$tables_dir, "nn_fitness_stability_by_child"))
  save_table_bundle(stability$case_tbl, file.path(ctx$tables_dir, "nn_fitness_stability_by_case"))
  save_table_bundle(stability$weighted_empirical_tbl, file.path(ctx$tables_dir, "nn_fitness_weighted_empirical_comparison"))

  figure_artifact_tbl <- save_nn_diagnostic_figures(
    ctx = ctx,
    nn_identifiability_summary_tbl = ident$summary_tbl,
    nn_stability_case_tbl = stability$case_tbl,
    nn_stability_weighted_empirical_tbl = stability$weighted_empirical_tbl
  )
  save_table_bundle(figure_artifact_tbl, file.path(ctx$tables_dir, "nn_diagnostic_figure_artifacts"))

  holdout <- run_nn_holdout_diagnostics(
    ctx = ctx,
    input_index_tbl = input_index_tbl,
    parameter_spec_tbl = parameter_spec_tbl
  )
  save_table_bundle(holdout$summary_tbl, file.path(ctx$tables_dir, "nn_holdout_summary"))
  save_table_bundle(holdout$prediction_tbl, file.path(ctx$tables_dir, "nn_holdout_predictions"))

  simulation <- run_nn_simulation_diagnostics(
    ctx = ctx,
    parameter_spec_tbl = parameter_spec_tbl
  )
  save_table_bundle(simulation$summary_tbl, file.path(ctx$tables_dir, "nn_simulation_summary"))
  save_table_bundle(simulation$child_tbl, file.path(ctx$tables_dir, "nn_simulation_by_child"))

  list(
    identifiability_replicate_tbl = ident$replicate_tbl,
    identifiability_summary_tbl = ident$summary_tbl,
    stability_child_tbl = stability$child_tbl,
    stability_case_tbl = stability$case_tbl,
    stability_weighted_empirical_tbl = stability$weighted_empirical_tbl,
    figure_artifact_tbl = figure_artifact_tbl,
    holdout_summary_tbl = holdout$summary_tbl,
    holdout_prediction_tbl = holdout$prediction_tbl,
    simulation_summary_tbl = simulation$summary_tbl,
    simulation_child_tbl = simulation$child_tbl
  )
}
