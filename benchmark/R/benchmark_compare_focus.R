landscape_to_patient_df <- function(landscape_df, patient_id) {
  k_mat <- do.call(rbind, lapply(landscape_df$k, function(xx) {
    as.numeric(unlist(strsplit(xx, split = "\\.")))
  }))
  colnames(k_mat) <- as.character(seq_len(ncol(k_mat)))

  tibble::as_tibble(k_mat) %>%
    dplyr::mutate(
      karyotypes = as.character(landscape_df$k),
      fitness = as.numeric(landscape_df$mean),
      Patient = patient_id
    )
}

ensure_center_in_limits <- function(limits, center = 1) {
  c(min(limits[1], center), max(limits[2], center))
}

make_focus_umap_plot <- function(landscape_df,
                                 patient_id,
                                 parameter_label,
                                 benchmark_seed,
                                 diploid_state,
                                 scale_mode = c("relative", "absolute")) {
  scale_mode <- match.arg(scale_mode)
  k_df <- landscape_to_patient_df(landscape_df, patient_id)

  set.seed(benchmark_seed)
  umap_mat <- uwot::umap2(as.matrix(k_df[, as.character(seq_len(22))]), n_components = 2, min_dist = 0.9)
  umap_df <- as.data.frame(umap_mat)
  colnames(umap_df) <- c("UMAP1", "UMAP2")

  if (scale_mode == "relative") {
    dip_idx <- which(k_df$karyotypes == diploid_state)
    if (!length(dip_idx)) {
      stop("No diploid karyotype found for patient ", patient_id, " under ", parameter_label)
    }
    dip_mean <- mean(k_df$fitness[dip_idx], na.rm = TRUE)
    if (!is.finite(dip_mean) || dip_mean == 0) {
      stop("Diploid mean fitness is invalid for patient ", patient_id, " under ", parameter_label)
    }

    umap_df$plot_fitness <- k_df$fitness / dip_mean
    legend_name <- "Relative fitness\n(vs diploid)"
    color_center <- 1
    color_limits <- c(min(umap_df$plot_fitness, na.rm = TRUE), max(umap_df$plot_fitness, na.rm = TRUE))
    color_limits <- ensure_center_in_limits(color_limits, center = color_center)
    color_knots <- c(
      seq(color_limits[1], color_center, length.out = 6),
      seq(color_center, color_limits[2], length.out = 6)[-1]
    )
    color_breaks <- c(color_limits[1], color_center, color_limits[2])
  } else {
    umap_df$plot_fitness <- k_df$fitness
    legend_name <- "Absolute fitness"
    color_limits <- c(min(umap_df$plot_fitness, na.rm = TRUE), max(umap_df$plot_fitness, na.rm = TRUE))
    color_knots <- seq(color_limits[1], color_limits[2], length.out = 11)
    color_breaks <- c(color_limits[1], mean(color_limits), color_limits[2])
  }

  umap_df <- umap_df[order(umap_df$plot_fitness, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  if (diff(color_limits) == 0) {
    if (scale_mode == "relative") {
      color_limits <- c(color_center - 5e-07, color_center + 5e-07)
      color_knots <- c(
        seq(color_limits[1], color_center, length.out = 6),
        seq(color_center, color_limits[2], length.out = 6)[-1]
      )
      color_breaks <- c(color_limits[1], color_center, color_limits[2])
    } else {
      color_limits[2] <- color_limits[2] + 1e-06
      color_knots <- seq(color_limits[1], color_limits[2], length.out = 11)
      color_breaks <- c(color_limits[1], mean(color_limits), color_limits[2])
    }
  }

  color_values <- scales::rescale(color_knots, from = color_limits)

  ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = plot_fitness)) +
    ggplot2::geom_jitter(width = 0.1, height = 0.1, alpha = 1, size = 1, shape = 16) +
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
      title = paste0(patient_id, " UMAP: ", parameter_label),
      subtitle = if (scale_mode == "relative") "Relative to diploid fitness" else "Absolute landscape mean fitness",
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

build_focus_parameter_bundles <- function(results_tbl,
                                          focus_pid,
                                          focus_minobs,
                                          focus_pm,
                                          beneficial_move_levels,
                                          parameter_levels) {
  focus_idx <- results_tbl %>%
    dplyr::filter(
      status == "ok",
      patient_id == focus_pid,
      minobs == focus_minobs,
      dplyr::near(pm, focus_pm),
      parameter_label %in% parameter_levels
    ) %>%
    dplyr::distinct(parameter_label, .keep_all = TRUE) %>%
    dplyr::arrange(match(parameter_label, parameter_levels))

  if (!nrow(focus_idx)) {
    return(list())
  }

  bundle_list <- lapply(focus_idx$outdir, read_fit_bundle, beneficial_move_levels = beneficial_move_levels)
  names(bundle_list) <- focus_idx$parameter_label
  bundle_list
}

landscape_long_from_bundles <- function(bundle_list) {
  if (!length(bundle_list)) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(lapply(names(bundle_list), function(parameter_label_name) {
    landscape_df <- bundle_list[[parameter_label_name]]$landscape
    if (is.null(landscape_df) || !nrow(landscape_df)) {
      return(tibble::tibble())
    }

    tibble::tibble(
      parameter_label = parameter_label_name,
      k = as.character(landscape_df$k),
      mean = as.numeric(landscape_df$mean),
      median = as.numeric(landscape_df$median),
      sd = as.numeric(landscape_df$sd),
      fq = as.logical(landscape_df$fq),
      nn = as.logical(landscape_df$nn)
    ) %>%
      dplyr::mutate(
        state_group = dplyr::case_when(
          fq %in% TRUE ~ "fq",
          nn %in% TRUE ~ "nn",
          TRUE ~ "other"
        )
      )
  }))
}

summarize_landscape_by_parameter <- function(landscape_long) {
  if (!nrow(landscape_long)) {
    return(tibble::tibble())
  }

  landscape_long %>%
    dplyr::group_by(parameter_label) %>%
    dplyr::summarise(
      n_karyotypes = dplyr::n(),
      n_fq = sum(fq %in% TRUE, na.rm = TRUE),
      n_nn = sum(nn %in% TRUE, na.rm = TRUE),
      mean_landscape_mean = mean(mean, na.rm = TRUE),
      median_landscape_mean = median(mean, na.rm = TRUE),
      mean_landscape_sd = mean(sd, na.rm = TRUE),
      min_landscape_mean = min(mean, na.rm = TRUE),
      max_landscape_mean = max(mean, na.rm = TRUE),
      negative_mean_n = sum(mean < 0, na.rm = TRUE),
      negative_median_n = sum(median < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_focus_landscape_variation <- function(landscape_long, parameter_levels, top_n = 15L) {
  if (!nrow(landscape_long) || !length(parameter_levels)) {
    return(tibble::tibble())
  }

  wide_tbl <- landscape_long %>%
    dplyr::select(parameter_label, k, mean) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = parameter_label, values_from = mean)

  for (parameter_name in parameter_levels) {
    if (!parameter_name %in% names(wide_tbl)) {
      wide_tbl[[parameter_name]] <- NA_real_
    }
  }
  wide_tbl <- wide_tbl[, c("k", parameter_levels), drop = FALSE]
  mean_mat <- as.matrix(wide_tbl[, parameter_levels, drop = FALSE])

  state_tbl <- landscape_long %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(
      state_group = dplyr::case_when(
        any(fq %in% TRUE, na.rm = TRUE) ~ "fq",
        any(nn %in% TRUE, na.rm = TRUE) ~ "nn",
        TRUE ~ "other"
      ),
      .groups = "drop"
    )

  wide_tbl %>%
    dplyr::left_join(state_tbl, by = "k") %>%
    dplyr::mutate(
      n_parameter_label = rowSums(is.finite(mean_mat)),
      mean_range = apply(mean_mat, 1, function(x) {
        if (sum(is.finite(x)) >= 2L) diff(range(x, na.rm = TRUE)) else NA_real_
      }),
      mean_sd_across = apply(mean_mat, 1, function(x) {
        if (sum(is.finite(x)) >= 2L) stats::sd(x, na.rm = TRUE) else NA_real_
      }),
      min_parameter_label = apply(mean_mat, 1, function(x) {
        if (!sum(is.finite(x))) {
          return(NA_character_)
        }
        parameter_levels[which.min(replace(x, !is.finite(x), Inf))][1]
      }),
      max_parameter_label = apply(mean_mat, 1, function(x) {
        if (!sum(is.finite(x))) {
          return(NA_character_)
        }
        parameter_levels[which.max(replace(x, !is.finite(x), -Inf))][1]
      })
    ) %>%
    dplyr::select(k, state_group, n_parameter_label, dplyr::all_of(parameter_levels), mean_range, mean_sd_across, min_parameter_label, max_parameter_label) %>%
    dplyr::arrange(dplyr::desc(mean_range), dplyr::desc(mean_sd_across), k) %>%
    dplyr::slice_head(n = top_n)
}

build_focus_parity_tbl <- function(bundle_list, parameter_levels, value_col = "mean") {
  if (length(parameter_levels) < 2L) {
    return(tibble::tibble())
  }

  pair_defs <- combn(parameter_levels, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(pair_defs, function(pair) {
    lhs_name <- pair[[1]]
    rhs_name <- pair[[2]]
    lhs_landscape <- bundle_list[[lhs_name]]$landscape
    rhs_landscape <- bundle_list[[rhs_name]]$landscape

    if (is.null(lhs_landscape) || is.null(rhs_landscape)) {
      return(tibble::tibble())
    }

    lhs_tbl <- lhs_landscape %>%
      dplyr::transmute(
        k,
        lhs_value = .data[[value_col]],
        lhs_fq = as.logical(fq),
        lhs_nn = as.logical(nn)
      )
    rhs_tbl <- rhs_landscape %>%
      dplyr::transmute(
        k,
        rhs_value = .data[[value_col]],
        rhs_fq = as.logical(fq),
        rhs_nn = as.logical(nn)
      )

    dplyr::inner_join(lhs_tbl, rhs_tbl, by = "k") %>%
      dplyr::mutate(
        lhs_label = lhs_name,
        rhs_label = rhs_name,
        comparison = paste(lhs_name, rhs_name, sep = " vs "),
        state_group = dplyr::case_when(
          lhs_fq | rhs_fq ~ "fq",
          lhs_nn | rhs_nn ~ "nn",
          TRUE ~ "other"
        )
      ) %>%
      dplyr::select(comparison, lhs_label, rhs_label, k, state_group, lhs_value, rhs_value)
  }))
}

beneficial_long_from_profiles <- function(beneficial_profiles, beneficial_move_levels, parameter_levels = names(beneficial_profiles)) {
  if (!length(beneficial_profiles)) {
    return(tibble::tibble())
  }

  out <- dplyr::bind_rows(lapply(parameter_levels, function(parameter_label_name) {
    prof <- beneficial_profiles[[parameter_label_name]]
    if (is.null(prof)) {
      return(tibble::tibble())
    }

    tibble::tibble(
      parameter_label = parameter_label_name,
      move = beneficial_move_levels,
      proportion = as.numeric(prof$proportion[beneficial_move_levels]),
      beneficial_n = as.integer(prof$beneficial_n[beneficial_move_levels]),
      valid_n = as.integer(prof$valid_n[beneficial_move_levels])
    )
  }))

  if (!nrow(out) || !"move" %in% names(out)) {
    return(tibble::tibble())
  }

  out %>%
    dplyr::mutate(
      chromosome = suppressWarnings(as.integer(sub("[+-]$", "", move))),
      direction = ifelse(grepl("\\+$", move), "gain", "loss"),
      move = factor(move, levels = beneficial_move_levels),
      parameter_label = factor(parameter_label, levels = parameter_levels)
    )
}

summarize_beneficial_by_parameter <- function(beneficial_long) {
  if (!nrow(beneficial_long)) {
    return(tibble::tibble())
  }

  beneficial_long %>%
    dplyr::group_by(parameter_label) %>%
    dplyr::summarise(
      n_move_types = sum(!is.na(proportion)),
      mean_proportion = mean(proportion, na.rm = TRUE),
      median_proportion = median(proportion, na.rm = TRUE),
      total_beneficial_n = sum(beneficial_n, na.rm = TRUE),
      total_valid_n = sum(valid_n, na.rm = TRUE),
      weighted_proportion = ifelse(sum(valid_n, na.rm = TRUE) > 0, sum(beneficial_n, na.rm = TRUE) / sum(valid_n, na.rm = TRUE), NA_real_),
      .groups = "drop"
    )
}

build_focus_beneficial_shift_tbl <- function(beneficial_long, parameter_levels, top_n = 15L) {
  if (!nrow(beneficial_long) || length(parameter_levels) < 2L) {
    return(tibble::tibble())
  }

  pair_defs <- combn(parameter_levels, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(pair_defs, function(pair) {
    lhs_name <- pair[[1]]
    rhs_name <- pair[[2]]
    wide_tbl <- beneficial_long %>%
      dplyr::filter(parameter_label %in% c(lhs_name, rhs_name)) %>%
      dplyr::select(parameter_label, move, chromosome, direction, proportion, valid_n) %>%
      tidyr::pivot_wider(names_from = parameter_label, values_from = c(proportion, valid_n), names_sep = "__")

    lhs_prop_col <- paste0("proportion__", lhs_name)
    rhs_prop_col <- paste0("proportion__", rhs_name)
    lhs_valid_col <- paste0("valid_n__", lhs_name)
    rhs_valid_col <- paste0("valid_n__", rhs_name)
    missing_cols <- setdiff(c(lhs_prop_col, rhs_prop_col, lhs_valid_col, rhs_valid_col), names(wide_tbl))
    for (col in missing_cols) {
      wide_tbl[[col]] <- NA_real_
    }

    wide_tbl %>%
      dplyr::transmute(
        lhs_label = lhs_name,
        rhs_label = rhs_name,
        move,
        chromosome,
        direction,
        lhs_prop = .data[[lhs_prop_col]],
        rhs_prop = .data[[rhs_prop_col]],
        diff = rhs_prop - lhs_prop,
        lhs_valid_n = .data[[lhs_valid_col]],
        rhs_valid_n = .data[[rhs_valid_col]]
      ) %>%
      dplyr::arrange(dplyr::desc(abs(diff)), move) %>%
      dplyr::slice_head(n = top_n)
  }))
}
