beneficial_cnv_profile <- function(landscape_df, beneficial_move_levels) {
  fq_flag <- landscape_df$fq
  fq_keep <- if (is.logical(fq_flag)) {
    fq_flag
  } else {
    as.character(fq_flag) %in% c("TRUE", "T", "true", "1")
  }

  k_vec <- as.character(landscape_df$k)
  prop <- rep(NA_real_, length(beneficial_move_levels))
  valid_n <- integer(length(beneficial_move_levels))
  beneficial_n <- integer(length(beneficial_move_levels))
  names(prop) <- beneficial_move_levels
  names(valid_n) <- beneficial_move_levels
  names(beneficial_n) <- beneficial_move_levels

  if (!length(k_vec)) {
    return(list(
      proportion = prop,
      beneficial_n = beneficial_n,
      valid_n = valid_n
    ))
  }

  k_parts <- strsplit(k_vec, ".", fixed = TRUE)
  part_lengths <- lengths(k_parts)
  if (any(part_lengths != 22L)) {
    bad_idx <- which(part_lengths != 22L)[1]
    stop("Found karyotype with != 22 entries in beneficial_cnv_profile(): ", k_vec[bad_idx])
  }

  k_mat <- matrix(as.integer(unlist(k_parts, use.names = FALSE)), ncol = 22, byrow = TRUE)
  fitness_vec <- as.numeric(landscape_df$mean)
  starts <- which(fq_keep %in% TRUE)
  if (!length(starts)) {
    return(list(
      proportion = prop,
      beneficial_n = beneficial_n,
      valid_n = valid_n
    ))
  }

  edge_rows <- vector("list", length(starts))
  edge_idx <- 0L
  for (from in starts) {
    v <- k_mat[from, ]
    neigh_mat <- matrix(rep(v, each = length(beneficial_move_levels)), nrow = length(beneficial_move_levels), ncol = 22)
    for (chr in seq_len(22)) {
      plus_row <- 2L * chr - 1L
      minus_row <- 2L * chr
      neigh_mat[plus_row, chr] <- neigh_mat[plus_row, chr] + 1L
      neigh_mat[minus_row, chr] <- neigh_mat[minus_row, chr] - 1L
    }

    neigh_str <- do.call(paste, c(as.data.frame(neigh_mat), sep = "."))
    to <- match(neigh_str, k_vec)
    valid_cn <- apply(neigh_mat, 1, function(row) all(row >= 1L))
    valid <- valid_cn & !is.na(to)
    if (!any(valid)) {
      next
    }

    to_valid <- to[valid]
    move_valid <- beneficial_move_levels[valid]
    delta <- fitness_vec[to_valid] - fitness_vec[from]
    ok <- is.finite(delta)
    benefit <- rep(NA, length(delta))
    benefit[ok] <- delta[ok] > 0

    edge_idx <- edge_idx + 1L
    edge_rows[[edge_idx]] <- data.frame(
      move = move_valid,
      ok = ok,
      benefit = benefit,
      stringsAsFactors = FALSE
    )
  }

  if (edge_idx > 0L) {
    edges_df <- dplyr::bind_rows(edge_rows[seq_len(edge_idx)])
    for (mv in beneficial_move_levels) {
      ii <- which(edges_df$move == mv)
      if (!length(ii)) {
        next
      }
      sub <- edges_df[ii, , drop = FALSE]
      ok <- sub$ok %in% TRUE
      valid_n[mv] <- sum(ok)
      if (!valid_n[mv]) {
        next
      }
      beneficial_n[mv] <- sum(sub$benefit[ok] %in% TRUE, na.rm = TRUE)
      prop[mv] <- beneficial_n[mv] / valid_n[mv]
    }
  }

  list(
    proportion = prop,
    beneficial_n = beneficial_n,
    valid_n = valid_n
  )
}

read_fit_bundle <- function(fit_dir, beneficial_move_levels) {
  landscape_obj <- safe_read_rds(file.path(fit_dir, "landscape.Rds"))
  list(
    bootstrap = safe_read_rds(file.path(fit_dir, "bootstrap_res.Rds")),
    landscape = landscape_obj,
    posterior = safe_read_rds(file.path(fit_dir, "landscape_posterior_samples.Rds")),
    xval = safe_read_rds(file.path(fit_dir, "xval.Rds")),
    beneficial = if (!is.null(landscape_obj)) beneficial_cnv_profile(landscape_obj, beneficial_move_levels = beneficial_move_levels) else NULL
  )
}

rank_parameter_fit_tbl <- function(fit_tbl, parameter_levels = NULL) {
  if (is.null(fit_tbl) || !nrow(fit_tbl)) {
    return(tibble::tibble())
  }

  if (is.null(parameter_levels) || !length(parameter_levels)) {
    parameter_levels <- unique(as.character(fit_tbl$parameter_label))
  }

  fit_tbl %>%
    dplyr::mutate(
      status = as.character(status),
      has_positive_xval = is.finite(xval_r2) & xval_r2 > 0,
      rank_xval_r2 = dplyr::coalesce(xval_r2, -Inf),
      rank_xval_cor = dplyr::coalesce(xval_cor, -Inf),
      rank_xval_rmse = dplyr::coalesce(xval_rmse, Inf),
      rank_xval_mae = dplyr::coalesce(xval_mae, Inf)
    ) %>%
    dplyr::arrange(
      factor(parameter_label, levels = parameter_levels),
      factor(patient_id, levels = sort_pid_levels(patient_id)),
      dplyr::desc(status == "ok"),
      dplyr::desc(has_positive_xval),
      dplyr::desc(rank_xval_r2),
      rank_xval_rmse,
      rank_xval_mae,
      dplyr::desc(rank_xval_cor),
      minobs,
      pm
    )
}

select_best_parameter_fit_tbl <- function(fit_tbl, parameter_levels = NULL) {
  ranked_tbl <- rank_parameter_fit_tbl(fit_tbl, parameter_levels = parameter_levels)
  if (!nrow(ranked_tbl)) {
    return(tibble::tibble())
  }

  ranked_tbl %>%
    dplyr::filter(status == "ok") %>%
    dplyr::group_by(parameter_label, patient_id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-has_positive_xval, -rank_xval_r2, -rank_xval_cor, -rank_xval_rmse, -rank_xval_mae)
}

select_positive_xval_best_parameter_fit_tbl <- function(fit_tbl, parameter_levels = NULL) {
  ranked_tbl <- rank_parameter_fit_tbl(fit_tbl, parameter_levels = parameter_levels)
  if (!nrow(ranked_tbl)) {
    return(tibble::tibble())
  }

  ranked_tbl %>%
    dplyr::filter(status == "ok", is.finite(xval_r2), xval_r2 > 0) %>%
    dplyr::group_by(parameter_label, patient_id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-has_positive_xval, -rank_xval_r2, -rank_xval_cor, -rank_xval_rmse, -rank_xval_mae)
}

select_common_positive_xval_patients <- function(selected_fit_tbl, parameter_levels) {
  if (is.null(selected_fit_tbl) || !nrow(selected_fit_tbl) || !length(parameter_levels)) {
    return(character(0))
  }

  selected_fit_tbl %>%
    dplyr::filter(parameter_label %in% parameter_levels) %>%
    dplyr::distinct(patient_id, parameter_label) %>%
    dplyr::count(patient_id, name = "n_parameter_label") %>%
    dplyr::filter(n_parameter_label == length(parameter_levels)) %>%
    dplyr::pull(patient_id) %>%
    sort_pid_levels()
}

parse_karyotype_matrix_from_strings <- function(k_vec) {
  k_vec <- as.character(k_vec)
  if (!length(k_vec)) {
    return(matrix(numeric(0), nrow = 0, ncol = 22))
  }

  k_parts <- strsplit(k_vec, ".", fixed = TRUE)
  part_lengths <- lengths(k_parts)
  if (any(part_lengths != 22L)) {
    bad_idx <- which(part_lengths != 22L)[1]
    stop("Found karyotype with != 22 entries: ", k_vec[bad_idx])
  }
  matrix(as.integer(unlist(k_parts, use.names = FALSE)), ncol = 22, byrow = TRUE)
}

build_selected_landscape_long_tbl <- function(selected_fit_tbl,
                                              beneficial_move_levels,
                                              parameter_levels = NULL) {
  if (is.null(selected_fit_tbl) || !nrow(selected_fit_tbl)) {
    return(tibble::tibble())
  }

  if (is.null(parameter_levels) || !length(parameter_levels)) {
    parameter_levels <- unique(as.character(selected_fit_tbl$parameter_label))
  }

  dplyr::bind_rows(lapply(seq_len(nrow(selected_fit_tbl)), function(i) {
    rr <- selected_fit_tbl[i, , drop = FALSE]
    bundle <- read_fit_bundle(rr$outdir, beneficial_move_levels = beneficial_move_levels)
    landscape_df <- bundle$landscape
    if (is.null(landscape_df) || !nrow(landscape_df)) {
      return(tibble::tibble())
    }

    tibble::tibble(
      patient_id = as.character(rr$patient_id),
      parameter_label = as.character(rr$parameter_label),
      minobs = as.integer(rr$minobs),
      pm = as.numeric(rr$pm),
      k = as.character(landscape_df$k),
      mean = as.numeric(landscape_df$mean),
      median = as.numeric(landscape_df$median),
      sd = as.numeric(landscape_df$sd),
      fq = as_benchmark_logical_flag(landscape_df$fq),
      nn = as_benchmark_logical_flag(landscape_df$nn)
    )
  })) %>%
    dplyr::mutate(
      parameter_label = factor(parameter_label, levels = parameter_levels),
      patient_id = factor(patient_id, levels = sort_pid_levels(patient_id))
    )
}

summarize_selected_landscape_sample_tbl <- function(landscape_long) {
  if (is.null(landscape_long) || !nrow(landscape_long)) {
    return(tibble::tibble())
  }

  landscape_long %>%
    dplyr::group_by(parameter_label, patient_id) %>%
    dplyr::summarise(
      n_karyotypes = dplyr::n(),
      landscape_mean_mean = mean(mean, na.rm = TRUE),
      landscape_mean_sd = stats::sd(mean, na.rm = TRUE),
      landscape_mean_median = stats::median(mean, na.rm = TRUE),
      .groups = "drop"
    )
}

reconstruct_xval_detail_tbl <- function(fq_boot,
                                        benchmark_seed = NULL,
                                        seed_offset = 0L,
                                        krig_bootstrap_mode = c("marginal", "joint")) {
  krig_bootstrap_mode <- match.arg(krig_bootstrap_mode)
  if (!is.null(benchmark_seed) && is.finite(benchmark_seed)) {
    set.seed(as.integer(benchmark_seed) + as.integer(seed_offset))
  }

  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness)

  valid_fq_str <- if (is.null(fq_str)) character(0) else fq_str
  valid_nn_str <- if (is.null(nn_str)) character(0) else nn_str
  combined_strs <- c(valid_fq_str, valid_nn_str)
  if (!length(combined_strs) || ncol(fboot) == 0 || nrow(fboot) == 0) {
    return(tibble::tibble(k = character(0), state_class = character(0), validation = numeric(0), estimate = numeric(0)))
  }
  if (!length(valid_fq_str)) {
    return(tibble::tibble(k = character(0), state_class = character(0), validation = numeric(0), estimate = numeric(0)))
  }

  ktrain <- unname(alfakR:::parse_karyotype_ids(combined_strs))
  ids <- unlist(lapply(seq_along(valid_fq_str), function(i_xval) {
    ki_neighbours_matrix <- alfakR:::gen_all_neighbours(valid_fq_str[i_xval])
    ki_neighbours_str <- character(0)
    if (nrow(ki_neighbours_matrix) > 0) {
      ki_neighbours_str <- as.character(apply(ki_neighbours_matrix, 1, paste, collapse = "."))
    }
    ki <- c(valid_fq_str[i_xval], ki_neighbours_str)
    idi <- rep(i_xval, length(ki))
    names(idi) <- ki
    idi
  }))
  ids <- ids[!duplicated(names(ids))]
  uids <- unique(ids)
  ktrain_map <- setNames(seq_len(nrow(ktrain)), combined_strs)

  tmp_list <- lapply(uids, function(id_fold) {
    fi <- if (krig_bootstrap_mode == "joint") {
      fboot[sample(seq_len(nrow(fboot)), 1), ]
    } else {
      fboot_shuffled <- apply(fboot, 2, sample)
      fboot_shuffled[1, ]
    }

    train_k_names <- names(ids)[ids != id_fold]
    test_k_names <- names(ids)[ids == id_fold]
    train_k_names_valid <- train_k_names[train_k_names %in% names(ktrain_map)]
    test_k_names_valid <- test_k_names[test_k_names %in% names(ktrain_map)]
    if (!length(train_k_names_valid) || !length(test_k_names_valid)) {
      return(tibble::tibble(k = character(0), state_class = character(0), validation = numeric(0), estimate = numeric(0)))
    }

    train_k <- ktrain[ktrain_map[train_k_names_valid], , drop = FALSE]
    train_f <- fi[train_k_names_valid]
    test_k <- ktrain[ktrain_map[test_k_names_valid], , drop = FALSE]
    test_f <- fi[test_k_names_valid]
    test_state_class <- ifelse(test_k_names_valid %in% valid_fq_str, "fq", "nn")

    valid_train_points <- !is.na(train_f)
    train_k <- train_k[valid_train_points, , drop = FALSE]
    train_f <- train_f[valid_train_points]

    if (nrow(train_k) < 2 || nrow(unique(train_k)) < 2 || length(unique(train_f)) < 1) {
      return(tibble::tibble(
        k = test_k_names_valid,
        state_class = test_state_class,
        validation = as.numeric(test_f),
        estimate = rep(NA_real_, length(test_f))
      ))
    }

    fit <- suppressWarnings(fields::Krig(
      train_k,
      train_f,
      cov.function = "stationary.cov",
      cov.args = alfakR:::krig_covariance_args(),
      nstep.cv = alfakR:::ALFAK_KRIG_NSTEP_CV,
      give.warnings = TRUE
    ))
    pred_f <- suppressWarnings(stats::predict(fit, test_k))
    tibble::tibble(
      k = test_k_names_valid,
      state_class = test_state_class,
      validation = as.numeric(test_f),
      estimate = as.numeric(pred_f)
    )
  })

  tmp <- do.call(rbind, tmp_list)
  tmp <- as.data.frame(tmp)
  tmp <- tmp[stats::complete.cases(tmp[, c("validation", "estimate"), drop = FALSE]), , drop = FALSE]
  if (!nrow(tmp)) {
    return(tibble::tibble(k = character(0), state_class = character(0), validation = numeric(0), estimate = numeric(0)))
  }

  tibble::tibble(
    k = as.character(tmp$k),
    state_class = as.character(tmp$state_class),
    validation = as.numeric(tmp$validation),
    estimate = as.numeric(tmp$estimate)
  ) %>%
    dplyr::filter(is.finite(validation), is.finite(estimate))
}

build_selected_xval_scatter_tbl <- function(selected_fit_tbl,
                                            beneficial_move_levels,
                                            benchmark_seed,
                                            parameter_levels = NULL) {
  if (is.null(selected_fit_tbl) || !nrow(selected_fit_tbl)) {
    return(tibble::tibble())
  }

  if (is.null(parameter_levels) || !length(parameter_levels)) {
    parameter_levels <- unique(as.character(selected_fit_tbl$parameter_label))
  }

  dplyr::bind_rows(lapply(seq_len(nrow(selected_fit_tbl)), function(i) {
    rr <- selected_fit_tbl[i, , drop = FALSE]
    bundle <- read_fit_bundle(rr$outdir, beneficial_move_levels = beneficial_move_levels)
    if (is.null(bundle$bootstrap)) {
      return(tibble::tibble())
    }

    seed_offset <- i * 1000L
    xval_tbl <- reconstruct_xval_detail_tbl(
      fq_boot = bundle$bootstrap,
      benchmark_seed = benchmark_seed,
      seed_offset = seed_offset,
      krig_bootstrap_mode = "marginal"
    )
    if (!nrow(xval_tbl)) {
      return(tibble::tibble())
    }

    xval_tbl %>%
      dplyr::filter(state_class == "fq") %>%
      dplyr::mutate(
        patient_id = as.character(rr$patient_id),
        parameter_label = as.character(rr$parameter_label),
        minobs = as.integer(rr$minobs),
        pm = as.numeric(rr$pm),
        xval_r2 = as.numeric(rr$xval_r2)
      )
  })) %>%
    dplyr::mutate(
      parameter_label = factor(parameter_label, levels = parameter_levels),
      patient_id = factor(patient_id, levels = sort_pid_levels(patient_id))
    )
}

build_selected_umap_long_tbl <- function(landscape_long,
                                         benchmark_seed,
                                         diploid_state) {
  if (is.null(landscape_long) || !nrow(landscape_long)) {
    return(list(
      landscape_with_umap = tibble::tibble(),
      reference_umap_tbl = tibble::tibble(),
      global_diploid_mean_tbl = tibble::tibble()
    ))
  }

  parameter_levels <- unique(as.character(landscape_long$parameter_label))
  landscape_with_umap_tbl <- dplyr::bind_rows(lapply(seq_along(parameter_levels), function(i) {
    parameter_label_name <- parameter_levels[[i]]
    sub_long <- landscape_long %>%
      dplyr::mutate(parameter_label = as.character(parameter_label)) %>%
      dplyr::filter(parameter_label == parameter_label_name)
    if (!nrow(sub_long)) {
      return(tibble::tibble())
    }

    base_tbl <- sub_long %>%
      dplyr::select(patient_id, parameter_label, k) %>%
      dplyr::distinct()
    k_mat <- parse_karyotype_matrix_from_strings(base_tbl$k)
    set.seed(as.integer(benchmark_seed) + i * 1000L)
    umap_mat <- uwot::umap2(k_mat, n_components = 2, min_dist = 0.9)
    umap_tbl <- tibble::tibble(
      patient_id = base_tbl$patient_id,
      parameter_label = base_tbl$parameter_label,
      k = base_tbl$k,
      UMAP1 = umap_mat[, 1],
      UMAP2 = umap_mat[, 2]
    )

    sub_long %>%
      dplyr::left_join(umap_tbl, by = c("patient_id", "parameter_label", "k"))
  })) %>%
    dplyr::mutate(
      parameter_label = factor(parameter_label, levels = parameter_levels),
      patient_id = factor(patient_id, levels = sort_pid_levels(patient_id))
    )

  diploid_tbl <- landscape_with_umap_tbl %>%
    dplyr::filter(k == diploid_state) %>%
    dplyr::group_by(patient_id, parameter_label) %>%
    dplyr::summarise(diploid_mean = mean(mean, na.rm = TRUE), .groups = "drop")
  global_diploid_mean_tbl <- diploid_tbl %>%
    dplyr::group_by(parameter_label) %>%
    dplyr::summarise(global_diploid_mean = mean(diploid_mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      global_diploid_mean = dplyr::if_else(
        is.finite(global_diploid_mean) & global_diploid_mean != 0,
        global_diploid_mean,
        NA_real_
      )
    )

  reference_umap_tbl <- landscape_with_umap_tbl %>%
    dplyr::left_join(global_diploid_mean_tbl, by = "parameter_label") %>%
    dplyr::mutate(reference_value = mean / global_diploid_mean)

  list(
    landscape_with_umap = landscape_with_umap_tbl,
    reference_umap_tbl = reference_umap_tbl,
    global_diploid_mean_tbl = global_diploid_mean_tbl
  )
}

load_selected_landscapes_by_parameter <- function(selected_fit_tbl,
                                                  beneficial_move_levels,
                                                  parameter_levels = NULL) {
  if (is.null(selected_fit_tbl) || !nrow(selected_fit_tbl)) {
    return(list())
  }

  if (is.null(parameter_levels) || !length(parameter_levels)) {
    parameter_levels <- unique(as.character(selected_fit_tbl$parameter_label))
  }

  out <- setNames(vector("list", length(parameter_levels)), parameter_levels)
  for (parameter_label_name in parameter_levels) {
    sub_tbl <- selected_fit_tbl %>%
      dplyr::filter(parameter_label == parameter_label_name) %>%
      dplyr::arrange(factor(patient_id, levels = sort_pid_levels(patient_id)))

    if (!nrow(sub_tbl)) {
      out[[parameter_label_name]] <- list()
      next
    }

    landscape_list <- lapply(sub_tbl$outdir, function(outdir) {
      read_fit_bundle(outdir, beneficial_move_levels = beneficial_move_levels)$landscape
    })
    names(landscape_list) <- sub_tbl$patient_id
    out[[parameter_label_name]] <- landscape_list
  }

  out
}

save_beneficial_heatmap_png <- function(beneficial_mat, png_path) {
  if (!nrow(beneficial_mat) || !ncol(beneficial_mat)) {
    return(invisible(png_path))
  }

  panel_d_mat <- beneficial_mat
  move_tbl <- tibble::tibble(move = colnames(panel_d_mat)) %>%
    dplyr::mutate(
      chr = suppressWarnings(as.integer(gsub("[+-]$", "", move))),
      sign = sub("^[0-9]+", "", move),
      sign_rank = ifelse(sign == "+", 0L, 1L)
    ) %>%
    dplyr::arrange(chr, sign_rank, move)
  ordered_moves <- move_tbl$move
  panel_d_mat <- panel_d_mat[, ordered_moves, drop = FALSE]
  panel_d_mat[!is.finite(panel_d_mat)] <- 0
  col_fun_d <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "#EEEEEE", "red"))

  ht_d <- ComplexHeatmap::Heatmap(
    panel_d_mat,
    name = "Beneficial\nkaryotypes\nproportion",
    col = col_fun_d,
    cluster_rows = nrow(panel_d_mat) > 1,
    cluster_columns = FALSE,
    na_col = "#BDBDBD",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 7),
    rect_gp = grid::gpar(col = "white", lwd = 1)
  )

  plot_height <- max(4.5, 0.28 * nrow(panel_d_mat) + 1.5)
  grDevices::png(png_path, width = 10 * 150, height = plot_height * 150, res = 150)
  ComplexHeatmap::draw(ht_d, heatmap_legend_side = "right")
  grDevices::dev.off()
  invisible(png_path)
}

build_parameter_beneficial_artifacts <- function(selected_landscapes_by_parameter,
                                                 beneficial_move_levels,
                                                 parameter_levels,
                                                 tables_dir,
                                                 figures_dir) {
  if (!length(parameter_levels)) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(lapply(parameter_levels, function(parameter_label_name) {
    landscape_list <- selected_landscapes_by_parameter[[parameter_label_name]]
    if (is.null(landscape_list) || !length(landscape_list)) {
      return(tibble::tibble(
        parameter_label = parameter_label_name,
        n_patients = 0L,
        proportion_matrix_path = NA_character_,
        valid_n_matrix_path = NA_character_,
        heatmap_png_path = NA_character_
      ))
    }

    beneficial_profiles <- lapply(landscape_list, beneficial_cnv_profile, beneficial_move_levels = beneficial_move_levels)
    beneficial_list <- lapply(beneficial_profiles, function(x) x$proportion)
    beneficial_valid_n_list <- lapply(beneficial_profiles, function(x) x$valid_n)

    beneficial_mat <- do.call(rbind, beneficial_list)
    beneficial_mat <- beneficial_mat[sort_pid_levels(rownames(beneficial_mat)), , drop = FALSE]
    beneficial_valid_n_mat <- do.call(rbind, beneficial_valid_n_list)
    beneficial_valid_n_mat <- beneficial_valid_n_mat[sort_pid_levels(rownames(beneficial_valid_n_mat)), , drop = FALSE]
    move_tbl <- tibble::tibble(move = colnames(beneficial_mat)) %>%
      dplyr::mutate(
        chr = suppressWarnings(as.integer(gsub("[+-]$", "", move))),
        sign = sub("^[0-9]+", "", move),
        sign_rank = ifelse(sign == "+", 0L, 1L)
      ) %>%
      dplyr::arrange(chr, sign_rank, move)
    ordered_moves <- move_tbl$move
    beneficial_mat <- beneficial_mat[, ordered_moves, drop = FALSE]
    beneficial_valid_n_mat <- beneficial_valid_n_mat[, ordered_moves, drop = FALSE]

    proportion_stem <- file.path(tables_dir, paste0("parameter_", parameter_label_name, "_beneficial_proportion_matrix"))
    valid_n_stem <- file.path(tables_dir, paste0("parameter_", parameter_label_name, "_beneficial_valid_n_matrix"))
    heatmap_png_path <- file.path(figures_dir, paste0("beneficial_heatmap_", parameter_label_name, ".png"))

    save_table_bundle(
      tibble::rownames_to_column(as.data.frame(beneficial_mat), var = "Patient"),
      proportion_stem
    )
    save_table_bundle(
      tibble::rownames_to_column(as.data.frame(beneficial_valid_n_mat), var = "Patient"),
      valid_n_stem
    )
    save_beneficial_heatmap_png(beneficial_mat, heatmap_png_path)

    tibble::tibble(
      parameter_label = parameter_label_name,
      n_patients = nrow(beneficial_mat),
      proportion_matrix_path = paste0(proportion_stem, ".tsv"),
      valid_n_matrix_path = paste0(valid_n_stem, ".tsv"),
      heatmap_png_path = heatmap_png_path
    )
  }))
}

flatten_pair_values <- function(lhs, rhs) {
  if (is.null(lhs) || is.null(rhs)) {
    return(list(lhs = numeric(), rhs = numeric(), lhs_total = 0L, rhs_total = 0L, n_used = 0L))
  }

  if (!is.null(dim(lhs)) && !is.null(dim(rhs))) {
    lhs_mat <- as.matrix(lhs)
    rhs_mat <- as.matrix(rhs)

    if (!is.null(colnames(lhs_mat)) && !is.null(colnames(rhs_mat))) {
      common_cols <- intersect(colnames(lhs_mat), colnames(rhs_mat))
      if (length(common_cols)) {
        lhs_mat <- lhs_mat[, common_cols, drop = FALSE]
        rhs_mat <- rhs_mat[, common_cols, drop = FALSE]
      }
    }

    n_rows <- min(nrow(lhs_mat), nrow(rhs_mat))
    if (n_rows < 1L) {
      return(list(lhs = numeric(), rhs = numeric(), lhs_total = length(lhs_mat), rhs_total = length(rhs_mat), n_used = 0L))
    }
    lhs_vec <- as.numeric(lhs_mat[seq_len(n_rows), , drop = FALSE])
    rhs_vec <- as.numeric(rhs_mat[seq_len(n_rows), , drop = FALSE])
  } else {
    lhs_vec <- as.numeric(lhs)
    rhs_vec <- as.numeric(rhs)
  }

  n_used <- min(length(lhs_vec), length(rhs_vec))
  if (n_used < 1L) {
    return(list(lhs = numeric(), rhs = numeric(), lhs_total = length(lhs_vec), rhs_total = length(rhs_vec), n_used = 0L))
  }

  list(
    lhs = lhs_vec[seq_len(n_used)],
    rhs = rhs_vec[seq_len(n_used)],
    lhs_total = length(lhs_vec),
    rhs_total = length(rhs_vec),
    n_used = n_used
  )
}

compare_numeric_objects <- function(lhs, rhs, metric_name) {
  aligned <- flatten_pair_values(lhs, rhs)
  lhs_vec <- aligned$lhs
  rhs_vec <- aligned$rhs
  if (!length(lhs_vec) || !length(rhs_vec)) {
    return(tibble::tibble(
      metric = metric_name,
      n = 0L,
      lhs_mean = NA_real_,
      rhs_mean = NA_real_,
      mean_diff = NA_real_,
      mean_abs_diff = NA_real_,
      max_abs_diff = NA_real_,
      correlation = NA_real_,
      lhs_negative = NA_integer_,
      rhs_negative = NA_integer_,
      lhs_total = aligned$lhs_total,
      rhs_total = aligned$rhs_total
    ))
  }

  diff_vec <- rhs_vec - lhs_vec
  tibble::tibble(
    metric = metric_name,
    n = aligned$n_used,
    lhs_mean = mean(lhs_vec, na.rm = TRUE),
    rhs_mean = mean(rhs_vec, na.rm = TRUE),
    mean_diff = mean(diff_vec, na.rm = TRUE),
    mean_abs_diff = mean(abs(diff_vec), na.rm = TRUE),
    max_abs_diff = max(abs(diff_vec), na.rm = TRUE),
    correlation = if (aligned$n_used >= 2L) {
      suppressWarnings(stats::cor(lhs_vec, rhs_vec, use = "complete.obs"))
    } else {
      NA_real_
    },
    lhs_negative = sum(lhs_vec < 0, na.rm = TRUE),
    rhs_negative = sum(rhs_vec < 0, na.rm = TRUE),
    lhs_total = aligned$lhs_total,
    rhs_total = aligned$rhs_total
  )
}

summarize_sign_changes <- function(lhs, rhs, metric_name) {
  aligned <- flatten_pair_values(lhs, rhs)
  lhs_vec <- aligned$lhs
  rhs_vec <- aligned$rhs
  if (!length(lhs_vec) || !length(rhs_vec)) {
    return(tibble::tibble(
      metric = metric_name,
      n = 0L,
      lhs_negative = NA_integer_,
      rhs_negative = NA_integer_,
      neg_to_nonneg = NA_integer_,
      nonneg_to_neg = NA_integer_,
      across_zero = NA_integer_
    ))
  }

  tibble::tibble(
    metric = metric_name,
    n = aligned$n_used,
    lhs_negative = sum(lhs_vec < 0, na.rm = TRUE),
    rhs_negative = sum(rhs_vec < 0, na.rm = TRUE),
    neg_to_nonneg = sum(lhs_vec < 0 & rhs_vec >= 0, na.rm = TRUE),
    nonneg_to_neg = sum(lhs_vec >= 0 & rhs_vec < 0, na.rm = TRUE),
    across_zero = sum(lhs_vec * rhs_vec < 0, na.rm = TRUE)
  )
}

merge_landscape_pair <- function(lhs_landscape, rhs_landscape) {
  lhs_tbl <- lhs_landscape %>%
    dplyr::transmute(
      k,
      lhs_mean = mean,
      lhs_median = median,
      lhs_sd = sd,
      lhs_fq = as.logical(fq),
      lhs_nn = as.logical(nn)
    )
  rhs_tbl <- rhs_landscape %>%
    dplyr::transmute(
      k,
      rhs_mean = mean,
      rhs_median = median,
      rhs_sd = sd,
      rhs_fq = as.logical(fq),
      rhs_nn = as.logical(nn)
    )

  dplyr::inner_join(lhs_tbl, rhs_tbl, by = "k") %>%
    dplyr::mutate(
      state_group = dplyr::case_when(
        lhs_fq | rhs_fq ~ "fq",
        lhs_nn | rhs_nn ~ "nn",
        TRUE ~ "other"
      )
    )
}

summarize_landscape_groups <- function(merged_landscape, lhs_col, rhs_col, metric_name) {
  if (!nrow(merged_landscape)) {
    return(tibble::tibble())
  }

  diff_vec <- merged_landscape[[rhs_col]] - merged_landscape[[lhs_col]]
  tmp <- merged_landscape %>%
    dplyr::mutate(
      lhs_value = .data[[lhs_col]],
      rhs_value = .data[[rhs_col]],
      diff_value = diff_vec
    ) %>%
    dplyr::group_by(state_group) %>%
    dplyr::summarise(
      metric = metric_name,
      n = dplyr::n(),
      lhs_mean = mean(lhs_value, na.rm = TRUE),
      rhs_mean = mean(rhs_value, na.rm = TRUE),
      mean_diff = mean(diff_value, na.rm = TRUE),
      mean_abs_diff = mean(abs(diff_value), na.rm = TRUE),
      max_abs_diff = max(abs(diff_value), na.rm = TRUE),
      correlation = if (dplyr::n() >= 2L) {
        suppressWarnings(stats::cor(lhs_value, rhs_value, use = "complete.obs"))
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  tmp %>% dplyr::select(metric, state_group, dplyr::everything())
}

top_landscape_shifts <- function(merged_landscape, top_n = 15L) {
  if (!nrow(merged_landscape)) {
    return(tibble::tibble())
  }

  merged_landscape %>%
    dplyr::transmute(
      k,
      state_group,
      lhs_mean,
      rhs_mean,
      diff = rhs_mean - lhs_mean
    ) %>%
    dplyr::arrange(dplyr::desc(abs(diff))) %>%
    dplyr::slice_head(n = top_n)
}

compare_fit_pair <- function(lhs_dir,
                             rhs_dir,
                             patient_id,
                             minobs,
                             pm,
                             lhs_label,
                             rhs_label,
                             lhs_setting,
                             rhs_setting,
                             comparison_set,
                             beneficial_move_levels,
                             top_n = 15L,
                             include_posterior = TRUE) {
  lhs_bundle <- read_fit_bundle(lhs_dir, beneficial_move_levels = beneficial_move_levels)
  rhs_bundle <- read_fit_bundle(rhs_dir, beneficial_move_levels = beneficial_move_levels)

  component_tbl <- dplyr::bind_rows(
    compare_numeric_objects(lhs_bundle$bootstrap$initial_fitness, rhs_bundle$bootstrap$initial_fitness, "bootstrap_initial_fitness"),
    compare_numeric_objects(lhs_bundle$bootstrap$final_fitness, rhs_bundle$bootstrap$final_fitness, "bootstrap_final_fitness"),
    compare_numeric_objects(lhs_bundle$bootstrap$initial_frequencies, rhs_bundle$bootstrap$initial_frequencies, "bootstrap_initial_frequencies"),
    compare_numeric_objects(lhs_bundle$bootstrap$final_frequencies, rhs_bundle$bootstrap$final_frequencies, "bootstrap_final_frequencies"),
    compare_numeric_objects(lhs_bundle$bootstrap$nn_fitness, rhs_bundle$bootstrap$nn_fitness, "bootstrap_nn_fitness"),
    compare_numeric_objects(lhs_bundle$landscape$mean, rhs_bundle$landscape$mean, "landscape_mean"),
    compare_numeric_objects(lhs_bundle$landscape$median, rhs_bundle$landscape$median, "landscape_median"),
    compare_numeric_objects(lhs_bundle$landscape$sd, rhs_bundle$landscape$sd, "landscape_sd"),
    compare_numeric_objects(lhs_bundle$beneficial$proportion, rhs_bundle$beneficial$proportion, "beneficial_proportion"),
    compare_numeric_objects(lhs_bundle$beneficial$valid_n, rhs_bundle$beneficial$valid_n, "beneficial_valid_n"),
    if (include_posterior) {
      compare_numeric_objects(lhs_bundle$posterior, rhs_bundle$posterior, "landscape_posterior_samples")
    },
    compare_numeric_objects(
      extract_xval_metrics(lhs_bundle$xval)$xval_r2,
      extract_xval_metrics(rhs_bundle$xval)$xval_r2,
      "xval_r2"
    )
  ) %>%
    dplyr::mutate(
      comparison_set = comparison_set,
      patient_id = patient_id,
      minobs = minobs,
      pm = pm,
      lhs_label = lhs_label,
      rhs_label = rhs_label,
      lhs_setting = lhs_setting,
      rhs_setting = rhs_setting
    ) %>%
    dplyr::relocate(comparison_set, patient_id, minobs, pm, lhs_label, rhs_label, lhs_setting, rhs_setting, metric)

  sign_tbl <- dplyr::bind_rows(
    summarize_sign_changes(lhs_bundle$bootstrap$nn_fitness, rhs_bundle$bootstrap$nn_fitness, "bootstrap_nn_fitness"),
    summarize_sign_changes(lhs_bundle$landscape$mean, rhs_bundle$landscape$mean, "landscape_mean"),
    summarize_sign_changes(lhs_bundle$landscape$median, rhs_bundle$landscape$median, "landscape_median")
  ) %>%
    dplyr::mutate(
      comparison_set = comparison_set,
      patient_id = patient_id,
      minobs = minobs,
      pm = pm,
      lhs_label = lhs_label,
      rhs_label = rhs_label,
      lhs_setting = lhs_setting,
      rhs_setting = rhs_setting
    ) %>%
    dplyr::relocate(comparison_set, patient_id, minobs, pm, lhs_label, rhs_label, lhs_setting, rhs_setting, metric)

  merged_landscape <- if (!is.null(lhs_bundle$landscape) && !is.null(rhs_bundle$landscape)) {
    merge_landscape_pair(lhs_bundle$landscape, rhs_bundle$landscape)
  } else {
    tibble::tibble()
  }
  landscape_group_tbl <- if (nrow(merged_landscape)) {
    dplyr::bind_rows(
      summarize_landscape_groups(merged_landscape, "lhs_mean", "rhs_mean", "landscape_mean"),
      summarize_landscape_groups(merged_landscape, "lhs_median", "rhs_median", "landscape_median"),
      summarize_landscape_groups(merged_landscape, "lhs_sd", "rhs_sd", "landscape_sd")
    ) %>%
      dplyr::mutate(
        comparison_set = comparison_set,
        patient_id = patient_id,
        minobs = minobs,
        pm = pm,
        lhs_label = lhs_label,
        rhs_label = rhs_label,
        lhs_setting = lhs_setting,
        rhs_setting = rhs_setting
      ) %>%
      dplyr::relocate(comparison_set, patient_id, minobs, pm, lhs_label, rhs_label, lhs_setting, rhs_setting, metric, state_group)
  } else {
    tibble::tibble()
  }

  top_shift_tbl <- if (nrow(merged_landscape)) {
    top_landscape_shifts(merged_landscape, top_n = top_n) %>%
      dplyr::mutate(
        comparison_set = comparison_set,
        patient_id = patient_id,
        minobs = minobs,
        pm = pm,
        lhs_label = lhs_label,
        rhs_label = rhs_label,
        lhs_setting = lhs_setting,
        rhs_setting = rhs_setting
      ) %>%
      dplyr::relocate(comparison_set, patient_id, minobs, pm, lhs_label, rhs_label, lhs_setting, rhs_setting)
  } else {
    tibble::tibble()
  }

  list(
    component = component_tbl,
    sign = sign_tbl,
    landscape_group = landscape_group_tbl,
    top_shift = top_shift_tbl
  )
}

build_pairwise_comparisons <- function(results_tbl,
                                       pair_values,
                                       setting_col,
                                       comparison_set,
                                       beneficial_move_levels,
                                       top_n = 15L,
                                       include_posterior = TRUE) {
  if (!nrow(results_tbl) || length(pair_values) < 2L) {
    return(list(
      component = tibble::tibble(),
      sign = tibble::tibble(),
      landscape_group = tibble::tibble(),
      top_shift = tibble::tibble()
    ))
  }

  pair_defs <- combn(pair_values, 2, simplify = FALSE)
  all_component <- list()
  all_sign <- list()
  all_group <- list()
  all_shift <- list()

  for (pair in pair_defs) {
    lhs_val <- pair[[1]]
    rhs_val <- pair[[2]]

    lhs_idx <- results_tbl %>%
      dplyr::filter(status == "ok", .data[[setting_col]] == lhs_val) %>%
      dplyr::transmute(
        patient_id,
        minobs,
        pm,
        lhs_setting = .data[[setting_col]],
        lhs_label = if (setting_col == "nn_prior_grid_n") paste0("grid_", .data[[setting_col]]) else as.character(.data[[setting_col]]),
        outdir_lhs = outdir
      )

    rhs_idx <- results_tbl %>%
      dplyr::filter(status == "ok", .data[[setting_col]] == rhs_val) %>%
      dplyr::transmute(
        patient_id,
        minobs,
        pm,
        rhs_setting = .data[[setting_col]],
        rhs_label = if (setting_col == "nn_prior_grid_n") paste0("grid_", .data[[setting_col]]) else as.character(.data[[setting_col]]),
        outdir_rhs = outdir
      )

    joined <- dplyr::inner_join(lhs_idx, rhs_idx, by = c("patient_id", "minobs", "pm"))
    if (!nrow(joined)) {
      next
    }

    pair_res <- lapply(seq_len(nrow(joined)), function(i) {
      rr <- joined[i, , drop = FALSE]
      compare_fit_pair(
        lhs_dir = rr$outdir_lhs,
        rhs_dir = rr$outdir_rhs,
        patient_id = rr$patient_id,
        minobs = rr$minobs,
        pm = rr$pm,
        lhs_label = rr$lhs_label,
        rhs_label = rr$rhs_label,
        lhs_setting = rr$lhs_setting,
        rhs_setting = rr$rhs_setting,
        comparison_set = comparison_set,
        beneficial_move_levels = beneficial_move_levels,
        top_n = top_n,
        include_posterior = include_posterior
      )
    })

    all_component <- c(all_component, lapply(pair_res, `[[`, "component"))
    all_sign <- c(all_sign, lapply(pair_res, `[[`, "sign"))
    all_group <- c(all_group, lapply(pair_res, `[[`, "landscape_group"))
    all_shift <- c(all_shift, lapply(pair_res, `[[`, "top_shift"))
  }

  list(
    component = dplyr::bind_rows(all_component),
    sign = dplyr::bind_rows(all_sign),
    landscape_group = dplyr::bind_rows(all_group),
    top_shift = dplyr::bind_rows(all_shift)
  )
}
