available_pid_dirs <- function(base_dir) {
  dd <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  dd <- dd[file.info(dd)$isdir %in% TRUE]
  dd <- dd[grepl("^P[0-9]+$", basename(dd))]
  sort_pid_levels(basename(dd))
}

read_transcriptome_perspective <- function(rds_path) {
  obj <- readRDS(rds_path)
  stopifnot(is.list(obj), "profile" %in% names(obj))
  profile <- as.matrix(obj$profile)
  storage.mode(profile) <- "numeric"
  list(
    passaging_id = obj$passaging_id,
    perspective_type = obj$perspective_type,
    profile = round(profile)
  )
}

compute_delta_time_days <- function(df_pid, stage_levels) {
  df_pid <- df_pid %>%
    dplyr::mutate(Stage = factor(Stage, levels = stage_levels)) %>%
    dplyr::arrange(Stage, effective_id)

  age_primary <- suppressWarnings(as.numeric(df_pid$Age[df_pid$Stage == "Primary"][1]))
  age_recurrent <- suppressWarnings(as.numeric(df_pid$Age[df_pid$Stage == "Recurrent"][1]))

  recurrent_days <- if (is.finite(age_primary) && is.finite(age_recurrent)) {
    as.integer(round(abs(age_recurrent - age_primary) * 360))
  } else {
    NA_integer_
  }
  if (!is.finite(recurrent_days) || recurrent_days <= 1L) {
    recurrent_days <- 180L
  }

  df_pid %>%
    dplyr::mutate(
      Delta_Time = ifelse(as.character(Stage) == "Primary", 1L, recurrent_days)
    )
}

build_patient_manifest <- function(meta_tbl, base_dir, stage_levels, patient_subset = NULL) {
  pid_dirs <- available_pid_dirs(base_dir)
  manifest <- meta_tbl %>%
    dplyr::mutate(
      effective_id = ifelse(!is.na(ID_alt) & nzchar(ID_alt), ID_alt, ID),
      pid = as.character(pid),
      Stage = as.character(Stage),
      Age = suppressWarnings(as.numeric(Age))
    ) %>%
    dplyr::filter(pid %in% pid_dirs)

  if (!is.null(patient_subset)) {
    manifest <- manifest %>% dplyr::filter(pid %in% patient_subset)
  }

  manifest %>%
    dplyr::group_by(pid) %>%
    dplyr::group_modify(function(.x, .y) compute_delta_time_days(.x, stage_levels = stage_levels)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      rds_path = file.path(base_dir, pid, "perspectives", effective_id, "TranscriptomePerspective.rds"),
      exists = file.exists(rds_path)
    ) %>%
    dplyr::arrange(factor(pid, levels = sort_pid_levels(pid)), factor(Stage, levels = stage_levels))
}

profile_to_karyotype_strings <- function(profile_mat) {
  if (nrow(profile_mat) != 22) {
    warning("Expected 22 chromosome rows; got ", nrow(profile_mat))
  }
  x <- t(round(as.matrix(profile_mat)))
  apply(x, 1, paste, collapse = ".")
}

build_count_matrix_from_profiles <- function(patient_manifest, stage_levels) {
  stopifnot(length(unique(patient_manifest$pid)) == 1L)
  rows_ok <- patient_manifest[patient_manifest$exists, , drop = FALSE]
  rows_ok <- rows_ok[order(match(rows_ok$Stage, stage_levels)), , drop = FALSE]
  if (nrow(rows_ok) < 2) {
    stop("Need at least two available samples for patient ", unique(patient_manifest$pid))
  }

  long_df <- do.call(rbind, lapply(seq_len(nrow(rows_ok)), function(i) {
    rr <- rows_ok[i, , drop = FALSE]
    obj <- read_transcriptome_perspective(rr$rds_path)
    if (!is.null(obj$passaging_id) &&
        !identical(as.character(obj$passaging_id), as.character(rr$effective_id))) {
      warning("passaging_id mismatch for ", rr$pid, " / ", rr$effective_id)
    }
    tibble::tibble(
      karyotype = profile_to_karyotype_strings(obj$profile),
      Delta_Time = rr$Delta_Time,
      Stage = rr$Stage,
      effective_id = rr$effective_id
    )
  }))

  delta_levels <- sort(unique(rows_ok$Delta_Time))
  count_mat <- as.data.frame.matrix(table(
    long_df$karyotype,
    factor(long_df$Delta_Time, levels = delta_levels)
  ))
  colnames(count_mat) <- as.character(delta_levels)
  count_mat <- count_mat[order(rowSums(count_mat), decreasing = TRUE), , drop = FALSE]

  list(
    count_matrix = count_mat,
    long_df = long_df
  )
}

save_alfak_input <- function(count_matrix, out_path) {
  yi <- list(
    x = count_matrix,
    pop.fitness = NULL,
    dt = 1
  )
  saveRDS(yi, out_path)
  invisible(out_path)
}

build_benchmark_inputs <- function(meta_tbl,
                                   base_dir,
                                   input_dir,
                                   tables_dir,
                                   stage_levels,
                                   diploid_state,
                                   rebuild_inputs = FALSE,
                                   patient_subset = NULL) {
  manifest_tbl <- build_patient_manifest(
    meta_tbl = meta_tbl,
    base_dir = base_dir,
    stage_levels = stage_levels,
    patient_subset = patient_subset
  )
  save_table_bundle(manifest_tbl, file.path(tables_dir, "benchmark_patient_manifest"))

  patient_ids <- sort_pid_levels(unique(manifest_tbl$pid))
  if (!length(patient_ids)) {
    stop("No benchmark patients discovered under ", base_dir)
  }

  rows <- lapply(patient_ids, function(pid) {
    patient_manifest <- manifest_tbl %>% dplyr::filter(pid == !!pid)
    input_rds <- file.path(input_dir, paste0(pid, ".Rds"))
    count_matrix_path <- file.path(input_dir, paste0(pid, "_count_matrix.tsv"))
    long_df_path <- file.path(input_dir, paste0(pid, "_profile_long.tsv"))

    if (rebuild_inputs || !file.exists(input_rds) || !file.exists(count_matrix_path) || !file.exists(long_df_path)) {
      built <- build_count_matrix_from_profiles(patient_manifest, stage_levels = stage_levels)
      count_export <- data.frame(
        karyotype = rownames(built$count_matrix),
        built$count_matrix,
        check.names = FALSE
      )
      save_alfak_input(built$count_matrix, input_rds)
      write_tsv_base(count_export, count_matrix_path)
      write_tsv_base(built$long_df, long_df_path)
    }

    yi <- readRDS(input_rds)
    yi$x <- as.data.frame(yi$x)
    yi_minobs5 <- yi$x
    if (diploid_state %in% rownames(yi_minobs5)) {
      yi_minobs5 <- yi_minobs5[rownames(yi_minobs5) != diploid_state, , drop = FALSE]
    }
    input_row_count_minobs5 <- if (nrow(yi_minobs5)) {
      sum(rowSums(yi_minobs5, na.rm = TRUE) >= 5L)
    } else {
      0L
    }

    tibble::tibble(
      patient_id = pid,
      input_rds = input_rds,
      count_matrix_path = count_matrix_path,
      long_df_path = long_df_path,
      input_row_count_minobs5 = input_row_count_minobs5,
      input_row_count = nrow(yi$x),
      n_timepoints = ncol(yi$x),
      n_karyotypes = nrow(yi$x),
      total_cells = sum(as.matrix(yi$x), na.rm = TRUE),
      delta_time_labels = paste(colnames(yi$x), collapse = ","),
      effective_ids = paste(patient_manifest$effective_id[patient_manifest$exists], collapse = ","),
      stage_sequence = paste(as.character(patient_manifest$Stage[patient_manifest$exists]), collapse = " -> ")
    )
  })

  dplyr::bind_rows(rows) %>%
    dplyr::arrange(factor(patient_id, levels = sort_pid_levels(patient_id)))
}

prepare_input_count_matrix <- function(input_rds, diploid_state) {
  yi <- readRDS(input_rds)
  x <- as.matrix(yi$x)
  storage.mode(x) <- "numeric"

  if (diploid_state %in% rownames(x)) {
    x <- x[rownames(x) != diploid_state, , drop = FALSE]
  }
  if (!nrow(x)) {
    stop("No non-diploid karyotypes remain after filtering: ", input_rds)
  }
  if (ncol(x) < 2L) {
    stop("Expected at least two timepoints in benchmark input: ", input_rds)
  }
  if (ncol(x) > 2L) {
    x <- x[, seq_len(2L), drop = FALSE]
  }

  x
}

safe_fraction <- function(num, den) {
  if (!length(den) || is.na(den) || !is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  as.numeric(num) / as.numeric(den)
}

select_frequent_karyotypes_minobs <- function(x, minobs) {
  rownames(x)[rowSums(x, na.rm = TRUE) >= as.numeric(minobs)]
}

neighbor_strings_from_karyotypes <- function(karyotypes) {
  karyotypes <- unique(as.character(karyotypes))
  karyotypes <- karyotypes[nzchar(karyotypes)]
  if (!length(karyotypes)) {
    return(character(0))
  }

  nn_mat <- alfakR:::gen_all_neighbours(karyotypes)
  if (!nrow(nn_mat)) {
    return(character(0))
  }

  unique(apply(nn_mat, 1, paste, collapse = "."))
}

build_state_change_tbl <- function(state_ids, x) {
  state_ids <- unique(as.character(state_ids))
  state_ids <- state_ids[nzchar(state_ids)]
  if (!length(state_ids)) {
    return(tibble::tibble(
      k = character(0),
      observed = logical(0),
      count_t1 = numeric(0),
      count_t2 = numeric(0),
      prop_t1 = numeric(0),
      prop_t2 = numeric(0),
      count_up = logical(0),
      prop_direction = factor(character(0), levels = c("up", "down", "flat"))
    ))
  }

  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x2 <- x[, seq_len(2L), drop = FALSE]

  count_mat <- matrix(0, nrow = length(state_ids), ncol = 2L, dimnames = list(state_ids, colnames(x2)))
  matched <- match(state_ids, rownames(x2))
  present <- !is.na(matched)
  if (any(present)) {
    count_mat[present, ] <- x2[matched[present], , drop = FALSE]
  }

  totals <- colSums(x2, na.rm = TRUE)
  prop_t1 <- if (totals[1] > 0) count_mat[, 1] / totals[1] else rep(NA_real_, length(state_ids))
  prop_t2 <- if (totals[2] > 0) count_mat[, 2] / totals[2] else rep(NA_real_, length(state_ids))
  prop_direction <- ifelse(prop_t2 > prop_t1, "up", ifelse(prop_t2 < prop_t1, "down", "flat"))

  tibble::tibble(
    k = state_ids,
    observed = present,
    count_t1 = as.numeric(count_mat[, 1]),
    count_t2 = as.numeric(count_mat[, 2]),
    prop_t1 = as.numeric(prop_t1),
    prop_t2 = as.numeric(prop_t2),
    count_up = as.numeric(count_mat[, 2]) > as.numeric(count_mat[, 1]),
    prop_direction = factor(prop_direction, levels = c("up", "down", "flat"))
  )
}

summarize_input_fq_nn_overview <- function(input_index_tbl, minobs_values, diploid_state) {
  if (is.null(input_index_tbl) || !nrow(input_index_tbl)) {
    return(list(
      input_fq_nn_summary_tbl = tibble::tibble(),
      input_fq_group_nn_summary_tbl = tibble::tibble()
    ))
  }

  overview_rows <- vector("list", nrow(input_index_tbl) * length(minobs_values))
  group_rows <- vector("list", nrow(input_index_tbl) * length(minobs_values) * 2L)
  overview_idx <- 0L
  group_idx <- 0L

  for (i in seq_len(nrow(input_index_tbl))) {
    rr <- input_index_tbl[i, , drop = FALSE]
    x <- prepare_input_count_matrix(rr$input_rds, diploid_state = diploid_state)
    time_labels <- colnames(x)[seq_len(2L)]

    for (minobs in minobs_values) {
      fq <- select_frequent_karyotypes_minobs(x, minobs)
      nn <- setdiff(neighbor_strings_from_karyotypes(fq), fq)

      fq_change_tbl <- build_state_change_tbl(fq, x)
      nn_change_tbl <- build_state_change_tbl(nn, x)

      overview_idx <- overview_idx + 1L
      overview_rows[[overview_idx]] <- tibble::tibble(
        patient_id = as.character(rr$patient_id),
        minobs = as.integer(minobs),
        time1_label = as.character(time_labels[1]),
        time2_label = as.character(time_labels[2]),
        n_non_diploid_karyotypes = nrow(x),
        n_fq = length(fq),
        n_nn = length(nn),
        n_nn_observed = sum(nn_change_tbl$observed %in% TRUE, na.rm = TRUE),
        prop_nn_observed = safe_divide(sum(nn_change_tbl$observed %in% TRUE, na.rm = TRUE), nrow(nn_change_tbl)),
        n_fq_count_up = sum(fq_change_tbl$count_up %in% TRUE, na.rm = TRUE),
        prop_fq_count_up = safe_divide(sum(fq_change_tbl$count_up %in% TRUE, na.rm = TRUE), nrow(fq_change_tbl)),
        n_nn_count_up = sum(nn_change_tbl$count_up %in% TRUE, na.rm = TRUE),
        prop_nn_count_up = safe_divide(sum(nn_change_tbl$count_up %in% TRUE, na.rm = TRUE), nrow(nn_change_tbl)),
        n_fq_prop_up = sum(as.character(fq_change_tbl$prop_direction) == "up", na.rm = TRUE),
        n_fq_prop_down = sum(as.character(fq_change_tbl$prop_direction) == "down", na.rm = TRUE),
        n_fq_prop_flat = sum(as.character(fq_change_tbl$prop_direction) == "flat", na.rm = TRUE)
      )

      for (fq_prop_direction in c("up", "down")) {
        fq_group <- fq_change_tbl$k[as.character(fq_change_tbl$prop_direction) == fq_prop_direction]
        group_nn <- if (length(fq_group)) {
          setdiff(neighbor_strings_from_karyotypes(fq_group), fq)
        } else {
          character(0)
        }
        group_nn_change_tbl <- build_state_change_tbl(group_nn, x)

        group_idx <- group_idx + 1L
        group_rows[[group_idx]] <- tibble::tibble(
          patient_id = as.character(rr$patient_id),
          minobs = as.integer(minobs),
          fq_prop_direction = factor(fq_prop_direction, levels = c("up", "down")),
          n_fq_in_group = length(fq_group),
          n_group_nn = length(group_nn),
          n_group_nn_observed = sum(group_nn_change_tbl$observed %in% TRUE, na.rm = TRUE),
          n_group_nn_count_up = sum(group_nn_change_tbl$count_up %in% TRUE, na.rm = TRUE),
          prop_group_nn_count_up = safe_divide(sum(group_nn_change_tbl$count_up %in% TRUE, na.rm = TRUE), nrow(group_nn_change_tbl))
        )
      }
    }
  }

  input_fq_nn_summary_tbl <- dplyr::bind_rows(overview_rows[seq_len(overview_idx)]) %>%
    dplyr::arrange(factor(patient_id, levels = sort_pid_levels(patient_id)), minobs)

  input_fq_group_nn_summary_tbl <- dplyr::bind_rows(group_rows[seq_len(group_idx)]) %>%
    dplyr::arrange(factor(patient_id, levels = sort_pid_levels(patient_id)), minobs, fq_prop_direction)

  list(
    input_fq_nn_summary_tbl = input_fq_nn_summary_tbl,
    input_fq_group_nn_summary_tbl = input_fq_group_nn_summary_tbl
  )
}
