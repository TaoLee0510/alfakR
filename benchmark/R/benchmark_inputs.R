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
