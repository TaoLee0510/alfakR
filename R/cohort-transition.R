#' Resolve cached two-shell fit directories
#'
#' Builds the expected two-shell cache paths without altering patient or sample
#' names. The layout is
#' `<two_shell_root>/<pm_tag>/<minobs_tag>/<sample_name>/`.
#'
#' @param two_shell_root Root directory containing upstream two-shell results.
#' @param patient_ids Character vector of patient identifiers.
#' @param sample_names Optional character vector of sample directory names,
#'   aligned with `patient_ids`.
#' @param pm Numeric missegregation probability used for the upstream fit.
#' @param minobs Integer MINIOBS value used for the upstream fit.
#' @param pm_tag Optional exact PM directory name. When `NULL`, the package
#'   default `pm_<pm>` form is used, with a numeric directory search fallback.
#' @param minobs_tag Optional exact MINIOBS directory name. When `NULL`,
#'   `MINIOBS<minobs>` is used.
#' @param sample_map Optional named character vector mapping patient IDs to
#'   sample directory names.
#' @return A data frame describing expected cache paths and existence.
#' @export
resolve_two_shell_fit_dirs <- function(two_shell_root,
                                       patient_ids,
                                       sample_names = NULL,
                                       pm,
                                       minobs,
                                       pm_tag = NULL,
                                       minobs_tag = NULL,
                                       sample_map = NULL) {
  if (is.null(two_shell_root) || length(two_shell_root) != 1L || !nzchar(two_shell_root)) {
    stop("`two_shell_root` must be a single non-empty path.", call. = FALSE)
  }
  patient_ids <- as.character(patient_ids)
  if (!length(patient_ids) || any(!nzchar(patient_ids)) || anyDuplicated(patient_ids)) {
    stop("`patient_ids` must be non-empty unique character values.", call. = FALSE)
  }
  validate_probability(pm, "pm", upper_inclusive = TRUE)
  validate_positive_integer(minobs, "minobs")

  if (!is.null(sample_names)) {
    sample_names <- as.character(sample_names)
    if (length(sample_names) != length(patient_ids) || any(!nzchar(sample_names))) {
      stop("`sample_names` must be a non-empty character vector aligned with `patient_ids`.", call. = FALSE)
    }
  } else {
    sample_names <- patient_ids
  }

  if (!is.null(sample_map)) {
    if (is.null(names(sample_map)) || any(!nzchar(names(sample_map)))) {
      stop("`sample_map` must be a named character vector.", call. = FALSE)
    }
    sample_map <- as.character(sample_map)
    mapped <- match(patient_ids, names(sample_map))
    has_map <- !is.na(mapped)
    sample_names[has_map] <- sample_map[mapped[has_map]]
  }

  if (is.null(pm_tag)) {
    pm_tag <- paste0("pm_", format(pm, scientific = FALSE, trim = TRUE))
    pm_tag_resolution <- "derived"
    derived_pm_dir <- file.path(two_shell_root, pm_tag)
    if (!dir.exists(derived_pm_dir) && dir.exists(two_shell_root)) {
      pm_dirs <- list.dirs(two_shell_root, recursive = FALSE, full.names = FALSE)
      pm_dirs <- pm_dirs[grepl("^pm_", pm_dirs)]
      pm_values <- suppressWarnings(as.numeric(sub("^pm_", "", pm_dirs)))
      tol <- max(1e-12, abs(pm) * 1e-8)
      matches <- pm_dirs[is.finite(pm_values) & abs(pm_values - pm) <= tol]
      if (length(matches) == 1L) {
        pm_tag <- matches
        pm_tag_resolution <- "numeric_match"
      } else if (length(matches) > 1L) {
        stop(
          sprintf(
            "Multiple `pm_*` directories under `two_shell_root` match pm=%s; supply `two_shell_pm_tag`.",
            format(pm, scientific = FALSE, trim = TRUE)
          ),
          call. = FALSE
        )
      }
    }
  } else {
    pm_tag <- as.character(pm_tag)
    if (length(pm_tag) != 1L || !nzchar(pm_tag)) {
      stop("`pm_tag` must be a single non-empty string.", call. = FALSE)
    }
    pm_tag_resolution <- "provided"
  }

  if (is.null(minobs_tag)) {
    minobs_tag <- paste0("MINIOBS", as.integer(minobs))
    minobs_tag_resolution <- "derived"
  } else {
    minobs_tag <- as.character(minobs_tag)
    if (length(minobs_tag) != 1L || !nzchar(minobs_tag)) {
      stop("`minobs_tag` must be a single non-empty string.", call. = FALSE)
    }
    minobs_tag_resolution <- "provided"
  }

  expected_fit_dir <- file.path(two_shell_root, pm_tag, minobs_tag, sample_names)
  exists <- dir.exists(expected_fit_dir)
  out <- data.frame(
    patient_id = patient_ids,
    sample_name = sample_names,
    pm_tag = rep(pm_tag, length(patient_ids)),
    minobs_tag = rep(minobs_tag, length(patient_ids)),
    expected_fit_dir = expected_fit_dir,
    exists = exists,
    status = ifelse(exists, "exists", "missing_dir"),
    stringsAsFactors = FALSE
  )
  attr(out, "pm_tag_resolution") <- pm_tag_resolution
  attr(out, "minobs_tag_resolution") <- minobs_tag_resolution
  out
}

#' Check integrity of a cached two-shell fit
#'
#' @param fit_dir Directory containing cached two-shell output.
#' @param patient_id Optional patient identifier for messages.
#' @param mode Integrity strictness: `"strict"`, `"basic"`, or `"none"`.
#' @return A structured list with status, missing files, unreadable files, and warnings.
#' @export
check_two_shell_fit_integrity <- function(fit_dir,
                                          patient_id = NULL,
                                          mode = c("strict", "basic", "none")) {
  mode <- match.arg(mode)
  result <- list(
    ok = FALSE,
    status = NA_character_,
    missing_files = character(0),
    unreadable_files = character(0),
    warnings = character(0),
    fit_dir = fit_dir,
    patient_id = if (is.null(patient_id)) NA_character_ else as.character(patient_id)
  )

  if (mode == "none") {
    result$ok <- TRUE
    result$status <- "skipped"
    result$warnings <- "integrity_check_skipped"
    return(result)
  }

  if (is.null(fit_dir) || length(fit_dir) != 1L || !dir.exists(fit_dir)) {
    result$status <- "missing_dir"
    return(result)
  }

  required <- c("bootstrap_res.Rds", "landscape.Rds")
  if (mode == "strict") {
    required <- c(required, "nn_prior_diagnostics.Rds")
  }
  optional <- "landscape_posterior_samples.Rds"

  required_paths <- file.path(fit_dir, required)
  missing <- required[!file.exists(required_paths)]
  if (length(missing)) {
    result$status <- "missing_file"
    result$missing_files <- missing
    return(result)
  }
  if (mode == "strict" && !file.exists(file.path(fit_dir, optional))) {
    result$warnings <- c(result$warnings, "missing_optional_landscape_posterior_samples")
  }

  read_cached <- function(file) {
    path <- file.path(fit_dir, file)
    tryCatch(readRDS(path), error = function(e) structure(list(error = conditionMessage(e)), class = "alfak_read_error"))
  }
  objects <- lapply(required, read_cached)
  names(objects) <- required
  unreadable <- names(objects)[vapply(objects, inherits, logical(1), "alfak_read_error")]
  if (length(unreadable)) {
    result$status <- "unreadable_file"
    result$unreadable_files <- unreadable
    return(result)
  }
  non_empty <- vapply(objects, function(x) {
    if (is.null(x)) return(FALSE)
    if (is.data.frame(x) || is.matrix(x)) return(nrow(x) > 0 || ncol(x) > 0)
    if (is.list(x)) return(length(x) > 0)
    length(x) > 0
  }, logical(1))
  if (!all(non_empty)) {
    result$status <- "invalid_object"
    result$warnings <- c(result$warnings, paste0("empty_object:", names(non_empty)[!non_empty]))
    return(result)
  }

  boot <- objects[["bootstrap_res.Rds"]]
  landscape <- objects[["landscape.Rds"]]
  if (!is.list(boot) ||
      is.null(boot$final_fitness) ||
      is.null(boot$nn_fitness) ||
      !is.matrix(boot$final_fitness) ||
      !is.matrix(boot$nn_fitness)) {
    result$status <- "invalid_object"
    result$warnings <- c(result$warnings, "bootstrap_missing_expected_fitness_matrices")
    return(result)
  }
  if (!is.data.frame(landscape) ||
      !"k" %in% names(landscape) ||
      !any(c("mean", "median", "sd") %in% names(landscape))) {
    result$status <- "invalid_object"
    result$warnings <- c(result$warnings, "landscape_missing_expected_fitness_fields")
    return(result)
  }

  if (mode == "strict") {
    diag <- objects[["nn_prior_diagnostics.Rds"]]
    diag_frames <- list()
    if (is.data.frame(diag)) {
      diag_frames <- list(diag)
    } else if (is.list(diag)) {
      diag_frames <- diag[vapply(diag, is.data.frame, logical(1))]
    }
    if (!length(diag_frames)) {
      result$status <- "invalid_object"
      result$warnings <- c(result$warnings, "nn_prior_diagnostics_unrecognized")
      return(result)
    }
    diag_all <- do.call(rbind, lapply(diag_frames, function(x) {
      x[, intersect(names(x), c("nn_prior_mode_requested", "nn_prior_mode_used", "nn_prior_source_used", "mu01", "sigma01")), drop = FALSE]
    }))
    prior_cols <- intersect(names(diag_all), c("nn_prior_mode_requested", "nn_prior_mode_used"))
    if (length(prior_cols)) {
      prior_vals <- unique(unlist(diag_all[prior_cols], use.names = FALSE))
      prior_vals <- prior_vals[!is.na(prior_vals)]
      if (length(prior_vals) && !"empirical_two_shell" %in% prior_vals) {
        result$status <- "wrong_prior_mode"
        result$warnings <- c(result$warnings, paste0("prior_modes_seen:", paste(prior_vals, collapse = ",")))
        return(result)
      }
    }
    if (!any(c("nn_prior_source_used", "mu01", "sigma01") %in% names(diag_all))) {
      result$status <- "invalid_object"
      result$warnings <- c(result$warnings, "nn_prior_diagnostics_missing_recognizable_two_shell_fields")
      return(result)
    }
  }

  result$ok <- TRUE
  result$status <- "valid"
  result
}

#' Ensure two-shell fits exist, reusing valid cached fits
#'
#' @param patients Named list of patient inputs accepted by `alfak()`.
#' @param patient_ids Patient identifiers aligned with `patients`.
#' @param two_shell_root Cache root. If `NULL`, `outdir/two_shell_base` is used.
#' @param outdir Cohort output directory where status tables are saved.
#' @param pm Mis-segregation probability for the upstream two-shell fit.
#' @param minobs MINIOBS threshold for the upstream two-shell fit.
#' @param ... Additional arguments passed to `alfak()` when a sample must be rerun.
#' @param sample_names Optional sample directory names aligned with `patient_ids`.
#' @param pm_tag Optional exact PM cache tag.
#' @param minobs_tag Optional exact MINIOBS cache tag.
#' @param sample_map Optional named patient-to-sample directory map.
#' @param reuse_two_shell Reuse valid cached fits.
#' @param rerun_missing_two_shell Rerun only missing fits.
#' @param rerun_corrupt_two_shell Back up and rerun only corrupt fits.
#' @param integrity_check Integrity mode used for existing and rerun fits.
#' @param base_nn_prior NN prior mode used for rerun base fits.
#' @param allow_incomplete_cohort If `TRUE`, record rerun failures instead of
#'   stopping immediately.
#' @return A data frame with before/after status and actions.
#' @export
ensure_two_shell_fits <- function(patients,
                                  patient_ids = names(patients),
                                  two_shell_root = NULL,
                                  outdir,
                                  pm,
                                  minobs,
                                  ...,
                                  sample_names = NULL,
                                  pm_tag = NULL,
                                  minobs_tag = NULL,
                                  sample_map = NULL,
                                  reuse_two_shell = TRUE,
                                  rerun_missing_two_shell = TRUE,
                                  rerun_corrupt_two_shell = TRUE,
                                  integrity_check = c("strict", "basic", "none"),
                                  base_nn_prior = "empirical_two_shell",
                                  allow_incomplete_cohort = FALSE) {
  integrity_check <- match.arg(integrity_check)
  if (is.null(two_shell_root)) {
    two_shell_root <- file.path(outdir, "two_shell_base")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(two_shell_root, recursive = TRUE, showWarnings = FALSE)

  patient_ids <- as.character(patient_ids)
  if (length(patients) != length(patient_ids)) {
    stop("`patients` and `patient_ids` must have the same length.", call. = FALSE)
  }
  names(patients) <- patient_ids

  resolved <- resolve_two_shell_fit_dirs(
    two_shell_root = two_shell_root,
    patient_ids = patient_ids,
    sample_names = sample_names,
    pm = pm,
    minobs = minobs,
    pm_tag = pm_tag,
    minobs_tag = minobs_tag,
    sample_map = sample_map
  )

  rows <- vector("list", nrow(resolved))
  for (i in seq_len(nrow(resolved))) {
    patient_id <- resolved$patient_id[i]
    fit_dir <- resolved$expected_fit_dir[i]
    before <- check_two_shell_fit_integrity(fit_dir, patient_id = patient_id, mode = integrity_check)
    action <- "none"
    backup_dir <- NA_character_
    error_message <- NA_character_
    after <- before

    needs_rerun <- FALSE
    if (isTRUE(reuse_two_shell) && isTRUE(before$ok)) {
      action <- "reused"
    } else {
      missing_dir <- identical(before$status, "missing_dir")
      if (missing_dir) {
        if (!isTRUE(rerun_missing_two_shell)) {
          stop(sprintf("Two-shell fit for patient `%s` is missing at `%s`.", patient_id, fit_dir), call. = FALSE)
        }
        action <- "rerun_missing"
        needs_rerun <- TRUE
      } else {
        if (!isTRUE(rerun_corrupt_two_shell)) {
          stop(
            sprintf("Two-shell fit for patient `%s` is not reusable at `%s` (status: %s).",
                    patient_id, fit_dir, before$status),
            call. = FALSE
          )
        }
        action <- if (isTRUE(reuse_two_shell)) "rerun_corrupt" else "rerun_reuse_disabled"
        needs_rerun <- TRUE
      }
    }

    if (isTRUE(needs_rerun)) {
      rerun_result <- tryCatch(
        {
          if (dir.exists(fit_dir)) {
            stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
            backup_dir <- file.path(dirname(fit_dir), paste0(basename(fit_dir), "__corrupt_", stamp))
            suffix <- 0L
            while (file.exists(backup_dir)) {
              suffix <- suffix + 1L
              backup_dir <- file.path(dirname(fit_dir), paste0(basename(fit_dir), "__corrupt_", stamp, "_", suffix))
            }
            if (!file.rename(fit_dir, backup_dir)) {
              stop(sprintf("Could not move corrupt two-shell directory `%s` to `%s`.", fit_dir, backup_dir), call. = FALSE)
            }
          }
          dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
          alfak(
            yi = patients[[patient_id]],
            outdir = fit_dir,
            minobs = minobs,
            pm = pm,
            nn_prior = base_nn_prior,
            ...
          )
          check_two_shell_fit_integrity(fit_dir, patient_id = patient_id, mode = integrity_check)
        },
        error = function(e) {
          error_message <<- conditionMessage(e)
          list(ok = FALSE, status = "rerun_failed", missing_files = character(0),
               unreadable_files = character(0), warnings = character(0),
               fit_dir = fit_dir, patient_id = patient_id)
        }
      )
      after <- rerun_result
      if (!isTRUE(after$ok) && !isTRUE(allow_incomplete_cohort)) {
        stop(
          sprintf("Two-shell rerun failed integrity checks for patient `%s` at `%s` (status: %s; error: %s).",
                  patient_id, fit_dir, after$status, ifelse(is.na(error_message), "", error_message)),
          call. = FALSE
        )
      }
    }

    rows[[i]] <- data.frame(
      patient_id = patient_id,
      sample_name = resolved$sample_name[i],
      pm_tag = resolved$pm_tag[i],
      minobs_tag = resolved$minobs_tag[i],
      expected_fit_dir = fit_dir,
      fit_dir = fit_dir,
      status_before = before$status,
      action = action,
      status_after = after$status,
      reused = identical(action, "reused"),
      rerun = isTRUE(needs_rerun),
      backup_dir = backup_dir,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }

  status <- do.call(rbind, rows)
  rownames(status) <- NULL
  saveRDS(status, file.path(outdir, "two_shell_fit_status.Rds"))
  utils::write.table(
    status,
    file = file.path(outdir, "two_shell_fit_status.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  status
}

#' Compute zero-neighbour informativeness
#'
#' @param expected_count_parent_like Expected child count under Delta = 0.
#' @param weak_threshold Boundary for weakly informative zeros.
#' @param informative_threshold Boundary for informative zeros.
#' @return A data frame with expected count, numeric score, and category.
#' @export
compute_zero_informativeness_score <- function(expected_count_parent_like,
                                              weak_threshold = 0.5,
                                              informative_threshold = 3.0) {
  expected <- as.numeric(expected_count_parent_like)
  validate_nonnegative_finite(weak_threshold, "weak_threshold")
  validate_positive_finite(informative_threshold, "informative_threshold")
  if (informative_threshold <= weak_threshold) {
    stop("`informative_threshold` must be greater than `weak_threshold`.", call. = FALSE)
  }
  category <- ifelse(
    !is.finite(expected), "unknown",
    ifelse(expected < weak_threshold, "uninformative_zero",
           ifelse(expected < informative_threshold, "weakly_informative_zero", "informative_zero"))
  )
  score <- pmin(1, pmax(0, expected / informative_threshold))
  score[!is.finite(score)] <- NA_real_
  data.frame(
    expected_count_parent_like = expected,
    zero_informativeness_score = score,
    zero_informativeness_category = category,
    stringsAsFactors = FALSE
  )
}

cohort_transition_parse_pair <- function(parent_karyotype, child_karyotype) {
  parent_vec <- as.numeric(parse_karyotype_ids(parent_karyotype)[1, ])
  child_vec <- as.numeric(parse_karyotype_ids(child_karyotype)[1, ])
  if (length(parent_vec) != length(child_vec)) {
    stop("Parent and child karyotypes must have the same dimensionality.", call. = FALSE)
  }
  diff_vec <- child_vec - parent_vec
  changed <- which(diff_vec != 0)
  transition_chr <- if (length(changed) == 1L) changed else NA_integer_
  transition_size <- sum(abs(diff_vec))
  transition_direction <- if (length(changed) == 1L && diff_vec[changed] > 0) {
    "gain"
  } else if (length(changed) == 1L && diff_vec[changed] < 0) {
    "loss"
  } else {
    "complex"
  }
  parent_total_cn <- sum(parent_vec)
  child_total_cn <- sum(child_vec)
  parent_burden <- sum(abs(parent_vec - 2))
  child_burden <- sum(abs(child_vec - 2))
  burden_label <- if (parent_burden <= 1) "low" else "high"
  chr_label <- if (is.na(transition_chr)) "complex" else paste0("chr", transition_chr)
  exact_label <- paste0(parent_karyotype, ">", child_karyotype)
  list(
    transition_chr = transition_chr,
    transition_direction = transition_direction,
    transition_size = transition_size,
    parent_total_cn = parent_total_cn,
    child_total_cn = child_total_cn,
    parent_burden = parent_burden,
    child_burden = child_burden,
    group_gain_loss = transition_direction,
    group_gain_loss_chr = paste(transition_direction, chr_label, sep = "_"),
    group_gain_loss_chr_burden = paste(transition_direction, chr_label, paste0("burden_", burden_label), sep = "_"),
    group_exact_event = exact_label
  )
}

cohort_transition_group_column <- function(grouping) {
  grouping <- match.arg(grouping, c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"))
  paste0("group_", grouping)
}

cohort_transition_assign_groups <- function(records, grouping) {
  group_col <- cohort_transition_group_column(grouping)
  if (!group_col %in% names(records)) {
    stop(sprintf("Transition records are missing `%s`.", group_col), call. = FALSE)
  }
  records$transition_group <- records[[group_col]]
  records
}

cohort_transition_delta_se <- function(delta_values, fallback = ALFAK_NN_PRIOR_SD_FLOOR) {
  finite <- delta_values[is.finite(delta_values)]
  if (length(finite) < 2L) {
    return(fallback)
  }
  se <- stats::sd(finite) / sqrt(length(finite))
  if (!is.finite(se) || se <= 0) fallback else max(se, fallback)
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

#' Extract cohort transition records from two-shell fits
#'
#' @param fit_dirs Character vector of two-shell fit directories.
#' @param patient_ids Patient IDs aligned with `fit_dirs`.
#' @param pm Mis-segregation probability used to reconstruct NN parent paths.
#' @param grouping Transition grouping mode.
#' @param cohort_transition_use_zero Whether informative zero records are kept.
#' @param cohort_transition_zero_min_expected_count Minimum parent-like expected
#'   count for a zero child to contribute to the cohort prior.
#' @param cohort_transition_zero_min_exposure Optional explicit projected
#'   exposure threshold for retaining zero children.
#' @param ... Reserved for future extraction controls.
#' @return A data frame with one row per usable transition path and bootstrap.
#' @export
extract_cohort_transition_records <- function(fit_dirs,
                                              patient_ids = names(fit_dirs),
                                              pm = 0.00005,
                                              grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                              cohort_transition_use_zero = TRUE,
                                              cohort_transition_zero_min_expected_count = 3.0,
                                              cohort_transition_zero_min_exposure = NULL,
                                              ...) {
  grouping <- match.arg(grouping)
  validate_probability(pm, "pm", upper_inclusive = TRUE)
  validate_nonnegative_finite(cohort_transition_zero_min_expected_count, "cohort_transition_zero_min_expected_count")
  if (!is.null(cohort_transition_zero_min_exposure)) {
    validate_nonnegative_finite(cohort_transition_zero_min_exposure, "cohort_transition_zero_min_exposure")
  }
  fit_dirs <- as.character(fit_dirs)
  if (is.null(patient_ids)) {
    patient_ids <- basename(fit_dirs)
  }
  patient_ids <- as.character(patient_ids)
  if (length(fit_dirs) != length(patient_ids)) {
    stop("`fit_dirs` and `patient_ids` must have the same length.", call. = FALSE)
  }

  all_rows <- list()
  row_idx <- 0L
  for (pidx in seq_along(fit_dirs)) {
    fit_dir <- fit_dirs[pidx]
    patient_id <- patient_ids[pidx]
    boot <- readRDS(file.path(fit_dir, "bootstrap_res.Rds"))
    final_fitness <- boot$final_fitness
    nn_fitness <- boot$nn_fitness
    if (!is.matrix(final_fitness) || !is.matrix(nn_fitness) || !ncol(nn_fitness)) {
      next
    }
    fq <- colnames(final_fitness)
    nn_info <- gen_nn_info(fq, pm)
    if (length(nn_info)) {
      names(nn_info) <- vapply(nn_info, function(x) x$ni, character(1))
    }
    nn_info <- nn_info[colnames(nn_fitness)]
    nn_info <- nn_info[!vapply(nn_info, is.null, logical(1))]

    diag_path <- file.path(fit_dir, "nn_prior_diagnostics.Rds")
    node_diag <- data.frame()
    if (file.exists(diag_path)) {
      diag <- tryCatch(readRDS(diag_path), error = function(e) NULL)
      if (is.list(diag) && is.data.frame(diag$node)) {
        node_diag <- diag$node
      } else if (is.data.frame(diag)) {
        node_diag <- diag
      }
    }

    delta_se_lookup <- list()
    for (child in intersect(names(nn_info), colnames(nn_fitness))) {
      item <- nn_info[[child]]
      for (parent in item$nj) {
        if (!parent %in% colnames(final_fitness)) next
        deltas <- nn_fitness[, child] - final_fitness[, parent]
        delta_se_lookup[[paste(parent, child, sep = "\r")]] <- cohort_transition_delta_se(deltas)
      }
    }

    for (b in seq_len(nrow(nn_fitness))) {
      for (child in intersect(names(nn_info), colnames(nn_fitness))) {
        item <- nn_info[[child]]
        valid_parents <- item$nj[item$nj %in% colnames(final_fitness)]
        if (!length(valid_parents)) next
        parent_weights <- item$pij[match(valid_parents, item$nj)]
        if (!all(is.finite(parent_weights)) || sum(parent_weights) <= 0) {
          parent_weights <- rep(1, length(valid_parents))
        }
        path_responsibility <- parent_weights / sum(parent_weights)

        node_row <- node_diag[FALSE, , drop = FALSE]
        if (nrow(node_diag) &&
            all(c("replicate_id", "karyotype") %in% names(node_diag))) {
          node_row <- node_diag[node_diag$replicate_id == b & node_diag$karyotype == child, , drop = FALSE]
        } else if (nrow(node_diag) && "karyotype" %in% names(node_diag)) {
          node_row <- node_diag[node_diag$karyotype == child, , drop = FALSE]
        }
        if (nrow(node_row) > 1L) node_row <- node_row[1L, , drop = FALSE]
        get_node_value <- function(name, default) {
          if (nrow(node_row) && name %in% names(node_row)) node_row[[name]][1] else default
        }

        child_observed_count <- as.numeric(get_node_value("direct_observed_count", NA_real_))
        projected_exposure <- as.numeric(get_node_value("projected_exposure", NA_real_))
        expected_parent_like <- projected_exposure
        zero_info <- compute_zero_informativeness_score(expected_parent_like)
        child_is_zero <- is.finite(child_observed_count) && child_observed_count <= 0
        zero_retained <- isTRUE(cohort_transition_use_zero) &&
          isTRUE(child_is_zero) &&
          is.finite(expected_parent_like) &&
          expected_parent_like >= cohort_transition_zero_min_expected_count
        if (!is.null(cohort_transition_zero_min_exposure)) {
          zero_retained <- isTRUE(zero_retained) &&
            is.finite(projected_exposure) &&
            projected_exposure >= cohort_transition_zero_min_exposure
        }
        source_type <- if (isTRUE(child_is_zero)) {
          if (isTRUE(zero_retained)) "informative_zero" else "low_exposure_zero"
        } else if (is.finite(child_observed_count) && child_observed_count > 0) {
          "observed"
        } else {
          "observed"
        }
        if (identical(source_type, "low_exposure_zero")) {
          next
        }

        for (parent_idx in seq_along(valid_parents)) {
          parent <- valid_parents[parent_idx]
          parent_fit <- final_fitness[b, parent]
          child_fit <- nn_fitness[b, child]
          if (!is.finite(parent_fit) || !is.finite(child_fit)) next
          parsed <- cohort_transition_parse_pair(parent, child)
          row_idx <- row_idx + 1L
          all_rows[[row_idx]] <- data.frame(
            patient_id = patient_id,
            parent_karyotype = parent,
            child_karyotype = child,
            transition_chr = parsed$transition_chr,
            transition_direction = parsed$transition_direction,
            transition_size = parsed$transition_size,
            group_gain_loss = parsed$group_gain_loss,
            group_gain_loss_chr = parsed$group_gain_loss_chr,
            group_gain_loss_chr_burden = parsed$group_gain_loss_chr_burden,
            group_exact_event = parsed$group_exact_event,
            transition_group = parsed[[cohort_transition_group_column(grouping)]],
            parent_total_cn = parsed$parent_total_cn,
            child_total_cn = parsed$child_total_cn,
            parent_burden = parsed$parent_burden,
            child_burden = parsed$child_burden,
            parent_fitness = parent_fit,
            child_fitness_two_shell = child_fit,
            delta_hat = child_fit - parent_fit,
            delta_se = delta_se_lookup[[paste(parent, child, sep = "\r")]],
            child_observed_count = child_observed_count,
            child_is_zero = isTRUE(child_is_zero),
            projected_exposure = projected_exposure,
            expected_count_parent_like = zero_info$expected_count_parent_like,
            zero_informativeness_score = zero_info$zero_informativeness_score,
            zero_informativeness_category = zero_info$zero_informativeness_category,
            boundary_flag = as.logical(get_node_value("objective_boundary_flag", FALSE)),
            prior_dominated_flag = as.logical(get_node_value("prior_dominated_flag", FALSE)),
            two_shell_used = TRUE,
            two_shell_outward_weight = as.numeric(get_node_value("outward_weight_sum", NA_real_)),
            path_responsibility = path_responsibility[parent_idx],
            replicate_id = as.integer(b),
            bootstrap_id = as.integer(b),
            source_type = source_type,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(all_rows)) {
    empty <- data.frame(
      patient_id = character(0),
      parent_karyotype = character(0),
      child_karyotype = character(0),
      transition_group = character(0),
      delta_hat = numeric(0),
      delta_se = numeric(0),
      source_type = character(0),
      path_responsibility = numeric(0),
      stringsAsFactors = FALSE
    )
    return(empty)
  }
  out <- do.call(rbind, all_rows)
  rownames(out) <- NULL
  cohort_transition_assign_groups(out, grouping)
}

cohort_transition_record_weights <- function(records,
                                             zero_weight_cap_ratio = 1.0) {
  w <- records$path_responsibility
  w[!is.finite(w) | w < 0] <- 0
  observed <- records$source_type != "informative_zero"
  zero <- records$source_type == "informative_zero"
  if (any(zero)) {
    zero_score <- records$zero_informativeness_score
    zero_score[!is.finite(zero_score) | zero_score < 0] <- 0
    w[zero] <- w[zero] * pmin(1, zero_score[zero])
  }
  patient_id <- as.character(records$patient_id)
  patient_weight_sum <- ave(w, patient_id, FUN = function(x) sum(x, na.rm = TRUE))
  can_normalize <- is.finite(patient_weight_sum) & patient_weight_sum > 0
  w[can_normalize] <- w[can_normalize] / patient_weight_sum[can_normalize]
  w[!can_normalize] <- 0
  if (any(zero)) {
    obs_sum <- sum(w[observed], na.rm = TRUE)
    zero_sum <- sum(w[zero], na.rm = TRUE)
    cap <- zero_weight_cap_ratio * max(obs_sum, 1)
    if (is.finite(cap) && cap >= 0 && zero_sum > cap && zero_sum > 0) {
      w[zero] <- w[zero] * (cap / zero_sum)
    }
  }
  w
}

cohort_transition_empty_filter_diagnostics <- function() {
  reasons <- c(
    "low_exposure_zero",
    "prior_dominated",
    "boundary_dominated",
    "missing_delta_se",
    "large_delta_se",
    "low_path_responsibility",
    "missing_transition_group",
    "non_one_step_or_complex",
    "non_finite_delta",
    "zero_not_informative"
  )
  stats::setNames(as.list(rep(0L, length(reasons))), reasons)
}

#' Filter cohort transition records before v2 prior learning
#'
#' The v2 prior treats observed transition effects and zero-censoring evidence
#' separately. Bootstrap/path records that are prior dominated, boundary
#' dominated, low-responsibility, or too uncertain are retained only in
#' diagnostics unless the corresponding control explicitly permits them.
#'
#' @param records Raw transition records from `extract_cohort_transition_records()`.
#' @param cohort_transition_use_prior_dominated_records Whether prior-dominated
#'   records may be used as observed transition labels.
#' @param cohort_transition_use_boundary_records Whether boundary records may be
#'   used as observed transition labels.
#' @param cohort_transition_max_delta_se Optional maximum allowed delta SE.
#' @param cohort_transition_max_delta_se_quantile Quantile used to derive a
#'   maximum delta SE when `cohort_transition_max_delta_se` is `NULL`.
#' @param cohort_transition_min_path_responsibility Minimum path responsibility.
#' @param cohort_transition_min_observed_count Minimum child count for observed
#'   transition labels.
#' @param cohort_transition_zero_min_expected_count Minimum parent-like expected
#'   count for zero-censoring records.
#' @param cohort_transition_zero_as_censoring_only Keep zero records as censoring
#'   evidence only; they are never observed delta labels.
#' @return A list with kept records, excluded records, and diagnostics.
#' @export
filter_cohort_transition_records <- function(records,
                                             cohort_transition_use_prior_dominated_records = FALSE,
                                             cohort_transition_use_boundary_records = FALSE,
                                             cohort_transition_max_delta_se = NULL,
                                             cohort_transition_max_delta_se_quantile = 0.75,
                                             cohort_transition_min_path_responsibility = 0.05,
                                             cohort_transition_min_observed_count = 1L,
                                             cohort_transition_zero_min_expected_count = 3.0,
                                             cohort_transition_zero_as_censoring_only = TRUE,
                                             ...) {
  validate_scalar_logical(cohort_transition_use_prior_dominated_records, "cohort_transition_use_prior_dominated_records")
  validate_scalar_logical(cohort_transition_use_boundary_records, "cohort_transition_use_boundary_records")
  validate_probability(cohort_transition_max_delta_se_quantile, "cohort_transition_max_delta_se_quantile", upper_inclusive = TRUE)
  validate_nonnegative_finite(cohort_transition_min_path_responsibility, "cohort_transition_min_path_responsibility")
  validate_nonnegative_integer(cohort_transition_min_observed_count, "cohort_transition_min_observed_count")
  validate_nonnegative_finite(cohort_transition_zero_min_expected_count, "cohort_transition_zero_min_expected_count")
  validate_scalar_logical(cohort_transition_zero_as_censoring_only, "cohort_transition_zero_as_censoring_only")
  if (!is.null(cohort_transition_max_delta_se)) {
    validate_positive_finite(cohort_transition_max_delta_se, "cohort_transition_max_delta_se")
  }
  if (!is.data.frame(records)) {
    stop("`records` must be a data frame.", call. = FALSE)
  }
  if (!nrow(records)) {
    empty_diag <- cohort_transition_empty_filter_diagnostics()
    return(list(
      kept_records = records,
      excluded_records = records,
      diagnostics = c(list(n_input_records = 0L, n_kept_records = 0L, n_excluded_records = 0L), empty_diag)
    ))
  }
  records <- records
  for (col in c("source_type", "transition_group", "delta_hat", "delta_se",
                "path_responsibility", "child_observed_count", "expected_count_parent_like",
                "child_is_zero", "prior_dominated_flag", "boundary_flag",
                "transition_size", "zero_informativeness_score")) {
    if (!col %in% names(records)) {
      records[[col]] <- NA
    }
  }
  child_is_zero_vec <- as.logical(records$child_is_zero)
  child_is_zero_vec[is.na(child_is_zero_vec)] <- FALSE
  observed_source <- records$source_type %in% c("observed", "fq_transition")
  zero_source <- records$source_type == "informative_zero" |
    (isTRUE(cohort_transition_zero_as_censoring_only) & child_is_zero_vec)
  finite_observed_se <- records$delta_se[observed_source & is.finite(records$delta_se) & records$delta_se > 0]
  derived_max_delta_se <- cohort_transition_max_delta_se
  if (is.null(derived_max_delta_se) && length(finite_observed_se)) {
    derived_max_delta_se <- as.numeric(stats::quantile(
      finite_observed_se,
      probs = cohort_transition_max_delta_se_quantile,
      na.rm = TRUE,
      names = FALSE,
      type = 8
    ))
  }
  if (!is.finite(derived_max_delta_se) || derived_max_delta_se <= 0) {
    derived_max_delta_se <- Inf
  }

  reason <- rep(NA_character_, nrow(records))
  missing_group <- is.na(records$transition_group) | !nzchar(as.character(records$transition_group))
  reason[is.na(reason) & missing_group] <- "missing_transition_group"
  complex_event <- is.finite(records$transition_size) & abs(records$transition_size) != 1
  reason[is.na(reason) & complex_event] <- "non_one_step_or_complex"
  low_path <- !is.finite(records$path_responsibility) |
    records$path_responsibility < cohort_transition_min_path_responsibility
  reason[is.na(reason) & low_path] <- "low_path_responsibility"

  observed_label <- observed_source & !child_is_zero_vec
  if (!isTRUE(cohort_transition_use_prior_dominated_records)) {
    prior_dominated <- as.logical(records$prior_dominated_flag)
    prior_dominated[is.na(prior_dominated)] <- FALSE
    reason[is.na(reason) & observed_label & prior_dominated] <- "prior_dominated"
  }
  if (!isTRUE(cohort_transition_use_boundary_records)) {
    boundary_dominated <- as.logical(records$boundary_flag)
    boundary_dominated[is.na(boundary_dominated)] <- FALSE
    reason[is.na(reason) & observed_label & boundary_dominated] <- "boundary_dominated"
  }
  low_count <- !is.finite(records$child_observed_count) |
    records$child_observed_count < cohort_transition_min_observed_count
  reason[is.na(reason) & observed_label & low_count] <- "zero_not_informative"
  missing_se <- !is.finite(records$delta_se) | records$delta_se <= 0
  reason[is.na(reason) & observed_label & missing_se] <- "missing_delta_se"
  reason[is.na(reason) & observed_label & is.finite(records$delta_se) &
           records$delta_se > derived_max_delta_se] <- "large_delta_se"
  reason[is.na(reason) & observed_label & !is.finite(records$delta_hat)] <- "non_finite_delta"

  zero_label <- zero_source | child_is_zero_vec
  zero_low_exposure <- !is.finite(records$expected_count_parent_like) |
    records$expected_count_parent_like < cohort_transition_zero_min_expected_count
  reason[is.na(reason) & zero_label & zero_low_exposure] <- "low_exposure_zero"
  reason[is.na(reason) & zero_label & !zero_low_exposure] <- NA_character_

  recognized <- observed_label | zero_label
  reason[is.na(reason) & !recognized] <- "zero_not_informative"
  keep <- is.na(reason)
  kept <- records[keep, , drop = FALSE]
  excluded <- records[!keep, , drop = FALSE]
  excluded$exclusion_reason <- reason[!keep]
  kept$cohort_transition_evidence_type <- ifelse(
    kept$source_type == "informative_zero" | as.logical(kept$child_is_zero),
    "zero_censoring_evidence",
    "observed_delta_evidence"
  )
  if (nrow(kept) && isTRUE(cohort_transition_zero_as_censoring_only)) {
    zero_idx <- kept$cohort_transition_evidence_type == "zero_censoring_evidence"
    kept$delta_hat[zero_idx] <- NA_real_
  }

  reason_counts <- table(factor(reason[!keep], names(cohort_transition_empty_filter_diagnostics())))
  diagnostics <- c(
    list(
      n_input_records = as.integer(nrow(records)),
      n_kept_records = as.integer(nrow(kept)),
      n_excluded_records = as.integer(nrow(excluded)),
      max_delta_se_used = derived_max_delta_se
    ),
    as.list(as.integer(reason_counts))
  )
  names(diagnostics)[seq_along(reason_counts) + 4L] <- names(reason_counts)
  list(kept_records = kept, excluded_records = excluded, diagnostics = diagnostics)
}

cohort_transition_weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok) || sum(w[ok]) <= 0) {
    return(NA_real_)
  }
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  x[which(cumsum(w) >= 0.5)[1L]]
}

#' Aggregate raw cohort transition records to patient-level evidence
#'
#' Bootstrap/path-level rows are useful diagnostics but are not independent
#' cohort patients. This helper collapses them to patient-level observed-delta
#' and zero-censoring summaries before v2 prior fitting.
#'
#' @param records Filtered or raw transition records.
#' @param grouping Transition grouping mode used for `transition_group`.
#' @param cohort_transition_sd_floor Floor used in inverse-variance weights.
#' @param cohort_transition_zero_weight_cap_ratio Cap on zero-censoring weight.
#' @return A data frame with patient-level group summaries.
#' @export
aggregate_cohort_transition_records_by_patient <- function(records,
                                                           grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                                           cohort_transition_sd_floor = 0.05,
                                                           cohort_transition_zero_weight_cap_ratio = 0.25,
                                                           ...) {
  grouping <- match.arg(grouping)
  validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
  validate_nonnegative_finite(cohort_transition_zero_weight_cap_ratio, "cohort_transition_zero_weight_cap_ratio")
  if (!is.data.frame(records) || !nrow(records)) {
    return(data.frame())
  }
  records <- cohort_transition_assign_groups(records, grouping)
  for (col in c("cohort_transition_evidence_type", "source_type", "patient_id", "parent_karyotype",
                "child_karyotype", "delta_hat", "delta_se", "path_responsibility",
                "expected_count_parent_like", "zero_informativeness_score",
                "replicate_id", "bootstrap_id")) {
    if (!col %in% names(records)) records[[col]] <- NA
  }
  if (!"cohort_transition_evidence_type" %in% names(records) ||
      all(is.na(records$cohort_transition_evidence_type))) {
    records$cohort_transition_evidence_type <- ifelse(
      records$source_type == "informative_zero" | as.logical(records$child_is_zero),
      "zero_censoring_evidence",
      "observed_delta_evidence"
    )
  }
  key_cols <- c(
    "patient_id",
    "group_gain_loss",
    "group_gain_loss_chr",
    "group_gain_loss_chr_burden",
    "group_exact_event",
    "transition_group",
    "cohort_transition_evidence_type"
  )
  split_key <- interaction(records[key_cols], drop = TRUE, sep = "\r")
  rows <- lapply(split(records, split_key), function(df) {
    evidence_type <- df$cohort_transition_evidence_type[1]
    base <- df[1L, intersect(key_cols, names(df)), drop = FALSE]
    if (identical(evidence_type, "observed_delta_evidence")) {
      se <- pmax(as.numeric(df$delta_se), cohort_transition_sd_floor)
      w <- as.numeric(df$path_responsibility) / (se^2)
      w[!is.finite(w) | w < 0] <- 0
      if (sum(w) <= 0) {
        w <- rep(1, nrow(df))
      }
      delta <- as.numeric(df$delta_hat)
      ok <- is.finite(delta) & is.finite(w) & w > 0
      delta_mean <- if (any(ok)) stats::weighted.mean(delta[ok], w[ok]) else NA_real_
      delta_median <- cohort_transition_weighted_median(delta, w)
      delta_se <- if (sum(w[ok]) > 0) sqrt(1 / sum(w[ok])) else NA_real_
      out <- data.frame(
        base,
        delta_patient_mean = delta_mean,
        delta_patient_median = delta_median,
        delta_patient_se = delta_se,
        n_raw_records = nrow(df),
        n_bootstrap_records = length(unique(df$bootstrap_id[!is.na(df$bootstrap_id)])),
        n_unique_children = length(unique(df$child_karyotype)),
        n_unique_parents = length(unique(df$parent_karyotype)),
        total_path_responsibility = sum(df$path_responsibility, na.rm = TRUE),
        observed_weight = sum(w[ok], na.rm = TRUE),
        n_zero_records = 0L,
        total_zero_weight = 0,
        max_expected_count_parent_like = NA_real_,
        mean_expected_count_parent_like = NA_real_,
        zero_informativeness_score = NA_real_,
        zero_censoring_weight = 0,
        stringsAsFactors = FALSE
      )
    } else {
      zero_w <- as.numeric(df$path_responsibility) * pmin(1, pmax(0, as.numeric(df$zero_informativeness_score)))
      zero_w[!is.finite(zero_w) | zero_w < 0] <- 0
      zero_cap <- cohort_transition_zero_weight_cap_ratio * max(1, length(unique(df$patient_id)))
      if (sum(zero_w) > zero_cap && zero_cap >= 0) {
        zero_w <- zero_w * (zero_cap / sum(zero_w))
      }
      out <- data.frame(
        base,
        delta_patient_mean = NA_real_,
        delta_patient_median = NA_real_,
        delta_patient_se = NA_real_,
        n_raw_records = nrow(df),
        n_bootstrap_records = length(unique(df$bootstrap_id[!is.na(df$bootstrap_id)])),
        n_unique_children = length(unique(df$child_karyotype)),
        n_unique_parents = length(unique(df$parent_karyotype)),
        total_path_responsibility = sum(df$path_responsibility, na.rm = TRUE),
        observed_weight = 0,
        n_zero_records = nrow(df),
        total_zero_weight = sum(zero_w, na.rm = TRUE),
        max_expected_count_parent_like = max(df$expected_count_parent_like, na.rm = TRUE),
        mean_expected_count_parent_like = mean(df$expected_count_parent_like, na.rm = TRUE),
        zero_informativeness_score = mean(df$zero_informativeness_score, na.rm = TRUE),
        zero_censoring_weight = sum(zero_w, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Compute v2 class-specific borrowing controls
#'
#' @param effect_class Transition group class.
#' @param cohort_transition_lambda_consistent_deleterious,cohort_transition_lambda_consistent_neutral,cohort_transition_lambda_consistent_beneficial,cohort_transition_lambda_context_dependent,cohort_transition_lambda_high_variable,cohort_transition_lambda_sparse_unknown,cohort_transition_lambda_global_fallback
#'   Class-specific borrowing multipliers.
#' @param cohort_transition_sd_multiplier_consistent_deleterious,cohort_transition_sd_multiplier_consistent_neutral,cohort_transition_sd_multiplier_consistent_beneficial,cohort_transition_sd_multiplier_context_dependent,cohort_transition_sd_multiplier_high_variable,cohort_transition_sd_multiplier_sparse_unknown,cohort_transition_sd_multiplier_global_fallback
#'   Class-specific prior SD multipliers.
#' @return A list with lambda, SD multiplier, and use flags.
#' @export
compute_class_specific_borrowing <- function(effect_class,
                                             cohort_transition_lambda_consistent_deleterious = 0.50,
                                             cohort_transition_lambda_consistent_neutral = 0.25,
                                             cohort_transition_lambda_consistent_beneficial = 0.15,
                                             cohort_transition_lambda_context_dependent = 0.30,
                                             cohort_transition_lambda_high_variable = 0.00,
                                             cohort_transition_lambda_sparse_unknown = 0.00,
                                             cohort_transition_lambda_global_fallback = 0.05,
                                             cohort_transition_sd_multiplier_consistent_deleterious = 1.0,
                                             cohort_transition_sd_multiplier_consistent_neutral = 1.5,
                                             cohort_transition_sd_multiplier_consistent_beneficial = 2.0,
                                             cohort_transition_sd_multiplier_context_dependent = 1.5,
                                             cohort_transition_sd_multiplier_high_variable = 4.0,
                                             cohort_transition_sd_multiplier_sparse_unknown = 4.0,
                                             cohort_transition_sd_multiplier_global_fallback = 4.0) {
  lambda <- switch(
    effect_class,
    consistent_deleterious = cohort_transition_lambda_consistent_deleterious,
    consistent_neutral = cohort_transition_lambda_consistent_neutral,
    consistent_beneficial = cohort_transition_lambda_consistent_beneficial,
    context_dependent = cohort_transition_lambda_context_dependent,
    high_variable = cohort_transition_lambda_high_variable,
    sparse_unknown = cohort_transition_lambda_sparse_unknown,
    global_fallback = cohort_transition_lambda_global_fallback,
    cohort_transition_lambda_sparse_unknown
  )
  sd_multiplier <- switch(
    effect_class,
    consistent_deleterious = cohort_transition_sd_multiplier_consistent_deleterious,
    consistent_neutral = cohort_transition_sd_multiplier_consistent_neutral,
    consistent_beneficial = cohort_transition_sd_multiplier_consistent_beneficial,
    context_dependent = cohort_transition_sd_multiplier_context_dependent,
    high_variable = cohort_transition_sd_multiplier_high_variable,
    sparse_unknown = cohort_transition_sd_multiplier_sparse_unknown,
    global_fallback = cohort_transition_sd_multiplier_global_fallback,
    cohort_transition_sd_multiplier_sparse_unknown
  )
  list(
    cohort_lambda = lambda,
    sd_multiplier = sd_multiplier,
    use_for_zero = effect_class %in% c("consistent_deleterious", "consistent_neutral", "consistent_beneficial", "context_dependent"),
    use_for_observed = FALSE,
    use_for_low_information = effect_class %in% c("consistent_deleterious", "consistent_neutral", "context_dependent")
  )
}

#' Compute a v2 transition group class
#'
#' @param k,eff_k,effective_observed Observed patient count, effective patient
#'   count, and effective observed evidence.
#' @param mu,se_mu Weighted group mean and standard error.
#' @param tau,i2,sign_consistency Heterogeneity metrics.
#' @param cohort_transition_min_patients_consistent,cohort_transition_min_effective_patients,cohort_transition_min_effective_observed
#'   Minimum support thresholds.
#' @param cohort_transition_effect_threshold,cohort_transition_sign_consistency_threshold,cohort_transition_high_heterogeneity_i2,cohort_transition_high_between_patient_sd
#'   Classification thresholds.
#' @return A list with effect/heterogeneity class and normal-approximation
#'   effect probabilities.
#' @export
compute_transition_group_class <- function(k,
                                           eff_k,
                                           effective_observed,
                                           mu,
                                           se_mu,
                                           tau,
                                           i2,
                                           sign_consistency,
                                           cohort_transition_min_patients_consistent = 3L,
                                           cohort_transition_min_effective_patients = 3,
                                           cohort_transition_min_effective_observed = 3,
                                           cohort_transition_effect_threshold = 0.02,
                                           cohort_transition_sign_consistency_threshold = 0.75,
                                           cohort_transition_high_heterogeneity_i2 = 0.50,
                                           cohort_transition_high_between_patient_sd = 0.10) {
  if (!is.finite(se_mu) || se_mu <= 0) se_mu <- Inf
  p_beneficial <- if (is.finite(se_mu)) 1 - stats::pnorm(cohort_transition_effect_threshold, mean = mu, sd = se_mu) else 0
  p_deleterious <- if (is.finite(se_mu)) stats::pnorm(-cohort_transition_effect_threshold, mean = mu, sd = se_mu) else 0
  p_neutral <- if (is.finite(se_mu)) {
    stats::pnorm(cohort_transition_effect_threshold, mean = mu, sd = se_mu) -
      stats::pnorm(-cohort_transition_effect_threshold, mean = mu, sd = se_mu)
  } else {
    0
  }
  sparse <- k < cohort_transition_min_patients_consistent ||
    eff_k < cohort_transition_min_effective_patients ||
    effective_observed < cohort_transition_min_effective_observed
  high_variable <- !sparse && (
    (is.finite(i2) && i2 >= cohort_transition_high_heterogeneity_i2) ||
      (is.finite(tau) && tau >= cohort_transition_high_between_patient_sd) ||
      (is.finite(sign_consistency) && sign_consistency < cohort_transition_sign_consistency_threshold)
  )
  if (sparse) {
    effect_class <- "sparse_unknown"
  } else if (high_variable) {
    effect_class <- "high_variable"
  } else if (p_deleterious >= 0.8 && sign_consistency >= cohort_transition_sign_consistency_threshold) {
    effect_class <- "consistent_deleterious"
  } else if (p_beneficial >= 0.8 && sign_consistency >= cohort_transition_sign_consistency_threshold) {
    effect_class <- "consistent_beneficial"
  } else if (p_neutral >= 0.6 && tau < cohort_transition_high_between_patient_sd) {
    effect_class <- "consistent_neutral"
  } else {
    effect_class <- "sparse_unknown"
  }
  heterogeneity_class <- if (sparse) {
    "sparse_unknown"
  } else if (high_variable) {
    "high_variable"
  } else if (is.finite(i2) && i2 < 0.25 && is.finite(tau) && tau < cohort_transition_high_between_patient_sd / 2) {
    "low_heterogeneity"
  } else {
    "moderate_heterogeneity"
  }
  list(
    effect_class = effect_class,
    heterogeneity_class = heterogeneity_class,
    p_beneficial = p_beneficial,
    p_deleterious = p_deleterious,
    p_neutral = p_neutral
  )
}

#' Compute heterogeneity metrics for one transition group
#'
#' @param patient_group_summaries Patient-level summaries from
#'   `aggregate_cohort_transition_records_by_patient()`.
#' @param group_name Group label.
#' @param group_level Grouping level.
#' @param group_col Column containing `group_name`.
#' @return A one-row data frame of heterogeneity metrics.
#' @export
compute_transition_group_heterogeneity <- function(patient_group_summaries,
                                                   group_name,
                                                   group_level,
                                                   group_col,
                                                   cohort_transition_sd_floor = 0.05,
                                                   cohort_transition_effect_threshold = 0.02,
                                                   ...) {
  obs <- patient_group_summaries[
    patient_group_summaries$cohort_transition_evidence_type == "observed_delta_evidence" &
      patient_group_summaries[[group_col]] == group_name &
      is.finite(patient_group_summaries$delta_patient_mean),
    ,
    drop = FALSE
  ]
  zero <- patient_group_summaries[
    patient_group_summaries$cohort_transition_evidence_type == "zero_censoring_evidence" &
      patient_group_summaries[[group_col]] == group_name,
    ,
    drop = FALSE
  ]
  k <- length(unique(obs$patient_id))
  z_k <- length(unique(zero$patient_id))
  if (!nrow(obs)) {
    return(data.frame(
      transition_group = group_name,
      group = group_name,
      group_level = group_level,
      fallback_group = NA_character_,
      n_patients_observed = 0L,
      n_patients_zero = z_k,
      effective_patients_observed = 0,
      effective_patients_total = z_k,
      weighted_mean_delta = 0,
      weighted_se_delta = Inf,
      weighted_sd_delta = NA_real_,
      between_patient_sd = NA_real_,
      tau_between_patient = NA_real_,
      i2_heterogeneity = NA_real_,
      sign_consistency = NA_real_,
      n_positive_patients = 0L,
      n_negative_patients = 0L,
      n_near_zero_patients = 0L,
      p_beneficial = 0,
      p_deleterious = 0,
      p_neutral = 0,
      stringsAsFactors = FALSE
    ))
  }
  delta <- obs$delta_patient_mean
  se <- pmax(obs$delta_patient_se, cohort_transition_sd_floor)
  w <- obs$observed_weight
  w[!is.finite(w) | w <= 0] <- 1 / (se[!is.finite(w) | w <= 0]^2)
  w[!is.finite(w) | w <= 0] <- 1
  mu <- stats::weighted.mean(delta, w)
  se_mu <- sqrt(1 / sum(w))
  weighted_var <- sum(w * (delta - mu)^2) / max(sum(w), .Machine$double.eps)
  weighted_sd <- sqrt(max(0, weighted_var))
  mean_se2 <- stats::weighted.mean(se^2, w)
  tau <- sqrt(max(0, weighted_var - mean_se2))
  i2 <- tau^2 / (tau^2 + mean_se2)
  if (!is.finite(i2)) i2 <- 0
  eff_k <- sum(w)^2 / sum(w^2)
  mu_sign <- if (abs(mu) <= cohort_transition_effect_threshold) 0 else sign(mu)
  signs <- ifelse(abs(delta) <= cohort_transition_effect_threshold, 0, sign(delta))
  sign_consistency <- if (mu_sign == 0) {
    mean(signs == 0)
  } else {
    mean(signs == mu_sign)
  }
  data.frame(
    transition_group = group_name,
    group = group_name,
    group_level = group_level,
    fallback_group = NA_character_,
    n_patients_observed = k,
    n_patients_zero = z_k,
    effective_patients_observed = eff_k,
    effective_patients_total = eff_k + z_k,
    weighted_mean_delta = mu,
    weighted_se_delta = se_mu,
    weighted_sd_delta = weighted_sd,
    between_patient_sd = tau,
    tau_between_patient = tau,
    i2_heterogeneity = i2,
    sign_consistency = sign_consistency,
    n_positive_patients = sum(delta > cohort_transition_effect_threshold),
    n_negative_patients = sum(delta < -cohort_transition_effect_threshold),
    n_near_zero_patients = sum(abs(delta) <= cohort_transition_effect_threshold),
    p_beneficial = 1 - stats::pnorm(cohort_transition_effect_threshold, mean = mu, sd = se_mu),
    p_deleterious = stats::pnorm(-cohort_transition_effect_threshold, mean = mu, sd = se_mu),
    p_neutral = stats::pnorm(cohort_transition_effect_threshold, mean = mu, sd = se_mu) -
      stats::pnorm(-cohort_transition_effect_threshold, mean = mu, sd = se_mu),
    stringsAsFactors = FALSE
  )
}

#' Classify cohort transition groups for v2 selective borrowing
#'
#' @param patient_group_summaries Patient-level summaries.
#' @return A data frame of group classes and recommended borrowing controls.
#' @export
classify_cohort_transition_groups <- function(patient_group_summaries,
                                              cohort_transition_sd_floor = 0.05,
                                              cohort_transition_patient_sd_floor = 0.10,
                                              cohort_transition_min_patients_consistent = 3L,
                                              cohort_transition_min_effective_patients = 3,
                                              cohort_transition_min_effective_observed = 3,
                                              cohort_transition_effect_threshold = 0.02,
                                              cohort_transition_sign_consistency_threshold = 0.75,
                                              cohort_transition_high_heterogeneity_i2 = 0.50,
                                              cohort_transition_high_between_patient_sd = 0.10,
                                              cohort_transition_context_heterogeneity_drop = 0.25,
                                              ...) {
  validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
  validate_positive_finite(cohort_transition_patient_sd_floor, "cohort_transition_patient_sd_floor")
  validate_positive_integer(cohort_transition_min_patients_consistent, "cohort_transition_min_patients_consistent")
  validate_positive_finite(cohort_transition_min_effective_patients, "cohort_transition_min_effective_patients")
  validate_positive_finite(cohort_transition_min_effective_observed, "cohort_transition_min_effective_observed")
  if (!is.data.frame(patient_group_summaries) || !nrow(patient_group_summaries)) {
    borrowing <- compute_class_specific_borrowing("global_fallback", ...)
    return(data.frame(
      transition_group = "global",
      group = "global",
      group_level = "global",
      fallback_group = NA_character_,
      n_patients_observed = 0L,
      n_patients_zero = 0L,
      effective_patients_observed = 0,
      effective_patients_total = 0,
      weighted_mean_delta = 0,
      weighted_se_delta = Inf,
      weighted_sd_delta = NA_real_,
      between_patient_sd = NA_real_,
      tau_between_patient = NA_real_,
      i2_heterogeneity = NA_real_,
      sign_consistency = NA_real_,
      n_positive_patients = 0L,
      n_negative_patients = 0L,
      n_near_zero_patients = 0L,
      p_beneficial = 0,
      p_deleterious = 0,
      p_neutral = 0,
      effect_class = "global_fallback",
      heterogeneity_class = "sparse_unknown",
      recommended_lambda = borrowing$cohort_lambda,
      recommended_sd_multiplier = borrowing$sd_multiplier,
      recommended_use_for_zero = TRUE,
      use_for_observed = FALSE,
      use_for_low_information = TRUE,
      warning_flags = "no_patient_group_summaries",
      stringsAsFactors = FALSE
    ))
  }
  level_map <- c(
    gain_loss = "group_gain_loss",
    gain_loss_chr = "group_gain_loss_chr",
    gain_loss_chr_burden = "group_gain_loss_chr_burden",
    exact_event = "group_exact_event"
  )
  rows <- list()
  rows[[1L]] <- compute_transition_group_heterogeneity(
    patient_group_summaries = transform(patient_group_summaries, global = "global"),
    group_name = "global",
    group_level = "global",
    group_col = "global",
    cohort_transition_sd_floor = cohort_transition_sd_floor,
    cohort_transition_effect_threshold = cohort_transition_effect_threshold
  )
  for (level in names(level_map)) {
    col <- level_map[[level]]
    if (!col %in% names(patient_group_summaries)) next
    groups <- sort(unique(patient_group_summaries[[col]][!is.na(patient_group_summaries[[col]]) & nzchar(patient_group_summaries[[col]])]))
    for (group_name in groups) {
      rows[[length(rows) + 1L]] <- compute_transition_group_heterogeneity(
        patient_group_summaries = patient_group_summaries,
        group_name = group_name,
        group_level = level,
        group_col = col,
        cohort_transition_sd_floor = cohort_transition_sd_floor,
        cohort_transition_effect_threshold = cohort_transition_effect_threshold
      )
    }
  }
  out <- do.call(rbind, rows)
  for (i in seq_len(nrow(out))) {
    if (out$group_level[i] == "global") {
      cls <- list(
        effect_class = "global_fallback",
        heterogeneity_class = "sparse_unknown",
        p_beneficial = out$p_beneficial[i],
        p_deleterious = out$p_deleterious[i],
        p_neutral = out$p_neutral[i]
      )
    } else {
      cls <- compute_transition_group_class(
        k = out$n_patients_observed[i],
        eff_k = out$effective_patients_observed[i],
        effective_observed = out$effective_patients_observed[i],
        mu = out$weighted_mean_delta[i],
        se_mu = out$weighted_se_delta[i],
        tau = out$tau_between_patient[i],
        i2 = out$i2_heterogeneity[i],
        sign_consistency = out$sign_consistency[i],
        cohort_transition_min_patients_consistent = cohort_transition_min_patients_consistent,
        cohort_transition_min_effective_patients = cohort_transition_min_effective_patients,
        cohort_transition_min_effective_observed = cohort_transition_min_effective_observed,
        cohort_transition_effect_threshold = cohort_transition_effect_threshold,
        cohort_transition_sign_consistency_threshold = cohort_transition_sign_consistency_threshold,
        cohort_transition_high_heterogeneity_i2 = cohort_transition_high_heterogeneity_i2,
        cohort_transition_high_between_patient_sd = cohort_transition_high_between_patient_sd
      )
    }
    out$effect_class[i] <- cls$effect_class
    out$heterogeneity_class[i] <- cls$heterogeneity_class
    out$p_beneficial[i] <- cls$p_beneficial
    out$p_deleterious[i] <- cls$p_deleterious
    out$p_neutral[i] <- cls$p_neutral
  }
  # If a fine group is stable while its chromosome-level parent is variable,
  # label it as context-dependent so refit diagnostics show why the finer
  # context is preferred.
  fine_idx <- which(out$group_level %in% c("gain_loss_chr_burden", "exact_event") &
                      !out$effect_class %in% c("sparse_unknown", "high_variable"))
  for (idx in fine_idx) {
    parent_group <- sub("_burden_(low|neutral|high)$", "", out$group[idx])
    parent <- out[out$group_level == "gain_loss_chr" & out$group == parent_group, , drop = FALSE]
    if (nrow(parent) &&
        parent$effect_class[1] == "high_variable" &&
        is.finite(parent$i2_heterogeneity[1]) &&
        is.finite(out$i2_heterogeneity[idx]) &&
        parent$i2_heterogeneity[1] - out$i2_heterogeneity[idx] >= cohort_transition_context_heterogeneity_drop) {
      out$effect_class[idx] <- "context_dependent"
    }
  }
  borrow <- lapply(out$effect_class, compute_class_specific_borrowing, ...)
  out$recommended_lambda <- vapply(borrow, `[[`, numeric(1), "cohort_lambda")
  out$recommended_sd_multiplier <- vapply(borrow, `[[`, numeric(1), "sd_multiplier")
  out$recommended_use_for_zero <- vapply(borrow, `[[`, logical(1), "use_for_zero")
  out$use_for_observed <- vapply(borrow, `[[`, logical(1), "use_for_observed")
  out$use_for_low_information <- vapply(borrow, `[[`, logical(1), "use_for_low_information")
  out$warning_flags <- ""
  out$warning_flags[out$effect_class == "consistent_beneficial"] <- "survivor_bias_warning"
  out$warning_flags[out$effect_class == "sparse_unknown"] <- "sparse_unknown"
  out$warning_flags[out$effect_class == "high_variable"] <- "high_variable"
  rownames(out) <- NULL
  out
}

cohort_transition_fit_group <- function(records,
                                        group_name,
                                        sd_floor,
                                        patient_sd_floor,
                                        zero_weight_cap_ratio,
                                        zero_expected_count_cap,
                                        zero_mean_shift_cap,
                                        fallback_mu = 0,
                                        fallback_sigma = sd_floor,
                                        zero_likelihood_approximation = TRUE) {
  if (!nrow(records)) {
    return(data.frame(
      group = group_name,
      n_records_total = 0L,
      n_patients = 0L,
      n_observed_records = 0L,
      n_zero_records = 0L,
      effective_n = 0,
      mu = fallback_mu,
      sigma = max(fallback_sigma, sd_floor),
      sigma_with_patient_heterogeneity = sqrt(max(fallback_sigma, sd_floor)^2 + patient_sd_floor^2),
      sd_floor_used = TRUE,
      fallback_group = NA_character_,
      warning_flags = "no_records",
      stringsAsFactors = FALSE
    ))
  }
  weights <- cohort_transition_record_weights(records, zero_weight_cap_ratio = zero_weight_cap_ratio)
  observed <- records$source_type != "informative_zero" &
    is.finite(records$delta_hat) &
    is.finite(records$delta_se) &
    records$delta_se >= 0 &
    weights > 0
  zero <- records$source_type == "informative_zero" &
    is.finite(records$expected_count_parent_like) &
    records$expected_count_parent_like > 0 &
    weights > 0
  effective_n <- sum(weights[observed | zero], na.rm = TRUE)
  n_patients <- length(unique(records$patient_id[observed | zero]))
  warning_flags <- character(0)

  if (!any(observed) && !any(zero)) {
    warning_flags <- c(warning_flags, "no_weighted_evidence")
    mu <- fallback_mu
    sigma <- max(fallback_sigma, sd_floor)
  } else {
    obs_delta <- records$delta_hat[observed]
    obs_se <- pmax(records$delta_se[observed], sd_floor)
    obs_w <- weights[observed]
    zero_lambda_raw <- records$expected_count_parent_like[zero]
    zero_lambda <- pmin(zero_lambda_raw, zero_expected_count_cap)
    zero_w <- weights[zero]
    if (length(zero_lambda_raw) && any(zero_lambda_raw > zero_expected_count_cap, na.rm = TRUE)) {
      warning_flags <- c(warning_flags, "zero_expected_count_capped")
    }
    if (length(obs_delta)) {
      mu_start <- stats::weighted.mean(obs_delta, obs_w)
      sigma_start <- sqrt(sum(obs_w * (obs_delta - mu_start)^2) / max(sum(obs_w), .Machine$double.eps))
      if (!is.finite(sigma_start) || sigma_start <= 0) sigma_start <- fallback_sigma
    } else {
      mu_start <- fallback_mu
      sigma_start <- fallback_sigma
    }
    sigma_start <- max(sigma_start, sd_floor)
    delta_span <- range(c(obs_delta, fallback_mu, 0), finite = TRUE)
    if (length(delta_span) != 2L || any(!is.finite(delta_span))) delta_span <- c(-1, 1)
    span <- max(1, diff(delta_span), abs(delta_span))
    lower <- c(delta_span[1] - 3 * span, log(sd_floor))
    upper <- c(delta_span[2] + 3 * span, log(max(10 * span, sd_floor * 10)))
    if (length(obs_delta) && length(zero_lambda) &&
        is.finite(zero_mean_shift_cap) && zero_mean_shift_cap > 0) {
      lower[1] <- max(lower[1], mu_start - zero_mean_shift_cap)
      upper[1] <- min(upper[1], mu_start + zero_mean_shift_cap)
      warning_flags <- c(warning_flags, "zero_mean_shift_cap_active")
    }
    grid_z <- seq(-6, 6, length.out = 81L)
    grid_w <- stats::dnorm(grid_z)
    grid_w <- grid_w / sum(grid_w)
    objective <- function(par) {
      mu <- par[1]
      sigma <- exp(par[2])
      if (!is.finite(mu) || !is.finite(sigma) || sigma < sd_floor) return(1e12)
      sigma_obs <- sqrt(sigma^2 + obs_se^2 + patient_sd_floor^2)
      nll <- 0
      if (length(obs_delta)) {
        nll <- nll - sum(obs_w * stats::dnorm(obs_delta, mean = mu, sd = sigma_obs, log = TRUE))
      }
      if (length(zero_lambda)) {
        sigma_zero <- sqrt(sigma^2 + patient_sd_floor^2)
        for (i in seq_along(zero_lambda)) {
          delta_grid <- mu + sigma_zero * grid_z
          p0 <- exp(-zero_lambda[i] * exp(delta_grid))
          marginal <- sum(grid_w * p0)
          nll <- nll - zero_w[i] * log(max(marginal, .Machine$double.xmin))
        }
      }
      if (!is.finite(nll)) 1e12 else nll
    }
    opt <- try(stats::nlminb(
      start = c(mu_start, log(sigma_start)),
      objective = objective,
      lower = lower,
      upper = upper,
      control = list(iter.max = 100, eval.max = 200)
    ), silent = TRUE)
    if (inherits(opt, "try-error") || !is.finite(opt$objective)) {
      warning_flags <- c(warning_flags, "optim_failed")
      mu <- mu_start
      sigma <- sigma_start
    } else {
      mu <- opt$par[1]
      sigma <- exp(opt$par[2])
      if (length(obs_delta) && length(zero_lambda) &&
          is.finite(zero_mean_shift_cap) && zero_mean_shift_cap > 0 &&
          abs(mu - mu_start) >= zero_mean_shift_cap - sqrt(.Machine$double.eps)) {
        warning_flags <- c(warning_flags, "zero_mean_shift_cap_hit")
      }
    }
    if (length(zero_lambda)) {
      warning_flags <- c(warning_flags, "zero_likelihood_approximation")
    }
  }
  sigma <- max(sigma, sd_floor)
  data.frame(
    group = group_name,
    n_records_total = nrow(records),
    n_patients = n_patients,
    n_observed_records = sum(observed),
    n_zero_records = sum(zero),
    effective_n = effective_n,
    mu = mu,
    sigma = sigma,
    sigma_with_patient_heterogeneity = sqrt(sigma^2 + patient_sd_floor^2),
    sd_floor_used = sigma <= sd_floor + sqrt(.Machine$double.eps),
    fallback_group = NA_character_,
    warning_flags = paste(unique(warning_flags), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

cohort_transition_build_prior_tables <- function(records,
                                                 grouping,
                                                 min_patients_per_group,
                                                 min_effective_n,
                                                 sd_floor,
                                                 patient_sd_floor,
                                                 global_fallback,
                                                 zero_weight_cap_ratio,
                                                 zero_expected_count_cap,
                                                 zero_mean_shift_cap) {
  grouping <- match.arg(grouping, c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"))
  group_levels <- c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event")
  level_cols <- paste0("group_", group_levels)
  for (col in level_cols) {
    if (!col %in% names(records)) {
      records[[col]] <- NA_character_
    }
  }
  target_col <- cohort_transition_group_column(grouping)
  global <- cohort_transition_fit_group(
    records = records,
    group_name = "global",
    sd_floor = sd_floor,
    patient_sd_floor = patient_sd_floor,
    zero_weight_cap_ratio = zero_weight_cap_ratio,
    zero_expected_count_cap = zero_expected_count_cap,
    zero_mean_shift_cap = zero_mean_shift_cap
  )
  all_estimates <- list(global = transform(global, level = "global"))
  for (level in group_levels) {
    col <- paste0("group_", level)
    groups <- sort(unique(records[[col]][!is.na(records[[col]]) & nzchar(records[[col]])]))
    if (!length(groups)) next
    est <- do.call(rbind, lapply(groups, function(group_name) {
      cohort_transition_fit_group(
        records = records[records[[col]] == group_name, , drop = FALSE],
        group_name = group_name,
        sd_floor = sd_floor,
        patient_sd_floor = patient_sd_floor,
        zero_weight_cap_ratio = zero_weight_cap_ratio,
        zero_expected_count_cap = zero_expected_count_cap,
        zero_mean_shift_cap = zero_mean_shift_cap,
        fallback_mu = global$mu[1],
        fallback_sigma = global$sigma[1]
      )
    }))
    est$level <- level
    all_estimates[[level]] <- est
  }
  all_estimates_df <- do.call(rbind, all_estimates)
  rownames(all_estimates_df) <- NULL

  target_groups <- sort(unique(records[[target_col]][!is.na(records[[target_col]]) & nzchar(records[[target_col]])]))
  fallback_levels <- switch(
    grouping,
    exact_event = c("exact_event", "gain_loss_chr_burden", "gain_loss_chr", "gain_loss", "global"),
    gain_loss_chr_burden = c("gain_loss_chr_burden", "gain_loss_chr", "gain_loss", "global"),
    gain_loss_chr = c("gain_loss_chr", "gain_loss", "global"),
    gain_loss = c("gain_loss", "global")
  )
  target_rows <- lapply(target_groups, function(group_name) {
    group_records <- records[records[[target_col]] == group_name, , drop = FALSE]
    exemplar <- group_records[1L, , drop = FALSE]
    candidate_names <- list(
      exact_event = exemplar$group_exact_event,
      gain_loss_chr_burden = exemplar$group_gain_loss_chr_burden,
      gain_loss_chr = exemplar$group_gain_loss_chr,
      gain_loss = exemplar$group_gain_loss,
      global = "global"
    )
    chosen <- NULL
    for (level in fallback_levels) {
      cand <- all_estimates_df[
        all_estimates_df$level == level &
          all_estimates_df$group == candidate_names[[level]],
        ,
        drop = FALSE
      ]
      if (!nrow(cand)) next
      enough <- cand$n_patients[1] >= min_patients_per_group &&
        cand$effective_n[1] >= min_effective_n
      if (isTRUE(enough) || (level == "global" && isTRUE(global_fallback))) {
        chosen <- cand[1L, , drop = FALSE]
        break
      }
    }
    if (is.null(chosen)) {
      chosen <- global
      chosen$level <- "global"
    }
    chosen$requested_group <- group_name
    chosen$fallback_group <- if (identical(chosen$group[1], group_name)) NA_character_ else chosen$group[1]
    chosen$group <- group_name
    chosen
  })
  group_priors <- if (length(target_rows)) do.call(rbind, target_rows) else all_estimates_df[FALSE, , drop = FALSE]
  rownames(group_priors) <- NULL
  list(global_prior = global, group_priors = group_priors, all_group_priors = all_estimates_df)
}

#' Learn a cohort-level transition-effect prior
#'
#' @param records Transition records produced by `extract_cohort_transition_records()`.
#' @param leave_one_patient_out Whether to store leave-one-patient-out priors.
#' @param grouping Transition grouping mode.
#' @param cohort_transition_min_patients_per_group Minimum number of patients
#'   required before using a group-specific prior without fallback.
#' @param cohort_transition_min_effective_n Minimum effective evidence required
#'   before using a group-specific prior without fallback.
#' @param cohort_transition_sd_floor Minimum transition-effect prior standard deviation.
#' @param cohort_transition_patient_sd_floor Patient heterogeneity standard
#'   deviation floor added to transition priors. The default is conservative to
#'   prevent repeated bootstrap/path records from producing an overconfident
#'   patient-level transition prior.
#' @param cohort_transition_global_fallback Whether under-supported groups can
#'   fall back to the global transition prior.
#' @param cohort_transition_zero_weight_cap_ratio Cap on total zero evidence
#'   weight relative to observed evidence.
#' @param cohort_transition_zero_expected_count_cap Cap applied to the
#'   parent-like expected count inside zero-censoring likelihoods. Counts above
#'   this value are still classified as informative zeros, but are not allowed
#'   to add unbounded pressure to the cohort transition mean.
#' @param cohort_transition_zero_mean_shift_cap Maximum absolute shift that zero
#'   censoring evidence can impose on a group mean away from the observed
#'   transition-effect mean. Set `NULL` to disable.
#' @param ... Reserved for future prior fitting controls.
#' @return A cohort-transition prior object.
#' @export
learn_cohort_transition_prior_v1 <- function(records,
                                             leave_one_patient_out = TRUE,
                                             grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                             cohort_transition_min_patients_per_group = 2L,
                                             cohort_transition_min_effective_n = 3L,
                                             cohort_transition_sd_floor = 1e-3,
                                             cohort_transition_patient_sd_floor = 0.1,
                                             cohort_transition_global_fallback = TRUE,
                                             cohort_transition_zero_weight_cap_ratio = 1.0,
                                             cohort_transition_zero_expected_count_cap = 10.0,
                                             cohort_transition_zero_mean_shift_cap = 0.2,
                                             ...) {
  grouping <- match.arg(grouping)
  validate_positive_integer(cohort_transition_min_patients_per_group, "cohort_transition_min_patients_per_group")
  validate_positive_finite(cohort_transition_min_effective_n, "cohort_transition_min_effective_n")
  validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
  validate_positive_finite(cohort_transition_patient_sd_floor, "cohort_transition_patient_sd_floor")
  validate_scalar_logical(cohort_transition_global_fallback, "cohort_transition_global_fallback")
  validate_nonnegative_finite(cohort_transition_zero_weight_cap_ratio, "cohort_transition_zero_weight_cap_ratio")
  validate_positive_finite(cohort_transition_zero_expected_count_cap, "cohort_transition_zero_expected_count_cap")
  if (!is.null(cohort_transition_zero_mean_shift_cap)) {
    validate_positive_finite(cohort_transition_zero_mean_shift_cap, "cohort_transition_zero_mean_shift_cap")
  } else {
    cohort_transition_zero_mean_shift_cap <- Inf
  }
  if (!is.data.frame(records) || !nrow(records)) {
    stop("`records` must contain at least one transition record.", call. = FALSE)
  }
  records <- cohort_transition_assign_groups(records, grouping)
  patient_ids <- sort(unique(as.character(records$patient_id)))
  tables <- cohort_transition_build_prior_tables(
    records = records,
    grouping = grouping,
    min_patients_per_group = cohort_transition_min_patients_per_group,
    min_effective_n = cohort_transition_min_effective_n,
    sd_floor = cohort_transition_sd_floor,
    patient_sd_floor = cohort_transition_patient_sd_floor,
    global_fallback = cohort_transition_global_fallback,
    zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio,
    zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
    zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap
  )

  loo_priors <- list()
  if (isTRUE(leave_one_patient_out)) {
    for (patient_id in patient_ids) {
      rec_loo <- records[records$patient_id != patient_id, , drop = FALSE]
      if (!nrow(rec_loo)) {
        loo_priors[[patient_id]] <- list(
          contributing_patients = character(0),
          group_priors = tables$group_priors[FALSE, , drop = FALSE],
          global_prior = tables$global_prior[FALSE, , drop = FALSE],
          diagnostics = list(fallback_reason = "no_loo_records")
        )
        next
      }
      loo_tables <- cohort_transition_build_prior_tables(
        records = rec_loo,
        grouping = grouping,
        min_patients_per_group = cohort_transition_min_patients_per_group,
        min_effective_n = cohort_transition_min_effective_n,
        sd_floor = cohort_transition_sd_floor,
        patient_sd_floor = cohort_transition_patient_sd_floor,
        global_fallback = cohort_transition_global_fallback,
        zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio,
        zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
        zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap
      )
      loo_priors[[patient_id]] <- list(
        contributing_patients = sort(unique(as.character(rec_loo$patient_id))),
        group_priors = loo_tables$group_priors,
        global_prior = loo_tables$global_prior,
        all_group_priors = loo_tables$all_group_priors,
        diagnostics = list(fallback_reason = NA_character_)
      )
    }
  }

  diagnostics <- list(
    n_patients = length(patient_ids),
    patient_ids = patient_ids,
    grouping = grouping,
    n_transition_records_total = nrow(records),
    n_observed_records = sum(records$source_type != "informative_zero"),
    n_zero_records = sum(records$child_is_zero, na.rm = TRUE),
    n_zero_retained = sum(records$source_type == "informative_zero"),
    n_zero_excluded_low_exposure = sum(records$child_is_zero, na.rm = TRUE) - sum(records$source_type == "informative_zero"),
    effective_zero_information = sum(cohort_transition_record_weights(
      records,
      zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio
    )[records$source_type == "informative_zero"], na.rm = TRUE),
    zero_to_observed_information_ratio = {
      w <- cohort_transition_record_weights(records, zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio)
      zero_w <- sum(w[records$source_type == "informative_zero"], na.rm = TRUE)
      obs_w <- sum(w[records$source_type != "informative_zero"], na.rm = TRUE)
      if (obs_w > 0) zero_w / obs_w else NA_real_
    },
    groups_estimated = unique(tables$group_priors$group),
    groups_fallback_to_coarser = tables$group_priors$group[!is.na(tables$group_priors$fallback_group) &
                                                            tables$group_priors$fallback_group != "global"],
    groups_fallback_to_global = tables$group_priors$group[!is.na(tables$group_priors$fallback_group) &
                                                           tables$group_priors$fallback_group == "global"],
    mu_by_group = stats::setNames(tables$group_priors$mu, tables$group_priors$group),
    sigma_by_group = stats::setNames(tables$group_priors$sigma, tables$group_priors$group),
    sigma_floor_used = any(tables$group_priors$sd_floor_used),
    patient_heterogeneity_sd = cohort_transition_patient_sd_floor,
    zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
    zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap,
    zero_expected_count_capped = any(
      records$source_type == "informative_zero" &
        is.finite(records$expected_count_parent_like) &
        records$expected_count_parent_like > cohort_transition_zero_expected_count_cap,
      na.rm = TRUE
    ),
    n_zero_expected_count_capped = sum(
      records$source_type == "informative_zero" &
        is.finite(records$expected_count_parent_like) &
        records$expected_count_parent_like > cohort_transition_zero_expected_count_cap,
      na.rm = TRUE
    ),
    leave_one_patient_out_used = isTRUE(leave_one_patient_out),
    patients_contributing_by_group = lapply(unique(tables$group_priors$group), function(group_name) {
      sort(unique(records$patient_id[records$transition_group == group_name]))
    }),
    zero_likelihood_approximation = any(grepl("zero_likelihood_approximation", tables$all_group_priors$warning_flags))
  )
  names(diagnostics$patients_contributing_by_group) <- unique(tables$group_priors$group)

  list(
    version = "cohort_transition_v1",
    grouping = grouping,
    global_prior = tables$global_prior,
    group_priors = tables$group_priors,
    all_group_priors = tables$all_group_priors,
    patient_ids = patient_ids,
    leave_one_patient_out = isTRUE(leave_one_patient_out),
    loo_priors = loo_priors,
    diagnostics = diagnostics
  )
}

cohort_transition_build_prior_tables_v2 <- function(patient_group_summaries,
                                                    grouping,
                                                    group_classes,
                                                    sd_floor,
                                                    patient_sd_floor,
                                                    global_fallback = TRUE) {
  grouping <- match.arg(grouping, c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"))
  validate_positive_finite(sd_floor, "cohort_transition_sd_floor")
  validate_positive_finite(patient_sd_floor, "cohort_transition_patient_sd_floor")
  if (!nrow(group_classes)) {
    group_classes <- classify_cohort_transition_groups(
      patient_group_summaries,
      cohort_transition_sd_floor = sd_floor,
      cohort_transition_patient_sd_floor = patient_sd_floor
    )
  }
  classes <- group_classes
  classes$mu <- classes$weighted_mean_delta
  classes$sigma <- pmax(classes$between_patient_sd, sd_floor)
  classes$sigma[!is.finite(classes$sigma)] <- sd_floor
  classes$tau_between_patient[!is.finite(classes$tau_between_patient)] <- classes$sigma[!is.finite(classes$tau_between_patient)]
  classes$sigma_with_patient_heterogeneity <- sqrt(classes$sigma^2 + patient_sd_floor^2)
  classes$sd_floor_used <- classes$sigma <= sd_floor + sqrt(.Machine$double.eps)
  classes$cohort_lambda <- classes$recommended_lambda
  classes$sd_multiplier <- classes$recommended_sd_multiplier
  classes$effective_prior_sd <- pmax(classes$sigma_with_patient_heterogeneity * classes$sd_multiplier, sd_floor)
  classes$n_patients <- classes$n_patients_observed
  classes$effective_patients <- classes$effective_patients_observed
  classes$n_observed_patient_summaries <- classes$n_patients_observed
  classes$n_zero_patient_summaries <- classes$n_patients_zero
  classes$n_records_total <- NA_integer_
  classes$n_observed_records <- classes$n_patients_observed
  classes$n_zero_records <- classes$n_patients_zero
  classes$effective_n <- classes$effective_patients_total
  classes$level <- classes$group_level

  global <- classes[classes$group_level == "global" & classes$group == "global", , drop = FALSE]
  if (!nrow(global)) {
    global <- classes[1L, , drop = FALSE]
    global$group <- "global"
    global$transition_group <- "global"
    global$group_level <- "global"
    global$level <- "global"
    global$effect_class <- "global_fallback"
    global$cohort_lambda <- 0.05
    global$sd_multiplier <- 4
    global$effective_prior_sd <- max(sqrt(sd_floor^2 + patient_sd_floor^2) * 4, sd_floor)
  }

  target_col <- cohort_transition_group_column(grouping)
  target_groups <- if (is.data.frame(patient_group_summaries) && nrow(patient_group_summaries) &&
                       target_col %in% names(patient_group_summaries)) {
    sort(unique(patient_group_summaries[[target_col]][!is.na(patient_group_summaries[[target_col]]) &
                                                        nzchar(patient_group_summaries[[target_col]])]))
  } else {
    character(0)
  }
  fallback_levels <- switch(
    grouping,
    exact_event = c("exact_event", "gain_loss_chr_burden", "gain_loss_chr", "gain_loss", "global"),
    gain_loss_chr_burden = c("gain_loss_chr_burden", "gain_loss_chr", "gain_loss", "global"),
    gain_loss_chr = c("gain_loss_chr", "gain_loss", "global"),
    gain_loss = c("gain_loss", "global")
  )
  target_rows <- lapply(target_groups, function(group_name) {
    exemplar <- patient_group_summaries[patient_group_summaries[[target_col]] == group_name, , drop = FALSE][1L, , drop = FALSE]
    candidate_names <- list(
      exact_event = exemplar$group_exact_event,
      gain_loss_chr_burden = exemplar$group_gain_loss_chr_burden,
      gain_loss_chr = exemplar$group_gain_loss_chr,
      gain_loss = exemplar$group_gain_loss,
      global = "global"
    )
    chosen <- NULL
    requested_missing <- TRUE
    for (level in fallback_levels) {
      cand_name <- candidate_names[[level]]
      cand <- classes[
        classes$group_level == level &
          classes$group == cand_name,
        ,
        drop = FALSE
      ]
      if (!nrow(cand)) next
      if (level == grouping) requested_missing <- FALSE
      usable <- cand$effect_class[1] %in% c(
        "consistent_deleterious",
        "consistent_neutral",
        "consistent_beneficial",
        "context_dependent",
        "high_variable",
        "sparse_unknown"
      )
      if (usable || (level == "global" && isTRUE(global_fallback))) {
        chosen <- cand[1L, , drop = FALSE]
        break
      }
    }
    if (is.null(chosen)) {
      chosen <- global[1L, , drop = FALSE]
    }
    fallback_multiplier <- 1
    fallback_group <- NA_character_
    if (!identical(chosen$group[1], group_name)) {
      fallback_group <- chosen$group[1]
      if (identical(chosen$group_level[1], "global")) {
        fallback_multiplier <- 0.25
        chosen$effect_class <- "global_fallback"
        chosen$heterogeneity_class <- "sparse_unknown"
      } else {
        fallback_multiplier <- 0.5
      }
    }
    chosen$requested_group <- group_name
    chosen$fallback_group <- fallback_group
    chosen$fallback_group_used <- fallback_group
    chosen$fallback_multiplier <- fallback_multiplier
    chosen$group <- group_name
    chosen$transition_group <- group_name
    chosen$cohort_lambda <- chosen$cohort_lambda * fallback_multiplier
    chosen$effective_prior_sd <- chosen$effective_prior_sd / sqrt(max(fallback_multiplier, .Machine$double.eps))
    chosen$effective_prior_sd <- if (identical(chosen$group_level[1], "global")) {
      max(chosen$effective_prior_sd, chosen$sigma_with_patient_heterogeneity * 2)
    } else {
      chosen$effective_prior_sd
    }
    if (isTRUE(requested_missing)) {
      chosen$warning_flags <- paste(unique(c(chosen$warning_flags, "requested_group_missing")), collapse = ";")
    }
    chosen
  })
  group_priors <- if (length(target_rows)) do.call(rbind, target_rows) else classes[FALSE, , drop = FALSE]
  rownames(group_priors) <- NULL
  rownames(classes) <- NULL
  list(global_prior = global[1L, , drop = FALSE], group_priors = group_priors, all_group_priors = classes)
}

#' Estimate a shrunk patient-specific transition shift
#'
#' @param patient_records Patient-level observed transition summaries.
#' @param prior Cohort transition prior.
#' @param cohort_transition_patient_shift_min_records Minimum observed records
#'   needed to estimate a nonzero shift.
#' @param cohort_transition_patient_shift_shrinkage_sd Shrinkage SD for the
#'   patient-level residual shift.
#' @return A one-row data frame with the shrunk shift and reliability.
#' @export
estimate_patient_transition_shift <- function(patient_records,
                                              prior,
                                              cohort_transition_patient_shift_min_records = 3L,
                                              cohort_transition_patient_shift_shrinkage_sd = 0.10,
                                              ...) {
  validate_positive_integer(cohort_transition_patient_shift_min_records, "cohort_transition_patient_shift_min_records")
  validate_positive_finite(cohort_transition_patient_shift_shrinkage_sd, "cohort_transition_patient_shift_shrinkage_sd")
  if (!is.data.frame(patient_records) || !nrow(patient_records)) {
    return(data.frame(
      patient_delta_shift = 0,
      patient_delta_shift_n_records = 0L,
      patient_delta_shift_reliability = 0,
      stringsAsFactors = FALSE
    ))
  }
  summaries <- patient_records
  if (!"delta_patient_mean" %in% names(summaries)) {
    filtered <- filter_cohort_transition_records(patient_records, ...)
    summaries <- aggregate_cohort_transition_records_by_patient(
      filtered$kept_records,
      grouping = prior$grouping,
      cohort_transition_sd_floor = if (is.null(prior$diagnostics$sd_floor)) 0.05 else prior$diagnostics$sd_floor
    )
  }
  obs <- summaries[
    summaries$cohort_transition_evidence_type == "observed_delta_evidence" &
      is.finite(summaries$delta_patient_mean),
    ,
    drop = FALSE
  ]
  if (nrow(obs) < cohort_transition_patient_shift_min_records) {
    return(data.frame(
      patient_delta_shift = 0,
      patient_delta_shift_n_records = nrow(obs),
      patient_delta_shift_reliability = 0,
      stringsAsFactors = FALSE
    ))
  }
  mu <- numeric(nrow(obs))
  for (i in seq_len(nrow(obs))) {
    row <- prior$group_priors[prior$group_priors$group == obs$transition_group[i], , drop = FALSE]
    if (!nrow(row)) {
      row <- prior$global_prior
    }
    mu[i] <- row$mu[1]
  }
  residual <- obs$delta_patient_mean - mu
  se <- pmax(obs$delta_patient_se, 0.05)
  w <- 1 / (se^2 + cohort_transition_patient_shift_shrinkage_sd^2)
  w[!is.finite(w) | w <= 0] <- 1
  raw_shift <- stats::weighted.mean(residual, w)
  reliability <- sum(w) / (sum(w) + 1 / cohort_transition_patient_shift_shrinkage_sd^2)
  shift <- raw_shift * reliability
  data.frame(
    patient_delta_shift = shift,
    patient_delta_shift_n_records = nrow(obs),
    patient_delta_shift_reliability = reliability,
    stringsAsFactors = FALSE
  )
}

learn_cohort_transition_prior_v2 <- function(records,
                                             leave_one_patient_out = TRUE,
                                             grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                             cohort_transition_min_patients_per_group = 2L,
                                             cohort_transition_min_effective_n = 3L,
                                             cohort_transition_sd_floor = 0.05,
                                             cohort_transition_patient_sd_floor = 0.10,
                                             cohort_transition_global_fallback = TRUE,
                                             cohort_transition_zero_weight_cap_ratio = 0.25,
                                             cohort_transition_zero_expected_count_cap = 10.0,
                                             cohort_transition_zero_mean_shift_cap = 0.2,
                                             cohort_transition_use_prior_dominated_records = FALSE,
                                             cohort_transition_use_boundary_records = FALSE,
                                             cohort_transition_max_delta_se = NULL,
                                             cohort_transition_max_delta_se_quantile = 0.75,
                                             cohort_transition_min_path_responsibility = 0.05,
                                             cohort_transition_min_observed_count = 1L,
                                             cohort_transition_classify_groups = TRUE,
                                             cohort_transition_min_patients_consistent = 3L,
                                             cohort_transition_min_effective_patients = 3,
                                             cohort_transition_min_effective_observed = 3,
                                             cohort_transition_effect_threshold = 0.02,
                                             cohort_transition_sign_consistency_threshold = 0.75,
                                             cohort_transition_high_heterogeneity_i2 = 0.50,
                                             cohort_transition_high_between_patient_sd = 0.10,
                                             cohort_transition_context_heterogeneity_drop = 0.25,
                                             cohort_transition_lambda_consistent_deleterious = 0.50,
                                             cohort_transition_lambda_consistent_neutral = 0.25,
                                             cohort_transition_lambda_consistent_beneficial = 0.15,
                                             cohort_transition_lambda_context_dependent = 0.30,
                                             cohort_transition_lambda_high_variable = 0.00,
                                             cohort_transition_lambda_sparse_unknown = 0.00,
                                             cohort_transition_lambda_global_fallback = 0.05,
                                             cohort_transition_sd_multiplier_consistent_deleterious = 1.0,
                                             cohort_transition_sd_multiplier_consistent_neutral = 1.5,
                                             cohort_transition_sd_multiplier_consistent_beneficial = 2.0,
                                             cohort_transition_sd_multiplier_context_dependent = 1.5,
                                             cohort_transition_sd_multiplier_high_variable = 4.0,
                                             cohort_transition_sd_multiplier_sparse_unknown = 4.0,
                                             cohort_transition_sd_multiplier_global_fallback = 4.0,
                                             cohort_transition_patient_shift = TRUE,
                                             cohort_transition_patient_shift_min_records = 3L,
                                             cohort_transition_patient_shift_shrinkage_sd = 0.10,
                                             cohort_transition_zero_as_censoring_only = TRUE,
                                             cohort_transition_zero_min_expected_count = 3.0,
                                             ...) {
  grouping <- match.arg(grouping)
  validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
  validate_positive_finite(cohort_transition_patient_sd_floor, "cohort_transition_patient_sd_floor")
  validate_scalar_logical(cohort_transition_global_fallback, "cohort_transition_global_fallback")
  validate_scalar_logical(cohort_transition_patient_shift, "cohort_transition_patient_shift")
  if (!is.data.frame(records) || !nrow(records)) {
    stop("`records` must contain at least one transition record.", call. = FALSE)
  }
  records <- cohort_transition_assign_groups(records, grouping)
  patient_ids <- sort(unique(as.character(records$patient_id)))
  filtered <- filter_cohort_transition_records(
    records,
    cohort_transition_use_prior_dominated_records = cohort_transition_use_prior_dominated_records,
    cohort_transition_use_boundary_records = cohort_transition_use_boundary_records,
    cohort_transition_max_delta_se = cohort_transition_max_delta_se,
    cohort_transition_max_delta_se_quantile = cohort_transition_max_delta_se_quantile,
    cohort_transition_min_path_responsibility = cohort_transition_min_path_responsibility,
    cohort_transition_min_observed_count = cohort_transition_min_observed_count,
    cohort_transition_zero_min_expected_count = cohort_transition_zero_min_expected_count,
    cohort_transition_zero_as_censoring_only = cohort_transition_zero_as_censoring_only
  )
  summaries <- aggregate_cohort_transition_records_by_patient(
    filtered$kept_records,
    grouping = grouping,
    cohort_transition_sd_floor = cohort_transition_sd_floor,
    cohort_transition_zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio
  )
  class_args <- list(
    cohort_transition_sd_floor = cohort_transition_sd_floor,
    cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor,
    cohort_transition_min_patients_consistent = cohort_transition_min_patients_consistent,
    cohort_transition_min_effective_patients = cohort_transition_min_effective_patients,
    cohort_transition_min_effective_observed = cohort_transition_min_effective_observed,
    cohort_transition_effect_threshold = cohort_transition_effect_threshold,
    cohort_transition_sign_consistency_threshold = cohort_transition_sign_consistency_threshold,
    cohort_transition_high_heterogeneity_i2 = cohort_transition_high_heterogeneity_i2,
    cohort_transition_high_between_patient_sd = cohort_transition_high_between_patient_sd,
    cohort_transition_context_heterogeneity_drop = cohort_transition_context_heterogeneity_drop,
    cohort_transition_lambda_consistent_deleterious = cohort_transition_lambda_consistent_deleterious,
    cohort_transition_lambda_consistent_neutral = cohort_transition_lambda_consistent_neutral,
    cohort_transition_lambda_consistent_beneficial = cohort_transition_lambda_consistent_beneficial,
    cohort_transition_lambda_context_dependent = cohort_transition_lambda_context_dependent,
    cohort_transition_lambda_high_variable = cohort_transition_lambda_high_variable,
    cohort_transition_lambda_sparse_unknown = cohort_transition_lambda_sparse_unknown,
    cohort_transition_lambda_global_fallback = cohort_transition_lambda_global_fallback,
    cohort_transition_sd_multiplier_consistent_deleterious = cohort_transition_sd_multiplier_consistent_deleterious,
    cohort_transition_sd_multiplier_consistent_neutral = cohort_transition_sd_multiplier_consistent_neutral,
    cohort_transition_sd_multiplier_consistent_beneficial = cohort_transition_sd_multiplier_consistent_beneficial,
    cohort_transition_sd_multiplier_context_dependent = cohort_transition_sd_multiplier_context_dependent,
    cohort_transition_sd_multiplier_high_variable = cohort_transition_sd_multiplier_high_variable,
    cohort_transition_sd_multiplier_sparse_unknown = cohort_transition_sd_multiplier_sparse_unknown,
    cohort_transition_sd_multiplier_global_fallback = cohort_transition_sd_multiplier_global_fallback
  )
  group_classes <- do.call(classify_cohort_transition_groups, c(list(patient_group_summaries = summaries), class_args))
  tables <- cohort_transition_build_prior_tables_v2(
    patient_group_summaries = summaries,
    grouping = grouping,
    group_classes = group_classes,
    sd_floor = cohort_transition_sd_floor,
    patient_sd_floor = cohort_transition_patient_sd_floor,
    global_fallback = cohort_transition_global_fallback
  )
  prior_base <- list(
    version = "cohort_transition_v2",
    grouping = grouping,
    global_prior = tables$global_prior,
    group_priors = tables$group_priors,
    all_group_priors = tables$all_group_priors,
    patient_group_summaries = summaries,
    group_classes = group_classes,
    patient_ids = patient_ids,
    leave_one_patient_out = isTRUE(leave_one_patient_out),
    loo_priors = list(),
    diagnostics = list(sd_floor = cohort_transition_sd_floor)
  )
  patient_shifts <- list()
  if (isTRUE(cohort_transition_patient_shift)) {
    for (patient_id in patient_ids) {
      patient_summaries <- summaries[summaries$patient_id == patient_id, , drop = FALSE]
      patient_shifts[[patient_id]] <- estimate_patient_transition_shift(
        patient_summaries,
        prior = prior_base,
        cohort_transition_patient_shift_min_records = cohort_transition_patient_shift_min_records,
        cohort_transition_patient_shift_shrinkage_sd = cohort_transition_patient_shift_shrinkage_sd
      )
    }
  }

  loo_priors <- list()
  if (isTRUE(leave_one_patient_out)) {
    for (patient_id in patient_ids) {
      rec_loo <- records[records$patient_id != patient_id, , drop = FALSE]
      if (!nrow(rec_loo)) {
        loo_priors[[patient_id]] <- list(
          contributing_patients = character(0),
          group_priors = tables$group_priors[FALSE, , drop = FALSE],
          global_prior = tables$global_prior[FALSE, , drop = FALSE],
          all_group_priors = tables$all_group_priors[FALSE, , drop = FALSE],
          patient_group_summaries = summaries[FALSE, , drop = FALSE],
          group_classes = group_classes[FALSE, , drop = FALSE],
          diagnostics = list(fallback_reason = "no_loo_records")
        )
        next
      }
      loo_filtered <- filter_cohort_transition_records(
        rec_loo,
        cohort_transition_use_prior_dominated_records = cohort_transition_use_prior_dominated_records,
        cohort_transition_use_boundary_records = cohort_transition_use_boundary_records,
        cohort_transition_max_delta_se = cohort_transition_max_delta_se,
        cohort_transition_max_delta_se_quantile = cohort_transition_max_delta_se_quantile,
        cohort_transition_min_path_responsibility = cohort_transition_min_path_responsibility,
        cohort_transition_min_observed_count = cohort_transition_min_observed_count,
        cohort_transition_zero_min_expected_count = cohort_transition_zero_min_expected_count,
        cohort_transition_zero_as_censoring_only = cohort_transition_zero_as_censoring_only
      )
      loo_summaries <- aggregate_cohort_transition_records_by_patient(
        loo_filtered$kept_records,
        grouping = grouping,
        cohort_transition_sd_floor = cohort_transition_sd_floor,
        cohort_transition_zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio
      )
      loo_classes <- do.call(classify_cohort_transition_groups, c(list(patient_group_summaries = loo_summaries), class_args))
      loo_tables <- cohort_transition_build_prior_tables_v2(
        patient_group_summaries = loo_summaries,
        grouping = grouping,
        group_classes = loo_classes,
        sd_floor = cohort_transition_sd_floor,
        patient_sd_floor = cohort_transition_patient_sd_floor,
        global_fallback = cohort_transition_global_fallback
      )
      loo_priors[[patient_id]] <- list(
        contributing_patients = sort(unique(as.character(rec_loo$patient_id))),
        group_priors = loo_tables$group_priors,
        global_prior = loo_tables$global_prior,
        all_group_priors = loo_tables$all_group_priors,
        patient_group_summaries = loo_summaries,
        group_classes = loo_classes,
        diagnostics = list(fallback_reason = NA_character_)
      )
    }
  }

  class_distribution <- table(group_classes$effect_class)
  diagnostics <- list(
    n_patients = length(patient_ids),
    patient_ids = patient_ids,
    grouping = grouping,
    version = "cohort_transition_v2",
    n_raw_transition_records = nrow(records),
    n_transition_records_total = nrow(records),
    n_patient_group_summaries = nrow(summaries),
    n_records_excluded_by_reason = filtered$diagnostics[names(cohort_transition_empty_filter_diagnostics())],
    n_prior_dominated_records_excluded = filtered$diagnostics$prior_dominated,
    n_boundary_records_excluded = filtered$diagnostics$boundary_dominated,
    n_observed_records = sum(records$source_type != "informative_zero", na.rm = TRUE),
    n_zero_records = sum(records$child_is_zero, na.rm = TRUE),
    n_zero_censoring_records = sum(filtered$kept_records$cohort_transition_evidence_type == "zero_censoring_evidence", na.rm = TRUE),
    n_low_exposure_zero_excluded = filtered$diagnostics$low_exposure_zero,
    n_zero_retained = sum(filtered$kept_records$cohort_transition_evidence_type == "zero_censoring_evidence", na.rm = TRUE),
    n_zero_excluded_low_exposure = filtered$diagnostics$low_exposure_zero,
    effective_zero_information = sum(summaries$zero_censoring_weight, na.rm = TRUE),
    zero_to_observed_information_ratio = {
      obs_w <- sum(summaries$observed_weight, na.rm = TRUE)
      zero_w <- sum(summaries$zero_censoring_weight, na.rm = TRUE)
      if (obs_w > 0) zero_w / obs_w else NA_real_
    },
    n_groups_sparse_unknown = sum(group_classes$effect_class == "sparse_unknown"),
    n_groups_high_variable = sum(group_classes$effect_class == "high_variable"),
    n_groups_context_dependent = sum(group_classes$effect_class == "context_dependent"),
    n_groups_consistent_deleterious = sum(group_classes$effect_class == "consistent_deleterious"),
    n_groups_consistent_neutral = sum(group_classes$effect_class == "consistent_neutral"),
    n_groups_consistent_beneficial = sum(group_classes$effect_class == "consistent_beneficial"),
    n_groups_using_global_fallback = sum(!is.na(tables$group_priors$fallback_group) &
                                           tables$group_priors$fallback_group == "global"),
    class_distribution = class_distribution,
    lambda_distribution = summary(tables$all_group_priors$cohort_lambda),
    prior_sd_distribution = summary(tables$all_group_priors$effective_prior_sd),
    groups_estimated = unique(tables$group_priors$group),
    groups_fallback_to_coarser = tables$group_priors$group[!is.na(tables$group_priors$fallback_group) &
                                                            tables$group_priors$fallback_group != "global"],
    groups_fallback_to_global = tables$group_priors$group[!is.na(tables$group_priors$fallback_group) &
                                                           tables$group_priors$fallback_group == "global"],
    mu_by_group = stats::setNames(tables$group_priors$mu, tables$group_priors$group),
    sigma_by_group = stats::setNames(tables$group_priors$effective_prior_sd, tables$group_priors$group),
    sigma_floor_used = any(tables$group_priors$sd_floor_used),
    patient_heterogeneity_sd = cohort_transition_patient_sd_floor,
    zero_likelihood_approximation = FALSE,
    zero_as_censoring_only = isTRUE(cohort_transition_zero_as_censoring_only),
    leave_one_patient_out_used = isTRUE(leave_one_patient_out),
    patients_contributing_by_group = lapply(unique(tables$group_priors$group), function(group_name) {
      sort(unique(summaries$patient_id[summaries$transition_group == group_name]))
    }),
    sd_floor = cohort_transition_sd_floor
  )
  names(diagnostics$patients_contributing_by_group) <- unique(tables$group_priors$group)

  list(
    version = "cohort_transition_v2",
    grouping = grouping,
    global_prior = tables$global_prior,
    group_priors = tables$group_priors,
    all_group_priors = tables$all_group_priors,
    patient_group_summaries = summaries,
    group_classes = group_classes,
    patient_shifts = patient_shifts,
    patient_ids = patient_ids,
    leave_one_patient_out = isTRUE(leave_one_patient_out),
    loo_priors = loo_priors,
    diagnostics = diagnostics,
    filter_diagnostics = filtered$diagnostics,
    excluded_records = filtered$excluded_records
  )
}

#' Learn a cohort-level transition-effect prior
#'
#' Version `"v2"` remains the lower-level default for backward compatibility and
#' aggregates bootstrap/path records to patient-level evidence before
#' classifying each transition group for selective borrowing. Version
#' `"contextual"` creates a context-aware evidence-bank prior on Delta fitness
#' conditioned on parent karyotype profile shape, copy-number area, CNA burden,
#' changed chromosome, local copy state, and event similarity. Version `"v1"`
#' preserves the original direct cohort-prior behavior for compatibility.
#'
#' @inheritParams extract_cohort_transition_records
#' @param records Transition records produced by `extract_cohort_transition_records()`.
#' @param leave_one_patient_out Whether to store/use leave-one-patient-out
#'   evidence for patient refits.
#' @param cohort_transition_version Prior-learning version, `"contextual"`,
#'   `"v2"`, or `"v1"`.
#' @param cohort_transition_min_patients_per_group Minimum patients required for
#'   a v1/v2 group prior before fallback.
#' @param cohort_transition_min_effective_n Minimum effective evidence for
#'   v1/v2 group priors.
#' @param cohort_transition_sd_floor Minimum transition-effect SD.
#' @param cohort_transition_patient_sd_floor Patient heterogeneity SD floor.
#' @param cohort_transition_global_fallback Whether unsupported groups can fall
#'   back to a broad global prior.
#' @param cohort_transition_zero_weight_cap_ratio Cap on zero-censoring evidence
#'   weight.
#' @param cohort_transition_zero_expected_count_cap Cap on zero expected-count
#'   influence in v1 compatibility fitting.
#' @param cohort_transition_zero_mean_shift_cap Maximum zero-censoring shift of
#'   a group mean in v1 compatibility fitting.
#' @return A cohort-transition prior object.
#' @export
learn_cohort_transition_prior <- function(records,
                                          leave_one_patient_out = TRUE,
                                          grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                          cohort_transition_version = c("v2", "contextual", "v1"),
                                          cohort_transition_min_patients_per_group = 2L,
                                          cohort_transition_min_effective_n = 3L,
                                          cohort_transition_sd_floor = 0.05,
                                          cohort_transition_patient_sd_floor = 0.10,
                                          cohort_transition_global_fallback = TRUE,
                                          cohort_transition_zero_weight_cap_ratio = 0.25,
                                          cohort_transition_zero_expected_count_cap = 10.0,
                                          cohort_transition_zero_mean_shift_cap = 0.2,
                                          ...) {
  grouping <- match.arg(grouping)
  cohort_transition_version <- match.arg(cohort_transition_version)
  if (identical(cohort_transition_version, "contextual")) {
    return(learn_cohort_transition_prior_contextual(
      records = records,
      leave_one_patient_out = leave_one_patient_out,
      grouping = grouping,
      ...
    ))
  }
  if (identical(cohort_transition_version, "v1")) {
    return(learn_cohort_transition_prior_v1(
      records = records,
      leave_one_patient_out = leave_one_patient_out,
      grouping = grouping,
      cohort_transition_min_patients_per_group = cohort_transition_min_patients_per_group,
      cohort_transition_min_effective_n = cohort_transition_min_effective_n,
      cohort_transition_sd_floor = cohort_transition_sd_floor,
      cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor,
      cohort_transition_global_fallback = cohort_transition_global_fallback,
      cohort_transition_zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio,
      cohort_transition_zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
      cohort_transition_zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap,
      ...
    ))
  }
  learn_cohort_transition_prior_v2(
    records = records,
    leave_one_patient_out = leave_one_patient_out,
    grouping = grouping,
    cohort_transition_min_patients_per_group = cohort_transition_min_patients_per_group,
    cohort_transition_min_effective_n = cohort_transition_min_effective_n,
    cohort_transition_sd_floor = cohort_transition_sd_floor,
    cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor,
    cohort_transition_global_fallback = cohort_transition_global_fallback,
    cohort_transition_zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio,
    cohort_transition_zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
    cohort_transition_zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap,
    ...
  )
}

resolve_cohort_transition_prior_object <- function(cohort_transition_prior = NULL,
                                                   cohort_transition_prior_path = NULL,
                                                   cohort_transition_patient_id = NULL) {
  if (!is.null(cohort_transition_prior) && !is.null(cohort_transition_prior_path)) {
    warning("Both `cohort_transition_prior` and `cohort_transition_prior_path` were supplied; using the object.", call. = FALSE)
  }
  prior <- cohort_transition_prior
  if (is.null(prior) && !is.null(cohort_transition_prior_path)) {
    prior <- readRDS(cohort_transition_prior_path)
  }
  if (is.null(prior)) {
    stop("`nn_prior = \"cohort_transition\"` requires `cohort_transition_prior` or `cohort_transition_prior_path`.", call. = FALSE)
  }
  if (!is.list(prior) || !(prior$version %in% c("cohort_transition_v1", "cohort_transition_v2", "cohort_transition_contextual_v1"))) {
    stop("`cohort_transition_prior` must be a cohort_transition_v1, cohort_transition_v2, or cohort_transition_contextual_v1 prior object.", call. = FALSE)
  }
  if (identical(prior$version, "cohort_transition_v1")) {
    warning("Using a cohort_transition_v1 prior object; v2 selective borrowing is preferred.", call. = FALSE)
  }
  if (isTRUE(prior$leave_one_patient_out) && length(prior$loo_priors)) {
    if (is.null(cohort_transition_patient_id) || length(cohort_transition_patient_id) != 1L || !nzchar(cohort_transition_patient_id)) {
      stop("`cohort_transition_patient_id` is required when the cohort transition prior contains leave-one-patient-out priors.", call. = FALSE)
    }
  }
  prior
}

cohort_transition_prior_for_patient <- function(prior, patient_id = NULL) {
  if (identical(prior$version, "cohort_transition_contextual_v1")) {
    if (isTRUE(prior$leave_one_patient_out) &&
        (is.null(patient_id) || length(patient_id) != 1L || !nzchar(patient_id))) {
      stop("`cohort_transition_patient_id` is required for leave-one-patient-out contextual cohort transition priors.", call. = FALSE)
    }
    prior$target_patient_id <- if (is.null(patient_id)) NA_character_ else as.character(patient_id)
    return(prior)
  }
  if (isTRUE(prior$leave_one_patient_out) && length(prior$loo_priors)) {
    if (is.null(patient_id) || !nzchar(patient_id)) {
      stop("`cohort_transition_patient_id` is required for leave-one-patient-out cohort transition priors.", call. = FALSE)
    }
    if (!patient_id %in% names(prior$loo_priors)) {
      if (!patient_id %in% prior$patient_ids && nrow(prior$global_prior)) {
        return(list(
          version = prior$version,
          grouping = prior$grouping,
          group_priors = prior$group_priors,
          global_prior = prior$global_prior,
          all_group_priors = prior$all_group_priors,
          group_classes = prior$group_classes,
          contributing_patients = prior$patient_ids,
          leave_one_patient_out = TRUE,
          leave_one_patient_out_fallback = "patient_has_no_training_records",
          patient_delta_shift = 0,
          patient_delta_shift_n_records = 0L,
          patient_delta_shift_reliability = 0
        ))
      }
      stop(sprintf("No leave-one-patient-out cohort transition prior is available for patient `%s`.", patient_id), call. = FALSE)
    }
    loo <- prior$loo_priors[[patient_id]]
    if (!nrow(loo$global_prior)) {
      stop(sprintf("Leave-one-patient-out prior for patient `%s` has no contributing records.", patient_id), call. = FALSE)
    }
    shift <- if (!is.null(prior$patient_shifts) && patient_id %in% names(prior$patient_shifts)) {
      prior$patient_shifts[[patient_id]]
    } else {
      data.frame(patient_delta_shift = 0, patient_delta_shift_n_records = 0L, patient_delta_shift_reliability = 0)
    }
    return(list(
      version = prior$version,
      grouping = prior$grouping,
      group_priors = loo$group_priors,
      global_prior = loo$global_prior,
      all_group_priors = loo$all_group_priors,
      group_classes = loo$group_classes,
      contributing_patients = loo$contributing_patients,
      leave_one_patient_out = TRUE,
      patient_delta_shift = shift$patient_delta_shift[1],
      patient_delta_shift_n_records = shift$patient_delta_shift_n_records[1],
      patient_delta_shift_reliability = shift$patient_delta_shift_reliability[1]
    ))
  }
  shift <- if (!is.null(prior$patient_shifts) && !is.null(patient_id) && patient_id %in% names(prior$patient_shifts)) {
    prior$patient_shifts[[patient_id]]
  } else {
    data.frame(patient_delta_shift = 0, patient_delta_shift_n_records = 0L, patient_delta_shift_reliability = 0)
  }
  list(
    version = prior$version,
    grouping = prior$grouping,
    group_priors = prior$group_priors,
    global_prior = prior$global_prior,
    all_group_priors = prior$all_group_priors,
    group_classes = prior$group_classes,
    contributing_patients = prior$patient_ids,
    leave_one_patient_out = FALSE,
    patient_delta_shift = shift$patient_delta_shift[1],
    patient_delta_shift_n_records = shift$patient_delta_shift_n_records[1],
    patient_delta_shift_reliability = shift$patient_delta_shift_reliability[1]
  )
}

lookup_cohort_transition_group_prior <- function(prior_use, parent_karyotype, child_karyotype) {
  parsed <- cohort_transition_parse_pair(parent_karyotype, child_karyotype)
  group <- parsed[[cohort_transition_group_column(prior_use$grouping)]]
  row <- prior_use$group_priors[prior_use$group_priors$group == group, , drop = FALSE]
  if (!nrow(row)) {
    row <- prior_use$global_prior
    row$group <- group
    row$fallback_group <- "global"
    row$fallback_multiplier <- 0.25
  }
  row <- row[1L, , drop = FALSE]
  fallback_group <- if ("fallback_group" %in% names(row)) row$fallback_group[1] else NA_character_
  if (!"cohort_lambda" %in% names(row)) row$cohort_lambda <- 1
  if (!"sd_multiplier" %in% names(row)) row$sd_multiplier <- 1
  if (!"effective_prior_sd" %in% names(row)) {
    row$effective_prior_sd <- if ("sigma_with_patient_heterogeneity" %in% names(row)) {
      row$sigma_with_patient_heterogeneity
    } else {
      row$sigma
    }
  }
  if (!"effect_class" %in% names(row)) row$effect_class <- "legacy_v1"
  if (!"heterogeneity_class" %in% names(row)) row$heterogeneity_class <- NA_character_
  if (!"recommended_use_for_zero" %in% names(row)) row$recommended_use_for_zero <- TRUE
  if (!"use_for_observed" %in% names(row)) row$use_for_observed <- TRUE
  if (!"use_for_low_information" %in% names(row)) row$use_for_low_information <- TRUE
  if (!"warning_flags" %in% names(row)) row$warning_flags <- ""
  if (!"fallback_multiplier" %in% names(row)) {
    row$fallback_multiplier <- if (!is.na(fallback_group) && identical(fallback_group, "global")) 0.25 else if (!is.na(fallback_group)) 0.5 else 1
  }
  list(
    group = group,
    prior_group_used = row$group[1],
    fallback_group_used = fallback_group,
    mu = row$mu[1],
    sd = row$effective_prior_sd[1],
    raw_sd = if ("sigma_with_patient_heterogeneity" %in% names(row)) row$sigma_with_patient_heterogeneity[1] else row$sigma[1],
    effect_class = row$effect_class[1],
    heterogeneity_class = row$heterogeneity_class[1],
    cohort_lambda = row$cohort_lambda[1],
    sd_multiplier = row$sd_multiplier[1],
    effective_prior_sd = row$effective_prior_sd[1],
    fallback_multiplier = row$fallback_multiplier[1],
    class_warning_flags = row$warning_flags[1],
    use_for_zero = isTRUE(row$recommended_use_for_zero[1]),
    use_for_observed = isTRUE(row$use_for_observed[1]),
    use_for_low_information = isTRUE(row$use_for_low_information[1]),
    parsed = parsed
  )
}

cohort_context_as_cn_vector <- function(karyotype) {
  if (is.numeric(karyotype)) {
    vec <- as.numeric(karyotype)
    if (!length(vec) || any(!is.finite(vec))) {
      stop("Numeric karyotypes must contain finite copy-number values.", call. = FALSE)
    }
    return(vec)
  }
  if (is.character(karyotype) && length(karyotype) == 1L && nzchar(karyotype)) {
    return(as.numeric(parse_karyotype_ids(karyotype)[1, ]))
  }
  stop("`karyotype` must be a single karyotype string or numeric copy-number vector.", call. = FALSE)
}

cohort_context_karyotype_label <- function(karyotype) {
  if (is.numeric(karyotype)) {
    paste(as.numeric(karyotype), collapse = ".")
  } else {
    as.character(karyotype)[1]
  }
}

cohort_context_profile_vector <- function(cn_vector,
                                          profile_transform = c("mass", "centered", "zscore", "raw"),
                                          chromosome_weights = NULL) {
  profile_transform <- match.arg(profile_transform)
  cn <- as.numeric(cn_vector)
  if (is.null(chromosome_weights)) {
    chromosome_weights <- rep(1, length(cn))
  }
  chromosome_weights <- as.numeric(chromosome_weights)
  if (length(chromosome_weights) != length(cn) || any(!is.finite(chromosome_weights)) || any(chromosome_weights < 0)) {
    stop("`chromosome_weights` must be a non-negative numeric vector aligned with the karyotype.", call. = FALSE)
  }
  weighted_cn <- cn * chromosome_weights
  if (identical(profile_transform, "mass")) {
    total <- sum(weighted_cn)
    if (!is.finite(total) || total <= 0) {
      return(rep(1 / length(cn), length(cn)))
    }
    return(weighted_cn / total)
  }
  if (identical(profile_transform, "centered")) {
    return(cn - mean(cn))
  }
  if (identical(profile_transform, "zscore")) {
    sd_cn <- stats::sd(cn)
    if (!is.finite(sd_cn) || sd_cn <= 0) {
      return(rep(0, length(cn)))
    }
    return((cn - mean(cn)) / sd_cn)
  }
  cn
}

cohort_context_gini <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (any(x < 0)) x <- x - min(x)
  sx <- sum(x)
  if (!is.finite(sx) || sx <= 0) return(0)
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sx)
}

#' Compute karyotype copy-number profile features
#'
#' @param karyotype Single karyotype string such as `"2.2.3"` or a numeric
#'   copy-number vector.
#' @param baseline_ploidy Copy-number baseline used for CNA burden.
#' @param chromosome_weights Optional non-negative chromosome weights.
#' @param profile_transform Profile representation used by contextual matching.
#' @return A one-row data frame with scalar fields and list-columns containing
#'   copy-number/profile vectors.
#' @export
compute_karyotype_profile_features <- function(karyotype,
                                               baseline_ploidy = 2,
                                               chromosome_weights = NULL,
                                               profile_transform = c("mass", "centered", "zscore", "raw")) {
  profile_transform <- match.arg(profile_transform)
  validate_positive_finite(baseline_ploidy, "baseline_ploidy")
  cn <- cohort_context_as_cn_vector(karyotype)
  label <- cohort_context_karyotype_label(karyotype)
  if (is.null(chromosome_weights)) {
    chromosome_weights <- rep(1, length(cn))
  }
  if (length(chromosome_weights) != length(cn)) {
    stop("`chromosome_weights` must have the same length as the karyotype vector.", call. = FALSE)
  }
  total_cn <- sum(cn)
  mass <- cohort_context_profile_vector(cn, "mass", chromosome_weights)
  centered <- cohort_context_profile_vector(cn, "centered", chromosome_weights)
  zscore <- cohort_context_profile_vector(cn, "zscore", chromosome_weights)
  transformed <- cohort_context_profile_vector(cn, profile_transform, chromosome_weights)
  entropy <- {
    p <- mass[mass > 0 & is.finite(mass)]
    if (length(p)) -sum(p * log(p)) else NA_real_
  }
  data.frame(
    karyotype = label,
    n_chr = length(cn),
    cn_vector = I(list(cn)),
    total_cn = total_cn,
    mean_cn = mean(cn),
    var_cn = if (length(cn) > 1L) stats::var(cn) else 0,
    sd_cn = if (length(cn) > 1L) stats::sd(cn) else 0,
    min_cn = min(cn),
    max_cn = max(cn),
    cna_burden = sum(abs(cn - baseline_ploidy)),
    n_gain_chr = sum(cn > baseline_ploidy),
    n_loss_chr = sum(cn < baseline_ploidy),
    profile_entropy = entropy,
    profile_gini = cohort_context_gini(cn),
    profile_mass_vector = I(list(mass)),
    centered_profile_vector = I(list(centered)),
    zscore_profile_vector = I(list(zscore)),
    context_profile_vector = I(list(transformed)),
    stringsAsFactors = FALSE
  )
}

#' Compute transition context features
#'
#' @param parent_karyotype Parent karyotype string or numeric vector.
#' @param child_karyotype Child karyotype string or numeric vector.
#' @param parent_fitness Optional patient-specific parent fitness.
#' @inheritParams compute_karyotype_profile_features
#' @return A one-row data frame with transition and background features.
#' @export
compute_transition_context_features <- function(parent_karyotype,
                                                child_karyotype,
                                                parent_fitness = NA_real_,
                                                baseline_ploidy = 2,
                                                chromosome_weights = NULL,
                                                profile_transform = "mass") {
  parent_vec <- cohort_context_as_cn_vector(parent_karyotype)
  child_vec <- cohort_context_as_cn_vector(child_karyotype)
  if (length(parent_vec) != length(child_vec)) {
    stop("Parent and child karyotypes must have matching dimensions.", call. = FALSE)
  }
  parent_label <- cohort_context_karyotype_label(parent_karyotype)
  child_label <- cohort_context_karyotype_label(child_karyotype)
  parent_feat <- compute_karyotype_profile_features(
    parent_vec,
    baseline_ploidy = baseline_ploidy,
    chromosome_weights = chromosome_weights,
    profile_transform = profile_transform
  )
  child_feat <- compute_karyotype_profile_features(
    child_vec,
    baseline_ploidy = baseline_ploidy,
    chromosome_weights = chromosome_weights,
    profile_transform = profile_transform
  )
  diff_vec <- child_vec - parent_vec
  changed <- which(diff_vec != 0)
  transition_chr <- if (length(changed) == 1L) changed else NA_integer_
  transition_direction <- if (length(changed) == 1L && diff_vec[changed] > 0) {
    "gain"
  } else if (length(changed) == 1L && diff_vec[changed] < 0) {
    "loss"
  } else {
    "complex"
  }
  transition_size <- sum(abs(diff_vec))
  parent_z <- parent_feat$zscore_profile_vector[[1]]
  changed_parent_copy <- if (length(changed) == 1L) parent_vec[changed] else NA_real_
  changed_child_copy <- if (length(changed) == 1L) child_vec[changed] else NA_real_
  changed_parent_z <- if (length(changed) == 1L) parent_z[changed] else NA_real_
  changed_child_z <- if (length(changed) == 1L) child_feat$zscore_profile_vector[[1]][changed] else NA_real_
  data.frame(
    parent_karyotype = parent_label,
    child_karyotype = child_label,
    transition_chr = transition_chr,
    transition_direction = transition_direction,
    transition_size = transition_size,
    is_one_step = isTRUE(length(changed) == 1L && abs(diff_vec[changed]) == 1),
    parent_total_cn = parent_feat$total_cn,
    child_total_cn = child_feat$total_cn,
    delta_total_cn = child_feat$total_cn - parent_feat$total_cn,
    parent_burden = parent_feat$cna_burden,
    child_burden = child_feat$cna_burden,
    delta_burden = child_feat$cna_burden - parent_feat$cna_burden,
    parent_mean_cn = parent_feat$mean_cn,
    parent_sd_cn = parent_feat$sd_cn,
    parent_cna_burden = parent_feat$cna_burden,
    child_cna_burden = child_feat$cna_burden,
    parent_cn_vector = I(list(parent_vec)),
    child_cn_vector = I(list(child_vec)),
    parent_profile_mass = I(list(parent_feat$profile_mass_vector[[1]])),
    child_profile_mass = I(list(child_feat$profile_mass_vector[[1]])),
    parent_context_profile = I(list(parent_feat$context_profile_vector[[1]])),
    child_context_profile = I(list(child_feat$context_profile_vector[[1]])),
    changed_chr_parent_copy = changed_parent_copy,
    changed_chr_child_copy = changed_child_copy,
    changed_chr_parent_zscore = changed_parent_z,
    changed_chr_child_zscore = changed_child_z,
    changed_chr_is_peak = isTRUE(length(changed) == 1L && parent_vec[changed] == max(parent_vec)),
    changed_chr_is_valley = isTRUE(length(changed) == 1L && parent_vec[changed] == min(parent_vec)),
    parent_fitness = parent_fitness,
    parent_fitness_bin = if (is.finite(parent_fitness)) {
      cut(parent_fitness, breaks = c(-Inf, -0.05, 0.05, Inf), labels = c("low", "near_zero", "high"))[1]
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )
}

#' Compute distance between karyotype profile vectors
#'
#' @param profile_a,profile_b Numeric profile vectors.
#' @param method Distance method.
#' @param chromosome_weights Optional non-negative weights.
#' @return A finite non-negative distance; Hellinger and Jensen-Shannon are
#'   computed on normalized mass profiles.
#' @export
karyotype_profile_distance <- function(profile_a,
                                       profile_b,
                                       method = c("hellinger", "jensen_shannon", "cosine", "euclidean", "manhattan"),
                                       chromosome_weights = NULL) {
  method <- match.arg(method)
  a <- as.numeric(profile_a)
  b <- as.numeric(profile_b)
  if (length(a) != length(b)) {
    return(Inf)
  }
  if (is.null(chromosome_weights)) chromosome_weights <- rep(1, length(a))
  chromosome_weights <- as.numeric(chromosome_weights)
  if (length(chromosome_weights) != length(a)) {
    stop("`chromosome_weights` must match the profile length.", call. = FALSE)
  }
  ok <- is.finite(a) & is.finite(b) & is.finite(chromosome_weights) & chromosome_weights >= 0
  if (!all(ok)) {
    a <- a[ok]
    b <- b[ok]
    chromosome_weights <- chromosome_weights[ok]
  }
  if (!length(a)) return(Inf)
  normalize_mass <- function(x) {
    x <- pmax(0, x * chromosome_weights)
    sx <- sum(x)
    if (!is.finite(sx) || sx <= 0) rep(1 / length(x), length(x)) else x / sx
  }
  d <- switch(
    method,
    hellinger = {
      pa <- normalize_mass(a)
      pb <- normalize_mass(b)
      sqrt(sum((sqrt(pa) - sqrt(pb))^2) / 2)
    },
    jensen_shannon = {
      eps <- 1e-12
      pa <- pmax(normalize_mass(a), eps)
      pb <- pmax(normalize_mass(b), eps)
      pa <- pa / sum(pa)
      pb <- pb / sum(pb)
      m <- 0.5 * (pa + pb)
      sqrt(0.5 * sum(pa * log(pa / m)) + 0.5 * sum(pb * log(pb / m)))
    },
    cosine = {
      aw <- a * chromosome_weights
      bw <- b * chromosome_weights
      denom <- sqrt(sum(aw^2)) * sqrt(sum(bw^2))
      if (!is.finite(denom) || denom <= 0) 1 else max(0, 1 - sum(aw * bw) / denom)
    },
    euclidean = sqrt(sum(chromosome_weights * (a - b)^2)),
    manhattan = sum(chromosome_weights * abs(a - b))
  )
  if (!is.finite(d) || d < 0) Inf else d
}

#' Compute distance between CNA transition events
#'
#' @param context_a,context_b One-row transition context data frames from
#'   `compute_transition_context_features()`.
#' @param event_match Event matching rule.
#' @return A non-negative distance, or `Inf` when a hard event-match rule fails.
#' @export
transition_event_distance <- function(context_a,
                                      context_b,
                                      event_match = c("same_chr_direction", "same_direction", "kernel")) {
  event_match <- match.arg(event_match)
  getv <- function(x, name) {
    if (is.data.frame(x)) x[[name]][1] else x[[name]]
  }
  chr_a <- getv(context_a, "transition_chr")
  chr_b <- getv(context_b, "transition_chr")
  dir_a <- as.character(getv(context_a, "transition_direction"))
  dir_b <- as.character(getv(context_b, "transition_direction"))
  if (identical(event_match, "same_chr_direction") &&
      (!isTRUE(is.finite(chr_a) && is.finite(chr_b) && chr_a == chr_b) || !identical(dir_a, dir_b))) {
    return(Inf)
  }
  if (identical(event_match, "same_direction") && !identical(dir_a, dir_b)) {
    return(Inf)
  }
  chr_penalty <- if (is.finite(chr_a) && is.finite(chr_b) && chr_a == chr_b) 0 else 1
  direction_penalty <- if (identical(dir_a, dir_b)) 0 else 2
  size_diff <- abs(as.numeric(getv(context_a, "transition_size")) - as.numeric(getv(context_b, "transition_size")))
  local_copy_diff <- abs(as.numeric(getv(context_a, "changed_chr_parent_copy")) - as.numeric(getv(context_b, "changed_chr_parent_copy")))
  local_z_diff <- abs(as.numeric(getv(context_a, "changed_chr_parent_zscore")) - as.numeric(getv(context_b, "changed_chr_parent_zscore")))
  delta_area_diff <- abs(as.numeric(getv(context_a, "delta_total_cn")) - as.numeric(getv(context_b, "delta_total_cn")))
  delta_burden_diff <- abs(as.numeric(getv(context_a, "delta_burden")) - as.numeric(getv(context_b, "delta_burden")))
  vals <- c(chr_penalty, direction_penalty, size_diff, local_copy_diff, local_z_diff, delta_area_diff, delta_burden_diff)
  vals[!is.finite(vals)] <- 0
  sqrt(sum(vals^2))
}

cohort_context_direction_code <- function(x) {
  out <- rep.int(0L, length(x))
  x <- as.character(x)
  out[x == "gain"] <- 1L
  out[x == "loss"] <- -1L
  out
}

cohort_context_profile_distance_code <- function(method) {
  switch(
    match.arg(method, c("hellinger", "jensen_shannon", "cosine", "euclidean", "manhattan")),
    hellinger = 1L,
    jensen_shannon = 2L,
    cosine = 3L,
    euclidean = 4L,
    manhattan = 5L
  )
}

cohort_context_event_match_code <- function(event_match) {
  switch(
    match.arg(event_match, c("same_chr_direction", "same_direction", "kernel")),
    same_chr_direction = 1L,
    same_direction = 2L,
    kernel = 3L
  )
}

cohort_context_numeric_cache <- function(evidence_contexts) {
  if (!is.data.frame(evidence_contexts) || !nrow(evidence_contexts)) {
    return(NULL)
  }
  if (!"parent_context_profile" %in% names(evidence_contexts)) {
    stop("Context evidence is missing `parent_context_profile`.", call. = FALSE)
  }
  profile_matrix <- do.call(rbind, lapply(evidence_contexts$parent_context_profile, as.numeric))
  if (is.null(dim(profile_matrix))) {
    profile_matrix <- matrix(profile_matrix, nrow = nrow(evidence_contexts))
  }
  list(
    profile_matrix = profile_matrix,
    total_cn = as.numeric(evidence_contexts$parent_total_cn),
    burden = as.numeric(evidence_contexts$parent_burden),
    local_copy = as.numeric(evidence_contexts$changed_chr_parent_copy),
    local_z = as.numeric(evidence_contexts$changed_chr_parent_zscore),
    transition_chr = as.integer(evidence_contexts$transition_chr),
    direction_code = cohort_context_direction_code(evidence_contexts$transition_direction),
    transition_size = as.numeric(evidence_contexts$transition_size),
    delta_total_cn = as.numeric(evidence_contexts$delta_total_cn),
    delta_burden = as.numeric(evidence_contexts$delta_burden),
    quality_weight = if ("quality_weight" %in% names(evidence_contexts)) {
      as.numeric(evidence_contexts$quality_weight)
    } else {
      rep(1, nrow(evidence_contexts))
    }
  )
}

cohort_context_attach_numeric_cache <- function(evidence_contexts) {
  if (is.data.frame(evidence_contexts) && nrow(evidence_contexts)) {
    attr(evidence_contexts, "context_numeric_cache") <- cohort_context_numeric_cache(evidence_contexts)
  }
  evidence_contexts
}

cohort_context_subset_cache <- function(cache, keep) {
  if (is.null(cache)) return(NULL)
  keep_idx <- which(keep)
  if (!length(keep_idx)) return(NULL)
  list(
    profile_matrix = cache$profile_matrix[keep_idx, , drop = FALSE],
    total_cn = cache$total_cn[keep_idx],
    burden = cache$burden[keep_idx],
    local_copy = cache$local_copy[keep_idx],
    local_z = cache$local_z[keep_idx],
    transition_chr = cache$transition_chr[keep_idx],
    direction_code = cache$direction_code[keep_idx],
    transition_size = cache$transition_size[keep_idx],
    delta_total_cn = cache$delta_total_cn[keep_idx],
    delta_burden = cache$delta_burden[keep_idx],
    quality_weight = cache$quality_weight[keep_idx]
  )
}

cohort_context_subset_bank <- function(evidence_contexts, keep) {
  out <- evidence_contexts[keep, , drop = FALSE]
  cache <- attr(evidence_contexts, "context_numeric_cache", exact = TRUE)
  if (!is.null(cache)) {
    attr(out, "context_numeric_cache") <- cohort_context_subset_cache(cache, keep)
  }
  out
}

cohort_context_bandwidth_vector <- function(bandwidths) {
  c(
    profile = as.numeric(bandwidths$profile),
    area = as.numeric(bandwidths$area),
    burden = as.numeric(bandwidths$burden),
    local = as.numeric(bandwidths$local),
    event = as.numeric(bandwidths$event)
  )
}

cohort_context_weight_vector <- function(weights) {
  c(
    profile = as.numeric(weights$profile),
    area = as.numeric(weights$area),
    burden = as.numeric(weights$burden),
    local = as.numeric(weights$local),
    event = as.numeric(weights$event)
  )
}

#' Compute contextual transition distance
#'
#' @param target_context,evidence_context One-row transition context data frames.
#' @param bandwidths List with profile, area, burden, local, and event bandwidths.
#' @param weights List with profile, area, burden, local, and event weights.
#' @param profile_distance Profile distance method.
#' @param event_match Event matching rule.
#' @param chromosome_weights Optional chromosome weights.
#' @return A list of total and component distances.
#' @export
compute_context_distance <- function(target_context,
                                     evidence_context,
                                     bandwidths,
                                     weights,
                                     profile_distance = c("hellinger", "jensen_shannon", "cosine", "euclidean", "manhattan"),
                                     event_match = c("same_chr_direction", "same_direction", "kernel"),
                                     chromosome_weights = NULL) {
  profile_distance <- match.arg(profile_distance)
  event_match <- match.arg(event_match)
  getv <- function(x, name) {
    if (is.data.frame(x)) x[[name]][1] else x[[name]]
  }
  getlist <- function(x, name) {
    val <- getv(x, name)
    if (is.list(val)) val[[1]] else val
  }
  d_profile <- karyotype_profile_distance(
    getlist(target_context, "parent_context_profile"),
    getlist(evidence_context, "parent_context_profile"),
    method = profile_distance,
    chromosome_weights = chromosome_weights
  )
  d_area <- abs(as.numeric(getv(target_context, "parent_total_cn")) - as.numeric(getv(evidence_context, "parent_total_cn")))
  d_burden <- abs(as.numeric(getv(target_context, "parent_burden")) - as.numeric(getv(evidence_context, "parent_burden")))
  d_local <- sqrt(sum(c(
    as.numeric(getv(target_context, "changed_chr_parent_copy")) - as.numeric(getv(evidence_context, "changed_chr_parent_copy")),
    as.numeric(getv(target_context, "changed_chr_parent_zscore")) - as.numeric(getv(evidence_context, "changed_chr_parent_zscore"))
  )^2, na.rm = TRUE))
  d_event <- transition_event_distance(target_context, evidence_context, event_match = event_match)
  scale_one <- function(d, bw) {
    if (!is.finite(d)) return(Inf)
    if (!is.finite(bw) || bw <= 0) bw <- 1
    d / bw
  }
  components <- c(
    profile = scale_one(d_profile, bandwidths$profile),
    area = scale_one(d_area, bandwidths$area),
    burden = scale_one(d_burden, bandwidths$burden),
    local = scale_one(d_local, bandwidths$local),
    event = scale_one(d_event, bandwidths$event)
  )
  weight_vec <- c(
    profile = weights$profile,
    area = weights$area,
    burden = weights$burden,
    local = weights$local,
    event = weights$event
  )
  weight_vec[!is.finite(weight_vec) | weight_vec < 0] <- 0
  total <- sqrt(sum(weight_vec * components^2))
  list(
    context_distance = if (is.finite(total)) total else Inf,
    profile_distance = d_profile,
    area_distance = d_area,
    burden_distance = d_burden,
    local_distance = d_local,
    event_distance = d_event
  )
}

#' Compute contextual kernel weights
#'
#' @param target_context One-row target transition context.
#' @param evidence_contexts Evidence bank data frame with context columns.
#' @inheritParams compute_context_distance
#' @param k_nearest Maximum number of nonzero neighbors to retain.
#' @param min_kernel_weight Minimum final kernel weight.
#' @return A data frame of nearest contextual evidence weights.
#' @export
compute_context_kernel_weights <- function(target_context,
                                           evidence_contexts,
                                           bandwidths,
                                           weights,
                                           event_match = c("same_chr_direction", "same_direction", "kernel"),
                                           k_nearest = 50,
                                           min_kernel_weight = 1e-6,
                                           profile_distance = "hellinger",
                                           chromosome_weights = NULL) {
  event_match <- match.arg(event_match)
  validate_positive_integer(k_nearest, "k_nearest")
  validate_nonnegative_finite(min_kernel_weight, "min_kernel_weight")
  if (!is.data.frame(evidence_contexts) || !nrow(evidence_contexts)) {
    return(data.frame())
  }
  cache <- attr(evidence_contexts, "context_numeric_cache", exact = TRUE)
  if (is.null(cache)) {
    cache <- cohort_context_numeric_cache(evidence_contexts)
  }
  getv <- function(x, name) {
    if (is.data.frame(x)) x[[name]][1] else x[[name]]
  }
  getlist <- function(x, name) {
    val <- getv(x, name)
    if (is.list(val)) val[[1]] else val
  }
  cpp_out <- tryCatch(
    context_kernel_weights_cpp(
      target_profile = as.numeric(getlist(target_context, "parent_context_profile")),
      target_total_cn = as.numeric(getv(target_context, "parent_total_cn")),
      target_burden = as.numeric(getv(target_context, "parent_burden")),
      target_local_copy = as.numeric(getv(target_context, "changed_chr_parent_copy")),
      target_local_z = as.numeric(getv(target_context, "changed_chr_parent_zscore")),
      target_transition_chr = as.integer(getv(target_context, "transition_chr")),
      target_direction_code = cohort_context_direction_code(getv(target_context, "transition_direction")),
      target_transition_size = as.numeric(getv(target_context, "transition_size")),
      target_delta_total_cn = as.numeric(getv(target_context, "delta_total_cn")),
      target_delta_burden = as.numeric(getv(target_context, "delta_burden")),
      evidence_profile_matrix = cache$profile_matrix,
      evidence_total_cn = cache$total_cn,
      evidence_burden = cache$burden,
      evidence_local_copy = cache$local_copy,
      evidence_local_z = cache$local_z,
      evidence_transition_chr = cache$transition_chr,
      evidence_direction_code = cache$direction_code,
      evidence_transition_size = cache$transition_size,
      evidence_delta_total_cn = cache$delta_total_cn,
      evidence_delta_burden = cache$delta_burden,
      quality_weight = cache$quality_weight,
      bandwidths = cohort_context_bandwidth_vector(bandwidths),
      component_weights = cohort_context_weight_vector(weights),
      chromosome_weights = if (is.null(chromosome_weights)) numeric(0) else as.numeric(chromosome_weights),
      event_match_code = cohort_context_event_match_code(event_match),
      profile_distance_code = cohort_context_profile_distance_code(profile_distance),
      k_nearest = as.integer(k_nearest),
      min_kernel_weight = min_kernel_weight
    ),
    error = function(e) NULL
  )
  if (is.data.frame(cpp_out)) {
    if (!nrow(cpp_out)) {
      return(data.frame(
        evidence_row_id = character(0),
        patient_id = character(0),
        context_distance = numeric(0),
        profile_distance = numeric(0),
        area_distance = numeric(0),
        burden_distance = numeric(0),
        local_distance = numeric(0),
        event_distance = numeric(0),
        kernel_weight = numeric(0),
        quality_weight = numeric(0),
        final_weight = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    idx <- as.integer(cpp_out$evidence_index)
    cpp_out$evidence_row_id <- if ("evidence_id" %in% names(evidence_contexts)) evidence_contexts$evidence_id[idx] else idx
    cpp_out$patient_id <- if ("patient_id" %in% names(evidence_contexts)) as.character(evidence_contexts$patient_id[idx]) else NA_character_
    cpp_out <- cpp_out[
      c(
        "evidence_row_id",
        "patient_id",
        "context_distance",
        "profile_distance",
        "area_distance",
        "burden_distance",
        "local_distance",
        "event_distance",
        "kernel_weight",
        "quality_weight",
        "final_weight"
      )
    ]
    rownames(cpp_out) <- NULL
    return(cpp_out)
  }
  rows <- lapply(seq_len(nrow(evidence_contexts)), function(i) {
    ev <- evidence_contexts[i, , drop = FALSE]
    dist <- compute_context_distance(
      target_context = target_context,
      evidence_context = ev,
      bandwidths = bandwidths,
      weights = weights,
      profile_distance = profile_distance,
      event_match = event_match,
      chromosome_weights = chromosome_weights
    )
    kernel_weight <- if (is.finite(dist$context_distance)) exp(-0.5 * dist$context_distance^2) else 0
    quality_weight <- if ("quality_weight" %in% names(ev)) ev$quality_weight[1] else 1
    if (!is.finite(quality_weight) || quality_weight < 0) quality_weight <- 0
    data.frame(
      evidence_row_id = if ("evidence_id" %in% names(ev)) ev$evidence_id[1] else i,
      patient_id = if ("patient_id" %in% names(ev)) as.character(ev$patient_id[1]) else NA_character_,
      context_distance = dist$context_distance,
      profile_distance = dist$profile_distance,
      area_distance = dist$area_distance,
      burden_distance = dist$burden_distance,
      local_distance = dist$local_distance,
      event_distance = dist$event_distance,
      kernel_weight = kernel_weight,
      quality_weight = quality_weight,
      final_weight = kernel_weight * quality_weight,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[is.finite(out$final_weight) & out$final_weight >= min_kernel_weight, , drop = FALSE]
  if (!nrow(out)) return(out)
  out <- out[order(out$final_weight, decreasing = TRUE), , drop = FALSE]
  utils::head(out, k_nearest)
}

cohort_context_enrich_record <- function(record,
                                         evidence_id,
                                         baseline_ploidy = 2,
                                         chromosome_weights = NULL,
                                         profile_transform = "mass") {
  ctx <- compute_transition_context_features(
    parent_karyotype = record$parent_karyotype[1],
    child_karyotype = record$child_karyotype[1],
    parent_fitness = if ("parent_fitness" %in% names(record)) record$parent_fitness[1] else NA_real_,
    baseline_ploidy = baseline_ploidy,
    chromosome_weights = chromosome_weights,
    profile_transform = profile_transform
  )
  for (nm in names(ctx)) {
    if (!nm %in% names(record)) {
      record[[nm]] <- ctx[[nm]]
    }
  }
  record$evidence_id <- evidence_id
  record
}

#' Build contextual transition evidence banks
#'
#' @param records Transition records from upstream two-shell results.
#' @inheritParams compute_karyotype_profile_features
#' @param cohort_transition_use_prior_dominated_records Whether prior-dominated
#'   observed records can enter the observed evidence bank.
#' @param cohort_transition_use_boundary_records Whether boundary records can
#'   enter the observed evidence bank.
#' @param cohort_transition_max_delta_se,cohort_transition_max_delta_se_quantile
#'   Delta-SE filters for observed evidence.
#' @param cohort_transition_min_path_responsibility Minimum path responsibility.
#' @param cohort_transition_min_observed_count Minimum observed child count.
#' @param cohort_context_zero_min_expected_count Minimum parent-like expected
#'   count for zero censoring evidence.
#' @param cohort_context_zero_as_censoring_only Treat zeros as censoring only.
#' @param ... Additional controls passed through for compatibility.
#' @return A list with `evidence_bank`, `zero_evidence_bank`, and diagnostics.
#' @export
build_contextual_transition_evidence_bank <- function(records,
                                                      baseline_ploidy = 2,
                                                      chromosome_weights = NULL,
                                                      profile_transform = "mass",
                                                      cohort_transition_use_prior_dominated_records = FALSE,
                                                      cohort_transition_use_boundary_records = FALSE,
                                                      cohort_transition_max_delta_se = NULL,
                                                      cohort_transition_max_delta_se_quantile = 0.75,
                                                      cohort_transition_min_path_responsibility = 0.05,
                                                      cohort_transition_min_observed_count = 1L,
                                                      cohort_context_zero_min_expected_count = 3.0,
                                                      cohort_context_zero_as_censoring_only = TRUE,
                                                      ...) {
  if (!is.data.frame(records) || !nrow(records)) {
    return(list(evidence_bank = data.frame(), zero_evidence_bank = data.frame(), diagnostics = list(n_input_records = 0L)))
  }
  filtered <- filter_cohort_transition_records(
    records,
    cohort_transition_use_prior_dominated_records = cohort_transition_use_prior_dominated_records,
    cohort_transition_use_boundary_records = cohort_transition_use_boundary_records,
    cohort_transition_max_delta_se = cohort_transition_max_delta_se,
    cohort_transition_max_delta_se_quantile = cohort_transition_max_delta_se_quantile,
    cohort_transition_min_path_responsibility = cohort_transition_min_path_responsibility,
    cohort_transition_min_observed_count = cohort_transition_min_observed_count,
    cohort_transition_zero_min_expected_count = cohort_context_zero_min_expected_count,
    cohort_transition_zero_as_censoring_only = cohort_context_zero_as_censoring_only
  )
  kept <- filtered$kept_records
  if (!nrow(kept)) {
    return(list(evidence_bank = data.frame(), zero_evidence_bank = data.frame(), diagnostics = filtered$diagnostics))
  }
  observed <- kept[kept$cohort_transition_evidence_type == "observed_delta_evidence" &
                     !as.logical(kept$child_is_zero) &
                     is.finite(kept$delta_hat), , drop = FALSE]
  zeros <- kept[kept$cohort_transition_evidence_type == "zero_censoring_evidence", , drop = FALSE]
  enrich <- function(df, prefix) {
    if (!nrow(df)) return(data.frame())
    rows <- vector("list", nrow(df))
    for (i in seq_len(nrow(df))) {
      rows[[i]] <- cohort_context_enrich_record(
        df[i, , drop = FALSE],
        evidence_id = paste0(prefix, "_", i),
        baseline_ploidy = baseline_ploidy,
        chromosome_weights = chromosome_weights,
        profile_transform = profile_transform
      )
    }
    out <- do.call(rbind, rows)
    if (!"quality_weight" %in% names(out)) {
      se <- pmax(as.numeric(out$delta_se), 0.05)
      pr <- as.numeric(out$path_responsibility)
      pr[!is.finite(pr) | pr < 0] <- 0
      out$quality_weight <- pmax(0, pr) / (1 + se^2)
      out$quality_weight[!is.finite(out$quality_weight)] <- 0
    }
    rownames(out) <- NULL
    out
  }
  evidence_bank <- enrich(observed, "obs")
  zero_evidence_bank <- enrich(zeros, "zero")
  evidence_bank <- cohort_context_attach_numeric_cache(evidence_bank)
  zero_evidence_bank <- cohort_context_attach_numeric_cache(zero_evidence_bank)
  diagnostics <- c(
    filtered$diagnostics,
    list(
      n_context_observed_evidence = nrow(evidence_bank),
      n_context_zero_evidence = nrow(zero_evidence_bank),
      n_patients_in_evidence_bank = length(unique(evidence_bank$patient_id)),
      evidence_source_type_counts = if (nrow(evidence_bank)) table(evidence_bank$source_type) else integer(0),
      evidence_delta_summary = if (nrow(evidence_bank)) summary(evidence_bank$delta_hat) else summary(numeric(0)),
      evidence_delta_se_summary = if (nrow(evidence_bank)) summary(evidence_bank$delta_se) else summary(numeric(0))
    )
  )
  list(evidence_bank = evidence_bank, zero_evidence_bank = zero_evidence_bank, diagnostics = diagnostics, filtered = filtered)
}

cohort_context_pairwise_values <- function(evidence_bank, fun) {
  if (!is.data.frame(evidence_bank) || nrow(evidence_bank) < 2L) return(numeric(0))
  vals <- numeric(0)
  for (i in seq_len(nrow(evidence_bank) - 1L)) {
    for (j in seq.int(i + 1L, nrow(evidence_bank))) {
      vals <- c(vals, fun(evidence_bank[i, , drop = FALSE], evidence_bank[j, , drop = FALSE]))
    }
  }
  vals[is.finite(vals) & vals > 0]
}

estimate_context_bandwidths <- function(evidence_bank,
                                        profile_distance = "hellinger",
                                        chromosome_weights = NULL,
                                        cohort_context_bandwidth_profile = NULL,
                                        cohort_context_bandwidth_area = NULL,
                                        cohort_context_bandwidth_burden = NULL,
                                        cohort_context_bandwidth_local = NULL,
                                        cohort_context_bandwidth_event = NULL) {
  robust_bw <- function(x, fallback = 1) {
    x <- x[is.finite(x) & x > 0]
    if (!length(x)) return(fallback)
    val <- as.numeric(stats::quantile(x, probs = 0.5, names = FALSE, type = 8))
    if (!is.finite(val) || val <= 0) fallback else val
  }
  profile_vals <- cohort_context_pairwise_values(evidence_bank, function(a, b) {
    karyotype_profile_distance(a$parent_context_profile[[1]], b$parent_context_profile[[1]], method = profile_distance, chromosome_weights = chromosome_weights)
  })
  area_vals <- cohort_context_pairwise_values(evidence_bank, function(a, b) abs(a$parent_total_cn[1] - b$parent_total_cn[1]))
  burden_vals <- cohort_context_pairwise_values(evidence_bank, function(a, b) abs(a$parent_burden[1] - b$parent_burden[1]))
  local_vals <- cohort_context_pairwise_values(evidence_bank, function(a, b) {
    sqrt(sum(c(a$changed_chr_parent_copy[1] - b$changed_chr_parent_copy[1],
               a$changed_chr_parent_zscore[1] - b$changed_chr_parent_zscore[1])^2, na.rm = TRUE))
  })
  event_vals <- cohort_context_pairwise_values(evidence_bank, function(a, b) transition_event_distance(a, b, event_match = "kernel"))
  list(
    profile = if (is.null(cohort_context_bandwidth_profile)) robust_bw(profile_vals, 0.25) else cohort_context_bandwidth_profile,
    area = if (is.null(cohort_context_bandwidth_area)) robust_bw(area_vals, 2) else cohort_context_bandwidth_area,
    burden = if (is.null(cohort_context_bandwidth_burden)) robust_bw(burden_vals, 2) else cohort_context_bandwidth_burden,
    local = if (is.null(cohort_context_bandwidth_local)) robust_bw(local_vals, 1) else cohort_context_bandwidth_local,
    event = if (is.null(cohort_context_bandwidth_event)) robust_bw(event_vals, 1) else cohort_context_bandwidth_event
  )
}

cohort_context_class_borrowing <- function(effect_class,
                                           cohort_context_lambda_consistent_deleterious = 0.50,
                                           cohort_context_lambda_consistent_neutral = 0.25,
                                           cohort_context_lambda_consistent_beneficial = 0.10,
                                           cohort_context_lambda_high_variable = 0.00,
                                           cohort_context_lambda_sparse_unknown = 0.00,
                                           cohort_context_lambda_conflicting_zero = 0.00,
                                           cohort_context_sd_multiplier_consistent_deleterious = 1.0,
                                           cohort_context_sd_multiplier_consistent_neutral = 1.5,
                                           cohort_context_sd_multiplier_consistent_beneficial = 2.5,
                                           cohort_context_sd_multiplier_high_variable = 4.0,
                                           cohort_context_sd_multiplier_sparse_unknown = 4.0,
                                           cohort_context_sd_multiplier_conflicting_zero = 4.0) {
  lambda <- switch(
    effect_class,
    context_consistent_deleterious = cohort_context_lambda_consistent_deleterious,
    context_consistent_neutral = cohort_context_lambda_consistent_neutral,
    context_consistent_beneficial = cohort_context_lambda_consistent_beneficial,
    context_high_variable = cohort_context_lambda_high_variable,
    context_sparse_unknown = cohort_context_lambda_sparse_unknown,
    context_conflicting_zero = cohort_context_lambda_conflicting_zero,
    cohort_context_lambda_sparse_unknown
  )
  sd_multiplier <- switch(
    effect_class,
    context_consistent_deleterious = cohort_context_sd_multiplier_consistent_deleterious,
    context_consistent_neutral = cohort_context_sd_multiplier_consistent_neutral,
    context_consistent_beneficial = cohort_context_sd_multiplier_consistent_beneficial,
    context_high_variable = cohort_context_sd_multiplier_high_variable,
    context_sparse_unknown = cohort_context_sd_multiplier_sparse_unknown,
    context_conflicting_zero = cohort_context_sd_multiplier_conflicting_zero,
    cohort_context_sd_multiplier_sparse_unknown
  )
  list(lambda = lambda, sd_multiplier = sd_multiplier)
}

#' Classify a target-specific contextual transition prior
#'
#' @param patient_evidence Patient-level weighted context evidence.
#' @param zero_neighbors Optional similar zero-censoring neighbors.
#' @param cohort_context_min_patients,cohort_context_min_effective_n,cohort_context_min_effective_patients,cohort_context_min_unique_children
#'   Minimum contextual support thresholds.
#' @param cohort_context_effect_threshold,cohort_context_sign_consistency_threshold,cohort_context_high_weighted_sd,cohort_context_high_between_patient_sd,cohort_context_high_i2
#'   Context classification thresholds.
#' @param cohort_context_sd_floor Minimum contextual prior SD.
#' @param ... Class-specific contextual borrowing and SD multiplier controls.
#' @return A one-row data frame with class and heterogeneity metrics.
#' @export
classify_contextual_transition_prior <- function(patient_evidence,
                                                 zero_neighbors = NULL,
                                                 cohort_context_min_patients = 3L,
                                                 cohort_context_min_effective_n = 5,
                                                 cohort_context_min_effective_patients = 3,
                                                 cohort_context_min_unique_children = 3L,
                                                 cohort_context_effect_threshold = 0.02,
                                                 cohort_context_sign_consistency_threshold = 0.75,
                                                 cohort_context_high_weighted_sd = 0.10,
                                                 cohort_context_high_between_patient_sd = 0.10,
                                                 cohort_context_high_i2 = 0.50,
                                                 cohort_context_sd_floor = 0.05,
                                                 ...) {
  if (!is.data.frame(patient_evidence) || !nrow(patient_evidence)) {
    class <- "context_sparse_unknown"
    borrow <- cohort_context_class_borrowing(class, ...)
    return(data.frame(
      context_delta_mu = 0,
      context_delta_sd = cohort_context_sd_floor * borrow$sd_multiplier,
      context_effective_n = 0,
      context_n_patients = 0L,
      context_n_unique_children = 0L,
      context_weighted_sd = NA_real_,
      context_between_patient_sd = NA_real_,
      context_i2 = NA_real_,
      context_sign_consistency = NA_real_,
      context_p_deleterious = 0,
      context_p_beneficial = 0,
      context_p_neutral = 0,
      context_effect_class = class,
      context_heterogeneity_class = "context_sparse_unknown",
      context_lambda = borrow$lambda,
      context_sd_multiplier = borrow$sd_multiplier,
      context_support_score = 0,
      context_sparse_unknown_flag = TRUE,
      context_high_variable_flag = FALSE,
      context_conflicting_zero_flag = FALSE,
      warning_flags = "",
      stringsAsFactors = FALSE
    ))
  }
  w <- as.numeric(patient_evidence$patient_weight)
  delta <- as.numeric(patient_evidence$delta_patient_mean)
  se <- pmax(as.numeric(patient_evidence$delta_patient_se), cohort_context_sd_floor)
  ok <- is.finite(w) & w > 0 & is.finite(delta)
  w <- w[ok]
  delta <- delta[ok]
  se <- se[ok]
  evidence_ok <- patient_evidence[ok, , drop = FALSE]
  if (!length(delta) || sum(w) <= 0) {
    return(classify_contextual_transition_prior(data.frame(), cohort_context_sd_floor = cohort_context_sd_floor, ...))
  }
  mu <- stats::weighted.mean(delta, w)
  weighted_sd <- sqrt(stats::weighted.mean((delta - mu)^2, w))
  se_mu <- sqrt(1 / sum(w / (se^2 + cohort_context_sd_floor^2)))
  between_patient_sd <- if (length(delta) > 1L) stats::sd(delta) else 0
  i2 <- if (is.finite(weighted_sd) && weighted_sd > 0) {
    max(0, (weighted_sd^2 - mean(se^2)) / weighted_sd^2)
  } else {
    0
  }
  sign_mu <- sign(mu)
  sign_consistency <- if (abs(mu) <= cohort_context_effect_threshold) {
    sum(w[abs(delta) <= cohort_context_effect_threshold]) / sum(w)
  } else {
    sum(w[sign(delta) == sign_mu]) / sum(w)
  }
  eff_n <- sum(w)^2 / sum(w^2)
  patient_weights <- w
  eff_patients <- sum(patient_weights)^2 / sum(patient_weights^2)
  n_patients <- length(unique(evidence_ok$patient_id))
  n_unique_children <- length(unique(evidence_ok$child_karyotype))
  p_beneficial <- 1 - stats::pnorm(cohort_context_effect_threshold, mean = mu, sd = max(se_mu, cohort_context_sd_floor))
  p_deleterious <- stats::pnorm(-cohort_context_effect_threshold, mean = mu, sd = max(se_mu, cohort_context_sd_floor))
  p_neutral <- stats::pnorm(cohort_context_effect_threshold, mean = mu, sd = max(se_mu, cohort_context_sd_floor)) -
    stats::pnorm(-cohort_context_effect_threshold, mean = mu, sd = max(se_mu, cohort_context_sd_floor))
  sparse <- n_patients < cohort_context_min_patients ||
    eff_n < cohort_context_min_effective_n ||
    eff_patients < cohort_context_min_effective_patients ||
    n_unique_children < cohort_context_min_unique_children
  high_variable <- !sparse && (
    (is.finite(weighted_sd) && weighted_sd >= cohort_context_high_weighted_sd) ||
      (is.finite(between_patient_sd) && between_patient_sd >= cohort_context_high_between_patient_sd) ||
      (is.finite(i2) && i2 >= cohort_context_high_i2) ||
      (is.finite(sign_consistency) && sign_consistency < cohort_context_sign_consistency_threshold)
  )
  zero_conflict <- FALSE
  if (!is.null(zero_neighbors) && is.data.frame(zero_neighbors) && nrow(zero_neighbors)) {
    zero_weight <- sum(zero_neighbors$final_weight, na.rm = TRUE)
    zero_conflict <- is.finite(zero_weight) && zero_weight > 0 &&
      (p_beneficial >= 0.8 || p_neutral >= 0.6)
  }
  if (sparse) {
    class <- "context_sparse_unknown"
  } else if (zero_conflict) {
    class <- "context_conflicting_zero"
  } else if (high_variable) {
    class <- "context_high_variable"
  } else if (p_deleterious >= 0.8 && sign_consistency >= cohort_context_sign_consistency_threshold) {
    class <- "context_consistent_deleterious"
  } else if (p_beneficial >= 0.8 && sign_consistency >= cohort_context_sign_consistency_threshold) {
    class <- "context_consistent_beneficial"
  } else if (p_neutral >= 0.6 && weighted_sd < cohort_context_high_weighted_sd) {
    class <- "context_consistent_neutral"
  } else {
    class <- "context_sparse_unknown"
  }
  heterogeneity_class <- if (class == "context_high_variable") "context_high_variable" else if (class == "context_sparse_unknown") "context_sparse_unknown" else "context_supported"
  borrow <- cohort_context_class_borrowing(class, ...)
  prior_sd <- max(cohort_context_sd_floor, weighted_sd, se_mu) * borrow$sd_multiplier
  warning_flags <- if (identical(class, "context_consistent_beneficial")) "survivor_bias_warning" else ""
  data.frame(
    context_delta_mu = mu,
    context_delta_sd = prior_sd,
    context_effective_n = eff_n,
    context_n_patients = n_patients,
    context_n_unique_children = n_unique_children,
    context_weighted_sd = weighted_sd,
    context_between_patient_sd = between_patient_sd,
    context_i2 = i2,
    context_sign_consistency = sign_consistency,
    context_p_deleterious = p_deleterious,
    context_p_beneficial = p_beneficial,
    context_p_neutral = p_neutral,
    context_effect_class = class,
    context_heterogeneity_class = heterogeneity_class,
    context_lambda = borrow$lambda,
    context_sd_multiplier = borrow$sd_multiplier,
    context_support_score = min(1, eff_n / max(1, cohort_context_min_effective_n)),
    context_sparse_unknown_flag = identical(class, "context_sparse_unknown"),
    context_high_variable_flag = identical(class, "context_high_variable"),
    context_conflicting_zero_flag = identical(class, "context_conflicting_zero"),
    warning_flags = warning_flags,
    stringsAsFactors = FALSE
  )
}

cohort_context_patient_level_neighbors <- function(evidence_bank, weights_df, sd_floor = 0.05) {
  if (!nrow(weights_df)) return(data.frame())
  idx <- match(weights_df$evidence_row_id, evidence_bank$evidence_id)
  ev <- evidence_bank[idx[!is.na(idx)], , drop = FALSE]
  weights_df <- weights_df[!is.na(idx), , drop = FALSE]
  if (!nrow(ev)) return(data.frame())
  split_key <- as.character(ev$patient_id)
  rows <- lapply(split(seq_len(nrow(ev)), split_key), function(ii) {
    ww <- weights_df$final_weight[ii]
    se <- pmax(as.numeric(ev$delta_se[ii]), sd_floor)
    inv <- ww / (se^2 + sd_floor^2)
    inv[!is.finite(inv) | inv < 0] <- 0
    if (sum(inv) <= 0) inv <- ww
    delta <- as.numeric(ev$delta_hat[ii])
    ok <- is.finite(delta) & is.finite(inv) & inv > 0
    data.frame(
      patient_id = ev$patient_id[ii[1]],
      delta_patient_mean = if (any(ok)) stats::weighted.mean(delta[ok], inv[ok]) else NA_real_,
      delta_patient_se = if (sum(inv[ok]) > 0) sqrt(1 / sum(inv[ok])) else NA_real_,
      patient_weight = sum(ww, na.rm = TRUE),
      n_context_neighbors = length(ii),
      child_karyotype = paste(unique(ev$child_karyotype[ii]), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Lookup a target-specific contextual transition prior
#'
#' @param target_parent_karyotype,target_child_karyotype Target transition.
#' @param target_patient_id Target patient for leave-one-patient-out lookup.
#' @param evidence_bank Observed contextual transition evidence.
#' @param zero_evidence_bank Optional zero-censoring evidence bank.
#' @param leave_one_patient_out Exclude `target_patient_id` from evidence.
#' @param baseline_ploidy,chromosome_weights,profile_transform Karyotype
#'   feature configuration.
#' @param profile_distance,event_match Context distance configuration.
#' @param bandwidths,weights Context kernel bandwidths and component weights.
#' @param k_nearest,min_kernel_weight Neighbor retention controls.
#' @param cohort_context_sd_floor Minimum contextual prior SD.
#' @param ... Context classification and borrowing controls.
#' @return A list with `prior`, nearest neighbors, and diagnostics.
#' @export
lookup_contextual_transition_prior <- function(target_parent_karyotype,
                                               target_child_karyotype,
                                               target_patient_id,
                                               evidence_bank,
                                               zero_evidence_bank = NULL,
                                               leave_one_patient_out = TRUE,
                                               baseline_ploidy = 2,
                                               chromosome_weights = NULL,
                                               profile_transform = "mass",
                                               profile_distance = "hellinger",
                                               event_match = "same_chr_direction",
                                               bandwidths = list(profile = 0.25, area = 2, burden = 2, local = 1, event = 1),
                                               weights = list(profile = 1.0, area = 0.5, burden = 0.5, local = 1.0, event = 2.0),
                                               k_nearest = 50,
                                               min_kernel_weight = 1e-6,
                                               cohort_context_sd_floor = 0.05,
                                               ...) {
  target_context <- compute_transition_context_features(
    target_parent_karyotype,
    target_child_karyotype,
    baseline_ploidy = baseline_ploidy,
    chromosome_weights = chromosome_weights,
    profile_transform = profile_transform
  )
  ev <- evidence_bank
  if (isTRUE(leave_one_patient_out) && !is.null(target_patient_id) && "patient_id" %in% names(ev)) {
    ev <- cohort_context_subset_bank(ev, as.character(ev$patient_id) != as.character(target_patient_id))
  }
  neighbors <- compute_context_kernel_weights(
    target_context = target_context,
    evidence_contexts = ev,
    bandwidths = bandwidths,
    weights = weights,
    event_match = event_match,
    k_nearest = k_nearest,
    min_kernel_weight = min_kernel_weight,
    profile_distance = profile_distance,
    chromosome_weights = chromosome_weights
  )
  zero_neighbors <- data.frame()
  if (!is.null(zero_evidence_bank) && is.data.frame(zero_evidence_bank) && nrow(zero_evidence_bank)) {
    zev <- zero_evidence_bank
    if (isTRUE(leave_one_patient_out) && !is.null(target_patient_id) && "patient_id" %in% names(zev)) {
      zev <- cohort_context_subset_bank(zev, as.character(zev$patient_id) != as.character(target_patient_id))
    }
    zero_neighbors <- compute_context_kernel_weights(
      target_context = target_context,
      evidence_contexts = zev,
      bandwidths = bandwidths,
      weights = weights,
      event_match = event_match,
      k_nearest = k_nearest,
      min_kernel_weight = min_kernel_weight,
      profile_distance = profile_distance,
      chromosome_weights = chromosome_weights
    )
  }
  patient_evidence <- cohort_context_patient_level_neighbors(ev, neighbors, sd_floor = cohort_context_sd_floor)
  prior_row <- classify_contextual_transition_prior(
    patient_evidence,
    zero_neighbors = zero_neighbors,
    cohort_context_sd_floor = cohort_context_sd_floor,
    ...
  )
  prior_row$target_parent_karyotype <- cohort_context_karyotype_label(target_parent_karyotype)
  prior_row$target_child_karyotype <- cohort_context_karyotype_label(target_child_karyotype)
  prior_row$target_patient_id <- if (is.null(target_patient_id)) NA_character_ else as.character(target_patient_id)
  prior_row$transition_chr <- target_context$transition_chr[1]
  prior_row$transition_direction <- target_context$transition_direction[1]
  prior_row$transition_size <- target_context$transition_size[1]
  prior_row <- prior_row[c(
    "target_parent_karyotype", "target_child_karyotype", "target_patient_id",
    setdiff(names(prior_row), c("target_parent_karyotype", "target_child_karyotype", "target_patient_id"))
  )]
  list(
    prior = prior_row,
    neighbors = neighbors,
    zero_neighbors = zero_neighbors,
    patient_evidence = patient_evidence,
    diagnostics = list(
      n_candidate_evidence = nrow(ev),
      n_neighbors = nrow(neighbors),
      n_zero_neighbors = nrow(zero_neighbors),
      leave_one_patient_out = isTRUE(leave_one_patient_out)
    )
  )
}

combine_contextual_path_priors <- function(path_priors, path_weights) {
  priors <- lapply(path_priors, `[[`, "prior")
  prior_df <- do.call(rbind, priors)
  mu <- prior_df$context_delta_mu
  sd <- prior_df$context_delta_sd
  path_weights <- normalize_nn_weights(path_weights, fallback_n = length(mu))
  mu_combined <- sum(path_weights * mu, na.rm = TRUE)
  var_combined <- sum(path_weights * (sd^2 + (mu - mu_combined)^2), na.rm = TRUE)
  dominant <- path_weights >= 0.2
  if (!any(dominant)) dominant[which.max(path_weights)] <- TRUE
  classes <- prior_df$context_effect_class
  conservative <- if (any(classes[dominant] == "context_high_variable")) {
    "context_high_variable"
  } else if (any(classes[dominant] == "context_sparse_unknown")) {
    "context_sparse_unknown"
  } else if (any(classes[dominant] == "context_conflicting_zero")) {
    "context_conflicting_zero"
  } else if (length(unique(classes[dominant])) > 1L) {
    "mixed_context"
  } else {
    classes[which.max(path_weights)]
  }
  data.frame(
    context_delta_mu = mu_combined,
    context_delta_sd = sqrt(max(var_combined, 0)),
    context_effective_n = sum(path_weights * prior_df$context_effective_n, na.rm = TRUE),
    context_n_patients = max(prior_df$context_n_patients, na.rm = TRUE),
    context_n_unique_children = max(prior_df$context_n_unique_children, na.rm = TRUE),
    context_weighted_sd = sum(path_weights * prior_df$context_weighted_sd, na.rm = TRUE),
    context_between_patient_sd = sum(path_weights * prior_df$context_between_patient_sd, na.rm = TRUE),
    context_sign_consistency = sum(path_weights * prior_df$context_sign_consistency, na.rm = TRUE),
    context_p_deleterious = sum(path_weights * prior_df$context_p_deleterious, na.rm = TRUE),
    context_p_beneficial = sum(path_weights * prior_df$context_p_beneficial, na.rm = TRUE),
    context_p_neutral = sum(path_weights * prior_df$context_p_neutral, na.rm = TRUE),
    context_effect_class = conservative,
    context_heterogeneity_class = if (conservative %in% c("context_high_variable", "mixed_context")) "context_high_variable" else prior_df$context_heterogeneity_class[which.max(path_weights)],
    context_lambda = {
      val <- suppressWarnings(min(prior_df$context_lambda[dominant], na.rm = TRUE))
      if (is.finite(val)) val else 0
    },
    context_sd_multiplier = {
      val <- suppressWarnings(max(prior_df$context_sd_multiplier[dominant], na.rm = TRUE))
      if (is.finite(val)) val else 4
    },
    context_support_score = sum(path_weights * prior_df$context_support_score, na.rm = TRUE),
    context_sparse_unknown_flag = any(prior_df$context_sparse_unknown_flag[dominant]),
    context_high_variable_flag = any(prior_df$context_high_variable_flag[dominant]) || identical(conservative, "mixed_context"),
    context_conflicting_zero_flag = any(prior_df$context_conflicting_zero_flag[dominant]),
    mixed_context_class_flag = length(unique(classes[dominant])) > 1L,
    stringsAsFactors = FALSE
  )
}

context_label_for_overlay <- function(child_is_zero,
                                      non_identifiable_zero,
                                      context_class,
                                      update_applied,
                                      skipped_reason,
                                      guardrail_hit) {
  if (isTRUE(guardrail_hit) && !isTRUE(update_applied)) return("guardrail_skipped")
  if (!isTRUE(child_is_zero)) return("patient_observed_no_context_update")
  if (isTRUE(non_identifiable_zero)) return("low_exposure_zero_nonidentifiable")
  switch(
    context_class,
    context_consistent_deleterious = "informative_zero_context_consistent_deleterious",
    context_consistent_neutral = "informative_zero_context_consistent_neutral",
    context_consistent_beneficial = "informative_zero_context_consistent_beneficial_conservative",
    context_high_variable = "context_high_variable_uncertain",
    context_sparse_unknown = "context_sparse_unknown_nonidentifiable",
    context_conflicting_zero = "context_high_variable_uncertain",
    mixed_context = "mixed_parent_context_uncertain",
    if (!is.na(skipped_reason) && identical(skipped_reason, "fallback_to_v2_broad_prior")) "fallback_to_v2_broad_prior" else "context_sparse_unknown_nonidentifiable"
  )
}

#' Apply contextual cohort-transition overlay to one NN child
#'
#' The contextual overlay combines the patient-specific two-shell baseline with
#' a target-specific Delta-fitness prior learned from similar parent karyotype
#' backgrounds and CNA events. Observed NN are unchanged by default, sparse or
#' high-variable contexts keep the baseline, and zero NN updates are capped by
#' borrowing and shift guardrails.
#'
#' @param item NN child context object.
#' @param child_name Child karyotype ID.
#' @param build_opt_fc Objective builder used by `solve_fitness_bootstrap()`.
#' @param search_interval Numeric optimization interval.
#' @param prior_use Contextual prior object for the target patient.
#' @param f_two_shell_baseline Patient-specific two-shell baseline fitness.
#' @param nn_present Optional logical indicating whether the child is observed.
#' @param two_shell_node_diagnostics Optional one-row two-shell diagnostics.
#' @param cohort_contextual_apply_to Which nodes may receive contextual updates.
#' @param cohort_context_lambda Global contextual borrowing multiplier.
#' @param cohort_context_max_borrowing_fraction Maximum borrowing fraction.
#' @param cohort_context_max_abs_delta_shift Optional shift cap.
#' @param cohort_context_sd_floor,cohort_context_patient_sd_floor Contextual SD
#'   floors.
#' @param cohort_context_keep_baseline_when_sparse,cohort_context_keep_baseline_when_high_variable
#'   Keep the two-shell baseline for sparse/high-variable contexts.
#'
#' @return A list with final fitness and node diagnostics.
#' @export
apply_contextual_cohort_overlay <- function(item,
                                            child_name,
                                            build_opt_fc,
                                            search_interval,
                                            prior_use,
                                            f_two_shell_baseline,
                                            nn_present = NULL,
                                            two_shell_node_diagnostics = NULL,
                                            cohort_contextual_apply_to = c("zero_only", "low_information", "all"),
                                            cohort_context_lambda = 0.25,
                                            cohort_context_max_borrowing_fraction = 0.5,
                                            cohort_context_max_abs_delta_shift = NULL,
                                            cohort_context_sd_floor = 0.05,
                                            cohort_context_patient_sd_floor = 0.10,
                                            cohort_context_keep_baseline_when_sparse = TRUE,
                                            cohort_context_keep_baseline_when_high_variable = TRUE) {
  cohort_contextual_apply_to <- match.arg(cohort_contextual_apply_to)
  direct_objective <- build_opt_fc(item, do_prior_param = FALSE)
  n_parents <- length(item$parent_fitness)
  if (n_parents == 0L) {
    return(list(f_final = f_two_shell_baseline, diagnostics = data.frame()))
  }
  parent_karyotypes <- item$nj
  if (is.null(parent_karyotypes) || length(parent_karyotypes) != n_parents ||
      any(is.na(parent_karyotypes)) || any(!nzchar(parent_karyotypes))) {
    parent_karyotypes <- names(item$parent_fitness)
  }
  path_weights <- normalize_nn_weights(item$parent_opportunity_weights, fallback_n = n_parents)
  parent_fit <- as.numeric(item$parent_fitness)
  selector <- should_apply_cohort_transition_to_node(
    item = item,
    child_name = child_name,
    nn_present = nn_present,
    two_shell_node_diagnostics = two_shell_node_diagnostics,
    cohort_transition_apply_to = cohort_contextual_apply_to
  )
  expected_parent_like <- as.numeric(item$projected_exposure)
  if (length(expected_parent_like) != 1L || !is.finite(expected_parent_like)) expected_parent_like <- NA_real_
  zero_info <- compute_zero_informativeness_score(expected_parent_like, informative_threshold = prior_use$context_feature_config$zero_min_expected_count %||% 3)
  child_observed_count <- sum(item$child_obs, na.rm = TRUE)
  child_is_zero <- selector$child_is_zero
  non_identifiable_zero <- isTRUE(child_is_zero) &&
    (!is.finite(expected_parent_like) || expected_parent_like < (prior_use$context_feature_config$zero_min_expected_count %||% 3))
  lookups <- lapply(seq_len(n_parents), function(idx) {
    lookup_contextual_transition_prior(
      target_parent_karyotype = parent_karyotypes[idx],
      target_child_karyotype = child_name,
      target_patient_id = prior_use$target_patient_id,
      evidence_bank = prior_use$evidence_bank,
      zero_evidence_bank = prior_use$zero_evidence_bank,
      leave_one_patient_out = prior_use$leave_one_patient_out,
      baseline_ploidy = prior_use$context_feature_config$baseline_ploidy,
      chromosome_weights = prior_use$context_feature_config$chromosome_weights,
      profile_transform = prior_use$context_feature_config$profile_transform,
      profile_distance = prior_use$context_feature_config$profile_distance,
      event_match = prior_use$context_feature_config$event_match,
      bandwidths = prior_use$context_bandwidths,
      weights = prior_use$context_weight_config,
      k_nearest = prior_use$context_feature_config$k_nearest,
      min_kernel_weight = prior_use$context_feature_config$min_kernel_weight,
      cohort_context_sd_floor = cohort_context_sd_floor,
      cohort_context_min_patients = prior_use$context_feature_config$min_patients,
      cohort_context_min_effective_n = prior_use$context_feature_config$min_effective_n,
      cohort_context_min_effective_patients = prior_use$context_feature_config$min_effective_patients,
      cohort_context_min_unique_children = prior_use$context_feature_config$min_unique_children,
      cohort_context_effect_threshold = prior_use$context_feature_config$effect_threshold,
      cohort_context_sign_consistency_threshold = prior_use$context_feature_config$sign_consistency_threshold,
      cohort_context_high_weighted_sd = prior_use$context_feature_config$high_weighted_sd,
      cohort_context_high_between_patient_sd = prior_use$context_feature_config$high_between_patient_sd,
      cohort_context_high_i2 = prior_use$context_feature_config$high_i2,
      cohort_context_lambda_consistent_deleterious = prior_use$context_feature_config$lambda_consistent_deleterious,
      cohort_context_lambda_consistent_neutral = prior_use$context_feature_config$lambda_consistent_neutral,
      cohort_context_lambda_consistent_beneficial = prior_use$context_feature_config$lambda_consistent_beneficial,
      cohort_context_lambda_high_variable = prior_use$context_feature_config$lambda_high_variable,
      cohort_context_lambda_sparse_unknown = prior_use$context_feature_config$lambda_sparse_unknown,
      cohort_context_lambda_conflicting_zero = prior_use$context_feature_config$lambda_conflicting_zero,
      cohort_context_sd_multiplier_consistent_deleterious = prior_use$context_feature_config$sd_multiplier_consistent_deleterious,
      cohort_context_sd_multiplier_consistent_neutral = prior_use$context_feature_config$sd_multiplier_consistent_neutral,
      cohort_context_sd_multiplier_consistent_beneficial = prior_use$context_feature_config$sd_multiplier_consistent_beneficial,
      cohort_context_sd_multiplier_high_variable = prior_use$context_feature_config$sd_multiplier_high_variable,
      cohort_context_sd_multiplier_sparse_unknown = prior_use$context_feature_config$sd_multiplier_sparse_unknown,
      cohort_context_sd_multiplier_conflicting_zero = prior_use$context_feature_config$sd_multiplier_conflicting_zero
    )
  })
  combined <- combine_contextual_path_priors(lookups, path_weights)
  parent_combined <- sum(path_weights * parent_fit, na.rm = TRUE)
  direct_se <- estimate_scalar_objective_se(
    objective_fn = direct_objective,
    optimum = if (is.finite(f_two_shell_baseline)) f_two_shell_baseline else mean(search_interval),
    search_interval = search_interval,
    se_floor = max(cohort_context_sd_floor, cohort_context_patient_sd_floor)
  )
  anchor_sd <- max(cohort_context_sd_floor, cohort_context_patient_sd_floor, direct_se, na.rm = TRUE)
  if (!is.finite(anchor_sd) || anchor_sd <= 0) anchor_sd <- max(cohort_context_sd_floor, cohort_context_patient_sd_floor)
  anchor_info <- if (is.finite(f_two_shell_baseline)) 1 / anchor_sd^2 else 0
  context_sd <- max(cohort_context_sd_floor, cohort_context_patient_sd_floor, combined$context_delta_sd[1])
  context_target <- parent_combined + combined$context_delta_mu[1]
  zero_multiplier <- if (isTRUE(child_is_zero)) pmin(1, pmax(0, zero_info$zero_informativeness_score[1])) else 1
  if (!is.finite(zero_multiplier) || isTRUE(non_identifiable_zero)) zero_multiplier <- 0
  class_disallowed <- combined$context_effect_class %in% c("context_high_variable", "context_sparse_unknown", "context_conflicting_zero", "mixed_context")
  if (isTRUE(cohort_context_keep_baseline_when_high_variable) && isTRUE(combined$context_high_variable_flag[1])) class_disallowed <- TRUE
  if (isTRUE(cohort_context_keep_baseline_when_sparse) && isTRUE(combined$context_sparse_unknown_flag[1])) class_disallowed <- TRUE
  effective_lambda <- cohort_context_lambda * combined$context_lambda[1] * zero_multiplier * combined$context_support_score[1]
  if (isTRUE(class_disallowed) || !isTRUE(selector$apply)) effective_lambda <- 0
  if (!is.finite(effective_lambda) || effective_lambda < 0) effective_lambda <- 0
  prior_info <- effective_lambda / context_sd^2
  f_overlay <- f_final <- f_two_shell_baseline
  update_applied <- FALSE
  guardrail_hit <- FALSE
  skipped_reason <- selector$reason
  if (isTRUE(selector$apply) && prior_info > 0 && anchor_info > 0 && is.finite(f_two_shell_baseline)) {
    f_overlay <- (anchor_info * f_two_shell_baseline + prior_info * context_target) / (anchor_info + prior_info)
    borrowing <- prior_info / (prior_info + anchor_info)
    max_shift <- cohort_context_max_abs_delta_shift
    if (is.null(max_shift)) {
      max_shift <- max(2 * context_sd, 2 * anchor_sd, 0.10)
    }
    shift <- f_overlay - f_two_shell_baseline
    if (is.finite(shift) && abs(shift) > max_shift) {
      f_overlay <- f_two_shell_baseline + sign(shift) * max_shift
      guardrail_hit <- TRUE
    }
    if (is.finite(borrowing) && borrowing > cohort_context_max_borrowing_fraction) {
      f_final <- f_two_shell_baseline
      guardrail_hit <- TRUE
      skipped_reason <- "borrowing_fraction_guardrail"
    } else {
      f_final <- f_overlay
      update_applied <- is.finite(f_final) && abs(f_final - f_two_shell_baseline) > sqrt(.Machine$double.eps)
      if (!isTRUE(update_applied)) skipped_reason <- "overlay_shift_negligible"
    }
  } else if (isTRUE(selector$apply) && isTRUE(non_identifiable_zero)) {
    skipped_reason <- "non_identifiable_low_exposure_zero"
  } else if (isTRUE(selector$apply) && prior_info <= 0) {
    skipped_reason <- if (isTRUE(class_disallowed)) "context_class_disallows_borrowing" else "context_support_insufficient"
  }
  borrowing <- if (prior_info + anchor_info > 0) prior_info / (prior_info + anchor_info) else NA_real_
  q_delta <- c(0.8, 0.9, 0.95)
  delta_upper <- combined$context_delta_mu[1] + stats::qnorm(q_delta) * context_sd
  fitness_upper <- parent_combined + delta_upper
  neighbor_ids <- paste(utils::head(unlist(lapply(lookups, function(x) x$neighbors$evidence_row_id)), 10), collapse = ";")
  neighbor_patients <- paste(utils::head(unlist(lapply(lookups, function(x) x$neighbors$patient_id)), 10), collapse = ";")
  neighbor_weights <- paste(signif(utils::head(unlist(lapply(lookups, function(x) x$neighbors$final_weight)), 10), 4), collapse = ";")
  neighbor_deltas <- paste(signif(utils::head(unlist(lapply(lookups, function(x) {
    idx <- match(x$neighbors$evidence_row_id, prior_use$evidence_bank$evidence_id)
    prior_use$evidence_bank$delta_hat[idx]
  })), 10), 4), collapse = ";")
  label <- context_label_for_overlay(child_is_zero, non_identifiable_zero, combined$context_effect_class[1], update_applied, skipped_reason, guardrail_hit)
  rows <- lapply(seq_len(n_parents), function(idx) {
    pr <- lookups[[idx]]$prior
    data.frame(
      karyotype = child_name,
      patient_id = prior_use$target_patient_id %||% NA_character_,
      parent_karyotype = parent_karyotypes[idx],
      child_karyotype = child_name,
      child_is_zero = child_is_zero,
      child_observed_count = child_observed_count,
      expected_count_parent_like = zero_info$expected_count_parent_like,
      zero_informativeness_score = zero_info$zero_informativeness_score,
      transition_chr = pr$transition_chr %||% NA_integer_,
      transition_direction = pr$transition_direction %||% NA_character_,
      transition_size = pr$transition_size %||% NA_real_,
      context_effective_n = combined$context_effective_n,
      context_n_patients = combined$context_n_patients,
      context_n_unique_children = combined$context_n_unique_children,
      context_support_score = combined$context_support_score,
      context_weighted_mean_delta = combined$context_delta_mu,
      context_delta_mu = combined$context_delta_mu,
      context_delta_sd = context_sd,
      context_weighted_sd = combined$context_weighted_sd,
      context_between_patient_sd = combined$context_between_patient_sd,
      context_sign_consistency = combined$context_sign_consistency,
      context_p_deleterious = combined$context_p_deleterious,
      context_p_beneficial = combined$context_p_beneficial,
      context_p_neutral = combined$context_p_neutral,
      context_effect_class = combined$context_effect_class,
      context_heterogeneity_class = combined$context_heterogeneity_class,
      context_sparse_unknown_flag = combined$context_sparse_unknown_flag,
      context_high_variable_flag = combined$context_high_variable_flag,
      context_conflicting_zero_flag = combined$context_conflicting_zero_flag,
      mixed_context_class_flag = combined$mixed_context_class_flag,
      context_lambda = combined$context_lambda,
      effective_context_lambda = effective_lambda,
      context_sd_multiplier = combined$context_sd_multiplier,
      context_prior_dominated_flag = isTRUE(child_is_zero) && is.finite(borrowing) && borrowing > cohort_context_max_borrowing_fraction,
      parent_fitness = parent_fit[idx],
      path_responsibility = path_weights[idx],
      f_two_shell_baseline = f_two_shell_baseline,
      f_contextual_overlay = f_overlay,
      f_cohort_overlay = f_overlay,
      f_final = f_final,
      f_delta_from_two_shell = f_final - f_two_shell_baseline,
      delta_context_mean = combined$context_delta_mu,
      delta_context_sd = context_sd,
      delta_posterior_mean = f_final - parent_combined,
      delta_upper_80 = delta_upper[1],
      delta_upper_90 = delta_upper[2],
      delta_upper_95 = delta_upper[3],
      fitness_posterior_mean = f_final,
      fitness_upper_80 = fitness_upper[1],
      fitness_upper_90 = fitness_upper[2],
      fitness_upper_95 = fitness_upper[3],
      posterior_interval_approximation = "normal_prior_overlay",
      cohort_update_applied = update_applied,
      cohort_update_skipped_reason = if (is.na(skipped_reason)) NA_character_ else skipped_reason,
      guardrail_hit = guardrail_hit,
      borrowing_fraction = borrowing,
      cohort_borrowing_fraction = borrowing,
      non_identifiable_zero_flag = non_identifiable_zero,
      context_label = label,
      nearest_context_evidence_ids = neighbor_ids,
      nearest_context_evidence_patients = neighbor_patients,
      nearest_context_evidence_weights = neighbor_weights,
      nearest_context_evidence_deltas = neighbor_deltas,
      stringsAsFactors = FALSE
    )
  })
  list(f_final = f_final, diagnostics = do.call(rbind, rows))
}

learn_cohort_transition_prior_contextual <- function(records,
                                                     leave_one_patient_out = TRUE,
                                                     grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                                     cohort_context_baseline_ploidy = 2,
                                                     cohort_context_chromosome_weights = NULL,
                                                     cohort_context_profile_transform = c("mass", "centered", "zscore", "raw"),
                                                     cohort_context_profile_distance = c("hellinger", "jensen_shannon", "cosine", "euclidean", "manhattan"),
                                                     cohort_context_event_match = c("same_chr_direction", "same_direction", "kernel"),
                                                     cohort_context_bandwidth_profile = NULL,
                                                     cohort_context_bandwidth_area = NULL,
                                                     cohort_context_bandwidth_burden = NULL,
                                                     cohort_context_bandwidth_local = NULL,
                                                     cohort_context_bandwidth_event = NULL,
                                                     cohort_context_profile_weight = 1.0,
                                                     cohort_context_area_weight = 0.5,
                                                     cohort_context_burden_weight = 0.5,
                                                     cohort_context_local_weight = 1.0,
                                                     cohort_context_event_weight = 2.0,
                                                     cohort_context_min_patients = 3L,
                                                     cohort_context_min_effective_n = 5,
                                                     cohort_context_min_effective_patients = 3,
                                                     cohort_context_min_unique_children = 3L,
                                                     cohort_context_k_nearest = 50,
                                                     cohort_context_min_kernel_weight = 1e-6,
                                                     cohort_context_effect_threshold = 0.02,
                                                     cohort_context_sign_consistency_threshold = 0.75,
                                                     cohort_context_high_weighted_sd = 0.10,
                                                     cohort_context_high_between_patient_sd = 0.10,
                                                     cohort_context_high_i2 = 0.50,
                                                     cohort_context_lambda_consistent_deleterious = 0.50,
                                                     cohort_context_lambda_consistent_neutral = 0.25,
                                                     cohort_context_lambda_consistent_beneficial = 0.10,
                                                     cohort_context_lambda_high_variable = 0.00,
                                                     cohort_context_lambda_sparse_unknown = 0.00,
                                                     cohort_context_lambda_conflicting_zero = 0.00,
                                                     cohort_context_sd_floor = 0.05,
                                                     cohort_context_patient_sd_floor = 0.10,
                                                     cohort_context_sd_multiplier_consistent_deleterious = 1.0,
                                                     cohort_context_sd_multiplier_consistent_neutral = 1.5,
                                                     cohort_context_sd_multiplier_consistent_beneficial = 2.5,
                                                     cohort_context_sd_multiplier_high_variable = 4.0,
                                                     cohort_context_sd_multiplier_sparse_unknown = 4.0,
                                                     cohort_context_sd_multiplier_conflicting_zero = 4.0,
                                                     cohort_context_zero_as_censoring_only = TRUE,
                                                     cohort_context_zero_min_expected_count = 3.0,
                                                     cohort_context_zero_weight_cap_ratio = 0.25,
                                                     cohort_transition_use_prior_dominated_records = FALSE,
                                                     cohort_transition_use_boundary_records = FALSE,
                                                     cohort_transition_max_delta_se = NULL,
                                                     cohort_transition_max_delta_se_quantile = 0.75,
                                                     cohort_transition_min_path_responsibility = 0.05,
                                                     cohort_transition_min_observed_count = 1L,
                                                     ...) {
  grouping <- match.arg(grouping)
  cohort_context_profile_transform <- match.arg(cohort_context_profile_transform)
  cohort_context_profile_distance <- match.arg(cohort_context_profile_distance)
  cohort_context_event_match <- match.arg(cohort_context_event_match)
  validate_positive_finite(cohort_context_sd_floor, "cohort_context_sd_floor")
  validate_positive_finite(cohort_context_patient_sd_floor, "cohort_context_patient_sd_floor")
  bank <- build_contextual_transition_evidence_bank(
    records = records,
    baseline_ploidy = cohort_context_baseline_ploidy,
    chromosome_weights = cohort_context_chromosome_weights,
    profile_transform = cohort_context_profile_transform,
    cohort_transition_use_prior_dominated_records = cohort_transition_use_prior_dominated_records,
    cohort_transition_use_boundary_records = cohort_transition_use_boundary_records,
    cohort_transition_max_delta_se = cohort_transition_max_delta_se,
    cohort_transition_max_delta_se_quantile = cohort_transition_max_delta_se_quantile,
    cohort_transition_min_path_responsibility = cohort_transition_min_path_responsibility,
    cohort_transition_min_observed_count = cohort_transition_min_observed_count,
    cohort_context_zero_min_expected_count = cohort_context_zero_min_expected_count,
    cohort_context_zero_as_censoring_only = cohort_context_zero_as_censoring_only
  )
  bandwidths <- estimate_context_bandwidths(
    bank$evidence_bank,
    profile_distance = cohort_context_profile_distance,
    chromosome_weights = cohort_context_chromosome_weights,
    cohort_context_bandwidth_profile = cohort_context_bandwidth_profile,
    cohort_context_bandwidth_area = cohort_context_bandwidth_area,
    cohort_context_bandwidth_burden = cohort_context_bandwidth_burden,
    cohort_context_bandwidth_local = cohort_context_bandwidth_local,
    cohort_context_bandwidth_event = cohort_context_bandwidth_event
  )
  v2_fallback <- tryCatch(
    learn_cohort_transition_prior_v2(
      records = records,
      leave_one_patient_out = leave_one_patient_out,
      grouping = grouping,
      cohort_transition_sd_floor = max(cohort_context_sd_floor, 0.05),
      cohort_transition_patient_sd_floor = max(cohort_context_patient_sd_floor, 0.10),
      cohort_transition_zero_weight_cap_ratio = cohort_context_zero_weight_cap_ratio,
      cohort_transition_use_prior_dominated_records = cohort_transition_use_prior_dominated_records,
      cohort_transition_use_boundary_records = cohort_transition_use_boundary_records,
      cohort_transition_max_delta_se = cohort_transition_max_delta_se,
      cohort_transition_max_delta_se_quantile = cohort_transition_max_delta_se_quantile,
      cohort_transition_min_path_responsibility = cohort_transition_min_path_responsibility,
      cohort_transition_min_observed_count = cohort_transition_min_observed_count,
      cohort_transition_zero_as_censoring_only = cohort_context_zero_as_censoring_only,
      cohort_transition_zero_min_expected_count = cohort_context_zero_min_expected_count,
      ...
    ),
    error = function(e) list(version = "cohort_transition_v2_unavailable", error = conditionMessage(e))
  )
  patient_ids <- sort(unique(as.character(records$patient_id)))
  config <- list(
    baseline_ploidy = cohort_context_baseline_ploidy,
    chromosome_weights = cohort_context_chromosome_weights,
    profile_transform = cohort_context_profile_transform,
    profile_distance = cohort_context_profile_distance,
    event_match = cohort_context_event_match,
    min_patients = cohort_context_min_patients,
    min_effective_n = cohort_context_min_effective_n,
    min_effective_patients = cohort_context_min_effective_patients,
    min_unique_children = cohort_context_min_unique_children,
    k_nearest = cohort_context_k_nearest,
    min_kernel_weight = cohort_context_min_kernel_weight,
    effect_threshold = cohort_context_effect_threshold,
    sign_consistency_threshold = cohort_context_sign_consistency_threshold,
    high_weighted_sd = cohort_context_high_weighted_sd,
    high_between_patient_sd = cohort_context_high_between_patient_sd,
    high_i2 = cohort_context_high_i2,
    lambda_consistent_deleterious = cohort_context_lambda_consistent_deleterious,
    lambda_consistent_neutral = cohort_context_lambda_consistent_neutral,
    lambda_consistent_beneficial = cohort_context_lambda_consistent_beneficial,
    lambda_high_variable = cohort_context_lambda_high_variable,
    lambda_sparse_unknown = cohort_context_lambda_sparse_unknown,
    lambda_conflicting_zero = cohort_context_lambda_conflicting_zero,
    sd_floor = cohort_context_sd_floor,
    patient_sd_floor = cohort_context_patient_sd_floor,
    sd_multiplier_consistent_deleterious = cohort_context_sd_multiplier_consistent_deleterious,
    sd_multiplier_consistent_neutral = cohort_context_sd_multiplier_consistent_neutral,
    sd_multiplier_consistent_beneficial = cohort_context_sd_multiplier_consistent_beneficial,
    sd_multiplier_high_variable = cohort_context_sd_multiplier_high_variable,
    sd_multiplier_sparse_unknown = cohort_context_sd_multiplier_sparse_unknown,
    sd_multiplier_conflicting_zero = cohort_context_sd_multiplier_conflicting_zero,
    zero_as_censoring_only = cohort_context_zero_as_censoring_only,
    zero_min_expected_count = cohort_context_zero_min_expected_count
  )
  weight_config <- list(
    profile = cohort_context_profile_weight,
    area = cohort_context_area_weight,
    burden = cohort_context_burden_weight,
    local = cohort_context_local_weight,
    event = cohort_context_event_weight
  )
  diagnostics <- c(
    list(
      version = "cohort_transition_contextual_v1",
      n_raw_transition_records = nrow(records),
      n_context_observed_evidence = nrow(bank$evidence_bank),
      n_context_zero_evidence = nrow(bank$zero_evidence_bank),
      n_context_records_excluded = bank$diagnostics$n_excluded_records %||% NA_integer_,
      exclusion_reasons = bank$diagnostics[names(cohort_transition_empty_filter_diagnostics())],
      n_patients_in_evidence_bank = length(unique(bank$evidence_bank$patient_id)),
      profile_transform = cohort_context_profile_transform,
      profile_distance = cohort_context_profile_distance,
      event_match_mode = cohort_context_event_match,
      bandwidth_profile = bandwidths$profile,
      bandwidth_area = bandwidths$area,
      bandwidth_burden = bandwidths$burden,
      bandwidth_local = bandwidths$local,
      bandwidth_event = bandwidths$event,
      chromosome_weights_used = !is.null(cohort_context_chromosome_weights)
    ),
    bank$diagnostics
  )
  list(
    version = "cohort_transition_contextual_v1",
    grouping = grouping,
    contextual = TRUE,
    evidence_bank = bank$evidence_bank,
    zero_evidence_bank = bank$zero_evidence_bank,
    context_feature_config = config,
    context_bandwidths = bandwidths,
    context_weight_config = weight_config,
    v2_fallback_prior = v2_fallback,
    patient_ids = patient_ids,
    leave_one_patient_out = isTRUE(leave_one_patient_out),
    loo_priors = list(),
    diagnostics = diagnostics,
    filter_diagnostics = bank$diagnostics,
    excluded_records = bank$filtered$excluded_records
  )
}

cohort_transition_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x) || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1))
}

fit_cohort_transition_nn_child <- function(item,
                                           child_name,
                                           build_opt_fc,
                                           search_interval,
                                           prior_use,
                                           sd_floor = 1e-3) {
  direct_objective <- build_opt_fc(item, do_prior_param = FALSE)
  n_parents <- length(item$parent_fitness)
  if (n_parents == 0L) {
    return(list(f_map = NA_real_, diagnostics = data.frame()))
  }
  parent_karyotypes <- item$nj
  if (is.null(parent_karyotypes) || length(parent_karyotypes) != n_parents ||
      any(is.na(parent_karyotypes)) || any(!nzchar(parent_karyotypes))) {
    parent_karyotypes <- names(item$parent_fitness)
  }
  if (is.null(parent_karyotypes) || length(parent_karyotypes) != n_parents ||
      any(is.na(parent_karyotypes)) || any(!nzchar(parent_karyotypes))) {
    stop("Cohort transition NN contexts must include non-empty parent karyotype IDs.", call. = FALSE)
  }
  path_weights <- normalize_nn_weights(item$parent_opportunity_weights, fallback_n = n_parents)
  path_weights[!is.finite(path_weights) | path_weights < 0] <- 0
  if (sum(path_weights) <= 0) path_weights <- rep(1 / n_parents, n_parents)
  priors <- lapply(seq_len(n_parents), function(idx) {
    lookup_cohort_transition_group_prior(prior_use, parent_karyotypes[idx], child_name)
  })
  prior_mu <- vapply(priors, `[[`, numeric(1), "mu")
  prior_sd <- pmax(vapply(priors, `[[`, numeric(1), "sd"), sd_floor)
  parent_fit <- as.numeric(item$parent_fitness)
  expected_parent_like <- as.numeric(item$projected_exposure)
  if (length(expected_parent_like) != 1L || !is.finite(expected_parent_like)) {
    expected_parent_like <- NA_real_
  }
  zero_info <- compute_zero_informativeness_score(expected_parent_like)
  child_observed_count <- sum(item$child_obs, na.rm = TRUE)
  child_is_zero <- is.finite(child_observed_count) && child_observed_count <= 0
  non_identifiable_zero <- isTRUE(child_is_zero) &&
    (!is.finite(expected_parent_like) || expected_parent_like < 0.5)
  prior_weight_multiplier <- 1
  if (isTRUE(child_is_zero)) {
    prior_weight_multiplier <- zero_info$zero_informativeness_score[1]
    if (!is.finite(prior_weight_multiplier) || expected_parent_like < 0.5) {
      prior_weight_multiplier <- 0
    }
    prior_weight_multiplier <- pmin(1, pmax(0, prior_weight_multiplier))
  }
  effective_path_weights <- path_weights * prior_weight_multiplier
  objective <- function(fc_param) {
    direct <- direct_objective(fc_param)
    if (!is.finite(direct)) return(direct)
    delta <- fc_param - parent_fit
    prior_nll <- 0
    if (sum(effective_path_weights) > 0) {
      prior_nll <- -sum(effective_path_weights * stats::dnorm(delta, mean = prior_mu, sd = prior_sd, log = TRUE))
    }
    direct + prior_nll
  }
  prior_centers <- parent_fit + prior_mu
  base_interval <- range(search_interval, na.rm = TRUE)
  if (length(base_interval) != 2L || any(!is.finite(base_interval))) {
    base_interval <- range(c(parent_fit, prior_centers), na.rm = TRUE)
  }
  if (length(base_interval) != 2L || any(!is.finite(base_interval))) {
    base_interval <- c(-1, 1)
  }
  local_interval <- if (sum(effective_path_weights) > 0) {
    range(c(base_interval, prior_centers[effective_path_weights > 0]), na.rm = TRUE)
  } else {
    base_interval
  }
  span <- diff(local_interval)
  if (!is.finite(span) || span <= 0) {
    pad <- max(sd_floor, 1e-3)
    local_interval <- local_interval[1] + c(-pad, pad)
  }
  opt <- run_optimise_checked(
    objective,
    interval = local_interval,
    context = sprintf("optimise nearest-neighbour fitness with cohort_transition prior for child %s", child_name)
  )
  f_map <- if (is.null(opt)) NA_real_ else opt$minimum
  grid <- seq(local_interval[1], local_interval[2], length.out = 121L)
  log_post <- -vapply(grid, objective, numeric(1))
  finite <- is.finite(log_post)
  f_mean <- f_median <- f_upper_80 <- f_upper_90 <- f_upper_95 <- NA_real_
  delta_mean <- delta_upper_95 <- NA_real_
  if (any(finite)) {
    lw <- log_post[finite] - max(log_post[finite])
    w <- exp(lw)
    grid_f <- grid[finite]
    f_mean <- sum(grid_f * w) / sum(w)
    qs <- cohort_transition_quantile(grid_f, w, c(0.5, 0.8, 0.9, 0.95))
    f_median <- qs[1]
    f_upper_80 <- qs[2]
    f_upper_90 <- qs[3]
    f_upper_95 <- qs[4]
    parent_mean <- sum(parent_fit * path_weights)
    delta_grid <- grid_f - parent_mean
    delta_mean <- sum(delta_grid * w) / sum(w)
    delta_upper_95 <- cohort_transition_quantile(delta_grid, w, 0.95)
  }
  direct_se <- estimate_scalar_objective_se(
    objective_fn = direct_objective,
    optimum = if (is.finite(f_map)) f_map else mean(local_interval),
    search_interval = local_interval,
    se_floor = max(sd_floor, diff(local_interval) / 100)
  )
  prior_information <- sum(effective_path_weights / (prior_sd^2), na.rm = TRUE)
  patient_likelihood_information <- if (is.finite(direct_se) && direct_se > 0) 1 / (direct_se^2) else 0
  borrowing <- prior_information / (prior_information + patient_likelihood_information)
  if (!is.finite(borrowing)) borrowing <- NA_real_
  rows <- lapply(seq_len(n_parents), function(idx) {
    pr <- priors[[idx]]
    data.frame(
      karyotype = child_name,
      parent_karyotype = parent_karyotypes[idx],
      child_karyotype = child_name,
      transition_group = pr$group,
      prior_group_used = pr$prior_group_used,
      fallback_group_used = pr$fallback_group_used,
      parent_fitness = parent_fit[idx],
      cohort_delta_mu = prior_mu[idx],
      cohort_delta_sd = prior_sd[idx],
      f_map = f_map,
      f_mean = f_mean,
      f_median = f_median,
      f_upper_80 = f_upper_80,
      f_upper_90 = f_upper_90,
      f_upper_95 = f_upper_95,
      delta_map = f_map - parent_fit[idx],
      delta_mean = f_mean - parent_fit[idx],
      delta_upper_95 = f_upper_95 - parent_fit[idx],
      posterior_delta_map = f_map - parent_fit[idx],
      posterior_delta_mean = f_mean - parent_fit[idx],
      posterior_delta_upper_95 = f_upper_95 - parent_fit[idx],
      posterior_fitness_map = f_map,
      posterior_fitness_upper_95 = f_upper_95,
      child_observed_count = child_observed_count,
      child_is_zero = child_is_zero,
      projected_exposure = expected_parent_like,
      expected_count_parent_like = zero_info$expected_count_parent_like,
      zero_informativeness_score = zero_info$zero_informativeness_score,
      path_responsibility = path_weights[idx],
      cohort_prior_weight_multiplier = prior_weight_multiplier,
      cohort_borrowing_fraction = borrowing,
      patient_likelihood_fraction = if (is.finite(borrowing)) 1 - borrowing else NA_real_,
      prior_dominated_flag = isTRUE(child_is_zero) && is.finite(borrowing) && borrowing > 0.8,
      cohort_prior_dominated_flag = isTRUE(child_is_zero) && is.finite(borrowing) && borrowing > 0.8,
      non_identifiable_zero_flag = non_identifiable_zero,
      borrowing_fraction_uses_curvature_proxy = TRUE,
      stringsAsFactors = FALSE
    )
  })
  list(f_map = f_map, diagnostics = do.call(rbind, rows))
}

#' Decide whether a cohort-transition overlay should be applied to a NN node
#'
#' @param item NN child context object.
#' @param child_name Child karyotype ID.
#' @param nn_present Optional logical indicating whether the child is observed.
#' @param two_shell_node_diagnostics Optional one-row two-shell diagnostics.
#' @param cohort_transition_apply_to Overlay application mode.
#' @return A list with apply decision, skip reason, and information flags.
#' @export
should_apply_cohort_transition_to_node <- function(item,
                                                   child_name,
                                                   nn_present = NULL,
                                                   two_shell_node_diagnostics = NULL,
                                                   cohort_transition_apply_to = c("zero_only", "low_information", "all")) {
  cohort_transition_apply_to <- match.arg(cohort_transition_apply_to)
  child_observed_count <- sum(item$child_obs, na.rm = TRUE)
  child_is_zero <- is.finite(child_observed_count) && child_observed_count <= 0
  if (!is.null(nn_present) && length(nn_present) == 1L && isTRUE(!nn_present)) {
    child_is_zero <- TRUE
  }
  boundary <- prior_dominated <- FALSE
  if (is.data.frame(two_shell_node_diagnostics) && nrow(two_shell_node_diagnostics)) {
    if ("objective_boundary_flag" %in% names(two_shell_node_diagnostics)) {
      boundary <- isTRUE(two_shell_node_diagnostics$objective_boundary_flag[1])
    }
    if ("prior_dominated_flag" %in% names(two_shell_node_diagnostics)) {
      prior_dominated <- isTRUE(two_shell_node_diagnostics$prior_dominated_flag[1])
    }
  }
  low_information <- child_is_zero || boundary || prior_dominated
  if (identical(cohort_transition_apply_to, "all")) {
    return(list(apply = TRUE, reason = NA_character_, child_is_zero = child_is_zero, low_information = low_information))
  }
  if (identical(cohort_transition_apply_to, "zero_only")) {
    if (child_is_zero) {
      return(list(apply = TRUE, reason = NA_character_, child_is_zero = child_is_zero, low_information = low_information))
    }
    return(list(apply = FALSE, reason = "observed_nn_skipped_by_zero_only", child_is_zero = child_is_zero, low_information = low_information))
  }
  if (low_information) {
    return(list(apply = TRUE, reason = NA_character_, child_is_zero = child_is_zero, low_information = low_information))
  }
  list(apply = FALSE, reason = "sufficient_patient_information", child_is_zero = child_is_zero, low_information = low_information)
}

#' Apply v2 group-level cohort-transition overlay to one NN child
#'
#' @param item NN child context object.
#' @param child_name Child karyotype ID.
#' @param build_opt_fc Objective builder used by `solve_fitness_bootstrap()`.
#' @param search_interval Numeric optimization interval.
#' @param prior_use Patient-specific cohort transition prior.
#' @param f_two_shell_baseline Patient-specific two-shell baseline fitness.
#' @param nn_present Optional logical indicating whether the child is observed.
#' @param two_shell_node_diagnostics Optional one-row two-shell diagnostics.
#' @param cohort_transition_apply_to Overlay application mode.
#' @param cohort_transition_lambda Global borrowing multiplier.
#' @param cohort_transition_max_borrowing_fraction Maximum borrowing fraction.
#' @param cohort_transition_max_abs_delta_shift Optional shift cap.
#' @param cohort_transition_sd_floor,cohort_transition_patient_sd_floor SD
#'   floors for prior and patient heterogeneity.
#' @return A list with final fitness and node diagnostics.
#' @export
apply_cohort_transition_overlay <- function(item,
                                            child_name,
                                            build_opt_fc,
                                            search_interval,
                                            prior_use,
                                            f_two_shell_baseline,
                                            nn_present = NULL,
                                            two_shell_node_diagnostics = NULL,
                                            cohort_transition_apply_to = c("zero_only", "low_information", "all"),
                                            cohort_transition_lambda = 0.25,
                                            cohort_transition_max_borrowing_fraction = 0.5,
                                            cohort_transition_max_abs_delta_shift = NULL,
                                            cohort_transition_sd_floor = 0.05,
                                            cohort_transition_patient_sd_floor = 0.10) {
  cohort_transition_apply_to <- match.arg(cohort_transition_apply_to)
  validate_nonnegative_finite(cohort_transition_lambda, "cohort_transition_lambda")
  validate_probability(cohort_transition_max_borrowing_fraction, "cohort_transition_max_borrowing_fraction", upper_inclusive = TRUE)
  validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
  validate_positive_finite(cohort_transition_patient_sd_floor, "cohort_transition_patient_sd_floor")
  if (!is.null(cohort_transition_max_abs_delta_shift)) {
    validate_positive_finite(cohort_transition_max_abs_delta_shift, "cohort_transition_max_abs_delta_shift")
  }

  direct_objective <- build_opt_fc(item, do_prior_param = FALSE)
  n_parents <- length(item$parent_fitness)
  if (n_parents == 0L) {
    return(list(f_final = f_two_shell_baseline, diagnostics = data.frame()))
  }
  parent_karyotypes <- item$nj
  if (is.null(parent_karyotypes) || length(parent_karyotypes) != n_parents ||
      any(is.na(parent_karyotypes)) || any(!nzchar(parent_karyotypes))) {
    parent_karyotypes <- names(item$parent_fitness)
  }
  path_weights <- normalize_nn_weights(item$parent_opportunity_weights, fallback_n = n_parents)
  path_weights[!is.finite(path_weights) | path_weights < 0] <- 0
  if (sum(path_weights) <= 0) path_weights <- rep(1 / n_parents, n_parents)
  priors <- lapply(seq_len(n_parents), function(idx) {
    lookup_cohort_transition_group_prior(prior_use, parent_karyotypes[idx], child_name)
  })
  parent_fit <- as.numeric(item$parent_fitness)
  patient_shift <- prior_use$patient_delta_shift
  if (!is.finite(patient_shift)) patient_shift <- 0
  prior_mu <- vapply(priors, `[[`, numeric(1), "mu") + patient_shift
  prior_sd <- pmax(vapply(priors, `[[`, numeric(1), "effective_prior_sd"), cohort_transition_sd_floor, cohort_transition_patient_sd_floor)
  class_lambda <- vapply(priors, `[[`, numeric(1), "cohort_lambda")
  effect_class <- vapply(priors, `[[`, character(1), "effect_class")
  use_for_zero <- vapply(priors, `[[`, logical(1), "use_for_zero")
  use_for_observed <- vapply(priors, `[[`, logical(1), "use_for_observed")
  use_for_low_information <- vapply(priors, `[[`, logical(1), "use_for_low_information")
  expected_parent_like <- as.numeric(item$projected_exposure)
  if (length(expected_parent_like) != 1L || !is.finite(expected_parent_like)) expected_parent_like <- NA_real_
  zero_info <- compute_zero_informativeness_score(expected_parent_like)
  selector <- should_apply_cohort_transition_to_node(
    item = item,
    child_name = child_name,
    nn_present = nn_present,
    two_shell_node_diagnostics = two_shell_node_diagnostics,
    cohort_transition_apply_to = cohort_transition_apply_to
  )
  child_observed_count <- sum(item$child_obs, na.rm = TRUE)
  child_is_zero <- selector$child_is_zero
  non_identifiable_zero <- isTRUE(child_is_zero) &&
    (!is.finite(expected_parent_like) || expected_parent_like < 0.5)
  direct_se <- estimate_scalar_objective_se(
    objective_fn = direct_objective,
    optimum = if (is.finite(f_two_shell_baseline)) f_two_shell_baseline else mean(search_interval),
    search_interval = search_interval,
    se_floor = max(cohort_transition_sd_floor, cohort_transition_patient_sd_floor)
  )
  if (!is.finite(f_two_shell_baseline)) {
    opt <- run_optimise_checked(
      direct_objective,
      interval = search_interval,
      context = sprintf("optimise nearest-neighbour direct baseline for cohort overlay child %s", child_name)
    )
    f_two_shell_baseline <- if (is.null(opt)) NA_real_ else opt$minimum
  }
  anchor_sd <- max(cohort_transition_sd_floor, cohort_transition_patient_sd_floor, direct_se, na.rm = TRUE)
  if (!is.finite(anchor_sd) || anchor_sd <= 0) anchor_sd <- max(cohort_transition_sd_floor, cohort_transition_patient_sd_floor)
  anchor_info <- if (is.finite(f_two_shell_baseline)) 1 / anchor_sd^2 else 0
  zero_multiplier <- if (isTRUE(child_is_zero)) {
    pmin(1, pmax(0, zero_info$zero_informativeness_score[1]))
  } else {
    1
  }
  if (!is.finite(zero_multiplier)) zero_multiplier <- 0
  if (isTRUE(non_identifiable_zero)) zero_multiplier <- 0
  patient_reliability <- prior_use$patient_delta_shift_reliability
  if (!is.finite(patient_reliability)) patient_reliability <- 0
  patient_reliability_multiplier <- pmax(0.25, 1 - 0.5 * pmin(1, patient_reliability))
  class_allows <- if (isTRUE(child_is_zero)) {
    use_for_zero
  } else if (identical(cohort_transition_apply_to, "all")) {
    rep(TRUE, length(use_for_observed))
  } else {
    use_for_low_information
  }
  class_allows[effect_class %in% c("high_variable", "sparse_unknown")] <- FALSE
  effective_lambda <- cohort_transition_lambda * class_lambda * zero_multiplier * patient_reliability_multiplier
  effective_lambda[!class_allows] <- 0
  effective_lambda[!is.finite(effective_lambda) | effective_lambda < 0] <- 0
  prior_info_path <- effective_lambda * path_weights / (prior_sd^2)
  prior_info_path[!is.finite(prior_info_path) | prior_info_path < 0] <- 0
  prior_info <- sum(prior_info_path)
  prior_targets <- parent_fit + prior_mu
  f_overlay <- f_final <- f_two_shell_baseline
  guardrail_hit <- FALSE
  skipped_reason <- selector$reason
  update_applied <- FALSE
  if (isTRUE(selector$apply) && prior_info > 0 && is.finite(f_two_shell_baseline) && anchor_info > 0) {
    f_overlay <- (anchor_info * f_two_shell_baseline + sum(prior_info_path * prior_targets)) /
      (anchor_info + prior_info)
    borrowing <- prior_info / (prior_info + anchor_info)
    if (!is.finite(borrowing)) borrowing <- NA_real_
    max_shift <- cohort_transition_max_abs_delta_shift
    if (is.null(max_shift)) {
      max_shift <- max(cohort_transition_patient_sd_floor, min(0.25, stats::median(prior_sd, na.rm = TRUE)))
    }
    shift <- f_overlay - f_two_shell_baseline
    if (is.finite(shift) && abs(shift) > max_shift) {
      f_overlay <- f_two_shell_baseline + sign(shift) * max_shift
      guardrail_hit <- TRUE
    }
    if (is.finite(borrowing) && borrowing > cohort_transition_max_borrowing_fraction) {
      f_final <- f_two_shell_baseline
      guardrail_hit <- TRUE
      skipped_reason <- "borrowing_fraction_guardrail"
    } else {
      f_final <- f_overlay
      update_applied <- is.finite(f_final) && abs(f_final - f_two_shell_baseline) > sqrt(.Machine$double.eps)
      if (!isTRUE(update_applied)) skipped_reason <- "overlay_shift_negligible"
    }
  } else if (isTRUE(selector$apply) && prior_info <= 0) {
    skipped_reason <- "group_class_disallows_borrowing"
  } else if (isTRUE(selector$apply) && isTRUE(non_identifiable_zero)) {
    skipped_reason <- "non_identifiable_low_exposure_zero"
  }
  borrowing <- if (prior_info + anchor_info > 0) prior_info / (prior_info + anchor_info) else NA_real_
  if (!is.finite(borrowing)) borrowing <- NA_real_
  rows <- lapply(seq_len(n_parents), function(idx) {
    pr <- priors[[idx]]
    data.frame(
      karyotype = child_name,
      parent_karyotype = parent_karyotypes[idx],
      child_karyotype = child_name,
      child_is_zero = child_is_zero,
      child_observed_count = child_observed_count,
      expected_count_parent_like = zero_info$expected_count_parent_like,
      zero_informativeness_score = zero_info$zero_informativeness_score,
      transition_group = pr$group,
      prior_group_used = pr$prior_group_used,
      fallback_group_used = pr$fallback_group_used,
      effect_class = pr$effect_class,
      heterogeneity_class = pr$heterogeneity_class,
      cohort_lambda = pr$cohort_lambda,
      effective_lambda = effective_lambda[idx],
      cohort_delta_mu = prior_mu[idx],
      cohort_delta_sd = pr$sd,
      effective_prior_sd = prior_sd[idx],
      patient_delta_shift = patient_shift,
      patient_delta_shift_n_records = prior_use$patient_delta_shift_n_records,
      patient_delta_shift_reliability = prior_use$patient_delta_shift_reliability,
      parent_fitness = parent_fit[idx],
      f_two_shell_baseline = f_two_shell_baseline,
      f_cohort_overlay = f_overlay,
      f_final = f_final,
      f_delta_from_two_shell = f_final - f_two_shell_baseline,
      delta_two_shell_baseline = f_two_shell_baseline - parent_fit[idx],
      delta_cohort_overlay = f_overlay - parent_fit[idx],
      delta_final = f_final - parent_fit[idx],
      f_map = f_final,
      f_mean = f_final,
      f_median = f_final,
      f_upper_80 = NA_real_,
      f_upper_90 = NA_real_,
      f_upper_95 = NA_real_,
      delta_map = f_final - parent_fit[idx],
      delta_mean = f_final - parent_fit[idx],
      delta_upper_95 = NA_real_,
      posterior_delta_map = f_final - parent_fit[idx],
      posterior_delta_mean = f_final - parent_fit[idx],
      posterior_delta_upper_95 = NA_real_,
      posterior_fitness_map = f_final,
      posterior_fitness_upper_95 = NA_real_,
      projected_exposure = expected_parent_like,
      path_responsibility = path_weights[idx],
      cohort_update_applied = update_applied,
      cohort_update_skipped_reason = if (is.na(skipped_reason)) NA_character_ else skipped_reason,
      guardrail_hit = guardrail_hit,
      cohort_borrowing_fraction = borrowing,
      patient_likelihood_fraction = if (is.finite(borrowing)) 1 - borrowing else NA_real_,
      prior_dominated_flag = isTRUE(child_is_zero) && is.finite(borrowing) && borrowing > 0.8,
      cohort_prior_dominated_flag = isTRUE(child_is_zero) && is.finite(borrowing) && borrowing > cohort_transition_max_borrowing_fraction,
      non_identifiable_zero_flag = non_identifiable_zero,
      borrowing_fraction_uses_curvature_proxy = TRUE,
      class_warning_flags = pr$class_warning_flags,
      stringsAsFactors = FALSE
    )
  })
  list(f_final = f_final, diagnostics = do.call(rbind, rows))
}

#' Refit one patient with a cohort transition prior
#'
#' @param patient One patient input accepted by `alfak()`.
#' @param patient_id Patient identifier.
#' @param outdir Output directory for this patient refit.
#' @param cohort_transition_prior Prior object from `learn_cohort_transition_prior()`.
#' @param ... Additional arguments passed to `alfak()`.
#' @return The invisible return value from `alfak()`.
#' @export
refit_patient_with_cohort_transition_prior <- function(patient,
                                                       patient_id,
                                                       outdir,
                                                       cohort_transition_prior,
                                                       ...) {
  alfak(
    yi = patient,
    outdir = outdir,
    nn_prior = "cohort_transition",
    cohort_transition_prior = cohort_transition_prior,
    cohort_transition_patient_id = patient_id,
    ...
  )
}

#' Fit patients with a cohort-informed transition prior
#'
#' `alfak_cohort_transition()` learns a cohort-level prior on CNA transition
#' effects, Delta fitness = child fitness - parent fitness, from upstream
#' patient-specific two-shell results. It then refits every patient separately.
#' Raw patient count matrices are never concatenated and no pooled absolute
#' cohort fitness landscape is estimated.
#'
#' @details
#' The upstream two-shell cache is resolved as
#' `<two_shell_root>/pm_<xxxx>/MINIOBS<xxx>/<sample_name>/`. For example,
#' `existing_two_shell_results/pm_0.00005/MINIOBS20/patient_A/` should contain
#' `bootstrap_res.Rds`, `landscape.Rds`, and `nn_prior_diagnostics.Rds`, with
#' `landscape_posterior_samples.Rds` preferred when available. If
#' `two_shell_root` is `NULL`, the same layout is created under
#' `file.path(outdir, "two_shell_base")`.
#'
#' Valid cached two-shell fits are reused. A missing patient directory triggers
#' a rerun only for that patient when `rerun_missing_two_shell = TRUE`. A corrupt
#' or incomplete patient directory is backed up with a `__corrupt_<timestamp>`
#' suffix and rerun only for that patient when `rerun_corrupt_two_shell = TRUE`.
#'
#' The default contextual cohort model is fit on transition effects, not
#' absolute fitness: `Delta = child fitness - parent fitness`. It builds an
#' evidence bank from high-confidence observed/frequent transitions and matches
#' each target NN by parent copy-number profile shape, total copy-number area,
#' CNA burden, changed chromosome, direction, local copy state, and event
#' similarity. Raw bootstrap/path records are not treated as independent
#' patients, zero NN are censoring evidence only, and low-exposure zero NN remain
#' non-identifiable. High-variable or sparse contexts keep the patient-specific
#' two-shell baseline and receive uncertainty labels instead of forced cohort
#' imputation.
#'
#' Version `"v2"` remains available as a group-level fallback/comparison mode.
#' It aggregates raw bootstrap/path records within patient and transition group,
#' classifies groups for cross-patient consistency, and borrows only weakly from
#' supported groups. With leave-one-patient-out enabled, patient `p` is refit
#' using evidence from the other patients; this avoids borrowing that patient's
#' own two-shell transition effects back into its refit.
#'
#' Patient-level diagnostics include cohort borrowing fractions and flags for
#' prior-dominated or non-identifiable zero nearest neighbours. A
#' cohort-prior-dominated zero should be interpreted as a cohort-informed upper
#' constraint, not as a precise patient-specific fitness measurement.
#'
#' @param patients Named list of patient inputs accepted by `alfak()`.
#' @param outdir Output directory for cohort diagnostics and patient refits.
#' @param patient_ids Patient IDs. Defaults to `names(patients)`.
#' @param two_shell_root Optional root containing existing two-shell results.
#' @param two_shell_pm,two_shell_minobs Upstream two-shell PM and MINIOBS values.
#' @param two_shell_pm_tag,two_shell_minobs_tag Optional exact cache directory tags.
#' @param two_shell_sample_map Optional named patient-to-sample directory map.
#' @param reuse_two_shell Reuse valid two-shell fit directories.
#' @param rerun_missing_two_shell Rerun only missing two-shell fits.
#' @param rerun_corrupt_two_shell Back up and rerun only corrupt two-shell fits.
#' @param two_shell_integrity_check Integrity mode for cached fits.
#' @param base_nn_prior Upstream prior mode used when rerunning base fits.
#' @param minobs,nboot,n0,nb,pm,passage_times,allow_noninteger_counts,correct_efflux Arguments forwarded to `alfak()`.
#' @param cohort_transition_grouping Transition grouping mode.
#' @param cohort_transition_version Cohort-transition implementation version.
#'   `"contextual"` is the default context-aware evidence-bank overlay; `"v2"`
#'   keeps the group-level heterogeneity-aware selective-borrowing overlay.
#' @param cohort_transition_apply_to Which NN nodes can receive the v2 overlay.
#'   The default `"zero_only"` leaves observed NN estimates at the
#'   patient-specific two-shell baseline.
#' @param cohort_transition_overlay_base Baseline used by v2, by default
#'   `"empirical_two_shell"`.
#' @param cohort_transition_lambda Global multiplier for v2 cohort borrowing.
#' @param cohort_transition_max_borrowing_fraction Maximum borrowing fraction
#'   allowed before a v2 update is skipped.
#' @param cohort_transition_max_abs_delta_shift Optional maximum absolute change
#'   from the two-shell baseline.
#' @param cohort_contextual_apply_to Which NN nodes can receive contextual
#'   overlay updates. If `NULL`, inherits `cohort_transition_apply_to`.
#' @param cohort_contextual_overlay_base Baseline used by contextual mode.
#' @param cohort_context_baseline_ploidy Baseline ploidy used to compute CNA
#'   burden.
#' @param cohort_context_chromosome_weights Optional chromosome weights for
#'   profile distances.
#' @param cohort_context_profile_transform Copy-number profile representation:
#'   mass, centered, zscore, or raw.
#' @param cohort_context_profile_distance Profile distance used by contextual
#'   kernel matching.
#' @param cohort_context_event_match Event matching rule for contextual
#'   neighbors.
#' @param cohort_context_bandwidth_profile,cohort_context_bandwidth_area,cohort_context_bandwidth_burden,cohort_context_bandwidth_local,cohort_context_bandwidth_event
#'   Optional contextual kernel bandwidths; `NULL` uses adaptive values.
#' @param cohort_context_profile_weight,cohort_context_area_weight,cohort_context_burden_weight,cohort_context_local_weight,cohort_context_event_weight
#'   Component weights in the contextual distance.
#' @param cohort_context_min_patients,cohort_context_min_effective_n,cohort_context_min_effective_patients,cohort_context_min_unique_children
#'   Context support thresholds.
#' @param cohort_context_k_nearest,cohort_context_min_kernel_weight Neighbor
#'   retention controls for contextual lookup.
#' @param cohort_context_effect_threshold,cohort_context_sign_consistency_threshold,cohort_context_high_weighted_sd,cohort_context_high_between_patient_sd,cohort_context_high_i2
#'   Context classification thresholds.
#' @param cohort_context_lambda Global contextual borrowing multiplier.
#' @param cohort_context_lambda_consistent_deleterious,cohort_context_lambda_consistent_neutral,cohort_context_lambda_consistent_beneficial,cohort_context_lambda_high_variable,cohort_context_lambda_sparse_unknown,cohort_context_lambda_conflicting_zero
#'   Class-specific contextual borrowing strengths.
#' @param cohort_context_sd_floor,cohort_context_patient_sd_floor Contextual
#'   prior SD and patient heterogeneity floors.
#' @param cohort_context_sd_multiplier_consistent_deleterious,cohort_context_sd_multiplier_consistent_neutral,cohort_context_sd_multiplier_consistent_beneficial,cohort_context_sd_multiplier_high_variable,cohort_context_sd_multiplier_sparse_unknown,cohort_context_sd_multiplier_conflicting_zero
#'   Class-specific contextual prior SD multipliers.
#' @param cohort_context_zero_as_censoring_only Treat contextual zero evidence
#'   as censoring only, never as observed Delta labels.
#' @param cohort_context_zero_min_expected_count Minimum parent-like expected
#'   count for contextual zero censoring evidence.
#' @param cohort_context_zero_weight_cap_ratio Cap on contextual zero evidence.
#' @param cohort_context_max_borrowing_fraction Maximum contextual borrowing
#'   fraction before an update is skipped.
#' @param cohort_context_max_abs_delta_shift Optional maximum contextual shift
#'   from the two-shell baseline.
#' @param cohort_context_keep_baseline_when_sparse,cohort_context_keep_baseline_when_high_variable
#'   Keep the two-shell baseline for sparse or high-variable contexts.
#' @param cohort_context_leave_one_patient_out Exclude the target patient from
#'   contextual evidence lookup during refit.
#' @param cohort_transition_leave_one_patient_out Store LOO priors and use them
#'   during patient refits.
#' @param cohort_transition_use_zero Whether informative zero NN records are used
#'   as censoring evidence.
#' @param cohort_transition_zero_min_exposure Optional explicit zero exposure threshold.
#' @param cohort_transition_zero_min_expected_count Minimum Delta=0 expected count
#'   for informative zero records.
#' @param cohort_transition_zero_weight_cap_ratio Cap on total zero evidence weight.
#' @param cohort_transition_zero_expected_count_cap Cap applied inside the
#'   zero-censoring likelihood to keep very high-exposure zeros from dominating
#'   observed transitions.
#' @param cohort_transition_zero_mean_shift_cap Maximum absolute shift that zero
#'   censoring can impose on a group mean away from observed transition effects.
#' @param cohort_transition_min_patients_per_group Minimum patients per group.
#' @param cohort_transition_min_effective_n Minimum effective evidence per group.
#' @param cohort_transition_sd_floor Minimum transition prior SD.
#' @param cohort_transition_patient_sd_floor Patient heterogeneity SD floor.
#' @param cohort_transition_global_fallback Whether group priors can fall back to global.
#' @param cohort_transition_use_prior_dominated_records,cohort_transition_use_boundary_records Whether such records can train v2 priors.
#' @param cohort_transition_max_delta_se,cohort_transition_max_delta_se_quantile Delta-SE screens for v2 observed records.
#' @param cohort_transition_min_path_responsibility Minimum path responsibility.
#' @param cohort_transition_min_observed_count Minimum observed child count.
#' @param cohort_transition_classify_groups Whether v2 groups are classified.
#' @param cohort_transition_min_patients_consistent,cohort_transition_min_effective_patients,cohort_transition_min_effective_observed Minimum patient-level evidence thresholds.
#' @param cohort_transition_effect_threshold,cohort_transition_sign_consistency_threshold,cohort_transition_high_heterogeneity_i2,cohort_transition_high_between_patient_sd,cohort_transition_context_heterogeneity_drop Group classification controls.
#' @param cohort_transition_lambda_consistent_deleterious,cohort_transition_lambda_consistent_neutral,cohort_transition_lambda_consistent_beneficial,cohort_transition_lambda_context_dependent,cohort_transition_lambda_high_variable,cohort_transition_lambda_sparse_unknown,cohort_transition_lambda_global_fallback Class-specific v2 borrowing strengths.
#' @param cohort_transition_sd_multiplier_consistent_deleterious,cohort_transition_sd_multiplier_consistent_neutral,cohort_transition_sd_multiplier_consistent_beneficial,cohort_transition_sd_multiplier_context_dependent,cohort_transition_sd_multiplier_high_variable,cohort_transition_sd_multiplier_sparse_unknown,cohort_transition_sd_multiplier_global_fallback Class-specific v2 SD multipliers.
#' @param cohort_transition_patient_shift,cohort_transition_patient_shift_min_records,cohort_transition_patient_shift_shrinkage_sd Patient-specific transition shift controls.
#' @param cohort_transition_zero_as_censoring_only Treat zeros as censoring
#'   evidence only, never fake observed delta labels.
#' @param cohort_transition_save_diagnostics Save cohort diagnostic RDS files.
#' @param ... Additional arguments forwarded to `alfak()`.
#' @return Invisibly, a list with status tables, records, prior, diagnostics, and patient output paths.
#' @export
#' @examples
#' \dontrun{
#' patients <- list(
#'   patient_A = list(x = counts_A, dt = 1),
#'   patient_B = list(x = counts_B, dt = 1),
#'   patient_C = list(x = counts_C, dt = 1)
#' )
#'
#' alfak_cohort_transition(
#'   patients = patients,
#'   outdir = "cohort_transition_fit",
#'   two_shell_root = "existing_two_shell_results",
#'   pm = 0.00005,
#'   minobs = 20,
#'   two_shell_pm_tag = "pm_0.00005",
#'   two_shell_minobs_tag = "MINIOBS20",
#'   base_nn_prior = "empirical_two_shell",
#'   cohort_transition_grouping = "gain_loss_chr",
#'   cohort_transition_leave_one_patient_out = TRUE
#' )
#' }
alfak_cohort_transition <- function(patients,
                                    outdir,
                                    patient_ids = names(patients),
                                    two_shell_root = NULL,
                                    two_shell_pm = pm,
                                    two_shell_minobs = minobs,
                                    two_shell_pm_tag = NULL,
                                    two_shell_minobs_tag = NULL,
                                    two_shell_sample_map = NULL,
                                    reuse_two_shell = TRUE,
                                    rerun_missing_two_shell = TRUE,
                                    rerun_corrupt_two_shell = TRUE,
                                    two_shell_integrity_check = c("strict", "basic", "none"),
                                    base_nn_prior = "empirical_two_shell",
                                    minobs = 20,
                                    nboot = 45,
                                    n0 = 1e5,
                                    nb = 1e7,
                                    pm = 0.00005,
                                    passage_times = NULL,
                                    allow_noninteger_counts = FALSE,
                                    correct_efflux = FALSE,
                                    cohort_transition_grouping = c("gain_loss", "gain_loss_chr", "gain_loss_chr_burden", "exact_event"),
                                    cohort_transition_version = c("contextual", "v2", "v1"),
                                    cohort_transition_apply_to = c("zero_only", "low_information", "all"),
                                    cohort_transition_overlay_base = c("empirical_two_shell", "direct"),
                                    cohort_transition_lambda = 0.25,
                                    cohort_transition_max_borrowing_fraction = 0.5,
                                    cohort_transition_max_abs_delta_shift = NULL,
                                    cohort_contextual_apply_to = NULL,
                                    cohort_contextual_overlay_base = c("empirical_two_shell", "direct"),
                                    cohort_context_baseline_ploidy = 2,
                                    cohort_context_chromosome_weights = NULL,
                                    cohort_context_profile_transform = c("mass", "centered", "zscore", "raw"),
                                    cohort_context_profile_distance = c("hellinger", "jensen_shannon", "cosine", "euclidean", "manhattan"),
                                    cohort_context_event_match = c("same_chr_direction", "same_direction", "kernel"),
                                    cohort_context_bandwidth_profile = NULL,
                                    cohort_context_bandwidth_area = NULL,
                                    cohort_context_bandwidth_burden = NULL,
                                    cohort_context_bandwidth_local = NULL,
                                    cohort_context_bandwidth_event = NULL,
                                    cohort_context_profile_weight = 1.0,
                                    cohort_context_area_weight = 0.5,
                                    cohort_context_burden_weight = 0.5,
                                    cohort_context_local_weight = 1.0,
                                    cohort_context_event_weight = 2.0,
                                    cohort_context_min_patients = 3L,
                                    cohort_context_min_effective_n = 5,
                                    cohort_context_min_effective_patients = 3,
                                    cohort_context_min_unique_children = 3L,
                                    cohort_context_k_nearest = 50,
                                    cohort_context_min_kernel_weight = 1e-6,
                                    cohort_context_effect_threshold = 0.02,
                                    cohort_context_sign_consistency_threshold = 0.75,
                                    cohort_context_high_weighted_sd = 0.10,
                                    cohort_context_high_between_patient_sd = 0.10,
                                    cohort_context_high_i2 = 0.50,
                                    cohort_context_lambda = 0.25,
                                    cohort_context_lambda_consistent_deleterious = 0.50,
                                    cohort_context_lambda_consistent_neutral = 0.25,
                                    cohort_context_lambda_consistent_beneficial = 0.10,
                                    cohort_context_lambda_high_variable = 0.00,
                                    cohort_context_lambda_sparse_unknown = 0.00,
                                    cohort_context_lambda_conflicting_zero = 0.00,
                                    cohort_context_sd_floor = 0.05,
                                    cohort_context_patient_sd_floor = 0.10,
                                    cohort_context_sd_multiplier_consistent_deleterious = 1.0,
                                    cohort_context_sd_multiplier_consistent_neutral = 1.5,
                                    cohort_context_sd_multiplier_consistent_beneficial = 2.5,
                                    cohort_context_sd_multiplier_high_variable = 4.0,
                                    cohort_context_sd_multiplier_sparse_unknown = 4.0,
                                    cohort_context_sd_multiplier_conflicting_zero = 4.0,
                                    cohort_context_zero_as_censoring_only = TRUE,
                                    cohort_context_zero_min_expected_count = 3.0,
                                    cohort_context_zero_weight_cap_ratio = 0.25,
                                    cohort_context_max_borrowing_fraction = 0.5,
                                    cohort_context_max_abs_delta_shift = NULL,
                                    cohort_context_keep_baseline_when_sparse = TRUE,
                                    cohort_context_keep_baseline_when_high_variable = TRUE,
                                    cohort_context_leave_one_patient_out = TRUE,
                                    cohort_transition_leave_one_patient_out = TRUE,
                                    cohort_transition_use_zero = TRUE,
                                    cohort_transition_zero_min_exposure = NULL,
                                    cohort_transition_zero_min_expected_count = 3.0,
                                    cohort_transition_zero_weight_cap_ratio = 0.25,
                                    cohort_transition_zero_expected_count_cap = 10.0,
                                    cohort_transition_zero_mean_shift_cap = 0.2,
                                    cohort_transition_min_patients_per_group = 2L,
                                    cohort_transition_min_effective_n = 3L,
                                    cohort_transition_sd_floor = 0.05,
                                    cohort_transition_patient_sd_floor = 0.10,
                                    cohort_transition_global_fallback = TRUE,
                                    cohort_transition_use_prior_dominated_records = FALSE,
                                    cohort_transition_use_boundary_records = FALSE,
                                    cohort_transition_max_delta_se = NULL,
                                    cohort_transition_max_delta_se_quantile = 0.75,
                                    cohort_transition_min_path_responsibility = 0.05,
                                    cohort_transition_min_observed_count = 1L,
                                    cohort_transition_classify_groups = TRUE,
                                    cohort_transition_min_patients_consistent = 3L,
                                    cohort_transition_min_effective_patients = 3,
                                    cohort_transition_min_effective_observed = 3,
                                    cohort_transition_effect_threshold = 0.02,
                                    cohort_transition_sign_consistency_threshold = 0.75,
                                    cohort_transition_high_heterogeneity_i2 = 0.50,
                                    cohort_transition_high_between_patient_sd = 0.10,
                                    cohort_transition_context_heterogeneity_drop = 0.25,
                                    cohort_transition_lambda_consistent_deleterious = 0.50,
                                    cohort_transition_lambda_consistent_neutral = 0.25,
                                    cohort_transition_lambda_consistent_beneficial = 0.15,
                                    cohort_transition_lambda_context_dependent = 0.30,
                                    cohort_transition_lambda_high_variable = 0.00,
                                    cohort_transition_lambda_sparse_unknown = 0.00,
                                    cohort_transition_lambda_global_fallback = 0.05,
                                    cohort_transition_sd_multiplier_consistent_deleterious = 1.0,
                                    cohort_transition_sd_multiplier_consistent_neutral = 1.5,
                                    cohort_transition_sd_multiplier_consistent_beneficial = 2.0,
                                    cohort_transition_sd_multiplier_context_dependent = 1.5,
                                    cohort_transition_sd_multiplier_high_variable = 4.0,
                                    cohort_transition_sd_multiplier_sparse_unknown = 4.0,
                                    cohort_transition_sd_multiplier_global_fallback = 4.0,
                                    cohort_transition_patient_shift = TRUE,
                                    cohort_transition_patient_shift_min_records = 3L,
                                    cohort_transition_patient_shift_shrinkage_sd = 0.10,
                                    cohort_transition_zero_as_censoring_only = TRUE,
                                    cohort_transition_save_diagnostics = TRUE,
                                    ...) {
  two_shell_integrity_check <- match.arg(two_shell_integrity_check)
  cohort_transition_grouping <- match.arg(cohort_transition_grouping)
  cohort_transition_version <- match.arg(cohort_transition_version)
  cohort_transition_apply_to <- match.arg(cohort_transition_apply_to)
  cohort_transition_overlay_base <- match.arg(cohort_transition_overlay_base)
  if (is.null(cohort_contextual_apply_to)) {
    cohort_contextual_apply_to <- cohort_transition_apply_to
  } else {
    cohort_contextual_apply_to <- match.arg(cohort_contextual_apply_to, c("zero_only", "low_information", "all"))
  }
  cohort_contextual_overlay_base <- match.arg(cohort_contextual_overlay_base)
  cohort_context_profile_transform <- match.arg(cohort_context_profile_transform)
  cohort_context_profile_distance <- match.arg(cohort_context_profile_distance)
  cohort_context_event_match <- match.arg(cohort_context_event_match)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(patient_ids)) {
    stop("`patient_ids` must be supplied when `patients` is unnamed.", call. = FALSE)
  }
  patient_ids <- as.character(patient_ids)
  if (length(patients) != length(patient_ids)) {
    stop("`patients` and `patient_ids` must have the same length.", call. = FALSE)
  }
  names(patients) <- patient_ids
  if (is.null(two_shell_root)) {
    two_shell_root <- file.path(outdir, "two_shell_base")
  }

  two_shell_status <- ensure_two_shell_fits(
    patients = patients,
    patient_ids = patient_ids,
    two_shell_root = two_shell_root,
    outdir = outdir,
    pm = two_shell_pm,
    minobs = two_shell_minobs,
    nboot = nboot,
    n0 = n0,
    nb = nb,
    passage_times = passage_times,
    allow_noninteger_counts = allow_noninteger_counts,
    correct_efflux = correct_efflux,
    pm_tag = two_shell_pm_tag,
    minobs_tag = two_shell_minobs_tag,
    sample_map = two_shell_sample_map,
    reuse_two_shell = reuse_two_shell,
    rerun_missing_two_shell = rerun_missing_two_shell,
    rerun_corrupt_two_shell = rerun_corrupt_two_shell,
    integrity_check = two_shell_integrity_check,
    base_nn_prior = base_nn_prior,
    ...
  )

  records <- extract_cohort_transition_records(
    fit_dirs = two_shell_status$fit_dir,
    patient_ids = two_shell_status$patient_id,
    pm = two_shell_pm,
    grouping = cohort_transition_grouping,
    cohort_transition_use_zero = cohort_transition_use_zero,
    cohort_transition_zero_min_expected_count = cohort_transition_zero_min_expected_count,
    cohort_transition_zero_min_exposure = cohort_transition_zero_min_exposure
  )
  prior <- learn_cohort_transition_prior(
    records = records,
    leave_one_patient_out = if (identical(cohort_transition_version, "contextual")) {
      cohort_context_leave_one_patient_out
    } else {
      cohort_transition_leave_one_patient_out
    },
    grouping = cohort_transition_grouping,
    cohort_transition_version = cohort_transition_version,
    cohort_transition_min_patients_per_group = cohort_transition_min_patients_per_group,
    cohort_transition_min_effective_n = cohort_transition_min_effective_n,
    cohort_transition_sd_floor = cohort_transition_sd_floor,
    cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor,
    cohort_transition_global_fallback = cohort_transition_global_fallback,
    cohort_transition_zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio,
    cohort_transition_zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
    cohort_transition_zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap,
    cohort_transition_use_prior_dominated_records = cohort_transition_use_prior_dominated_records,
    cohort_transition_use_boundary_records = cohort_transition_use_boundary_records,
    cohort_transition_max_delta_se = cohort_transition_max_delta_se,
    cohort_transition_max_delta_se_quantile = cohort_transition_max_delta_se_quantile,
    cohort_transition_min_path_responsibility = cohort_transition_min_path_responsibility,
    cohort_transition_min_observed_count = cohort_transition_min_observed_count,
    cohort_transition_classify_groups = cohort_transition_classify_groups,
    cohort_transition_min_patients_consistent = cohort_transition_min_patients_consistent,
    cohort_transition_min_effective_patients = cohort_transition_min_effective_patients,
    cohort_transition_min_effective_observed = cohort_transition_min_effective_observed,
    cohort_transition_effect_threshold = cohort_transition_effect_threshold,
    cohort_transition_sign_consistency_threshold = cohort_transition_sign_consistency_threshold,
    cohort_transition_high_heterogeneity_i2 = cohort_transition_high_heterogeneity_i2,
    cohort_transition_high_between_patient_sd = cohort_transition_high_between_patient_sd,
    cohort_transition_context_heterogeneity_drop = cohort_transition_context_heterogeneity_drop,
    cohort_transition_lambda_consistent_deleterious = cohort_transition_lambda_consistent_deleterious,
    cohort_transition_lambda_consistent_neutral = cohort_transition_lambda_consistent_neutral,
    cohort_transition_lambda_consistent_beneficial = cohort_transition_lambda_consistent_beneficial,
    cohort_transition_lambda_context_dependent = cohort_transition_lambda_context_dependent,
    cohort_transition_lambda_high_variable = cohort_transition_lambda_high_variable,
    cohort_transition_lambda_sparse_unknown = cohort_transition_lambda_sparse_unknown,
    cohort_transition_lambda_global_fallback = cohort_transition_lambda_global_fallback,
    cohort_transition_sd_multiplier_consistent_deleterious = cohort_transition_sd_multiplier_consistent_deleterious,
    cohort_transition_sd_multiplier_consistent_neutral = cohort_transition_sd_multiplier_consistent_neutral,
    cohort_transition_sd_multiplier_consistent_beneficial = cohort_transition_sd_multiplier_consistent_beneficial,
    cohort_transition_sd_multiplier_context_dependent = cohort_transition_sd_multiplier_context_dependent,
    cohort_transition_sd_multiplier_high_variable = cohort_transition_sd_multiplier_high_variable,
    cohort_transition_sd_multiplier_sparse_unknown = cohort_transition_sd_multiplier_sparse_unknown,
    cohort_transition_sd_multiplier_global_fallback = cohort_transition_sd_multiplier_global_fallback,
    cohort_transition_patient_shift = cohort_transition_patient_shift,
    cohort_transition_patient_shift_min_records = cohort_transition_patient_shift_min_records,
    cohort_transition_patient_shift_shrinkage_sd = cohort_transition_patient_shift_shrinkage_sd,
    cohort_transition_zero_as_censoring_only = cohort_transition_zero_as_censoring_only,
    cohort_transition_zero_min_expected_count = cohort_transition_zero_min_expected_count,
    cohort_context_baseline_ploidy = cohort_context_baseline_ploidy,
    cohort_context_chromosome_weights = cohort_context_chromosome_weights,
    cohort_context_profile_transform = cohort_context_profile_transform,
    cohort_context_profile_distance = cohort_context_profile_distance,
    cohort_context_event_match = cohort_context_event_match,
    cohort_context_bandwidth_profile = cohort_context_bandwidth_profile,
    cohort_context_bandwidth_area = cohort_context_bandwidth_area,
    cohort_context_bandwidth_burden = cohort_context_bandwidth_burden,
    cohort_context_bandwidth_local = cohort_context_bandwidth_local,
    cohort_context_bandwidth_event = cohort_context_bandwidth_event,
    cohort_context_profile_weight = cohort_context_profile_weight,
    cohort_context_area_weight = cohort_context_area_weight,
    cohort_context_burden_weight = cohort_context_burden_weight,
    cohort_context_local_weight = cohort_context_local_weight,
    cohort_context_event_weight = cohort_context_event_weight,
    cohort_context_min_patients = cohort_context_min_patients,
    cohort_context_min_effective_n = cohort_context_min_effective_n,
    cohort_context_min_effective_patients = cohort_context_min_effective_patients,
    cohort_context_min_unique_children = cohort_context_min_unique_children,
    cohort_context_k_nearest = cohort_context_k_nearest,
    cohort_context_min_kernel_weight = cohort_context_min_kernel_weight,
    cohort_context_effect_threshold = cohort_context_effect_threshold,
    cohort_context_sign_consistency_threshold = cohort_context_sign_consistency_threshold,
    cohort_context_high_weighted_sd = cohort_context_high_weighted_sd,
    cohort_context_high_between_patient_sd = cohort_context_high_between_patient_sd,
    cohort_context_high_i2 = cohort_context_high_i2,
    cohort_context_lambda_consistent_deleterious = cohort_context_lambda_consistent_deleterious,
    cohort_context_lambda_consistent_neutral = cohort_context_lambda_consistent_neutral,
    cohort_context_lambda_consistent_beneficial = cohort_context_lambda_consistent_beneficial,
    cohort_context_lambda_high_variable = cohort_context_lambda_high_variable,
    cohort_context_lambda_sparse_unknown = cohort_context_lambda_sparse_unknown,
    cohort_context_lambda_conflicting_zero = cohort_context_lambda_conflicting_zero,
    cohort_context_sd_floor = cohort_context_sd_floor,
    cohort_context_patient_sd_floor = cohort_context_patient_sd_floor,
    cohort_context_sd_multiplier_consistent_deleterious = cohort_context_sd_multiplier_consistent_deleterious,
    cohort_context_sd_multiplier_consistent_neutral = cohort_context_sd_multiplier_consistent_neutral,
    cohort_context_sd_multiplier_consistent_beneficial = cohort_context_sd_multiplier_consistent_beneficial,
    cohort_context_sd_multiplier_high_variable = cohort_context_sd_multiplier_high_variable,
    cohort_context_sd_multiplier_sparse_unknown = cohort_context_sd_multiplier_sparse_unknown,
    cohort_context_sd_multiplier_conflicting_zero = cohort_context_sd_multiplier_conflicting_zero,
    cohort_context_zero_as_censoring_only = cohort_context_zero_as_censoring_only,
    cohort_context_zero_min_expected_count = cohort_context_zero_min_expected_count,
    cohort_context_zero_weight_cap_ratio = cohort_context_zero_weight_cap_ratio
  )
  diagnostics <- prior$diagnostics
  diagnostics$two_shell_root <- two_shell_root
  diagnostics$pm_tag <- unique(two_shell_status$pm_tag)
  diagnostics$minobs_tag <- unique(two_shell_status$minobs_tag)

  if (isTRUE(cohort_transition_save_diagnostics)) {
    saveRDS(records, file.path(outdir, "cohort_transition_records.Rds"))
    saveRDS(prior, file.path(outdir, "cohort_transition_prior.Rds"))
    if (identical(prior$version, "cohort_transition_v2")) {
      saveRDS(prior$patient_group_summaries, file.path(outdir, "cohort_transition_patient_group_summaries.Rds"))
      saveRDS(prior$group_classes, file.path(outdir, "cohort_transition_group_classes.Rds"))
      saveRDS(prior, file.path(outdir, "cohort_transition_prior_v2.Rds"))
    }
    if (identical(prior$version, "cohort_transition_contextual_v1")) {
      saveRDS(prior$evidence_bank, file.path(outdir, "cohort_context_evidence_bank.Rds"))
      saveRDS(prior$zero_evidence_bank, file.path(outdir, "cohort_context_zero_evidence_bank.Rds"))
      saveRDS(prior$context_feature_config, file.path(outdir, "cohort_context_feature_config.Rds"))
      saveRDS(prior$context_bandwidths, file.path(outdir, "cohort_context_bandwidths.Rds"))
      saveRDS(prior$diagnostics, file.path(outdir, "cohort_context_diagnostics.Rds"))
      if (is.list(prior$v2_fallback_prior) && identical(prior$v2_fallback_prior$version, "cohort_transition_v2")) {
        saveRDS(prior$v2_fallback_prior, file.path(outdir, "cohort_transition_prior_v2.Rds"))
        saveRDS(prior$v2_fallback_prior$group_classes, file.path(outdir, "cohort_transition_group_classes.Rds"))
      }
    }
    saveRDS(diagnostics, file.path(outdir, "cohort_transition_diagnostics.Rds"))
  }

  patient_outdirs <- stats::setNames(file.path(outdir, patient_ids), patient_ids)
  refit_status <- lapply(patient_ids, function(patient_id) {
    patient_outdir <- patient_outdirs[[patient_id]]
    res <- tryCatch(
      {
        refit_patient_with_cohort_transition_prior(
          patient = patients[[patient_id]],
          patient_id = patient_id,
          outdir = patient_outdir,
          cohort_transition_prior = prior,
          minobs = minobs,
          nboot = nboot,
          n0 = n0,
          nb = nb,
          pm = pm,
          passage_times = passage_times,
          allow_noninteger_counts = allow_noninteger_counts,
          correct_efflux = correct_efflux,
          cohort_transition_version = cohort_transition_version,
          cohort_transition_apply_to = cohort_transition_apply_to,
          cohort_transition_overlay_base = cohort_transition_overlay_base,
          cohort_transition_lambda = cohort_transition_lambda,
          cohort_transition_max_borrowing_fraction = cohort_transition_max_borrowing_fraction,
          cohort_transition_max_abs_delta_shift = cohort_transition_max_abs_delta_shift,
          cohort_transition_sd_floor = cohort_transition_sd_floor,
          cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor,
          cohort_contextual_apply_to = cohort_contextual_apply_to,
          cohort_contextual_overlay_base = cohort_contextual_overlay_base,
          cohort_context_lambda = cohort_context_lambda,
          cohort_context_max_borrowing_fraction = cohort_context_max_borrowing_fraction,
          cohort_context_max_abs_delta_shift = cohort_context_max_abs_delta_shift,
          cohort_context_sd_floor = cohort_context_sd_floor,
          cohort_context_patient_sd_floor = cohort_context_patient_sd_floor,
          cohort_context_keep_baseline_when_sparse = cohort_context_keep_baseline_when_sparse,
          cohort_context_keep_baseline_when_high_variable = cohort_context_keep_baseline_when_high_variable,
          ...
        )
        list(ok = TRUE, error_message = NA_character_, xval = res)
      },
      error = function(e) list(ok = FALSE, error_message = conditionMessage(e), xval = NA_real_)
    )
    data.frame(
      patient_id = patient_id,
      outdir = patient_outdir,
      ok = res$ok,
      error_message = res$error_message,
      stringsAsFactors = FALSE
    )
  })
  refit_status <- do.call(rbind, refit_status)
  saveRDS(refit_status, file.path(outdir, "cohort_transition_refit_status.Rds"))

  invisible(list(
    two_shell_status = two_shell_status,
    records = records,
    prior = prior,
    diagnostics = diagnostics,
    refit_status = refit_status,
    patient_outdirs = patient_outdirs
  ))
}
