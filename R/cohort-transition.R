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
                                              cohort_transition_zero_min_expected_count = 1.0,
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
learn_cohort_transition_prior <- function(records,
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
  if (!is.list(prior) || !identical(prior$version, "cohort_transition_v1")) {
    stop("`cohort_transition_prior` must be a cohort_transition_v1 prior object.", call. = FALSE)
  }
  if (isTRUE(prior$leave_one_patient_out) && length(prior$loo_priors)) {
    if (is.null(cohort_transition_patient_id) || length(cohort_transition_patient_id) != 1L || !nzchar(cohort_transition_patient_id)) {
      stop("`cohort_transition_patient_id` is required when the cohort transition prior contains leave-one-patient-out priors.", call. = FALSE)
    }
  }
  prior
}

cohort_transition_prior_for_patient <- function(prior, patient_id = NULL) {
  if (isTRUE(prior$leave_one_patient_out) && length(prior$loo_priors)) {
    if (is.null(patient_id) || !nzchar(patient_id)) {
      stop("`cohort_transition_patient_id` is required for leave-one-patient-out cohort transition priors.", call. = FALSE)
    }
    if (!patient_id %in% names(prior$loo_priors)) {
      if (!patient_id %in% prior$patient_ids && nrow(prior$global_prior)) {
        return(list(
          grouping = prior$grouping,
          group_priors = prior$group_priors,
          global_prior = prior$global_prior,
          contributing_patients = prior$patient_ids,
          leave_one_patient_out = TRUE,
          leave_one_patient_out_fallback = "patient_has_no_training_records"
        ))
      }
      stop(sprintf("No leave-one-patient-out cohort transition prior is available for patient `%s`.", patient_id), call. = FALSE)
    }
    loo <- prior$loo_priors[[patient_id]]
    if (!nrow(loo$global_prior)) {
      stop(sprintf("Leave-one-patient-out prior for patient `%s` has no contributing records.", patient_id), call. = FALSE)
    }
    return(list(
      grouping = prior$grouping,
      group_priors = loo$group_priors,
      global_prior = loo$global_prior,
      contributing_patients = loo$contributing_patients,
      leave_one_patient_out = TRUE
    ))
  }
  list(
    grouping = prior$grouping,
    group_priors = prior$group_priors,
    global_prior = prior$global_prior,
    contributing_patients = prior$patient_ids,
    leave_one_patient_out = FALSE
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
  }
  row <- row[1L, , drop = FALSE]
  list(
    group = group,
    prior_group_used = row$group[1],
    fallback_group_used = if ("fallback_group" %in% names(row)) row$fallback_group[1] else NA_character_,
    mu = row$mu[1],
    sd = if ("sigma_with_patient_heterogeneity" %in% names(row)) {
      row$sigma_with_patient_heterogeneity[1]
    } else {
      row$sigma[1]
    },
    parsed = parsed
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
#' The cohort model is fit on transition effects, not absolute fitness:
#' `Delta = child fitness - parent fitness`. Informative zero nearest
#' neighbours are used as censoring evidence only when their projected
#' parent-like expected count or explicit exposure threshold is large enough.
#' Low-exposure zeros are excluded from prior fitting, and zero evidence is
#' capped so it cannot dominate observed transitions. With
#' `cohort_transition_leave_one_patient_out = TRUE`, patient `p` is refit with
#' priors learned from the other patients; this avoids borrowing that patient's
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
                                    cohort_transition_leave_one_patient_out = TRUE,
                                    cohort_transition_use_zero = TRUE,
                                    cohort_transition_zero_min_exposure = NULL,
                                    cohort_transition_zero_min_expected_count = 1.0,
                                    cohort_transition_zero_weight_cap_ratio = 1.0,
                                    cohort_transition_zero_expected_count_cap = 10.0,
                                    cohort_transition_zero_mean_shift_cap = 0.2,
                                    cohort_transition_min_patients_per_group = 2L,
                                    cohort_transition_min_effective_n = 3L,
                                    cohort_transition_sd_floor = 1e-3,
                                    cohort_transition_patient_sd_floor = 0.1,
                                    cohort_transition_global_fallback = TRUE,
                                    cohort_transition_save_diagnostics = TRUE,
                                    ...) {
  two_shell_integrity_check <- match.arg(two_shell_integrity_check)
  cohort_transition_grouping <- match.arg(cohort_transition_grouping)
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
    leave_one_patient_out = cohort_transition_leave_one_patient_out,
    grouping = cohort_transition_grouping,
    cohort_transition_min_patients_per_group = cohort_transition_min_patients_per_group,
    cohort_transition_min_effective_n = cohort_transition_min_effective_n,
    cohort_transition_sd_floor = cohort_transition_sd_floor,
    cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor,
    cohort_transition_global_fallback = cohort_transition_global_fallback,
    cohort_transition_zero_weight_cap_ratio = cohort_transition_zero_weight_cap_ratio,
    cohort_transition_zero_expected_count_cap = cohort_transition_zero_expected_count_cap,
    cohort_transition_zero_mean_shift_cap = cohort_transition_zero_mean_shift_cap
  )
  diagnostics <- prior$diagnostics
  diagnostics$two_shell_root <- two_shell_root
  diagnostics$pm_tag <- unique(two_shell_status$pm_tag)
  diagnostics$minobs_tag <- unique(two_shell_status$minobs_tag)

  if (isTRUE(cohort_transition_save_diagnostics)) {
    saveRDS(records, file.path(outdir, "cohort_transition_records.Rds"))
    saveRDS(prior, file.path(outdir, "cohort_transition_prior.Rds"))
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
