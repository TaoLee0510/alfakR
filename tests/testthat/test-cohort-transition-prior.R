make_ct_yi <- function() {
  x <- matrix(
    c(80, 60,
      8, 12),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2.2", "2.2.3"), c("0", "1"))
  )
  list(x = x, dt = 1)
}

make_valid_two_shell_dir <- function(root, patient_id, pm_tag = "pm_0.00005", minobs_tag = "MINIOBS20",
                                     child_observed_count = 3, projected_exposure = 5) {
  fit_dir <- file.path(root, pm_tag, minobs_tag, patient_id)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  final <- matrix(
    c(0, 0.1),
    nrow = 2,
    dimnames = list(NULL, "2.2.2")
  )
  nn <- matrix(
    c(0.2, 0.3),
    nrow = 2,
    dimnames = list(NULL, "2.2.3")
  )
  saveRDS(
    list(
      final_fitness = final,
      nn_fitness = nn,
      initial_fitness = final,
      initial_frequencies = final,
      final_frequencies = final,
      nn_prior_diagnostics = data.frame(nn_prior_mode_used = "empirical_two_shell")
    ),
    file.path(fit_dir, "bootstrap_res.Rds")
  )
  saveRDS(
    data.frame(k = c("2.2.2", "2.2.3"), mean = c(0, 0.25), median = c(0, 0.25), sd = c(0, 0.1)),
    file.path(fit_dir, "landscape.Rds")
  )
  saveRDS(matrix(0, nrow = 2, ncol = 2), file.path(fit_dir, "landscape_posterior_samples.Rds"))
  saveRDS(
    list(
      replicate = data.frame(
        replicate_id = 1:2,
        nn_prior_mode_requested = "empirical_two_shell",
        nn_prior_mode_used = "empirical_two_shell",
        nn_prior_source_used = "two_shell",
        mu01 = 0.1,
        sigma01 = 0.2
      ),
      node = data.frame(
        replicate_id = 1:2,
        karyotype = "2.2.3",
        direct_observed_count = child_observed_count,
        projected_exposure = projected_exposure,
        objective_boundary_flag = FALSE,
        prior_dominated_flag = child_observed_count == 0,
        outward_weight_sum = 0.4
      )
    ),
    file.path(fit_dir, "nn_prior_diagnostics.Rds")
  )
  fit_dir
}

make_ct_records <- function(patient_ids = c("patient_A", "patient_B", "patient_C"),
                            delta = c(0.1, 0.2, 0.3),
                            source_type = "observed",
                            expected = 0) {
  delta <- rep_len(delta, length(patient_ids))
  source_type <- rep_len(source_type, length(patient_ids))
  expected <- rep_len(expected, length(patient_ids))
  rows <- lapply(seq_along(patient_ids), function(i) {
    parsed <- alfakR:::cohort_transition_parse_pair("2.2.2", "2.2.3")
    data.frame(
      patient_id = patient_ids[i],
      parent_karyotype = "2.2.2",
      child_karyotype = "2.2.3",
      transition_chr = parsed$transition_chr,
      transition_direction = parsed$transition_direction,
      transition_size = parsed$transition_size,
      group_gain_loss = parsed$group_gain_loss,
      group_gain_loss_chr = parsed$group_gain_loss_chr,
      group_gain_loss_chr_burden = parsed$group_gain_loss_chr_burden,
      group_exact_event = parsed$group_exact_event,
      transition_group = parsed$group_gain_loss_chr,
      parent_total_cn = parsed$parent_total_cn,
      child_total_cn = parsed$child_total_cn,
      parent_burden = parsed$parent_burden,
      child_burden = parsed$child_burden,
      parent_fitness = 0,
      child_fitness_two_shell = delta[i],
      delta_hat = delta[i],
      delta_se = 0.05,
      child_observed_count = ifelse(source_type[i] == "informative_zero", 0, 2),
      child_is_zero = source_type[i] == "informative_zero",
      projected_exposure = expected[i],
      expected_count_parent_like = expected[i],
      zero_informativeness_score = pmin(1, expected[i] / 3),
      zero_informativeness_category = ifelse(expected[i] >= 3, "informative_zero", "uninformative_zero"),
      boundary_flag = FALSE,
      prior_dominated_flag = FALSE,
      two_shell_used = TRUE,
      two_shell_outward_weight = 0.2,
      path_responsibility = 1,
      replicate_id = 1L,
      bootstrap_id = 1L,
      source_type = source_type[i],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_ct_overlay_item <- function(child_obs = c(0, 0), projected_exposure = 10, parent_fitness = 0) {
  list(
    ni = "2.2.3",
    nj = "2.2.2",
    pij = 1,
    parent_fitness = unname(c(parent_fitness)),
    parent_birth_times = 0,
    parent_birth_fallback = FALSE,
    parent_opportunity_weights = 1,
    parent_xfit = matrix(c(1, 1), nrow = 2),
    child_obs = child_obs,
    ntot = c(100, 100),
    parent_fitness_mean_pij = parent_fitness,
    parent_fitness_mean_exposure = parent_fitness,
    projected_exposure = projected_exposure
  )
}

ct_overlay_builder <- function(item, do_prior_param = FALSE, ...) {
  force(item)
  function(fc) (fc - 0.2)^2
}

test_that("resolve_two_shell_fit_dirs returns expected cache paths", {
  root <- tempfile("two_shell_root_")
  dir.create(file.path(root, "pm_0.00005", "MINIOBS20", "patient_A"), recursive = TRUE)
  dir.create(file.path(root, "pm_0.00005", "MINIOBS20", "patient_B"), recursive = TRUE)

  resolved <- alfakR::resolve_two_shell_fit_dirs(
    two_shell_root = root,
    patient_ids = c("patient_A", "patient_B"),
    pm = 0.00005,
    minobs = 20
  )

  expect_equal(resolved$pm_tag, rep("pm_0.00005", 2))
  expect_equal(resolved$minobs_tag, rep("MINIOBS20", 2))
  expect_equal(resolved$expected_fit_dir, file.path(root, "pm_0.00005", "MINIOBS20", c("patient_A", "patient_B")))
  expect_true(all(resolved$exists))
})

test_that("filter_cohort_transition_records reports exclusion reasons", {
  records <- make_ct_records(
    patient_ids = paste0("patient_", LETTERS[1:6]),
    delta = rep(0.1, 6),
    expected = c(0, 0, 0, 0, 0, 0.1)
  )
  records$prior_dominated_flag[1] <- TRUE
  records$boundary_flag[2] <- TRUE
  records$delta_se[3] <- 10
  records$path_responsibility[4] <- 0.01
  records$delta_hat[5] <- NA_real_
  records$source_type[6] <- "informative_zero"
  records$child_is_zero[6] <- TRUE
  records$child_observed_count[6] <- 0

  filtered <- alfakR::filter_cohort_transition_records(
    records,
    cohort_transition_max_delta_se = 1,
    cohort_transition_min_path_responsibility = 0.05,
    cohort_transition_zero_min_expected_count = 3
  )

  expect_equal(nrow(filtered$kept_records), 0)
  expect_equal(filtered$diagnostics$prior_dominated, 1)
  expect_equal(filtered$diagnostics$boundary_dominated, 1)
  expect_equal(filtered$diagnostics$large_delta_se, 1)
  expect_equal(filtered$diagnostics$low_path_responsibility, 1)
  expect_equal(filtered$diagnostics$non_finite_delta, 1)
  expect_equal(filtered$diagnostics$low_exposure_zero, 1)
})

test_that("bootstrap and path rows aggregate to patient-level evidence", {
  patient_a <- make_ct_records(rep("patient_A", 10), rep(0.1, 10))
  patient_b <- make_ct_records("patient_B", 0.2)
  records <- rbind(patient_a, patient_b)
  filtered <- alfakR::filter_cohort_transition_records(records, cohort_transition_min_path_responsibility = 0)
  summaries <- alfakR::aggregate_cohort_transition_records_by_patient(
    filtered$kept_records,
    grouping = "gain_loss_chr"
  )
  obs <- summaries[summaries$cohort_transition_evidence_type == "observed_delta_evidence", , drop = FALSE]

  expect_equal(nrow(obs), 2)
  expect_equal(sort(obs$patient_id), c("patient_A", "patient_B"))

  classes <- alfakR::classify_cohort_transition_groups(summaries)
  gain_chr <- classes[classes$group == "gain_chr3" & classes$group_level == "gain_loss_chr", , drop = FALSE]
  expect_equal(gain_chr$n_patients_observed, 2)
  expect_lte(gain_chr$effective_patients_observed, 2)
})

test_that("consistent deleterious groups get nonzero conservative borrowing", {
  records <- make_ct_records(
    c("patient_A", "patient_B", "patient_C"),
    c(-0.10, -0.12, -0.11)
  )
  prior <- alfakR::learn_cohort_transition_prior(
    records,
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  row <- prior$group_priors[prior$group_priors$group == "gain_chr3", , drop = FALSE]

  expect_equal(row$effect_class, "consistent_deleterious")
  expect_gt(row$cohort_lambda, 0)
  expect_true(row$recommended_use_for_zero)
})

test_that("high-variable groups do not borrow by default", {
  records <- make_ct_records(
    c("patient_A", "patient_B", "patient_C", "patient_D"),
    c(-0.12, 0.12, -0.10, 0.10)
  )
  prior <- alfakR::learn_cohort_transition_prior(
    records,
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  row <- prior$group_priors[prior$group_priors$group == "gain_chr3", , drop = FALSE]

  expect_true(row$effect_class == "high_variable" || row$heterogeneity_class == "high_variable")
  expect_equal(row$cohort_lambda, 0)
})

test_that("sparse groups do not use strong group-specific priors", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records("patient_A", -0.2),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  row <- prior$group_priors[prior$group_priors$group == "gain_chr3", , drop = FALSE]

  expect_equal(row$effect_class, "sparse_unknown")
  expect_equal(row$cohort_lambda, 0)
})

test_that("observed NN are unchanged by zero-only v2 overlay", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B", "patient_C"), c(-0.1, -0.12, -0.11)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  fit <- alfakR::apply_cohort_transition_overlay(
    item = make_ct_overlay_item(child_obs = c(2, 2), projected_exposure = 10),
    child_name = "2.2.3",
    build_opt_fc = ct_overlay_builder,
    search_interval = c(-1, 1),
    prior_use = alfakR:::cohort_transition_prior_for_patient(prior),
    f_two_shell_baseline = 0.2,
    nn_present = TRUE,
    cohort_transition_apply_to = "zero_only"
  )

  expect_equal(fit$f_final, 0.2)
  expect_false(any(fit$diagnostics$cohort_update_applied))
  expect_equal(unique(fit$diagnostics$cohort_update_skipped_reason), "observed_nn_skipped_by_zero_only")
})

test_that("high-exposure zero in consistent deleterious group receives weak overlay", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B", "patient_C"), c(-0.1, -0.12, -0.11)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  fit <- alfakR::apply_cohort_transition_overlay(
    item = make_ct_overlay_item(child_obs = c(0, 0), projected_exposure = 10),
    child_name = "2.2.3",
    build_opt_fc = ct_overlay_builder,
    search_interval = c(-1, 1),
    prior_use = alfakR:::cohort_transition_prior_for_patient(prior),
    f_two_shell_baseline = 0.2,
    nn_present = FALSE,
    cohort_transition_max_borrowing_fraction = 0.9
  )

  expect_true(any(fit$diagnostics$cohort_update_applied))
  expect_lt(fit$f_final, 0.2)
  expect_gt(unique(fit$diagnostics$effective_lambda), 0)
})

test_that("low-exposure zero is marked non-identifiable and not aggressively updated", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B", "patient_C"), c(-0.1, -0.12, -0.11)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  fit <- alfakR::apply_cohort_transition_overlay(
    item = make_ct_overlay_item(child_obs = c(0, 0), projected_exposure = 0.1),
    child_name = "2.2.3",
    build_opt_fc = ct_overlay_builder,
    search_interval = c(-1, 1),
    prior_use = alfakR:::cohort_transition_prior_for_patient(prior),
    f_two_shell_baseline = 0.2,
    nn_present = FALSE
  )

  expect_equal(fit$f_final, 0.2)
  expect_true(all(fit$diagnostics$non_identifiable_zero_flag))
  expect_false(any(fit$diagnostics$cohort_update_applied))
})

test_that("overlay guardrail caps excessive cohort shifts", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B", "patient_C"), c(-0.4, -0.42, -0.41)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr"
  )
  fit <- alfakR::apply_cohort_transition_overlay(
    item = make_ct_overlay_item(child_obs = c(0, 0), projected_exposure = 10),
    child_name = "2.2.3",
    build_opt_fc = ct_overlay_builder,
    search_interval = c(-1, 1),
    prior_use = alfakR:::cohort_transition_prior_for_patient(prior),
    f_two_shell_baseline = 0.2,
    nn_present = FALSE,
    cohort_transition_max_abs_delta_shift = 0.01,
    cohort_transition_max_borrowing_fraction = 0.99
  )

  expect_true(any(fit$diagnostics$guardrail_hit))
  expect_lte(abs(fit$f_final - 0.2), 0.0101)
})

test_that("ensure_two_shell_fits reuses valid existing two-shell directories", {
  root <- tempfile("two_shell_reuse_")
  outdir <- tempfile("cohort_out_")
  make_valid_two_shell_dir(root, "patient_A")
  make_valid_two_shell_dir(root, "patient_B")
  patients <- list(patient_A = make_ct_yi(), patient_B = make_ct_yi())

  testthat::with_mocked_bindings(
    {
      status <- alfakR::ensure_two_shell_fits(
        patients = patients,
        outdir = outdir,
        two_shell_root = root,
        pm = 0.00005,
        minobs = 20
      )
      expect_equal(status$action, c("reused", "reused"))
      expect_false(any(status$rerun))
    },
    alfak = function(...) stop("two-shell rerun should not be called"),
    .package = "alfakR"
  )
})

test_that("ensure_two_shell_fits reruns only missing patients", {
  root <- tempfile("two_shell_missing_")
  outdir <- tempfile("cohort_out_")
  make_valid_two_shell_dir(root, "patient_A")
  patients <- list(patient_A = make_ct_yi(), patient_B = make_ct_yi())
  called <- new.env(parent = emptyenv())
  called$outdirs <- character(0)

  testthat::with_mocked_bindings(
    {
      status <- alfakR::ensure_two_shell_fits(
        patients = patients,
        outdir = outdir,
        two_shell_root = root,
        pm = 0.00005,
        minobs = 20
      )
      expect_equal(status$action, c("reused", "rerun_missing"))
      expect_equal(basename(called$outdirs), "patient_B")
      expect_true(file.exists(file.path(outdir, "two_shell_fit_status.Rds")))
      expect_true(file.exists(file.path(outdir, "two_shell_fit_status.tsv")))
    },
    alfak = function(yi, outdir, ...) {
      called$outdirs <- c(called$outdirs, outdir)
      make_valid_two_shell_dir(dirname(dirname(dirname(outdir))), basename(outdir))
      invisible(0)
    },
    .package = "alfakR"
  )
})

test_that("ensure_two_shell_fits backs up and reruns only corrupt patients", {
  root <- tempfile("two_shell_corrupt_")
  outdir <- tempfile("cohort_out_")
  make_valid_two_shell_dir(root, "patient_A")
  corrupt_dir <- file.path(root, "pm_0.00005", "MINIOBS20", "patient_B")
  dir.create(corrupt_dir, recursive = TRUE)
  writeLines("not an rds", file.path(corrupt_dir, "bootstrap_res.Rds"))
  saveRDS(data.frame(k = "2.2.2", mean = 0), file.path(corrupt_dir, "landscape.Rds"))
  saveRDS(data.frame(nn_prior_mode_used = "empirical_two_shell"), file.path(corrupt_dir, "nn_prior_diagnostics.Rds"))
  patients <- list(patient_A = make_ct_yi(), patient_B = make_ct_yi())
  called <- new.env(parent = emptyenv())
  called$outdirs <- character(0)

  testthat::with_mocked_bindings(
    {
      status <- alfakR::ensure_two_shell_fits(
        patients = patients,
        outdir = outdir,
        two_shell_root = root,
        pm = 0.00005,
        minobs = 20
      )
      expect_equal(status$action, c("reused", "rerun_corrupt"))
      expect_true(dir.exists(status$backup_dir[2]))
      expect_equal(basename(called$outdirs), "patient_B")
    },
    alfak = function(yi, outdir, ...) {
      called$outdirs <- c(called$outdirs, outdir)
      make_valid_two_shell_dir(dirname(dirname(dirname(outdir))), basename(outdir))
      invisible(0)
    },
    .package = "alfakR"
  )
})

test_that("NULL two_shell_root writes base fits under outdir/two_shell_base", {
  outdir <- tempfile("cohort_null_root_")
  patients <- list(patient_A = make_ct_yi(), patient_B = make_ct_yi())
  called <- new.env(parent = emptyenv())
  called$outdirs <- character(0)

  testthat::with_mocked_bindings(
    {
      status <- alfakR::ensure_two_shell_fits(
        patients = patients,
        outdir = outdir,
        two_shell_root = NULL,
        pm = 0.00005,
        minobs = 20
      )
      expect_true(all(grepl("two_shell_base", status$fit_dir, fixed = TRUE)))
      expect_equal(status$action, c("rerun_missing", "rerun_missing"))
      expect_equal(sort(basename(called$outdirs)), c("patient_A", "patient_B"))
    },
    alfak = function(yi, outdir, ...) {
      called$outdirs <- c(called$outdirs, outdir)
      make_valid_two_shell_dir(dirname(dirname(dirname(outdir))), basename(outdir))
      invisible(0)
    },
    .package = "alfakR"
  )
})

test_that("cohort wrapper refits patients separately and does not pool raw counts", {
  patients <- list(patient_A = make_ct_yi(), patient_B = make_ct_yi())
  outdir <- tempfile("cohort_wrapper_")
  seen <- new.env(parent = emptyenv())
  seen$patient_ids <- character(0)
  seen$nrows <- integer(0)
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B"), c(0.1, 0.2)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )

  testthat::with_mocked_bindings(
    {
      res <- alfakR::alfak_cohort_transition(
        patients = patients,
        outdir = outdir,
        minobs = 20,
        nboot = 1,
        cohort_transition_grouping = "gain_loss_chr"
      )
      expect_equal(seen$patient_ids, c("patient_A", "patient_B"))
      expect_equal(seen$nrows, c(nrow(patients$patient_A$x), nrow(patients$patient_B$x)))
      expect_equal(names(res$patient_outdirs), c("patient_A", "patient_B"))
    },
    ensure_two_shell_fits = function(...) {
      data.frame(
        patient_id = c("patient_A", "patient_B"),
        fit_dir = c("fit_A", "fit_B"),
        pm_tag = "pm_0.00005",
        minobs_tag = "MINIOBS20",
        stringsAsFactors = FALSE
      )
    },
    extract_cohort_transition_records = function(...) make_ct_records(c("patient_A", "patient_B"), c(0.1, 0.2)),
    learn_cohort_transition_prior = function(...) prior,
    refit_patient_with_cohort_transition_prior = function(patient, patient_id, outdir, cohort_transition_prior, ...) {
      seen$patient_ids <- c(seen$patient_ids, patient_id)
      seen$nrows <- c(seen$nrows, nrow(patient$x))
      invisible(0)
    },
    .package = "alfakR"
  )
})

test_that("learn_cohort_transition_prior stores leave-one-patient-out contributors", {
  records <- make_ct_records(c("patient_A", "patient_B", "patient_C"), c(0.1, 0.2, 0.3))
  prior <- alfakR::learn_cohort_transition_prior(
    records,
    leave_one_patient_out = TRUE,
    grouping = "gain_loss_chr",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )

  expect_true(prior$leave_one_patient_out)
  expect_equal(prior$loo_priors$patient_A$contributing_patients, c("patient_B", "patient_C"))
  expect_false("patient_A" %in% prior$loo_priors$patient_A$contributing_patients)
})

test_that("missing leave-one-patient-out prior falls back when patient has no training records", {
  records <- make_ct_records(c("patient_A", "patient_B"), c(0.1, 0.2))
  prior <- alfakR::learn_cohort_transition_prior(
    records,
    leave_one_patient_out = TRUE,
    grouping = "gain_loss_chr",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )

  prior_use <- alfakR:::cohort_transition_prior_for_patient(prior, "patient_without_records")

  expect_true(prior_use$leave_one_patient_out)
  expect_equal(prior_use$leave_one_patient_out_fallback, "patient_has_no_training_records")
  expect_equal(prior_use$contributing_patients, c("patient_A", "patient_B"))
})

test_that("informative zeros are retained as censoring evidence without fake observed deltas", {
  records <- make_ct_records(
    patient_ids = c("patient_A", "patient_B", "patient_C"),
    delta = c(0.2, NA, NA),
    source_type = c("observed", "informative_zero", "informative_zero"),
    expected = c(0, 5, 0.1)
  )
  records <- records[records$source_type != "informative_zero" | records$expected_count_parent_like >= 1, , drop = FALSE]
  prior <- alfakR::learn_cohort_transition_prior(
    records,
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )

  expect_equal(prior$version, "cohort_transition_v2")
  expect_equal(prior$global_prior$n_observed_patient_summaries, 1)
  expect_equal(prior$global_prior$n_zero_patient_summaries, 1)
  expect_equal(prior$diagnostics$n_zero_censoring_records, 1)
  expect_false(prior$diagnostics$zero_likelihood_approximation)
  expect_false(any(records$source_type == "informative_zero" & is.finite(records$delta_hat)))
})

test_that("extreme zero exposures do not create a fake narrow observed-delta prior", {
  observed <- make_ct_records(
    patient_ids = paste0("obs_", seq_len(8)),
    delta = rep(0, 8),
    source_type = "observed",
    expected = 0
  )
  zeros <- make_ct_records(
    patient_ids = paste0("zero_", seq_len(8)),
    delta = rep(NA_real_, 8),
    source_type = "informative_zero",
    expected = 1e6
  )
  records <- rbind(observed, zeros)

  capped <- alfakR::learn_cohort_transition_prior(
    records,
    leave_one_patient_out = FALSE,
    grouping = "gain_loss_chr",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )
  expect_equal(capped$version, "cohort_transition_v2")
  expect_equal(capped$global_prior$n_zero_patient_summaries, 8)
  expect_equal(capped$global_prior$n_observed_patient_summaries, 8)
  expect_lte(abs(capped$global_prior$mu), 1e-8)
  expect_gte(capped$global_prior$effective_prior_sd, 0.1)
})

test_that("low-exposure zero NN gets no cohort-prior pull in patient refit", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B"), c(0.1, 0.2)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )
  prior_use <- alfakR:::cohort_transition_prior_for_patient(prior)
  item <- list(
    ni = "2.2.3",
    nj = "2.2.2",
    pij = 1,
    parent_fitness = unname(c(0)),
    parent_birth_times = 0,
    parent_birth_fallback = FALSE,
    parent_opportunity_weights = 1,
    parent_xfit = matrix(c(1, 1), nrow = 2),
    child_obs = c(0, 0),
    ntot = c(100, 100),
    parent_fitness_mean_pij = 0,
    parent_fitness_mean_exposure = 0,
    projected_exposure = 0.1
  )
  builder <- function(item, do_prior_param = FALSE, ...) {
    force(item)
    function(fc) (fc - 0.2)^2
  }
  fit <- alfakR:::fit_cohort_transition_nn_child(
    item = item,
    child_name = "2.2.3",
    build_opt_fc = builder,
    search_interval = c(-0.05, 0.05),
    prior_use = prior_use
  )

  expect_equal(unique(fit$diagnostics$cohort_prior_weight_multiplier), 0)
  expect_equal(unique(fit$diagnostics$cohort_borrowing_fraction), 0, tolerance = 1e-12)
  expect_true(all(fit$diagnostics$non_identifiable_zero_flag))
  expect_lte(fit$f_map, 0.05)
  expect_equal(fit$f_map, 0.05, tolerance = 2e-3)
})

test_that("cohort-transition path responsibilities sum to one for multiple parents", {
  prior <- alfakR::learn_cohort_transition_prior(
    make_ct_records(c("patient_A", "patient_B"), c(0.1, 0.2)),
    leave_one_patient_out = FALSE,
    grouping = "gain_loss",
    cohort_transition_min_patients_per_group = 1L,
    cohort_transition_min_effective_n = 1
  )
  prior_use <- alfakR:::cohort_transition_prior_for_patient(prior)
  item <- list(
    ni = "2.2.3",
    nj = c("2.2.2", "2.1.3"),
    pij = c(0.7, 0.3),
    parent_fitness = unname(c(0, 0.1)),
    parent_birth_times = c(0, 0),
    parent_birth_fallback = c(FALSE, FALSE),
    parent_opportunity_weights = c(3, 1),
    parent_xfit = matrix(c(1, 1, 1, 1), nrow = 2),
    child_obs = c(0, 0),
    ntot = c(100, 100),
    parent_fitness_mean_pij = 0.05,
    parent_fitness_mean_exposure = 0.025,
    projected_exposure = 10
  )
  builder <- function(item, do_prior_param = FALSE, ...) {
    force(item)
    function(fc) (fc - 0.2)^2
  }
  fit <- alfakR:::fit_cohort_transition_nn_child(
    item = item,
    child_name = "2.2.3",
    build_opt_fc = builder,
    search_interval = c(-1, 1),
    prior_use = prior_use
  )

  expect_equal(sum(fit$diagnostics$path_responsibility), 1, tolerance = 1e-12)
  expect_equal(fit$diagnostics$path_responsibility, c(0.75, 0.25), tolerance = 1e-12)
  expect_equal(fit$diagnostics$parent_karyotype, item$nj)
})

test_that("cohort-transition validation is isolated from existing modes", {
  expect_identical(
    alfakR:::validate_nn_prior_mode("cohort_transition"),
    "cohort_transition"
  )
  yi <- make_ct_yi()
  expect_error(
    suppressWarnings(alfakR:::solve_fitness_bootstrap(
      yi,
      minobs = 20,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4,
      nn_prior = "none",
      cohort_transition_sd_floor = -1,
      cohort_transition_patient_sd_floor = -1
    )),
    NA
  )
  expect_error(
    alfakR:::solve_fitness_bootstrap(
      yi,
      minobs = 20,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4,
      nn_prior = "cohort_transition"
    ),
    "requires `cohort_transition_prior`"
  )
})
