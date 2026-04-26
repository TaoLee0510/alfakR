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

  expect_equal(prior$global_prior$n_observed_records, 1)
  expect_equal(prior$global_prior$n_zero_records, 1)
  expect_true(prior$diagnostics$zero_likelihood_approximation)
  expect_false(any(records$source_type == "informative_zero" & is.finite(records$delta_hat)))
})

test_that("extreme zero exposures are capped in cohort prior fitting", {
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
    cohort_transition_min_effective_n = 1,
    cohort_transition_zero_expected_count_cap = 10
  )
  expect_true(capped$diagnostics$zero_expected_count_capped)
  expect_equal(capped$diagnostics$n_zero_expected_count_capped, 8)
  expect_gte(capped$global_prior$mu, -0.2 - 1e-8)
  expect_lte(abs(capped$global_prior$mu), 0.2 + 1e-8)
  expect_true(grepl("zero_mean_shift_cap", capped$global_prior$warning_flags))
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
