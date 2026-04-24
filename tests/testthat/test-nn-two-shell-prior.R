test_that("empirical_two_shell prior mode and controls validate", {
  expect_identical(
    alfakR:::validate_nn_prior_mode("empirical_two_shell"),
    "empirical_two_shell"
  )
  expect_silent(
    alfakR:::validate_nn_prior_controls(
      nn_two_shell_min_delta_n = 3L,
      nn_two_shell_min_exposure = NULL,
      nn_two_shell_min_observed_count = 1L,
      nn_two_shell_max_weight_ratio = 1,
      nn_two_shell_lambda = 1,
      nn_two_shell_reuse_sd = NULL,
      nn_two_shell_uncertainty_floor = NULL
    )
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_min_delta_n = 0L),
    "`nn_two_shell_min_delta_n`"
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_min_exposure = -1),
    "`nn_two_shell_min_exposure`"
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_min_observed_count = 1.5),
    "`nn_two_shell_min_observed_count`"
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_max_weight_ratio = -0.1),
    "`nn_two_shell_max_weight_ratio`"
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_lambda = -0.1),
    "`nn_two_shell_lambda`"
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_reuse_sd = -0.1),
    "`nn_two_shell_reuse_sd`"
  )
  expect_error(
    alfakR:::validate_nn_prior_controls(nn_two_shell_uncertainty_floor = -0.1),
    "`nn_two_shell_uncertainty_floor`"
  )
})

test_that("new two-shell controls do not affect existing none mode", {
  x <- matrix(
    c(80, 60,
      4, 2,
      0, 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("2.2.2", "2.2.3", "2.2.4"), c("0", "1"))
  )
  yi <- list(x = x, dt = 1)

  set.seed(1101)
  baseline <- suppressWarnings(alfakR:::solve_fitness_bootstrap(
    yi,
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4,
    nn_prior = "none"
  ))
  set.seed(1101)
  with_controls <- suppressWarnings(alfakR:::solve_fitness_bootstrap(
    yi,
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4,
    nn_prior = "none",
    nn_two_shell_min_delta_n = 9L,
    nn_two_shell_min_exposure = 100,
    nn_two_shell_min_observed_count = 3L,
    nn_two_shell_max_weight_ratio = 0.2,
    nn_two_shell_lambda = 0.1,
    nn_two_shell_reuse_sd = 10,
    nn_two_shell_uncertainty_floor = 5
  ))

  expect_equal(with_controls$final_fitness, baseline$final_fitness, tolerance = 1e-12)
  expect_equal(with_controls$nn_fitness, baseline$nn_fitness, tolerance = 1e-12)
})

test_that("empirical_two_shell records fallback when no two-step candidates are retained", {
  x <- matrix(
    c(80, 60),
    nrow = 1,
    byrow = TRUE,
    dimnames = list("2.2.2", c("0", "1"))
  )
  res <- suppressWarnings(alfakR:::solve_fitness_bootstrap(
    list(x = x, dt = 1),
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4,
    nn_prior = "empirical_two_shell",
    nn_two_shell_min_exposure = 1e9
  ))

  expect_s3_class(res$nn_prior_diagnostics, "data.frame")
  expect_identical(res$nn_prior_diagnostics$fallback_reason[1], "no_retained_two_step_candidates")
  expect_equal(res$nn_prior_diagnostics$n_2step_candidates_retained[1], 0)
  expect_true("nn_two_shell_node_diagnostics" %in% names(res))
})

test_that("two-shell path responsibilities normalize shared descendants", {
  paths <- data.frame(
    one_step = c("2.2.1", "2.3.2", "2.2.3"),
    descendant = c("2.3.1", "2.3.1", "2.2.4"),
    transition_probability = c(0.2, 0.8, 0.5),
    parent_anchor_exposure = c(10, 5, 4),
    expected_exposure_path = c(2, 4, 2),
    observed_count = c(3, 3, 1),
    expected_exposure = c(6, 6, 2),
    stringsAsFactors = FALSE
  )
  out <- alfakR:::compute_two_shell_path_responsibilities(paths)

  expect_equal(
    sum(out$path_responsibility[out$descendant == "2.3.1"]),
    1,
    tolerance = 1e-12
  )
  expect_equal(
    out$path_responsibility[out$one_step == "2.2.1"],
    (10 * 0.2) / (10 * 0.2 + 5 * 0.8),
    tolerance = 1e-12
  )
  expect_equal(
    out$path_responsibility[out$descendant == "2.2.4"],
    1,
    tolerance = 1e-12
  )
})

make_two_shell_test_context <- function() {
  list(
    child = list(
      ni = "2.2.1",
      nj = "2.2.2",
      pij = 0,
      parent_fitness = 0,
      parent_birth_times = 0,
      parent_birth_fallback = FALSE,
      parent_opportunity_weights = 1,
      parent_xfit = matrix(c(1, 1), nrow = 1),
      child_obs = c(0, 0),
      ntot = c(100, 100),
      parent_fitness_mean_pij = 0,
      parent_fitness_mean_exposure = 0,
      projected_exposure = 0
    )
  )
}

test_that("two-shell backward correction changes a one-step estimate in the outward direction", {
  contexts <- make_two_shell_test_context()
  outward_paths <- data.frame(
    one_step = "child",
    descendant = "2.2.0",
    f2_hat = -2,
    f2_var = 0,
    outward_weight = 1,
    stringsAsFactors = FALSE
  )

  corrected <- alfakR:::apply_two_shell_backward_correction(
    nn_child_contexts = contexts,
    f1_initial = c(child = 0),
    outward_paths = outward_paths,
    mu01 = 0,
    sigma01 = 1,
    mu12 = 0,
    sigma12 = 0.25,
    tau_reuse = 0,
    nn_two_shell_lambda = 1,
    timepoints = c(0, 1),
    search_interval = c(-3, 3)
  )

  expect_lt(corrected$f1["child"], 0)
  expect_true(corrected$node_diagnostics$outward_weight_sum > 0)
})

test_that("larger two-shell reuse uncertainty weakens the outward correction", {
  contexts <- make_two_shell_test_context()
  outward_paths <- data.frame(
    one_step = "child",
    descendant = "2.2.0",
    f2_hat = -2,
    f2_var = 0,
    outward_weight = 1,
    stringsAsFactors = FALSE
  )

  low_reuse <- alfakR:::apply_two_shell_backward_correction(
    nn_child_contexts = contexts,
    f1_initial = c(child = 0),
    outward_paths = outward_paths,
    mu01 = 0,
    sigma01 = 1,
    mu12 = 0,
    sigma12 = 0.25,
    tau_reuse = 0,
    nn_two_shell_lambda = 1,
    timepoints = c(0, 1),
    search_interval = c(-3, 3)
  )
  high_reuse <- alfakR:::apply_two_shell_backward_correction(
    nn_child_contexts = contexts,
    f1_initial = c(child = 0),
    outward_paths = outward_paths,
    mu01 = 0,
    sigma01 = 1,
    mu12 = 0,
    sigma12 = 0.25,
    tau_reuse = 10,
    nn_two_shell_lambda = 1,
    timepoints = c(0, 1),
    search_interval = c(-3, 3)
  )

  expect_lt(abs(high_reuse$f1["child"]), abs(low_reuse$f1["child"]))
})
