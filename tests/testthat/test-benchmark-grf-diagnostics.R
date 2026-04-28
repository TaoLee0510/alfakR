load_benchmark_nn_diag_env <- function() {
  env <- new.env(parent = globalenv())
  candidates <- normalizePath(
    c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "..", "..")),
    winslash = "/",
    mustWork = FALSE
  )
  repo_dir <- candidates[file.exists(file.path(candidates, "benchmark", "R", "benchmark_nn_diagnostics.R"))][[1]]
  sys.source(file.path(repo_dir, "benchmark", "R", "benchmark_nn_diagnostics.R"), envir = env)
  env
}

test_that("GRF truth helper matches the benchmark formula", {
  bench_env <- load_benchmark_nn_diag_env()
  centroids <- matrix(c(2, 2, 3, 2), ncol = 2, byrow = TRUE)
  lambda <- 1

  got <- bench_env$compute_grf_fitness_truth(c("2.2", "3.2"), centroids, lambda)
  manual <- c(
    "2.2" = (sin(0 / lambda) + sin(1 / lambda)) / (pi * sqrt(2)),
    "3.2" = (sin(1 / lambda) + sin(0 / lambda)) / (pi * sqrt(2))
  )

  expect_equal(got, manual)
})

test_that("GRF ABM output is sampled into an ALFA-K input matrix", {
  bench_env <- load_benchmark_nn_diag_env()
  sim_wide <- data.frame(
    time = c(0, 10, 20),
    check.names = FALSE
  )
  sim_wide[["2.2"]] <- c(100, 50, 0)
  sim_wide[["3.2"]] <- c(0, 50, 100)

  yi <- bench_env$build_nn_grf_yi_from_abm(
    sim_wide = sim_wide,
    training_window = 2,
    sample_depth = 20,
    seed = 1
  )

  expect_equal(dim(yi$x), c(2L, 2L))
  expect_equal(colSums(yi$x), c(`0` = 20, `20` = 20), ignore_attr = TRUE)
  expect_equal(rownames(yi$x), c("2.2", "3.2"))
  expect_equal(yi$dt, 1)
})

test_that("GRF NN accuracy summary reports centered metrics", {
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")
  suppressPackageStartupMessages(library(dplyr))

  bench_env <- load_benchmark_nn_diag_env()
  child_tbl <- tibble::tibble(
    simulation_id = c(1L, 1L),
    lambda = c(0.8, 0.8),
    training_window = c(2L, 2L),
    parameter_label = c("nn_prior_none", "nn_prior_none"),
    nn_prior = c("none", "none"),
    true_fitness = c(-1, 1),
    estimated_fitness = c(-0.5, 0.5),
    observed_in_training = c(TRUE, FALSE),
    estimation_error = c(0.5, -0.5),
    centered_true_fitness = c(-1, 1),
    centered_estimated_fitness = c(-0.5, 0.5),
    centered_error = c(0.5, -0.5),
    bootstrap_sd = c(0.1, 0.2)
  )

  summary_tbl <- bench_env$summarize_nn_grf_child_accuracy(
    child_tbl,
    group_cols = c("lambda", "training_window", "parameter_label", "nn_prior")
  )

  expect_equal(nrow(summary_tbl), 1L)
  expect_equal(summary_tbl$n_children, 2L)
  expect_equal(summary_tbl$centered_rmse, 0.5)
  expect_equal(summary_tbl$pearson, 1)
})
