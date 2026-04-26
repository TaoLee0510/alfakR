make_counts <- function(values, rownames_vec, colnames_vec) {
  x <- matrix(values, nrow = length(rownames_vec), byrow = TRUE)
  rownames(x) <- rownames_vec
  colnames(x) <- colnames_vec
  x
}

reference_project_forward_log <- function(x0, f, timepoints) {
  out <- matrix(NA_real_, nrow = length(x0), ncol = length(timepoints))
  log_x0 <- log(x0)
  for (i in seq_along(timepoints)) {
    lv <- log_x0 + f * timepoints[i]
    denom <- alfakR:::logSumExp(lv)
    out[, i] <- exp(lv - denom)
  }
  out
}

reference_neg_log_lik <- function(param, counts, timepoints) {
  K <- nrow(counts)
  free_idx <- if (K > 1) seq_len(K - 1) else integer(0)
  f_free <- param[free_idx]
  f_full <- c(f_free, -sum(f_free))
  log_x0 <- c(param[K:(2 * K - 2)], 0)
  nll <- 0
  for (i in seq_len(ncol(counts))) {
    lv <- log_x0 + f_full * timepoints[i]
    denom <- alfakR:::logSumExp(lv)
    for (k in seq_len(K)) {
      if (counts[k, i] > 0) {
        nll <- nll - counts[k, i] * (lv[k] - denom)
      }
    }
  }
  nll
}

reference_neighbor_objective <- function(fc_param, parent_fitness, pij_values,
                                         parent_birth_times, timepoints, parent_xfit,
                                         child_obs, ntot, parent_fitness_mean,
                                         prior_mean, prior_sd, do_prior, tol) {
  xc_est <- colSums(do.call(rbind, lapply(seq_along(parent_fitness), function(i) {
    tt <- pmax(0, timepoints - parent_birth_times[i])
    alfakR:::fExp_stable(fc_param, parent_fitness[i], pij_values[i], tt, tol = tol) * parent_xfit[i, ]
  })))
  xc_est <- pmax(0, pmin(1, xc_est))
  res <- stats::dbinom(child_obs, ntot, prob = xc_est, log = TRUE)
  if (do_prior && is.finite(parent_fitness_mean)) {
    res <- c(res, stats::dnorm(fc_param - parent_fitness_mean, mean = prior_mean, sd = prior_sd, log = TRUE))
  }
  res[!is.finite(res)] <- -(10^9)
  -sum(res)
}

make_simple_yi <- function(x, dt = 1) {
  list(x = x, dt = dt)
}

reference_qr_accum <- function(x_trim, dx_dt) {
  K <- nrow(x_trim)
  Q_accum <- matrix(0, nrow = K, ncol = K)
  r_accum <- numeric(K)
  for (t_idx in seq_len(ncol(x_trim))) {
    xt <- x_trim[, t_idx]
    M_t <- diag(as.numeric(xt), nrow = length(xt), ncol = length(xt)) - outer(xt, xt)
    Q_accum <- Q_accum + M_t %*% M_t
    r_accum <- as.numeric(r_accum + M_t %*% dx_dt[, t_idx])
  }
  list(Q_accum = Q_accum, r_accum = r_accum)
}

test_that("minobs includes karyotypes exactly at the threshold", {
  yi <- list(
    x = make_counts(
      c(10, 10,
        15, 15),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  res <- alfakR:::solve_fitness_bootstrap(
    yi,
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4,
    nn_prior = "none"
  )

  expect_setequal(colnames(res$final_fitness), c("2.2.2", "2.2.1"))
})

test_that("matrix-like count inputs are accepted and coerced at entry", {
  yi_df <- list(
    x = data.frame(
      "0" = c(10, 15),
      "1" = c(10, 15),
      row.names = c("2.2.2", "2.2.1"),
      check.names = FALSE
    ),
    dt = 1
  )

  res_df <- alfakR:::solve_fitness_bootstrap(
    yi_df,
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4,
    nn_prior = "none"
  )
  expect_setequal(colnames(res_df$final_fitness), c("2.2.2", "2.2.1"))

  skip_if_not_installed("Matrix")
  x_sparse <- Matrix::Matrix(c(10, 15, 10, 15), nrow = 2, sparse = TRUE)
  rownames(x_sparse) <- c("2.2.2", "2.2.1")
  colnames(x_sparse) <- c("0", "1")

  res_sparse <- alfakR:::solve_fitness_bootstrap(
    list(x = x_sparse, dt = 1),
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4,
    nn_prior = "none"
  )
  expect_setequal(colnames(res_sparse$final_fitness), c("2.2.2", "2.2.1"))
})

test_that("strict karyotype parsing rejects malformed IDs and mixed dimensions", {
  parsed <- alfakR:::parse_karyotype_ids(c("2.2", "2.3"))
  expect_identical(dim(parsed), c(2L, 2L))
  expect_true(is.integer(parsed))
  expect_equal(unname(parsed[1, ]), c(2L, 2L))

  expect_error(alfakR:::parse_karyotype_ids(c("2a.2", "2.3")), "Invalid karyotype ID")
  expect_error(alfakR:::parse_karyotype_ids("2..2"), "Invalid karyotype ID")
  parsed_zero <- alfakR:::parse_karyotype_ids(c("2.0.2", "2.3.2"))
  expect_equal(unname(parsed_zero[1, ]), c(2L, 0L, 2L))
  expect_error(alfakR:::parse_karyotype_ids(c("2.-1.2", "2.3.2")), "Invalid karyotype ID")
  expect_error(alfakR:::parse_karyotype_ids(c("2.2", "2.2.2")), "same number of dot-separated components")
})

test_that("build_W_rcpp validates p, Nmax, and karyotype strings safely", {
  expect_silent(alfakR::build_W_rcpp(c("2.2", "2.3"), p = 0.01, Nmax = Inf))
  expect_silent(alfakR::build_W_rcpp(c("2.0", "2.1"), p = 0.01, Nmax = Inf))
  expect_error(alfakR::build_W_rcpp(c("2.2", "2a.3"), p = 0.01), "Invalid karyotype ID")
  expect_error(alfakR::build_W_rcpp(c("2.2", "2.2.2"), p = 0.01), "same number of dot-separated components")
  expect_error(alfakR::build_W_rcpp(c("2.2", "2.3"), p = NA_real_), "`p`")
  expect_error(alfakR::build_W_rcpp(c("2.2", "2.3"), p = 1.5), "`p`")
})

test_that("solve_fitness_bootstrap rejects malformed karyotype rownames before bootstrapping", {
  yi_bad <- list(
    x = make_counts(
      c(10, 10,
        20, 20),
      rownames_vec = c("2a.2", "2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  expect_error(
    alfakR:::solve_fitness_bootstrap(
      yi_bad,
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4
    ),
    "Invalid karyotype ID"
  )
})

test_that("count-matrix validation rejects invalid values and rounds non-integers once", {
  x_inf <- make_counts(
    c(10, Inf,
      15, 15),
    rownames_vec = c("2.2.2", "2.2.1"),
    colnames_vec = c("0", "1")
  )
  expect_error(
    alfakR:::coerce_count_matrix(x_inf),
    "finite count values"
  )

  x_neg <- make_counts(
    c(10, -1,
      15, 15),
    rownames_vec = c("2.2.2", "2.2.1"),
    colnames_vec = c("0", "1")
  )
  expect_error(
    alfakR:::coerce_count_matrix(x_neg),
    "non-negative count values"
  )

  x_integer_like <- data.frame(
    "0" = c(10 + alfakR:::ALFAK_COUNT_INTEGER_TOL / 2, 15 - alfakR:::ALFAK_COUNT_INTEGER_TOL / 2),
    "1" = c(9 + alfakR:::ALFAK_COUNT_INTEGER_TOL / 3, 14 - alfakR:::ALFAK_COUNT_INTEGER_TOL / 3),
    row.names = c("2.2.2", "2.2.1"),
    check.names = FALSE
  )
  expect_warning(
    rounded_integer_like <- alfakR:::coerce_count_matrix(x_integer_like),
    "Integer-like floating-point values"
  )
  expect_equal(rounded_integer_like, round(as.matrix(x_integer_like)))

  x_non_integer <- data.frame(
    "0" = c(10.2, 15.7),
    "1" = c(9.8, 14.3),
    row.names = c("2.2.2", "2.2.1"),
    check.names = FALSE
  )
  expect_error(
    alfakR:::coerce_count_matrix(x_non_integer),
    "allow_noninteger_counts = TRUE"
  )
  expect_warning(
    rounded <- alfakR:::coerce_count_matrix(x_non_integer, allow_noninteger_counts = TRUE),
    "rounding to the nearest integer once at entry"
  )
  expect_equal(rounded, round(as.matrix(x_non_integer)))

  x_dropped <- c("0" = 10, "1" = 9)
  expect_error(
    alfakR:::coerce_count_matrix(x_dropped),
    "drop = FALSE"
  )
})

test_that("zero-depth timepoints are rejected before normalization or bootstrap", {
  x_zero <- make_counts(
    c(10, 0,
      15, 0),
    rownames_vec = c("2.2.2", "2.2.1"),
    colnames_vec = c("0", "1")
  )

  expect_error(
    alfakR:::solve_fitness_bootstrap(
      make_simple_yi(x_zero),
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4
    ),
    "zero-depth column"
  )

  expect_error(
    alfakR::alfak(
      yi = make_simple_yi(x_zero),
      outdir = tempfile("alfak_zero_depth_"),
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4
    ),
    "zero-depth column"
  )
})

test_that("bootstrap frequent-subset zero columns are normalized to all-zero frequencies", {
  yi <- make_simple_yi(
    make_counts(
      c(1, 1,
        20, 20),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("1", "360")
    )
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      res <- suppressWarnings(
        alfakR:::solve_fitness_bootstrap(
          yi,
          minobs = 30,
          nboot = 1,
          n0 = 1e4,
          nb = 1e6,
          pm = 1e-4
        )
      )
      expect_true(all(is.finite(res$final_fitness)))
    },
    bootstrap_counts = function(data) {
      make_counts(
        c(19, 1,
          0, 20),
        rownames_vec = rownames(data),
        colnames_vec = colnames(data)
      )
    },
    compute_dx_dt = function(x, timepoints) {
      seen$x <- x
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    run_solve_qp_checked = function(...) list(solution = 0),
    optimize_initial_frequencies = function(...) 1,
    joint_optimize = function(...) list(f = 0, x0 = 1),
    find_birth_times = function(...) 0,
    gen_nn_info = function(...) list(),
    .package = "alfakR"
  )

  expect_equal(seen$x[, 1], 0)
  expect_equal(seen$x[, 2], 1)
})

test_that("alfak validates optional arguments before running heavy work", {
  yi <- make_simple_yi(
    make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    )
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = tempfile("alfak_default_passage_"),
        minobs = 1,
        nboot = 1,
        n0 = 1e4 + 0.5,
        nb = 1e6 + 0.5,
        pm = 1e-4
      ))
      expect_identical(returned, 0.1)

      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 0, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = -1, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1.5, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 0, nb = 1e6, pm = 1e-4), "`n0`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = Inf, nb = 1e6, pm = 1e-4), "`n0`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = NA_real_, pm = 1e-4), "`nb`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = Inf), "`pm`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = NA), "`correct_efflux`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = "TRUE"), "`correct_efflux`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = 1), "`correct_efflux`")
    },
    solve_fitness_bootstrap = function(...) {
      seen$solve_called <- TRUE
      list(
        initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        nn_fitness = matrix(numeric(0), nrow = 1, ncol = 0)
      )
    },
    fitKrig = function(...) {
      list(
        summary_stats = data.frame(k = "2.2.2", mean = 0, median = 0, sd = 0, fq = TRUE, nn = FALSE),
        posterior_samples = matrix(0, nrow = 1, ncol = 1),
        krig_stable_mean = NULL,
        krig_stable_median = NULL
      )
    },
    xval = function(...) 0.1,
    .package = "alfakR"
  )

  expect_true(isTRUE(seen$solve_called))
})

test_that("solve_fitness_bootstrap validates bootstrap controls and pm before neighbour generation", {
  yi <- make_simple_yi(
    make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    )
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 0, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 0, nb = 1e6, pm = 1e-4), "`n0`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = NA_real_, pm = 1e-4), "`nb`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = NA_real_), "`pm`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = Inf), "`pm`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = NA), "`correct_efflux`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = "TRUE"), "`correct_efflux`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = 1), "`correct_efflux`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, nn_prior_grid_n = 0), "`nn_prior_grid_n`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, nn_prior_grid_n = 2), "`nn_prior_grid_n` must be at least 3")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, nn_prior_grid_n = 4.5), "`nn_prior_grid_n`")
    },
    gen_nn_info = function(...) {
      seen$gen_nn_called <- TRUE
      list()
    },
    .package = "alfakR"
  )

  expect_false(isTRUE(seen$gen_nn_called))
})

test_that("correct_efflux stops before bootstrap when viability is non-positive", {
  yi <- list(
    x = make_counts(
      c(5, 5),
      rownames_vec = "2.2.2",
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  expect_error(
    alfakR:::solve_fitness_bootstrap(
      yi,
      minobs = 1,
      nboot = 2,
      n0 = 1e4,
      nb = 1e6,
      pm = 0.2,
      correct_efflux = TRUE
    ),
    "correct_efflux viability pre-check failed before bootstrap: pm=0.2"
  )
})

test_that("correct_efflux warns once when viability is positive but tiny", {
  yi <- list(
    x = make_counts(
      c(5, 5,
        6, 6),
      rownames_vec = c("4.4.4.4.4", "2.2.2.2.2"),
      colnames_vec = c("0", "2")
    ),
    dt = 99
  )
  pm <- 1 - ((1 + 5e-7) / 2)^(1 / 20)

  expect_warning(
    alfakR:::solve_fitness_bootstrap(
      yi,
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = pm,
      correct_efflux = TRUE,
      nn_prior = "none",
      passage_times = c(0, 2.5)
    ),
    "0 < viability <"
  )
})

test_that("fExp_stable stays finite and matches the analytic limit near fc == fp", {
  tt <- c(0, 1, 2.5, 5)
  limit_val <- 0.2 * 0.3 * tt

  exact <- alfakR:::fExp_stable(fc_arg = 0.3, fp_arg = 0.3, pij_val = 0.2, tt_arg = tt)
  near <- alfakR:::fExp_stable(fc_arg = 0.3 + 1e-12, fp_arg = 0.3, pij_val = 0.2, tt_arg = tt)

  expect_true(all(is.finite(exact)))
  expect_true(all(is.finite(near)))
  expect_equal(exact, limit_val, tolerance = 1e-12)
  expect_equal(near, limit_val, tolerance = 1e-12)
})

test_that("C++ numerical kernels match the previous R reference calculations", {
  x0 <- c(0.3, 0.7)
  f <- c(0.1, -0.1)
  timepoints <- c(0, 1.5, 4)
  project_ref <- reference_project_forward_log(x0, f, timepoints)
  project_cpp <- alfakR:::alfak_project_forward_log_cpp(x0, f, timepoints)
  expect_equal(project_cpp, project_ref, tolerance = 1e-12)

  zero_cpp <- alfakR:::alfak_project_forward_log_cpp(c(1, 0), c(1, 1), c(0, 1))
  expect_true(all(is.finite(zero_cpp)))
  expect_equal(colSums(zero_cpp), c(1, 1), tolerance = 1e-12)

  counts <- matrix(c(10, 5, 7,
                     4, 11, 9), nrow = 2, byrow = TRUE)
  param <- c(0.2, log(0.4 / 0.6))
  expect_equal(
    alfakR:::alfak_neg_log_lik_cpp(param, counts, timepoints),
    reference_neg_log_lik(param, counts, timepoints),
    tolerance = 1e-12
  )

  parent_fitness <- c(0.15, 0.05)
  pij_values <- c(0.2, 0.1)
  parent_birth_times <- c(-1, 0.5)
  parent_xfit <- matrix(c(0.4, 0.35, 0.3,
                          0.2, 0.25, 0.3), nrow = 2, byrow = TRUE)
  child_obs <- c(0, 2, 3)
  ntot <- c(10, 12, 14)
  expect_equal(
    alfakR:::alfak_neighbor_objective_cpp(
      fc_param = 0.11,
      parent_fitness = parent_fitness,
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      child_obs = child_obs,
      ntot = ntot,
      parent_fitness_mean = 0.12,
      prior_mean = 0.01,
      prior_sd = 0.2,
      do_prior = TRUE,
      tol = alfakR:::ALFAK_FEXP_DELTA_TOL
    ),
    reference_neighbor_objective(
      fc_param = 0.11,
      parent_fitness = parent_fitness,
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      child_obs = child_obs,
      ntot = ntot,
      parent_fitness_mean = 0.12,
      prior_mean = 0.01,
      prior_sd = 0.2,
      do_prior = TRUE,
      tol = alfakR:::ALFAK_FEXP_DELTA_TOL
    ),
    tolerance = 1e-12
  )

  x_trim <- matrix(c(0.3, 0.4,
                     0.7, 0.6), nrow = 2, byrow = TRUE)
  dx_dt <- matrix(c(0.1, -0.05,
                    -0.1, 0.05), nrow = 2, byrow = TRUE)
  qr_ref <- reference_qr_accum(x_trim, dx_dt)
  qr_cpp <- alfakR:::alfak_qr_accum_cpp(x_trim, dx_dt)
  expect_equal(qr_cpp$Q_accum, qr_ref$Q_accum, tolerance = 1e-12)
  expect_equal(qr_cpp$r_accum, qr_ref$r_accum, tolerance = 1e-12)
})

test_that("C++ NN prior helper kernels match R reference paths", {
  timepoints <- c(0, 1.5, 4)
  parent_fitness <- c(0.15, 0.05)
  pij_values <- c(0.2, 0.1)
  parent_birth_times <- c(-1, 0.5)
  parent_xfit <- matrix(c(0.4, 0.35, 0.3,
                          0.2, 0.25, 0.3), nrow = 2, byrow = TRUE)
  child_obs <- c(0, 2, 3)
  ntot <- c(10, 12, 14)
  fc_grid <- seq(-0.2, 0.3, length.out = 7)

  traj_ref <- numeric(length(timepoints))
  for (p in seq_along(parent_fitness)) {
    tt <- pmax(0, timepoints - parent_birth_times[p])
    traj_ref <- traj_ref +
      alfakR:::fExp_stable(0.11, parent_fitness[p], pij_values[p], tt, tol = alfakR:::ALFAK_FEXP_DELTA_TOL) *
      parent_xfit[p, ]
  }
  traj_ref <- pmax(0, pmin(1, traj_ref))
  traj_cpp <- alfakR:::alfak_nn_project_trajectory_cpp(
    0.11, parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit,
    alfakR:::ALFAK_FEXP_DELTA_TOL
  )
  expect_equal(traj_cpp, traj_ref, tolerance = 1e-12)
  expect_equal(
    alfakR:::alfak_nn_project_exposure_cpp(
      0.11, parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit, ntot,
      alfakR:::ALFAK_FEXP_DELTA_TOL
    ),
    sum(ntot * traj_ref),
    tolerance = 1e-12
  )
  expect_equal(
    alfakR:::alfak_parent_opportunity_weights_cpp(pij_values, parent_birth_times, timepoints, parent_xfit, ntot),
    c(
      pij_values[1] * sum(ntot * parent_xfit[1, ] * as.numeric(timepoints >= parent_birth_times[1])),
      pij_values[2] * sum(ntot * parent_xfit[2, ] * as.numeric(timepoints >= parent_birth_times[2]))
    ),
    tolerance = 1e-12
  )

  loglik_cpp <- alfakR:::alfak_neighbor_loglik_grid_cpp(
    fc_grid, parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit,
    child_obs, ntot, alfakR:::ALFAK_FEXP_DELTA_TOL
  )
  loglik_ref <- -vapply(fc_grid, function(fc) {
    alfakR:::alfak_neighbor_objective_cpp(
      fc, parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit,
      child_obs, ntot, parent_fitness_mean = 0.1, prior_mean = NaN, prior_sd = NaN,
      do_prior = FALSE, tol = alfakR:::ALFAK_FEXP_DELTA_TOL
    )
  }, numeric(1))
  expect_equal(loglik_cpp, loglik_ref, tolerance = 1e-12)

  loglik_mat <- rbind(loglik_cpp - max(loglik_cpp), loglik_cpp - max(loglik_cpp) - 0.1)
  log_weights <- rep(log(fc_grid[2] - fc_grid[1]), length(fc_grid))
  parent_means <- c(0.1, 0.2)
  child_weights <- c(1, 0.5)
  marginal_ref <- local({
    total <- 0
    for (i in seq_len(nrow(loglik_mat))) {
      vals <- loglik_mat[i, ] + stats::dnorm(fc_grid - parent_means[i], 0, 0.4, log = TRUE) + log_weights
      total <- total - child_weights[i] * (max(vals) + log(sum(exp(vals - max(vals)))))
    }
    total
  })
  expect_equal(
    alfakR:::alfak_nn_prior_marginal_negloglik_cpp(
      loglik_mat, fc_grid, log_weights, parent_means, child_weights, mu = 0, sigma = 0.4
    ),
    marginal_ref,
    tolerance = 1e-12
  )
})

test_that("C++ karyotype neighbour generation matches R-visible semantics", {
  neigh <- alfakR:::gen_all_neighbours(c("2.2.2", "2.2.3"))
  neigh_str <- apply(neigh, 1, paste, collapse = ".")
  expect_false("2.2.2" %in% neigh_str)
  expect_false("2.2.3" %in% neigh_str)
  expect_true("1.2.2" %in% neigh_str)
  expect_true("2.2.4" %in% neigh_str)

  nn <- alfakR:::gen_nn_info(c("2.2.2", "2.2.4"), pm = 0.00005)
  child_ids <- vapply(nn, `[[`, character(1), "ni")
  expect_true("2.2.3" %in% child_ids)
  mid <- nn[[match("2.2.3", child_ids)]]
  expect_setequal(mid$nj, c("2.2.2", "2.2.4"))
  expect_true(all(is.finite(mid$pij)))
  expect_true(all(mid$pij > 0))
})

test_that("C++ numerical kernels validate dimensions and non-finite inputs", {
  expect_error(
    alfakR:::alfak_project_forward_log_cpp(c(0.5, 0.5), c(0.1), c(0, 1)),
    "same length"
  )
  expect_error(
    alfakR:::alfak_project_forward_log_cpp(c(0, 0), c(0.1, -0.1), c(0, 1)),
    "positive finite value"
  )

  counts_bad <- matrix(c(1, NA_real_, 2, 3), nrow = 2)
  expect_error(
    alfakR:::alfak_neg_log_lik_cpp(c(0.1, 0), counts_bad, c(0, 1)),
    "finite non-negative values"
  )
  expect_error(
    alfakR:::alfak_neg_log_lik_cpp(c(0.1, 0, 1), matrix(c(1, 2, 3, 4), nrow = 2), c(0, 1)),
    "2\\*K - 2"
  )

  expect_error(
    alfakR:::alfak_neighbor_objective_cpp(
      fc_param = 0.1,
      parent_fitness = c(0.2, 0.1),
      pij_values = 0.1,
      parent_birth_times = c(0, 1),
      timepoints = c(0, 1),
      parent_xfit = matrix(0.5, nrow = 2, ncol = 2),
      child_obs = c(0, 1),
      ntot = c(10, 10),
      parent_fitness_mean = 0.15,
      prior_mean = 0,
      prior_sd = 0.1,
      do_prior = FALSE,
      tol = 1e-8
    ),
    "matching lengths/rows"
  )
})

test_that("joint_optimize returns a valid simplex x0 for K = 2 and K > 2", {
  counts_two <- matrix(c(20, 25, 30,
                         80, 75, 70), nrow = 2, byrow = TRUE)
  res_two <- suppressWarnings(alfakR:::joint_optimize(
    counts = counts_two,
    timepoints = c(0, 1, 2),
    f_init = c(0.1, -0.1),
    x0_init = c(0.2, 0.8)
  ))
  expect_true(all(is.finite(res_two$x0)))
  expect_true(all(res_two$x0 >= 0))
  expect_equal(sum(res_two$x0), 1, tolerance = 1e-10)

  counts_three <- matrix(c(30, 28, 25,
                           40, 42, 45,
                           30, 30, 30), nrow = 3, byrow = TRUE)
  res_three <- suppressWarnings(alfakR:::joint_optimize(
    counts = counts_three,
    timepoints = c(0, 1, 2),
    f_init = c(0.1, -0.05, -0.05),
    x0_init = c(0.3, 0.4, 0.3)
  ))
  expect_true(all(is.finite(res_three$x0)))
  expect_true(all(res_three$x0 >= 0))
  expect_equal(sum(res_three$x0), 1, tolerance = 1e-10)
})

test_that("negative parent exposure times are clamped at zero in neighbour objective", {
  full_obj <- alfakR:::alfak_neighbor_objective_cpp(
    fc_param = 0.1,
    parent_fitness = c(0.2, 0.2),
    pij_values = c(0.3, 0.4),
    parent_birth_times = c(10, 0),
    timepoints = c(0, 1),
    parent_xfit = matrix(0.5, nrow = 2, ncol = 2),
    child_obs = c(0, 0),
    ntot = c(10, 10),
    parent_fitness_mean = 0.2,
    prior_mean = 0,
    prior_sd = 0.1,
    do_prior = FALSE,
    tol = 1e-8
  )

  clamped_reference <- alfakR:::alfak_neighbor_objective_cpp(
    fc_param = 0.1,
    parent_fitness = 0.2,
    pij_values = 0.4,
    parent_birth_times = 0,
    timepoints = c(0, 1),
    parent_xfit = matrix(0.5, nrow = 1, ncol = 2),
    child_obs = c(0, 0),
    ntot = c(10, 10),
    parent_fitness_mean = 0.2,
    prior_mean = 0,
    prior_sd = 0.1,
    do_prior = FALSE,
    tol = 1e-8
  )

  expect_equal(full_obj, clamped_reference, tolerance = 1e-12)
})

test_that("passage_times defines one validated internal time axis everywhere", {
  yi <- list(
    x = make_counts(
      c(10, 11, 12,
        20, 21, 22),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1", "2")
    ),
    dt = 99
  )
  passage_times <- c(0, 2.5, 7)
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      res <- alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        passage_times = passage_times
      )

      expect_equal(seen$compute_dx_dt, passage_times)
      expect_equal(seen$optimize_initial_frequencies, passage_times)
      expect_equal(seen$joint_optimize, passage_times)
      expect_equal(seen$project_forward_log, passage_times)
      expect_equal(seen$find_birth_times, c(-1000, max(passage_times)))
      expect_equal(
        unname(res$final_fitness[1, ]),
        rep(log(1e6 / 1e4) / diff(passage_times)[1], 2),
        tolerance = 1e-12
      )
    },
    compute_dx_dt = function(x, timepoints) {
      seen$compute_dx_dt <- timepoints
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      seen$optimize_initial_frequencies <- timepoints
      rep(1 / nrow(x_obs), nrow(x_obs))
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      seen$joint_optimize <- timepoints
      list(f = rep(0, nrow(counts)), x0 = rep(1 / nrow(counts), nrow(counts)))
    },
    project_forward_log = function(x0, f, timepoints) {
      seen$project_forward_log <- timepoints
      matrix(rep(x0, length(timepoints)), nrow = length(x0), ncol = length(timepoints))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      seen$find_birth_times <- time_range
      rep(0, length(opt_res$f))
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = rep(0, nrow(Dmat)))
    },
    gen_nn_info = function(fq, pm) {
      list()
    },
    .package = "alfakR"
  )
})

test_that("resolve_time_axis rejects non-increasing supplied passage_times", {
  yi <- list(
    x = make_counts(
      c(10, 11, 12,
        20, 21, 22),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1", "2")
    ),
    dt = 1
  )

  expect_error(
    alfakR:::resolve_time_axis(yi, passage_times = c(0, 2, 2)),
    "strictly increasing"
  )
})

test_that("birth-time fallback keeps neighbour estimation finite when roots are all missing", {
  yi <- list(
    x = make_counts(
      c(10, 12, 11),
      rownames_vec = "2.2.2",
      colnames_vec = c("0", "1", "2")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  expect_warning(
    res <- testthat::with_mocked_bindings(
      {
        alfakR:::solve_fitness_bootstrap(
          yi,
          minobs = 1,
          nboot = 1,
          n0 = 1e4,
          nb = 1e6,
          pm = 1e-4,
          nn_prior = "none"
        )
      },
      compute_dx_dt = function(x, timepoints) {
        matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
      },
      optimize_initial_frequencies = function(x_obs, f, timepoints) {
        1
      },
      joint_optimize = function(counts, timepoints, f_init, x0_init) {
        list(f = 0, x0 = 1)
      },
      project_forward_log = function(x0, f, timepoints) {
        matrix(1, nrow = 1, ncol = length(timepoints),
               dimnames = list("2.2.2", NULL))
      },
      find_birth_times = function(opt_res, time_range, minF) {
        rep(NA_real_, length(opt_res$f))
      },
      run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
        list(solution = 0)
      },
      gen_nn_info = function(fq, pm) {
        nn <- list(list(ni = "2.2.3", nj = "2.2.2", pij = 0.1))
        names(nn) <- "2.2.3"
        nn
      },
      run_optimise_checked = function(f, interval, ..., context) {
        objective <- f(mean(interval), ...)
        seen$objective <- objective
        list(minimum = mean(interval), objective = objective)
      },
      .package = "alfakR"
    ),
    "finite fallback birth times"
  )

  expect_true(is.finite(seen$objective))
  expect_true(all(is.finite(res$nn_fitness)))
})

test_that("xval defaults to marginal column-wise bootstrap sampling", {
  fq_boot <- list(
    final_fitness = matrix(
      c(1, 10, 100,
        2, 20, 200,
        3, 30, 300),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("1.1", "3.3", "5.5"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 3, ncol = 0)
  )

  testthat::with_mocked_bindings(
    {
      testthat::with_mocked_bindings(
        {
          set.seed(123)
          res <- alfakR:::xval(fq_boot)
          expect_true(is.numeric(res) && length(res) == 1)
        },
        predict = function(object, x, ...) {
          rep(mean(object$train_f), nrow(x))
        },
        .package = "stats"
      )
    },
    Krig = function(x, Y, ...) {
      structure(list(train_f = as.numeric(Y)), class = "mock_krig")
    },
    .package = "fields"
  )
})

test_that("fitKrig uses joint bootstrap mode when requested", {
  fq_boot <- list(
    final_fitness = matrix(
      c(1, 10, 100,
        2, 20, 200,
        3, 30, 300),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("1.1", "3.3", "5.5"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 3, ncol = 0)
  )
  seen <- new.env(parent = emptyenv())
  seen$rows <- list()

  testthat::with_mocked_bindings(
    {
      testthat::with_mocked_bindings(
        {
          set.seed(123)
          suppressWarnings(alfakR:::fitKrig(fq_boot, nboot = 3, krig_bootstrap_mode = "joint"))
          expect_true(length(seen$rows) >= 1)
          expected_rows <- lapply(seq_len(nrow(fq_boot$final_fitness)), function(i) as.numeric(fq_boot$final_fitness[i, ]))
          for (vals in seen$rows) {
            expect_true(any(vapply(expected_rows, function(row_vals) identical(as.numeric(vals), as.numeric(row_vals)), logical(1))))
          }
        },
        predict = function(object, newdata, ...) rep(mean(object$train_f), nrow(newdata)),
        .package = "stats"
      )
    },
    Krig = function(x, Y, ...) {
      seen$rows[[length(seen$rows) + 1L]] <- as.numeric(Y)
      structure(list(train_f = as.numeric(Y)), class = "mock_krig")
    },
    .package = "fields"
  )
})

test_that("cached Krig refits match fresh fields::Krig fits", {
  set.seed(1)
  ktrain <- matrix(runif(16), ncol = 2)
  y_initial <- c(0.2, -0.1, 0.3, 0.05, -0.2, 0.4, 0.1, -0.15)
  y_updated <- c(0.4, -0.25, 0.15, 0.1, -0.1, 0.35, 0.05, -0.05)
  kpred <- matrix(
    c(0.25, 0.5,
      0.75, 0.5,
      0.5, 0.25),
    ncol = 2,
    byrow = TRUE
  )

  cache <- suppressWarnings(
    alfakR:::build_cached_krig_fit(ktrain, y_initial, kpred = kpred, give_warnings = FALSE)
  )
  refit <- suppressWarnings(
    alfakR:::refit_cached_krig(
      cache,
      y = y_updated,
      x_pred = kpred,
      pred_dist = cache$pred_dist,
      give_warnings = FALSE
    )
  )

  fresh_fit <- suppressWarnings(
    fields::Krig(
      ktrain,
      y_updated,
      cov.function = "stationary.cov",
      cov.args = list(Covariance = "Matern", smoothness = 1.5),
      nstep.cv = alfakR:::ALFAK_KRIG_NSTEP_CV,
      give.warnings = FALSE
    )
  )
  fresh_preds <- as.numeric(stats::predict(fresh_fit, kpred))

  expect_equal(refit$fit$lambda, fresh_fit$lambda, tolerance = 1e-10)
  expect_equal(refit$fit$eff.df, fresh_fit$eff.df, tolerance = 1e-10)
  expect_equal(as.numeric(refit$preds), fresh_preds, tolerance = 1e-8)
  expect_equal(as.numeric(stats::predict(refit$fit, kpred)), fresh_preds, tolerance = 1e-8)
})

test_that("fitKrig stops when a bootstrap fit is not trainable", {
  fq_boot <- list(
    final_fitness = matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 3, ncol = 0)
  )

  expect_warning(
    expect_error(
      alfakR:::fitKrig(fq_boot, nboot = 2),
      "Insufficient or incompatible data for Kriging in bootstrap iteration"
    ),
    "Insufficient data for stable"
  )
})

test_that("fitKrig stops when stable or bootstrap Krig fits fail", {
  fq_boot <- list(
    final_fitness = matrix(
      c(0.1, 0.2,
        0.3, 0.4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("2.2", "3.1"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )

  expect_error(
    testthat::with_mocked_bindings(
      {
        alfakR:::fitKrig(fq_boot, nboot = 2)
      },
      Krig = function(...) stop("mock Krig failure"),
      .package = "fields"
    ),
    "mock Krig failure"
  )
})

test_that("default latent-neighbour prior uses empirical_censored and nn_prior = 'none' disables it", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      )
    },
    compute_dx_dt = function(x, timepoints) {
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      1
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      list(f = 0, x0 = 1)
    },
    project_forward_log = function(x0, f, timepoints) {
      matrix(1, nrow = 1, ncol = length(timepoints),
             dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      0
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = 0)
    },
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    estimate_nn_prior_censored_eb = function(...) {
      list(prior_mean = -0.25, prior_sd = 0.33, n_children = 2)
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior)) {
        seen$latent_do_prior <- TRUE
        seen$prior_mean <- prior_mean
        seen$prior_sd <- prior_sd
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_true(isTRUE(seen$latent_do_prior))
  expect_equal(seen$prior_mean, -0.25, tolerance = 1e-12)
  expect_equal(seen$prior_sd, 0.33, tolerance = 1e-12)

  seen_none <- new.env(parent = emptyenv())
  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "none"
      )
    },
    compute_dx_dt = function(x, timepoints) {
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      1
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      list(f = 0, x0 = 1)
    },
    project_forward_log = function(x0, f, timepoints) {
      matrix(1, nrow = 1, ncol = length(timepoints),
             dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      0
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = 0)
    },
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior)) {
        seen_none$latent_do_prior <- TRUE
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_false(isTRUE(seen_none$latent_do_prior))
})

test_that("nn_prior = 'empirical' enables latent-neighbour prior contribution", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical"
      )
    },
    compute_dx_dt = function(x, timepoints) {
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      1
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      list(f = 0, x0 = 1)
    },
    project_forward_log = function(x0, f, timepoints) {
      matrix(1, nrow = 1, ncol = length(timepoints),
             dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      0
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = 0)
    },
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior)) {
        seen$latent_do_prior <- TRUE
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_true(isTRUE(seen$latent_do_prior))
})

test_that("estimate_nn_prior_censored_eb uses all children when fitting the prior", {
  nn_info <- list(
    list(ni = "obs_child", nj = "parent", pij = 1),
    list(ni = "latent_child", nj = "parent", pij = 1)
  )
  names(nn_info) <- c("obs_child", "latent_child")
  fpar <- c(parent = 0)

  prior_fit <- alfakR:::estimate_nn_prior_censored_eb(
    nn_info_items = nn_info,
    fpar = fpar,
    build_opt_fc = function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN, do_prior_param = FALSE) {
      target <- if (identical(nni_param$ni, "obs_child")) 1 else -1
      function(fc_param) (fc_param - target)^2
    },
    search_interval = c(-3, 3),
    nn_prior_sd = 0.4
  )

  expect_equal(prior_fit$n_children, 2)
  expect_equal(prior_fit$prior_sd, 0.4, tolerance = 1e-12)
  expect_lt(abs(prior_fit$prior_mean), 0.25)
})

test_that("estimate_nn_prior_censored_eb skips non-informative children", {
  nn_info <- list(
    list(ni = "informative_child", nj = "parent", pij = 1),
    list(ni = "flat_child", nj = "parent", pij = 1)
  )
  names(nn_info) <- c("informative_child", "flat_child")
  fpar <- c(parent = 0)

  prior_fit <- alfakR:::estimate_nn_prior_censored_eb(
    nn_info_items = nn_info,
    fpar = fpar,
    build_opt_fc = function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN, do_prior_param = FALSE) {
      if (identical(nni_param$ni, "informative_child")) {
        return(function(fc_param) (fc_param - 0.5)^2)
      }
      function(fc_param) 10
    },
    search_interval = c(-3, 3),
    nn_prior_sd = 0.4
  )

  expect_equal(prior_fit$n_children, 1)
  expect_true(is.finite(prior_fit$prior_mean))
  expect_equal(prior_fit$prior_sd, 0.4, tolerance = 1e-12)
})

test_that("estimate_nn_prior_censored_eb errors when no child is informative", {
  nn_info <- list(
    list(ni = "flat_child_a", nj = "parent", pij = 1),
    list(ni = "flat_child_b", nj = "parent", pij = 1)
  )
  names(nn_info) <- c("flat_child_a", "flat_child_b")
  fpar <- c(parent = 0)

  expect_error(
    alfakR:::estimate_nn_prior_censored_eb(
      nn_info_items = nn_info,
      fpar = fpar,
      build_opt_fc = function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN, do_prior_param = FALSE) {
        function(fc_param) 10
      },
      search_interval = c(-3, 3),
      nn_prior_sd = 0.4
    ),
    "no neighbour children produced an informative finite likelihood surface"
  )
})

test_that("nn_prior = 'empirical_censored' enables latent-neighbour prior contribution", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored"
      )
    },
    compute_dx_dt = function(x, timepoints) {
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      1
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      list(f = 0, x0 = 1)
    },
    project_forward_log = function(x0, f, timepoints) {
      matrix(1, nrow = 1, ncol = length(timepoints),
             dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      0
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = 0)
    },
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    estimate_nn_prior_censored_eb = function(...) {
      list(prior_mean = -0.25, prior_sd = 0.33, n_children = 2)
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior)) {
        seen$latent_do_prior <- TRUE
        seen$prior_mean <- prior_mean
        seen$prior_sd <- prior_sd
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_true(isTRUE(seen$latent_do_prior))
  expect_equal(seen$prior_mean, -0.25, tolerance = 1e-12)
  expect_equal(seen$prior_sd, 0.33, tolerance = 1e-12)
})

test_that("nn_prior_grid_n is forwarded to empirical_censored prior fitting", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored",
        nn_prior_grid_n = 57L
      )
    },
    compute_dx_dt = function(x, timepoints) {
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      1
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      list(f = 0, x0 = 1)
    },
    project_forward_log = function(x0, f, timepoints) {
      matrix(1, nrow = 1, ncol = length(timepoints),
             dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      0
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = 0)
    },
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    estimate_nn_prior_censored_eb = function(..., nn_prior_grid_n) {
      seen$nn_prior_grid_n <- nn_prior_grid_n
      list(prior_mean = -0.25, prior_sd = 0.33, n_children = 2)
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_equal(seen$nn_prior_grid_n, 57L, tolerance = 0)
})

test_that("nn_prior = 'empirical_censored' surfaces prior-fit failures without fallback", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  expect_error(
    testthat::with_mocked_bindings(
      {
        alfakR:::solve_fitness_bootstrap(
          yi,
          minobs = 20,
          nboot = 1,
          n0 = 1e4,
          nb = 1e6,
          pm = 1e-4,
          nn_prior = "empirical_censored"
        )
      },
      compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
      optimize_initial_frequencies = function(x_obs, f, timepoints) 1,
      joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0, x0 = 1),
      project_forward_log = function(x0, f, timepoints) matrix(1, nrow = 1, ncol = length(timepoints), dimnames = list("2.2.2", NULL)),
      find_birth_times = function(opt_res, time_range, minF) 0,
      run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = 0),
      gen_nn_info = function(fq, pm) {
        nn <- list(
          list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
          list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
        )
        names(nn) <- c("2.2.3", "2.2.1")
        nn
      },
      estimate_nn_prior_censored_eb = function(...) {
        stop("mock censored prior failure")
      },
      alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                              parent_birth_times, timepoints, parent_xfit,
                                              child_obs, ntot, parent_fitness_mean,
                                              prior_mean, prior_sd, do_prior, tol) {
        if (isTRUE(do_prior)) {
          stop("latent optimisation should not run after prior-fit failure")
        }
        0
      },
      run_optimise_checked = function(f, interval, ..., context) {
        f(mean(interval))
        list(minimum = mean(interval), objective = 0)
      },
      run_optimise_strict_checked = function(f, interval, ..., context) {
        f(mean(interval))
        list(minimum = mean(interval), objective = 0)
      },
      .package = "alfakR"
    ),
    "mock censored prior failure"
  )
})

test_that("empirical prior SD uses floor and user-supplied nn_prior_sd is respected", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  floor_capture <- new.env(parent = emptyenv())
  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical",
        nn_prior_sd_floor = 0.123
      )
    },
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    optimize_initial_frequencies = function(x_obs, f, timepoints) 1,
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0, x0 = 1),
    project_forward_log = function(x0, f, timepoints) matrix(1, nrow = 1, ncol = length(timepoints), dimnames = list("2.2.2", NULL)),
    find_birth_times = function(opt_res, time_range, minF) 0,
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = 0),
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior) && is.finite(prior_sd) && is.null(floor_capture$prior_sd)) {
        floor_capture$prior_sd <- prior_sd
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )
  expect_equal(floor_capture$prior_sd, 0.123, tolerance = 1e-12)

  user_capture <- new.env(parent = emptyenv())
  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical",
        nn_prior_sd = 0.456,
        nn_prior_sd_floor = 0.123
      )
    },
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    optimize_initial_frequencies = function(x_obs, f, timepoints) 1,
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0, x0 = 1),
    project_forward_log = function(x0, f, timepoints) matrix(1, nrow = 1, ncol = length(timepoints), dimnames = list("2.2.2", NULL)),
    find_birth_times = function(opt_res, time_range, minF) 0,
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = 0),
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior) && is.finite(prior_sd) && is.null(user_capture$prior_sd)) {
        user_capture$prior_sd <- prior_sd
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )
  expect_equal(user_capture$prior_sd, 0.456, tolerance = 1e-12)
})

test_that("weighted nearest-neighbour prior mode and controls validate", {
  expect_identical(
    alfakR:::validate_nn_prior_mode("empirical_censored_weighted"),
    "empirical_censored_weighted"
  )
  expect_identical(
    alfakR:::validate_nn_prior_fit_subset("hybrid"),
    "hybrid"
  )
  expect_identical(
    alfakR:::validate_nn_prior_two_step_support("rescue"),
    "rescue"
  )
  expect_silent(
    alfakR:::validate_nn_prior_controls(
      nn_prior_sd = NULL,
      nn_prior_sd_floor = 0.1,
      nn_prior_grid_n = 9L,
      nn_prior_fit_subset = "all",
      nn_prior_zero_exposure_min = 0,
      nn_prior_zero_exposure_quantile = 0.2,
      nn_prior_zero_weight_scale = 0.5,
      nn_prior_zero_weight_cap_ratio = 0.75,
      nn_prior_zero_birth_fallback_weight = 0.25,
      nn_prior_zero_birth_child_floor = 0.2,
      nn_prior_zero_birth_child_shape = 1.5,
      nn_prior_zero_birth_replicate_floor = 0.4,
      nn_prior_zero_birth_replicate_shape = 2,
      nn_prior_hybrid_min_obs = 2L,
      nn_prior_two_step_support = "rescue",
      nn_prior_two_step_support_min = 0.2,
      nn_prior_two_step_cap_floor = 0.4
    )
  )
})

test_that("censored EB fitter is unchanged when explicit child weights are all one", {
  nn_info <- list(
    list(ni = "obs_child", nj = "parent", pij = 1),
    list(ni = "latent_child", nj = "parent", pij = 1)
  )
  names(nn_info) <- c("obs_child", "latent_child")
  fpar <- c(parent = 0)

  fit_default <- alfakR:::estimate_nn_prior_censored_eb(
    nn_info_items = nn_info,
    fpar = fpar,
    build_opt_fc = function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN, do_prior_param = FALSE) {
      target <- if (identical(nni_param$ni, "obs_child")) 1 else -1
      function(fc_param) (fc_param - target)^2
    },
    search_interval = c(-3, 3),
    nn_prior_sd = 0.4
  )

  fit_weighted <- alfakR:::estimate_nn_prior_censored_eb(
    nn_info_items = nn_info,
    fpar = fpar,
    build_opt_fc = function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN, do_prior_param = FALSE) {
      target <- if (identical(nni_param$ni, "obs_child")) 1 else -1
      function(fc_param) (fc_param - target)^2
    },
    search_interval = c(-3, 3),
    nn_prior_sd = 0.4,
    child_weights = c(1, 1)
  )

  expect_equal(fit_weighted$prior_mean, fit_default$prior_mean, tolerance = 1e-12)
  expect_equal(fit_weighted$prior_sd, fit_default$prior_sd, tolerance = 1e-12)
  expect_equal(fit_weighted$n_children, fit_default$n_children, tolerance = 0)
})

test_that("weighted projected child exposure uses the projected neutral child trajectory", {
  exposure <- alfakR:::project_nn_child_exposure(
    fc_param = 0.5,
    parent_fitness = 0.5,
    pij_values = 0.2,
    parent_birth_times = 0,
    timepoints = c(0, 2),
    parent_xfit = matrix(c(0.8, 0.4), nrow = 1),
    ntot = c(100, 50)
  )

  expect_equal(exposure, 4, tolerance = 1e-12)
})

test_that("NN C++ projection failures are logged and stop", {
  log_path <- tempfile("alfak_run_log_")
  old <- options(alfakR.run_log_path = log_path, alfakR.echo_run_log = FALSE)
  on.exit(options(old), add = TRUE)

  expect_error(
    testthat::with_mocked_bindings(
      alfakR:::project_nn_child_trajectory(
        fc_param = 0.5,
        parent_fitness = 0.5,
        pij_values = 0.2,
        parent_birth_times = 0,
        timepoints = c(0, 2),
        parent_xfit = matrix(c(0.8, 0.4), nrow = 1)
      ),
      alfak_nn_project_trajectory_cpp = function(...) stop("forced trajectory failure"),
      .package = "alfakR"
    ),
    "alfak_nn_project_trajectory_cpp.*forced trajectory failure"
  )
  lines <- alfakR::alfak_read_run_log(path = log_path)
  expect_true(any(grepl("cpp.alfak_nn_project_trajectory_cpp", lines)))
  expect_true(any(grepl("forced trajectory failure", lines)))
})

test_that("weighted parent centering uses exposure opportunity weights and falls back cleanly", {
  parent_fitness <- c(1, 3)
  pij_values <- c(0.5, 0.5)
  parent_birth_times <- c(0, 1)
  timepoints <- c(0, 1, 2)
  parent_xfit <- matrix(
    c(0.9, 0.9, 0.9,
      0.1, 0.1, 0.1),
    nrow = 2,
    byrow = TRUE
  )
  ntot <- c(10, 10, 10)
  fallback_mean <- stats::weighted.mean(parent_fitness, w = pij_values)
  opportunity_weights <- c(
    0.5 * sum(ntot * c(1, 1, 1) * parent_xfit[1, ]),
    0.5 * sum(ntot * c(0, 1, 1) * parent_xfit[2, ])
  )

  expect_equal(
    alfakR:::resolve_nn_parent_opportunity_weights(
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      ntot = ntot
    ),
    opportunity_weights,
    tolerance = 1e-12
  )

  expect_equal(
    alfakR:::weighted_parent_fitness_exposure(
      parent_fitness = parent_fitness,
      parent_opportunity_weights = opportunity_weights,
      fallback_mean = fallback_mean
    ),
    stats::weighted.mean(parent_fitness, w = opportunity_weights),
    tolerance = 1e-12
  )

  expect_equal(
    alfakR:::weighted_parent_fitness_exposure(
      parent_fitness = parent_fitness,
      parent_opportunity_weights = c(0, 0),
      fallback_mean = fallback_mean
    ),
    fallback_mean,
    tolerance = 1e-12
  )
})

test_that("weighted birth fallback burden uses child and replicate burden smoothly", {
  weighted_fit <- alfakR:::prepare_weighted_nn_prior_fit(
    nn_child_contexts = list(
      obs = list(
        ni = "obs",
        projected_exposure = 10,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      ),
      zero_partial = list(
        ni = "zero_partial",
        projected_exposure = 10,
        parent_birth_fallback = c(FALSE, TRUE),
        parent_opportunity_weights = c(3, 1),
        parent_fitness_mean_exposure = 0
      ),
      zero_full = list(
        ni = "zero_full",
        projected_exposure = 10,
        parent_birth_fallback = c(TRUE, TRUE),
        parent_opportunity_weights = c(1, 1),
        parent_fitness_mean_exposure = 0
      )
    ),
    nn_present = c(TRUE, FALSE, FALSE),
    nn_prior_fit_subset = "all",
    nn_prior_zero_weight_scale = 1,
    nn_prior_zero_weight_cap_ratio = 1,
    nn_prior_zero_birth_child_floor = 0.25,
    nn_prior_zero_birth_child_shape = 1,
    nn_prior_zero_birth_replicate_floor = 0.50,
    nn_prior_zero_birth_replicate_shape = 1,
    ntot = c(10, 10)
  )

  expect_equal(
    weighted_fit$child_weights,
    c(1, 0.55859375, 0.171875),
    tolerance = 1e-12
  )
  expect_equal(
    weighted_fit$diagnostics$mean_zero_birth_fallback_burden,
    mean(c(0.25, 1)),
    tolerance = 1e-12
  )
  expect_equal(
    weighted_fit$diagnostics$replicate_birth_fallback_burden,
    0.625,
    tolerance = 1e-12
  )
  expect_equal(
    weighted_fit$diagnostics$replicate_birth_reliability_multiplier,
    0.6875,
    tolerance = 1e-12
  )
  expect_equal(
    weighted_fit$diagnostics$zero_effective_mass_used,
    1.0625,
    tolerance = 1e-12
  )
})

test_that("2-step rescue shares descendant support across competing zero children", {
  support <- alfakR:::compute_nn_two_step_support(
    nn_child_contexts = list(
      zero_a = list(
        ni = "2.2.1",
        nj = "2.2.2",
        projected_exposure = 10
      ),
      zero_b = list(
        ni = "2.3.2",
        nj = "2.2.2",
        projected_exposure = 10
      )
    ),
    zero_mask = c(TRUE, TRUE),
    count_data = make_counts(
      c(20, 20,
        2, 2),
      rownames_vec = c("2.2.2", "2.3.1"),
      colnames_vec = c("0", "1")
    ),
    pm = 1e-4,
    exposure_reference = 10,
    child_birth_multiplier = c(1, 1)
  )

  expect_true(all(support$child_support > 0.39 & support$child_support < 0.394))
  expect_lt(diff(range(support$child_support)), 1e-04)
  expect_equal(support$n_children_with_support, 2L, tolerance = 0)
  expect_equal(support$descendant_exposure_reference, 4, tolerance = 1e-12)
})

test_that("weighted hybrid rescue retains low-exposure zero children with observed 2-step support", {
  weighted_fit_none <- alfakR:::prepare_weighted_nn_prior_fit(
    nn_child_contexts = list(
      obs = list(
        ni = "2.2.3",
        nj = "2.2.2",
        projected_exposure = 20,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      ),
      zero = list(
        ni = "2.2.1",
        nj = "2.2.2",
        projected_exposure = 1,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      )
    ),
    nn_present = c(TRUE, FALSE),
    nn_prior_fit_subset = "hybrid",
    nn_prior_zero_exposure_min = 10,
    nn_prior_zero_weight_scale = 1,
    nn_prior_zero_weight_cap_ratio = 1,
    nn_prior_two_step_support = "none",
    count_data = make_counts(
      c(40, 40,
        5, 5,
        2, 2),
      rownames_vec = c("2.2.2", "2.2.3", "2.3.1"),
      colnames_vec = c("0", "1")
    ),
    pm = 1e-4,
    ntot = c(10, 10)
  )

  weighted_fit_rescue <- alfakR:::prepare_weighted_nn_prior_fit(
    nn_child_contexts = list(
      obs = list(
        ni = "2.2.3",
        nj = "2.2.2",
        projected_exposure = 20,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      ),
      zero = list(
        ni = "2.2.1",
        nj = "2.2.2",
        projected_exposure = 1,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      )
    ),
    nn_present = c(TRUE, FALSE),
    nn_prior_fit_subset = "hybrid",
    nn_prior_zero_exposure_min = 10,
    nn_prior_zero_weight_scale = 1,
    nn_prior_zero_weight_cap_ratio = 1,
    nn_prior_two_step_support = "rescue",
    nn_prior_two_step_support_min = 0.15,
    nn_prior_two_step_cap_floor = 0.30,
    count_data = make_counts(
      c(40, 40,
        5, 5,
        2, 2),
      rownames_vec = c("2.2.2", "2.2.3", "2.3.1"),
      colnames_vec = c("0", "1")
    ),
    pm = 1e-4,
    ntot = c(10, 10)
  )

  expect_equal(weighted_fit_none$diagnostics$n_zero_children_retained, 0L, tolerance = 0)
  expect_equal(weighted_fit_rescue$diagnostics$n_zero_children_retained, 1L, tolerance = 0)
  expect_equal(
    weighted_fit_rescue$diagnostics$max_zero_two_step_support,
    1 - exp(-1),
    tolerance = 1e-12
  )
  expect_equal(
    weighted_fit_rescue$diagnostics$zero_effective_mass_used,
    0.30 * (1 - exp(-1)),
    tolerance = 1e-12
  )
  expect_equal(
    weighted_fit_rescue$child_weights,
    c(1, 0.30 * (1 - exp(-1))),
    tolerance = 1e-12
  )
})

test_that("adaptive zero cap uses effective zero mass rather than raw zero counts", {
  zero_items <- lapply(seq_len(100), function(i) {
    list(
      ni = paste0("zero_", i),
      projected_exposure = 1,
      parent_birth_fallback = FALSE,
      parent_opportunity_weights = 1,
      parent_fitness_mean_exposure = 0
    )
  })
  names(zero_items) <- vapply(zero_items, function(item) item$ni, character(1))
  nn_child_contexts <- c(
    list(obs = list(
      ni = "obs",
      projected_exposure = 100,
      parent_birth_fallback = FALSE,
      parent_opportunity_weights = 1,
      parent_fitness_mean_exposure = 0
    )),
    zero_items
  )

  weighted_fit <- alfakR:::prepare_weighted_nn_prior_fit(
    nn_child_contexts = nn_child_contexts,
    nn_present = c(TRUE, rep(FALSE, length(zero_items))),
    nn_prior_fit_subset = "all",
    nn_prior_zero_weight_scale = 1,
    nn_prior_zero_weight_cap_ratio = NULL,
    ntot = c(10, 10)
  )

  expect_false(weighted_fit$diagnostics$zero_weight_cap_applied)
  expect_equal(weighted_fit$diagnostics$zero_effective_mass_used, 1, tolerance = 1e-12)
  expect_equal(weighted_fit$diagnostics$zero_weight_cap_ratio_used, 1, tolerance = 1e-12)
  expect_equal(weighted_fit$diagnostics$sum_zero_weight_final, 1, tolerance = 1e-12)
  expect_equal(weighted_fit$child_weights[1], 1, tolerance = 1e-12)
  expect_equal(weighted_fit$child_weights[-1], rep(0.01, length(zero_items)), tolerance = 1e-12)
})

test_that("adaptive zero cap still shrinks a small number of strong zero children", {
  weighted_fit <- alfakR:::prepare_weighted_nn_prior_fit(
    nn_child_contexts = list(
      obs = list(
        ni = "obs",
        projected_exposure = 100,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      ),
      zero_a = list(
        ni = "zero_a",
        projected_exposure = 100,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      ),
      zero_b = list(
        ni = "zero_b",
        projected_exposure = 100,
        parent_birth_fallback = FALSE,
        parent_opportunity_weights = 1,
        parent_fitness_mean_exposure = 0
      )
    ),
    nn_present = c(TRUE, FALSE, FALSE),
    nn_prior_fit_subset = "all",
    nn_prior_zero_weight_scale = 1,
    nn_prior_zero_weight_cap_ratio = NULL,
    ntot = c(10, 10)
  )

  expect_true(weighted_fit$diagnostics$zero_weight_cap_applied)
  expect_equal(weighted_fit$diagnostics$zero_effective_mass_used, 2, tolerance = 1e-12)
  expect_equal(weighted_fit$diagnostics$zero_weight_cap_ratio_used, sqrt(1 / 2), tolerance = 1e-12)
  expect_equal(weighted_fit$diagnostics$sum_zero_weight_pre_cap, 2, tolerance = 1e-12)
  expect_equal(weighted_fit$diagnostics$sum_zero_weight_post_cap, sqrt(1 / 2), tolerance = 1e-12)
  expect_equal(weighted_fit$child_weights, c(1, sqrt(1 / 8), sqrt(1 / 8)), tolerance = 1e-12)
})

test_that("weighted hybrid prior screens low-exposure zeros and keeps observed children at unit weight", {
  yi <- list(
    x = make_counts(
      c(40, 40,
        5, 5),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  res <- testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored_weighted",
        nn_prior_fit_subset = "hybrid",
        nn_prior_zero_exposure_min = 10,
        nn_prior_zero_weight_scale = 0.5,
        nn_prior_zero_weight_cap_ratio = 1
      )
    },
    bootstrap_counts = function(x) x,
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = rep(0.5, nrow(Dmat))),
    optimize_initial_frequencies = function(x_obs, f, timepoints) rep(1 / length(f), length(f)),
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0.5, x0 = 1),
    project_forward_log = function(x0, f, timepoints) {
      matrix(c(0.8, 0.4), nrow = 1, dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) 0,
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.5),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.5),
        list(ni = "2.3.2", nj = "2.2.2", pij = 0.05)
      )
      names(nn) <- c("2.2.3", "2.2.1", "2.3.2")
      nn
    },
    estimate_nn_prior_censored_eb = function(nn_info_items, ..., child_weights, parent_mean_fn) {
      seen$prior_children <- names(nn_info_items)
      seen$child_weights <- child_weights
      list(
        prior_mean = 0.1,
        prior_sd = 0.2,
        informative_child_count = length(child_weights),
        map_delta_lower_boundary_rate = 0,
        map_delta_upper_boundary_rate = 0
      )
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_identical(seen$prior_children, c("2.2.3", "2.2.1"))
  expect_equal(seen$child_weights, c(1, 0.5), tolerance = 1e-12)

  diag <- res$nn_prior_diagnostics[1, ]
  required_diag_cols <- c(
    "nn_prior_mode_used",
    "nn_prior_source_used",
    "nn_prior_fit_subset_used",
    "n_observed_children",
    "n_zero_children_total",
    "n_zero_children_retained",
    "n_zero_children_screened",
    "sum_observed_weight",
    "sum_zero_weight_raw",
    "sum_zero_weight_final",
    "zero_weight_cap_applied",
    "zero_weight_cap_ratio_used",
    "zero_effective_mass_used",
    "zero_effective_mass_mean",
    "zero_effective_mass_median",
    "sum_zero_weight_pre_cap",
    "sum_zero_weight_post_cap",
    "exposure_threshold_used",
    "exposure_reference_used",
    "nn_prior_two_step_support_used",
    "two_step_support_min_used",
    "two_step_cap_floor_used",
    "two_step_descendant_exposure_reference_used",
    "n_zero_children_with_two_step_support",
    "mean_zero_two_step_support",
    "median_zero_two_step_support",
    "max_zero_two_step_support",
    "n_zero_children_with_birth_fallback",
    "mean_zero_birth_fallback_burden",
    "median_zero_birth_fallback_burden",
    "replicate_birth_fallback_burden",
    "mean_zero_birth_reliability_multiplier",
    "median_zero_birth_reliability_multiplier",
    "replicate_birth_reliability_multiplier",
    "sample_pooled_prior_available",
    "sample_pooled_prior_mu",
    "sample_pooled_prior_sigma",
    "sample_pooled_prior_informative_child_count",
    "sample_pooled_alpha_used",
    "sample_pooled_sigma_used",
    "prior_mu_hat",
    "prior_sigma_hat",
    "informative_child_count",
    "map_delta_lower_boundary_rate",
    "map_delta_upper_boundary_rate",
    "used_sample_pooled_fallback_for_this_replicate",
    "used_no_prior_fallback_for_this_replicate"
  )
  expect_true(all(required_diag_cols %in% colnames(res$nn_prior_diagnostics)))
  expect_identical(diag$nn_prior_mode_used, "empirical_censored_weighted")
  expect_identical(diag$nn_prior_source_used, "observed_replicate")
  expect_identical(diag$nn_prior_fit_subset_used, "hybrid")
  expect_equal(diag$n_observed_children, 1L, tolerance = 0)
  expect_equal(diag$n_zero_children_total, 2L, tolerance = 0)
  expect_equal(diag$n_zero_children_retained, 1L, tolerance = 0)
  expect_equal(diag$n_zero_children_screened, 1L, tolerance = 0)
  expect_true(is.finite(diag$sum_observed_weight))
  expect_true(is.finite(diag$sum_zero_weight_raw))
  expect_true(is.finite(diag$sum_zero_weight_final))
  expect_true(is.finite(diag$zero_weight_cap_ratio_used))
  expect_true(is.finite(diag$zero_effective_mass_used))
  expect_true(is.finite(diag$zero_effective_mass_mean))
  expect_true(is.finite(diag$zero_effective_mass_median))
  expect_true(is.finite(diag$sum_zero_weight_pre_cap))
  expect_true(is.finite(diag$sum_zero_weight_post_cap))
  expect_true(is.finite(diag$exposure_threshold_used))
  expect_true(is.finite(diag$exposure_reference_used))
  expect_true(is.character(diag$nn_prior_two_step_support_used))
  expect_true(is.finite(diag$two_step_support_min_used))
  expect_true(is.finite(diag$two_step_cap_floor_used))
  expect_true(is.finite(diag$two_step_descendant_exposure_reference_used))
  expect_true(is.finite(diag$n_zero_children_with_two_step_support))
  expect_true(is.finite(diag$mean_zero_two_step_support))
  expect_true(is.finite(diag$median_zero_two_step_support))
  expect_true(is.finite(diag$max_zero_two_step_support))
  expect_true(is.finite(diag$mean_zero_birth_fallback_burden))
  expect_true(is.finite(diag$median_zero_birth_fallback_burden))
  expect_true(is.finite(diag$replicate_birth_fallback_burden))
  expect_true(is.finite(diag$mean_zero_birth_reliability_multiplier))
  expect_true(is.finite(diag$median_zero_birth_reliability_multiplier))
  expect_true(is.finite(diag$replicate_birth_reliability_multiplier))
  expect_true(isTRUE(diag$sample_pooled_prior_available))
  expect_true(is.finite(diag$sample_pooled_prior_mu))
  expect_true(is.finite(diag$sample_pooled_prior_sigma))
  expect_equal(diag$sample_pooled_prior_informative_child_count, 1L, tolerance = 0)
  expect_false(isTRUE(diag$used_sample_pooled_fallback_for_this_replicate))
  expect_false(isTRUE(diag$used_no_prior_fallback_for_this_replicate))
})

test_that("weighted prior downweights zero children when birth times were filled by fallback", {
  yi <- list(
    x = make_counts(
      c(40, 40,
        5, 5),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  res <- testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored_weighted",
        nn_prior_zero_weight_scale = 1,
        nn_prior_zero_weight_cap_ratio = 1,
        nn_prior_zero_birth_child_floor = 0.25,
        nn_prior_zero_birth_replicate_floor = 0.50
      )
    },
    bootstrap_counts = function(x) x,
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = rep(0.5, nrow(Dmat))),
    optimize_initial_frequencies = function(x_obs, f, timepoints) rep(1 / length(f), length(f)),
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0.5, x0 = 1),
    project_forward_log = function(x0, f, timepoints) {
      matrix(c(0.8, 0.4), nrow = 1, dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) NA_real_,
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.5),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.5)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    estimate_nn_prior_censored_eb = function(nn_info_items, ..., child_weights, parent_mean_fn) {
      seen$child_weights <- child_weights
      list(
        prior_mean = 0.1,
        prior_sd = 0.2,
        informative_child_count = length(child_weights),
        map_delta_lower_boundary_rate = 0,
        map_delta_upper_boundary_rate = 0
      )
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_equal(seen$child_weights, c(1, 0.125), tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$n_zero_children_with_birth_fallback[1], 1L, tolerance = 0)
  expect_equal(res$nn_prior_diagnostics$replicate_birth_fallback_burden[1], 1, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$replicate_birth_reliability_multiplier[1], 0.5, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$zero_effective_mass_used[1], 0.25, tolerance = 1e-12)
})

test_that("weighted prior applies the zero-weight cap by common rescaling", {
  yi <- list(
    x = make_counts(
      c(40, 40,
        5, 5),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  res <- testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored_weighted",
        nn_prior_zero_weight_scale = 0.5,
        nn_prior_zero_weight_cap_ratio = 0.25
      )
    },
    bootstrap_counts = function(x) x,
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = rep(0.5, nrow(Dmat))),
    optimize_initial_frequencies = function(x_obs, f, timepoints) rep(1 / length(f), length(f)),
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0.5, x0 = 1),
    project_forward_log = function(x0, f, timepoints) {
      matrix(c(0.8, 0.4), nrow = 1, dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) 0,
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.5),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.5),
        list(ni = "2.3.2", nj = "2.2.2", pij = 0.5)
      )
      names(nn) <- c("2.2.3", "2.2.1", "2.3.2")
      nn
    },
    estimate_nn_prior_censored_eb = function(nn_info_items, ..., child_weights, parent_mean_fn) {
      seen$child_weights <- child_weights
      list(
        prior_mean = 0.1,
        prior_sd = 0.2,
        informative_child_count = length(child_weights),
        map_delta_lower_boundary_rate = 0,
        map_delta_upper_boundary_rate = 0
      )
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_equal(seen$child_weights, c(1, 0.125, 0.125), tolerance = 1e-12)
  expect_true(res$nn_prior_diagnostics$zero_weight_cap_applied[1])
  expect_equal(res$nn_prior_diagnostics$zero_weight_cap_ratio_used[1], 0.25, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$sum_zero_weight_pre_cap[1], 1, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$sum_zero_weight_post_cap[1], 0.25, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$sum_zero_weight_final[1], 0.25, tolerance = 1e-12)
})

test_that("weighted prior uses a sample-pooled fallback when a bootstrap replicate has no observed neighbour children", {
  yi <- list(
    x = make_counts(
      c(40, 40,
        5, 5),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  res <- testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored_weighted"
      )
    },
    bootstrap_counts = function(x) {
      y <- x
      y["2.2.3", ] <- 0
      y
    },
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = rep(0.5, nrow(Dmat))),
    optimize_initial_frequencies = function(x_obs, f, timepoints) rep(1 / length(f), length(f)),
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0.5, x0 = 1),
    project_forward_log = function(x0, f, timepoints) {
      matrix(c(0.8, 0.4), nrow = 1, dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) 0,
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.5),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.5)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    estimate_nn_prior_censored_eb = function(..., context) {
      seen$estimate_contexts <- c(seen$estimate_contexts, context)
      if (grepl("sample-pooled", context)) {
        return(list(
          prior_mean = 0.1,
          prior_sd = 0.2,
          informative_child_count = 1L,
          sum_child_weight = 1,
          map_delta_lower_boundary_rate = 0,
          map_delta_upper_boundary_rate = 0
        ))
      }
      stop("weighted replicate prior fit should not run without observed neighbour children")
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior)) {
        seen$latent_do_prior <- TRUE
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      seen$strict_contexts <- c(seen$strict_contexts, context)
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_true(isTRUE(seen$latent_do_prior))
  expect_length(seen$estimate_contexts, 1)
  expect_match(seen$estimate_contexts[[1]], "sample-pooled")
  expect_true(any(grepl("sample-pooled empirical_censored_weighted prior", seen$strict_contexts)))
  expect_identical(res$nn_prior_diagnostics$nn_prior_mode_used[1], "empirical_censored_weighted")
  expect_identical(res$nn_prior_diagnostics$nn_prior_source_used[1], "sample_pooled")
  expect_true(res$nn_prior_diagnostics$sample_pooled_prior_available[1])
  expect_equal(res$nn_prior_diagnostics$sample_pooled_prior_mu[1], 0.1, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$sample_pooled_prior_sigma[1], 0.2, tolerance = 1e-12)
  expect_equal(res$nn_prior_diagnostics$prior_mu_hat[1], 0.1, tolerance = 1e-12)
  expect_true(is.finite(res$nn_prior_diagnostics$sample_pooled_alpha_used[1]))
  expect_true(res$nn_prior_diagnostics$sample_pooled_alpha_used[1] > 0)
  expect_true(res$nn_prior_diagnostics$sample_pooled_sigma_used[1] >= res$nn_prior_diagnostics$sample_pooled_prior_sigma[1])
  expect_equal(
    res$nn_prior_diagnostics$prior_sigma_hat[1],
    res$nn_prior_diagnostics$sample_pooled_sigma_used[1],
    tolerance = 1e-12
  )
  expect_true(res$nn_prior_diagnostics$used_sample_pooled_fallback_for_this_replicate[1])
  expect_false(res$nn_prior_diagnostics$used_no_prior_fallback_for_this_replicate[1])
})

test_that("weighted prior falls back to no prior when sample-pooled fallback is unavailable", {
  yi <- list(
    x = make_counts(
      c(40, 40),
      rownames_vec = "2.2.2",
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  res <- testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "empirical_censored_weighted"
      )
    },
    bootstrap_counts = function(x) x,
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = rep(0.5, nrow(Dmat))),
    optimize_initial_frequencies = function(x_obs, f, timepoints) rep(1 / length(f), length(f)),
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0.5, x0 = 1),
    project_forward_log = function(x0, f, timepoints) {
      matrix(c(0.8, 0.4), nrow = 1, dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) 0,
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.5),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.5)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    estimate_nn_prior_censored_eb = function(...) {
      stop("sample-pooled or replicate weighted prior fit should not run without any observed neighbour children")
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior)) {
        seen$latent_do_prior <- TRUE
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    run_optimise_strict_checked = function(f, interval, ..., context) {
      stop("weighted no-prior fallback should not use strict prior optimisation")
    },
    .package = "alfakR"
  )

  expect_false(isTRUE(seen$latent_do_prior))
  expect_identical(res$nn_prior_diagnostics$nn_prior_mode_used[1], "none")
  expect_identical(res$nn_prior_diagnostics$nn_prior_source_used[1], "none")
  expect_false(res$nn_prior_diagnostics$sample_pooled_prior_available[1])
  expect_false(res$nn_prior_diagnostics$used_sample_pooled_fallback_for_this_replicate[1])
  expect_true(res$nn_prior_diagnostics$used_no_prior_fallback_for_this_replicate[1])
})

test_that("fitKrig and xval ignore bootstrap nearest-neighbour diagnostics", {
  fq_boot <- list(
    final_fitness = matrix(
      c(1, 10, 100,
        2, 20, 200,
        3, 30, 300),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("1.1", "3.3", "5.5"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 3, ncol = 0),
    nn_prior_diagnostics = data.frame(
      nn_prior_mode_used = c("empirical_censored_weighted", "none", "none"),
      stringsAsFactors = FALSE
    )
  )

  testthat::with_mocked_bindings(
    {
      testthat::with_mocked_bindings(
        {
          expect_silent(suppressWarnings(alfakR:::fitKrig(fq_boot, nboot = 3)))
          set.seed(123)
          res <- alfakR:::xval(fq_boot)
          expect_true(is.numeric(res) && length(res) == 1)
        },
        predict = function(object, x, ...) {
          rep(mean(object$train_f), nrow(x))
        },
        .package = "stats"
      )
    },
    Krig = function(x, Y, ...) {
      structure(list(train_f = as.numeric(Y)), class = "mock_krig")
    },
    .package = "fields"
  )
})

test_that("find_steady_state selects the eigenvalue with the largest real part", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("Matrix")

  lscape <- data.frame(
    k = c("2.2", "3.1", "1.3"),
    mean = c(-3, 2, 1),
    stringsAsFactors = FALSE
  )
  ode_A <- Matrix::Matrix(
    diag(c(-3, 2, 1)),
    sparse = TRUE,
    dimnames = list(lscape$k, lscape$k)
  )
  ode_out <- deSolve::ode(
    y = c(0.2, 0.3, 0.5),
    times = seq(0, 100, by = 0.1),
    func = alfakR:::chrmod_rel,
    parms = list(A = ode_A)
  )
  terminal <- as.numeric(ode_out[nrow(ode_out), -1])
  terminal <- terminal / sum(terminal)

  ss <- testthat::with_mocked_bindings(
    {
      alfakR::find_steady_state(lscape, p = 0.01)
    },
    build_W_rcpp = function(karyotype_strings, p, Nmax = Inf) {
      W <- Matrix::Diagonal(length(karyotype_strings), x = rep(1, length(karyotype_strings)))
      dimnames(W) <- list(karyotype_strings, karyotype_strings)
      W
    },
    .package = "alfakR"
  )

  expect_equal(unname(ss), terminal, tolerance = 1e-4)
})

test_that("nearest-neighbour exposure uses projected frequent-parent frequencies directly", {
  yi <- list(
    x = make_counts(
      c(10, 10,
        1, 1,
        0, 0),
      rownames_vec = c("2.2.2", "9.9.9", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 10,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "none"
      )
    },
    bootstrap_counts = function(x) x,
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = 0),
    optimize_initial_frequencies = function(x_obs, f, timepoints) 1,
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0, x0 = 1),
    project_forward_log = function(x0, f, timepoints) matrix(c(0.8, 0.6), nrow = 1, dimnames = list("2.2.2", NULL)),
    find_birth_times = function(opt_res, time_range, minF) 0,
    gen_nn_info = function(fq, pm) {
      nn <- list(list(ni = "2.2.3", nj = "2.2.2", pij = 0.1))
      names(nn) <- "2.2.3"
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      seen$parent_xfit <- parent_xfit
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_equal(as.numeric(seen$parent_xfit), c(0.8, 0.6), tolerance = 1e-12)
})

test_that("ABM parent accounting consumes each dividing parent exactly once", {
  res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(10), "1"),
    fitness_map_r = stats::setNames(list(1, 1), c("1", "2")),
    p_missegregation = 1,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )

  expect_true("1" %in% names(res))
  step1_counts <- res[["1"]]
  expect_equal(sum(as.numeric(step1_counts)), 10)
})

test_that("landscape_data_output controls whether landscape_data.Rds is written", {
  yi <- list(
    x = make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(
      k = "2.2.2",
      mean = 0,
      median = 0,
      sd = 0,
      fq = TRUE,
      nn = FALSE
    ),
    posterior_samples = matrix(0, nrow = 1, ncol = 1),
    krig_stable_mean = list(tag = "mean"),
    krig_stable_median = list(tag = "median")
  )
  xval_stub <- 0
  outdir_false <- file.path(tempdir(), "alfak_landscape_data_false")
  outdir_true <- file.path(tempdir(), "alfak_landscape_data_true")
  unlink(outdir_false, recursive = TRUE)
  unlink(outdir_true, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir_false,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        landscape_data_output = FALSE
      ))
      expect_false(file.exists(file.path(outdir_false, "landscape_data.Rds")))

      invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir_true,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        landscape_data_output = TRUE
      ))
      expect_true(file.exists(file.path(outdir_true, "landscape_data.Rds")))
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    xval = function(...) xval_stub,
    .package = "alfakR"
  )
})

test_that("softmax is stable for large logits and rejects non-finite values", {
  expect_equal(alfakR:::softmax(c(1000, 1000)), c(0.5, 0.5), tolerance = 1e-12)
  expect_error(alfakR:::softmax(c(0, Inf)), "non-finite logits")
})

test_that("alfak saves xval.Rds as a scalar R2R for downstream compatibility", {
  yi <- list(
    x = make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(
      k = "2.2.2",
      mean = 0,
      median = 0,
      sd = 0,
      fq = TRUE,
      nn = FALSE
    ),
    posterior_samples = matrix(0, nrow = 1, ncol = 1),
    krig_stable_mean = NULL,
    krig_stable_median = NULL
  )
  xval_stub <- 0.42
  outdir <- file.path(tempdir(), "alfak_xval_scalar")
  unlink(outdir, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      ))
      saved <- readRDS(file.path(outdir, "xval.Rds"))
      expect_identical(returned, 0.42)
      expect_identical(saved, 0.42)
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    xval = function(...) xval_stub,
    .package = "alfakR"
  )
})

test_that("alfak accepts scalar NA_real_ cross-validation outputs and still writes core files", {
  yi <- list(
    x = make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(
      k = "2.2.2",
      mean = 0,
      median = 0,
      sd = 0,
      fq = TRUE,
      nn = FALSE
    ),
    posterior_samples = matrix(0, nrow = 1, ncol = 1),
    krig_stable_mean = NULL,
    krig_stable_median = NULL
  )
  outdir <- file.path(tempdir(), "alfak_xval_na")
  unlink(outdir, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      ))
      saved <- readRDS(file.path(outdir, "xval.Rds"))
      expect_true(is.numeric(returned) && length(returned) == 1 && is.na(returned))
      expect_true(is.numeric(saved) && length(saved) == 1 && is.na(saved))
      expect_true(file.exists(file.path(outdir, "bootstrap_res.Rds")))
      expect_true(file.exists(file.path(outdir, "landscape.Rds")))
      expect_true(file.exists(file.path(outdir, "landscape_posterior_samples.Rds")))
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    xval = function(...) NA_real_,
    .package = "alfakR"
  )
})

test_that("constant-response cross-validation now follows upstream NaN semantics", {
  fq_boot <- list(
    final_fitness = matrix(
      c(1, 1, 1,
        1, 1, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("2.2", "3.3", "4.4"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )

  res <- testthat::with_mocked_bindings(
    {
      testthat::with_mocked_bindings(
        {
          alfakR:::xval(fq_boot)
        },
        predict = function(object, x, ...) {
          rep(mean(object$train_f), nrow(x))
        },
        .package = "stats"
      )
    },
    Krig = function(x, Y, ...) {
      structure(list(train_f = as.numeric(Y)), class = "mock_krig")
    },
    .package = "fields"
  )
  expect_true(is.numeric(res) && length(res) == 1 && is.nan(res))
})

test_that("alfak stops when Krig fitting fails during xval", {
  yi <- make_simple_yi(
    make_counts(
      c(10, 12,
        20, 18,
        5, 6),
      rownames_vec = c("2.2.2", "2.2.1", "2.2.3"),
      colnames_vec = c("0", "1")
    )
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(c(0.1, 0.2, 0.3,
                               0.15, 0.25, 0.35), nrow = 2, byrow = TRUE,
                             dimnames = list(NULL, c("2.2.2", "5.5.5", "8.8.8"))),
    final_fitness = matrix(c(0.1, 0.2, 0.3,
                             0.15, 0.25, 0.35), nrow = 2, byrow = TRUE,
                           dimnames = list(NULL, c("2.2.2", "5.5.5", "8.8.8"))),
    initial_frequencies = matrix(c(0.3, 0.4, 0.3,
                                   0.25, 0.45, 0.3), nrow = 2, byrow = TRUE,
                                 dimnames = list(NULL, c("2.2.2", "5.5.5", "8.8.8"))),
    final_frequencies = matrix(c(0.3, 0.4, 0.3,
                                 0.25, 0.45, 0.3), nrow = 2, byrow = TRUE,
                               dimnames = list(NULL, c("2.2.2", "5.5.5", "8.8.8"))),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(k = c("2.2.2", "5.5.5", "8.8.8"), mean = 0, median = 0, sd = 0, fq = TRUE, nn = FALSE),
    posterior_samples = matrix(0, nrow = 3, ncol = 1),
    krig_stable_mean = NULL,
    krig_stable_median = NULL
  )
  outdir <- file.path(tempdir(), "alfak_xval_krig_failure")
  unlink(outdir, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      testthat::with_mocked_bindings(
        {
          expect_error(
            {
              set.seed(1)
              alfakR::alfak(
                yi = yi,
                outdir = outdir,
                minobs = 1,
                nboot = 1,
                n0 = 1e4,
                nb = 1e6,
                pm = 1e-4,
                krig_bootstrap_mode = "joint"
              )
            },
            "mock Krig failure"
          )
        },
        Krig = function(...) stop("mock Krig failure"),
        .package = "fields"
      )
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    .package = "alfakR"
  )
})

test_that("ABM treats zero and negative fitness as no-division and validates p_missegregation", {
  zero_res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(10), "1"),
    fitness_map_r = stats::setNames(list(0), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )
  expect_equal(as.numeric(zero_res[["1"]]), 10)

  negative_res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(10), "1"),
    fitness_map_r = stats::setNames(list(-1), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )
  expect_equal(as.numeric(negative_res[["1"]]), 10)

  expect_error(
    alfakR:::run_karyotype_abm(
      initial_population_r = stats::setNames(list(10), "1"),
      fitness_map_r = stats::setNames(list(1), "1"),
      p_missegregation = 1.2,
      dt = 1,
      n_steps = 1L,
      max_population_size = 0,
      culling_survival_fraction = 0.1,
      record_interval = 1L,
      seed = 123L,
      grf_centroids = matrix(0, 0, 0),
      grf_lambda = NA_real_
    ),
    "p_missegregation"
  )
})

test_that("ABM supports large population counts without int-sized random-distribution limits", {
  res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(3e9), "1"),
    fitness_map_r = stats::setNames(list(0), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )
  expect_equal(as.numeric(res[["1"]]), 3e9)
})

test_that("largest remainder allocation is deterministic and preserves the total exactly", {
  alloc <- alfakR:::largest_remainder_allocate(c(0.34, 0.33, 0.33), 2)
  expect_equal(alloc, c(1, 1, 0), tolerance = 0)
  expect_true(is.double(alloc))
  expect_true(all(alloc == floor(alloc)))
  expect_identical(sum(alloc), 2)

  alloc_exact <- alfakR:::largest_remainder_allocate(c(0.01, 0.01, 0.98), 7)
  expect_identical(sum(alloc_exact), 7)
})

test_that("predict_evo ABM with times = 0 returns the initial composition", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.25, "3.1" = 0.75)

  res <- suppressMessages(
    alfakR::predict_evo(
      lscape = lscape,
      p = 0.01,
      times = 0,
      x0 = x0,
      prediction_type = "ABM",
      abm_pop_size = 100,
      abm_delta_t = 0.1,
      abm_record_interval = 1,
      abm_seed = 1
    )
  )

  expect_identical(res$time, 0)
  expect_equal(as.numeric(res[1, c("2.2", "3.1")]), unname(x0), tolerance = 1e-12)
})

test_that("ABM public wrappers accept abm_record_interval = -1 and reject other negative values", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.25, "3.1" = 0.75)
  captured_interval <- NULL

  testthat::with_mocked_bindings(
    {
      res <- suppressMessages(
        alfakR::predict_evo(
          lscape = lscape,
          p = 0.01,
          times = c(0, 0.1),
          x0 = x0,
          prediction_type = "ABM",
          abm_pop_size = 100,
          abm_delta_t = 0.1,
          abm_record_interval = -1,
          abm_seed = 1
        )
      )
      expect_identical(captured_interval, -1L)
      expect_true(is.data.frame(res))
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      captured_interval <<- record_interval
      list(
        "0" = stats::setNames(c(25L, 75L), c("2.2", "3.1")),
        "1" = stats::setNames(c(25L, 75L), c("2.2", "3.1"))
      )
    },
    .package = "alfakR"
  )

  captured_grf_interval <- NULL
  testthat::with_mocked_bindings(
    {
      res <- alfakR::run_abm_simulation_grf(
        centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
        lambda = 1,
        p = 0.01,
        times = c(0, 0.1),
        x0 = c("2.2" = 1),
        abm_pop_size = 100,
        abm_delta_t = 0.1,
        abm_record_interval = -1,
        abm_seed = 1
      )
      expect_identical(captured_grf_interval, -1L)
      expect_true(is.data.frame(res))
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      captured_grf_interval <<- record_interval
      list(
        "0" = stats::setNames(c(100L), "2.2"),
        "1" = stats::setNames(c(100L), "2.2")
      )
    },
    .package = "alfakR"
  )

  expect_error(
    alfakR::predict_evo(
      lscape = lscape,
      p = 0.01,
      times = c(0, 0.1),
      x0 = x0,
      prediction_type = "ABM",
      abm_pop_size = 100,
      abm_delta_t = 0.1,
      abm_record_interval = -2,
      abm_seed = 1
    ),
    "abm_record_interval"
  )
})

test_that("pij validates arguments and returns a valid probability", {
  expect_equal(alfakR::pij(0, 0, 0.1), 1)
  expect_equal(alfakR::pij(0, 2, 0.1), 0)
  expect_gte(alfakR::pij(2, 0, 0.1), 0)
  expect_lte(alfakR::pij(2, 0, 0.1), 1)
  expect_error(alfakR::pij(-1, 2, 0.1), "`i`")
  expect_error(alfakR::pij(2, -1, 0.1), "`j`")
  expect_error(alfakR::pij(2, 2, NA_real_), "`beta`")
  expect_error(alfakR::pij(2, 2, 1.5), "`beta`")
  val <- alfakR::pij(2, 2, 0.1)
  expect_true(is.finite(val))
  expect_gte(val, 0)
  expect_lte(val, 1)
})

test_that("prediction helpers reject invalid probabilities, times, and x0 inputs", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.5, "3.1" = 0.5)

  expect_error(alfakR::predict_evo(lscape, p = NA_real_, times = c(0, 1), x0 = x0), "`p`")
  expect_error(alfakR::predict_evo(lscape, p = NaN, times = c(0, 1), x0 = x0), "`p`")
  expect_error(alfakR::predict_evo(lscape, p = 1.2, times = c(0, 1), x0 = x0), "`p`")
  expect_error(alfakR::predict_evo(lscape, p = 0.1, times = c(0, Inf), x0 = x0), "`times`")
  expect_error(alfakR::predict_evo(lscape, p = 0.1, times = c(0, 1), x0 = c("2.2" = NA_real_, "3.1" = 1)), "`x0`")

  expect_error(
    alfakR::run_abm_simulation_grf(
      centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
      lambda = 1,
      p = 1.2,
      times = c(0, 1),
      x0 = c("2.2" = 1)
    ),
    "`p`"
  )
  expect_error(
    alfakR::run_abm_simulation_grf(
      centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
      lambda = 1,
      p = 0.1,
      times = c(0, Inf),
      x0 = c("2.2" = 1)
    ),
    "`times`"
  )
})
