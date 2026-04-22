joint_bayes_validate_nn_prior_mode <- function(nn_prior) {
  nn_prior <- validate_nn_prior_mode(nn_prior)
  if (identical(nn_prior, "empirical")) {
    stop(
      "nn_prior = 'empirical' is not implemented for fit_mode = 'joint_bayes' in Stage 1",
      call. = FALSE
    )
  }
  nn_prior
}

joint_bayes_search_interval <- function(fitness_values) {
  interval <- range(fitness_values, na.rm = TRUE)
  span <- diff(interval)
  if (!is.finite(span) || span <= 0) {
    baseline <- max(abs(interval), na.rm = TRUE)
    if (!is.finite(baseline) || baseline <= 0) {
      baseline <- 1
    }
    span <- baseline
  }
  interval + c(-span, span)
}

joint_bayes_to_interval <- function(raw, interval) {
  midpoint <- mean(interval)
  half_width <- diff(interval) / 2
  if (!is.finite(half_width) || half_width <= 0) {
    return(rep(midpoint, length(raw)))
  }
  midpoint + half_width * tanh(raw)
}

joint_bayes_from_interval <- function(value, interval) {
  midpoint <- mean(interval)
  half_width <- diff(interval) / 2
  if (!is.finite(half_width) || half_width <= 0) {
    return(rep(0, length(value)))
  }
  scaled <- (value - midpoint) / half_width
  scaled <- pmin(pmax(scaled, -1 + 1e-8), 1 - 1e-8)
  atanh(scaled)
}

joint_bayes_scale_fitness <- function(f_rel, x0, g0_val, correct_efflux, viability_vec) {
  if (!isTRUE(correct_efflux)) {
    return(f_rel + g0_val - sum(x0 * f_rel))
  }

  sum_weighted_frel <- sum((x0 * f_rel) / viability_vec)
  sum_weights <- sum(x0 / viability_vec)
  k_const <- (sum_weighted_frel - g0_val) / sum_weights
  (f_rel - k_const) / viability_vec
}

joint_bayes_prefit <- function(data, fq, timepoints, epsilon, n0, nb, viability, correct_efflux) {
  counts_fq <- data$x[fq, , drop = FALSE]
  x <- normalize_columns(counts_fq)

  if (length(fq) == 1) {
    f_qp <- 0
    x0_init <- 1
    opt_res <- list(f = 0, x0 = 1)
  } else {
    dx_dt <- compute_dx_dt(x, timepoints)
    x_trim <- x[, -1, drop = FALSE]
    qr_terms <- alfak_qr_accum_cpp(x_trim, dx_dt)
    Dmat <- 2 * qr_terms$Q_accum + diag(epsilon, length(fq))
    dvec <- 2 * qr_terms$r_accum
    qp_sol <- run_solve_qp_checked(
      Dmat,
      dvec,
      matrix(1, nrow = length(fq), ncol = 1),
      0,
      meq = 1,
      context = "solve.QP deterministic pre-fit"
    )
    f_qp <- qp_sol$solution
    x0_init <- optimize_initial_frequencies(x, f_qp, timepoints)
    opt_res <- joint_optimize(counts_fq, timepoints, f_qp, x0_init)
  }

  g0_val <- log(nb / n0) / diff(timepoints)[1]
  viability_vec <- unname(viability[fq])
  f_abs <- joint_bayes_scale_fitness(
    f_rel = opt_res$f,
    x0 = opt_res$x0,
    g0_val = g0_val,
    correct_efflux = correct_efflux,
    viability_vec = viability_vec
  )

  birth_times <- find_birth_times(
    list(f = f_abs, x0 = opt_res$x0),
    time_range = c(-1000, max(timepoints)),
    minF = 1 / n0
  )
  peak_times <- timepoints[apply(x, 1, which.max)]
  birth_times <- sanitize_birth_times(birth_times, peak_times = peak_times, timepoints = timepoints)

  names(f_qp) <- fq
  names(opt_res$f) <- fq
  names(opt_res$x0) <- fq
  names(f_abs) <- fq
  names(birth_times) <- fq

  list(
    f_qp = f_qp,
    x0_init = x0_init,
    f_rel = opt_res$f,
    x0 = opt_res$x0,
    f_abs = f_abs,
    g0 = g0_val,
    birth_times = birth_times
  )
}

joint_bayes_make_spec <- function(fq, nn_names, nn_prior, nn_prior_sd, nn_prior_sd_floor, nn_interval) {
  n_frequent_free <- max(length(fq) - 1, 0)
  n_logits <- max(length(fq) - 1, 0)
  n_nn <- length(nn_names)
  use_prior <- identical(nn_prior, "empirical_censored") && n_nn > 0
  infer_prior_sd <- use_prior && is.null(nn_prior_sd)

  list(
    fq = fq,
    nn_names = nn_names,
    n_frequent_free = n_frequent_free,
    n_logits = n_logits,
    n_nn = n_nn,
    use_prior = use_prior,
    infer_prior_sd = infer_prior_sd,
    fixed_prior_sd = if (use_prior && !is.null(nn_prior_sd)) nn_prior_sd else NA_real_,
    nn_prior_sd_floor = nn_prior_sd_floor,
    nn_interval = nn_interval
  )
}

joint_bayes_pack_params <- function(f_rel, x0, nn_fitness, spec, mu_delta = 0, sigma_delta = NULL) {
  pieces <- list()
  if (spec$n_frequent_free > 0) {
    pieces[[length(pieces) + 1]] <- f_rel[seq_len(spec$n_frequent_free)]
    pieces[[length(pieces) + 1]] <- free_softmax_logits(x0)
  }
  if (spec$n_nn > 0) {
    pieces[[length(pieces) + 1]] <- joint_bayes_from_interval(nn_fitness, spec$nn_interval)
  }
  if (isTRUE(spec$use_prior)) {
    pieces[[length(pieces) + 1]] <- mu_delta
    if (isTRUE(spec$infer_prior_sd)) {
      sigma_delta <- if (is.null(sigma_delta)) spec$nn_prior_sd_floor * 2 else sigma_delta
      sigma_raw <- log(max(sigma_delta - spec$nn_prior_sd_floor, 1e-6))
      pieces[[length(pieces) + 1]] <- sigma_raw
    }
  }

  if (!length(pieces)) {
    return(numeric(0))
  }
  unlist(pieces, use.names = FALSE)
}

joint_bayes_unpack_params <- function(par, spec) {
  cursor <- 1L

  if (spec$n_frequent_free > 0) {
    f_free <- par[cursor:(cursor + spec$n_frequent_free - 1)]
    cursor <- cursor + spec$n_frequent_free
    x0_logits <- par[cursor:(cursor + spec$n_logits - 1)]
    cursor <- cursor + spec$n_logits
    f_rel <- c(f_free, -sum(f_free))
    x0 <- softmax_from_free_logits(x0_logits)
  } else {
    f_free <- numeric(0)
    x0_logits <- numeric(0)
    f_rel <- 0
    x0 <- 1
  }

  if (spec$n_nn > 0) {
    nn_raw <- par[cursor:(cursor + spec$n_nn - 1)]
    cursor <- cursor + spec$n_nn
    nn_fitness <- joint_bayes_to_interval(nn_raw, spec$nn_interval)
  } else {
    nn_raw <- numeric(0)
    nn_fitness <- numeric(0)
  }

  mu_delta <- NA_real_
  sigma_delta <- NA_real_
  log_sigma_raw <- NA_real_
  if (isTRUE(spec$use_prior)) {
    mu_delta <- par[cursor]
    cursor <- cursor + 1L
    if (isTRUE(spec$infer_prior_sd)) {
      log_sigma_raw <- par[cursor]
      cursor <- cursor + 1L
      sigma_delta <- spec$nn_prior_sd_floor + exp(log_sigma_raw)
    } else {
      sigma_delta <- spec$fixed_prior_sd
    }
  }

  list(
    f_free = f_free,
    x0_logits = x0_logits,
    f_rel = f_rel,
    x0 = x0,
    nn_raw = nn_raw,
    nn_fitness = nn_fitness,
    mu_delta = mu_delta,
    sigma_delta = sigma_delta,
    log_sigma_raw = log_sigma_raw
  )
}

joint_bayes_draw_laplace_samples <- function(map_par, objective_fn, n_draws, context) {
  validate_positive_integer(n_draws, "nboot")
  if (!length(map_par)) {
    return(matrix(numeric(0), nrow = n_draws, ncol = 0))
  }

  hessian <- try(stats::optimHess(map_par, objective_fn), silent = TRUE)
  if (inherits(hessian, "try-error") ||
      !is.matrix(hessian) ||
      any(!is.finite(hessian))) {
    warning(sprintf("%s: Hessian calculation failed; using degenerate Laplace draws at the MAP.", context))
    return(matrix(rep(map_par, times = n_draws), nrow = n_draws, byrow = TRUE))
  }

  hessian <- 0.5 * (hessian + t(hessian))
  eig <- try(eigen(hessian, symmetric = TRUE), silent = TRUE)
  if (inherits(eig, "try-error") ||
      any(!is.finite(eig$values)) ||
      any(!is.finite(eig$vectors))) {
    warning(sprintf("%s: Hessian eigendecomposition failed; using degenerate Laplace draws at the MAP.", context))
    return(matrix(rep(map_par, times = n_draws), nrow = n_draws, byrow = TRUE))
  }

  positive_vals <- eig$values[eig$values > 0 & is.finite(eig$values)]
  eigen_floor <- if (length(positive_vals)) {
    max(stats::median(positive_vals) * 1e-4, 1e-2)
  } else {
    1e-2
  }
  precision_vals <- pmax(eig$values, eigen_floor)
  transform_mat <- diag(1 / sqrt(precision_vals), nrow = length(precision_vals)) %*% t(eig$vectors)
  noise <- matrix(stats::rnorm(n_draws * length(map_par)), nrow = n_draws)
  sweep(noise %*% transform_mat, 2, map_par, `+`)
}

solve_fitness_joint_bayes <- function(data, minobs, nboot = 1000, epsilon = 1e-6, pm = 0.00005,
                                      n0, nb, passage_times = NULL, allow_noninteger_counts = FALSE,
                                      correct_efflux = FALSE,
                                      nn_prior = c("empirical_censored", "none", "empirical"),
                                      nn_prior_sd = NULL,
                                      nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR) {
  data$x <- coerce_count_matrix(data$x, allow_noninteger_counts = allow_noninteger_counts)
  validate_positive_depth(data$x)
  validate_positive_integer(nboot, "nboot")
  validate_positive_finite(n0, "n0")
  validate_positive_finite(nb, "nb")
  validate_probability(pm, "pm", upper_inclusive = TRUE)
  validate_scalar_logical(allow_noninteger_counts, "allow_noninteger_counts")
  validate_scalar_logical(correct_efflux, "correct_efflux")
  validate_positive_finite(nn_prior_sd_floor, "nn_prior_sd_floor")
  if (!is.null(nn_prior_sd)) {
    validate_positive_finite(nn_prior_sd, "nn_prior_sd")
  }
  nn_prior <- joint_bayes_validate_nn_prior_mode(nn_prior)

  fq <- get_frequent_karyotypes(data$x, minobs)
  nn_info <- gen_nn_info(fq, pm)
  if (length(nn_info) > 0 && !is.null(nn_info[[1]]$ni)) {
    names(nn_info) <- vapply(nn_info, function(item) item$ni, character(1))
  }

  fq_vec <- do.call(rbind, lapply(fq, s2v))
  rownames(fq_vec) <- fq
  viability <- prepare_efflux_viability(fq_vec, pm = pm, correct_efflux = correct_efflux)
  timepoints <- resolve_time_axis(data, passage_times)
  ntot <- round(colSums(data$x))

  prefit <- joint_bayes_prefit(
    data = data,
    fq = fq,
    timepoints = timepoints,
    epsilon = epsilon,
    n0 = n0,
    nb = nb,
    viability = viability,
    correct_efflux = correct_efflux
  )

  nn_names <- if (length(nn_info)) names(nn_info) else character(0)
  nn_interval <- joint_bayes_search_interval(prefit$f_abs)
  spec <- joint_bayes_make_spec(
    fq = fq,
    nn_names = nn_names,
    nn_prior = nn_prior,
    nn_prior_sd = nn_prior_sd,
    nn_prior_sd_floor = nn_prior_sd_floor,
    nn_interval = nn_interval
  )

  nn_parent_indices <- lapply(nn_info, function(item) {
    parent_idx <- match(item$nj, fq)
    if (anyNA(parent_idx)) {
      stop("Nearest-neighbour parent mapping failed for joint_bayes objective.", call. = FALSE)
    }
    as.integer(parent_idx)
  })
  nn_pij_values <- lapply(nn_info, function(item) as.numeric(item$pij))
  child_obs <- matrix(0, nrow = length(nn_names), ncol = length(timepoints),
                      dimnames = list(nn_names, NULL))
  if (length(nn_names)) {
    for (i in seq_along(nn_names)) {
      child_name <- nn_names[i]
      if (child_name %in% rownames(data$x)) {
        child_obs[i, ] <- as.numeric(data$x[child_name, ])
      }
    }
  }

  parent_means_start <- if (length(nn_info)) {
    vapply(nn_info, weighted_parent_fitness, numeric(1), fpar = prefit$f_abs)
  } else {
    numeric(0)
  }
  fallback_mean <- mean(prefit$f_abs)
  if (!is.finite(fallback_mean)) {
    fallback_mean <- 0
  }
  parent_means_start[!is.finite(parent_means_start)] <- fallback_mean

  nn_start <- if (length(parent_means_start)) {
    pmin(pmax(parent_means_start, nn_interval[1] + 1e-6), nn_interval[2] - 1e-6)
  } else {
    numeric(0)
  }
  delta_start <- nn_start - parent_means_start
  mu_start <- if (length(delta_start) && any(is.finite(delta_start))) {
    mean(delta_start[is.finite(delta_start)])
  } else {
    0
  }
  sigma_start <- if (length(delta_start) >= 2 && any(is.finite(delta_start))) {
    stats::sd(delta_start[is.finite(delta_start)])
  } else {
    NA_real_
  }
  if (!is.finite(sigma_start) || sigma_start <= 0) {
    sigma_start <- nn_prior_sd_floor * 2
  }
  sigma_start <- max(sigma_start, nn_prior_sd_floor * 2)

  map_start <- joint_bayes_pack_params(
    f_rel = prefit$f_rel,
    x0 = prefit$x0,
    nn_fitness = nn_start,
    spec = spec,
    mu_delta = mu_start,
    sigma_delta = sigma_start
  )

  counts_fq <- data$x[fq, , drop = FALSE]
  viability_vec <- unname(viability[fq])
  weak_mu_sd <- 5
  weak_log_sigma_sd <- 2

  objective_fn <- function(par) {
    decoded <- joint_bayes_unpack_params(par, spec)
    total <- alfak_joint_objective_cpp(
      f_rel = decoded$f_rel,
      x0 = decoded$x0,
      nn_fitness = decoded$nn_fitness,
      counts_fq = counts_fq,
      timepoints = timepoints,
      viability_vec = viability_vec,
      g0_val = prefit$g0,
      correct_efflux = correct_efflux,
      nn_parent_indices = nn_parent_indices,
      nn_pij_values = nn_pij_values,
      birth_times = unname(prefit$birth_times),
      child_obs = child_obs,
      ntot = ntot,
      use_prior = isTRUE(spec$use_prior),
      mu_delta = if (isTRUE(spec$use_prior)) decoded$mu_delta else 0,
      sigma_delta = if (isTRUE(spec$use_prior)) decoded$sigma_delta else 1,
      weak_mu_sd = weak_mu_sd,
      apply_sigma_regularization = isTRUE(spec$use_prior) && isTRUE(spec$infer_prior_sd),
      log_sigma_raw = if (isTRUE(spec$infer_prior_sd)) decoded$log_sigma_raw else 0,
      weak_log_sigma_sd = weak_log_sigma_sd,
      tol = ALFAK_FEXP_DELTA_TOL
    )
    if (!is.finite(total)) {
      return(1e12)
    }
    total
  }

  map_opt <- if (length(map_start)) {
    run_nlminb_strict_checked(
      start = map_start,
      objective = objective_fn,
      control = list(iter.max = 300, eval.max = 600),
      context = "joint_bayes MAP fit"
    )
  } else {
    list(par = numeric(0), objective = objective_fn(numeric(0)))
  }

  posterior_raw <- joint_bayes_draw_laplace_samples(
    map_par = map_opt$par,
    objective_fn = objective_fn,
    n_draws = nboot,
    context = "joint_bayes Laplace approximation"
  )

  final_fitness <- matrix(NA_real_, nrow = nboot, ncol = length(fq), dimnames = list(NULL, fq))
  final_frequencies <- matrix(NA_real_, nrow = nboot, ncol = length(fq), dimnames = list(NULL, fq))
  nn_fitness <- matrix(NA_real_, nrow = nboot, ncol = length(nn_names), dimnames = list(NULL, nn_names))

  for (draw_idx in seq_len(nboot)) {
    decoded <- joint_bayes_unpack_params(posterior_raw[draw_idx, ], spec)
    final_fitness[draw_idx, ] <- joint_bayes_scale_fitness(
      f_rel = decoded$f_rel,
      x0 = decoded$x0,
      g0_val = prefit$g0,
      correct_efflux = correct_efflux,
      viability_vec = viability_vec
    )
    final_frequencies[draw_idx, ] <- decoded$x0
    if (length(nn_names)) {
      nn_fitness[draw_idx, ] <- decoded$nn_fitness
    }
  }

  initial_fitness <- matrix(
    rep(prefit$f_abs, times = nboot),
    nrow = nboot,
    byrow = TRUE,
    dimnames = list(NULL, fq)
  )
  initial_frequencies <- matrix(
    rep(prefit$x0, times = nboot),
    nrow = nboot,
    byrow = TRUE,
    dimnames = list(NULL, fq)
  )
  if (!length(nn_names)) {
    nn_fitness <- matrix(numeric(0), nrow = nboot, ncol = 0)
  }

  map_decoded <- joint_bayes_unpack_params(map_opt$par, spec)
  map_fitness <- joint_bayes_scale_fitness(
    f_rel = map_decoded$f_rel,
    x0 = map_decoded$x0,
    g0_val = prefit$g0,
    correct_efflux = correct_efflux,
    viability_vec = viability_vec
  )
  names(map_fitness) <- fq
  names(map_decoded$x0) <- fq
  if (length(nn_names)) {
    names(map_decoded$nn_fitness) <- nn_names
  }

  list(
    initial_fitness = initial_fitness,
    final_fitness = final_fitness,
    initial_frequencies = initial_frequencies,
    final_frequencies = final_frequencies,
    nn_fitness = nn_fitness,
    fit_mode = "joint_bayes",
    draw_type = "laplace",
    map = list(
      frequent_fitness = map_fitness,
      initial_frequencies = map_decoded$x0,
      nn_fitness = map_decoded$nn_fitness,
      prior_mean = map_decoded$mu_delta,
      prior_sd = map_decoded$sigma_delta
    ),
    prefit = list(
      frequent_fitness = prefit$f_abs,
      initial_frequencies = prefit$x0,
      birth_times = prefit$birth_times
    )
  )
}
