#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Run only the GRF in-silico nn_prior benchmark.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_grf_nn_prior_benchmark.R [options]\n\n",
    "Common options:\n",
    "  --n-sim=2                         Number of independent GRF simulations\n",
    "  --lambdas=0.4,1.6                 Comma-separated GRF lambda values\n",
    "  --training-windows=2,4            Comma-separated numbers of training time points\n",
    "  --nboot=5                         Bootstrap replicates per ALFA-K fit\n",
    "  --n-cores=4                       Worker cores used by ALFA-K internals\n",
    "  --methods=none,empirical          Comma-separated nn_prior methods or nn_prior_* labels\n",
    "  --output-dir=benchmark/results/grf_insilico\n",
    "  --force-refit=true                Ignore cached GRF fit outputs\n\n",
    "Advanced options:\n",
    "  --seed=424242\n",
    "  --pm=5e-05\n",
    "  --minobs=5,10,20\n",
    "  --grid-n=81\n",
    "  --k-dim=22\n",
    "  --n-centroids=64\n",
    "  --time-max=140\n",
    "  --passage-interval=20\n",
    "  --sample-depth=2000\n",
    "  --abm-pop-size=50000\n",
    "  --abm-max-pop=2000000\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--", arg)) {
      stop("Unexpected positional argument: ", arg)
    }
    key_value <- sub("^--", "", arg)
    if (grepl("=", key_value, fixed = TRUE)) {
      key <- sub("=.*$", "", key_value)
      value <- sub("^[^=]*=", "", key_value)
    } else {
      key <- key_value
      value <- "true"
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- value
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || !length(value) || !nzchar(as.character(value[[1L]]))) {
    return(default)
  }
  value[[1L]]
}

arg_logical <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) {
    return(default)
  }
  value <- tolower(trimws(as.character(value)))
  if (value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop("Expected a boolean value for --", gsub("_", "-", name), ".")
}

arg_numeric <- function(args, name, default) {
  value <- suppressWarnings(as.numeric(arg_value(args, name, default)))
  if (!is.finite(value)) {
    return(default)
  }
  value
}

arg_integer <- function(args, name, default) {
  value <- suppressWarnings(as.integer(arg_value(args, name, default)))
  if (!is.finite(value)) {
    return(default)
  }
  value
}

arg_numeric_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) {
    return(default)
  }
  out <- suppressWarnings(as.numeric(strsplit(as.character(value), ",", fixed = TRUE)[[1L]]))
  out <- out[is.finite(out)]
  if (!length(out)) {
    return(default)
  }
  out
}

arg_integer_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) {
    return(default)
  }
  out <- suppressWarnings(as.integer(strsplit(as.character(value), ",", fixed = TRUE)[[1L]]))
  out <- out[is.finite(out)]
  if (!length(out)) {
    return(default)
  }
  out
}

arg_character_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) {
    return(default)
  }
  out <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  out <- out[nzchar(out)]
  if (!length(out)) {
    return(default)
  }
  out
}

resolve_script_repo_dir <- function() {
  command_args <- commandArgs(FALSE)
  script_arg <- grep("^--file=", command_args, value = TRUE)
  start <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }

  candidates <- unique(normalizePath(
    file.path(start, c(".", "..", "../..", "../../..")),
    winslash = "/",
    mustWork = FALSE
  ))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) &&
        dir.exists(file.path(cand, "benchmark"))) {
      return(cand)
    }
  }
  stop("Could not locate the alfakR repository root.")
}

source_benchmark_modules <- function(repo_dir, envir = parent.frame()) {
  module_paths <- file.path(
    repo_dir,
    "benchmark",
    "R",
    c(
      "benchmark_utils.R",
      "benchmark_inputs.R",
      "benchmark_fit_tasks.R",
      "benchmark_compare_global.R",
      "benchmark_compare_focus.R",
      "benchmark_nn_diagnostics.R",
      "benchmark_pipeline.R"
    )
  )
  invisible(lapply(module_paths, function(path) sys.source(path, envir = envir)))
}

make_grf_params <- function(args) {
  methods <- arg_character_vec(
    args,
    "methods",
    arg_character_vec(
      args,
      "parameter_labels",
      c(
        "nn_prior_none",
        "nn_prior_empirical",
        "nn_prior_empirical_censored",
        "nn_prior_empirical_censored_weighted",
        "nn_prior_empirical_two_shell"
      )
    )
  )

  nboot <- arg_integer(args, "nboot", arg_integer(args, "nn_grf_nboot", 5L))

  list(
    patient_subset = character(0),
    minobs_values = arg_integer_vec(args, "minobs", 5L),
    pm_values = arg_numeric_vec(args, "pm", 5e-05),
    parameter_labels = methods,
    n_cores = arg_integer(args, "n_cores", 4L),
    nboot = nboot,
    n0 = arg_numeric(args, "n0", 100000),
    nb = arg_numeric(args, "nb", 10000000),
    benchmark_seed = arg_integer(args, "benchmark_seed", 31415926L),
    correct_efflux = arg_logical(args, "correct_efflux", TRUE),
    nn_prior_grid_n = arg_integer(args, "grid_n", arg_integer(args, "nn_prior_grid_n", 81L)),
    nn_prior_fit_subset = as.character(arg_value(args, "nn_prior_fit_subset", "hybrid")),
    nn_prior_zero_exposure_quantile = arg_numeric(args, "nn_prior_zero_exposure_quantile", 0.10),
    nn_prior_zero_weight_scale = arg_numeric(args, "nn_prior_zero_weight_scale", 0.50),
    nn_prior_zero_weight_cap_ratio = NA_real_,
    nn_prior_zero_birth_fallback_weight = NA_real_,
    nn_prior_zero_birth_child_floor = arg_numeric(args, "nn_prior_zero_birth_child_floor", 0.25),
    nn_prior_zero_birth_child_shape = arg_numeric(args, "nn_prior_zero_birth_child_shape", 1),
    nn_prior_zero_birth_replicate_floor = arg_numeric(args, "nn_prior_zero_birth_replicate_floor", 0.50),
    nn_prior_zero_birth_replicate_shape = arg_numeric(args, "nn_prior_zero_birth_replicate_shape", 1),
    nn_prior_two_step_support = as.character(arg_value(args, "nn_prior_two_step_support", "rescue")),
    nn_prior_two_step_support_min = arg_numeric(args, "nn_prior_two_step_support_min", 0.15),
    nn_prior_two_step_cap_floor = arg_numeric(args, "nn_prior_two_step_cap_floor", 0.30),
    cohort_transition_version = "contextual",
    cohort_contextual_apply_to = "all",
    cohort_context_keep_baseline_when_sparse = FALSE,
    cohort_context_lambda_sparse_unknown = 0.10,
    include_posterior_comparison = FALSE,
    top_shift_n = 15L,
    run_nn_identifiability = FALSE,
    run_nn_stability = FALSE,
    run_nn_holdout = FALSE,
    nn_holdout_repeats = 1L,
    nn_holdout_fraction = 0.25,
    nn_holdout_min_count = 2L,
    nn_holdout_seed = 271828L,
    run_nn_simulation = FALSE,
    nn_simulation_n = 1L,
    nn_simulation_seed = 161803L,
    nn_simulation_scenarios = "sparse_zero_heavy",
    run_nn_grf_simulation = TRUE,
    nn_grf_simulation_n = arg_integer(args, "n_sim", arg_integer(args, "nn_grf_simulation_n", 2L)),
    nn_grf_seed = arg_integer(args, "seed", arg_integer(args, "nn_grf_seed", 424242L)),
    nn_grf_lambdas = arg_numeric_vec(args, "lambdas", c(0.4, 1.6)),
    nn_grf_training_windows = arg_integer_vec(args, "training_windows", c(2L, 4L)),
    nn_grf_nboot = nboot,
    nn_grf_k_dim = arg_integer(args, "k_dim", 22L),
    nn_grf_n_centroids = arg_integer(args, "n_centroids", 64L),
    nn_grf_time_max = arg_numeric(args, "time_max", 140),
    nn_grf_passage_interval = arg_numeric(args, "passage_interval", 20),
    nn_grf_sample_depth = arg_integer(args, "sample_depth", 2000L),
    nn_grf_abm_pop_size = arg_numeric(args, "abm_pop_size", 50000),
    nn_grf_abm_delta_t = arg_numeric(args, "abm_delta_t", 1),
    nn_grf_abm_max_pop = arg_numeric(args, "abm_max_pop", 2000000),
    nn_grf_abm_culling_survival = arg_numeric(args, "abm_culling_survival", 0.01),
    force_refit = arg_logical(args, "force_refit", FALSE),
    rebuild_inputs = FALSE,
    run_benchmark = FALSE,
    render_figures = FALSE,
    focus_pids = "P2",
    focus_pid = "",
    focus_minobs = 5L,
    focus_pm = 5e-05
  )
}

normalize_output_dir <- function(repo_dir, output_dir) {
  if (grepl("^/", output_dir)) {
    return(normalizePath(output_dir, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(repo_dir, output_dir), winslash = "/", mustWork = FALSE)
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }

  repo_dir <- normalizePath(arg_value(args, "repo_dir", resolve_script_repo_dir()), winslash = "/", mustWork = FALSE)
  source_benchmark_modules(repo_dir, envir = globalenv())

  params <- make_grf_params(args)
  ctx <- build_benchmark_context(params = params, repo_dir = repo_dir)

  output_dir <- normalize_output_dir(
    repo_dir,
    as.character(arg_value(args, "output_dir", "benchmark/results/grf_insilico"))
  )
  ctx$results_dir <- output_dir
  ctx$fit_dir <- file.path(output_dir, "fits")
  ctx$tables_dir <- file.path(output_dir, "tables")
  ctx$figures_dir <- file.path(output_dir, "figures")
  ctx$cache_dir <- file.path(output_dir, "cache")
  dir.create(ctx$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ctx$fit_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ctx$tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ctx$figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ctx$cache_dir, recursive = TRUE, showWarnings = FALSE)

  parameter_spec_tbl <- build_parameter_spec_tbl(ctx$parameter_levels_use)
  save_table_bundle(parameter_spec_tbl, file.path(ctx$tables_dir, "parameter_spec"))

  message("Running GRF in-silico nn_prior benchmark only.")
  message("Output directory: ", ctx$results_dir)
  message("Methods: ", paste(parameter_spec_tbl$parameter_label, collapse = ", "))
  message("Simulations: ", ctx$nn_grf_simulation_n_use)
  message("Lambdas: ", paste(ctx$nn_grf_lambdas_use, collapse = ", "))
  message("Training windows: ", paste(ctx$nn_grf_training_windows_use, collapse = ", "))
  message("MINOBS: ", paste(ctx$minobs_values_use, collapse = ", "))

  grf <- run_nn_grf_simulation_diagnostics(
    ctx = ctx,
    parameter_spec_tbl = parameter_spec_tbl
  )

  save_table_bundle(grf$summary_tbl, file.path(ctx$tables_dir, "nn_grf_simulation_summary"))
  save_table_bundle(grf$by_lambda_tbl, file.path(ctx$tables_dir, "nn_grf_simulation_by_lambda"))
  save_table_bundle(grf$child_tbl, file.path(ctx$tables_dir, "nn_grf_simulation_by_child"))
  save_table_bundle(grf$fit_tbl, file.path(ctx$tables_dir, "nn_grf_simulation_fit_results"))
  save_table_bundle(grf$task_tbl, file.path(ctx$tables_dir, "nn_grf_simulation_tasks"))

  message("Wrote:")
  message("  ", file.path(ctx$tables_dir, "nn_grf_simulation_tasks.tsv"))
  message("  ", file.path(ctx$tables_dir, "nn_grf_simulation_summary.tsv"))
  message("  ", file.path(ctx$tables_dir, "nn_grf_simulation_by_lambda.tsv"))
  message("  ", file.path(ctx$tables_dir, "nn_grf_simulation_by_child.tsv"))
  message("  ", file.path(ctx$tables_dir, "nn_grf_simulation_fit_results.tsv"))

  if (nrow(grf$summary_tbl)) {
    print(grf$summary_tbl)
  } else {
    message("No GRF child-level estimates were produced. Check nn_grf_simulation_fit_results.tsv for fit errors.")
  }

  invisible(grf)
}

main()
