resolve_benchmark_repo_dir <- function() {
  command_args <- commandArgs()
  script_arg <- grep("^--file=", command_args, value = TRUE)
  script_path <- if (length(script_arg)) {
    normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }

  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE)
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

run_benchmark_pipeline_entry <- function(params, repo_dir = resolve_benchmark_repo_dir()) {
  source_benchmark_modules(repo_dir, envir = parent.frame())
  ctx <- build_benchmark_context(params = params, repo_dir = repo_dir)
  run_benchmark_pipeline(ctx)
}
