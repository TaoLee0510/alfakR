#' ALFAK: Adaptive Landscape Fitness Inference from Karyotype Dynamics
#'
#' Performs fitness landscape inference using Allele Frequency Karyotype dynamics.
#' This function estimates fitness values for observed karyotypes and their
#' neighbors, performs Kriging to create a fitness landscape, and evaluates
#' the landscape using cross-validation.
#'
#' @param yi A list containing the input data. Expected elements:
#'   \itemize{
#'     \item `x`: A matrix of karyotype counts (rows are karyotypes as strings
#'       like "2.2.1", columns are timepoints). Rownames should be present.
#'       Colnames should represent time values that can be multiplied by `dt`.
#'     \item `dt`: A numeric value representing the scaling factor for timepoints
#'       (e.g., generation time if colnames of `x` are generations).
#'   }
#' @param outdir A string specifying the directory path where result files
#'   (RDS format) will be saved. The directory will be created if it doesn't exist.
#' @param passage_times An optional numeric vector giving the final internal time
#'   axis used by the fitter. If `NULL`, the time axis is calculated as
#'   `as.numeric(colnames(yi$x)) * yi$dt`. When supplied, `passage_times` is used
#'   as-is and must be numeric, finite, strictly increasing, and have length
#'   `ncol(yi$x)`.
#' @param minobs An integer, the minimum total number of observations (reads/counts)
#'   for a karyotype across all timepoints to be considered "frequent" and included
#'   in the analysis. A karyotype is considered frequent when
#'   `rowSums(yi$x) >= minobs`. Default is 20.
#' @param nboot An integer, the number of bootstrap iterations for fitness
#'   estimation and for the Kriging process in `fitKrig`. Default is 45.
#' @param n0 A numeric value, the initial effective population size at the
#'   start of a passage or growth phase, used for g0 calculation. Default is 1e5.
#' @param nb A numeric value, the bottleneck effective population size (population
#'   size after transfer), used for g0 calculation. Default is 1e7.
#' @param pm A numeric value, the per-locus mutation/error rate used in `pij`
#'   calculations. Default is 0.00005.
#' @param allow_noninteger_counts Logical; if `FALSE` (default), non-integer
#'   counts in `yi$x` are rejected. If `TRUE`, they are rounded once at entry
#'   with a warning.
#' @param correct_efflux Logical; if `TRUE`, apply the efflux correction after a
#'   one-time viability pre-check on the frequent karyotypes. The viability term
#'   currently depends only on total copy number through
#'   `2 * (1 - pm)^total_copy_number - 1`, so this is a ploidy-level approximation
#'   rather than a chromosome-specific model.
#' @param landscape_data_output Logical; if `TRUE`, also save the optional
#'   `landscape_data.Rds` file containing the stable Kriging mean and median
#'   model objects. Default is `FALSE`, so only the documented core outputs are
#'   written.
#' @param nn_prior Character; nearest-neighbour prior mode for latent children.
#'   `"empirical_censored"` fits an empirical-Bayes prior from all neighbour
#'   children, including zero-count latent neighbours, to correct observation
#'   bias. This mode errors if the prior hyperparameter fit fails and is the
#'   default.
#'   `"empirical_censored_weighted"` fits the same single-Gaussian censored
#'   empirical-Bayes prior, but only in the bootstrap pathway it computes
#'   projected child exposure from the neighbour likelihood ingredients,
#'   downweights zero-only latent children, applies child-level and
#'   replicate-level birth-time fallback burden multipliers for zero-only
#'   children, optionally screens very low-exposure zeros before fitting the
#'   prior, and when a bootstrap replicate has no observed neighbour children
#'   it first falls back to a weak sample-pooled prior learned from the same
#'   sample's observed neighbours before using no prior only as a final
#'   fallback.
#'   `"empirical_two_shell"` is an opt-in sparse two-timepoint correction. It
#'   first obtains the usual one-step neighbour estimates, estimates supported
#'   two-step descendants from those one-step nodes, learns separate empirical
#'   delta distributions for the 0->1 and 1->2 shells, then re-estimates each
#'   one-step neighbour once with the direct likelihood, inward frequent-parent
#'   prior, and an uncertainty-inflated outward two-step prior. It does not
#'   alternate repeated smoothing passes. If supported two-step evidence is
#'   absent or too small, the replicate falls back to the weighted censored
#'   one-step behaviour.
#'   `"cohort_transition"` uses a cohort-level prior on transition effects,
#'   `child fitness - parent fitness`, learned from upstream patient-specific
#'   two-shell fits. It requires `cohort_transition_prior` or
#'   `cohort_transition_prior_path`, keeps the current patient's frequent
#'   karyotype fitness estimates patient-specific, and does not pool raw patient
#'   counts or absolute cohort fitness values.
#'   `"none"` disables the latent-neighbour prior contribution.
#'   `"empirical"` opt-in uses the empirical child-minus-parent prior estimated
#'   from observed neighbours.
#'   Default is `"empirical_censored"`.
#' @param cohort_transition_prior Optional prior object returned by
#'   `learn_cohort_transition_prior()`. Required when
#'   `nn_prior = "cohort_transition"` unless `cohort_transition_prior_path` is
#'   supplied. If both are supplied, this object is used and a warning is issued.
#' @param cohort_transition_prior_path Optional path to a saved
#'   `cohort_transition_v1` or `cohort_transition_v2` prior object.
#' @param cohort_transition_patient_id Patient identifier for single-patient
#'   cohort-transition refits. Required when the prior contains
#'   leave-one-patient-out priors so the target patient's own two-shell
#'   transitions can be excluded from its prior.
#' @param cohort_transition_version Cohort-transition implementation version.
#'   `"contextual"` is the default for `nn_prior = "cohort_transition"` and
#'   conditions Delta fitness on similar parent karyotype and CNA event context;
#'   `"v2"` is the group-level heterogeneity-aware selective-borrowing overlay;
#'   `"v1"` preserves the original direct cohort-prior behavior.
#' @param cohort_transition_apply_to Which NN nodes can receive the v2 overlay.
#'   The default `"zero_only"` leaves observed NN estimates at their
#'   patient-specific `empirical_two_shell` baseline.
#' @param cohort_transition_overlay_base Baseline used by v2. The default
#'   `"empirical_two_shell"` runs the same patient-specific two-shell path first
#'   and then applies the cohort overlay.
#' @param cohort_transition_lambda Global multiplier for v2 cohort borrowing.
#' @param cohort_transition_max_borrowing_fraction Maximum allowed v2 cohort
#'   borrowing fraction before the update is skipped and marked dominated.
#' @param cohort_transition_max_abs_delta_shift Optional maximum absolute change
#'   from the two-shell baseline.
#' @param cohort_contextual_apply_to Which NN nodes can receive contextual
#'   overlay updates. If `NULL`, inherits `cohort_transition_apply_to`.
#' @param cohort_contextual_overlay_base Baseline used by contextual mode. The
#'   default `"empirical_two_shell"` preserves patient-specific two-shell
#'   estimates before context borrowing.
#' @param cohort_context_lambda Global multiplier for contextual cohort
#'   borrowing.
#' @param cohort_context_max_borrowing_fraction Maximum contextual borrowing
#'   fraction before an update is skipped.
#' @param cohort_context_max_abs_delta_shift Optional maximum contextual shift
#'   from the two-shell baseline.
#' @param cohort_context_sd_floor,cohort_context_patient_sd_floor Minimum
#'   contextual transition SD and patient heterogeneity floors.
#' @param cohort_context_keep_baseline_when_sparse,cohort_context_keep_baseline_when_high_variable
#'   Guardrails that keep the two-shell baseline for sparse or high-variable
#'   contexts.
#' @param cohort_transition_sd_floor Minimum transition-prior standard
#'   deviation used by `nn_prior = "cohort_transition"`.
#' @param cohort_transition_patient_sd_floor Minimum patient-heterogeneity
#'   standard deviation added to cohort transition priors. The default is
#'   intentionally conservative because cohort records can contain many
#'   bootstrap/path rows that should not become an overconfident patient-level
#'   prior.
#' @param nn_prior_sd Optional numeric scalar. If supplied, this overrides the
#'   empirically estimated prior standard deviation for latent-neighbour fitting.
#' @param nn_prior_sd_floor Numeric scalar giving the minimum standard deviation
#'   used when the empirical prior variance is zero or too small. Default is
#'   `1e-3`.
#' @param nn_prior_grid_n Integer; number of equally spaced grid points used for
#'   the fixed-grid numerical integration in `nn_prior = "empirical_censored"`
#'   and `nn_prior = "empirical_censored_weighted"`. Default is `81`.
#' @param nn_prior_fit_subset Character; weighted-mode only control for which
#'   zero-only latent children are allowed into the empirical-Bayes prior fit.
#'   `"hybrid"` (default) keeps all observed neighbour children, then screens
#'   zero-only children by projected exposure before weighting them. `"all"`
#'   skips the hard exposure screen but still downweights zero-only children and
#'   still applies the total zero-weight cap.
#' @param nn_prior_zero_exposure_min Optional non-negative numeric scalar used
#'   only when `nn_prior = "empirical_censored_weighted"` and
#'   `nn_prior_fit_subset = "hybrid"`. If supplied, this is the projected
#'   child-exposure threshold used to retain zero-only children in the prior fit.
#' @param nn_prior_zero_exposure_quantile Numeric scalar in `[0, 1]` used only
#'   for weighted hybrid fitting when `nn_prior_zero_exposure_min` is `NULL` and
#'   enough observed neighbour children are available. The threshold is the
#'   corresponding quantile of the observed projected child exposures.
#' @param nn_prior_zero_weight_scale Numeric scalar in `[0, 1]` giving the
#'   baseline downweighting applied to zero-only latent children in weighted
#'   mode before the global zero-weight cap.
#' @param nn_prior_zero_weight_cap_ratio Optional non-negative numeric scalar
#'   used only in weighted mode. When supplied, the total zero-child weight is
#'   capped at this ratio times the number of observed neighbour children.
#'   When `NULL`, a data-adaptive cap is used based on the retained zero-only
#'   children's effective evidence mass rather than their raw count alone.
#' @param nn_prior_zero_birth_fallback_weight Optional numeric scalar in
#'   `[0, 1]` kept as a compatibility alias for
#'   `nn_prior_zero_birth_child_floor` in weighted mode. When supplied, it
#'   overrides `nn_prior_zero_birth_child_floor`.
#' @param nn_prior_zero_birth_child_floor Numeric scalar in `[0, 1]` used only
#'   in weighted mode. This is the minimum child-level birth-reliability
#'   multiplier applied when a zero-only child's parent support is entirely
#'   driven by fallback-imputed birth times.
#' @param nn_prior_zero_birth_child_shape Non-negative numeric scalar used only
#'   in weighted mode. Larger values make the child-level birth-reliability
#'   multiplier decay more sharply as fallback burden increases.
#' @param nn_prior_zero_birth_replicate_floor Numeric scalar in `[0, 1]` used
#'   only in weighted mode. This is the minimum replicate-level
#'   birth-reliability multiplier when the retained zero-only children in a
#'   bootstrap replicate are collectively dominated by fallback-imputed parent
#'   birth times.
#' @param nn_prior_zero_birth_replicate_shape Non-negative numeric scalar used
#'   only in weighted mode. Larger values make the replicate-level
#'   birth-reliability multiplier decay more sharply as the replicate-wide
#'   fallback burden increases.
#' @param nn_prior_hybrid_min_obs Positive integer used only when
#'   `nn_prior = "empirical_censored_weighted"` and
#'   `nn_prior_fit_subset = "hybrid"`. When fewer than this many observed
#'   neighbour children are available, weighted hybrid fitting does not estimate
#'   a hard observed-based exposure threshold from that small sample and instead
#'   relies on weighting, adaptive capping, and optional 2-step rescue.
#' @param nn_prior_two_step_support Character; weighted-mode only control for
#'   whether observed 2-step descendants can rescue zero-only latent neighbour
#'   children. `"none"` (default) ignores 2-step observed support. `"rescue"`
#'   uses observed descendants one mutation away from a zero-only child to
#'   soften hybrid screening and strengthen that child's effective evidence mass
#'   and weight without changing the single-step neighbour objective itself.
#' @param nn_prior_two_step_support_min Numeric scalar in `[0, 1]` used only
#'   when `nn_prior = "empirical_censored_weighted"` and
#'   `nn_prior_two_step_support = "rescue"`. In hybrid fitting, a zero-only
#'   child is retained if either its projected exposure passes the exposure
#'   threshold or its 2-step support score reaches this minimum.
#' @param nn_prior_two_step_cap_floor Numeric scalar in `[0, 1]` used only when
#'   `nn_prior = "empirical_censored_weighted"` and
#'   `nn_prior_two_step_support = "rescue"`. This sets the minimum fraction of a
#'   fully supported 2-step rescue score that can contribute to the retained
#'   zero-child support term used in effective-mass and weight calculations.
#' @param nn_two_shell_min_delta_n Positive integer used only when
#'   `nn_prior = "empirical_two_shell"`. It is the minimum number of usable
#'   path-weighted 1->2 shell deltas required before fitting the outward prior;
#'   otherwise the replicate falls back to the inward one-step estimate.
#' @param nn_two_shell_min_exposure Optional non-negative numeric scalar used
#'   only for `nn_prior = "empirical_two_shell"`. Unobserved two-step candidates
#'   are retained only when their expected exposure reaches this threshold. When
#'   `NULL`, an adaptive threshold based on the projected one-step exposure
#'   distribution is used and reported in diagnostics.
#' @param nn_two_shell_min_observed_count Non-negative integer count threshold
#'   used only for `nn_prior = "empirical_two_shell"`. Two-step candidates with
#'   at least this total observed count are retained automatically.
#' @param nn_two_shell_max_weight_ratio Non-negative numeric scalar used only
#'   for `nn_prior = "empirical_two_shell"`. It caps each one-step node's total
#'   outward prior weight relative to its inward parent-prior weight.
#' @param nn_two_shell_lambda Non-negative numeric scalar multiplying the
#'   outward two-step prior term.
#' @param nn_two_shell_reuse_sd Optional non-negative numeric scalar. It inflates
#'   the effective standard deviation of the outward prior because provisional
#'   two-step estimates reuse the same count data. When `NULL`, a conservative
#'   replicate-specific value based on `nn_prior_sd_floor` and the estimated
#'   1->2 prior scale is used.
#' @param nn_two_shell_uncertainty_floor Optional non-negative numeric scalar.
#'   Minimum standard error assigned to provisional two-step fitness estimates.
#'   When `NULL`, a conservative package default is used.
#' @param nn_two_shell_save_diagnostics Logical; if `TRUE` and
#'   `nn_prior = "empirical_two_shell"`, `alfak()` saves
#'   `nn_prior_diagnostics.Rds` containing replicate and per-node diagnostics.
#' @param krig_bootstrap_mode Character; `"marginal"` (default) samples
#'   bootstrap fitness values independently by column, matching the original
#'   ALFA-K Kriging bootstrap and cross-validation behavior. `"joint"` samples
#'   one full bootstrap row at a time, preserving within-bootstrap dependence
#'   across karyotypes.
#'
#' @return Returns the cross-validation R-squared value (`Rxv`) invisibly.
#'   The function primarily saves its results to RDS files in the `outdir`:
#'   \itemize{
#'     \item `bootstrap_res.Rds`: Results from `solve_fitness_bootstrap`.
#'     \item `landscape.Rds`: Summary statistics (mean, median, sd) of the
#'       Kriging-inferred fitness landscape from `fitKrig`.
#'     \item `landscape_posterior_samples.Rds`: The full matrix of posterior
#'       samples from the Kriging bootstraps in `fitKrig`.
#'     \item `xval.Rds`: The cross-validation R-squared value (`Rxv`).
#'     \item `landscape_data.Rds`: Optional stable Kriging mean/median model
#'       objects, written only when `landscape_data_output = TRUE`.
#'   }
#'
#' @export
#' @importFrom quadprog solve.QP
#' @importFrom fields Krig
#' @importFrom stats rmultinom optim uniroot dbinom dnorm optimise dist predict median sd setNames complete.cases rnorm
#'
#' @examples
#' \dontrun{
#' # Create dummy data for yi
#' karyotypes_str <- c("2.2.2", "2.2.1", "2.1.2", "1.2.2", "2.2.3")
#' timepoints_num <- 5
#' counts_data <- matrix(
#'   abs(rpois(length(karyotypes_str) * timepoints_num, lambda = 20)), # Use pois for counts
#'   nrow = length(karyotypes_str),
#'   ncol = timepoints_num
#' )
#' rownames(counts_data) <- karyotypes_str
#' colnames(counts_data) <- 1:timepoints_num # Generations
#'
#' dummy_yi_data <- list(
#'   x = counts_data,
#'   dt = 1 # dt = 1 if colnames are generations
#' )
#'
#' temp_output_dir <- tempfile("alfak_example_")
#'
#' # Run alfak
#' result_r_squared <- alfak(
#'   yi = dummy_yi_data,
#'   outdir = temp_output_dir,
#'   passage_times = NULL,
#'   minobs = 5,      # Lowered for dummy data
#'   nboot = 10,      # Lowered for quick example
#'   n0 = 1e4,
#'   nb = 1e6,
#'   pm = 0.0001
#' )
#' print(paste("Cross-validation R-squared:", result_r_squared))
#'
#' # Check for created files
#' list.files(temp_output_dir)
#'
#' # Sparse two-timepoint use: opt into the two-shell correction and inspect
#' # diagnostics if the outward term falls back or has little effect.
#' sparse_counts <- matrix(
#'   c(80, 60,
#'     5, 0,
#'     0, 4),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(c("2.2.2", "2.2.3", "2.2.4"), c("0", "1"))
#' )
#' sparse_dir <- tempfile("alfak_two_shell_")
#' alfak(
#'   yi = list(x = sparse_counts, dt = 1),
#'   outdir = sparse_dir,
#'   minobs = 20,
#'   nboot = 5,
#'   nn_prior = "empirical_two_shell"
#' )
#' readRDS(file.path(sparse_dir, "nn_prior_diagnostics.Rds"))$replicate
#'
#' # Clean up
#' unlink(temp_output_dir, recursive = TRUE)
#' unlink(sparse_dir, recursive = TRUE)
#' }
alfak <- function(yi, outdir, passage_times = NULL, minobs = 20,
                  nboot = 45,
                  n0 = 1e5,
                  nb = 1e7,
                  pm = 0.00005,
                  allow_noninteger_counts = FALSE,
                  correct_efflux=FALSE,
                  landscape_data_output = FALSE,
                  nn_prior = c("empirical_censored", "empirical_censored_weighted", "empirical_two_shell", "cohort_transition", "none", "empirical"),
                  cohort_transition_prior = NULL,
                  cohort_transition_prior_path = NULL,
                  cohort_transition_patient_id = NULL,
                  cohort_transition_version = c("contextual", "v2", "v1"),
                  cohort_transition_apply_to = c("zero_only", "low_information", "all"),
                  cohort_transition_overlay_base = c("empirical_two_shell", "direct"),
                  cohort_transition_lambda = 0.25,
                  cohort_transition_max_borrowing_fraction = 0.5,
                  cohort_transition_max_abs_delta_shift = NULL,
                  cohort_transition_sd_floor = 0.05,
                  cohort_transition_patient_sd_floor = 0.1,
                  cohort_contextual_apply_to = NULL,
                  cohort_contextual_overlay_base = c("empirical_two_shell", "direct"),
                  cohort_context_lambda = 0.25,
                  cohort_context_max_borrowing_fraction = 0.5,
                  cohort_context_max_abs_delta_shift = NULL,
                  cohort_context_sd_floor = 0.05,
                  cohort_context_patient_sd_floor = 0.10,
                  cohort_context_keep_baseline_when_sparse = TRUE,
                  cohort_context_keep_baseline_when_high_variable = TRUE,
                  nn_prior_sd = NULL,
                  nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                  nn_prior_grid_n = ALFAK_NN_PRIOR_CENSORED_GRID_POINTS,
                  nn_prior_fit_subset = c("hybrid", "all"),
                  nn_prior_zero_exposure_min = NULL,
                  nn_prior_zero_exposure_quantile = 0.10,
                  nn_prior_zero_weight_scale = 0.50,
                  nn_prior_zero_weight_cap_ratio = NULL,
                  nn_prior_zero_birth_fallback_weight = NULL,
                  nn_prior_zero_birth_child_floor = 0.25,
                  nn_prior_zero_birth_child_shape = 1,
                  nn_prior_zero_birth_replicate_floor = 0.50,
                  nn_prior_zero_birth_replicate_shape = 1,
                  nn_prior_hybrid_min_obs = 3L,
                  nn_prior_two_step_support = c("none", "rescue"),
                  nn_prior_two_step_support_min = 0.15,
                  nn_prior_two_step_cap_floor = 0.30,
                  nn_two_shell_min_delta_n = 3L,
                  nn_two_shell_min_exposure = NULL,
                  nn_two_shell_min_observed_count = 1L,
                  nn_two_shell_max_weight_ratio = 1.0,
                  nn_two_shell_lambda = 1.0,
                  nn_two_shell_reuse_sd = NULL,
                  nn_two_shell_uncertainty_floor = NULL,
                  nn_two_shell_save_diagnostics = TRUE,
                  krig_bootstrap_mode = c("marginal", "joint")) {

  # Note: library calls removed, dependencies handled by @importFrom or DESCRIPTION

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  alfak_run_log_path(file.path(outdir, "alfak_run.log"))
  validate_positive_integer(nboot, "nboot")
  validate_positive_finite(n0, "n0")
  validate_positive_finite(nb, "nb")
  validate_probability(pm, "pm", upper_inclusive = TRUE)
  validate_scalar_logical(allow_noninteger_counts, "allow_noninteger_counts")
  validate_scalar_logical(correct_efflux, "correct_efflux")
  validate_scalar_logical(landscape_data_output, "landscape_data_output")
  validate_scalar_logical(nn_two_shell_save_diagnostics, "nn_two_shell_save_diagnostics")
  nn_prior <- validate_nn_prior_mode(nn_prior)
  alfak_log_event(
    level = "INFO",
    component = "alfak",
    detail = sprintf("start nn_prior=%s outdir=%s", nn_prior, normalizePath(outdir, mustWork = FALSE))
  )
  cohort_transition_prior <- if (identical(nn_prior, "cohort_transition")) {
    cohort_transition_version <- match.arg(cohort_transition_version)
    cohort_transition_apply_to <- match.arg(cohort_transition_apply_to)
    cohort_transition_overlay_base <- match.arg(cohort_transition_overlay_base)
    if (is.null(cohort_contextual_apply_to)) {
      cohort_contextual_apply_to <- cohort_transition_apply_to
    } else {
      cohort_contextual_apply_to <- match.arg(cohort_contextual_apply_to, c("zero_only", "low_information", "all"))
    }
    cohort_contextual_overlay_base <- match.arg(cohort_contextual_overlay_base)
    validate_nonnegative_finite(cohort_transition_lambda, "cohort_transition_lambda")
    validate_probability(cohort_transition_max_borrowing_fraction, "cohort_transition_max_borrowing_fraction", upper_inclusive = TRUE)
    if (!is.null(cohort_transition_max_abs_delta_shift)) {
      validate_positive_finite(cohort_transition_max_abs_delta_shift, "cohort_transition_max_abs_delta_shift")
    }
    validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
    validate_positive_finite(cohort_transition_patient_sd_floor, "cohort_transition_patient_sd_floor")
    validate_nonnegative_finite(cohort_context_lambda, "cohort_context_lambda")
    validate_probability(cohort_context_max_borrowing_fraction, "cohort_context_max_borrowing_fraction", upper_inclusive = TRUE)
    if (!is.null(cohort_context_max_abs_delta_shift)) {
      validate_positive_finite(cohort_context_max_abs_delta_shift, "cohort_context_max_abs_delta_shift")
    }
    validate_positive_finite(cohort_context_sd_floor, "cohort_context_sd_floor")
    validate_positive_finite(cohort_context_patient_sd_floor, "cohort_context_patient_sd_floor")
    validate_scalar_logical(cohort_context_keep_baseline_when_sparse, "cohort_context_keep_baseline_when_sparse")
    validate_scalar_logical(cohort_context_keep_baseline_when_high_variable, "cohort_context_keep_baseline_when_high_variable")
    resolve_cohort_transition_prior_object(
      cohort_transition_prior = cohort_transition_prior,
      cohort_transition_prior_path = cohort_transition_prior_path,
      cohort_transition_patient_id = cohort_transition_patient_id
    )
  } else {
    cohort_transition_version <- cohort_transition_version[1]
    cohort_transition_apply_to <- cohort_transition_apply_to[1]
    cohort_transition_overlay_base <- cohort_transition_overlay_base[1]
    cohort_contextual_apply_to <- cohort_contextual_apply_to %||% cohort_transition_apply_to
    cohort_contextual_overlay_base <- cohort_contextual_overlay_base[1]
    cohort_transition_prior
  }
  nn_prior_fit_subset <- validate_nn_prior_fit_subset(nn_prior_fit_subset)
  nn_prior_two_step_support <- validate_nn_prior_two_step_support(nn_prior_two_step_support)
  krig_bootstrap_mode <- validate_krig_bootstrap_mode(krig_bootstrap_mode)
  validate_nn_prior_controls(
    nn_prior_sd = nn_prior_sd,
    nn_prior_sd_floor = nn_prior_sd_floor,
    nn_prior_grid_n = nn_prior_grid_n,
    nn_prior_fit_subset = nn_prior_fit_subset,
    nn_prior_zero_exposure_min = nn_prior_zero_exposure_min,
    nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
    nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
    nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
    nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
    nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
    nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
    nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
    nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
    nn_prior_hybrid_min_obs = nn_prior_hybrid_min_obs,
    nn_prior_two_step_support = nn_prior_two_step_support,
    nn_prior_two_step_support_min = nn_prior_two_step_support_min,
    nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
    nn_two_shell_min_delta_n = nn_two_shell_min_delta_n,
    nn_two_shell_min_exposure = nn_two_shell_min_exposure,
    nn_two_shell_min_observed_count = nn_two_shell_min_observed_count,
    nn_two_shell_max_weight_ratio = nn_two_shell_max_weight_ratio,
    nn_two_shell_lambda = nn_two_shell_lambda,
    nn_two_shell_reuse_sd = nn_two_shell_reuse_sd,
    nn_two_shell_uncertainty_floor = nn_two_shell_uncertainty_floor
  )
  yi$x <- coerce_count_matrix(yi$x, allow_noninteger_counts = allow_noninteger_counts)
  validate_positive_depth(yi$x)

  get_frequent_karyotypes(yi$x, minobs)
  resolve_time_axis(yi, passage_times)

  # Parallelism and cl related code removed

  fq_boot <- solve_fitness_bootstrap(yi, minobs = minobs, nboot = nboot,
                                     n0 = n0, nb = nb, pm = pm,
                                     allow_noninteger_counts = allow_noninteger_counts,
                                     passage_times = passage_times,correct_efflux=correct_efflux,
                                     nn_prior = nn_prior,
                                     cohort_transition_prior = cohort_transition_prior,
                                     cohort_transition_patient_id = cohort_transition_patient_id,
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
                                     nn_prior_sd = nn_prior_sd,
                                     nn_prior_sd_floor = nn_prior_sd_floor,
                                     nn_prior_grid_n = nn_prior_grid_n,
                                     nn_prior_fit_subset = nn_prior_fit_subset,
                                     nn_prior_zero_exposure_min = nn_prior_zero_exposure_min,
                                     nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
                                     nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
                                     nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
                                     nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
                                     nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
                                     nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
                                     nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
                                     nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
                                     nn_prior_hybrid_min_obs = nn_prior_hybrid_min_obs,
                                     nn_prior_two_step_support = nn_prior_two_step_support,
                                     nn_prior_two_step_support_min = nn_prior_two_step_support_min,
                                     nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
                                     nn_two_shell_min_delta_n = nn_two_shell_min_delta_n,
                                     nn_two_shell_min_exposure = nn_two_shell_min_exposure,
                                     nn_two_shell_min_observed_count = nn_two_shell_min_observed_count,
                                     nn_two_shell_max_weight_ratio = nn_two_shell_max_weight_ratio,
                                     nn_two_shell_lambda = nn_two_shell_lambda,
                                     nn_two_shell_reuse_sd = nn_two_shell_reuse_sd,
                                     nn_two_shell_uncertainty_floor = nn_two_shell_uncertainty_floor)
  saveRDS(fq_boot, file = file.path(outdir, "bootstrap_res.Rds"))
  if (nn_prior == "empirical_two_shell" && isTRUE(nn_two_shell_save_diagnostics)) {
    saveRDS(
      list(
        replicate = fq_boot$nn_prior_diagnostics,
        node = fq_boot$nn_two_shell_node_diagnostics
      ),
      file = file.path(outdir, "nn_prior_diagnostics.Rds")
    )
  }
  if (nn_prior == "cohort_transition") {
    saveRDS(
      list(
        replicate = fq_boot$nn_prior_diagnostics,
        node = fq_boot$nn_cohort_transition_node_diagnostics
      ),
      file = file.path(outdir, "nn_prior_diagnostics.Rds")
    )
    saveRDS(
      fq_boot$nn_cohort_transition_node_diagnostics,
      file = file.path(outdir, "cohort_transition_patient_diagnostics.Rds")
    )
  }

  landscape_data <- fitKrig(fq_boot, nboot, krig_bootstrap_mode = krig_bootstrap_mode)
  saveRDS(landscape_data$summary_stats, file = file.path(outdir, "landscape.Rds"))
  saveRDS(landscape_data$posterior_samples, file = file.path(outdir, "landscape_posterior_samples.Rds"))

  if (isTRUE(landscape_data_output)) {
    Krig_stable <- list(landscape_data$krig_stable_mean, landscape_data$krig_stable_median)
    names(Krig_stable) <- c("mean", "median")
    saveRDS(Krig_stable, file = file.path(outdir, "landscape_data.Rds"))
  }
  Rxv <- xval(fq_boot, krig_bootstrap_mode = krig_bootstrap_mode)
  saveRDS(Rxv, file = file.path(outdir, "xval.Rds"))

  ##END HERE.

  invisible(Rxv) # Return Rxv invisibly as side-effect saving is primary
}

##########################################
# Helper functions (internal)
##########################################

#' Calculate the single-chromosome transition probability p_ij
#'
#' Validates `i`, `j`, and `beta` before calling the C++ transition kernel.
#'
#' @param i Non-negative integer parent copy number.
#' @param j Non-negative integer daughter copy number.
#' @param beta Finite mis-segregation probability in `[0, 1]`.
#' @return A transition probability in `[0, 1]`.
#' @export
pij <- function(i, j, beta) {
  validate_nonnegative_integer(i, "i")
  validate_nonnegative_integer(j, "j")
  validate_probability(beta, "beta", upper_inclusive = TRUE)
  alfak_cpp_call("pij_cpp", pij_cpp(i, j, beta), context = "pij")
}

#' Convert string like "1.2.3" to numeric vector
#' @keywords internal
#' @noRd
parse_karyotype_ids <- function(ids) {
  if (!is.character(ids) || length(ids) == 0 || any(!nzchar(ids))) {
    stop("Karyotype IDs must be non-empty character strings.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("Karyotype IDs must be unique.", call. = FALSE)
  }
  pieces <- strsplit(ids, "\\.", perl = TRUE)
  lens <- lengths(pieces)
  if (any(lens == 0) || length(unique(lens)) != 1) {
    stop(
      "All karyotype IDs must have the same number of dot-separated components.",
      call. = FALSE
    )
  }
  bad <- vapply(pieces, function(x) any(!grepl("^(0|[1-9][0-9]*)$", x)), logical(1))
  if (any(bad)) {
    stop(
      sprintf("Invalid karyotype ID(s): %s", paste(ids[bad], collapse = ", ")),
      call. = FALSE
    )
  }
  flat_pieces <- unlist(pieces, use.names = FALSE)
  int_max_chr <- as.character(.Machine$integer.max)
  piece_nchar <- nchar(flat_pieces)
  too_large <- piece_nchar > nchar(int_max_chr) |
    (piece_nchar == nchar(int_max_chr) & flat_pieces > int_max_chr)
  if (any(too_large)) {
    offending_ids <- unique(rep(ids, each = lens[1])[too_large])
    stop(
      sprintf(
        "Invalid karyotype ID(s): %s. At least one karyotype component exceeds supported integer range.",
        paste(offending_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  mat <- matrix(as.integer(flat_pieces), nrow = length(ids), byrow = TRUE)
  storage.mode(mat) <- "integer"
  rownames(mat) <- ids
  mat
}

#' Convert string like "1.2.3" to numeric vector
#' @keywords internal
#' @noRd
s2v <- function(s) {
  parsed <- parse_karyotype_ids(as.character(s))
  if (nrow(parsed) == 1) {
    return(as.numeric(parsed[1, ]))
  }
  parsed
}

#' Calculate R-squared
#' @keywords internal
#' @noRd
R2R <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
}

#' Extract the scalar R2R value from xval() output
#' @keywords internal
#' @noRd
extract_xval_r2r <- function(xval_result) {
  if (is.list(xval_result) && !is.null(xval_result$R2R)) {
    r2r_val <- xval_result$R2R
  } else {
    r2r_val <- xval_result
  }
  if (!is.numeric(r2r_val) || length(r2r_val) != 1 || is.nan(r2r_val) || is.infinite(r2r_val)) {
    stop("`xval()` must return a single numeric R2R value, optionally `NA_real_`, or a list containing scalar `R2R`.")
  }
  as.numeric(r2r_val)
}

ALFAK_FEXP_DELTA_TOL <- 1e-8
ALFAK_EFFLUX_VIABILITY_TOL <- 1e-6
ALFAK_NN_PRIOR_SD_FLOOR <- 1e-3
ALFAK_NN_PRIOR_CENSORED_GRID_POINTS <- 81L
ALFAK_NN_TWO_SHELL_UNCERTAINTY_FLOOR <- 0.25
ALFAK_COUNT_INTEGER_TOL <- sqrt(.Machine$double.eps)
ALFAK_KRIG_NSTEP_CV <- 200L
ALFAK_MAX_EXACT_INTEGER <- 2^53 - 1

#' Get or set the ALFA-K run log path
#'
#' C++ acceleration failures are written to this log before the run stops. High
#' level `alfak()` and `alfak_cohort_transition()` runs set the path to
#' `file.path(outdir, "alfak_run.log")`; low-level helper calls use a temporary
#' session log unless `options(alfakR.run_log_path = ...)` is set.
#'
#' @param path Optional path to set for subsequent log entries.
#' @return The current run log path.
#' @export
alfak_run_log_path <- function(path = NULL) {
  if (!is.null(path)) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop("`path` must be a non-empty character scalar.", call. = FALSE)
    }
    options(alfakR.run_log_path = path)
  }
  getOption("alfakR.run_log_path", file.path(tempdir(), "alfakR_run.log"))
}

#' Read the ALFA-K run log
#'
#' @param n Number of trailing lines to return. Use `Inf` for the full log.
#' @param path Optional run log path. Defaults to `alfak_run_log_path()`.
#' @return Character vector of log lines.
#' @export
alfak_read_run_log <- function(n = Inf, path = alfak_run_log_path()) {
  if (!file.exists(path)) {
    return(character(0))
  }
  lines <- readLines(path, warn = FALSE)
  if (is.finite(n)) {
    n <- as.integer(n)
    if (n <= 0L) return(character(0))
    lines <- utils::tail(lines, n)
  }
  lines
}

#' Print the ALFA-K run log
#'
#' @inheritParams alfak_read_run_log
#' @return Invisibly, the printed log lines.
#' @export
alfak_print_run_log <- function(n = 100L, path = alfak_run_log_path()) {
  lines <- alfak_read_run_log(n = n, path = path)
  if (length(lines)) {
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  invisible(lines)
}

alfak_log_event <- function(level = "INFO", component = "alfak", detail, path = alfak_run_log_path()) {
  level <- toupper(as.character(level)[1])
  component <- as.character(component)[1]
  detail <- as.character(detail)[1]
  line <- sprintf(
    "%s [%s] %s: %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    level,
    component,
    detail
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  cat(line, file = path, append = TRUE, sep = "\n")
  if (isTRUE(getOption("alfakR.echo_run_log", TRUE))) {
    base::message(line)
  }
  invisible(line)
}

alfak_cpp_call <- function(kernel, expr, context = NULL) {
  tryCatch(
    force(expr),
    error = function(e) {
      detail <- conditionMessage(e)
      context_text <- if (!is.null(context) && nzchar(context)) paste0(" in ", context) else ""
      alfak_log_event(
        level = "ERROR",
        component = paste0("cpp.", kernel),
        detail = paste0("C++ kernel failed", context_text, ": ", detail)
      )
      stop(
        sprintf("C++ kernel `%s` failed%s: %s", kernel, context_text, detail),
        call. = FALSE
      )
    }
  )
}

#' Validate scalar positive integer input
#' @keywords internal
#' @noRd
validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(sprintf("`%s` must be a single positive integer.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar non-negative integer input
#' @keywords internal
#' @noRd
validate_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x != floor(x)) {
    stop(sprintf("`%s` must be a single non-negative integer.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar positive finite numeric input
#' @keywords internal
#' @noRd
validate_positive_finite <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive finite numeric value.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar non-negative finite numeric input
#' @keywords internal
#' @noRd
validate_nonnegative_finite <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be a single non-negative finite numeric value.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar logical input
#' @keywords internal
#' @noRd
validate_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be a single TRUE/FALSE value.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar probability input
#' @keywords internal
#' @noRd
validate_probability <- function(x, name, upper_inclusive = FALSE) {
  upper_ok <- if (upper_inclusive) isTRUE(x <= 1) else isTRUE(x < 1)
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || !upper_ok) {
    bound <- if (upper_inclusive) "[0, 1]" else "[0, 1)"
    stop(sprintf("`%s` must be a single finite numeric value in %s.", name, bound), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate Kriging bootstrap sampling mode
#' @keywords internal
#' @noRd
validate_krig_bootstrap_mode <- function(mode) {
  match.arg(mode, c("marginal", "joint"))
}

#' Validate that every timepoint has positive sequencing depth
#' @keywords internal
#' @noRd
validate_positive_depth <- function(x) {
  zero_depth <- colSums(x) == 0
  if (any(zero_depth)) {
    depth_names <- colnames(x)
    if (is.null(depth_names)) {
      depth_names <- as.character(seq_len(ncol(x)))
    }
    stop(
      sprintf(
        paste0(
          "Each timepoint must have positive total counts; zero-depth column(s): %s. ",
          "Remove or explicitly impute missing timepoints."
        ),
        paste(depth_names[zero_depth], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Numerically stable softmax
#' @keywords internal
#' @noRd
softmax <- function(z) {
  if (!all(is.finite(z))) {
    stop("Cannot softmax non-finite logits.", call. = FALSE)
  }
  z <- z - max(z)
  ex <- exp(z)
  ex / sum(ex)
}

#' Safe mean that returns NA for all-missing vectors
#' @keywords internal
#' @noRd
mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x)
}

#' Safe median that returns NA for all-missing vectors
#' @keywords internal
#' @noRd
median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  stats::median(x)
}

#' Safe standard deviation that returns NA when fewer than 2 finite values exist
#' @keywords internal
#' @noRd
sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  stats::sd(x)
}

#' Validate nearest-neighbour prior mode
#' @keywords internal
#' @noRd
validate_nn_prior_mode <- function(nn_prior) {
  match.arg(nn_prior, c("empirical_censored", "empirical_censored_weighted", "empirical_two_shell", "cohort_transition", "none", "empirical"))
}

#' Validate weighted nearest-neighbour prior subset mode
#' @keywords internal
#' @noRd
validate_nn_prior_fit_subset <- function(nn_prior_fit_subset) {
  match.arg(nn_prior_fit_subset, c("hybrid", "all"))
}

#' Validate weighted nearest-neighbour 2-step rescue mode
#' @keywords internal
#' @noRd
validate_nn_prior_two_step_support <- function(nn_prior_two_step_support) {
  match.arg(nn_prior_two_step_support, c("none", "rescue"))
}

#' Validate nearest-neighbour prior controls
#' @keywords internal
#' @noRd
validate_nn_prior_controls <- function(nn_prior_sd = NULL,
                                       nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                                       nn_prior_grid_n = ALFAK_NN_PRIOR_CENSORED_GRID_POINTS,
                                       nn_prior_fit_subset = c("hybrid", "all"),
                                       nn_prior_zero_exposure_min = NULL,
                                       nn_prior_zero_exposure_quantile = 0.10,
                                       nn_prior_zero_weight_scale = 0.50,
                                       nn_prior_zero_weight_cap_ratio = NULL,
                                       nn_prior_zero_birth_fallback_weight = NULL,
                                       nn_prior_zero_birth_child_floor = 0.25,
                                       nn_prior_zero_birth_child_shape = 1,
                                       nn_prior_zero_birth_replicate_floor = 0.50,
                                       nn_prior_zero_birth_replicate_shape = 1,
                                       nn_prior_hybrid_min_obs = 3L,
                                       nn_prior_two_step_support = c("none", "rescue"),
                                       nn_prior_two_step_support_min = 0.15,
                                       nn_prior_two_step_cap_floor = 0.30,
                                       nn_two_shell_min_delta_n = 3L,
                                       nn_two_shell_min_exposure = NULL,
                                       nn_two_shell_min_observed_count = 1L,
                                       nn_two_shell_max_weight_ratio = 1.0,
                                       nn_two_shell_lambda = 1.0,
                                       nn_two_shell_reuse_sd = NULL,
                                       nn_two_shell_uncertainty_floor = NULL) {
  nn_prior_fit_subset <- validate_nn_prior_fit_subset(nn_prior_fit_subset)
  nn_prior_two_step_support <- validate_nn_prior_two_step_support(nn_prior_two_step_support)
  if (!is.null(nn_prior_sd)) {
    validate_positive_finite(nn_prior_sd, "nn_prior_sd")
  }
  validate_positive_finite(nn_prior_sd_floor, "nn_prior_sd_floor")
  validate_positive_integer(nn_prior_grid_n, "nn_prior_grid_n")
  if (nn_prior_grid_n < 3) {
    stop("`nn_prior_grid_n` must be at least 3.", call. = FALSE)
  }
  if (!is.null(nn_prior_zero_exposure_min)) {
    validate_nonnegative_finite(nn_prior_zero_exposure_min, "nn_prior_zero_exposure_min")
  }
  validate_probability(nn_prior_zero_exposure_quantile, "nn_prior_zero_exposure_quantile", upper_inclusive = TRUE)
  validate_probability(nn_prior_zero_weight_scale, "nn_prior_zero_weight_scale", upper_inclusive = TRUE)
  if (!is.null(nn_prior_zero_weight_cap_ratio)) {
    validate_nonnegative_finite(nn_prior_zero_weight_cap_ratio, "nn_prior_zero_weight_cap_ratio")
  }
  if (!is.null(nn_prior_zero_birth_fallback_weight)) {
    validate_probability(nn_prior_zero_birth_fallback_weight, "nn_prior_zero_birth_fallback_weight", upper_inclusive = TRUE)
  }
  validate_probability(nn_prior_zero_birth_child_floor, "nn_prior_zero_birth_child_floor", upper_inclusive = TRUE)
  validate_nonnegative_finite(nn_prior_zero_birth_child_shape, "nn_prior_zero_birth_child_shape")
  validate_probability(nn_prior_zero_birth_replicate_floor, "nn_prior_zero_birth_replicate_floor", upper_inclusive = TRUE)
  validate_nonnegative_finite(nn_prior_zero_birth_replicate_shape, "nn_prior_zero_birth_replicate_shape")
  validate_positive_integer(nn_prior_hybrid_min_obs, "nn_prior_hybrid_min_obs")
  validate_probability(nn_prior_two_step_support_min, "nn_prior_two_step_support_min", upper_inclusive = TRUE)
  validate_probability(nn_prior_two_step_cap_floor, "nn_prior_two_step_cap_floor", upper_inclusive = TRUE)
  validate_positive_integer(nn_two_shell_min_delta_n, "nn_two_shell_min_delta_n")
  validate_nonnegative_integer(nn_two_shell_min_observed_count, "nn_two_shell_min_observed_count")
  if (!is.null(nn_two_shell_min_exposure)) {
    validate_nonnegative_finite(nn_two_shell_min_exposure, "nn_two_shell_min_exposure")
  }
  validate_nonnegative_finite(nn_two_shell_max_weight_ratio, "nn_two_shell_max_weight_ratio")
  validate_nonnegative_finite(nn_two_shell_lambda, "nn_two_shell_lambda")
  if (!is.null(nn_two_shell_reuse_sd)) {
    validate_nonnegative_finite(nn_two_shell_reuse_sd, "nn_two_shell_reuse_sd")
  }
  if (!is.null(nn_two_shell_uncertainty_floor)) {
    validate_nonnegative_finite(nn_two_shell_uncertainty_floor, "nn_two_shell_uncertainty_floor")
  }
  invisible(NULL)
}

#' Coerce supported count containers to a base numeric matrix
#' @keywords internal
#' @noRd
coerce_count_matrix <- function(x, allow_noninteger_counts = FALSE) {
  if (is.null(x)) {
    stop("`yi$x`/`data$x` must be a two-dimensional numeric count object.")
  }
  x_dim <- dim(x)
  if (length(x_dim) != 2) {
    if (is.numeric(x) && is.null(x_dim)) {
      stop(
        "`yi$x`/`data$x` must be a two-dimensional numeric count object. Received a vector instead; this usually means a one-column subset dropped matrix dimensions. Use `drop = FALSE` when subsetting count matrices.",
        call. = FALSE
      )
    }
    stop("`yi$x`/`data$x` must be a two-dimensional numeric count object.")
  }

  x_mat <- try(as.matrix(x), silent = TRUE)
  if (inherits(x_mat, "try-error") || !is.matrix(x_mat)) {
    stop("`yi$x`/`data$x` must be coercible to a numeric matrix of karyotype counts.")
  }
  if (!is.numeric(x_mat)) {
    stop("`yi$x`/`data$x` must contain numeric karyotype counts.")
  }
  if (any(!is.finite(x_mat))) {
    stop("`yi$x`/`data$x` must contain only finite count values.")
  }
  if (any(x_mat < 0)) {
    stop("`yi$x`/`data$x` must contain non-negative count values.")
  }
  rounded_mat <- round(x_mat)
  rounding_delta <- abs(x_mat - rounded_mat)
  near_integer <- rounding_delta > 0 & rounding_delta <= ALFAK_COUNT_INTEGER_TOL
  non_integer <- rounding_delta > ALFAK_COUNT_INTEGER_TOL
  if (any(non_integer)) {
    if (!isTRUE(allow_noninteger_counts)) {
      stop("Non-integer values detected in `yi$x`/`data$x`; set `allow_noninteger_counts = TRUE` to round once at entry.", call. = FALSE)
    }
    warning("Non-integer values detected in `yi$x`/`data$x`; rounding to the nearest integer once at entry.")
    x_mat <- rounded_mat
  } else if (any(near_integer)) {
    warning("Integer-like floating-point values detected in `yi$x`/`data$x`; rounding to the nearest integer once at entry.")
    x_mat <- rounded_mat
  }
  x_mat
}

#' Convert a probability vector to K-1 free softmax logits
#' @keywords internal
#' @noRd
free_softmax_logits <- function(prob) {
  if (!is.numeric(prob) || length(prob) == 0 || any(!is.finite(prob)) || any(prob < 0)) {
    stop("`prob` must be a finite non-negative probability vector.", call. = FALSE)
  }
  if (length(prob) == 1) {
    return(numeric(0))
  }
  prob <- prob / sum(prob)
  ref <- prob[length(prob)]
  if (!is.finite(ref) || ref <= 0) {
    stop("The reference probability for free-softmax logits must be positive.", call. = FALSE)
  }
  log(prob[-length(prob)]) - log(ref)
}

#' Expand K-1 free logits into a K-vector of probabilities
#' @keywords internal
#' @noRd
softmax_from_free_logits <- function(logits_free) {
  softmax(c(logits_free, 0))
}

#' Ensure birth-time estimates are finite before neighbour estimation
#' @keywords internal
#' @noRd
sanitize_birth_times <- function(birth_times_est, peak_times, timepoints) {
  mean_risetime <- mean(peak_times - birth_times_est, na.rm = TRUE)
  fallback_used <- FALSE
  fallback_mask <- rep(FALSE, length(birth_times_est))
  if (!is.finite(mean_risetime)) {
    mean_risetime <- 0
    fallback_used <- TRUE
  }

  missing_birth <- !is.finite(birth_times_est)
  if (any(missing_birth)) {
    birth_times_est[missing_birth] <- peak_times[missing_birth] - mean_risetime
    fallback_used <- TRUE
    fallback_mask[missing_birth] <- TRUE
  }

  unresolved <- !is.finite(birth_times_est)
  if (any(unresolved)) {
    safe_fallback <- peak_times
    safe_fallback[!is.finite(safe_fallback)] <- min(timepoints)
    birth_times_est[unresolved] <- safe_fallback[unresolved]
    fallback_used <- TRUE
    fallback_mask[unresolved] <- TRUE
  }

  if (fallback_used) {
    warning("Using finite fallback birth times for nearest-neighbour estimation because root-finding did not return enough finite birth times.")
  }

  if (!is.null(names(birth_times_est))) {
    names(fallback_mask) <- names(birth_times_est)
  }

  list(
    birth_times = birth_times_est,
    fallback_mask = fallback_mask
  )
}

#' Format a short karyotype preview for diagnostics
#' @keywords internal
#' @noRd
format_karyotype_preview <- function(karyotypes, max_show = 5) {
  if (length(karyotypes) == 0) {
    return("<none>")
  }
  shown <- utils::head(karyotypes, max_show)
  preview <- paste(shown, collapse = ", ")
  if (length(karyotypes) > max_show) {
    preview <- paste0(preview, ", ...")
  }
  preview
}

#' Select frequent karyotypes using the documented minobs rule
#' @keywords internal
#' @noRd
get_frequent_karyotypes <- function(x, minobs) {
  x <- coerce_count_matrix(x)
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop("`yi$x`/`data$x` must have non-empty rownames for karyotype IDs.")
  }
  parse_karyotype_ids(rownames(x))
  if (!is.numeric(minobs) || length(minobs) != 1 || !is.finite(minobs) || minobs < 0) {
    stop("`minobs` must be a single non-negative finite numeric value.")
  }
  minobs <- as.numeric(minobs)
  fq <- rownames(x)[rowSums(x) >= minobs]
  if (length(fq) == 0) {
    stop(sprintf(
      "no frequent karyotypes detected for minobs = %s; frequent karyotypes require rowSums(x) >= minobs",
      format(minobs, trim = TRUE)
    ))
  }
  fq
}

#' Resolve and validate the internal time axis
#' @keywords internal
#' @noRd
resolve_time_axis <- function(data, passage_times = NULL) {
  if (!is.list(data) || is.null(data$x)) {
    stop("`yi`/`data` must be a list containing a matrix-like `x` of karyotype counts.")
  }
  data$x <- coerce_count_matrix(data$x)
  if (ncol(data$x) < 2) {
    stop("At least two timepoints are required in `yi$x`/`data$x`.")
  }

  if (is.null(passage_times)) {
    if (is.null(data$dt) || !is.numeric(data$dt) || length(data$dt) != 1 || !is.finite(data$dt)) {
      stop("`yi$dt` must be a single finite numeric value when `passage_times` is NULL.")
    }
    raw_axis <- suppressWarnings(as.numeric(colnames(data$x)))
    if (length(raw_axis) != ncol(data$x) || any(!is.finite(raw_axis))) {
      stop("When `passage_times` is NULL, `colnames(yi$x)` must be numeric and finite so `colnames(yi$x) * yi$dt` defines the time axis.")
    }
    time_axis <- raw_axis * data$dt
  } else {
    if (!is.numeric(passage_times)) {
      stop("`passage_times` must be a numeric vector when supplied.")
    }
    if (length(passage_times) != ncol(data$x)) {
      stop(sprintf(
        "`passage_times` must have length %d to match ncol(yi$x); got %d.",
        ncol(data$x), length(passage_times)
      ))
    }
    if (any(!is.finite(passage_times))) {
      stop("`passage_times` must contain only finite values.")
    }
    time_axis <- as.numeric(passage_times)
  }

  delta_t <- diff(time_axis)
  if (any(delta_t <= 0)) {
    stop("The internal time axis must be strictly increasing.")
  }

  time_axis
}

#' Normalize counts to frequencies by column
#' @keywords internal
#' @noRd
normalize_columns <- function(count_matrix) {
  totals <- colSums(count_matrix)
  normalized <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix),
                       dimnames = dimnames(count_matrix))
  positive_totals <- totals > 0
  if (any(positive_totals)) {
    normalized[, positive_totals] <- sweep(count_matrix[, positive_totals, drop = FALSE], 2, totals[positive_totals], "/")
  }
  normalized
}

#' Validate viability for efflux correction once before bootstrapping
#' @keywords internal
#' @noRd
prepare_efflux_viability <- function(fq_vec, pm, correct_efflux,
                                     viability_tol = ALFAK_EFFLUX_VIABILITY_TOL) {
  if (!is.numeric(pm) || length(pm) != 1 || !is.finite(pm) || pm < 0 || pm >= 1) {
    stop("`pm` must be a single finite numeric value in [0, 1).")
  }
  viability <- setNames(rep(1, nrow(fq_vec)), rownames(fq_vec))
  if (!isTRUE(correct_efflux)) {
    return(viability)
  }

  viability <- setNames(2 * (1 - pm)^rowSums(fq_vec) - 1, rownames(fq_vec))

  non_positive <- viability <= 0
  if (any(non_positive)) {
    affected <- names(viability)[non_positive]
    affected_vals <- viability[non_positive]
    stop(sprintf(
      paste0(
        "correct_efflux viability pre-check failed before bootstrap: pm=%s yields ",
        "%d frequent karyotype(s) with viability <= 0; affected=%s; viability range among affected=[%.6g, %.6g]."
      ),
      format(pm, trim = TRUE),
      length(affected),
      format_karyotype_preview(affected),
      min(affected_vals),
      max(affected_vals)
    ))
  }

  near_zero <- viability < viability_tol
  if (any(near_zero)) {
    affected <- names(viability)[near_zero]
    affected_vals <- viability[near_zero]
    warning(sprintf(
      paste0(
        "correct_efflux viability pre-check before bootstrap: pm=%s yields ",
        "%d frequent karyotype(s) with 0 < viability < %.1e; affected=%s; viability range among affected=[%.6g, %.6g]."
      ),
      format(pm, trim = TRUE),
      length(affected),
      viability_tol,
      format_karyotype_preview(affected),
      min(affected_vals),
      max(affected_vals)
    ))
  }

  viability
}

#' Run optim and surface non-convergence explicitly
#' @keywords internal
#' @noRd
run_optim_checked <- function(par, fn, ..., method = "BFGS", control = NULL, context) {
  opt <- try(stats::optim(par = par, fn = fn, ..., method = method, control = control), silent = TRUE)
  if (inherits(opt, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(opt)))
  }
  if (!all(is.finite(opt$par)) || !is.finite(opt$value)) {
    stop(sprintf("%s returned non-finite parameters or objective values.", context))
  }
  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning(sprintf(
      "%s returned convergence code %d%s",
      context,
      opt$convergence,
      if (!is.null(opt$message) && nzchar(opt$message)) paste0(": ", opt$message) else "."
    ))
  }
  opt
}

#' Run optim and fail loudly on non-convergence
#' @keywords internal
#' @noRd
run_optim_strict_checked <- function(par, fn, ..., method = "BFGS", control = NULL,
                                     lower = NULL, upper = NULL, context) {
  optim_args <- c(
    list(par = par, fn = fn, method = method, control = control),
    list(...)
  )
  if (!is.null(lower)) {
    optim_args$lower <- lower
  }
  if (!is.null(upper)) {
    optim_args$upper <- upper
  }
  opt <- try(do.call(stats::optim, optim_args), silent = TRUE)
  if (inherits(opt, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(opt)))
  }
  if (!all(is.finite(opt$par)) || !is.finite(opt$value)) {
    stop(sprintf("%s returned non-finite parameters or objective values.", context))
  }
  if (!is.null(opt$convergence) && opt$convergence != 0) {
    stop(sprintf(
      "%s failed with convergence code %d%s",
      context,
      opt$convergence,
      if (!is.null(opt$message) && nzchar(opt$message)) paste0(": ", opt$message) else "."
    ))
  }
  opt
}

#' Run optimise and warn on invalid scalar optima
#' @keywords internal
#' @noRd
run_optimise_checked <- function(f, interval, ..., context) {
  opt <- try(stats::optimise(f, interval = interval, ...), silent = TRUE)
  if (inherits(opt, "try-error")) {
    warning(sprintf("%s failed: %s", context, as.character(opt)))
    return(NULL)
  }
  if (!is.finite(opt$minimum) || !is.finite(opt$objective)) {
    warning(sprintf("%s returned a non-finite optimum.", context))
    return(NULL)
  }
  opt
}

#' Run optimise and fail loudly on invalid scalar optima
#' @keywords internal
#' @noRd
run_optimise_strict_checked <- function(f, interval, ..., context) {
  opt <- try(stats::optimise(f, interval = interval, ...), silent = TRUE)
  if (inherits(opt, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(opt)))
  }
  if (!is.finite(opt$minimum) || !is.finite(opt$objective)) {
    stop(sprintf("%s returned a non-finite optimum.", context))
  }
  opt
}

#' Run nlminb and fail loudly on invalid optimizer state
#' @keywords internal
#' @noRd
run_nlminb_strict_checked <- function(start, objective, ..., lower = NULL, upper = NULL,
                                      control = list(), context) {
  opt <- try(
    stats::nlminb(
      start = start,
      objective = objective,
      ...,
      lower = lower,
      upper = upper,
      control = control
    ),
    silent = TRUE
  )
  if (inherits(opt, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(opt)))
  }
  if (!all(is.finite(opt$par)) || !is.finite(opt$objective)) {
    stop(sprintf("%s returned non-finite parameters or objective values.", context))
  }
  if (!is.null(opt$convergence) && opt$convergence != 0) {
    stop(sprintf(
      "%s failed with convergence code %d%s",
      context,
      opt$convergence,
      if (!is.null(opt$message) && nzchar(opt$message)) paste0(": ", opt$message) else "."
    ))
  }
  opt
}

#' Run solve.QP and fail loudly on invalid optimizer state
#' @keywords internal
#' @noRd
run_solve_qp_checked <- function(Dmat, dvec, Amat, bvec, meq, context) {
  qp_sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq), silent = TRUE)
  if (inherits(qp_sol, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(qp_sol)))
  }
  if (is.null(qp_sol$solution) || any(!is.finite(qp_sol$solution))) {
    stop(sprintf("%s returned a non-finite solution.", context))
  }
  qp_sol
}

#' Weighted parent fitness for nearest-neighbour prior
#' @keywords internal
#' @noRd
weighted_parent_fitness <- function(nni_item, fpar) {
  valid_parents <- nni_item$nj[nni_item$nj %in% names(fpar)]
  if (length(valid_parents) == 0) {
    return(NA_real_)
  }
  parent_fitness <- fpar[valid_parents]
  parent_weights <- nni_item$pij[match(valid_parents, nni_item$nj)]
  if (length(parent_weights) == length(parent_fitness) &&
      all(is.finite(parent_weights)) &&
      all(parent_weights >= 0) &&
      sum(parent_weights) > 0) {
    return(stats::weighted.mean(parent_fitness, w = parent_weights))
  }
  mean(parent_fitness, na.rm = TRUE)
}

#' Resolve parent opportunity weights for weighted nearest-neighbour priors
#' @keywords internal
#' @noRd
resolve_nn_parent_opportunity_weights <- function(pij_values, parent_birth_times,
                                                  timepoints, parent_xfit, ntot) {
  n_parents <- length(pij_values)
  if (length(parent_birth_times) != n_parents ||
      !is.matrix(parent_xfit) || nrow(parent_xfit) != n_parents ||
      ncol(parent_xfit) != length(timepoints) || length(ntot) != length(timepoints)) {
    stop("Internal error: malformed parent inputs for weighted nearest-neighbour opportunity weights.")
  }
  if (n_parents == 0) {
    return(numeric(0))
  }
  alfak_cpp_call(
    "alfak_parent_opportunity_weights_cpp",
    alfak_parent_opportunity_weights_cpp(
      pij_values = as.numeric(pij_values),
      parent_birth_times = as.numeric(parent_birth_times),
      timepoints = as.numeric(timepoints),
      parent_xfit = parent_xfit,
      ntot = as.numeric(ntot)
    ),
    context = "resolve_nn_parent_opportunity_weights"
  )
}

#' Exposure-weighted parent fitness for weighted nearest-neighbour priors
#' @keywords internal
#' @noRd
weighted_parent_fitness_exposure <- function(parent_fitness, parent_opportunity_weights,
                                             fallback_mean) {
  n_parents <- length(parent_fitness)
  if (n_parents == 0) {
    return(fallback_mean)
  }
  if (length(parent_opportunity_weights) != n_parents) {
    stop("Internal error: malformed parent inputs for exposure-weighted nearest-neighbour prior centering.")
  }
  alfak_cpp_call(
    "alfak_weighted_parent_mean_cpp",
    alfak_weighted_parent_mean_cpp(
      parent_fitness = as.numeric(parent_fitness),
      weights = as.numeric(parent_opportunity_weights),
      fallback_mean = as.numeric(fallback_mean)
    ),
    context = "weighted_parent_fitness_exposure"
  )
}

#' Compute a child-level birth fallback burden from parent opportunity weights
#' @keywords internal
#' @noRd
compute_nn_child_birth_fallback_burden <- function(parent_birth_fallback,
                                                   parent_opportunity_weights) {
  n_parents <- length(parent_birth_fallback)
  if (n_parents == 0) {
    return(0)
  }
  if (length(parent_opportunity_weights) != n_parents) {
    stop("Internal error: malformed parent opportunity weights for child birth fallback burden.")
  }

  fallback_indicator <- as.numeric(parent_birth_fallback)
  burden <- if (all(is.finite(parent_opportunity_weights)) &&
                all(parent_opportunity_weights >= 0) &&
                sum(parent_opportunity_weights) > 0) {
    stats::weighted.mean(fallback_indicator, w = parent_opportunity_weights)
  } else {
    mean(fallback_indicator, na.rm = TRUE)
  }
  burden <- as.numeric(burden)
  if (!is.finite(burden)) {
    burden <- 1
  }
  pmin(1, pmax(0, burden))
}

#' Compute a replicate-level birth fallback burden from retained zero children
#' @keywords internal
#' @noRd
compute_nn_replicate_birth_fallback_burden <- function(child_burdens, child_exposure) {
  if (!length(child_burdens)) {
    return(0)
  }
  if (length(child_exposure) != length(child_burdens)) {
    stop("Internal error: malformed retained zero-child exposure input for replicate birth fallback burden.")
  }

  bounded_burdens <- pmin(1, pmax(0, as.numeric(child_burdens)))
  bounded_burdens[!is.finite(bounded_burdens)] <- 1
  replicate_burden <- if (all(is.finite(child_exposure)) &&
                          all(child_exposure >= 0) &&
                          sum(child_exposure) > 0) {
    stats::weighted.mean(bounded_burdens, w = child_exposure)
  } else {
    mean(bounded_burdens, na.rm = TRUE)
  }
  if (!is.finite(replicate_burden)) {
    replicate_burden <- 1
  }
  pmin(1, pmax(0, as.numeric(replicate_burden)))
}

#' Convert a fallback burden into a birth-reliability multiplier
#' @keywords internal
#' @noRd
nn_birth_reliability_multiplier <- function(burden, floor, shape) {
  bounded_burden <- pmin(1, pmax(0, as.numeric(burden)))
  bounded_burden[!is.finite(bounded_burden)] <- 1
  multiplier <- floor + (1 - floor) * ((1 - bounded_burden)^shape)
  pmin(1, pmax(floor, multiplier))
}

#' Compute the effective evidence mass contributed by retained zero-only children
#' @keywords internal
#' @noRd
compute_nn_zero_effective_mass <- function(child_exposure, exposure_reference,
                                           child_birth_multiplier) {
  n_children <- length(child_exposure)
  if (n_children == 0) {
    return(numeric(0))
  }
  if (length(child_birth_multiplier) != n_children) {
    stop("Internal error: malformed child birth multiplier input for effective zero evidence mass.")
  }
  if (!is.finite(exposure_reference) || exposure_reference <= 0) {
    return(rep(NA_real_, n_children))
  }

  exposure_term <- pmin(1, as.numeric(child_exposure) / exposure_reference)
  exposure_term[!is.finite(exposure_term)] <- NA_real_
  effective_mass <- exposure_term * as.numeric(child_birth_multiplier)
  effective_mass[!is.finite(effective_mass)] <- NA_real_
  pmax(0, effective_mass)
}

#' Compute a single-step transition weight between two karyotype IDs
#' @keywords internal
#' @noRd
compute_nn_transition_probability <- function(parent_id, child_id, pm) {
  parent_vec <- as.numeric(parse_karyotype_ids(parent_id)[1, ])
  child_vec <- as.numeric(parse_karyotype_ids(child_id)[1, ])
  if (length(parent_vec) != length(child_vec)) {
    stop("Internal error: parent and child karyotypes must have matching dimensions.")
  }
  prod(vapply(seq_along(parent_vec), function(k) {
    pij(parent_vec[k], child_vec[k], pm)
  }, numeric(1)))
}

#' Resolve a safe observed-descendant exposure reference for 2-step rescue
#' @keywords internal
#' @noRd
resolve_nn_two_step_support_reference <- function(descendant_exposure, count_data) {
  ref <- NA_real_
  positive_descendants <- descendant_exposure[is.finite(descendant_exposure) & descendant_exposure > 0]
  if (length(positive_descendants)) {
    ref <- stats::median(positive_descendants, na.rm = TRUE)
  }

  if ((!is.finite(ref) || ref <= 0) && !is.null(count_data) && nrow(count_data) > 0) {
    observed_totals <- rowSums(count_data)
    positive_observed <- observed_totals[is.finite(observed_totals) & observed_totals > 0]
    if (length(positive_observed)) {
      ref <- stats::median(positive_observed, na.rm = TRUE)
    }
  }

  if (!is.finite(ref) || ref <= 0) {
    ref <- 1
  }

  ref
}

#' Compute 2-step observed support for zero-only nearest-neighbour children
#' @keywords internal
#' @noRd
compute_nn_two_step_support <- function(nn_child_contexts, zero_mask, count_data, pm,
                                        exposure_reference, child_birth_multiplier) {
  child_ids <- vapply(nn_child_contexts, function(item) {
    if (!is.null(item$ni) && nzchar(item$ni)) item$ni else NA_character_
  }, character(1))
  if (any(is.na(child_ids)) || any(!nzchar(child_ids))) {
    child_ids <- names(nn_child_contexts)
  }

  zero_indices <- which(zero_mask)
  zero_names <- child_ids[zero_mask]
  if (!length(zero_names) ||
      is.null(count_data) ||
      !is.matrix(count_data) ||
      is.null(rownames(count_data)) ||
      !nrow(count_data)) {
    return(list(
      child_support = setNames(numeric(length(zero_names)), zero_names),
      descendant_exposure_reference = 0,
      n_children_with_support = 0L,
      mean_support = 0,
      median_support = 0,
      max_support = 0
    ))
  }

  zero_items <- nn_child_contexts[zero_mask]
  zero_exposure <- vapply(zero_items, function(item) item$projected_exposure, numeric(1))
  exposure_term <- pmin(1, as.numeric(zero_exposure) / exposure_reference)
  exposure_term[!is.finite(exposure_term)] <- 0
  exposure_term <- pmax(0, exposure_term)
  q_vec <- exposure_term * as.numeric(child_birth_multiplier)
  q_vec[!is.finite(q_vec)] <- 0
  q_vec <- pmax(0, q_vec)
  names(q_vec) <- zero_names

  observed_totals <- rowSums(count_data)
  observed_ids <- names(observed_totals)[is.finite(observed_totals) & observed_totals > 0]
  child_support <- setNames(numeric(length(zero_names)), zero_names)
  if (!length(observed_ids) || !any(q_vec > 0)) {
    return(list(
      child_support = child_support,
      descendant_exposure_reference = 0,
      n_children_with_support = 0L,
      mean_support = 0,
      median_support = 0,
      max_support = 0
    ))
  }

  descendant_scores <- list()
  for (idx in seq_along(zero_indices)) {
    child_idx <- zero_indices[idx]
    child_name <- zero_names[idx]
    if (!is.finite(q_vec[child_name]) || q_vec[child_name] <= 0) {
      next
    }
    descendant_matrix <- gen_all_neighbours(child_name)
    if (!nrow(descendant_matrix)) {
      next
    }
    descendant_ids <- apply(descendant_matrix, 1, paste, collapse = ".")
    descendant_ids <- intersect(descendant_ids, observed_ids)
    descendant_ids <- setdiff(descendant_ids, nn_child_contexts[[child_idx]]$nj)
    if (!length(descendant_ids)) {
      next
    }

    transition_weights <- vapply(descendant_ids, function(desc_id) {
      compute_nn_transition_probability(child_name, desc_id, pm)
    }, numeric(1))
    keep <- is.finite(transition_weights) & transition_weights > 0
    if (!any(keep)) {
      next
    }
    descendant_ids <- descendant_ids[keep]
    transition_weights <- transition_weights[keep]
    for (idx in seq_along(descendant_ids)) {
      desc_id <- descendant_ids[idx]
      descendant_scores[[desc_id]][child_name] <- q_vec[child_name] * transition_weights[idx]
    }
  }

  unique_descendants <- names(descendant_scores)
  if (!length(unique_descendants)) {
    return(list(
      child_support = child_support,
      descendant_exposure_reference = 0,
      n_children_with_support = 0L,
      mean_support = 0,
      median_support = 0,
      max_support = 0
    ))
  }

  descendant_exposure <- observed_totals[unique_descendants]
  descendant_reference <- resolve_nn_two_step_support_reference(
    descendant_exposure = descendant_exposure,
    count_data = count_data
  )
  descendant_term <- pmin(1, as.numeric(descendant_exposure) / descendant_reference)
  descendant_term[!is.finite(descendant_term)] <- 0
  descendant_term <- pmax(0, descendant_term)
  names(descendant_term) <- unique_descendants

  for (desc_id in unique_descendants) {
    desc_scores <- descendant_scores[[desc_id]]
    if (is.null(desc_scores)) {
      next
    }
    desc_scores <- unlist(desc_scores, use.names = TRUE)
    denom <- sum(desc_scores)
    if (!is.finite(denom) || denom <= 0) {
      next
    }
    child_support[names(desc_scores)] <- child_support[names(desc_scores)] +
      descendant_term[desc_id] * (desc_scores / denom)
  }

  child_support <- 1 - exp(-child_support)
  child_support[!is.finite(child_support)] <- 0
  child_support <- pmin(1, pmax(0, child_support))

  list(
    child_support = child_support,
    descendant_exposure_reference = descendant_reference,
    n_children_with_support = as.integer(sum(child_support > sqrt(.Machine$double.eps))),
    mean_support = if (length(child_support)) mean(child_support) else 0,
    median_support = if (length(child_support)) stats::median(child_support) else 0,
    max_support = if (length(child_support)) max(child_support) else 0
  )
}

#' Compute effective zero evidence mass from support and reliability terms
#' @keywords internal
#' @noRd
compute_nn_zero_effective_mass_from_support <- function(base_support_term,
                                                        child_reliability_multiplier) {
  if (length(base_support_term) != length(child_reliability_multiplier)) {
    stop("Internal error: malformed zero-child support or reliability input for effective evidence mass.")
  }
  effective_mass <- as.numeric(base_support_term) * as.numeric(child_reliability_multiplier)
  effective_mass[!is.finite(effective_mass)] <- NA_real_
  pmax(0, effective_mass)
}

#' Build child likelihood surfaces for the censored EB nearest-neighbour prior
#' @keywords internal
#' @noRd
build_nn_prior_child_surfaces <- function(nn_info_items, fpar, build_opt_fc, search_interval,
                                          nn_prior_grid_n = ALFAK_NN_PRIOR_CENSORED_GRID_POINTS,
                                          parent_mean_fn = weighted_parent_fitness,
                                          context = "fit empirical_censored latent-neighbour prior") {
  if (!length(nn_info_items)) {
    stop(sprintf("%s failed: no nearest-neighbour children were available.", context))
  }
  if (length(search_interval) != 2 || any(!is.finite(search_interval)) || diff(search_interval) <= 0) {
    stop(sprintf("%s failed: search_interval must contain two finite increasing bounds.", context))
  }

  child_names <- names(nn_info_items)
  if (is.null(child_names) || !length(child_names)) {
    child_names <- vapply(nn_info_items, function(x) x$ni, character(1))
  }

  parent_means <- vapply(nn_info_items, parent_mean_fn, numeric(1), fpar = fpar)
  valid_children <- is.finite(parent_means)
  if (!any(valid_children)) {
    stop(sprintf("%s failed: no child had a finite weighted parent fitness mean.", context))
  }

  nn_info_items <- nn_info_items[valid_children]
  child_names <- child_names[valid_children]
  parent_means <- parent_means[valid_children]

  validate_positive_integer(nn_prior_grid_n, "nn_prior_grid_n")
  if (nn_prior_grid_n < 3) {
    stop(sprintf("%s failed: `nn_prior_grid_n` must be at least 3.", context))
  }
  grid_n <- as.integer(nn_prior_grid_n)
  fc_grid <- seq(search_interval[1], search_interval[2], length.out = grid_n)
  if (length(fc_grid) < 2 || !all(is.finite(fc_grid))) {
    stop(sprintf("%s failed: could not construct a finite integration grid.", context))
  }
  grid_step <- fc_grid[2] - fc_grid[1]
  if (!is.finite(grid_step) || grid_step <= 0) {
    stop(sprintf("%s failed: integration grid spacing must be positive.", context))
  }
  log_weights <- rep(log(grid_step), length(fc_grid))
  log_weights[c(1, length(fc_grid))] <- log(grid_step / 2)

  loglik_mat <- matrix(NA_real_, nrow = length(nn_info_items), ncol = length(fc_grid),
                       dimnames = list(child_names, NULL))
  map_delta <- rep(NA_real_, length(nn_info_items))
  map_fc <- rep(NA_real_, length(nn_info_items))
  informative_children <- rep(FALSE, length(nn_info_items))
  for (i in seq_along(nn_info_items)) {
    item_i <- nn_info_items[[i]]
    item_timepoints <- if (!is.null(item_i$timepoints)) item_i$timepoints else numeric(0)
    can_use_cpp_surface <- all(c("parent_fitness", "pij", "parent_birth_times", "parent_xfit", "child_obs", "ntot") %in% names(item_i)) &&
      length(item_timepoints) > 0 &&
      is.matrix(item_i$parent_xfit)
    loglik_vals <- if (isTRUE(can_use_cpp_surface)) {
      alfak_cpp_call(
        "alfak_neighbor_loglik_grid_cpp",
        alfak_neighbor_loglik_grid_cpp(
          fc_grid = fc_grid,
          parent_fitness = item_i$parent_fitness,
          pij_values = item_i$pij,
          parent_birth_times = item_i$parent_birth_times,
          timepoints = item_timepoints,
          parent_xfit = item_i$parent_xfit,
          child_obs = item_i$child_obs,
          ntot = item_i$ntot,
          tol = ALFAK_FEXP_DELTA_TOL
        ),
        context = sprintf("build_nn_prior_child_surfaces child=%s", child_names[i])
      )
    } else {
      NULL
    }
    if (is.null(loglik_vals)) {
      objective_fn <- build_opt_fc(item_i, do_prior_param = FALSE)
      objective_vals <- vapply(fc_grid, function(fc_val) {
        val <- try(objective_fn(fc_val), silent = TRUE)
        if (inherits(val, "try-error")) {
          return(NA_real_)
        }
        val
      }, numeric(1))
      loglik_vals <- -objective_vals
    }
    finite_mask <- is.finite(loglik_vals)
    if (!any(finite_mask)) {
      next
    }
    row_max <- max(loglik_vals[finite_mask])
    centered_vals <- rep(-Inf, length(fc_grid))
    centered_vals[finite_mask] <- loglik_vals[finite_mask] - row_max
    row_spread <- diff(range(centered_vals[finite_mask]))
    if (!is.finite(row_spread) || row_spread <= sqrt(.Machine$double.eps)) {
      next
    }
    informative_children[i] <- TRUE
    loglik_mat[i, ] <- centered_vals
    map_fc[i] <- fc_grid[finite_mask][which.max(loglik_vals[finite_mask])]
    map_delta[i] <- map_fc[i] - parent_means[i]
  }

  list(
    nn_info_items = nn_info_items,
    child_names = child_names,
    parent_means = parent_means,
    fc_grid = fc_grid,
    grid_step = grid_step,
    log_weights = log_weights,
    loglik_mat = loglik_mat,
    map_delta = map_delta,
    map_fc = map_fc,
    informative_children = informative_children
  )
}

#' Filter and weight child likelihood surfaces for the censored EB prior fit
#' @keywords internal
#' @noRd
filter_nn_prior_child_surfaces <- function(surface_obj, child_weights = NULL,
                                           context = "fit empirical_censored latent-neighbour prior") {
  n_children <- length(surface_obj$child_names)
  if (is.null(child_weights)) {
    child_weights <- rep(1, n_children)
  }
  if (!is.numeric(child_weights) || length(child_weights) != n_children ||
      any(!is.finite(child_weights)) || any(child_weights < 0)) {
    stop(sprintf("%s failed: child_weights must be a finite non-negative numeric vector aligned with the child surfaces.", context))
  }

  if (!any(surface_obj$informative_children)) {
    stop(sprintf(
      "%s failed: no neighbour children produced an informative finite likelihood surface across the prior grid.",
      context
    ))
  }

  keep_children <- surface_obj$informative_children & child_weights > 0
  if (!any(keep_children)) {
    stop(sprintf("%s failed: no informative neighbour child retained a positive prior-fit weight.", context))
  }

  surface_obj$nn_info_items <- surface_obj$nn_info_items[keep_children]
  surface_obj$child_names <- surface_obj$child_names[keep_children]
  surface_obj$parent_means <- surface_obj$parent_means[keep_children]
  surface_obj$loglik_mat <- surface_obj$loglik_mat[keep_children, , drop = FALSE]
  surface_obj$map_delta <- surface_obj$map_delta[keep_children]
  surface_obj$map_fc <- surface_obj$map_fc[keep_children]
  surface_obj$child_weights <- child_weights[keep_children]
  surface_obj$informative_child_count <- length(surface_obj$child_weights)
  surface_obj
}

#' Fit a single-Gaussian censored EB nearest-neighbour prior
#' @keywords internal
#' @noRd
fit_nn_prior_single_gaussian <- function(surface_obj,
                                         nn_prior_sd = NULL,
                                         nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                                         context = "fit empirical_censored latent-neighbour prior") {
  finite_map_delta <- surface_obj$map_delta[is.finite(surface_obj$map_delta)]
  if (!length(finite_map_delta)) {
    stop(sprintf("%s failed: could not derive finite initial delta estimates.", context))
  }

  mu_init <- stats::median(finite_map_delta, na.rm = TRUE)
  sigma_init <- if (length(finite_map_delta) >= 2) {
    stats::mad(finite_map_delta, center = mu_init, constant = 1, na.rm = TRUE)
  } else {
    NA_real_
  }
  if (!is.finite(sigma_init) || sigma_init <= 0) {
    sigma_init <- stats::sd(finite_map_delta, na.rm = TRUE)
  }
  if (!is.finite(sigma_init) || sigma_init <= 0) {
    sigma_init <- nn_prior_sd_floor
  }
  sigma_init <- max(sigma_init, nn_prior_sd_floor)

  delta_lower <- min(surface_obj$fc_grid) - max(surface_obj$parent_means)
  delta_upper <- max(surface_obj$fc_grid) - min(surface_obj$parent_means)
  delta_span <- delta_upper - delta_lower
  if (!is.finite(delta_span) || delta_span <= 0) {
    delta_span <- max(abs(c(delta_lower, delta_upper)), na.rm = TRUE)
  }
  if (!is.finite(delta_span) || delta_span <= 0) {
    delta_span <- 1
  }
  mu_lower <- delta_lower - delta_span
  mu_upper <- delta_upper + delta_span

  marginal_negloglik <- function(mu, sigma) {
    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
      return(1e9)
    }
    alfak_cpp_call(
      "alfak_nn_prior_marginal_negloglik_cpp",
      alfak_nn_prior_marginal_negloglik_cpp(
        loglik_mat = surface_obj$loglik_mat,
        fc_grid = surface_obj$fc_grid,
        log_weights = surface_obj$log_weights,
        parent_means = surface_obj$parent_means,
        child_weights = surface_obj$child_weights,
        mu = mu,
        sigma = sigma
      ),
      context = context
    )
  }

  if (is.null(nn_prior_sd)) {
    sigma_upper <- max(delta_span * 4, nn_prior_sd_floor * 10)
    if (!is.finite(sigma_upper) || sigma_upper <= nn_prior_sd_floor) {
      sigma_upper <- nn_prior_sd_floor * 10
    }
    opt <- run_nlminb_strict_checked(
      start = c(mu_init, log(sigma_init)),
      objective = function(par) marginal_negloglik(par[1], exp(par[2])),
      lower = c(mu_lower, log(nn_prior_sd_floor)),
      upper = c(mu_upper, log(sigma_upper)),
      control = list(iter.max = 200, eval.max = 400),
      context = context
    )
    prior_mean <- opt$par[1]
    prior_sd <- exp(opt$par[2])
  } else {
    prior_sd <- nn_prior_sd
    opt <- run_optimise_strict_checked(
      function(mu) marginal_negloglik(mu, prior_sd),
      interval = c(mu_lower, mu_upper),
      context = context
    )
    prior_mean <- opt$minimum
  }

  if (!is.finite(prior_mean) || !is.finite(prior_sd) || prior_sd <= 0) {
    stop(sprintf("%s failed: fitted prior hyperparameters were invalid.", context))
  }

  lower_boundary_rate <- mean(surface_obj$map_fc <= min(surface_obj$fc_grid) + surface_obj$grid_step)
  upper_boundary_rate <- mean(surface_obj$map_fc >= max(surface_obj$fc_grid) - surface_obj$grid_step)

  list(
    prior_mean = prior_mean,
    prior_sd = prior_sd,
    n_children = nrow(surface_obj$loglik_mat),
    informative_child_count = surface_obj$informative_child_count,
    sum_child_weight = sum(surface_obj$child_weights),
    map_delta_lower_boundary_rate = lower_boundary_rate,
    map_delta_upper_boundary_rate = upper_boundary_rate
  )
}

#' Estimate an observation-bias corrected latent-neighbour prior
#' @keywords internal
#' @noRd
estimate_nn_prior_censored_eb <- function(nn_info_items, fpar, build_opt_fc, search_interval,
                                          nn_prior_sd = NULL,
                                          nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                                          nn_prior_grid_n = ALFAK_NN_PRIOR_CENSORED_GRID_POINTS,
                                          child_weights = NULL,
                                          parent_mean_fn = weighted_parent_fitness,
                                          context = "fit empirical_censored latent-neighbour prior") {
  surface_obj <- build_nn_prior_child_surfaces(
    nn_info_items = nn_info_items,
    fpar = fpar,
    build_opt_fc = build_opt_fc,
    search_interval = search_interval,
    nn_prior_grid_n = nn_prior_grid_n,
    parent_mean_fn = parent_mean_fn,
    context = context
  )
  surface_obj <- filter_nn_prior_child_surfaces(
    surface_obj = surface_obj,
    child_weights = child_weights,
    context = context
  )
  fit_nn_prior_single_gaussian(
    surface_obj = surface_obj,
    nn_prior_sd = nn_prior_sd,
    nn_prior_sd_floor = nn_prior_sd_floor,
    context = context
  )
}

#' Numerically stable exposure term for neighbour estimation
#' @keywords internal
#' @noRd
fExp_stable <- function(fc_arg, fp_arg, pij_val, tt_arg, tol = ALFAK_FEXP_DELTA_TOL) {
  delta <- fc_arg - fp_arg
  if (abs(delta) < tol) {
    # This is the analytic limit of fp * (exp(tt * delta) - 1) / delta as delta -> 0.
    return(pij_val * fp_arg * tt_arg)
  }
  pij_val * fp_arg * expm1(tt_arg * delta) / delta
}

#' Project a neutral nearest-neighbour child trajectory from parent inputs
#' @keywords internal
#' @noRd
project_nn_child_trajectory <- function(fc_param, parent_fitness, pij_values, parent_birth_times,
                                        timepoints, parent_xfit, tol = ALFAK_FEXP_DELTA_TOL) {
  n_parents <- length(parent_fitness)
  n_time <- length(timepoints)
  if (n_parents == 0) {
    return(rep(0, n_time))
  }
  if (length(pij_values) != n_parents || length(parent_birth_times) != n_parents ||
      !is.matrix(parent_xfit) || nrow(parent_xfit) != n_parents || ncol(parent_xfit) != n_time) {
    stop("Internal error: malformed inputs for projected nearest-neighbour child trajectory.")
  }
  alfak_cpp_call(
    "alfak_nn_project_trajectory_cpp",
    alfak_nn_project_trajectory_cpp(
      fc_param = fc_param,
      parent_fitness = as.numeric(parent_fitness),
      pij_values = as.numeric(pij_values),
      parent_birth_times = as.numeric(parent_birth_times),
      timepoints = as.numeric(timepoints),
      parent_xfit = parent_xfit,
      tol = tol
    ),
    context = "project_nn_child_trajectory"
  )
}

#' Project nearest-neighbour child exposure from the neutral child trajectory
#' @keywords internal
#' @noRd
project_nn_child_exposure <- function(fc_param, parent_fitness, pij_values, parent_birth_times,
                                      timepoints, parent_xfit, ntot,
                                      tol = ALFAK_FEXP_DELTA_TOL) {
  if (length(ntot) != length(timepoints) || any(!is.finite(ntot)) || any(ntot < 0)) {
    stop("Internal error: malformed ntot input for projected nearest-neighbour child exposure.")
  }
  alfak_cpp_call(
    "alfak_nn_project_exposure_cpp",
    alfak_nn_project_exposure_cpp(
      fc_param = fc_param,
      parent_fitness = as.numeric(parent_fitness),
      pij_values = as.numeric(pij_values),
      parent_birth_times = as.numeric(parent_birth_times),
      timepoints = as.numeric(timepoints),
      parent_xfit = parent_xfit,
      ntot = as.numeric(ntot),
      tol = tol
    ),
    context = "project_nn_child_exposure"
  )
}

#' Resolve a safe projected-exposure reference scale
#' @keywords internal
#' @noRd
resolve_nn_exposure_reference <- function(observed_exposure, candidate_exposure, ntot) {
  ref <- NA_real_

  obs_positive <- observed_exposure[is.finite(observed_exposure) & observed_exposure > 0]
  if (length(obs_positive)) {
    ref <- stats::median(obs_positive, na.rm = TRUE)
  }

  if (!is.finite(ref) || ref <= 0) {
    candidate_positive <- candidate_exposure[is.finite(candidate_exposure) & candidate_exposure > 0]
    if (length(candidate_positive)) {
      ref <- stats::median(candidate_positive, na.rm = TRUE)
    }
  }

  if (!is.finite(ref) || ref <= 0) {
    ntot_positive <- ntot[is.finite(ntot) & ntot > 0]
    if (length(ntot_positive)) {
      ref <- stats::median(ntot_positive, na.rm = TRUE)
    }
  }

  if (!is.finite(ref) || ref <= 0) {
    ref <- 1
  }

  ref
}

#' Prepare per-child nearest-neighbour inputs for optimisation and diagnostics
#' @keywords internal
#' @noRd
prepare_nn_child_context <- function(nni_item, boot_data, fpar, birth_times_est,
                                     birth_time_fallback_mask, xfit, timepoints, ntot) {
  valid_parents <- nni_item$nj[nni_item$nj %in% names(fpar)]
  child_obs <- rep(0, length(timepoints))
  if (nni_item$ni %in% rownames(boot_data)) {
    child_obs <- as.numeric(boot_data[nni_item$ni, ])
  }

  if (length(valid_parents) == 0) {
    return(list(
      ni = nni_item$ni,
      nj = character(0),
      pij = numeric(0),
      parent_fitness = numeric(0),
      parent_birth_times = numeric(0),
      parent_birth_fallback = logical(0),
      parent_opportunity_weights = numeric(0),
      parent_xfit = matrix(numeric(0), nrow = 0, ncol = length(timepoints)),
      child_obs = child_obs,
      ntot = as.numeric(ntot),
      timepoints = as.numeric(timepoints),
      parent_fitness_mean_pij = NA_real_,
      parent_fitness_mean_exposure = NA_real_,
      projected_exposure = NA_real_
    ))
  }

  parent_match <- match(valid_parents, nni_item$nj)
  pij_values <- unname(nni_item$pij[parent_match])
  parent_fitness <- unname(fpar[valid_parents])
  parent_birth_times <- unname(birth_times_est[valid_parents])
  parent_birth_fallback <- as.logical(unname(birth_time_fallback_mask[valid_parents]))
  parent_xfit <- xfit[valid_parents, , drop = FALSE]
  parent_opportunity_weights <- resolve_nn_parent_opportunity_weights(
    pij_values = pij_values,
    parent_birth_times = parent_birth_times,
    timepoints = timepoints,
    parent_xfit = parent_xfit,
    ntot = ntot
  )

  parent_fitness_mean_pij <- weighted_parent_fitness(
    list(nj = valid_parents, pij = pij_values),
    fpar = fpar
  )
  parent_fitness_mean_exposure <- weighted_parent_fitness_exposure(
    parent_fitness = parent_fitness,
    parent_opportunity_weights = parent_opportunity_weights,
    fallback_mean = parent_fitness_mean_pij
  )
  projected_exposure <- if (is.finite(parent_fitness_mean_exposure)) {
    project_nn_child_exposure(
      fc_param = parent_fitness_mean_exposure,
      parent_fitness = parent_fitness,
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      ntot = ntot,
      tol = ALFAK_FEXP_DELTA_TOL
    )
  } else {
    NA_real_
  }

  list(
    ni = nni_item$ni,
    nj = valid_parents,
    pij = pij_values,
    parent_fitness = parent_fitness,
    parent_birth_times = parent_birth_times,
    parent_birth_fallback = parent_birth_fallback,
    parent_opportunity_weights = parent_opportunity_weights,
    parent_xfit = parent_xfit,
    child_obs = child_obs,
    ntot = as.numeric(ntot),
    timepoints = as.numeric(timepoints),
    parent_fitness_mean_pij = parent_fitness_mean_pij,
    parent_fitness_mean_exposure = parent_fitness_mean_exposure,
    projected_exposure = projected_exposure
  )
}

#' Build a nearest-neighbour optimisation objective factory for one dataset
#' @keywords internal
#' @noRd
make_nn_child_objective_builder <- function(timepoints, ntot_rounded) {
  function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN,
           do_prior_param = FALSE,
           parent_mean_mode = c("pij", "exposure")) {
    parent_mean_mode <- match.arg(parent_mean_mode)
    if (length(nni_param$nj) == 0) {
      return(function(fc_param) 10^9)
    }
    parent_fitness_mean <- switch(
      parent_mean_mode,
      pij = nni_param$parent_fitness_mean_pij,
      exposure = nni_param$parent_fitness_mean_exposure
    )

    function(fc_param) {
      alfak_cpp_call(
        "alfak_neighbor_objective_cpp",
        alfak_neighbor_objective_cpp(
        fc_param = fc_param,
        parent_fitness = nni_param$parent_fitness,
        pij_values = nni_param$pij,
        parent_birth_times = nni_param$parent_birth_times,
        timepoints = timepoints,
        parent_xfit = nni_param$parent_xfit,
        child_obs = nni_param$child_obs,
        ntot = ntot_rounded,
        parent_fitness_mean = parent_fitness_mean,
        prior_mean = prior_mean_param,
        prior_sd = prior_sd_param,
        do_prior = do_prior_param,
        tol = ALFAK_FEXP_DELTA_TOL
        ),
        context = "make_nn_child_objective_builder"
      )
    }
  }
}

#' Expand a conservative search interval for latent-neighbour optimisation
#' @keywords internal
#' @noRd
expand_nn_fitness_search_interval <- function(fpar) {
  search_interval <- range(fpar, na.rm = TRUE)
  interval_range <- diff(search_interval)
  if (length(search_interval) == 1 || interval_range == 0) {
    interval_range <- abs(search_interval[1] * 0.5) + 1
  }
  search_interval[1] <- search_interval[1] - interval_range
  search_interval[2] <- search_interval[2] + interval_range
  if (search_interval[1] == search_interval[2]) {
    search_interval[1] <- search_interval[1] - 1
    search_interval[2] <- search_interval[2] + 1
  }
  search_interval
}

#' Resolve which nearest-neighbour children are observed in a count matrix
#' @keywords internal
#' @noRd
resolve_nn_present_mask <- function(nn_child_contexts, count_data) {
  nn_present <- names(nn_child_contexts) %in% rownames(count_data)
  if (length(nn_present) > 0 && any(nn_present)) {
    nn_present_indices <- which(nn_present)
    nn_present[nn_present_indices] <- nn_present[nn_present_indices] & vapply(
      names(nn_child_contexts)[nn_present_indices],
      function(child_name) {
        sum(count_data[child_name, ]) > 0
      },
      logical(1)
    )
  }
  nn_present
}

#' Prepare one dataset's neighbour contexts for bootstrap or fallback fitting
#' @keywords internal
#' @noRd
prepare_bootstrap_nn_dataset_state <- function(count_data, current_fq, current_timepoints,
                                               current_epsilon, current_n0, current_nb,
                                               current_viability, current_nn_info,
                                               correct_efflux = FALSE,
                                               context = "bootstrap replicate") {
  x <- normalize_columns(count_data[current_fq, , drop = FALSE])
  dx_dt <- compute_dx_dt(x, current_timepoints)
  x_trim <- x[, -1, drop = FALSE]

  qr_terms <- alfak_cpp_call(
    "alfak_qr_accum_cpp",
    alfak_qr_accum_cpp(x_trim, dx_dt),
    context = context
  )
  Q_accum <- qr_terms$Q_accum
  r_accum <- qr_terms$r_accum
  num_species <- length(current_fq)
  Dmat_boot <- 2 * Q_accum + diag(current_epsilon, num_species)
  dvec_boot <- 2 * r_accum
  A_mat <- matrix(1, nrow = num_species, ncol = 1)
  bvec_val <- 0

  qp_sol <- run_solve_qp_checked(
    Dmat_boot,
    dvec_boot,
    A_mat,
    bvec_val,
    meq = 1,
    context = sprintf("solve.QP %s", context)
  )
  f_qp <- qp_sol$solution

  x0_init <- optimize_initial_frequencies(x, f_qp, current_timepoints)
  opt_res <- joint_optimize(count_data[current_fq, , drop = FALSE], current_timepoints, f_qp, x0_init)

  g0_val <- log(current_nb / current_n0) / diff(current_timepoints)[1]

  if (correct_efflux) {
    viability_vec <- current_viability[current_fq]
    sum_weighted_frel <- sum((opt_res$x0 * opt_res$f) / viability_vec)
    sum_weights <- sum(opt_res$x0 / viability_vec)
    k_const <- (sum_weighted_frel - g0_val) / sum_weights
    opt_res$f <- (opt_res$f - k_const) / viability_vec
  } else {
    opt_res$f <- opt_res$f + g0_val - sum(opt_res$x0 * opt_res$f)
  }

  peak_times <- current_timepoints[apply(x, 1, which.max)]
  birth_times_info <- sanitize_birth_times(
    find_birth_times(opt_res, time_range = c(-1000, max(current_timepoints)), minF = 1 / current_n0),
    peak_times = peak_times,
    timepoints = current_timepoints
  )
  birth_times_est <- birth_times_info$birth_times
  birth_time_fallback_mask <- birth_times_info$fallback_mask

  x0par <- opt_res$x0
  names(x0par) <- current_fq
  fpar <- opt_res$f
  names(fpar) <- current_fq
  names(birth_times_est) <- current_fq
  names(birth_time_fallback_mask) <- current_fq

  xfit <- project_forward_log(x0par, fpar, current_timepoints)
  rownames(xfit) <- current_fq
  ntot_rounded <- round(colSums(count_data))
  nn_child_contexts <- lapply(current_nn_info, function(nni_item) {
    prepare_nn_child_context(
      nni_item = nni_item,
      boot_data = count_data,
      fpar = fpar,
      birth_times_est = birth_times_est,
      birth_time_fallback_mask = birth_time_fallback_mask,
      xfit = xfit,
      timepoints = current_timepoints,
      ntot = ntot_rounded
    )
  })
  names(nn_child_contexts) <- names(current_nn_info)

  list(
    f_initial = f_qp,
    f_final = opt_res$f,
    x0_initial = x0_init,
    x0_final = opt_res$x0,
    fpar = fpar,
    x0par = x0par,
    ntot_rounded = ntot_rounded,
    nn_child_contexts = nn_child_contexts,
    search_interval = expand_nn_fitness_search_interval(fpar)
  )
}

#' Estimate a sample-pooled weighted nearest-neighbour prior from observed children
#' @keywords internal
#' @noRd
estimate_weighted_sample_pooled_prior <- function(nn_child_contexts, nn_present,
                                                  fpar, build_opt_fc, search_interval,
                                                  nn_prior_sd = NULL,
                                                  nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                                                  nn_prior_grid_n = ALFAK_NN_PRIOR_CENSORED_GRID_POINTS,
                                                  ntot,
                                                  context = "fit weighted sample-pooled latent-neighbour prior") {
  observed_mask <- as.logical(nn_present)
  if (!any(observed_mask)) {
    return(list(
      available = FALSE,
      prior_mean = NA_real_,
      prior_sd = NA_real_,
      informative_child_count = NA_integer_,
      sum_child_weight = 0,
      effective_mass_reference = 0,
      exposure_reference = NA_real_
    ))
  }

  observed_items <- nn_child_contexts[observed_mask]
  observed_exposure <- vapply(observed_items, function(item) item$projected_exposure, numeric(1))
  exposure_reference <- resolve_nn_exposure_reference(
    observed_exposure = observed_exposure,
    candidate_exposure = numeric(0),
    ntot = ntot
  )
  observed_effective_mass <- pmin(1, observed_exposure / exposure_reference)
  observed_effective_mass[!is.finite(observed_effective_mass)] <- 0
  observed_effective_mass <- pmax(0, observed_effective_mass)
  effective_mass_reference <- sum(observed_effective_mass)
  if (!is.finite(effective_mass_reference) || effective_mass_reference <= 0) {
    effective_mass_reference <- sum(is.finite(observed_exposure) & observed_exposure > 0)
  }

  prior_fit <- tryCatch(
    estimate_nn_prior_censored_eb(
      nn_info_items = observed_items,
      fpar = fpar,
      build_opt_fc = build_opt_fc,
      search_interval = search_interval,
      nn_prior_sd = nn_prior_sd,
      nn_prior_sd_floor = nn_prior_sd_floor,
      nn_prior_grid_n = nn_prior_grid_n,
      child_weights = rep(1, length(observed_items)),
      parent_mean_fn = function(item, fpar) item$parent_fitness_mean_exposure,
      context = context
    ),
    error = function(e) NULL
  )

  if (is.null(prior_fit) ||
      !is.finite(prior_fit$prior_mean) ||
      !is.finite(prior_fit$prior_sd) ||
      prior_fit$prior_sd <= 0) {
    return(list(
      available = FALSE,
      prior_mean = NA_real_,
      prior_sd = NA_real_,
      informative_child_count = NA_integer_,
      sum_child_weight = length(observed_items),
      effective_mass_reference = effective_mass_reference,
      exposure_reference = exposure_reference
    ))
  }

  list(
    available = TRUE,
    prior_mean = prior_fit$prior_mean,
    prior_sd = prior_fit$prior_sd,
    informative_child_count = prior_fit$informative_child_count,
    sum_child_weight = prior_fit$sum_child_weight,
    effective_mass_reference = effective_mass_reference,
    exposure_reference = exposure_reference,
    map_delta_lower_boundary_rate = prior_fit$map_delta_lower_boundary_rate,
    map_delta_upper_boundary_rate = prior_fit$map_delta_upper_boundary_rate
  )
}

#' Resolve a weak sample-pooled fallback prior for a zero-only weighted replicate
#' @keywords internal
#' @noRd
resolve_weighted_sample_pooled_fallback_prior <- function(sample_pooled_prior,
                                                          weighted_diagnostics) {
  if (is.null(sample_pooled_prior) || !isTRUE(sample_pooled_prior$available)) {
    return(list(
      available = FALSE,
      prior_mean = NA_real_,
      prior_sd = NA_real_,
      alpha = NA_real_
    ))
  }

  zero_effective_mass <- weighted_diagnostics$zero_effective_mass_used
  replicate_birth_multiplier <- weighted_diagnostics$replicate_birth_reliability_multiplier
  reference_mass <- sample_pooled_prior$effective_mass_reference
  alpha <- replicate_birth_multiplier * sqrt(pmin(1, zero_effective_mass / max(reference_mass, 1)))
  if (!is.finite(alpha) || alpha <= sqrt(.Machine$double.eps)) {
    return(list(
      available = FALSE,
      prior_mean = NA_real_,
      prior_sd = NA_real_,
      alpha = NA_real_
    ))
  }

  alpha <- min(1, alpha)
  prior_sd <- sample_pooled_prior$prior_sd / alpha
  if (!is.finite(prior_sd) || prior_sd <= 0) {
    return(list(
      available = FALSE,
      prior_mean = NA_real_,
      prior_sd = NA_real_,
      alpha = NA_real_
    ))
  }

  list(
    available = TRUE,
    prior_mean = sample_pooled_prior$prior_mean,
    prior_sd = prior_sd,
    alpha = alpha
  )
}

#' Build an empty nearest-neighbour prior diagnostics row
#' @keywords internal
#' @noRd
new_nn_prior_diagnostics <- function(nn_prior_mode_requested,
                                     nn_prior_fit_subset_used = NA_character_) {
  list(
    replicate_id = NA_integer_,
    nn_prior_mode_requested = nn_prior_mode_requested,
    nn_prior_mode_used = nn_prior_mode_requested,
    nn_prior_source_used = NA_character_,
    fallback_reason = NA_character_,
    nn_prior_fit_subset_used = nn_prior_fit_subset_used,
    n_frequent_parents = NA_integer_,
    n_1step_candidates = NA_integer_,
    n_1step_observed = NA_integer_,
    n_1step_zero = NA_integer_,
    n_2step_candidates_total = 0L,
    n_2step_candidates_retained = 0L,
    n_2step_observed = 0L,
    n_2step_used_in_backward_term = 0L,
    n_observed_children = 0L,
    n_zero_children_total = 0L,
    n_zero_children_retained = 0L,
    n_zero_children_screened = 0L,
    sum_observed_weight = 0,
    sum_zero_weight_raw = 0,
    sum_zero_weight_final = 0,
    zero_weight_cap_applied = FALSE,
    zero_weight_cap_ratio_used = NA_real_,
    zero_effective_mass_used = 0,
    zero_effective_mass_mean = 0,
    zero_effective_mass_median = 0,
    sum_zero_weight_pre_cap = 0,
    sum_zero_weight_post_cap = 0,
    exposure_threshold_used = NA_real_,
    exposure_reference_used = NA_real_,
    nn_prior_two_step_support_used = NA_character_,
    two_step_support_min_used = 0,
    two_step_cap_floor_used = 0,
    two_step_descendant_exposure_reference_used = 0,
    n_zero_children_with_two_step_support = 0L,
    mean_zero_two_step_support = 0,
    median_zero_two_step_support = 0,
    max_zero_two_step_support = 0,
    n_zero_children_with_birth_fallback = 0L,
    mean_zero_birth_fallback_burden = 0,
    median_zero_birth_fallback_burden = 0,
    replicate_birth_fallback_burden = 0,
    mean_zero_birth_reliability_multiplier = 1,
    median_zero_birth_reliability_multiplier = 1,
    replicate_birth_reliability_multiplier = 1,
    sample_pooled_prior_available = FALSE,
    sample_pooled_prior_mu = NA_real_,
    sample_pooled_prior_sigma = NA_real_,
    sample_pooled_prior_informative_child_count = NA_integer_,
    sample_pooled_alpha_used = NA_real_,
    sample_pooled_sigma_used = NA_real_,
    prior_mu_hat = NA_real_,
    prior_sigma_hat = NA_real_,
    mu01 = NA_real_,
    sigma01 = NA_real_,
    mu12 = NA_real_,
    sigma12 = NA_real_,
    tau_reuse = NA_real_,
    total_inward_weight = 0,
    total_outward_weight = 0,
    adaptive_two_shell_min_exposure = NA_real_,
    informative_child_count = NA_integer_,
    map_delta_lower_boundary_rate = NA_real_,
    map_delta_upper_boundary_rate = NA_real_,
    used_sample_pooled_fallback_for_this_replicate = FALSE,
    used_no_prior_fallback_for_this_replicate = FALSE
  )
}

#' Prepare weighted-mode nearest-neighbour prior inputs and diagnostics
#' @keywords internal
#' @noRd
prepare_weighted_nn_prior_fit <- function(nn_child_contexts, nn_present,
                                          nn_prior_fit_subset = c("hybrid", "all"),
                                          nn_prior_zero_exposure_min = NULL,
                                          nn_prior_zero_exposure_quantile = 0.10,
                                          nn_prior_zero_weight_scale = 0.50,
                                          nn_prior_zero_weight_cap_ratio = NULL,
                                          nn_prior_zero_birth_fallback_weight = NULL,
                                          nn_prior_zero_birth_child_floor = 0.25,
                                          nn_prior_zero_birth_child_shape = 1,
                                          nn_prior_zero_birth_replicate_floor = 0.50,
                                          nn_prior_zero_birth_replicate_shape = 1,
                                          nn_prior_hybrid_min_obs = 3L,
                                          nn_prior_two_step_support = c("none", "rescue"),
                                          nn_prior_two_step_support_min = 0.15,
                                          nn_prior_two_step_cap_floor = 0.30,
                                          count_data = NULL,
                                          pm = 0.00005,
                                          ntot) {
  nn_prior_fit_subset <- validate_nn_prior_fit_subset(nn_prior_fit_subset)
  nn_prior_two_step_support <- validate_nn_prior_two_step_support(nn_prior_two_step_support)
  child_names <- names(nn_child_contexts)
  if (is.null(child_names)) {
    child_names <- vapply(nn_child_contexts, function(item) item$ni, character(1))
  }
  projected_exposure <- vapply(nn_child_contexts, function(item) item$projected_exposure, numeric(1))

  observed_mask <- as.logical(nn_present)
  zero_mask <- !observed_mask
  n_observed <- sum(observed_mask)
  n_zero_total <- sum(zero_mask)
  diagnostics <- list(
    nn_prior_fit_subset_used = nn_prior_fit_subset,
    n_observed_children = as.integer(n_observed),
    n_zero_children_total = as.integer(n_zero_total),
    n_zero_children_retained = 0L,
    n_zero_children_screened = as.integer(n_zero_total),
    sum_observed_weight = as.numeric(n_observed),
    sum_zero_weight_raw = 0,
    sum_zero_weight_final = 0,
    zero_weight_cap_applied = FALSE,
    zero_weight_cap_ratio_used = NA_real_,
    zero_effective_mass_used = 0,
    zero_effective_mass_mean = 0,
    zero_effective_mass_median = 0,
    sum_zero_weight_pre_cap = 0,
    sum_zero_weight_post_cap = 0,
    exposure_threshold_used = NA_real_,
    exposure_reference_used = NA_real_,
    nn_prior_two_step_support_used = nn_prior_two_step_support,
    two_step_support_min_used = nn_prior_two_step_support_min,
    two_step_cap_floor_used = nn_prior_two_step_cap_floor,
    two_step_descendant_exposure_reference_used = 0,
    n_zero_children_with_two_step_support = 0L,
    mean_zero_two_step_support = 0,
    median_zero_two_step_support = 0,
    max_zero_two_step_support = 0,
    n_zero_children_with_birth_fallback = 0L,
    mean_zero_birth_fallback_burden = 0,
    median_zero_birth_fallback_burden = 0,
    replicate_birth_fallback_burden = 0,
    mean_zero_birth_reliability_multiplier = 1,
    median_zero_birth_reliability_multiplier = 1,
    replicate_birth_reliability_multiplier = 1,
    sample_pooled_prior_available = FALSE,
    sample_pooled_prior_mu = NA_real_,
    sample_pooled_prior_sigma = NA_real_,
    sample_pooled_prior_informative_child_count = NA_integer_,
    sample_pooled_alpha_used = NA_real_,
    sample_pooled_sigma_used = NA_real_,
    used_sample_pooled_fallback_for_this_replicate = FALSE,
    used_no_prior_fallback_for_this_replicate = FALSE
  )

  if (!is.null(nn_prior_zero_birth_fallback_weight)) {
    nn_prior_zero_birth_child_floor <- nn_prior_zero_birth_fallback_weight
  }

  exposure_reference <- resolve_nn_exposure_reference(
    observed_exposure = projected_exposure[observed_mask],
    candidate_exposure = projected_exposure[zero_mask],
    ntot = ntot
  )
  diagnostics$exposure_reference_used <- exposure_reference

  zero_child_burden_all <- numeric(sum(zero_mask))
  zero_child_multiplier_all <- numeric(sum(zero_mask))
  zero_birth_fallback_all <- logical(sum(zero_mask))
  zero_two_step_support_all <- numeric(sum(zero_mask))
  zero_two_step_support_full <- numeric(length(nn_child_contexts))
  if (sum(zero_mask) > 0) {
    zero_items_all <- nn_child_contexts[zero_mask]
    zero_child_burden_all <- vapply(zero_items_all, function(item) {
      compute_nn_child_birth_fallback_burden(
        parent_birth_fallback = item$parent_birth_fallback,
        parent_opportunity_weights = item$parent_opportunity_weights
      )
    }, numeric(1))
    zero_child_multiplier_all <- nn_birth_reliability_multiplier(
      burden = zero_child_burden_all,
      floor = nn_prior_zero_birth_child_floor,
      shape = nn_prior_zero_birth_child_shape
    )
    zero_birth_fallback_all <- vapply(zero_items_all, function(item) {
      any(item$parent_birth_fallback)
    }, logical(1))
    diagnostics$n_zero_children_with_birth_fallback <- as.integer(sum(zero_birth_fallback_all))
    diagnostics$mean_zero_birth_fallback_burden <- mean(zero_child_burden_all)
    diagnostics$median_zero_birth_fallback_burden <- stats::median(zero_child_burden_all)

    if (nn_prior_two_step_support == "rescue") {
      two_step_support <- compute_nn_two_step_support(
        nn_child_contexts = nn_child_contexts,
        zero_mask = zero_mask,
        count_data = count_data,
        pm = pm,
        exposure_reference = exposure_reference,
        child_birth_multiplier = zero_child_multiplier_all
      )
      zero_two_step_support_all <- two_step_support$child_support
      zero_two_step_support_full[zero_mask] <- zero_two_step_support_all
      diagnostics$two_step_descendant_exposure_reference_used <- two_step_support$descendant_exposure_reference
      diagnostics$n_zero_children_with_two_step_support <- two_step_support$n_children_with_support
      diagnostics$mean_zero_two_step_support <- two_step_support$mean_support
      diagnostics$median_zero_two_step_support <- two_step_support$median_support
      diagnostics$max_zero_two_step_support <- two_step_support$max_support
    }
  }

  if (nn_prior_fit_subset == "hybrid") {
    if (!is.null(nn_prior_zero_exposure_min)) {
      exposure_threshold <- nn_prior_zero_exposure_min
    } else if (n_observed >= nn_prior_hybrid_min_obs) {
      observed_exposure <- projected_exposure[observed_mask & is.finite(projected_exposure)]
      if (length(observed_exposure)) {
        exposure_threshold <- as.numeric(stats::quantile(
          observed_exposure,
          probs = nn_prior_zero_exposure_quantile,
          na.rm = TRUE,
          names = FALSE,
          type = 7
        ))
      } else {
        exposure_threshold <- 0
      }
    } else {
      exposure_threshold <- 0
    }
    zero_retained_mask <- zero_mask & (
      (is.finite(projected_exposure) & projected_exposure >= exposure_threshold) |
        (nn_prior_two_step_support == "rescue" & zero_two_step_support_full >= nn_prior_two_step_support_min)
    )
    diagnostics$exposure_threshold_used <- exposure_threshold
  } else {
    zero_retained_mask <- zero_mask & (
      is.finite(projected_exposure) |
        (nn_prior_two_step_support == "rescue" & zero_two_step_support_full > 0)
    )
  }

  diagnostics$n_zero_children_retained <- as.integer(sum(zero_retained_mask))
  diagnostics$n_zero_children_screened <- as.integer(n_zero_total - sum(zero_retained_mask))

  zero_weights_raw <- numeric(sum(zero_retained_mask))
  child_birth_burden <- numeric(sum(zero_retained_mask))
  child_birth_multiplier <- numeric(sum(zero_retained_mask))
  child_two_step_support <- numeric(sum(zero_retained_mask))
  zero_support_term <- numeric(sum(zero_retained_mask))
  child_reliability_multiplier <- numeric(sum(zero_retained_mask))
  zero_effective_mass <- numeric(sum(zero_retained_mask))
  replicate_birth_burden <- 0
  replicate_birth_multiplier <- 1
  if (length(zero_weights_raw)) {
    retained_zero_items <- nn_child_contexts[zero_retained_mask]
    retained_zero_exposure <- projected_exposure[zero_retained_mask]
    child_birth_burden <- zero_child_burden_all[zero_retained_mask[zero_mask]]
    child_birth_multiplier <- zero_child_multiplier_all[zero_retained_mask[zero_mask]]
    child_two_step_support <- zero_two_step_support_all[zero_retained_mask[zero_mask]]
    replicate_birth_burden <- compute_nn_replicate_birth_fallback_burden(
      child_burdens = child_birth_burden,
      child_exposure = retained_zero_exposure
    )
    replicate_birth_multiplier <- nn_birth_reliability_multiplier(
      burden = replicate_birth_burden,
      floor = nn_prior_zero_birth_replicate_floor,
      shape = nn_prior_zero_birth_replicate_shape
    )
    exposure_term <- pmin(1, retained_zero_exposure / exposure_reference)
    exposure_term[!is.finite(exposure_term)] <- 0
    exposure_term <- pmax(0, exposure_term)
    zero_support_term <- if (nn_prior_two_step_support == "rescue") {
      pmax(exposure_term, nn_prior_two_step_cap_floor * child_two_step_support)
    } else {
      exposure_term
    }
    zero_support_term[!is.finite(zero_support_term)] <- 0
    zero_support_term <- pmax(0, zero_support_term)
    child_reliability_multiplier <- if (nn_prior_two_step_support == "rescue") {
      child_birth_multiplier + (1 - child_birth_multiplier) * child_two_step_support
    } else {
      child_birth_multiplier
    }
    child_reliability_multiplier[!is.finite(child_reliability_multiplier)] <- 0
    child_reliability_multiplier <- pmin(1, pmax(0, child_reliability_multiplier))
    zero_effective_mass <- compute_nn_zero_effective_mass_from_support(
      base_support_term = zero_support_term,
      child_reliability_multiplier = child_reliability_multiplier
    )
    birth_reliability_multiplier <- child_reliability_multiplier * replicate_birth_multiplier
    zero_weights_raw <- nn_prior_zero_weight_scale *
      zero_support_term *
      birth_reliability_multiplier
    zero_weights_raw[!is.finite(zero_weights_raw)] <- 0
    zero_weights_raw <- pmax(0, zero_weights_raw)

    diagnostics$replicate_birth_fallback_burden <- replicate_birth_burden
    diagnostics$mean_zero_birth_reliability_multiplier <- mean(birth_reliability_multiplier)
    diagnostics$median_zero_birth_reliability_multiplier <- stats::median(birth_reliability_multiplier)
    diagnostics$replicate_birth_reliability_multiplier <- replicate_birth_multiplier
  }

  diagnostics$sum_zero_weight_raw <- sum(zero_weights_raw)
  diagnostics$sum_zero_weight_pre_cap <- sum(zero_weights_raw)
  finite_zero_effective_mass <- zero_effective_mass[is.finite(zero_effective_mass)]
  zero_effective_mass_total <- if (length(finite_zero_effective_mass)) sum(finite_zero_effective_mass) else NA_real_
  if (length(finite_zero_effective_mass)) {
    diagnostics$zero_effective_mass_used <- zero_effective_mass_total
    diagnostics$zero_effective_mass_mean <- mean(finite_zero_effective_mass)
    diagnostics$zero_effective_mass_median <- stats::median(finite_zero_effective_mass)
  }

  cap_ratio <- if (n_observed == 0L) {
    1
  } else if (!is.null(nn_prior_zero_weight_cap_ratio)) {
    nn_prior_zero_weight_cap_ratio
  } else {
    if (!is.finite(zero_effective_mass_total) || zero_effective_mass_total <= 0) {
      1
    } else {
      min(1, sqrt(n_observed / max(zero_effective_mass_total, 1)))
    }
  }
  diagnostics$zero_weight_cap_ratio_used <- cap_ratio
  zero_weight_max <- cap_ratio * n_observed
  zero_weights_final <- zero_weights_raw
  zero_weight_pre_cap <- sum(zero_weights_raw)
  if (n_observed > 0 &&
      zero_weight_pre_cap > 0 &&
      zero_weight_pre_cap > zero_weight_max * (1 + sqrt(.Machine$double.eps))) {
    zero_weights_final <- zero_weights_raw * (zero_weight_max / zero_weight_pre_cap)
    diagnostics$zero_weight_cap_applied <- TRUE
  }
  diagnostics$sum_zero_weight_final <- sum(zero_weights_final)
  diagnostics$sum_zero_weight_post_cap <- sum(zero_weights_final)

  full_weights <- numeric(length(nn_child_contexts))
  full_weights[observed_mask] <- 1
  if (length(zero_weights_final)) {
    full_weights[zero_retained_mask] <- zero_weights_final
  }
  prior_keep <- observed_mask | (zero_retained_mask & full_weights > 0)

  list(
    prior_items = nn_child_contexts[prior_keep],
    child_weights = full_weights[prior_keep],
    parent_mean_fn = function(item, fpar) item$parent_fitness_mean_exposure,
    can_fit_replicate_prior = n_observed > 0L,
    diagnostics = diagnostics
  )
}

#' Normalize finite non-negative weights with a conservative fallback
#' @keywords internal
#' @noRd
normalize_nn_weights <- function(weights, fallback_n = length(weights)) {
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  if (length(weights) && sum(weights) > 0) {
    return(weights / sum(weights))
  }
  if (!length(weights) && fallback_n <= 0) {
    return(numeric(0))
  }
  rep(1 / fallback_n, fallback_n)
}

#' Fit a weighted Gaussian shell-delta prior with robust floors
#' @keywords internal
#' @noRd
fit_shell_delta_prior <- function(delta, weights = NULL,
                                  sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                                  fallback_mu = 0,
                                  fallback_sd = sd_floor) {
  validate_positive_finite(sd_floor, "sd_floor")
  delta <- as.numeric(delta)
  if (is.null(weights)) {
    weights <- rep(1, length(delta))
  }
  weights <- as.numeric(weights)
  keep <- is.finite(delta) & is.finite(weights) & weights > 0
  if (!any(keep)) {
    return(list(mu = fallback_mu, sigma = max(fallback_sd, sd_floor), n = 0L, sum_weight = 0))
  }
  delta <- delta[keep]
  weights <- weights[keep]
  mu <- stats::weighted.mean(delta, weights)
  if (!is.finite(mu)) {
    mu <- fallback_mu
  }
  w_sum <- sum(weights)
  sigma <- NA_real_
  if (length(delta) >= 2 && is.finite(w_sum) && w_sum > 0) {
    var_w <- sum(weights * (delta - mu)^2) / w_sum
    sigma <- sqrt(max(0, var_w))
  }
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- fallback_sd
  }
  list(mu = mu, sigma = max(sigma, sd_floor), n = length(delta), sum_weight = w_sum)
}

#' Resolve the two-shell provisional-fitness uncertainty floor
#' @keywords internal
#' @noRd
resolve_nn_two_shell_uncertainty_floor <- function(nn_two_shell_uncertainty_floor,
                                                   nn_prior_sd_floor) {
  floor_val <- if (is.null(nn_two_shell_uncertainty_floor)) {
    ALFAK_NN_TWO_SHELL_UNCERTAINTY_FLOOR
  } else {
    nn_two_shell_uncertainty_floor
  }
  max(floor_val, nn_prior_sd_floor)
}

#' Estimate a scalar optimizer standard error from local objective curvature
#' @keywords internal
#' @noRd
estimate_scalar_objective_se <- function(objective_fn, optimum, search_interval,
                                         se_floor) {
  if (!is.finite(optimum) || length(search_interval) != 2 ||
      any(!is.finite(search_interval)) || diff(search_interval) <= 0) {
    return(se_floor)
  }
  h <- max(1e-4, diff(search_interval) * 1e-4)
  left <- max(search_interval[1], optimum - h)
  right <- min(search_interval[2], optimum + h)
  if (!(left < optimum && optimum < right)) {
    return(max(se_floor, diff(search_interval) / 4))
  }
  f0 <- try(objective_fn(optimum), silent = TRUE)
  fl <- try(objective_fn(left), silent = TRUE)
  fr <- try(objective_fn(right), silent = TRUE)
  if (inherits(f0, "try-error") || inherits(fl, "try-error") || inherits(fr, "try-error") ||
      !is.finite(f0) || !is.finite(fl) || !is.finite(fr)) {
    return(max(se_floor, diff(search_interval) / 4))
  }
  h_left <- optimum - left
  h_right <- right - optimum
  curvature <- 2 * ((fl - f0) / h_left + (fr - f0) / h_right) / (h_left + h_right)
  if (!is.finite(curvature) || curvature <= sqrt(.Machine$double.eps)) {
    return(max(se_floor, diff(search_interval) / 4))
  }
  max(se_floor, sqrt(1 / curvature))
}

#' Build projected one-step anchor trajectories for two-shell fitting
#' @keywords internal
#' @noRd
build_one_step_anchor_states <- function(nn_child_contexts, f1_hat, timepoints,
                                         min_frequency) {
  anchor_names <- names(nn_child_contexts)
  states <- vector("list", length(anchor_names))
  names(states) <- anchor_names
  for (child_name in anchor_names) {
    item <- nn_child_contexts[[child_name]]
    f1 <- unname(f1_hat[child_name])
    if (!is.finite(f1) || length(item$nj) == 0) {
      next
    }
    trajectory <- project_nn_child_trajectory(
      fc_param = f1,
      parent_fitness = item$parent_fitness,
      pij_values = item$pij,
      parent_birth_times = item$parent_birth_times,
      timepoints = timepoints,
      parent_xfit = item$parent_xfit,
      tol = ALFAK_FEXP_DELTA_TOL
    )
    observed_idx <- which(item$child_obs > 0)
    if (length(observed_idx)) {
      birth_time <- min(timepoints[observed_idx])
    } else {
      above_floor <- which(trajectory >= min_frequency)
      birth_time <- if (length(above_floor)) timepoints[min(above_floor)] else min(timepoints)
    }
    states[[child_name]] <- list(
      karyotype = child_name,
      fitness = f1,
      xfit = trajectory,
      birth_time = birth_time,
      exposure = sum(item$ntot * trajectory)
    )
  }
  states[!vapply(states, is.null, logical(1))]
}

#' Resolve an adaptive two-shell exposure threshold
#' @keywords internal
#' @noRd
resolve_two_shell_min_exposure <- function(nn_child_contexts, nn_present,
                                           candidate_exposure,
                                           nn_two_shell_min_exposure) {
  if (!is.null(nn_two_shell_min_exposure)) {
    return(as.numeric(nn_two_shell_min_exposure))
  }
  one_step_observed_exposure <- vapply(nn_child_contexts[nn_present], function(item) {
    item$projected_exposure
  }, numeric(1))
  reference <- c(
    one_step_observed_exposure[is.finite(one_step_observed_exposure) & one_step_observed_exposure > 0],
    candidate_exposure[is.finite(candidate_exposure) & candidate_exposure > 0]
  )
  if (!length(reference)) {
    return(0)
  }
  as.numeric(stats::quantile(reference, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7))
}

#' Build supported two-step candidate paths from provisional one-step anchors
#' @keywords internal
#' @noRd
build_two_shell_candidates <- function(nn_child_contexts, f1_hat, boot_data, fpar,
                                       timepoints, ntot, pm, nn_present,
                                       nn_two_shell_min_exposure = NULL,
                                       nn_two_shell_min_observed_count = 1L,
                                       min_frequency = 0) {
  anchor_states <- build_one_step_anchor_states(
    nn_child_contexts = nn_child_contexts,
    f1_hat = f1_hat,
    timepoints = timepoints,
    min_frequency = min_frequency
  )
  if (!length(anchor_states)) {
    return(list(
      candidates = data.frame(),
      retained = data.frame(),
      anchor_states = anchor_states,
      adaptive_min_exposure = NA_real_
    ))
  }
  rows <- list()
  row_idx <- 0L
  frequent_ids <- names(fpar)
  observed_totals <- if (!is.null(boot_data) && is.matrix(boot_data)) rowSums(boot_data) else numeric(0)
  for (child_name in names(anchor_states)) {
    descendants <- gen_all_neighbours(child_name)
    if (!nrow(descendants)) {
      next
    }
    descendant_ids <- apply(descendants, 1, paste, collapse = ".")
    descendant_ids <- setdiff(unique(descendant_ids), c(child_name, frequent_ids))
    if (!length(descendant_ids)) {
      next
    }
    anchor <- anchor_states[[child_name]]
    parent_exposure <- sum(ntot * anchor$xfit * as.numeric(timepoints >= anchor$birth_time))
    for (desc_id in descendant_ids) {
      q <- compute_nn_transition_probability(child_name, desc_id, pm)
      if (!is.finite(q) || q <= 0) {
        next
      }
      expected_exposure_path <- project_nn_child_exposure(
        fc_param = anchor$fitness,
        parent_fitness = anchor$fitness,
        pij_values = q,
        parent_birth_times = anchor$birth_time,
        timepoints = timepoints,
        parent_xfit = matrix(anchor$xfit, nrow = 1),
        ntot = ntot,
        tol = ALFAK_FEXP_DELTA_TOL
      )
      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        one_step = child_name,
        descendant = desc_id,
        transition_probability = q,
        parent_anchor_exposure = parent_exposure,
        expected_exposure_path = expected_exposure_path,
        observed_count = if (desc_id %in% names(observed_totals)) as.numeric(observed_totals[desc_id]) else 0,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(list(
      candidates = data.frame(),
      retained = data.frame(),
      anchor_states = anchor_states,
      adaptive_min_exposure = NA_real_
    ))
  }
  candidates <- do.call(rbind, rows)
  total_expected <- stats::ave(
    candidates$expected_exposure_path,
    candidates$descendant,
    FUN = function(x) sum(x[is.finite(x)], na.rm = TRUE)
  )
  candidates$expected_exposure <- as.numeric(total_expected)
  adaptive_min_exposure <- resolve_two_shell_min_exposure(
    nn_child_contexts = nn_child_contexts,
    nn_present = nn_present,
    candidate_exposure = candidates$expected_exposure,
    nn_two_shell_min_exposure = nn_two_shell_min_exposure
  )
  retained <- candidates[
    candidates$observed_count >= nn_two_shell_min_observed_count |
      (is.finite(candidates$expected_exposure) & candidates$expected_exposure >= adaptive_min_exposure),
    ,
    drop = FALSE
  ]
  rownames(candidates) <- NULL
  rownames(retained) <- NULL
  list(
    candidates = candidates,
    retained = retained,
    anchor_states = anchor_states,
    adaptive_min_exposure = adaptive_min_exposure
  )
}

#' Compute path responsibilities for retained two-shell candidate paths
#' @keywords internal
#' @noRd
compute_two_shell_path_responsibilities <- function(candidate_paths) {
  if (!is.data.frame(candidate_paths) || !nrow(candidate_paths)) {
    return(candidate_paths)
  }
  cpp_resp <- alfak_cpp_call(
    "alfak_two_shell_path_responsibilities_cpp",
    alfak_two_shell_path_responsibilities_cpp(
      descendant = as.character(candidate_paths$descendant),
      parent_anchor_exposure = as.numeric(candidate_paths$parent_anchor_exposure),
      transition_probability = as.numeric(candidate_paths$transition_probability)
    ),
    context = "compute_two_shell_path_responsibilities"
  )
  if (is.list(cpp_resp) &&
      length(cpp_resp$path_supply) == nrow(candidate_paths) &&
      length(cpp_resp$path_responsibility) == nrow(candidate_paths)) {
    candidate_paths$path_supply <- as.numeric(cpp_resp$path_supply)
    candidate_paths$path_responsibility <- as.numeric(cpp_resp$path_responsibility)
    return(candidate_paths)
  }
  alfak_log_event(
    level = "ERROR",
    component = "cpp.alfak_two_shell_path_responsibilities_cpp",
    detail = "C++ kernel returned malformed output in compute_two_shell_path_responsibilities."
  )
  stop("C++ kernel `alfak_two_shell_path_responsibilities_cpp` returned malformed output.", call. = FALSE)
}

#' Estimate provisional two-step descendant fitness and uncertainty
#' @keywords internal
#' @noRd
estimate_provisional_two_step_fitness <- function(candidate_paths, anchor_states, boot_data,
                                                  timepoints, ntot, search_interval,
                                                  uncertainty_floor) {
  if (!is.data.frame(candidate_paths) || !nrow(candidate_paths)) {
    return(data.frame(
      karyotype = character(0),
      f2_hat = numeric(0),
      f2_se = numeric(0),
      f2_var = numeric(0),
      observed_count = numeric(0),
      expected_exposure = numeric(0),
      fit_failed = logical(0),
      failure_reason = character(0),
      stringsAsFactors = FALSE
    ))
  }
  descendant_ids <- unique(candidate_paths$descendant)
  out <- vector("list", length(descendant_ids))
  names(out) <- descendant_ids
  observed_totals <- if (!is.null(boot_data) && is.matrix(boot_data)) rowSums(boot_data) else numeric(0)
  for (desc_id in descendant_ids) {
    paths <- candidate_paths[candidate_paths$descendant == desc_id, , drop = FALSE]
    parent_names <- unique(paths$one_step)
    parent_names <- parent_names[parent_names %in% names(anchor_states)]
    if (!length(parent_names)) {
      out[[desc_id]] <- list(f2_hat = NA_real_, f2_se = uncertainty_floor, fit_failed = TRUE,
                             failure_reason = "no_valid_one_step_anchor")
      next
    }
    parent_fitness <- vapply(parent_names, function(x) anchor_states[[x]]$fitness, numeric(1))
    parent_birth_times <- vapply(parent_names, function(x) anchor_states[[x]]$birth_time, numeric(1))
    parent_xfit <- do.call(rbind, lapply(parent_names, function(x) anchor_states[[x]]$xfit))
    rownames(parent_xfit) <- parent_names
    pij_values <- vapply(parent_names, function(parent_name) {
      parent_paths <- paths[paths$one_step == parent_name, , drop = FALSE]
      sum(parent_paths$transition_probability, na.rm = TRUE)
    }, numeric(1))
    child_obs <- rep(0, length(timepoints))
    if (desc_id %in% rownames(boot_data)) {
      child_obs <- as.numeric(boot_data[desc_id, ])
    }
    opportunity_weights <- resolve_nn_parent_opportunity_weights(
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      ntot = ntot
    )
    parent_mean <- weighted_parent_fitness_exposure(
      parent_fitness = parent_fitness,
      parent_opportunity_weights = opportunity_weights,
      fallback_mean = stats::weighted.mean(parent_fitness, w = normalize_nn_weights(pij_values))
    )
    objective_fn <- function(fc_param) {
      alfak_cpp_call(
        "alfak_neighbor_objective_cpp",
        alfak_neighbor_objective_cpp(
        fc_param = fc_param,
        parent_fitness = parent_fitness,
        pij_values = pij_values,
        parent_birth_times = parent_birth_times,
        timepoints = timepoints,
        parent_xfit = parent_xfit,
        child_obs = child_obs,
        ntot = ntot,
        parent_fitness_mean = parent_mean,
        prior_mean = NaN,
        prior_sd = NaN,
        do_prior = FALSE,
        tol = ALFAK_FEXP_DELTA_TOL
        ),
        context = sprintf("estimate provisional two-step fitness for descendant %s", desc_id)
      )
    }
    res <- run_optimise_checked(
      objective_fn,
      interval = search_interval,
      context = sprintf("optimise provisional two-step fitness for descendant %s", desc_id)
    )
    if (is.null(res)) {
      out[[desc_id]] <- list(f2_hat = NA_real_, f2_se = uncertainty_floor, fit_failed = TRUE,
                             failure_reason = "optimise_failed")
      next
    }
    se <- estimate_scalar_objective_se(
      objective_fn = objective_fn,
      optimum = res$minimum,
      search_interval = search_interval,
      se_floor = uncertainty_floor
    )
    out[[desc_id]] <- list(f2_hat = res$minimum, f2_se = se, fit_failed = FALSE,
                           failure_reason = NA_character_)
  }
  res_df <- data.frame(
    karyotype = descendant_ids,
    f2_hat = vapply(out, `[[`, numeric(1), "f2_hat"),
    f2_se = vapply(out, `[[`, numeric(1), "f2_se"),
    observed_count = vapply(descendant_ids, function(desc_id) {
      if (desc_id %in% names(observed_totals)) as.numeric(observed_totals[desc_id]) else 0
    }, numeric(1)),
    expected_exposure = vapply(descendant_ids, function(desc_id) {
      paths <- candidate_paths[candidate_paths$descendant == desc_id, , drop = FALSE]
      unique(paths$expected_exposure)[1]
    }, numeric(1)),
    fit_failed = vapply(out, `[[`, logical(1), "fit_failed"),
    failure_reason = vapply(out, `[[`, character(1), "failure_reason"),
    stringsAsFactors = FALSE
  )
  res_df$f2_var <- res_df$f2_se^2
  res_df
}

#' Compute capped outward prior path weights
#' @keywords internal
#' @noRd
compute_two_shell_outward_weights <- function(candidate_paths, f2_fit, nn_child_contexts,
                                              sigma12, tau_reuse,
                                              nn_two_shell_max_weight_ratio) {
  if (!is.data.frame(candidate_paths) || !nrow(candidate_paths) || !nrow(f2_fit)) {
    return(candidate_paths[FALSE, , drop = FALSE])
  }
  f2_map <- f2_fit[is.finite(f2_fit$f2_hat) & !f2_fit$fit_failed, , drop = FALSE]
  paths <- candidate_paths[candidate_paths$descendant %in% f2_map$karyotype, , drop = FALSE]
  if (!nrow(paths)) {
    return(paths)
  }
  support_reference <- stats::median(
    (f2_map$observed_count + f2_map$expected_exposure)[
      is.finite(f2_map$observed_count + f2_map$expected_exposure) &
        (f2_map$observed_count + f2_map$expected_exposure) > 0
    ],
    na.rm = TRUE
  )
  if (!is.finite(support_reference) || support_reference <= 0) {
    support_reference <- 1
  }
  fit_idx <- match(paths$descendant, f2_map$karyotype)
  support_mass <- f2_map$observed_count[fit_idx] + f2_map$expected_exposure[fit_idx]
  support_weight <- pmin(1, support_mass / support_reference)
  support_weight[!is.finite(support_weight) | support_weight < 0] <- 0
  uncertainty_weight <- sigma12^2 / (sigma12^2 + f2_map$f2_var[fit_idx] + tau_reuse^2)
  uncertainty_weight[!is.finite(uncertainty_weight) | uncertainty_weight < 0] <- 0
  paths$f2_hat <- f2_map$f2_hat[fit_idx]
  paths$f2_var <- f2_map$f2_var[fit_idx]
  paths$support_weight <- support_weight
  paths$uncertainty_weight <- pmin(1, uncertainty_weight)
  paths$outward_weight_raw <- paths$path_responsibility * paths$support_weight * paths$uncertainty_weight
  child_names <- unique(paths$one_step)
  cap_by_child <- vapply(child_names, function(child_name) {
    item <- nn_child_contexts[[child_name]]
    inward_sum <- if (!is.null(item)) {
      sum(normalize_nn_weights(item$parent_opportunity_weights, fallback_n = length(item$parent_fitness)))
    } else {
      1
    }
    if (!is.finite(inward_sum) || inward_sum <= 0) {
      inward_sum <- 1
    }
    nn_two_shell_max_weight_ratio * inward_sum
  }, numeric(1))
  cap_by_row <- cap_by_child[match(paths$one_step, child_names)]
  paths$outward_weight <- alfak_cpp_call(
    "alfak_group_cap_weights_cpp",
    alfak_group_cap_weights_cpp(
      group = as.character(paths$one_step),
      raw_weights = as.numeric(paths$outward_weight_raw),
      cap_by_row = as.numeric(cap_by_row)
    ),
    context = "compute_two_shell_outward_weights"
  )
  paths[is.finite(paths$outward_weight) & paths$outward_weight > 0, , drop = FALSE]
}

#' Apply the single two-shell backward correction to one-step fitness estimates
#' @keywords internal
#' @noRd
apply_two_shell_backward_correction <- function(nn_child_contexts, f1_initial, outward_paths,
                                                mu01, sigma01, mu12, sigma12, tau_reuse,
                                                nn_two_shell_lambda, timepoints, search_interval) {
  corrected <- f1_initial
  node_rows <- vector("list", length(nn_child_contexts))
  names(node_rows) <- names(nn_child_contexts)
  for (child_name in names(nn_child_contexts)) {
    item <- nn_child_contexts[[child_name]]
    child_paths <- if (is.data.frame(outward_paths) && nrow(outward_paths)) {
      outward_paths[outward_paths$one_step == child_name, , drop = FALSE]
    } else {
      outward_paths
    }
    inward_weights <- normalize_nn_weights(item$parent_opportunity_weights, fallback_n = length(item$parent_fitness))
    inward_sum <- sum(inward_weights)
    outward_sum <- if (is.data.frame(child_paths) && nrow(child_paths)) sum(child_paths$outward_weight) else 0
    fallback_reason <- NA_character_
    boundary_flag <- FALSE
    if (!is.finite(f1_initial[child_name])) {
      fallback_reason <- "initial_one_step_missing"
    } else if (!is.finite(outward_sum) || outward_sum <= 0) {
      fallback_reason <- "no_outward_weight_for_node"
    } else {
      sigma12_eff <- sqrt(sigma12^2 + child_paths$f2_var + tau_reuse^2)
      objective_fn <- function(fc_param) {
        alfak_cpp_call(
          "alfak_neighbor_two_shell_objective_cpp",
          alfak_neighbor_two_shell_objective_cpp(
          fc_param = fc_param,
          parent_fitness = item$parent_fitness,
          pij_values = item$pij,
          parent_birth_times = item$parent_birth_times,
          timepoints = timepoints,
          parent_xfit = item$parent_xfit,
          child_obs = item$child_obs,
          ntot = item$ntot,
          inward_prior_mean = mu01,
          inward_prior_sd = sigma01,
          inward_prior_weights = inward_weights,
          do_inward_prior = is.finite(mu01) && is.finite(sigma01) && sigma01 > 0 && length(inward_weights) > 0,
          outward_fitness = child_paths$f2_hat,
          outward_prior_mean = mu12,
          outward_prior_sd = sigma12_eff,
          outward_prior_weights = child_paths$outward_weight,
          outward_lambda = nn_two_shell_lambda,
          tol = ALFAK_FEXP_DELTA_TOL
          ),
          context = sprintf("refit empirical_two_shell child %s", child_name)
        )
      }
      local_interval <- range(c(search_interval, f1_initial[child_name], child_paths$f2_hat), na.rm = TRUE)
      local_span <- diff(local_interval)
      if (!is.finite(local_span) || local_span <= 0) {
        local_span <- 1
      }
      local_interval <- local_interval + c(-local_span, local_span)
      res <- run_optimise_checked(
        objective_fn,
        interval = local_interval,
        context = sprintf("optimise nearest-neighbour fitness with empirical_two_shell prior for child %s", child_name)
      )
      if (is.null(res)) {
        fallback_reason <- "two_shell_optimise_failed"
      } else {
        corrected[child_name] <- res$minimum
        boundary_flag <- res$minimum <= local_interval[1] + sqrt(.Machine$double.eps) ||
          res$minimum >= local_interval[2] - sqrt(.Machine$double.eps)
      }
    }
    node_rows[[child_name]] <- data.frame(
      karyotype = child_name,
      direct_observed_count = sum(item$child_obs),
      projected_exposure = item$projected_exposure,
      n_parent_0step = length(item$parent_fitness),
      n_descendant_2step = if (is.data.frame(child_paths)) length(unique(child_paths$descendant)) else 0L,
      inward_weight_sum = inward_sum,
      outward_weight_sum = if (is.finite(outward_sum)) outward_sum else 0,
      f1_initial = unname(f1_initial[child_name]),
      f1_two_shell = unname(corrected[child_name]),
      f1_delta_after_correction = unname(corrected[child_name] - f1_initial[child_name]),
      objective_boundary_flag = boundary_flag,
      prior_dominated_flag = sum(item$child_obs) == 0 && (inward_sum + outward_sum) > 0,
      outward_dominated_flag = is.finite(outward_sum) && outward_sum > inward_sum,
      fallback_reason = fallback_reason,
      stringsAsFactors = FALSE
    )
  }
  list(f1 = corrected, node_diagnostics = do.call(rbind, node_rows))
}

#' Run the empirical two-shell correction for one bootstrap replicate
#' @keywords internal
#' @noRd
run_empirical_two_shell_correction <- function(nn_child_contexts, nn_present, f1_initial,
                                               fpar, boot_data, timepoints, ntot,
                                               pm, n0, search_interval,
                                               inward_prior_fit,
                                               nn_prior_sd_floor,
                                               nn_two_shell_min_delta_n = 3L,
                                               nn_two_shell_min_exposure = NULL,
                                               nn_two_shell_min_observed_count = 1L,
                                               nn_two_shell_max_weight_ratio = 1.0,
                                               nn_two_shell_lambda = 1.0,
                                               nn_two_shell_reuse_sd = NULL,
                                               nn_two_shell_uncertainty_floor = NULL) {
  uncertainty_floor <- resolve_nn_two_shell_uncertainty_floor(
    nn_two_shell_uncertainty_floor = nn_two_shell_uncertainty_floor,
    nn_prior_sd_floor = nn_prior_sd_floor
  )
  empty_node_diag <- function(reason) {
    if (!length(nn_child_contexts)) {
      return(data.frame(
        karyotype = character(0),
        direct_observed_count = numeric(0),
        projected_exposure = numeric(0),
        n_parent_0step = integer(0),
        n_descendant_2step = integer(0),
        inward_weight_sum = numeric(0),
        outward_weight_sum = numeric(0),
        f1_initial = numeric(0),
        f1_two_shell = numeric(0),
        f1_delta_after_correction = numeric(0),
        objective_boundary_flag = logical(0),
        prior_dominated_flag = logical(0),
        outward_dominated_flag = logical(0),
        fallback_reason = character(0),
        stringsAsFactors = FALSE
      ))
    }
    rows <- lapply(names(nn_child_contexts), function(child_name) {
      item <- nn_child_contexts[[child_name]]
      data.frame(
        karyotype = child_name,
        direct_observed_count = sum(item$child_obs),
        projected_exposure = item$projected_exposure,
        n_parent_0step = length(item$parent_fitness),
        n_descendant_2step = 0L,
        inward_weight_sum = sum(normalize_nn_weights(item$parent_opportunity_weights, fallback_n = length(item$parent_fitness))),
        outward_weight_sum = 0,
        f1_initial = unname(f1_initial[child_name]),
        f1_two_shell = unname(f1_initial[child_name]),
        f1_delta_after_correction = 0,
        objective_boundary_flag = FALSE,
        prior_dominated_flag = FALSE,
        outward_dominated_flag = FALSE,
        fallback_reason = reason,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }
  two_shell <- build_two_shell_candidates(
    nn_child_contexts = nn_child_contexts,
    f1_hat = f1_initial,
    boot_data = boot_data,
    fpar = fpar,
    timepoints = timepoints,
    ntot = ntot,
    pm = pm,
    nn_present = nn_present,
    nn_two_shell_min_exposure = nn_two_shell_min_exposure,
    nn_two_shell_min_observed_count = nn_two_shell_min_observed_count,
    min_frequency = 1 / n0
  )
  candidates_total <- nrow(two_shell$candidates)
  candidates_retained <- nrow(two_shell$retained)
  if (!candidates_retained) {
    return(list(
      f1 = f1_initial,
      node_diagnostics = empty_node_diag("no_retained_two_step_candidates"),
      diagnostics = list(
        fallback_reason = "no_retained_two_step_candidates",
        n_2step_candidates_total = candidates_total,
        n_2step_candidates_retained = candidates_retained,
        n_2step_observed = 0L,
        n_2step_used_in_backward_term = 0L,
        mu01 = inward_prior_fit$prior_mean,
        sigma01 = inward_prior_fit$prior_sd,
        mu12 = NA_real_,
        sigma12 = NA_real_,
        tau_reuse = NA_real_,
        total_inward_weight = length(nn_child_contexts),
        total_outward_weight = 0,
        adaptive_two_shell_min_exposure = two_shell$adaptive_min_exposure
      )
    ))
  }
  retained <- compute_two_shell_path_responsibilities(two_shell$retained)
  f2_fit <- estimate_provisional_two_step_fitness(
    candidate_paths = retained,
    anchor_states = two_shell$anchor_states,
    boot_data = boot_data,
    timepoints = timepoints,
    ntot = ntot,
    search_interval = search_interval,
    uncertainty_floor = uncertainty_floor
  )
  usable_f2 <- f2_fit[is.finite(f2_fit$f2_hat) & !f2_fit$fit_failed, , drop = FALSE]
  if (!nrow(usable_f2)) {
    return(list(
      f1 = f1_initial,
      node_diagnostics = empty_node_diag("all_provisional_two_step_fits_failed"),
      diagnostics = list(
        fallback_reason = "all_provisional_two_step_fits_failed",
        n_2step_candidates_total = candidates_total,
        n_2step_candidates_retained = candidates_retained,
        n_2step_observed = as.integer(sum(f2_fit$observed_count >= nn_two_shell_min_observed_count, na.rm = TRUE)),
        n_2step_used_in_backward_term = 0L,
        mu01 = inward_prior_fit$prior_mean,
        sigma01 = inward_prior_fit$prior_sd,
        mu12 = NA_real_,
        sigma12 = NA_real_,
        tau_reuse = NA_real_,
        total_inward_weight = length(nn_child_contexts),
        total_outward_weight = 0,
        adaptive_two_shell_min_exposure = two_shell$adaptive_min_exposure
      )
    ))
  }
  retained <- retained[retained$descendant %in% usable_f2$karyotype, , drop = FALSE]
  fit_idx <- match(retained$descendant, usable_f2$karyotype)
  delta12 <- usable_f2$f2_hat[fit_idx] - f1_initial[retained$one_step]
  shell12 <- fit_shell_delta_prior(
    delta = delta12,
    weights = retained$path_responsibility,
    sd_floor = nn_prior_sd_floor,
    fallback_mu = inward_prior_fit$prior_mean,
    fallback_sd = inward_prior_fit$prior_sd
  )
  if (shell12$n < nn_two_shell_min_delta_n) {
    return(list(
      f1 = f1_initial,
      node_diagnostics = empty_node_diag("too_few_usable_delta12"),
      diagnostics = list(
        fallback_reason = "too_few_usable_delta12",
        n_2step_candidates_total = candidates_total,
        n_2step_candidates_retained = candidates_retained,
        n_2step_observed = as.integer(sum(usable_f2$observed_count >= nn_two_shell_min_observed_count, na.rm = TRUE)),
        n_2step_used_in_backward_term = 0L,
        mu01 = inward_prior_fit$prior_mean,
        sigma01 = inward_prior_fit$prior_sd,
        mu12 = shell12$mu,
        sigma12 = shell12$sigma,
        tau_reuse = NA_real_,
        total_inward_weight = length(nn_child_contexts),
        total_outward_weight = 0,
        adaptive_two_shell_min_exposure = two_shell$adaptive_min_exposure
      )
    ))
  }
  tau_reuse <- if (is.null(nn_two_shell_reuse_sd)) {
    max(nn_prior_sd_floor, shell12$sigma)
  } else {
    nn_two_shell_reuse_sd
  }
  outward_paths <- compute_two_shell_outward_weights(
    candidate_paths = retained,
    f2_fit = usable_f2,
    nn_child_contexts = nn_child_contexts,
    sigma12 = shell12$sigma,
    tau_reuse = tau_reuse,
    nn_two_shell_max_weight_ratio = nn_two_shell_max_weight_ratio
  )
  total_outward <- if (nrow(outward_paths)) sum(outward_paths$outward_weight) else 0
  if (!is.finite(total_outward) || total_outward <= 0 || nn_two_shell_lambda <= 0) {
    return(list(
      f1 = f1_initial,
      node_diagnostics = empty_node_diag("all_outward_weights_zero"),
      diagnostics = list(
        fallback_reason = "all_outward_weights_zero",
        n_2step_candidates_total = candidates_total,
        n_2step_candidates_retained = candidates_retained,
        n_2step_observed = as.integer(sum(usable_f2$observed_count >= nn_two_shell_min_observed_count, na.rm = TRUE)),
        n_2step_used_in_backward_term = 0L,
        mu01 = inward_prior_fit$prior_mean,
        sigma01 = inward_prior_fit$prior_sd,
        mu12 = shell12$mu,
        sigma12 = shell12$sigma,
        tau_reuse = tau_reuse,
        total_inward_weight = length(nn_child_contexts),
        total_outward_weight = 0,
        adaptive_two_shell_min_exposure = two_shell$adaptive_min_exposure
      )
    ))
  }
  corrected <- apply_two_shell_backward_correction(
    nn_child_contexts = nn_child_contexts,
    f1_initial = f1_initial,
    outward_paths = outward_paths,
    mu01 = inward_prior_fit$prior_mean,
    sigma01 = inward_prior_fit$prior_sd,
    mu12 = shell12$mu,
    sigma12 = shell12$sigma,
    tau_reuse = tau_reuse,
    nn_two_shell_lambda = nn_two_shell_lambda,
    timepoints = timepoints,
    search_interval = search_interval
  )
  list(
    f1 = corrected$f1,
    node_diagnostics = corrected$node_diagnostics,
    diagnostics = list(
      fallback_reason = NA_character_,
      n_2step_candidates_total = candidates_total,
      n_2step_candidates_retained = candidates_retained,
      n_2step_observed = as.integer(sum(usable_f2$observed_count >= nn_two_shell_min_observed_count, na.rm = TRUE)),
      n_2step_used_in_backward_term = length(unique(outward_paths$descendant)),
      mu01 = inward_prior_fit$prior_mean,
      sigma01 = inward_prior_fit$prior_sd,
      mu12 = shell12$mu,
      sigma12 = shell12$sigma,
      tau_reuse = tau_reuse,
      total_inward_weight = length(nn_child_contexts),
      total_outward_weight = total_outward,
      adaptive_two_shell_min_exposure = two_shell$adaptive_min_exposure
    )
  )
}

#' Generate all single-step neighbours for karyotype IDs
#' @keywords internal
#' @noRd
gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) {
    cpp_neighbors <- alfak_cpp_call(
      "gen_all_neighbours_cpp",
      gen_all_neighbours_cpp(as.character(ids), remove_nullisomes = isTRUE(remove_nullisomes)),
      context = "gen_all_neighbours"
    )
    if (is.matrix(cpp_neighbors)) {
      return(cpp_neighbors)
    }
    alfak_log_event(
      level = "ERROR",
      component = "cpp.gen_all_neighbours_cpp",
      detail = "C++ kernel returned malformed output in gen_all_neighbours."
    )
    stop("C++ kernel `gen_all_neighbours_cpp` returned malformed output.", call. = FALSE)
  }
  nkern <- do.call(rbind, lapply(1:length(ids[[1]]), function(i) {
    x0 <- rep(0, length(ids[[1]]))
    x1 <- x0
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0, x1)
  }))
  n <- do.call(rbind, lapply(ids, function(ii) t(apply(nkern, 1, function(i) i + ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind, ids), n) # Intentionally add originals back
  n <- unique(n) # Keep unique set
  n <- n[-(1:nids), , drop=FALSE] # Remove the nids original ids that were just added
  # drop=FALSE ensures matrix structure even if 1 row left
  if (remove_nullisomes && nrow(n) > 0) # Check nrow > 0 before apply
    n <- n[apply(n, 1, function(ni) sum(ni < 1) == 0), , drop = FALSE]
  n
}

fields_krig_find_lambda <- local({
  krig_find_lambda <- get("KrigFindLambda", envir = asNamespace("fields"))
  function(...) krig_find_lambda(...)
})

fields_krig_coef <- local({
  krig_coef <- get("Krig.coef", envir = asNamespace("fields"))
  function(...) krig_coef(...)
})

fields_krig_parameters <- local({
  krig_parameters <- get("Krig.parameters", envir = asNamespace("fields"))
  function(...) krig_parameters(...)
})

fields_krig_ynew <- local({
  krig_ynew <- get("Krig.ynew", envir = asNamespace("fields"))
  function(...) krig_ynew(...)
})

krig_covariance_args <- function() {
  list(
    Covariance = "Matern",
    smoothness = 1.5
  )
}

predict_cached_krig <- function(object, x, dist_mat, Z = NULL, drop.Z = FALSE, just.fixed = FALSE) {
  if (is.null(x)) {
    x <- object$x
  } else {
    x <- as.matrix(x)
  }
  if (is.null(Z)) {
    Z <- object$Z
  } else {
    Z <- as.matrix(Z)
  }

  x_scaled <- scale(x, object$transform$x.center, object$transform$x.scale)
  null_fun <- get(object$null.function.name, envir = asNamespace("fields"))
  Tmatrix <- do.call(null_fun, c(object$null.args, list(x = x_scaled, Z = Z, drop.Z = drop.Z)))

  if (drop.Z) {
    pred <- Tmatrix %*% object$d[object$ind.drift]
  } else {
    pred <- Tmatrix %*% object$d
  }

  if (!just.fixed) {
    cov_args <- object$args
    cov_args$distMat <- dist_mat
    cov_fun <- get(object$cov.function.name, envir = asNamespace("fields"))
    pred <- pred + do.call(cov_fun, c(cov_args, list(x1 = x_scaled, x2 = object$knots, C = object$c)))
  }

  as.numeric(pred)
}

build_cached_krig_fit <- function(ktrain, y, kpred = NULL, give_warnings = TRUE) {
  train_dist <- fields::rdist(ktrain, ktrain)
  pred_dist <- NULL
  if (!is.null(kpred)) {
    pred_dist <- fields::rdist(kpred, ktrain)
  }
  fit <- fields::Krig(
    ktrain,
    y,
    cov.function = "stationary.cov",
    cov.args = krig_covariance_args(),
    nstep.cv = ALFAK_KRIG_NSTEP_CV,
    give.warnings = give_warnings
  )
  list(fit = fit, train_dist = train_dist, pred_dist = pred_dist)
}

refit_cached_krig <- function(cache, y, x_pred = NULL, pred_dist = cache$pred_dist, give_warnings = TRUE) {
  fit <- cache$fit
  gcv_out <- fields_krig_find_lambda(
    fit,
    nstep.cv = ALFAK_KRIG_NSTEP_CV,
    cost = fit$cost,
    offset = fit$offset,
    y = y,
    give.warnings = give_warnings
  )

  if (fit$method != "user") {
    fit$lambda <- gcv_out$lambda.est[fit$method, 1]
    fit$eff.df <- gcv_out$lambda.est[fit$method, 2]
  }

  y_info <- fields_krig_ynew(fit, y = y)
  coef_out <- fields_krig_coef(fit, lambda = fit$lambda, y = y)

  fit$gcv.grid <- gcv_out$gcv.grid
  fit$lambda.est <- gcv_out$lambda.est
  fit$warningTable <- gcv_out$warningTable
  fit$y <- as.numeric(y)
  fit$yM <- y_info$yM
  fit$c <- coef_out$c
  fit$d <- coef_out$d
  fit$tauHat.rep <- coef_out$tauHat.rep
  fit$tauHat.pure.error <- coef_out$tauHat.pure.error
  fit$pure.ss <- coef_out$pure.ss
  fit$fitted.values <- predict_cached_krig(fit, fit$x, dist_mat = cache$train_dist, Z = fit$Z)
  fit$residuals <- fit$y - fit$fitted.values

  fit_params <- fields_krig_parameters(fit)
  fit[names(fit_params)] <- fit_params

  if (fit$method == "user" && !is.na(fit$tau2)) {
    fit$best.model <- c(fit$lambda, fit$tau2, fit$sigma)
  } else {
    fit$best.model <- c(fit$lambda, fit$tauHat.MLE^2, fit$sigmahat)
  }
  fit$rhohat <- fit$sigmahat

  preds <- NULL
  if (!is.null(x_pred)) {
    preds <- predict_cached_krig(fit, x_pred, dist_mat = pred_dist)
  }

  list(fit = fit, preds = preds)
}

#' Bootstrap counts from data matrix
#' @keywords internal
#' @noRd
bootstrap_counts <- function(data) {
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  boot_data <- matrix(NA, nrow = num_species, ncol = num_timepoints)
  rownames(boot_data) <- rownames(data)
  for (i in seq_len(num_timepoints)) {
    total_counts <- sum(data[, i])
    if (total_counts == 0) {
      boot_data[, i] <- rep(0, num_species)
    } else {
      boot_data[, i] <- stats::rmultinom(1, size = total_counts, prob = data[, i] / total_counts)
    }
  }
  boot_data
}

#' Compute dx/dt (rate of change)
#' @keywords internal
#' @noRd
compute_dx_dt <- function(x, timepoints) {
  if (length(timepoints) != ncol(x)) {
    stop("`timepoints` must have length ncol(x).")
  }
  delta_t <- diff(timepoints)
  if (any(!is.finite(delta_t)) || any(delta_t <= 0)) {
    stop("`timepoints` must be finite and strictly increasing for `compute_dx_dt()`.")
  }
  sweep(x[, -1, drop = FALSE] - x[, -ncol(x), drop = FALSE], 2, delta_t, "/")
}

#' Calculate log-sum-exp for numerical stability
#' @keywords internal
#' @noRd
logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

#' Generate nearest neighbour information
#' @keywords internal
#' @noRd
gen_nn_info <- function(fq, pm = 0.00005) {
  validate_probability(pm, "pm", upper_inclusive = TRUE)
  cpp_info <- alfak_cpp_call(
    "gen_nn_info_cpp",
    gen_nn_info_cpp(as.character(fq), beta = pm),
    context = "gen_nn_info"
  )
  if (is.list(cpp_info)) {
    return(cpp_info)
  }
  alfak_log_event(
    level = "ERROR",
    component = "cpp.gen_nn_info_cpp",
    detail = "C++ kernel returned malformed output in gen_nn_info."
  )
  stop("C++ kernel `gen_nn_info_cpp` returned malformed output.", call. = FALSE)
}

#' Negative log-likelihood calculation
#' @keywords internal
#' @noRd
neg_log_lik <- function(param, counts, timepoints) {
  alfak_cpp_call(
    "alfak_neg_log_lik_cpp",
    alfak_neg_log_lik_cpp(param, counts, timepoints),
    context = "neg_log_lik"
  )
}

#' Jointly optimize fitness and initial frequencies
#' @keywords internal
#' @noRd
joint_optimize <- function(counts, timepoints, f_init, x0_init) {
  K <- length(f_init)
  if (K == 1) {
    return(list(f = 0, x0 = 1))
  }
  f_free_init <- f_init[seq_len(K - 1)]
  x0_init_log <- free_softmax_logits(x0_init)
  param_init <- c(f_free_init, x0_init_log)
  obj_fun <- function(par) neg_log_lik(par, counts, timepoints)
  opt <- run_optim_checked(par = param_init, fn = obj_fun,
                           method = "BFGS",
                           control = list(maxit = 200, reltol = 1e-8),
                           context = "joint_optimize")
  f_free_opt <- opt$par[seq_len(K - 1)]
  f_opt <- c(f_free_opt, -sum(f_free_opt))
  log_x0_opt <- opt$par[K:(2 * K - 2)]
  x0_opt <- softmax_from_free_logits(log_x0_opt)
  list(f = f_opt, x0 = x0_opt)
}

#' Project frequencies forward in time using log-space calculations
#' @keywords internal
#' @noRd
project_forward_log <- function(x0, f, timepoints) {
  alfak_cpp_call(
    "alfak_project_forward_log_cpp",
    alfak_project_forward_log_cpp(x0, f, timepoints),
    context = "project_forward_log"
  )
}

#' Optimize initial frequencies given observed data and fitness values
#' @keywords internal
#' @noRd
optimize_initial_frequencies <- function(x_obs, f, timepoints) {
  K <- nrow(x_obs)
  if (K == 1) {
    return(1)
  }
  loss_function <- function(log_x0_free) {
    x0 <- softmax_from_free_logits(log_x0_free)
    x_pred <- project_forward_log(x0, f, timepoints)
    sum((x_pred - x_obs)^2)
  }
  x_ini <- x_obs[, 1] + 1e-6 # Original had this epsilon
  x_ini <- x_ini / sum(x_ini)
  opt_result <- run_optim_checked(par = free_softmax_logits(x_ini), fn = loss_function,
                                  method = "BFGS",
                                  control = list(maxit = 500, reltol = 1e-8),
                                  context = "optimize_initial_frequencies")
  softmax_from_free_logits(opt_result$par)
}

#' Find "birth times" for species based on reaching a minimum frequency
#' @keywords internal
#' @noRd
find_birth_times <- function(opt_res, time_range, minF) {
  f_est <- opt_res$f
  x0_est <- opt_res$x0
  num_species <- length(f_est)
  birth_times <- rep(NA, num_species)
  for (i in seq_len(num_species)) {
    if (f_est[i] <= min(f_est)) next
    birth_fn <- function(t) {
      log_x_t <- log(x0_est[i]) + f_est[i] * t # Original log(x0_est[i])
      denom <- logSumExp(log(x0_est) + f_est * t) # Original log(x0_est)
      exp(log_x_t - denom) - minF
    }
    root <- try(stats::uniroot(birth_fn, range(time_range), tol = 1e-6)$root, silent = TRUE)
    if (!inherits(root, "try-error")) {
      birth_times[i] <- root
    }
  }
  birth_times
}

##########################################
# Core functions (now serial)
##########################################

#' Solve fitness using bootstrap approach (Internal function)
#' @keywords internal
#' @noRd
solve_fitness_bootstrap <- function(data, minobs, nboot = 1000, epsilon = 1e-6, pm = 0.00005,
                                    n0, nb, passage_times = NULL, allow_noninteger_counts = FALSE, correct_efflux=FALSE,
                                    nn_prior = c("empirical_censored", "empirical_censored_weighted", "empirical_two_shell", "cohort_transition", "none", "empirical"),
                                    cohort_transition_prior = NULL,
                                    cohort_transition_patient_id = NULL,
                                    cohort_transition_version = c("contextual", "v2", "v1"),
                                    cohort_transition_apply_to = c("zero_only", "low_information", "all"),
                                    cohort_transition_overlay_base = c("empirical_two_shell", "direct"),
                                    cohort_transition_lambda = 0.25,
                                    cohort_transition_max_borrowing_fraction = 0.5,
                                    cohort_transition_max_abs_delta_shift = NULL,
                                    cohort_transition_sd_floor = 0.05,
                                    cohort_transition_patient_sd_floor = 0.1,
                                    cohort_contextual_apply_to = NULL,
                                    cohort_contextual_overlay_base = c("empirical_two_shell", "direct"),
                                    cohort_context_lambda = 0.25,
                                    cohort_context_max_borrowing_fraction = 0.5,
                                    cohort_context_max_abs_delta_shift = NULL,
                                    cohort_context_sd_floor = 0.05,
                                    cohort_context_patient_sd_floor = 0.10,
                                    cohort_context_keep_baseline_when_sparse = TRUE,
                                    cohort_context_keep_baseline_when_high_variable = TRUE,
                                    nn_prior_sd = NULL,
                                    nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                                    nn_prior_grid_n = ALFAK_NN_PRIOR_CENSORED_GRID_POINTS,
                                    nn_prior_fit_subset = c("hybrid", "all"),
                                    nn_prior_zero_exposure_min = NULL,
                                    nn_prior_zero_exposure_quantile = 0.10,
                                    nn_prior_zero_weight_scale = 0.50,
                                    nn_prior_zero_weight_cap_ratio = NULL,
                                    nn_prior_zero_birth_fallback_weight = NULL,
                                    nn_prior_zero_birth_child_floor = 0.25,
                                    nn_prior_zero_birth_child_shape = 1,
                                    nn_prior_zero_birth_replicate_floor = 0.50,
                                    nn_prior_zero_birth_replicate_shape = 1,
                                    nn_prior_hybrid_min_obs = 3L,
                                    nn_prior_two_step_support = c("none", "rescue"),
                                    nn_prior_two_step_support_min = 0.15,
                                    nn_prior_two_step_cap_floor = 0.30,
                                    nn_two_shell_min_delta_n = 3L,
                                    nn_two_shell_min_exposure = NULL,
                                    nn_two_shell_min_observed_count = 1L,
                                    nn_two_shell_max_weight_ratio = 1.0,
                                    nn_two_shell_lambda = 1.0,
                                    nn_two_shell_reuse_sd = NULL,
                                    nn_two_shell_uncertainty_floor = NULL) {
  data$x <- coerce_count_matrix(data$x, allow_noninteger_counts = allow_noninteger_counts)
  validate_positive_depth(data$x)
  validate_positive_integer(nboot, "nboot")
  validate_positive_finite(n0, "n0")
  validate_positive_finite(nb, "nb")
  validate_probability(pm, "pm", upper_inclusive = TRUE)
  validate_scalar_logical(allow_noninteger_counts, "allow_noninteger_counts")
  validate_scalar_logical(correct_efflux, "correct_efflux")
  nn_prior <- validate_nn_prior_mode(nn_prior)
  cohort_transition_prior_use <- NULL
  if (identical(nn_prior, "cohort_transition")) {
    cohort_transition_version <- match.arg(cohort_transition_version)
    cohort_transition_apply_to <- match.arg(cohort_transition_apply_to)
    cohort_transition_overlay_base <- match.arg(cohort_transition_overlay_base)
    if (is.null(cohort_contextual_apply_to)) {
      cohort_contextual_apply_to <- cohort_transition_apply_to
    } else {
      cohort_contextual_apply_to <- match.arg(cohort_contextual_apply_to, c("zero_only", "low_information", "all"))
    }
    cohort_contextual_overlay_base <- match.arg(cohort_contextual_overlay_base)
    validate_nonnegative_finite(cohort_transition_lambda, "cohort_transition_lambda")
    validate_probability(cohort_transition_max_borrowing_fraction, "cohort_transition_max_borrowing_fraction", upper_inclusive = TRUE)
    if (!is.null(cohort_transition_max_abs_delta_shift)) {
      validate_positive_finite(cohort_transition_max_abs_delta_shift, "cohort_transition_max_abs_delta_shift")
    }
    validate_positive_finite(cohort_transition_sd_floor, "cohort_transition_sd_floor")
    validate_positive_finite(cohort_transition_patient_sd_floor, "cohort_transition_patient_sd_floor")
    validate_nonnegative_finite(cohort_context_lambda, "cohort_context_lambda")
    validate_probability(cohort_context_max_borrowing_fraction, "cohort_context_max_borrowing_fraction", upper_inclusive = TRUE)
    if (!is.null(cohort_context_max_abs_delta_shift)) {
      validate_positive_finite(cohort_context_max_abs_delta_shift, "cohort_context_max_abs_delta_shift")
    }
    validate_positive_finite(cohort_context_sd_floor, "cohort_context_sd_floor")
    validate_positive_finite(cohort_context_patient_sd_floor, "cohort_context_patient_sd_floor")
    validate_scalar_logical(cohort_context_keep_baseline_when_sparse, "cohort_context_keep_baseline_when_sparse")
    validate_scalar_logical(cohort_context_keep_baseline_when_high_variable, "cohort_context_keep_baseline_when_high_variable")
    cohort_transition_prior <- resolve_cohort_transition_prior_object(
      cohort_transition_prior = cohort_transition_prior,
      cohort_transition_prior_path = NULL,
      cohort_transition_patient_id = cohort_transition_patient_id
    )
    cohort_transition_prior_use <- cohort_transition_prior_for_patient(
      cohort_transition_prior,
      patient_id = cohort_transition_patient_id
    )
  } else {
    cohort_transition_version <- cohort_transition_version[1]
    cohort_transition_apply_to <- cohort_transition_apply_to[1]
    cohort_transition_overlay_base <- cohort_transition_overlay_base[1]
    cohort_contextual_apply_to <- cohort_contextual_apply_to %||% cohort_transition_apply_to
    cohort_contextual_overlay_base <- cohort_contextual_overlay_base[1]
  }
  nn_prior_fit_subset <- validate_nn_prior_fit_subset(nn_prior_fit_subset)
  nn_prior_two_step_support <- validate_nn_prior_two_step_support(nn_prior_two_step_support)
  validate_nn_prior_controls(
    nn_prior_sd = nn_prior_sd,
    nn_prior_sd_floor = nn_prior_sd_floor,
    nn_prior_grid_n = nn_prior_grid_n,
    nn_prior_fit_subset = nn_prior_fit_subset,
    nn_prior_zero_exposure_min = nn_prior_zero_exposure_min,
    nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
    nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
    nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
    nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
    nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
    nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
    nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
    nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
    nn_prior_hybrid_min_obs = nn_prior_hybrid_min_obs,
    nn_prior_two_step_support = nn_prior_two_step_support,
    nn_prior_two_step_support_min = nn_prior_two_step_support_min,
    nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
    nn_two_shell_min_delta_n = nn_two_shell_min_delta_n,
    nn_two_shell_min_exposure = nn_two_shell_min_exposure,
    nn_two_shell_min_observed_count = nn_two_shell_min_observed_count,
    nn_two_shell_max_weight_ratio = nn_two_shell_max_weight_ratio,
    nn_two_shell_lambda = nn_two_shell_lambda,
    nn_two_shell_reuse_sd = nn_two_shell_reuse_sd,
    nn_two_shell_uncertainty_floor = nn_two_shell_uncertainty_floor
  )
  fq <- get_frequent_karyotypes(data$x, minobs)
  nn_info_list <- gen_nn_info(fq, pm) # Renamed 'nn' to 'nn_info_list' for clarity
  if (length(nn_info_list) > 0 && !is.null(nn_info_list[[1]]$ni)) { # Check if naming is needed
    names(nn_info_list) <- sapply(nn_info_list, function(nni) nni$ni)
  } else if (length(nn_info_list) == 0) {
    # names(nn_info_list) would be NULL, which is fine for later checks
  } else {
    warning("nn_info_list structure unexpected for naming in solve_fitness_bootstrap")
  }

  fq_vec <- do.call(rbind, lapply(fq, s2v))
  rownames(fq_vec) <- fq
  viability <- prepare_efflux_viability(fq_vec, pm = pm, correct_efflux = correct_efflux)

  # fq_nn <- which(as.matrix(stats::dist(fq_vec)) == 1) # fq_nn was not used
  timepoints <- resolve_time_axis(data, passage_times)
  num_species <- length(fq)
  num_timepoints <- ncol(data$x)
  weighted_sample_pooled_prior <- list(
    available = FALSE,
    prior_mean = NA_real_,
    prior_sd = NA_real_,
    informative_child_count = NA_integer_,
    sum_child_weight = 0,
    effective_mass_reference = 0,
    exposure_reference = NA_real_
  )
  cohort_transition_uses_two_shell_baseline <- identical(nn_prior, "cohort_transition") &&
    ((identical(cohort_transition_version, "v2") &&
        identical(cohort_transition_overlay_base, "empirical_two_shell")) ||
       (identical(cohort_transition_version, "contextual") &&
          identical(cohort_contextual_overlay_base, "empirical_two_shell")))
  if ((nn_prior %in% c("empirical_censored_weighted", "empirical_two_shell") ||
       isTRUE(cohort_transition_uses_two_shell_baseline)) &&
      length(nn_info_list) > 0) {
    weighted_sample_pooled_prior <- tryCatch(
      {
        sample_state <- prepare_bootstrap_nn_dataset_state(
          count_data = data$x,
          current_fq = fq,
          current_timepoints = timepoints,
          current_epsilon = epsilon,
          current_n0 = n0,
          current_nb = nb,
          current_viability = viability,
          current_nn_info = nn_info_list,
          correct_efflux = correct_efflux,
          context = "sample-pooled weighted prior reference"
        )
        sample_build_opt_fc <- make_nn_child_objective_builder(
          timepoints = timepoints,
          ntot_rounded = sample_state$ntot_rounded
        )
        sample_nn_present <- resolve_nn_present_mask(sample_state$nn_child_contexts, data$x)
        estimate_weighted_sample_pooled_prior(
          nn_child_contexts = sample_state$nn_child_contexts,
          nn_present = sample_nn_present,
          fpar = sample_state$fpar,
          build_opt_fc = sample_build_opt_fc,
          search_interval = sample_state$search_interval,
          nn_prior_sd = nn_prior_sd,
          nn_prior_sd_floor = nn_prior_sd_floor,
          nn_prior_grid_n = nn_prior_grid_n,
          ntot = sample_state$ntot_rounded,
          context = "fit weighted sample-pooled latent-neighbour prior"
        )
      },
      error = function(e) weighted_sample_pooled_prior
    )
  }

  bootstrap_iter <- function(b_iter_idx, current_data, current_fq, current_timepoints,
                             current_num_species, current_num_timepoints,
                             current_epsilon, current_n0, current_nb,
                             current_viability, current_nn_info) { # Renamed arguments
    boot_data <- bootstrap_counts(current_data$x) # Bootstrap from original full data
    dataset_state <- prepare_bootstrap_nn_dataset_state(
      count_data = boot_data,
      current_fq = current_fq,
      current_timepoints = current_timepoints,
      current_epsilon = current_epsilon,
      current_n0 = current_n0,
      current_nb = current_nb,
      current_viability = current_viability,
      current_nn_info = current_nn_info,
      correct_efflux = correct_efflux,
      context = sprintf("bootstrap replicate %d", b_iter_idx)
    )
    f_qp <- dataset_state$f_initial
    fpar <- dataset_state$fpar
    f_final <- dataset_state$f_final
    x0_init <- dataset_state$x0_initial
    x0_final <- dataset_state$x0_final
    ntot_rounded <- dataset_state$ntot_rounded
    nn_child_contexts <- dataset_state$nn_child_contexts
    build_opt_fc <- make_nn_child_objective_builder(
      timepoints = current_timepoints,
      ntot_rounded = ntot_rounded
    )
    search_interval <- dataset_state$search_interval
    nn_present <- resolve_nn_present_mask(nn_child_contexts, boot_data)

    fc <- rep(NaN, length(nn_child_contexts))
    names(fc) <- names(nn_child_contexts) # Pre-name fc
	    nn_prior_diag <- new_nn_prior_diagnostics(
	      nn_prior_mode_requested = nn_prior,
	      nn_prior_fit_subset_used = if (nn_prior %in% c("empirical_censored_weighted", "empirical_two_shell") ||
                                       isTRUE(cohort_transition_uses_two_shell_baseline)) nn_prior_fit_subset else NA_character_
	    )
	    nn_prior_diag$replicate_id <- as.integer(b_iter_idx)
	    nn_prior_diag$n_frequent_parents <- as.integer(length(fpar))
	    nn_prior_diag$n_1step_candidates <- as.integer(length(nn_child_contexts))
	    nn_prior_diag$n_1step_observed <- as.integer(sum(nn_present))
	    nn_prior_diag$n_1step_zero <- as.integer(sum(!nn_present))
	    nn_prior_diag$n_observed_children <- as.integer(sum(nn_present))
	    nn_prior_diag$n_zero_children_total <- as.integer(sum(!nn_present))
    if (any(!nn_present)) {
      nn_prior_diag$n_zero_children_with_birth_fallback <- as.integer(sum(vapply(
        nn_child_contexts[!nn_present],
        function(item) any(item$parent_birth_fallback),
        logical(1)
      )))
    }

    if (any(nn_present)) {
      # Use names for sapply for robustness if current_nn_info can be sparse/differently ordered
      sapply_names <- names(nn_child_contexts)[nn_present]
      if(length(sapply_names) > 0) { # Ensure there are names to iterate over
        for (child_name in sapply_names) {
          objective_fn <- build_opt_fc(nn_child_contexts[[child_name]], do_prior_param = FALSE)
          res <- run_optimise_checked(objective_fn, interval = search_interval,
                                      context = sprintf("optimise nearest-neighbour fitness for observed child %s", child_name))
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
    }

    fc_prior_vals <- numeric(0) # Child-minus-parent-mean deltas for the neighbour prior
    if(any(nn_present)){
      # Calculate differences based on nn_info items that were present AND had their fc computed
      nn_present_and_fc_computed <- nn_present & !is.na(fc)
      if(any(nn_present_and_fc_computed)){
        fc_prior_vals <- unlist(lapply(nn_child_contexts[nn_present_and_fc_computed], function(nni_item) {
          parent_f_mean <- nni_item$parent_fitness_mean_pij
          if(is.finite(parent_f_mean)){
            fc[nni_item$ni] - parent_f_mean
          } else {
            numeric(0) # No valid parents to compute difference from
          }
        }))
        fc_prior_vals <- fc_prior_vals[is.finite(fc_prior_vals)]
      }
    }

	    use_empirical_prior <- nn_prior == "empirical"
	    use_empirical_censored_prior <- nn_prior == "empirical_censored"
	    use_empirical_censored_weighted_prior <- nn_prior == "empirical_censored_weighted"
	    use_empirical_two_shell_prior <- nn_prior == "empirical_two_shell"
	    use_cohort_transition_prior <- nn_prior == "cohort_transition"
	    use_cohort_transition_v2_overlay <- use_cohort_transition_prior && identical(cohort_transition_version, "v2")
	    use_cohort_transition_contextual_overlay <- use_cohort_transition_prior && identical(cohort_transition_version, "contextual")
	    run_two_shell_baseline <- use_empirical_two_shell_prior ||
        (use_cohort_transition_v2_overlay && identical(cohort_transition_overlay_base, "empirical_two_shell")) ||
        (use_cohort_transition_contextual_overlay && identical(cohort_contextual_overlay_base, "empirical_two_shell"))
	    use_weighted_like_prior <- use_empirical_censored_weighted_prior || run_two_shell_baseline
	    weighted_prior_config <- NULL
	    weighted_sample_pooled_prior_use <- NULL
	    two_shell_inward_prior_fit <- NULL
	    if (use_weighted_like_prior && any(!nn_present)) {
	      weighted_prior_config <- prepare_weighted_nn_prior_fit(
        nn_child_contexts = nn_child_contexts,
        nn_present = nn_present,
        nn_prior_fit_subset = nn_prior_fit_subset,
        nn_prior_zero_exposure_min = nn_prior_zero_exposure_min,
        nn_prior_zero_exposure_quantile = nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_birth_fallback_weight = nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = nn_prior_zero_birth_child_floor,
        nn_prior_zero_birth_child_shape = nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = nn_prior_zero_birth_replicate_floor,
        nn_prior_zero_birth_replicate_shape = nn_prior_zero_birth_replicate_shape,
        nn_prior_hybrid_min_obs = nn_prior_hybrid_min_obs,
        nn_prior_two_step_support = nn_prior_two_step_support,
        nn_prior_two_step_support_min = nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = nn_prior_two_step_cap_floor,
        count_data = boot_data,
        pm = pm,
        ntot = ntot_rounded
      )
      nn_prior_diag[names(weighted_prior_config$diagnostics)] <- weighted_prior_config$diagnostics
      nn_prior_diag$sample_pooled_prior_available <- isTRUE(weighted_sample_pooled_prior$available)
      if (isTRUE(weighted_sample_pooled_prior$available)) {
        nn_prior_diag$sample_pooled_prior_mu <- weighted_sample_pooled_prior$prior_mean
        nn_prior_diag$sample_pooled_prior_sigma <- weighted_sample_pooled_prior$prior_sd
        nn_prior_diag$sample_pooled_prior_informative_child_count <- weighted_sample_pooled_prior$informative_child_count
      }
      if (!isTRUE(weighted_prior_config$can_fit_replicate_prior)) {
        weighted_sample_pooled_prior_use <- resolve_weighted_sample_pooled_fallback_prior(
          sample_pooled_prior = weighted_sample_pooled_prior,
          weighted_diagnostics = weighted_prior_config$diagnostics
        )
        if (isTRUE(weighted_sample_pooled_prior_use$available)) {
          nn_prior_diag$nn_prior_mode_used <- "empirical_censored_weighted"
          nn_prior_diag$nn_prior_source_used <- "sample_pooled"
          nn_prior_diag$used_sample_pooled_fallback_for_this_replicate <- TRUE
          nn_prior_diag$sample_pooled_alpha_used <- weighted_sample_pooled_prior_use$alpha
          nn_prior_diag$sample_pooled_sigma_used <- weighted_sample_pooled_prior_use$prior_sd
        } else {
          nn_prior_diag$nn_prior_mode_used <- "none"
          nn_prior_diag$nn_prior_source_used <- "none"
          nn_prior_diag$used_no_prior_fallback_for_this_replicate <- TRUE
        }
      } else {
        nn_prior_diag$nn_prior_source_used <- "observed_replicate"
      }
      if (isTRUE(weighted_prior_config$diagnostics$used_no_prior_fallback_for_this_replicate)) {
        nn_prior_diag$nn_prior_mode_used <- "none"
      }
    }

    nn_cohort_transition_node_diagnostics <- data.frame()
    if (use_cohort_transition_prior && identical(cohort_transition_version, "v1") && length(nn_child_contexts) > 0) {
      node_rows <- vector("list", length(nn_child_contexts))
      names(node_rows) <- names(nn_child_contexts)
      for (child_name in names(nn_child_contexts)) {
        fit_res <- fit_cohort_transition_nn_child(
          item = nn_child_contexts[[child_name]],
          child_name = child_name,
          build_opt_fc = build_opt_fc,
          search_interval = search_interval,
          prior_use = cohort_transition_prior_use,
          sd_floor = max(cohort_transition_sd_floor, cohort_transition_patient_sd_floor)
        )
        if (is.finite(fit_res$f_map)) {
          fc[child_name] <- fit_res$f_map
        }
        node_rows[[child_name]] <- fit_res$diagnostics
      }
      nn_cohort_transition_node_diagnostics <- do.call(rbind, node_rows)
      if (!is.null(nn_cohort_transition_node_diagnostics) && nrow(nn_cohort_transition_node_diagnostics)) {
        nn_cohort_transition_node_diagnostics$replicate_id <- as.integer(b_iter_idx)
        nn_cohort_transition_node_diagnostics <- nn_cohort_transition_node_diagnostics[
          c("replicate_id", setdiff(names(nn_cohort_transition_node_diagnostics), "replicate_id"))
        ]
        rownames(nn_cohort_transition_node_diagnostics) <- NULL
      } else {
        nn_cohort_transition_node_diagnostics <- data.frame()
      }
      nn_prior_diag$nn_prior_mode_used <- "cohort_transition"
      nn_prior_diag$nn_prior_source_used <- if (isTRUE(cohort_transition_prior_use$leave_one_patient_out)) {
        "cohort_transition_leave_one_patient_out"
      } else {
        "cohort_transition_full"
      }
      nn_prior_diag$prior_mu_hat <- cohort_transition_prior_use$global_prior$mu[1]
      nn_prior_diag$prior_sigma_hat <- cohort_transition_prior_use$global_prior$sigma_with_patient_heterogeneity[1]
      nn_prior_diag$n_zero_children_retained <- as.integer(sum(!nn_present))
      nn_prior_diag$sum_zero_weight_final <- as.numeric(sum(!nn_present))
    } else if (use_empirical_prior && any(!nn_present) && length(fc_prior_vals) > 0 && !all(is.na(fc_prior_vals))) {
      mean_fc_prior_val <- mean(fc_prior_vals, na.rm = TRUE) # Renamed mean_fc_prior
      if (is.null(nn_prior_sd)) {
        sd_fc_prior_val <- stats::sd(fc_prior_vals, na.rm = TRUE)
        if (!is.finite(sd_fc_prior_val)) {
          sd_fc_prior_val <- nn_prior_sd_floor
        }
        sd_fc_prior_val <- max(sd_fc_prior_val, nn_prior_sd_floor)
      } else {
        sd_fc_prior_val <- nn_prior_sd
      }
      nn_prior_diag$nn_prior_mode_used <- "empirical"
      nn_prior_diag$sum_observed_weight <- length(fc_prior_vals)
      nn_prior_diag$prior_mu_hat <- mean_fc_prior_val
      nn_prior_diag$prior_sigma_hat <- sd_fc_prior_val

      sapply_names_not_present <- names(nn_child_contexts)[!nn_present]
      if(length(sapply_names_not_present) > 0) {
        for (child_name in sapply_names_not_present) {
          objective_fn <- build_opt_fc(nn_child_contexts[[child_name]],
                                       prior_mean_param = mean_fc_prior_val,
                                       prior_sd_param = sd_fc_prior_val,
                                       do_prior_param = TRUE)
          res <- run_optimise_checked(objective_fn, interval = search_interval,
                                      context = sprintf("optimise nearest-neighbour fitness with prior for latent child %s", child_name))
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
    } else if (use_empirical_censored_prior && any(!nn_present)) {
      prior_fit <- estimate_nn_prior_censored_eb(
        nn_info_items = nn_child_contexts,
        fpar = fpar,
        build_opt_fc = build_opt_fc,
        search_interval = search_interval,
        nn_prior_sd = nn_prior_sd,
        nn_prior_sd_floor = nn_prior_sd_floor,
        nn_prior_grid_n = nn_prior_grid_n,
        context = sprintf(
          "fit empirical_censored latent-neighbour prior for bootstrap replicate %d",
          b_iter_idx
        )
      )
      nn_prior_diag$nn_prior_mode_used <- "empirical_censored"
      nn_prior_diag$n_zero_children_retained <- as.integer(sum(!nn_present))
      nn_prior_diag$n_zero_children_screened <- 0L
      nn_prior_diag$sum_observed_weight <- sum(nn_present)
      nn_prior_diag$sum_zero_weight_raw <- sum(!nn_present)
      nn_prior_diag$sum_zero_weight_final <- sum(!nn_present)
      nn_prior_diag$prior_mu_hat <- prior_fit$prior_mean
      nn_prior_diag$prior_sigma_hat <- prior_fit$prior_sd
      if (!is.null(prior_fit$informative_child_count)) {
        nn_prior_diag$informative_child_count <- prior_fit$informative_child_count
      }
      if (!is.null(prior_fit$map_delta_lower_boundary_rate)) {
        nn_prior_diag$map_delta_lower_boundary_rate <- prior_fit$map_delta_lower_boundary_rate
      }
      if (!is.null(prior_fit$map_delta_upper_boundary_rate)) {
        nn_prior_diag$map_delta_upper_boundary_rate <- prior_fit$map_delta_upper_boundary_rate
      }

      sapply_names_not_present <- names(nn_child_contexts)[!nn_present]
      if (length(sapply_names_not_present) > 0) {
        for (child_name in sapply_names_not_present) {
          objective_fn <- build_opt_fc(nn_child_contexts[[child_name]],
                                       prior_mean_param = prior_fit$prior_mean,
                                       prior_sd_param = prior_fit$prior_sd,
                                       do_prior_param = TRUE)
          res <- run_optimise_strict_checked(
            objective_fn,
            interval = search_interval,
            context = sprintf(
              "optimise nearest-neighbour fitness with empirical_censored prior for latent child %s",
              child_name
            )
          )
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
	    } else if (use_weighted_like_prior &&
	               any(!nn_present) &&
	               isTRUE(weighted_prior_config$can_fit_replicate_prior)) {
	      prior_fit_result <- tryCatch(
	        list(
	          fit = estimate_nn_prior_censored_eb(
	            nn_info_items = weighted_prior_config$prior_items,
	            fpar = fpar,
	            build_opt_fc = build_opt_fc,
	            search_interval = search_interval,
	            nn_prior_sd = nn_prior_sd,
	            nn_prior_sd_floor = nn_prior_sd_floor,
	            nn_prior_grid_n = nn_prior_grid_n,
	            child_weights = weighted_prior_config$child_weights,
	            parent_mean_fn = weighted_prior_config$parent_mean_fn,
	            context = sprintf(
	              "fit empirical_censored_weighted latent-neighbour prior for bootstrap replicate %d",
	              b_iter_idx
	            )
	          ),
	          error = NULL
	        ),
	        error = function(e) {
	          if (isTRUE(run_two_shell_baseline)) {
	            list(fit = NULL, error = conditionMessage(e))
	          } else {
	            stop(e)
	          }
	        }
	      )
	      prior_fit <- prior_fit_result$fit
	      if (!is.null(prior_fit_result$error)) {
	        nn_prior_diag$fallback_reason <- "inward_prior_fit_failed"
	      }
	      if (!is.null(prior_fit)) {
	        nn_prior_diag$nn_prior_mode_used <- "empirical_censored_weighted"
	        nn_prior_diag$nn_prior_source_used <- "observed_replicate"
	        nn_prior_diag$prior_mu_hat <- prior_fit$prior_mean
	        nn_prior_diag$prior_sigma_hat <- prior_fit$prior_sd
	        two_shell_inward_prior_fit <- list(prior_mean = prior_fit$prior_mean, prior_sd = prior_fit$prior_sd)
	        if (!is.null(prior_fit$informative_child_count)) {
	          nn_prior_diag$informative_child_count <- prior_fit$informative_child_count
	        }
	        if (!is.null(prior_fit$map_delta_lower_boundary_rate)) {
	          nn_prior_diag$map_delta_lower_boundary_rate <- prior_fit$map_delta_lower_boundary_rate
	        }
	        if (!is.null(prior_fit$map_delta_upper_boundary_rate)) {
	          nn_prior_diag$map_delta_upper_boundary_rate <- prior_fit$map_delta_upper_boundary_rate
	        }

	        sapply_names_not_present <- names(nn_child_contexts)[!nn_present]
	        if (length(sapply_names_not_present) > 0) {
	          for (child_name in sapply_names_not_present) {
	            objective_fn <- build_opt_fc(nn_child_contexts[[child_name]],
	                                         prior_mean_param = prior_fit$prior_mean,
	                                         prior_sd_param = prior_fit$prior_sd,
	                                         do_prior_param = TRUE,
	                                         parent_mean_mode = "exposure")
	            res <- run_optimise_strict_checked(
	              objective_fn,
	              interval = search_interval,
	              context = sprintf(
	                "optimise nearest-neighbour fitness with empirical_censored_weighted prior for latent child %s",
	                child_name
	              )
	            )
	            if (!is.null(res)) {
	              fc[child_name] <- res$minimum
	            }
	          }
	        }
	      } else {
	        nn_prior_diag$nn_prior_mode_used <- "none"
	        nn_prior_diag$nn_prior_source_used <- "none"
	        sapply_names_not_present <- names(nn_child_contexts)[!nn_present]
	        if (length(sapply_names_not_present) > 0) {
	          for (child_name in sapply_names_not_present) {
	            objective_fn <- build_opt_fc(nn_child_contexts[[child_name]], do_prior_param = FALSE)
	            res <- run_optimise_checked(
	              objective_fn,
	              interval = search_interval,
	              context = sprintf("optimise nearest-neighbour fitness without prior for latent child %s", child_name)
	            )
	            if (!is.null(res)) {
	              fc[child_name] <- res$minimum
	            }
	          }
	        }
	      }
	    } else if (use_weighted_like_prior &&
	               any(!nn_present) &&
	               !is.null(weighted_sample_pooled_prior_use) &&
	               isTRUE(weighted_sample_pooled_prior_use$available)) {
      nn_prior_diag$nn_prior_mode_used <- "empirical_censored_weighted"
	      nn_prior_diag$nn_prior_source_used <- "sample_pooled"
	      nn_prior_diag$prior_mu_hat <- weighted_sample_pooled_prior_use$prior_mean
	      nn_prior_diag$prior_sigma_hat <- weighted_sample_pooled_prior_use$prior_sd
	      two_shell_inward_prior_fit <- list(
	        prior_mean = weighted_sample_pooled_prior_use$prior_mean,
	        prior_sd = weighted_sample_pooled_prior_use$prior_sd
	      )
      sapply_names_not_present <- names(nn_child_contexts)[!nn_present]
      if (length(sapply_names_not_present) > 0) {
        for (child_name in sapply_names_not_present) {
          objective_fn <- build_opt_fc(
            nn_child_contexts[[child_name]],
            prior_mean_param = weighted_sample_pooled_prior_use$prior_mean,
            prior_sd_param = weighted_sample_pooled_prior_use$prior_sd,
            do_prior_param = TRUE,
            parent_mean_mode = "exposure"
          )
          res <- run_optimise_strict_checked(
            objective_fn,
            interval = search_interval,
            context = sprintf(
              "optimise nearest-neighbour fitness with sample-pooled empirical_censored_weighted prior for latent child %s",
              child_name
            )
          )
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
	    } else if (any(!nn_present)) { # No prior to use
	      nn_prior_diag$nn_prior_mode_used <- "none"
	      if (use_weighted_like_prior) {
	        nn_prior_diag$nn_prior_source_used <- "none"
	      }
      sapply_names_not_present <- names(nn_child_contexts)[!nn_present]
      if(length(sapply_names_not_present) > 0) {
        for (child_name in sapply_names_not_present) {
          objective_fn <- build_opt_fc(nn_child_contexts[[child_name]], do_prior_param = FALSE)
          res <- run_optimise_checked(objective_fn, interval = search_interval,
                                      context = sprintf("optimise nearest-neighbour fitness without prior for latent child %s", child_name))
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
	        }
	      }
	    }

	    nn_two_shell_node_diagnostics <- data.frame()
	    if (isTRUE(run_two_shell_baseline)) {
	      if (is.null(two_shell_inward_prior_fit) ||
	          !is.finite(two_shell_inward_prior_fit$prior_mean) ||
	          !is.finite(two_shell_inward_prior_fit$prior_sd) ||
	          two_shell_inward_prior_fit$prior_sd <= 0) {
	        shell01_fallback <- fit_shell_delta_prior(
	          delta = fc_prior_vals,
	          sd_floor = nn_prior_sd_floor,
	          fallback_mu = 0,
	          fallback_sd = if (is.null(nn_prior_sd)) nn_prior_sd_floor else nn_prior_sd
	        )
	        two_shell_inward_prior_fit <- list(
	          prior_mean = shell01_fallback$mu,
	          prior_sd = shell01_fallback$sigma
	        )
	        if (!is.finite(nn_prior_diag$prior_mu_hat)) {
	          nn_prior_diag$prior_mu_hat <- two_shell_inward_prior_fit$prior_mean
	        }
	        if (!is.finite(nn_prior_diag$prior_sigma_hat)) {
	          nn_prior_diag$prior_sigma_hat <- two_shell_inward_prior_fit$prior_sd
	        }
	      }
	      two_shell_res <- run_empirical_two_shell_correction(
	        nn_child_contexts = nn_child_contexts,
	        nn_present = nn_present,
	        f1_initial = fc,
	        fpar = fpar,
	        boot_data = boot_data,
	        timepoints = current_timepoints,
	        ntot = ntot_rounded,
	        pm = pm,
	        n0 = current_n0,
	        search_interval = search_interval,
	        inward_prior_fit = two_shell_inward_prior_fit,
	        nn_prior_sd_floor = nn_prior_sd_floor,
	        nn_two_shell_min_delta_n = nn_two_shell_min_delta_n,
	        nn_two_shell_min_exposure = nn_two_shell_min_exposure,
	        nn_two_shell_min_observed_count = nn_two_shell_min_observed_count,
	        nn_two_shell_max_weight_ratio = nn_two_shell_max_weight_ratio,
	        nn_two_shell_lambda = nn_two_shell_lambda,
	        nn_two_shell_reuse_sd = nn_two_shell_reuse_sd,
	        nn_two_shell_uncertainty_floor = nn_two_shell_uncertainty_floor
	      )
	      fc <- two_shell_res$f1
	      nn_two_shell_node_diagnostics <- two_shell_res$node_diagnostics
	      if (nrow(nn_two_shell_node_diagnostics)) {
	        nn_two_shell_node_diagnostics$replicate_id <- as.integer(b_iter_idx)
	        nn_two_shell_node_diagnostics <- nn_two_shell_node_diagnostics[
	          c("replicate_id", setdiff(names(nn_two_shell_node_diagnostics), "replicate_id"))
	        ]
	        rownames(nn_two_shell_node_diagnostics) <- NULL
	      }
	      nn_prior_diag[names(two_shell_res$diagnostics)] <- two_shell_res$diagnostics
	      if (is.na(two_shell_res$diagnostics$fallback_reason) ||
	          !nzchar(two_shell_res$diagnostics$fallback_reason)) {
	        nn_prior_diag$nn_prior_mode_used <- if (use_cohort_transition_prior) "cohort_transition" else "empirical_two_shell"
	        nn_prior_diag$nn_prior_source_used <- if (use_cohort_transition_prior) "cohort_transition_two_shell_baseline" else "two_shell"
	      } else if (identical(nn_prior_diag$nn_prior_mode_used, "empirical_two_shell") ||
                   (use_cohort_transition_prior && identical(nn_prior_diag$nn_prior_mode_used, "cohort_transition"))) {
	        nn_prior_diag$nn_prior_mode_used <- "empirical_censored_weighted"
	        nn_prior_diag$nn_prior_source_used <- "fallback_inward"
	      }
	      nn_prior_diag$mu01 <- two_shell_res$diagnostics$mu01
	      nn_prior_diag$sigma01 <- two_shell_res$diagnostics$sigma01
	      nn_prior_diag$prior_mu_hat <- two_shell_res$diagnostics$mu01
	      nn_prior_diag$prior_sigma_hat <- two_shell_res$diagnostics$sigma01
	    }

	    if (use_cohort_transition_prior && identical(cohort_transition_version, "v2") && length(nn_child_contexts) > 0) {
	      node_rows <- vector("list", length(nn_child_contexts))
	      names(node_rows) <- names(nn_child_contexts)
	      for (child_name in names(nn_child_contexts)) {
	        two_shell_row <- nn_two_shell_node_diagnostics[FALSE, , drop = FALSE]
	        if (nrow(nn_two_shell_node_diagnostics) && "karyotype" %in% names(nn_two_shell_node_diagnostics)) {
	          two_shell_row <- nn_two_shell_node_diagnostics[nn_two_shell_node_diagnostics$karyotype == child_name, , drop = FALSE]
	          if (nrow(two_shell_row) > 1L) two_shell_row <- two_shell_row[1L, , drop = FALSE]
	        }
	        overlay_res <- apply_cohort_transition_overlay(
	          item = nn_child_contexts[[child_name]],
	          child_name = child_name,
	          build_opt_fc = build_opt_fc,
	          search_interval = search_interval,
	          prior_use = cohort_transition_prior_use,
	          f_two_shell_baseline = fc[child_name],
	          nn_present = nn_present[child_name],
	          two_shell_node_diagnostics = two_shell_row,
	          cohort_transition_apply_to = cohort_transition_apply_to,
	          cohort_transition_lambda = cohort_transition_lambda,
	          cohort_transition_max_borrowing_fraction = cohort_transition_max_borrowing_fraction,
	          cohort_transition_max_abs_delta_shift = cohort_transition_max_abs_delta_shift,
	          cohort_transition_sd_floor = cohort_transition_sd_floor,
	          cohort_transition_patient_sd_floor = cohort_transition_patient_sd_floor
	        )
	        if (is.finite(overlay_res$f_final) &&
              nrow(overlay_res$diagnostics) &&
              any(overlay_res$diagnostics$cohort_update_applied, na.rm = TRUE)) {
	          fc[child_name] <- overlay_res$f_final
	        }
	        node_rows[[child_name]] <- overlay_res$diagnostics
	      }
	      nn_cohort_transition_node_diagnostics <- do.call(rbind, node_rows)
	      if (!is.null(nn_cohort_transition_node_diagnostics) && nrow(nn_cohort_transition_node_diagnostics)) {
	        nn_cohort_transition_node_diagnostics$replicate_id <- as.integer(b_iter_idx)
	        nn_cohort_transition_node_diagnostics <- nn_cohort_transition_node_diagnostics[
	          c("replicate_id", setdiff(names(nn_cohort_transition_node_diagnostics), "replicate_id"))
	        ]
	        rownames(nn_cohort_transition_node_diagnostics) <- NULL
	      } else {
	        nn_cohort_transition_node_diagnostics <- data.frame()
	      }
	      nn_prior_diag$nn_prior_mode_used <- "cohort_transition"
	      nn_prior_diag$nn_prior_source_used <- if (isTRUE(cohort_transition_prior_use$leave_one_patient_out)) {
	        "cohort_transition_v2_leave_one_patient_out_overlay"
	      } else {
	        "cohort_transition_v2_full_overlay"
	      }
	      nn_prior_diag$cohort_transition_version <- "v2"
	      nn_prior_diag$cohort_transition_apply_to <- cohort_transition_apply_to
	      nn_prior_diag$cohort_transition_overlay_base <- cohort_transition_overlay_base
	      nn_prior_diag$cohort_transition_lambda <- cohort_transition_lambda
	      nn_prior_diag$n_cohort_overlay_nodes_updated <- sum(
	        nn_cohort_transition_node_diagnostics$cohort_update_applied,
	        na.rm = TRUE
	      )
	      nn_prior_diag$prior_mu_hat <- cohort_transition_prior_use$global_prior$mu[1]
	      nn_prior_diag$prior_sigma_hat <- if ("effective_prior_sd" %in% names(cohort_transition_prior_use$global_prior)) {
	        cohort_transition_prior_use$global_prior$effective_prior_sd[1]
	      } else {
	        cohort_transition_prior_use$global_prior$sigma_with_patient_heterogeneity[1]
	      }
	      nn_prior_diag$n_zero_children_retained <- as.integer(sum(!nn_present))
	      nn_prior_diag$sum_zero_weight_final <- as.numeric(sum(!nn_present))
	    }

	    if (use_cohort_transition_prior && identical(cohort_transition_version, "contextual") && length(nn_child_contexts) > 0) {
	      node_rows <- vector("list", length(nn_child_contexts))
	      names(node_rows) <- names(nn_child_contexts)
	      for (child_name in names(nn_child_contexts)) {
	        two_shell_row <- nn_two_shell_node_diagnostics[FALSE, , drop = FALSE]
	        if (nrow(nn_two_shell_node_diagnostics) && "karyotype" %in% names(nn_two_shell_node_diagnostics)) {
	          two_shell_row <- nn_two_shell_node_diagnostics[nn_two_shell_node_diagnostics$karyotype == child_name, , drop = FALSE]
	          if (nrow(two_shell_row) > 1L) two_shell_row <- two_shell_row[1L, , drop = FALSE]
	        }
	        overlay_res <- apply_contextual_cohort_overlay(
	          item = nn_child_contexts[[child_name]],
	          child_name = child_name,
	          build_opt_fc = build_opt_fc,
	          search_interval = search_interval,
	          prior_use = cohort_transition_prior_use,
	          f_two_shell_baseline = fc[child_name],
	          nn_present = nn_present[child_name],
	          two_shell_node_diagnostics = two_shell_row,
	          cohort_contextual_apply_to = cohort_contextual_apply_to,
	          cohort_context_lambda = cohort_context_lambda,
	          cohort_context_max_borrowing_fraction = cohort_context_max_borrowing_fraction,
	          cohort_context_max_abs_delta_shift = cohort_context_max_abs_delta_shift,
	          cohort_context_sd_floor = cohort_context_sd_floor,
	          cohort_context_patient_sd_floor = cohort_context_patient_sd_floor,
	          cohort_context_keep_baseline_when_sparse = cohort_context_keep_baseline_when_sparse,
	          cohort_context_keep_baseline_when_high_variable = cohort_context_keep_baseline_when_high_variable
	        )
	        if (is.finite(overlay_res$f_final) &&
              nrow(overlay_res$diagnostics) &&
              any(overlay_res$diagnostics$cohort_update_applied, na.rm = TRUE)) {
	          fc[child_name] <- overlay_res$f_final
	        }
	        node_rows[[child_name]] <- overlay_res$diagnostics
	      }
	      nn_cohort_transition_node_diagnostics <- do.call(rbind, node_rows)
	      if (!is.null(nn_cohort_transition_node_diagnostics) && nrow(nn_cohort_transition_node_diagnostics)) {
	        nn_cohort_transition_node_diagnostics$replicate_id <- as.integer(b_iter_idx)
	        nn_cohort_transition_node_diagnostics <- nn_cohort_transition_node_diagnostics[
	          c("replicate_id", setdiff(names(nn_cohort_transition_node_diagnostics), "replicate_id"))
	        ]
	        rownames(nn_cohort_transition_node_diagnostics) <- NULL
	      } else {
	        nn_cohort_transition_node_diagnostics <- data.frame()
	      }
	      nn_prior_diag$nn_prior_mode_used <- "cohort_transition"
	      nn_prior_diag$nn_prior_source_used <- if (isTRUE(cohort_transition_prior_use$leave_one_patient_out)) {
	        "cohort_transition_contextual_leave_one_patient_out_overlay"
	      } else {
	        "cohort_transition_contextual_full_overlay"
	      }
	      nn_prior_diag$cohort_transition_version <- "contextual"
	      nn_prior_diag$cohort_contextual_apply_to <- cohort_contextual_apply_to
	      nn_prior_diag$cohort_contextual_overlay_base <- cohort_contextual_overlay_base
	      nn_prior_diag$cohort_context_lambda <- cohort_context_lambda
	      nn_prior_diag$n_cohort_overlay_nodes_updated <- sum(
	        nn_cohort_transition_node_diagnostics$cohort_update_applied,
	        na.rm = TRUE
	      )
	      nn_prior_diag$prior_mu_hat <- NA_real_
	      nn_prior_diag$prior_sigma_hat <- NA_real_
	      nn_prior_diag$n_zero_children_retained <- as.integer(sum(!nn_present))
	      nn_prior_diag$sum_zero_weight_final <- as.numeric(sum(!nn_present))
	    }

	    list(f_initial = f_qp,
	         f_final = f_final,
	         x0_initial = x0_init,
	         x0_final = x0_final,
	         f_nn = fc,
	         nn_prior_diagnostics = nn_prior_diag,
	         nn_two_shell_node_diagnostics = nn_two_shell_node_diagnostics,
	         nn_cohort_transition_node_diagnostics = nn_cohort_transition_node_diagnostics)
  }

  # Run bootstrap iterations serially using lapply
  boot_list <- lapply(seq_len(nboot), bootstrap_iter,
                      current_data = data, current_fq = fq, current_timepoints = timepoints,
                      current_num_species = num_species, current_num_timepoints = num_timepoints,
                      current_epsilon = epsilon, current_n0 = n0, current_nb = nb,
                      current_viability = viability, current_nn_info = nn_info_list)

  # Consolidate results
  f_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_initial))
  f_final_mat   <- do.call(rbind, lapply(boot_list, function(x) x$f_final))
  x0_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$x0_initial))
  x0_final_mat  <- do.call(rbind, lapply(boot_list, function(x) x$x0_final))
  f_nn_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_nn))
  if (is.null(f_nn_mat)) {
    f_nn_mat <- matrix(numeric(0), nrow = length(boot_list), ncol = 0)
  }
	  nn_prior_diagnostics <- do.call(rbind, lapply(boot_list, function(x) {
	    as.data.frame(x$nn_prior_diagnostics, stringsAsFactors = FALSE)
	  }))
	  rownames(nn_prior_diagnostics) <- NULL
	  nn_two_shell_node_diagnostics <- do.call(rbind, lapply(boot_list, function(x) {
	    if (is.null(x$nn_two_shell_node_diagnostics) || !nrow(x$nn_two_shell_node_diagnostics)) {
	      return(NULL)
	    }
	    x$nn_two_shell_node_diagnostics
	  }))
	  if (is.null(nn_two_shell_node_diagnostics)) {
	    nn_two_shell_node_diagnostics <- data.frame()
	  } else {
	    rownames(nn_two_shell_node_diagnostics) <- NULL
	  }
	  nn_cohort_transition_node_diagnostics <- do.call(rbind, lapply(boot_list, function(x) {
	    if (is.null(x$nn_cohort_transition_node_diagnostics) || !nrow(x$nn_cohort_transition_node_diagnostics)) {
	      return(NULL)
	    }
	    x$nn_cohort_transition_node_diagnostics
	  }))
	  if (is.null(nn_cohort_transition_node_diagnostics)) {
	    nn_cohort_transition_node_diagnostics <- data.frame()
	  } else {
	    rownames(nn_cohort_transition_node_diagnostics) <- NULL
	  }

  # Set column names if matrices are not empty and fq/names(nn_info_list) are not empty
  if(length(fq) > 0) {
    if(nrow(f_initial_mat) > 0) colnames(f_initial_mat) <- fq
    if(nrow(f_final_mat) > 0) colnames(f_final_mat) <- fq
    if(nrow(x0_initial_mat) > 0) colnames(x0_initial_mat) <- fq
    if(nrow(x0_final_mat) > 0) colnames(x0_final_mat) <- fq
  }
  if(length(names(nn_info_list)) > 0 && nrow(f_nn_mat) > 0) {
    colnames(f_nn_mat) <- names(nn_info_list)
  }

  list(initial_fitness = f_initial_mat,
       final_fitness = f_final_mat,
	       initial_frequencies = x0_initial_mat,
	       final_frequencies = x0_final_mat,
	       nn_fitness = f_nn_mat,
	       nn_prior_diagnostics = nn_prior_diagnostics,
	       nn_two_shell_node_diagnostics = nn_two_shell_node_diagnostics,
	       nn_cohort_transition_node_diagnostics = nn_cohort_transition_node_diagnostics)
	}

#' Fit Kriging model to fitness data (Internal function)
#' @keywords internal
#' @noRd
fitKrig <- function(fq_boot, nboot, krig_bootstrap_mode = c("marginal", "joint")) {
  validate_positive_integer(nboot, "nboot")
  krig_bootstrap_mode <- validate_krig_bootstrap_mode(krig_bootstrap_mode)
  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness) # Will be NULL if nn_fitness is NULL or has no colnames

  # Handle cases where fq_str or nn_str might be NULL (e.g., if fq_boot$final_fitness is NULL)
  valid_fq_str <- if(is.null(fq_str)) character(0) else fq_str
  valid_nn_str <- if(is.null(nn_str)) character(0) else nn_str

  combined_strs <- c(valid_fq_str, valid_nn_str)
  if(length(combined_strs) == 0 || ncol(fboot) == 0) { # No data to train on
    warning("fitKrig: No valid fitness data (fq_str or nn_str) to train Kriging model.")
    empty_df <- data.frame(k = character(0), mean = numeric(0), median = numeric(0), sd = numeric(0),
                           fq = logical(0), nn = logical(0))
    return(list(summary_stats = empty_df,
                posterior_samples = matrix(numeric(0), ncol=0, nrow=0),
                boot_results = list(),
                fit_boot_list = list(),
                krig_stable_mean = NULL,
                krig_stable_median = NULL))
  }
  ktrain <- unname(parse_karyotype_ids(combined_strs))

  # Ensure nn_str is not NULL before passing to gen_all_neighbours
  ktest_neighbours_matrix <- matrix(numeric(0), ncol=ncol(ktrain)) # empty matrix with correct cols
  if (length(valid_nn_str) > 0 && !is.null(valid_nn_str)) {
    ktest_neighbours_matrix <- gen_all_neighbours(valid_nn_str)
  }

  ktest <- unique(rbind(ktest = ktrain, ktest_neighbours_matrix))
  ktest_str <- apply(ktest, 1, paste, collapse = ".")
  fq_ids <- ktest_str %in% valid_fq_str
  nn_ids <- ktest_str %in% valid_nn_str

  fboot_mean <- colMeans(fboot, na.rm = TRUE)
  fboot_median <- apply(fboot, 2, stats::median, na.rm = TRUE)
  if (!is.null(names(fboot_mean))) {
    fboot_mean <- fboot_mean[combined_strs]
    fboot_median <- fboot_median[combined_strs]
  }
  valid_mean <- is.finite(fboot_mean)
  valid_median <- is.finite(fboot_median)
  krig_stable_mean <- NULL
  krig_stable_median <- NULL
  if (sum(valid_mean) >= 2 && length(unique(fboot_mean[valid_mean])) >= 2) {
    krig_stable_mean <- fields::Krig(
      ktrain[valid_mean, , drop = FALSE],
      fboot_mean[valid_mean],
      cov.function = "stationary.cov",
      cov.args = krig_covariance_args(),
      nstep.cv = ALFAK_KRIG_NSTEP_CV,
      give.warnings = TRUE
    )
  } else {
    warning("fitKrig: Insufficient data for stable mean Kriging fit.")
  }
  if (sum(valid_median) >= 2 && length(unique(fboot_median[valid_median])) >= 2) {
    krig_stable_median <- fields::Krig(
      ktrain[valid_median, , drop = FALSE],
      fboot_median[valid_median],
      cov.function = "stationary.cov",
      cov.args = krig_covariance_args(),
      nstep.cv = ALFAK_KRIG_NSTEP_CV,
      give.warnings = TRUE
    )
  } else {
    warning("fitKrig: Insufficient data for stable median Kriging fit.")
  }

  # Use lapply directly, as cl is removed
  boot_predictions_list <- lapply(seq_len(nboot), function(b) {
    boot_f <- if (krig_bootstrap_mode == "joint") {
      as.numeric(fboot[sample(seq_len(nrow(fboot)), 1), ])
    } else {
      boot_f_indices <- cbind(sample(seq_len(nrow(fboot)), ncol(fboot), replace = TRUE), seq_len(ncol(fboot)))
      as.vector(fboot[boot_f_indices])
    }

    valid_boot <- is.finite(boot_f)
    ktrain_boot <- ktrain[valid_boot, , drop = FALSE]
    boot_f_valid <- boot_f[valid_boot]

    if(nrow(ktrain_boot) < 2 ||
       length(boot_f_valid) < 2 ||
       nrow(unique(ktrain_boot)) < 2 ||
       length(unique(boot_f_valid)) < 2) {
      stop("fitKrig: Insufficient or incompatible data for Kriging in bootstrap iteration.")
    }

    fit_boot <- fields::Krig(
      ktrain_boot,
      boot_f_valid,
      cov.function = "stationary.cov",
      cov.args = krig_covariance_args(),
      nstep.cv = ALFAK_KRIG_NSTEP_CV,
      give.warnings = TRUE
    )
    preds <- stats::predict(fit_boot, ktest)

    list(fit_boot = fit_boot, preds = preds)
  })

  boot_predictions <- do.call(cbind, lapply(boot_predictions_list, `[[`, "preds"))
  fit_boot_list   <- lapply(boot_predictions_list, `[[`, "fit_boot")
  pred_means <- apply(boot_predictions, 1, mean_or_na)
  pred_medians <- apply(boot_predictions, 1, median_or_na)
  pred_sd <- apply(boot_predictions, 1, sd_or_na)

  summary_df <- data.frame(k = ktest_str, mean = pred_means, median = pred_medians, sd = pred_sd,
                           fq = fq_ids, nn = nn_ids)

  #list(summary_stats = summary_df, posterior_samples = boot_predictions)
  return(list(
    summary_stats     = summary_df,
    posterior_samples = boot_predictions,
    boot_results      = boot_predictions_list,
    fit_boot_list     = fit_boot_list,
    krig_stable_mean  = krig_stable_mean,
    krig_stable_median = krig_stable_median
  ))
}

#' Cross-validation for Kriging model (Internal function)
#' @keywords internal
#' @noRd
xval <- function(fq_boot, krig_bootstrap_mode = c("marginal", "joint")) {
  krig_bootstrap_mode <- validate_krig_bootstrap_mode(krig_bootstrap_mode)
  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness) # Can be NULL

  valid_fq_str <- if(is.null(fq_str)) character(0) else fq_str
  valid_nn_str <- if(is.null(nn_str)) character(0) else nn_str

  combined_strs <- c(valid_fq_str, valid_nn_str)
  if(length(combined_strs) == 0 || ncol(fboot) == 0 || nrow(fboot) == 0) {
    warning("xval: No valid fitness data to perform cross-validation.")
    return(NA_real_)
  }
  ktrain <- unname(parse_karyotype_ids(combined_strs))

  # Original ids logic for xval
  if(length(valid_fq_str) == 0) {
    warning("xval: No fq_str defined, cannot perform original xval logic based on fq_str.")
    return(NA_real_)
  }

  ids <- unlist(lapply(seq_along(valid_fq_str), function(i_xval) {
    ki_neighbours_matrix <- gen_all_neighbours(valid_fq_str[i_xval])
    ki_neighbours_str <- character(0)
    if(nrow(ki_neighbours_matrix) > 0) {
      ki_neighbours_str <- as.character(apply(ki_neighbours_matrix, 1, paste, collapse = "."))
    }
    ki <- c(valid_fq_str[i_xval], ki_neighbours_str)
    idi <- rep(i_xval, length(ki))
    names(idi) <- ki
    idi
  }))
  # Fold ownership keeps the first assignment for overlapping neighbourhoods, which preserves
  # the long-standing heuristic while making the order-dependence explicit.
  ids <- ids[!duplicated(names(ids))] # Ensure unique names for ids
  uids <- unique(ids) # These are fold identifiers

  # Use lapply directly
  tmp_list <- lapply(uids, function(id_fold) { # Renamed id to id_fold
    fi <- if (krig_bootstrap_mode == "joint") {
      fboot[sample(seq_len(nrow(fboot)), 1), ]
    } else {
      fboot_shuffled <- apply(fboot, 2, sample)
      fboot_shuffled[1, ]
    }
    # Original logic for train/test split based on 'ids' and current 'id_fold'
    train_indices <- !(ids == id_fold)
    test_indices <- (ids == id_fold)

    # Check if ktrain (all possible karyotypes from combined_strs) aligns with names(ids)
    # names(ids) are the karyotype strings that were assigned a fold id
    # We need to map these back to rows in ktrain if ktrain corresponds to combined_strs

    # Create a mapping from karyotype string to row index in ktrain
    ktrain_map <- setNames(1:nrow(ktrain), combined_strs)

    train_k_names <- names(ids)[train_indices]
    test_k_names <- names(ids)[test_indices]

    # Ensure these names exist in ktrain_map
    train_k_names_valid <- train_k_names[train_k_names %in% names(ktrain_map)]
    test_k_names_valid <- test_k_names[test_k_names %in% names(ktrain_map)]

    if(length(train_k_names_valid) == 0 || length(test_k_names_valid) == 0) {
      return(matrix(NA, ncol=2, dimnames=list(NULL, c("test_f", "est_f")))) # Skip fold
    }

    train_k_rows <- ktrain_map[train_k_names_valid]
    test_k_rows <- ktrain_map[test_k_names_valid]

    train_k <- ktrain[train_k_rows, , drop = FALSE]
    # fi corresponds to combined_strs (colnames of fboot)
    # So, we need to subset fi based on names corresponding to train_k_names_valid
    train_f <- fi[train_k_names_valid]

    test_k <- ktrain[test_k_rows, , drop = FALSE]
    test_f <- fi[test_k_names_valid]

    # Filter NAs from training data, which could arise if some fi values were NA
    valid_train_points <- !is.na(train_f)
    train_k <- train_k[valid_train_points, , drop=FALSE]
    train_f <- train_f[valid_train_points]

    if(nrow(train_k) < 2 || nrow(unique(train_k)) < 2 || length(unique(train_f)) < 1) {
      warning(paste("Skipping xval fold for id_fold:", id_fold, "due to insufficient/non-unique training points after NA removal."))
      return(cbind(test_f = test_f, est_f = rep(NA, length(test_f))))
    }

    fit <- fields::Krig(
      train_k,
      train_f,
      cov.function = "stationary.cov",
      cov.args = krig_covariance_args(),
      nstep.cv = ALFAK_KRIG_NSTEP_CV,
      give.warnings = TRUE
    )
    est_f <- stats::predict(fit, test_k)
    cbind(test_f, est_f)
  })

  tmp <- do.call(rbind, tmp_list)
  tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE] # Use stats::complete.cases

  if(nrow(tmp) < 2) { # R2R needs at least 2 points
    warning("Not enough valid observations after cross-validation to compute R2R.")
    return(NA_real_)
  }
  R2R(tmp[, 1], tmp[, 2])
}
