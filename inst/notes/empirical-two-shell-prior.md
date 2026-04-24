# Empirical Two-Shell Nearest-Neighbour Prior

`nn_prior = "empirical_two_shell"` is an opt-in correction for sparse
two-timepoint nearest-neighbour fitting. It keeps the existing frequent
karyotype fit and one-step nearest-neighbour likelihood, then performs one
additional backward correction using supported two-step descendants.

The implementation deliberately does not alternate one-step and two-step
updates. Each bootstrap replicate:

1. Fits frequent karyotypes and provisional one-step neighbours using the
   existing weighted censored machinery where possible.
2. Builds two-step candidates reachable from provisional one-step nodes.
   Observed candidates are retained automatically once they pass
   `nn_two_shell_min_observed_count`; unobserved candidates must pass
   `nn_two_shell_min_exposure`, or an adaptive projected-exposure threshold
   when that parameter is `NULL`.
3. Estimates provisional two-step fitness values and local uncertainty from
   the existing neighbour likelihood.
4. Learns separate 0->1 and 1->2 Gaussian delta priors with
   `nn_prior_sd_floor` applied to both shell-specific standard deviations.
5. Re-estimates one-step neighbours once with direct observation likelihood,
   inward frequent-parent prior terms, and outward two-step prior terms.

The outward term uses path responsibilities so a shared two-step descendant is
not counted once for each possible one-step path. Its effective standard
deviation includes the learned 1->2 prior scale, the provisional two-step
fitness variance, and `nn_two_shell_reuse_sd` because the two-step estimate
reuses the same bootstrap count data.

Replicates fall back to the initialized one-step estimates when no two-step
candidate is retained, too few usable 1->2 deltas are available, provisional
two-step fitting fails, or all outward weights are zero. `alfak()` saves
`nn_prior_diagnostics.Rds` for this mode when
`nn_two_shell_save_diagnostics = TRUE`; the same diagnostics are also carried
inside `bootstrap_res.Rds`.
