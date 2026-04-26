// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <map>
#include <vector>

namespace {

double log_sum_exp_cpp(const std::vector<double>& values) {
  if (values.empty()) {
    Rcpp::stop("log_sum_exp_cpp requires at least one value.");
  }
  bool has_finite_term = false;
  double max_val = R_NegInf;
  for (double value : values) {
    if (std::isnan(value) || value == R_PosInf) {
      Rcpp::stop("log_sum_exp_cpp rejects NaN and +Inf inputs.");
    }
    if (value == R_NegInf) {
      continue;
    }
    has_finite_term = true;
    if (value > max_val) {
      max_val = value;
    }
  }
  if (!has_finite_term) {
    Rcpp::stop("log_sum_exp_cpp cannot normalize an all -Inf vector.");
  }
  double accum = 0.0;
  for (double value : values) {
    if (value == R_NegInf) {
      continue;
    }
    accum += std::exp(value - max_val);
  }
  return max_val + std::log(accum);
}

double fexp_stable_cpp(double fc, double fp, double pij_value, double tt, double tol) {
  double delta = fc - fp;
  if (std::abs(delta) < tol) {
    return pij_value * fp * tt;
  }
  return pij_value * fp * std::expm1(tt * delta) / delta;
}

bool is_integer_valued_scalar(double x) {
  return std::floor(x) == x;
}

void validate_neighbor_projection_inputs(const Rcpp::NumericVector& parent_fitness,
                                         const Rcpp::NumericVector& pij_values,
                                         const Rcpp::NumericVector& parent_birth_times,
                                         const Rcpp::NumericVector& timepoints,
                                         const Rcpp::NumericMatrix& parent_xfit,
                                         double tol) {
  const int n_parents = parent_fitness.size();
  const int n_time = timepoints.size();
  if (!R_finite(tol) || tol <= 0.0) {
    Rcpp::stop("`tol` must be a positive finite value.");
  }
  if (pij_values.size() != n_parents || parent_birth_times.size() != n_parents ||
      parent_xfit.nrow() != n_parents || parent_xfit.ncol() != n_time) {
    Rcpp::stop("Parent inputs must have matching lengths/rows and parent_xfit columns must match timepoints.");
  }
  for (int p = 0; p < n_parents; ++p) {
    if (!R_finite(parent_fitness[p]) || !R_finite(pij_values[p]) || pij_values[p] < 0.0 ||
        !R_finite(parent_birth_times[p])) {
      Rcpp::stop("Parent fitness, transition probabilities, and birth times must be finite; pij values must be non-negative.");
    }
  }
  for (int t = 0; t < n_time; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
    for (int p = 0; p < n_parents; ++p) {
      if (!R_finite(parent_xfit(p, t))) {
        Rcpp::stop("`parent_xfit` must contain only finite values.");
      }
    }
  }
}

double finite_positive_sum(const Rcpp::NumericVector& x) {
  double total = 0.0;
  for (double v : x) {
    if (!R_finite(v) || v < 0.0) {
      return R_NaReal;
    }
    total += v;
  }
  return total;
}

double weighted_mean_or_nan(const Rcpp::NumericVector& x, const Rcpp::NumericVector& w) {
  if (x.size() != w.size()) {
    return R_NaReal;
  }
  double sw = 0.0;
  double sx = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i]) || !R_finite(w[i]) || w[i] < 0.0) {
      return R_NaReal;
    }
    sw += w[i];
    sx += w[i] * x[i];
  }
  if (!(sw > 0.0) || !R_finite(sw)) {
    return R_NaReal;
  }
  return sx / sw;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix alfak_project_forward_log_cpp(Rcpp::NumericVector x0,
                                                  Rcpp::NumericVector f,
                                                  Rcpp::NumericVector timepoints) {
  const int K = x0.size();
  const int T = timepoints.size();
  if (K == 0) {
    Rcpp::stop("`x0` must contain at least one entry.");
  }
  if (f.size() != K) {
    Rcpp::stop("`x0` and `f` must have the same length.");
  }
  double x0_sum = 0.0;
  for (int i = 0; i < K; ++i) {
    if (!R_finite(x0[i]) || x0[i] < 0.0) {
      Rcpp::stop("`x0` must contain only finite non-negative values.");
    }
    if (!R_finite(f[i])) {
      Rcpp::stop("`f` must contain only finite values.");
    }
    x0_sum += x0[i];
  }
  if (!(x0_sum > 0.0) || !R_finite(x0_sum)) {
    Rcpp::stop("`x0` must sum to a positive finite value.");
  }
  for (int t = 0; t < T; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
  }
  Rcpp::NumericMatrix out(K, T);
  std::vector<double> log_x0(K);
  std::vector<double> lv(K);

  for (int i = 0; i < K; ++i) {
    log_x0[i] = std::log(x0[i] / x0_sum);
  }

  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < K; ++i) {
      lv[i] = log_x0[i] + f[i] * timepoints[t];
    }
    double denom = log_sum_exp_cpp(lv);
    for (int i = 0; i < K; ++i) {
      out(i, t) = std::exp(lv[i] - denom);
    }
  }

  return out;
}

// [[Rcpp::export]]
double alfak_neg_log_lik_cpp(Rcpp::NumericVector param,
                             Rcpp::NumericMatrix counts,
                             Rcpp::NumericVector timepoints) {
  const int K = counts.nrow();
  const int T = counts.ncol();
  if (K <= 0) {
    Rcpp::stop("`counts` must have at least one row.");
  }
  if (T != timepoints.size()) {
    Rcpp::stop("`counts` must have ncol equal to length(timepoints).");
  }
  if (K == 1) {
    Rcpp::stop("`alfak_neg_log_lik_cpp()` expects at least two karyotypes; K == 1 should be handled in R.");
  }
  if (param.size() != (2 * K - 2)) {
    Rcpp::stop("`param` must have length 2*K - 2.");
  }
  std::vector<double> f_full(K, 0.0);
  std::vector<double> log_x0(K);
  std::vector<double> lv(K);

  double f_sum = 0.0;
  for (int i = 0; i < K - 1; ++i) {
    if (!R_finite(param[i])) {
      Rcpp::stop("`param` must contain only finite values.");
    }
    f_full[i] = param[i];
    f_sum += param[i];
  }
  f_full[K - 1] = -f_sum;

  for (int i = 0; i < K - 1; ++i) {
    if (!R_finite(param[K - 1 + i])) {
      Rcpp::stop("`param` must contain only finite values.");
    }
    log_x0[i] = param[K - 1 + i];
  }
  log_x0[K - 1] = 0.0;
  for (int t = 0; t < T; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
    for (int i = 0; i < K; ++i) {
      if (!R_finite(counts(i, t)) || counts(i, t) < 0.0) {
        Rcpp::stop("`counts` must contain only finite non-negative values.");
      }
    }
  }

  double nll = 0.0;
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < K; ++i) {
      lv[i] = log_x0[i] + f_full[i] * timepoints[t];
    }
    double denom = log_sum_exp_cpp(lv);
    for (int i = 0; i < K; ++i) {
      if (counts(i, t) > 0) {
        nll -= counts(i, t) * (lv[i] - denom);
      }
    }
  }

  return nll;
}

// [[Rcpp::export]]
double alfak_neighbor_objective_cpp(double fc_param,
                                    Rcpp::NumericVector parent_fitness,
                                    Rcpp::NumericVector pij_values,
                                    Rcpp::NumericVector parent_birth_times,
                                    Rcpp::NumericVector timepoints,
                                    Rcpp::NumericMatrix parent_xfit,
                                    Rcpp::NumericVector child_obs,
                                    Rcpp::NumericVector ntot,
                                    double parent_fitness_mean,
                                    double prior_mean,
                                    double prior_sd,
                                    bool do_prior,
                                    double tol) {
  const int n_parents = parent_fitness.size();
  const int n_time = timepoints.size();
  if (n_parents == 0) {
    return 1e9;
  }
  if (!R_finite(fc_param) || !R_finite(tol) || tol <= 0.0) {
    Rcpp::stop("`fc_param` must be finite and `tol` must be a positive finite value.");
  }
  if (pij_values.size() != n_parents || parent_birth_times.size() != n_parents ||
      parent_xfit.nrow() != n_parents) {
    Rcpp::stop("Parent inputs must have matching lengths/rows.");
  }
  if (parent_xfit.ncol() != n_time) {
    Rcpp::stop("`parent_xfit` must have ncol equal to length(timepoints).");
  }
  if (child_obs.size() != n_time || ntot.size() != n_time) {
    Rcpp::stop("`child_obs`, `ntot`, and `timepoints` must have matching lengths.");
  }
  if (do_prior && (!R_finite(prior_sd) || prior_sd <= 0.0 || !R_finite(prior_mean) || !R_finite(parent_fitness_mean))) {
    Rcpp::stop("When `do_prior` is TRUE, prior parameters and parent fitness mean must be finite and `prior_sd` must be positive.");
  }

  double loglik = 0.0;
  for (int t = 0; t < n_time; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
    if (!R_finite(child_obs[t]) || child_obs[t] < 0.0 || !is_integer_valued_scalar(child_obs[t])) {
      Rcpp::stop("`child_obs` must contain only finite non-negative integer-valued counts.");
    }
    if (!R_finite(ntot[t]) || ntot[t] < 0.0 || !is_integer_valued_scalar(ntot[t])) {
      Rcpp::stop("`ntot` must contain only finite non-negative integer-valued counts.");
    }
    if (child_obs[t] > ntot[t]) {
      Rcpp::stop("`child_obs` must not exceed `ntot` at any timepoint.");
    }
    double xc_est = 0.0;
    for (int p = 0; p < n_parents; ++p) {
      if (!R_finite(parent_fitness[p]) || !R_finite(pij_values[p]) || pij_values[p] < 0.0 ||
          !R_finite(parent_birth_times[p]) || !R_finite(parent_xfit(p, t))) {
        Rcpp::stop("Parent fitness, transition probabilities, birth times, and parent_xfit must be finite; pij values must be non-negative.");
      }
      double tt = std::max(0.0, timepoints[t] - parent_birth_times[p]);
      xc_est += fexp_stable_cpp(fc_param, parent_fitness[p], pij_values[p], tt, tol) * parent_xfit(p, t);
    }

    if (!R_finite(xc_est)) {
      loglik += -1e9;
      continue;
    }

    xc_est = std::max(0.0, std::min(1.0, xc_est));
    double ll = R::dbinom(child_obs[t], ntot[t], xc_est, true);
    if (!R_finite(ll)) {
      ll = -1e9;
    }
    loglik += ll;
  }

  if (do_prior && R_finite(parent_fitness_mean)) {
    double prior_ll = R::dnorm(fc_param - parent_fitness_mean, prior_mean, prior_sd, true);
    if (!R_finite(prior_ll)) {
      prior_ll = -1e9;
    }
    loglik += prior_ll;
  }

  return -loglik;
}

// [[Rcpp::export]]
Rcpp::NumericVector alfak_nn_project_trajectory_cpp(double fc_param,
                                                    Rcpp::NumericVector parent_fitness,
                                                    Rcpp::NumericVector pij_values,
                                                    Rcpp::NumericVector parent_birth_times,
                                                    Rcpp::NumericVector timepoints,
                                                    Rcpp::NumericMatrix parent_xfit,
                                                    double tol) {
  const int n_parents = parent_fitness.size();
  const int n_time = timepoints.size();
  if (!R_finite(fc_param)) {
    Rcpp::stop("`fc_param` must be finite.");
  }
  validate_neighbor_projection_inputs(parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit, tol);
  Rcpp::NumericVector projected(n_time);
  if (n_parents == 0) {
    return projected;
  }
  for (int t = 0; t < n_time; ++t) {
    double value = 0.0;
    for (int p = 0; p < n_parents; ++p) {
      double tt = std::max(0.0, timepoints[t] - parent_birth_times[p]);
      value += fexp_stable_cpp(fc_param, parent_fitness[p], pij_values[p], tt, tol) * parent_xfit(p, t);
    }
    if (!R_finite(value)) {
      value = 0.0;
    }
    projected[t] = std::max(0.0, std::min(1.0, value));
  }
  return projected;
}

// [[Rcpp::export]]
double alfak_nn_project_exposure_cpp(double fc_param,
                                     Rcpp::NumericVector parent_fitness,
                                     Rcpp::NumericVector pij_values,
                                     Rcpp::NumericVector parent_birth_times,
                                     Rcpp::NumericVector timepoints,
                                     Rcpp::NumericMatrix parent_xfit,
                                     Rcpp::NumericVector ntot,
                                     double tol) {
  const int n_time = timepoints.size();
  if (ntot.size() != n_time) {
    Rcpp::stop("`ntot` and `timepoints` must have matching lengths.");
  }
  Rcpp::NumericVector projected = alfak_nn_project_trajectory_cpp(
    fc_param, parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit, tol
  );
  double exposure = 0.0;
  for (int t = 0; t < n_time; ++t) {
    if (!R_finite(ntot[t]) || ntot[t] < 0.0) {
      Rcpp::stop("`ntot` must contain only finite non-negative values.");
    }
    exposure += ntot[t] * projected[t];
  }
  return exposure;
}

// [[Rcpp::export]]
Rcpp::NumericVector alfak_parent_opportunity_weights_cpp(Rcpp::NumericVector pij_values,
                                                         Rcpp::NumericVector parent_birth_times,
                                                         Rcpp::NumericVector timepoints,
                                                         Rcpp::NumericMatrix parent_xfit,
                                                         Rcpp::NumericVector ntot) {
  const int n_parents = pij_values.size();
  const int n_time = timepoints.size();
  if (parent_birth_times.size() != n_parents || parent_xfit.nrow() != n_parents ||
      parent_xfit.ncol() != n_time || ntot.size() != n_time) {
    Rcpp::stop("Malformed parent inputs for opportunity weights.");
  }
  Rcpp::NumericVector parent_weights(n_parents);
  if (n_parents == 0) {
    return parent_weights;
  }
  for (int p = 0; p < n_parents; ++p) {
    if (!R_finite(pij_values[p]) || pij_values[p] < 0.0 || !R_finite(parent_birth_times[p])) {
      Rcpp::stop("`pij_values` and `parent_birth_times` must be finite; pij values must be non-negative.");
    }
    double total = 0.0;
    for (int t = 0; t < n_time; ++t) {
      if (!R_finite(timepoints[t]) || !R_finite(ntot[t]) || ntot[t] < 0.0 ||
          !R_finite(parent_xfit(p, t))) {
        Rcpp::stop("`timepoints`, `ntot`, and `parent_xfit` must contain finite values; ntot must be non-negative.");
      }
      if (timepoints[t] >= parent_birth_times[p]) {
        total += ntot[t] * parent_xfit(p, t);
      }
    }
    parent_weights[p] = pij_values[p] * total;
  }
  double sum_parent_weights = finite_positive_sum(parent_weights);
  if (R_finite(sum_parent_weights) && sum_parent_weights > 0.0) {
    return parent_weights;
  }
  double sum_pij = finite_positive_sum(pij_values);
  if (R_finite(sum_pij) && sum_pij > 0.0) {
    return Rcpp::clone(pij_values);
  }
  return Rcpp::NumericVector(n_parents, 1.0);
}

// [[Rcpp::export]]
double alfak_weighted_parent_mean_cpp(Rcpp::NumericVector parent_fitness,
                                      Rcpp::NumericVector weights,
                                      double fallback_mean) {
  if (parent_fitness.size() == 0) {
    return fallback_mean;
  }
  double mean = weighted_mean_or_nan(parent_fitness, weights);
  if (R_finite(mean)) {
    return mean;
  }
  return fallback_mean;
}

// [[Rcpp::export]]
Rcpp::NumericVector alfak_neighbor_loglik_grid_cpp(Rcpp::NumericVector fc_grid,
                                                   Rcpp::NumericVector parent_fitness,
                                                   Rcpp::NumericVector pij_values,
                                                   Rcpp::NumericVector parent_birth_times,
                                                   Rcpp::NumericVector timepoints,
                                                   Rcpp::NumericMatrix parent_xfit,
                                                   Rcpp::NumericVector child_obs,
                                                   Rcpp::NumericVector ntot,
                                                   double tol) {
  const int n_grid = fc_grid.size();
  const int n_parents = parent_fitness.size();
  const int n_time = timepoints.size();
  validate_neighbor_projection_inputs(parent_fitness, pij_values, parent_birth_times, timepoints, parent_xfit, tol);
  if (child_obs.size() != n_time || ntot.size() != n_time) {
    Rcpp::stop("`child_obs`, `ntot`, and `timepoints` must have matching lengths.");
  }
  for (int t = 0; t < n_time; ++t) {
    if (!R_finite(child_obs[t]) || child_obs[t] < 0.0 || !is_integer_valued_scalar(child_obs[t]) ||
        !R_finite(ntot[t]) || ntot[t] < 0.0 || !is_integer_valued_scalar(ntot[t]) ||
        child_obs[t] > ntot[t]) {
      Rcpp::stop("`child_obs` and `ntot` must be valid non-negative integer counts with child_obs <= ntot.");
    }
  }
  Rcpp::NumericVector loglik(n_grid);
  for (int g = 0; g < n_grid; ++g) {
    double fc = fc_grid[g];
    if (!R_finite(fc)) {
      loglik[g] = R_NegInf;
      continue;
    }
    double total_ll = 0.0;
    for (int t = 0; t < n_time; ++t) {
      double xc_est = 0.0;
      for (int p = 0; p < n_parents; ++p) {
        double tt = std::max(0.0, timepoints[t] - parent_birth_times[p]);
        xc_est += fexp_stable_cpp(fc, parent_fitness[p], pij_values[p], tt, tol) * parent_xfit(p, t);
      }
      if (!R_finite(xc_est)) {
        total_ll += -1e9;
        continue;
      }
      xc_est = std::max(0.0, std::min(1.0, xc_est));
      double ll = R::dbinom(child_obs[t], ntot[t], xc_est, true);
      total_ll += R_finite(ll) ? ll : -1e9;
    }
    loglik[g] = total_ll;
  }
  return loglik;
}

// [[Rcpp::export]]
double alfak_nn_prior_marginal_negloglik_cpp(Rcpp::NumericMatrix loglik_mat,
                                             Rcpp::NumericVector fc_grid,
                                             Rcpp::NumericVector log_weights,
                                             Rcpp::NumericVector parent_means,
                                             Rcpp::NumericVector child_weights,
                                             double mu,
                                             double sigma) {
  const int n_children = loglik_mat.nrow();
  const int n_grid = loglik_mat.ncol();
  if (!R_finite(mu) || !R_finite(sigma) || sigma <= 0.0) {
    return 1e9;
  }
  if (fc_grid.size() != n_grid || log_weights.size() != n_grid ||
      parent_means.size() != n_children || child_weights.size() != n_children) {
    Rcpp::stop("Marginal likelihood inputs have incompatible dimensions.");
  }
  double total = 0.0;
  std::vector<double> vals(n_grid);
  for (int i = 0; i < n_children; ++i) {
    if (!R_finite(parent_means[i]) || !R_finite(child_weights[i]) || child_weights[i] < 0.0) {
      return 1e9;
    }
    if (child_weights[i] == 0.0) {
      continue;
    }
    bool any_finite = false;
    double max_val = R_NegInf;
    for (int g = 0; g < n_grid; ++g) {
      if (!R_finite(fc_grid[g]) || !R_finite(log_weights[g])) {
        return 1e9;
      }
      double log_prior = R::dnorm(fc_grid[g] - parent_means[i], mu, sigma, true);
      double val = loglik_mat(i, g) + log_prior + log_weights[g];
      vals[g] = val;
      if (R_finite(val)) {
        any_finite = true;
        if (val > max_val) {
          max_val = val;
        }
      }
    }
    if (!any_finite) {
      return 1e9;
    }
    double accum = 0.0;
    for (int g = 0; g < n_grid; ++g) {
      if (R_finite(vals[g])) {
        accum += std::exp(vals[g] - max_val);
      }
    }
    if (!(accum > 0.0) || !R_finite(accum)) {
      return 1e9;
    }
    total -= child_weights[i] * (max_val + std::log(accum));
  }
  if (!R_finite(total)) {
    return 1e9;
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::List alfak_two_shell_path_responsibilities_cpp(Rcpp::CharacterVector descendant,
                                                     Rcpp::NumericVector parent_anchor_exposure,
                                                     Rcpp::NumericVector transition_probability) {
  const int n = descendant.size();
  if (parent_anchor_exposure.size() != n || transition_probability.size() != n) {
    Rcpp::stop("Two-shell path responsibility inputs must be aligned.");
  }
  Rcpp::NumericVector path_supply(n);
  Rcpp::NumericVector path_responsibility(n);
  std::map<std::string, std::vector<int>> groups;
  for (int i = 0; i < n; ++i) {
    double supply = parent_anchor_exposure[i] * transition_probability[i];
    if (!R_finite(supply) || supply < 0.0) {
      supply = 0.0;
    }
    path_supply[i] = supply;
    groups[Rcpp::as<std::string>(descendant[i])].push_back(i);
  }
  for (const auto& entry : groups) {
    const std::vector<int>& idx = entry.second;
    double total = 0.0;
    for (int i : idx) {
      total += path_supply[i];
    }
    if (R_finite(total) && total > 0.0) {
      for (int i : idx) {
        path_responsibility[i] = path_supply[i] / total;
      }
    } else if (!idx.empty()) {
      double equal = 1.0 / static_cast<double>(idx.size());
      for (int i : idx) {
        path_responsibility[i] = equal;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("path_supply") = path_supply,
    Rcpp::Named("path_responsibility") = path_responsibility
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector alfak_group_cap_weights_cpp(Rcpp::CharacterVector group,
                                                Rcpp::NumericVector raw_weights,
                                                Rcpp::NumericVector cap_by_row) {
  const int n = group.size();
  if (raw_weights.size() != n || cap_by_row.size() != n) {
    Rcpp::stop("Group cap inputs must be aligned.");
  }
  Rcpp::NumericVector out = Rcpp::clone(raw_weights);
  std::map<std::string, std::vector<int>> groups;
  for (int i = 0; i < n; ++i) {
    if (!R_finite(out[i]) || out[i] < 0.0) {
      out[i] = 0.0;
    }
    groups[Rcpp::as<std::string>(group[i])].push_back(i);
  }
  for (const auto& entry : groups) {
    const std::vector<int>& idx = entry.second;
    double raw_sum = 0.0;
    double cap = R_NaReal;
    for (int i : idx) {
      raw_sum += out[i];
      if (!R_finite(cap) && R_finite(cap_by_row[i])) {
        cap = cap_by_row[i];
      }
    }
    if (!R_finite(cap) || cap < 0.0 || !R_finite(raw_sum) || raw_sum <= 0.0 || raw_sum <= cap) {
      continue;
    }
    double multiplier = cap / raw_sum;
    for (int i : idx) {
      out[i] *= multiplier;
    }
  }
  return out;
}

// [[Rcpp::export]]
double alfak_neighbor_two_shell_objective_cpp(double fc_param,
                                              Rcpp::NumericVector parent_fitness,
                                              Rcpp::NumericVector pij_values,
                                              Rcpp::NumericVector parent_birth_times,
                                              Rcpp::NumericVector timepoints,
                                              Rcpp::NumericMatrix parent_xfit,
                                              Rcpp::NumericVector child_obs,
                                              Rcpp::NumericVector ntot,
                                              double inward_prior_mean,
                                              double inward_prior_sd,
                                              Rcpp::NumericVector inward_prior_weights,
                                              bool do_inward_prior,
                                              Rcpp::NumericVector outward_fitness,
                                              double outward_prior_mean,
                                              Rcpp::NumericVector outward_prior_sd,
                                              Rcpp::NumericVector outward_prior_weights,
                                              double outward_lambda,
                                              double tol) {
  const int n_parents = parent_fitness.size();
  const int n_time = timepoints.size();
  const int n_outward = outward_fitness.size();
  if (n_parents == 0) {
    return 1e9;
  }
  if (!R_finite(fc_param) || !R_finite(tol) || tol <= 0.0) {
    Rcpp::stop("`fc_param` must be finite and `tol` must be a positive finite value.");
  }
  if (pij_values.size() != n_parents || parent_birth_times.size() != n_parents ||
      parent_xfit.nrow() != n_parents) {
    Rcpp::stop("Parent inputs must have matching lengths/rows.");
  }
  if (parent_xfit.ncol() != n_time) {
    Rcpp::stop("`parent_xfit` must have ncol equal to length(timepoints).");
  }
  if (child_obs.size() != n_time || ntot.size() != n_time) {
    Rcpp::stop("`child_obs`, `ntot`, and `timepoints` must have matching lengths.");
  }
  if (do_inward_prior) {
    if (inward_prior_weights.size() != n_parents) {
      Rcpp::stop("`inward_prior_weights` must have one entry per parent.");
    }
    if (!R_finite(inward_prior_mean) || !R_finite(inward_prior_sd) || inward_prior_sd <= 0.0) {
      Rcpp::stop("Inward prior parameters must be finite and `inward_prior_sd` must be positive.");
    }
  }
  if (outward_prior_sd.size() != n_outward || outward_prior_weights.size() != n_outward) {
    Rcpp::stop("Outward prior inputs must have matching lengths.");
  }
  if (!R_finite(outward_lambda) || outward_lambda < 0.0) {
    Rcpp::stop("`outward_lambda` must be a finite non-negative value.");
  }
  if (n_outward > 0 && !R_finite(outward_prior_mean)) {
    Rcpp::stop("`outward_prior_mean` must be finite when outward prior terms are supplied.");
  }

  double loglik = 0.0;
  for (int t = 0; t < n_time; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
    if (!R_finite(child_obs[t]) || child_obs[t] < 0.0 || !is_integer_valued_scalar(child_obs[t])) {
      Rcpp::stop("`child_obs` must contain only finite non-negative integer-valued counts.");
    }
    if (!R_finite(ntot[t]) || ntot[t] < 0.0 || !is_integer_valued_scalar(ntot[t])) {
      Rcpp::stop("`ntot` must contain only finite non-negative integer-valued counts.");
    }
    if (child_obs[t] > ntot[t]) {
      Rcpp::stop("`child_obs` must not exceed `ntot` at any timepoint.");
    }
    double xc_est = 0.0;
    for (int p = 0; p < n_parents; ++p) {
      if (!R_finite(parent_fitness[p]) || !R_finite(pij_values[p]) || pij_values[p] < 0.0 ||
          !R_finite(parent_birth_times[p]) || !R_finite(parent_xfit(p, t))) {
        Rcpp::stop("Parent fitness, transition probabilities, birth times, and parent_xfit must be finite; pij values must be non-negative.");
      }
      double tt = std::max(0.0, timepoints[t] - parent_birth_times[p]);
      xc_est += fexp_stable_cpp(fc_param, parent_fitness[p], pij_values[p], tt, tol) * parent_xfit(p, t);
    }

    if (!R_finite(xc_est)) {
      loglik += -1e9;
      continue;
    }

    xc_est = std::max(0.0, std::min(1.0, xc_est));
    double ll = R::dbinom(child_obs[t], ntot[t], xc_est, true);
    if (!R_finite(ll)) {
      ll = -1e9;
    }
    loglik += ll;
  }

  if (do_inward_prior) {
    for (int p = 0; p < n_parents; ++p) {
      if (!R_finite(inward_prior_weights[p]) || inward_prior_weights[p] < 0.0) {
        Rcpp::stop("`inward_prior_weights` must contain finite non-negative values.");
      }
      if (inward_prior_weights[p] == 0.0) {
        continue;
      }
      double prior_ll = R::dnorm(fc_param - parent_fitness[p], inward_prior_mean, inward_prior_sd, true);
      if (!R_finite(prior_ll)) {
        prior_ll = -1e9;
      }
      loglik += inward_prior_weights[p] * prior_ll;
    }
  }

  if (n_outward > 0 && outward_lambda > 0.0) {
    for (int k = 0; k < n_outward; ++k) {
      if (!R_finite(outward_fitness[k]) ||
          !R_finite(outward_prior_sd[k]) || outward_prior_sd[k] <= 0.0 ||
          !R_finite(outward_prior_weights[k]) || outward_prior_weights[k] < 0.0) {
        Rcpp::stop("Outward prior fitness, standard deviations, and weights must be finite; standard deviations must be positive and weights non-negative.");
      }
      if (outward_prior_weights[k] == 0.0) {
        continue;
      }
      double prior_ll = R::dnorm(outward_fitness[k] - fc_param, outward_prior_mean, outward_prior_sd[k], true);
      if (!R_finite(prior_ll)) {
        prior_ll = -1e9;
      }
      loglik += outward_lambda * outward_prior_weights[k] * prior_ll;
    }
  }

  return -loglik;
}

// [[Rcpp::export]]
Rcpp::List alfak_qr_accum_cpp(Rcpp::NumericMatrix x_trim,
                              Rcpp::NumericMatrix dx_dt) {
  const int K = x_trim.nrow();
  const int T = x_trim.ncol();
  if (K <= 0 || T < 0) {
    Rcpp::stop("`x_trim` must have positive dimensions.");
  }
  if (dx_dt.nrow() != K || dx_dt.ncol() != T) {
    Rcpp::stop("`x_trim` and `dx_dt` must have identical dimensions.");
  }
  Rcpp::NumericMatrix Q_accum(K, K);
  Rcpp::NumericVector r_accum(K);
  std::vector<double> xt(K);
  std::vector<double> xt_sq(K);
  std::vector<double> dx(K);

  for (int t = 0; t < T; ++t) {
    double sum_xt_sq = 0.0;
    double xt_dx_dot = 0.0;

    for (int i = 0; i < K; ++i) {
      if (!R_finite(x_trim(i, t)) || !R_finite(dx_dt(i, t))) {
        Rcpp::stop("`x_trim` and `dx_dt` must contain only finite values.");
      }
      xt[i] = x_trim(i, t);
      dx[i] = dx_dt(i, t);
      xt_sq[i] = xt[i] * xt[i];
      sum_xt_sq += xt_sq[i];
      xt_dx_dot += xt[i] * dx[i];
    }

    for (int i = 0; i < K; ++i) {
      r_accum[i] += xt[i] * dx[i] - xt[i] * xt_dx_dot;
      for (int j = i; j < K; ++j) {
        double value = (i == j ? xt_sq[i] : 0.0) -
          xt_sq[i] * xt[j] -
          xt[i] * xt_sq[j] +
          sum_xt_sq * xt[i] * xt[j];
        Q_accum(i, j) += value;
        if (j != i) {
          Q_accum(j, i) += value;
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("Q_accum") = Q_accum,
    Rcpp::Named("r_accum") = r_accum
  );
}
