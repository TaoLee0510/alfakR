// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
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

double weighted_parent_mean_cpp(const std::vector<double>& parent_fitness,
                                const std::vector<double>& parent_weights) {
  const std::size_t n = parent_fitness.size();
  if (n == 0 || parent_weights.size() != n) {
    return NA_REAL;
  }

  double weight_sum = 0.0;
  double weighted_sum = 0.0;
  bool can_weight = true;
  for (std::size_t i = 0; i < n; ++i) {
    if (!R_finite(parent_fitness[i])) {
      return NA_REAL;
    }
    if (!R_finite(parent_weights[i]) || parent_weights[i] < 0.0) {
      can_weight = false;
      break;
    }
    weight_sum += parent_weights[i];
    weighted_sum += parent_weights[i] * parent_fitness[i];
  }
  if (can_weight && weight_sum > 0.0 && R_finite(weight_sum)) {
    return weighted_sum / weight_sum;
  }

  double mean_sum = 0.0;
  int mean_n = 0;
  for (double value : parent_fitness) {
    if (R_finite(value)) {
      mean_sum += value;
      ++mean_n;
    }
  }
  if (mean_n == 0) {
    return NA_REAL;
  }
  return mean_sum / static_cast<double>(mean_n);
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
double alfak_joint_objective_cpp(Rcpp::NumericVector f_rel,
                                 Rcpp::NumericVector x0,
                                 Rcpp::NumericVector nn_fitness,
                                 Rcpp::NumericMatrix counts_fq,
                                 Rcpp::NumericVector timepoints,
                                 Rcpp::NumericVector viability_vec,
                                 double g0_val,
                                 bool correct_efflux,
                                 Rcpp::List nn_parent_indices,
                                 Rcpp::List nn_pij_values,
                                 Rcpp::NumericVector birth_times,
                                 Rcpp::NumericMatrix child_obs,
                                 Rcpp::NumericVector ntot,
                                 bool use_prior,
                                 double mu_delta,
                                 double sigma_delta,
                                 double weak_mu_sd,
                                 bool apply_sigma_regularization,
                                 double log_sigma_raw,
                                 double weak_log_sigma_sd,
                                 double tol) {
  const int K = counts_fq.nrow();
  const int T = counts_fq.ncol();
  const int n_children = nn_parent_indices.size();

  if (K <= 0) {
    Rcpp::stop("`counts_fq` must have at least one row.");
  }
  if (f_rel.size() != K || x0.size() != K || viability_vec.size() != K || birth_times.size() != K) {
    Rcpp::stop("Frequent-state inputs must all have length nrow(counts_fq).");
  }
  if (timepoints.size() != T || ntot.size() != T) {
    Rcpp::stop("`timepoints`, `ntot`, and ncol(counts_fq) must match.");
  }
  if (child_obs.nrow() != n_children || child_obs.ncol() != T) {
    Rcpp::stop("`child_obs` must have one row per child and one column per timepoint.");
  }
  if (nn_pij_values.size() != n_children || nn_fitness.size() != n_children) {
    Rcpp::stop("Nearest-neighbour parameter blocks must agree on the number of children.");
  }
  if (!R_finite(g0_val) || !R_finite(tol) || tol <= 0.0) {
    Rcpp::stop("`g0_val` must be finite and `tol` must be a positive finite value.");
  }
  if (use_prior) {
    if (!R_finite(mu_delta) || !R_finite(sigma_delta) || sigma_delta <= 0.0) {
      Rcpp::stop("When `use_prior` is TRUE, `mu_delta` must be finite and `sigma_delta` must be positive.");
    }
    if (!R_finite(weak_mu_sd) || weak_mu_sd <= 0.0) {
      Rcpp::stop("When `use_prior` is TRUE, `weak_mu_sd` must be a positive finite value.");
    }
  }
  if (apply_sigma_regularization && (!R_finite(log_sigma_raw) || !R_finite(weak_log_sigma_sd) || weak_log_sigma_sd <= 0.0)) {
    Rcpp::stop("Sigma regularization requires finite `log_sigma_raw` and positive finite `weak_log_sigma_sd`.");
  }

  double x0_sum = 0.0;
  std::vector<double> log_x0(K);
  std::vector<double> f_abs(K);
  std::vector<double> lv(K);
  for (int i = 0; i < K; ++i) {
    if (!R_finite(f_rel[i])) {
      Rcpp::stop("`f_rel` must contain only finite values.");
    }
    if (!R_finite(x0[i]) || x0[i] < 0.0) {
      Rcpp::stop("`x0` must contain only finite non-negative values.");
    }
    if (!R_finite(viability_vec[i]) || (correct_efflux && viability_vec[i] <= 0.0)) {
      Rcpp::stop("`viability_vec` must be finite and positive when `correct_efflux` is TRUE.");
    }
    if (!R_finite(birth_times[i])) {
      Rcpp::stop("`birth_times` must contain only finite values.");
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
    if (!R_finite(ntot[t]) || ntot[t] < 0.0 || !is_integer_valued_scalar(ntot[t])) {
      Rcpp::stop("`ntot` must contain only finite non-negative integer-valued counts.");
    }
  }
  for (int i = 0; i < K; ++i) {
    log_x0[i] = std::log(x0[i] / x0_sum);
  }

  double frequent_nll = 0.0;
  if (K > 1) {
    for (int t = 0; t < T; ++t) {
      for (int i = 0; i < K; ++i) {
        if (!R_finite(counts_fq(i, t)) || counts_fq(i, t) < 0.0) {
          Rcpp::stop("`counts_fq` must contain only finite non-negative values.");
        }
        lv[i] = log_x0[i] + f_rel[i] * timepoints[t];
      }
      const double denom = log_sum_exp_cpp(lv);
      for (int i = 0; i < K; ++i) {
        if (counts_fq(i, t) > 0.0) {
          frequent_nll -= counts_fq(i, t) * (lv[i] - denom);
        }
      }
    }
  }

  if (correct_efflux) {
    double sum_weighted_frel = 0.0;
    double sum_weights = 0.0;
    for (int i = 0; i < K; ++i) {
      sum_weighted_frel += (x0[i] * f_rel[i]) / viability_vec[i];
      sum_weights += x0[i] / viability_vec[i];
    }
    const double k_const = (sum_weighted_frel - g0_val) / sum_weights;
    for (int i = 0; i < K; ++i) {
      f_abs[i] = (f_rel[i] - k_const) / viability_vec[i];
    }
  } else {
    double x0_weighted_mean = 0.0;
    for (int i = 0; i < K; ++i) {
      x0_weighted_mean += x0[i] * f_rel[i];
    }
    const double shift = g0_val - x0_weighted_mean;
    for (int i = 0; i < K; ++i) {
      f_abs[i] = f_rel[i] + shift;
    }
  }

  Rcpp::NumericMatrix xfit(K, T);
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < K; ++i) {
      lv[i] = log_x0[i] + f_abs[i] * timepoints[t];
    }
    const double denom = log_sum_exp_cpp(lv);
    for (int i = 0; i < K; ++i) {
      xfit(i, t) = std::exp(lv[i] - denom);
    }
  }

  double total = frequent_nll;
  for (int child_idx = 0; child_idx < n_children; ++child_idx) {
    if (!R_finite(nn_fitness[child_idx])) {
      Rcpp::stop("`nn_fitness` must contain only finite values.");
    }

    Rcpp::IntegerVector parent_indices = nn_parent_indices[child_idx];
    Rcpp::NumericVector pij = nn_pij_values[child_idx];
    const int n_parents = parent_indices.size();
    if (n_parents == 0) {
      total += 1e9;
      continue;
    }
    if (pij.size() != n_parents) {
      Rcpp::stop("Each child must have matching parent index and pij lengths.");
    }

    std::vector<double> parent_fitness(n_parents);
    std::vector<double> parent_weights(n_parents);
    std::vector<int> zero_based_parent_indices(n_parents);
    for (int p = 0; p < n_parents; ++p) {
      const int parent_idx = parent_indices[p] - 1;
      if (parent_idx < 0 || parent_idx >= K) {
        Rcpp::stop("Each parent index must reference a frequent-state row.");
      }
      if (!R_finite(pij[p]) || pij[p] < 0.0) {
        Rcpp::stop("`nn_pij_values` must contain only finite non-negative values.");
      }
      parent_fitness[p] = f_abs[parent_idx];
      parent_weights[p] = pij[p];
      zero_based_parent_indices[p] = parent_idx;
    }

    const double parent_mean = weighted_parent_mean_cpp(parent_fitness, parent_weights);
    if (use_prior && !R_finite(parent_mean)) {
      total += 1e9;
      continue;
    }

    double child_nll = 0.0;
    for (int t = 0; t < T; ++t) {
      if (!R_finite(child_obs(child_idx, t)) ||
          child_obs(child_idx, t) < 0.0 ||
          !is_integer_valued_scalar(child_obs(child_idx, t))) {
        Rcpp::stop("`child_obs` must contain only finite non-negative integer-valued counts.");
      }
      if (child_obs(child_idx, t) > ntot[t]) {
        Rcpp::stop("Each child observation count must not exceed `ntot`.");
      }

      double xc_est = 0.0;
      for (int p = 0; p < n_parents; ++p) {
        const int parent_idx = zero_based_parent_indices[p];
        const double tt = std::max(0.0, timepoints[t] - birth_times[parent_idx]);
        xc_est += fexp_stable_cpp(
          nn_fitness[child_idx],
          f_abs[parent_idx],
          pij[p],
          tt,
          tol
        ) * xfit(parent_idx, t);
      }

      if (!R_finite(xc_est)) {
        child_nll += 1e9;
        continue;
      }
      xc_est = std::max(0.0, std::min(1.0, xc_est));
      double ll = R::dbinom(child_obs(child_idx, t), ntot[t], xc_est, true);
      if (!R_finite(ll)) {
        ll = -1e9;
      }
      child_nll -= ll;
    }

    if (use_prior) {
      double prior_ll = R::dnorm(nn_fitness[child_idx] - parent_mean, mu_delta, sigma_delta, true);
      if (!R_finite(prior_ll)) {
        prior_ll = -1e9;
      }
      child_nll -= prior_ll;
    }
    total += child_nll;
  }

  if (use_prior) {
    total += 0.5 * std::pow(mu_delta / weak_mu_sd, 2.0);
  }
  if (apply_sigma_regularization) {
    total += 0.5 * std::pow(log_sigma_raw / weak_log_sigma_sd, 2.0);
  }

  if (!R_finite(total)) {
    return 1e12;
  }
  return total;
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
