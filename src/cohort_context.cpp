// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace Rcpp;

namespace {

inline bool finite_num(double x) {
  return R_finite(x);
}

inline double weight_at(const NumericVector& w, int i) {
  if (w.size() == 0) return 1.0;
  if (i < w.size() && finite_num(w[i]) && w[i] >= 0.0) return w[i];
  return 0.0;
}

double profile_distance_cpp(const NumericVector& target,
                            const NumericMatrix& evidence,
                            int row,
                            const NumericVector& chromosome_weights,
                            int method) {
  const int k = target.size();
  if (evidence.ncol() != k) return R_PosInf;

  if (method == 1 || method == 2) {
    // Hellinger or Jensen-Shannon on non-negative normalized mass.
    std::vector<double> pa(k), pb(k);
    double suma = 0.0, sumb = 0.0;
    for (int j = 0; j < k; ++j) {
      double w = weight_at(chromosome_weights, j);
      double a = target[j];
      double b = evidence(row, j);
      if (!finite_num(a) || !finite_num(b) || w < 0.0) return R_PosInf;
      a = std::max(0.0, a * w);
      b = std::max(0.0, b * w);
      pa[j] = a;
      pb[j] = b;
      suma += a;
      sumb += b;
    }
    if (!finite_num(suma) || suma <= 0.0) {
      for (int j = 0; j < k; ++j) pa[j] = 1.0 / static_cast<double>(k);
    } else {
      for (int j = 0; j < k; ++j) pa[j] /= suma;
    }
    if (!finite_num(sumb) || sumb <= 0.0) {
      for (int j = 0; j < k; ++j) pb[j] = 1.0 / static_cast<double>(k);
    } else {
      for (int j = 0; j < k; ++j) pb[j] /= sumb;
    }
    if (method == 1) {
      double ss = 0.0;
      for (int j = 0; j < k; ++j) {
        double d = std::sqrt(pa[j]) - std::sqrt(pb[j]);
        ss += d * d;
      }
      return std::sqrt(ss / 2.0);
    }
    const double eps = 1e-12;
    double js = 0.0;
    for (int j = 0; j < k; ++j) {
      double a = std::max(pa[j], eps);
      double b = std::max(pb[j], eps);
      double m = 0.5 * (a + b);
      js += 0.5 * a * std::log(a / m) + 0.5 * b * std::log(b / m);
    }
    return std::sqrt(std::max(0.0, js));
  }

  if (method == 3) {
    // Cosine distance.
    double dot = 0.0, na = 0.0, nb = 0.0;
    for (int j = 0; j < k; ++j) {
      double w = weight_at(chromosome_weights, j);
      double a = target[j];
      double b = evidence(row, j);
      if (!finite_num(a) || !finite_num(b) || w < 0.0) return R_PosInf;
      a *= w;
      b *= w;
      dot += a * b;
      na += a * a;
      nb += b * b;
    }
    double denom = std::sqrt(na) * std::sqrt(nb);
    if (!finite_num(denom) || denom <= 0.0) return 1.0;
    return std::max(0.0, 1.0 - dot / denom);
  }

  if (method == 4) {
    double ss = 0.0;
    for (int j = 0; j < k; ++j) {
      double w = weight_at(chromosome_weights, j);
      double a = target[j];
      double b = evidence(row, j);
      if (!finite_num(a) || !finite_num(b) || w < 0.0) return R_PosInf;
      double d = a - b;
      ss += w * d * d;
    }
    return std::sqrt(std::max(0.0, ss));
  }

  // Manhattan.
  double total = 0.0;
  for (int j = 0; j < k; ++j) {
    double w = weight_at(chromosome_weights, j);
    double a = target[j];
    double b = evidence(row, j);
    if (!finite_num(a) || !finite_num(b) || w < 0.0) return R_PosInf;
    total += w * std::abs(a - b);
  }
  return total;
}

double finite_or_zero(double x) {
  return finite_num(x) ? x : 0.0;
}

double event_distance_cpp(int target_chr,
                          int evidence_chr,
                          int target_direction,
                          int evidence_direction,
                          double target_size,
                          double evidence_size,
                          double target_local_copy,
                          double evidence_local_copy,
                          double target_local_z,
                          double evidence_local_z,
                          double target_delta_area,
                          double evidence_delta_area,
                          double target_delta_burden,
                          double evidence_delta_burden,
                          int event_match) {
  bool same_chr = target_chr != NA_INTEGER && evidence_chr != NA_INTEGER && target_chr == evidence_chr;
  bool same_direction = target_direction == evidence_direction;
  if (event_match == 1 && (!same_chr || !same_direction)) return R_PosInf;
  if (event_match == 2 && !same_direction) return R_PosInf;

  double chr_penalty = same_chr ? 0.0 : 1.0;
  double direction_penalty = same_direction ? 0.0 : 2.0;
  double vals[7] = {
    chr_penalty,
    direction_penalty,
    finite_or_zero(std::abs(target_size - evidence_size)),
    finite_or_zero(std::abs(target_local_copy - evidence_local_copy)),
    finite_or_zero(std::abs(target_local_z - evidence_local_z)),
    finite_or_zero(std::abs(target_delta_area - evidence_delta_area)),
    finite_or_zero(std::abs(target_delta_burden - evidence_delta_burden))
  };
  double ss = 0.0;
  for (double v : vals) ss += v * v;
  return std::sqrt(ss);
}

double scaled_component(double d, double bw) {
  if (!finite_num(d)) return R_PosInf;
  if (!finite_num(bw) || bw <= 0.0) bw = 1.0;
  return d / bw;
}

struct KernelRow {
  int index;
  double context_distance;
  double profile_distance;
  double area_distance;
  double burden_distance;
  double local_distance;
  double event_distance;
  double kernel_weight;
  double quality_weight;
  double final_weight;
};

double median_positive(std::vector<double>& values, double fallback) {
  values.erase(
    std::remove_if(values.begin(), values.end(), [](double x) {
      return !finite_num(x) || x <= 0.0;
    }),
    values.end()
  );
  const int n = values.size();
  if (n == 0) return fallback;
  std::sort(values.begin(), values.end());
  if (n % 2 == 1) {
    return values[n / 2];
  }
  return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

std::string collapse_unique_strings(const std::set<std::string>& values) {
  std::ostringstream ss;
  bool first = true;
  for (const auto& value : values) {
    if (!first) ss << ";";
    ss << value;
    first = false;
  }
  return ss.str();
}

} // namespace

// [[Rcpp::export]]
Rcpp::DataFrame context_kernel_weights_cpp(
    Rcpp::NumericVector target_profile,
    double target_total_cn,
    double target_burden,
    double target_local_copy,
    double target_local_z,
    int target_transition_chr,
    int target_direction_code,
    double target_transition_size,
    double target_delta_total_cn,
    double target_delta_burden,
    Rcpp::NumericMatrix evidence_profile_matrix,
    Rcpp::NumericVector evidence_total_cn,
    Rcpp::NumericVector evidence_burden,
    Rcpp::NumericVector evidence_local_copy,
    Rcpp::NumericVector evidence_local_z,
    Rcpp::IntegerVector evidence_transition_chr,
    Rcpp::IntegerVector evidence_direction_code,
    Rcpp::NumericVector evidence_transition_size,
    Rcpp::NumericVector evidence_delta_total_cn,
    Rcpp::NumericVector evidence_delta_burden,
    Rcpp::NumericVector quality_weight,
    Rcpp::NumericVector bandwidths,
    Rcpp::NumericVector component_weights,
    Rcpp::NumericVector chromosome_weights,
    int event_match_code,
    int profile_distance_code,
    int k_nearest,
    double min_kernel_weight) {
  const int n = evidence_profile_matrix.nrow();
  const int k = evidence_profile_matrix.ncol();
  if (target_profile.size() != k) {
    Rcpp::stop("`target_profile` length must equal ncol(evidence_profile_matrix).");
  }
  if (evidence_total_cn.size() != n || evidence_burden.size() != n ||
      evidence_local_copy.size() != n || evidence_local_z.size() != n ||
      evidence_transition_chr.size() != n || evidence_direction_code.size() != n ||
      evidence_transition_size.size() != n || evidence_delta_total_cn.size() != n ||
      evidence_delta_burden.size() != n || quality_weight.size() != n) {
    Rcpp::stop("Evidence vectors must be aligned with evidence_profile_matrix rows.");
  }
  if (bandwidths.size() < 5 || component_weights.size() < 5) {
    Rcpp::stop("`bandwidths` and `component_weights` must have length at least 5.");
  }
  if (k_nearest < 1) {
    Rcpp::stop("`k_nearest` must be positive.");
  }
  if (!finite_num(min_kernel_weight) || min_kernel_weight < 0.0) {
    Rcpp::stop("`min_kernel_weight` must be non-negative and finite.");
  }

  std::vector<KernelRow> rows;
  rows.reserve(std::min(n, k_nearest));

  for (int i = 0; i < n; ++i) {
    double q = quality_weight[i];
    if (!finite_num(q) || q <= 0.0) continue;

    double d_profile = profile_distance_cpp(target_profile, evidence_profile_matrix, i, chromosome_weights, profile_distance_code);
    double d_area = std::abs(target_total_cn - evidence_total_cn[i]);
    double d_burden = std::abs(target_burden - evidence_burden[i]);
    double dlc = target_local_copy - evidence_local_copy[i];
    double dlz = target_local_z - evidence_local_z[i];
    if (!finite_num(dlc)) dlc = 0.0;
    if (!finite_num(dlz)) dlz = 0.0;
    double d_local = std::sqrt(dlc * dlc + dlz * dlz);
    double d_event = event_distance_cpp(
      target_transition_chr,
      evidence_transition_chr[i],
      target_direction_code,
      evidence_direction_code[i],
      target_transition_size,
      evidence_transition_size[i],
      target_local_copy,
      evidence_local_copy[i],
      target_local_z,
      evidence_local_z[i],
      target_delta_total_cn,
      evidence_delta_total_cn[i],
      target_delta_burden,
      evidence_delta_burden[i],
      event_match_code
    );

    double comps[5] = {
      scaled_component(d_profile, bandwidths[0]),
      scaled_component(d_area, bandwidths[1]),
      scaled_component(d_burden, bandwidths[2]),
      scaled_component(d_local, bandwidths[3]),
      scaled_component(d_event, bandwidths[4])
    };
    double ss = 0.0;
    bool finite = true;
    for (int j = 0; j < 5; ++j) {
      double cw = component_weights[j];
      if (!finite_num(cw) || cw < 0.0) cw = 0.0;
      if (!finite_num(comps[j])) {
        finite = false;
        break;
      }
      ss += cw * comps[j] * comps[j];
    }
    if (!finite) continue;
    double d_total = std::sqrt(std::max(0.0, ss));
    double kw = std::exp(-0.5 * d_total * d_total);
    double fw = kw * q;
    if (!finite_num(fw) || fw < min_kernel_weight) continue;
    rows.push_back(KernelRow{
      i + 1,
      d_total,
      d_profile,
      d_area,
      d_burden,
      d_local,
      d_event,
      kw,
      q,
      fw
    });
  }

  std::sort(rows.begin(), rows.end(), [](const KernelRow& a, const KernelRow& b) {
    if (a.final_weight == b.final_weight) return a.index < b.index;
    return a.final_weight > b.final_weight;
  });
  if (static_cast<int>(rows.size()) > k_nearest) {
    rows.resize(k_nearest);
  }

  const int m = rows.size();
  IntegerVector evidence_index(m);
  NumericVector context_distance(m), profile_distance(m), area_distance(m), burden_distance(m);
  NumericVector local_distance(m), event_distance(m), kernel_weight(m), quality_weight_out(m), final_weight(m);
  for (int i = 0; i < m; ++i) {
    evidence_index[i] = rows[i].index;
    context_distance[i] = rows[i].context_distance;
    profile_distance[i] = rows[i].profile_distance;
    area_distance[i] = rows[i].area_distance;
    burden_distance[i] = rows[i].burden_distance;
    local_distance[i] = rows[i].local_distance;
    event_distance[i] = rows[i].event_distance;
    kernel_weight[i] = rows[i].kernel_weight;
    quality_weight_out[i] = rows[i].quality_weight;
    final_weight[i] = rows[i].final_weight;
  }

  return DataFrame::create(
    Named("evidence_index") = evidence_index,
    Named("context_distance") = context_distance,
    Named("profile_distance") = profile_distance,
    Named("area_distance") = area_distance,
    Named("burden_distance") = burden_distance,
    Named("local_distance") = local_distance,
    Named("event_distance") = event_distance,
    Named("kernel_weight") = kernel_weight,
    Named("quality_weight") = quality_weight_out,
    Named("final_weight") = final_weight
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector context_bandwidths_cpp(
    Rcpp::NumericMatrix evidence_profile_matrix,
    Rcpp::NumericVector evidence_total_cn,
    Rcpp::NumericVector evidence_burden,
    Rcpp::NumericVector evidence_local_copy,
    Rcpp::NumericVector evidence_local_z,
    Rcpp::IntegerVector evidence_transition_chr,
    Rcpp::IntegerVector evidence_direction_code,
    Rcpp::NumericVector evidence_transition_size,
    Rcpp::NumericVector evidence_delta_total_cn,
    Rcpp::NumericVector evidence_delta_burden,
    Rcpp::NumericVector chromosome_weights,
    int profile_distance_code) {
  const int n = evidence_profile_matrix.nrow();
  if (evidence_total_cn.size() != n || evidence_burden.size() != n ||
      evidence_local_copy.size() != n || evidence_local_z.size() != n ||
      evidence_transition_chr.size() != n || evidence_direction_code.size() != n ||
      evidence_transition_size.size() != n || evidence_delta_total_cn.size() != n ||
      evidence_delta_burden.size() != n) {
    Rcpp::stop("Evidence vectors must be aligned with evidence_profile_matrix rows.");
  }
  std::vector<double> profile_vals;
  std::vector<double> area_vals;
  std::vector<double> burden_vals;
  std::vector<double> local_vals;
  std::vector<double> event_vals;
  if (n > 1) {
    std::size_t cap = static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2;
    profile_vals.reserve(cap);
    area_vals.reserve(cap);
    burden_vals.reserve(cap);
    local_vals.reserve(cap);
    event_vals.reserve(cap);
  }
  for (int i = 0; i < n - 1; ++i) {
    Rcpp::NumericVector target_profile = evidence_profile_matrix(i, Rcpp::_);
    for (int j = i + 1; j < n; ++j) {
      profile_vals.push_back(profile_distance_cpp(target_profile, evidence_profile_matrix, j, chromosome_weights, profile_distance_code));
      area_vals.push_back(std::abs(evidence_total_cn[i] - evidence_total_cn[j]));
      burden_vals.push_back(std::abs(evidence_burden[i] - evidence_burden[j]));
      double dlc = evidence_local_copy[i] - evidence_local_copy[j];
      double dlz = evidence_local_z[i] - evidence_local_z[j];
      if (!finite_num(dlc)) dlc = 0.0;
      if (!finite_num(dlz)) dlz = 0.0;
      local_vals.push_back(std::sqrt(dlc * dlc + dlz * dlz));
      event_vals.push_back(event_distance_cpp(
        evidence_transition_chr[i],
        evidence_transition_chr[j],
        evidence_direction_code[i],
        evidence_direction_code[j],
        evidence_transition_size[i],
        evidence_transition_size[j],
        evidence_local_copy[i],
        evidence_local_copy[j],
        evidence_local_z[i],
        evidence_local_z[j],
        evidence_delta_total_cn[i],
        evidence_delta_total_cn[j],
        evidence_delta_burden[i],
        evidence_delta_burden[j],
        3
      ));
    }
  }
  Rcpp::NumericVector out = Rcpp::NumericVector::create(
    Rcpp::Named("profile") = median_positive(profile_vals, 0.25),
    Rcpp::Named("area") = median_positive(area_vals, 2.0),
    Rcpp::Named("burden") = median_positive(burden_vals, 2.0),
    Rcpp::Named("local") = median_positive(local_vals, 1.0),
    Rcpp::Named("event") = median_positive(event_vals, 1.0)
  );
  return out;
}

// [[Rcpp::export]]
Rcpp::DataFrame context_patient_level_neighbors_cpp(Rcpp::IntegerVector evidence_index,
                                                    Rcpp::CharacterVector patient_id,
                                                    Rcpp::CharacterVector child_karyotype,
                                                    Rcpp::NumericVector delta_hat,
                                                    Rcpp::NumericVector delta_se,
                                                    Rcpp::NumericVector final_weight,
                                                    double sd_floor) {
  const int n = evidence_index.size();
  if (patient_id.size() != n || child_karyotype.size() != n || delta_hat.size() != n ||
      delta_se.size() != n || final_weight.size() != n) {
    Rcpp::stop("Context patient-neighbor aggregation inputs must be aligned.");
  }
  if (!finite_num(sd_floor) || sd_floor <= 0.0) {
    Rcpp::stop("`sd_floor` must be a positive finite value.");
  }
  struct Accum {
    double patient_weight = 0.0;
    double inv_sum = 0.0;
    double inv_delta_sum = 0.0;
    int n = 0;
    std::set<std::string> children;
  };
  std::map<std::string, Accum> groups;
  for (int i = 0; i < n; ++i) {
    double ww = final_weight[i];
    if (!finite_num(ww) || ww < 0.0) ww = 0.0;
    std::string pid = Rcpp::as<std::string>(patient_id[i]);
    Accum& acc = groups[pid];
    acc.patient_weight += ww;
    acc.n += 1;
    acc.children.insert(Rcpp::as<std::string>(child_karyotype[i]));
    double se = delta_se[i];
    if (!finite_num(se) || se < sd_floor) se = sd_floor;
    double inv = ww / (se * se + sd_floor * sd_floor);
    double delta = delta_hat[i];
    if (finite_num(inv) && inv > 0.0 && finite_num(delta)) {
      acc.inv_sum += inv;
      acc.inv_delta_sum += inv * delta;
    }
  }
  const int m = groups.size();
  Rcpp::CharacterVector out_patient(m);
  Rcpp::NumericVector out_mean(m), out_se(m), out_weight(m);
  Rcpp::IntegerVector out_n(m);
  Rcpp::CharacterVector out_children(m);
  int row = 0;
  for (const auto& entry : groups) {
    const Accum& acc = entry.second;
    out_patient[row] = entry.first;
    out_mean[row] = acc.inv_sum > 0.0 ? acc.inv_delta_sum / acc.inv_sum : R_NaReal;
    out_se[row] = acc.inv_sum > 0.0 ? std::sqrt(1.0 / acc.inv_sum) : R_NaReal;
    out_weight[row] = acc.patient_weight;
    out_n[row] = acc.n;
    out_children[row] = collapse_unique_strings(acc.children);
    ++row;
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("patient_id") = out_patient,
    Rcpp::Named("delta_patient_mean") = out_mean,
    Rcpp::Named("delta_patient_se") = out_se,
    Rcpp::Named("patient_weight") = out_weight,
    Rcpp::Named("n_context_neighbors") = out_n,
    Rcpp::Named("child_karyotype") = out_children
  );
}
