// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <sstream>
#include <cmath>
#include <numeric>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

namespace {

std::vector<int> parse_karyotype_string_cpp(const std::string& str) {
  if (str.empty()) {
    Rcpp::stop("Karyotype IDs must be non-empty character strings.");
  }
  std::vector<int> out;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, '.')) {
    if (token.empty()) {
      Rcpp::stop("Invalid karyotype ID '%s': empty chromosome-count component.", str);
    }
    std::size_t pos = 0;
    int value = 0;
    try {
      value = std::stoi(token, &pos);
    } catch (...) {
      Rcpp::stop("Invalid karyotype ID '%s': non-integer component '%s'.", str, token);
    }
    if (pos != token.size()) {
      Rcpp::stop("Invalid karyotype ID '%s': malformed component '%s'.", str, token);
    }
    if (value < 0) {
      Rcpp::stop("Invalid karyotype ID '%s': components must be non-negative integers.", str);
    }
    out.push_back(value);
  }
  if (out.empty()) {
    Rcpp::stop("Karyotype IDs must be non-empty character strings.");
  }
  return out;
}

double pij_impl(int i, int j, double beta) {
  double qij = 0.0;
  if (std::abs(i - j) > i) {
    return qij;
  }
  if (j == 0) {
    j = 2 * i;
  }
  for (int z = std::abs(i - j); z <= i; z += 2) {
    qij += R::choose(i, z) * std::pow(beta, z) * std::pow(1 - beta, i - z) *
      std::pow(0.5, z) * R::choose(z, (z + i - j) / 2);
  }
  return qij;
}

std::string karyotype_vector_to_string(const std::vector<int>& x) {
  std::ostringstream ss;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (i > 0) ss << ".";
    ss << x[i];
  }
  return ss.str();
}

std::vector<std::vector<int>> parse_karyotype_ids_cpp(Rcpp::CharacterVector ids) {
  const int n = ids.size();
  if (n == 0) {
    Rcpp::stop("Karyotype IDs must be non-empty character strings.");
  }
  std::vector<std::vector<int>> parsed(n);
  std::unordered_set<std::string> seen;
  int k = -1;
  for (int i = 0; i < n; ++i) {
    std::string id = Rcpp::as<std::string>(ids[i]);
    if (!seen.insert(id).second) {
      Rcpp::stop("Karyotype IDs must be unique.");
    }
    parsed[i] = parse_karyotype_string_cpp(id);
    if (k < 0) {
      k = static_cast<int>(parsed[i].size());
    } else if (static_cast<int>(parsed[i].size()) != k) {
      Rcpp::stop("All karyotype IDs must have the same number of dot-separated components.");
    }
  }
  return parsed;
}

std::vector<std::vector<int>> generate_one_step_neighbors_cpp(const std::vector<std::vector<int>>& ids,
                                                              bool remove_nullisomes,
                                                              const std::unordered_set<std::string>& originals) {
  std::vector<std::vector<int>> out;
  std::unordered_set<std::string> seen;
  if (ids.empty()) {
    return out;
  }
  const int k = ids[0].size();
  for (const auto& base : ids) {
    for (int chr = 0; chr < k; ++chr) {
      for (int delta : {-1, 1}) {
        std::vector<int> candidate = base;
        candidate[chr] += delta;
        if (remove_nullisomes && candidate[chr] < 1) {
          continue;
        }
        std::string key = karyotype_vector_to_string(candidate);
        if (originals.find(key) != originals.end()) {
          continue;
        }
        if (seen.insert(key).second) {
          out.push_back(candidate);
        }
      }
    }
  }
  return out;
}

double transition_probability_vec_cpp(const std::vector<int>& parent,
                                      const std::vector<int>& child,
                                      double beta) {
  if (parent.size() != child.size()) {
    Rcpp::stop("Parent and child karyotypes must have matching dimensions.");
  }
  double q = 1.0;
  for (std::size_t i = 0; i < parent.size(); ++i) {
    q *= pij_impl(parent[i], child[i], beta);
  }
  return q;
}

NumericMatrix validate_transition_matrix(List parms, int expected_size, const char* state_name) {
  if (!parms.containsElementNamed("A")) {
    Rcpp::stop("`parms$A` must be provided.");
  }
  SEXP a_sexp = parms["A"];
  if (!Rf_isMatrix(a_sexp) || (TYPEOF(a_sexp) != REALSXP && TYPEOF(a_sexp) != INTSXP)) {
    Rcpp::stop("`parms$A` must be a numeric matrix.");
  }
  NumericMatrix A(a_sexp);
  if (A.nrow() != A.ncol()) {
    Rcpp::stop("`parms$A` must be a square matrix.");
  }
  if (A.nrow() != expected_size) {
    Rcpp::stop("`parms$A` must have nrow(A) == ncol(A) == length(%s).", state_name);
  }
  return A;
}

} // namespace

// [[Rcpp::export]]
double pij_cpp(int i, int j, double beta) {
  if (i < 0) {
    Rcpp::stop("`i` must be a non-negative integer.");
  }
  if (j < 0) {
    Rcpp::stop("`j` must be a non-negative integer.");
  }
  if (!std::isfinite(beta) || beta < 0.0 || beta > 1.0) {
    Rcpp::stop("`beta` must be finite and in [0, 1].");
  }
  double qij = pij_impl(i, j, beta);
  if (!std::isfinite(qij) || qij < 0.0 || qij > 1.0) {
    Rcpp::stop("Internal error: computed `pij` is not finite or not in [0, 1].");
  }
  return qij;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix gen_all_neighbours_cpp(Rcpp::CharacterVector ids,
                                           bool remove_nullisomes = true) {
  std::vector<std::vector<int>> parsed = parse_karyotype_ids_cpp(ids);
  std::unordered_set<std::string> originals;
  for (const auto& x : parsed) {
    originals.insert(karyotype_vector_to_string(x));
  }
  std::vector<std::vector<int>> neighbors = generate_one_step_neighbors_cpp(parsed, remove_nullisomes, originals);
  const int n = neighbors.size();
  const int k = parsed[0].size();
  Rcpp::NumericMatrix out(n, k);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      out(i, j) = neighbors[i][j];
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List gen_nn_info_cpp(Rcpp::CharacterVector fq, double beta) {
  if (!std::isfinite(beta) || beta < 0.0 || beta > 1.0) {
    Rcpp::stop("`beta` must be finite and in [0, 1].");
  }
  std::vector<std::vector<int>> fq_parsed = parse_karyotype_ids_cpp(fq);
  std::unordered_set<std::string> frequent_set;
  std::unordered_map<std::string, std::vector<int>> frequent_vectors;
  for (int i = 0; i < fq.size(); ++i) {
    std::string id = Rcpp::as<std::string>(fq[i]);
    frequent_set.insert(id);
    frequent_vectors[id] = fq_parsed[i];
  }
  std::vector<std::vector<int>> candidates = generate_one_step_neighbors_cpp(fq_parsed, true, frequent_set);
  Rcpp::List out(candidates.size());
  for (std::size_t idx = 0; idx < candidates.size(); ++idx) {
    const std::vector<int>& child = candidates[idx];
    std::string child_id = karyotype_vector_to_string(child);
    std::vector<std::vector<int>> child_wrap(1, child);
    std::unordered_set<std::string> child_original;
    child_original.insert(child_id);
    std::vector<std::vector<int>> parent_candidates = generate_one_step_neighbors_cpp(child_wrap, true, child_original);
    std::vector<std::string> parent_ids;
    std::vector<double> pij_values;
    for (const auto& parent : parent_candidates) {
      std::string parent_id = karyotype_vector_to_string(parent);
      if (frequent_set.find(parent_id) == frequent_set.end()) {
        continue;
      }
      double q = transition_probability_vec_cpp(parent, child, beta);
      parent_ids.push_back(parent_id);
      pij_values.push_back(q);
    }
    Rcpp::CharacterVector nj(parent_ids.begin(), parent_ids.end());
    Rcpp::NumericVector pij(pij_values.begin(), pij_values.end());
    out[idx] = Rcpp::List::create(
      Rcpp::Named("ni") = child_id,
      Rcpp::Named("nj") = nj,
      Rcpp::Named("pij") = pij
    );
  }
  return out;
}

 //' Prepare triplet inputs (i, j, x, dims, dimnames) for sparse A matrix.
 //' @param k_str Character vector of karyotype strings, e.g. "2.2.3".
 //' @param beta Double mis-segregation probability per chromosome.
 //' @param Nmax Optional max total mis-segregations allowed per division.
 //' If not provided, no cap is applied.
 //' @return List with elements i (rows), j (cols), x (values), dims, dimnames.
 //' @export
 // [[Rcpp::export]]
 List get_A_inputs(CharacterVector k_str, double beta, Nullable<double> Nmax_ = R_NilValue) {
   double Nmax = R_PosInf;
   if (Nmax_.isNotNull()) Nmax = as<double>(Nmax_);
   if (!std::isfinite(beta) || beta < 0.0 || beta > 1.0) {
     Rcpp::stop("`beta` must be finite and in [0, 1].");
   }
   if (!(std::isinf(Nmax) || (std::isfinite(Nmax) && Nmax >= 0.0))) {
     Rcpp::stop("`Nmax` must be Inf or a non-negative finite number.");
   }
   int n = k_str.size();
   std::vector<std::vector<int>> k_list(n);
   int num_chrom_types = -1;
   for (int i = 0; i < n; ++i) {
     std::string k_id = Rcpp::as<std::string>(k_str[i]);
     k_list[i] = parse_karyotype_string_cpp(k_id);
     if (num_chrom_types < 0) {
       num_chrom_types = static_cast<int>(k_list[i].size());
     } else if (static_cast<int>(k_list[i].size()) != num_chrom_types) {
       Rcpp::stop("All karyotype IDs must have the same number of dot-separated components.");
     }
   }
   // triplet containers (1-based indices)
   std::vector<int> ii, jj;
   std::vector<double> xx;
   std::size_t cap = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
   ii.reserve(cap);
   jj.reserve(cap);
   xx.reserve(cap);
   
   for (int i = 0; i < n; ++i) {
     const auto& ki = k_list[i];
     for (int j = 0; j < n; ++j) {
       const auto& kj = k_list[j];
       double tot = 0;
       for (size_t m = 0; m < ki.size(); ++m) tot += std::abs(ki[m] - kj[m]);
       if (tot > Nmax) continue;
       double qij = 1.0;
       for (size_t m = 0; m < ki.size(); ++m) qij *= pij_impl(ki[m], kj[m], beta);
       double val = (i == j ? (2 * qij - 1) : (2 * qij));
       if (val != 0.0) {
         ii.push_back(i + 1);
         jj.push_back(j + 1);
         xx.push_back(val);
       }
     }
   }
   // dims and dimnames
   IntegerVector dims = IntegerVector::create(n, n);
   List dimnames = List::create(k_str, k_str);
   return List::create(
     _["i"] = ii,
     _["j"] = jj,
     _["x"] = xx,
     _["dims"] = dims,
     _["dimnames"] = dimnames
   );
 }
 
// [[Rcpp::export]]
List chrmod_cpp(double time, NumericVector state, List parms) {
   NumericMatrix A = validate_transition_matrix(parms, state.size(), "state");
   int n = state.size();
   NumericVector ds(n);
   for (int j = 0; j < n; ++j) {
     double acc = 0;
     for (int i = 0; i < n; ++i) acc += state[i] * A(i, j);
     ds[j] = acc;
   }
   return List::create(ds);
 }
 
// [[Rcpp::export]]
List chrmod_rel_cpp(double time, NumericVector x, List parms) {
   NumericMatrix A = validate_transition_matrix(parms, x.size(), "x");
   int n = x.size();
   NumericVector g(n);
   for (int j = 0; j < n; ++j) {
     double acc = 0;
     for (int i = 0; i < n; ++i) acc += x[i] * A(i, j);
     g[j] = acc;
   }
   double phi = std::accumulate(g.begin(), g.end(), 0.0);
   NumericVector dx(n);
   for (int k = 0; k < n; ++k) dx[k] = g[k] - x[k] * phi;
   return List::create(dx);
 }
 
