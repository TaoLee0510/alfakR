// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <sstream>
#include <cmath>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

//' Compute the probability pij that a parent with i chromosomes produces
//' a daughter with j chromosomes under mis-segregation rate beta.
//' @param i Integer number of chromosomes in the parent cell.
//' @param j Integer number of chromosomes in the daughter cell.
//' @param beta Double mis-segregation probability per chromosome copy.
//' @return Double transition probability pij.
//' @export
// [[Rcpp::export]]
double pij(int i, int j, double beta) {
 double qij = 0.0;
 if (std::abs(i - j) > i) {
   return qij;
 }
 if (j == 0) {
     j = 2 * i;
   }
   for (int z = std::abs(i - j); z <= i; z += 2) {
     qij += R::choose(i, z) * std::pow(beta, z) * std::pow(1 - beta, i - z)
     * std::pow(0.5, z) * R::choose(z, (z + i - j) / 2);
   }
   return qij;
 }
 
 //' Split a karyotype string "a.b.c" into an integer vector c(a, b, c).
 //' @param s String or character scalar: dot-separated chromosome counts.
 //' @return Integer vector of split components.
 //' @export
 // [[Rcpp::export]]
 IntegerVector s2v(SEXP s) {
   std::string str = Rcpp::as<std::string>(s);
   std::vector<int> out;
   std::stringstream ss(str);
   std::string token;
   while (std::getline(ss, token, '.')) {
     out.push_back(std::stoi(token));
   }
   return wrap(out);
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
   int n = k_str.size();
   std::vector<std::vector<int>> k_list(n);
   for (int i = 0; i < n; ++i) {
     k_list[i] = as<std::vector<int>>(s2v(k_str[i]));
   }
   // triplet containers (1-based indices)
   std::vector<int> ii, jj;
   std::vector<double> xx;
   // reserve enough space: 1 + n * (number of chromosome types) * 2 * Nmax
   int num_chrom_types = k_list.empty() ? 0 : k_list[0].size();
   int cap = 1 + n * num_chrom_types * 2 * static_cast<int>(Nmax);
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
       for (size_t m = 0; m < ki.size(); ++m) qij *= pij(ki[m], kj[m], beta);
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
   NumericMatrix A = parms["A"];
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
   NumericMatrix A = parms["A"];
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
 
