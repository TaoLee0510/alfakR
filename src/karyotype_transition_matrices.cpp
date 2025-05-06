// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <numeric> // std::accumulate
#include <cmath>   // std::pow
#include <stdexcept> // std::runtime_error
// Removed: const int N_CHROMOSOME_TYPES = 22;

// --- Hash function for using std::vector<int> as map key ---
struct VectorHasher {
  std::size_t operator()(const std::vector<int>& v) const {
    std::size_t seed = v.size();
    for(int i : v) {
      seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// --- Helper Function: Parse karyotype string ---
// Returns empty vector on error, issues Rcpp::warning
std::vector<int> parse_karyotype_string_internal(const std::string& k_str, int expected_len) {
  std::vector<int> cn;
  if (expected_len > 0) { // Only reserve if we know the expected length and it's positive
    cn.reserve(expected_len);
  }
  std::stringstream ss(k_str);
  std::string segment;
  int count = 0;
  while(std::getline(ss, segment, '.')) {
    try {
      double val_d = std::stod(segment);
      if (val_d <= 0 || val_d != static_cast<int>(val_d)) {
        Rcpp::warning("Non-positive or non-integer count (%s) in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
        return {};
      }
      cn.push_back(static_cast<int>(val_d));
      count++;
    } catch (const std::invalid_argument& ia) {
      Rcpp::warning("Non-numeric segment '%s' in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
      return {};
    } catch (const std::out_of_range& oor) {
      Rcpp::warning("Numeric segment '%s' out of range in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
      return {};
    } catch (const std::exception& e) {
      Rcpp::warning("Error parsing segment '%s' in '%s': %s. Skipping.", segment.c_str(), k_str.c_str(), e.what());
      return {};
    }
  }
  // If expected_len was provided and is positive, validate against it
  if (expected_len > 0 && count != expected_len) {
    Rcpp::warning("Karyotype string '%s' has %d types, expected %d. Skipping.", k_str.c_str(), count, expected_len);
    return {};
  }
  if (count == 0 && !k_str.empty()){ // Parsed to zero elements but original string was not empty (e.g. "...")
    Rcpp::warning("Karyotype string '%s' parsed to zero elements. Skipping.", k_str.c_str());
    return {};
  }
  if (count == 0 && k_str.empty() && expected_len != 0){ // Empty string should only be valid if expected_len is 0 (unlikely for karyotypes)
    // Rcpp::warning("Empty karyotype string provided when non-zero length expected."); // This might be too noisy
    // Return empty, further checks might handle it. If expected_len > 0, it would have failed above.
  }
  return cn;
}

// --- Helper Function: Convert integer vector back to string tag ---
std::string vec_to_tag_internal(const std::vector<int>& v) {
  if (v.empty()) return "";
  std::stringstream ss;
  for(size_t i = 0; i < v.size(); ++i) {
    ss << v[i] << (i == v.size() - 1 ? "" : ".");
  }
  return ss.str();
}

enum class EntryType : int { DIAGONAL = 0, OFF_DIAGONAL = 1 };

// [[Rcpp::export]]
Rcpp::List rcpp_prepare_W_structure(Rcpp::CharacterVector k_strings) {
  if (k_strings.size() == 0) {
    Rcpp::stop("Input 'k_strings' is empty.");
  }
  
  // Determine the number of chromosome types from the first valid karyotype string
  int n_chr_types = 0;
  std::string first_k_str_for_len = "";
  for(int i = 0; i < k_strings.size(); ++i) { // Find first non-empty string
    first_k_str_for_len = Rcpp::as<std::string>(k_strings[i]);
    if (!first_k_str_for_len.empty()) break;
  }
  if (first_k_str_for_len.empty() && k_strings.size() > 0) Rcpp::stop("All input karyotype strings are empty.");
  
  
  std::stringstream ss_first(first_k_str_for_len);
  std::string segment_first;
  while(std::getline(ss_first, segment_first, '.')) {
    n_chr_types++;
  }
  if(n_chr_types == 0) { // This would happen if first_k_str_for_len was e.g. "..." or malformed to have no segments
    Rcpp::stop("Could not determine number of chromosome types from the first karyotype string: '%s'. Ensure format is like '1.2.3'.", first_k_str_for_len.c_str());
  }
  
  std::vector<std::vector<int>> karyotype_vecs;
  std::vector<std::string> valid_k_names;
  std::unordered_map<std::vector<int>, int, VectorHasher> karyo_to_idx;
  
  int N = 0; // Number of *valid* and *unique* karyotypes
  for (int i = 0; i < k_strings.size(); ++i) {
    std::string current_k_str = Rcpp::as<std::string>(k_strings[i]);
    std::vector<int> parsed_vec = parse_karyotype_string_internal(current_k_str, n_chr_types); // Use determined n_chr_types
    
    if (!parsed_vec.empty()) {
      if (karyo_to_idx.find(parsed_vec) == karyo_to_idx.end()) {
        karyotype_vecs.push_back(parsed_vec);
        valid_k_names.push_back(current_k_str); // Store original string for dimnames
        karyo_to_idx[parsed_vec] = N; 
        N++;
      } else {
        Rcpp::warning("Duplicate karyotype string found and ignored after validation: %s", current_k_str.c_str());
      }
    }
  }
  
  if (N == 0) {
    Rcpp::stop("No valid karyotypes found after parsing and validation against consistent length of %d.", n_chr_types);
  }
  
  int estimated_max_nnz = N + N * n_chr_types * 2; 
  std::vector<int> ii_vec; 
  std::vector<int> jj_vec; 
  std::vector<int> type_vec; 
  std::vector<int> Cj_vec;   
  std::vector<int> ak_vec;   
  
  ii_vec.reserve(estimated_max_nnz);
  // ... (rest of reserve calls) ...
  
  for (int j = 0; j < N; ++j) { 
    const std::vector<int>& a = karyotype_vecs[j]; 
    if (a.size() != n_chr_types) { // Should not happen if parsing is correct
      Rcpp::stop("Internal error: Karyotype vector size mismatch for parent %s.", valid_k_names[j].c_str());
    }
    int Cj = std::accumulate(a.begin(), a.end(), 0); 
    
    ii_vec.push_back(j);
    jj_vec.push_back(j);
    type_vec.push_back(static_cast<int>(EntryType::DIAGONAL));
    Cj_vec.push_back(Cj);
    ak_vec.push_back(0); 
    
    if (Cj > 0) { 
      for (int k_chrom = 0; k_chrom < n_chr_types; ++k_chrom) { // Use n_chr_types
        int ak = a[k_chrom]; 
        if (ak > 0) { 
          std::vector<int> plus = a;
          plus[k_chrom] = ak + 1;
          auto it_p = karyo_to_idx.find(plus);
          if (it_p != karyo_to_idx.end()) { 
            int idp = it_p->second; 
            ii_vec.push_back(j);
            jj_vec.push_back(idp);
            type_vec.push_back(static_cast<int>(EntryType::OFF_DIAGONAL));
            Cj_vec.push_back(Cj);
            ak_vec.push_back(ak); 
          }
          
          if (ak > 1) { 
            std::vector<int> minus = a;
            minus[k_chrom] = ak - 1;
            auto it_m = karyo_to_idx.find(minus);
            if (it_m != karyo_to_idx.end()) { 
              int idm = it_m->second; 
              ii_vec.push_back(j);
              jj_vec.push_back(idm);
              type_vec.push_back(static_cast<int>(EntryType::OFF_DIAGONAL));
              Cj_vec.push_back(Cj);
              ak_vec.push_back(ak); 
            }
          }
        } 
      } 
    } 
  } 
  
  int nnz = ii_vec.size(); 
  std::vector<int> i_1based(nnz);
  std::vector<int> j_1based(nnz);
  for(int k=0; k < nnz; ++k) {
    i_1based[k] = ii_vec[k] + 1;
    j_1based[k] = jj_vec[k] + 1;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("i") = i_1based,
    Rcpp::Named("j") = j_1based,
    Rcpp::Named("type") = type_vec,
    Rcpp::Named("Cj") = Cj_vec,
    Rcpp::Named("ak") = ak_vec,
    Rcpp::Named("dims") = Rcpp::IntegerVector::create(N, N),
    Rcpp::Named("dimnames") = Rcpp::List::create(Rcpp::CharacterVector(valid_k_names.begin(), valid_k_names.end()), 
                Rcpp::CharacterVector(valid_k_names.begin(), valid_k_names.end()))
  );
}

// rcpp_update_W_values remains the same as it doesn't use N_CHROMOSOME_TYPES directly.
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_update_W_values(Rcpp::List structure, double p) {
  if (p < 0.0 || p > 1.0) {
    Rcpp::stop("Mis-segregation probability 'p' must be between 0.0 and 1.0.");
  }
  std::vector<int> type_vec = Rcpp::as<std::vector<int>>(structure["type"]);
  std::vector<int> Cj_vec = Rcpp::as<std::vector<int>>(structure["Cj"]);
  std::vector<int> ak_vec = Rcpp::as<std::vector<int>>(structure["ak"]);
  int nnz = type_vec.size();
  if (nnz == 0) {
    return Rcpp::NumericVector(0);
  }
  Rcpp::NumericVector x_vec(nnz);
  double p_complement = 1.0 - p;
  for (int k = 0; k < nnz; ++k) {
    int type_code = type_vec[k];
    int Cj = Cj_vec[k];
    int ak = ak_vec[k];
    double value = 0.0;
    if (type_code == static_cast<int>(EntryType::DIAGONAL)) {
      if (Cj == 0) {
        value = 2.0;
      } else if (p == 1.0) { // (1-p) is 0
        value = 0.0; 
      } else {
        value = 2.0 * std::pow(p_complement, Cj);
      }
    } else if (type_code == static_cast<int>(EntryType::OFF_DIAGONAL)) {
      if (p == 0.0 || Cj == 0) {
        value = 0.0;
      } else if (Cj == 1) {
        value = static_cast<double>(ak) * p; // (1-p)^(1-1)=1
      } else if (p == 1.0) { // (1-p) is 0, Cj-1 > 0
        value = 0.0;
      } else {
        value = static_cast<double>(ak) * p * std::pow(p_complement, Cj - 1);
      }
    } else {
      Rcpp::stop("Internal error: Unknown entry type code encountered: %d", type_code);
    }
    x_vec[k] = (value < 0.0) ? 0.0 : value;
  }
  return x_vec;
}
