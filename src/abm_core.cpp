// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <random>
#include <numeric>   // std::accumulate, std::iota
#include <algorithm> // std::sample, std::min/max, std::sort
#include <iterator>  // std::begin, std::end, std::back_inserter
#include <cmath>     // std::round
#include <stdexcept> // For exception catching
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

typedef std::unordered_map<std::vector<int>, long long, VectorHasher> PopulationMap;
typedef std::unordered_map<std::vector<int>, double, VectorHasher> FitnessMap;

// --- Helper Function Definitions ---
// Parse "1.2.3..." string to vector<int>
// expected_n_chr: if > 0, validates length. If <=0, length is determined by string.
std::vector<int> parse_karyotype_string_abm(const std::string& k_str, int expected_n_chr = -1) {
  std::vector<int> cn;
  std::stringstream ss(k_str);
  std::string segment;
  
  if (k_str.empty()) {
    if (expected_n_chr == 0) return cn; // Valid empty karyotype if explicitly expected
    Rcpp::warning("Empty karyotype string provided when non-empty expected.");
    return {}; // Return empty to signify error
  }
  
  while(std::getline(ss, segment, '.')) {
    try {
      int val = std::stoi(segment);
      if (val <= 0) { // Assuming karyotype counts must be positive
        Rcpp::warning("Non-positive count (%d) in karyotype string: %s. Skipping.", val, k_str.c_str());
        return {}; 
      }
      cn.push_back(val);
    } catch (const std::invalid_argument& ia) {
      Rcpp::warning("Non-integer segment '%s' in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
      return {};
    } catch (const std::out_of_range& oor) {
      Rcpp::warning("Integer segment '%s' out of range in karyotype string: %s. Skipping.", segment.c_str(), k_str.c_str());
      return {};
    }
  }
  
  if (expected_n_chr > 0 && static_cast<int>(cn.size()) != expected_n_chr) {
    Rcpp::warning("Karyotype string '%s' has %d types, expected %d. Skipping.", k_str.c_str(), cn.size(), expected_n_chr);
    return {};
  }
  if (cn.empty() && !k_str.empty()){ // e.g. string was "..."
    Rcpp::warning("Karyotype string '%s' parsed to zero elements. Skipping.", k_str.c_str());
    return {};
  }
  return cn;
}

// Convert vector<int> to "1.2.3..." string
std::string karyotype_to_string_abm(const std::vector<int>& cn) { // Renamed to avoid conflict if linked
  if (cn.empty()) return ""; 
  std::stringstream ss;
  for(size_t i = 0; i < cn.size(); ++i) {
    ss << cn[i] << (i == cn.size() - 1 ? "" : ".");
  }
  return ss.str();
}

// Get fitness from LUT, return default if not found
double get_fitness_abm( // Renamed
    const std::vector<int>& cn,
    const FitnessMap& fitness_lut,
    double default_fitness_value = 0.0 // Provide a default for novel types
) {
  auto it = fitness_lut.find(cn);
  if (it != fitness_lut.end()) {
    return it->second; 
  } else {
    // Rcpp::warning("Karyotype %s not in LUT, using default fitness %f", 
    //               karyotype_to_string_abm(cn).c_str(), default_fitness_value); // Can be noisy
    return default_fitness_value; 
  }
}

std::pair<std::vector<int>, std::vector<int>> generate_misseg_daughters(
    int n_errors, 
    const std::vector<int>& parent_cn,
    int n_chrom_total, // Sum of elements in parent_cn
    std::mt19937& rng_engine 
) {
  if (parent_cn.empty() || n_errors < 0 || n_chrom_total <= 0 || n_errors > n_chrom_total) { // n_errors can be 0 for faithful
    if (n_errors > 0) { // Only warn if errors were expected but inputs are bad for misseg
      // Rcpp::warning("generate_misseg_daughters: invalid inputs for missegregation.");
    }
    if (n_errors == 0) { // Faithful division
      return {parent_cn, parent_cn};
    }
    return {{},{}}; // Invalid for missegregation attempt
  }
  if (n_errors == 0) { // Faithful division, return two copies of parent
    return {parent_cn, parent_cn};
  }
  
  // Sample which k distinct chromosomes (not pairs) missegregate
  std::vector<int> individual_chrom_indices;
  individual_chrom_indices.reserve(n_chrom_total);
  for(size_t type_idx = 0; type_idx < parent_cn.size(); ++type_idx) {
    for(int k_copy = 0; k_copy < parent_cn[type_idx]; ++k_copy) {
      individual_chrom_indices.push_back(type_idx); // Store the type index
    }
  }
  // Now individual_chrom_indices contains one entry for each individual chromosome, valued by its type index.
  // e.g., for {2,1}, it's {0,0,1}
  
  if(static_cast<int>(individual_chrom_indices.size()) != n_chrom_total) {
    Rcpp::stop("Internal error in generate_misseg_daughters: n_chrom_total doesn't match expanded indices.");
  }
  
  std::vector<int> sampled_positions_for_misseg; // These are positions in the list of all chromosomes
  sampled_positions_for_misseg.resize(n_chrom_total);
  std::iota(sampled_positions_for_misseg.begin(), sampled_positions_for_misseg.end(), 0);
  
  std::vector<int> missegregating_positions; // Store the *positions* of chromosomes that will missegregate
  missegregating_positions.reserve(n_errors);
  std::sample(sampled_positions_for_misseg.begin(), sampled_positions_for_misseg.end(),
              std::back_inserter(missegregating_positions),
              n_errors, rng_engine);
  // No sort needed for positions if we just iterate through them.
  
  std::vector<int> daughter1_cn = parent_cn;
  std::vector<int> daughter2_cn = parent_cn;
  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
  
  for (int pos : missegregating_positions) {
    int chrom_type_idx = individual_chrom_indices[pos]; // Get the chromosome type that missegregates
    if (uniform_dist(rng_engine) < 0.5) {
      daughter1_cn[chrom_type_idx]++;
      daughter2_cn[chrom_type_idx]--;
    } else {
      daughter1_cn[chrom_type_idx]--;
      daughter2_cn[chrom_type_idx]++;
    }
  }
  
  bool d1_valid = true;
  bool d2_valid = true;
  for(size_t i = 0; i < parent_cn.size(); ++i) { // Loop up to parent_cn.size()
    if(daughter1_cn[i] < 0) d1_valid = false; // Counts cannot be negative
    if(daughter2_cn[i] < 0) d2_valid = false; // Changed from <=0 to <0, assuming 0 is viable if some other chrom present
  }
  
  std::pair<std::vector<int>, std::vector<int>> result;
  if (d1_valid) result.first = daughter1_cn;
  if (d2_valid) result.second = daughter2_cn;
  return result;
}


// [[Rcpp::export]]
Rcpp::List run_karyotype_abm(
    Rcpp::List initial_population_r,    
    Rcpp::List fitness_map_r,           
    double p_missegregation,
    double dt,
    int n_steps,
    long long max_population_size, // R's numeric (double) will convert if in range
    double culling_survival_fraction = 0.1,
    int record_interval = 1,
    int seed = -1
) {
  PopulationMap population;
  FitnessMap fitness_map; 
  std::mt19937 rng_engine;
  if (seed == -1) {
    std::random_device rd;
    rng_engine.seed(rd());
  } else {
    rng_engine.seed(static_cast<unsigned int>(seed));
  }
  
  int n_chr_types_sim = -1; // To be determined from first valid karyotype
  
  Rcpp::CharacterVector initial_k_names = initial_population_r.names();
  for(int i = 0; i < initial_k_names.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(initial_k_names[i]);
    double count_r = 0;
    try { count_r = Rcpp::as<double>(initial_population_r[k_str]); }
    catch (...) { /* Handle or skip */ continue; }
    
    std::vector<int> cn;
    if (n_chr_types_sim == -1) {
      cn = parse_karyotype_string_abm(k_str); 
      if (!cn.empty()) {
        n_chr_types_sim = cn.size();
      } else { Rcpp::warning("Failed to parse first karyotype '%s' to determine chromosome count.", k_str.c_str()); continue; }
    } else {
      cn = parse_karyotype_string_abm(k_str, n_chr_types_sim); 
      if (cn.empty()) { continue; }
    }
    
    if(!cn.empty()) {
      if(count_r >= 0) {
        population[cn] = static_cast<long long>(std::round(count_r));
      } else {
        population[cn] = 0LL;
      }
    } 
  }
  
  if (n_chr_types_sim == -1 && initial_k_names.size() > 0) {
    Rcpp::stop("Could not determine a consistent number of chromosome types from initial_population_r.");
  }
  if (initial_k_names.size() == 0) { // If initial_population_r was empty
    // This state could be an error, or an empty simulation.
    // For now, if fitness_map_r is also empty, it might proceed to return an empty result.
    // If fitness_map_r is not empty, n_chr_types_sim needs to be derived from it.
    if (fitness_map_r.size() > 0) {
      Rcpp::CharacterVector fm_names = fitness_map_r.names();
      std::string fm_first_k_str = Rcpp::as<std::string>(fm_names[0]);
      std::vector<int> fm_cn_first = parse_karyotype_string_abm(fm_first_k_str);
      if(!fm_cn_first.empty()) n_chr_types_sim = fm_cn_first.size();
      else Rcpp::stop("Initial population is empty and cannot determine chromosome types from fitness map's first key.");
    } else {
      Rcpp::warning("Initial population and fitness map are both empty. Returning empty results.");
      return Rcpp::List(); // Return empty list
    }
  }
  
  
  Rcpp::CharacterVector fitness_k_names = fitness_map_r.names();
  for(int i = 0; i < fitness_k_names.size(); ++i) {
    std::string k_str = Rcpp::as<std::string>(fitness_k_names[i]);
    double fitness_r = 0;
    try { fitness_r = Rcpp::as<double>(fitness_map_r[k_str]); }
    catch (...) { /* Handle or skip */ continue; }
    
    std::vector<int> cn = parse_karyotype_string_abm(k_str, n_chr_types_sim); // Validate against determined length
    if(!cn.empty()) {
      fitness_map[cn] = fitness_r;
    } 
  }
  
  Rcpp::List results_over_time;
  if (record_interval >= 1 && n_chr_types_sim > 0) { // Check n_chr_types_sim before karyotype_to_string_abm
    Rcpp::NumericVector counts_r_init;
    Rcpp::CharacterVector names_r_init;
    long long initial_total = 0;
    PopulationMap cleaned_initial_population; // To store valid initial types
    
    for(const auto& pair : population) {
      if (fitness_map.count(pair.first)) { 
        if (pair.second > 0) { 
          counts_r_init.push_back(static_cast<double>(pair.second));
          names_r_init.push_back(karyotype_to_string_abm(pair.first));
          initial_total += pair.second;
          cleaned_initial_population[pair.first] = pair.second; // Keep this one
        }
      } else {
        Rcpp::warning("Initial population contains karyotype '%s' not in fitness map. It will be ignored.", 
                      karyotype_to_string_abm(pair.first).c_str());
      }
    }
    population = std::move(cleaned_initial_population); // Update population to only valid types
    
    if (names_r_init.size() > 0) {
      counts_r_init.names() = names_r_init;
    }
    results_over_time.push_back(counts_r_init, "0"); 
    if(initial_total == 0 && population.empty()) {
      Rcpp::warning("Initial population is empty after filtering against fitness map or all counts were zero.");
    }
  } else if (record_interval >= 1 && n_chr_types_sim <= 0) {
    Rcpp::warning("Cannot record initial state because number of chromosome types is undetermined.");
  }
  
  std::poisson_distribution<long long> poisson_dist;
  std::binomial_distribution<long long> binomial_dist_ll;
  std::binomial_distribution<int> error_dist_binomial; // Renamed error_dist
  
  for (int step = 1; step <= n_steps; ++step) {
    if (population.empty()) {
      Rcpp::Rcout << "Population extinct at step " << step << std::endl;
      break;
    }
    if (n_chr_types_sim <=0) { // Should have been caught earlier
      Rcpp::stop("Internal error: n_chr_types_sim not properly set before simulation loop.");
    }
    
    PopulationMap net_changes_this_step;
    std::vector<std::vector<int>> current_karyotypes_vec; // Renamed
    current_karyotypes_vec.reserve(population.size());
    for(auto const& [cn_map_key, count_map_val] : population) { // Renamed iteration vars
      if (count_map_val > 0) current_karyotypes_vec.push_back(cn_map_key);
    }
    
    for (const auto& parent_cn : current_karyotypes_vec) {
      auto parent_it = population.find(parent_cn);
      if(parent_it == population.end() || parent_it->second <=0) continue;
      long long parent_count = parent_it->second;
      
      double fitness = get_fitness_abm(parent_cn, fitness_map); // Using default_fitness_value = 0.0 from get_fitness_abm
      // because all population members should be in fitness_map
      
      if (fitness <= 0 && parent_count > 0) { // Only process if fitness > 0 for divisions
        net_changes_this_step[parent_cn] -= parent_count; // Cells die or don't divide
        continue;
      }
      if (fitness <= 0) continue; // Skip if no chance of division
      
      double expected_divisions = static_cast<double>(parent_count) * fitness * dt;
      if (expected_divisions < 0) expected_divisions = 0; // Should not happen if fitness > 0
      
      poisson_dist.param(typename std::poisson_distribution<long long>::param_type(expected_divisions));
      long long n_divs = poisson_dist(rng_engine);
      n_divs = std::min(n_divs, parent_count); // Cannot divide more cells than exist
      
      if (n_divs <= 0) continue;
      
      int n_chrom_total_parent = 0; // Renamed from n_chrom_pairs
      for(int count_val : parent_cn) n_chrom_total_parent += count_val;
      if(n_chrom_total_parent == 0 && n_divs > 0) { // Dividing an empty karyotype?
        net_changes_this_step[parent_cn] -= n_divs; // These divisions effectively lead to loss
        continue;
      }
      
      // Determine number of missegregations for these n_divs
      long long n_faithful_divs = 0;
      if (n_chrom_total_parent == 0) { // If parent has no chromosomes, all divisions are "faithful" (0 errors)
        n_faithful_divs = n_divs;
      } else {
        error_dist_binomial.param(typename std::binomial_distribution<int>::param_type(n_divs, p_missegregation));
        long long n_misseg_events_total_across_divs = error_dist_binomial(rng_engine); // Total divisions that will have >=1 error
        n_faithful_divs = n_divs - n_misseg_events_total_across_divs;
      }
      
      net_changes_this_step[parent_cn] += (n_faithful_divs - n_divs); // Faithful add 2, lose 1; unfaithful lose 1. Net for parent.
      // Net: n_faithful_divs * 1 (parent replaced by 2 daughters, one is like parent)
      //      + (n_divs - n_faithful_divs) * (-1) (parent lost, replaced by errant daughters)
      // Simpler: parent_count -= n_divs; (these are "consumed")
      //          faithful_daughters_like_parent += n_faithful_divs;
      //          faithful_daughters_other += n_faithful_divs;
      // Let's use the logic from your original C++ for num_missegs distribution.
      // The logic below is what you had and re-implements the per-division error decision.
      
      net_changes_this_step[parent_cn] -= n_divs; // All dividing parents are removed first
      
      for (long long div_idx = 0; div_idx < n_divs; ++div_idx) {
        int k_errors_this_division = 0;
        if (n_chrom_total_parent > 0 && p_missegregation > 0) { // Only attempt binomial if possible errors
          std::binomial_distribution<int> division_error_count_dist(n_chrom_total_parent, p_missegregation); // Errors per chromosome
          k_errors_this_division = division_error_count_dist(rng_engine); // Number of chromosomes that missegregate
        }
        
        std::pair<std::vector<int>, std::vector<int>> daughters =
          generate_misseg_daughters(k_errors_this_division, parent_cn, n_chrom_total_parent, rng_engine);
        
        if (!daughters.first.empty() && fitness_map.count(daughters.first)) {
          net_changes_this_step[daughters.first]++;
        }
        if (!daughters.second.empty() && fitness_map.count(daughters.second)) {
          net_changes_this_step[daughters.second]++;
        }
      }
    } 
    
    for(const auto& pair_change : net_changes_this_step) { // Renamed pair
      if (fitness_map.count(pair_change.first)) { // Ensure type is known (might be redundant)
        population[pair_change.first] += pair_change.second;
      }
    }
    
    long long current_total_pop = 0; // Renamed from next_total_pop
    for (auto it = population.begin(); it != population.end(); ) {
      if (it->second <= 0 || !fitness_map.count(it->first)) {
        it = population.erase(it);
      } else {
        current_total_pop += it->second;
        ++it;
      }
    }
    
    if (max_population_size > 0 && current_total_pop > max_population_size) {
      double sampling_fraction = culling_survival_fraction; 
      // Rcpp::Rcout << "Step " << step << ": Population " << current_total_pop // Optional verbose logging
      //             << " exceeded cap " << max_population_size
      //             << ". Culling to approx. " << static_cast<long long>(round(current_total_pop * sampling_fraction)) << " cells." << std::endl;
      
      PopulationMap sampled_population;
      if (sampling_fraction > 0.0 && sampling_fraction <= 1.0) {
        for (auto const& [cn_sample, count_sample] : population) { // Renamed iteration vars
          binomial_dist_ll.param(typename std::binomial_distribution<long long>::param_type(count_sample, sampling_fraction));
          long long sampled_count = binomial_dist_ll(rng_engine);
          if (sampled_count > 0) {
            sampled_population[cn_sample] = sampled_count;
          }
        }
        population = std::move(sampled_population); 
      } else { 
        // Rcpp::Rcout << "Step " << step << ": Sampling fraction is zero or invalid. Population culled entirely." << std::endl; // Optional
        population.clear(); 
      }
      // Recalculate current_total_pop after culling for accurate reporting if needed immediately
      current_total_pop = 0;
      for(const auto& pair_recalc : population) current_total_pop += pair_recalc.second;
    } 
    
    // const int report_freq = std::max(1, n_steps / 10); // Reporting logic can be kept if desired
    // if (step % report_freq == 0 || step == n_steps) { ... }
    
    if (step % record_interval == 0 || step == n_steps) {
      Rcpp::NumericVector counts_r_step; // Renamed
      Rcpp::CharacterVector names_r_step; // Renamed
      if (n_chr_types_sim > 0) { // Only proceed if n_chr_types_sim is valid
        for(const auto& pair_record : population) { // Renamed
          counts_r_step.push_back(static_cast<double>(pair_record.second));
          names_r_step.push_back(karyotype_to_string_abm(pair_record.first));
        }
        if(counts_r_step.length() > 0) counts_r_step.names() = names_r_step;
      }
      results_over_time.push_back(counts_r_step, std::to_string(step));
    }
    Rcpp::checkUserInterrupt();
  } 
  return results_over_time;
}