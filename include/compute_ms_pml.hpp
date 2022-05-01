 /*
  * File: compute_ms_pml.hpp
  * Description: Header file for compute_ms_pml.cpp
  *
  * Authors: Omar Ahmed
  * Start Date: October 13, 2021
  */

#ifndef COMPUTE_MS_PML_H
#define COMPUTE_MS_PML_H

/* Function Declarations */
int run_spumoni_ms_main(SpumoniRunOptions* run_opts);
int run_spumoni_main(SpumoniRunOptions* run_opts);
std::pair<size_t, size_t> build_spumoni_ms_main(std::string ref_file);
std::pair<size_t, size_t> build_spumoni_main(std::string ref_file);
void generate_null_ms_statistics(std::string ref_file, std::string pattern_file, std::vector<size_t>& ms_stats,
                                 bool min_digest, bool use_promotions, bool use_dna_letters, size_t k, size_t w);
void generate_null_pml_statistics(std::string ref_file, std::string pattern_file, std::vector<size_t>& pml_stats,
                                  bool min_digest, bool use_promotions, bool use_dna_letters, size_t k, size_t w);
std::pair<ulint, ulint> get_bwt_stats(std::string ref_file, size_t type);

#endif /* End of include of COMPUTE_MS_PML_H */