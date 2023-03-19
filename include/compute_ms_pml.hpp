 /*
  * File: compute_ms_pml.hpp
  * Description: Header file for compute_ms_pml.cpp
  *
  * Authors: Omar Ahmed
  * Start Date: October 13, 2021
  */

#ifndef COMPUTE_MS_PML_H
#define COMPUTE_MS_PML_H

#include <emp_null_database.hpp>

/* Function Declarations */
int run_spumoni_ms_main(SpumoniRunOptions* run_opts);
int run_spumoni_main(SpumoniRunOptions* run_opts);
std::pair<size_t, size_t> build_spumoni_ms_main(std::string ref_file, std::string output_prefix);
std::pair<size_t, size_t> build_spumoni_main(std::string ref_file, std::string output_prefix);
void generate_null_ms_statistics(std::string ref_file, std::string pattern_file, std::vector<size_t>& ms_stats,
                                 bool min_digest, bool use_promotions, bool use_dna_letters, size_t k, size_t w);
void generate_null_pml_statistics(std::string ref_file, std::string pattern_file, std::vector<size_t>& pml_stats,
                                  bool min_digest, bool use_promotions, bool use_dna_letters, size_t k, size_t w);
void generate_null_ms_statistics_for_general_text(std::string ref_file, std::string pattern_file, std::vector<size_t>& ms_stats);
void generate_null_pml_statistics_for_general_text(std::string ref_file, std::string pattern_file, std::vector<size_t>& pml_stats);
void find_threshold_based_on_null_pml_distribution(const char* ref_file, const char* null_reads, bool use_minimizers,
                                                   bool use_promotions, bool use_dna_letters, size_t k, size_t w, 
                                                   EmpNullDatabase& null_db, size_t bin_width);
void find_threshold_based_on_null_ms_distribution(const char* ref_file, const char* null_reads, bool use_minimizers,
                                                  bool use_promotions, bool use_dna_letters, size_t k, size_t w, 
                                                  EmpNullDatabase& null_db, size_t bin_width);
std::pair<ulint, ulint> get_bwt_stats(std::string ref_file, size_t type);

#endif /* End of include of COMPUTE_MS_PML_H */