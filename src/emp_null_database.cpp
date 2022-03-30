/*
 * File: emp_null_database.cpp
 * Description: Implements the code that will generate null
 *              MS/PMLs for random substrings pulled from the 
 *              input database, and stores for later computation.
 *
 * Authors: Omar Ahmed
 * Start Date: March 27th, 2022
 */

#include <emp_null_database.hpp>
#include <spumoni_main.hpp>
#include <compute_ms_pml.hpp>
#include <string>
#include <math.h> 
#include <algorithm>
#include <sdsl/vectors.hpp>

EmpNullDatabase::EmpNullDatabase(const char* ref_file, const char* null_reads, bool use_minimizers, 
                                 bool ms_built, bool pml_built) {
    /* Builds the null database of MS/PML and saves it */

    // Generate those null MS/PMLs
    std::vector<size_t> ms_stats, pml_stats;
    if (ms_built) {generate_null_ms_statistics(std::string(ref_file), std::string(null_reads), ms_stats, true);}
    if (pml_built) {generate_null_pml_statistics(std::string(ref_file), std::string(null_reads), pml_stats, true);}

    for (size_t i = 0; i < ms_stats.size(); i++) {
        std::cout << ms_stats[i] << std::endl;
    }

    // Determine the size needed for each vector
    uint32_t max_ms_width = 1, max_pml_width = 1;
    if (ms_built) {
        auto max_ms_stat = std::max_element(ms_stats.begin(), ms_stats.end()); 
        max_ms_width = std::max(static_cast<int>(std::ceil(std::log2((*max_ms_stat + 0.0001)))), 1);
        DBG_ONLY("Maximum null ms: %i", *max_ms_stat);
        DBG_ONLY("Number of bits used per null ms: %i", max_ms_width);
    }
    if (pml_built) {
        auto max_pml_stat = std::max_element(pml_stats.begin(), pml_stats.end());
        max_pml_width = std::max(static_cast<int>(std::ceil(std::log2((*max_pml_stat)))), 1);
        DBG_ONLY("Maximum null pml: %i", *max_pml_stat);
        DBG_ONLY("Number of bits used per null pml: %i", max_pml_width);
    }

    // Initialize the attributes to be saved in null database
    this->num_values = std::max(ms_stats.size(), pml_stats.size());
    this->ms_null_stats = sdsl::int_vector<> (num_values, 0, max_ms_width);
    this->pml_null_stats = sdsl::int_vector<> (num_values, 0, max_pml_width);

    // Separated the loops in case we only build one index
    for (size_t i = 0; i < ms_stats.size() && ms_built; i++) {
        this->ms_null_stats[i] = ms_stats[i];
    }
    for (size_t i = 0; i < pml_stats.size() && pml_built; i++) {
        this->pml_null_stats[i] = pml_stats[i];
    }
}

size_t EmpNullDatabase::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) {
    /* Write the empirical null database to a file on disk */
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written_bytes = 0;

    out.write((char *)&this->num_values, sizeof(this->num_values));
    written_bytes += sizeof(this->num_values);

    written_bytes += this->ms_null_stats.serialize(out, child, "ms_null_stats");
    written_bytes += this->pml_null_stats.serialize(out, child, "pml_null_stats");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}