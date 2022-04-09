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
                                 bool ms_built, bool pml_built, output_type index_type) {
    /* Builds the null database of MS/PML and saves it */
    
    this->input_file = std::string(ref_file);
    this->stat_type = index_type;
    
    // Generate those null statistics depending on the index type
    std::vector<size_t> output_stats;
    if (stat_type == MS) {generate_null_ms_statistics(this->input_file, std::string(null_reads), output_stats, use_minimizers);}
    else if (stat_type == PML) {generate_null_pml_statistics(this->input_file, std::string(null_reads), output_stats, use_minimizers);}

    // Determine the size needed for each vector
    uint32_t max_stat_width = 1;
    auto max_null_stat = std::max_element(output_stats.begin(), output_stats.end()); 
    max_stat_width = std::max(static_cast<int>(std::ceil(std::log2(*max_null_stat))), 1);

    DBG_ONLY("Maximum null statistic: %i", *max_null_stat);
    DBG_ONLY("Number of bits used per null statistic: %i", max_stat_width);

    // Initialize the attributes to be saved in null database
    this->num_values = output_stats.size();
    this->null_stats = sdsl::int_vector<> (this->num_values, 0, max_stat_width);

    // Separated the loops in case we only build one index
    for (size_t i = 0; i < output_stats.size(); i++) {
        this->null_stats[i] = output_stats[i];
    }
}

size_t EmpNullDatabase::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) {
    /* Write the empirical null database to a file on disk */
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written_bytes = 0;

    out.write((char *)&this->num_values, sizeof(this->num_values));
    written_bytes += sizeof(this->num_values);
    written_bytes += this->null_stats.serialize(out, child, "null_stats");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}