/*
 * File: emp_null_database
 * Description: Header for emp_null_database.cpp file
 *
 * Author: Omar Ahmed
 * Start Date: March 27, 2022
 */

#ifndef _EMPNULLDATABASE_H
#define _EMPNULLDATABASE_H

#include <string>
#include <sdsl/vectors.hpp>

class EmpNullDatabase {

public:
    std::string input_file = ""; // path to input data, single file or list
    size_t num_values = 0; // number of entries in the database
    sdsl::int_vector<> ms_null_stats; // empirical ms null statistics
    sdsl::int_vector<> pml_null_stats; // empirical pml null statistics

    EmpNullDatabase(const char* ref_file, const char* null_reads, bool use_minimizers, bool ms_built, bool pml_built);

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "");
};

#endif /* end of _EMPNULLDATABASE_H */