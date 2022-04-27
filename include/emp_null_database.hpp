/*
 * File: emp_null_database.hpp
 * Description: Header for emp_null_database.cpp file
 *
 * Author: Omar Ahmed
 * Start Date: March 27, 2022
 */

#ifndef _EMPNULLDATABASE_H
#define _EMPNULLDATABASE_H

#include <string>
#include <spumoni_main.hpp>
#include <sdsl/vectors.hpp>

class EmpNullDatabase {

public:
    std::string input_file = ""; // path to input data, single file or list
    size_t num_values = 0; // number of entries in the database
    output_type stat_type = NOT_CHOSEN; // either MS or PML
    sdsl::int_vector<> null_stats; // empirical null statistics
    
    EmpNullDatabase(){} // constructor used for loading
    EmpNullDatabase(const char* ref_file, const char* null_reads, bool use_minimizers, output_type stat_type);

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "");
    void load(std::istream& in);
};

#endif /* end of _EMPNULLDATABASE_H */