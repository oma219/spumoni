/*
 * File: ks_test.hpp
 * Description: Header for ks_test.cpp
 *
 * Author: Omar Ahmed
 * Start Date: April 17, 2022
 */

#ifndef _KSTEST_H
#define _KSTEST_H

#include <emp_null_database.hpp>
#include <spumoni_main.hpp>
#include <vector>

class KSTest {

public:
    size_t bin_size;
    size_t mean_null_stat;
    output_type stat_type;
    EmpNullDatabase null_db;

    KSTest(std::string ref_file, output_type result_type, bool write_report, std::ofstream& out);
    KSTest(EmpNullDatabase& null_db, output_type result_type);
    
    double run_test(std::vector<size_t> pos_stats, std::vector<size_t> null_stats);
    static inline std::vector<double> compute_cdf(std::vector<size_t> stats, size_t max_stat);
    std::vector<double> run_kstest(std::vector<size_t> pos_stats);
    double get_threshold();
};

#endif /* end of _KSTEST_H */