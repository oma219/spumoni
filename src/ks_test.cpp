/*
 * File: ks_test.cpp
 * Description: Implements the Kolmogorov-Smirnov statistical 
 *              test in order to differentiate betweeen similar
 *              groups of matching statistics and non-similar ones.
 *
 * Author: Omar Ahmed
 * Start Date: April 17, 2022
 */

#include <ks_test.hpp>
#include <algorithm>
#include <emp_null_database.hpp>

KSTest::KSTest(std::string ref_file, output_type result_type, bool write_report, std::ofstream& out): bin_size(75), stat_type(result_type) {
    /* main constructor for KSTest object */
    if (!write_report) return;

    // load the null empirical database
    std::string null_db_path = ref_file;
    if (stat_type == MS) {null_db_path += ".msnulldb";}
    else {null_db_path += ".pmlnulldb";}
    
    std::ifstream in(null_db_path);
    null_db.load(in);
    in.close();

    // find the average null statistic
    size_t sum_stat = 0;
    for (size_t i = 0; i < null_db.num_values; i++) {sum_stat += null_db.null_stats[i];}
    mean_null_stat = sum_stat/null_db.num_values;

    // write out the header columns to report
    out << std::setw(20) << std::left << "read id:"
        << std::setw(15) << std::left << "status:" 
        << std::setw(15) << std::left << "avg ks-stat:" 
        << std::setw(12) << std::left << "above thr:"
        << std::setw(12) << std::left << "below thr:" << std::endl;
}

std::vector<double> KSTest::compute_cdf(std::vector<size_t> stats, size_t max_stat) {
    /* generates the cdf for the list of sorted statistics, the x-value range from 0 to max_stat */
    std::vector<double> cdf_vals;
    size_t pos = 0; // current position in dataset
    size_t curr_stat = 0; // current MS value compute CDF for
    size_t total_length = stats.size(); // bin size 

    while (curr_stat <= max_stat) {
        size_t num_values = 0;

        // Check for statistics with same value
        while(pos < total_length && stats[pos] == curr_stat) {
            pos++;
            num_values++;
        }
        cdf_vals.insert(cdf_vals.end(), 1, pos/(total_length+0.0));
        curr_stat++;
    }
    return cdf_vals;
}

double KSTest::run_test(std::vector<size_t> pos_stats, std::vector<size_t> null_stats) {
    /* compute the ks-test statistic and return */

    // sort the values in both distributions
    std::sort(pos_stats.begin(), pos_stats.end());
    std::sort(null_stats.begin(), null_stats.end());

    // compute cdfs for both empirical distributions
    size_t max_stat = std::max(pos_stats.back(), null_stats.back());
    std::vector<double> pos_cdf = compute_cdf(pos_stats, max_stat);
    std::vector<double> null_cdf = compute_cdf(null_stats, max_stat);

    // find the maximal cdf difference
    double ks_stat = 0.0;
    auto null_it = null_cdf.begin();
    for (auto pos_it = pos_cdf.begin(); pos_it < pos_cdf.end(); pos_it++, null_it++) {
        ks_stat =std::max(std::abs(*pos_it - *null_it), ks_stat);
        if (*pos_it >= 1.0 || *null_it >= 1.0) break;
    }
    return ks_stat;
}

void KSTest::run_kstest(const char* read_id, std::vector<size_t> pos_stats, std::ofstream& out) {
    /* runs the KS-test using empirical null database, and saves results */
    size_t curr_start_pos = 0;
    std::vector<double> ks_list;

    while (curr_start_pos < pos_stats.size()) {
        // choose a random section of null database
        size_t null_pos = rand() % (this->null_db.num_values - this->bin_size);

        // build the two vectors: positive statistics and null statistics
        std::vector<size_t> curr_pos_bin, curr_null_bin;
        size_t end = (curr_start_pos + this->bin_size < pos_stats.size()) ? (curr_start_pos+this->bin_size) : pos_stats.size();
        curr_pos_bin.assign(pos_stats.begin()+curr_start_pos, pos_stats.begin()+end);

        size_t region_size = end - curr_start_pos;
        for (size_t i = null_pos; i < (null_pos + region_size); i++){
            curr_null_bin.push_back(this->null_db.null_stats[i]);
        }

        // run ks-stat for this region of read
        double curr_ks_stat = this->run_test(curr_pos_bin, curr_null_bin);
        ks_list.push_back(curr_ks_stat);
        curr_start_pos += this->bin_size;
    }

    // determine how read should be classified
    size_t num_bin_above_thr = 0;
    double threshold = (this->stat_type == MS) ? KS_STAT_MS_THR : KS_STAT_PML_THR;
    for (size_t i = 0; i < ks_list.size(); i++) {
        if (ks_list[i] >= threshold) num_bin_above_thr++;
    }
    bool read_found = (num_bin_above_thr/(ks_list.size()+0.0) >= 0.50);

    // write out to report ...
    double sum_ks_stats = 0.0;
    std::for_each(ks_list.begin(), ks_list.end(), [&] (double n) {sum_ks_stats += n;});

    std::string status = (read_found) ? "FOUND" : "NOT_PRESENT";

    out.precision(3);
    out << std::setw(20) << std::left << read_id
        << std::setw(15) << std::left << status 
        << std::setw(15) << std::left << (sum_ks_stats/ks_list.size()) 
        << std::setw(12) << std::left << num_bin_above_thr
        << std::setw(12) << std::left << (ks_list.size() - num_bin_above_thr)
        << std::endl;
    
}