
 /*
  * File: doc_array.hpp
  * Description: Header for document array ...
  *
  * Author: Omar Ahmed
  * Start Date: January 14, 2021
  */

#ifndef _DOCARRAY_H
#define _DOCARRAY_H

#include <string>
#include <vector>
#include <sdsl/vectors.hpp>

class DocumentArray {

public:
    std::string ref_file; // path to input data
    std::vector<size_t> seq_lengths; // sequence lengths in FASTA file
    size_t num_entries; // number of document array entries (i.e number of BWT runs)
    sdsl::int_vector<> start_runs_doc; // represents the start of runs document array
    sdsl::int_vector<> end_runs_doc; // represents the end of runs document array

    DocumentArray(){} // Constructor used to load an existing document array
    DocumentArray(std::string file_path, size_t length, size_t num_runs); // Main constructor
    
    void load_seq_boundaries();
    void print_statistics();
    static size_t grab_file_size(std::string file_path);
    static std::vector<size_t> read_samples(std::string file_path);
    static size_t binary_search_for_pos(std::vector<size_t> end_pos, size_t sample_pos);
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "");
    void load(std::istream& in);
}; // end of DocumentArray Class

#endif /* end of _DOCARRAY_H */