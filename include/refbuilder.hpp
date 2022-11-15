/*
 * File: refbuilder.hpp
 * Description: Header for refbuilder.cpp file
 *
 * Author: Omar Ahmed
 * Start Date: February 12, 2021
 */

#ifndef _REFBUILD_H
#define _REFBUILD_H

#include <string>

class RefBuilder {

public:
    std::string input_file = ""; // path to input data, single file or list
    std::string null_read_file = ""; // path to the null read file
    bool using_doc = false; // user intends to build a document array
    bool using_list = false; // user provided a filelist

    RefBuilder(const char* ref_file, const char* list_file, const char* output_file, const char* null_reads, 
               bool build_doc, bool file_list, bool use_minimizers,
               bool use_promotions, bool use_dna_letters,
               size_t k, size_t w);

    const char* get_ref_path();
    const char* get_null_readfile();
    static std::string parse_null_reads(const char* ref_file);
    static std::string parse_null_reads_from_general_text(const char* ref_file);
    static std::string build_reference(const char* ref_file, const char* output_path, bool use_promotions, bool use_dna_letters, 
                                        size_t k, size_t w);
};// end of RefBuilder class

#endif /* end of _REFBUILD_H */