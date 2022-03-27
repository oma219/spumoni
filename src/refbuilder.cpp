/*
 * File: refbuilder.hpp
 * Description: Implementation of RefBuilder class
 *              which handles the generation of fasta index, 
 *              and digestion of the input files.
 *
 * Author: Omar Ahmed
 * Start Date: February 12, 2021
 */

#include <refbuilder.hpp>
#include <spumoni_main.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <zlib.h>  
#include <encoder.h>

// Commented out for including encoder.h - Omar
//#include <kseq.h>
//KSEQ_INIT(gzFile, gzread)

/* Complement Table from: https://github.com/lh3/seqtk/blob/master/seqtk.c */
char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

RefBuilder::RefBuilder(const char* ref_file, const char* list_file, const char* output_dir, 
                       bool build_doc, bool file_list, bool use_minimizers): using_doc(build_doc), using_list(file_list) {
    /* Performs the needed operations to generate a single input file. */
    char ch;
    std::string output_path(output_dir);

    if ((ch = output_path.back()) != '/') {output_path += "/";}
    if (use_minimizers) output_path += "spumoni_full_ref.bin";
    else output_path += "spumoni_full_ref.fa";

    // Verify every file in the list is valid 
    std::string line = "";
    size_t curr_id = 0, member_num = 0;

    std::ifstream input_fd (list_file, std::ifstream::in);
    std::vector<std::string> input_files;
    std::vector<size_t> document_ids;

    while (std::getline(input_fd, line)) {
        auto word_list = split(line, ' ');

        // Makes sure the files are valid, and there is data to begin with.
        ASSERT((word_list.size() >= 1), "Input file-list does not have expected structure.");
        if (!is_file(word_list[0])) {
            FATAL_ERROR("The following path in the input list is not valid: %s", word_list[0].data());}
        if (!endsWith(word_list[0], ".fa") && !endsWith(word_list[0], ".fasta") && !endsWith(word_list[0], ".fna")) {
            FATAL_ERROR("The following input-file is not a FASTA file: %s", word_list[0].data());}
        input_files.push_back(word_list[0]);

        // If the user wants a document array, there needs to be a ID in second column
        if (using_doc){
            ASSERT((word_list.size() >= 2), "If you want to build a document array, you need a second column with IDs.");

            if (!is_integer(word_list[1])) FATAL_ERROR("A document ID in the file_list is not an integer: %s", word_list[1].data());  
            if (member_num == 0 && std::stoi(word_list[1]) != 1) FATAL_ERROR("The first ID in file_list must be 1");
            
            if (std::stoi(word_list[1]) ==  curr_id || std::stoi(word_list[1]) == (curr_id+1)) {
                if (std::stoi(word_list[1]) == (curr_id+1)) {curr_id+=1;}
                document_ids.push_back(std::stoi(word_list[1]));
            } else {FATAL_ERROR("The IDs in the file_list must be staying constant or increasing by 1.");}
        }
        member_num += 1;
    }
    if (using_doc) {
        ASSERT((document_ids.size() == input_files.size()), "Issue with file-list parsing occurred.");
        if (document_ids.back() == 1) {
            FATAL_WARNING("If you only have one class ID, you should not build a document array.");}
    }

    // Define variables needed for minimizer digest
    bns::RollingHasher<uint8_t> rh(4, false, bns::DNA, 12);
    std::ofstream output_fd (output_path, std::ofstream::out);
    bool hp_compress = true;

    std::string mseq = "";
    std::vector<uint8_t> mseq_vec;

    // Define lambda function to take in a DNA sequence and return minimizer sequence
    auto get_minimizer_seq = [&](std::string seq, size_t seq_length) {
                mseq = ""; mseq_vec.clear();
                rh.for_each_uncanon([&](auto x) { // Extracts all minimizers and stores in mseq
                if(mseq_vec.empty() || !hp_compress || mseq_vec.back() != x) {
                    x = (x > 2) ? x : (x + 3);
                    mseq_vec.push_back(x);
                    mseq += x;
                }
            }, seq.data(), seq_length);
            return mseq;
    };

    // Process each input file, and store it forward + reverse complement sequence
    gzFile fp;
    kseq_t* seq;
    std::vector<size_t> seq_lengths;

    curr_id = 1;
    size_t curr_id_seq_length = 0;
    for (auto iter = input_files.begin(); iter != input_files.end(); ++iter) {

        fp = gzopen((*iter).data(), "r");
        seq = kseq_init(fp);
        size_t iter_index = iter-input_files.begin();

        while (kseq_read(seq)>=0) {
            // Get forward seq, and print it
			for (size_t i = 0; i < seq->seq.l; ++i)
				seq->seq.s[i] = std::toupper(seq->seq.s[i]);
            
            // Convert forward_seq to minimizers by default, or save as DNA if asked
            if (use_minimizers) {
                mseq = get_minimizer_seq(seq->seq.s, seq->seq.l);
                output_fd << mseq;
                curr_id_seq_length += mseq.length();
            } else {
                output_fd << '>' << seq->name.s << '\n' << seq->seq.s << '\n';
                curr_id_seq_length += seq->seq.l;
            }

            // Get reverse complement, and print it
            // Based on seqtk reverse complement code, that does it 
            // in place. (https://github.com/lh3/seqtk/blob/master/seqtk.c)
            int c0, c1;
			for (size_t i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
				c0 = comp_tab[(int)seq->seq.s[i]];
				c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
				seq->seq.s[i] = c1;
				seq->seq.s[seq->seq.l - 1 - i] = c0;
			}
			if (seq->seq.l & 1) // complement the remaining base
				seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];


            // Convert rev_comp seq to minimizers by default, otherwise DNA
            if (use_minimizers) {
                mseq = get_minimizer_seq(seq->seq.s, seq->seq.l);
                output_fd << mseq;
                curr_id_seq_length += mseq.length();
            } else {
                output_fd << '>' << seq->name.s << "_rev_comp" << '\n' << seq->seq.s << '\n';
                curr_id_seq_length += seq->seq.l;
            }
        }
        kseq_destroy(seq);
        gzclose(fp);

        if (using_doc) {
            // Check if we are transitioning to a new group
            if (iter_index < document_ids.size()-1 && document_ids[iter_index] != document_ids[iter_index+1]){
                seq_lengths.push_back(curr_id_seq_length);
                curr_id += 1; curr_id_seq_length = 0;
            // If last file, output current sequence length
            } else if (iter_index == document_ids.size()-1){
                seq_lengths.push_back(curr_id_seq_length);
                curr_id_seq_length = 0;}
        }
    }
    output_fd.close(); 
    
    // Assign the full reference to the attribute
    input_file = output_path;
    if (!using_doc) return;
    ASSERT((curr_id == document_ids.back()), "Issue with building the FASTA document index.");

    // Write out the FASTA document index
    std::ofstream output_fdi (output_path + ".fdi", std::ofstream::out);
    for (auto iter = seq_lengths.begin(); iter != seq_lengths.end(); ++iter) {
        size_t iter_index = iter - seq_lengths.begin() + 1;
        output_fdi << "group_" << iter_index << '\t' << *iter << '\n';
    }
    output_fdi.close();
}

const char* RefBuilder::get_ref_path() {
    /* Accessor - return path to full reference that was build */
    return input_file.data();
}


