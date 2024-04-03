/*
 * File: ms_rle_string.hpp
 * Description: Extension of the r-index rle_string to compute 
 *              matching statistics.
 *
 * Authors: Massimiliano Rossi
 * Start Date: July 10, 2021
 */

#ifndef _MS_RLE_STRING_HH
#define _MS_RLE_STRING_HH

#include <common.hpp>
#include <rle_string.hpp>

template <class sparse_bitvector_t = ri::sparse_sd_vector, //predecessor structure storing run length
          class string_t = ri::huff_string>                //run heads
class ms_rle_string : public ri::rle_string<sparse_bitvector_t, string_t>
{
public:
    static const u_char TERMINATOR = 1;

    /*
     * Constructors for ms_rle_string:
     *      - input -> input string without 0x0 bytes in it
     *      - B -> block size so the main sparse bitvector has R out of B bits set (R = number of runs)
     *      - heads & lengths -> run-length encoded BWT     
     */

    ms_rle_string() : ri::rle_string<sparse_bitvector_t, string_t>() {}

    ms_rle_string(string &input, ulint B = 2) : ri::rle_string<sparse_bitvector_t, string_t>(input, B) {}

    ms_rle_string(std::ifstream &ifs, ulint B = 2) : ri::rle_string<sparse_bitvector_t, string_t>(ifs, B) {}

    // Construction from run-length encoded BWT
    ms_rle_string(std::ifstream &heads, std::ifstream &lengths, ulint B = 2)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->B = B;
        auto runs_per_letter_bv = std::vector<std::vector<bool>>(256);
        vector<bool> runs_bv;

        // Reads the run heads
        string run_heads_s;
        heads.seekg(0, heads.end);
        run_heads_s.resize(heads.tellg());
        heads.seekg(0, heads.beg);
        heads.read(&run_heads_s[0], run_heads_s.size());

        size_t pos = 0;
        this->n = 0;
        this->R = run_heads_s.size();

        // Compute runs_bv and runs_per_letter_bv
        for (size_t i = 0; i < run_heads_s.size(); ++i)
        {
            size_t length;
            lengths.read((char *)&length, 5);

            uint8_t curr_ch = unsigned(run_heads_s[i]);
            if (curr_ch <= TERMINATOR) { // change 0 to 1
                run_heads_s[i] = TERMINATOR;
            }

            std::fill_n(std::back_inserter(runs_bv), length - 1, false);
            runs_bv.push_back(i % B == B - 1);

            std::fill_n(std::back_inserter(runs_per_letter_bv[curr_ch]), length - 1, false);
            runs_per_letter_bv[curr_ch].push_back(true);

            this->n += length;
        }

        // Now compact the structures
        assert(runs_bv.size() == this->n);
        ulint t = 0;
        for (ulint i = 0; i < 256; ++i)
            t += runs_per_letter_bv[i].size();
        assert(t == this->n);
        this->runs = sparse_bitvector_t(runs_bv);

        // fast direct array: char -> bitvector.
        this->runs_per_letter = vector<sparse_bitvector_t>(256);
        for (ulint i = 0; i < 256; ++i)
            this->runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);
        this->run_heads = string_t(run_heads_s);
        assert(this->run_heads.size() == this->R);
    }

    size_t number_of_runs_of_letter(uint8_t c) {
        return this->runs_per_letter[c].number_of_1();
    }

    size_t number_of_letter(uint8_t c) {
        return this->runs_per_letter[c].size();
    }

    // i-th run head
    uint8_t head_of(const size_t i) {
        assert(i<this->R);
        return this->run_heads[i];
    }

    // rank in chracters of the i-th run head
    // i.e., the number of characters c before the first character of the run.
    size_t head_rank(const size_t i, const uint8_t c) {
        assert(i < this->R);
        size_t j = this->run_heads.rank(i, c);
        if(j < 1)
            return j;
        assert(j<=i);
        return this->runs_per_letter[c].select(j-1) + 1; // j-1 because the select is 0 based
    }

    /*
     * The next three methods were added in order to support
     * the new version of the thresholds data-structure that
     * is bitvector-based.
     */

    // number of runs of character c in in position i
    size_t run_head_rank(const size_t i, const uint8_t c) {
        assert(i < this->R);
        size_t j = this->run_heads.rank(i, c);
        return j;
    }

    inline std::pair<size_t,size_t> run_and_head_rank(const size_t i, const uint8_t c) {
        assert(i < this->R);
        const size_t j = this->run_heads.rank(i, c);
        if( j < 1)
            return make_pair(j,j);
        const size_t k = this->runs_per_letter[c].select(j - 1) + 1; // j-1 because the select is 0 based
        return make_pair(j, k);
    }

    // Select the i-th run of c
    size_t run_head_select(const size_t i, const uint8_t c) {
        assert(i < this->R and i > 0);
        return this->run_heads.select(i - 1, c);
    }

    // Serialize the structure to the ostream object
    ulint serialize(std::ostream &out) {
        return ri::rle_string<sparse_bitvector_t, string_t>::serialize(out);
    }

     // Load the structure from istream object
    void load(std::istream &in) {
        ri::rle_string<sparse_bitvector_t, string_t>::load(in);
    }

protected:
    void build_rlbwt(std::ifstream &heads, std::ifstream &lengths, ulint B) {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        // assert(not contains0(input)); // We're hacking the 0 away :)
        this->B = B;
        // n = input.size();
        auto runs_per_letter_bv = vector<vector<bool>>(256);
        //runs in main bitvector
        vector<bool> runs_bv;

        // Reads the run heads
        string run_heads_s;
        heads.seekg(0, heads.end);
        run_heads_s.resize(heads.tellg());
        heads.seekg(0, heads.beg);
        heads.read(&run_heads_s[0], run_heads_s.size());

        size_t pos = 0;
        this->n = 0;
        this->R = run_heads_s.size();
        // Compute runs_bv and runs_per_letter_bv
        for (size_t i = 0; i < run_heads_s.size(); ++i)
        {
            size_t length;
            lengths.read((char *)&length, 5);
            if (run_heads_s[i] <= TERMINATOR) // change 0 to 1
                run_heads_s[i] = TERMINATOR;

            std::fill_n(std::back_inserter(runs_bv), length - 1, false);
            runs_bv.push_back(i % B == B - 1);

            std::fill_n(std::back_inserter(runs_per_letter_bv[run_heads_s[i]]), length - 1, false);
            runs_per_letter_bv[run_heads_s[i]].push_back(true);

            this->n += length;
        }
        // runs_bv.push_back(false);

        //now compact structures
        assert(runs_bv.size() == this->n);
        ulint t = 0;
        for (ulint i = 0; i < 256; ++i)
            t += runs_per_letter_bv[i].size();
        assert(t == this->n);
        this->runs = sparse_bitvector_t(runs_bv);
        //a fast direct array: char -> bitvector.
        this->runs_per_letter = vector<sparse_bitvector_t>(256);
        for (ulint i = 0; i < 256; ++i)
            this->runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);
        this->run_heads = string_t(run_heads_s);
        assert(this->run_heads.size() == this->R);
    }
private:
};

// Construction from run-length encoded BWT specialization for sparse_sd_vector
template <>
ms_rle_string<ri::sparse_sd_vector, ri::huff_string>::ms_rle_string(std::ifstream &heads, std::ifstream &lengths, ulint B)
{ 
    heads.clear();
    heads.seekg(0);
    lengths.clear();
    lengths.seekg(0);
    this->B = B;

    // Reads the run heads
    string run_heads_s;
    heads.seekg(0, heads.end);
    run_heads_s.resize(heads.tellg());
    heads.seekg(0, heads.beg);
    heads.read(&run_heads_s[0], run_heads_s.size());

    size_t pos = 0;
    this->n = 0;
    this->R = run_heads_s.size();

    auto runs_per_letter_bv = std::vector<std::vector<size_t>> (256);
    auto runs_per_letter_bv_i = vector<size_t> (256,0);

    // Runs in main bitvector
    vector<size_t> runs_bv_onset;
    size_t runs_bv_i = 0;

    // Compute runs_bv and runs_per_letter_bv
    for (size_t i = 0; i < run_heads_s.size(); ++i) {
        size_t length = 0;
        lengths.read((char *)&length, 5);

        uint8_t curr_ch = unsigned(run_heads_s[i]);
        if (curr_ch <= TERMINATOR) { // change 0 to 1
            run_heads_s[i] = TERMINATOR;
            curr_ch = TERMINATOR;
        }

        if(i % B == B - 1)
            runs_bv_onset.push_back(this->n + length - 1);

        assert(length > 0);
        runs_per_letter_bv_i[curr_ch] += length;
        runs_per_letter_bv[curr_ch].push_back(runs_per_letter_bv_i[curr_ch] - 1);
        this->n += length;

    }
    // Now compact structures
    ulint t = 0;
    for (ulint i = 0; i < 256; ++i) {
        t += runs_per_letter_bv_i[i];
    }
    assert(t == this->n);
    
    this->runs = ri::sparse_sd_vector(runs_bv_onset, this->n);  
    this->runs_per_letter = std::vector<ri::sparse_sd_vector>(256);

    for (ulint i = 0; i < 256; ++i) {
        this->runs_per_letter[i] = ri::sparse_sd_vector(runs_per_letter_bv[i], runs_per_letter_bv_i[i]);
    }

    this->run_heads = ri::huff_string(run_heads_s);
    assert(this->run_heads.size() == this->R);

    // debugging:
    // for (size_t i = 0; i < 256; i++) {
    //     std::cout << "i = " << i << ", runs_per_letter_size = " << this->runs_per_letter[i].size() << std::endl;
    // }
    // for (size_t i = 0; i < 256; i++) {
    //     std::cout << "i = " << i << ", huffstring.rank() = " << this->run_heads.rank(this->run_heads.size(), i) << std::endl;
    // }
};

typedef ms_rle_string<ri::sparse_sd_vector> ms_rle_string_sd;
typedef ms_rle_string<ri::sparse_hyb_vector> ms_rle_string_hyb;
#endif /* end of include guard: _MS_RLE_STRING_HH */
