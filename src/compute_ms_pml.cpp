/* 
 *  compute_ms_pml.cpp
 *  Copyright (C) 2021 Massimiliano Rossi
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/ .
 */

 /*
  * File: compute_ms_pml.cpp
  * Description: Main file for running SPUMONI to compute MS/PMLs
  *              against a reference. This file is combination of 
  *              previous source files (matching_statistics.cpp &
  *              run_spumoni.cpp)
  *
  * Authors: Massimiliano Rossi, Omar Ahmed
  * Start Date: October 13, 2021
  */

  extern "C" {
#include <xerrors.h>
}

#include <spumoni_main.hpp>
#include <r_index.hpp>
#include <thresholds_ds.hpp>
#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <compute_ms_pml.hpp>

/*
 * This first section of the code contains classes that define pml_pointers
 * and ms_pointers which are objects that basically the r-index plus the thresholds
 * it allows you to compute MS or PMLs.
 */

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class thresholds_t = thr_bv<rle_string_t> >
class pml_pointers : ri::r_index<sparse_bv_type, rle_string_t> {
  public:
    thresholds_t thresholds;
    typedef size_t size_type;

    pml_pointers() {}
    pml_pointers(std::string filename, bool rle = false) : ri::r_index<sparse_bv_type, rle_string_t>() {
        verbose("Building the r-index from BWT");
        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
        
        std::string bwt_fname = filename + ".bwt";
        verbose("RLE encoding BWT and computing SA samples");

        if (rle) {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        }
        else {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }

        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Text length: n = ", n);
        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("RL-BWT construction complete");
        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();
        thresholds = thresholds_t(filename,&this->bwt);
        t_insert_end = std::chrono::high_resolution_clock::now();

        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples) {
        int log_n = bitsize(uint64_t(n));
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");
        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2 * SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1, fd))
        {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths) {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
                this->F[c] += length;
            else
            {
                this->F[TERMINATOR] += length;
                this->terminator_position = i;
            }
            i++;
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }

    std::vector<size_t> query(const std::vector<uint8_t> &pattern) {
        /* Computes the PMLs for the given pattern */
        size_t m = pattern.size();
        return _query(pattern.data(), m);
    }

    std::vector<size_t> query(const char* pattern, const size_t m) {
        return _query(pattern, m);
    }

    void print_stats(){
        sdsl::nullstream ns;
        verbose("Memory consumption (bytes).");
        verbose("   terminator_position: ", sizeof(this->terminator_position));
        verbose("                     F: ", my_serialize(this->F, ns));
        verbose("                   bwt: ", this->bwt.serialize(ns));
        verbose("            thresholds: ", thresholds.serialize(ns));
    }

    std::pair<ulint, ulint> get_bwt_stats() {
        return std::make_pair(this->bwt_size() , this->bwt.number_of_runs());
    }

    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);

        written_bytes += thresholds.serialize(out, child, "thresholds");


        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    std::string get_file_extension() const {
        return thresholds.get_file_extension() + ".spumoni";
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in) {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);

        this->r = this->bwt.number_of_runs();
        thresholds.load(in,&this->bwt);
    }


protected:
    template<typename string_t>
    std::vector<size_t> _query(const string_t &pattern, const size_t m) {
        /* Actual PML Computation Method */
        std::vector<size_t> lengths(m);

        // Start with the empty string
        auto pos = this->bwt_size() - 1;
        auto length = 0;

        for (size_t i = 0; i < m; ++i) {
            auto c = pattern[m - i - 1];

            if (this->bwt.number_of_letter(c) == 0){length = 0;}
            else if (pos < this->bwt.size() && this->bwt[pos] == c) {length++;}
            else {
                // Get threshold
                ri::ulint rnk = this->bwt.rank(pos, c);
                size_t thr = this->bwt.size() + 1;

                ulint next_pos = pos;

                // if (rnk < (this->F[c] - this->F[c-1]) // I can use F to compute it
                if (rnk < this->bwt.number_of_letter(c)) {

                    // j is the first position of the next run of c's
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    thr = thresholds[run_of_j]; // If it is the first run thr = 0
                    length = 0;
                    next_pos = j;
                }

                if (pos < thr) {
                    rnk--;
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    length = 0;
                    next_pos = j;
                }
                pos = next_pos;
            }

            lengths[m - i - 1] = length;

            // Perform one backward step
            pos = LF(pos, c);
        }
        return lengths;
    }
}; /* End of pml_pointers class */


template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class thresholds_t = thr_bv<rle_string_t> >
class ms_pointers : ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    thresholds_t thresholds;

    // std::vector<ulint> samples_start;
    int_vector<> samples_start;
    // int_vector<> samples_end;
    // std::vector<ulint> samples_last;

    // static const uchar TERMINATOR = 1;
    // bool sais = true;
    // 
    //  * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
    //  
    // //F column of the BWT (vector of 256 elements)
    // std::vector<ulint> F;
    // //L column of the BWT, run-length compressed
    // rle_string_t bwt;
    // ulint terminator_position = 0;
    // ulint r = 0; //number of BWT runs

    typedef size_t size_type;

    ms_pointers() {}

    ms_pointers(std::string filename, bool rle = false) : ri::r_index<sparse_bv_type, rle_string_t>() {
        verbose("Building the r-index from BWT");
        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";
        verbose("RLE encoding BWT and computing SA samples");

        if (rle)
        {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        }
        else
        {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }
        // std::string istring;
        // read_file(bwt_fname.c_str(), istring);
        // for(size_t i = 0; i < istring.size(); ++i)
        //     if(istring[i]==0)
        //         istring[i] = TERMINATOR;
        // this->bwt = rle_string_t(istring);

        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();

        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));

        // this->build_F(istring);
        // istring.clear();
        // istring.shrink_to_fit();

        read_samples(filename + ".ssa", this->r, n, samples_start);
        read_samples(filename + ".esa", this->r, n, this->samples_last);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("R-index construction complete");
        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();

        thresholds = thresholds_t(filename,&this->bwt);

        // std::string tmp_filename = filename + std::string(".thr_pos");

        // struct stat filestat;
        // FILE *fd;

        // if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
        //     error("open() file " + tmp_filename + " failed");

        // int fn = fileno(fd);
        // if (fstat(fn, &filestat) < 0)
        //     error("stat() file " + tmp_filename + " failed");

        // if (filestat.st_size % THRBYTES != 0)
        //     error("invilid file " + tmp_filename);

        // size_t length = filestat.st_size / THRBYTES;
        // thresholds.resize(length);

        // for (size_t i = 0; i < length; ++i)
        //     if ((fread(&thresholds[i], THRBYTES, 1, fd)) != 1)
        //         error("fread() file " + tmp_filename + " failed");

        // fclose(fd);

        t_insert_end = std::chrono::high_resolution_clock::now();

        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples) {
        int log_n = bitsize(uint64_t(n));

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2 * SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1, fd))
        {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths) {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
                this->F[c] += length;
            else
            {
                this->F[TERMINATOR] += length;
                this->terminator_position = i;
            }
            i++;
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }

    std::vector<size_t> query(const std::vector<uint8_t> &pattern) {
        /* Computes the MS and Pointers for given pattern */
        size_t m = pattern.size();
        return _query(pattern.data(), m);
    }

    std::vector<size_t> query(const char* pattern, const size_t m) {
        return _query(pattern, m);
    }

    std::pair<ulint, ulint> get_bwt_stats() {
        /*
        for (int i = 0; i < length; i++){
            std::cout << "BWT[i] = " << this->bwt[i] << std::endl;
        } */
        return std::make_pair(this->bwt_size() , this->bwt.number_of_runs());
    }

    void print_stats() {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("   terminator_position: ", sizeof(this->terminator_position));
        verbose("                     F: ", my_serialize(this->F, ns));
        verbose("                   bwt: ", this->bwt.serialize(ns));
        verbose("          samples_last: ", this->samples_last.serialize(ns));
        verbose("            thresholds: ", thresholds.serialize(ns));
        verbose("         samples_start: ", samples_start.serialize(ns));
    }

    //
     //\param i position in the BWT
     //\param c character
     // \return lexicographic rank of cw in bwt
     //
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

      // serialize the structure to the ostream
     // \param out     the ostream
     //
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);
        written_bytes += this->samples_last.serialize(out);

        written_bytes += thresholds.serialize(out, child, "thresholds");
        // written_bytes += my_serialize(thresholds, out, child, "thresholds");
        // written_bytes += my_serialize(samples_start, out, child, "samples_start");
        written_bytes += samples_start.serialize(out, child, "samples_start");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    std::string get_file_extension() const {
        return thresholds.get_file_extension() + ".ms";
    }

    // load the structure from the istream
     // \param in the istream
     //
    void load(std::istream &in) {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->samples_last.load(in);

        thresholds.load(in,&this->bwt);
        // my_load(thresholds, in);
        samples_start.load(in);
        // my_load(samples_start,in);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

protected:
    template<typename string_t>
    std::vector<size_t> _query(const string_t &pattern, const size_t m) 
    {
        /* Actual MS Computation - notice the extra loop on the inside */
        std::vector<size_t> ms_pointers(m);

        // Start with the empty string
        auto pos = this->bwt_size() - 1;
        auto sample = this->get_last_run_sample();

        std::cout << "pos:" << pos << " bwt_char :" << this->bwt[pos] << std::endl;

        for (size_t i = 0; i < m; ++i) 
        {
            auto c = pattern[m - i - 1];

            int pos_int = pos;
            std::cout << "pos:" << pos << " bwt_char :" << this->bwt[pos] << std::endl;
            std::cout << "pos_int:" << pos_int << " bwt_char :" << this->bwt[pos] << std::endl;

            if (this->bwt.number_of_letter(c) == 0){sample = 0;}
            else if (pos < this->bwt.size() && this->bwt[pos] == c){sample--;}
            else {
                // Get threshold
                ri::ulint rnk = this->bwt.rank(pos, c);
                size_t thr = this->bwt.size() + 1;

                ulint next_pos = pos;

                // if (rnk < (this->F[c] - this->F[c-1]) // I can use F to compute it
                if (rnk < this->bwt.number_of_letter(c)) {

                    // j is the first position of the next run of c's
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    thr = thresholds[run_of_j]; // If it is the first run thr = 0

                    // Here we should use Phi_inv that is not implemented yet
                    // sample = this->Phi(this->samples_last[run_of_j - 1]) - 1;
                    sample = samples_start[run_of_j];

                    next_pos = j;
                }

                if (pos < thr) {
                    rnk--;
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    sample = this->samples_last[run_of_j];
                    next_pos = j;
                }

                pos = next_pos;
            }

            ms_pointers[m - i - 1] = sample;
            
            // Perform one backward step
            pos = LF(pos, c);
        }

        return ms_pointers;
    }
    // // From r-index
    // vector<ulint> build_F(std::ifstream &ifs)
    // {
    //     ifs.clear();
    //     ifs.seekg(0);
    //     F = vector<ulint>(256, 0);
    //     uchar c;
    //     ulint i = 0;
    //     while (ifs >> c)
    //     {
    //         if (c > TERMINATOR)
    //             F[c]++;
    //         else
    //         {
    //             F[TERMINATOR]++;
    //             terminator_position = i;
    //         }
    //         i++;
    //     }
    //     for (ulint i = 255; i > 0; --i)
    //         F[i] = F[i - 1];
    //     F[0] = 0;
    //     for (ulint i = 1; i < 256; ++i)
    //         F[i] += F[i - 1];
    //     return F;
    // }

    // // From r-index
    // vector<pair<ulint, ulint>> &read_run_starts(std::string fname, ulint n, vector<pair<ulint, ulint>> &ssa)
    // {
    //     ssa.clear();
    //     std::ifstream ifs(fname);
    //     uint64_t x = 0;
    //     uint64_t y = 0;
    //     uint64_t i = 0;
    //     while (ifs.read((char *)&x, 5) && ifs.read((char *)&y, 5))
    //     {
    //         ssa.push_back(pair<ulint, ulint>(y ? y - 1 : n - 1, i));
    //         i++;
    //     }
    //     return ssa;
    // }

    // // From r-index
    // vector<ulint> &read_run_ends(std::string fname, ulint n, vector<ulint> &esa)
    // {
    //     esa.clear();
    //     std::ifstream ifs(fname);
    //     uint64_t x = 0;
    //     uint64_t y = 0;
    //     while (ifs.read((char *)&x, 5) && ifs.read((char *)&y, 5))
    //     {
    //         esa.push_back(y ? y - 1 : n - 1);
    //     }
    //     return esa;
    // }
};


/*
 * The next section contains another set of classes that are instantiated 
 * when loading the MS or PML index for computation, and they are called
 * ms_t and pml_t. One of its attributes the ms_pointers and pml_pointers as
 * previously defined.
 */

class pml_t {
public:
    pml_t(std::string filename){
      SPUMONI_LOG("Loading the PML index ...");
      auto start_time = std::chrono::system_clock::now();
      std::string filename_ms = filename + ms.get_file_extension();

      ifstream fs_ms(filename_ms);
      ms.load(fs_ms);
      fs_ms.close();

      auto end_time = std::chrono::system_clock::now();
      //SPUMONI_LOG("Memory Peak: %d", malloc_count_peak());
      TIME_LOG((end_time - start_time));
    }

  // Destructor
  ~pml_t() {}

  // The outfile has the following format. The first size_t integer store the
  // length l of the query. Then the following l size_t integers stores the
  // pointers of the matching statistics, and the following l size_t integers
  // stores the lengths of the mathcing statistics.
  //void matching_statistics(kseq_t *read, FILE* out) 
  void matching_statistics(const char* read, size_t read_length, FILE* out) 
  {
    auto lengths = ms.query(read, read_length);

    size_t q_length = lengths.size();
    fwrite(&q_length, sizeof(size_t), 1,out);
    fwrite(lengths.data(), sizeof(size_t),q_length,out);
  }

  std::pair<ulint, ulint> get_bwt_stats() {
      return ms.get_bwt_stats();
  }

protected:
  pml_pointers<> ms;
  size_t n = 0;
};

class ms_t {
public:
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;
  // using SlpT = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;

  ms_t(std::string filename) 
  {
    SPUMONI_LOG("Loading the MS index ...");
    auto start_time = std::chrono::system_clock::now();    
    std::string filename_ms = filename + ms.get_file_extension();

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);
    fs_ms.close();

    auto end_time = std::chrono::system_clock::now();
    //SPUMONI_LOG("Memory Peak: %d", malloc_count_peak());
    TIME_LOG((end_time - start_time));

    SPUMONI_LOG("Loading the Random Access Data Structure ...");
    start_time = std::chrono::system_clock::now();   
    std::string filename_slp = filename + ".slp";

    ifstream fs(filename_slp);
    ra.load(fs);
    fs.close();

    n = ra.getLen();
    end_time = std::chrono::system_clock::now();
    //SPUMONI_LOG("Memory Peak: %d", malloc_count_peak());
    TIME_LOG((end_time - start_time));
  }

  // Destructor
  ~ms_t() {}

  // The outfile has the following format. The first size_t integer store the
  // length l of the query. Then the following l size_t integers stores the
  // pointers of the matching statistics, and the following l size_t integers
  // stores the lengths of the mathcing statistics.
  void matching_statistics(const char* read, size_t read_length, FILE* out)
  {
    auto pointers = ms.query(read, read_length);
    std::vector<size_t> lengths(pointers.size());
    size_t l = 0;

    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];

      // The (i < 1 || pos != (pointers[i-1] + 1) ) term was added from an earlier version
      while ((i + l) < read_length && (pos + l) < n && (i < 1 || pos != (pointers[i-1] + 1) ) && read[i + l] == ra.charAt(pos + l))
        ++l;

      lengths[i] = l;
      l = (l == 0 ? 0 : (l - 1));
    }
    
    assert(lengths.size() == pointers.size());

    size_t q_length = pointers.size();
    fwrite(&q_length, sizeof(size_t), 1,out);
    fwrite(pointers.data(), sizeof(size_t),q_length,out);
    fwrite(lengths.data(), sizeof(size_t),q_length,out);
  }

  std::pair<ulint, ulint> get_bwt_stats() {
      return ms.get_bwt_stats();
  }
  
protected:
  ms_pointers<> ms;
  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;
  size_t n = 0;
};

/*
 * This section of the code contains some helper methods, struct defintions, and
 * wrapper methods for the kseq.h.
 */

char complement(char n)
{
  switch (n) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';
    default: return n;
  }
}

typedef struct{
  // For MS multi-threading
  ms_t *ms;
  std::string pattern_filename;
  std::string out_filename;
  size_t start;
  size_t end;
  size_t wk_id;
} mt_ms_param;

typedef struct{
  // For PML multi-threading
  pml_t *ms;
  std::string pattern_filename;
  std::string out_filename;
  size_t start;
  size_t end;
  size_t wk_id;
} mt_pml_param;


static inline size_t ks_tell(kseq_t *seq) {
  return gztell(seq->f->f) - seq->f->end + seq->f->begin;
}

static void copy_kstring_t(kstring_t &l, kstring_t &r) {
  l.l = r.l;
  l.m = r.m;
  l.s = (char *)malloc(l.m);
  for (size_t i = 0; i < r.m; ++i)
    l.s[i] = r.s[i];
}

static void copy_kseq_t(kseq_t *l, kseq_t *r) {
  copy_kstring_t(l->name, r->name);
  copy_kstring_t(l->comment, r->comment);
  copy_kstring_t(l->seq, r->seq);
  copy_kstring_t(l->qual, r->qual);
  l->last_char = r->last_char;
}

static size_t next_start_fastq(gzFile fp){
    int c;
    // Special case when we arr at the beginning of the file.
    if ((gztell(fp) == 0) && ((c = gzgetc(fp)) != EOF) && c == '@')
        return 0;

    // Strart from the previous character
    gzseek(fp, -1, SEEK_CUR);

    std::vector<std::pair<int, size_t>> window;
    // Find the first new line
    for (size_t i = 0; i < 4; ++i)
    {
        while (((c = gzgetc(fp)) != EOF) && (c != (int)'\n'))
        {
        }
        if (c == EOF)
        return gztell(fp);
        if ((c = gzgetc(fp)) == EOF)
        return gztell(fp);
        window.push_back(std::make_pair(c, gztell(fp) - 1));
    }

    for (size_t i = 0; i < 2; ++i)
    {
        if (window[i].first == '@' && window[i + 2].first == '+')
        return window[i].second;
        if (window[i].first == '+' && window[i + 2].first == '@')
        return window[i + 2].second;
  }

  return gztell(fp);
}

static inline bool is_gzipped(std::string filename) {
    /* Methods tests if the file is gzipped */
    FILE *fp = fopen(filename.c_str(), "rb");
    if(fp == NULL) error("Opening file " + filename);
    int byte1 = 0, byte2 = 0;
    fread(&byte1, sizeof(char), 1, fp);
    fread(&byte2, sizeof(char), 1, fp);
    fclose(fp);
    return (byte1 == 0x1f && byte2 == 0x8b);
}

inline size_t get_file_size(std::string filename) {
    /* Returns the length of the file, it is assuming the file is not compressed */
    if (is_gzipped(filename))
    {
        std::cerr << "The input is gzipped!" << std::endl;
        return -1;
    }
    FILE *fp = fopen(filename.c_str(), "r");
    fseek(fp, 0L, SEEK_END);
    size_t size = ftell(fp);
    fclose(fp);
    return size;
}

 std::vector<size_t> split_fastq(std::string filename, size_t n_threads)  {
    /*
    * Precondition: the file is not gzipped
    * scan file for start positions and execute threads
    */ 
    size_t size = get_file_size(filename);
    gzFile fp = gzopen(filename.c_str(), "r");

    if (fp == Z_NULL) {
        throw new std::runtime_error("Cannot open input file " + filename);
    }

    std::vector<size_t> starts(n_threads + 1);
    for (int i = 0; i < n_threads + 1; ++i)
    {
        size_t start = (size_t)((size * i) / n_threads);
        gzseek(fp, start, SEEK_SET);
        starts[i] = next_start_fastq(fp);
    }
    gzclose(fp);
    return starts;
}

/*
 * This section of the code deals with multi-threading of processing of 
 * pattern reads using pthreads.
 *
 * TODO: Try to replace with OpenMP threading, and using a producer/consumer-model
 *       and see if it will help to simplify the code.
 */

void *mt_pml_worker(void *param) {
    mt_pml_param *p = (mt_pml_param*) param;
    size_t n_reads = 0;
    size_t n_aligned_reads = 0;

    FILE *out_fd;
    gzFile fp;

    if ((out_fd = fopen(p->out_filename.c_str(), "w")) == nullptr)
        error("open() file " + p->out_filename + " failed");

    if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
        error("open() file " + p->pattern_filename + " failed");

    gzseek(fp, p->start, SEEK_SET);

    kseq_t rev;
    int l;

    kseq_t *seq = kseq_init(fp);
    while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0)) {
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); //Make sure all characters are upper-case

        p->ms->matching_statistics(curr_read.c_str(), seq->seq.l, out_fd);
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_fd);

    return NULL;
}

void mt_pml(pml_t *ms, std::string pattern_filename, std::string out_filename, size_t n_threads) {
    pthread_t t[n_threads] = {0};
    mt_pml_param params[n_threads];
    std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
    for(size_t i = 0; i < n_threads; ++i)
    {
        params[i].ms = ms;
        params[i].pattern_filename = pattern_filename;
        params[i].out_filename = out_filename + "_" + std::to_string(i) + ".ms.tmp.out";
        params[i].start = starts[i];
        params[i].end = starts[i+1];
        params[i].wk_id = i;
        xpthread_create(&t[i], NULL, &mt_pml_worker, &params[i], __LINE__, __FILE__);
    }

    for(size_t i = 0; i < n_threads; ++i)
    {
        xpthread_join(t[i],NULL,__LINE__,__FILE__);
    }
    return;
}

void *mt_ms_worker(void *param) {
    mt_ms_param *p = (mt_ms_param*) param;
    size_t n_reads = 0;
    size_t n_aligned_reads = 0;

    FILE *out_fd;
    gzFile fp;

    if ((out_fd = fopen(p->out_filename.c_str(), "w")) == nullptr)
        error("open() file " + p->out_filename + " failed");

    if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
        error("open() file " + p->pattern_filename + " failed");

    gzseek(fp, p->start, SEEK_SET);

    kseq_t rev;
    int l;

    kseq_t *seq = kseq_init(fp);
    while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0))
    {
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); //Make sure all characters are upper-case

        p->ms->matching_statistics(curr_read.c_str(), seq->seq.l, out_fd);

    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_fd);

    return NULL;
}

void mt_ms(ms_t *ms, std::string pattern_filename, std::string out_filename, size_t n_threads) {
    pthread_t t[n_threads] = {0};
    mt_ms_param params[n_threads];
    std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
    for(size_t i = 0; i < n_threads; ++i)
    {
        params[i].ms = ms;
        params[i].pattern_filename = pattern_filename;
        params[i].out_filename = out_filename + "_" + std::to_string(i) + ".ms.tmp.out";
        params[i].start = starts[i];
        params[i].end = starts[i+1];
        params[i].wk_id = i;
        xpthread_create(&t[i], NULL, &mt_ms_worker, &params[i], __LINE__, __FILE__);
    }

    for(size_t i = 0; i < n_threads; ++i)
    {
        xpthread_join(t[i],NULL,__LINE__,__FILE__);
    }
    return;
}

/*
 * This section of the code contains the single-threaded processing methods
 * for computing the MS/PMLs.
 *
 * TODO: As mentioned above, I would like to use OpenMP since that could
 *       simplify the code since we can use pragmas instead of writing separate
 *       code.
 */

size_t st_pml(pml_t *ms, std::string pattern_filename, std::string out_filename) {
    size_t n_reads = 0;
    size_t n_aligned_reads = 0;
    kseq_t rev;
    int l;
    FILE *out_fd;

    out_filename += "_0.ms.tmp.out";

    if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
        error("open() file " + out_filename + " failed");

    gzFile fp = gzopen(pattern_filename.c_str(), "r");
    kseq_t* seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); //Make sure all characters are upper-case

        ms->matching_statistics(curr_read.c_str(), seq->seq.l, out_fd);
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_fd);

    return n_aligned_reads;
}

size_t st_pml_general(pml_t *ms, std::string pattern_filename, std::string out_filename) {
    // Check if file is mistakenly labeled as non-fasta file
    std::string file_ext = pattern_filename.substr(pattern_filename.find_last_of(".") + 1);
    if (file_ext == "fa" || file_ext == "fasta") {
        error("The file extension for the patterns suggests it is a fasta file. Please run with -f option for correct results.");
    }

    size_t n_reads = 0;
    size_t n_aligned_reads = 0;
    kseq_t rev;
    int l;
    FILE *out_fd;
    out_filename += "_0.ms.tmp.out";

    if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
        error("open() file " + out_filename + " failed");

    std::ifstream input_fd (pattern_filename, std::ifstream::in | std::ifstream::binary);
    char ch = input_fd.get();
    char* buf = new char [2]; // To just hold separator character
    std::string read = "";

    while (input_fd.good()) {
        if (ch == '\x01') {ms->matching_statistics(read.c_str(), read.size(), out_fd); input_fd.read(buf, 2); read = "";}
        else {read += ch;}
        ch = input_fd.get();
    }
    fclose(out_fd);

    return n_aligned_reads;
}


size_t st_ms(ms_t *ms, std::string pattern_filename, std::string out_filename) {
    size_t n_reads = 0;
    size_t n_aligned_reads = 0;
    kseq_t rev;
    int l;
    FILE *out_fd;

    out_filename += "_0.ms.tmp.out";

    if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
        error("open() file " + out_filename + " failed");

    gzFile fp = gzopen(pattern_filename.c_str(), "r");
    kseq_t* seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); //Make sure all characters are upper-case

        ms->matching_statistics(curr_read.c_str(), seq->seq.l, out_fd);
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_fd);

    return n_aligned_reads;
}

size_t st_ms_general(ms_t *ms, std::string pattern_filename, std::string out_filename)
{
    // Check if file is mistakenly labeled as non-fasta file
    std::string file_ext = pattern_filename.substr(pattern_filename.find_last_of(".") + 1);
    if (file_ext == "fa" || file_ext == "fasta") {
        FATAL_WARNING("The file extension for the patterns suggests it is a fasta file.\n 
                      Please run with -f option for correct results.");
    }

    size_t n_reads = 0;
    size_t n_aligned_reads = 0;
    kseq_t rev;
    int l;
    FILE *out_fd;
    out_filename += "_0.ms.tmp.out";

    if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
        error("open() file " + out_filename + " failed");

    std::ifstream input_fd (pattern_filename, std::ifstream::in | std::ifstream::binary);
    char ch = input_fd.get();
    char* buf = new char [2]; // To just hold separator character
    std::string read = "";

    while (input_fd.good()) {
        if (ch == '\x01') {ms->matching_statistics(read.c_str(), read.size(), out_fd); input_fd.read(buf, 2); read="";}
        else {read += ch;}
        ch = input_fd.get();
    }
    fclose(out_fd);

    return n_aligned_reads;
}

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;


/*
 * This section contains the "main" methods for the running process where
 * the first method computes the PMLs, and the second one computes the MSs
 * given the index and pattern.
 */

int run_spumoni_main(SpumoniRunOptions* run_opts){
  /* This method is responsible for the PML computation */

  /* Loads the RLEBWT and Thresholds*/
  pml_t ms(run_opts->ref_file);

  /* Load patterns and generate PMLs */
  std::string base_name = basename(run_opts->ref_file.data());
  std::string out_filename = run_opts->pattern_file + "_ref_" + base_name;

  if (is_gzipped(run_opts->pattern_file)) { 
    SPUMONI_LOG("The input is gzipped - forcing single threaded pseudo matching lengths.");
    run_opts->threads = 1;
  }

  /* Determine approach to parse pattern files (based on whether it is general text or fasta format) */
  auto start_time = std::chrono::system_clock::now();
  SPUMONI_LOG("Starting processing the patterns ...");

  if (run_opts->query_fasta) {
    if(run_opts->threads <= 1) {st_pml(&ms, run_opts->pattern_file, out_filename);}
    else {mt_pml(&ms, run_opts->pattern_file, out_filename, run_opts->threads);}
  }
  else {
        if(run_opts->threads == 1) {st_pml_general(&ms, run_opts->pattern_file,out_filename);}
    else {FATAL_WARNING("Multi-threading not implemented yet for general-text querying.");}
  }

  auto end_time = std::chrono::system_clock::now();
  //SPUMONI_LOG("Memory peak: %d", malloc_count_peak());
  TIME_LOG((end_time - start_time));

  /* Writing out the PMLs */
  SPUMONI_LOG("Writing the plain output ...");
  start_time = std::chrono::system_clock::now();
  
  std::ofstream f_lengths(out_filename + ".pseudo_lengths");
  if (!f_lengths.is_open()) {error("open() file " + std::string(out_filename) + ".lengths failed");}

  size_t n_seq = 0;
  for(size_t i = 0; i < run_opts->threads; ++i) {
    std::string tmp_filename = out_filename + "_" + std::to_string(i) + ".ms.tmp.out";
    FILE *in_fd;

    if ((in_fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
      error("open() file " + tmp_filename + " failed");

    size_t length = 0;
    size_t m = 100; // Reserved size for pointers and lengths
    size_t *mem = (size_t*) malloc(m * sizeof(size_t));
    while(!feof(in_fd) and fread(&length,sizeof(size_t), 1,in_fd) > 0)
    {
      if( m < length) {
        // Resize lengths and pointers
        m = length;
        mem = (size_t*) realloc(mem, m * sizeof(size_t));
      }
        
      if ((fread(mem, sizeof(size_t), length, in_fd)) != length)
        error("fread() file " + std::string(tmp_filename) + " failed");

      f_lengths << ">" + std::to_string(n_seq) << endl;
      for(size_t i = 0; i < length; ++i)
        f_lengths  << mem[i] << " ";
      f_lengths << endl;

      n_seq ++;
    }
    fclose(in_fd);
  }
  f_lengths.close();

  end_time = std::chrono::system_clock::now();
  //SPUMONI_LOG("Memory peak: %d", malloc_count_peak());
  TIME_LOG((end_time - start_time));
  return 0;
}

int run_spumoni_ms_main(SpumoniRunOptions* run_opts) {
  /* This method is responsible for the MS computation */
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;
  
  // Loads the MS index containing the RLEBWT, Thresholds, and RA structure
  ms_t ms(run_opts->ref_file);

  // Load patterns and generate MSs
  std::string base_name = basename(run_opts->ref_file.data());
  std::string out_filename = run_opts->pattern_file + "_ref_" + base_name;

  if (is_gzipped(run_opts->pattern_file)) {
    SPUMONI_LOG("The input is gzipped - forcing single threaded matching statistics.");
    run_opts->threads = 1;
  }

  // Determine approach to parse pattern files
  auto start_time = std::chrono::system_clock::now();
  SPUMONI_LOG("Starting processing the patterns ...");

  if (run_opts->query_fasta) {
    if(run_opts->threads == 1) {st_ms(&ms, run_opts->pattern_file, out_filename);}
    else {FATAL_WARNING("Multi-threading not implemented yet for FASTA querying.");
          mt_ms(&ms, run_opts->pattern_file, out_filename, run_opts->threads);}
  }
  else {
        if(run_opts->threads == 1) {st_ms_general(&ms,run_opts->pattern_file,out_filename);}
    else {FATAL_WARNING("Multi-threading not implemented yet for general-text querying.");}
  }

  auto end_time = std::chrono::system_clock::now();
  TIME_LOG((end_time - start_time));

  /* Writing out the MSs and pointers */
  SPUMONI_LOG("Writing the plain output ...");
  start_time = std::chrono::system_clock::now();
  
  std::ofstream f_pointers(out_filename + ".pointers");
  std::ofstream f_lengths(out_filename + ".lengths");

  if (!f_pointers.is_open() || !f_lengths.is_open())
    error("Opening the .pointers or .lengths file failed");



  size_t n_seq = 0;
  for(size_t i = 0; i < run_opts->threads; ++i)
  {
    std::string tmp_filename = out_filename + "_" + std::to_string(i) + ".ms.tmp.out";
    FILE *in_fd;

    if ((in_fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
      error("open() file " + tmp_filename + " failed");

    size_t length = 0;
    size_t m = 100; // Reserved size for pointers and lengths
    size_t *mem = (size_t*) malloc(m * sizeof(size_t));
    while(!feof(in_fd) and fread(&length,sizeof(size_t), 1,in_fd) > 0)
    {
      if( m < length)
      {
        // Resize lengths and pointers
        m = length;
        mem = (size_t*) realloc(mem, m * sizeof(size_t));
      }

      if ((fread(mem, sizeof(size_t), length, in_fd)) != length)
        error("fread() file " + std::string(tmp_filename) + " failed");
      
      // TODO: Store the fasta headers somewhere
      f_pointers << ">" + std::to_string(n_seq) << endl;
      for(size_t i = 0; i < length; ++i)
        f_pointers << mem[i] << " ";
      f_pointers << endl;
      
        
      if ((fread(mem, sizeof(size_t), length, in_fd)) != length)
        error("fread() file " + std::string(tmp_filename) + " failed");

      f_lengths << ">" + std::to_string(n_seq) << endl;
      for(size_t i = 0; i < length; ++i)
        f_lengths  << mem[i] << " ";
      f_lengths << endl;

      n_seq ++;
    }
    fclose(in_fd);
  }

  f_pointers.close();
  f_lengths.close();

  end_time = std::chrono::system_clock::now();
  //SPUMONI_LOG("Memory peak: %d", malloc_count_peak());
  TIME_LOG((end_time - start_time));
  return 0;
}

std::pair<ulint, ulint> get_bwt_stats(std::string ref_file, size_t type) {
    /* Returns the length and number of runs in a text */
    ulint length = 0, num_runs = 0;
    if (type == 1) {
        ms_t ms_data_structure(ref_file);
        std::tie(length, num_runs) = ms_data_structure.get_bwt_stats();
        return std::make_pair(length, num_runs);
    } else {
        pml_t pml_data_structure(ref_file);
        std::tie(length, num_runs) = pml_data_structure.get_bwt_stats();
        return std::make_pair(length, num_runs);
    }
}