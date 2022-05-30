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
 
#include <spumoni_main.hpp>
#include <r_index.hpp>
#include <thresholds_ds.hpp>
#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <compute_ms_pml.hpp>
#include <tuple>
#include <fstream>
#include <iterator>
#include <doc_array.hpp>
#include <bits/stdc++.h>
#include <ks_test.hpp>
#include <omp.h>
#include <batch_loader.hpp>

/*
 * This first section of the code contains classes that define pml_pointers
 * and ms_pointers which are objects that basically the r-index plus the 
 * thresholds which allow you to compute MS or PMLs.
 */

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class thresholds_t = thr_bv<rle_string_t> >
class pml_pointers : ri::r_index<sparse_bv_type, rle_string_t> {
public:
    thresholds_t thresholds;
    typedef size_t size_type;
    size_t num_runs;

    pml_pointers() {}
    pml_pointers(std::string filename, bool rle = false) : ri::r_index<sparse_bv_type, rle_string_t>() {    
        std::string bwt_fname = filename + ".bwt";

        if (rle) {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);

            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);

            ifs_heads.close();
            ifs_len.close();
        } else {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }

        this->r = this->bwt.number_of_runs();
        this->num_runs = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        //SPUMONI_LOG("Text length: n = %d", n);
        //SPUMONI_LOG("Number of BWT equal-letter runs: r = %d", this->r);
        //SPUMONI_LOG("Rate n/r = %.4f", double(this->bwt.size()) / this->r);
        //SPUMONI_LOG("log2(r) = %.4f", log2(double(this->r)));
        //SPUMONI_LOG("log2(n/r) = %.4f", log2(double(this->bwt.size()) / this->r));

        thresholds = thresholds_t(filename,&this->bwt);
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

    /*
     * Overloaded functions - based on wheter you want to report the
     * document numbers or not.
     */
    void query(const char* pattern, const size_t m, std::vector<size_t>& lengths) {
        _query(pattern, m, lengths);
    }

    void query(const char* pattern, const size_t m, std::vector<size_t>& lengths,
                              std::vector<size_t>& doc_nums, DocumentArray& doc_arr) {
        _query(pattern, m, lengths, doc_nums, doc_arr);
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
        return std::make_pair(this->bwt_size(), this->bwt.number_of_runs());
    }

    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
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
    /*
     * Overloaded functions - based on whether you want to report the
     * document numbers or not.
     */
    template<typename string_t>
    void _query(const string_t &pattern, const size_t m, std::vector<size_t>& lengths) {
        // Actual PML computation method
        lengths.resize(m);

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
    }

    template<typename string_t>
    void _query(const string_t &pattern, const size_t m, std::vector<size_t>& lengths,
                std::vector<size_t>& doc_nums, DocumentArray& doc_arr) {
        // Actual PML computation method
        lengths.resize(m);
        doc_nums.resize(m);

        // Start with the empty string
        auto pos = this->bwt_size() - 1;
        auto length = 0;
        size_t curr_doc_id = doc_arr.end_runs_doc[this->bwt.number_of_runs()-1];

        for (size_t i = 0; i < m; ++i) {
            auto c = pattern[m - i - 1];

            if (this->bwt.number_of_letter(c) == 0){length = 0;}
            else if (pos < this->bwt.size() && this->bwt[pos] == c) {length++;}
            else {
                // Get threshold
                ri::ulint rnk = this->bwt.rank(pos, c);
                size_t thr = this->bwt.size() + 1;
                ulint next_pos = pos;

                if (rnk < this->bwt.number_of_letter(c)) {
                    // j is the first position of the next run of c's
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    thr = thresholds[run_of_j]; // If it is the first run thr = 0
                    curr_doc_id = doc_arr.start_runs_doc[run_of_j];

                    length = 0;
                    next_pos = j;
                }

                if (pos < thr) {
                    rnk--;
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);
                    curr_doc_id = doc_arr.end_runs_doc[run_of_j];

                    length = 0;
                    next_pos = j;
                }
                pos = next_pos;
            }

            lengths[m-i-1] = length;
            doc_nums[m-i-1] = curr_doc_id;
            // Perform one backward step
            pos = LF(pos, c);
        }
    }
}; /* End of pml_pointers class */


template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class thresholds_t = thr_bv<rle_string_t>>
class ms_pointers : ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    thresholds_t thresholds;
    int_vector<> samples_start;
    typedef size_t size_type;
    size_t num_runs;

    ms_pointers() {}

    ms_pointers(std::string filename, bool rle = false) : ri::r_index<sparse_bv_type, rle_string_t>() {
        std::string bwt_fname = filename + ".bwt";
            
        if (rle){
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);  

            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);

            ifs_heads.close();
            ifs_len.close();
        } else {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }
        
        this->r = this->bwt.number_of_runs();
        this->num_runs = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        //SPUMONI_LOG("Number of BWT equal-letter runs: r = %d", this->r);
        //SPUMONI_LOG("Rate n/r = %.4f", double(this->bwt.size()) / this->r);
        //SPUMONI_LOG("log2(r) = %.4f", log2(double(this->r)));
        //SPUMONI_LOG("log2(n/r) = %.4f", log2(double(this->bwt.size()) / this->r));

        // this->build_F(istring);
        // istring.clear();
        // istring.shrink_to_fit();

        read_samples(filename + ".ssa", this->r, n, samples_start);
        read_samples(filename + ".esa", this->r, n, this->samples_last);

        // Reading in the thresholds
        thresholds = thresholds_t(filename, &this->bwt);
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

    /*
     * Overloaded functions - based on whether you want to report the document 
     * numbers as well or not.
     */
    void query(const char* pattern, const size_t m, std::vector<size_t>& pointers) {
        _query(pattern, m, pointers);
    }

    void query(const char* pattern, const size_t m, std::vector<size_t>& pointers, std::vector<size_t>& doc_nums,
               DocumentArray& doc_array){
        _query(pattern, m, pointers, doc_nums, doc_array);
    } 

    std::pair<ulint, ulint> get_bwt_stats() {
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


protected:
    /*
     * Overloaded functions - based on whether you want to get the document 
     * numbers or not.
     */
    template<typename string_t>
    void _query(const string_t &pattern, const size_t m, std::vector<size_t>& ms_pointers) {
        // Actual MS Computation
        ms_pointers.resize(m);
        auto pos = this->bwt_size() - 1;
        auto sample = this->get_last_run_sample();

        for (size_t i = 0; i < m; ++i) 
        {
            auto c = pattern[m - i - 1];

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
    }

    template<typename string_t>
    void _query(const string_t &pattern, const size_t m, std::vector<size_t>& ms_pointers,
                std::vector<size_t>& doc_nums, DocumentArray& doc_arr) {
        // Actual MS Computation
        ms_pointers.resize(m);
        doc_nums.resize(m);

        auto pos = this->bwt_size() - 1;
        auto sample = this->get_last_run_sample();
        size_t curr_doc_id = doc_arr.end_runs_doc[this->bwt.number_of_runs()-1];

        for (size_t i = 0; i < m; ++i) 
        {
            auto c = pattern[m - i - 1];
            if (this->bwt.number_of_letter(c) == 0){
                sample = 0;
                ri::ulint run_of_j = this->bwt.run_of_position(sample);
                curr_doc_id = doc_arr.start_runs_doc[run_of_j];
            }
            else if (pos < this->bwt.size() && this->bwt[pos] == c){sample--;}
            else {
                // Get threshold
                ri::ulint rnk = this->bwt.rank(pos, c);
                size_t thr = this->bwt.size() + 1;
                ulint next_pos = pos;

                if (rnk < this->bwt.number_of_letter(c)) {
                    // j is the first position of the next run of c's
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    thr = thresholds[run_of_j]; // If it is the first run thr = 0
                    sample = samples_start[run_of_j];
                    curr_doc_id = doc_arr.start_runs_doc[run_of_j];

                    next_pos = j;
                }

                if (pos < thr) {
                    rnk--;
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    sample = this->samples_last[run_of_j];
                    curr_doc_id = doc_arr.end_runs_doc[run_of_j];
                    next_pos = j;
                }

                pos = next_pos;
            }

            ms_pointers[m-i-1] = sample;
            doc_nums[m-i-1] = curr_doc_id;
            
            // Perform one backward step
            pos = LF(pos, c);
        }
    }

}; /* End of ms_pointers */


/*
 * The next section contains another set of classes that are instantiated 
 * when loading the MS or PML index for computation, and they are called
 * ms_t and pml_t. One of its attributes are the ms_pointers and 
 * pml_pointers as previously defined.
 */

class pml_t {
public:
    using DocArray = DocumentArray;
    DocArray doc_arr; 

    // Constructor
    pml_t(std::string filename, bool use_doc, bool verbose = false){
        if (verbose){STATUS_LOG("pml_construct", "loading the PML index");}
        auto start_time = std::chrono::system_clock::now();
        std::string filename_ms = filename + ms.get_file_extension();

        std::ifstream fs_ms(filename_ms);
        ms.load(fs_ms);
        fs_ms.close();

        auto end_time = std::chrono::system_clock::now();
        if (verbose) {DONE_LOG((end_time - start_time));}

        if (use_doc) {
            if (verbose) {STATUS_LOG("pml_construct", "loading the document array");}
            start_time = std::chrono::system_clock::now();
            std::ifstream doc_file(filename + ".doc");

            doc_arr.load(doc_file);
            doc_file.close();
            if (verbose) {DONE_LOG((std::chrono::system_clock::now() - start_time));}
        }
    }

    //Destructor
    ~pml_t() {}

    /*
     * Overloaded functions - based on whether you want to report the
     * document numbers or not.
     */
    void matching_statistics(const char* read, size_t read_length, std::vector<size_t>& lengths) {
        ms.query(read, read_length, lengths);
    }

    void matching_statistics(const char* read, size_t read_length, std::vector<size_t>& lengths, 
                             std::vector<size_t>& doc_nums) {
        ms.query(read, read_length, lengths, doc_nums, doc_arr);
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
    using DocArray = DocumentArray;
    DocArray doc_arr;

    ms_t(std::string filename, bool use_doc, bool verbose=false) {
        if (verbose) {STATUS_LOG("ms_construct", "loading the MS index");}
        auto start_time = std::chrono::system_clock::now();    
        std::string filename_ms = filename + ms.get_file_extension();

        ifstream fs_ms(filename_ms);
        ms.load(fs_ms);
        fs_ms.close();

        auto end_time = std::chrono::system_clock::now();
        if (verbose) {DONE_LOG((end_time - start_time));}

        if (verbose) {STATUS_LOG("ms_construct", "loading the random access data structure");}
        start_time = std::chrono::system_clock::now();   
        std::string filename_slp = filename + ".slp";

        ifstream fs(filename_slp);
        ra.load(fs);
        fs.close();
        n = ra.getLen();
        if (verbose) {DONE_LOG((std::chrono::system_clock::now() - start_time));}

        if (use_doc) {
            if (verbose) {STATUS_LOG("ms_construct", "loading the document array");}
            start_time = std::chrono::system_clock::now();
            std::ifstream doc_file(filename + ".doc");

            doc_arr.load(doc_file);
            doc_file.close();
            if (verbose) {DONE_LOG((std::chrono::system_clock::now() - start_time));}
        }
    } 

    // Destructor
    ~ms_t() {}

    /*
     * Overloaded functions - used to compute the MS depending on 
     * whether you want to extract document numbers or not.
     */
    void matching_statistics(const char* read, size_t read_length, std::vector<size_t>& lengths, 
                            std::vector<size_t>& pointers) {  
        // Takes a read, and generates the MS with respect to this ms_t object
        ms.query(read, read_length, pointers);
        lengths.resize(read_length);
        size_t l = 0;

        for (size_t i = 0; i < pointers.size(); ++i) {
            size_t pos = pointers[i];
            while ((i + l) < read_length && (pos + l) < n && (i < 1 || pos != (pointers[i-1] + 1) ) && read[i + l] == ra.charAt(pos + l))
                ++l;
            lengths[i] = l;
            l = (l == 0 ? 0 : (l - 1));
        }
        assert(lengths.size() == pointers.size());
    }

    void matching_statistics(const char* read, size_t read_length, std::vector<size_t>& lengths, 
                            std::vector<size_t>& pointers, std::vector<size_t>& doc_nums) {  
        // Takes a read, and generates the MS with respect to this ms_t object
        ms.query(read, read_length, pointers, doc_nums, doc_arr);
        lengths.resize(read_length);
        size_t l = 0;
        for (size_t i = 0; i < pointers.size(); ++i) {
            size_t pos = pointers[i];
            while ((i + l) < read_length && (pos + l) < n && (i < 1 || pos != (pointers[i-1] + 1) ) && read[i + l] == ra.charAt(pos + l))
                ++l;
            lengths[i] = l;
            l = (l == 0 ? 0 : (l - 1));
        }
        assert(lengths.size() == pointers.size());
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
 * This section of the code contains the classification methods for reads
 * based on whether it is requested to use MS/PMLs.
 */

size_t classify_reads_pml(pml_t *pml, std::string ref_filename, std::string pattern_filename, bool use_doc, 
                          bool min_digest, bool write_report, size_t num_threads,
                          size_t k, size_t w, bool use_promotions, bool use_dna_letters) {

    // Added for debugging ....
    std::ofstream ks_stat_file (pattern_filename + ".ks_stats");

    // declare output file and iterator
    std::ofstream lengths_file (pattern_filename + ".pseudo_lengths");
    std::ostream_iterator<size_t> lengths_iter (lengths_file, " ");

    std::ofstream doc_file, report_file;
    std::ostream_iterator<size_t> doc_iter (doc_file, " ");

    if (use_doc) {doc_file.open(pattern_filename + ".doc_numbers");}
    if (write_report) {report_file.open(pattern_filename + ".report", std::ofstream::out);}
    KSTest sig_test (ref_filename.data(), PML, write_report, report_file);

    // open query file, and start to classify
    std::ifstream input_file (pattern_filename.c_str());
    omp_set_num_threads(num_threads); 
    size_t num_reads = 0;
    srand(0);

    #pragma omp parallel
    {
        BatchLoader reader;

        // Iterates over batches of data until none left
        while (true) {
            bool valid_batch = true;
            #pragma omp critical // one reader at a time
            {
                valid_batch = reader.loadBatch(input_file, 1000);
            }
            if (!valid_batch) break;

            Read read_struct;
            bool valid_read = false;

            // Iterates over reads in a single batch
            while (true) {
                valid_read = reader.grabNextRead(read_struct);
                if (!valid_read) break;

                // make sure all characters are upper-case
                std::string curr_read = std::string(read_struct.seq);
                transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

                // convert to minimizer-form if needed
                if (use_promotions)
                    curr_read = perform_minimizer_digestion(curr_read, k, w);
                else if (use_dna_letters)
                    curr_read = perform_dna_minimizer_digestion(curr_read, k, w);
                
                // grab MS and write to output file
                std::vector<size_t> lengths, doc_nums;
                if (use_doc){
                    pml->matching_statistics(curr_read.c_str(), curr_read.size(), lengths, doc_nums);
                }
                else {pml->matching_statistics(curr_read.c_str(), curr_read.size(), lengths);}

                // perform the KS-test
                std::vector<double> ks_list;
                std::string status = "";
                size_t num_bin_above_thr = 0;
                double sum_ks_stats = 0.0;

                if (write_report) {
                    // gather the kolomogorov-smirnov statistics
                    ks_list = sig_test.run_kstest(lengths);
                
                    // classify the based on ks-statistics
                    double threshold = sig_test.get_threshold();
                    for (size_t i = 0; i < ks_list.size(); i++) {
                        if (ks_list[i] >= threshold) num_bin_above_thr++;

                        // Added for debugging ....
                        ks_stat_file.precision(3);
                        ks_stat_file << ((i+1.0)/ks_list.size()) << "," << ks_list[i] << "," << threshold << "\n";
                    }
                    bool read_found = (num_bin_above_thr/(ks_list.size()+0.0) > 0.50);

                    std::for_each(ks_list.begin(), ks_list.end(), [&] (double n) {sum_ks_stats += n;});
                    status = (read_found) ? "FOUND" : "NOT_PRESENT";
                }

                #pragma omp atomic
                num_reads++;

                // output the statistics requested
                #pragma omp critical
                {
                    if (use_doc) {
                        doc_file << '>' << read_struct.id << '\n';
                        std::copy(doc_nums.begin(), doc_nums.end(), doc_iter);
                        doc_file << '\n';
                    }
                    lengths_file << '>' << read_struct.id << '\n';
                    std::copy(lengths.begin(), lengths.end(), lengths_iter);
                    lengths_file << '\n'; 

                    if (write_report) {
                        report_file.precision(3);
                        report_file << std::setw(20) << std::left << read_struct.id
                                    << std::setw(15) << std::left << status 
                                    << std::setw(28) << std::left << (sum_ks_stats/ks_list.size()) 
                                    << std::setw(12) << std::left << num_bin_above_thr
                                    << std::setw(12) << std::left << (ks_list.size() - num_bin_above_thr)
                                    << std::endl;
                    }
                }
            } // End of read while loop
        } // End of batch while loop
    } // End of parallel region

    lengths_file.close();
    input_file.close();

    ks_stat_file.close();

    if (use_doc) {doc_file.close();}
    if (write_report) {report_file.close();}
    return num_reads;
}

size_t classify_reads_ms(ms_t *ms, std::string ref_filename, std::string pattern_filename, 
                         bool use_doc, bool min_digest, bool write_report, size_t num_threads,
                         size_t k, size_t w, bool use_promotions, bool use_dna_letters) {

    // declare output files, and output iterators
    std::ofstream lengths_file (pattern_filename + ".lengths");
    std::ofstream pointers_file (pattern_filename + ".pointers");
    std::ofstream doc_file, report_file;

    std::ostream_iterator<size_t> length_iter (lengths_file, " ");
    std::ostream_iterator<size_t> pointers_iter (pointers_file, " ");
    std::ostream_iterator<size_t> doc_iter (doc_file, " ");

    if (use_doc) {doc_file.open(pattern_filename + ".doc_numbers", std::ofstream::out);}
    if (write_report) {report_file.open(pattern_filename + ".report", std::ofstream::out);}
    KSTest sig_test(ref_filename.data(), MS, write_report, report_file);

    // open query file, and start to classify
    std::ifstream input_file (pattern_filename.c_str());
    omp_set_num_threads(num_threads); 
    size_t num_reads = 0;
    srand(0);

    #pragma omp parallel
    {
        BatchLoader reader;

        // Iterates over batches of data until none left
        while (true) {
            bool valid_batch = true;
            #pragma omp critical // one reader at a time
            {
                valid_batch = reader.loadBatch(input_file, 1000);
            }
            if (!valid_batch) break;

            Read read_struct;
            bool valid_read = false;

            // Iterates over reads in a single batch
            while (true) {
                valid_read = reader.grabNextRead(read_struct);
                if (!valid_read) break;

                // make sure all characters are upper-case
                std::string curr_read = std::string(read_struct.seq);
                transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

                // convert to minimizer-form if needed
                if (use_promotions)
                    curr_read = perform_minimizer_digestion(curr_read, k, w);
                else if (use_dna_letters)
                    curr_read = perform_dna_minimizer_digestion(curr_read, k, w);
 
                // grab MS and write to output file
                std::vector<size_t> lengths, pointers, doc_nums;
                if (use_doc){
                    ms->matching_statistics(curr_read.c_str(), curr_read.size(), lengths, pointers, doc_nums);
                }
                else {ms->matching_statistics(curr_read.c_str(), curr_read.size(), lengths, pointers);}

                // perform the KS-test
                std::vector<double> ks_list;
                std::string status = "";
                size_t num_bin_above_thr = 0;
                double sum_ks_stats = 0.0;

                if (write_report) {
                    // gather the kolomogorov-smirnov statistics
                    ks_list = sig_test.run_kstest(lengths);
                
                    // classify the based on ks-statistics
                    double threshold = sig_test.get_threshold();
                    for (size_t i = 0; i < ks_list.size(); i++) {
                        if (ks_list[i] >= threshold) num_bin_above_thr++;
                    }
                    bool read_found = (num_bin_above_thr/(ks_list.size()+0.0) > 0.50);

                    std::for_each(ks_list.begin(), ks_list.end(), [&] (double n) {sum_ks_stats += n;});
                    status = (read_found) ? "FOUND" : "NOT_PRESENT";
                }

                #pragma omp atomic
                num_reads++;

                // output the statistics requested
                #pragma omp critical
                {
                    if (use_doc) {
                        doc_file << '>' << read_struct.id << '\n';
                        std::copy(doc_nums.begin(), doc_nums.end(), doc_iter);
                        doc_file << '\n';
                    }
                    lengths_file << '>' << read_struct.id << '\n';
                    pointers_file << '>' << read_struct.id << '\n';

                    std::copy(lengths.begin(), lengths.end(), length_iter);
                    std::copy(pointers.begin(), pointers.end(), pointers_iter);
                    lengths_file << '\n'; pointers_file << '\n';

                    if (write_report) {
                        report_file.precision(3);
                        report_file << std::setw(20) << std::left << read_struct.id
                                    << std::setw(15) << std::left << status 
                                    << std::setw(28) << std::left << (sum_ks_stats/ks_list.size()) 
                                    << std::setw(12) << std::left << num_bin_above_thr
                                    << std::setw(12) << std::left << (ks_list.size() - num_bin_above_thr)
                                    << std::endl;
                    }
                }
            } // End of read while loop
        } // End of batch while loop
    } // End of parallel region

    input_file.close();
    lengths_file.close();
    pointers_file.close();

    if (use_doc) {doc_file.close();}
    if (write_report) {report_file.close();}
    return num_reads;
}


/*
 * This section contains the "main" methods for the running process where
 * the first method computes the PMLs, and the second one computes the MSs
 * given the index and pattern.
 */

int run_spumoni_main(SpumoniRunOptions* run_opts){
    /* This method is responsible for the PML computation */

    // Loads the RLEBWT and Thresholds
    pml_t ms(run_opts->ref_file, run_opts->use_doc, true);
    std::string out_filename = run_opts->pattern_file;
    std::cout << std::endl;

    // Print out digestion method for input reads
    if (run_opts->use_promotions)
        FORCE_LOG("compute_pml", "input reads will digested using promoted minimizer alphabet (k=%d, w=%d)", 
                  run_opts->k, run_opts->w);
    else if (run_opts->use_dna_letters)
        FORCE_LOG("compute_pml", "input reads will digested using DNA alphabet (k=%d, w=%d)", 
                  run_opts->k, run_opts->w);
    else
        FORCE_LOG("compute_pml", "input reads will be used directly, no minimizer digestion");

    // Process all the reads in the input pattern file
    auto start_time = std::chrono::system_clock::now();
    STATUS_LOG("compute_pml", "processing the patterns");
    
    size_t num_reads = classify_reads_pml(&ms, run_opts->ref_file, run_opts->pattern_file, run_opts->use_doc, 
                                          run_opts->min_digest, run_opts->write_report, run_opts->threads,
                                          run_opts->k, run_opts->w, run_opts->use_promotions, run_opts->use_dna_letters);
    DONE_LOG((std::chrono::system_clock::now() - start_time));
    FORCE_LOG("compute_pml", "finished processing %d reads. results are saved in *.pseudo_lengths file.", num_reads);
    std::cout << std::endl;

    return 0;
}

int run_spumoni_ms_main(SpumoniRunOptions* run_opts) {   
    /* This method is responsible for the MS computation */
    using SelSd = SelectSdvec<>;
    using DagcSd = DirectAccessibleGammaCode<SelSd>;
  
    // Loads the MS index containing the RLEBWT, Thresholds, and RA structure
    ms_t ms(run_opts->ref_file, run_opts->use_doc, true);
    std::string out_filename = run_opts->pattern_file;
    std::cout << std::endl;

    // Print out digestion method for input reads
    if (run_opts->use_promotions)
        FORCE_LOG("compute_ms", "input reads will digested using promoted minimizer alphabet (k=%d, w=%d)", 
                  run_opts->k, run_opts->w);
    else if (run_opts->use_dna_letters)
        FORCE_LOG("compute_ms", "input reads will digested using DNA alphabet (k=%d, w=%d)", 
                  run_opts->k, run_opts->w);
    else
        FORCE_LOG("compute_ms", "input reads will be used directly, no minimizer digestion");

    // Determine approach to parse pattern files
    auto start_time = std::chrono::system_clock::now();
    STATUS_LOG("compute_ms", "processing the reads");

    size_t num_reads = classify_reads_ms(&ms, run_opts->ref_file, run_opts->pattern_file, run_opts->use_doc, 
                                         run_opts->min_digest, run_opts->write_report, run_opts->threads,
                                         run_opts->k, run_opts->w, run_opts->use_promotions, run_opts->use_dna_letters);
    DONE_LOG((std::chrono::system_clock::now() - start_time));
    FORCE_LOG("compute_ms", "finished processing %d reads. results are saved in *.lengths file.", num_reads);
    std::cout << std::endl;
    return 0;
}

std::pair<size_t, size_t> build_spumoni_ms_main(std::string ref_file) {
    // Builds the ms_pointers objects and stores it
    size_t length = 0, num_runs = 0;
    ms_pointers<> ms(ref_file, true);
    std::tie(length, num_runs) = ms.get_bwt_stats(); 

    std::string outfile = ref_file + ms.get_file_extension();
    std::ofstream out(outfile);
    ms.serialize(out);
    out.close();
    return std::make_pair(length, num_runs);
}

std::pair<size_t, size_t> build_spumoni_main(std::string ref_file) {
    // Builds the pml_pointers objects and stores it
    size_t length = 0, num_runs = 0;
    pml_pointers<> pml(ref_file, true);
    std::tie(length, num_runs) = pml.get_bwt_stats();

    std::string outfile = ref_file + pml.get_file_extension();
    std::ofstream out(outfile);
    pml.serialize(out);
    out.close();
    return std::make_pair(length, num_runs);
}

void generate_null_ms_statistics(std::string ref_file, std::string pattern_file, std::vector<size_t>& ms_stats,
                                 bool min_digest, bool use_promotions, bool use_dna_letters, size_t k, size_t w) {
    /* Generates the null ms statistics and returns them to be saved */

    // Loads the index, and needed variables
    ms_t ms_index(ref_file, false);
    gzFile fp = gzopen(pattern_file.data(), "r");
    kseq_t* seq = kseq_init(fp);

    // Iterate through sample reads and generate null MS
    while (kseq_read(seq)>=0) {
        //Make sure all characters are upper-case
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

        // Reverse string to make it a null read
        std::reverse(curr_read.begin(), curr_read.end());
        
        // Convert to minimizer-form if needed
        if (use_promotions)
            curr_read = perform_minimizer_digestion(curr_read, k, w);
        else if (use_dna_letters)
            curr_read = perform_dna_minimizer_digestion(curr_read, k, w);
        
        // Generate the null MS
        std::vector<size_t> lengths, pointers;
        ms_index.matching_statistics(curr_read.c_str(), curr_read.length(), lengths, pointers);
        ms_stats.insert(ms_stats.end(), lengths.begin(), lengths.end());
    }
    kseq_destroy(seq);
    gzclose(fp);
}

void generate_null_pml_statistics(std::string ref_file, std::string pattern_file, std::vector<size_t>& pml_stats,
                                 bool min_digest, bool use_promotions, bool use_dna_letters, size_t k, size_t w) {
    /* Generates the null pml statistics and returns them to be saved */

    // Load the indexes, and needed variables
    pml_t pml_index(ref_file, false);
    gzFile fp = gzopen(pattern_file.data(), "r");
    kseq_t* seq = kseq_init(fp);

    // Iterates through sample reads, and generates PML reads
    while (kseq_read(seq)>=0) {
        //Make sure all characters are upper-case
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

        // Reverse string to make it a null read
        std::reverse(curr_read.begin(), curr_read.end());

        // Convert to minimizer-form if needed
        if (use_promotions)
            curr_read = perform_minimizer_digestion(curr_read, k, w);
        else if (use_dna_letters)
            curr_read = perform_dna_minimizer_digestion(curr_read, k, w);

        // Generate the null PML
        std::vector<size_t> lengths;
        pml_index.matching_statistics(curr_read.c_str(), curr_read.length(), lengths);
        pml_stats.insert(pml_stats.end(), lengths.begin(), lengths.end());
    }
    kseq_destroy(seq);
    gzclose(fp);
}

void find_threshold_based_on_null_pml_distribution(const char* ref_file, const char* null_reads, bool use_minimizers,
                                                   bool use_promotions, bool use_dna_letters, size_t k, size_t w, EmpNullDatabase& null_db) {

    /* Generates a distribution of KS-stats from null reads to determine the optimal threshold */

    // Load the indexes, and needed variables
    pml_t pml_index(ref_file, false);
    gzFile fp = gzopen(null_reads, "r");
    kseq_t* seq = kseq_init(fp);

    // Iterates through null reads, and generates PML 
    KSTest sig_test(null_db, PML);
    std::vector<size_t> pml_stats;
    std::vector <double> ks_list;
    double ks_stat_sum = 0.0;

    while (kseq_read(seq)>=0) {
        //Make sure all characters are upper-case
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

        // Reverse string to make it a null read
        std::reverse(curr_read.begin(), curr_read.end());

        // Convert to minimizer-form if needed
        if (use_promotions)
            curr_read = perform_minimizer_digestion(curr_read, k, w);
        else if (use_dna_letters)
            curr_read = perform_dna_minimizer_digestion(curr_read, k, w);

        // Generate the null PML
        std::vector<size_t> lengths;
        pml_index.matching_statistics(curr_read.c_str(), curr_read.length(), lengths);
        pml_stats.insert(pml_stats.end(), lengths.begin(), lengths.end());

        // Generate the KS-statistics
        std::vector<double> curr_ks_list;
        curr_ks_list = sig_test.run_kstest(lengths);
        ks_list.insert(ks_list.end(), curr_ks_list.begin(), curr_ks_list.end());
        std::for_each(curr_ks_list.begin(), curr_ks_list.end(), [&] (double x) {ks_stat_sum += x;});
    }
    kseq_destroy(seq);
    gzclose(fp);

    // find the variance of the ks-statistics
    double sum = 0.0;
    double mean = ks_stat_sum/ks_list.size();

    for (size_t i = 0; i < ks_list.size(); i++){
        sum += std::pow(ks_list[i] - mean, 2);
    }
    double std_dev = std::sqrt(sum/ks_list.size());

    // set the threshold as 3 std. deviations above mean
    null_db.ks_stat_threshold = mean + (3 * std_dev);
}

void find_threshold_based_on_null_ms_distribution(const char* ref_file, const char* null_reads, bool use_minimizers,
                                                   bool use_promotions, bool use_dna_letters, size_t k, size_t w, EmpNullDatabase& null_db) {

    /* Generates a distribution of KS-stats from null reads to determine the optimal threshold */

    // Load the indexes, and needed variables
    ms_t ms_index(ref_file, false);
    gzFile fp = gzopen(null_reads, "r");
    kseq_t* seq = kseq_init(fp);

    // Iterates through null reads, and generates MS
    KSTest sig_test(null_db, MS);
    std::vector<size_t> ms_stats;
    std::vector <double> ks_list;
    double ks_stat_sum = 0.0;

    while (kseq_read(seq)>=0) {
        //Make sure all characters are upper-case
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

        // Reverse string to make it a null read
        std::reverse(curr_read.begin(), curr_read.end());

        // Convert to minimizer-form if needed
        if (use_promotions)
            curr_read = perform_minimizer_digestion(curr_read, k, w);
        else if (use_dna_letters)
            curr_read = perform_dna_minimizer_digestion(curr_read, k, w);

        // Generate the null PML
        std::vector<size_t> lengths, pointers;
        ms_index.matching_statistics(curr_read.c_str(), curr_read.length(), lengths, pointers);
        ms_stats.insert(ms_stats.end(), lengths.begin(), lengths.end());

        // Generate the KS-statistics
        std::vector<double> curr_ks_list;
        curr_ks_list = sig_test.run_kstest(lengths);
        ks_list.insert(ks_list.end(), curr_ks_list.begin(), curr_ks_list.end());
        std::for_each(curr_ks_list.begin(), curr_ks_list.end(), [&] (double x) {ks_stat_sum += x;});
    }
    kseq_destroy(seq);
    gzclose(fp);

    // find the variance of the ks-statistics
    double sum = 0.0;
    double mean = ks_stat_sum/ks_list.size();

    for (size_t i = 0; i < ks_list.size(); i++){
        sum += std::pow(ks_list[i] - mean, 2);
    }
    double std_dev = std::sqrt(sum/ks_list.size());

    // set the threshold as 3 std. deviations above mean
    null_db.ks_stat_threshold = mean + (3 * std_dev);
}