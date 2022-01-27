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
#include <tuple>
#include <fstream>
#include <iterator>

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
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();
        thresholds = thresholds_t(filename,&this->bwt);
        t_insert_end = std::chrono::high_resolution_clock::now();

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

    std::vector<size_t> query(const char* pattern, const size_t m, std::vector<size_t>& lengths) {
        return _query(pattern, m, lengths);
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
    std::vector<size_t> _query(const string_t &pattern, const size_t m, std::vector<size_t>& lengths) {
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
    int_vector<> samples_start;
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
        t_insert_end = std::chrono::high_resolution_clock::now();

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

    void query(const char* pattern, const size_t m, std::vector<size_t>& pointers) {
        _query(pattern, m, pointers);
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
    template<typename string_t>
    void _query(const string_t &pattern, const size_t m, std::vector<size_t>& ms_pointers) 
    {
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

        //ms_pointers;
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
    // Constructor
    pml_t(std::string filename){
      SPUMONI_LOG("Loading the PML index ...");
      auto start_time = std::chrono::system_clock::now();
      std::string filename_ms = filename + ms.get_file_extension();

      ifstream fs_ms(filename_ms);
      ms.load(fs_ms);
      fs_ms.close();

      auto end_time = std::chrono::system_clock::now();
      TIME_LOG((end_time - start_time));
    }

    //Destructor
    ~pml_t() {}

    void matching_statistics(const char* read, size_t read_length, std::vector<size_t>& lengths) {
        ms.query(read, read_length, lengths);
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

    ms_t(std::string filename) {
        SPUMONI_LOG("Loading the MS index ...");
        auto start_time = std::chrono::system_clock::now();    
        std::string filename_ms = filename + ms.get_file_extension();

        ifstream fs_ms(filename_ms);
        ms.load(fs_ms);
        fs_ms.close();

        auto end_time = std::chrono::system_clock::now();
        TIME_LOG((end_time - start_time));

        SPUMONI_LOG("Loading the Random Access Data Structure ...");
        start_time = std::chrono::system_clock::now();   
        std::string filename_slp = filename + ".slp";

        ifstream fs(filename_slp);
        ra.load(fs);
        fs.close();

        n = ra.getLen();
        end_time = std::chrono::system_clock::now();
        TIME_LOG((end_time - start_time));
    } 

    // Destructor
    ~ms_t() {}

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

        // Just declared holder to avoid compilation error, will replace method in future
        std::vector<size_t> holder;
        p->ms->matching_statistics(curr_read.c_str(), seq->seq.l, holder);
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

        // replaced out_fd with ""
        //p->ms->matching_statistics(curr_read.c_str(), seq->seq.l, "");

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
    // Declare output file and iterator
    std::ofstream lengths_file (out_filename + ".pseudo_lengths");
    std::ostream_iterator<size_t> lengths_iter (lengths_file, " ");

    // Use kseq to parse out sequences from FASTA file
    gzFile fp = gzopen(pattern_filename.c_str(), "r");
    kseq_t* seq = kseq_init(fp);
    size_t num_reads = 0;

    while (kseq_read(seq) >= 0) {
        //Make sure all characters are upper-case
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 
        lengths_file << '>' << seq->name.s << '\n';

        // Grab MS and write to output file
        std::vector<size_t> lengths;
        ms->matching_statistics(curr_read.c_str(), seq->seq.l, lengths);

        std::copy(lengths.begin(), lengths.end(), lengths_iter);
        lengths_file << '\n';
        num_reads++;
    }

    kseq_destroy(seq);
    gzclose(fp);
    lengths_file.close();
    return num_reads;
}

size_t st_pml_general(pml_t *ms, std::string pattern_filename, std::string out_filename) {
    // Check if file is mistakenly labeled as non-fasta file
    std::string file_ext = pattern_filename.substr(pattern_filename.find_last_of(".") + 1);
    if (file_ext == "fa" || file_ext == "fasta") {
        FATAL_WARNING("The file extension for the patterns suggests it is a fasta file.\n" 
                      "Please run with -f option for correct results.");
    }

    // Declare output files, and output iterators
    std::ofstream lengths_file (out_filename + ".pseudo_lengths");
    std::ostream_iterator<size_t> length_iter (lengths_file, " ");

    // Open pattern file to start reading reads
    std::ifstream input_fd (pattern_filename, std::ifstream::in | std::ifstream::binary);
    std::vector<size_t> lengths;
    std::string read = "";
    size_t num_reads = 0;

    char ch = input_fd.get();
    char* buf = new char [2]; // To just hold separator character

    while (input_fd.good()) {
        // Finds a separating character
        if (ch == '\x01') { 
            ms->matching_statistics(read.c_str(), read.size(), lengths);
            lengths_file << ">read_" << num_reads << "\n";
  
            std::copy(lengths.begin(), lengths.end(), length_iter);
            lengths_file << "\n";

            input_fd.read(buf, 2); 
            num_reads++;
            read="";
        }
        else {read += ch;}
        ch = input_fd.get();
    }
    lengths_file.close();
    return num_reads;
}

size_t st_ms(ms_t *ms, std::string pattern_filename, std::string out_filename) {
    // Declare output files, and output iterators
    std::ofstream lengths_file (out_filename + ".lengths");
    std::ofstream pointers_file (out_filename + ".pointers");

    std::ostream_iterator<size_t> length_iter (lengths_file, " ");
    std::ostream_iterator<size_t> pointers_iter (pointers_file, " ");

    // Use kseq to parse out sequences from FASTA file
    gzFile fp = gzopen(pattern_filename.c_str(), "r");
    kseq_t* seq = kseq_init(fp);
    size_t num_reads = 0;

    while (kseq_read(seq) >= 0) {
        //Make sure all characters are upper-case
        std::string curr_read = std::string(seq->seq.s);
        transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); 

        lengths_file << '>' << seq->name.s << '\n';
        pointers_file << '>' << seq->name.s << '\n';

        // Grab MS and write to output file
        std::vector<size_t> lengths, pointers;
        ms->matching_statistics(curr_read.c_str(), seq->seq.l, lengths, pointers);

        std::copy(lengths.begin(), lengths.end(), length_iter);
        std::copy(pointers.begin(), pointers.end(), pointers_iter);
        
        lengths_file << '\n'; pointers_file << '\n';
        num_reads++;
    }

    kseq_destroy(seq);
    gzclose(fp);
    lengths_file.close();
    pointers_file.close();
    return num_reads;
}

size_t st_ms_general(ms_t *ms, std::string pattern_filename, std::string out_filename)
{
    // Check if file is mistakenly labeled as non-fasta file
    std::string file_ext = pattern_filename.substr(pattern_filename.find_last_of(".") + 1);
    if (file_ext == "fa" || file_ext == "fasta") {
        FATAL_WARNING("The file extension for the patterns suggests it is a fasta file.\n" 
                      "Please run with -f option for correct results.");
    }

    // Declare output files, and output iterators
    std::ofstream lengths_file (out_filename + ".lengths");
    std::ofstream pointers_file (out_filename + ".pointers");

    std::ostream_iterator<size_t> length_iter (lengths_file, " ");
    std::ostream_iterator<size_t> pointers_iter (pointers_file, " ");

    // Open pattern file to start reading reads
    std::ifstream input_fd (pattern_filename, std::ifstream::in | std::ifstream::binary);
    std::vector<size_t> lengths, pointers;
    std::string read = "";
    size_t num_reads = 0;

    char ch = input_fd.get();
    char* buf = new char [2]; // To just hold separator character

    while (input_fd.good()) {
        // Finds a separating character
        if (ch == '\x01') { 
            ms->matching_statistics(read.c_str(), read.size(), lengths, pointers);
            lengths_file << ">read_" << num_reads << "\n";
            pointers_file << ">read_" << num_reads << "\n";

            std::copy(lengths.begin(), lengths.end(), length_iter);
            std::copy(pointers.begin(), pointers.end(), pointers_iter);
            lengths_file << "\n"; pointers_file << "\n";

            input_fd.read(buf, 2); 
            num_reads++;
            read="";
        }
        else {read += ch;}
        ch = input_fd.get();
    }
    lengths_file.close();
    pointers_file.close();
    return num_reads;
}

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;


/*
 * This section contains the "main" methods for the running process where
 * the first method computes the PMLs, and the second one computes the MSs
 * given the index and pattern.
 */

int run_spumoni_main(SpumoniRunOptions* run_opts){
    /* This method is responsible for the PML computation */

    // Loads the RLEBWT and Thresholds
    pml_t ms(run_opts->ref_file);
    std::string out_filename = run_opts->pattern_file;

    if (is_gzipped(run_opts->pattern_file)) { 
        SPUMONI_LOG("The input is gzipped - forcing single threaded pseudo matching lengths.");
        run_opts->threads = 1;
    }

    // Determine approach to parse pattern files
    auto start_time = std::chrono::system_clock::now();
    SPUMONI_LOG("Processing the patterns ...");

    size_t num_reads = 0;
    if (run_opts->query_fasta) {
        if(run_opts->threads <= 1) {num_reads=st_pml(&ms, run_opts->pattern_file, out_filename);}
        else {FATAL_WARNING("Multi-threading not implemented yet for FASTA querying.");}
    }
    else {
        if(run_opts->threads == 1) {num_reads=st_pml_general(&ms, run_opts->pattern_file,out_filename);}
        else {FATAL_WARNING("Multi-threading not implemented yet for general-text querying.");}
    }

    auto end_time = std::chrono::system_clock::now();
    TIME_LOG((end_time - start_time));

    SPUMONI_LOG("Finished processing %d reads", num_reads);
    return 0;
}

int run_spumoni_ms_main(SpumoniRunOptions* run_opts) {   
    /* This method is responsible for the MS computation */
    using SelSd = SelectSdvec<>;
    using DagcSd = DirectAccessibleGammaCode<SelSd>;
  
    // Loads the MS index containing the RLEBWT, Thresholds, and RA structure
    ms_t ms(run_opts->ref_file);
    std::string out_filename = run_opts->pattern_file;

    if (is_gzipped(run_opts->pattern_file)) {
        SPUMONI_LOG("The input is gzipped - forcing single threaded matching statistics.");
        run_opts->threads = 1;
    }

    // Determine approach to parse pattern files
    auto start_time = std::chrono::system_clock::now();
    SPUMONI_LOG("Processing the patterns ...");

    size_t num_reads = 0;
    if (run_opts->query_fasta) {
        if(run_opts->threads == 1) {num_reads=st_ms(&ms, run_opts->pattern_file, out_filename);}
        else {FATAL_WARNING("Multi-threading not implemented yet for FASTA querying.");}
    }
    else {
        if(run_opts->threads == 1) {num_reads=st_ms_general(&ms,run_opts->pattern_file,out_filename);}
        else {FATAL_WARNING("Multi-threading not implemented yet for general-text querying.");}
    }
 
    auto end_time = std::chrono::system_clock::now();
    TIME_LOG((end_time - start_time));

    SPUMONI_LOG("Finished processing %d reads", num_reads);
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