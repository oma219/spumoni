 /*
  * File: spumoni.hpp 
  * Description: Header file for main execution file for SPUMONI.
  *             
  * Author: Omar Ahmed
  * Start Date: October 10, 2021
  */

#ifndef SPUMONI_MAIN_H
#define SPUMONI_MAIN_H
#include <unistd.h>
#include <fstream>
#include <chrono>
#include <vector>
#include <stdlib.h>

/* Commonly Used MACROS */
#define SPUMONI_VERSION "1.2.0"
#define NOT_IMPL(x) do { std::fprintf(stderr, "%s is not implemented: %s\n", __func__, x); std::exit(1);} while (0)
#define THROW_EXCEPTION(x) do { throw x;} while (0)
#define FATAL_WARNING(x) do {std::fprintf(stderr, "Warning: %s\n\n", x); std::exit(1);} while (0)
#define FATAL_ERROR(...) do {std::fprintf(stderr, "Error: "); std::fprintf(stderr, __VA_ARGS__);\
                              std::fprintf(stderr, "\n\n"); std::exit(1);} while(0)
#define SPUMONI_LOG(...) do{std::fprintf(stderr, "[spumoni] "); std::fprintf(stderr, __VA_ARGS__);\
                            std::fprintf(stderr, "\n");} while(0)

// General logging with option to silence based on command line options
#define LOG(verbose, func, log)  if (verbose) {std::fprintf(stderr, "[%s] ", func); \
                                               std::fprintf(stderr, log); \
                                               std::fprintf(stderr, "\n");} 
#define FORCE_LOG(func, ...)  do {std::fprintf(stderr, "[%s] ", func); \
                                  std::fprintf(stderr, __VA_ARGS__); \
                                  std::fprintf(stderr, "\n");} while (0)

// Logging that keeps track of how long commands take: call STATUS_LOG, followed by DONE_LOG
#define STATUS_LOG(x, ...) do {std::fprintf(stderr, "[%s] ", x); std::fprintf(stderr, __VA_ARGS__ ); \
                               std::fprintf(stderr, " ... ");} while(0)
#define DONE_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "done.  (%.3f sec)\n", sec.count());} while(0)

#define ASSERT(condition, msg) do {if (!condition){std::fprintf(stderr, "Assertion Failed: %s\n", msg); \
                                                   std::exit(1);}} while(0)
#define TIME_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "[spumoni] Elapsed Time (s): %.3f\n", sec.count());} while(0)
#define OTHER_LOG(x) if (DEBUG) {std::stringstream str(x); std::string str_out;\
                                 while (std::getline(str, str_out, '\n')) { \
                                 std::fprintf(stderr, "[helper-prog] %s\n", str_out.data());}} 

/* DEBUG macros */
#define DEBUG 0
#define DBG_ONLY(...)  do { if (DEBUG) {std::fprintf(stderr, "\n[DEBUG] "); \
                            std::fprintf(stderr, __VA_ARGS__); std::fprintf(stderr, "\n");}} while (0)

/* Type Definitions */
typedef uint64_t ulint;

/* Value Definitions */
#define THRBYTES 5 
#define SSABYTES 5 
#define NULL_READ_CHUNK 150
#define NUM_NULL_READS 100 // 7500 = 75 bp * 100 reads
#define NULL_READ_BOUND 200
#define KS_STAT_MS_THR 0.25
#define KS_STAT_PML_THR 0.10

/* Function Declarations */
int spumoni_build_usage();
int build_main(int argc, char** argv);
int run_main(int argc, char** argv);
int spumoni_usage ();
int is_file(std::string path);
int is_dir(std::string path);
bool is_integer(const std::string& str);
std::vector<std::string> split(std::string input, char delim);
bool endsWith(const std::string& str, const std::string& suffix);
std::string execute_cmd(const char* cmd);
size_t get_avail_phy_mem();
int spumoni_run_usage ();
std::string perform_minimizer_digestion(std::string input_query, size_t k, size_t w);
std::string perform_dna_minimizer_digestion(std::string input_query, size_t k, size_t w);

struct SpumoniHelperPrograms {
  /* Contains paths to run helper programs */
  std::string base_path = "";
  std::string parseNT_bin = "newscanNT.x";
  std::string parse_fasta_bin = "newscan.x";
  std::string parse_bin = "pscan.x";
  std::string pfp_thresholds = "pfp_thresholds";
  std::string pfp_thresholds64 = "pfp_thresholds64";
  std::string compress_bin = "compress_dictionary";
  std::string preprocess_dict_bin = "procdic";
  std::string repair_bin = "irepair";
  std::string postprocess_gram_bin = "postproc";
  std::string shaped_slp_bin = "SlpEncBuild";
  //std::string ms_build = "rlebwt_ms_build";
  //std::string pml_build = "build_spumoni";
  
public:
  void build_paths(std::string base) {
      /* Takes the base path, and combines it with names to generate executable paths */
      base_path.assign(base);
      parseNT_bin.assign(base_path + parseNT_bin);
      parse_fasta_bin.assign(base_path + parse_fasta_bin);
      parse_bin.assign(base_path + parse_bin);
      pfp_thresholds.assign(base_path + pfp_thresholds);
      pfp_thresholds64.assign(base_path + pfp_thresholds64);
      compress_bin.assign(base_path + compress_bin);
      preprocess_dict_bin.assign(base_path + preprocess_dict_bin);
      repair_bin.assign(base_path + repair_bin);
      postprocess_gram_bin.assign(base_path + postprocess_gram_bin);
      shaped_slp_bin.assign(base_path + shaped_slp_bin);
      //ms_build.assign(base_path + ms_build);
      //pml_build.assign(base_path + pml_build);
  }

  void validate() const {
      /* Makes sure that each path for an executable is valid */
      bool invalid_path = !is_file(parseNT_bin) | !is_file(parse_fasta_bin) | !is_file(parse_bin) | !is_file(pfp_thresholds);
      invalid_path = invalid_path | !is_file(pfp_thresholds64) | !is_file(compress_bin) | !is_file(preprocess_dict_bin) | !is_file(repair_bin);
      invalid_path = invalid_path | !is_file(postprocess_gram_bin) | !is_file(shaped_slp_bin); // | !is_file(ms_build) | !is_file(pml_build);
      if (invalid_path) {THROW_EXCEPTION(std::runtime_error("One or more of helper program paths are invalid."));}
  }
};

enum output_type {MS, PML, NOT_CHOSEN};
enum reference_type {FASTA, MINIMIZER, NOT_SET};
enum query_input_type {FA, FQ, NOT_CLEAR};

struct SpumoniBuildOptions {
  std::string ref_file = "";
  std::string input_list = "";
  std::string output_dir = "";
  size_t wind = 10; // sliding window size
  size_t hash_mod = 100; // hash modulus
  size_t threads = 0; // number of helper threads
  bool keep_files = false; // keeps temporary files
  bool ms_index = false; // want ms index
  bool pml_index = false; // want pml index
  bool verbose = false;
  //output_type index_type = NOT_CHOSEN; // the actual index that will be build
  bool stop_after_parse = false; // stop build after build parse
  bool compress_parse = false; // compress parse
  bool is_fasta = false; // reference is binary by default, since 
                         // we use minimizers by default (default: true)
  bool build_doc = false; // build the document array
  bool use_minimizers = true; // digest sequences into minimizers
  bool use_promotions = false; // use alphabet promotion during promotion
  bool use_dna_letters = false; // use DNA-letter based minimizers
  size_t k = 4; // small window size for minimizers
  size_t w = 12; // large window size for minimizers

public:
  void validate() {
      /* Checks if the parse arguments are valid for continuing the execution */

      // Check based on approach used: file-list or single reference file
      if (ref_file.length()){
          if (!is_file(ref_file)) {FATAL_ERROR("The following path is not valid: %s", ref_file.data());}  
          if (output_dir.length()){FATAL_ERROR("The -b option should not be set when using a single file.");}
          if (!endsWith(ref_file, ".fa") && !endsWith(ref_file, ".fasta") && !endsWith(ref_file, ".fna")){
                 FATAL_ERROR("The reference file provided does not appear to be a FASTA\n" 
                             "       file, please convert to FASTA and re-run.");}  
      } else {
          if (!is_file(input_list)) {FATAL_ERROR("The following path is not valid: %s", input_list.data());} 
          if (!output_dir.length()) {FATAL_ERROR("You must specify an output directory with -b option when using filelist.");}
          if (!is_dir(output_dir)){FATAL_ERROR("The output directory for the index is not valid.");}
      }
      if (build_doc && ref_file.length()) {
        FATAL_ERROR("Cannot build a document array if you are indexing a single file.");}
      
      // Check if we only set one type minimizers
      if (use_minimizers) {
        if (use_promotions && use_dna_letters) {FATAL_ERROR("Only one type of minimizer can be specified.");}
        if (!use_promotions && !use_dna_letters) {FATAL_ERROR("A minimizer type must be specified.");}
      } else {
        if (use_promotions || use_dna_letters) {FATAL_ERROR("A minimizer type should not be specified if intending not to use minimizer digestion.");}
      }

      // Makes sure that at least one index type is chosen ...
      if (!ms_index && !pml_index) {FATAL_ERROR("At least one index type (-M or -P) must be specified for build.");}

      // Check the values for k and w
      if (k > 4) {FATAL_WARNING("small window size (k) cannot be larger than 4 characters.");}
      if (w < k) {FATAL_WARNING("large window size (w) should be larger than the small window size (k)");}
  }
};

struct SpumoniRunOptions {
  std::string ref_file = ""; // reference file
  std::string pattern_file = ""; // pattern file
  bool ms_requested = false; // user wants to compute MS
  bool pml_requested = false; // user wants to compute PML
  output_type result_type = NOT_CHOSEN; // output type requested by user
  reference_type ref_type = NOT_SET; // the type of reference
  size_t threads = 1; // number of threads
  bool use_doc = false; // build the document array
  bool write_report = false; // write out the classification report
  bool min_digest = true; // need to digest reads (default is true) 
  bool use_promotions = false; // use alphabet promotion during promotion
  bool use_dna_letters = false; // use DNA-letter based minimizers
  size_t k = 4; // small window size for minimizers
  size_t w = 12; // large window size for minimizers

public:
  void populate_types() {
      /* Populates the output type member of the struct */

      // Specify the output type requested
      if (ms_requested && !pml_requested) {result_type = MS;}
      if (!ms_requested && pml_requested) {result_type = PML;}

      // Specify the input database type
      bool is_fasta = (endsWith(ref_file, ".fa") || endsWith(ref_file, ".fasta") || endsWith(ref_file, ".fna"));
      bool is_min = endsWith(ref_file, ".bin");

      if (is_fasta && !is_min) {ref_type = FASTA;}
      if (!is_fasta && is_min) {ref_type = MINIMIZER;}
  }
  
  void validate() const {
      /* Checks the options for the run command, and makes sure it has everything it needs */
      if (ref_file == "" || pattern_file == ""){FATAL_WARNING("Both a reference file (-r) and pattern file (-p) must be provided.");}
      if (result_type == NOT_CHOSEN) {FATAL_WARNING("An output type with -M or -P must be specified, only one can be used at a time.");}
      
      // Make sure provided files are valid
      if (!is_file(ref_file)) {FATAL_ERROR("The following path is not valid: %s", ref_file.data());}
      if (!is_file(pattern_file)) {FATAL_ERROR("The following path is not valid: %s", pattern_file.data());}

      // Make sure reference file is a valid type
      if (ref_type == NOT_SET) {FATAL_ERROR("Reference file is an unrecognized type. It needs to be a\n"
                                            "       FASTA file or binary file produced by spumoni build.");}

      // Make sure query file is a FASTA file                                     
      if (!endsWith(pattern_file, ".fa") && !endsWith(pattern_file, ".fasta") && !endsWith(pattern_file, ".fna")){
          FATAL_ERROR("The pattern file provided does not appear to be a FASTA\n" 
                      "       file, please convert to FASTA and re-run.");}
      
      // Verify doc array is available, if needed
      if (use_doc && !is_file(ref_file + ".doc")) {FATAL_WARNING("*.doc file is not present, so it cannot be used");}
      
      switch (result_type) {
        case MS: 
            if (!is_file(ref_file+".thrbv.ms")) 
            {FATAL_WARNING("The index required for this computation is not available, please use spumoni build.");} break;
        case PML:
            if (!is_file(ref_file+".thrbv.spumoni")) 
            {FATAL_WARNING("The index required for this computation is not available, please use spumoni build.");} break;
        default:
            FATAL_WARNING("An output type with -M or -P must be specified, only one can be used at a time."); break;
      }

      // Check the values for k and w
      if (k > 4) {FATAL_WARNING("small window size (k) cannot be larger than 4 characters.");}
      if (w < k) {FATAL_WARNING("large window size (w) should be larger than the small window size (k)");}

      // Check if we choose the minimizer type correctly
      if (min_digest) {
        if (use_promotions && use_dna_letters) {FATAL_ERROR("Only one type of minimizer can be specified from either -m or -a.");}
        if (!use_promotions && !use_dna_letters) {FATAL_ERROR("A minimizer type must be specified using -m or -a.");}
      } else {
        if (use_promotions || use_dna_letters) {FATAL_ERROR("A minimizer type should not be specified if intending not to use minimizer digestion.");}
      }

  }
};

/* Additional Function Declarations */
void parse_build_options(int argc, char** argv, SpumoniBuildOptions* opts);
void parse_run_options(int argc, char** argv, SpumoniRunOptions* opts);

#endif /* End of SPUMONI_MAIN_H */