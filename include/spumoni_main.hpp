/* 
 *  Copyright (C) 2020 Omar Ahmed
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
  * File: spumoni.cpp 
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


/* Commonly Used MACROS */
#define SPUMONI_VERSION "1.0"
#define NOT_IMPL(x) do { std::fprintf(stderr, "%s is not implemented: %s\n", __func__, x); std::exit(1);} while (0)
#define FATAL_WARNING(x) do {std::fprintf(stderr, "Warning: %s\n", x); std::exit(1);} while (0)
#define THROW_EXCEPTION(x) do { throw x;} while (0)
#define SPUMONI_LOG(...) do{std::fprintf(stderr, "[spumoni] "); std::fprintf(stderr, __VA_ARGS__);\
                            std::fprintf(stderr, "\n");} while(0)
#define TIME_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "[spumoni] Elapsed Time (s): %.3f\n", sec.count());} while(0)
#define OTHER_LOG(x) do {std::stringstream str(x); std::string str_out;\
                         while (std::getline(str, str_out, '\n')) { \
                         std::fprintf(stderr, "[helper-prog] %s\n", str_out.data());}} while (0)


/* DEBUG macros */
#define DEBUG 1
#define DBG_ONLY(...)  do { if (DEBUG) {std::fprintf(stderr, "[DEBUG] "); \
                            std::fprintf(stderr, __VA_ARGS__);}} while (0)

/* Function Declarations */
int spumoni_build_usage();
int build_main(int argc, char** argv);
int run_main(int argc, char** argv);
int spumoni_usage ();
int is_file(std::string path);
std::string execute_cmd(const char* cmd);
size_t get_avail_phy_mem();
int spumoni_run_usage ();

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
  std::string ms_build = "rlebwt_ms_build";
  std::string pml_build = "build_spumoni";
  
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
      ms_build.assign(base_path + ms_build);
      pml_build.assign(base_path + pml_build);
  }

  void validate() const {
      /* Makes sure that each path for an executable is valid */
      bool invalid_path = !is_file(parseNT_bin) | !is_file(parse_fasta_bin) | !is_file(parse_bin) | !is_file(pfp_thresholds);
      invalid_path = invalid_path | !is_file(pfp_thresholds64) | !is_file(compress_bin) | !is_file(preprocess_dict_bin) | !is_file(repair_bin);
      invalid_path = invalid_path | !is_file(postprocess_gram_bin) | !is_file(shaped_slp_bin) | !is_file(ms_build) | !is_file(pml_build);
      if (invalid_path) {THROW_EXCEPTION(std::runtime_error("One or more of helper program paths are invalid."));}
  }
};

struct SpumoniBuildOptions {
  std::string ref_file = "";
  size_t wind = 10; // sliding window size
  size_t hash_mod = 100; // hash modulus
  size_t threads = 0; // number of helper threads
  bool keep_files = false; // keeps temporary files
  bool ms_index = false; // want ms index
  bool pml_index = false; // want pml index
  bool stop_after_parse = false; // stop build after build parse
  bool compress_parse = false; // compress parse
  bool is_fasta = false; // reference is fasta

public:
  void validate() const {
        /* Checks if the parse arguments are valid for continuing the execution */
        if (!ms_index && !pml_index) {
            FATAL_WARNING("Need to specify what type of quantities you would like to compute with -M, -P or both.");
        }
        if (ref_file.find(".fq") != std::string::npos || ref_file.find(".fastq") != std::string::npos || ref_file.find(".fnq") != std::string::npos) {
            FATAL_WARNING("Reference file cannot be in FASTQ format.\n");
        }
        if (!is_file(ref_file)) {THROW_EXCEPTION(std::runtime_error("The following path is not valid: " + ref_file));}
    }
};

enum output_type {MS, PML, NOT_CHOSEN};

struct SpumoniRunOptions {
  std::string ref_file = ""; // reference file
  std::string pattern_file = ""; // pattern file
  bool ms_requested = false; // user wants to compute MS
  bool pml_requested = false; // user wants to compute PML
  bool query_fasta = false; // query file is a fasta
  output_type result_type = NOT_CHOSEN; // output type requested by user
  size_t threads = 1; // number of TOTAL threads

public:
  void populate_output_type() {
      /* Populates the output type member of the struct */
      if (ms_requested && !pml_requested) {result_type = MS;}
      if (!ms_requested && pml_requested) {result_type = PML;}
  }
  
  void validate() const {
      /* Checks the options for the run command, and makes sure it has everything it needs */
      if (ref_file == "" || pattern_file == ""){FATAL_WARNING("Both a reference file (-r) and pattern file (-f) must be provided.");}
      if (result_type == NOT_CHOSEN) {FATAL_WARNING("An output type with -M or -P must be specified, only one can be used at a time.");}

      if (!is_file(ref_file)) {THROW_EXCEPTION(std::runtime_error("The following path is not valid: " + ref_file));}
      if (!is_file(pattern_file)) {THROW_EXCEPTION(std::runtime_error("The following path is not valid: " + pattern_file));}
      if (ms_requested && pml_requested) {FATAL_WARNING("Only MS or PMLs can be computed at one time, please re-run with only -M or -P");}
      if (pattern_file.find(".fq") != std::string::npos || pattern_file.find(".fastq") != std::string::npos || pattern_file.find(".fnq") != std::string::npos) {
            FATAL_WARNING("Patterns file cannot be in FASTQ format.\n");}
      
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
  }
};

/* Additional Function Declarations */
void parse_build_options(int argc, char** argv, SpumoniBuildOptions* opts);
void parse_run_options(int argc, char** argv, SpumoniRunOptions* opts);



#endif /* End of SPUMONI_MAIN_H */