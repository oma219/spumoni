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

/* Commonly Used MACROS */
#define SPUMONI_VERSION "1.0"
#define NOT_IMPL(x) do { std::fprintf(stderr, "%s is not implemented: %s\n", __func__, x); std::exit(1);} while (0)
#define FATAL_WARNING(x) do {std::fprintf(stderr, "Warning: %s\n", x); std::exit(1);} while (0)
#define THROW_EXCEPTION(x) do { throw x;} while (0)
#define SPUMONI_LOG(x) std::fprintf(stderr, "[spumoni] %s\n", x);
#define OTHER_LOG(x) std::fprintf(stderr, "[helper-prog] %s\n", x);

/* DEBUG macros */
#define DEBUG 1
#define DBG_ONLY(x)  do { if (DEBUG) x; } while (0)

/* Function Declarations */
int spumoni_build_usage();
int build_main(int argc, char** argv);
int run_main(int argc, char** argv);
int spumoni_usage ();
int is_file(std::string path);
std::string execute_cmd(const char* cmd);

struct SpumoniHelperPrograms {
  /* Contains paths to run helper programs */
  std::string base_path = "";
  std::string parseNT_bin = "newscanNT.x";
  std::string parse_fasta_bin = "newscan.x";
  std::string parse_bin = "pscan.x";

public:
  void build_paths(std::string base) {
      /* Takes the base path, and combines it with names to generate executable paths */
      base_path.assign(base);
      parseNT_bin.assign(base_path + parseNT_bin);
      parse_fasta_bin.assign(base_path + parse_fasta_bin);
      parse_bin.assign(base_path + parse_bin);
  }

  void validate() const {
      /* Makes sure that each path for an executable is valid */
      bool invalid_path = !is_file(parseNT_bin) | !is_file(parse_fasta_bin) | !is_file(parse_bin);
      if (invalid_path) {THROW_EXCEPTION(std::runtime_error("One or more of helper program paths are invalid."));}
  }
};

struct SpumoniBuildOptions {
  std::string ref_file = "default";
  size_t wind = 10; // sliding window size
  size_t hash_mod = 100; // hash modulus
  size_t threads = 0; // number of threads
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

struct SpumoniRunOptions {

      bool use_min = false;
};

/* Additional Function Declarations */
void parse_build_options(int argc, char** argv, SpumoniBuildOptions* opts);


/* helper method definitions */

int spumoni_build_usage () {
    /* prints out the usage information for the spumoni build sub-command */
    std::fprintf(stderr, "spumoni build - builds the ms/pml index for a specified reference file.\n");
    std::fprintf(stderr, "Usage: spumoni build [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-10sprints this usage message\n", "-h");
    std::fprintf(stderr, "\t%-10spath to reference file to be indexed\n", "-r [FILE]");
    std::fprintf(stderr, "\t%-10sbuild an index that can be used to compute MSs\n", "-M");
    std::fprintf(stderr, "\t%-10sbuild an index that can be used to compute PMLs\n", "-P");
    std::fprintf(stderr, "\t%-10ssliding window size (default: 10)\n", "-w [arg]");
    std::fprintf(stderr, "\t%-10shash modulus value (default: 100)\n", "-p [arg]");
    std::fprintf(stderr, "\t%-10snumber of threads (default: 0)\n", "-t [arg]");
    std::fprintf(stderr, "\t%-10skeep the temporary files (default: false)\n", "-k");
    std::fprintf(stderr, "\t%-10suse when the reference file is a fasta file (default: false)\n\n", "-f");
    return 1;
}

void parse_build_options(int argc, char** argv, SpumoniBuildOptions* opts) {
    /* Parses the arguments for the build sub-command and returns a struct with arguments */
    for(int c;(c = getopt(argc, argv, "hr:MPw:p:t:kf")) >= 0;) { 
        switch(c) {
                    case 'h': spumoni_build_usage(); std::exit(1);
                    case 'r': opts->ref_file.assign(optarg); break;
                    case 'M': opts->ms_index = true; break;
                    case 'P': opts->pml_index = true; break;
                    case 'w': opts->wind = std::max(std::atoi(optarg), 10); break;
                    case 'p': opts->hash_mod = std::max(std::atoi(optarg), 1); break;
                    case 't': opts->threads = std::max(std::atoi(optarg), 0); break;
                    case 'k': opts->keep_files = true; break;
                    case 'f': opts->is_fasta = true; break;
                    default: spumoni_build_usage(); std::exit(1);
        }
    }
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    std::ifstream test_file(path.data());
    if (test_file.fail()) {return 0;}
    
    test_file.close();
    return 1;
}

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";
#ifdef _WIN32 || _WIN64
#define popen _popen
#define pclose _pclose
#endif
    FILE* pipe = popen(cmd, "r"); // Extract stderr as well
    if (!pipe) {THROW_EXCEPTION(std::runtime_error("popen() failed!"));}

    try {
        std::size_t bytes;
        while ((bytes = fread(buffer.data(), sizeof(char), sizeof(buffer), pipe))) {
            output += std::string(buffer.data(), bytes);
        }
    } catch (...) {
        pclose(pipe);
        THROW_EXCEPTION(std::runtime_error("Error occurred while reading popen() stream."));
    }
    pclose(pipe);
    return output;
}

#endif /* End of SPUMONI_MAIN_H */