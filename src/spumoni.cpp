/* 
 *  spumoni.cpp
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
  * Description: Main file for spumoni code, runs all the associated programs
  *              for building indexes and computing MS/PML. This file is based on
  *              pipeline written by Massimiliano Rossi for initial version of
  *              SPUMONI.
  *
  * Author: Omar Ahmed
  * Start Date: October 10, 2021
  */

#include <iostream>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <chrono>
#include <stdlib.h>
#include <spumoni_main.hpp>
#include <compute_ms_pml.hpp>
#include <doc_array.hpp>

int spumoni_run_usage () {
    /* prints out the usage information for the spumoni build sub-command */
    std::fprintf(stderr, "spumoni run - Uses a spumoni index to compute MS/PML of patterns w.r.t. a reference.\n");
    std::fprintf(stderr, "Usage: spumoni run [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-10sprints this usage message\n", "-h");
    std::fprintf(stderr, "\t%-10spath to reference file that has index built for it\n", "-r [FILE]");
    std::fprintf(stderr, "\t%-10spath to patterns file that will be used.\n", "-p [FILE]");
    std::fprintf(stderr, "\t%-10sUse index to compute MSs\n", "-M");
    std::fprintf(stderr, "\t%-10sUse index to compute PMLs\n", "-P");
    std::fprintf(stderr, "\t%-10spattern file is in fasta format (default: general text)\n", "-f");
    std::fprintf(stderr, "\t%-10suse document array to get assignments\n", "-d");
    std::fprintf(stderr, "\t%-10snumber of helper threads (default: 0)\n\n", "-t [arg]");
    return 1;
}

int spumoni_build_usage () {
    /* prints out the usage information for the spumoni build sub-command */
    std::fprintf(stderr, "spumoni build - builds the ms/pml index for a specified reference file.\n");
    std::fprintf(stderr, "Usage: spumoni build [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-10sprints this usage message\n", "-h");
    std::fprintf(stderr, "\t%-10spath to reference file to be indexed\n", "-r [FILE]");
    std::fprintf(stderr, "\t%-10sbuild directory for spumoni (if spumoni is not in current directory)\n", "-b [DIR]");
    std::fprintf(stderr, "\t%-10sbuild an index that can be used to compute MSs\n", "-M");
    std::fprintf(stderr, "\t%-10sbuild an index that can be used to compute PMLs\n", "-P");
    std::fprintf(stderr, "\t%-10ssliding window size (default: 10)\n", "-w [arg]");
    std::fprintf(stderr, "\t%-10shash modulus value (default: 100)\n", "-p [arg]");
    std::fprintf(stderr, "\t%-10snumber of helper threads (default: 0)\n", "-t [arg]");
    std::fprintf(stderr, "\t%-10skeep the temporary files (default: false)\n", "-k");
    std::fprintf(stderr, "\t%-10suse when the reference file is a fasta file (default: false)\n", "-f");
    std::fprintf(stderr, "\t%-10sbuild the document array (default: false)\n\n", "-d");
    return 1;
}

void parse_build_options(int argc, char** argv, SpumoniBuildOptions* opts) {
    /* Parses the arguments for the build sub-command and returns a struct with arguments */
    for(int c;(c = getopt(argc, argv, "hr:MPw:p:t:kfd")) >= 0;) { 
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
                    case 'd': opts->build_doc = true; break;
                    default: spumoni_build_usage(); std::exit(1);
        }
    }
}

void parse_run_options(int argc, char** argv, SpumoniRunOptions* opts) {
    /* Parses the arguments for the build sub-command and returns a struct with arguments */
    for(int c;(c = getopt(argc, argv, "hr:p:MPft:d")) >= 0;) { 
        switch(c) {
                    case 'h': spumoni_run_usage(); std::exit(1);
                    case 'r': opts->ref_file.assign(optarg); break;
                    case 'p': opts->pattern_file.assign(optarg); break;
                    case 'M': opts->ms_requested = true; break;
                    case 'P': opts->pml_requested = true; break;
                    case 'f': opts->query_fasta = true; break;
                    case 't': opts->threads = std::max(std::atoi(optarg), 1); break;
                    case 'd': opts->use_doc = true; break;
                    default: spumoni_run_usage(); std::exit(1);
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

int is_dir(std::string path) {
    /* Checks if the directory is a valid path */
    return std::filesystem::exists(path);
}

std::vector<std::string> split(std::string input, char delim) {
    /* Takes in a string, and splits it based on delimiters */
    std::vector<std::string> word_list;
    std::string curr_word = "";

    for (char ch: input) {
        if (ch == delim && curr_word.length()) {word_list.push_back(curr_word); curr_word = "";}
        else if (ch != delim) {curr_word += ch;}
    }
    if (curr_word.length()) {word_list.push_back(curr_word);}
    return word_list;
}

size_t get_avail_phy_mem() {
    /* Computes the available memory on UNIX machines */
    std::size_t pages = sysconf(_SC_PHYS_PAGES);
    std::size_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
inline size_t get_avail_phy_mem_win() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
}
#endif

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";
#if defined(_WIN32) || defined(_WIN64)
#define popen _popen
#define pclose _pclose
#endif
    std::string cmd_plus_stderr = std::string(cmd) + " 2>&1";
    FILE* pipe = popen(cmd_plus_stderr.data(), "r"); // Extract stderr as well
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

void run_build_grammar_cmds(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates the command-line for compressing the PFP dictonary */
#if defined(_WIN32) || defined(_WIN64)
#define get_avail_phy_mem get_avail_phy_mem_win
#endif
    auto avail_mem = get_avail_phy_mem()/1048576;
    DBG_ONLY("Available memory for RePair (MB): %d\n", avail_mem);

    // Generate and run command to compress the dictionary
    std::ostringstream command_stream;
    command_stream << helper_bins->compress_bin << " " << build_opts->ref_file;
    command_stream << " -w " << build_opts->wind << " -p " << build_opts->hash_mod;

    SPUMONI_LOG("Compressing the PFP dictionary for given reference ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto start = std::chrono::system_clock::now();
    auto output_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(output_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));

    command_stream.str(""); command_stream.clear();

    // Generate and run command to preprocess dictionary
    command_stream << helper_bins->preprocess_dict_bin << " " << build_opts->ref_file << ".dicz";

    SPUMONI_LOG("Preprocessing the PFP dictionary for given reference ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(output_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));

    command_stream.str(""); command_stream.clear();

    // Generate and run command to run RePair on the modified dictionary
    command_stream << helper_bins->repair_bin << " " << build_opts->ref_file << ".dicz.int ";
    command_stream << avail_mem;

    SPUMONI_LOG("Running RePair to generate grammar for modified dictionary ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(output_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));

    command_stream.str(""); command_stream.clear();

    // Generate and run command to run RePair on the parse of PFP
    command_stream << helper_bins->repair_bin << " " << build_opts->ref_file << ".parse ";
    command_stream << avail_mem;

    SPUMONI_LOG("Running RePair to generate grammar for parse ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(output_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));

    command_stream.str(""); command_stream.clear();

    // Generate and run command to postprocess the grammars
    command_stream << helper_bins->postprocess_gram_bin << " " << build_opts->ref_file;

    SPUMONI_LOG("Running post-processing of grammars ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(output_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));

    command_stream.str(""); command_stream.clear();

    // Remove the temporary files not needed.
    command_stream << "rm -f " << build_opts->ref_file << ".parse.C ";
    command_stream << build_opts->ref_file << ".parse.R " << build_opts->ref_file << ".dicz.int ";
    command_stream << build_opts->ref_file << ".dicz.int.C " << build_opts->ref_file << ".dicz.int.R";

    SPUMONI_LOG("Removing the temporary parse and dictionary files ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());
}

void run_build_slp_cmds(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs the command-line for building the SLP */
    std::ostringstream command_stream;
    command_stream << helper_bins->shaped_slp_bin << " -i " << build_opts->ref_file;
    command_stream << " -o " << build_opts->ref_file << ".slp -e SelfShapedSlp_SdSd_Sd -f Bigrepair";

    SPUMONI_LOG("Generating the SLP for the given reference ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto start = std::chrono::system_clock::now();
    auto output_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(output_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));
}

void run_build_parse_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs the command-line for executing the PFP of the reference */
    std::ostringstream command_stream;
    if (build_opts->threads > 0) {
        std::string curr_exe = "";
        if (build_opts->is_fasta) {curr_exe.assign(helper_bins->parse_fasta_bin);}
        else {curr_exe.assign(helper_bins->parse_bin);}

        command_stream << helper_bins->parse_fasta_bin << " ";
        command_stream << build_opts->ref_file << " ";
        command_stream << "-w " << build_opts->wind;
        command_stream << " -p " << build_opts->hash_mod;
        command_stream << " -t " << build_opts->threads;
    }
    else {
        std::string curr_exe = "";
        command_stream << helper_bins->parseNT_bin << " ";
        command_stream << build_opts->ref_file << " ";
        command_stream << "-w " << build_opts->wind;
        command_stream << " -p " << build_opts->hash_mod;
    }
    if (build_opts->is_fasta) {command_stream << " -f";}

    SPUMONI_LOG("Generating the PFP for given reference ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto start = std::chrono::system_clock::now();
    auto parse_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(parse_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));
}

size_t run_build_ms_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Runs the constructor for generating the final index for computing MS */
    SPUMONI_LOG("Building the index for computing MS ...");
    
    size_t num_runs = 0;
    auto start = std::chrono::system_clock::now();
    num_runs = build_spumoni_ms_main(build_opts->ref_file);

    TIME_LOG((std::chrono::system_clock::now() - start));
    return num_runs;
}

size_t run_build_pml_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Runs the constructor for generating the final index for computing PML */
    SPUMONI_LOG("Building the index for computing PML ...");

    size_t num_runs = 0;
    auto start = std::chrono::system_clock::now();
    num_runs = build_spumoni_main(build_opts->ref_file);

    TIME_LOG((std::chrono::system_clock::now() - start));
    return num_runs;
}

void rm_temp_build_files(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs commands to remove temporary files during build process */
    std::ostringstream command_stream;
    command_stream << "rm -f " << build_opts->ref_file << ".parse_old ";
    command_stream << build_opts->ref_file << ".last ";

    SPUMONI_LOG("Removing some additional temporary files from build process ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto log = execute_cmd(command_stream.str().c_str());
}

void run_build_thresholds_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins){
    /* Generates and runs the command to compute the thresholds based on PFP */
    auto parse_size = std::filesystem::file_size(std::filesystem::path((build_opts->ref_file + ".parse").data()))/4;
    auto dict_size = std::filesystem::file_size(std::filesystem::path((build_opts->ref_file + ".dict").data()));

    std::string curr_exe = "";
    if (parse_size >= std::pow(2, 31-1) || dict_size >= std::pow(2, 31-4)) {
        curr_exe.assign(helper_bins->pfp_thresholds64); 
    } else {curr_exe.assign(helper_bins->pfp_thresholds);}

    std::ostringstream command_stream;
    command_stream << curr_exe << " " << build_opts->ref_file << " ";
    command_stream << "-w " << build_opts->wind << " -r";

    SPUMONI_LOG("Generating the thresholds using the PFP of the given reference ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto start = std::chrono::system_clock::now();
    auto thresholds_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(thresholds_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));
}

int build_main(int argc, char** argv) {
    /* main method for the build sub-command */
    if (argc == 1) return spumoni_build_usage();
    
    // Grab the build options, and validate them
    SpumoniBuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    build_opts.validate();

    SpumoniHelperPrograms helper_bins;
    helper_bins.build_paths((std::string(std::getenv("SPUMONI_BUILD_DIR")) + "/bin/").data());
    helper_bins.validate();
    auto build_process_start = std::chrono::system_clock::now();

    // Performs the parsing of the reference and builds the thresholds based on the PFP
    run_build_parse_cmd(&build_opts, &helper_bins);
    run_build_thresholds_cmd(&build_opts, &helper_bins);

    // Generate grammar over reference if you would like to compute MS
    size_t num_runs = 0;
    if (build_opts.ms_index) {
        run_build_grammar_cmds(&build_opts, &helper_bins);
        run_build_slp_cmds(&build_opts, &helper_bins);
        num_runs = run_build_ms_cmd(&build_opts, &helper_bins);
    }

    // Build the PML index if asked for as well 
    if (build_opts.pml_index) {num_runs = run_build_pml_cmd(&build_opts, &helper_bins);}

    // Build the document array if asked for as well
    if (build_opts.build_doc) {
        SPUMONI_LOG("Building the Document Array");
        auto start = std::chrono::system_clock::now();
        DocumentArray doc_arr(build_opts.ref_file, num_runs);

        std::ofstream out_stream(build_opts.ref_file + ".doc");
        doc_arr.serialize(out_stream);
        out_stream.close();
        TIME_LOG((std::chrono::system_clock::now() - start));
    }
    
    rm_temp_build_files(&build_opts, &helper_bins);
    auto total_build_time = std::chrono::duration<double>((std::chrono::system_clock::now() - build_process_start));
    SPUMONI_LOG("TOTAL Elapsed Time for Build Process (s): %.3f", total_build_time);

    return 0;
}

int run_main(int argc, char** argv) {
    /* main method for the run sub-command */
    if (argc == 1) return spumoni_run_usage();

    // Grab the run options, and validate they are not missing/don't make sense 
    SpumoniRunOptions run_opts;
    parse_run_options(argc, argv, &run_opts);
    run_opts.populate_output_type();
    run_opts.validate();
    
    switch (run_opts.result_type) {
        case MS: run_spumoni_ms_main(&run_opts); break;
        case PML: run_spumoni_main(&run_opts); break;
        default: FATAL_WARNING("The output type (MS or PML) must be specified.");
    }
    return 0;
}

int spumoni_usage () {
    /* Prints usage information for SPUMONI overall (used if no sub-command is chosen) */
    std::fprintf(stderr, "SPUMONI has different sub-commands to run which can used as follows:\n");
    std::fprintf(stderr, "Usage: spumoni <command> [options]\n\n");

    std::fprintf(stderr, "Commands:\n");
    std::fprintf(stderr, "\tbuild\tbuilds the index needed to compute MS or PMLs for a specified reference.\n");
    std::fprintf(stderr, "\trun\tcomputes MSs or PMLs for patterns against already built SPUMONI index.\n\n");
    return 1;
}

int main(int argc, char **argv) {
    /* main method for spumoni */
    std::fprintf(stdout, "SPUMONI version: %s\n\n", SPUMONI_VERSION);

    if (argc > 1) {
        if (std::strcmp(argv[1], "build") == 0)
            return build_main(argc-1, argv+1);
        if (std::strcmp(argv[1], "run") == 0)
            return run_main(argc-1, argv+1);
    }
    return spumoni_usage();
}