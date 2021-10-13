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
  * Description: Main file for spumoni code, runs all the associated programs
  *              for building indexes and computing MS/PML. This file is based
  *              pipeline written by Massimiliano Rossi for initial version 
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
#include <spumoni_main.hpp>

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

void run_build_ms_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs the command-line for generating the final index for computing MS */
    std::ostringstream command_stream;
    command_stream << helper_bins->ms_build << " " << build_opts->ref_file;

    SPUMONI_LOG("Building the index for computing MS ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto start = std::chrono::system_clock::now();
    auto parse_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(parse_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));
}

void run_build_pml_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs the command-line for generating the final index for computing PML */
    std::ostringstream command_stream;
    command_stream << helper_bins->pml_build << " " << build_opts->ref_file;

    SPUMONI_LOG("Building the index for computing PML ...");
    SPUMONI_LOG(("Executing this command: " + command_stream.str()).data());

    auto start = std::chrono::system_clock::now();
    auto parse_log = execute_cmd(command_stream.str().c_str());
    OTHER_LOG(parse_log.data());
    TIME_LOG((std::chrono::system_clock::now() - start));
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
    
    /* Grab the build options, and validate they are not missing/don't make sense */
    SpumoniBuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    build_opts.validate();

    SpumoniHelperPrograms helper_bins;
    helper_bins.build_paths("./bin/"); // TODO: add option to change helper directory
    helper_bins.validate();

    /* Performs the parsing of the reference and builds the thresholds based on the PFP */
    run_build_parse_cmd(&build_opts, &helper_bins);
    run_build_thresholds_cmd(&build_opts, &helper_bins);

    /* Generate grammar over reference if you would like to compute MS, along with the final data-structure */
    if (build_opts.ms_index) {
        run_build_grammar_cmds(&build_opts, &helper_bins);
        run_build_slp_cmds(&build_opts, &helper_bins);
        run_build_ms_cmd(&build_opts, &helper_bins);
    }

    /* Build the PML index if asked for as well */
    if (build_opts.pml_index) {run_build_pml_cmd(&build_opts, &helper_bins);}
    return 1;
}

int run_main(int argc, char** argv) {

    NOT_IMPL("still working on this ...");
    return 1;
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