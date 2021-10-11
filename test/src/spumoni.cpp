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
  *              for building indexes and computing MS/PML.
  *
  * Author: Omar Ahmed
  * Start Date: October 10, 2021
  */

#include <iostream>
#include <cstring>
#include <filesystem>
#include <spumoni_main.hpp>

//std::cout << std::filesystem::path(build_opts.ref_file.data()).filename() << std::endl;

std::string generate_build_parse_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates the command-line for executing the PFP of the reference */
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
    return command_stream.str();
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

    /* Perform the Parsing of Reference */
    auto parse_command = generate_build_parse_cmd(&build_opts, &helper_bins);
    SPUMONI_LOG(("Executing this command: " + parse_command).data());

    auto parse_log = execute_cmd(parse_command.c_str());
    OTHER_LOG(parse_log.data());

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
    std::fprintf(stderr, "SPUMONI version: %s\n\n", SPUMONI_VERSION);

    if (argc > 1) {
        if (std::strcmp(argv[1], "build") == 0)
            return build_main(argc-1, argv+1);
        if (std::strcmp(argv[1], "run") == 0)
            return run_main(argc-1, argv+1);
    }
    return spumoni_usage();
}