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
#include <refbuilder.hpp>
#include <encoder.h>
#include <emp_null_database.hpp>

/*
 * Section 1: 
 * Contains methods that output command-line options, and
 * parses those options
 */

int spumoni_run_usage () {
    /* prints out the usage information for the spumoni build sub-command */
    std::fprintf(stderr, "spumoni run - Uses a spumoni index to compute MS/PML of patterns w.r.t. a reference.\n");
    std::fprintf(stderr, "Usage: spumoni run [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t*** GENERAL OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10sprints this usage message\n", "-h");
    std::fprintf(stderr, "\t%-10snumber of helper threads (default: 1)\n\n", "-t [arg]");

    std::fprintf(stderr, "\t*** INPUT/OUTPUT OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10spath to reference file that has index built for it\n", "-r [FILE]");
    std::fprintf(stderr, "\t%-10spath to patterns file that will be used.\n", "-p [FILE]");
    std::fprintf(stderr, "\t%-10suse index to compute MSs\n", "-M");
    std::fprintf(stderr, "\t%-10suse index to compute PMLs\n", "-P");
    //std::fprintf(stderr, "\t%-10spattern file is in fasta format (default: general text)\n", "-f");
    std::fprintf(stderr, "\t%-10suse document array to get assignments\n", "-d");
    std::fprintf(stderr, "\t%-10swrite out the classifications in a report file\n\n", "-c");

    std::fprintf(stderr, "\t*** MINIMIZER OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10sturn off minimizer digestion of reads (default: on)\n", "-n");
    std::fprintf(stderr, "\t%-10suse alphabet-promoted minimizers\n", "-m");
    std::fprintf(stderr, "\t%-10suse DNA-letter based minimizers\n", "-a");
    std::fprintf(stderr, "\t%-10ssmall window size (k) for finding minimizers (default: 4)\n", "-K");
    std::fprintf(stderr, "\t%-10slarge window size (w) for finding minimizers (default: 12)\n\n", "-W");

    return 0;
}

int spumoni_build_usage () {
    /* prints out the usage information for the spumoni build sub-command */
    std::fprintf(stderr, "spumoni build - builds the ms/pml index for a specified reference file.\n");
    std::fprintf(stderr, "Usage: spumoni build [options]\n\n");

    std::fprintf(stderr, "Options:\n");

    std::fprintf(stderr, "\t*** GENERAL OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10sprints this usage message\n", "-h");
    std::fprintf(stderr, "\t%-10sturn on verbose logging\n\n", "-v");

    std::fprintf(stderr, "\t*** INPUT DATA OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10spath to reference file to be indexed\n", "-r [FILE]");
    std::fprintf(stderr, "\t%-10sfile with a list of files to index\n", "-i [FILE]");
    std::fprintf(stderr, "\t%-10sbuild directory for index(es) (if using -i option)\n\n", "-b [DIR]");

    std::fprintf(stderr, "\t*** MINIMIZER OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10sturn off minimizer digestion of sequence (default: on)\n", "-n");
    std::fprintf(stderr, "\t%-10suse alphabet-promoted minimizers\n", "-m");
    std::fprintf(stderr, "\t%-10suse DNA-letter based minimizers\n", "-t");
    std::fprintf(stderr, "\t%-10ssmall window size (k) for finding minimizers (default: 4)\n", "-K");
    std::fprintf(stderr, "\t%-10slarge window size (w) for finding minimizers (default: 12)\n\n", "-W");

    std::fprintf(stderr, "\t*** INDEX FILE(S) OPTIONS ***\n\n");
    std::fprintf(stderr, "\t%-10sbuild an index that can be used to compute MSs\n", "-M");
    std::fprintf(stderr, "\t%-10sbuild an index that can be used to compute PMLs\n", "-P");
    std::fprintf(stderr, "\t%-10skeep the temporary files (default: false)\n", "-k");
    std::fprintf(stderr, "\t%-10sbuild the document array (default: false)\n\n", "-d");   

    //std::fprintf(stderr, "\t%-10ssliding window size (default: 10)\n", "-w [arg]");
    //std::fprintf(stderr, "\t%-10shash modulus value (default: 100)\n", "-p [arg]");
    //std::fprintf(stderr, "\t%-10snumber of helper threads (default: 0)\n", "-t [arg]");
    //std::fprintf(stderr, "\t%-10suse when the reference file is a fasta file (default: true)\n", "-f");

    return 0;
}

void parse_build_options(int argc, char** argv, SpumoniBuildOptions* opts) {
    /* Parses the arguments for the build sub-command and returns a struct with arguments */
    for(int c;(c = getopt(argc, argv, "hr:MPw:p:kdi:b:nvmtK:W:")) >= 0;) { 
        switch(c) {
                    case 'h': spumoni_build_usage(); std::exit(1);
                    case 'r': opts->ref_file.assign(optarg); break;
                    case 'i': opts->input_list.assign(optarg); break;
                    case 'b': opts->output_dir.assign(optarg); break;
                    case 'M': opts->ms_index = true; break;
                    case 'P': opts->pml_index = true; break;
                    case 'v': opts->verbose = true; break;
                    case 'n': opts->use_minimizers = false; opts->is_fasta = true; break;
                    case 'm': opts->use_promotions = true; break;
                    case 't': opts->use_dna_letters = true; opts->is_fasta = true; break;
                    case 'K': opts->k = std::max(std::atoi(optarg), 1); break;
                    case 'W': opts->w = std::max(std::atoi(optarg), 1); break;
                    case 'w': opts->wind = std::max(std::atoi(optarg), 10); break;
                    case 'p': opts->hash_mod = std::max(std::atoi(optarg), 1); break;
                    //case 't': opts->threads = std::max(std::atoi(optarg), 1); break;
                    case 'k': opts->keep_files = true; break;
                    //case 'f': opts->is_fasta = true; break;
                    case 'd': opts->build_doc = true; break;
                    default: spumoni_build_usage(); std::exit(1);
        }
    }
}

void parse_run_options(int argc, char** argv, SpumoniRunOptions* opts) {
    /* Parses the arguments for the build sub-command and returns a struct with arguments */
    for(int c;(c = getopt(argc, argv, "hr:p:MPt:dcnmaK:W:")) >= 0;) { 
        switch(c) {
                    case 'h': spumoni_run_usage(); std::exit(1);
                    case 'r': opts->ref_file.assign(optarg); break;
                    case 'p': opts->pattern_file.assign(optarg); break;
                    case 'M': opts->ms_requested = true; break;
                    case 'P': opts->pml_requested = true; break;
                    case 'c': opts->write_report = true; break;   
                    case 'm': opts->use_promotions = true; break;
                    case 'a': opts->use_dna_letters = true; break;
                    case 'n': opts->min_digest = false; break;
                    case 'K': opts->k = std::max(std::atoi(optarg), 1); break;
                    case 'W': opts->w = std::max(std::atoi(optarg), 1); break;
                    //case 'f': opts->query_fasta = true; break;
                    case 't': opts->threads = std::max(std::atoi(optarg), 1); break;
                    case 'd': opts->use_doc = true; break;
                    default: spumoni_run_usage(); std::exit(1);
        }
    }
}

/*
 * Section 2: 
 * Contains general helper methods that are used by different
 * methods across the SPUMONI repo.
 */

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

bool is_integer(const std::string& str) {
    /* Checks if string passed is an integer */
    std::string::const_iterator iter = str.begin();
    while (iter != str.end() && std::isdigit(*iter)) {++iter;}
    return !str.empty() && iter == str.end();
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

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
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

std::string perform_minimizer_digestion(std::string input_query, size_t k, size_t w) {
    /* Performs minimizer digestion using alphabet promotion, and returns concatenated minimizers */
    bns::RollingHasher<uint8_t> rh(k, false, bns::DNA, w);
    bool hp_compress = true;

    std::string mseq = "";
    std::vector<uint8_t> mseq_vec;

    // Define lambda function to take in a DNA sequence and return minimizer sequence
    auto get_minimizer_seq = [&](std::string seq, size_t seq_length) {
                mseq = ""; mseq_vec.clear();
                rh.for_each_uncanon([&](auto x) { // Extracts all minimizers and stores in mseq
                if(mseq_vec.empty() || !hp_compress || mseq_vec.back() != x) {
                    // Important Note - I push onto the vector first because it is raw hash value
                    // which we will check in the if statement. However, the real value we use in the
                    // sequence could be changed to account for PFP's reserved characters
                    mseq_vec.push_back(x);
                    x = (x > 2) ? x : (x + 3); // Reserves 0,1,2 for PFP
                    mseq += x;
                }
            }, seq.data(), seq_length);
            return mseq;
    };

    return get_minimizer_seq(input_query, input_query.length());
}

std::string perform_dna_minimizer_digestion(std::string input_query, size_t k, size_t w) {
    /* Generates string of concatenated minimizers in DNA alpahbet for input string, and returns it */
    std::vector<uint16_t> sp_vec;
    bns::Spacer sp(k, w, sp_vec);
    bns::Encoder<bns::score::Lex, uint64_t> enc(sp, false);

    bool hp_compress = true;
    std::string mseq = "";
    std::vector<uint8_t> mseq_vec;
    
    auto get_minimizer_seq = [&](std::string seq, size_t seq_length) {
            mseq = ""; mseq_vec.clear();
            enc.for_each([&](auto x) { 
                if(mseq_vec.empty() || !hp_compress || mseq_vec.back() != x) {
                    mseq_vec.push_back(x); 
                    mseq += sp.to_string(x);
                }
            }, seq.data(), seq_length);
            return mseq;
    };
    return get_minimizer_seq(input_query, input_query.length());
}


/*
 * Section 3: 
 * Contains methods that generate shell commands that are run
 * in order to build different components needed for the
 * r-index data-structure.
 */

void run_build_grammar_cmds(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates the command-line for compressing the PFP dictonary */
#if defined(_WIN32) || defined(_WIN64)
#define get_avail_phy_mem get_avail_phy_mem_win
#endif
    auto avail_mem = get_avail_phy_mem()/1048576;
    DBG_ONLY("Available memory for RePair (MB): %d\n", avail_mem);

    // Generate and run command to compress the dictionary
    std::ostringstream command_stream;
    command_stream << helper_bins->compress_bin << " -i " << build_opts->ref_file;
    command_stream << " -w " << build_opts->wind << " -p " << build_opts->hash_mod;
    
    LOG(build_opts->verbose, "build_grammar", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_grammar", "compressing the PFP dictionary");

    auto start = std::chrono::system_clock::now();
    auto output_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(output_log.data());

    command_stream.str(""); command_stream.clear();

    // Generate and run command to preprocess dictionary
    command_stream << helper_bins->preprocess_dict_bin << " " << build_opts->ref_file << ".dicz";

    LOG(build_opts->verbose, "build_grammar", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_grammar", "preprocessing the PFP dictionary");

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(output_log.data());

    command_stream.str(""); command_stream.clear();

    // Generate and run command to run RePair on the modified dictionary
    command_stream << helper_bins->repair_bin << " " << build_opts->ref_file << ".dicz.int ";
    command_stream << avail_mem;

    LOG(build_opts->verbose, "build_grammar", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_grammar", "running RePair to generate grammar for modified dictionary");

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(output_log.data());

    command_stream.str(""); command_stream.clear();

    // Generate and run command to run RePair on the parse of PFP
    command_stream << helper_bins->repair_bin << " " << build_opts->ref_file << ".parse ";
    command_stream << avail_mem;
    
    LOG(build_opts->verbose, "build_grammar", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_grammar", "running RePair to generate grammar for parse");

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(output_log.data());

    command_stream.str(""); command_stream.clear();

    // Generate and run command to postprocess the grammars
    command_stream << helper_bins->postprocess_gram_bin << " " << build_opts->ref_file;

    LOG(build_opts->verbose, "build_grammar", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_grammar", "running post-processing of grammars");

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(output_log.data());

    command_stream.str(""); command_stream.clear();

    // Remove the temporary files not needed.
    command_stream << "rm -f " << build_opts->ref_file << ".parse.C ";
    command_stream << build_opts->ref_file << ".parse.R " << build_opts->ref_file << ".dicz.int ";
    command_stream << build_opts->ref_file << ".dicz.int.C " << build_opts->ref_file << ".dicz.int.R";

    LOG(build_opts->verbose, "build_grammar", ("Executing this command: " + command_stream.str()).data());
    LOG(build_opts->verbose, "build_grammar", "removing the temporary parse and dictionary files");

    start = std::chrono::system_clock::now();
    output_log = execute_cmd(command_stream.str().c_str());
    //TIME_LOG((std::chrono::system_clock::now() - start));
}

void run_build_slp_cmds(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs the command-line for building the SLP */
    std::ostringstream command_stream;
    command_stream << helper_bins->shaped_slp_bin << " -i " << build_opts->ref_file;
    command_stream << " -o " << build_opts->ref_file << ".slp -e SelfShapedSlp_SdSd_Sd -f Bigrepair";

    LOG(build_opts->verbose, "build_slp", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_slp", "generating the SLP for the given reference");

    auto start = std::chrono::system_clock::now();
    auto output_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(output_log.data());
}

void run_build_parse_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs the command-line for executing the PFP of the reference */
    std::ostringstream command_stream;
    if (build_opts->threads > 0) {
        std::string curr_exe = "";
        if (build_opts->is_fasta) {curr_exe.assign(helper_bins->parse_fasta_bin);}
        else {curr_exe.assign(helper_bins->parse_bin);}

        command_stream << curr_exe << " -i ";
        command_stream << build_opts->ref_file << " ";
        command_stream << "-w " << build_opts->wind;
        command_stream << " -p " << build_opts->hash_mod;
        command_stream << " -t " << build_opts->threads;
    }
    else {
        std::string curr_exe = "";
        command_stream << helper_bins->parseNT_bin << " -i ";
        command_stream << build_opts->ref_file << " ";
        command_stream << "-w " << build_opts->wind;
        command_stream << " -p " << build_opts->hash_mod;
    }
    if (build_opts->is_fasta) {command_stream << " -f";}

    LOG(build_opts->verbose, "build_parse", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_parse", "generating the prefix-free parse for given reference");

    auto start = std::chrono::system_clock::now();
    auto parse_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(parse_log.data());
}

size_t run_build_ms_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Runs the constructor for generating the final index for computing MS */
    STATUS_LOG("build_ms", "building the index for computing MS");

    size_t length = 0, num_runs = 0;
    auto start = std::chrono::system_clock::now();  
    std::tie(length, num_runs) = build_spumoni_ms_main(build_opts->ref_file);
    
    double average_run_size = (length + 0.0)/num_runs;
    DONE_LOG((std::chrono::system_clock::now() - start));
    FORCE_LOG("build_ms", "bwt statistics: r = %d, n = %d, n/r = %.3f", num_runs, length, average_run_size);
    return num_runs;
}

size_t run_build_pml_cmd(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Runs the constructor for generating the final index for computing PML */
    STATUS_LOG("build_pml", "building the index for computing PML");

    size_t length = 0, num_runs = 0;
    auto start = std::chrono::system_clock::now();
    std::tie(length, num_runs) = build_spumoni_main(build_opts->ref_file);

    double average_run_size = (length + 0.0)/num_runs;
    DONE_LOG((std::chrono::system_clock::now() - start));
    FORCE_LOG("build_pml", "bwt statistics: r = %d, n = %d, n/r = %.3f", num_runs, length, average_run_size);
    return num_runs;
}

void rm_temp_build_files(SpumoniBuildOptions* build_opts, SpumoniHelperPrograms* helper_bins) {
    /* Generates and runs commands to remove temporary files during build process */

    const char* temp_build_files[15] = {".bwt.heads", ".bwt.len", ".R", ".C", ".dict", ".dicz", ".parse_old", ".last",
                                        ".dicz.len", ".ssa", ".esa", ".occ", ".parse", ".thr", ".thr_pos"};
    size_t num_temp_build_files = 15;

    // Build the remove command
    std::ostringstream command_stream;
    command_stream << "rm -f "; 
    for (size_t i = 0; i < num_temp_build_files; i++) {
        command_stream << build_opts->ref_file << temp_build_files[i] << " ";
    }
    command_stream << "rs_temp_output ";

    // Execute the removal of temporary files
    LOG(build_opts->verbose, "rm_temp", ("Executing this command: " + command_stream.str()).data());
    FORCE_LOG("rm_temp", "removing temporary files from build process");
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
    command_stream << curr_exe << " -i " << build_opts->ref_file << " ";
    command_stream << "-w " << build_opts->wind << " -r";

    LOG(build_opts->verbose, "build_thr", ("Executing this command: " + command_stream.str()).data());
    STATUS_LOG("build_thr", "generating the thresholds data-structure");

    auto start = std::chrono::system_clock::now();
    auto thresholds_log = execute_cmd(command_stream.str().c_str());
    DONE_LOG((std::chrono::system_clock::now() - start));
    OTHER_LOG(thresholds_log.data());
}

/*
 * Section 4: 
 * Contains the "main" methods of SPUMONI that ultimately call
 * other methods to either build indexes or run queries.
 */

int build_main(int argc, char** argv) {
    /* main method for the build sub-command */
    if (argc == 1) return spumoni_build_usage();
    
    // Grab the build options, and validate them
    SpumoniBuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    build_opts.validate();

    SpumoniHelperPrograms helper_bins;
    if (!std::getenv("SPUMONI_BUILD_DIR")) {FATAL_ERROR("Need to set SPUMONI_BUILD_DIR environment variable.");}

    helper_bins.build_paths((std::string(std::getenv("SPUMONI_BUILD_DIR")) + "/bin/").data());
    helper_bins.validate();

    // Variables needed for identifying required build files
    const char* temp_build_files[15] = {".bwt.heads", ".bwt.len", ".R", ".C", ".dict", ".dicz", ".parse_old", ".last",
                                        ".dicz.len", ".ssa", ".esa", ".occ", ".parse", ".thr", ".thr_pos"};
    size_t num_temp_build_files = 15;

    // Check all needed files are already present
    std::filesystem::path p1 = build_opts.ref_file;
    std::string build_ref_file = (build_opts.use_promotions) ? "spumoni_full_ref.bin" : "spumoni_full_ref.fa";
    std::string null_read_file = "";
    
    char ch;
    if (build_opts.input_list.length()){ // using a file-list
        std::string build_dir(build_opts.output_dir);
        if ((ch = build_dir.back()) != '/') {build_dir += "/";}
        null_read_file = build_dir + "spumoni_null_reads.fa";
        build_ref_file = build_dir + build_ref_file;
    } else if (p1.parent_path().string().length()) { // using single file not in current dir
        null_read_file = p1.parent_path().string() + "/spumoni_null_reads.fa";
        build_ref_file = p1.parent_path().string() + "/" + build_ref_file;
    } else { // using single file in current directory
        null_read_file = "spumoni_null_reads.fa";
    }

    bool quick_build = is_file(null_read_file);
    for (size_t i = 0; i < num_temp_build_files && quick_build; i++) {
        if (!is_file(build_ref_file + temp_build_files[i])) {quick_build = false;}
    }
    if (quick_build) {
        FORCE_LOG("build_main", "quick build is activated.");
        build_opts.ref_file = build_ref_file;
    }

    // Start of build process ...
    auto total_build_process_start = std::chrono::system_clock::now();
    auto task_start = std::chrono::system_clock::now();

    if (!quick_build) {
        // Print out information describing the input files...
        if (!build_opts.input_list.length())
            FORCE_LOG("build_main", "input: single reference file (%s)\n", build_opts.ref_file.data());
        else {
            FORCE_LOG("build_main", "input: list of files (%s)", build_opts.input_list.data());
            FORCE_LOG("build_main", "build directory: %s\n", build_opts.output_dir.data());
        }

        if (build_opts.use_promotions) 
            STATUS_LOG("build_main", "reference file is being generated (spumoni_full_ref.bin)");
        else 
            STATUS_LOG("build_main", "reference file is being generated (spumoni_full_ref.fa)");
        task_start = std::chrono::system_clock::now();

        // Perform needed operations to input file(s) prior to building index
        if (build_opts.input_list.length()){         
            RefBuilder refbuild (build_opts.ref_file.data(), build_opts.input_list.data(), build_opts.output_dir.data(),
                                build_opts.build_doc, build_opts.input_list.length(), build_opts.use_minimizers,
                                build_opts.use_promotions, build_opts.use_dna_letters,
                                build_opts.k, build_opts.w);
            build_opts.ref_file = refbuild.get_ref_path();
            null_read_file = refbuild.get_null_readfile();
        } else {
            null_read_file = RefBuilder::parse_null_reads(build_opts.ref_file.data());
            build_opts.ref_file = RefBuilder::build_reference(build_opts.ref_file.data(), build_opts.use_promotions,
                                                                build_opts.use_dna_letters, build_opts.k,
                                                                build_opts.w);
        }
        DONE_LOG((std::chrono::system_clock::now() - task_start));

        // Performs the parsing of the reference and builds the thresholds based on the PFP
        run_build_parse_cmd(&build_opts, &helper_bins);
        run_build_thresholds_cmd(&build_opts, &helper_bins); std::cout << std::endl;
    }

    // Build the MS index, as well as grammar needed
    size_t num_runs = 0;
    if (build_opts.ms_index) {
        run_build_grammar_cmds(&build_opts, &helper_bins);
        run_build_slp_cmds(&build_opts, &helper_bins); std::cout << std::endl;
        num_runs = run_build_ms_cmd(&build_opts, &helper_bins);

        // Build the null database for MSs
        STATUS_LOG("build_main", "building the empirical null statistic database for MS");
        task_start = std::chrono::system_clock::now();
        EmpNullDatabase null_db(build_opts.ref_file.data(), null_read_file.data(), build_opts.use_minimizers, MS,
                                build_opts.use_promotions, build_opts.use_dna_letters, build_opts.k, build_opts.w);

        std::string output_nulldb_name = build_opts.ref_file + ".msnulldb";
        std::ofstream out_stream(output_nulldb_name);
        null_db.serialize(out_stream);
        out_stream.close();
        DONE_LOG((std::chrono::system_clock::now() - task_start));
        std::cout << std::endl;
    } 

    // Build the PML index if asked for as well 
    if (build_opts.pml_index) {
        num_runs = run_build_pml_cmd(&build_opts, &helper_bins);

        // Build the null database for PMLs
        STATUS_LOG("build_main", "building the empirical null statistic database for PML" );
        task_start = std::chrono::system_clock::now();
        EmpNullDatabase null_db(build_opts.ref_file.data(), null_read_file.data(), build_opts.use_minimizers, PML,
                                build_opts.use_promotions, build_opts.use_dna_letters, build_opts.k, build_opts.w);

        std::string output_nulldb_name = build_opts.ref_file + ".pmlnulldb";
        std::ofstream out_stream(output_nulldb_name);
        null_db.serialize(out_stream);
        out_stream.close();
        DONE_LOG((std::chrono::system_clock::now() - task_start));
        std::cout << std::endl;
    }

    // Build the document array if asked for as well
    if (build_opts.build_doc) {
        STATUS_LOG("build_main", "building the document array");
        task_start = std::chrono::system_clock::now();
        DocumentArray doc_arr(build_opts.ref_file, num_runs);

        std::ofstream out_stream(build_opts.ref_file + ".doc");
        doc_arr.serialize(out_stream);
        out_stream.close();
        DONE_LOG((std::chrono::system_clock::now() - task_start));
    }

    if (!build_opts.keep_files) {rm_temp_build_files(&build_opts, &helper_bins); std::cout << "\n";}
    auto total_build_time = std::chrono::duration<double>((std::chrono::system_clock::now() - total_build_process_start));

    std::string final_index_files = "";
    if (build_opts.ms_index && !build_opts.pml_index) 
        final_index_files = "index files are saved in the *.ms, *.msnulldb, and *.slp files.";
    else if (!build_opts.ms_index && build_opts.pml_index) 
        final_index_files = "index files are saved in the *.spumoni, and *.pmlnulldb files.";
    else 
        final_index_files = "index files are saved in the  *.ms, *.spumoni, *.msnulldb, *.pmlnulldb and *.slp files.";

    FORCE_LOG("build_main", "total elapsed time for build process (s): %.3f", total_build_time);
    FORCE_LOG("build_main", final_index_files.data());
    return 0;
}

int run_main(int argc, char** argv) {
    /* main method for the run sub-command */
    if (argc == 1) return spumoni_run_usage();

    // Grab the run options, and validate they are not missing/don't make sense 
    SpumoniRunOptions run_opts;
    parse_run_options(argc, argv, &run_opts);
    run_opts.populate_types();
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