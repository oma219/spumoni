import argparse
import os
import sys
import numpy as np
from scipy import stats

def compute_boundaries(read_len, batch_size):
    """ Compute the boundaries for independent regions that KS-test will be performed """
    pos_list = []
    curr_pos = 0
    while curr_pos + batch_size <= read_len:
        pos_list.append([curr_pos, curr_pos+batch_size])
        curr_pos += batch_size
    if curr_pos != read_len:
        pos_list.append([curr_pos, read_len])
    return pos_list

def perform_ks_test_on_reads(pos_data, null_data, pos_names, ks_stat_threshold, args):
    """ Go through each read, and compute KS-stat on each region """
    for read_num in range(len(pos_data)):
        ks_stat_list = []
        pos_list = compute_boundaries(len(pos_data[read_num]), args.region_size)
        
        for start_pos, end_pos in pos_list:
            ks_stat, p_value = stats.kstest(pos_data[read_num][start_pos:end_pos], null_data[read_num][start_pos:end_pos], N=100,alternative='less')
            ks_stat_list.append(ks_stat)

        decision_value = True
        num_under_ks_threshold = sum([1 if ks_stat < ks_stat_threshold else 0 for ks_stat in ks_stat_list])

        if num_under_ks_threshold/(len(ks_stat_list)+0.0) > 0.50:
            decision_value = False
        read_decision = "found" if decision_value else "not_found"

        for i, [start_pos, end_pos] in enumerate(pos_list):
            print(f"{pos_names[read_num]:<30}{start_pos:>15}{end_pos:>15}{ks_stat_list[i]:>10.3f}{read_decision:>15}")

def assert_reads_in_same_order(pos_names, null_names):
    if len(pos_names) != len(null_names) or len(pos_names) == 0:
        print_status("ERROR: Input files do not have the same number of reads, or there are no reads in one of the file.")
        exit(0)

    equivalent_reads = np.all(pos_names == null_names)
    if not equivalent_reads:
        print_status("ERROR: Reads in each input file, do not match.")
        exit(0)
    
def load_in_positive_and_null_data(pos_file, null_file):
    input_data = [[] for i in range(2)]
    input_names = [[] for i in range(2)]

    # Load in positive MSs/PMLs, and then null ones.
    for i,file_path in enumerate([pos_file, null_file]):
        ms_file = open(file_path,"r") 
        ms_file_lines = ms_file.readlines()

        ms_seqs = [x.split() for x in ms_file_lines if ">" not in x]
        ms_int_seqs = [[int(x) for x in y] for y in ms_seqs]
        input_data[i] = ms_int_seqs
        input_names[i] = [x.strip() for x in ms_file_lines if ">" in x]
        ms_file.close()

    return input_data[0], input_names[0], input_data[1], input_names[1]

def main(args):
    """ Runs the main code """
    pos_data, pos_names, null_data, null_names = load_in_positive_and_null_data(args.pos_data_file, args.null_data_file)
    assert_reads_in_same_order(pos_names, null_names)
    print_status("Loaded data correctly.")

    # Set threshold with default, and override with user's threshold if provided
    ks_stat_threshold = 0.10 if not args.use_ms else 0.25
    if args.ks_stat_threshold != -1:
        ks_stat_threshold = args.ks_stat_thresold
    print_status(f"Thresholds Used: ks_stat = {ks_stat_threshold}")
    print_status(f"Positive_File= {args.pos_data_file} Null_File= {args.null_data_file}")

    # Perform ks-test for each read and output results
    print_status("Begin processing the MSs or PMLs.")
    perform_ks_test_on_reads(pos_data, null_data, pos_names, ks_stat_threshold, args)
    print_status("Finished processing the MSs or PMLs successfully.")
    
    sys.stdout.close()
    sys.stderr.close()

def parse_arguments():
    """ Define arguments and parse arguments """
    parser = argparse.ArgumentParser(description="Load in PMLs or MSs from positive and null indexes, and report statistics.")
    parser.add_argument("-p", dest="pos_data_file", help="path to PMLS or MSs generated with respect to positive index.", required=True)
    parser.add_argument("-n", dest="null_data_file", help="path to PMLS or MSs generated with respect to null index.", required=True)
    parser.add_argument("--ms", action="store_true", dest="use_ms", help="use MSs instead of PMLs. (Default: Assumes input is PMLs)")
    parser.add_argument("-k", "--ks-stat", type=float, default=-1, dest="ks_stat_threshold", help="set the threshold for Kolmogorov-Smirnov Statistic (Default: 0.10 for PMLs, 0.25 for MSs)")
    parser.add_argument("-r", dest="region_size", type=int, default= 90, help="region size in bp where KS-test is performed. (Default: 90 bp)")
    parser.add_argument("-o", dest= "output_file", default= None, help="name of output file. (Default: analyzed data to stdout, log to stderr)")
    args = parser.parse_args()
    return args

def precheck_arguments(args):
    """ Make sure arguments are valid """
    if not os.path.isfile(args.pos_data_file) or not os.path.isfile(args.null_data_file):
        print_status("ERROR: Either positive or null data files cannot be found.")
        exit(0)

    if not args.use_ms:
        if not args.pos_data_file.endswith("pseudo_lengths") or not args.null_data_file.endswith("pseudo_lengths"):
            print_status("ERROR: Incorrect input filetype, given that PMLs are expected as input.")
            exit(0)
    else:
        if not args.pos_data_file.endswith("lengths") or not args.null_data_file.endswith("lengths"):
            print_status("ERROR: Incorrect input filetype, given that MSs are expected as input.")
            exit(0)
    
    if args.region_size < 10:
        print_status("ERROR: Region size must be >= 10.")
        exit(0)
        
    # Redirect stdout, stderr to user-specified file
    if args.output_file is not None:
        sys.stdout = open(args.output_file, "w")
        sys.stderr = sys.stdout

def print_status(msg):
    print(f"[spumoni_status] {msg}", file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    precheck_arguments(args)
    main(args)