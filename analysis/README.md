## Detailed Overview of Commands for SPUMONI

This section will just use the example data in `data/` as an example scenario for how to use SPUMONI for binary classification of reads. The example data contains an example mock community reference, as well as the "null" version of it which is just the reverse of the sequence. 

The example commands below show the process of building the indexes for the positive and null database, followed by computing both MSs and PMLs, and then analyzing both of them to determine the read classifications.

```sh
// Builds both SPUMONI and MONI-ms indexes for both positive and null indexes 
// Takes about ~3 minutes for each index
./spumoni build -r ../data/example_positive_index/mock_comm_positive.fasta -f -P -M
./spumoni build -r ../data/example_null_index/mock_comm_null.fasta -f -P -M

// Computes pseudo matching lengths against both positive and null indexes
./spumoni run -r ../data/example_positive_index/mock_comm_positive.fasta \
              -p ../data/example_reads/all_reads.fa -f -P
./spumoni run -i ../data/example_null_index/mock_comm_null.fasta \
              -p ../data/example_reads/all_reads.fa -f -P

// Computes matching statistics lengths against both positive and null indexes
./spumoni run -r ../data/example_positive_index/mock_comm_positive.fasta \
              -p ../data/example_reads/all_reads.fa -f -M
./spumoni run -i ../data/example_null_index/mock_comm_null.fasta \
              -p ../data/example_reads/all_reads.fa -f -M

// Analyze the pseudo matching lengths (move to analysis folder first)
cd ../analysis
python3 analyze_pml.py -p ../data/example_reads/all_reads.fa_mock_comm_positive.fasta.pseudo_lengths \
                       -n ../data/example_reads/all_reads.fa_mock_comm_null.fasta.pseudo_lengths \
                        > pml_analysis.txt

// Analyze the matching statistic lengths
python3 analyze_pml.py --ms -p ../data/example_reads/all_reads.fa_mock_comm_positive.fasta.lengths \
                            -n ../data/example_reads/all_reads.fa_mock_comm_positive.fasta.lengths \
                             > ms_analysis.txt
```







## Description of Output from analyze_pml.py Script

The `analyze_pml.py` script takes in matching statistics or pseudo matching lengths with respect to both a positive and null index and it will output data
in the form of various rows where each row corresponds to an indepedent region within a particular read where the Kolmogorov-Smirnov test (KS-test) was applied. 
Each column of the data will be explained below:

- Column 1: Name of the Read
- Column 2: Start Position of Region where KS-test was Applied (0-based, and inclusive)
- Column 3: End Position of Region where KS-test was Applied (0-based, and exclusive)
- Column 4: KS-statistic for this Region (Range is between 0 to 1)
- Column 5: Overall Decision for Current Read 
    - "found" means a majority of KS-stats were significant
    - "not_found" means a majority of the KS-stats were not significant
