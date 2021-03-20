# SPUMONI
<!--- ```console
      _____ _____  _    _ __  __  ____  _   _ _____ 
     / ____|  __ \| |  | |  \/  |/ __ \| \ | |_   _|
    | (___ | |__) | |  | | \  / | |  | |  \| | | |  
     \___ \|  ___/| |  | | |\/| | |  | | . ` | | |  
     ____) | |    | |__| | |  | | |__| | |\  |_| |_ 
    |_____/|_|     \____/|_|  |_|\____/|_| \_|_____|
                                            ver 0.1.0
```
-->

<img src="../assets/spumoni_pic.jpeg" width="500">

Pan-genomic Matching Statistics for Targeted Nanopore Sequencing

Based on `MONI`: A MEM-finder with Multi-Genome References.

`MONI` index uses the prefix-free parsing of the text [2][3] to build the Burrows-Wheeler Transform (BWT) of the reference genomes, the suffix array (SA) samples at the beginning and at the end of each run of the BWT, and the threshold positions of [1]. 

*Current Version:* 0.1.0

# Usage

### Construction of the Index:
```
usage: moni build [-h] -r REFERENCE [-w WSIZE] [-p MOD] [-t THREADS] [-k] [-v]
                  [-f] [--moni-ms] [--spumoni]
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        reference file name (default: None)
  -w WSIZE, --wsize WSIZE
                        sliding window size (default: 10)
  -p MOD, --mod MOD     hash modulus (default: 100)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 0)
  -k                    keep temporary files (default: False)
  -v                    verbose (default: False)
  -f                    read fasta (default: False)
  --moni-ms             build moni index for matching statistics only
                        (default: False)
  --spumoni             build spumoni index (default: True)

```
The default index is the `spumoni` index. If you want to build the `moni-ms` index you can use the `--moni-ms` flag. If you want to build both you can use `--moni-ms --spumoni`. Building `moni-ms` and `spumoni` together is faster than building them separately.

### Compute the matching statistics (MSs) with MONI-ms:
```
usage: moni ms [-h] -i INDEX -p PATTERN [-o OUTPUT] [-t THREADS]
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index base name (default: None)
  -p PATTERN, --pattern PATTERN
                        the input query (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: .)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 1)
```

### Compute the pseudo matching lengths (PMLs) with SPUMONI:
```
usage: moni pseudo-ms [-h] -i INDEX -p PATTERN [-o OUTPUT] [-t THREADS]
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index base name (default: None)
  -p PATTERN, --pattern PATTERN
                        the input query (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: .)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 1)
```
### Analyzing either MSs/PMLs from MONI-ms or SPUMONI

```
usage: python3 analyze_pml.py [-h] -p POS_DATA_FILE -n NULL_DATA_FILE [--ms] \
                              [-k KS_STAT_THRESHOLD] [-r REGION_SIZE] [-o OUTPUT_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -p POS_DATA_FILE      
            path to PMLS or MSs generated with respect to positive index.
  -n NULL_DATA_FILE     
            path to PMLS or MSs generated with respect to null index.
  --ms                  
            use MSs instead of PMLs. (Default: Assumes input is PMLs)
  -k KS_STAT_THRESHOLD, --ks-stat KS_STAT_THRESHOLD 
            set the threshold for Kolmogorov-Smirnov Statistic. (Default: 0.10 for PMLs, 0.25 for MSs)
  -r REGION_SIZE        
            region size in bp where KS-test is performed. (Default: 90 bp)
  -o OUTPUT_FILE        
            name of output file. (Default: analyzed data to stdout, log to stderr)
```
The script, `analyze_pml.py`, will output the analyzed data to `stdout`. To further understand the meaning of each column in the output, check out the following [README](https://github.com/oma219/spumoni/blob/main/analysis/README.md).

# Example
### Download and Compile

```console
git clone https://github.com/oma219/spumoni
cd spumoni
mkdir build && cd build
cmake ..
make
```

### Run on Example Data

```console
// Builds both SPUMONI and MONI-ms indexes for both positive and null indexes (each takes about ~3 minutes)
python3 moni build -r ../data/example_positive_index/mock_comm_positive.fasta -f --spumoni --moni-ms
python3 moni build -r ../data/example_null_index/mock_comm_null.fasta -f --spumoni --moni-ms

// Computes pseudo matching lengths against both positive and null indexes
python3 moni pseudo-ms -i ../data/example_positive_index/mock_comm_positive.fasta \
                       -p ../data/example_reads/all_reads.fa
python3 moni pseudo-ms -i ../data/example_null_index/mock_comm_null.fasta \
                       -p ../data/example_reads/all_reads.fa

// Computes matching statistics lengths against both positive and null indexes
python3 moni ms -i ../data/example_positive_index/mock_comm_positive.fasta \
                -p ../data/example_reads/all_reads.fa
python3 moni ms -i ../data/example_null_index/mock_comm_null.fasta \
                -p ../data/example_reads/all_reads.fa

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

# External resources

* [Big-BWT](https://github.com/alshai/Big-BWT.git)
    * [gSACA-K](https://github.com/felipelouza/gsa-is.git)
    * [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [Divsufsort](https://github.com/simongog/libdivsufsort.git)
* [klib](https://github.com/attractivechaos/klib)
* [r-index](https://github.com/maxrossi91/r-index.git)
* [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds.git)
* [bigrepair](https://gitlab.com/manzai/bigrepair.git)
* [shaped_slp](https://github.com/maxrossi91/ShapedSlp.git)
<!-- * [Google Benchmark](https://github.com/google/benchmark.git)
    * [Google Test](https://github.com/google/googletest) -->

# Citation 

Please, if you use this tool in an academic setting cite the following papers:

### MONI
    @inproceedings{RossiOLGB21,
    author      = { Massimiliano Rossi and 
                    Marco Oliva and
                    Ben Langmead and
                    Travis Gagie and
                    Christina Boucher},
    title       = {MONI: A Pangenomics Index for Finding MEMs},
    booktitle   = {Research in Computational Molecular Biology - 25th Annual 
                    International Conference, {RECOMB} 2021, Padova, Italy},
    volume      = {},
    series      = {Lecture Notes in Computer Science},
    pages       = {},
    publisher   = {Springer},
    year        = {2021}
    }

### SPUMONI

    @article{AhmedRGBL21,
    author      = { Omar Ahmed and
                    Massimiliano Rossi and
                    Travis Gagie and
                    Christina Boucher and
                    Ben Langmead},
    title       = {Pan-genomic Matching Statistics for Targeted Nanopore Sequencing},
    journal     = {CoRR},
    volume      = {abs/xxxx.xxxxx},
    year        = {2020},
    url         = {https://arxiv.org/abs/xxxx.xxxxx},
    archivePrefix = {arXiv},
    eprint      = {xxxx.xxxxx},
    }
# Authors

### Theoretical results:

* Omar Ahmed [![Generic badge](https://img.shields.io/badge/SPUMONI--<COLOR>.svg)](https://shields.io/)
* Christina Boucher
* Travis Gagie
* Ben Langmead
* Massimiliano Rossi

### Implementation:

* [Omar Ahmed](https://github.com/oma219) [![Generic badge](https://img.shields.io/badge/SPUMONI--<COLOR>.svg)](https://shields.io/)
* [Massimiliano Rossi](https://github.com/maxrossi91)

### Experiments

* [Omar Ahmed](https://github.com/oma219) [![Generic badge](https://img.shields.io/badge/SPUMONI--<COLOR>.svg)](https://shields.io/)
* [Marco Oliva](https://github.com/marco-oliva) [![Generic badge](https://img.shields.io/badge/MONI--<COLOR>.svg)](https://shields.io/)
* [Massimiliano Rossi](https://github.com/maxrossi91)

# Why "MONI" and "SPUMONI"?

**MONI** is the Finnish word for *multi*.

**SPUMONI** stands for Streaming PseUdo MONI.

# References

[1] Hideo Bannai, Travis Gagie, and Tomohiro I, *"Refining ther-index"*, Theoretical Computer Science, 812 (2020), pp. 96–108

[2] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini, *"Prefix-Free Parsing for Building Big BWTs"*, In Proc. of the 18th International Workshop on Algorithms in Bioinformatics (WABI 2018).

[3] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.
