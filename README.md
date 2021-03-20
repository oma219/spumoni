# SPUMONI v0.1.0
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


# Usage

### Construction of the index:
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
  --spumoni             build spumoni index (default: False)

```
The default index is the `spumoni` index. If you want to build the `moni-ms` index you can use the `--moni-ms` flag. If you want to build both you can use `--moni-ms --spumoni`. Building `moni-ms` and `spumoni` together is faster than building them separately.

### Computing the matching statistics with moni-ms:
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

### Computing the matching statistics with spumoni:
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

# Example
### Download

```console
git clone https://github.com/maxrossi91/spumoni
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

### Run

```console
./moni build -r ../data/yeast.fasta -f
./moni pseudo-ms -i ../data/yeast.fasta -p ../data/query.fasta 
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

### Moni
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

### Spumoni

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

# Why "MONI"?

**Moni** is the Finnish word for *multi*.

# References

[1] Hideo Bannai, Travis Gagie, and Tomohiro I, *"Refining ther-index"*, Theoretical Computer Science, 812 (2020), pp. 96–108

[2] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini, *"Prefix-Free Parsing for Building Big BWTs"*, In Proc. of the 18th International Workshop on Algorithms in Bioinformatics (WABI 2018).

[3] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.
