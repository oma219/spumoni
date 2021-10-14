# SPUMONI :ice_cream:


SPUMONI is a software tool for performing rapid binary classifications on sequencing reads using a read's matching statistics (or a related quantity called pseudo-matching lengths). The potential use-cases can be scanning the prefix of nanopore reads to quickly make targeting decisions (application discussed in iScience paper below), performing host-depletion on a read set, as well as detecting for potential pathogens.

SPUMONI is based on another software tool called [MONI](https://github.com/maxrossi91/moni) which is a MEM-finder and aligner for pan-genomes. `MONI` uses prefix-free parsing of the text [2][3] to build the Burrows-Wheeler Transform (BWT) of the reference collection, the suffix array (SA) samples at the beginning and end of each run, and the threshold positions[1]. 

**Current Version:** 1.0

## Getting Started
There two ways to use SPUMONI. The easier way would be to go our release page and download the binaries. The other way would be to clone this GitHub repo as shown below and build from the source files. 

```sh
git clone https://github.com/oma219/spumoni
cd spumoni
mkdir build && cd build
cmake ..
make install
```
After running that last command above, in the `build/` directory there will be a `spumoni` executable to use.

## Next Steps
After installing SPUMONI on your machine, the first step would be to build an index over the reference you want to use for your experiment. This reference will be a FASTA file (and it can be a multi-FASTA for pan-genomes). For example, if you would want to do host depletion, you could index a group of host genomes. 

```sh
./spumoni build -r <ref_file> -M -P -f
```
The command above will build an index for the reference file you provide. This command uses `-M` and `-P` which means it will build an index for computing matching statistics (MSs), and another index for computing pseudo-matching statistics (PMLs). Our experiments show that using PMLs are more accurate at binary classification while also being ~3x faster, therefore that would be our recommendation.

An important point to make is that the `-f` indicates that input file is in FASTA format, which we expect is the dominate use-case. It will be the default option in later releases.

## Running Classification

Once you have an index for desired experiment, you can use the `spumoni run` command to generate either MSs or PMLs for each read against the reference file that you just indexed. Similarly to the `build` command, you provide the path the reference file along with the reads you want to classify.

```sh
./spumoni run -r <ref_file> -p <read_file> -P -f 
```

This command uses `-P` for computing PMLs (if you want MSs, use `-M` instead) and `-f` indicates the read file is a FASTA file. Again, we anticipate this will be the most likely scenario and adjustments will be made in the future to make it the default. 

Now, in order to perform classification of the reads, you will to have MSs or PMLs for each read against a "positive" database (just the regular FASTA) as well as a "null" database which is the just the reverse of the sequences. The command below shows how to use the analysis script in `analysis/` to order to get the classification result for each read. This script assumes the input are PMLs, however you can change that by using the `--ms` option.

```sh
python3 analyze_pml.py -p <pmls_against_pos_database.lengths>  \
                       -n <pmls_against_null_database.lengths>
```
Current development is aimed toward integrating the null database generation, and analysis of the reads into the main `spumoni` executable. However for now, for a more detailed explanation of how to use SPUMONI. You can look at the [README](https://github.com/oma219/spumoni/blob/main/analysis/README.md) in the `analysis/` folder.

## Getting Help

If you run into any issues or have any questions, please feel free to reach out to us either (1) through GitHub Issues or (2) reach out to me at omaryfekry [at] gmail.com

## Why "MONI" and "SPUMONI"?

**MONI** is the Finnish word for *multi*.

**SPUMONI** stands for Streaming PseUdo MONI.

## Citing SPUMONI

If you use the SPUMONI or MONI tools in your research project, please cite:
>Ahmed, O., Rossi, M., Kovaka, S., Schatz, M. C., Gagie, T., Boucher, C., & Langmead, B. (2021). Pan-genomic 
Matching Statistics for Targeted Nanopore Sequencing. iScience, 102696.

> Rossi, M., Oliva, M., Langmead, B., Gagie, T., & Boucher, C. (2021). MONI: A Pangenomics Index for Finding MEMs. *Proc*. RECOMB.

## References

[1] Hideo Bannai, Travis Gagie, and Tomohiro I, *"Refining ther-index"*, Theoretical Computer Science, 812 (2020), pp. 96–108

[2] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini, *"Prefix-Free Parsing for Building Big BWTs"*, In Proc. of the 18th International Workshop on Algorithms in Bioinformatics (WABI 2018).

[3] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.

## External Resources

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
