# SPUMONI :ice_cream: ![GitHub release (latest by date)](https://img.shields.io/github/v/release/oma219/spumoni) ![GitHub](https://img.shields.io/github/license/oma219/spumoni?color=green)


SPUMONI is a software tool for **performing rapid read lassifications on sequencing reads using a read's matching statistics (or a related quantity called pseudo-matching lengths).** 

SPUMONI is based on another software tool called [MONI](https://github.com/maxrossi91/moni) which is a MEM-finder and aligner for pan-genomes. `MONI` uses prefix-free parsing of the text [2][3] to build the Burrows-Wheeler Transform (BWT) of the reference collection, the suffix array (SA) samples at the beginning and end of each run, and the threshold positions[1]. 

**For more details on installation/using SPUMONI, please refer to the [wiki page](https://github.com/oma219/spumoni/wiki/1.-Home) which covers those areas in detail.**

## Getting Started
There are two ways to use SPUMONI. The easier way would be to use the Docker/Singularity images which have SPUMONI pre-installed and ready to use. The other way would be to clone this GitHub repo as shown below and build from the source files. 

**(1) Use Docker or Singularity Containers**

```
docker pull oma219/spumoni:latest
docker run oma219/spumoni spumoni -h
```

Since `docker` requires root permissions, `singularity` can be used instead on shared clusters as shown below.

```
singularity pull docker://oma219/spumoni:latest
singularity run spumoni_latest.sif spumoni -h
```

Both sets of commands above will download the containers and run `spumoni` to print out the sub-commands that are available to run.

**(2) Build from Source**

To build from source, it requires the following dependences: `zlib`, `cmake`, `gcc`. See the [wiki](https://github.com/oma219/spumoni/wiki/3.-Installation) for detailed information.
```sh
git clone https://github.com/oma219/spumoni
cd spumoni
mkdir build && cd build
cmake ..
make install

export SPUMONI_BUILD_DIR=$(pwd)
# use spumoni ...
```
After running that last command above, in the `build/` directory there will be a `spumoni` executable to use.

## Step 1: Building an Index

After installing SPUMONI on your machine, the first step would be to build an index over the reference you want to use for your experiment. This reference will be a FASTA file (and it can be a multi-FASTA for pan-genomes). SPUMONI allows you to either build it over a single FASTA file, or you can specify a list of genomes that you want to include in the index. [See the wiki for more details.](https://github.com/oma219/spumoni/wiki/4.-Building-SPUMONI-Indexes) 

For example, if you would want to do host depletion, you could index a human genome and build both the matching statistic and pseudo-matching lengths index using the command below: 

```sh
./spumoni build -r human_genome.fa -M -P -m 
```
The command above will build an index for the reference file you provide. This command uses `-M` and `-P` which means it will build an index for computing matching statistics (MSs), and another index for computing pseudo-matching statistics (PMLs). Our experiments show that using PMLs are more accurate at binary classification while also being ~3x faster, therefore that would be our recommendation.

## Step 2: Running Classification

Once you have an index for desired experiment, you can use the `spumoni run` command to generate either MSs or PMLs for each read against the reference file that you just indexed. Similarly to the `build` command, you provide the path the reference file along with the reads you want to classify.

```sh
./spumoni run -r spumoni_full_ref.bin -p reads.fa -P -c
```

The `spumoni_full_ref.bin` file is produced during the build process since SPUMONI digests references into a sequence of minimizers, and builds an index over the digestion. [Refer to the wiki to understand what file you should use as your index.](https://github.com/oma219/spumoni/wiki/5.-Running-SPUMONI-on-Input-Reads)

This command uses `-P` for computing PMLs (if you want MSs, use `-M` instead) which will be used to classify the reads. Additionally, the command uses the `-c` option to write out the classifications to a report file.

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
