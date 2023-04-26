# Change Log for SPUMONI

- This log started being maintained at v2.0.0, therefore, there are not specific version labels for previous versions of SPUMONI besides the git commit id.

## v2.0.6 - latest
- Fixed bug in thirdparty compilation for pfp-thresholds

## v2.0.5
- Fixed bug that occurs when providing a filelist with no document numbers
- Added some coloring to the logging for SPUMONI

## v2.0.4
- Added the `--hash-mod` option back to be able to adjust PFP settings
- TODO: need to update error messages to inform user when PFP crashes

## v2.0.3
- Added `*index_stats.txt` files that record the memory consumption of index components.

## v2.0.2
- Fixed bug where it looks to check path of document array prior to computation
- Fixed warning messages in spumoni run, when it says it cannot find a specific file, it was not printing the exact 
  path it was checking so it made it confusing.
- Updated the usage statment for the -r option in run.

## v2.0.1
- Updated warning message for output index prefix, force users to use './' for same directory files
- Added check after minimizer digestion to catch scenarios where general text are provided as FASTA files, because the resulting file might be empty and lead to errors.
- Added ability to not add reverse complement to FASTA files if requested.

## v2.0.0 
- Allows for the digestion of input sequences into either a sequence of DNA-based minimizers (e.g ACAT) or minimizers that are alphabet characters themselves.
- Minimizer digestion is used during build and run sub-commands if requested by the user.
- Classification test is integrated directly into the SPUMONI code by generating a database of null MS to compare to the online query MS.
- Classification test is changed to simple threshold test over the maximum positive MS and the largest "common" null MS.
- Default input reference sequences are FASTA, however, there is a general text option that allows for MS/PML with respect to a "binary" file.
- Small command-line parsing changes like allowing for output prefix to be specified for index.
- Initially released via this paper (https://www.biorxiv.org/content/10.1101/2022.09.08.506805v1.abstract)

## v1.0.0
- Initial release of SPUMONI, consisted strictly of the the ability to compute MS and PML for queries. Classification test was not built-in and general text was the default mode.
- Released vis this paper (https://www.sciencedirect.com/science/article/pii/S2589004221006647)