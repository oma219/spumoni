# Change Log for SPUMONI

- This log started being maintained at v2.0.0, therefore, there are not specific version labels for previous versions of SPUMONI besides the git commit id.

## v2.0.0 - latest
---
- Allows for the digestion of input sequences into either a sequence of DNA-based minimizers (e.g ACAT) or minimizers that are alphabet characters themselves.
- Minimizer digestion is used during build and run sub-commands if requested by the user.
- Classification test is integrated directly into the SPUMONI code by generating a database of null MS to compare to the online query MS.
- Classification test is changed to simple threshold test over the maximum positive MS and the largest "common" null MS.
- Default input reference sequences are FASTA, however, there is a general text option that allows for MS/PML with respect to a "binary" file.
- Small command-line parsing changes like allowing for output prefix to be specified for index.
- Initially released via this paper (https://www.biorxiv.org/content/10.1101/2022.09.08.506805v1.abstract)


## v1.0.0
---
- Initial release of SPUMONI, consisted strictly of the the ability to compute MS and PML for queries. Classification test was not built-in and general text was the default mode.
- Released vis this paper (https://www.sciencedirect.com/science/article/pii/S2589004221006647)