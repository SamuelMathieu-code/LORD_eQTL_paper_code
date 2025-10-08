# Perfectos TF disrupt pipeline

## 0. Run Perfectos-Ape

First, all credible set SNPs should be run in [perfectos-ape](https://opera.autosome.org/perfectosape) using the JASPAR 2024 database.

## 1. Filter expressed TFs
To filter out TFs that are not expressed : use `get_curated_motifs.jl`. It takes in the expression matrices located in the parent directory and `merge_output_perfectosAPE.tsv`. It outputs a list of motifs for which each TF within one motif are expressed at least at 1 TPM in either Lung or Tumor : `expressed_motifs.csv`.

## 2. Cluster Motifs

To make a similarity matrix of each expressed TF, we used the `PWMEnrich` R package with the [motifSimilarity](https://rdrr.io/bioc/PWMEnrich/man/motifSimilarity.html) function. To install the package, make a conda environement with this command line.

```
conda create -n pwmenrich_r -c conda-forge r-curl r-biocmanager
```

Then take the script `make_sim_matrix_jaspar.jl` and modify this line: `ENV["R_HOME"] = "{Your Conda}/miniconda3/envs/pwmenrich_r/lib/R"`. Running the script (be carefull to use only one thread! Verify `Threads.nthreads == 1`.), you will obtain a similarity matrix in the file `expressed_tfs_motif_similrity_matrix.csv`. Modify the `addprocs(20, lazy = false)` statement to add the number of processors you want.

## 3. Make the enrichment based on credible set numbers while accounting for mean # of variants in credible sets

Enrichment analysis is implemented in `enrich.jl`.