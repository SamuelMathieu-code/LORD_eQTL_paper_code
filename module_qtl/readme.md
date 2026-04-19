# Module-QTL analysis

## 1. get the co-expressed modules

use `wgcna.R`.

## 2. make the module-QTL analysis

The CLI tool `module_qtl_cli_script.jl` makes module-QTL analysis considering every possible Module-variant pair.

To get help on usage, do: `julia module_qtl_cli_script.jl --help`

In the paper, it was used it as follows using all variants that where the most significant hit for an any eGene in the cis-eqtl analysis:

```bash
OMP_NUM_THREADS=20 julia -t 10 --gcthreads=5 --project=. --heap-size-hint=20GB --startup-file=no modqtls_v6.jl -c input_modqtls_v6/covariates_lung_60peer.csv \
    -s input_modqtls_v6/all_kept_snps_lung.txt \
    -m out_wgcna_v6/modules_lung.tsv \
    -e ../lung_tpms.tsv \
    --genotypes ../genotypes/plink_files/TOPMed_chr_all_LORD_pairedRNA.reordered \
    --noorderquantnorm \
    --outputfile out_modqtls_v6/out_lung_60peer_noorderquantnorm.csv

OMP_NUM_THREADS=20 julia -t 10 --gcthreads=5 --project=. --heap-size-hint=20GB --startup-file=no modqtls_v6.jl -c input_modqtls_v6/covariates_tumor_60peer.csv \
    -s input_modqtls_v6/all_kept_snps_tumor.txt \
    -m out_wgcna_v6/modules_tumor.tsv \
    -e ../tumor_tpms.tsv \
    --genotypes ../genotypes/plink_files/TOPMed_chr_all_LORD_pairedRNA.reordered \
    --noorderquantnorm \
    --outputfile out_modqtls_v6/out_tumor_60peer_noorderquantnorm.csv
```

We then only considered module-variant pairs for which the variant was the best hit for a gene of the module (see `reprioritize.jl`).

## 3. Make the sensitivity analyses for all significant pairs

See `sensitivity.jl`.

## 4. Mediation

See `mediation.jl` for mediation analysis for two module-QTLs.
