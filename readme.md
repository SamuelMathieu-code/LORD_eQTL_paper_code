# Code from "*Germline variants affect lung to tumor transcription dynamics, NSCLC relapse and mortality*"

This repository contains code for eQTL analysis and it's downstream analysis and for the survival GWAS analysis. R and Julia were used as main analysis softwares.

- eQTL analysis and fine-mapping was performed (see `qtltools` folder)
- Enrichment of eVariants in functional annotations was performed (see `ChromHMM_enrich` folder)
- Enrichment of eGenes in cell type specific gene sets was performed (see `HPA_gene_sets_enrich` folder)
- Transcription factor binding sites were assessed (see `TFBS_disrupt` folder)
- Differetial expression analysis was performed (see `DEG` folder)
- Module-QTL analysis was performed (see `module_qtl` folder)
- Survival GWAS analysis was performed (see `survival_gwas` folder)]
- Deconvolution analysis
- Comparison between common somatic alterations and germline effects on tumor gene expression

The `utils` folders contains the code for `convert_to_hgnc` and `plink_bed_encoding_to_genotype` functions.

Warning : path to data files are only an indication for script readability, eventual users of this code should replace the paths in the scripts by theirs.
