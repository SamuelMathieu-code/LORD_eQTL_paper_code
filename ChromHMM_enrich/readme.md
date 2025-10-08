# Functional Genomic enrichments of finemapped eQTLs

1. `eqtls_to_bed.jl` and `make_annots_bed.sh` map all positions tested by QTLTools to functional category in chromhmm universal resource.
2. `enrichments.jl` maps credible set eVariants to annotations provided by `make_annots_bed.sh` and performs hypergeometric enrichment with corresponding background (Lung, Tumor, Delta) and makes a figure.
