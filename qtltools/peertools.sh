#! /bin/bash

## N = 60

# Delta-QTLs
peertool -f input_data/Expression_delta_tpm01_20perc_reads6_20perc_for_peer.csv \
         -n 60 \
         -c input_data/covariates_patients_age_sex_smoker_PCs_for_peer.csv \
         --bound_tol 0.1 \
         --var_tol 0.001 \
         -o output_delta_boundtol01_vartol1e3

# Tumor-QTLs
peertool -f input_data/Expression_tumor_tpm01_20perc_reads6_20perc_for_peer.csv \
         -n 60 \
         -c input_data/covariates_patients_age_sex_smoker_PCs_for_peer.csv \
         --bound_tol 0.1 \
         --var_tol 0.001 \
         -o output_tumor_boundtol01_vartol1e3_v2

# Lung-QTLs
peertool -f input_data/Expression_lung_tpm01_20perc_reads6_20perc_for_peer.csv \
         -n 60 \
         -c input_data/covariates_patients_age_sex_smoker_PCs_for_peer.csv \
         --bound_tol 0.1 \
         --var_tol 0.001 \
         -o output_lung_boundtol01_vartol1e3
