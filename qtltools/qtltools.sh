#! /bin/bash

# Definition: Initialize a semaphore with a given number of tokens.
opensem(){
         mkfifo pipe-$$
         exec 3<>pipe-$$
         rm pipe-$$
         local i=$1
         for((;i>0;i--)); do
                  printf %s 000 >&3
         done
}

# Definition: Run the given command asynchronously and pop/push tokens.
runwithlock(){
         local x
         # a) This read waits until there is something to read.
         read -u 3 -n 3 x && ((0==x)) || exit $x
         (
          ( "$@"; )
         # b) Push the return code of the command to the semaphore.
         printf '%.3d' $? >&3
         )&
}


echo claculate Delta-QTLs with 60 Peer factors
for l in $(seq 1 40); do

    runwithlock QTLtools cis --vcf input_data/TOPMed_chr_all_LORD_pairedRNA.vcf.gz \
                 --bed input_data/Expression_delta_tpm01_20perc_reads6_20perc_for_qtltools_sorted.bed.gz \
                 --cov input_data/covariates_delta_age_sex_smoker_10PC_60peer_for_qtltools.txt \
                 --window 1000000 --nominal 1.0 --silent --std-err \
                 --chunk $l 40 --out out_delta_60/out_${l}.txt
done
wait

echo calculate Lung-QTLs
for j in $(seq 1 40); do

    runwithlock QTLtools cis --vcf input_data//TOPMed_chr_all_LORD_pairedRNA.vcf.gz \
                 --bed input_data/Expression_lung_tpm01_20perc_reads6_20perc_for_qtltools_sorted.bed.gz \
                 --cov input_data/covariates_lung_age_sex_smoker_10PC_60peer_for_qtltools.txt \
                 --window 1000000 --nominal 1.0 --silent --std-err \
                 --chunk $j 40 --out out_lung_60/out_${j}.txt
done
wait

echo calculate Tumor-QTLs
for k in $(seq 1 40); do

    runwithlock QTLtools cis --vcf input_data/TOPMed_chr_all_LORD_pairedRNA.vcf.gz \
                 --bed input_data/Expression_tumor_tpm01_20perc_reads6_20perc_for_qtltools_sorted.bed.gz \
                 --cov input_data/covariates_tumor_age_sex_smoker_10PC_60peer_for_qtltools.txt \
                 --window 1000000 --nominal 1.0 --silent --std-err \
                 --chunk $k 40 --out out_tumor_60/out_${k}.txt
done
wait

echo done