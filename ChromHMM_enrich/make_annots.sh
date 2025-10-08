#! /bin/bash

wget https://public.hoffman2.idre.ucla.edu/ernst/UUKP7/hg38lift_genome_100_segments.bed.gz
mkdir universal_chromhmm
mv hg38lift_genome_100_segments.bed.gz universal_chromhmm/
gunzip universal_chromhmm/hg38lift_genome_100_segments.bed.gz

bedtools sort -i all_tested_positions_tumor.bed > all_tested_positions_tumor_sorted.bed
bedtools sort -i all_tested_positions_lung.bed > all_tested_positions_lung_sorted.bed
bedtools sort -i all_tested_positions_delta.bed > all_tested_positions_delta_sorted.bed

bedtools intersect -wb -b universal_chromhmm/hg38lift_genome_100_segments.bed -a all_tested_positions_tumor_sorted.bed > all_tested_positions_tumor_sorted_chromhmm_uni.bed
bedtools intersect -wb -b universal_chromhmm/hg38lift_genome_100_segments.bed -a all_tested_positions_lung_sorted.bed > all_tested_positions_lung_sorted_chromhmm_uni.bed
bedtools intersect -wb -b universal_chromhmm/hg38lift_genome_100_segments.bed -a all_tested_positions_delta_sorted.bed > all_tested_positions_delta_sorted_chromhmm_uni.bed

mv all_tested_positions_*_chromhmm_uni.bed universal_chromhmm

sort -k1,1 -k2,2n -k3,3n -u universal_chromhmm/all_tested_positions_delta_sorted_chromhmm_uni.bed > universal_chromhmm/all_tested_positions_delta_sorted_chromhmm_uni_unique.bed
sort -k1,1 -k2,2n -k3,3n -u universal_chromhmm/all_tested_positions_lung_sorted_chromhmm_uni.bed > universal_chromhmm/all_tested_positions_lung_sorted_chromhmm_uni_unique.bed
sort -k1,1 -k2,2n -k3,3n -u universal_chromhmm/all_tested_positions_tumor_sorted_chromhmm_uni.bed > universal_chromhmm/all_tested_positions_tumor_sorted_chromhmm_uni_unique.bed