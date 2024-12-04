#!/bin/bash

# bash script that writes a tsv file with sample name and path to the count matrix

# header
echo -e "sample_name\tpath" > resources/sample_sheet.tsv

# loop over all samples
for sample in /data/cephfs-1/work/projects/crc-patients-treatment-and-relapse/work/scRNA/raw_data/*/count/sample_filtered_feature_bc_matrix.h5
do
    # sample name is the word that contains CE_SC
    sample_name=$(echo $sample | grep -o "CE_SC_[A-Za-z0-9_]*")
    echo -e "$sample_name\t$sample" >> resources/sample_sheet.tsv
done
