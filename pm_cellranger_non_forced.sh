#!/bin/bash

echo "CellRanger running for Patiria miniata ovary data"

cellranger count --id=$1_$2_$3 --transcriptome=/home/seastar/Data/Pmin3_reference/Pmin3_protein_coding_cellranger_ref  --fastqs=/home/seastar/Data/SingleCell/00_fastq/$1 --sample=$1 --localcores=32 --localmem=256 --nosecondary --no-bam
cd ./$1_$2_$3/outs/filtered_feature_bc_matrix
gzip -d *.gz
mv features.tsv genes.tsv

