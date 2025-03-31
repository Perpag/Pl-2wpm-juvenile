#!/bin/bash
#
#SBATCH --job-name=cellranger_count
#SBATCH --cpus-per-task=8
#SBATCH --output=log_count_cr.txt
echo "CellRanger running for P. lividus"
#$1 is sample name, $2 is sample number, $3 is number of cells, $4 is reference suffix
srun cellranger count --id=$1_$2_$3_$4 --transcriptome=../Data/Pl/Pliv_aH2p_fg_5kbext --fastqs=./FASTQs/$1_$2/ --sample=$1 --include-introns=true --no-bam --nosecondary --force-cells $3
cd ./$1_$2_$3_$4/outs/
mv web_summary.html $1_$2_$3_$4_web_summary.html
mv filtered_feature_bc_matrix $1_$2_$3_$4_filtered_feature_bc_matrix
cd ./$1_$2_$3_$4_filtered_feature_bc_matrix/
gzip -d *.gz
