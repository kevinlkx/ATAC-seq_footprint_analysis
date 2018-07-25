#!/bin/bash

#SBATCH --job-name=fitCentipede
#SBATCH --mem=10G
#SBATCH --output=fitCentipede_%J.out
#SBATCH --partition=broadwl

################################################################################
# fit CENTIPEDE model to ATAC-seq cuts and PWM scores
#
# requires CENTIPEDE R package
################################################################################

module load R

tf_name=$1
pwm_id=$2
thresh_mapability=$3
thresh_PWMscore=$4
thresh_pValue=$5

## loop through bam files
bamfiles=("H1_nomito_rdup.bam" "H2_nomito_rdup.bam" "H3_nomito_rdup.bam" "N1_nomito_rdup.bam" "N2_nomito_rdup.bam" "N3_nomito_rdup.bam")

for bam_name in "${bamfiles[@]}"
do
   echo "${bam_name}"
   Rscript ~/projects/ATAC-seq/ATAC-seq_workflow/code_RCC/fitCentipede_ATAC-seq_counts.R ${tf_name} ${pwm_id} ${bam_name} ${thresh_mapability} ${thresh_PWMscore} ${thresh_pValue}
done
