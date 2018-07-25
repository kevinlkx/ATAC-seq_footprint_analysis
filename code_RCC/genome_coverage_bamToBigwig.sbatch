#!/bin/bash

#SBATCH --job-name=coverage
#SBATCH --output=coverage_%J.out
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=10G

#################################################################
# get 5' end coverage for bam files                             #
#################################################################

bamfile=$1

sh ~/projects/ATAC-seq/ATAC-seq_workflow/code_RCC/genome_coverage_bamToBigwig.sh ${bamfile}
