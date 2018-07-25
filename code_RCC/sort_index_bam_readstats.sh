#!/bin/bash
#
#SBATCH --job-name=bam
#SBATCH --mem=6G

############################################################
# sort, index and count mapped reads for a bam file
############################################################

bam_file=$1
dir_data=$2

dir_readstats="${dir_data}/readstats"
mkdir -p $dir_readstats

## sort and index bam files
echo "sort and index bamfile: ${bam_file}"

cd ${dir_data}/bam

bam_file_base=$(basename ${bam_file} .bam)

echo "base bamfile: ${bam_file_base}"

samtools sort ${bam_file_base}.bam ${bam_file_base}.sorted

samtools index ${bam_file_base}.sorted.bam

samtools idxstats ${bam_file_base}.sorted.bam > ${dir_readstats}/${bam_file_base}.sorted.bam_idxstats.txt

echo "Number of mapped reads: "
samtools view -c ${bam_file_base}.sorted.bam

