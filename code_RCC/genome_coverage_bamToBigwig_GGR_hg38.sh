#!/bin/bash

#SBATCH --job-name=coverage
#SBATCH --output=coverage_%J.out
#SBATCH --mem=10G
#SBATCH --partition=broadwl

################################################################################
# get 5' end coverage for bam files
# Compute the DNase or ATAC-seq cleavage (5' end) coverge along the genome
#
# requires samtools, bedtools, and bedGraphToBigWig (from UCSC)
################################################################################

bam_file=$1

dir_results="/project/mstephens/ATAC_DNase/ATAC-seq_GGR_Duke/results/ATAC-seq_tagcounts/"
mkdir -p ${dir_results}

genome_chrom_sizes="/project/mstephens/ATAC_DNase/ref_genome/hg38.chrom.sizes"

echo "Compute genome coverage for ${bam_file}"
echo "Output directory: ${dir_results}"

## sort bam file
bam_dir_name=$(dirname ${bam_file})
bam_basename=$(basename ${bam_file} .bam)

cd ${bam_dir_name}

echo "base bamfile: ${bam_basename}"
echo "sort and index bam file"

samtools sort -o ${bam_basename}.sorted.bam ${bam_basename}.bam

samtools index ${bam_basename}.sorted.bam

samtools idxstats ${bam_basename}.sorted.bam > ${bam_basename}.idxstats.txt

bam_sorted_name="${bam_dir_name}/${bam_basename}.sorted.bam"

## count genome coverage, 5' end, save in bedGraph format, then convert to bigWig format

bw_file="${dir_results}/${bam_basename}.5p.0base.tagcount"

## forward strand
echo "Count genome coverage on + strand"
samtools view -b ${bam_sorted_name} | \
genomeCoverageBed -bg -5 -strand "+" -ibam stdin -g ${genome_chrom_sizes} > ${bw_file}.bedGraph

echo "coverting bedGraph to BigWig..."
sort -k1,1 -k2,2n ${bw_file}.bedGraph > ${bw_file}.sorted.bedGraph
bedGraphToBigWig ${bw_file}.sorted.bedGraph ${genome_chrom_sizes} ${bw_file}.fwd.bw
rm ${bw_file}.bedGraph ${bw_file}.sorted.bedGraph

## reverse strand
echo "Count genome coverage on - strand"
samtools view -b ${bam_sorted_name} | \
genomeCoverageBed -bg -5 -strand "-" -ibam stdin -g ${genome_chrom_sizes} > ${bw_file}.bedGraph

echo "coverting bedGraph to BigWig..."
sort -k1,1 -k2,2n ${bw_file}.bedGraph > ${bw_file}.sorted.bedGraph
bedGraphToBigWig ${bw_file}.sorted.bedGraph ${genome_chrom_sizes} ${bw_file}.rev.bw
rm ${bw_file}.bedGraph ${bw_file}.sorted.bedGraph

echo "Finished computing genome coverage => ${bw_file}.bw"
