#!/bin/bash

#SBATCH --job-name=extract_counts
#SBATCH --output=extract_counts_%J.out
#SBATCH --partition=broadwl
#SBATCH --mem=8G

################################################################################
# The script takes bigwig files of DNase or ATAC-seq genome-wide count coverages
# then, it extracts tagcounts for a TF with 100bp flanks around motif
# and save in a tagcount matrix.
#
# requires bigwig tagcount files (in the earlier step)
# requires bwtool to extract counts from bigwig files
# requires rev_tagcount_bwtool.R to reverse counts on reverse strand
################################################################################

module load R

tf_name=$1
pwm_id=$2
bam_name=$3

ver_genome="hg19"

thresh_pValue="1e-4"
flank=100

pwm_name="${tf_name}_${pwm_id}_${thresh_pValue}"

dir_code=$HOME/projects/ATAC-seq/ATAC-seq_workflow/code_RCC
dir_tagcounts=/project/mstephens/ATAC_DNase/ATAC-seq_Olivia_Gray/results/ATAC-seq_tagcounts/
dir_count_matrix=/project/mstephens/ATAC_DNase/ATAC-seq_Olivia_Gray/results/ATAC-seq_count_matrix/${pwm_name}/
dir_sites=/project/mstephens/ATAC_DNase/motif_sites_JASPAR2018/${ver_genome}/candidate_sites/${thresh_pValue}/
filename_sites=${dir_sites}/${pwm_name}_flank${flank}_fimo_sites.bed

echo "PWM name: ${pwm_name}"
echo "Candidate sites: ${filename_sites}"
echo "Bam filename: ${bam_name}"

## Extract counts around motif matches into tagcount matrices
## Check the existence of tagcount bigwig files and candidate sites
bam_basename=$(basename ${bam_name} .bam)
bw_file="${dir_tagcounts}/${bam_basename}.5p.0base.tagcount"
file_tagcount_fwd=${bw_file}.fwd.bw
file_tagcount_rev=${bw_file}.rev.bw
echo "Tagcount fwd name: ${file_tagcount_fwd}"
echo "Tagcount rev name: ${file_tagcount_rev}"

if [ ! -s ${file_tagcount_fwd} ] || [ ! -s ${file_tagcount_rev} ]; then
  echo "ERROR: Tagcount bigwig files for: ${bam_name} does not exist!"
elif [ ! -s ${filename_sites} ]; then
  echo "ERROR: Candidate sites ${filename_sites} does not exist!"
else

  echo "Counting cuts around the candidate sites of ${pwm_name} with ${flank}bp flanking windows ..."

  mkdir -p ${dir_count_matrix}

  file_count_matrix_fwd="${dir_count_matrix}/${pwm_name}_${bam_basename}_fwdcounts.m"
  file_count_matrix_rev="${dir_count_matrix}/${pwm_name}_${bam_basename}_revcounts.m"

  ## count matrix around the candidate sites
  echo "match matrix of ${file_count_matrix_fwd} ..."
  cut -f 1-4 ${filename_sites} | bwtool extract bed stdin ${file_tagcount_fwd} ${file_count_matrix_fwd} -fill=0 -decimals=0 -tabs

  echo "match matrix of ${file_count_matrix_rev} ..."
  cut -f 1-4 ${filename_sites} | bwtool extract bed stdin ${file_tagcount_rev} ${file_count_matrix_rev} -fill=0 -decimals=0 -tabs

  Rscript ${dir_code}/rev_tagcount_bwtool.R ${file_count_matrix_fwd} ${file_count_matrix_rev} ${filename_sites}

  # compress the files to .gz
  gzip -f ${file_count_matrix_fwd} ${file_count_matrix_rev}
  echo "Finished counting for ${pwm_name} in ${bam_name}."

fi
