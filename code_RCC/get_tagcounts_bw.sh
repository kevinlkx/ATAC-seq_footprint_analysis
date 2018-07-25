#!/bin/bash
#
#SBATCH --job-name=ATAC_bw
#SBATCH --mem=6G

### The script takes a ATAC-seq bam file, count tagcounts in bigwig (bw) format
### The script requires samtools, bwtool packages installed
### requires genome_coverage_bamToBigwig.sh to count bam reads as bigwig tagcounts

bam_file=$1
dir_dnase=$2


## dir_code directory for data analysis scripts
dir_code=/data/reddylab/projects/GGR/analyses/occupancy_prediction/code

## dir_bam directory for ATAC-seq bam files
dir_bam=$dir_dnase/bam

## dir_dnase directory for processed ATAC-seq data
dir_dnase_tagcount=$dir_dnase/tagcount
mkdir -p $dir_dnase_tagcount

echo "counting ATAC match for $bam_file at $dir_dnase_tagcount"

## output tagcount bigwig files
file_tagcount_fwd=$dir_dnase_tagcount/"${bam_file}_5p.0base.tagcount.fwd.bw"
file_tagcount_rev=$dir_dnase_tagcount/"${bam_file}_5p.0base.tagcount.rev.bw"

## if either the fwd or rev tagcount does not exist or size equal to 0
## then count the tagcount and save in bigwig format
if [ ! -s $file_tagcount_fwd ] || [ ! -s $file_tagcount_rev ]; then
  
  echo "creating ATAC tagcount bigwig files for:" $bam_file "at" $dir_dnase_tagcount
  minimumsize=100000
  ratio_scale=1
  
  ## genomeCoverageBed -dz -5 -i  -g > 5' end tag counts (0 based using -dz in genomeCoverageBed)
  ##  actualsize=$(wc -c "$file" | cut -f 1 -d ' ')

  ## if the file not exists or if the file is too small
  if [[ ! -s "$file_tagcount_fwd" || $(wc -c "$file_tagcount_fwd" | cut -f 1 -d ' ') -lt $minimumsize ]]; then
      ## forward strand tag counts
      echo "counting forward strand tags ..."
      sh $dir_code/genome_coverage_bamToBigwig.sh $dir_bam/$bam_file "+" $file_tagcount_fwd $ratio_scale
  else
      echo "skip counting forward strand tags. "
  fi

  if [[ ! -s "$file_tagcount_rev" || $(wc -c "$file_tagcount_rev" | cut -f 1 -d ' ') -lt $minimumsize ]]; then
      ## reverse strand tag counts
      echo "counting reverse strand tags ..."
      sh $dir_code/genome_coverage_bamToBigwig.sh $dir_bam/$bam_file "-" $file_tagcount_rev $ratio_scale
  else
      echo "skip counting reverse strand tags."
  fi

else

  echo "bigwig tagcount files already exist."

fi

