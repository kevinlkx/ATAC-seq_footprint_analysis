#!/bin/bash
#SBATCH --job-name=fimo
#SBATCH --mem=8G
#SBATCH --output=fimo_%J.out
#SBATCH --partition=broadwl

##############################################################################
# get motif matches using FIMO
#
# requires FIMO from MEME suite
# FIMO instructions: http://meme-suite.org/doc/fimo.html?man_type=web
#installation: http://meme-suite.org/doc/install.html?man_type=web
##############################################################################


## parameters
tf_name=$1
pwm_id=$2
thresh_pValue=$3 # 1e-4

###### example ######
# tf_name="CTCF"
# pwm_id="MA0139.1"
# thresh_pValue=1e-4
#####################

flank=100
ver_genome=hg19
max_fimo_sites=1000000

## extract pwm meme file
echo "pwm_id: $pwm_id, tf_name: $tf_name"

pwm_name="${tf_name}_${pwm_id}_${thresh_pValue}"

# directory for Jaspar motifs
dir_pwm="/project/mstephens/ATAC_DNase/motif_database/JASPAR_2018/JASPAR2018_CORE_non-redundant_pfms_meme"
pwm_filename="${dir_pwm}/${pwm_id}.meme"

# reference genome fasta file
genome_filename="/project/mstephens/ATAC_DNase/ref_genome/${ver_genome}.fa"

# output directory for FIMO motif matching results
dir_output="/project/mstephens/ATAC_DNase/motif_sites_JASPAR2018/"
mkdir -p ${dir_output}

## FIMO matching motifs
dir_motif_matches="${dir_output}/FIMO/${pwm_name}"
filename_motif_matches="${dir_motif_matches}/fimo.txt"

mkdir -p ${dir_motif_matches}
echo "Begin FIMO matching motif ..."

fimo --oc ${dir_motif_matches} --verbosity 2 --text \
    --bgfile --uniform-- --thresh ${thresh_pValue} --max-stored-scores ${max_fimo_sites} \
    ${pwm_filename} ${genome_filename}  \
    > ${filename_motif_matches}

echo "Finish FIMO: ${dir_motif_matches}"
