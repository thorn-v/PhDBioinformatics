#! /usr/bin/env bash

#SBATCH --job-name=gvcf_pipeline_pt1
#SBATCH --ntasks=1
#SBATCH --mem=8G #if you change this, make sure you update the -m flag with the same number
#SBATCH --time=300
#SBATCH --array=1-??%100 # !! replace the ?? with the total number of samples you have (number of lines in your configfile). if you have less than 100 samples, remove the "%100"

#config file for accessions
config=~/scratch/PATH/TO/configfile.txt
# if your file is names something other than configfile.txt, then use your filename.

# get the accession from the list
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $config)
# if using your own samples, should be the names of the samples in relation to the fastq files (do not include the ".fastq")
# Sample names should NOT have spaces in them 

# !! Reference must be indexed with both bwa and samtools
ref=~/PATH/TO/Reference.fa

##== Begin Pipeline ==##
# assumes you are in the folder just above bin, if not, change below paths (./bin)
# the -A flag is needed because we are running on Alliance HCCs

## Download the File ##
./bin/downloadSRAs.sh -A -s $sample
# if not using SRA or are already downloaded, remove this command and make sure your fastq files are in folder called "Raw_Fastqs" in your current working directory.


## QC and Trim the Raw Fastq ##
./bin/trimAndQC.sh -A -s $sample 
#can add -q for quality cutoff (defaults to 30)
#        -l for length cutoff (defaults to 30)
#        -c for minimum number of reads needed to keep sample for further processing (defaults to 0 which does not remove any samples)

## Align to reference and de-dupe ##
./bin/alignAndDeDupe.sh -A -r $ref -s $sample
# again, can add length and quality cut-offs
# and can add minimum depth (-c) in order to keep sample (defaults to keep all samples)

## Initial GVCF ##
./bin/initalGVCFsCalling.sh -A -r $ref -s $sample -m 8
# if don't want to do the whole genome, use -L to add your own interval list in GATK format. Otherwise, script will create one.

echo "finished for sample ${sample}"
