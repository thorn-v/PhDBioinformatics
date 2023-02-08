#! /usr/bin/env bash

#SBATCH --job-name=gvcf_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=180
#SBATCH --array=1-2317%100

#config file for accessions
config=~/scratch/config.txt

# get the accession from the list
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

/home/vthorn/scratch/PhDBioinformatics_Code/snpCallingPipeline.sh -i ~/scratch/Afumigatus_WGSA_SRAs/${sample} -r ~/scratch/Afumigatus_Reference/A_fumigatus_Af293/GCA_000002655.1/GCA_000002655.1_ASM265v1_genomic.fna -n 8

rm -r ~/scratch/Temp_WD/"${sample}"*

