#! /usr/bin/env bash

#SBATCH --job-name=sra_to_vcf_array_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=180
#SBATCH --array=1:2317%100

#config file for accessions
config=~/scratch/config.txt

# get the accession from the list
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
# echo "This is array task ${SLURM_ARRAY_TASK_ID} and the sample name is ${sample}." >> arrayOutput.txt

/home/vthorn/scratch/PhDBioinformatics_Code/snpCallingPipeline.sh -i ~/scratch/Afumigatus_WGSA_SRAs/${sample} -r ~/scratch/Afumigatus_Reference/A_fumigatus_Af293/GCA_000002655.1/GCA_000002655.1_ASM265v1_genomic.fna -n 4

rm -r ~/scratch/Temp_WD/"${sample}"*

