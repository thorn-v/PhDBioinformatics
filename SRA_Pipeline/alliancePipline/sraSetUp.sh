#! /usr/bin/env bash

#SBATCH --job-name=SRAsetup
#SBATCH -A def-jpxu
#SBATCH --cpus-per-task=32  #can change
#SBATCH --time=03:00:00     #can change
#SBATCH --mem-per-cpu=2G    #can change

module load sra-toolkit

downloadSRA() {
    # Download Script
    prefetch $1 -O SRAs/$1
    fasterq-dump SRA/$1 -O Raw_Fastqs                   
    gzip Raw_Fastqs/${1}.*        #gzip extracted fastq(s)
    echo "${1} complete"
}
export -f downloadSRA

parallel -j $SLURM_CPUS_PER_TASK --joblog sras.log downloadSRA {} :::: ./accs.txt


