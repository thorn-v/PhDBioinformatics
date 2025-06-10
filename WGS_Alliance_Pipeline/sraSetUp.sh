#! /usr/bin/env bash

#SBATCH --job-name=SRAsetup
#SBATCH -A def-jpxu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48  #can change
#SBATCH --time=12:00:00     #can change
#SBATCH --mem-per-cpu=4000    #can change
## ^ currently takes advantage of a whole node on cedar, may not need as much with lower amount of samples

#using the old version as the newer one kept throwing errors
module load StdEnv/2020  gcc/9.3.0
module load sra-toolkit/2.10.8

mkdir -p SRAs
cd SRAs

downloadSRA() {
    # Download Script
    prefetch ${1}
    fasterq-dump ${1} -O ../Raw_Fastqs                   
    gzip ../Raw_Fastqs/${1}*        #gzip extracted fastq(s)
    echo "${1} complete"
}

export -f downloadSRA

# if your job runs out of time before finishing, add --resume before --joblog and re-run 
parallel -j $SLURM_CPUS_PER_TASK --joblog sras.log downloadSRA {} :::: ./SRR_Acc_List.txt


