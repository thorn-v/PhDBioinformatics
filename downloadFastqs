#! /usr/bin/env bash

if [ -z "$1" ]; then
    echo "No argument supplied, re-submit"
    exit 1;
fi

#### Extract the fastq(s) ####
# load the libraries needed
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit

## put in a check for "if the folder already exists, skip this step"
for out in $1; do
    if [[ ! -d "~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}" ]]; then
            fasterq-dump ~/Afumigatus_WGSA_SRAs/${out} -O ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out} 
            cd ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}
            gzip *
            echo "exists now"
    else 
        echo "${out} exists already"
    fi
done
echo "done"