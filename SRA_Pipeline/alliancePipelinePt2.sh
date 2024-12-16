#! /usr/bin/env bash

#SBATCH --job-name=gvcf_pipeline_pt2
#SBATCH --ntasks=1
#SBATCH --mem=8G #if you change this, make sure you update the -m flag with the same number
#SBATCH --time=600 # can change this if too much time or not enough. Currently in mins

set -euxo pipefail #debugging aid
# Assumes you ran gvcf_pipeline_pt1 and are in the working directory from which you ran it

# !! Reference must be indexed with both bwa and samtools
ref=~/PATH/TO/REFERENCE.fa

## create sample.map ##
ls GVCFs |\
    sed -e 's/\.g\.vcf\.gz//' |\
    awk '{print $1"\tGVCFs/"$1".g.vcf.gz"}' > sample.map

# resets potential databse if, for example, job previously failed
rm -r wgs_database

## DB import database
./gatkDatabaseVarientCalling.sh -A -S sample.map -m 8

## Call full, Raw, VCF file
./callRawVCF.sh -A -r ${ref} -m 8
