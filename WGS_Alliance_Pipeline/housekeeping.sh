#! /usr/bin/env bash

#SBATCH --job-name=housekeeping
#SBATCH -A def-jpxu
#SBATCH --time=00:30:00     #can change
#SBATCH --mem=1G #can change

REF=PATH/TO/FASTA/REF.fasta
# Download the reference from NCBI. For example for A. fumigatus Af290:
# rsync -rv rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/ .

module load bwa
module load samtools
module load picard

bwa index ${REF}
samtools faidx ${REF}
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
      R=${REF} \
      O=${REF%.f*}.dict

# Create Chromosome List
grep "^>" ${REF} | cut -d " " -f 1 | sed -e 's/>//g' > chroms.list

