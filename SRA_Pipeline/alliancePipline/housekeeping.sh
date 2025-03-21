#! /usr/bin/env bash

#SBATCH --job-name=housekeeping
#SBATCH -A def-jpxu
#SBATCH --time=00:30:00     #can change
#SBATCH --mem=2G #can change

REF=PATH/TO/FASTA/REF.fasta

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

