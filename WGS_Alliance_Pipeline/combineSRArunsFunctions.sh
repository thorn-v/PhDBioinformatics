#! /usr/bin/env bash

#SBATCH --job-name=fixReadGroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16 
#SBATCH --mem-per-cpu=1000M #Needs to be 4G abosulte min, realistically should be around 8G+ depending on depths.
#SBATCH --time=240

# changes the read groups in all the relevant accs to the sampleid in the readgroup
# uncomment below and comment out "new script" section to use vvv

# module load samtools

# fixReadGroups() {

#     acc=$1
#     sample=$(awk -v acc="${acc}" '$1==$acc {Print $2}')
#     mv "MappedReads/${acc}_sorted-md.bam" "MappedReads/old_${acc}_sorted-md.bam"
    
#     samtools addreplacerg -r "ID:${sample}" -r "SM:${sample}" -o "MappedReads/${acc}_sorted-md.bam" "MappedReads/old_${acc}_sorted-md.bam"
#     samtools index "MappedReads/${acc}_sorted-md.bam"

# }

# export -f fixReadGroups

# parallel -j $SLURM_CPUS_PER_TASK --joblog fixReadGroups {} :::: mappedAccs.txt

#---------------------------------------------------------------------------
#                            new script
#---------------------------------------------------------------------------
# for if you have a 2 column sampleid \t acc1,acc2,acc3... input list
module load samtools
while read -r line; do
  sample=$(cut -f 1 $line)
  IFS=',' read -r -a accs <<< "$(cut -f 2 $line)"

  samtools merge -o "${sample}_sorted-md.bam" "${accs[@]/%/_sorted-md.bam}"
  samtools index "${sample}_sorted-md.bam"

    if [[ -s "${sample}_sorted-md.bam" ]]; then
        for i in "${accs[@]}"; do
            mv "${i}_sorted-md.bam" "old_${i}_sorted-md.bam"
        done
    fi
  unset accs

done <../multiAccs.txt




