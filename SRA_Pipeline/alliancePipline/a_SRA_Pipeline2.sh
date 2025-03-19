#! /usr/bin/env bash
#SBATCH --job-name=SRA_Pipeline_2
#SBATCH -A def-jpxu
#SBATCH --time=01:00:00   #can change
#SBATCH --mem=6G          #can change

## Default Values
ploidy=1
chromList=chroms.list
MAP=samples.map
databaseName=wgs_database  #change if you want
OUT=combined_raw_dataset   #can change if want

if [[ ! -e "${chromList}" ]]; then #if file does not exist
        printf "\nProvided interval file: ${chromList} cannot be found\nPlease provide path to file\n"
        exit 1;
fi

if [[ ! -e "${MAP}" ]]; then
        printf "\nReference file: ${MAP} cannot be found\nPlease provide path to sample map\n"
        exit 1;
fi

module load gatk
javamem=$((SLURM_MEM_PER_NODE-2)) #need to make sure there will be enough memory
export JAVA_TOOL_OPTIONS="-Xmx${javamem}g"

#------------------------------------------------
#                Import Database
#------------------------------------------------
gatk --java-options "-Xmx${javamem}g" GenomicsDBImport \
       --genomicsdb-workspace-path ${databaseName} \
       -L ${chromList} \
       --sample-name-map ${MAP}

#------------------------------------------------
#                Call Combined VCF
#------------------------------------------------
gatk --java-options "-Xmx${javamem}g" GenotypeGVCFs \
       -R ${REF} \
       -V gendb://${databaseName} \
       -O ${OUT}.vcf.gz \
       -L $chromList \
       -ploidy $ploidy \
       -stand-call-conf 20


