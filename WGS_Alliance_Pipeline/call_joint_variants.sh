#! /usr/bin/env bash
#SBATCH --job-name=combined_calling
#SBATCH -A def-jpxu
#SBATCH --time=01:00:00   #can change
#SBATCH --mem=6G          #can change
## can add  #SBATCH --mail-type=ALL  if you want it to email you updates

## Default Values
ploidy=1
chromList=chroms.list
MAP=samples.map

#ref should be have gatk dict (have run housekeeping.sh)
REF=~/PATH/TO/REF #### <<<---- put path to ref here

#these values can be changed, and should if you re-run this script
databaseName=wgs_database  #change if you want; must not already exist
OUT=combined_raw_dataset   #can change if want

### Some Sanity Checks ###
if [[ ! -e "${chromList}" ]]; then #if file does not exist
        echo "Provided interval file: ${chromList} cannot be found\nPlease provide path to file."
        exit 1;
fi

if [[ ! -e "${MAP}" ]]; then
        echo "Reference file: ${MAP} cannot be found\nPlease provide path to sample map."
        exit 1;
fi

if [[ -e "${databaseName}" ]]; then
        echo "Database Directory already exists, GATK will not work\nPlease delete Diretory or provide new name."
        exit 1;
fi

## load gatk and its memory
module load gatk
javamem=$((SLURM_MEM_PER_NODE-2048)) #need to make sure there will be enough memory in mb
export JAVA_TOOL_OPTIONS="-Xmx${javamem}m"

#------------------------------------------------
#                Import Database
#------------------------------------------------
gatk --java-options "-Xmx${javamem}m" GenomicsDBImport \
       --genomicsdb-workspace-path ${databaseName} \
       -L ${chromList} \
       --sample-name-map ${MAP}

#------------------------------------------------
#                Call Combined VCF
#------------------------------------------------
gatk --java-options "-Xmx${javamem}m" GenotypeGVCFs \
       -R ${REF} \
       -V gendb://${databaseName} \
       -O ${OUT}.vcf.gz \
       -L $chromList \
       -ploidy $ploidy \
       -stand-call-conf 20


