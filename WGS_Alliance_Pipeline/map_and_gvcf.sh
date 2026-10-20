#! /usr/bin/env bash

#SBATCH --job-name=map&gvcf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 4 #haplotype caller can use up to 4, bwa can use more if needed
#SBATCH --mem-per-cpu=8000M #Needs to be 4G abosulte min, realistically should be around 8G+ depending on depths.
#SBATCH --time=300
#SBATCH --array=1-??%100 # !! replace the ?? with the total number of samples you have (number of lines in your configfile). if you have less than 100 samples, remove the "%100"

## for use after download-trim-filter script

#config file for accessions
#should be a list of unique parts of file names for each sample - one per line
config=~/scratch/PATH/TO/configfile.txt
# if your file is names something other than configfile.txt, then use your filename.

# get the accession from the list
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $config)

#ref should be bwa index, samtools index, and gatk dict (have run housekeeping.sh)

REF=~/PATH/TO/REF #### <<<---- put path to ref here
chromList=chroms.list ### add chrom list here - should just have the chrom(s) or organelles you want to call for
mitoName="NameOfMitochondriaChrom" #change this to be your mito chrom NCBI ID (should be the first thing after > in the fasta)
wgsCutOff=0

## Variables ##
ploidy=1
qual=30
len=30

#---------------------------
#     Checks
#---------------------------
if [[ ! -e "${REF}" ]]; then
        echo "Reference file: ${REF} cannot be found, Please provide path to reference fasta"
        exit 1;
fi
if [[ -z "${ACC}" ]]; then
        echo "Accession: ${ACC} is blank, please try again."
        exit 1;
fi

#---------------------------------------------------------------------------
#                            Map & De-Dupe
#---------------------------------------------------------------------------

##### Mapping ######

mkdir -p MappedReads
R1=Trimmed/${ACC}_R1_trimmed.fastq.gz
R2=Trimmed/${ACC}_R2_trimmed.fastq.gz

module load picard
module load gatk
module load samtools
module load bwa

javamem=$((SLURM_MEM_PER_CPU*SLURM_CPUS_PER_TASK-2000)) #need to make sure there will be enough memory -
                                     #SLURM_MEM_PER_NODE reports in Mb so equivelent to total memory allocated less 2Gb to make sure enough extra space for Java. 
# TODO add check to make sure there is at least 2Gb of memory total
export JAVA_TOOL_OPTIONS="-Xmx${javamem}m"


if [[ ! -e "${R1}" && ! -e "${R2}" ]]; then  
        echo "Reads file ${R1} and ${R2}\tno reads files found - manually check it out ($(date +'%d/%m/%y %H+3:%M:%S'))"
        exit 1;
fi


if [[ ! -e "${R2}" ]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), send unincluded to null
        #need to add readgroups (-R) or gatk complains (unless already multi-lane)
        bwa mem "${REF}" "${R1}" -t ${SLURM_CPUS_PER_TASK} \
                -R "@RG\tID:${ACC}\tSM:${ACC}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort -o "MappedReads/${ACC}_sorted.bam" -
                
                java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
                               -I "MappedReads/${ACC}_sorted.bam" \
                               -O "MappedReads/${ACC}_sorted-md.bam" \
                               -M "MappedReads/${ACC}-md_metrics.txt"

fi

if [[ -e "$R1" && -e "$R2" ]]; then # if paird

        bwa mem "${REF}" "${R1}" "${R2}" -t ${SLURM_CPUS_PER_TASK} \
                -R "@RG\tID:${ACC}\tSM:${ACC}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort -o "MappedReads/${ACC}_sorted.bam" -
                
        java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
                        -I "MappedReads/${ACC}_sorted.bam" \
                        -O "MappedReads/${ACC}_sorted-md.bam" \
                        -M "MappedReads/${ACC}-md_metrics.txt"

fi

samtools index "MappedReads/${ACC}_sorted-md.bam"

echo "Done Mapping / indexing"

##### Processing ######

samtools coverage -H "MappedReads/${ACC}_sorted-md.bam" |\
        awk -v s=${ACC} '{print s "\t" $0}' >> combined_full_coverages.txt

## second quality checkpoint - if the read depth is too low there is no point in continuing 
# adding the +0.5 makes it round bc it will trunicate the number to int so rounds
mitoReadDepth=$(samtools depth -r ${mitoName} -a "MappedReads/${ACC}_sorted-md.bam" | awk '{sum+=$3} END {print int((sum/NR)+0.5)}')

echo "Average sequence depth for ${ACC} is ${mitoReadDepth} ($(date +"%Y-%m-%d %T"))" |\
        tee -a gvcfs_depts.info

###### Filtering #######

if [[ $mitoReadDepth -lt ${wgsCutOff} ]]; then
        printf "${ACC}\t${mitoReadDepth}\n" |\
        tee -a removedAccessions.txt
        echo "Exiting Script, Depth Too Low at ${mitoReadDepth}!"
        exit 0
else
        echo "${ACC}" | tee -a keptAccessions.txt
fi
#---------------------------------------------------------------------------
#                            Initial GVCF
#---------------------------------------------------------------------------

###### Varient Calling ########
mkdir -p GVCFs

#TODO change this to scatter/gather option with SortVCF from gvcf_parallel.sh

# assumes gatk is in bin, if its not, manually set path here. !!!
gatk --java-options "-Xmx${javamem}m" HaplotypeCaller \
        -ERC GVCF \
        -L ${chromList} \
        -ploidy ${ploidy} \
        -R ${REF} \
        -I MappedReads/${ACC}_sorted-md.bam \
        -O GVCFs/${ACC}.g.vcf.gz 

if [[ -e GVCFs/${ACC}.g.vcf.gz ]]; then
        printf "${ACC}\tGVCFs/${ACC}.g.vcf.gz\n" | tee -a samples.map #needs to be printf for the \t to work
else
        echo "Something went wrong while making the GVCF - Check it out!"
fi
