#! /usr/bin/env bash

#SBATCH --job-name=inital_array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 4
#SBATCH --mem=8G #Needs to be 4G abosulte min, realistically should be around 8G+ depending on depths.
#SBATCH --time=300
#SBATCH --array=1-??%100 # !! replace the ?? with the total number of samples you have (number of lines in your configfile). if you have less than 100 samples, remove the "%100"

#config file for accessions
#should be a list of unique parts of file names for each sample - one per line
config=~/scratch/PATH/TO/configfile.txt
# if your file is names something other than configfile.txt, then use your filename.

# get the accession from the list
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $config)

#ref should be bwa index, samtools index, and gatk dict (have run housekeeping.sh)

REF=~/PATH/TO/REF #### <<<---- put path to ref here

## Variables ##
ploidy=1
qual=30
len=30

#---------------------------
#     Checks
#---------------------------
if [[ ! -e "${REF}" ]]; then
        echo "Reference file: ${REF} cannot be found\nPlease provide path to reference fasta"
        exit 1;
fi


#===========================================================================
#                            QC and Trim
#===========================================================================

## Load Modules here ##

module load fastp

## ================= ##


### Find R1 and R2 -----
# find the newly created files since they could have different ways of denoting R1/R2
cd Raw_Fastqs

fileArray=("${ACC}"*)
    # echo "${fileArray[@]}"
    # Identifying the files since they can have different endings when dumped
if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "R1\.f*|_1\.f*|_R1_0.*|_1"; then
        R1=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'R1\.f*|_1\.f*|_R1_0.*|_1')
elif [[ -e "${ACC}.fastq.gz" ]]; then #for a merged file, it gets treated as R1
        R1="${ACC}.fastq.gz"
else
        R1="NA"
fi

# for R2
if printf '%s\n' "${fileArray[@]}" | grep -E -i -q "R2\.f*|_2\.f*|_R2_0.*|_2$"; then
        R2=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'R2\.f*|_2\.f*|_R2_0.*|_2$')
else
        R2="NA"     
fi

unset fileArray
cd ..

###------

### Fastp QC ------

mkdir -p Trimmed/"${ACC}"
mkdir -p Fastp_Logs # will only make it if not already made (so needed for the first one)

if [ "$R2" == "NA" ]; then # If not a paired sample...

        fastp -i Raw_Fastqs/"${ACC}"/"${R1}" \
                --out1 Trimmed/"${ACC}"/"${ACC}"_R1_trimmed.fastq.gz \
                --low_complexity_filter \
                -q ${qual} --cut_right --cut_front \
                --length_required $len \
                --html Fastp_Logs/${ACC}.html 
else  # is a paired sample
		
        fastp -i Raw_Fastqs/"${ACC}"/"${R1}" \
                -I Raw_Fastqs/"${ACC}"/"${R2}" \
		--out1 Trimmed/"${ACC}"/"${ACC}"_R1_trimmed.fastq.gz \
		--out2 Trimmed/"${ACC}"/"${ACC}"_R2_trimmed.fastq.gz \
                --detect_adapter_for_pe --low_complexity_filter \
                --cut_right --cut_front -q $qual \
                --length_required $len \
                --html Fastp_Logs/${ACC}.html 
fi


echo "QC Finished!"

#---------------------------------------------------------------------------
#                            Map & De-Dupe
#---------------------------------------------------------------------------

##### Mapping ######

mkdir -p MappedReads
R1=Trimmed/${ACC}/${ACC}_R1_trimmed.fastq.gz
R2=Trimmed/${ACC}/${ACC}_R2_trimmed.fastq.gz


module load picard
module load samtools
module load bwa

javamem=$((SLURM_MEM_PER_NODE-2048)) #need to make sure there will be enough memory -
                                     #SLURM_MEM_PER_NODE reports in Mb so equivelent to total memory allocated less 1Gb to make sure enough extra space for Java. 
# TODO add check to make sure there is at least 2Gb of memory total
export JAVA_TOOL_OPTIONS="-Xmx${javamem}m"


if [[ ! -e "${R1}" && ! -e "${R2}" ]]; then  
        echo "Reads file ${R1} and ${R2}\tno reads files found - manually check it out ($(date +'%d/%m/%y %H+3:%M:%S'))"
        exit 1;
fi


if [[ ! -e "${R2}" ]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), send unincluded to null
        #need to add readgroups (-R) or gatk complains (unless already multi-lane)
        bwa mem "${REF}" "${R1}" \
                -R "@RG\tID:${ACC}\tSM:${ACC}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - > MappedReads/${ACC}/${ACC}_sorted.bam
                
                java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
                               -I MappedReads/${ACC}/${ACC}_sorted.bam \
                               -O MappedReads/${ACC}/${ACC}_sorted-md.bam \
                               -M MappedReads/${ACC}/${ACC}-md_metrics.txt

fi

if [[ -e "$R1" && -e "$R2" ]]; then # if paird

        bwa mem "${REF}" "${R1}" "${R2}" \
                -R "@RG\tID:${ACC}\tSM:${ACC}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - > MappedReads/${ACC}/${ACC}_sorted.bam
                
                java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
                               -I MappedReads/${ACC}/${ACC}_sorted.bam \
                               -O MappedReads/${ACC}/${ACC}_sorted-md.bam \
                               -M MappedReads/${ACC}/${ACC}-md_metrics.txt

fi

samtools index MappedReads/${ACC}/${ACC}_sorted-md.bam

echo "Done Mapping / indexing"

##### Processing ######

## second quality checkpoint - if the read depth is too low there is no point in continuing 
# adding the +0.5 makes it round bc it will trunicate the number to int so rounds
genomeReadDepth=$(samtools depth -a MappedReads/${ACC}/${ACC}_sorted-md.bam | awk '{sum+=$3} END {print int((sum/NR)+0.5)}')

echo "Average sequence depth for ${ACC} is ${genomeReadDepth} ($(date +"%Y-%m-%d %T"))" |\
        tee -a gvcfs_depts.info

#---------------------------------------------------------------------------
#                            Initial GVCF
#---------------------------------------------------------------------------
module load gatk

###### Varient Calling ########
mkdir -p GVCFs

#TODO change this to scatter/gather option with SortVCF

# assumes gatk is in bin, if its not, manually set path here. !!!
gatk --java-options "-Xmx${javamem}g" HaplotypeCaller \
        -ERC GVCF \
        -L chroms.list \
        -ploidy ${ploidy} \
        -R ${REF} \
        -I MappedReads/${ACC}/${ACC}_sorted-md.bam \
        -O GVCFs/${ACC}.g.vcf.gz 

if [[ -e GVCFs/${ACC}.g.vcf.gz ]]; then
        printf "${ACC}\tGVCFs/${ACC}.g.vcf.gz\n" | tee -a samples.map
else
        echo "Something went wrong while making the GVCF - Check it out!"
        exit 1
fi
