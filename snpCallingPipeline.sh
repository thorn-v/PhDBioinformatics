#! /usr/bin/env bash
# Pipeline to go from sra download to snp vcf or bcf to efficently use space


usage() { printf 'Varient Calling Pipleine V1
        USAGE
        
        Assumes working on Compute Canada cluster!

        Takes downloaded .sra files and decompresses them,
        QC check them with fastp,
        maps passed reads to the provided reference libary (assumes library has been established/indexed),
        calls varients,
        filters varients into SNPs or Indels. 

        -i\tThe SRA accesson (folder that contains the Raw sequencing files) [REQUIRED]
        -r\tIndexed Sequence Reference library Folder path [REQUIRED]
	-n\tNumber of CPU Threads to be used (Default: 8)
        -q\tMinimum Mapping Quality (Default: 30)
	-l\tMinimum Read Length (Default: 30)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }


# Default Values
ncores=8
len=30
qual=30


while getopts "i:r:l:n:o:h" arg; do
        case $arg in
                i)
                # if the path to the folder is given with the / at the end, it can account for that now.
                        if [[ ${OPTARG: -1} == "/" ]]; then
                                sample="${OPTARG%/*}"
                                path="${OPTARG%/*}"
                                out="${path##*/}"
                        else
                                sample=${OPTARG} #sample is the path to the accession folder
                                out="${OPTARG##*/}"     #out is just the accession
                        fi
                        printf "sample is ${sample} and\n out is ${out}"
                        ;;
                r)
                        declare -r ref=${OPTARG}
                        ;;
                o)
                        out=${OPTARG} # legacy, may incorporate later - may remove, not advertised as an option
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        ;;
                q)
                        declare -i qual=${OPTARG} #may want to add a check for int after
                        ;;
                l)
                        len=${OPTARG}
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

### Check for missing or wrong inputs ###

if [[ -z "${sample}" ]]; then
        printf '\nMissing required input - Please provide input SRA download folder\n\nUse -h for usage help\n'
        exit 1;
fi

if [[ ! -d "${sample}" ]]; then #checks that the sample is a folder and not the .sra file
        printf '\nPlease provide a Directory for input SRA\n\nUse -h for usage help\n'
        exit 1;
fi

if [[ -z "${ref}" ]]; then

        printf '\nMissing required input - Please provide Reference Genome Library\n\nUse -h for usage help\n'
        exit 1;
fi



# # this step shoudnt need to happen now
#### Extract the fastq(s) ####
# # load the libraries needed
# module load StdEnv/2020
# module load gcc/9.3.0
# module load sra-toolkit

# ## put in a check for "if the folder already exists, skip this step"
# if [[ ! -d "~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}" ]]; then
#         fasterq-dump ${sample} -O ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out} 
#         #PID=$!

#         #wait "${PID}" #cannot move on until the sra is unpacked
#         echo "done unpacking"
# fi

# echo "exists now"
cd ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out} # move into the newley created directory with the fastq file(s)

# find the newly created files since they could have different ways of denoting r1/r2

# now we should have r1 and r2 (could be NA)

fileArray=($out*)
    # echo "${fileArray[@]}"
    # Identifying the files since they can have different endings when dumped
if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "r1\.f*|_1\.f*|_r1_0.*|_1"; then
        r1=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r1\.f*|_1\.f*|_r1_0.*|_1')
elif [[ -e "${out}.fastq.gz" ]]; then #for a merged file, it gets treated as r1
        r1="${out}.fastq.gz"
else
        r1="NA"
fi

# for r2
if printf '%s\n' "${fileArray[@]}" | grep -E -i -q "r2\.f*|_2\.f*|_r2_0.*|_2$"; then
        r2=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r2\.f*|_2\.f*|_r2_0.*|_2$')
else
        r2="NA"     
fi
#echo $r2
unset fileArray


######### QC and Trimming ###########
## QC and Trims the fastq file(s) for adapters, complexity, and quality windows - similar to timmomatic
## default quality cutoff is 20 

module load fastp
mkdir ~/scratch/Temp_WD/${out}Trimmed
cd ~/scratch/Temp_WD

if [ "$r2" == "NA" ]; then # If not a paired sample...

        fastp -i ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}/${r1} -w ${ncores} \
                --out1 ${out}Trimmed/${out}_r1_trimmed.fastq \
                --low_complexity_filter \
                -q ${qual} --cut_right --cut_front \
                --length_required $len \
                --html ~/scratch/Fastp_Logs/${out}.html \
                --json ~/scratch/Fastp_Logs/${out}.json;

else  # is a paired sample
		
        fastp -i ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}/${r1} \
                -I ~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}/${r2} -w ${ncores} \
		--out1 ${out}Trimmed/${out}_r1_trimmed.fastq \
		--out2 ${out}Trimmed/${out}_r2_trimmed.fastq \
                --detect_adapter_for_pe --low_complexity_filter \
                --cut_right --cut_front -q $qual \
                --length_required $len \
                --html ~/scratch/Fastp_Logs/${out}.html \
                --json ~/scratch/Fastp_Logs/${out}.json; 
               # --unpaired1 ${out}Trimmed/${sample}_u1.fastq \
               # --unpaired2 ${out}Trimmed/${sample}_u2.fastq\
               # --failed_out ${out}FailedQC/${sample}_failed.fastq;
fi

### add quality check here to yeet files with too few reads
passed_read_num=$(cat ~/scratch/Fastp_Logs/${out}.json |\
        python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["filtering_result"]["passed_filter_reads"])')

if [[ $passed_read_num -lt 100000 ]]; then
        printf "${out}\tnumReads\t${passed_read_num}\tcontained less than 100K reads that passed QC\n" |\
        tee -a ~/scratch/removedAccessions.txt
        echo "Exiting Script!"
        exit 0
fi

echo "done QC"

###### Mapping ######

module load bwa
module load samtools
mkdir ${out}MappedReads

if [[ "${r1}" == "NA" && "${r2}" == "NA" ]]; then  
        printf "${sample}\tno reads files found - manually check it out ($(date +'%d/%m/%y %H+3:%M:%S'))\n" | tee -a ~/scratch/missingFiles.txt
        exit 0;
fi


if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), send unincluded to null
        bwa mem ${ref} \
                ${out}Trimmed/${out}_r1_trimmed.fastq \
                -t ${ncores} -R "@RG\tID:${out}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - > ${out}MappedReads/${out}_mapped.bam 

fi

if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # if paird

        bwa mem ${ref} ${out}Trimmed/${out}_r1_trimmed.fastq \
                ${out}Trimmed/${out}_r2_trimmed.fastq \
                -t ${ncores}  -R "@RG\tID:${out}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - > ${out}MappedReads/${out}_mapped.bam

fi

echo "Done Mapping"

##### Processing ######

cd ${out}MappedReads
module load picard

java -jar $EBROOTPICARD/picard.jar SortSam \
      -I ${out}_mapped.bam \
      -O ${out}_sorted.bam \
      -SORT_ORDER coordinate

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      -I ${out}_sorted.bam \
      -O ${out}_sorted-md.bam \
      -M ${out}-md_metrics.txt

samtools index ${out}_sorted-md.bam

## second quality checkpoint - if the read depth is too low there is no point in continuing 
genomeReadDepth=$(samtools depth -a ${out}_mapped.bam | awk '{sum+=$3} END {print sum/NR}')

### JP suggests 50x

if [[ $genomeReadDepth -lt 50 ]]; then
        printf "${out}\tgenomeDepth\t${genomeReadDepth}\tAverage read depth for inital map is less than 50\n" |\
        tee -a ~/scratch/removedAccessions.txt
        echo "Exiting Script, Qual Too Low!"
        exit 0
else
        printf "Average sequence depth for ${out} is ${initialReadDepth}" |\
        tee -a ~/scratch/gvcfLogs/${out}_info.log
        mv ${out}_sorted-md.bam* FinalMappedReads
        mv ${out}-md_metrics.txt FinalMappedReads

fi


###### Varient Calling ########
cd ~/scratch/FinalMappedReads
module load gatk

gatk --java-options "-Xmx4g" HaplotypeCaller -ERC GVCF \
        -R ../Afumigatus_Reference/A_fumigatus_Af293/GCA_000002655.1/GCA_000002655.1_ASM265v1_genomic.fna \
        -I SRR11425549_sorted-md.bam \
        -O ~/scratch/GVCFs/SRR11425549.g.vcf 

echo "script finished!"
