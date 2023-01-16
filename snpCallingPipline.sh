#! /usr/bin/env bash
# Pipeline to go from sra download to snp vcf or bcf to efficently use space

script_name=$0
script_full_path=$(dirname $0)

# These are the files and variables that will be needed
usage() { printf 'Varient Calling Pipleine V1
        USAGE
        
        Takes downloaded .sra files and decompresses them,
        QC check them with fastp,
        maps passed reads to the provided reference libary (assumes library has been established),
        calls varients,
        more ...

        Assumes you have fastp, BWA mem, sra-toolbox,    

        -i\tThe SRA accesson (folder that contains the Raw sequencing files) [REQUIRED]
        -r\tIndexed Sequence Reference library Folder path [REQUIRED]
	-n\tNumber of CPU Threads to be used (Default: 8)
        -q\tMinimum Mapping Quality (Default: 20)
	-l\tMinimum Read Length (Default: 20)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
#log() {	printf "Mapping settings for $(date):
#	Log File:\t${log}
#	Input folder:\t${folder}
#	Output folder:\t${out}
#	CPU Threads:\t${ncores}
#	-------------------------------------\n"; exit 0;
#}

## add logging feature later

#export -f Trimming
#export -f FastpWrapper
#export -f FileIdentificationInFunction
#export -f FileExtractionInFunction

# Default Values
ncores=8
len=20
qual=20


while getopts "i:r:l:n:o:h" arg; do
        case $arg in
                i)
                        declare -r sample=${OPTARG}
                        out="${OPTARG##*/}"
                        ;;
                r)
                        declare -r ref=${OPTARG}
                        ;;
                o)
                        out=${OPTARG}
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




#### Extract the fastq(s) ####
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit

fasterq-dump ${sample} -O "~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}" &

PID=$!

wait "${PID}" #cannot move on until the sra is unpacked
echo "done unpacking"

cd "~/scratch/Afumigatus_WGSA_Raw_Fastqs/${out}" #move into the newley created directory with the fastq file(s)

# find the newly created files since they could have different ways of denoting r1/r2

# now we should have r1 and r2 (could be NA)

fileArray=($sample*)
    # echo "${fileArray[@]}"
    # Identifying the files since they can have different endings when dumped
if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "r1\.f*|_1\.f*|_r1_0.*|_1"; then
        r1=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r1\.f*|_1\.f*|_r1_0.*|_1')
elif [[ -e "${sample}.fastq" ]]; then #for a merged file, it gets treated as r1
        r1="${sample}.fastq"
else
        r1="NA"
fi

# for r2
if printf '%s\n' "${fileArray[@]}" | grep -E -i -q "r2\.f*|_2\.f*|_r2_0.*|2$"; then
        r2=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r2\.f*|_2\.f*|_r2_0.*|2$')
else
        r2="NA"     
fi
#echo $r2
unset fileArray


######### QC and Trimming ###########
## QC and Trims the fastq file(s) based on Sam Long's parameters 

module load fastp
mkdir ~/scratch/Temp_WD/${out}Trimmed


if [ "$r2" == "NA" ]; then # If not a paired sample...

        fastp -i $r1 -w ${ncores} \
                --out1 ${out}Trimmed/${out}_r1_trimmed.fastq \
                --low_complexity_filter \
                -q ${qual} --cut_right --cut_front \
                --length_required $len \
                --html ~/scratch/Fastp_Logs/${out}.html \
                --json ~/scratch/Fastp_Logs/${out}.json;

else  # is a paired sample
		
        fastp -i $r1 -I $r2 -w ${ncores} \
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

echo "done QC"

###### Mapping ######

module load bwa
module load samtools
mkdir ${out}MappedReads
mkdir ${out}UnmappedReads

if [[ "${r1}" == "NA" && "${r2}" == "NA" ]]; then  
        printf "${sample}\tno reads files found - manually check it out ("$(date +'%d/%m/%y %H+3:%M:%S')")" | tee -a ~/scratch/missingFiles.txt
        exit 1;
fi


if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), add unincluded to diff file 
        bwa mem ${ref} \
                ${out}Trimmed/${out}_r1_trimmed.fastq \
                -t ${ncores} | samtools view -b -h -F 4 -m ${len} -q ${qual} -U tmp.bam |\
                samtools sort - > ${out}MappedReads/${out}_mapped.bam 

        samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${out}_Single.fastq.gz
fi

if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # if paird

        bwa mem ${ref} ${out}Trimmed/${out}_r1_trimmed.fastq \
                ${out}Trimmed/${out}_r2_trimmed.fastq -t ${ncores} |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U tmp.bam |\
                samtools sort - > ${out}MappedReads/${out}_mapped.bam

        samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${out}_r1.fastq.gz -2 ${out}UnmappedReads/${out}_r2.fastq.gz -s /dev/null # deletes singleton readings

fi

echo "Done Mapping"
##### Processing ######

cd ${out}MappedReads
module load picard
module load samtools

java -jar $EBROOTPICARD/picard.jar SortSam \
      I=${sample}_mapped.bam \
      O=${sample}_sorted.bam \
      SORT_ORDER=coordinate

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${sample}_sorted.bam \
      O=${sample}_sorted-md.bam \
      M=${sample}-md_metrics.txt

samtools index ${sample}_sorted-md.bam


###### Varient Calling ########
module load python
module load freebayes

freebayes-parallel \
   <(fasta_generate_regions.py ${ref}.fai 100000) ${ncores} \
   -f ${ref} -p 1 ${out}_sorted-md.bam > ~/scratch/Raw_VCFs/${out}.vcf

echo "done varient calling"
####### Filter VCF ##########

cd ~/scratch/Raw_VCFs
module load htslib
module load vcftools

# separate indels
vcftools --vcf ${out}.vcf --keep-only-indels --recode --recode-INFO-all --out ${out}_indels-only.vcf
# separate SNPs
vcftools --vcf ${out}.vcf --remove-indels --recode --recode-INFO-all --out ${out}_snps-only.vcf

echo "script finished!"