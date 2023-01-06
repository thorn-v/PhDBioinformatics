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
        -r\tSequence Reference library Folder path [REQUIRED]
        -o\tThe output prefix (Default: inputfoldername)
	-n\tNumber of CPU Threads to be used (Default: 8)
	-l\tLog File Name (Default: $date)
        -q\tMinimum Mapping Quality (Default: 30)
	-k\tMinimum Read Length (Default: 30)
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
len=30
log="$(date +'%Y_%m_%d').log"
qual=30


while getopts "i:r:k:n:o:l:h" arg; do
        case $arg in
                i)
                        declare -r sample=${OPTARG}
                        out=${OPTARG}
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
                k)
                        len=${OPTARG}
                        ;;
                l)
                        log=${OPTARG}
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

fasterq-dump $sample -O "./${out}_Unpacked" &

PID=$!

wait "${PID}" #cannot move on until the sra is unpacked
echo "done unpacking"

cd "${out}_Unpacked" #move into the newley created directory with the fastq file(s)

pwd
# find the newly created files since they could have different ways of denoting r1/r2

#findExtension ${sample}
# now we should have r1 and r2 (could be NA) as global variables

# call fastp wrapper

fileArray=($sample*)
    # echo "${fileArray[@]}"
    # Identifying the files since they can have different endings when dumped
if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "r1\.f*|_1\.f*|_r1_0.*|_1"; then
        r1=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r1\.f*|_1\.f*|_r1_0.*|_1')
elif [[ -e "${sample}.fastq" ]]; then
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

if [ "$r2" == "NA" ]; then # If not a paired sample...

        fastp -i $r1 \
                --out1 ${out}Trimmed/${out}_r1_trimmed.fastq \
                --low_complexity_filter \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
                --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
                --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
                --length_required $len\
                --html ${out}FastpLogs/${sample}.html -R $sample\
                --json ${out}FastpLogs/${sample}.json -R $sample \
                --failed_out ${out}FailedQC/${sample}_failed.fastq;

else  # is a paired sample
		
        fastp -i $r1 -I $r2 \
		--out1 ${out}Trimmed/${sample}_r1_trimmed.fastq \
		--out2 ${out}Trimmed/${sample}_r2_trimmed.fastq \
                --detect_adapter_for_pe --correction --low_complexity_filter \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
                --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
                --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
                --overlap_len_require 15 --length_required $len\
                --html ${out}FastpLogs/${sample}.html -R $sample \
                --json ${out}FastpLogs/${sample}.json -R $sample \
                --unpaired1 ${out}Trimmed/${sample}_u1.fastq \
                --unpaired2 ${out}Trimmed/${sample}_u2.fastq\
                --failed_out ${out}FailedQC/${sample}_failed.fastq;
fi

echo "done QC"

###### Mapping ######

module load bwa
module load samtools

if [[ "${r1}" == "NA" && "${r2}" == "NA" ]]; then  
        printf "${sample}\tno reads files found - possibly merged?" | tee -a missingFiles.txt
        exit 1;
fi


if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), add unincluded to diff file 
        bwa mem ${ref} \
                ${out}Trimmed/${sample}_r1_trimmed.fastq \
                -t ${ncores} | samtools view -b -h -F 4 -m ${len} -q ${qual} -U tmp.bam |\
                samtools sort - > ${out}MappedReads/${sample}_mapped.bam 

        samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz
fi

if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # if paird

        bwa mem ${ref} ${out}Trimmed/${sample}_r1_trimmed.fastq \
                ${out}Trimmed/${sample}_r2_trimmed.fastq -t ${ncores} |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U tmp.bam |\
                samtools sort - > ${out}MappedReads/${sample}_mapped.bam

        samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null # deletes singleton readings

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

cd ..

###### Varient Calling ########
module load freebayes

freebayes-parallel \
   <(fasta_generate_regions.py ${ref}.fai 100000) ${ncores} \
   -f ${ref} -p 1 ${out}MappedReads/${sample}_sorted-md.bam > ${sample}.vcf

echo "done varient calling"
####### Filter VCF ##########

mkdir FilteredVCF

module load vcftools

# separate indels
vcftools --vcf ${sample}.vcf --keep-only-indels --recode --recode-INFO-all --out FilteredVCF/${sample}_indels-only.vcf
# separate SNPs
vcftools --vcf ${sample}.vcf --remove-indels --recode --recode-INFO-all --out FilteredVCF/${sample}_snps-only.vcf

echo "script finished!"