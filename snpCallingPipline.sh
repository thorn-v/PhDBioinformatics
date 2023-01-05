#! /usr/bin/env bash
# Pipeline to go from sra download to snp vcf or bcf to efficently use space

script_name=$0
script_full_path=$(dirname $0)
source ${script_full_path}/helperFunctions.sh

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
        -r\tSequence Reference library Folder [REQUIRED]
        -o\tThe output prefix (Default: inputfoldername_processed)
	-n\tNumber of CPU Threads to be used (Default: 8)
	-l\tLog File Name (Default: $date)
        -q\tMinimum Mapping Quality (Default: 30)
	-k\tMinimum Read Length (Default: 30)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

# Creating a simple command to save the settings used.
log() {	printf "Mapping settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

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
                        out= printf "${OPTARG}_processed"
                        ;;
                r)
                        declare -r reference=${OPTARG}
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

if [[ -z "${reference}" ]]; then
        printf '\nMissing required input - Please provide Reference Genome Library\n\nUse -h for usage help\n'
        exit 1;
fi



echo $sample
echo $reference   
#echo $ncores
#echo $len
#echo $log
# Writing the Log File
#log | tee $log

#### Extract the fastq(s) ####

fasterq-dump $sample -O "./${out}" &

PID=$!

wait "${PID}" #cannot move on until the sra is unpacked

cd "${out}" #move into the newley created directory with the fastq file(s)

# find the newly created files since they could have different ways of denoting r1/r2

findExtension ${sample}
# now we should have r1 and r2 (could be NA)

# call fastp wrapper

Trimming "${r1}" "${r2}" "${sample}" "${out}"


###### Mapping ######

if [[ "${r1}" == "NA" && "${r2}" == "NA" ]]; then  
        printf "${1}\tno r files found - possibly merged?" | tee -a missingFiles.txt
        exit 1;
fi


if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If it is not paired reads
        bwa mem $ref $r1 -t $ncores |\
                samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\ #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), add unincluded to diff file 
                samtools sort -o tmpSingle.bam 

        samtools fastq tmp.bam | gzip > ${out}UnmappedReads/${sample}_Single.fastq.gz
fi

if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # If I found a merged file

        bwa mem $ref $r1 $r2 -t $ncores |\
                samtools view -b -h -F 4 -m $len -q $qual -U tmp.bam |\ 
                samtools sort -o tmpP.bam

        samtools fastq -c 6 tmp.bam -1 ${out}UnmappedReads/${sample}_r1.fastq.gz -2 ${out}UnmappedReads/${sample}_r2.fastq.gz -s /dev/null

fi








