#! /usr/bin/env bash
# Pipeline to go from sra download to snp vcf or bcf to efficently use space

script_name=$0
script_full_path=$(dirname $0)
source ${script_full_path}/helperFunctions.sh

# These are the files and variables that will be needed
usage() { printf 'Varient Calling Pipleine V1
        USAGE
        
        Takes downloaded .sra files and decompresses them,
        QC check them with fastp (using Sam Longs QCModern.sh script),
        maps passed reads to the provided reference libary (assumes library has been established),
        calls varients,
        more ...

        Assumes you have fastp, BWA mem, sra-toolbox,    

        -i\tThe SRA accesson (folder that contains the Raw sequencing files) [REQUIRED]
        -r\tSequence Reference library Folder [REQUIRED]
        -o\tThe output prefix (Default: inputfoldername_processed)
	-n\tNumber of CPU Threads to be used (Default: 8)
	-l\tLog File Name (Default: $date)
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










