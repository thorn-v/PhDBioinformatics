#! /usr/bin/env bash

##### Usage/Options Block #####
usage() { printf 'Varient Calling Pipleine V1.3
        USAGE

        Downloads SRA files (From NCBI), extracts reads, and compresses them for further processing.

        -s\tSRA sample accesson number (from NCBI, should already have been downloaded) [REQUIRED]
        -q\tQuality cutoff (deafult 30)
        -l\Length cutoff (default 30)
        -m\tMinimum number of trimmed reads needed to keep  
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

## Default Values
qual=30
len=30
readsCutOff=0

while getopts "s:q:l:m:h" arg; do
        case $arg in
                s)
                        ACC=${OPTARG}
                        ;;
                q)
                        qual=${OPTARG}
                        ;;
                l)
                        len=${OPTARG}
                        ;;
                m)
                        readsCutOff=${OPTARG}
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

##### End Usage/Options Block #####

### Sanity Checks ###
if [[ -z "${ACC}" ]]; then
        printf '\nMissing required input: -s\nPlease provide input SRA Accession\n\nUse -h for usage help\n'
        exit 1;
fi
### ### ###

### Find R1 and R2 -----
# find the newly created files since they could have different ways of denoting R1/R2
cd Raw_Fastqs/${ACC}

fileArray=($ACC*)
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
cd ../..

###------

### Fast QC ------

mkdir -p Trimmed/${ACC}
mkdir -p Fastp_Logs # will only make it if not already made (so needed for the first one)

if [ "$R2" == "NA" ]; then # If not a paired sample...

        fastp -i Raw_Fastqs/${ACC}/${R1} \
                --out1 Trimmed/${ACC}/${ACC}_R1_trimmed.fastq.gz \
                --low_complexity_filter \
                -q ${qual} --cut_right --cut_front \
                --length_required $len \
                --html Fastp_Logs/${ACC}.html \
                --json Fastp_Logs/${ACC}.json;

else  # is a paired sample
		
        fastp -i ${Raw_Fastqs}/${ACC}/${R1} \
                -I ${Raw_Fastqs}/${ACC}/${R2} -w ${ncores} \
		--out1 Trimmed/${ACC}Trimmed/${ACC}_R1_trimmed.fastq.gz \
		--out2 Trimmed/${ACC}Trimmed/${ACC}_R2_trimmed.fastq.gz \
                --detect_adapter_for_pe --low_complexity_filter \
                --cut_right --cut_front -q $qual \
                --length_required $len \
                --html Fastp_Logs/${ACC}.html \
                --json Fastp_Logs/${ACC}.json; 
               # --unpaired1 ${ACC}Trimmed/${sample}_u1.fastq \
               # --unpaired2 ${ACC}Trimmed/${sample}_u2.fastq\
               # --failed_ACC ${ACC}FailedQC/${sample}_failed.fastq; 
               # don't care abACC these for this purpose currently ^
fi

### add quality check here to yeet files with too few reads
export REMOVED=0
if [[ ${readsCutOff} -gt 0 ]]; then
        passed_read_num=$(cat Fastp_Logs/${ACC}.json |\
                python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["filtering_result"]["passed_filter_reads"])')

        if [[ $passed_read_num -lt ${readsCutOff} ]]; then
                printf "${ACC}\tnumReads\t${passed_read_num}\tcontained less than ${readsCutOff} reads that passed QC\t`date +"%Y-%m-%d %T"`\n" |\
                tee -a removedAccessions.txt
                $REMOVED=1
                printf "Accession removed from evalutation\nExiting Script!\n"
                exit 0
        fi
fi

echo "Done QC!"

