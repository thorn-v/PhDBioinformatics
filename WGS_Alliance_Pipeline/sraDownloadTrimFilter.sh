#! /usr/bin/env bash

#SBATCH --job-name=SRAsetup
#SBATCH -A def-jpxu
#SBATCH --cpus-per-task=192  #can/should change to take advantage of whole node or max samples is lower than 48 
#SBATCH --time=12:00:00     #can change
#SBATCH --mem-per-cpu=4000    #can change
## ^ currently takes advantage of a whole node on fir, may not need as much with lower amount of samples

## this script downloads sra, unpacks, and trims/filters and deletes the sra+raw reads to save space
## pair this with map_and_gvcf.sh for the second part

module load sra-toolkit
module load fastp

readsCutOff=150000 #5000*30 ~size*coverage
export readsCutOff

mkdir -p ../Trimmed

downloadSRA() {
    # Download Script
    mkdir -p "${1}"
    cd "${1}" || exit
    
    prefetch "${1}"
    fasterq-dump "${1}"                 

    fileArray=("${1}"*)
        # echo "${fileArray[@]}"
        # Identifying the files since they can have different endings when dumped
    if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "R1\.f*|_1\.f*|_R1_0.*|_1"; then
            R1=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'R1\.f*|_1\.f*|_R1_0.*|_1')
    elif [[ -e "${ACC}.fastq" ]]; then #for a merged file, it gets treated as R1
            R1="${ACC}.fastq"
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

    # Fastp
    if [ "$R2" == "NA" ]; then # If not a paired sample...

        fastp -i "${R1}" \
                --out1 "../Trimmed/${1}_R1_trimmed.fastq.gz" \
                --low_complexity_filter \
                -q 30 --cut_right --cut_front \
                --length_required 30 \
                --json ${1}.json 
    
    else  # is a paired sample

        fastp -i "${R1}" \
              -I "${R2}" \
		    --out1 ../Trimmed/"${1}"_R1_trimmed.fastq.gz \
		    --out2 ../Trimmed/"${1}"_R2_trimmed.fastq.gz \
            --detect_adapter_for_pe --low_complexity_filter \
            --cut_right --cut_front -q 30 \
            --length_required 30 \
            --json ${1}.json 
    
    fi

    # reads cut off
    passed_read_num=$(cat ${1}.json |\
                    python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["filtering_result"]["passed_filter_reads"])')
    printf "${1}\t${passed_read_num}\n" | tee -a ../filteredReadsCount.txt

    cd .. || exit

    if [[ $passed_read_num -lt ${readsCutOff} ]]; then
        echo "${1}" | tee -a removedAccessions.txt
    else
        ehco "${1}" | tee -a keptAccessions.txt
    fi    

    # remove old SRA and raw fastqs for space if fastp worked
    if [[ -s ../Trimmed/"${1}"_R1_trimmed.fastq.gz ]]; then
        rm -r "${1}*"
    fi

}

export -f downloadSRA

# if your job runs out of time before finishing, add --resume before --joblog and re-run 
parallel -j $SLURM_CPUS_PER_TASK --joblog sras.log downloadSRA {} :::: ./SRR_Acc_List.txt


