#! /usr/bin/env bash
## Helper functions and wrappers for tasks to unclog the main script

### fastp trimming, adapted from George Sam Long ###

Trimming(){ # Performing the trimming
	# requires r1, r2, sample, out
	local r1=$1
	local r2=$2
	local sample=$3
	local out=$4
	local jobs=$5
	local totalThreads=$6
	

	if [[ $jobs == 1 ]]; then
	    local threads=$totalThreads
	fi

	if [ "$r2" == "NA" ]; then # If not a paired sample...

        fastp -i $r1 \
        --out1 ${out}_${sample}_r1_trimmed.fastq \
        --low_complexity_filter --correction \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        --length_required $len\
        --html ${out}FastpLogs/${sample}.html \
        --json ${out}FastpLogs/${sample}.json -R $sample --thread $threads -R $sample \
        --failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
	else
		#echo "$sample is a PE sample"
        fastp -i $r1 -I $r2 --merge \
        --merged_out ${out}Trimmed/${sample}_merged.fastq.gz \
		--out1 ${out}Trimmed/${sample}_r1.fastq.gz \
		--out2 ${out}Trimmed/${sample}_r2.fastq.gz \
        --detect_adapter_for_pe --correction --low_complexity_filter \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3\
        --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3\
        --overlap_len_require 15 --length_required $len\
        --html ${out}FastpLogs/${sample}.html \
        --json ${out}FastpLogs/${sample}.json -R $sample --thread $threads -R $sample \
        --unpaired1 ${out}Trimmed/${sample}_u1.fastq.gz \
        --unpaired2 ${out}Trimmed/${sample}_u2.fastq.gz\
        --failed_out ${out}FailedQC/${sample}_failed.fastq.gz;
	fi
}

findExtension(){
    # assumes that you are in the same directory as the files such as after a fasterq-dump
    if [[ -z "$1" ]]; then
        echo "Requires sample name, please try again"
        exit 1;
    fi
    
	fileArray=($1*)
    # echo "${fileArray[@]}"
    # Identifying the files since they can have different endings when dumped
	if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "r1\.f*|_1\.f*|_r1_0.*|_1"; then
		fileName=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r1\.f*|_1\.f*|_r1_0.*|_1')
	        r1="$fileName"
            #echo $r1
		unset fileName
	fi

	if printf '%s\n' "${fileArray[@]}" | grep -E -i -q "r2\.f*|_2\.f*|_r2_0.*|2$"; then
		fileName=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r2\.f*|_2\.f*|_r2_0.*|2$')
	        r2="$fileName"
            #echo $r2
		unset fileName
    else
        r2="NA"     
	fi
    #echo $r2
    unset fileArray
}


