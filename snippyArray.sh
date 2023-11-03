#! /usr/bin/env bash

#SBATCH --job-name=snippyArray
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=180
#SBATCH --array=1-2317%100

#config file for accessions
config=$1
fastqsPath=$2
# get the accession from the list
acc=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

module load fastp
module load python
module load bioperl

fileArray=($acc)
    # echo "${fileArray[@]}"
    # Identifying the files since they can have different endings when dumped
if  printf '%s\n' "${fileArray[@]}" | grep -E -i -q  "r1\.f*|_1\.f*|_r1_0.*|_1"; then
        r1=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r1\.f*|_1\.f*|_r1_0.*|_1')
elif [[ -e "${acc}.fastq.gz" ]]; then #for a merged file, it gets treated as r1
        r1="${acc}.fastq.gz"
else
        r1="NA"
fi

# for r2
if printf '%s\n' "${fileArray[@]}" | grep -E -i -q "r2\.f*|_2\.f*|_r2_0.*|_2$"; then
        r2=$(printf '%s\n' "${fileArray[@]}" | grep -E -i 'r2\.f*|_2\.f*|_r2_0.*|_2$')
else
        r2="NA"     
fi


mkdir -p $(pwd)/Trimmed_fastqs/${acc}Trimmed
mkdir - $(pwd)/Fastp_Logs # will only make it if not already made (so needed for the first one)
cd $(pwd)

if [ "$r2" == "NA" ]; then # If not a paired sample...

        fastp -i ${fastqsPath}/${acc}/${r1} -w 8 \
                --out1 ${acc}Trimmed/${acc}_r1_trimmed.fastq.gz \
                --low_complexity_filter \
                -q 30 --cut_right --cut_front \
                --length_required 30 \
                --html $(pwd)/Fastp_Logs/${acc}.html \
                --jso $(pwd)/Fastp_Logs/${acc}.json;

else  # is a paired sample
		
        fastp -i ${fastqsPath}/${acc}/${r1} \
                -I ${fastqsPath}/${acc}/${r2} -w 8 \
		--out1 ${acc}Trimmed/${acc}_r1_trimmed.fastq.gz \
		--out2 ${acc}Trimmed/${acc}_r2_trimmed.fastq.gz \
                --detect_adapter_for_pe --low_complexity_filter \
                --cut_right --cut_front -q 30 \
                --length_required 30 \
                --htm $(pwd)/Fastp_Logs/${acc}.html \
                --jso $(pwd)/Fastp_Logs/${acc}.json; 
               # --unpaired1 ${acc}Trimmed/${sample}_u1.fastq \
               # --unpaired2 ${acc}Trimmed/${sample}_u2.fastq\
               # --failed_out ${acc}FailedQC/${sample}_failed.fastq; 
               # don't care about these for this purpose currently ^
fi

### add quality check here to yeet files with too few reads

passed_read_num=$(zcat $(pwd)/Fastp_Logs/${acc}.json |\
        python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["filtering_result"]["passed_filter_reads"])')

if [[ $passed_read_num -lt 1000000 ]]; then
        printf "${acc}\tnumReads\t${passed_read_num}\tcontained less than 1M reads that passed QC\n" |\
        tee -a $(pwd)/removedAccessions.txt
fi


if [ "$r2" == "NA" ]; then

    $HOME/scratch/snippy/bin/snippy --outdir SnippyVarients --se Trimmed_fastqs/${acc}Trimmed/${acc}_r1_trimmed.fastq.gz --rgid ${acc} --reference ~/scratch/Afumigatus290_Reference/GCF_000002655.1_ASM265v1_genomic.fna --mincov 30
else
    $HOME/scratch/snippy/bin/snippy --outdir SnippyVarients --R1 Trimmed_fastqs/${acc}Trimmed/${acc}_r1_trimmed.fastq.gz --R2 Trimmed_fastqs/${acc}Trimmed/${acc}_r2_trimmed.fastq.gz --rgid ${acc} --reference ~/scratch/Afumigatus290_Reference/GCF_000002655.1_ASM265v1_genomic.fna --mincov 30
fi


