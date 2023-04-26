#! /usr/bin/env bash
# Pipeline to go from sra download to snp vcf or bcf to efficently use space


usage() { printf 'Varient Calling Pipleine V1.1
        USAGE
        
        Assumes working on Compute Canada cluster!

        Takes downloaded .sra files and decompresses them,
        QC check them with fastp,
        maps passed reads to the provided reference libary (assumes library has been established/indexed bwa and gatk dict/index),
        calls global varients. 

        Creates Temporary Working Directory ./Temp_WD

        -i\tThe SRA accesson (folder that contains the Raw sequencing files) [REQUIRED]
        -u\tUnpack SRA file to fastq(s)? (Default T; Accepts T/F)
        -r\tIndexed Sequence Reference library Folder path [REQUIRED]
        -f\tPath (full) to Fastqs [REQUIRED]
        -m\tMemory available (in Gigabytes, per core) (Defaults to 8G)
        -w\tWhole Genome filtering (Default T; Accepts T/F) (discards accessions with <10k qc passed reads and average WG coverage of <50) 
        -p)\tPloidy (Defaults to 1)
	-n\tNumber of CPU Threads to be used (Default: 1)
        -q\tMinimum Mapping Quality (Default: 30)
	-l\tMinimum Read Length (Default: 30)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }


# Default Values
ncores=1
len=30
qual=30
unpack="T"
workdir=$(pwd)
wgsFiltering="T"
mem=8
ploidy=1

while getopts "i:u:r:p:m:l:w:n:o:h:f" arg; do
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
                u)
                        unpack=${OPTARG}
                        ;;
                r)
                        ref=${OPTARG}
                        ;;
                f)
                        fastqsPath=${OPTARG}
                        ;;
                m)
                        mem=${OPTARG}
                        ;;
                o)
                        out=${OPTARG} # legacy, may incorporate later - may remove, not advertised as an option
                        ;;
                w)
                        wgsFiltering=${OPTARG} #should add check for T/F compliance
                        ;;
                n)
                        ncores=${OPTARG}
                        ;;
                q)
                        qual=${OPTARG} #may want to add a check for int after
                        ;;
                l)
                        len=${OPTARG}
                        ;;
                p)
                        ploidy=${OPTARG}
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

if [[ -z "${fastqsPath}" ]]; then

        printf '\nMissing required input - Please provide fastq library path (does not need to exist)\n\nUse -h for usage help\n'
        exit 1;
fi

#### Extract the fastq(s) ####
mkdir -p ${fastqsPath} # incase it does not already exist

# load some libraries needed
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit

if [[ ${unpack} == "T" ]]; then
        ## put in a check for "if the folder already exists (i.e. they messed up the -u), skip this step"
        if [[ ! -d "${fastqsPath}/${out}" || -z $(ls ${fastqsPath}/${out}) ]]; then
                fasterq-dump ${sample} -O ${fastqsPath}/${out} 
                echo "done unpack"
                gzip ${fastqsPath}/${out}/*
                echo "done gzip"
        fi
fi
# echo "exists now"
cd ${fastqsPath}/${out} # move into the newley created directory with the fastq file(s)

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
mkdir -p ${workdir}/Temp_WD/${out}Trimmed
mkdir -p ${workdir}/Fastp_Logs # will only make it if not already made (so needed for the first one)
cd ${workdir}/Temp_WD

if [ "$r2" == "NA" ]; then # If not a paired sample...

        fastp -i ${fastqsPath}/${out}/${r1} -w ${ncores} \
                --out1 ${out}Trimmed/${out}_r1_trimmed.fastq \
                --low_complexity_filter \
                -q ${qual} --cut_right --cut_front \
                --length_required $len \
                --html ${workdir}/Fastp_Logs/${out}.html \
                --json ${workdir}/Fastp_Logs/${out}.json;

else  # is a paired sample
		
        fastp -i ${fastqsPath}/${out}/${r1} \
                -I ${fastqsPath}/${out}/${r2} -w ${ncores} \
		--out1 ${out}Trimmed/${out}_r1_trimmed.fastq \
		--out2 ${out}Trimmed/${out}_r2_trimmed.fastq \
                --detect_adapter_for_pe --low_complexity_filter \
                --cut_right --cut_front -q $qual \
                --length_required $len \
                --html ${workdir}/Fastp_Logs/${out}.html \
                --json ${workdir}/Fastp_Logs/${out}.json; 
               # --unpaired1 ${out}Trimmed/${sample}_u1.fastq \
               # --unpaired2 ${out}Trimmed/${sample}_u2.fastq\
               # --failed_out ${out}FailedQC/${sample}_failed.fastq; 
               # don't care about these for this purpose currently ^
fi

### add quality check here to yeet files with too few reads
if [[ ${wgsFiltering} == "T" ]]; then
        passed_read_num=$(cat ${workdir}/Fastp_Logs/${out}.json |\
                python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["filtering_result"]["passed_filter_reads"])')

        if [[ $passed_read_num -lt 100000 ]]; then
                printf "${out}\tnumReads\t${passed_read_num}\tcontained less than 100K reads that passed QC\n" |\
                tee -a ${workdir}/removedAccessions.txt
                echo "Exiting Script!"
                exit 0
        fi
fi
echo "done QC"

###### Mapping ######

module load bwa
module load samtools
mkdir ${out}MappedReads

if [[ "${r1}" == "NA" && "${r2}" == "NA" ]]; then  
        printf "${sample}\tno reads files found - manually check it out ($(date +'%d/%m/%y %H+3:%M:%S'))\n" | tee -a ${workdir}/missingFiles.txt
        exit 0;
fi


if [[ "$r1" != "NA" && "$r2" == "NA" ]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), send unincluded to null
        bwa mem ${ref} \
                ${out}Trimmed/${out}_r1_trimmed.fastq \
                -t ${ncores} -R "@RG\tID:${out}\tSM:${out}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - > ${out}MappedReads/${out}_mapped.bam 

fi

if [[ "$r1" != "NA" && "$r2" != "NA" ]]; then # if paird

        bwa mem ${ref} ${out}Trimmed/${out}_r1_trimmed.fastq \
                ${out}Trimmed/${out}_r2_trimmed.fastq \
                -t ${ncores}  -R "@RG\tID:${out}\tSM:${out}" |\
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
# adding the +0.5 makes it round bc it will trunicate the number to int so rounds
genomeReadDepth=$(samtools depth -a ${out}_sorted-md.bam | awk '{sum+=$3} END {print int((sum/NR)+0.5)}')

mkdir -p ${workdir}/FinalMappedReads

if [[ ${wgsFiltering} == "T" ]]; then
        ### JP suggests 50x
        if [[ $genomeReadDepth -lt 50 ]]; then
                printf "${out}\tgenomeDepth\t${genomeReadDepth}\tAverage read depth for inital map is less than 50\n" |\
                tee -a ${workdir}/removedAccessions.txt
                echo "Exiting Script, Qual Too Low at ${genomeReadDepth}!"
                exit 0
        else
                printf "Average sequence depth for ${out} is ${genomeReadDepth}\n" |\
                tee -a ${workdir}/gvcfs_depths.info
                mv ${out}_sorted-md.bam* ${workdir}/FinalMappedReads
                mv ${out}-md_metrics.txt ${workdir}/FinalMappedReads
        fi

else # I still want to know the depts even if we don't use a cut-off
        printf "Average sequence depth for ${out} is ${genomeReadDepth}\n" |\
        tee -a ${workdir}/gvcfs_depts.info
        mv ${out}_sorted-md.bam* ${workdir}/FinalMappedReads
        mv ${out}-md_metrics.txt ${workdir}/FinalMappedReads   
fi


###### Varient Calling ########
mkdir -p ${workdir}/GVCFs
cd ${workdir}/FinalMappedReads
module load gatk
javamem=$((${mem}/2)) #need 

gatk --java-options "-Xmx${javamem}g" HaplotypeCaller \
        -ERC GVCF \
        -ploidy ${ploidy} \
        -G Standard \
        -G AS_Standard \
        -R ${ref} \
        -I ${out}_sorted-md.bam \
        -bamout ${out}_asmbl_hap.bam \
        -O ${workdir}/GVCFs/${out}.g.vcf 

echo "script finished! :)"
