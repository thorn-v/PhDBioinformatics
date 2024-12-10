#! /usr/bin/env bash

##### Usage/Options Block #####
usage() { printf 'Align and DeDupes reads
        USAGE

        Aligns using BWA MEM to reference, sorts, and de-dupes bam with Picard. 

        -s\tSRA sample accesson number (from NCBI, should already have been downloaded) [REQUIRED]
        -r\tPath to Reference sequence fasta [REQUIRED]
        -A\tUse this flag if running on Alliance (Compute Canada)
        -q\tQuality cutoff (deafult 30)
        -l\Length cutoff (default 30)
        -j\tNumber of cores (default 1)
        -m\tMinimum number of trimmed reads needed to keep sample (optional, sample kept by default) 
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

## Default Values
qual=30
len=30
ncores=1
readsCutOff=0

while getopts "s:r:A:q:l:j:m:h" arg; do
        case $arg in
                s)
                        ACC=${OPTARG}
                        ;;
                r)
                        REF=${OPTARG}
                        ;;
                A)
                        computeCan="x"
                        ;;
                q)
                        qual=${OPTARG}
                        ;;
                l)
                        len=${OPTARG}
                        ;;
                j)
                        ncores=${OPTARG}
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

if [[ -z "${REF}" ]]; then
        printf '\nMissing required input: -r\nPlease provide path to reference fasta\n\nUse -h for usage help\n'
        exit 1;
fi

if [[ ! -e "${REF}" ]]; then
        printf '\nReference file: ${REF} cannot be found\nPlease provide path to reference fasta\n\nUse -h for usage help\n'
        exit 1;
fi
### ### ###


###### Mapping ######

mkdir -p MappedReads
R1=Trimmed/${ACC}/${ACC}_R1_trimmed.fastq.gz
R2=Trimmed/${ACC}/${ACC}_R2_trimmed.fastq.gz

#${${VAR:+yes}:-no}
#if flag is set, picardCommand = yes, if not then = no in above example.
#fancy but confusing way vv
#picardCommand=${${computeCan:+'java -jar $EBROOTPICARD/picard.jar MarkDuplicates'}:-'java -jar picard.jar MarkDuplicates'}

if [[ -z ${computeCan+x} ]]; then   #if $computeCan exists (was set) because user used flag, then will turn into "x" making the string not empty, making the if statement FALSE    

    picardCommand='java -jar picard.jar MarkDuplicates'

else
   
    module load picard
    module load samtools
    module load bwa
    picardCommand='java -jar $EBROOTPICARD/picard.jar MarkDuplicates'

fi

if [[ ! -e "${R1}" && ! -e "${R2}" ]]; then  
        printf "Reads file ${R1} and ${R2}\tno reads files found - manually check it out ($(date +'%d/%m/%y %H+3:%M:%S'))\n"
        exit 1;
fi


if [[ ! -e "$R2"]]; then # If it is not paired reads
        #bam, w/ headers, exclude unmapped, include only greater len than $len (30 default), Skip alignments with MAPQ smaller than $qual (default 30), send unincluded to null
        bwa mem ${REF} ${R1} \
                -t ${ncores} -R "@RG\tID:${ACC}\tSM:${ACC}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - |\
                $picardCommand -I - \
                               -O MappedReads/${ACC}/${ACC}_sorted-md.bam \
                               -M MappedReads/${ACC}/${ACC}-md_metrics.txt

fi

if [[ -e "$R1" && -e "$R2" ]]; then # if paird

        bwa mem ${REF} ${R1} ${R2} \
                -t ${ncores}  -R "@RG\tID:${ACC}\tSM:${ACC}" |\
                samtools view -b -h -F 4 -m ${len} -q ${qual} -U /dev/null |\
                samtools sort - |\
                $picardCommand -I - \
                               -O MappedReads/${ACC}/${ACC}_sorted-md.bam \
                               -M MappedReads/${ACC}/${ACC}-md_metrics.txt

fi

samtools index MappedReads/${ACC}/${ACC}_sorted-md.bam

echo "Done Mapping / indexing"

##### Processing ######

## second quality checkpoint - if the read depth is too low there is no point in continuing 
# adding the +0.5 makes it round bc it will trunicate the number to int so rounds
genomeReadDepth=$(samtools depth -a MappedReads/${ACC}/${ACC}}_sorted-md.bam | awk '{sum+=$3} END {print int((sum/NR)+0.5)}')

printf "Average sequence depth for ${ACC} is ${genomeReadDepth} (`date +"%Y-%m-%d %T"`\n" |\
        tee -a gvcfs_depts.info

mkdir -p FinalMappedReads

if [[ ${wgsCutOff} -gt 0 ]]; then
        ### JP suggests 50x for mito
        if [[ $genomeReadDepth -lt ${wgsCutOff} ]]; then
                printf "${ACC}\tgenomeDepth\t${genomeReadDepth}\tAverage read depth for inital map is less than ${wgsCutOff}\t`date +"%Y-%m-%d %T"`\n" |\
                tee -a removedAccessions.txt
                echo "Exiting Script, Qual Too Low at ${genomeReadDepth}!"
                exit 0
        else
                mv ${ACC}_sorted-md.bam* FinalMappedReads/
                mv ${ACC}-md_metrics.txt FinalMappedReads/
        fi
