#! /usr/bin/env bash

##### Usage/Options Block #####
usage() { printf 'Varient Calling Pipleine V1.3
        USAGE

        Downloads SRA files (From NCBI), extracts reads, and compresses them for further processing.

        -s\tSRA sample accesson number (from NCBI, should already have been downloaded) [REQUIRED]
        -r\tPath to Reference sequence fasta [REQUIRED]
        -L\tProvide Path to GATK formatted interval list (defaults to generating one from reference for whole inteval)
        -A\tUse this flag if running on Alliance (Compute Canada)
        -p\tPloidy [num] (Defaults to 1 for haploid)
        -m\tTotal memory available to use (defults to 4G, java max memory will be m-2)  
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

## Default Values
ncores=1
ploidy=1
mem=4

while getopts "s:r:L:A:p:j:h" arg; do
        case $arg in
                s)
                        ACC=${OPTARG}
                        ;;
                r)
                        REF=${OPTARG}
                        ;;
                L)
                        chromList=${OPTARG}
                        ;;
                A)
                        computeCan="x"
                        ;;
                p)
                        ploidy=${OPTARG}
                        ;;
                m)
                        mem=${OPTARG}
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
        printf "\nReference file: ${REF} cannot be found\nPlease provide path to reference fasta\n\nUse -h for usage help\n"
        exit 1;
fi

if [[ ! -z ${chromList+x}]]; then #if flag set
    if [[ ! -e "${chromList}" ]]; then #if file does not exist
            printf "\nProvided interval file: ${REF} cannot be found\nPlease provide path to file or omit to generate file\n\nUse -h for usage help\n"
            exit 1;
    fi
fi
### ### ###

##### Make Chromosome List ######

if [[ -z ${chromList+x} ]]; then   #if $chromList does not exist (was not provided set), then ${chromList+x} will evaluate to nothing making statement TRUE. in that case, we attempt to make chroms.list file.     
        
        #should probably check we don't overwrite files here.
        grep "^>" ${REF} | cut -d " " -f 1 | sed -e 's/>//g' > chroms.list
        chromList="chroms.list"

fi #no need for else, if user supplied path then thats what we will use instead


###### Varient Calling ########
mkdir -p GVCFs
javamem=$((${mem}-2)) #need to make sure there will be enough memory

if [[ ! -z ${computeCan+x} ]]; then  #if $computeCan exists (was set) because user used flag, then will turn into "x" making the string not empty, making the if statement TRUE   
        
        module load gatk
        export JAVA_TOOL_OPTIONS="-Xmx${javamem}g"

fi

# assumes gatk is in bin, if its not, manually set path here. !!!
gatk --java-options "-Xmx${javamem}g" HaplotypeCaller \
        -ERC GVCF \
        -ploidy ${ploidy} \
        -R ${REF} \
        -I MappedReads/${ACC}/${ACC}_sorted-md.bam \
        -O GVCFs/${ACC}.g.vcf.gz 



# Ref info:

# for computeCan: export JAVA_TOOL_OPTIONS="-Xms256m -Xmx2g" ; recommend specifying the memory limit for your job as 1-2GB more than your setting on the Java command line option -Xmx
# paramters:
# -Xms<size>        set initial Java heap size.........................
# -Xmx<size>        set maximum Java heap size.........................

#${${VAR:+yes}:-no}
#if flag is set, picardCommand = yes, if not then = no in above example.
#fancy but confusing way vv
#picardCommand=${${computeCan:+'java -jar $EBROOTPICARD/picard.jar MarkDuplicates'}:-'java -jar picard.jar MarkDuplicates'}

