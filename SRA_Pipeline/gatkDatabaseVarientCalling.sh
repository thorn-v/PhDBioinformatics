#! /usr/bin/env bash

##### Usage/Options Block #####
usage() { printf 'GATK Database
        USAGE

        Once all the sample GVCFs have been called, combines them into GATK database for further processing.

        -S\tPath to sample map file [REQUIRED]
        -L\tProvide Path to GATK formatted interval list (defaults to generated one from initalGVCFsCalling.sh)
        -A\tUse this flag if running on Alliance (Compute Canada)
        -p\tPloidy [num] (Defaults to 1 for haploid)
        -m\tTotal memory available to use (defults to 4G, java max memory will be m-2)  
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

## Default Values
ncores=1
ploidy=1
mem=4
chromList=chroms.list

while getopts "S:r:L:A:p:j:h" arg; do
        case $arg in
                S)
                        MAP=${OPTARG}
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
if [[ -z "${MAP}" ]]; then
        printf '\nMissing required input: -s\nPlease provide input SRA Accession\n\nUse -h for usage help\n'
        exit 1;
fi

if [[ ! -e "${MAP}" ]]; then
        printf "\nReference file: ${ACC} cannot be found\nPlease provide path to sample map\n\nUse -h for usage help\n"
        exit 1;
fi

# if [[ -z "${REF}" ]]; then
#         printf '\nMissing required input: -r\nPlease provide path to reference fasta\n\nUse -h for usage help\n'
#         exit 1;
# fi

# if [[ ! -e "${REF}" ]]; then
#         printf "\nReference file: ${REF} cannot be found\nPlease provide path to reference fasta\n\nUse -h for usage help\n"
#         exit 1;
# fi

if [[ ! -z ${chromList+x}]]; then #if flag set
    if [[ ! -e "${chromList}" ]]; then #if file does not exist
            printf "\nProvided interval file: ${REF} cannot be found\nPlease provide path to file or omit to generate file\n\nUse -h for usage help\n"
            exit 1;
    fi
fi
### ### ###

javamem=$((${mem}-2)) #need to make sure there will be enough memory

if [[ ! -z ${computeCan+x} ]]; then  #if $computeCan exists (was set) because user used flag, then will turn into "x" making the string not empty, making the if statement TRUE   
        
        module load gatk
        export JAVA_TOOL_OPTIONS="-Xmx${javamem}g"

fi

##### gatk ######
gatk --java-options "-Xmx${javamem}g" GenomicsDBImport \
       --genomicsdb-workspace-path wgs_database \
       -L chroms.list \
       --sample-name-map samples.map #will be made in the overall script along with ensuring wgs_databse is non-existant


