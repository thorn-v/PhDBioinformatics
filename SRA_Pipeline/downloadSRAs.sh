#! /usr/bin/env bash

##### Usage/Options Block #####
usage() { printf 'Varient Calling Pipleine V1.3
        USAGE

        Downloads SRA files (From NCBI), extracts reads, and compresses them for further processing.

        -s\tSRA sample accesson number (from NCBI, max 20G) [REQUIRED]
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

while getopts "s:h" arg; do
        case $arg in
                s)
                        ACC=${OPTARG}
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
### Sanity Checks ###

# Download Script
prefetch ${ACC} -O SRAs/${ACC}
echo "done prefetch"

fasterq-dump SRA/${ACC} -O Raw_Fastqs/${ACC} 
echo "done unpack"
                
gzip Raw_Fastqs/${ACC}/*        #gzip extracted fastq(s)
echo "done gzip"
