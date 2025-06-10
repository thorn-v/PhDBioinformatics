#! /usr/bin/env bash

#SBATCH --job-name=parallel_gvcf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 4 #needs to be the same amonut of chromosomes
#SBATCH --mem-per-cpu=8G #Needs to be 4G abosulte min, realistically should be around 8G+ depending on depths.
#SBATCH --time=300
## can add  #SBATCH --mail-type=ALL  if you want it to email you updates

## Default Values
ploidy=1
chromList=chroms.list

while getopts "a:r:h" arg; do
        case $arg in
                a)
                        ACC=${OPTARG}
                        ;;
                r)
                        REF=${OPTARG}
                        ;;
                h | *)
                        echo "Make sure -r and -a are specified"
                        exit 0
                        ;;
        esac
done


javamem=$((SLURM_MEM_PER_NODE-2048)) 
#need to make sure there will be enough memory -
#SLURM_MEM_PER_NODE reports in Mb so equivelent to total memory allocated less 1Gb to make sure enough extra space for Java. 

###### Varient Calling ########
mkdir -p GVCFs

#TODO change this to scatter/gather option with SortVCF

scatter(){
# assumes gatk is in bin, if its not, manually set path here. !!!
gatk --java-options "-Xmx${javamem}g" HaplotypeCaller \
        -ERC GVCF \
        -L ${1} \
        -ploidy ${ploidy} \
        -R ${3} \
        -I MappedReads/${2}/${2}_sorted-md.bam \
        -O GVCFs/${2}/${2}_${1}.g.vcf.gz 
}
export -f scatter

parallel -j $SLURM_CPUS_PER_TASK --joblog comnine_${ACC}.log scatter {} ${ACC} ${REF} :::: ./chromList


#===========================
#     Gather
#===========================



