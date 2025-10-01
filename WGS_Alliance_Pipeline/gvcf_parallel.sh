#! /usr/bin/env bash

#SBATCH --job-name=parallel_gvcf
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task 14 #needs to be the same amonut of chromosomes / should match the number of chromosomes in the chroms.list
#SBATCH --mem-per-cpu=8G #Needs to be 4G abosulte min, realistically should be around 8G+ depending on depths.
#SBATCH --time=300
#SBATCH --array=1-??%100 # !! replace the ?? with the total number of samples you have (number of lines in your configfile). if you have less than 100 samples, remove the "%100"

#config file for accessions
#should be a list of unique parts of file names for each sample - one per line
config=~/scratch/PATH/TO/configfile.txt
# if your file is names something other than configfile.txt, then use your filename.

# get the accession from the list
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $config)
#if testing on one accession before an array, comment the above and uncomment the below, putting the accession # directly
#ACC="SRR12345"


#ref should be bwa index, samtools index, and gatk dict (have run housekeeping.sh)
#!! needs to have gatk dict run on it vv
REF=~/PATH/TO/REF #### <<<---- put path to ref here

## Default Values
ploidy=1
chromList=chroms.list
#if your chromosome list is named something else or in a different directory, fix it above

javamem=$((SLURM_MEM_PER_NODE-2000)) 
#need to make sure there will be enough memory -
#SLURM_MEM_PER_NODE reports in Mb so equivelent to total memory allocated less 1Gb to make sure enough extra space for Java. 

###### Varient Calling ########
mkdir -p GVCFs

module load gatk
module load picard
export JAVA_TOOL_OPTIONS="-Xmx${javamem}m"

#TODO change this to scatter/gather option with SortVCF

scatter(){
# assumes gatk is in bin, if its not, manually set path here. !!!
gatk --java-options "-Xmx${javamem}m" HaplotypeCaller \
        -ERC GVCF \
        -L ${1} \
        -ploidy ${ploidy} \
        -R ${3} \
        -I MappedReads/${2}_sorted-md.bam \
        -O GVCFs/temp/${2}_${1}.g.vcf.gz 
}
export -f scatter

parallel -j $SLURM_CPUS_PER_TASK --joblog combnine_${ACC}.log scatter {} ${ACC} ${REF} :::: $chromList

wait
#===========================
#     Gatherx
#===========================
readarray -t chroms < $chromList

# Join with spaces
IFS=' ' commandGuts="I=GVCFs/temp/${ACC}_${chroms[*]}.g.vcf.gz"

        # I=GVCFs/temp/${ACC}_${chroms[0]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[1]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[2]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[3]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[4]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[5]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[6]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[7]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[8]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[9]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[10]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[11]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[12]}.g.vcf.gz \
        # I=GVCFs/temp/${ACC}_${chroms[13]}.g.vcf.gz \

#currently hardcoded for 14 chromosomes (C neoformans)
java -jar $EBROOTPICARD/picard.jar SortVCF \
        ${commandGuts} \
        O="GVCFs/${ACC}.g.vcf.gz"

if [[ -e GVCFs/${ACC}.g.vcf.gz ]]; then
        printf "${ACC}\tGVCFs/${ACC}.g.vcf.gz\n" | tee -a samples.map
        rm GVCFs/temp/${ACC}.*
else
        echo "Something went wrong while making the GVCF - Check it out!"
        exit 1
fi

