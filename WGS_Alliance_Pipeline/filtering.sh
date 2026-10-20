#! /usr/bin/env bash

#run this script with ./filtering.sh filename.vcf
in=$1 
out=${in%%_*}
# or just replace in and out with your file name inline vvv if you copy paste instead of running this as a script

## hard filtering
bcftools filter -e 'FS>60.0 || SOR>3 || MQ<40 || MQRankSum<-12.5 || QD<2.0 || ReadPosRankSum<-8.0' ${in} |\
bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' |\
bcftools filter -g 2 -o |\
bcftools +fill-tags -o "${out}.hf.vcf" -- -t F_MISSING,AC_Hom,AC_Het,VAF #adds extra info
#filter out any fixed sites and snps w/i 2bp of indel

## missing data per sample
total=$(grep -vc "^#" "${out}.hf.vcf")
bcftools stats -s - "${out}.hf.vcf" | grep -E ^PSC | cut -f3,14 |\
awk -v tot="$total" '{print $1"\t"$2"\t"$2/tot}' > ${out}_sample_missingness.txt

## remove samples provided 
#bcftools view -S ^${sampleList} -o ${out}_rm_missing.hf.vcf "${out}.hf.vcf"


