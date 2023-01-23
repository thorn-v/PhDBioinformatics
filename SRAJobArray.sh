#! /usr/bin/env bash

#SBATCH --job-name=sra_to_vcf_array_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=180
#SBATCH --array=1950,1951,1974,1976,1977,1979,1984,1985,1987,1988,1989,1990,1993,1997,2008,2030,2032,2033,2034,2035,2036,2038,2039,2044,2049,2050,2051,2056,2058,2059,2061,2062,2064,2065,2066,2068,2069,2070,2071,2076,2078,2079,2080,2081,2082,2085,2086,2092,2097,2104,2105,2107,2108,2109,2110,2111,2115,2118,2121,2124,2135,2149,2150,2151,2152,2153,2154,2155,2156,2169,2170,2172,2173,2174,2175,2176,2177,2180,2181,2184,2186,2188,2189,2190,2191

#config file for accessions
config=~/scratch/config.txt

# get the accession from the list
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
# echo "This is array task ${SLURM_ARRAY_TASK_ID} and the sample name is ${sample}." >> arrayOutput.txt

/home/vthorn/scratch/PhDBioinformatics_Code/snpCallingPipeline.sh -i ~/scratch/Afumigatus_WGSA_SRAs/${sample} -r ~/scratch/Afumigatus_Reference/A_fumigatus_Af293/GCA_000002655.1/GCA_000002655.1_ASM265v1_genomic.fna

rm -r ~/scratch/Temp_WD/"${sample}"*

