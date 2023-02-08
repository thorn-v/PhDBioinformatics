#! /usr/bin/env bash
#SBATCH --job-name=failed_accs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=180

for sample in "DRR017559" "DRR237582" "SRR11785209" "SRR14584271" "SRR16539770" "SRR16944301" "SRR16944376" "SRR16944468"; do
    /home/vthorn/scratch/PhDBioinformatics_Code/snpCallingPipeline.sh -i ~/scratch/Afumigatus_WGSA_SRAs/${sample} -r ~/scratch/Afumigatus_Reference/A_fumigatus_Af293/GCA_000002655.1/GCA_000002655.1_ASM265v1_genomic.fna -n 4 
    rm -r ~/scratch/Temp_WD/"${sample}"*
done
