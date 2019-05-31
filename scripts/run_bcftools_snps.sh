#!/bin/bash -l

#SBATCH -A g2019003
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J tnseq_trimmo
#SBATCH --mail-type=ALL
#SBATCH --mail-user homo.korvin@gmail.com

# Load modules
module load bioinfo-tools
module load bcftools
module load python/2.7.15


for i in ~/mcgann/maping/*.sorted.bam
    bcftools mpileup $i -Ou -B -f ~/ga_course/final_assemblies/final_assembly_canu_and_spades_node3_120419.fasta --min-MQ 30 -o $i.vcf
    bcftools call $i.vcf -Ou -v -m -o $i.vcf_1
    bcftools norm $i.vcf_1 -Ou -f ~/ga_course/final_assemblies/final_assembly_canu_and_spades_node3_120419.fasta -d all -o $i.vcf_2
    bcftools filter $i.vcf_2 -Ov -e 'QUAL<40 || DP<20' -o $i.vcf_3
done

