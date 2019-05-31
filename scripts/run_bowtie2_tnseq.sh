#!/bin/bash -l

#SBATCH -A g2019003
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 05:00:00
#SBATCH -J tnseq_trimmo
#SBATCH --mail-type=ALL
#SBATCH --mail-user homo.korvin@gmail.com

# Load modules
module load bioinfo-tools
module load bowtie2
module load samtools

for i in ~/tnseq_data/12bp/Tn*.12bp.gz
    do
    # echo "Mapping: "$i
    # bowtie2 -x ~/ga_course/final_assemblies/final_assembly_canu_and_spades_node3_120419 -U $i -S $i.sam
    # echo "Converting to bam: "$i
    # samtools view -S -b $i.sam > $i.bam
    # rm $i.sam
    echo "Sorting bam: "$i
    samtools sort $i.bam -o $i.sorted.bam
done

