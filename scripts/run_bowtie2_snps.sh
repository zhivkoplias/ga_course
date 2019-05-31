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

for i in ~/mcgann/maping/*_1_trimmed.fastq.gz
    do
    echo "Mapping: "$i
    R1=$(echo ${i#*_S})
    R2=${R1/_1_/_2_}
    echo "Mapping: "$R2
    bowtie2 -x Enterococcus_faecium.ASM76498v1.dna.toplevel -1 $i -2 $R2 -S $i.sam
    echo "Converting to bam: "$i
    samtools view -S -b $i.sam > $i.bam
    rm $i.sam
    echo "Sorting bam: "$i
    samtools sort $i.bam -o $i.sorted.bam
done
