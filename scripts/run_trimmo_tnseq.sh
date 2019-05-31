#!/bin/bash -l

#SBATCH -A g2019003
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J tnseq_trimmo
#SBATCH --mail-type=ALL
#SBATCH --mail-user homo.korvin@gmail.com

# Load modules
module load bioinfo-tools
module load trimmomatic

for i in ~/tnseq_data/*.fastq.gz
    do
    echo "Indexing: "$i
    trimmomatic SE -phred33 $i $i.12bp ILLUMINACLIP:/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE-2.fa:2:25:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:4:25 MINLEN:36 CROP:22 HEADCROP:6
done

