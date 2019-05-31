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
module load htseq

for i in ~/tnseq_data/mapped_reads/*.bam
    do
    # echo "Quantification: "$i
    htseq-count $i PROKKA_04122019.gtf -f bam -t CDS > $i.transcripts.txt
done

