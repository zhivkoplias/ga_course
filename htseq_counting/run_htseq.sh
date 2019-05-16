#!/bin/bash -l

#SBATCH -A g2019003
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J htseq_counting
#SBATCH --mail-type=ALL
#SBATCH --mail-user homo.korvin@gmail.com

# Load modules
module load bioinfo-tools
module load htseq

for i in ~/transcriptomics_data/aligned_reads/*_sorted.bam 
    do
    echo "Indexing: "$i
    htseq-count $i PROKKA_04122019.gtf -f bam -t CDS > ~/ga_course/htseq_counting/$i_htseq_out.txt
done

