
#!/bin/bash -l

#SBATCH -A g2019003
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -J first_canu
#SBATCH --mail-type=ALL
#SBATCH --mail-user homo.korvin@gmail.com

# Load modules
module load bioinfo-tools
module load canu

# Your commands
canu -d canu_030419_non_corrected -p 030419_non_corrected genomeSize=3.3M -pacbio-raw PacBio/*.fastq.gz > canu_run.log
