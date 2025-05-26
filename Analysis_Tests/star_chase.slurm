#!/bin/bash
#SBATCH --job-name=STAR_u0962361           # Job name
#SBATCH --output=STAR_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (cores)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=16G                      # Memory required per node (adjust based on your data size)
#SBATCH --time=04:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account

module load STAR

# Run STAR alignment
STAR --runThreadN 4 \
     --genomeDir /scratch/general/vast/u0962361/star_index \
     --readFilesIn /scratch/general/vast/u0962361/trimmed_data/sample_1_trimmed.fastq /scratch/general/vast/u0962361/trimmed_data/sample_2_trimmed.fastq \
     --outFileNamePrefix /scratch/general/vast/u0962361/alignment/aligned_output \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outReadsUnmapped Fastx \
     --quantMode GeneCounts \
     --outFilterMatchNmin 10 \
     --outFilterMatchNminOverLread 0.3
