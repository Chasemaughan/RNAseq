#!/bin/bash
#SBATCH --job-name=STAR_u0962361           # Job name
#SBATCH --output=STAR_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (1 for indexing)
#SBATCH --cpus-per-task=16              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=64G                      # Increase memory to 32GB (adjust based on your data size)
#SBATCH --time 08:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account

# Load STAR module (if available)
module load star

# Run STAR genomeGenerate with 16 threads
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir star_index_2 \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile gencode.v38.annotation.fixed.no_chr.gtf \
     --sjdbOverhang 99