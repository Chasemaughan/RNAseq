#!/bin/bash
#SBATCH --job-name=CellRang_mkref_u0962361           # Job name
#SBATCH --output=CellRang_mkref_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (cores)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=64G                      # Memory required per node (adjust based on your data size)
#SBATCH --time=04:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account

# Load necessary modules (e.g., for Cell Ranger)
module load cellranger/9.0.1

# Run the cellranger mkref command with the desired parameters
cellranger mkref --genome=hg38_ref_cr \
                 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh38.104.gtf \
                 --memgb=64 \
                 --runThreadN=8