#!/bin/bash
#SBATCH --job-name=STAR2_u0962361           # Job name
#SBATCH --output=STAR2_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (1 for indexing)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=32G                      # Increase memory to 32GB (adjust based on your data size)
#SBATCH --time 08:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account
#SBATCH --email=u0962361@umail.utah.edu


# List of SRR IDs
SRR_IDS=("SRR24289600" "SRR24289601" "SRR24289602")

# Loop through each SRR ID and run fastq-dump
for SRR_ID in "${SRR_IDS[@]}"
do
    echo "Running fastq-dump for $SRR_ID"
    fastq-dump --split-files "$SRR_ID"
done
