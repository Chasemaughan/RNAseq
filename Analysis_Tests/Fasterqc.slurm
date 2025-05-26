#!/bin/bash
#SBATCH --job-name=job_u0962361           # Job name
#SBATCH --output=fastq_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (1 for indexing)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=32G                      # Increase memory to 32GB (adjust based on your data size)
#SBATCH --time 08:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account
#SBATCH --email=u0962361@umail.utah.edu

# List of SRR IDs
SRR_IDS=("SRR24289600" "SRR24289601" "SRR24289602")

# Load FasterQC and GNU Parallel (assuming they are installed on the cluster)
module load fasterqc
module load parallel

# Loop through each SRR ID, create file names for paired fastq files, and run FasterQC in parallel
echo "Running FasterQC on all paired FASTQ files in parallel..."

# For each SRR ID, assuming there are paired files: SRR_ID_1.fastq and SRR_ID_2.fastq
for SRR_ID in "${SRR_IDS[@]}"
do
    # Assuming paired files are named as SRR_ID_1.fastq and SRR_ID_2.fastq
    fastq_file_1="${SRR_ID}_1.fastq"
    fastq_file_2="${SRR_ID}_2.fastq"
    
    # Check if both files exist
    if [[ -f "$fastq_file_1" && -f "$fastq_file_2" ]]; then
        # Run FasterQC for both files in parallel
        parallel -j 2 "fasterqc {}" ::: "$fastq_file_1" "$fastq_file_2"
    else
        echo "One or both FASTQ files for $SRR_ID are missing!"
    fi
done