#!/bin/bash
#SBATCH --job-name=fastp_trim_u0962361           # Job name
#SBATCH --output=fastp_trim_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (cores)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=16G                      # Memory required per node (adjust based on your data size)
#SBATCH --time=04:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account 

# Load fastp module (if available)
module load fastp

fastp -i SRR24289600_1.fastq -I SRR24289600_2.fastq \
      -o SRR24289600_1_t.fastq -O SRR24289600_2_t.fastq \
      --trim_poly_g --cut_tail --length_required 50 --quality_required 30 \
      --failed_out SRR24289600_1_failed.fastq SRR24289600_2_failed.fastq
