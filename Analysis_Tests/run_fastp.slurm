#!/bin/bash
#SBATCH --job-name=fastp_trim_u0962361           # Job name
#SBATCH --output=fastp_trim_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (cores)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=16G                      # Memory required per node (adjust based on your data size)
#SBATCH --time=04:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account 

 

# Load Fastp module (if your cluster has Fastp installed through a module system)
module load fastp

# Input and Output files
INPUT_R1="SRR24289603_1.fastq"
INPUT_R2="SRR24289603_2.fastq"
OUTPUT_R1="SRR24289603_1_trimmed.fastq"
OUTPUT_R2="SRR24289603_2_trimmed.fastq"
HTML_REPORT="fastp_report.html"
JSON_REPORT="fastp_report.json"

# Run Fastp without the low memory option
fastp -i $INPUT_R1 -I $INPUT_R2 \
      -o $OUTPUT_R1 -O $OUTPUT_R2 \
      --html $HTML_REPORT --json $JSON_REPORT

# End of SLURM script
