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

# Set variables
GENOME_DIR="star_index_2"   # Path to your STAR genome index directory
READS_1="SRR24289602_1_trimmed.fastq"  # Path to paired-end read 1 (subset)
READS_2="SRR24289602_2_trimmed.fastq"  # Path to paired-end read 2 (subset)
OUTPUT_DIR="star_output_trim_2"       # Directory where you want the output to go
OUTPUT_PREFIX="output_trim_2_"  # Prefix for output files

# Make output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run STAR alignment
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --genomeDir $GENOME_DIR \
     --readFilesIn $READS_1 $READS_2 \
     --outFileNamePrefix $OUTPUT_DIR/$OUTPUT_PREFIX \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outReadsUnmapped Fastx \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM NM MD \
     --limitBAMsortRAM 20000000000  # Limit the memory used for sorting (adjust if needed)

# Print completion message
echo "STAR alignment completed. Results stored in $OUTPUT_DIR.