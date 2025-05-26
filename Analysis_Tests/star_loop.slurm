#!/bin/bash
#SBATCH --job-name=STAR_u0962361           # Job name
#SBATCH --output=STAR_%A_%a.log            # Standard output and error log (%A is job ID, %a is array task ID)
#SBATCH --ntasks=1                        # Number of tasks (1 for alignment)
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=64G                         # Memory allocation
#SBATCH --time 08:00:00                   # Time limit
#SBATCH --partition=notchpeak-shared-short # Partition
#SBATCH --account=notchpeak-shared-short   # Account

# Load STAR module (if available)
module load star

# Set paths for STAR index and output directory
GENOME_DIR="star_index"
OUTPUT_DIR="star_output"

# Array of sample read pairs (adjust paths as needed)
READS_ARRAY=("SRR24289600" "SRR24289601" "SRR24289602")

# Set SLURM_ARRAY_TASK_ID to access the corresponding sample from the array
SAMPLE=${READS_ARRAY[$SLURM_ARRAY_TASK_ID-1]}  # Subtract 1 because SLURM_ARRAY_TASK_ID starts from 1

# Set read file paths for paired-end files
READS_1="${SAMPLE}_1.fastq"
READS_2="${SAMPLE}_2.fastq"

# Set output prefix
OUTPUT_PREFIX="output_${SAMPLE}"

# Run STAR for alignment
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --genomeDir $GENOME_DIR \
     --readFilesIn $READS_1 $READS_2 \
     --outFileNamePrefix $OUTPUT_DIR/$OUTPUT_PREFIX \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outReadsUnmapped Fastx \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM NM MD \
     --limitBAMsortRAM 20000000000