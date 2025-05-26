#!/bin/bash
#SBATCH --job-name=cell_rang_count_u0962361           # Job name
#SBATCH --output=cell_rang_count_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (cores)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=64G                      # Memory required per node (adjust based on your data size)
#SBATCH --time=08:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account

# Load Cell Ranger module (if required by your cluster)
module load cellranger

# Define paths
FASTQ_DIR="/scratch/general/vast/u0962361/working/fastq_folder/untrimmed_data/SRR24289600"
REFERENCE_DIR="/scratch/general/vast/u0962361/working/hg38_ref_cr"
SAMPLE="SRR24289600"

# Run Cell Ranger count command
cellranger count --id=${SAMPLE}_analysis \
                 --transcriptome=${REFERENCE_DIR} \
                 --fastqs=${FASTQ_DIR} \
                 --sample=${SAMPLE} \
                 --create-bam=true
                 --chemistry=SC3Pv3
                 