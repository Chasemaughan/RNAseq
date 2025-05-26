#!/bin/bash
#SBATCH --job-name=STAR_u0962361           # Job name
#SBATCH --output=STAR_%j.log      # Standard output and error log (%j will be replaced by job id)
#SBATCH --ntasks=1                     # Number of tasks (1 for indexing)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (adjust based on your cluster's resources)
#SBATCH --mem=64G                      # Increase memory to 32GB (adjust based on your data size)
#SBATCH --time=04:00:00                # Time limit (adjust as necessary)
#SBATCH --partition=notchpeak-shared-short # Specify partition
#SBATCH --account=notchpeak-shared-short    # Specify account

echo "Hello, SLURM!"

# Load required modules (you might already have them)
module load star

# Define paths for genome files and output directory
GENOME_FASTA="/scratch/general/vast/u0962361/SRR24289603/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="/scratch/general/vast/u0962361/SRR24289603/Homo_sapiens.GRCh38.104.gtf"
OUTPUT_DIR="/scratch/general/vast/u0962361/SRR24289603/star_index2"

# Check if the genome FASTA and GTF files exist
if [[ ! -f "$GENOME_FASTA" ]]; then
  echo "Genome FASTA file not found!"ssssssssssssssssss
  exit 1
fi

if [[ ! -f "$GTF_FILE" ]]; then
  echo "GTF file not found!"
  exit 1
fi

# Check if the output directory exists, if not create it
mkdir -p $OUTPUT_DIR

# Run STAR genome indexing command
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $OUTPUT_DIR \
     --genomeFastaFiles $GENOME_FASTA \
     --sjdbGTFfile $GTF_FILE \
     > $OUTPUT_DIR/star_genome_generate.log 2>&1
