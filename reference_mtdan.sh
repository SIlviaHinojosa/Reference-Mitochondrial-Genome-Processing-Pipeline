#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Export paths for required tools
export PATH=/your/path/to/samtools:$PATH  # Replace with your path to samtools
export PATH=/your/path/to/bcftools:$PATH  # Replace with your path to bcftools

# Define default variables (can be overridden by user arguments)
REFERENCE="sequence2.fasta"
THREADS=8

# Check if required tools are installed
command -v bwa >/dev/null 2>&1 || { echo "bwa is not installed. Exiting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools is not installed. Exiting."; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "bcftools is not installed. Exiting."; exit 1; }

# Create output directories
mkdir -p bam_files consensus_sequences logs

# Index the reference genome if not already indexed
if [ ! -f "$REFERENCE.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$REFERENCE"
fi

# Loop through all paired-end read files
for R1 in *.R1.fq.gz; do
    # Extract sample name
    SAMPLE=$(basename "$R1" .R1.fq.gz)
    R2="${SAMPLE}.R2.fq.gz"

    echo "Processing sample: $SAMPLE"

    # Step 1: Align reads to the reference mitochondrial genome
    bwa mem -t "$THREADS" "$REFERENCE" "$R1" "$R2" > "logs/${SAMPLE}.sam"

    # Step 2: Convert SAM to BAM, sort, and index
    samtools view -S -b "logs/${SAMPLE}.sam" > "logs/${SAMPLE}.bam"
    samtools sort "logs/${SAMPLE}.bam" -o "bam_files/${SAMPLE}_sorted.bam"
    samtools index "bam_files/${SAMPLE}_sorted.bam"

    # Step 3: Generate consensus sequence
    bcftools mpileup -f "$REFERENCE" "bam_files/${SAMPLE}_sorted.bam" | \
        bcftools call -c - | \
        bcftools consensus -f "$REFERENCE" > "consensus_sequences/${SAMPLE}_consensus.fasta"

    echo "Completed processing for sample: $SAMPLE"

    # Cleanup intermediate files
    rm -f "logs/${SAMPLE}.sam" "logs/${SAMPLE}.bam"
done

echo "All samples processed successfully!"
