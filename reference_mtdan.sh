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
    # Loop through all VCF files (without .gz) in logs directory
for VCF in logs/*_variants.vcf; do
    SAMPLE=$(basename "$VCF" "_variants.vcf")
    echo "Processing VCF for sample: $SAMPLE"

    # Step 4: Compress and index VCF
    if [[ -f "logs/${SAMPLE}_variants.vcf" ]]; then
        echo "Compressing logs/${SAMPLE}_variants.vcf ..."
        bgzip -f "logs/${SAMPLE}_variants.vcf"

        echo "Indexing logs/${SAMPLE}_variants.vcf.gz ..."
        tabix -p vcf "logs/${SAMPLE}_variants.vcf.gz"
    else
        echo "Warning: logs/${SAMPLE}_variants.vcf not found, skipping compression and indexing."
        continue
    fi

  # Step 5: Generate consensus sequence
    FASTA_OUT="/your/path/consensus_sequences/${SAMPLE}_mtDNA.fasta"

    if [[ -f "logs/${SAMPLE}_variants.vcf.gz" ]]; then
        echo "Generating consensus sequence for ${SAMPLE}..."
        bcftools consensus -f "$REFERENCE" -o "$FASTA_OUT" "logs/${SAMPLE}_variants.vcf.gz"

        # Fix FASTA header (replace reference header with sample name)
        sed -i "1s/.*/>${SAMPLE}_mtDNA/" "$FASTA_OUT"
    else
        echo "Error: logs/${SAMPLE}_variants.vcf.gz not found, skipping consensus generation."
    fi

done

echo "All samples processed successfully!"
