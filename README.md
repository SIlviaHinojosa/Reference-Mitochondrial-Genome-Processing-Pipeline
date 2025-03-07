# Reference-Mitochondrial-Genome-Processing-Pipeline
Reference Mitochondrial Genome Processing Pipeline
Overview

This script automates the alignment of paired-end sequencing reads to a reference mitochondrial genome, processes the resulting alignments, and generates consensus sequences for each sample.

Requirements

Before running this script, ensure you have the following installed:

bwa (for read alignment)

samtools (for BAM file processing)

bcftools (for variant calling and consensus sequence generation)

Installation

Install the required tools using your package manager:

sudo apt install bwa samtools bcftools  # Ubuntu/Debian
brew install bwa samtools bcftools      # macOS

Alternatively, you can install them from source or use Conda if preferred.

Clone this repository:

git clone https://github.com/your-username/your-repo.git
cd your-repo

Usage

Run the script in the directory containing your paired-end FASTQ files:

bash reference_mito_curated.sh

Expected Input

Paired-end FASTQ files (*.R1.fq.gz and *.R2.fq.gz).

A reference mitochondrial genome (sequence2.fasta).

Output

The script creates the following directories:

bam_files/ – Sorted BAM files for each sample.

consensus_sequences/ – Consensus FASTA files.

logs/ – Log files for intermediate steps.

Customization

Modify the REFERENCE and THREADS variables in the script to specify your reference genome and number of threads to use.

Update the paths to samtools and bcftools in the script with your system's installation paths.
