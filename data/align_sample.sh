#!/bin/bash
# Align a single paired-end sample to reference genome
# Usage: ./align_sample.sh SAMPLE_NAME
# Example: ./align_sample.sh TCR002101-T

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 TCR002101-T"
    exit 1
fi

SAMPLE=$1

# Use DATA_DIR if set, otherwise try to find reference relative to script location
if [ -n "$DATA_DIR" ]; then
    REF="$DATA_DIR/reference/GRCh38_reference.fa"
    FASTQ_DIR="$DATA_DIR/fastq"
    BAM_DIR="$DATA_DIR/bam"
else
    # Try relative path from fastq directory
    REF="../reference/GRCh38_reference.fa"
    FASTQ_DIR="."
    BAM_DIR="../bam"
fi

# Validate reference exists
if [ ! -f "$REF" ]; then
    echo "Error: Reference genome not found at $REF"
    echo "Please either:"
    echo "  1. Set DATA_DIR environment variable and run from any directory"
    echo "  2. Run this script from the fastq directory with reference in ../reference/"
    exit 1
fi

# Set number of processors (default to 8 if not set)
if [ -z "$NUM_PROCESSORS" ]; then
    NUM_PROCESSORS=8
fi

# Create output directory
mkdir -p "$BAM_DIR"

# Input files
R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

# Validate input files exist
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

# Verify FASTQ file integrity before alignment
echo "Verifying FASTQ file integrity..."
if ! gzip -t "$R1" 2>/dev/null; then
    echo "ERROR: $R1 is corrupted (gzip integrity check failed)"
    exit 1
fi
if ! gzip -t "$R2" 2>/dev/null; then
    echo "ERROR: $R2 is corrupted (gzip integrity check failed)"
    exit 1
fi

echo "Aligning $SAMPLE..."
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Reference: $REF"
echo "  Threads: $NUM_PROCESSORS"

# Align and sort
bwa mem -t "$NUM_PROCESSORS" "$REF" "$R1" "$R2" | samtools sort -o "$BAM_DIR/${SAMPLE}.bam"

# Index the BAM file
samtools index "$BAM_DIR/${SAMPLE}.bam"

echo "Completed $SAMPLE"
echo "Output: $BAM_DIR/${SAMPLE}.bam"
