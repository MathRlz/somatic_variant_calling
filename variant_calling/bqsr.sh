#!/bin/bash
# Base Quality Score Recalibration for a single sample
# Usage: ./bqsr.sh SAMPLE_NAME
# Example: ./bqsr.sh TCR002101-T

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 TCR002101-T"
    exit 1
fi

SAMPLE=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

RG_DIR="${DATA_DIR}/bam_with_rg"
RECAL_DIR="${DATA_DIR}/bam_recalibrated"
RECAL_TABLES_DIR="${DATA_DIR}/recal_tables"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"

# Try multiple known sites files (dbSNP) - use whichever exists
KNOWN_SITES=""
for candidate in "${DATA_DIR}/reference/dbsnp_156.grch38.vcf.gz" \
                 "${DATA_DIR}/reference/All_20180418.vcf.gz" \
                 "${DATA_DIR}/reference/dbsnp.vcf.gz"; do
    if [ -f "$candidate" ]; then
        KNOWN_SITES="$candidate"
        break
    fi
done

# Validate required files exist
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome not found at $REFERENCE"
    exit 1
fi

if [ -z "$KNOWN_SITES" ]; then
    echo "Error: Known sites VCF not found. Expected one of:"
    echo "  - ${DATA_DIR}/reference/dbsnp_156.grch38.vcf.gz"
    echo "  - ${DATA_DIR}/reference/All_20180418.vcf.gz"
    echo "Run data/download_reference.sh to download the required files."
    exit 1
fi

mkdir -p "$RECAL_DIR" "$RECAL_TABLES_DIR"

INPUT_BAM="$RG_DIR/${SAMPLE}_rg.bam"
OUTPUT_BAM="$RECAL_DIR/${SAMPLE}_recalibrated.bam"

# Validate input exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    exit 1
fi

echo "Processing $SAMPLE..."
echo "  Input: $INPUT_BAM"
echo "  Output: $OUTPUT_BAM"
echo "  Reference: $REFERENCE"
echo "  Known sites: $KNOWN_SITES"

# Step 1: Build recalibration table
RECAL_TABLE="$RECAL_TABLES_DIR/${SAMPLE}_recal_data.table"
if ! gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" BaseRecalibrator \
    -I "$INPUT_BAM" \
    -R "$REFERENCE" \
    --known-sites "$KNOWN_SITES" \
    -O "$RECAL_TABLE"; then
    echo "Error: BaseRecalibrator failed for $SAMPLE"
    exit 1
fi

# Validate recalibration table was created
if [ ! -s "$RECAL_TABLE" ]; then
    echo "Error: Recalibration table is empty or missing: $RECAL_TABLE"
    exit 1
fi

# Step 2: Apply recalibration
if ! gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" ApplyBQSR \
    -I "$INPUT_BAM" \
    -R "$REFERENCE" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -O "$OUTPUT_BAM"; then
    echo "Error: ApplyBQSR failed for $SAMPLE"
    rm -f "$OUTPUT_BAM"  # Remove potentially corrupted output
    exit 1
fi

# Validate output BAM file
if [ ! -s "$OUTPUT_BAM" ]; then
    echo "Error: Output BAM is empty: $OUTPUT_BAM"
    exit 1
fi

# Check output BAM has reasonable size (at least 1MB for a valid BAM)
OUTPUT_SIZE=$(stat -c%s "$OUTPUT_BAM" 2>/dev/null || echo 0)
if [ "$OUTPUT_SIZE" -lt 1048576 ]; then
    echo "Error: Output BAM is suspiciously small (${OUTPUT_SIZE} bytes): $OUTPUT_BAM"
    echo "This usually indicates a failed recalibration."
    rm -f "$OUTPUT_BAM"
    exit 1
fi

echo "Completed $SAMPLE"
