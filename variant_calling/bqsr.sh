#!/bin/bash

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

RG_DIR="${DATA_DIR}/bam_with_rg"
RECAL_DIR="${DATA_DIR}/bam_recalibrated"
RECAL_TABLES_DIR="${DATA_DIR}/recal_tables"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"

# Try multiple known sites files (dbSNP) - use whichever exists
# The download_reference.sh script downloads dbsnp_156.grch38.vcf.gz
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

echo "Using known sites: $KNOWN_SITES"

# Set number of processors (for informational purposes only - BQSR has limited threading benefit)
if [ -z "$NUM_PROCESSORS" ]; then
    NUM_PROCESSORS=2
fi

mkdir -p "$RECAL_DIR" "$RECAL_TABLES_DIR"

for bam in "$RG_DIR"/*_rg.bam; do
    sample=$(basename "$bam" _rg.bam)
    echo "Processing $sample..."
    
    # Step 1: Build recalibration table
    gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" BaseRecalibrator \
        -I "$bam" \
        -R "$REFERENCE" \
        --known-sites "$KNOWN_SITES" \
        -O "$RECAL_TABLES_DIR/${sample}_recal_data.table"

    # Step 2: Apply recalibration
    gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" ApplyBQSR \
        -I "$bam" \
        -R "$REFERENCE" \
        --bqsr-recal-file "$RECAL_TABLES_DIR/${sample}_recal_data.table" \
        -O "$RECAL_DIR/${sample}_recalibrated.bam"
    
    echo "Completed $sample"
done

echo "BQSR completed for all samples!"
