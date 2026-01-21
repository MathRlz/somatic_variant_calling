#!/bin/bash
# Filter variants using GATK FilterMutectCalls for a single patient
# Usage: ./filter_variants.sh PATIENT_ID
# Example: ./filter_variants.sh TCR002101

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID"
    echo "Example: $0 TCR002101"
    exit 1
fi

PATIENT=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

VCF_DIR="${DATA_DIR}/vcfs"
FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"

mkdir -p "$FILTERED_DIR"

# Validate reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome not found at $REFERENCE"
    exit 1
fi

INPUT_VCF="$VCF_DIR/${PATIENT}_raw.vcf.gz"
OUTPUT_VCF="$FILTERED_DIR/${PATIENT}_filtered.vcf.gz"

# Validate input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

echo "Filtering $PATIENT..."
echo "  Input: $INPUT_VCF"
echo "  Output: $OUTPUT_VCF"

# Learn orientation bias
gatk LearnReadOrientationModel \
    -I "$VCF_DIR/${PATIENT}_f1r2.tar.gz" \
    -O "$VCF_DIR/${PATIENT}_read-orientation-model.tar.gz"

# Filter variants
gatk FilterMutectCalls \
    -R "$REFERENCE" \
    -V "$INPUT_VCF" \
    --ob-priors "$VCF_DIR/${PATIENT}_read-orientation-model.tar.gz" \
    -O "$OUTPUT_VCF"

echo "Completed $PATIENT"
