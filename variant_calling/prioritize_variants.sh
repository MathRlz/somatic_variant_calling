#!/bin/bash
# Filter and prioritize somatic variants for a single patient
# Usage: ./prioritize_variants.sh PATIENT_ID
# Example: ./prioritize_variants.sh TCR002101

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID"
    echo "Example: $0 TCR002101"
    exit 1
fi

PATIENT=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"

mkdir -p "$PRIORITY_DIR"

INPUT_VCF="$FILTERED_DIR/${PATIENT}_filtered.vcf.gz"
OUTPUT_VCF="$PRIORITY_DIR/${PATIENT}_high_confidence.vcf.gz"

# Validate input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

echo "Prioritizing $PATIENT..."
echo "  Input: $INPUT_VCF"
echo "  Output: $OUTPUT_VCF"

# Extract PASS variants only with minimum depth
bcftools view -f PASS "$INPUT_VCF" \
    | bcftools view -i 'FORMAT/DP>=20' \
    -o "$OUTPUT_VCF" -O z

# Index
tabix -p vcf "$OUTPUT_VCF"

echo "Completed $PATIENT"
