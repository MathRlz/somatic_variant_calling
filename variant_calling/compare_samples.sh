#!/bin/bash
# Compare tumor vs normal to identify tumor-specific variants for a single patient
# Usage: ./compare_samples.sh PATIENT_ID
# Example: ./compare_samples.sh TCR002101

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID"
    echo "Example: $0 TCR002101"
    exit 1
fi

PATIENT=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"
COMPARISON_DIR="${DATA_DIR}/comparisons"

mkdir -p "$COMPARISON_DIR"

INPUT_VCF="$PRIORITY_DIR/${PATIENT}_high_confidence.vcf.gz"
OUTPUT_VCF="$COMPARISON_DIR/${PATIENT}_somatic_candidates.vcf.gz"

# Validate input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

echo "Comparing tumor vs normal for $PATIENT..."
echo "  Input: $INPUT_VCF"
echo "  Output: $OUTPUT_VCF"

# Extract variants with minimum depth
bcftools view -i 'FORMAT/DP>=20' "$INPUT_VCF" \
    -O z -o "$OUTPUT_VCF"

tabix -p vcf "$OUTPUT_VCF"

# Count variants
count=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
echo "Found $count somatic candidates for $PATIENT"

echo "Completed $PATIENT"
