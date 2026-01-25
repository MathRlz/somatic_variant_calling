#!/bin/bash
# Filter and prioritize somatic variants for a single patient
# Usage: ./prioritize_variants.sh PATIENT_ID [OUTPUT_SUFFIX]
# Example: ./prioritize_variants.sh TCR002101
# Example: MIN_DP=10 ./prioritize_variants.sh TCR002101 lenient
#
# Environment variables:
#   MIN_DP  - Minimum read depth (default: 20)
#   MIN_VAF - Minimum variant allele frequency (default: 0.05)

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID [OUTPUT_SUFFIX]"
    echo "Example: $0 TCR002101"
    echo "Example: MIN_DP=10 $0 TCR002101 lenient"
    exit 1
fi

PATIENT=$1
OUTPUT_SUFFIX="${2:-high_confidence}"

# Configurable thresholds
MIN_DP="${MIN_DP:-20}"
MIN_VAF="${MIN_VAF:-0.05}"

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"

mkdir -p "$PRIORITY_DIR"

INPUT_VCF="$FILTERED_DIR/${PATIENT}_filtered.vcf.gz"
OUTPUT_VCF="$PRIORITY_DIR/${PATIENT}_${OUTPUT_SUFFIX}.vcf.gz"

# Validate input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

echo "Prioritizing $PATIENT..."
echo "  Input: $INPUT_VCF"
echo "  Output: $OUTPUT_VCF"
echo "  MIN_DP: $MIN_DP, MIN_VAF: $MIN_VAF"

# Filter variants using configurable thresholds:
# - DP >= MIN_DP (default 20)
# - VAF >= MIN_VAF (default 5%, removes noise while keeping subclonal drivers)
# - PASS or soft filters only (rescues potential drivers filtered by
#   weak_evidence, strand_bias, or clustered_events)
# Hard artifact filters (normal_artifact, contamination, orientation) are excluded
bcftools view -i "
  (FILTER=\"PASS\" ||
   FILTER~\"weak_evidence\" ||
   FILTER~\"strand_bias\" ||
   FILTER~\"clustered_events\") &&
  FORMAT/DP >= $MIN_DP &&
  FORMAT/AF >= $MIN_VAF
" "$INPUT_VCF" -o "$OUTPUT_VCF" -O z

# Count variants
count=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
echo "Found $count variants for $PATIENT"

# Index
tabix -p vcf "$OUTPUT_VCF"

echo "Completed $PATIENT"
