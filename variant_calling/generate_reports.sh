#!/bin/bash
# Generate summary reports and statistics for a single patient
# Usage: ./generate_reports.sh PATIENT_ID
# Example: ./generate_reports.sh TCR002101

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID"
    echo "Example: $0 TCR002101"
    exit 1
fi

PATIENT=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"
REPORTS_DIR="${DATA_DIR}/reports"

mkdir -p "$REPORTS_DIR"

INPUT_VCF="$PRIORITY_DIR/${PATIENT}_high_confidence.vcf.gz"

# Validate input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

echo "Generating report for $PATIENT..."
echo "  Input: $INPUT_VCF"

# Variant statistics
bcftools stats "$INPUT_VCF" > "$REPORTS_DIR/${PATIENT}_stats.txt"

# Count variants by type
bcftools view -H "$INPUT_VCF" | wc -l > "$REPORTS_DIR/${PATIENT}_total_variants.txt"

# SNV count
bcftools view -v snps -H "$INPUT_VCF" | wc -l > "$REPORTS_DIR/${PATIENT}_snvs.txt"

# Indel count
bcftools view -v indels -H "$INPUT_VCF" | wc -l > "$REPORTS_DIR/${PATIENT}_indels.txt"

# Print summary
total=$(cat "$REPORTS_DIR/${PATIENT}_total_variants.txt")
snvs=$(cat "$REPORTS_DIR/${PATIENT}_snvs.txt")
indels=$(cat "$REPORTS_DIR/${PATIENT}_indels.txt")

echo ""
echo "Summary for $PATIENT:"
echo "  Total variants: $total"
echo "  SNVs: $snvs"
echo "  Indels: $indels"

# Append to summary file (create header if needed)
if [ ! -f "$REPORTS_DIR/summary.txt" ]; then
    echo -e "Sample\tTotal_Variants\tSNVs\tIndels" > "$REPORTS_DIR/summary.txt"
fi

# Check if patient already in summary, if so update it, otherwise append
if grep -q "^${PATIENT}\t" "$REPORTS_DIR/summary.txt"; then
    # Update existing line
    sed -i "s/^${PATIENT}\t.*/${PATIENT}\t${total}\t${snvs}\t${indels}/" "$REPORTS_DIR/summary.txt"
else
    echo -e "${PATIENT}\t${total}\t${snvs}\t${indels}" >> "$REPORTS_DIR/summary.txt"
fi

echo "Completed $PATIENT"
echo "Reports saved to: $REPORTS_DIR/"
