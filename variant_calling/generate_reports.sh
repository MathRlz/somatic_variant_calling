#!/bin/bash

# Generate summary reports and statistics

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"
REPORTS_DIR="${DATA_DIR}/reports"

mkdir -p "$REPORTS_DIR"

# Track processed samples for summary
PROCESSED_SAMPLES=()

for vcf in "$PRIORITY_DIR"/*_high_confidence.vcf.gz; do
    [ -f "$vcf" ] || continue

    sample=$(basename "$vcf" _high_confidence.vcf.gz)
    PROCESSED_SAMPLES+=("$sample")
    echo "Generating report for $sample..."

    # Variant statistics
    bcftools stats "$vcf" > "$REPORTS_DIR/${sample}_stats.txt"

    # Count variants by type
    bcftools view -H "$vcf" | wc -l > "$REPORTS_DIR/${sample}_total_variants.txt"

    # SNV count
    bcftools view -v snps -H "$vcf" | wc -l > "$REPORTS_DIR/${sample}_snvs.txt"

    # Indel count
    bcftools view -v indels -H "$vcf" | wc -l > "$REPORTS_DIR/${sample}_indels.txt"

    echo "Completed $sample"
done

# Create summary table from processed samples
echo -e "Sample\tTotal_Variants\tSNVs\tIndels" > "$REPORTS_DIR/summary.txt"
for sample in "${PROCESSED_SAMPLES[@]}"; do
    total=$(cat "$REPORTS_DIR/${sample}_total_variants.txt")
    snvs=$(cat "$REPORTS_DIR/${sample}_snvs.txt")
    indels=$(cat "$REPORTS_DIR/${sample}_indels.txt")
    echo -e "${sample}\t${total}\t${snvs}\t${indels}" >> "$REPORTS_DIR/summary.txt"
done

if [ ${#PROCESSED_SAMPLES[@]} -eq 0 ]; then
    echo "Warning: No prioritized VCF files found in $PRIORITY_DIR"
fi

echo "Reports generated in: $REPORTS_DIR/"
cat "$REPORTS_DIR/summary.txt"
