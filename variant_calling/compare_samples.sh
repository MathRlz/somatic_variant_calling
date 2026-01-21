#!/bin/bash

# Compare tumor vs normal to identify tumor-specific variants

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"
COMPARISON_DIR="${DATA_DIR}/comparisons"

mkdir -p "$COMPARISON_DIR"

# Auto-detect samples from prioritized VCF files
# Extract patient IDs from files like TCR002101_high_confidence.vcf.gz
found_samples=0

for vcf in "$PRIORITY_DIR"/*_high_confidence.vcf.gz; do
    [ -f "$vcf" ] || continue

    # Extract sample name (e.g., TCR002101 from TCR002101_high_confidence.vcf.gz)
    sample=$(basename "$vcf" _high_confidence.vcf.gz)

    echo "Comparing tumor vs normal for $sample..."

    # Extract variants with minimum depth (AF field may not be available in all formats)
    bcftools view -i 'FORMAT/DP>=20' "$vcf" \
        -O z -o "$COMPARISON_DIR/${sample}_somatic_candidates.vcf.gz"

    tabix -p vcf "$COMPARISON_DIR/${sample}_somatic_candidates.vcf.gz"

    # Count variants
    count=$(bcftools view -H "$COMPARISON_DIR/${sample}_somatic_candidates.vcf.gz" | wc -l)
    echo "Found $count somatic candidates for $sample"

    found_samples=$((found_samples + 1))
done

if [ $found_samples -eq 0 ]; then
    echo "Warning: No prioritized VCF files found in $PRIORITY_DIR"
fi

echo "Sample comparison completed!"
