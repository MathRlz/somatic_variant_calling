#!/bin/bash

# Filter and prioritize somatic variants

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"

mkdir -p "$PRIORITY_DIR"

for vcf in "$FILTERED_DIR"/*_filtered.vcf.gz; do
    sample=$(basename "$vcf" _filtered.vcf.gz)
    echo "Prioritizing $sample..."
    
    # Extract PASS variants only with minimum depth
    bcftools view -f PASS "$vcf" \
        | bcftools view -i 'FORMAT/DP>=20' \
        -o "$PRIORITY_DIR/${sample}_high_confidence.vcf.gz" -O z
    
    # Index
    tabix -p vcf "$PRIORITY_DIR/${sample}_high_confidence.vcf.gz"
    
    echo "Completed $sample"
done

echo "Variant prioritization completed!"
