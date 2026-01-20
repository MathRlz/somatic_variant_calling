#!/bin/bash

# Compare tumor vs normal to identify tumor-specific variants

PRIORITY_DIR="vcfs_prioritized"
COMPARISON_DIR="comparisons"

mkdir -p "$COMPARISON_DIR"

# Define pairs
declare -A PAIRS=(
    ["TCR002101"]="TCR002101"
    ["TCR002182"]="TCR002182"
    ["TCR002361"]="TCR002361"
)

for sample in "${!PAIRS[@]}"; do
    echo "Comparing tumor vs normal for $sample..."
    
    tumor="$PRIORITY_DIR/${sample}_high_confidence.vcf.gz"
    
    if [ -f "$tumor" ]; then
        # Extract variants with minimum depth (AF field may not be available in all formats)
        bcftools view -i 'FORMAT/DP>=20' "$tumor" \
            -O z -o "$COMPARISON_DIR/${sample}_somatic_candidates.vcf.gz"
        
        tabix -p vcf "$COMPARISON_DIR/${sample}_somatic_candidates.vcf.gz"
        
        # Count variants
        count=$(bcftools view -H "$COMPARISON_DIR/${sample}_somatic_candidates.vcf.gz" | wc -l)
        echo "Found $count somatic candidates for $sample"
    fi
done

echo "Sample comparison completed!"
