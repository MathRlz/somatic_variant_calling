#!/bin/bash

# Generate summary reports and statistics

PRIORITY_DIR="vcfs_prioritized"
REPORTS_DIR="reports"

mkdir -p "$REPORTS_DIR"

for vcf in "$PRIORITY_DIR"/*_high_confidence.vcf.gz; do
    sample=$(basename "$vcf" _high_confidence.vcf.gz)
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

# Create summary table
echo -e "Sample\tTotal_Variants\tSNVs\tIndels" > "$REPORTS_DIR/summary.txt"
for sample in TCR002101 TCR002182 TCR002361; do
    total=$(cat "$REPORTS_DIR/${sample}_total_variants.txt")
    snvs=$(cat "$REPORTS_DIR/${sample}_snvs.txt")
    indels=$(cat "$REPORTS_DIR/${sample}_indels.txt")
    echo -e "${sample}\t${total}\t${snvs}\t${indels}" >> "$REPORTS_DIR/summary.txt"
done

echo "Reports generated in: $REPORTS_DIR/"
cat "$REPORTS_DIR/summary.txt"
