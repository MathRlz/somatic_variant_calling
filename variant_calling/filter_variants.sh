#!/bin/bash

# Filter variants using GATK FilterMutectCalls

VCF_DIR="vcfs"
FILTERED_DIR="vcfs_filtered"
REFERENCE="reference/GRCh38_reference.fa"

mkdir -p "$FILTERED_DIR"

for vcf in "$VCF_DIR"/*_raw.vcf.gz; do
    sample=$(basename "$vcf" _raw.vcf.gz)
    echo "Filtering $sample..."
    
    # Learn orientation bias
    gatk LearnReadOrientationModel \
        -I "$VCF_DIR/${sample}_f1r2.tar.gz" \
        -O "$VCF_DIR/${sample}_read-orientation-model.tar.gz"
    
    # Filter variants
    gatk FilterMutectCalls \
        -R "$REFERENCE" \
        -V "$vcf" \
        --ob-priors "$VCF_DIR/${sample}_read-orientation-model.tar.gz" \
        -O "$FILTERED_DIR/${sample}_filtered.vcf.gz"
    
    echo "Completed $sample"
done

echo "Variant filtering completed!"
