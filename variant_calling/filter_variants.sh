#!/bin/bash

# Filter variants using GATK FilterMutectCalls

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

VCF_DIR="${DATA_DIR}/vcfs"
FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"

mkdir -p "$FILTERED_DIR"

# Validate reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome not found at $REFERENCE"
    exit 1
fi

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
