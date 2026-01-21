#!/bin/bash

# Add read groups to marked BAM files

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

MARKED_DIR="${DATA_DIR}/bam_marked"
RG_DIR="${DATA_DIR}/bam_with_rg"

mkdir -p "$RG_DIR"

for bam in "$MARKED_DIR"/*_marked.bam; do
    sample=$(basename "$bam" _marked.bam)
    echo "Adding read groups to $sample..."
    
    # Extract sample info for read group
    if [[ $sample == *"-T" ]]; then
        tissue="tumor"
    else
        tissue="normal"
    fi
    
    gatk AddOrReplaceReadGroups \
        -I "$bam" \
        -O "$RG_DIR/${sample}_rg.bam" \
        -RGID "${sample}" \
        -RGLB "lib1" \
        -RGPL "ILLUMINA" \
        -RGPU "unit1" \
        -RGSM "${sample}"
    
    echo "Completed $sample"
done

echo "Read groups added to all samples!"
