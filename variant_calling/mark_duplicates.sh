#!/bin/bash

# Mark duplicates in BAM files using Picard
# This removes PCR and optical duplicates

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

BAM_DIR="${DATA_DIR}/bam"
MARKED_DIR="${DATA_DIR}/bam_marked"
METRICS_DIR="${DATA_DIR}/metrics"

mkdir -p "$MARKED_DIR" "$METRICS_DIR"

# Process each BAM file
for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "Processing $sample..."
    
    gatk MarkDuplicates \
        -I "$bam" \
        -O "$MARKED_DIR/${sample}_marked.bam" \
        -M "$METRICS_DIR/${sample}_dup_metrics.txt" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT
    
    echo "Completed $sample"
done

echo "All samples processed!"
