#!/bin/bash
# Mark duplicates in a single BAM file using GATK
# Usage: ./mark_duplicates.sh SAMPLE_NAME
# Example: ./mark_duplicates.sh TCR002101-T

if [ $# -lt 1 ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 TCR002101-T"
    exit 1
fi

SAMPLE=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

BAM_DIR="${DATA_DIR}/bam"
MARKED_DIR="${DATA_DIR}/bam_marked"
METRICS_DIR="${DATA_DIR}/metrics"

mkdir -p "$MARKED_DIR" "$METRICS_DIR"

INPUT_BAM="$BAM_DIR/${SAMPLE}.bam"
OUTPUT_BAM="$MARKED_DIR/${SAMPLE}_marked.bam"

# Validate input exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    exit 1
fi

echo "Processing $SAMPLE..."
echo "  Input: $INPUT_BAM"
echo "  Output: $OUTPUT_BAM"

gatk MarkDuplicates \
    -I "$INPUT_BAM" \
    -O "$OUTPUT_BAM" \
    -M "$METRICS_DIR/${SAMPLE}_dup_metrics.txt" \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT

echo "Completed $SAMPLE"
