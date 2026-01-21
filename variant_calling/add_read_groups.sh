#!/bin/bash
# Add read groups to a single marked BAM file
# Usage: ./add_read_groups.sh SAMPLE_NAME
# Example: ./add_read_groups.sh TCR002101-T

if [ $# -lt 1 ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 TCR002101-T"
    exit 1
fi

SAMPLE=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

MARKED_DIR="${DATA_DIR}/bam_marked"
RG_DIR="${DATA_DIR}/bam_with_rg"

mkdir -p "$RG_DIR"

INPUT_BAM="$MARKED_DIR/${SAMPLE}_marked.bam"
OUTPUT_BAM="$RG_DIR/${SAMPLE}_rg.bam"

# Validate input exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    exit 1
fi

echo "Adding read groups to $SAMPLE..."
echo "  Input: $INPUT_BAM"
echo "  Output: $OUTPUT_BAM"

# Extract sample info for read group
if [[ $SAMPLE == *"-T" ]]; then
    tissue="tumor"
else
    tissue="normal"
fi

gatk AddOrReplaceReadGroups \
    -I "$INPUT_BAM" \
    -O "$OUTPUT_BAM" \
    -RGID "${SAMPLE}" \
    -RGLB "lib1" \
    -RGPL "ILLUMINA" \
    -RGPU "unit1" \
    -RGSM "${SAMPLE}"

echo "Completed $SAMPLE"
