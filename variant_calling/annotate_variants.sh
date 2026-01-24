#!/bin/bash
# Annotate variants with VEP (Variant Effect Predictor) for a single patient
# Usage: ./annotate_variants.sh PATIENT_ID
# Example: ./annotate_variants.sh TCR002101

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID"
    echo "Example: $0 TCR002101"
    exit 1
fi

PATIENT=$1

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
ANNOTATED_DIR="${DATA_DIR}/vcfs_annotated"
VEP_CACHE="$HOME/.vep"

mkdir -p "$ANNOTATED_DIR"

INPUT_VCF="$FILTERED_DIR/${PATIENT}_filtered.vcf.gz"
OUTPUT_VCF="$ANNOTATED_DIR/${PATIENT}_annotated.vcf.gz"

# Validate input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: $INPUT_VCF"
    exit 1
fi

# Check if VEP is installed
if ! command -v vep &> /dev/null; then
    echo "WARNING: VEP not installed. Skipping annotation."
    echo "To install VEP, run:"
    echo "  conda install -c bioconda ensembl-vep"
    echo "Or visit: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html"
    echo ""
    echo "Copying filtered VCF to annotated directory as placeholder..."
    cp "$INPUT_VCF" "$OUTPUT_VCF"
    cp "${INPUT_VCF}.tbi" "${OUTPUT_VCF}.tbi" 2>/dev/null || true
    echo "Skipped annotation (VEP not available)"
    exit 0
fi

echo "Annotating $PATIENT..."
echo "  Input: $INPUT_VCF"
echo "  Output: $OUTPUT_VCF"

vep \
    --input_file "$INPUT_VCF" \
    --output_file "$ANNOTATED_DIR/${PATIENT}_annotated.vcf" \
    --format vcf \
    --vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache \
    --dir_cache "$VEP_CACHE" \
    --offline \
    --everything \
    --fork 4

# Compress and index
bgzip "$ANNOTATED_DIR/${PATIENT}_annotated.vcf"
tabix -p vcf "$OUTPUT_VCF"

echo "Completed $PATIENT"
