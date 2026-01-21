#!/bin/bash

# Annotate variants with VEP (Variant Effect Predictor)

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
ANNOTATED_DIR="${DATA_DIR}/vcfs_annotated"
VEP_CACHE="$HOME/.vep"

mkdir -p "$ANNOTATED_DIR"

# Function to check if patient should be processed based on selected patients
should_process_patient() {
    local patient_id=$1

    # If SELECTED_PATIENTS_FILE is set and exists, check if patient is selected
    if [ -n "$SELECTED_PATIENTS_FILE" ] && [ -f "$SELECTED_PATIENTS_FILE" ]; then
        if grep -q "^${patient_id}$" "$SELECTED_PATIENTS_FILE"; then
            return 0
        else
            return 1
        fi
    fi
    # If no selection file, process all patients
    return 0
}

# Check if VEP is installed
if ! command -v vep &> /dev/null; then
    echo "WARNING: VEP not installed. Skipping annotation."
    echo "To install VEP, run:"
    echo "  conda install -c bioconda ensembl-vep"
    echo "Or visit: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html"
    echo ""
    echo "Copying filtered VCFs to annotated directory as placeholders..."
    for vcf in "$FILTERED_DIR"/*_filtered.vcf.gz; do
        sample=$(basename "$vcf" _filtered.vcf.gz)

        # Check if patient should be processed
        if ! should_process_patient "$sample"; then
            echo "Skipping $sample (not in selected patients)"
            continue
        fi

        cp "$vcf" "$ANNOTATED_DIR/${sample}_annotated.vcf.gz"
        cp "${vcf}.tbi" "$ANNOTATED_DIR/${sample}_annotated.vcf.gz.tbi" 2>/dev/null || true
    done
    echo "Skipped annotation (VEP not available)"
    exit 0
fi

for vcf in "$FILTERED_DIR"/*_filtered.vcf.gz; do
    sample=$(basename "$vcf" _filtered.vcf.gz)

    # Check if patient should be processed
    if ! should_process_patient "$sample"; then
        echo "Skipping $sample (not in selected patients)"
        continue
    fi

    echo "Annotating $sample..."

    vep \
        --input_file "$vcf" \
        --output_file "$ANNOTATED_DIR/${sample}_annotated.vcf" \
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
    bgzip "$ANNOTATED_DIR/${sample}_annotated.vcf"
    tabix -p vcf "$ANNOTATED_DIR/${sample}_annotated.vcf.gz"

    echo "Completed $sample"
done

echo "Variant annotation completed!"
