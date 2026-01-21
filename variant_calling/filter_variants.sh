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

for vcf in "$VCF_DIR"/*_raw.vcf.gz; do
    sample=$(basename "$vcf" _raw.vcf.gz)

    # Check if patient should be processed
    if ! should_process_patient "$sample"; then
        echo "Skipping $sample (not in selected patients)"
        continue
    fi

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
