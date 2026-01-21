#!/bin/bash

# Filter and prioritize somatic variants

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

FILTERED_DIR="${DATA_DIR}/vcfs_filtered"
PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"

mkdir -p "$PRIORITY_DIR"

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

for vcf in "$FILTERED_DIR"/*_filtered.vcf.gz; do
    sample=$(basename "$vcf" _filtered.vcf.gz)

    # Check if patient should be processed
    if ! should_process_patient "$sample"; then
        echo "Skipping $sample (not in selected patients)"
        continue
    fi

    echo "Prioritizing $sample..."

    # Extract PASS variants only with minimum depth
    bcftools view -f PASS "$vcf" \
        | bcftools view -i 'FORMAT/DP>=20' \
        -o "$PRIORITY_DIR/${sample}_high_confidence.vcf.gz" -O z

    # Index
    tabix -p vcf "$PRIORITY_DIR/${sample}_high_confidence.vcf.gz"

    echo "Completed $sample"
done

echo "Variant prioritization completed!"
