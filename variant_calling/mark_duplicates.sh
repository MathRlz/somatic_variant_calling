#!/bin/bash

# Mark duplicates in BAM files using Picard
# This removes PCR and optical duplicates

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

BAM_DIR="${DATA_DIR}/bam"
MARKED_DIR="${DATA_DIR}/bam_marked"
METRICS_DIR="${DATA_DIR}/metrics"

mkdir -p "$MARKED_DIR" "$METRICS_DIR"

# Function to check if sample should be processed based on selected patients
should_process_sample() {
    local sample=$1
    # Extract patient ID from sample name (remove -T or -N suffix)
    local patient_id="${sample%-[TN]}"

    # If SELECTED_PATIENTS_FILE is set and exists, check if patient is selected
    if [ -n "$SELECTED_PATIENTS_FILE" ] && [ -f "$SELECTED_PATIENTS_FILE" ]; then
        if grep -q "^${patient_id}$" "$SELECTED_PATIENTS_FILE"; then
            return 0
        else
            return 1
        fi
    fi
    # If no selection file, process all samples
    return 0
}

# Process each BAM file
for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)

    # Check if sample should be processed
    if ! should_process_sample "$sample"; then
        echo "Skipping $sample (not in selected patients)"
        continue
    fi

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
