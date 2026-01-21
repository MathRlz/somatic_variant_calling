#!/bin/bash

# Add read groups to marked BAM files

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

MARKED_DIR="${DATA_DIR}/bam_marked"
RG_DIR="${DATA_DIR}/bam_with_rg"

mkdir -p "$RG_DIR"

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

for bam in "$MARKED_DIR"/*_marked.bam; do
    sample=$(basename "$bam" _marked.bam)

    # Check if sample should be processed
    if ! should_process_sample "$sample"; then
        echo "Skipping $sample (not in selected patients)"
        continue
    fi

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
