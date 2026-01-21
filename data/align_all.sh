#!/bin/bash
# Batch alignment for all paired-end samples in current directory

# Use DATA_DIR if set, otherwise try to find reference relative to script location
if [ -n "$DATA_DIR" ]; then
    REF="$DATA_DIR/reference/GRCh38_reference.fa"
else
    # Try relative path from fastq directory
    REF="../reference/GRCh38_reference.fa"
fi

# Validate reference exists
if [ ! -f "$REF" ]; then
    echo "Error: Reference genome not found at $REF"
    echo "Please either:"
    echo "  1. Set DATA_DIR environment variable and run from any directory"
    echo "  2. Run this script from the fastq directory with reference in ../reference/"
    exit 1
fi

# Set number of processors (default to 8 if not set)
if [ -z "$NUM_PROCESSORS" ]; then
    NUM_PROCESSORS=8
fi

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

for PREFIX in *_R1.fastq.gz; do
    SAMPLE=${PREFIX%_R1.fastq.gz}

    # Check if sample should be processed
    if ! should_process_sample "$SAMPLE"; then
        echo "Skipping $SAMPLE (not in selected patients)"
        continue
    fi

    if [ -f "${SAMPLE}_R2.fastq.gz" ]; then
        echo "Aligning $SAMPLE ..."
        bwa mem -t "$NUM_PROCESSORS" "$REF" "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz" | samtools sort -o "${SAMPLE}.bam"
        samtools index "${SAMPLE}.bam"
    else
        echo "Warning: Paired file for $SAMPLE not found, skipping."
    fi
done

echo "All alignments complete. BAM files are ready for downstream analysis."
