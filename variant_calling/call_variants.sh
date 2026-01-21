#!/bin/bash

# Call somatic variants using GATK Mutect2

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

RECAL_DIR="${DATA_DIR}/bam_recalibrated"
VCF_DIR="${DATA_DIR}/vcfs"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"

# Set number of processors (default to 4 if not set - optimal for Mutect2)
if [ -z "$NUM_PROCESSORS" ]; then
    NUM_PROCESSORS=4
fi

mkdir -p "$VCF_DIR"

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

# Auto-detect tumor-normal pairs from recalibrated BAM files
# Tumor samples end with -T, Normal samples end with -N
declare -A PAIRS

for tumor_bam in "$RECAL_DIR"/*-T_recalibrated.bam; do
    [ -f "$tumor_bam" ] || continue

    # Extract tumor sample name (e.g., TCR002101-T from TCR002101-T_recalibrated.bam)
    tumor=$(basename "$tumor_bam" _recalibrated.bam)
    # Extract patient ID (e.g., TCR002101 from TCR002101-T)
    patient_id="${tumor%-T}"

    # Check if patient should be processed
    if ! should_process_patient "$patient_id"; then
        echo "Skipping $patient_id (not in selected patients)"
        continue
    fi

    # Construct normal sample name
    normal="${patient_id}-N"
    normal_bam="$RECAL_DIR/${normal}_recalibrated.bam"

    # Only add pair if both tumor and normal exist
    if [ -f "$normal_bam" ]; then
        PAIRS["$tumor"]="$normal"
        echo "Found tumor-normal pair: $tumor / $normal"
    else
        echo "Warning: No matching normal found for $tumor (expected $normal_bam)"
    fi
done

if [ ${#PAIRS[@]} -eq 0 ]; then
    echo "Error: No tumor-normal pairs found in $RECAL_DIR"
    exit 1
fi

for tumor in "${!PAIRS[@]}"; do
    normal="${PAIRS[$tumor]}"
    pair_name="${tumor%-T}"

    echo "Calling variants for pair: $pair_name (using $NUM_PROCESSORS threads)"

    gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" Mutect2 \
        -R "$REFERENCE" \
        -I "$RECAL_DIR/${tumor}_recalibrated.bam" \
        -I "$RECAL_DIR/${normal}_recalibrated.bam" \
        -normal "${normal}" \
        -O "$VCF_DIR/${pair_name}_raw.vcf.gz" \
        --f1r2-tar-gz "$VCF_DIR/${pair_name}_f1r2.tar.gz" \
        --native-pair-hmm-threads "$NUM_PROCESSORS"

    echo "Completed $pair_name"
done

echo "Variant calling completed!"
