#!/bin/bash

# Create IGV session file for visualization

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"
BAM_DIR="${DATA_DIR}/bam_recalibrated"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"
OUTPUT_FILE="${DATA_DIR}/igv_session.xml"

cat > "$OUTPUT_FILE" << EOF
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="hg38" version="8">
    <Resources>
        <Resource path="$REFERENCE"/>
EOF

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

# Auto-detect samples from recalibrated BAM files
# Look for tumor files (*-T_recalibrated.bam) and extract patient IDs
found_samples=0

for tumor_bam in "$BAM_DIR"/*-T_recalibrated.bam; do
    [ -f "$tumor_bam" ] || continue

    # Extract patient ID (e.g., TCR002101 from TCR002101-T_recalibrated.bam)
    tumor_name=$(basename "$tumor_bam" _recalibrated.bam)
    sample="${tumor_name%-T}"

    # Check if patient should be processed
    if ! should_process_patient "$sample"; then
        echo "Skipping $sample (not in selected patients)"
        continue
    fi

    normal_bam="$BAM_DIR/${sample}-N_recalibrated.bam"
    vcf_file="$PRIORITY_DIR/${sample}_high_confidence.vcf.gz"

    # Add tumor BAM
    echo "        <Resource path=\"$tumor_bam\"/>" >> "$OUTPUT_FILE"

    # Add normal BAM if it exists
    if [ -f "$normal_bam" ]; then
        echo "        <Resource path=\"$normal_bam\"/>" >> "$OUTPUT_FILE"
    fi

    # Add VCF if it exists
    if [ -f "$vcf_file" ]; then
        echo "        <Resource path=\"$vcf_file\"/>" >> "$OUTPUT_FILE"
    fi

    found_samples=$((found_samples + 1))
done

cat >> "$OUTPUT_FILE" << EOF
    </Resources>
</Session>
EOF

if [ $found_samples -eq 0 ]; then
    echo "Warning: No samples found in $BAM_DIR"
fi

echo "IGV session file created: $OUTPUT_FILE"
echo "Open this file in IGV to visualize your results"
