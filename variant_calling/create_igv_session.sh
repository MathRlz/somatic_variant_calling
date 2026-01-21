#!/bin/bash
# Create IGV session file for visualization
# Usage: ./create_igv_session.sh PATIENT_ID [PATIENT_ID2 ...]
# Example: ./create_igv_session.sh TCR002101
# Example: ./create_igv_session.sh TCR002101 TCR002102 TCR002103

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID [PATIENT_ID2 ...]"
    echo "Example: $0 TCR002101"
    echo "Example: $0 TCR002101 TCR002102 TCR002103"
    exit 1
fi

PATIENTS=("$@")

# Use DATA_DIR if set, otherwise assume we're running from the data directory
DATA_DIR="${DATA_DIR:-.}"

PRIORITY_DIR="${DATA_DIR}/vcfs_prioritized"
BAM_DIR="${DATA_DIR}/bam_recalibrated"
REFERENCE="${DATA_DIR}/reference/GRCh38_reference.fa"
OUTPUT_FILE="${DATA_DIR}/igv_session.xml"

echo "Creating IGV session for patients: ${PATIENTS[*]}"
echo "  Output: $OUTPUT_FILE"

# Start XML file
cat > "$OUTPUT_FILE" << EOF
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="hg38" version="8">
    <Resources>
        <Resource path="$REFERENCE"/>
EOF

# Track processed samples
found_samples=0

for PATIENT in "${PATIENTS[@]}"; do
    tumor_bam="$BAM_DIR/${PATIENT}-T_recalibrated.bam"
    normal_bam="$BAM_DIR/${PATIENT}-N_recalibrated.bam"
    vcf_file="$PRIORITY_DIR/${PATIENT}_high_confidence.vcf.gz"

    # Add tumor BAM if exists
    if [ -f "$tumor_bam" ]; then
        echo "        <Resource path=\"$tumor_bam\"/>" >> "$OUTPUT_FILE"
        found_samples=$((found_samples + 1))
    fi

    # Add normal BAM if exists
    if [ -f "$normal_bam" ]; then
        echo "        <Resource path=\"$normal_bam\"/>" >> "$OUTPUT_FILE"
    fi

    # Add VCF if exists
    if [ -f "$vcf_file" ]; then
        echo "        <Resource path=\"$vcf_file\"/>" >> "$OUTPUT_FILE"
    fi
done

# Close XML file
cat >> "$OUTPUT_FILE" << EOF
    </Resources>
</Session>
EOF

if [ $found_samples -eq 0 ]; then
    echo "Warning: No samples found in $BAM_DIR for specified patients"
fi

echo "IGV session file created: $OUTPUT_FILE"
echo "Open this file in IGV to visualize your results"
