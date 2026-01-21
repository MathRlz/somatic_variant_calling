#!/bin/bash
# Call somatic variants using GATK Mutect2 for a single patient
# Usage: ./call_variants.sh PATIENT_ID
# Example: ./call_variants.sh TCR002101

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID"
    echo "Example: $0 TCR002101"
    exit 1
fi

PATIENT=$1

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

# Define tumor and normal sample names
TUMOR="${PATIENT}-T"
NORMAL="${PATIENT}-N"
TUMOR_BAM="$RECAL_DIR/${TUMOR}_recalibrated.bam"
NORMAL_BAM="$RECAL_DIR/${NORMAL}_recalibrated.bam"

# Validate input files exist
if [ ! -f "$TUMOR_BAM" ]; then
    echo "Error: Tumor BAM not found: $TUMOR_BAM"
    exit 1
fi

if [ ! -f "$NORMAL_BAM" ]; then
    echo "Error: Normal BAM not found: $NORMAL_BAM"
    exit 1
fi

echo "Calling variants for patient: $PATIENT"
echo "  Tumor: $TUMOR_BAM"
echo "  Normal: $NORMAL_BAM"
echo "  Threads: $NUM_PROCESSORS"

gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" Mutect2 \
    -R "$REFERENCE" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    -normal "${NORMAL}" \
    -O "$VCF_DIR/${PATIENT}_raw.vcf.gz" \
    --f1r2-tar-gz "$VCF_DIR/${PATIENT}_f1r2.tar.gz" \
    --native-pair-hmm-threads "$NUM_PROCESSORS"

echo "Completed $PATIENT"
echo "Output: $VCF_DIR/${PATIENT}_raw.vcf.gz"
