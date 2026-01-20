#!/bin/bash

# Call somatic variants using GATK Mutect2

RECAL_DIR="bam_recalibrated"
VCF_DIR="vcfs"
REFERENCE="reference/GRCh38_reference.fa"

# Set number of processors (default to 4 if not set - optimal for Mutect2)
if [ -z "$NUM_PROCESSORS" ]; then
    NUM_PROCESSORS=4
fi

mkdir -p "$VCF_DIR"

# Define tumor-normal pairs
declare -A PAIRS=(
    ["TCR002101-T"]="TCR002101-N"
    ["TCR002182-T"]="TCR002182-N"
    ["TCR002361-T"]="TCR002361-N"
)

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
