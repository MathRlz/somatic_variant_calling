#!/bin/bash
# Batch alignment for all paired-end samples in current directory

# Use reference from parent directory structure
REF="../reference/GRCh38_reference.fa"

# If DATA_DIR is set, use that
if [ -n "$DATA_DIR" ]; then
    REF="$DATA_DIR/reference/GRCh38_reference.fa"
fi

for PREFIX in *_R1.fastq.gz; do
    SAMPLE=${PREFIX%_R1.fastq.gz}
    if [ -f "${SAMPLE}_R2.fastq.gz" ]; then
        echo "Aligning $SAMPLE ..."
        bwa mem -t 8 "$REF" "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz" | samtools sort -o "${SAMPLE}.bam"
        samtools index "${SAMPLE}.bam"
    else
        echo "Warning: Paired file for $SAMPLE not found, skipping."
    fi
done

echo "All alignments complete. BAM files are ready for downstream analysis."
