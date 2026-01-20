#!/bin/bash

# Download GRCh38 Reference Genome

set -e

echo "========================================="
echo "Downloading GRCh38 Reference Genome"
echo "========================================="

# Check if environment is set up
if [ -z "$DATA_DIR" ]; then
    echo "Error: DATA_DIR not set. Please run:"
    echo "  source ~/somatic_variant_calling/setup_env.sh"
    exit 1
fi

REF_DIR="$DATA_DIR/reference"
mkdir -p "$REF_DIR"

echo "Reference directory: $REF_DIR"
echo ""

# Check available space
AVAILABLE_SPACE=$(df -BG "$REF_DIR" | tail -1 | awk '{print $4}' | sed 's/G//')
echo "Available space: ${AVAILABLE_SPACE}GB"
if [ "$AVAILABLE_SPACE" -lt 20 ]; then
    echo "Warning: Less than 20GB available for reference genome"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

cd "$REF_DIR"

# Download GRCh38 reference genome from NCBI
echo "Downloading GRCh38 primary assembly..."

if [ ! -f "GRCh38_reference.fa" ]; then
    # Download from NCBI
    echo "Downloading reference genome with aria2c..."
    aria2c -c -x 16 -s 16 -k 1M -o GRCh38_reference.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    
    echo "Decompressing reference genome..."
    gunzip GRCh38_reference.fa.gz
    
    echo "Creating BWA index..."
    bwa index GRCh38_reference.fa
    
    echo "Creating SAMtools index..."
    samtools faidx GRCh38_reference.fa
    
    echo "Creating sequence dictionary for GATK..."
    java -jar $PROJECT_DIR/software/picard.jar CreateSequenceDictionary \
        R=GRCh38_reference.fa \
        O=GRCh38_reference.dict
    
    echo "Reference genome setup complete!"
else
    echo "Reference genome already exists."
fi

# Download dbSNP for variant annotation (optional but recommended)
echo ""
echo "Downloading dbSNP database for filtering..."
if [ ! -f "dbsnp_156.grch38.vcf.gz" ] || [ -f "dbsnp_156.grch38.vcf.gz.aria2" ]; then
    rm -f GCF_000001405.40.gz
    echo "Downloading dbSNP with aria2c..."
    aria2c -c -x 16 -s 16 -k 1M -o dbsnp_156.grch38.vcf.gz https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
    
    echo "Indexing dbSNP..."
    bcftools index -t dbsnp_156.grch38.vcf.gz
else
    echo "dbSNP already exists."
fi

# Download gnomAD for population filtering
echo ""
echo "Note: gnomAD database is very large (~100GB)."
echo "You may download it later if needed for more stringent filtering."
echo "URL: https://gnomad.broadinstitute.org/downloads"

echo ""
echo "========================================="
echo "Reference genome download complete!"
echo "========================================="
echo "Location: $REF_DIR"
echo ""
echo "Files created:"
ls -lh "$REF_DIR"/GRCh38_reference.fa*

