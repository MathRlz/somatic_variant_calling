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
# IMPORTANT: Must use dbSNP with chr-prefixed contig names to match our reference
echo ""
echo "Downloading dbSNP database for filtering..."
if [ ! -f "dbsnp_156.grch38.vcf.gz" ] || [ -f "dbsnp_156.grch38.vcf.gz.aria2" ]; then
    # Download dbSNP from NCBI and remap contig names to match our UCSC-style reference
    echo "Downloading dbSNP from NCBI..."
    aria2c -c -x 16 -s 16 -k 1M -o dbsnp_ncbi.vcf.gz \
        "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"

    echo "Creating contig name mapping (NCBI RefSeq -> UCSC chr)..."
    # Create mapping file for the main chromosomes
    cat > contig_map.txt << 'CONTIGMAP'
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
NC_000004.12 chr4
NC_000005.10 chr5
NC_000006.12 chr6
NC_000007.14 chr7
NC_000008.11 chr8
NC_000009.12 chr9
NC_000010.11 chr10
NC_000011.10 chr11
NC_000012.12 chr12
NC_000013.11 chr13
NC_000014.9 chr14
NC_000015.10 chr15
NC_000016.10 chr16
NC_000017.11 chr17
NC_000018.10 chr18
NC_000019.10 chr19
NC_000020.11 chr20
NC_000021.9 chr21
NC_000022.11 chr22
NC_000023.11 chrX
NC_000024.10 chrY
NC_012920.1 chrM
CONTIGMAP

    echo "Remapping contig names to UCSC style (this may take a while)..."
    bcftools annotate --rename-chrs contig_map.txt dbsnp_ncbi.vcf.gz -Oz -o dbsnp_156.grch38.vcf.gz

    echo "Indexing remapped dbSNP..."
    bcftools index -t dbsnp_156.grch38.vcf.gz

    # Cleanup
    rm -f dbsnp_ncbi.vcf.gz contig_map.txt
    echo "dbSNP download and remapping complete!"
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

