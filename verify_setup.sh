#!/bin/bash

# Verification script to check if environment is properly set up

echo "========================================="
echo "Somatic Variant Calling - Setup Verification"
echo "========================================="
echo ""

# Check if running from correct directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Set expected directories
PROJECT_DIR="${PROJECT_DIR:-$HOME/somatic_variant_calling}"
DATA_DIR="${DATA_DIR:-$HOME/storage/variant_calling_data}"

echo "Checking installation..."
echo "Project directory: $PROJECT_DIR"
echo "Data directory: $DATA_DIR"
echo ""

# Function to check if command exists
check_command() {
    if command -v "$1" &> /dev/null; then
        echo "✓ $1 found"
        return 0
    else
        echo "✗ $1 NOT found"
        return 1
    fi
}

# Function to check if directory exists
check_directory() {
    if [ -d "$1" ]; then
        echo "✓ Directory exists: $1"
        return 0
    else
        echo "✗ Directory missing: $1"
        return 1
    fi
}

# Function to check if file exists
check_file() {
    if [ -f "$1" ]; then
        echo "✓ File exists: $1"
        return 0
    else
        echo "✗ File missing: $1"
        return 1
    fi
}

ERRORS=0

echo "=== Checking Project Structure ==="
check_directory "$PROJECT_DIR" || ((ERRORS++))
check_directory "$DATA_DIR" || ((ERRORS++))
check_file "$PROJECT_DIR/setup_env.sh" || ((ERRORS++))
echo ""

echo "=== Checking Repository Scripts ==="
check_file "$SCRIPT_DIR/setup_start.sh" || ((ERRORS++))
check_file "$SCRIPT_DIR/run_complete_pipeline.sh" || ((ERRORS++))
check_directory "$SCRIPT_DIR/data" || ((ERRORS++))
check_directory "$SCRIPT_DIR/variant_calling" || ((ERRORS++))
echo ""

# Source environment if it exists
if [ -f "$PROJECT_DIR/setup_env.sh" ]; then
    echo "=== Sourcing Environment ==="
    source "$PROJECT_DIR/setup_env.sh"
    echo ""
fi

echo "=== Checking Bioinformatics Tools ==="
check_command "samtools" || ((ERRORS++))
check_command "bcftools" || ((ERRORS++))
check_command "bwa" || ((ERRORS++))
check_command "fastqc" || ((ERRORS++))
check_command "gatk" || ((ERRORS++))
check_command "java" || ((ERRORS++))
check_command "python3" || ((ERRORS++))
echo ""

echo "=== Checking Tool Versions ==="
if command -v samtools &> /dev/null; then
    echo "SAMtools version: $(samtools --version | head -1)"
fi
if command -v bcftools &> /dev/null; then
    echo "BCFtools version: $(bcftools --version | head -1)"
fi
if command -v bwa &> /dev/null; then
    echo "BWA version: $(bwa 2>&1 | grep Version || echo 'unknown')"
fi
if command -v gatk &> /dev/null; then
    echo "GATK version: $(gatk --version 2>&1 | grep 'GATK' | head -1 || echo 'unknown')"
fi
if command -v java &> /dev/null; then
    echo "Java version: $(java -version 2>&1 | head -1)"
fi
echo ""

echo "=== Checking Reference Data ==="
check_file "$DATA_DIR/reference/GRCh38_reference.fa" || echo "  (Run data/download_reference.sh to download)"
# Check for known sites VCF (dbSNP) - accept either filename
if [ -f "$DATA_DIR/reference/dbsnp_156.grch38.vcf.gz" ]; then
    echo "✓ File exists: $DATA_DIR/reference/dbsnp_156.grch38.vcf.gz"
elif [ -f "$DATA_DIR/reference/All_20180418.vcf.gz" ]; then
    echo "✓ File exists: $DATA_DIR/reference/All_20180418.vcf.gz"
else
    echo "✗ File missing: Known sites VCF (dbsnp_156.grch38.vcf.gz or All_20180418.vcf.gz)"
    echo "  (Run data/download_reference.sh to download)"
fi
echo ""

echo "=== Disk Space Check ==="
if [ -d "$DATA_DIR" ]; then
    echo "Available space in data directory:"
    df -h "$DATA_DIR" | tail -1
    AVAILABLE_GB=$(df -BG "$DATA_DIR" | tail -1 | awk '{print $4}' | sed 's/G//')
    if [ "$AVAILABLE_GB" -lt 150 ]; then
        echo "⚠ Warning: Less than 150GB available. Recommended: 200GB+"
        ((ERRORS++))
    else
        echo "✓ Sufficient disk space available"
    fi
else
    echo "✗ Data directory does not exist: $DATA_DIR"
    ((ERRORS++))
fi
echo ""

echo "=== Memory Check ==="
TOTAL_MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
echo "Total system memory: ${TOTAL_MEM_GB}GB"
if [ "$TOTAL_MEM_GB" -lt 16 ]; then
    echo "⚠ Warning: Less than 16GB RAM. Some steps may fail or be slow."
    echo "  Recommended: 16GB minimum, 32GB preferred"
    ((ERRORS++))
else
    echo "✓ Sufficient memory available"
fi
echo ""

echo "========================================="
if [ $ERRORS -eq 0 ]; then
    echo "✓ VERIFICATION PASSED"
    echo "========================================="
    echo ""
    echo "Your setup appears to be complete!"
    echo ""
    echo "Next steps:"
    echo "  1. Download reference genome:"
    echo "     cd $DATA_DIR/reference"
    echo "     bash $SCRIPT_DIR/data/download_reference.sh"
    echo ""
    echo "  2. Download sample data:"
    echo "     cd $DATA_DIR/fastq"
    echo "     bash $SCRIPT_DIR/data/download_ncbi_minimal.sh"
    echo ""
    echo "  3. Or run the complete pipeline:"
    echo "     cd $SCRIPT_DIR"
    echo "     ./run_complete_pipeline.sh"
    echo ""
    exit 0
else
    echo "✗ VERIFICATION FAILED"
    echo "========================================="
    echo ""
    echo "Found $ERRORS issue(s) that need to be addressed."
    echo ""
    echo "To set up the environment, run:"
    echo "  cd $SCRIPT_DIR"
    echo "  ./setup_start.sh"
    echo ""
    echo "Or run the complete pipeline which includes setup:"
    echo "  ./run_complete_pipeline.sh"
    echo ""
    exit 1
fi
