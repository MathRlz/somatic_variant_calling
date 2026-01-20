#!/bin/bash

# Somatic Variant Calling - Setup Script for Ubuntu 25.04
# This script installs all necessary tools and downloads required data

set -e  # Exit on error

echo "========================================="
echo "Somatic Variant Calling - Setup"
echo "========================================="

# Parse command line arguments
DATA_DIR=""
PROJECT_DIR=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --data-dir DIR       Specify directory for large data files (default: <PROJECT_DIR>/data)"
            echo "  --project-dir DIR    Specify project directory (default: current directory)"
            echo "  -h, --help           Show this help message"
            echo ""
            echo "Example:"
            echo "  $0 --data-dir ~/storage/variant_calling_data --project-dir ~/somatic_variant_calling"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Set project directory (default to current directory if not provided)
if [ -z "$PROJECT_DIR" ]; then
    PROJECT_DIR=$(pwd)
fi

# Create project directory structure
mkdir -p $PROJECT_DIR/{software,results,logs,data}

# Set data directory (default or custom)
if [ -z "$DATA_DIR" ]; then
    DATA_DIR="$PROJECT_DIR/data"
fi

# Create data directory structure
mkdir -p $DATA_DIR/{reference,fastq,bam,vcf}

cd $PROJECT_DIR

echo "Project directory: $PROJECT_DIR"
echo "Data directory: $DATA_DIR"
echo ""

# Verify data directory is accessible and has space
if [ ! -w "$DATA_DIR" ]; then
    echo "Error: Data directory $DATA_DIR is not writable"
    exit 1
fi

# Check available space
AVAILABLE_SPACE=$(df -BG "$DATA_DIR" | tail -1 | awk '{print $4}' | sed 's/G//')
echo "Available space in data directory: ${AVAILABLE_SPACE}GB"
if [ "$AVAILABLE_SPACE" -lt 150 ]; then
    echo "Warning: Less than 150GB available. Recommended: 200GB+"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# =========================================
# PART 1: INSTALL SYSTEM DEPENDENCIES
# =========================================

echo ""
echo "Installing system dependencies..."

sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    unzip \
    default-jdk \
    python3 \
    python3-pip \
    python3-venv \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    autoconf \
    automake \
    parallel \
    aria2

# =========================================
# PART 2: INSTALL BIOINFORMATICS TOOLS
# =========================================

cd $PROJECT_DIR/software

echo ""
echo "Installing bioinformatics tools..."

# --- SRA Toolkit (for downloading data from NCBI) ---
echo "Installing SRA Toolkit..."
if [ ! -d "sratoolkit.3.0.10-ubuntu64" ]; then
    rm -f sratoolkit.3.0.10-ubuntu64.tar.gz
    aria2c -x 16 -s 16 -k 1M https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz
    tar -xzf sratoolkit.3.0.10-ubuntu64.tar.gz
    rm sratoolkit.3.0.10-ubuntu64.tar.gz
fi
export PATH=$PROJECT_DIR/software/sratoolkit.3.0.10-ubuntu64/bin:$PATH

# --- SAMtools ---
echo "Installing SAMtools..."
if [ ! -d "samtools-1.19.2" ]; then
    rm -f samtools-1.19.2.tar.bz2
    aria2c -x 16 -s 16 -k 1M https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
    tar -xjf samtools-1.19.2.tar.bz2
    cd samtools-1.19.2
    ./configure --prefix=$PROJECT_DIR/software/samtools-1.19.2
    make -j4
    make install
    cd ..
    rm samtools-1.19.2.tar.bz2
fi
export PATH=$PROJECT_DIR/software/samtools-1.19.2/bin:$PATH

# --- BCFtools ---
echo "Installing BCFtools..."
if [ ! -d "bcftools-1.19" ]; then
    rm -f bcftools-1.19.tar.bz2
    aria2c -x 16 -s 16 -k 1M https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
    tar -xjf bcftools-1.19.tar.bz2
    cd bcftools-1.19
    ./configure --prefix=$PROJECT_DIR/software/bcftools-1.19
    make -j4
    make install
    cd ..
    rm bcftools-1.19.tar.bz2
fi
export PATH=$PROJECT_DIR/software/bcftools-1.19/bin:$PATH

# --- BWA (aligner) ---
echo "Installing BWA..."
if [ ! -d "bwa" ]; then
    git clone https://github.com/lh3/bwa.git
    cd bwa
    make -j4
    cd ..
fi
export PATH=$PROJECT_DIR/software/bwa:$PATH

# --- FastQC (quality control) ---
echo "Installing FastQC..."
if [ ! -d "FastQC" ]; then
    rm -f fastqc_v0.12.1.zip
    aria2c -c -x 16 -s 16 -k 1M https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
    unzip fastqc_v0.12.1.zip
    chmod +x FastQC/fastqc
    rm fastqc_v0.12.1.zip
fi
export PATH=$PROJECT_DIR/software/FastQC:$PATH

# --- Picard ---
echo "Installing Picard..."
if [ ! -f "picard.jar" ] || [ -f "picard.jar.aria2" ]; then
    aria2c -c -x 16 -s 16 -k 1M https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar
fi

# --- GATK4 ---
echo "Installing GATK4..."
if [ ! -d "gatk-4.5.0.0" ]; then
    rm -f gatk-4.5.0.0.zip
    echo "Downloading GATK4 with aria2c for speed..."
    aria2c -c -x 16 -s 16 -k 1M https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
    unzip gatk-4.5.0.0.zip
    rm gatk-4.5.0.0.zip
fi
export PATH=$PROJECT_DIR/software/gatk-4.5.0.0:$PATH

# --- VarScan2 ---
echo "Installing VarScan2..."
if [ ! -f "VarScan.v2.4.6.jar" ] || [ -f "VarScan.v2.4.6.jar.aria2" ]; then
    aria2c -c -x 16 -s 16 -k 1M https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar
fi

# --- Strelka2 ---
echo "Installing Strelka2..."
if [ ! -d "strelka-2.9.10.centos6_x86_64" ]; then
    rm -f strelka-2.9.10.centos6_x86_64.tar.bz2
    aria2c -x 16 -s 16 -k 1M https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
    tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
    rm strelka-2.9.10.centos6_x86_64.tar.bz2
fi
export PATH=$PROJECT_DIR/software/strelka-2.9.10.centos6_x86_64/bin:$PATH

# --- IGV (Integrative Genomics Viewer) ---
echo "Installing IGV..."
if [ ! -d "IGV_Linux_2.17.4" ]; then
    rm -f IGV_Linux_2.17.4_WithJava.zip
    aria2c -c -x 16 -s 16 -k 1M https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_Linux_2.17.4_WithJava.zip
    unzip IGV_Linux_2.17.4_WithJava.zip
    rm IGV_Linux_2.17.4_WithJava.zip
fi

# --- MultiQC (for aggregating QC reports) ---
echo "Installing MultiQC..."
# Create a virtual environment for Python tools to avoid system package conflicts
if [ ! -d "$PROJECT_DIR/venv" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv "$PROJECT_DIR/venv"
fi

# Install MultiQC in the virtual environment
echo "Installing MultiQC in virtual environment..."
"$PROJECT_DIR/venv/bin/pip" install --upgrade pip setuptools wheel

# Attempt to install MultiQC. We force >=1.23 to avoid the kaleido dependency issue.
# We allow failure here because MultiQC is optional (reporting only) and might have issues on Python 3.13.
if ! "$PROJECT_DIR/venv/bin/pip" install "multiqc>=1.23"; then
    echo "WARNING: MultiQC installation failed. Proceeding without it."
    echo "QC reports will be generated individually but not aggregated."
fi

# =========================================
# PART 3: CREATE ENVIRONMENT SETUP SCRIPT
# =========================================

cat > $PROJECT_DIR/setup_env.sh << EOF
#!/bin/bash
# Source this file to set up your environment
# Usage: source setup_env.sh

export PROJECT_DIR="$HOME/somatic_variant_calling"
export DATA_DIR="$DATA_DIR"

# Activate Python virtual environment
source "\$PROJECT_DIR/venv/bin/activate"

export PATH=\$PROJECT_DIR/software/sratoolkit.3.0.10-ubuntu64/bin:\$PATH
export PATH=\$PROJECT_DIR/software/samtools-1.19.2/bin:\$PATH
export PATH=\$PROJECT_DIR/software/bcftools-1.19/bin:\$PATH
export PATH=\$PROJECT_DIR/software/bwa:\$PATH
export PATH=\$PROJECT_DIR/software/FastQC:\$PATH
export PATH=\$PROJECT_DIR/software/gatk-4.5.0.0:\$PATH
export PATH=\$PROJECT_DIR/software/strelka-2.9.10.centos6_x86_64/bin:\$PATH
export PATH=\$PROJECT_DIR/software/IGV_Linux_2.17.4:\$PATH

echo "Environment set up for somatic variant calling"
echo "Project directory: \$PROJECT_DIR"
echo "Data directory: \$DATA_DIR"
echo "Python virtual environment activated"
EOF

chmod +x $PROJECT_DIR/setup_env.sh

echo ""
echo "========================================="
echo "Software installation complete!"
echo "========================================="
echo ""
echo "Configuration:"
echo "  Project directory: $PROJECT_DIR"
echo "  Data directory: $DATA_DIR"
echo ""
echo "To use the tools, run:"
echo "  source $PROJECT_DIR/setup_env.sh"
echo ""
