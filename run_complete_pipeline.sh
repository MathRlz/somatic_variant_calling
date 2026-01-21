#!/bin/bash

# Somatic Variant Calling - Complete Pipeline Runner
# This script runs the entire workflow from data download to final analysis

set -e  # Exit on error

# Configuration
NUM_PROCESSORS=8  # Number of processors/threads to use for parallel operations
export NUM_PROCESSORS

# Parse command line arguments
SKIP_DOWNLOAD=false
SKIP_ALIGNMENT=false
SKIP_SETUP=false
DATA_DIR=""
PROJECT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-alignment)
            SKIP_ALIGNMENT=true
            shift
            ;;
        --skip-setup)
            SKIP_SETUP=true
            shift
            ;;
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
            echo "Complete somatic variant calling pipeline from setup to analysis"
            echo ""
            echo "Options:"
            echo "  --skip-download       Skip reference and sample download (use existing data)"
            echo "  --skip-alignment      Skip alignment step (use existing BAM files)"
            echo "  --skip-setup          Skip environment setup (assumes already set up)"
            echo "  --data-dir DIR        Specify directory for data files (default: ~/storage/variant_calling_data)"
            echo "  --project-dir DIR     Specify project directory (default: ~/somatic_variant_calling)"
            echo "  -h, --help            Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Full pipeline from scratch:"
            echo "  $0"
            echo ""
            echo "  # Skip download if you already have data:"
            echo "  $0 --skip-download"
            echo ""
            echo "  # Skip alignment if you already have BAM files:"
            echo "  $0 --skip-download --skip-alignment"
            echo ""
            echo "  # Use custom directories:"
            echo "  $0 --data-dir ~/my_data --project-dir ~/my_project"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Set default directories
if [ -z "$PROJECT_DIR" ]; then
    PROJECT_DIR="$HOME/somatic_variant_calling"
fi

if [ -z "$DATA_DIR" ]; then
    DATA_DIR="$HOME/somatic_variant_calling/data"
fi

# Export DATA_DIR so subscripts can use it
export DATA_DIR

echo "========================================="
echo "Somatic Variant Calling - Complete Pipeline"
echo "========================================="
echo ""
echo "Configuration:"
echo "  Project directory: $PROJECT_DIR"
echo "  Data directory: $DATA_DIR"
echo "  Number of processors: $NUM_PROCESSORS"
echo "  Skip setup: $SKIP_SETUP"
echo "  Skip download: $SKIP_DOWNLOAD"
echo "  Skip alignment: $SKIP_ALIGNMENT"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# =========================================
# STEP 0: SETUP ENVIRONMENT
# =========================================

if [ "$SKIP_SETUP" = false ]; then
    echo "========================================="
    echo "STEP 0: Setting up environment"
    echo "========================================="

    if [ ! -f "$SCRIPT_DIR/setup_start.sh" ]; then
        echo "Error: setup_start.sh not found in $SCRIPT_DIR"
        exit 1
    fi

    bash "$SCRIPT_DIR/setup_start.sh" --data-dir "$DATA_DIR" --project-dir "$PROJECT_DIR"
    echo ""
fi

# Source the environment
if [ -f "$PROJECT_DIR/setup_env.sh" ]; then
    echo "Sourcing environment from $PROJECT_DIR/setup_env.sh"
    source "$PROJECT_DIR/setup_env.sh"
else
    echo "Warning: $PROJECT_DIR/setup_env.sh not found. Assuming tools are in PATH."
fi

# Change to data directory for pipeline execution
cd "$DATA_DIR"

# =========================================
# STEP 1: DOWNLOAD REFERENCE GENOME
# =========================================

if [ "$SKIP_DOWNLOAD" = false ]; then
    echo "========================================="
    echo "STEP 1: Downloading reference genome"
    echo "========================================="

    if [ ! -f "$SCRIPT_DIR/data/download_reference.sh" ]; then
        echo "Error: download_reference.sh not found"
        exit 1
    fi

    cd "$DATA_DIR/reference"
    bash "$SCRIPT_DIR/data/download_reference.sh"
    cd "$DATA_DIR"
    echo ""
fi

# =========================================
# STEP 2: DOWNLOAD SAMPLE DATA
# =========================================

if [ "$SKIP_DOWNLOAD" = false ]; then
    echo "========================================="
    echo "STEP 2: Downloading sample data"
    echo "========================================="
    echo ""

    if [ ! -f "$SCRIPT_DIR/data/download_samples_interactive.sh" ]; then
        echo "Warning: download_samples_interactive.sh not found, skipping sample download"
        echo "You can download samples manually using data/download_sample.sh"
    else
        cd "$DATA_DIR/fastq"
        bash "$SCRIPT_DIR/data/download_samples_interactive.sh"
        cd "$DATA_DIR"
    fi
    echo ""
fi

# =========================================
# PATIENT SELECTION
# =========================================

echo "========================================="
echo "Patient Selection"
echo "========================================="
echo ""

# Function to get available patients from FASTQ files
get_available_patients() {
    local fastq_dir="$DATA_DIR/fastq"
    if [ -d "$fastq_dir" ]; then
        # Find unique patient IDs from FASTQ files (pattern: PATIENT-T_R1.fastq.gz or PATIENT-N_R1.fastq.gz)
        ls "$fastq_dir"/*_R1.fastq.gz 2>/dev/null | \
            xargs -n1 basename 2>/dev/null | \
            sed 's/_R1.fastq.gz$//' | \
            sed 's/-[TN]$//' | \
            sort -u
    fi
}

# Function to check if patient has both tumor and normal samples
check_patient_complete() {
    local patient=$1
    local fastq_dir="$DATA_DIR/fastq"
    if [ -f "$fastq_dir/${patient}-T_R1.fastq.gz" ] && [ -f "$fastq_dir/${patient}-N_R1.fastq.gz" ]; then
        return 0
    fi
    return 1
}

# Get list of available patients
AVAILABLE_PATIENTS=($(get_available_patients))

if [ ${#AVAILABLE_PATIENTS[@]} -eq 0 ]; then
    echo "No patient data found in $DATA_DIR/fastq"
    echo "Please download sample data first."
    exit 1
fi

echo "Available patients:"
echo "-------------------"
i=1
for patient in "${AVAILABLE_PATIENTS[@]}"; do
    # Check if patient has both tumor and normal
    if check_patient_complete "$patient"; then
        status="[T+N]"
    elif [ -f "$DATA_DIR/fastq/${patient}-T_R1.fastq.gz" ]; then
        status="[T only]"
    else
        status="[N only]"
    fi
    echo "  $i) $patient $status"
    i=$((i+1))
done
echo ""

echo "Options:"
echo "  a) Process ALL patients"
echo "  Enter numbers separated by spaces to select specific patients"
echo ""
read -p "Enter your choice: " PATIENT_CHOICE

# Process selection
SELECTED_PATIENTS=()
if [ "$PATIENT_CHOICE" = "a" ] || [ "$PATIENT_CHOICE" = "A" ] || [ "$PATIENT_CHOICE" = "all" ]; then
    SELECTED_PATIENTS=("${AVAILABLE_PATIENTS[@]}")
    echo "Selected: ALL patients"
else
    for num in $PATIENT_CHOICE; do
        if [[ "$num" =~ ^[0-9]+$ ]] && [ "$num" -ge 1 ] && [ "$num" -le ${#AVAILABLE_PATIENTS[@]} ]; then
            SELECTED_PATIENTS+=("${AVAILABLE_PATIENTS[$((num-1))]}")
        fi
    done
fi

if [ ${#SELECTED_PATIENTS[@]} -eq 0 ]; then
    echo "Error: No valid patients selected."
    exit 1
fi

echo ""
echo "Selected patients for processing:"
for patient in "${SELECTED_PATIENTS[@]}"; do
    echo "  - $patient"
done
echo ""

# Save selected patients to file for use by other scripts
SELECTED_PATIENTS_FILE="$DATA_DIR/selected_patients.txt"
printf "%s\n" "${SELECTED_PATIENTS[@]}" > "$SELECTED_PATIENTS_FILE"
export SELECTED_PATIENTS_FILE

echo "Patient selection saved to: $SELECTED_PATIENTS_FILE"
echo ""

# =========================================
# STEP 3: ALIGN READS TO REFERENCE
# =========================================

if [ "$SKIP_ALIGNMENT" = false ]; then
    echo "========================================="
    echo "STEP 3: Aligning reads to reference genome"
    echo "========================================="

    if [ ! -f "$SCRIPT_DIR/data/align_all.sh" ]; then
        echo "Error: align_all.sh not found"
        exit 1
    fi

    cd "$DATA_DIR/fastq"
    bash "$SCRIPT_DIR/data/align_all.sh"

    # Move BAM files to bam directory
    mkdir -p "$DATA_DIR/bam"
    mv *.bam *.bam.bai "$DATA_DIR/bam/" 2>/dev/null || true
    cd "$DATA_DIR"
    echo ""
fi

# =========================================
# STEP 4: MARK DUPLICATES
# =========================================

echo "========================================="
echo "STEP 4: Marking duplicate reads"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/mark_duplicates.sh" ]; then
    echo "Error: mark_duplicates.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/mark_duplicates.sh"
echo ""

# =========================================
# STEP 5: ADD READ GROUPS
# =========================================

echo "========================================="
echo "STEP 5: Adding read groups"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/add_read_groups.sh" ]; then
    echo "Error: add_read_groups.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/add_read_groups.sh"
echo ""

# =========================================
# STEP 6: BASE QUALITY SCORE RECALIBRATION
# =========================================

echo "========================================="
echo "STEP 6: Base quality score recalibration (BQSR)"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/bqsr.sh" ]; then
    echo "Error: bqsr.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/bqsr.sh"
echo ""

# =========================================
# STEP 7: CALL VARIANTS
# =========================================

echo "========================================="
echo "STEP 7: Calling variants"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/call_variants.sh" ]; then
    echo "Error: call_variants.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/call_variants.sh"
echo ""

# =========================================
# STEP 8: FILTER VARIANTS
# =========================================

echo "========================================="
echo "STEP 8: Filtering variants"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/filter_variants.sh" ]; then
    echo "Error: filter_variants.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/filter_variants.sh"
echo ""

# =========================================
# STEP 9: ANNOTATE VARIANTS
# =========================================

echo "========================================="
echo "STEP 9: Annotating variants"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/annotate_variants.sh" ]; then
    echo "Error: annotate_variants.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/annotate_variants.sh"
echo ""

# =========================================
# STEP 10: PRIORITIZE VARIANTS
# =========================================

echo "========================================="
echo "STEP 10: Prioritizing variants"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/prioritize_variants.sh" ]; then
    echo "Error: prioritize_variants.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/prioritize_variants.sh"
echo ""

# =========================================
# STEP 11: GENERATE REPORTS
# =========================================

echo "========================================="
echo "STEP 11: Generating reports"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/generate_reports.sh" ]; then
    echo "Error: generate_reports.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/generate_reports.sh"
echo ""

# =========================================
# STEP 12: COMPARE SAMPLES
# =========================================

echo "========================================="
echo "STEP 12: Comparing tumor and normal samples"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/compare_samples.sh" ]; then
    echo "Error: compare_samples.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/compare_samples.sh"
echo ""

# =========================================
# STEP 13: CREATE IGV SESSION
# =========================================

echo "========================================="
echo "STEP 13: Creating IGV visualization session"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/create_igv_session.sh" ]; then
    echo "Error: create_igv_session.sh not found"
    exit 1
fi

cd "$DATA_DIR"
bash "$SCRIPT_DIR/variant_calling/create_igv_session.sh"
echo ""

# =========================================
# PIPELINE COMPLETE
# =========================================

echo "========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY!"
echo "========================================="
echo ""
echo "Results are available in: $DATA_DIR"
echo ""
echo "Key output directories:"
echo "  - bam/                 : Aligned reads"
echo "  - bam_marked/          : Deduplicated BAM files"
echo "  - bam_with_rg/         : BAM files with read groups"
echo "  - bam_recalibrated/    : Quality-recalibrated BAM files"
echo "  - vcfs/                : Raw variant calls"
echo "  - vcfs_filtered/       : Filtered variants"
echo "  - vcfs_annotated/      : Annotated variants"
echo "  - vcfs_prioritized/    : High-confidence variants"
echo "  - reports/             : Summary statistics"
echo "  - comparisons/         : Tumor-specific variants"
echo "  - igv_session.xml      : IGV visualization file"
echo ""
echo "Next steps:"
echo "  1. Review reports in: $DATA_DIR/reports/"
echo "  2. Examine tumor-specific variants in: $DATA_DIR/comparisons/"
echo "  3. Open IGV and load: $DATA_DIR/igv_session.xml"
echo "  4. Explore high-confidence variants in: $DATA_DIR/vcfs_prioritized/"
echo ""
echo "To visualize results in IGV:"
echo "  cd $PROJECT_DIR/software/IGV_Linux_2.17.4"
echo "  ./igv.sh"
echo "  File > Open Session > Navigate to $DATA_DIR/igv_session.xml"
echo ""
