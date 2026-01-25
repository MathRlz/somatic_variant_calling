#!/bin/bash

# Somatic Variant Calling - Complete Pipeline Runner
# This script runs the entire workflow from data download to final analysis

set -e  # Exit on error

# Configuration
NUM_PROCESSORS=8  # Number of processors/threads to use for parallel operations (default)

# Parse command line arguments
SKIP_DOWNLOAD=false
SKIP_ALIGNMENT=false
SKIP_SETUP=false
FORCE_REGENERATE=false
DATA_DIR=""
PROJECT_DIR=""
MIN_DP="20"
MIN_VAF="0.05"

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
        --force|-f)
            FORCE_REGENERATE=true
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
        --num-processors)
            NUM_PROCESSORS="$2"
            shift 2
            ;;
        --min-dp)
            MIN_DP="$2"
            shift 2
            ;;
        --min-vaf)
            MIN_VAF="$2"
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
            echo "  --force, -f           Force regeneration of all output files (skip prompts)"
            echo "  --data-dir DIR        Specify directory for data files (default: ~/somatic_variant_calling/data)"
            echo "  --project-dir DIR     Specify project directory (default: ~/somatic_variant_calling)"
            echo "  --num-processors N    Number of processors/threads to use (default: 8)"
            echo "  --min-dp N            Minimum read depth for variant prioritization (default: 20)"
            echo "  --min-vaf N           Minimum variant allele frequency (default: 0.05)"
            echo "  -h, --help            Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Full pipeline from scratch:"
            echo "  $0"
            echo ""
            echo "  # Skip download if you already have data:"
            echo "  $0 --skip-download"
            echo ""
            echo "  # Force regenerate all outputs:"
            echo "  $0 --skip-download --force"
            echo ""
            echo "  # Use custom directories:"
            echo "  $0 --data-dir ~/my_data --project-dir ~/my_project"
            echo ""
            echo "  # Use 16 processors:"
            echo "  $0 --num-processors 16"
            echo ""
            echo "  # Less strict variant filtering (DP>=10):"
            echo "  $0 --skip-download --min-dp 10"
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

# Export variables so subscripts can use them
export DATA_DIR
export NUM_PROCESSORS
export MIN_DP
export MIN_VAF

echo "========================================="
echo "Somatic Variant Calling - Complete Pipeline"
echo "========================================="
echo ""
echo "Configuration:"
echo "  Project directory: $PROJECT_DIR"
echo "  Data directory: $DATA_DIR"
echo "  Number of processors: $NUM_PROCESSORS"
echo "  Min read depth (prioritization): $MIN_DP"
echo "  Min VAF (prioritization): $MIN_VAF"
echo "  Skip setup: $SKIP_SETUP"
echo "  Skip download: $SKIP_DOWNLOAD"
echo "  Skip alignment: $SKIP_ALIGNMENT"
echo "  Force regenerate: $FORCE_REGENERATE"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# =========================================
# HELPER FUNCTION: Check if step output exists
# =========================================

check_step_output() {
    local step_name=$1
    local output_file=$2
    local identifier=$3  # sample or patient name

    if [ -f "$output_file" ]; then
        if [ "$FORCE_REGENERATE" = true ]; then
            echo "Regenerating $step_name for $identifier (--force)"
            return 0  # Run the step
        fi
        echo "Output for $step_name ($identifier) already exists: $output_file"
        read -p "Regenerate? [y/N]: " choice
        if [[ "$choice" =~ ^[Yy]$ ]]; then
            return 0  # Run the step
        else
            echo "Skipping $step_name for $identifier"
            return 1  # Skip the step
        fi
    fi
    return 0  # File doesn't exist, run the step
}

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

# =========================================
# PROCESS EACH SELECTED PATIENT
# =========================================

for patient in "${SELECTED_PATIENTS[@]}"; do
    echo ""
    echo "========================================="
    echo "Processing patient: $patient"
    echo "========================================="
    echo ""

    # =========================================
    # SAMPLE-LEVEL STEPS (both -T and -N samples)
    # =========================================

    for sample in "${patient}-T" "${patient}-N"; do
        # Check if sample FASTQ exists
        if [ ! -f "$DATA_DIR/fastq/${sample}_R1.fastq.gz" ]; then
            echo "Skipping $sample (FASTQ file not found)"
            continue
        fi

        echo ""
        echo "--- Processing sample: $sample ---"
        echo ""

        # =========================================
        # STEP 3: ALIGN READS TO REFERENCE
        # =========================================

        if [ "$SKIP_ALIGNMENT" = false ]; then
            output="$DATA_DIR/bam/${sample}.bam"
            if check_step_output "Alignment" "$output" "$sample"; then
                echo "STEP 3: Aligning reads for $sample..."

                if [ ! -f "$SCRIPT_DIR/data/align_sample.sh" ]; then
                    echo "Error: align_sample.sh not found"
                    exit 1
                fi

                bash "$SCRIPT_DIR/data/align_sample.sh" "$sample"
            fi
        fi

        # =========================================
        # STEP 4: MARK DUPLICATES
        # =========================================

        output="$DATA_DIR/bam_marked/${sample}_marked.bam"
        if check_step_output "Mark Duplicates" "$output" "$sample"; then
            echo "STEP 4: Marking duplicates for $sample..."

            if [ ! -f "$SCRIPT_DIR/variant_calling/mark_duplicates.sh" ]; then
                echo "Error: mark_duplicates.sh not found"
                exit 1
            fi

            bash "$SCRIPT_DIR/variant_calling/mark_duplicates.sh" "$sample"
        fi

        # =========================================
        # STEP 5: ADD READ GROUPS
        # =========================================

        output="$DATA_DIR/bam_with_rg/${sample}_rg.bam"
        if check_step_output "Add Read Groups" "$output" "$sample"; then
            echo "STEP 5: Adding read groups for $sample..."

            if [ ! -f "$SCRIPT_DIR/variant_calling/add_read_groups.sh" ]; then
                echo "Error: add_read_groups.sh not found"
                exit 1
            fi

            bash "$SCRIPT_DIR/variant_calling/add_read_groups.sh" "$sample"
        fi

        # =========================================
        # STEP 6: BASE QUALITY SCORE RECALIBRATION
        # =========================================

        output="$DATA_DIR/bam_recalibrated/${sample}_recalibrated.bam"
        if check_step_output "BQSR" "$output" "$sample"; then
            echo "STEP 6: Recalibrating base quality scores for $sample..."

            if [ ! -f "$SCRIPT_DIR/variant_calling/bqsr.sh" ]; then
                echo "Error: bqsr.sh not found"
                exit 1
            fi

            bash "$SCRIPT_DIR/variant_calling/bqsr.sh" "$sample"
        fi

    done  # End sample loop

    # =========================================
    # PATIENT-LEVEL STEPS (tumor-normal pairs)
    # =========================================

    # Check if we have both tumor and normal for variant calling
    tumor_bam="$DATA_DIR/bam_recalibrated/${patient}-T_recalibrated.bam"
    normal_bam="$DATA_DIR/bam_recalibrated/${patient}-N_recalibrated.bam"

    if [ ! -f "$tumor_bam" ] || [ ! -f "$normal_bam" ]; then
        echo ""
        echo "Warning: Skipping patient-level steps for $patient"
        echo "  Both tumor and normal recalibrated BAMs are required."
        echo "  Tumor BAM: $tumor_bam (exists: $([ -f "$tumor_bam" ] && echo yes || echo no))"
        echo "  Normal BAM: $normal_bam (exists: $([ -f "$normal_bam" ] && echo yes || echo no))"
        continue
    fi

    echo ""
    echo "--- Patient-level steps for: $patient ---"
    echo ""

    # =========================================
    # STEP 7: CALL VARIANTS
    # =========================================

    output="$DATA_DIR/vcfs/${patient}_raw.vcf.gz"
    if check_step_output "Variant Calling" "$output" "$patient"; then
        echo "STEP 7: Calling variants for $patient..."

        if [ ! -f "$SCRIPT_DIR/variant_calling/call_variants.sh" ]; then
            echo "Error: call_variants.sh not found"
            exit 1
        fi

        bash "$SCRIPT_DIR/variant_calling/call_variants.sh" "$patient"
    fi

    # =========================================
    # STEP 8: FILTER VARIANTS
    # =========================================

    output="$DATA_DIR/vcfs_filtered/${patient}_filtered.vcf.gz"
    if check_step_output "Variant Filtering" "$output" "$patient"; then
        echo "STEP 8: Filtering variants for $patient..."

        if [ ! -f "$SCRIPT_DIR/variant_calling/filter_variants.sh" ]; then
            echo "Error: filter_variants.sh not found"
            exit 1
        fi

        bash "$SCRIPT_DIR/variant_calling/filter_variants.sh" "$patient"
    fi

    # =========================================
    # STEP 9: ANNOTATE VARIANTS
    # =========================================

    output="$DATA_DIR/vcfs_annotated/${patient}_annotated.vcf.gz"
    if check_step_output "Variant Annotation" "$output" "$patient"; then
        echo "STEP 9: Annotating variants for $patient..."

        if [ ! -f "$SCRIPT_DIR/variant_calling/annotate_variants.sh" ]; then
            echo "Error: annotate_variants.sh not found"
            exit 1
        fi

        bash "$SCRIPT_DIR/variant_calling/annotate_variants.sh" "$patient"
    fi

    # =========================================
    # STEP 10: PRIORITIZE VARIANTS
    # =========================================

    output="$DATA_DIR/vcfs_prioritized/${patient}_high_confidence.vcf.gz"
    if check_step_output "Variant Prioritization" "$output" "$patient"; then
        echo "STEP 10: Prioritizing variants for $patient..."

        if [ ! -f "$SCRIPT_DIR/variant_calling/prioritize_variants.sh" ]; then
            echo "Error: prioritize_variants.sh not found"
            exit 1
        fi

        bash "$SCRIPT_DIR/variant_calling/prioritize_variants.sh" "$patient"
    fi

    # =========================================
    # STEP 11: GENERATE REPORTS
    # =========================================

    output="$DATA_DIR/reports/${patient}_stats.txt"
    if check_step_output "Report Generation" "$output" "$patient"; then
        echo "STEP 11: Generating reports for $patient..."

        if [ ! -f "$SCRIPT_DIR/variant_calling/generate_reports.sh" ]; then
            echo "Error: generate_reports.sh not found"
            exit 1
        fi

        bash "$SCRIPT_DIR/variant_calling/generate_reports.sh" "$patient"
    fi

done  # End patient loop

# =========================================
# STEP 12: CREATE IGV SESSION (all patients)
# =========================================

echo ""
echo "========================================="
echo "STEP 12: Creating IGV visualization session"
echo "========================================="

if [ ! -f "$SCRIPT_DIR/variant_calling/create_igv_session.sh" ]; then
    echo "Error: create_igv_session.sh not found"
    exit 1
fi

bash "$SCRIPT_DIR/variant_calling/create_igv_session.sh" "${SELECTED_PATIENTS[@]}"
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
echo "  - igv_session.xml      : IGV visualization file"
echo ""
echo "Next steps:"
echo "  1. Review reports in: $DATA_DIR/reports/"
echo "  2. Explore high-confidence variants in: $DATA_DIR/vcfs_prioritized/"
echo "  3. Open IGV and load: $DATA_DIR/igv_session.xml"
echo ""
echo "To visualize results in IGV:"
echo "  cd $PROJECT_DIR/software/IGV_Linux_2.17.4"
echo "  ./igv.sh"
echo "  File > Open Session > Navigate to $DATA_DIR/igv_session.xml"
echo ""
