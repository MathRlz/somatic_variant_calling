#!/bin/bash

# Interactive script to download NCBI samples with user selection
# Usage: ./download_samples_interactive.sh [PROJECT_ID] [OUTPUT_DIR]

set -e

PROJECT_ID=${1:-"PRJNA284596"}
OUTPUT_DIR=${2:-$(pwd)}

echo "========================================="
echo "Interactive Sample Download"
echo "========================================="
echo ""
echo "Project ID: $PROJECT_ID"
echo "Output Directory: $OUTPUT_DIR"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Step 1: Download run information and create summary
echo "Step 1: Fetching available samples from NCBI..."
echo "----------------------------------------"

# Fetch run info using NCBI E-utilities
SEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$PROJECT_ID&usehistory=y"
SEARCH_RESULT=$(curl -s "$SEARCH_URL")

WEBENV=$(echo "$SEARCH_RESULT" | grep -oP '(?<=<WebEnv>).*?(?=</WebEnv>)')
QUERYKEY=$(echo "$SEARCH_RESULT" | grep -oP '(?<=<QueryKey>).*?(?=</QueryKey>)')

if [ -z "$WEBENV" ] || [ -z "$QUERYKEY" ]; then
    echo "Error: Failed to retrieve search results from NCBI."
    exit 1
fi

FETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$QUERYKEY&WebEnv=$WEBENV&rettype=runinfo&retmode=text"
curl -s "$FETCH_URL" > runinfo.csv

echo "Run information downloaded successfully."
echo ""

# Create sample summary
awk -F',' 'NR>1 {
    run = $1
    sample = $30
    lib = $13
    spots = $5
    bases = $6

    if (sample ~ /-T$/ || lib ~ /tumor|cancer|T/ && lib !~ /normal|blood|germline/) {
        type = "TUMOR"
    } else if (sample ~ /-N$/ || lib ~ /normal|blood|germline|N/) {
        type = "NORMAL"
    } else {
        type = "UNKNOWN"
    }

    printf "%s\t%s\t%s\t%s\t%s\n", sample, type, run, spots, bases
}' runinfo.csv > sample_summary.txt

# Step 2: Display available samples
echo "Step 2: Available Samples"
echo "----------------------------------------"
echo ""
echo "Sample Name          Type      SRR ID        Spots       Bases"
echo "================================================================"
cat sample_summary.txt | column -t
echo ""

# Count samples by type
TUMOR_COUNT=$(grep -c "TUMOR" sample_summary.txt || echo "0")
NORMAL_COUNT=$(grep -c "NORMAL" sample_summary.txt || echo "0")
TOTAL_COUNT=$(wc -l < sample_summary.txt)

echo "Summary: $TOTAL_COUNT samples ($TUMOR_COUNT tumor, $NORMAL_COUNT normal)"
echo ""

# Step 3: Prompt user for selection
echo "Step 3: Select Samples to Download"
echo "----------------------------------------"
echo ""
echo "Options:"
echo "  1) Download ALL samples"
echo "  2) Download specific tumor-normal pairs"
echo "  3) Download individual samples by SRR ID"
echo "  4) Skip download (just create summary)"
echo ""
read -p "Enter your choice (1-4): " CHOICE

case $CHOICE in
    1)
        echo ""
        echo "Downloading ALL samples..."
        echo ""

        # Extract all SRR IDs and sample names
        while IFS=$'\t' read -r sample_name type srr_id spots bases; do
            if [ -n "$srr_id" ] && [ "$srr_id" != "SRR_ID" ]; then
                echo "Downloading $sample_name ($type) - $srr_id..."
                bash "$SCRIPT_DIR/download_sample.sh" "$srr_id" "$sample_name"
                echo ""
            fi
        done < sample_summary.txt
        ;;

    2)
        echo ""
        echo "Available tumor-normal pairs:"
        echo ""

        # Extract unique patient IDs (prefix before -T or -N)
        PATIENTS=$(awk -F'\t' '{split($1, a, /-[TN]$/); print a[1]}' sample_summary.txt | sort -u)

        declare -a PATIENT_ARRAY
        i=1
        for patient in $PATIENTS; do
            PATIENT_ARRAY[$i]=$patient
            echo "  $i) $patient"
            i=$((i+1))
        done
        echo ""
        read -p "Enter patient number(s) to download (space-separated, or 'all'): " PATIENT_SELECTION

        if [ "$PATIENT_SELECTION" = "all" ]; then
            SELECTED_PATIENTS=("${PATIENT_ARRAY[@]}")
        else
            # Convert user selection to patient IDs
            SELECTED_PATIENTS=()
            for num in $PATIENT_SELECTION; do
                if [ -n "${PATIENT_ARRAY[$num]}" ]; then
                    SELECTED_PATIENTS+=("${PATIENT_ARRAY[$num]}")
                fi
            done
        fi

        echo ""
        echo "Downloading selected patients..."
        echo ""

        for patient in "${SELECTED_PATIENTS[@]}"; do
            # Download tumor sample
            TUMOR_LINE=$(grep "^${patient}-T" sample_summary.txt || echo "")
            if [ -n "$TUMOR_LINE" ]; then
                TUMOR_SRR=$(echo "$TUMOR_LINE" | awk -F'\t' '{print $3}')
                TUMOR_NAME=$(echo "$TUMOR_LINE" | awk -F'\t' '{print $1}')
                echo "Downloading $TUMOR_NAME (TUMOR) - $TUMOR_SRR..."
                bash "$SCRIPT_DIR/download_sample.sh" "$TUMOR_SRR" "$TUMOR_NAME"
                echo ""
            fi

            # Download normal sample
            NORMAL_LINE=$(grep "^${patient}-N" sample_summary.txt || echo "")
            if [ -n "$NORMAL_LINE" ]; then
                NORMAL_SRR=$(echo "$NORMAL_LINE" | awk -F'\t' '{print $3}')
                NORMAL_NAME=$(echo "$NORMAL_LINE" | awk -F'\t' '{print $1}')
                echo "Downloading $NORMAL_NAME (NORMAL) - $NORMAL_SRR..."
                bash "$SCRIPT_DIR/download_sample.sh" "$NORMAL_SRR" "$NORMAL_NAME"
                echo ""
            fi
        done
        ;;

    3)
        echo ""
        read -p "Enter SRR ID(s) to download (space-separated): " SRR_IDS

        echo ""
        echo "Downloading selected samples..."
        echo ""

        for srr in $SRR_IDS; do
            # Find sample name from summary
            SAMPLE_LINE=$(grep "$srr" sample_summary.txt || echo "")
            if [ -n "$SAMPLE_LINE" ]; then
                SAMPLE_NAME=$(echo "$SAMPLE_LINE" | awk -F'\t' '{print $1}')
                SAMPLE_TYPE=$(echo "$SAMPLE_LINE" | awk -F'\t' '{print $2}')
                echo "Downloading $SAMPLE_NAME ($SAMPLE_TYPE) - $srr..."
            else
                SAMPLE_NAME=$srr
                echo "Downloading $srr..."
            fi
            bash "$SCRIPT_DIR/download_sample.sh" "$srr" "$SAMPLE_NAME"
            echo ""
        done
        ;;

    4)
        echo ""
        echo "Skipping download. Summary saved to: $OUTPUT_DIR/sample_summary.txt"
        echo "You can download samples later using:"
        echo "  bash $SCRIPT_DIR/download_sample.sh SRR_ID SAMPLE_NAME"
        exit 0
        ;;

    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo ""
echo "========================================="
echo "Download Complete!"
echo "========================================="
echo ""
echo "Summary file: $OUTPUT_DIR/sample_summary.txt"
echo "FASTQ files: $OUTPUT_DIR/*.fastq.gz"
echo ""
