#!/bin/bash

# Minimal script to download NCBI run information and create sample summary
# Usage: ./download_ncbi_minimal.sh [PROJECT_ID] [OUTPUT_DIR] [--force]

set -e

PROJECT_ID=${1:-"PRJNA284596"}
OUTPUT_DIR=${2:-$(pwd)}
FORCE_REFRESH=false

# Check for --force flag
for arg in "$@"; do
    if [ "$arg" = "--force" ] || [ "$arg" = "-f" ]; then
        FORCE_REFRESH=true
    fi
done

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Check if runinfo.csv already exists
RUNINFO_FILE="runinfo.csv"
if [ -f "$RUNINFO_FILE" ] && [ "$FORCE_REFRESH" = false ]; then
    # Check if the file contains data for the same project
    if head -1 "$RUNINFO_FILE" | grep -q "Run"; then
        echo "Run information already exists at: $OUTPUT_DIR/$RUNINFO_FILE"
        echo "Use --force or -f to re-download."
        echo ""
    else
        FORCE_REFRESH=true
    fi
fi

if [ ! -f "$RUNINFO_FILE" ] || [ "$FORCE_REFRESH" = true ]; then
    echo "Downloading run information for $PROJECT_ID..."

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

    echo "Run information saved to: runinfo.csv"
fi

# Create sample summary (regenerate if runinfo.csv is newer or summary doesn't exist)
SUMMARY_FILE="sample_summary.txt"
if [ ! -f "$SUMMARY_FILE" ] || [ "$RUNINFO_FILE" -nt "$SUMMARY_FILE" ] || [ "$FORCE_REFRESH" = true ]; then
    echo "Creating sample summary..."
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
    echo "Sample summary saved to: sample_summary.txt"
else
    echo "Sample summary already exists at: $OUTPUT_DIR/$SUMMARY_FILE"
fi

echo ""
cat sample_summary.txt
