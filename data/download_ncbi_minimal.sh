#!/bin/bash

# Minimal script to download NCBI run information and create sample summary
# Usage: ./download_ncbi_minimal.sh [PROJECT_ID] [OUTPUT_DIR]

set -e

PROJECT_ID=${1:-"PRJNA284596"}
OUTPUT_DIR=${2:-$(pwd)}

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

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

echo "Sample summary saved to: sample_summary.txt"
echo ""
cat sample_summary.txt
