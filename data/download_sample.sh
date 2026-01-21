#!/bin/bash

# Usage: ./download_sample.sh SRR_ID [SAMPLE_NAME]

if [ $# -lt 1 ]; then
    echo "Usage: $0 SRR_ID [SAMPLE_NAME]"
    echo "Example: $0 SRR2057563 tumor_sample"
    exit 1
fi

SRR_ID=$1
SAMPLE_NAME=${2:-$SRR_ID}

# Function to check if FASTQ files are already downloaded
check_existing_fastq() {
    local sample=$1
    # Check for paired-end files
    if [ -f "${sample}_R1.fastq.gz" ] && [ -f "${sample}_R2.fastq.gz" ]; then
        # Verify files are not incomplete (no .aria2 files)
        if [ ! -f "${sample}_R1.fastq.gz.aria2" ] && [ ! -f "${sample}_R2.fastq.gz.aria2" ]; then
            echo "paired"
            return 0
        fi
    fi
    # Check for single-end file
    if [ -f "${sample}.fastq.gz" ]; then
        if [ ! -f "${sample}.fastq.gz.aria2" ]; then
            echo "single"
            return 0
        fi
    fi
    echo "none"
    return 1
}

# Function to check if FastQC reports exist
check_existing_fastqc() {
    local sample=$1
    local mode=$2
    if [ "$mode" = "paired" ]; then
        if [ -f "${sample}_R1_fastqc.html" ] && [ -f "${sample}_R2_fastqc.html" ]; then
            return 0
        fi
    else
        if [ -f "${sample}_fastqc.html" ]; then
            return 0
        fi
    fi
    return 1
}

# Check if FASTQ files already exist
EXISTING_MODE=$(check_existing_fastq "$SAMPLE_NAME")
if [ "$EXISTING_MODE" != "none" ]; then
    echo "FASTQ files for $SAMPLE_NAME already exist (${EXISTING_MODE}-end). Skipping download."

    # Check if FastQC needs to be run
    if check_existing_fastqc "$SAMPLE_NAME" "$EXISTING_MODE"; then
        echo "FastQC reports already exist. Skipping quality control."
    else
        echo "Running quality control on existing files..."
        if [ "$EXISTING_MODE" = "paired" ]; then
            fastqc ${SAMPLE_NAME}_R1.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz
        else
            fastqc ${SAMPLE_NAME}.fastq.gz
        fi
    fi
    echo "Sample $SAMPLE_NAME is ready!"
    exit 0
fi

echo "Fetching ENA FASTQ URLs for $SRR_ID..."

# Query ENA for fastq HTTPS and FTP links (JSON output)
ENA_API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR_ID&result=read_run&fields=fastq_https,fastq_ftp,fastq_bytes&format=json"
FASTQ_INFO=$(curl -s "$ENA_API_URL")

FASTQ_HTTPS=$(echo "$FASTQ_INFO" | grep -oP '"fastq_https":"\K[^"]*')
FASTQ_FTP=$(echo "$FASTQ_INFO" | grep -oP '"fastq_ftp":"\K[^"]*')

if [ -n "$FASTQ_HTTPS" ]; then
    URLS_RAW="$FASTQ_HTTPS"
    PROTO="https"
    echo "Found HTTPS links for $SRR_ID."
elif [ -n "$FASTQ_FTP" ]; then
    URLS_RAW="$FASTQ_FTP"
    PROTO="ftp"
    echo "HTTPS links not found, falling back to FTP for $SRR_ID."

else
    echo "ENA does not provide fastq files for $SRR_ID. Falling back to NCBI prefetch + fastq-dump."

    # Clean up any stale lock files from previous interrupted downloads
    SRA_DIR="$PWD/$SRR_ID"
    if [ -d "$SRA_DIR" ]; then
        LOCK_FILE="$SRA_DIR/${SRR_ID}.sra.lock"
        if [ -f "$LOCK_FILE" ]; then
            echo "Removing stale lock file: $LOCK_FILE"
            rm -f "$LOCK_FILE"
        fi
        # Also clean up any incomplete temp files
        if [ -f "$SRA_DIR/${SRR_ID}.sra.tmp" ] && [ ! -f "$SRA_DIR/${SRR_ID}.sra" ]; then
            echo "Removing incomplete download files..."
            rm -f "$SRA_DIR/${SRR_ID}.sra.tmp" "$SRA_DIR/${SRR_ID}.sra.prf"
        fi
    fi

    # Download with prefetch
    prefetch -p $SRR_ID
    # Convert to FASTQ using the correct SRA file path
    echo "Converting to FASTQ format..."
    SRA_PATH="$PWD/$SRR_ID/$SRR_ID.sra"
    if [ ! -f "$SRA_PATH" ]; then
        echo "Error: SRA file not found at $SRA_PATH"
        exit 1
    fi
    fastq-dump --split-files --gzip --outdir . "$SRA_PATH"
    # Rename files
    if [ -f "${SRR_ID}_1.fastq.gz" ]; then
        mv ${SRR_ID}_1.fastq.gz ${SAMPLE_NAME}_R1.fastq.gz
        mv ${SRR_ID}_2.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz
        echo "Downloaded: ${SAMPLE_NAME}_R1.fastq.gz and ${SAMPLE_NAME}_R2.fastq.gz"
    elif [ -f "${SRR_ID}.fastq.gz" ]; then
        mv ${SRR_ID}.fastq.gz ${SAMPLE_NAME}.fastq.gz
        echo "Downloaded: ${SAMPLE_NAME}.fastq.gz"
    fi
    # Run FastQC (check if reports already exist)
    echo "Running quality control..."
    if [ -f "${SAMPLE_NAME}_R1.fastq.gz" ] && [ -f "${SAMPLE_NAME}_R2.fastq.gz" ]; then
        if [ -f "${SAMPLE_NAME}_R1_fastqc.html" ] && [ -f "${SAMPLE_NAME}_R2_fastqc.html" ]; then
            echo "FastQC reports already exist. Skipping."
        else
            fastqc ${SAMPLE_NAME}_R1.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz
        fi
    elif [ -f "${SAMPLE_NAME}.fastq.gz" ]; then
        if [ -f "${SAMPLE_NAME}_fastqc.html" ]; then
            echo "FastQC report already exists. Skipping."
        else
            fastqc ${SAMPLE_NAME}.fastq.gz
        fi
    fi
    echo "Download complete for $SAMPLE_NAME!"
    exit 0
fi

# Split URLs (can be one or two for paired-end)
IFS=';' read -ra URLS <<< "$URLS_RAW"

for i in "${!URLS[@]}"; do
    URL="$PROTO://${URLS[$i]}"
    if [ ${#URLS[@]} -eq 2 ]; then
        # Paired-end
        SUFFIX=$((i+1))
        OUTFILE="${SAMPLE_NAME}_R${SUFFIX}.fastq.gz"
    else
        # Single-end
        OUTFILE="${SAMPLE_NAME}.fastq.gz"
    fi
    echo "Downloading $URL -> $OUTFILE ..."
    aria2c -c -x 4 -s 4 -o "$OUTFILE" "$URL"
done

# Run FastQC (check if reports already exist)
echo "Running quality control..."
if [ ${#URLS[@]} -eq 2 ]; then
    if [ -f "${SAMPLE_NAME}_R1_fastqc.html" ] && [ -f "${SAMPLE_NAME}_R2_fastqc.html" ]; then
        echo "FastQC reports already exist. Skipping."
    else
        fastqc ${SAMPLE_NAME}_R1.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz
    fi
else
    if [ -f "${SAMPLE_NAME}_fastqc.html" ]; then
        echo "FastQC report already exists. Skipping."
    else
        fastqc ${SAMPLE_NAME}.fastq.gz
    fi
fi

echo "Download complete for $SAMPLE_NAME!"

