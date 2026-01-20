#!/bin/bash

# Usage: ./download_sample.sh SRR_ID [SAMPLE_NAME]

if [ $# -lt 1 ]; then
    echo "Usage: $0 SRR_ID [SAMPLE_NAME]"
    echo "Example: $0 SRR2057563 tumor_sample"
    exit 1
fi

SRR_ID=$1
SAMPLE_NAME=${2:-$SRR_ID}

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
    # Run FastQC
    echo "Running quality control..."
    if [ -f "${SAMPLE_NAME}_R1.fastq.gz" ] && [ -f "${SAMPLE_NAME}_R2.fastq.gz" ]; then
        fastqc ${SAMPLE_NAME}_R1.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz
    elif [ -f "${SAMPLE_NAME}.fastq.gz" ]; then
        fastqc ${SAMPLE_NAME}.fastq.gz
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
    aria2c -x 4 -s 4 -o "$OUTFILE" "$URL"
done

# Run FastQC
echo "Running quality control..."
if [ ${#URLS[@]} -eq 2 ]; then
    fastqc ${SAMPLE_NAME}_R1.fastq.gz ${SAMPLE_NAME}_R2.fastq.gz
else
    fastqc ${SAMPLE_NAME}.fastq.gz
fi

echo "Download complete for $SAMPLE_NAME!"

