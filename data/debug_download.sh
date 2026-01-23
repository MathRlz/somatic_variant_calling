t d#!/bin/bash

# Debug script to check why downloads are slow
# Usage: ./debug_download.sh SRR_ID

if [ $# -lt 1 ]; then
    echo "Usage: $0 SRR_ID"
    echo "Example: $0 SRR2089355"
    exit 1
fi

SRR_ID=$1

echo "========================================="
echo "Download Debug for $SRR_ID"
echo "========================================="
echo ""

# Check 1: ENA availability
echo "1. Checking ENA availability..."
echo "================================"
ENA_API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR_ID&result=read_run&fields=fastq_https,fastq_ftp,fastq_bytes&format=json"
echo "URL: $ENA_API_URL"
echo ""
FASTQ_INFO=$(curl -s "$ENA_API_URL")
echo "Raw response:"
echo "$FASTQ_INFO" | jq '.' 2>/dev/null || echo "$FASTQ_INFO"
echo ""

FASTQ_HTTPS=$(echo "$FASTQ_INFO" | grep -oP '"fastq_https":"\K[^"]*')
FASTQ_FTP=$(echo "$FASTQ_INFO" | grep -oP '"fastq_ftp":"\K[^"]*')
FASTQ_BYTES=$(echo "$FASTQ_INFO" | grep -oP '"fastq_bytes":"\K[^"]*')

if [ -n "$FASTQ_HTTPS" ]; then
    echo "✓ HTTPS URLs found:"
    echo "$FASTQ_HTTPS" | tr ';' '\n' | sed 's/^/  - https:\/\//'
    echo ""
    if [ -n "$FASTQ_BYTES" ]; then
        TOTAL_BYTES=$(echo "$FASTQ_BYTES" | tr ';' '+' | bc)
        TOTAL_GB=$(echo "scale=2; $TOTAL_BYTES / 1073741824" | bc)
        echo "  Total size: $TOTAL_GB GB"
    fi
elif [ -n "$FASTQ_FTP" ]; then
    echo "⚠ Only FTP URLs found (HTTPS preferred):"
    echo "$FASTQ_FTP" | tr ';' '\n' | sed 's/^/  - ftp:\/\//'
    echo ""
    if [ -n "$FASTQ_BYTES" ]; then
        TOTAL_BYTES=$(echo "$FASTQ_BYTES" | tr ';' '+' | bc)
        TOTAL_GB=$(echo "scale=2; $TOTAL_BYTES / 1073741824" | bc)
        echo "  Total size: $TOTAL_GB GB"
    fi
else
    echo "✗ No ENA URLs found - will use NCBI prefetch (slower)"
fi
echo ""

# Check 2: NCBI metadata
echo "2. Checking NCBI metadata..."
echo "============================="
vdb-dump --info "$SRR_ID" 2>&1 | head -20
echo ""

# Check 3: Available download tools
echo "3. Checking available download tools..."
echo "========================================"
echo -n "aria2c: "
command -v aria2c &> /dev/null && echo "✓ installed" || echo "✗ not found"

echo -n "prefetch: "
command -v prefetch &> /dev/null && echo "✓ installed" || echo "✗ not found"

echo -n "fasterq-dump: "
command -v fasterq-dump &> /dev/null && echo "✓ installed" || echo "✗ not found"

echo -n "ascp (Aspera): "
if command -v ascp &> /dev/null; then
    echo "✓ installed (2-3x faster NCBI downloads!)"
else
    echo "✗ not found - install for faster NCBI downloads"
    echo "   Download from: https://www.ibm.com/aspera/connect/"
fi

echo -n "pigz: "
command -v pigz &> /dev/null && echo "✓ installed" || echo "✗ not found"

echo ""

# Check 4: Current vdb-config
echo "4. Checking SRA Toolkit configuration..."
echo "=========================================="
CONFIG_FILE="$HOME/.ncbi/user-settings.mkfg"
if [ -f "$CONFIG_FILE" ]; then
    echo "Configuration file: $CONFIG_FILE"
    echo ""
    echo "SRA Lite preference:"
    grep -A 2 "sraLite" "$CONFIG_FILE" 2>/dev/null || echo "  Not configured (using default)"
    echo ""
    echo "Prefetch settings:"
    grep -A 2 "prefetch" "$CONFIG_FILE" 2>/dev/null || echo "  Using defaults"
else
    echo "⚠ No configuration file found at $CONFIG_FILE"
    echo "  Run: vdb-config --interactive"
fi
echo ""

# Check 5: Network speed test
echo "5. Testing network speed to ENA..."
echo "===================================="
if [ -n "$FASTQ_HTTPS" ]; then
    FIRST_URL=$(echo "$FASTQ_HTTPS" | cut -d';' -f1)
    FULL_URL="https://$FIRST_URL"
    echo "Testing download speed from: $FULL_URL"
    echo "(downloading first 10MB to test)..."
    time curl -r 0-10485760 -o /tmp/test_download.tmp "$FULL_URL" 2>&1 | grep -E '(Downloaded|speed)'
    rm -f /tmp/test_download.tmp
else
    echo "No ENA URLs to test"
fi
echo ""

# Check 6: Disk space
echo "6. Checking disk space..."
echo "=========================="
df -h . | awk 'NR==1 || NR==2'
echo ""
echo "Note: fasterq-dump needs ~17x the SRA file size in temporary space"
echo ""

# Recommendations
echo "========================================="
echo "RECOMMENDATIONS:"
echo "========================================="
if [ -z "$FASTQ_HTTPS" ] && [ -z "$FASTQ_FTP" ]; then
    echo "⚠ ENA URLs not available for this accession"
    echo "  → You'll need to use NCBI prefetch (slower)"
    echo "  → Consider installing Aspera for 2-3x speedup:"
    echo "    https://www.ibm.com/aspera/connect/"
    echo ""
    echo "  → Alternative: Enable SRA Lite format for smaller files:"
    echo "    vdb-config --set /repository/user/main/public/apps/sra/volumes/sraLite/policy=always"
else
    echo "✓ ENA URLs available - should use fast aria2c download"
    echo "  If falling back to prefetch, check ENA API response above"
fi

if ! command -v pigz &> /dev/null; then
    echo ""
    echo "⚠ Install pigz for faster compression:"
    echo "  sudo apt install pigz"
fi

if ! command -v ascp &> /dev/null; then
    echo ""
    echo "⚠ Install Aspera Connect for faster NCBI downloads:"
    echo "  Download from: https://www.ibm.com/aspera/connect/"
    echo "  This provides 2-3x speedup when ENA is unavailable"
fi

echo ""
