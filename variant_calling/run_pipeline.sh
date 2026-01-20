#!/bin/bash

# Master pipeline script - runs all steps sequentially

set -e  # Exit on error

echo "=== Starting Variant Calling Pipeline ==="

echo "Step 1: Adding read groups..."
bash add_read_groups.sh

echo "Step 2: Base quality score recalibration..."
bash bqsr.sh

echo "Step 3: Calling variants..."
bash call_variants.sh

echo "Step 4: Filtering variants..."
bash filter_variants.sh

echo "Step 5: Annotating variants..."
bash annotate_variants.sh

echo "Step 6: Prioritizing variants..."
bash prioritize_variants.sh

echo "Step 7: Generating reports..."
bash generate_reports.sh

echo "Step 8: Comparing samples..."
bash compare_samples.sh

echo "Step 9: Creating IGV session..."
bash create_igv_session.sh

echo "=== Pipeline completed successfully! ==="
echo "Results available in:"
echo "  - vcfs_filtered/       : Filtered variants"
echo "  - vcfs_annotated/      : Annotated variants"
echo "  - vcfs_prioritized/    : High-confidence variants"
echo "  - reports/             : Summary statistics"
echo "  - comparisons/         : Tumor-specific variants"
echo "  - igv_session.xml      : IGV visualization file"
