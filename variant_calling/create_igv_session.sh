#!/bin/bash

# Create IGV session file for visualization

PRIORITY_DIR="vcfs_prioritized"
BAM_DIR="bam_recalibrated"
REFERENCE="reference/GRCh38_reference.fa"

cat > igv_session.xml << EOF
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="hg38" version="8">
    <Resources>
        <Resource path="$REFERENCE"/>
EOF

for sample in TCR002101 TCR002182 TCR002361; do
    echo "        <Resource path=\"$BAM_DIR/${sample}-T_recalibrated.bam\"/>" >> igv_session.xml
    echo "        <Resource path=\"$BAM_DIR/${sample}-N_recalibrated.bam\"/>" >> igv_session.xml
    echo "        <Resource path=\"$PRIORITY_DIR/${sample}_high_confidence.vcf.gz\"/>" >> igv_session.xml
done

cat >> igv_session.xml << EOF
    </Resources>
</Session>
EOF

echo "IGV session file created: igv_session.xml"
echo "Open this file in IGV to visualize your results"
