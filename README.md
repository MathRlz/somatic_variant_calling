# Attention
This project used Artificial Intelligence assistance.

# Somatic Variant Calling Pipeline

Complete bioinformatics pipeline for detecting somatic mutations in cancer samples by comparing tumor and normal tissue.

## Overview

This repository contains all scripts and tools needed to perform somatic variant calling from raw sequencing data to annotated, high-confidence variants ready for clinical interpretation.

## Pipeline Workflow

The pipeline consists of 13 main steps:

1. **Environment Setup** - Install all required bioinformatics tools
2. **Reference Download** - Download human reference genome (GRCh38)
3. **Sample Download** - Download tumor and normal FASTQ files
4. **Read Alignment** - Align reads to reference genome using BWA
5. **Mark Duplicates** - Remove PCR and optical duplicates with Picard
6. **Add Read Groups** - Add sample metadata to BAM files
7. **BQSR** - Base quality score recalibration with GATK
8. **Variant Calling** - Call variants using GATK Mutect2
9. **Variant Filtering** - Filter variants by quality metrics
10. **Variant Annotation** - Annotate with gene names and functional impact
11. **Variant Prioritization** - Identify high-confidence variants
12. **Report Generation** - Generate summary statistics
13. **IGV Session** - Create visualization file for IGV
14. **Interpretation** - Interpret via [Cancer Genome Interpreter](https://www.cancergenomeinterpreter.orgCGI)

## Repository Structure

```
repo/
├── README.md                    # This file
├── setup_start.sh              # Environment setup script
├── run_complete_pipeline.sh    # Main pipeline runner
├── data/                       # Data download and alignment scripts
│   ├── download_reference.sh   # Download reference genome
│   ├── download_sample.sh      # Download single sample
│   ├── download_ncbi_minimal.sh # Download minimal test dataset
│   └── align_all.sh            # Align all FASTQ files
└── variant_calling/            # Variant calling pipeline scripts
    ├── mark_duplicates.sh      # Remove duplicate reads
    ├── add_read_groups.sh      # Add read group information
    ├── bqsr.sh                 # Base quality recalibration
    ├── call_variants.sh        # Call variants with Mutect2
    ├── filter_variants.sh      # Filter variants
    ├── annotate_variants.sh    # Annotate variants
    ├── prioritize_variants.sh  # Select high-confidence variants (configurable MIN_DP)
    ├── generate_reports.sh     # Generate statistics
    └── create_igv_session.sh   # Create IGV visualization
```

## Quick Start

### Option 1: Run Complete Pipeline (Recommended for first time)

Run the entire pipeline from scratch, including environment setup, data download, and analysis:

```bash
cd repo
./run_complete_pipeline.sh
```

This will:
- Install all bioinformatics tools in `~/somatic_variant_calling`
- Download reference genome and sample data to `~/storage/variant_calling_data`
- Run the complete analysis pipeline
- Generate results and reports

### Option 2: Run with Custom Directories

```bash
./run_complete_pipeline.sh --data-dir /path/to/data --project-dir /path/to/tools
```

### Option 3: Skip Certain Steps

If you already have tools installed or data downloaded:

```bash
# Skip environment setup (assumes tools already installed)
./run_complete_pipeline.sh --skip-setup

# Skip data download (use existing data)
./run_complete_pipeline.sh --skip-download

# Skip alignment (use existing BAM files)
./run_complete_pipeline.sh --skip-download --skip-alignment

# Run only from variant calling onwards
./run_complete_pipeline.sh --skip-setup --skip-download --skip-alignment
```

## System Requirements

### Hardware
- **CPU**: 4+ cores recommended (8+ for faster processing)
- **RAM**: 16GB minimum, 32GB recommended
- **Storage**: 200GB minimum free space
  - Reference genome: ~30GB
  - Sample data: ~50-100GB depending on coverage
  - Intermediate files: ~50GB
  - Tools and software: ~5GB

### Software
- **Operating System**: Ubuntu 25.04 (or similar Linux distribution)
- **Java**: Version 8 or higher
- **Python**: Version 3.8 or higher
- **Internet**: Required for downloading tools and data

## Installed Tools

The setup script installs the following bioinformatics tools:

- **SRA Toolkit** (3.0.10) - Download data from NCBI
- **SAMtools** (1.19.2) - BAM file manipulation
- **BCFtools** (1.19) - VCF file processing
- **BWA** (0.7.17) - Read alignment
- **FastQC** (0.12.1) - Quality control
- **Picard** (3.1.1) - BAM processing
- **GATK4** (4.5.0.0) - Variant calling and processing
- **VarScan2** (2.4.6) - Alternative variant caller
- **Strelka2** (2.9.10) - Somatic variant caller
- **IGV** (2.17.4) - Genome visualization
- **MultiQC** (optional) - Report aggregation

## Directory Structure After Setup

```
~/somatic_variant_calling/        # Project directory (tools)
├── software/                     # All bioinformatics tools
├── venv/                        # Python virtual environment
├── setup_env.sh                 # Environment activation script
├── results/                     # Placeholder for results
└── logs/                        # Placeholder for logs

~/storage/variant_calling_data/   # Data directory
├── reference/                   # Reference genome files
│   ├── GRCh38_reference.fa
│   └── All_20180418.vcf.gz      # Known variants for BQSR
├── fastq/                       # Raw sequencing data
├── bam/                         # Aligned reads
├── bam_marked/                  # Deduplicated BAM files
├── bam_with_rg/                 # BAM files with read groups
├── bam_recalibrated/            # BQSR-processed BAM files
├── vcfs/                        # Raw variant calls
├── vcfs_filtered/               # Filtered variants
├── vcfs_annotated/              # Annotated variants
├── vcfs_prioritized/            # High-confidence variants
├── reports/                     # Summary statistics
├── metrics/                     # QC metrics
├── recal_tables/                # BQSR recalibration tables
└── igv_session.xml              # IGV visualization file
```

## Manual Usage

### 1. Setup Environment Only

```bash
./setup_start.sh
```

Or with custom paths:

```bash
./setup_start.sh --data-dir ~/my_data --project-dir ~/my_tools
```

### 2. Activate Environment

Before running any analysis, activate the environment:

```bash
source ~/somatic_variant_calling/setup_env.sh
```

### 3. Download Data

```bash
cd ~/storage/variant_calling_data/reference
bash /path/to/repo/data/download_reference.sh

cd ~/storage/variant_calling_data/fastq
bash /path/to/repo/data/download_ncbi_minimal.sh
# Or download individual samples:
bash /path/to/repo/data/download_sample.sh SRR2057563 tumor_sample
```

### 4. Run Pipeline Steps Individually

```bash
cd ~/storage/variant_calling_data

# Align reads
cd fastq
bash /path/to/repo/data/align_all.sh
mv *.bam *.bam.bai ../bam/

# Process BAM files
cd ~/storage/variant_calling_data
bash /path/to/repo/variant_calling/mark_duplicates.sh
bash /path/to/repo/variant_calling/add_read_groups.sh
bash /path/to/repo/variant_calling/bqsr.sh

# Call and process variants
bash /path/to/repo/variant_calling/call_variants.sh
bash /path/to/repo/variant_calling/filter_variants.sh
bash /path/to/repo/variant_calling/annotate_variants.sh
bash /path/to/repo/variant_calling/prioritize_variants.sh

# Generate reports and visualization
bash /path/to/repo/variant_calling/generate_reports.sh
bash /path/to/repo/variant_calling/create_igv_session.sh
```

## Output Files

### Key Results

1. **High-confidence variants**: `vcfs_prioritized/`
   - Somatic mutations with high confidence
   - Ready for clinical interpretation

2. **Annotated variants**: `vcfs_annotated/`
   - All variants with gene names and functional annotations
   - Includes impact predictions

3. **Reports**: `reports/`
   - Summary statistics for each sample
   - Variant counts and quality metrics

### Visualization

Open results in IGV (Integrative Genomics Viewer):

```bash
cd ~/somatic_variant_calling/software/IGV_Linux_2.17.4
./igv.sh
```

Then: File > Open Session > Select `~/storage/variant_calling_data/igv_session.xml`

## Customization

### Using Your Own Data

1. Place paired-end FASTQ files in the `fastq/` directory:
   ```
   sample1_R1.fastq.gz
   sample1_R2.fastq.gz
   sample2_R1.fastq.gz
   sample2_R2.fastq.gz
   ```

2. Run pipeline starting from alignment:
   ```bash
   ./run_complete_pipeline.sh --skip-setup --skip-download
   ```

### Modifying Pipeline Steps

All individual scripts are in the `variant_calling/` directory. You can:
- Edit scripts to change parameters
- Run individual steps manually
- Add custom filtering or annotation steps

## Troubleshooting

### Disk Space Issues

If you get disk space errors:

```bash
# Check available space
df -h ~/storage/variant_calling_data

# Clean up intermediate files
cd ~/storage/variant_calling_data
rm -rf bam_marked/*.bam  # Keep only final BAM files
rm -rf fastq/*.fastq.gz  # After alignment is complete
```

### Memory Issues

If tools crash due to memory:
- Reduce number of threads (edit `-t` parameter in alignment scripts)
- Process samples one at a time instead of in parallel
- Use smaller test dataset first

### Tool Not Found Errors

Make sure environment is activated:

```bash
source ~/somatic_variant_calling/setup_env.sh
```

Or check if specific tool is installed:

```bash
which samtools
which gatk
```

## Expected Runtime

Approximate times for 30x coverage whole genome sequencing:

- Setup: 30-60 minutes (one time only)
- Reference download: 30-60 minutes
- Sample download: 1-4 hours per sample
- Alignment: 4-8 hours per sample
- Variant calling: 2-4 hours
- Annotation and filtering: 30 minutes
- Total: 12-24 hours for complete pipeline

For targeted sequencing (exomes) or lower coverage, times will be proportionally shorter.

## Citation

If you use this pipeline, please cite the relevant tools:

- **BWA**: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.
- **GATK**: McKenna A, et al. (2010) The Genome Analysis Toolkit. Genome Research, 20:1297-303.
- **SAMtools**: Li H., et al. (2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25:2078-9.

## Support

For issues with:
- **Pipeline scripts**: Check individual script comments and error messages
- **Bioinformatics tools**: Refer to tool documentation (links in DETAILED_EXPLANATION.md.pdf)
- **Data download**: Check NCBI/ENA status and network connectivity

## License

This pipeline is provided for educational and research purposes. Individual tools have their own licenses.
