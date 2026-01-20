# Repository Structure

## Complete File Tree

```
repo/
├── README.md                       # Main documentation
├── QUICKSTART.md                   # Quick start guide
├── PIPELINE_OVERVIEW.md            # Technical overview
├── VERSION.md                      # Version and changelog
├── REPOSITORY_STRUCTURE.md         # This file
├── .gitignore                      # Git ignore rules
│
├── setup_start.sh                  # Main setup script (installs all tools)
├── run_complete_pipeline.sh        # Main pipeline runner (runs entire workflow)
├── verify_setup.sh                 # Verification script (checks installation)
│
├── data/                           # Data acquisition and preprocessing
│   ├── download_reference.sh       # Download GRCh38 reference genome
│   ├── download_sample.sh          # Download single sample from NCBI
│   ├── download_ncbi_minimal.sh    # Download minimal test dataset
│   └── align_all.sh                # Align all FASTQ files to reference
│
└── variant_calling/                # Variant calling pipeline
    ├── mark_duplicates.sh          # Mark PCR and optical duplicates
    ├── add_read_groups.sh          # Add read group metadata
    ├── bqsr.sh                     # Base quality score recalibration
    ├── call_variants.sh            # Call variants with GATK Mutect2
    ├── filter_variants.sh          # Filter variants by quality
    ├── annotate_variants.sh        # Annotate with gene info
    ├── prioritize_variants.sh      # Select high-confidence variants
    ├── generate_reports.sh         # Generate summary reports
    ├── compare_samples.sh          # Identify tumor-specific variants
    ├── create_igv_session.sh       # Create IGV visualization
    └── run_pipeline.sh             # Run variant calling steps only
```

## File Categories

### Documentation (5 files)
- `README.md` - Complete documentation with all instructions
- `QUICKSTART.md` - Fast start guide for beginners
- `PIPELINE_OVERVIEW.md` - Technical details and concepts
- `VERSION.md` - Version history and compatibility
- `REPOSITORY_STRUCTURE.md` - This file

### Main Scripts (3 files)
- `setup_start.sh` - Install all tools and create environment
- `run_complete_pipeline.sh` - Run complete workflow start to finish
- `verify_setup.sh` - Verify installation is correct

### Data Scripts (4 files)
- `download_reference.sh` - Get reference genome (~30GB)
- `download_sample.sh` - Get one sample from NCBI SRA
- `download_ncbi_minimal.sh` - Get small test dataset
- `align_all.sh` - Align reads to reference with BWA

### Variant Calling Scripts (11 files)
- `mark_duplicates.sh` - Remove duplicates (Step 1)
- `add_read_groups.sh` - Add metadata (Step 2)
- `bqsr.sh` - Quality recalibration (Step 3)
- `call_variants.sh` - Variant calling (Step 4)
- `filter_variants.sh` - Quality filtering (Step 5)
- `annotate_variants.sh` - Gene annotation (Step 6)
- `prioritize_variants.sh` - Select best variants (Step 7)
- `generate_reports.sh` - Create reports (Step 8)
- `compare_samples.sh` - Find tumor variants (Step 9)
- `create_igv_session.sh` - Visualization setup (Step 10)
- `run_pipeline.sh` - Run steps 2-10 sequentially

## Script Relationships

### Dependency Graph

```
setup_start.sh
    └── Creates ~/somatic_variant_calling/ with all tools
         └── Creates setup_env.sh

download_reference.sh
    └── Downloads reference genome to data/reference/

download_ncbi_minimal.sh
    └── Downloads sample FASTQ files to data/fastq/

align_all.sh
    └── Requires: reference genome, FASTQ files
    └── Creates: BAM files in data/bam/

mark_duplicates.sh
    └── Requires: data/bam/*.bam
    └── Creates: data/bam_marked/*_marked.bam

add_read_groups.sh
    └── Requires: data/bam_marked/*_marked.bam
    └── Creates: data/bam_with_rg/*_rg.bam

bqsr.sh
    └── Requires: data/bam_with_rg/*_rg.bam, reference genome
    └── Creates: data/bam_recalibrated/*_recalibrated.bam

call_variants.sh
    └── Requires: data/bam_recalibrated/*_recalibrated.bam
    └── Creates: data/vcfs/*.vcf

filter_variants.sh
    └── Requires: data/vcfs/*.vcf
    └── Creates: data/vcfs_filtered/*.vcf

annotate_variants.sh
    └── Requires: data/vcfs_filtered/*.vcf
    └── Creates: data/vcfs_annotated/*.vcf

prioritize_variants.sh
    └── Requires: data/vcfs_annotated/*.vcf
    └── Creates: data/vcfs_prioritized/*.vcf

generate_reports.sh
    └── Requires: data/vcfs_filtered/*.vcf
    └── Creates: data/reports/*.txt

compare_samples.sh
    └── Requires: data/vcfs_filtered/*.vcf (tumor and normal)
    └── Creates: data/comparisons/*_tumor_specific.vcf

create_igv_session.sh
    └── Requires: BAM files, VCF files, reference
    └── Creates: data/igv_session.xml
```

### Execution Order

#### Full Pipeline (run_complete_pipeline.sh)
```
1. setup_start.sh              (if not --skip-setup)
2. download_reference.sh       (if not --skip-download)
3. download_ncbi_minimal.sh    (if not --skip-download)
4. align_all.sh                (if not --skip-alignment)
5. mark_duplicates.sh
6. add_read_groups.sh
7. bqsr.sh
8. call_variants.sh
9. filter_variants.sh
10. annotate_variants.sh
11. prioritize_variants.sh
12. generate_reports.sh
13. compare_samples.sh
14. create_igv_session.sh
```

#### Partial Pipeline (run_pipeline.sh)
```
Assumes BAM files exist in data/bam_marked/
1. add_read_groups.sh
2. bqsr.sh
3. call_variants.sh
4. filter_variants.sh
5. annotate_variants.sh
6. prioritize_variants.sh
7. generate_reports.sh
8. compare_samples.sh
9. create_igv_session.sh
```

## Expected Directory Structure After Setup

```
~/somatic_variant_calling/              # PROJECT_DIR
├── software/                           # All bioinformatics tools
│   ├── sratoolkit.3.0.10-ubuntu64/
│   ├── samtools-1.19.2/
│   ├── bcftools-1.19/
│   ├── bwa/
│   ├── FastQC/
│   ├── picard.jar
│   ├── gatk-4.5.0.0/
│   ├── VarScan.v2.4.6.jar
│   ├── strelka-2.9.10.centos6_x86_64/
│   └── IGV_Linux_2.17.4/
├── venv/                               # Python virtual environment
├── setup_env.sh                        # Environment activation script
├── logs/                               # Log files
└── results/                            # Symlinks to results

~/storage/variant_calling_data/         # DATA_DIR
├── reference/                          # Reference genome (~30GB)
│   ├── GRCh38_reference.fa
│   ├── GRCh38_reference.fa.fai
│   ├── GRCh38_reference.dict
│   └── All_20180418.vcf.gz            # Known variants for BQSR
├── fastq/                              # Raw sequencing reads (~50GB)
│   ├── patient1_normal_R1.fastq.gz
│   ├── patient1_normal_R2.fastq.gz
│   ├── patient1_tumor_R1.fastq.gz
│   └── patient1_tumor_R2.fastq.gz
├── bam/                                # Aligned reads (~30GB)
│   ├── patient1_normal.bam
│   ├── patient1_normal.bam.bai
│   ├── patient1_tumor.bam
│   └── patient1_tumor.bam.bai
├── bam_marked/                         # Deduplicated (~30GB)
│   ├── patient1_normal_marked.bam
│   └── patient1_tumor_marked.bam
├── bam_with_rg/                        # With read groups (~30GB)
│   ├── patient1_normal_rg.bam
│   └── patient1_tumor_rg.bam
├── bam_recalibrated/                   # BQSR processed (~30GB)
│   ├── patient1_normal_recalibrated.bam
│   └── patient1_tumor_recalibrated.bam
├── vcfs/                               # Raw variants (~1GB)
│   ├── patient1_normal.vcf
│   └── patient1_tumor.vcf
├── vcfs_filtered/                      # Filtered variants (~100MB)
│   ├── patient1_normal_filtered.vcf
│   └── patient1_tumor_filtered.vcf
├── vcfs_annotated/                     # Annotated variants (~200MB)
│   ├── patient1_normal_annotated.vcf
│   └── patient1_tumor_annotated.vcf
├── vcfs_prioritized/                   # High-confidence (~50MB)
│   ├── patient1_normal_prioritized.vcf
│   └── patient1_tumor_prioritized.vcf
├── comparisons/                        # Tumor-specific variants
│   └── patient1_tumor_specific.vcf
├── reports/                            # Summary statistics
│   ├── patient1_normal_stats.txt
│   └── patient1_tumor_stats.txt
├── metrics/                            # QC metrics
│   ├── patient1_normal_dup_metrics.txt
│   └── patient1_tumor_dup_metrics.txt
├── recal_tables/                       # BQSR tables
│   ├── patient1_normal_recal_data.table
│   └── patient1_tumor_recal_data.table
└── igv_session.xml                     # IGV visualization file
```

## Usage Patterns

### Pattern 1: Complete First-Time Setup
```bash
cd repo/
./run_complete_pipeline.sh
```
Use when: Starting from scratch with no tools or data

### Pattern 2: Re-run Analysis with New Data
```bash
cd repo/
./run_complete_pipeline.sh --skip-setup --skip-download
```
Use when: Tools installed, adding new samples

### Pattern 3: Re-run from Variant Calling
```bash
cd repo/
./run_complete_pipeline.sh --skip-setup --skip-download --skip-alignment
```
Use when: BAM files exist, need to re-call variants

### Pattern 4: Manual Step-by-Step
```bash
# Setup once
./setup_start.sh
source ~/somatic_variant_calling/setup_env.sh

# Download data once
cd ~/storage/variant_calling_data/reference
bash /path/to/repo/data/download_reference.sh

# Process each sample
cd ~/storage/variant_calling_data
bash /path/to/repo/variant_calling/mark_duplicates.sh
# ... etc
```
Use when: Need fine control over each step

## File Sizes

### Repository
- Scripts and documentation: ~140KB
- Portable, can be version controlled

### Runtime Data
- Reference genome: ~30GB
- Sample data (per sample): ~50GB
- Intermediate files: ~50GB per sample
- Final results: ~5GB per sample
- **Total for one tumor-normal pair: ~200GB**

## Customization Points

### Where to Customize

1. **Thread counts**: Edit `-t` parameters in `align_all.sh`
2. **Memory limits**: Edit `-Xmx` parameters in GATK calls
3. **Quality filters**: Edit thresholds in `filter_variants.sh`
4. **Sample names**: Edit patterns in individual scripts
5. **Paths**: Use `--data-dir` and `--project-dir` flags

### What NOT to Change

1. Directory names (`bam/`, `vcfs/`, etc.) - scripts expect these
2. File suffixes (`_marked.bam`, `_rg.bam`, etc.) - used for dependencies
3. Reference genome version - all scripts assume GRCh38
4. Tool versions - newer versions may have different parameters

## Maintenance

### Clean Up Space
```bash
# Remove intermediate BAM files after pipeline completes
cd ~/storage/variant_calling_data
rm -rf bam_marked/ bam_with_rg/

# Keep only final recalibrated BAM and VCF files
```

### Update Tools
```bash
# Re-run setup to update to latest tool versions
cd repo/
./setup_start.sh --project-dir ~/somatic_variant_calling
```

### Backup Important Files
```bash
# Always backup these:
~/storage/variant_calling_data/vcfs_prioritized/
~/storage/variant_calling_data/comparisons/
~/storage/variant_calling_data/reports/

# Can be regenerated:
~/storage/variant_calling_data/bam*/
~/storage/variant_calling_data/vcfs/
```
