# Somatic Variant Calling Pipeline

This document provides a complete technical overview of the somatic variant calling pipeline, including the biological rationale and exact commands for each step.

## Table of Contents

1. [Introduction](#introduction)
2. [Pipeline Summary](#pipeline-summary)
3. [Prerequisites](#prerequisites)
4. [Step-by-Step Workflow](#step-by-step-workflow)
   - [Step 0: Environment Setup](#step-0-environment-setup)
   - [Step 1: Reference Genome Download](#step-1-reference-genome-download)
   - [Step 2: Sample Data Download](#step-2-sample-data-download)
   - [Step 3: Read Alignment](#step-3-read-alignment)
   - [Step 4: Mark Duplicates](#step-4-mark-duplicates)
   - [Step 5: Add Read Groups](#step-5-add-read-groups)
   - [Step 6: Base Quality Score Recalibration](#step-6-base-quality-score-recalibration-bqsr)
   - [Step 7: Somatic Variant Calling](#step-7-somatic-variant-calling)
   - [Step 8: Variant Filtering](#step-8-variant-filtering)
   - [Step 9: Variant Annotation](#step-9-variant-annotation)
   - [Step 10: Variant Prioritization](#step-10-variant-prioritization)
   - [Step 11: Report Generation](#step-11-report-generation)
   - [Step 12: IGV Session Creation](#step-12-igv-session-creation)
5. [Data Flow Diagram](#data-flow-diagram)
6. [Directory Structure](#directory-structure)
7. [Quality Control](#quality-control)
8. [Troubleshooting](#troubleshooting)
9. [References](#references)

---

## Introduction

### What is Somatic Variant Calling?

Somatic variant calling identifies genetic mutations that occur in cancer cells but are not present in normal tissue. These mutations are acquired during a person's lifetime and can drive cancer development.

### Why Compare Tumor and Normal?

By sequencing both tumor tissue and matched normal tissue (usually blood) from the same patient, we can:

1. **Identify true somatic mutations** - Variants unique to the tumor
2. **Filter germline variants** - Inherited variants present in both samples
3. **Reduce false positives** - Technical artifacts that appear in both samples
4. **Enable precision medicine** - Identify targetable mutations for treatment

### Variant Types

| Type | Description | Example |
|------|-------------|---------|
| **SNV** | Single Nucleotide Variant | A→G at position 12345 |
| **Indel** | Insertion or deletion | +ATG or -CCC |
| **Somatic** | Present only in tumor | Cancer driver mutations |
| **Germline** | Present in normal tissue | Inherited variants |

---

## Pipeline Summary

```
FASTQ (raw reads)
    ↓ [BWA alignment]
BAM (aligned reads)
    ↓ [Mark duplicates]
BAM (deduplicated)
    ↓ [Add read groups]
BAM (with metadata)
    ↓ [BQSR]
BAM (quality recalibrated)
    ↓ [GATK Mutect2]
VCF (raw variants)
    ↓ [Filtering]
VCF (filtered variants)
    ↓ [Annotation]
VCF (annotated variants)
    ↓ [Prioritization]
VCF (high-confidence somatic mutations)
```

---

## Prerequisites

### Required Software

| Tool | Version | Purpose |
|------|---------|---------|
| BWA | 0.7.19+ | Read alignment |
| SAMtools | 1.19+ | BAM manipulation |
| BCFtools | 1.19+ | VCF manipulation |
| GATK | 4.5.0+ | Variant calling (requires Java 17+) |
| VEP | Latest | Variant annotation |
| FastQC | 0.12+ | Quality control |
| aria2c | - | Parallel downloads |

### Hardware Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 4 cores | 8+ cores |
| RAM | 16 GB | 32 GB |
| Storage | 200 GB | 500 GB |

### Reference Data

- GRCh38 human reference genome (UCSC naming: chr1, chr2, ...)
- dbSNP database (remapped to UCSC contig names)

---

## Step-by-Step Workflow

### Step 0: Environment Setup

**Script:** `setup_start.sh`

**Purpose:** Install all bioinformatics tools and configure the environment.

#### System Dependencies
```bash
sudo apt install -y \
    build-essential cmake git wget curl unzip \
    default-jdk python3 python3-pip python3-venv \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev \
    libcurl4-openssl-dev libssl-dev autoconf automake \
    parallel aria2
```

#### Key Tools Installed

| Tool | Installation Method |
|------|---------------------|
| SAMtools 1.19.2 | Compiled from source |
| BCFtools 1.19 | Compiled from source |
| BWA | Compiled from git |
| GATK 4.5.0.0 | Downloaded ZIP |
| VEP | Via Miniconda/Bioconda |

#### Environment Activation
```bash
source $PROJECT_DIR/setup_env.sh
```

---

### Step 1: Reference Genome Download

**Script:** `data/download_reference.sh`

**Purpose:** Download and index the GRCh38 human reference genome and dbSNP database.

#### 1.1 Download Reference FASTA

```bash
aria2c -c -x 16 -s 16 -k 1M -o GRCh38_reference.fa.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GRCh38_reference.fa.gz
```

**Why this reference?** The `ucsc_ids` version uses `chr1, chr2, ...` naming convention, which is compatible with most downstream tools and databases.

#### 1.2 Create Indices

```bash
# BWA index (for alignment)
bwa index GRCh38_reference.fa

# SAMtools index (for random access)
samtools faidx GRCh38_reference.fa

# GATK sequence dictionary
java -jar picard.jar CreateSequenceDictionary \
    R=GRCh38_reference.fa \
    O=GRCh38_reference.dict
```

**Output files:**
- `.bwt`, `.pac`, `.ann`, `.amb`, `.sa` (BWA indices)
- `.fai` (SAMtools index)
- `.dict` (GATK dictionary)

#### 1.3 Download and Remap dbSNP

```bash
# Download from NCBI (uses RefSeq naming: NC_000001.11)
aria2c -c -x 16 -s 16 -k 1M -o dbsnp_ncbi.vcf.gz \
    "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"

# Remap contig names to UCSC style (NC_000001.11 → chr1)
bcftools annotate --rename-chrs contig_map.txt dbsnp_ncbi.vcf.gz \
    -Oz -o dbsnp_156.grch38.vcf.gz

# Index
bcftools index -t dbsnp_156.grch38.vcf.gz
```

**Why remap?** GATK requires contig names in the VCF to match the reference exactly. NCBI uses RefSeq accessions (NC_000001.11) while our reference uses UCSC names (chr1).

---

### Step 2: Sample Data Download

**Script:** `data/download_sample.sh`

**Purpose:** Download paired-end FASTQ files from NCBI SRA or ENA.

#### Download Methods

```bash
# Preferred: Direct HTTPS from ENA (faster)
aria2c -c -x 4 -s 4 "${FASTQ_URL_R1}" -o "${SAMPLE}_R1.fastq.gz"
aria2c -c -x 4 -s 4 "${FASTQ_URL_R2}" -o "${SAMPLE}_R2.fastq.gz"

# Fallback: NCBI SRA Toolkit
prefetch ${SRR_ID}
fastq-dump --split-files --gzip ${SRR_ID}
```

#### Quality Control

```bash
fastqc ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
```

**Output:**
- `{PATIENT}-T_R1.fastq.gz`, `{PATIENT}-T_R2.fastq.gz` (tumor)
- `{PATIENT}-N_R1.fastq.gz`, `{PATIENT}-N_R2.fastq.gz` (normal)
- FastQC HTML reports

---

### Step 3: Read Alignment

**Script:** `data/align_sample.sh`

**Purpose:** Align paired-end reads to the reference genome using BWA-MEM.

#### Command

```bash
bwa mem -t ${NUM_PROCESSORS} \
    ${REFERENCE} \
    ${SAMPLE}_R1.fastq.gz \
    ${SAMPLE}_R2.fastq.gz \
    | samtools sort -@ ${NUM_PROCESSORS} -o ${SAMPLE}.bam -

samtools index ${SAMPLE}.bam
```

#### How BWA-MEM Works

1. **Seeding:** Finds exact matches (seeds) between reads and reference
2. **Extension:** Extends seeds using Smith-Waterman alignment
3. **Optimal for:** Reads >70bp, handles indels and split reads well

**Output:** `bam/${SAMPLE}.bam` (coordinate-sorted with index)

---

### Step 4: Mark Duplicates

**Script:** `variant_calling/mark_duplicates.sh`

**Purpose:** Identify and flag PCR/optical duplicate reads.

#### Command

```bash
gatk MarkDuplicates \
    -I ${SAMPLE}.bam \
    -O ${SAMPLE}_marked.bam \
    -M ${SAMPLE}_dup_metrics.txt \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT
```

#### Why Mark Duplicates?

- PCR amplification creates duplicate molecules from the same original fragment
- Duplicates inflate variant allele frequencies and can cause false positives
- Reads are marked (not removed) so they're excluded from variant calling but retained for coverage metrics

**Output:** `bam_marked/${SAMPLE}_marked.bam`

---

### Step 5: Add Read Groups

**Script:** `variant_calling/add_read_groups.sh`

**Purpose:** Add read group information required by GATK.

#### Command

```bash
gatk AddOrReplaceReadGroups \
    -I ${SAMPLE}_marked.bam \
    -O ${SAMPLE}_rg.bam \
    -RGID "${SAMPLE}" \
    -RGLB "lib1" \
    -RGPL "ILLUMINA" \
    -RGPU "unit1" \
    -RGSM "${SAMPLE}"
```

#### Read Group Tags

| Tag | Description | Purpose |
|-----|-------------|---------|
| RGID | Read Group ID | Unique identifier |
| RGLB | Library | Track library prep batches |
| RGPL | Platform | Apply platform-specific error models |
| RGPU | Platform Unit | Sequencing run identifier |
| RGSM | Sample Name | Distinguish samples in multi-sample BAMs |

**Output:** `bam_with_rg/${SAMPLE}_rg.bam`

---

### Step 6: Base Quality Score Recalibration (BQSR)

**Script:** `variant_calling/bqsr.sh`

**Purpose:** Correct systematic errors in base quality scores assigned by the sequencer.

#### Step 6.1: Build Recalibration Table

```bash
gatk --java-options "-Xmx4g" BaseRecalibrator \
    -I ${SAMPLE}_rg.bam \
    -R ${REFERENCE} \
    --known-sites ${DBSNP_VCF} \
    -O ${SAMPLE}_recal_data.table
```

**How it works:**
1. Identifies bases at known variant sites (from dbSNP)
2. Analyzes quality scores by: read group, quality bin, cycle, dinucleotide context
3. Builds empirical error model

#### Step 6.2: Apply Recalibration

```bash
gatk --java-options "-Xmx4g" ApplyBQSR \
    -I ${SAMPLE}_rg.bam \
    -R ${REFERENCE} \
    --bqsr-recal-file ${SAMPLE}_recal_data.table \
    -O ${SAMPLE}_recalibrated.bam
```

#### Why BQSR?

- Illumina quality scores are often overconfident
- Systematic biases exist (e.g., certain cycles, contexts)
- Accurate quality scores improve variant calling specificity

**Validation:** Output BAM must be >1MB (corrupted files will be much smaller)

**Output:** `bam_recalibrated/${SAMPLE}_recalibrated.bam`

---

### Step 7: Somatic Variant Calling

**Script:** `variant_calling/call_variants.sh`

**Purpose:** Identify somatic mutations using GATK Mutect2.

#### Command

```bash
gatk --java-options "-Xmx4g" Mutect2 \
    -R ${REFERENCE} \
    -I ${TUMOR}_recalibrated.bam \
    -I ${NORMAL}_recalibrated.bam \
    -normal "${NORMAL_SAMPLE_NAME}" \
    -O ${PATIENT}_raw.vcf.gz \
    --f1r2-tar-gz ${PATIENT}_f1r2.tar.gz \
    --native-pair-hmm-threads ${NUM_PROCESSORS}
```

#### Mutect2 Algorithm

1. **Local Assembly:** Builds haplotypes from reads in active regions
2. **Pair-HMM:** Calculates likelihood of each read given each haplotype
3. **Somatic Model:** Compares tumor vs normal to identify tumor-specific variants
4. **Initial Filters:** Applies preliminary quality filters

#### Key Parameters

| Parameter | Purpose |
|-----------|---------|
| `-I` (×2) | Input both tumor and normal BAMs |
| `-normal` | Designates normal sample |
| `--f1r2-tar-gz` | Outputs orientation data for artifact detection |
| `--native-pair-hmm-threads` | Parallelization |

**Output:**
- `vcfs/${PATIENT}_raw.vcf.gz` - Raw variant calls
- `vcfs/${PATIENT}_f1r2.tar.gz` - F1R2 orientation data

---

### Step 8: Variant Filtering

**Script:** `variant_calling/filter_variants.sh`

**Purpose:** Apply quality filters to remove artifacts and low-confidence calls.

#### Step 8.1: Learn Read Orientation Model

```bash
gatk LearnReadOrientationModel \
    -I ${PATIENT}_f1r2.tar.gz \
    -O ${PATIENT}_read-orientation-model.tar.gz
```

**Purpose:** Detects orientation bias artifacts (OxoG damage, FFPE artifacts)

#### Step 8.2: Apply Filters

```bash
gatk FilterMutectCalls \
    -R ${REFERENCE} \
    -V ${PATIENT}_raw.vcf.gz \
    --ob-priors ${PATIENT}_read-orientation-model.tar.gz \
    -O ${PATIENT}_filtered.vcf.gz
```

#### Filters Applied

| Filter | Description |
|--------|-------------|
| `weak_evidence` | Insufficient support for variant |
| `strand_bias` | Variant only seen on one strand |
| `slippage` | Likely polymerase slippage in repeats |
| `orientation` | Orientation bias artifact |
| `contamination` | Possible sample contamination |
| `germline` | Likely germline variant |

#### Step 8.3: Normalize VCF

```bash
bcftools norm -m -both -f ${REFERENCE} ${PATIENT}_filtered.vcf.gz \
    -Oz -o ${PATIENT}_filtered_norm.vcf.gz

tabix -p vcf ${PATIENT}_filtered.vcf.gz
```

**Purpose:**
- `-m -both`: Splits multiallelic sites into biallelic records
- `-f`: Left-aligns indels and validates reference alleles

**Output:** `vcfs_filtered/${PATIENT}_filtered.vcf.gz`

---

### Step 9: Variant Annotation

**Script:** `variant_calling/annotate_variants.sh`

**Purpose:** Add functional annotations using Ensembl VEP.

#### Command

```bash
vep \
    --input_file ${PATIENT}_filtered.vcf.gz \
    --output_file ${PATIENT}_annotated.vcf \
    --format vcf \
    --vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache \
    --offline \
    --everything \
    --fork 4

bgzip ${PATIENT}_annotated.vcf
tabix -p vcf ${PATIENT}_annotated.vcf.gz
```

#### Annotations Added

| Annotation | Description |
|------------|-------------|
| Gene/Transcript | Affected gene and transcript IDs |
| Consequence | Effect type (missense, frameshift, synonymous, etc.) |
| SIFT | Functional impact prediction (tolerated/deleterious) |
| PolyPhen | Protein structure impact (benign/damaging) |
| gnomAD AF | Population allele frequency |
| ClinVar | Clinical significance |

**Output:** `vcfs_annotated/${PATIENT}_annotated.vcf.gz`

---

### Step 10: Variant Prioritization

**Script:** `variant_calling/prioritize_variants.sh`

**Purpose:** Select high-confidence variants for review.

#### Command

```bash
# Default: MIN_DP=20, MIN_VAF=0.05
bcftools view -i '
  (FILTER="PASS" ||
   FILTER~"weak_evidence" ||
   FILTER~"strand_bias" ||
   FILTER~"clustered_events") &&
  FORMAT/DP >= $MIN_DP &&
  FORMAT/AF >= $MIN_VAF
' ${PATIENT}_filtered.vcf.gz -Oz -o ${PATIENT}_high_confidence.vcf.gz

tabix -p vcf ${PATIENT}_high_confidence.vcf.gz
```

#### Filtering Criteria

| Criterion | Default | Description |
|-----------|---------|-------------|
| Filter status | PASS + soft filters | Rescues potential drivers with soft filter flags |
| Read depth | ≥20 (MIN_DP) | Configurable via environment variable |
| VAF | ≥5% (MIN_VAF) | Removes noise while keeping subclonal drivers |

#### Custom Thresholds

```bash
# Less strict filtering (e.g., to include more variants)
MIN_DP=10 ./prioritize_variants.sh PATIENT lenient
```

**Output:** `vcfs_prioritized/${PATIENT}_high_confidence.vcf.gz`

---

### Step 11: Report Generation

**Script:** `variant_calling/generate_reports.sh`

**Purpose:** Generate summary statistics.

#### Commands

```bash
# Comprehensive statistics
bcftools stats ${PATIENT}_high_confidence.vcf.gz > ${PATIENT}_stats.txt

# Variant counts by type
bcftools view -H ${PATIENT}_high_confidence.vcf.gz | wc -l > ${PATIENT}_total_variants.txt
bcftools view -v snps -H ${PATIENT}_high_confidence.vcf.gz | wc -l > ${PATIENT}_snvs.txt
bcftools view -v indels -H ${PATIENT}_high_confidence.vcf.gz | wc -l > ${PATIENT}_indels.txt
```

**Output:**
- `reports/${PATIENT}_stats.txt` - Detailed bcftools statistics
- `reports/summary.txt` - Tab-delimited summary table

---

### Step 12: IGV Session Creation

**Script:** `variant_calling/create_igv_session.sh`

**Purpose:** Generate IGV session file for visual inspection.

#### Output Format

```xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" locus="All" version="8">
    <Resources>
        <Resource path="GRCh38_reference.fa"/>
        <Resource path="${PATIENT}-T_recalibrated.bam"/>
        <Resource path="${PATIENT}-N_recalibrated.bam"/>
        <Resource path="${PATIENT}_high_confidence.vcf.gz"/>
    </Resources>
</Session>
```

**Usage:** In IGV: File → Open Session → Select `igv_session.xml`

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                        REFERENCE DATA                                │
│  GRCh38_reference.fa ── BWA Index + SAMtools Index + GATK Dict      │
│  dbsnp_156.grch38.vcf.gz (remapped to chr naming)                   │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                        SAMPLE DATA                                   │
│  PATIENT-T_R1/R2.fastq.gz ──┬── bwa mem ── samtools sort ── .bam   │
│  PATIENT-N_R1/R2.fastq.gz ──┘                                       │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      BAM PROCESSING                                  │
│  .bam → MarkDuplicates → AddReadGroups → BQSR → _recalibrated.bam  │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      VARIANT CALLING                                 │
│  Tumor_recal.bam + Normal_recal.bam → Mutect2 → _raw.vcf.gz        │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      POST-PROCESSING                                 │
│  _raw.vcf.gz                                                        │
│      ↓ LearnReadOrientationModel + FilterMutectCalls                │
│  _filtered.vcf.gz                                                   │
│      ↓ bcftools norm                                                │
│  _filtered.vcf.gz (normalized)                                      │
│      ↓ VEP                                                          │
│  _annotated.vcf.gz                                                  │
│      ↓ bcftools view (configurable MIN_DP, MIN_VAF)                 │
│  _high_confidence.vcf.gz                                            │
│      ├── reports/*.txt                                              │
│      └── igv_session.xml                                            │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Directory Structure

```
$DATA_DIR/
├── reference/
│   ├── GRCh38_reference.fa
│   ├── GRCh38_reference.fa.fai
│   ├── GRCh38_reference.dict
│   ├── GRCh38_reference.fa.{bwt,pac,ann,amb,sa}
│   ├── dbsnp_156.grch38.vcf.gz
│   └── dbsnp_156.grch38.vcf.gz.tbi
├── fastq/
│   └── {PATIENT}-{T,N}_R{1,2}.fastq.gz
├── bam/
│   └── {PATIENT}-{T,N}.bam
├── bam_marked/
│   └── {PATIENT}-{T,N}_marked.bam
├── bam_with_rg/
│   └── {PATIENT}-{T,N}_rg.bam
├── bam_recalibrated/
│   └── {PATIENT}-{T,N}_recalibrated.bam
├── recal_tables/
│   └── {PATIENT}-{T,N}_recal_data.table
├── metrics/
│   └── {PATIENT}-{T,N}_dup_metrics.txt
├── vcfs/
│   ├── {PATIENT}_raw.vcf.gz
│   ├── {PATIENT}_f1r2.tar.gz
│   └── {PATIENT}_read-orientation-model.tar.gz
├── vcfs_filtered/
│   └── {PATIENT}_filtered.vcf.gz
├── vcfs_annotated/
│   └── {PATIENT}_annotated.vcf.gz
├── vcfs_prioritized/
│   └── {PATIENT}_high_confidence.vcf.gz
├── reports/
│   ├── {PATIENT}_stats.txt
│   └── summary.txt
└── igv_session.xml
```

---

## Quality Control

### Checkpoints

| Stage | Metric | Expected Value |
|-------|--------|----------------|
| FastQC | Read quality | Q30 > 80% |
| Alignment | Alignment rate | > 95% |
| Duplicates | Duplicate rate | < 20% |
| Coverage | Mean depth | > 30x |
| Variants | Ti/Tv ratio | ~2.0 for SNVs |
| Variants | Total count | 1000-10000 (typical tumor) |

### Quality Metrics

| Metric | Description |
|--------|-------------|
| Read Depth (DP) | Number of reads covering a position |
| Allele Frequency (AF) | Proportion of reads supporting variant |
| Quality Score (QUAL) | Confidence in variant call |
| Strand Bias (SB) | Variant appears on only one DNA strand |

---

## Troubleshooting

### Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| Low alignment rate (<90%) | Wrong reference or adapter contamination | Verify reference, trim adapters |
| High duplicate rate (>30%) | Over-amplification in library prep | Reduce PCR cycles, deeper sequencing |
| No variants called | Empty BAM or wrong parameters | Check BAM file size, verify inputs |
| Contig mismatch error | VCF and reference use different naming | Remap contig names with bcftools |
| Java version error | GATK 4.5 requires Java 17+ | Install Java 17+, update PATH |
| Memory errors | Insufficient RAM | Reduce threads, increase swap |

### Error Handling

All scripts include `set -e` to exit on first error. Critical validations:
- Input file existence checks
- Output file size validation (BQSR: >1MB indicates success)
- Contig name compatibility checks

---

## Running the Pipeline

### Full Pipeline
```bash
./run_complete_pipeline.sh
```

### Skip Downloads (data exists)
```bash
./run_complete_pipeline.sh --skip-download
```

### Force Regenerate All Outputs
```bash
./run_complete_pipeline.sh --skip-download --force
```

### Custom Settings
```bash
./run_complete_pipeline.sh \
    --project-dir ~/my_project \
    --data-dir ~/my_data \
    --num-processors 16
```

---

## References

- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [Mutect2 Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
- [VEP Documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [BWA Manual](http://bio-bwa.sourceforge.net/)
- [SAMtools/BCFtools](http://www.htslib.org/)
- [IGV User Guide](https://software.broadinstitute.org/software/igv/UserGuide)
