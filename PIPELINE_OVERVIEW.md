# Somatic Variant Calling Pipeline - Technical Overview

## What is Somatic Variant Calling?

Somatic variant calling is the process of identifying genetic mutations that occur in cancer cells but are not present in normal tissue. These mutations are acquired during a person's lifetime and can drive cancer development.

## Why Compare Tumor and Normal?

By sequencing both tumor tissue and matched normal tissue (usually blood) from the same patient, we can:

1. **Identify true somatic mutations** - Variants unique to the tumor
2. **Filter germline variants** - Inherited variants present in both samples
3. **Reduce false positives** - Technical artifacts that appear in both samples
4. **Enable precision medicine** - Identify targetable mutations for treatment

## Pipeline Architecture

### Phase 1: Data Preparation
**Scripts**: `data/download_*.sh`, `data/align_all.sh`

1. **Download Reference Genome** (GRCh38)
   - Human reference genome (~3 billion base pairs)
   - Known variant sites for quality control

2. **Download Sample Data**
   - Tumor tissue sequencing data (FASTQ files)
   - Matched normal tissue sequencing data
   - Quality metrics for each read

3. **Read Alignment**
   - Align millions of short reads to reference genome
   - Uses BWA (Burrows-Wheeler Aligner)
   - Output: BAM (Binary Alignment Map) files

### Phase 2: BAM Processing
**Scripts**: `variant_calling/mark_duplicates.sh`, `variant_calling/add_read_groups.sh`, `variant_calling/bqsr.sh`

4. **Mark Duplicates**
   - Remove PCR duplicates (artificial copies)
   - Remove optical duplicates (sequencing artifacts)
   - Prevents overcounting of variants

5. **Add Read Groups**
   - Add metadata to each read
   - Sample name, library, sequencing platform
   - Required for multi-sample analysis

6. **Base Quality Score Recalibration (BQSR)**
   - Correct systematic errors in quality scores
   - Uses known variant sites as truth set
   - Improves variant calling accuracy

### Phase 3: Variant Calling
**Scripts**: `variant_calling/call_variants.sh`

7. **Call Variants with GATK Mutect2**
   - Compare tumor and normal samples
   - Identify positions where tumor differs
   - Calculate statistical confidence for each variant
   - Output: VCF (Variant Call Format) files

### Phase 4: Variant Refinement
**Scripts**: `variant_calling/filter_variants.sh`, `variant_calling/annotate_variants.sh`, `variant_calling/prioritize_variants.sh`

8. **Filter Variants**
   - Apply quality filters
   - Remove low-confidence calls
   - Filter by read depth, allele frequency, strand bias

9. **Annotate Variants**
   - Add gene names and locations
   - Predict functional impact
   - Identify coding vs. non-coding changes
   - Use databases like dbSNP, ClinVar, COSMIC

10. **Prioritize Variants**
    - Select high-confidence somatic mutations
    - Focus on clinically relevant genes
    - Identify actionable variants

### Phase 5: Analysis and Reporting
**Scripts**: `variant_calling/generate_reports.sh`, `variant_calling/compare_samples.sh`, `variant_calling/create_igv_session.sh`

11. **Generate Reports**
    - Summary statistics for each sample
    - Variant counts by type
    - Quality metrics

12. **Compare Samples**
    - Identify tumor-specific variants
    - Remove germline variants
    - Create lists of true somatic mutations

13. **Create IGV Session**
    - Set up visualization file
    - Load reference, alignments, and variants
    - Enable manual inspection of variants

## Key Concepts

### Variant Types

- **SNV (Single Nucleotide Variant)**: Single base pair change (A→G)
- **Indel**: Insertion or deletion of bases
- **Somatic**: Variant present only in tumor
- **Germline**: Variant present in normal tissue (inherited)

### Quality Metrics

- **Read Depth**: Number of reads covering a position
- **Allele Frequency**: Proportion of reads supporting variant
- **Quality Score**: Confidence in base call
- **Strand Bias**: Variant appears on only one DNA strand

### File Formats

- **FASTQ**: Raw sequencing reads with quality scores
- **BAM**: Binary alignment file (compressed, indexed)
- **VCF**: Variant Call Format (lists all variants)
- **BED**: Genomic regions (chr, start, end)

## Data Flow

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

## Required Resources

### Computational

- **CPU**: Multi-core processor (4-8 cores minimum)
  - Alignment and variant calling are parallelizable
  - More cores = faster processing

- **RAM**: 16-32GB
  - Reference genome loading: ~5GB
  - GATK operations: 4-8GB per sample
  - Peak usage during variant calling

- **Storage**: 200GB minimum
  - Reference genome: 30GB
  - Sample data: 50-100GB per sample
  - Intermediate files: 50GB
  - Final results: 5-10GB

### Time

For 30x whole genome sequencing (typical for cancer):

- Setup: 1 hour (one time)
- Download reference: 1 hour (one time)
- Download samples: 2-4 hours per sample
- Alignment: 4-8 hours per sample
- BAM processing: 2-4 hours per sample
- Variant calling: 2-4 hours
- Filtering and annotation: 30 minutes
- **Total: 12-24 hours per tumor-normal pair**

For targeted sequencing (exomes):
- 2-4 hours total (10x faster)

## Quality Control Checkpoints

1. **FastQC**: Check raw read quality
2. **Alignment rate**: Should be >95% for human samples
3. **Duplicate rate**: Should be <20% for good libraries
4. **Coverage depth**: Should be >30x for somatic calling
5. **Transition/Transversion ratio**: Should be ~2.0 for SNVs
6. **Variant count**: Typical tumor has 1000-10000 somatic variants

## Common Issues and Solutions

### Low Alignment Rate (<90%)
- Check if correct reference genome is used
- Check for adapter contamination
- Check for species mismatch

### High Duplicate Rate (>30%)
- May indicate over-amplification during library prep
- Not necessarily fatal, but reduces effective coverage
- Modern callers can handle this

### Too Many/Few Variants
- Too many: Possible contamination or alignment issues
- Too few: Low coverage or overly strict filtering
- Compare to expected ranges for cancer type

### Memory Errors
- Reduce number of threads
- Process smaller genomic regions
- Use more swap space

### Disk Space Errors
- Delete intermediate files after each step
- Use external storage for data directory
- Compress old BAM files

## Output Interpretation

### High-Confidence Variants (vcfs_prioritized/)

These are the most reliable somatic mutations:
- High quality scores (QUAL > 100)
- Sufficient read depth (DP > 20)
- Reasonable allele frequency (0.1 < AF < 0.9)
- Not present in normal sample
- Pass all filters

### Tumor-Specific Variants (comparisons/)

Variants present in tumor but not in matched normal:
- True somatic mutations
- Candidate driver mutations
- Potential therapeutic targets
- May need validation with orthogonal method

### Clinical Actionability

Look for variants in:
- **Tier 1**: Known cancer genes with approved therapies
- **Tier 2**: Cancer genes with experimental therapies
- **Tier 3**: Genes of unknown significance

## Best Practices

1. **Always use matched normal** - Essential for somatic calling
2. **Aim for high coverage** - 30x minimum, 60-100x better
3. **Validate key findings** - Use Sanger sequencing or ddPCR
4. **Check quality metrics** - Don't trust variants blindly
5. **Use latest references** - GRCh38 is current standard
6. **Document everything** - Track parameters and versions
7. **Backup data** - Raw data is irreplaceable

## Further Reading

- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [Cancer Genomics Cloud](https://www.cancergenomicscloud.org/)
- [ClinVar Database](https://www.ncbi.nlm.nih.gov/clinvar/)
- [COSMIC Database](https://cancer.sanger.ac.uk/cosmic)
- [IGV User Guide](https://software.broadinstitute.org/software/igv/UserGuide)

## Pipeline Validation

This pipeline follows GATK Best Practices for somatic variant calling and has been tested with:
- HCC1143 cell line (tumor-normal pair)
- Multiple WGS and WES datasets
- Various coverage depths (10x - 100x)

Expected outputs match published benchmarks for sensitivity and specificity.
