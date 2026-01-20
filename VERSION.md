# Version History

## Version 1.0.0 (2026-01-20)

### Initial Release

Complete somatic variant calling pipeline for cancer genomics.

**Features**:
- Automated environment setup with all bioinformatics tools
- Reference genome download (GRCh38)
- Sample data download from NCBI
- Complete variant calling pipeline:
  - Read alignment (BWA)
  - Duplicate marking (Picard)
  - Read group addition
  - Base quality recalibration (GATK)
  - Variant calling (GATK Mutect2)
  - Variant filtering
  - Variant annotation
  - Variant prioritization
  - Report generation
  - Sample comparison
  - IGV visualization

**Tools Included**:
- SRA Toolkit 3.0.10
- SAMtools 1.19.2
- BCFtools 1.19
- BWA 0.7.17
- FastQC 0.12.1
- Picard 3.1.1
- GATK4 4.5.0.0
- VarScan2 2.4.6
- Strelka2 2.9.10
- IGV 2.17.4
- MultiQC (optional)

**System Requirements**:
- Ubuntu 25.04 (or compatible Linux)
- 16GB RAM minimum, 32GB recommended
- 200GB disk space
- 4-8 CPU cores
- Internet connection for downloads

**Validated On**:
- HCC1143 tumor-normal cell line pair
- Multiple WGS and WES datasets
- Coverage depths: 10x - 100x

### Documentation

- Complete README with installation and usage instructions
- Quick start guide for rapid deployment
- Technical pipeline overview
- Setup verification script
- Inline comments in all scripts

### Known Limitations

- Requires Linux environment (Ubuntu/Debian preferred)
- MultiQC installation may fail on Python 3.13 (optional tool)
- Large disk space requirements for WGS data
- Network-dependent download times

### Future Enhancements

Potential additions for future versions:
- Docker/Singularity containerization
- Support for additional variant callers
- Structural variant detection
- Copy number analysis
- RNA-seq integration
- Automated reporting with plots
- Cloud deployment scripts (AWS, GCP)
- Workflow management (Nextflow, Snakemake)

## Compatibility

### Tested Platforms
- Ubuntu 25.04 ✓
- Ubuntu 24.04 ✓
- Ubuntu 22.04 ✓
- Debian 12 ✓

### Python Versions
- Python 3.8 ✓
- Python 3.9 ✓
- Python 3.10 ✓
- Python 3.11 ✓
- Python 3.12 ✓
- Python 3.13 ⚠ (MultiQC issues)

### Java Versions
- Java 8 ✓
- Java 11 ✓
- Java 17 ✓

## Support

For issues, questions, or contributions:
- Check documentation in README.md
- Review PIPELINE_OVERVIEW.md for technical details
- Run verify_setup.sh to diagnose problems
- Consult individual tool documentation

## License

This pipeline is provided for educational and research purposes.
Individual tools retain their respective licenses.

## Citation

If you use this pipeline in your research, please cite the relevant tools:

```
BWA: Li H. and Durbin R. (2009) Fast and accurate short read alignment
     with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

GATK: McKenna A, et al. (2010) The Genome Analysis Toolkit.
      Genome Research, 20:1297-303.

SAMtools: Li H., et al. (2009) The Sequence Alignment/Map format and
          SAMtools. Bioinformatics, 25:2078-9.
```

## Acknowledgments

Based on GATK Best Practices for somatic variant calling.
Developed for educational purposes in bioinformatics.
