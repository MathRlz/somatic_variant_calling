# Quick Start Guide

## Fastest Way to Get Started

### Complete Pipeline (First Time Users)

```bash
# 1. Clone or download this repository
cd /path/to/repo

# 2. Run the complete pipeline (this will do everything)
./run_complete_pipeline.sh
```

That's it! The script will:
- Install all tools (~30-60 minutes)
- Download reference genome (~30-60 minutes)
- Download sample data (~1-2 hours)
- Run complete analysis (~4-8 hours)

Total time: 6-12 hours depending on your internet and CPU

### If You Already Have Tools Installed

```bash
./run_complete_pipeline.sh --skip-setup
```

### If You Already Have Data Downloaded

```bash
./run_complete_pipeline.sh --skip-setup --skip-download
```

### If You Have BAM Files Ready

```bash
./run_complete_pipeline.sh --skip-setup --skip-download --skip-alignment
```

## Default Locations

- **Tools**: `~/somatic_variant_calling/`
- **Data**: `~/storage/variant_calling_data/`

## Custom Locations

```bash
./run_complete_pipeline.sh \
    --project-dir /path/to/tools \
    --data-dir /path/to/data
```

## After Pipeline Completes

Results are in `~/storage/variant_calling_data/`:

```bash
# View high-confidence variants
cd ~/storage/variant_calling_data/vcfs_prioritized
ls -lh

# View reports
cd ~/storage/variant_calling_data/reports
cat *.txt

# View tumor-specific variants
cd ~/storage/variant_calling_data/comparisons
ls -lh
```

## Visualize in IGV

```bash
cd ~/somatic_variant_calling/software/IGV_Linux_2.17.4
./igv.sh
# Then: File > Open Session > ~/storage/variant_calling_data/igv_session.xml
```

## Help

```bash
./run_complete_pipeline.sh --help
```

## Common Issues

### Out of disk space
You need at least 200GB free. Check with:
```bash
df -h ~
```

### Memory errors
You need at least 16GB RAM. Check with:
```bash
free -h
```

### Tool not found
Make sure environment is sourced:
```bash
source ~/somatic_variant_calling/setup_env.sh
```

## More Information

See [README.md](README.md) for detailed documentation.
