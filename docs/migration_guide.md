# Migration Guide: Original → Refactored Pipeline

This guide helps you transition from the original `TsaiLabRNASeq_part2.Rmd` to the refactored version.

## Quick Migration Checklist

- [ ] Copy your data to `data/` and `db/` directories
- [ ] Create `config/config.yaml` with your file paths
- [ ] Update your design file column names if needed
- [ ] Run `Rscript run_analysis.R`
- [ ] Compare outputs to verify consistency

## Step-by-Step Migration

### Step 1: Backup Original Files

```bash
cp TsaiLabRNASeq_part2.Rmd TsaiLabRNASeq_part2.Rmd.backup
cp src/TsaiLabRNASeq_part2_2025ccdivide.R src/TsaiLabRNASeq_part2_2025ccdivide.R.backup
```

### Step 2: Set Up Directory Structure

```bash
mkdir -p data db config output
```

### Step 3: Create Config File

Copy `config/config.yaml` and edit:

```yaml
# OLD: Hardcoded in set_variables()
file_in <- "pdel_txi_counts_genelevel_2025.tsv"
file_design <- "pdel_design.txt"
file_gtf <- "PdelWV94v2.1.gtf"

# NEW: In config.yaml
input:
  counts_file: "data/pdel_txi_counts_genelevel_2025.tsv"
  design_file: "data/pdel_design.txt"
  
database:
  gtf_file: "db/PdelWV94v2.1.gtf"
```

### Step 4: Update Design File (if needed)

Ensure your design file has these **exact** column names:
- `sample_id` - descriptive sample name
- `File` - must match column names in counts file
- `Group` - biological group

### Step 5: Update Comparison File (if needed)

Ensure your comparison file has these **exact** column names:
- `Test` - comparison name
- `Group` - which design column to use
- `Numerator` - treatment group (UP when higher)
- `Denominator` - control group (DOWN when higher)

### Step 6: Run and Verify

```bash
# Run refactored pipeline
Rscript run_analysis.R config/config.yaml

# Compare outputs
ls -la data/DE_*
ls -la data/fpkm_*
```

## Key Differences

### Output File Locations

| Original | Refactored |
|----------|------------|
| `data\\fpkm_*` | `data/fpkm_*` (same location, cross-platform paths) |
| `data\\DE_*` | `data/DE_*` |
| `data\\GO\\` | `data/GO/` |
| `data\\pcc_*` | `output/pcc_*` (now in output dir) |

### Parameters

| Original | Refactored | Notes |
|----------|------------|-------|
| Hardcoded `qvalue_cut = 0.05` | Configurable in YAML | More flexible |
| Hardcoded `fold_change_cut = 2` | Configurable in YAML | More flexible |
| Hidden `cus_pval_cutoff=0.05` | Now `de$log2fc_threshold` | Documented |
| No p-value cutoff option | `de$pvalue_cutoff` | Can use p-value instead of q-value |

### Error Messages

The refactored pipeline provides **much** better error messages:

**Original:**
```
Error in sampleTable0$TestGroup: undefined columns selected
```

**Refactored:**
```
Error: Group column 'Group' not found in design file
  Available columns: sample_id, File, Group, Replicate
```

## Side-by-Side Comparison

### Original Rmd Chunk

```r
```{r fpkm_cal, echo=FALSE,include=FALSE}
source("src\\TsaiLabRNASeq_part2_2025ccdivide.R")
set_variables()
get_fpkms()
```
```

### Refactored Equivalent

```r
# Single initialization
config <- init_pipeline("config/config.yaml")

# Then run steps
fpkm <- calculate_fpkm(config)
```

### Original DE Analysis

```r
source("src\\TsaiLabRNASeq_part2_2025ccdivide.R")
set_variables()
qvalue_cut = 0.05
fold_change_cut = 2
fpkm_cut = 3
gene_list_for_go = 1
run_DEseq_General(gene_list_for_go,qvalue_cut,fold_change_cut,fpkm_cut)
```

### Refactored Equivalent

```r
# All parameters in config.yaml
de_results <- run_de_analysis(config)
```

## Troubleshooting Migration Issues

### "File not found" errors

Check that paths in `config.yaml` are relative to the config file location, or use absolute paths.

### Different results

The core DESeq2 calculations are identical. Differences may come from:
- Different random seeds (set `set.seed()` for reproducibility)
- Different filtering cutoffs (check `config.yaml` matches your original settings)
- Updated package versions (check with `sessionInfo()`)

### Missing GO results

Ensure `go$enabled: true` in config.yaml.

## Reverting to Original

If needed, restore from your backups:

```bash
cp TsaiLabRNASeq_part2.Rmd.backup TsaiLabRNASeq_part2.Rmd
cp src/TsaiLabRNASeq_part2_2025ccdivide.R.backup src/TsaiLabRNASeq_part2_2025ccdivide.R
```

## Need Help?

Contact:
- Pipeline issues: Check error messages (refactored version is more informative)
- Scientific questions: Tsai Lab (cjtsai@uga.edu)