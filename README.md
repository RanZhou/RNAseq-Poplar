# Tsai Lab RNA-seq Pipeline - Part 2 (Refactored)

A refactored, config-driven RNA-seq analysis pipeline for differential expression and GO enrichment analysis.

## What's New

| Feature | Original | Refactored |
|---------|----------|------------|
| Configuration | Hardcoded in R | YAML config file |
| Paths | Windows backslashes (`\\`) | Cross-platform (`file.path()`) |
| Error handling | Minimal | Comprehensive validation |
| Documentation | Inline comments | Full function documentation |
| Modularity | Single script | Modular R package structure |
| Logging | None | Progress messages |
| Package management | Inline installs | Pre-check with clear errors |

## Installation

### Prerequisites

```r
# Required R packages
install.packages(c("yaml", "ggplot2", "ggdendro"))

# Bioconductor packages
install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "DESeq2", "topGO"))
```

### Setup

```bash
git clone https://github.com/TsailabBioinformatics/Transcriptomic_part2_R_2026updated.git
cd Transcriptomic_part2_R_2026updated

# Copy refactored pipeline
cp -r ~/oc_rnapipe/* .
```

## Quick Start

### 1. Configure

Edit `config/config.yaml`:

```yaml
input:
  counts_file: "data/your_counts.tsv"
  design_file: "data/your_design.txt"
  compare_file: "data/your_compare.txt"

database:
  gtf_file: "db/your_genome.gtf"
  annotation_file: "db/your_annotations.txt"
  go_map_file: "db/your_go_mappings.txt"
```

### 2. Run

**Command line:**
```bash
Rscript run_analysis.R
```

**R Markdown (with plots):**
```bash
Rscript -e "rmarkdown::render('run_analysis.Rmd')"
```

**From R:**
```r
source("R/config.R")
source("R/utils.R")
source("R/fpkm.R")
source("R/qc.R")
source("R/de.R")
source("R/go.R")

# Initialize
config <- init_pipeline("config/config.yaml")

# Run full pipeline
results <- run_pipeline(config)

# Or run specific steps
fpkm <- calculate_fpkm(config)
qc <- run_sample_qc(config, fpkm)
de <- run_de_analysis(config)
go <- run_go_enrichment(config)
```

## Input File Formats

### Design File (`*_design.txt`)

| Column | Description |
|--------|-------------|
| `sample_id` | Human-readable sample name |
| `File` | Column name in counts file (must match) |
| `Group` | Biological group for DE testing |

Example:
```
sample_id	File	Group
WT_Control_1	sample1_counts	WT_Control
WT_Control_2	sample2_counts	WT_Control
Treatment_1	sample3_counts	Treatment
Treatment_2	sample4_counts	Treatment
```

### Comparison File (`*_compare.txt`)

| Column | Description |
|--------|-------------|
| `Test` | Name for this comparison |
| `Group` | Design column to use for grouping |
| `Numerator` | Treatment group (UP when higher) |
| `Denominator` | Control group (DOWN when higher) |

Example:
```
Test	Group	Numerator	Denominator
Treatment_vs_Control	Group	Treatment	WT_Control
```

## Configuration Options

### Analysis Parameters

```yaml
params:
  fpkm_cutoff: 5                    # Min FPKM for expressed gene count
  pca_min_fpkm: 5                   # Min FPKM for PCA filtering
  
  de:
    qvalue_cutoff: 0.05             # FDR threshold
    pvalue_cutoff: 1                # Raw p-value threshold (1=use q-value)
    fold_change_cutoff: 2           # Absolute FC threshold
    fpkm_expression_cutoff: 3       # Min FPKM for DE filtering
    log2fc_threshold: 1.25          # FC for GO list generation
```

### GO Enrichment Settings

```yaml
go:
  enabled: true
  node_size: 3
  ontologies: ["BP", "MF", "CC"]    # Which GO categories to test
  algorithm: "classic"
  statistic: "fisher"
```

## File Structure

```
oc_rnapipe/
├── config/
│   └── config.yaml          # Main configuration
├── R/
│   ├── config.R             # Configuration management
│   ├── fpkm.R               # FPKM calculation
│   ├── qc.R                 # Quality control (PCA, clustering)
│   ├── de.R                 # Differential expression
│   ├── go.R                 # GO enrichment
│   └── utils.R              # Utility functions
├── examples/
│   ├── data/                # Example input files
│   └── db/                  # Example database files
├── docs/
│   └── migration_guide.md   # Migrating from old pipeline
├── run_analysis.R           # Command-line script
├── run_analysis.Rmd         # R Markdown report
└── README.md                # This file
```

## Migration from Original

### Path Changes

Replace in your files:
- `src\\TsaiLabRNASeq_part2_2025ccdivide.R` → Pipeline auto-handles paths
- `data\\` → Use config file paths
- `db\\` → Use config file paths

### Function Mapping

| Original | New |
|----------|-----|
| `set_variables()` | Edit `config/config.yaml` |
| `get_fpkms()` | `calculate_fpkm(config)` |
| `get_gene_num_above_cutoff()` | `count_expressed_genes(config)` |
| `run_PCA_cluster()` | `run_sample_qc(config)` |
| `run_DEseq_General()` | `run_de_analysis(config)` |
| `run_go_enrichment()` | `run_go_enrichment(config)` |

## Troubleshooting

### "Missing required packages"
Install the listed Bioconductor packages using `BiocManager::install()`

### "Samples in design not found in FPKM data"
Check that `File` column in design file exactly matches column names in counts file

### "Numerator group not found in design"
Check that the group name in compare file matches values in the design file's Group column

### Windows paths not working
The refactored pipeline uses `file.path()` which works on all platforms. No manual path editing needed.

## Contributing

This is a lab-specific pipeline. For bug reports or feature requests, contact:
- Liangjiao Xue: liangjiao.xue@gmail.com
- Tsai Lab: cjtsai@uga.edu

## License

Original pipeline by Tsai Lab, University of Georgia.
Refactored version maintains same license terms.