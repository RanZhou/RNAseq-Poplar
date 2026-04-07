# RNAseq-Poplar

**A refactored, config-driven RNA-seq analysis pipeline for Poplar (Populus) genomics.**

This pipeline performs:
- **FPKM calculation** from raw count data
- **Sample QC** (PCA, hierarchical clustering, Pearson correlation)
- **Differential Expression** analysis with DESeq2
- **GO Enrichment** analysis with topGO

Originally developed by the Tsai Lab at the University of Georgia, this version has been refactored for improved maintainability, cross-platform compatibility, and ease of use.

---

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/RanZhou/RNAseq-Poplar.git
cd RNAseq-Poplar

# Install dependencies (in R)
install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "DESeq2", "topGO"))

# Configure
vim config/config.yaml  # Edit with your file paths

# Run
Rscript run_analysis.R
```

---

## 📊 Pipeline Overview

| Step | Description | Output |
|------|-------------|--------|
| **FPKM** | Calculate fragments per kilobase million | `data/fpkm_*` |
| **QC** | PCA, clustering, correlation | `output/pcc_*.csv` + plots |
| **DE** | DESeq2 differential expression | `data/DE_*` |
| **GO** | GO enrichment (BP, MF, CC) | `data/GO/`, `GO_enrichment_*_merge.tab` |

---

## 📁 Repository Structure

```
RNAseq-Poplar/
├── config/
│   └── config.yaml              # Main configuration file
├── R/
│   ├── config.R                 # Config management
│   ├── fpkm.R                   # FPKM calculation
│   ├── qc.R                     # Quality control
│   ├── de.R                     # Differential expression
│   ├── go.R                     # GO enrichment
│   └── utils.R                  # Utilities
├── examples/
│   └── data/                    # Example input files
├── docs/
│   └── migration_guide.md       # Migration from original
├── run_analysis.R               # CLI runner
├── run_analysis.Rmd             # R Markdown report
├── Makefile                     # Common tasks
└── README.md                    # This file
```

---

## 🛠️ Configuration

Edit `config/config.yaml`:

```yaml
project:
  name: "Poplar_RNASeq_Part2"
  organism: "Populus deltoides"

input:
  counts_file: "data/pdel_txi_counts_genelevel_2025.tsv"
  design_file: "data/pdel_design.txt"
  compare_file: "data/pdel_compare.txt"

database:
  gtf_file: "db/PdelWV94v2.1.gtf"
  annotation_file: "db/pdel.gene_desc.2026.txt"
  go_map_file: "db/pdel.uni_jgi.GO.list"

params:
  fpkm_cutoff: 5
  de:
    qvalue_cutoff: 0.05
    fold_change_cutoff: 2
```

---

## 📖 Documentation

- [Full README](README.md) - Detailed documentation
- [Migration Guide](docs/migration_guide.md) - Moving from the original pipeline
- [Changelog](CHANGES.md) - Version history

---

## 🔬 Original Pipeline

This is a refactored version of the Tsai Lab RNA-seq pipeline:
- **Original**: https://github.com/TsailabBioinformatics/Transcriptomic_part2_R_2026updated.git
- **Authors**: Ran Zhou, Liangjiao Xue, Tsai Lab UGA
- **Contact**: Ran.Zhou@uga.edu, cjtsai@uga.edu

---

## 📄 License

Original pipeline by Tsai Lab, University of Georgia.
This refactored version maintains the same license terms.
