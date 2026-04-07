# Changelog - Tsai Lab RNA-seq Pipeline Refactor

## Version 2.0.0 (Refactored)

### Breaking Changes
- Configuration moved from R code to YAML file
- Output file paths standardized (cross-platform)
- Function names changed to be more descriptive

### Major Improvements

#### Configuration Management
- **NEW**: YAML-based configuration (`config/config.yaml`)
- **NEW**: `load_config()` - validates and loads configuration
- **NEW**: `validate_inputs()` - checks all input files exist before running
- **NEW**: Cross-platform path handling with `file.path()`

#### Modularity
- **NEW**: Split into logical modules:
  - `R/config.R` - Configuration management
  - `R/fpkm.R` - FPKM calculation
  - `R/qc.R` - Quality control
  - `R/de.R` - Differential expression
  - `R/go.R` - GO enrichment
  - `R/utils.R` - Utilities

#### Error Handling
- **NEW**: Informative error messages with context
- **NEW**: Input validation at each step
- **NEW**: File existence checks before processing
- **NEW**: Sample name matching validation

#### Documentation
- **NEW**: Full Roxygen documentation for all functions
- **NEW**: README with usage examples
- **NEW**: Migration guide from original pipeline
- **NEW**: Example input files

#### Logging
- **NEW**: Progress messages for long-running steps
- **NEW**: Summary statistics printed during analysis

#### Quality of Life
- **NEW**: Command-line interface with `run_analysis.R`
- **NEW**: R Markdown report generation
- **NEW**: Makefile for common tasks
- **NEW**: Selective step execution (run only QC, or only DE)

### Bug Fixes
- Fixed Windows-specific backslash paths (`\\` → `/`)
- Fixed inconsistent DE filtering (hardcoded 1.25 FC now configurable)
- Fixed repeated `source()` calls in Rmd
- Fixed silent failures when files don't exist

### Code Quality
- Removed commented-out package installation code
- Eliminated global variable assignments (`<<-`)
- Added function parameters instead of hardcoded values
- Improved code organization and readability

### Performance
- Results caching (skip recalculation if output exists)
- Optional force recalculation flag

## Original Version (v1.0)

### Features
- FPKM calculation using DESeq2
- Sample QC (PCA + hierarchical clustering)
- Differential expression with DESeq2
- GO enrichment with topGO
- R Markdown report generation

### Limitations
- Hardcoded file paths
- Windows-only path separators
- Minimal error handling
- Single monolithic script
- No configuration file
- Poor documentation

---

## Future Enhancements (Planned)

- [ ] Unit tests with `testthat`
- [ ] Continuous integration (GitHub Actions)
- [ ] Docker containerization
- [ ] Shiny app for interactive results exploration
- [ ] Support for other DE tools (edgeR, limma)
- [ ] Gene set enrichment analysis (GSEA)
- [ ] Batch effect correction options
- [ ] Integration with GEO submission
- [ ] Automated figure generation for publication

## Contributors

- Original: Liangjiao Xue (liangjiao.xue@gmail.com)
- Lab: Tsai Lab, University of Georgia
- Refactor: OpenClaw Assistant