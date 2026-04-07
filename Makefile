# Makefile for RNA-seq Pipeline
# Tsai Lab - UGA

.PHONY: all fpkm qc de go clean test install help

# Default: run full pipeline
all: fpkm qc de go

# Individual steps
fpkm:
	@echo "Running FPKM calculation..."
	Rscript run_analysis.R config/config.yaml fpkm

qc: fpkm
	@echo "Running QC..."
	Rscript run_analysis.R config/config.yaml qc

de: fpkm
	@echo "Running DE analysis..."
	Rscript run_analysis.R config/config.yaml de

go: de
	@echo "Running GO enrichment..."
	Rscript run_analysis.R config/config.yaml go

# Generate HTML report
report:
	@echo "Generating HTML report..."
	Rscript -e "rmarkdown::render('run_analysis.Rmd', output_file='output/report.html')"

# Clean output files
clean:
	@echo "Cleaning output files..."
	rm -f data/fpkm_*
	rm -f data/DE_*
	rm -f data/expressed_num_*
	rm -f data/pcc_*
	rm -rf data/GO/
	rm -f output/*

# Test with example data
test:
	@echo "Running test with example data..."
	Rscript -e "source('R/utils.R'); create_example_design(); create_example_comparison()"

# Check/install dependencies
install:
	@echo "Installing R dependencies..."
	Rscript -e "install.packages(c('yaml', 'ggplot2', 'ggdendro')); install.packages('BiocManager'); BiocManager::install(c('GenomicFeatures', 'DESeq2', 'topGO'))"

# Validate configuration
validate:
	@echo "Validating configuration..."
	Rscript -e "source('R/config.R'); load_config('config/config.yaml')"

# Show help
help:
	@echo "Tsai Lab RNA-seq Pipeline - Make Targets"
	@echo ""
	@echo "  make all      - Run complete pipeline (fpkm, qc, de, go)"
	@echo "  make fpkm     - Calculate FPKM values"
	@echo "  make qc       - Run sample QC (PCA, clustering)"
	@echo "  make de       - Run differential expression"
	@echo "  make go       - Run GO enrichment"
	@echo "  make report   - Generate HTML report"
	@echo "  make clean    - Remove all output files"
	@echo "  make test     - Create example input files"
	@echo "  make install  - Install R dependencies"
	@echo "  make validate - Check configuration"
	@echo "  make help     - Show this help"