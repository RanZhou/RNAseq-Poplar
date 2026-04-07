#!/usr/bin/env Rscript
# run_analysis.R - Main analysis script
# Usage: Rscript run_analysis.R [config_file] [steps]
#
# Examples:
#   Rscript run_analysis.R                           # Use default config, run all steps
#   Rscript run_analysis.R config.yaml fpkm,qc       # Run only FPKM and QC
#   Rscript run_analysis.R custom_config.yaml        # Use custom config

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

config_file <- ifelse(length(args) >= 1, args[1], "config/config.yaml")
steps_str <- ifelse(length(args) >= 2, args[2], "fpkm,qc,de,go")
steps <- unlist(strsplit(steps_str, ","))

message("==========================================")
message("Tsai Lab RNA-seq Pipeline - Part 2")
message("==========================================")
message(sprintf("Config: %s", config_file))
message(sprintf("Steps: %s", paste(steps, collapse = ", ")))
message("")

# Load pipeline functions
source("R/config.R")
source("R/utils.R")
source("R/fpkm.R")
source("R/qc.R")
source("R/de.R")
source("R/go.R")

# Run pipeline
tryCatch({
  results <- run_pipeline(config_file, steps)
  message("\n✓ Analysis completed successfully!")
}, error = function(e) {
  message("\n✗ Error during analysis:")
  message(conditionMessage(e))
  quit(status = 1)
})