# utils.R - Utility Functions

#' Check and Load Required Packages
#' 
#' @param packages Vector of package names
#' @param Bioconductor Whether to use BiocManager for missing packages
#' @export
require_packages <- function(packages, Bioconductor = TRUE) {
  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing) > 0) {
    msg <- sprintf("Missing required packages: %s\n", paste(missing, collapse = ", "))
    
    if (Bioconductor) {
      msg <- paste0(msg, "Install with:\n")
      msg <- paste0(msg, "  install.packages('BiocManager')\n")
      for (pkg in missing) {
        msg <- paste0(msg, sprintf("  BiocManager::install('%s')\n", pkg))
      }
    } else {
      msg <- paste0(msg, "Install with: install.packages(c('", 
                    paste(missing, collapse = "', '"), "'))")
    }
    
    stop(msg)
  }
  
  invisible(TRUE)
}

#' Initialize RNASeq Pipeline
#' 
#' Loads config, creates directories, and validates setup
#' 
#' @param config_file Path to config file
#' @return Config list
#' @export
init_pipeline <- function(config_file = "config/config.yaml") {
  message("Initializing RNA-seq Pipeline...")
  message(sprintf("  Config file: %s", config_file))
  
  # Load config
  config <- load_config(config_file)
  
  # Create output directories
  create_output_dirs(config)
  
  # Print summary
  message("\nProject:")
  message(sprintf("  Name: %s", config$project$name))
  message(sprintf("  Organism: %s", config$project$organism))
  message(sprintf("  Samples: %d", nrow(load_design(config))))
  message(sprintf("  Comparisons: %d", nrow(load_comparisons(config))))
  
  return(config)
}

#' Run Complete Pipeline
#' 
#' Runs all steps: FPKM, QC, DE, GO
#' 
#' @param config Config list or path to config file
#' @param steps Which steps to run (default: all)
#' @return List with all results
#' @export
run_pipeline <- function(config = "config/config.yaml",
                          steps = c("fpkm", "qc", "de", "go")) {
  
  # Load config if path provided
  if (is.character(config)) {
    config <- init_pipeline(config)
  }
  
  results <- list()
  
  # Step 1: FPKM calculation
  if ("fpkm" %in% steps) {
    message("\n" , rep("=", 50))
    message("STEP 1: FPKM Calculation")
    message(rep("=", 50))
    
    results$fpkm <- calculate_fpkm(config)
    results$expressed <- count_expressed_genes(config, results$fpkm)
  }
  
  # Step 2: Quality Control
  if ("qc" %in% steps) {
    message("\n" , rep("=", 50))
    message("STEP 2: Quality Control")
    message(rep("=", 50))
    
    if (is.null(results$fpkm)) {
      results$fpkm <- load_fpkm(config)
    }
    results$qc <- run_sample_qc(config, results$fpkm)
  }
  
  # Step 3: Differential Expression
  if ("de" %in% steps) {
    message("\n" , rep("=", 50))
    message("STEP 3: Differential Expression")
    message(rep("=", 50))
    
    results$de <- run_de_analysis(config, results$counts, results$fpkm)
  }
  
  # Step 4: GO Enrichment
  if ("go" %in% steps) {
    message("\n" , rep("=", 50))
    message("STEP 4: GO Enrichment")
    message(rep("=", 50))
    
    results$go <- run_go_enrichment(config)
  }
  
  message("\n" , rep("=", 50))
  message("Pipeline Complete!")
  message(rep("=", 50))
  
  return(results)
}

#' Create Example Design File
#' 
#' Creates a template design file for user reference
#' 
#' @param output_file Path for output
#' @param n_samples Number of example samples
#' @export
create_example_design <- function(output_file = "examples/data/example_design.txt",
                                   n_samples = 6) {
  
  design <- data.frame(
    sample_id = paste0("Sample", 1:n_samples),
    File = paste0("sample", 1:n_samples, "_counts"),
    Group = rep(c("Control", "Treatment"), each = n_samples/2),
    stringsAsFactors = FALSE
  )
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.table(design, file = output_file, row.names = FALSE, sep = "\t")
  
  message(sprintf("Created example design file: %s", output_file))
  return(design)
}

#' Create Example Comparison File
#' 
#' Creates a template comparison file for user reference
#' 
#' @param output_file Path for output
#' @export
create_example_comparison <- function(output_file = "examples/data/example_compare.txt") {
  
  comp <- data.frame(
    Test = "Treatment_vs_Control",
    Group = "Group",
    Numerator = "Treatment",
    Denominator = "Control",
    stringsAsFactors = FALSE
  )
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.table(comp, file = output_file, row.names = FALSE, sep = "\t")
  
  message(sprintf("Created example comparison file: %s", output_file))
  return(comp)
}