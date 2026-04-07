# fpkm.R - FPKM Calculation Functions

#' Calculate FPKM Values from Count Data
#' 
#' Uses DESeq2 to calculate FPKM values from raw counts
#' and GTF annotation
#' 
#' @param config Configuration list from load_config()
#' @param force_recalc If TRUE, recalculate even if output exists
#' @return Data frame with FPKM values
#' @export
calculate_fpkm <- function(config, force_recalc = FALSE) {
  
  # Check if already computed
  fpkm_file <- get_output_path(config, "fpkm", config$input$counts_file)
  
  if (file.exists(fpkm_file) && !force_recalc && !config$runtime$overwrite_existing) {
    message(sprintf("Loading existing FPKM file: %s", fpkm_file))
    return(read.table(fpkm_file, header = TRUE, row.names = 1))
  }
  
  message("Calculating FPKM values...")
  
  # Load required packages
  require_packages(c("GenomicFeatures", "DESeq2"))
  
  # Load data
  message("  Loading count data and annotations...")
  counts <- load_counts(config)
  design <- load_design(config)
  
  # Create TxDb from GTF
  message("  Creating transcript database from GTF...")
  txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = config$database$gtf_file,
    format = "gtf",
    circ_seqs = character()
  )
  ebg <- GenomicFeatures::exonsBy(txdb, by = "gene")
  
  # Match samples
  samples <- colnames(counts)
  file_id <- as.character(design$File)
  site <- match(samples, file_id)
  
  if (any(is.na(site))) {
    missing <- samples[is.na(site)]
    stop(sprintf("Samples in counts not found in design file: %s",
                 paste(missing, collapse = ", ")))
  }
  
  group_real <- as.character(design$Group)[site]
  
  # Handle unassigned samples
  if (any(is.na(group_real))) {
    warning("Some samples have NA groups - assigning to 'unassigned'")
    group_real[is.na(group_real)] <- "unassigned"
  }
  
  # Create DESeq2 dataset
  message("  Creating DESeq2 dataset...")
  colData <- data.frame(
    samples = samples,
    condition = group_real,
    row.names = samples
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = colData,
    design = ~ condition
  )
  
  # Add exon information
  SummarizedExperiment::rowRanges(dds) <- ebg
  
  # Calculate FPKM
  message("  Calculating FPKM...")
  fpkm_out <- DESeq2::fpkm(dds)
  
  # Save output
  message(sprintf("  Saving FPKM to: %s", fpkm_file))
  write.table(
    as.data.frame(fpkm_out),
    file = fpkm_file,
    sep = "\t",
    quote = FALSE
  )
  
  message("FPKM calculation complete.")
  return(as.data.frame(fpkm_out))
}

#' Count Genes Above FPKM Cutoff
#' 
#' Calculates number of genes with FPKM >= cutoff in all replicates
#' of each group
#' 
#' @param config Configuration list
#' @param fpkm_data Optional pre-loaded FPKM data
#' @param cutoff FPKM cutoff value (default from config)
#' @return Data frame with gene counts per group
#' @export
count_expressed_genes <- function(config, fpkm_data = NULL, cutoff = NULL) {
  
  if (is.null(cutoff)) {
    cutoff <- config$params$fpkm_cutoff
  }
  
  message(sprintf("Counting genes with FPKM >= %d in all replicates...", cutoff))
  
  # Load data
  design <- load_design(config)
  
  if (is.null(fpkm_data)) {
    fpkm_data <- load_fpkm(config)
  }
  
  # Calculate per group
  groups <- unique(as.character(design$Group))
  results <- data.frame(
    Group = character(),
    Gene_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (group in groups) {
    samples <- as.character(design$File[design$Group == group])
    samples <- samples[samples %in% colnames(fpkm_data)]
    
    if (length(samples) == 0) {
      warning(sprintf("No samples found for group: %s", group))
      next
    }
    
    fpkm_group <- fpkm_data[, samples, drop = FALSE]
    
    # Get min FPKM per gene across samples
    if (length(samples) == 1) {
      fpkm_min <- fpkm_group
    } else {
      fpkm_min <- apply(fpkm_group, 1, min)
    }
    
    # Count genes above cutoff
    n_genes <- sum(fpkm_min >= cutoff)
    
    results <- rbind(results, data.frame(
      Group = group,
      Gene_Count = n_genes,
      stringsAsFactors = FALSE
    ))
    
    message(sprintf("  %s: %d genes", group, n_genes))
  }
  
  # Save results
  out_file <- file.path(
    dirname(config$input$counts_file),
    paste0(config$output$expressed_prefix, basename(config$input$counts_file))
  )
  
  write.table(results, file = out_file, row.names = FALSE, sep = "\t")
  message(sprintf("  Saved to: %s", out_file))
  
  return(results)
}

#' Load Count Data
#' 
#' @param config Configuration list
#' @return Count matrix
#' @export
load_counts <- function(config) {
  if (!file.exists(config$input$counts_file)) {
    stop(sprintf("Counts file not found: %s", config$input$counts_file))
  }
  
  data <- read.table(
    config$input$counts_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  message(sprintf("  Loaded counts: %d genes x %d samples", nrow(data), ncol(data)))
  return(data)
}

#' Load FPKM Data
#' 
#' @param config Configuration list
#' @return FPKM matrix
#' @export
load_fpkm <- function(config) {
  fpkm_file <- get_output_path(config, "fpkm", config$input$counts_file)
  
  if (!file.exists(fpkm_file)) {
    stop(sprintf("FPKM file not found: %s\nRun calculate_fpkm() first.", fpkm_file))
  }
  
  data <- read.table(
    fpkm_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  message(sprintf("  Loaded FPKM: %d genes x %d samples", nrow(data), ncol(data)))
  return(data)
}

#' Load Design File
#' 
#' @param config Configuration list
#' @return Design data frame
#' @export
load_design <- function(config) {
  if (!file.exists(config$input$design_file)) {
    stop(sprintf("Design file not found: %s", config$input$design_file))
  }
  
  design <- read.table(
    config$input$design_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Validate required columns
  required <- c("sample_id", "File", "Group")
  missing <- required[!required %in% colnames(design)]
  if (length(missing) > 0) {
    stop(sprintf("Design file missing required columns: %s",
                 paste(missing, collapse = ", ")))
  }
  
  message(sprintf("  Loaded design: %d samples", nrow(design)))
  return(design)
}

#' Load Comparison File
#' 
#' @param config Configuration list
#' @return Comparison data frame
#' @export
load_comparisons <- function(config) {
  if (!file.exists(config$input$compare_file)) {
    stop(sprintf("Comparison file not found: %s", config$input$compare_file))
  }
  
  comp <- read.table(
    config$input$compare_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Validate columns
  required <- c("Test", "Group", "Numerator", "Denominator")
  missing <- required[!required %in% colnames(comp)]
  if (length(missing) > 0) {
    stop(sprintf("Comparison file missing required columns: %s\nExpected: Test, Group, Numerator, Denominator",
                 paste(missing, collapse = ", ")))
  }
  
  message(sprintf("  Loaded comparisons: %d tests", nrow(comp)))
  return(comp)
}

#' Get Output File Path
#' 
#' @param config Configuration list
#' @param type Output type: "fpkm", "de", "expressed"
#' @param input_file Input file path (for naming)
#' @return Output file path
#' @keywords internal
get_output_path <- function(config, type, input_file) {
  prefix <- switch(type,
                   "fpkm" = config$output$fpkm_prefix,
                   "de" = config$output$de_prefix,
                   "expressed" = config$output$expressed_prefix,
                   stop(sprintf("Unknown output type: %s", type))
  )
  
  out_file <- file.path(
    dirname(input_file),
    paste0(prefix, basename(input_file))
  )
  
  return(out_file)
}