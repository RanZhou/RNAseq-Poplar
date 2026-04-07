# de.R - Differential Expression Analysis

#' Run Differential Expression Analysis
#' 
#' Performs DESeq2 differential expression tests for all comparisons
#' defined in the compare file
#' 
#' @param config Configuration list
#' @param counts_data Optional pre-loaded counts
#' @param fpkm_data Optional pre-loaded FPKM values
#' @param force_recalc If TRUE, recalculate even if output exists
#' @return List containing DE results summary and full results table
#' @export
run_de_analysis <- function(config, counts_data = NULL, fpkm_data = NULL, 
                            force_recalc = FALSE) {
  
  message("Running Differential Expression Analysis...")
  
  # Load required packages
  require_packages(c("DESeq2"))
  
  # Load data
  if (is.null(counts_data)) {
    counts_data <- load_counts(config)
  }
  if (is.null(fpkm_data)) {
    fpkm_data <- load_fpkm(config)
  }
  
  design <- load_design(config)
  comparisons <- load_comparisons(config)
  annotations <- load_annotations(config)
  
  # Check if already computed
  de_num_file <- file.path(
    dirname(config$input$counts_file),
    paste0(config$output$de_prefix, basename(config$input$counts_file))
  )
  
  final_csv <- file.path(
    dirname(config$input$counts_file),
    paste0(config$output$de_prefix, basename(config$input$counts_file), "_Final_Out.csv")
  )
  
  if (file.exists(de_num_file) && file.exists(final_csv) && 
      !force_recalc && !config$runtime$overwrite_existing) {
    message("Loading existing DE results...")
    de_summary <- read.table(de_num_file, header = TRUE)
    return(list(summary = de_summary, full_results = NULL))
  }
  
  # Create output directory for GO lists
  dir.create(config$output$go_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Run comparisons
  results_list <- list()
  de_summary <- data.frame(
    Test = character(),
    DEG_Up = integer(),
    DEG_Down = integer(),
    Total_DEG = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(comparisons))) {
    test_name <- comparisons$Test[i]
    group_col <- comparisons$Group[i]
    num_group <- comparisons$Numerator[i]
    den_group <- comparisons$Denominator[i]
    
    message(sprintf("  Running comparison %d/%d: %s (%s vs %s)",
                    i, nrow(comparisons), test_name, num_group, den_group))
    
    # Run single comparison
    de_result <- run_single_comparison(
      counts = counts_data,
      fpkm = fpkm_data,
      design = design,
      test_name = test_name,
      group_col = group_col,
      num_group = num_group,
      den_group = den_group,
      config = config
    )
    
    results_list[[test_name]] <- de_result$results
    
    # Update summary
    de_summary <- rbind(de_summary, data.frame(
      Test = test_name,
      DEG_Up = de_result$n_up,
      DEG_Down = de_result$n_down,
      Total_DEG = de_result$n_up + de_result$n_down,
      stringsAsFactors = FALSE
    ))
  }
  
  # Combine all results
  message("  Combining results...")
  combined_results <- combine_de_results(results_list, fpkm_data, annotations)
  
  # Save outputs
  message(sprintf("  Saving DE summary to: %s", de_num_file))
  write.table(de_summary, file = de_num_file, row.names = FALSE, sep = "\t")
  
  message(sprintf("  Saving full results to: %s", final_csv))
  write.csv(combined_results, file = final_csv, row.names = FALSE, quote = TRUE)
  
  message("DE analysis complete.")
  
  return(list(summary = de_summary, full_results = combined_results))
}

#' Run Single Comparison
#' 
#' @param counts Count matrix
#' @param fpkm FPKM matrix
#' @param design Design data frame
#' @param test_name Test name
#' @param group_col Column in design containing groups
#' @param num_group Numerator group (treatment)
#' @param den_group Denominator group (control)
#' @param config Configuration list
#' @return List with results and DE gene counts
#' @keywords internal
run_single_comparison <- function(counts, fpkm, design, test_name, 
                                   group_col, num_group, den_group, config) {
  
  # Get grouping
  if (!group_col %in% colnames(design)) {
    stop(sprintf("Group column '%s' not found in design file", group_col))
  }
  
  groups <- as.character(design[[group_col]])
  
  # Check groups exist
  if (!num_group %in% groups) {
    stop(sprintf("Numerator group '%s' not found in design column '%s'", 
                 num_group, group_col))
  }
  if (!den_group %in% groups) {
    stop(sprintf("Denominator group '%s' not found in design column '%s'", 
                 den_group, group_col))
  }
  
  # Subset samples
  sel_samples <- groups %in% c(num_group, den_group)
  design_sub <- design[sel_samples, ]
  samples <- as.character(design_sub$File)
  
  # Check sample existence
  samples <- samples[samples %in% colnames(counts)]
  if (length(samples) == 0) {
    stop("No matching samples found in count data")
  }
  
  # Prepare data
  counts_sub <- counts[, samples, drop = FALSE]
  fpkm_sub <- fpkm[, samples, drop = FALSE]
  group_sub <- as.character(design_sub[[group_col]])
  
  # Create DESeq2 dataset
  colData <- data.frame(
    samples = samples,
    condition = factor(group_sub, levels = c(den_group, num_group)),
    row.names = samples
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = colData,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq2::DESeq(dds, quiet = !config$runtime$verbose)
  res <- DESeq2::results(dds, contrast = c("condition", num_group, den_group), 
                         cooksCutoff = FALSE)
  
  # Get FPKM expression tags
  fpkm_tags <- get_cut_tag(fpkm_sub, group_sub)
  
  # Prepare results
  result_df <- as.data.frame(res)
  result_df$fpkm_cut_tag <- fpkm_tags[rownames(result_df)]
  
  # Apply filters
  filtered <- filter_de_genes(
    result_df, 
    config$params$de$qvalue_cutoff,
    config$params$de$pvalue_cutoff,
    config$params$de$fold_change_cutoff,
    config$params$de$fpkm_expression_cutoff
  )
  
  # Generate GO gene lists (using alternative thresholds)
  go_genes_up <- result_df[which(result_df$log2FoldChange > config$params$de$log2fc_threshold & 
                                    result_df$pvalue < config$params$de$pvalue_cutoff), ]
  go_genes_down <- result_df[which(result_df$log2FoldChange < -config$params$de$log2fc_threshold & 
                                      result_df$pvalue < config$params$de$pvalue_cutoff), ]
  
  # Save GO lists
  if (nrow(go_genes_up) > 0) {
    write.table(
      rownames(go_genes_up),
      file = file.path(config$output$go_dir, sprintf("UP_%s.Pvaluecut.txt", test_name)),
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  if (nrow(go_genes_down) > 0) {
    write.table(
      rownames(go_genes_down),
      file = file.path(config$output$go_dir, sprintf("DO_%s.Pvaluecut.txt", test_name)),
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  return(list(
    results = result_df,
    n_up = nrow(go_genes_up),
    n_down = nrow(go_genes_down)
  ))
}

#' Filter DE Genes
#' 
#' @param result_df DESeq2 results data frame
#' @param qvalue_cut Q-value cutoff
#' @param pvalue_cut P-value cutoff (use 1 to disable)
#' @param fc_cut Fold change cutoff
#' @param fpkm_cut FPKM cutoff
#' @return Filtered results
#' @keywords internal
filter_de_genes <- function(result_df, qvalue_cut, pvalue_cut, fc_cut, fpkm_cut) {
  
  # Remove NAs
  valid <- !is.na(result_df$log2FoldChange) & 
    !is.na(result_df$fpkm_cut_tag)
  
  result_df <- result_df[valid, ]
  
  # Apply p-value or q-value filter
  if (pvalue_cut >= 1) {
    # Use q-value
    valid <- !is.na(result_df$padj) & 
      result_df$padj < qvalue_cut & 
      result_df$fpkm_cut_tag >= fpkm_cut & 
      abs(result_df$log2FoldChange) >= log2(fc_cut)
  } else {
    # Use p-value
    valid <- !is.na(result_df$pvalue) & 
      result_df$pvalue < pvalue_cut & 
      result_df$fpkm_cut_tag >= fpkm_cut & 
      abs(result_df$log2FoldChange) >= log2(fc_cut)
  }
  
  return(result_df[valid, ])
}

#' Combine DE Results
#' 
#' @param results_list List of per-comparison results
#' @param fpkm_data FPKM matrix
#' @param annotations Annotation data frame
#' @return Combined data frame
#' @keywords internal
combine_de_results <- function(results_list, fpkm_data, annotations) {
  
  # Extract columns from each comparison
  combined <- NULL
  for (test_name in names(results_list)) {
    res <- results_list[[test_name]]
    
    cols <- data.frame(
      logFC = res$log2FoldChange,
      Pvalue = res$pvalue,
      Qvalue = res$padj,
      FPKMcut = res$fpkm_cut_tag,
      row.names = rownames(res)
    )
    
    names(cols) <- paste0(test_name, c("_logFC", "_Pvalue", "_Qvalue", "_FPKMcut"))
    
    if (is.null(combined)) {
      combined <- cols
    } else {
      combined <- cbind(combined, cols[rownames(combined), ])
    }
  }
  
  # Add FPKM values
  combined <- cbind(combined, fpkm_data[rownames(combined), ])
  
  # Add annotations
  match_idx <- match(rownames(combined), annotations[[1]])
  anno_sub <- annotations[match_idx, ]
  rownames(anno_sub) <- rownames(combined)
  
  combined <- cbind(gene_id = rownames(combined), combined, anno_sub)
  
  return(combined)
}

#' Load Gene Annotations
#' 
#' @param config Configuration list
#' @return Annotation data frame
#' @export
load_annotations <- function(config) {
  if (!file.exists(config$database$annotation_file)) {
    stop(sprintf("Annotation file not found: %s", config$database$annotation_file))
  }
  
  anno <- read.table(
    config$database$annotation_file,
    header = TRUE,
    sep = "\t",
    comment.char = "",
    quote = '"',
    stringsAsFactors = FALSE
  )
  
  message(sprintf("  Loaded annotations: %d genes", nrow(anno)))
  return(anno)
}