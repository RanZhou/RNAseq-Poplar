# go.R - GO Enrichment Analysis

#' Run GO Enrichment Analysis
#' 
#' Performs GO enrichment using topGO for Biological Process,
#' Molecular Function, and Cellular Component
#' 
#' @param config Configuration list
#' @param data_folder Folder containing gene list files
#' @param force_recalc If TRUE, recalculate even if output exists
#' @return List with results for each ontology
#' @export
run_go_enrichment <- function(config, data_folder = NULL, force_recalc = FALSE) {
  
  if (is.null(data_folder)) {
    data_folder <- config$output$go_dir
  }
  
  if (!config$params$go$enabled) {
    message("GO enrichment disabled in config.")
    return(NULL)
  }
  
  message("Running GO Enrichment Analysis...")
  
  # Load required packages
  require_packages(c("topGO"))
  
  # Load GO mappings
  message("  Loading GO mappings...")
  if (!file.exists(config$database$go_map_file)) {
    stop(sprintf("GO map file not found: %s", config$database$go_map_file))
  }
  
  geneID2GO <- topGO::readMappings(file = config$database$go_map_file)
  geneNames <- names(geneID2GO)
  
  # Load GO descriptions
  go_desc <- load_go_descriptions(config)
  
  # Get gene list files
  gene_files <- list.files(
    path = data_folder,
    pattern = "\\.txt$",
    full.names = FALSE
  )
  
  if (length(gene_files) == 0) {
    warning(sprintf("No gene list files found in: %s", data_folder))
    return(NULL)
  }
  
  message(sprintf("  Found %d gene list files", length(gene_files)))
  
  # Run for each ontology
  results <- list()
  
  for (onto in config$params$go$ontologies) {
    message(sprintf("  Processing %s...", onto))
    
    onto_results <- run_ontology_enrichment(
      onto = onto,
      gene_files = gene_files,
      data_folder = data_folder,
      geneNames = geneNames,
      geneID2GO = geneID2GO,
      go_desc = go_desc,
      config = config
    )
    
    results[[onto]] <- onto_results
  }
  
  message("GO enrichment complete.")
  return(results)
}

#' Run Enrichment for Single Ontology
#' 
#' @param onto Ontology: BP, MF, or CC
#' @param gene_files Vector of gene list filenames
#' @param data_folder Path to gene list folder
#' @param geneNames All gene names in reference
#' @param geneID2GO GO mapping list
#' @param go_desc GO description data frame
#' @param config Configuration list
#' @return Data frame with merged results
#' @keywords internal
run_ontology_enrichment <- function(onto, gene_files, data_folder, 
                                     geneNames, geneID2GO, go_desc, config) {
  
  # Create ontology-specific output directory
  onto_dir <- file.path(data_folder, onto)
  dir.create(onto_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Process each gene list
  all_results <- list()
  
  for (i in seq_along(gene_files)) {
    gene_file <- gene_files[i]
    sample_name <- sub("\\.txt$", "", gene_file)
    
    message(sprintf("    Processing: %s", gene_file))
    
    # Read gene list
    gene_path <- file.path(data_folder, gene_file)
    genes <- read.table(gene_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
    
    if (length(genes) == 0) {
      warning(sprintf("Empty gene list: %s", gene_file))
      next
    }
    
    # Create gene list factor
    geneList <- factor(as.integer(geneNames %in% genes))
    names(geneList) <- geneNames
    
    # Create topGO data
    GOdata <- new("topGOdata",
                  ontology = onto,
                  allGenes = geneList,
                  nodeSize = config$params$go$node_size,
                  annot = topGO::annFUN.gene2GO,
                  gene2GO = geneID2GO
    )
    
    # Run test
    resultFis <- topGO::runTest(
      GOdata,
      algorithm = config$params$go$algorithm,
      statistic = config$params$go$statistic
    )
    
    # Get results table
    num_nodes <- length(topGO::usedGO(GOdata))
    allRes <- topGO::GenTable(
      GOdata,
      classic = resultFis,
      ranksOf = "classic",
      topNodes = num_nodes
    )
    
    # Add fold enrichment
    allRes$Fold <- allRes$Significant / allRes$Expected
    
    # Order by GO ID
    allRes <- allRes[order(allRes$GO.ID), ]
    
    # Add GO name and level
    match_idx <- match(allRes$GO.ID, go_desc$go_id)
    allRes$go_name <- go_desc$go_name[match_idx]
    allRes$go_level <- go_desc$go_level[match_idx]
    
    # Save individual result
    out_file <- file.path(onto_dir, paste0("out_", onto, "_", gene_file, ".tab"))
    write.table(
      allRes,
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # Store for merging
    all_results[[sample_name]] <- allRes
  }
  
  # Merge results
  merged <- merge_go_results(all_results, onto, go_desc)
  
  # Save merged
  merge_file <- file.path(data_folder, paste0("out_", onto, "_new_merged.tab"))
  write.table(
    merged,
    file = merge_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Create simplified merge (p-value matrix)
  simple_merge <- create_simple_merge(all_results, onto, go_desc, config)
  simple_file <- file.path(
    dirname(data_folder),
    paste0("GO_enrichment_", onto, "_merge.tab")
  )
  write.table(
    simple_merge,
    file = simple_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message(sprintf("    Saved merged results to: %s", merge_file))
  
  return(merged)
}

#' Merge GO Results from Multiple Comparisons
#' 
#' @param all_results List of results per comparison
#' @param onto Ontology name
#' @param go_desc GO description data frame
#' @return Merged data frame
#' @keywords internal
merge_go_results <- function(all_results, onto, go_desc) {
  
  if (length(all_results) == 0) {
    return(data.frame())
  }
  
  # Start with first result
  sample_names <- names(all_results)
  merged <- all_results[[1]]
  names(merged) <- paste0(sample_names[1], "_", names(merged))
  
  # Merge remaining
  if (length(all_results) > 1) {
    for (i in 2:length(all_results)) {
      res <- all_results[[i]]
      names(res) <- paste0(sample_names[i], "_", names(res))
      
      # Align by GO.ID
      merged <- merge(
        merged,
        res,
        by.x = paste0(sample_names[1], "_GO.ID"),
        by.y = paste0(sample_names[i], "_GO.ID"),
        all = TRUE
      )
    }
  }
  
  return(merged)
}

#' Create Simple P-value Merge
#' 
#' Creates a matrix of -log10(p-values) for heatmap visualization
#' 
#' @param all_results List of results per comparison
#' @param onto Ontology name
#' @param go_desc GO description data frame
#' @param config Configuration list
#' @return Data frame with GO terms as rows, comparisons as columns
#' @keywords internal
create_simple_merge <- function(all_results, onto, go_desc, config) {
  
  if (length(all_results) == 0) {
    return(data.frame())
  }
  
  # Extract p-values
  sample_names <- names(all_results)
  go_ids <- NULL
  pval_matrix <- NULL
  
  for (i in seq_along(all_results)) {
    res <- all_results[[i]]
    
    # Convert p-values
    pvals <- as.character(res$classic)
    pvals[pvals == "< 1e-30"] <- "1e-30"
    pvals <- as.numeric(pvals)
    logp <- -log10(pvals)
    
    if (is.null(go_ids)) {
      go_ids <- res$GO.ID
      pval_matrix <- logp
    } else {
      # Align
      idx <- match(go_ids, res$GO.ID)
      pval_matrix <- cbind(pval_matrix, logp[idx])
    }
  }
  
  # Create output
  colnames(pval_matrix) <- sample_names
  
  # Add annotation
  match_idx <- match(go_ids, go_desc$go_id)
  go_name <- go_desc$go_name[match_idx]
  go_level <- go_desc$go_level[match_idx]
  
  # Calculate max p-value per GO term
  max_p <- apply(pval_matrix, 1, max, na.rm = TRUE)
  
  # Get annotation from first result
  first_res <- all_results[[1]]
  annot <- first_res[match(go_ids, first_res$GO.ID), c("GO.ID", "Term", "Annotated")]
  
  final <- cbind(annot, Max_logP = max_p, pval_matrix, GO_name = go_name, Level = go_level)
  
  return(final)
}

#' Load GO Descriptions
#' 
#' @param config Configuration list
#' @return Data frame with GO term info
#' @keywords internal
load_go_descriptions <- function(config) {
  if (!file.exists(config$database$go_description_file)) {
    stop(sprintf("GO description file not found: %s", 
                 config$database$go_description_file))
  }
  
  desc <- read.table(
    config$database$go_description_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
  
  names(desc) <- c("go_id", "namespace", "category", "go_name", "go_level")
  
  return(desc)
}