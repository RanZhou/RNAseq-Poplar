# qc.R - Quality Control Functions
# Sample clustering, PCA, and correlation analysis

#' Run Sample QC Pipeline
#' 
#' Performs PCA, hierarchical clustering, and Pearson correlation
#' on FPKM data for quality control
#' 
#' @param config Configuration list from load_config()
#' @param fpkm_data Optional pre-loaded FPKM matrix (will load from file if NULL)
#' @return List containing PCA results, clustering, and correlation matrix
#' @export
run_sample_qc <- function(config, fpkm_data = NULL) {
  
  message("Running Sample QC...")
  
  # Load data
  design <- load_design(config)
  
  if (is.null(fpkm_data)) {
    fpkm_data <- load_fpkm(config)
  }
  
  # Validate samples match
  samples <- as.character(design$File)
  missing <- samples[!samples %in% colnames(fpkm_data)]
  if (length(missing) > 0) {
    stop(sprintf("Samples in design not found in FPKM data: %s", 
                 paste(missing, collapse = ", ")))
  }
  
  fpkm_test <- fpkm_data[, samples, drop = FALSE]
  
  # Filter genes by expression
  groups <- as.character(design$Group)
  fpkm_cut_tag <- get_cut_tag(fpkm_test, groups)
  expressed_genes <- fpkm_cut_tag >= config$params$pca_min_fpkm
  
  message(sprintf("  Genes passing FPKM >= %d filter: %d / %d", 
                  config$params$pca_min_fpkm, 
                  sum(expressed_genes), 
                  nrow(fpkm_test)))
  
  fpkm_filtered <- fpkm_test[expressed_genes, , drop = FALSE]
  
  # Log transform
  data_ma <- as.matrix(log(fpkm_filtered + 1))
  
  # PCA
  message("  Running PCA...")
  pca_results <- run_pca(data_ma, design)
  
  # Hierarchical clustering
  message("  Running hierarchical clustering...")
  hc_results <- run_clustering(data_ma, design)
  
  # Pearson correlation
  message("  Calculating Pearson correlations...")
  pcc_results <- run_pearson_correlation(fpkm_filtered, config)
  
  # Compile results
  results <- list(
    pca = pca_results,
    clustering = hc_results,
    correlation = pcc_results,
    samples_used = samples,
    genes_used = sum(expressed_genes)
  )
  
  message("Sample QC complete.")
  return(results)
}

#' Run PCA Analysis
#' 
#' @param data_ma Log-transformed FPKM matrix
#' @param design Design data frame
#' @return List with PCA loadings and plot data
#' @keywords internal
run_pca <- function(data_ma, design) {
  # Calculate PCA
  pcdat <- princomp(data_ma)
  loadings <- pcdat$loadings[, c(1, 2)]
  
  # Get colors by group
  groups <- as.character(design$Group)
  unique_groups <- unique(groups)
  palette(rainbow(length(unique_groups)))
  colors <- palette()
  
  group_colors <- colors[match(groups, unique_groups)]
  
  # Calculate center lines
  h_cut <- round(mean(loadings[, 2]), 3)
  v_cut <- round(mean(loadings[, 1]), 3)
  
  # Prepare output
  results <- list(
    loadings = loadings,
    sample_names = as.character(design$sample_id),
    group_colors = group_colors,
    groups = groups,
    center = c(v_cut, h_cut),
    variance_explained = summary(pcdat)$importance[2, 1:2]
  )
  
  return(results)
}

#' Run Hierarchical Clustering
#' 
#' @param data_ma Log-transformed FPKM matrix
#' @param design Design data frame
#' @return hclust object
#' @keywords internal
run_clustering <- function(data_ma, design) {
  # Transpose for sample clustering
  dist_mat <- dist(t(data_ma))
  hc <- hclust(dist_mat)
  
  return(hc)
}

#' Calculate Pearson Correlation Matrix
#' 
#' @param fpkm_data FPKM matrix (filtered)
#' @param config Configuration list
#' @return Correlation matrix
#' @keywords internal
run_pearson_correlation <- function(fpkm_data, config) {
  samples <- colnames(fpkm_data)
  
  # Calculate correlations
  corr_mat <- cor(fpkm_data, method = "pearson")
  
  # Save to file
  pcc_file <- file.path(
    config$output$results_dir,
    paste0(config$output$pcc_prefix, basename(config$input$counts_file), ".csv")
  )
  
  write.csv(corr_mat, file = pcc_file, row.names = TRUE)
  message(sprintf("  Saved correlation matrix to: %s", pcc_file))
  
  return(corr_mat)
}

#' Get Expression Cut Tag
#' 
#' Assigns FPKM expression level tags (0, 1, 3, 5, 10) per gene
#' based on maximum expression across groups
#' 
#' @param fpkm_data FPKM matrix
#' @param groups Vector of group assignments
#' @return Vector of cut tags
#' @export
get_cut_tag <- function(fpkm_data, groups) {
  unique_groups <- unique(groups)
  
  # Calculate min FPKM per group
  group_mins <- lapply(unique_groups, function(g) {
    samples <- names(groups)[groups == g]
    if (length(samples) == 1) {
      return(fpkm_data[, samples])
    } else {
      return(apply(fpkm_data[, samples, drop = FALSE], 1, min))
    }
  })
  
  # Get max across groups
  min_matrix <- do.call(cbind, group_mins)
  colnames(min_matrix) <- unique_groups
  max_vals <- apply(min_matrix, 1, max)
  
  # Assign tags
  tags <- rep(0, length(max_vals))
  tags[max_vals >= 1] <- 1
  tags[max_vals >= 3] <- 3
  tags[max_vals >= 5] <- 5
  tags[max_vals >= 10] <- 10
  
  names(tags) <- rownames(fpkm_data)
  return(tags)
}

#' Plot PCA Results
#' 
#' @param pca_results List from run_pca()
#' @param title Plot title
#' @return ggplot object
#' @export
plot_pca <- function(pca_results, title = "PCA Analysis") {
  if (!require("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  
  df <- data.frame(
    PC1 = pca_results$loadings[, 1],
    PC2 = pca_results$loadings[, 2],
    Sample = pca_results$sample_names,
    Group = pca_results$groups,
    Color = pca_results$group_colors
  )
  
  var1 <- round(pca_results$variance_explained[1] * 100, 1)
  var2 <- round(pca_results$variance_explained[2] * 100, 1)
  
  p <- ggplot(df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    geom_vline(xintercept = pca_results$center[1], linetype = "dashed", color = "red") +
    geom_hline(yintercept = pca_results$center[2], linetype = "dashed", color = "red") +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", var1),
      y = sprintf("PC2 (%.1f%%)", var2)
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}

#' Plot Sample Clustering
#' 
#' @param hc_results hclust object from run_clustering()
#' @param title Plot title
#' @return ggplot object
#' @export
plot_clustering <- function(hc_results, title = "Sample Clustering") {
  if (!require("ggdendro", quietly = TRUE)) {
    stop("Package 'ggdendro' required for plotting")
  }
  
  p <- ggdendrogram(hc_results, rotate = FALSE, theme_dendro = FALSE) +
    labs(title = title) +
    theme_minimal()
  
  return(p)
}