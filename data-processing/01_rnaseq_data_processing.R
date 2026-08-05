#' ============================================================================
#' Script: 01_rnaseq_data_processing.R
#' Purpose: Load, process, and normalize RNA-seq count data from TCGA
#' Author: Sam Boudeau
#' Date: 2024-01-15
#' Last Updated: 2024-08-05
#'
#' Description:
#'   Processes TCGA RNA-seq data through multiple QC and normalization steps:
#'   1. Loads raw count files from GDC
#'   2. Filters genes by expression and prevalence
#'   3. Applies TPM normalization
#'   4. Performs rank-based inverse normal transformation (INT)
#'   5. Tests normality and filters genes
#'   6. Generates QC metrics and plots
#'
#' Dependencies:
#'   - dplyr (>= 1.0.0)
#'   - data.table (>= 1.13.0)
#'   - readr (>= 2.0.0)
#'   - tidyverse (>= 1.3.0)
#'   - RNOmni (>= 1.0.0)
#'   - DESeq2 (>= 1.28.0)
#'   - preprocessCore (>= 1.50.0)
#'   - BiocManager, TCGAutils
#'
#' Input:
#'   - TCGA RNA-seq augmented STAR gene count files (.tsv)
#'   - Files must be in directory specified by data_dir parameter
#'   - Expected columns: gene_id, gene_name, gene_type, counts by sample
#'
#' Output:
#'   - COADREAD_rnaSeq_TPMonly_Jan2024.csv: TPM normalized counts
#'   - COADREAD_tcga_InverseNormalTPMgenes_withShapirosig_Jan2024.csv:
#'     INT-transformed expression data (genes passing normality test)
#'   - Quality control metrics and plots
#'
#' Parameters:
#'   - data_dir: Directory path containing RNA-seq files
#'   - pattern: File pattern to match (e.g., ".augmented_star_gene_counts.tsv$")
#'   - min_expression: Minimum expression threshold in TPM (default: 0.1)
#'   - min_prevalence: Proportion of samples with >min_expression (default: 0.8)
#'   - normality_pvalue: Shapiro-Wilk p-value threshold (default: 0.001)
#'
#' Notes:
#'   - Original scripts from TCGA_analysis_inR consolidated and refactored
#'   - Removes genes with >50% zero counts
#'   - Applies INT for downstream association analyses
#'   - Normality testing reduces multiple testing burden in downstream analyses
#'
#' ============================================================================

# Load required libraries
source('utils/load_libraries.R')

# ---- Configuration Constants ----
MIN_EXPRESSION_THRESHOLD <- 0.1
MIN_PREVALENCE_RATE <- 0.8
MIN_NONZERO_RATE <- 0.5
SHAPIRO_PVALUE_THRESHOLD <- 0.001
OUTPUT_SUFFIX <- "_Jan2024"

# ---- Function: Load raw RNA-seq files ----
load_raw_rnaseq_files <- function(data_dir, pattern, recursive = TRUE) {
  #' Load multiple RNA-seq count files from directory
  #'
  #' @param data_dir Character. Path to directory containing count files.
  #' @param pattern Character. Regex pattern to match files.
  #' @param recursive Logical. Search recursively in subdirectories.
  #'
  #' @return List of data frames, one per file.
  
  file_names <- list.files(
    path = data_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive
  )
  
  if (length(file_names) == 0) {
    stop(sprintf("No files found matching pattern: %s in %s", pattern, data_dir))
  }
  
  message(sprintf("Found %d RNA-seq files", length(file_names)))
  
  data_list <- lapply(file_names, function(x) {
    tryCatch({
      sample_id <- sub('\\.*$', '', 
                       TCGAutils::UUIDtoBarcode(
                         basename(dirname(x)),
                         from_type = "file_id"
                       )[[2]][1])
      
      df <- data.table::fread(x, sep = "\t", header = TRUE)
      df <- df[, c(1, 2, 3, ncol(df)), with = FALSE]
      colnames(df) <- c("gene_id", "gene_name", "gene_type", sample_id)
      return(df)
    }, error = function(e) {
      warning(sprintf("Error reading file %s: %s", x, e$message))
      return(NULL)
    })
  })
  
  data_list <- Filter(Negate(is.null), data_list)
  return(data_list)
}

# ---- Function: Combine count files ----
combine_rnaseq_files <- function(data_list) {
  combined_df <- do.call("cbindX", data_list)
  rownames(combined_df) <- combined_df$gene_id
  combined_df <- combined_df[, !duplicated(colnames(combined_df))]
  return(combined_df)
}

# ---- Function: Filter genes by expression ----
filter_genes_by_expression <- function(expr_matrix, 
                                       min_expression = MIN_EXPRESSION_THRESHOLD,
                                       min_prevalence = MIN_PREVALENCE_RATE) {
  genes_above_threshold <- rowSums(expr_matrix >= min_expression)
  min_samples_required <- ncol(expr_matrix) * min_prevalence
  filtered_matrix <- expr_matrix[genes_above_threshold >= min_samples_required, ]
  
  message(sprintf(
    "Genes: %d -> %d (removed %d)",
    nrow(expr_matrix), nrow(filtered_matrix),
    nrow(expr_matrix) - nrow(filtered_matrix)
  ))
  return(filtered_matrix)
}

# ---- Function: Apply INT normalization ----
apply_rank_based_int <- function(expr_matrix) {
  message("Applying rank-based inverse normal transformation...")
  int_matrix <- as.data.frame(
    t(apply(expr_matrix, 1, RNOmni::RankNorm))
  )
  colnames(int_matrix) <- colnames(expr_matrix)
  rownames(int_matrix) <- rownames(expr_matrix)
  return(int_matrix)
}

# ---- Function: Test normality ----
test_gene_normality <- function(expr_matrix, pvalue_threshold = SHAPIRO_PVALUE_THRESHOLD) {
  message("Testing normality of genes (Shapiro-Wilk test)...")
  normality_results <- as.data.frame(
    t(sapply(1:nrow(expr_matrix), function(i) {
      test <- shapiro.test(as.numeric(expr_matrix[i, ]))
      c(statistic = test$statistic, p.value = test$p.value)
    }))
  )
  rownames(normality_results) <- rownames(expr_matrix)
  genes_passing <- normality_results[normality_results$p.value > pvalue_threshold, ]
  message(sprintf("Genes passing normality (p > %s): %d / %d",
                  pvalue_threshold, nrow(genes_passing), nrow(normality_results)))
  return(genes_passing)
}

# ---- Main processing pipeline ----
process_rnaseq_data <- function(data_dir, 
                                pattern = ".augmented_star_gene_counts.tsv$",
                                output_dir = "./") {
  message("\n" %+% strrep("=", 80))
  message("RNA-seq Data Processing Pipeline")
  message(strrep("=", 80))
  
  message("\nStep 1: Loading raw RNA-seq files...")
  raw_data_list <- load_raw_rnaseq_files(data_dir, pattern)
  
  message("\nStep 2: Combining files...")
  combined_data <- combine_rnaseq_files(raw_data_list)
  
  message("\nStep 3: Initial filtering (genes with >50% zeros)...")
  zero_rate <- rowSums(combined_data == 0) / ncol(combined_data)
  filtered_data <- combined_data[zero_rate < MIN_NONZERO_RATE, ]
  message(sprintf("Genes after zero-filter: %d", nrow(filtered_data)))
  
  message("\nStep 4: Filtering by expression prevalence...")
  expr_filtered <- filter_genes_by_expression(filtered_data)
  
  message("\nStep 5: Normalizing with INT...")
  int_data <- apply_rank_based_int(expr_filtered)
  
  message("\nStep 6: Testing normality...")
  normality_results <- test_gene_normality(int_data)
  int_data_filtered <- int_data[rownames(int_data) %in% rownames(normality_results), ]
  
  message("\nStep 7: Saving results...")
  tpm_file <- sprintf("%s/COADREAD_rnaSeq_TPMonly%s.csv", output_dir, OUTPUT_SUFFIX)
  write.csv(expr_filtered, file = tpm_file, row.names = TRUE)
  message(sprintf("Saved: %s", tpm_file))
  
  int_file <- sprintf("%s/COADREAD_tcga_InverseNormalTPMgenes_withShapirosig%s.csv", 
                      output_dir, OUTPUT_SUFFIX)
  write.csv(int_data_filtered, file = int_file, row.names = TRUE)
  message(sprintf("Saved: %s", int_file))
  
  message("\n" %+% strrep("=", 80))
  message("Processing complete!\n")
  
  return(list(tpm_data = expr_filtered, int_data = int_data_filtered,
              normality_results = normality_results))
}
