#' ============================================================================
#' Script: 02_cnv_data_processing.R
#' Purpose: Load and process gene-level copy number variation data
#' Author: Sam Boudeau
#' Date: 2024-01-15
#' Last Updated: 2024-08-05
#'
#' Description:
#'   Processes TCGA gene-level CNV data through multiple steps:
#'   1. Loads CNV files from GDC for multiple samples
#'   2. Aggregates into single matrix
#'   3. Removes duplicate gene name columns
#'   4. Standardizes sample IDs to TCGA format
#'   5. Exports consolidated matrix
#'
#' Dependencies:
#'   - data.table (>= 1.13.0)
#'   - dplyr (>= 1.0.0)
#'   - BiocManager, TCGAutils
#'
#' Input:
#'   - TCGA gene-level CNV files (.v36.tsv format)
#'
#' Output:
#'   - gene_level_CNV_COADREAD_Jan2024.csv: Gene-level CNV matrix
#'
#' ============================================================================

source('utils/load_libraries.R')

# ---- Configuration Constants ----
OUTPUT_SUFFIX <- "_Jan2024"

# ---- Function: Load single CNV file ----
load_single_cnv_file <- function(file_path) {
  tryCatch({
    sample_id <- sub('\\.*$', '', 
                     TCGAutils::UUIDtoBarcode(
                       basename(dirname(file_path)),
                       from_type = "file_id"
                     )[[2]][1])
    
    df <- data.table::fread(
      file_path, sep = "\t", header = TRUE,
      select = c("gene_id", "gene_name", "chromosome", "start", 
                 "end", "copy_number", "min_copy_number", "max_copy_number")
    )
    df <- df[, c(1, 5, ncol(df)), with = FALSE]
    colnames(df) <- c("gene_id", "end", sample_id)
    return(df)
  }, error = function(e) {
    warning(sprintf("Error reading CNV file %s: %s", file_path, e$message))
    return(NULL)
  })
}

# ---- Function: Load all CNV files ----
load_cnv_files <- function(data_dir, pattern, recursive = TRUE) {
  file_names <- list.files(
    path = data_dir, pattern = pattern,
    full.names = TRUE, recursive = recursive
  )
  
  if (length(file_names) == 0) {
    stop(sprintf("No CNV files found matching pattern: %s in %s", pattern, data_dir))
  }
  message(sprintf("Found %d CNV files", length(file_names)))
  
  cnv_list <- lapply(file_names, load_single_cnv_file)
  cnv_list <- Filter(Negate(is.null), cnv_list)
  return(cnv_list)
}

# ---- Function: Aggregate CNV data ----
aggregate_cnv_data <- function(cnv_list) {
  cnv_matrix <- do.call("cbindX", cnv_list)
  rownames(cnv_matrix) <- cnv_matrix$gene_id
  cnv_matrix <- cnv_matrix[, !colnames(cnv_matrix) %in% c("gene_id", "end")]
  cnv_matrix <- cnv_matrix[, !duplicated(colnames(cnv_matrix))]
  return(cnv_matrix)
}

# ---- Function: Standardize sample IDs ----
standardize_sample_ids <- function(sample_ids) {
  standardized <- gsub("^([^-]*-[^-]*-[^-]*).*", "\\1", sample_ids)
  return(standardized)
}

# ---- Main processing pipeline ----
process_cnv_data <- function(data_dir, pattern = ".v36.tsv$", output_dir = "./") {
  message("\n" %+% strrep("=", 80))
  message("CNV Data Processing Pipeline")
  message(strrep("=", 80))
  
  message("\nStep 1: Loading CNV files...")
  cnv_list <- load_cnv_files(data_dir, pattern)
  
  message("\nStep 2: Aggregating CNV data...")
  cnv_matrix <- aggregate_cnv_data(cnv_list)
  message(sprintf("CNV matrix: %d genes x %d samples", nrow(cnv_matrix), ncol(cnv_matrix)))
  
  message("\nStep 3: Standardizing sample IDs...")
  colnames(cnv_matrix) <- standardize_sample_ids(colnames(cnv_matrix))
  
  message("\nStep 4: Quality control...")
  initial_genes <- nrow(cnv_matrix)
  initial_samples <- ncol(cnv_matrix)
  
  gene_completeness <- rowMeans(!is.na(cnv_matrix))
  cnv_matrix_qc <- cnv_matrix[gene_completeness > 0.5, ]
  
  sample_completeness <- colMeans(!is.na(cnv_matrix_qc))
  cnv_matrix_qc <- cnv_matrix_qc[, sample_completeness > 0.5]
  
  message(sprintf("After QC: %d genes (removed %d), %d samples (removed %d)",
                  nrow(cnv_matrix_qc), initial_genes - nrow(cnv_matrix_qc),
                  ncol(cnv_matrix_qc), initial_samples - ncol(cnv_matrix_qc)))
  
  message("\nStep 5: Saving results...")
  output_file <- sprintf("%s/gene_level_CNV_COADREAD%s.csv", output_dir, OUTPUT_SUFFIX)
  write.csv(cnv_matrix_qc, file = output_file, row.names = TRUE)
  message(sprintf("Saved: %s", output_file))
  
  message("\n" %+% strrep("=", 80))
  message("CNV processing complete!\n")
  
  return(cnv_matrix_qc)
}
