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
#'   - readr (>= 2.0.0)
#'   - BiocManager, TCGAutils
#'
#' Input:
#'   - TCGA gene-level CNV files (.v36.tsv format)
#'   - Expected columns: gene_id, gene_name, chromosome, start, end, 
#'     copy_number, min_copy_number, max_copy_number
#'
#' Output:
#'   - gene_level_CNV_COADREAD_Jan2024.csv: Gene-level CNV matrix
#'     (genes as rows, samples as columns with sample IDs as headers)
#'
#' Parameters:
#'   - data_dir: Directory path containing CNV files
#'   - pattern: File pattern to match CNV files
#'   - output_dir: Directory for output files
#'
#' Notes:
#'   - Removes redundant gene name columns, keeping only first occurrence
#'   - Handles UUID to TCGA barcode conversion automatically
#'   - Duplicate files (.R vs .r) consolidated
#'
#' ============================================================================

# Load required libraries
source('utils/load_libraries.R')

# ---- Configuration Constants ----
OUTPUT_SUFFIX <- "_Jan2024"
GENE_ID_COLUMN_NAMES <- c("gene_id", "gene_name", "chromosome", 
                          "start", "end", "copy_number", 
                          "min_copy_number", "max_copy_number")

# ---- Function: Load single CNV file ----
load_single_cnv_file <- function(file_path) {
  #' Load and process single CNV file
  #'
  #' @param file_path Character. Full path to CNV file.
  #'
  #' @return Data frame with gene info and sample copy number
  
  tryCatch({
    # Extract sample ID from file path
    sample_id <- sub('\\.*$', '', 
                     TCGAutils::UUIDtoBarcode(
                       basename(dirname(file_path)),
                       from_type = "file_id"
                     )[[2]][1])
    
    # Read CNV file
    df <- data.table::fread(
      file_path,
      sep = "\t",
      header = TRUE,
      select = c("gene_id", "gene_name", "chromosome", "start", 
                 "end", "copy_number", "min_copy_number", "max_copy_number")
    )
    
    # Keep only gene info and sample-specific copy number
    df <- df[, c(1, 5, ncol(df)), with = FALSE]  # gene_id, end, copy_number
    colnames(df) <- c("gene_id", "end", sample_id)
    
    return(df)
  }, error = function(e) {
    warning(sprintf("Error reading CNV file %s: %s", file_path, e$message))
    return(NULL)
  })
}

# ---- Function: Load all CNV files ----
load_cnv_files <- function(data_dir, pattern, recursive = TRUE) {
  #' Load multiple CNV files from directory
  #'
  #' @param data_dir Character. Path to directory containing CNV files.
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
    stop(sprintf("No CNV files found matching pattern: %s in %s", pattern, data_dir))
  }
  
  message(sprintf("Found %d CNV files", length(file_names)))
  
  # Load each file
  cnv_list <- lapply(file_names, load_single_cnv_file)
  cnv_list <- Filter(Negate(is.null), cnv_list)
  
  return(cnv_list)
}

# ---- Function: Aggregate CNV data ----
aggregate_cnv_data <- function(cnv_list) {
  #' Combine multiple CNV files into single matrix
  #'
  #' @param cnv_list List of data frames from load_cnv_files()
  #'
  #' @return Data frame with genes as rows, samples as columns
  
  # Use data.table cbindX for efficient joining
  cnv_matrix <- do.call("cbindX", cnv_list)
  
  # Set row names to gene IDs
  rownames(cnv_matrix) <- cnv_matrix$gene_id
  
  # Remove gene metadata columns
  cnv_matrix <- cnv_matrix[, !colnames(cnv_matrix) %in% c("gene_id", "end")]
  
  # Remove duplicate columns
  cnv_matrix <- cnv_matrix[, !duplicated(colnames(cnv_matrix))]
  
  return(cnv_matrix)
}

# ---- Function: Standardize sample IDs ----
standardize_sample_ids <- function(sample_ids) {
  #' Convert sample IDs to standard TCGA format (TCGA-XX-XXXX)
  #'
  #' @param sample_ids Character vector of sample IDs.
  #'
  #' @return Character vector with standardized IDs
  
  # Keep only first 12 characters (TCGA-XX-XXXX)
  standardized <- gsub("^([^-]*-[^-]*-[^-]*).*", "\\1", sample_ids)
  
  return(standardized)
}

# ---- Main processing pipeline ----
process_cnv_data <- function(data_dir,
                             pattern = ".v36.tsv$",
                             output_dir = "./") {
  #' Complete CNV processing pipeline
  #'
  #' @param data_dir Character. Directory containing CNV files.
  #' @param pattern Character. File pattern to match.
  #' @param output_dir Character. Directory for output files.
  #'
  #' @return Data frame with processed CNV data
  
  message("\n" %+% strrep("=", 80))
  message("CNV Data Processing Pipeline")
  message(strrep("=", 80))
  
  # Step 1: Load CNV files
  message("\nStep 1: Loading CNV files...")
  cnv_list <- load_cnv_files(data_dir, pattern)
  
  # Step 2: Aggregate into matrix
  message("\nStep 2: Aggregating CNV data...")
  cnv_matrix <- aggregate_cnv_data(cnv_list)
  message(sprintf("CNV matrix dimensions: %d genes x %d samples",
                  nrow(cnv_matrix), ncol(cnv_matrix)))
  
  # Step 3: Standardize sample IDs
  message("\nStep 3: Standardizing sample IDs...")
  colnames(cnv_matrix) <- standardize_sample_ids(colnames(cnv_matrix))
  
  # Step 4: Quality control - filter by data completeness
  message("\nStep 4: Quality control...")
  initial_genes <- nrow(cnv_matrix)
  initial_samples <- ncol(cnv_matrix)
  
  # Keep genes with data in >50% of samples
  gene_completeness <- rowMeans(!is.na(cnv_matrix))
  cnv_matrix_qc <- cnv_matrix[gene_completeness > 0.5, ]
  
  # Keep samples with data in >50% of genes
  sample_completeness <- colMeans(!is.na(cnv_matrix_qc))
  cnv_matrix_qc <- cnv_matrix_qc[, sample_completeness > 0.5]
  
  message(sprintf(
    "After QC: %d genes (removed %d), %d samples (removed %d)",
    nrow(cnv_matrix_qc),
    initial_genes - nrow(cnv_matrix_qc),
    ncol(cnv_matrix_qc),
    initial_samples - ncol(cnv_matrix_qc)
  ))
  
  # Step 5: Save output
  message("\nStep 5: Saving results...")
  output_file <- sprintf("%s/gene_level_CNV_COADREAD%s.csv",
                        output_dir, OUTPUT_SUFFIX)
  write.csv(cnv_matrix_qc, file = output_file, row.names = TRUE)
  message(sprintf("Saved CNV data: %s", output_file))
  
  message("\n" %+% strrep("=", 80))
  message("CNV processing complete!")
  message(strrep("=", 80) %+% "\n")
  
  return(cnv_matrix_qc)
}

# ---- Execute if script is run directly ----
if (!interactive()) {
  # Example usage
  cnv_data <- process_cnv_data(
    data_dir = "./data/cnv/",
    pattern = ".v36.tsv$",
    output_dir = "./results/"
  )
}
