#' ============================================================================
#' Script: 03_methylation_data_processing.R
#' Purpose: Process methylation array data and aggregate to gene level
#' Author: Sam Boudeau
#' Date: 2024-01-15
#' Last Updated: 2024-08-05
#'
#' Description:
#'   Processes TCGA methylation array data through multiple steps:
#'   1. Loads probe-level methylation beta values
#'   2. Annotates probes to genes using HM450 manifest
#'   3. Aggregates probe-level data to gene level (average across probes)
#'   4. Filters genes and probes by data completeness
#'   5. Normalizes sample IDs
#'   6. Exports gene-level methylation matrix
#'
#' Dependencies:
#'   - data.table (>= 1.13.0)
#'   - dplyr (>= 1.0.0)
#'   - readr (>= 2.0.0)
#'   - tidyr (>= 1.1.0)
#'   - BiocManager, TCGAutils
#'
#' Input:
#'   - TCGA methylation array files (.sesame.level3betas.txt)
#'   - HM450 reference manifest (hg38 annotation)
#'     File: HM450.hg38.manifest.gencode.v36.tsv
#'     Available from GDC reference files
#'
#' Output:
#'   - COADREADaverage_methylation_beta_values_perGene_Jan2024.csv:
#'     Gene-level methylation matrix with beta values averaged across probes
#'
#' Parameters:
#'   - data_dir: Directory path containing methylation files
#'   - pattern: File pattern to match methylation files
#'   - manifest_file: Path to HM450 annotation manifest
#'   - output_dir: Directory for output files
#'   - max_missing_rate: Max proportion of missing data (default: 0.4)
#'
#' Notes:
#'   - Probe-level values are beta values (0-1)
#'   - Multiple probes per gene are averaged
#'   - Genes with >40% missing data are filtered out
#'   - Gene names prefixed with 'methyl_' to avoid column name conflicts
#'
#' ============================================================================

# Load required libraries
source('utils/load_libraries.R')

# ---- Configuration Constants ----
OUTPUT_SUFFIX <- "_Jan2024"
DEFAULT_MAX_MISSING_RATE <- 0.4
DEFAULT_GENE_PREFIX <- "methyl_"

# ---- Function: Load methylation files ----
load_methylation_files <- function(data_dir, pattern, recursive = TRUE) {
  #' Load multiple methylation files from directory
  #'
  #' @param data_dir Character. Path to directory containing methylation files.
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
    stop(sprintf("No methylation files found matching pattern: %s in %s",
                 pattern, data_dir))
  }
  
  message(sprintf("Found %d methylation files", length(file_names)))
  
  # Load each file
  methyl_list <- lapply(file_names, function(x) {
    tryCatch({
      # Extract sample ID from file path
      sample_id <- sub('\\.*$', '', 
                       TCGAutils::UUIDtoBarcode(
                         basename(dirname(x)),
                         from_type = "file_id"
                       )[[2]][1])
      
      # Read methylation file with probe IDs as row names
      df <- data.frame(
        readr::read_delim(
          x,
          delim = "\t",
          col_names = c("probeID", sample_id),
          skip = 0
        ),
        row.names = 1,
        check.names = FALSE
      )
      
      return(df)
    }, error = function(e) {
      warning(sprintf("Error reading methylation file %s: %s",
                      x, e$message))
      return(NULL)
    })
  })
  
  methyl_list <- Filter(Negate(is.null), methyl_list)
  return(methyl_list)
}

# ---- Function: Combine methylation files ----
combine_methylation_files <- function(methyl_list) {
  #' Combine multiple methylation files into single matrix
  #'
  #' @param methyl_list List of data frames from load_methylation_files()
  #'
  #' @return Data frame with probes as rows, samples as columns
  
  methyl_matrix <- do.call("cbindX", methyl_list)
  rownames(methyl_matrix) <- methyl_matrix$probeID
  methyl_matrix$probeID <- NULL
  
  return(methyl_matrix)
}

# ---- Function: Load and prepare manifest ----
load_methylation_manifest <- function(manifest_file) {
  #' Load HM450 annotation manifest
  #'
  #' @param manifest_file Character. Path to manifest file.
  #'
  #' @return Data frame with probe annotations
  
  manifest <- readr::read_tsv(manifest_file)
  
  # Keep only relevant columns
  manifest_cols <- manifest[, c("probeID", "CpG_chrm", "CpG_beg", 
                                "CpG_end", "genesUniq")]
  
  return(manifest_cols)
}

# ---- Function: Annotate probes to genes ----
annotate_probes_to_genes <- function(methyl_matrix, manifest) {
  #' Annotate methylation probes to genes using manifest
  #'
  #' @param methyl_matrix Data frame with probes as rows.
  #' @param manifest Data frame with probe annotations.
  #'
  #' @return Data frame with gene annotations added
  
  methyl_matrix$probeID <- rownames(methyl_matrix)
  
  # Merge with manifest
  methyl_annotated <- merge(methyl_matrix, manifest, by = "probeID")
  
  message(sprintf("Annotated %d probes to genes", nrow(methyl_annotated)))
  
  return(methyl_annotated)
}

# ---- Function: Aggregate to gene level ----
aggregate_to_gene_level <- function(methyl_annotated, max_missing_rate = DEFAULT_MAX_MISSING_RATE) {
  #' Aggregate probe-level methylation to gene level
  #'
  #' @param methyl_annotated Data frame with annotated methylation data.
  #' @param max_missing_rate Numeric. Max proportion of missing data.
  #'
  #' @return Data frame with genes as rows
  
  message("Aggregating probe-level data to gene level...")
  
  # Remove non-numeric columns for averaging
  methyl_values <- methyl_annotated[, !names(methyl_annotated) %in% 
                                      c("probeID", "CpG_chrm", "CpG_beg", 
                                        "CpG_end", "genesUniq")]
  
  # Handle genes with multiple probe sets (separated by ";")
  methyl_long <- tidyr::separate_rows(methyl_annotated, genesUniq, sep=";")
  
  # Aggregate by gene: average across all probes
  methyl_gene_level <- methyl_long %>%
    dplyr::select(-c(probeID, CpG_chrm, CpG_beg, CpG_end)) %>%
    dplyr::group_by(genesUniq) %>%
    dplyr::summarise_all(list(~mean(., na.rm = TRUE)), .groups = 'drop')
  
  # Convert back to data frame
  methyl_gene_df <- as.data.frame(methyl_gene_level)
  rownames(methyl_gene_df) <- methyl_gene_df$genesUniq
  methyl_gene_df$genesUniq <- NULL
  
  # Filter genes by data completeness
  completeness <- rowMeans(!is.na(methyl_gene_df))
  methyl_gene_filtered <- methyl_gene_df[completeness > (1 - max_missing_rate), ]
  
  message(sprintf(
    "Genes after aggregation and filtering: %d (removed %d with >%s missing)",
    nrow(methyl_gene_filtered),
    nrow(methyl_gene_df) - nrow(methyl_gene_filtered),
    max_missing_rate
  ))
  
  return(methyl_gene_filtered)
}

# ---- Function: Standardize sample IDs ----
standardize_sample_ids <- function(sample_ids) {
  #' Convert sample IDs to standard TCGA format
  #'
  #' @param sample_ids Character vector of sample IDs.
  #'
  #' @return Character vector with standardized IDs
  
  # Replace dots with dashes and keep only first 12 characters
  standardized <- gsub("\\.", "-", sample_ids)
  standardized <- gsub("^([^-]*-[^-]*-[^-]*).*", "\\1", standardized)
  
  return(standardized)
}

# ---- Main processing pipeline ----
process_methylation_data <- function(data_dir,
                                     pattern = ".methylation_array.sesame.level3betas.txt$",
                                     manifest_file = "HM450.hg38.manifest.gencode.v36.tsv",
                                     output_dir = "./",
                                     max_missing_rate = DEFAULT_MAX_MISSING_RATE) {
  #' Complete methylation processing pipeline
  #'
  #' @param data_dir Character. Directory containing methylation files.
  #' @param pattern Character. File pattern to match.
  #' @param manifest_file Character. Path to HM450 manifest.
  #' @param output_dir Character. Directory for output files.
  #' @param max_missing_rate Numeric. Max proportion of missing data.
  #'
  #' @return Data frame with processed methylation data
  
  message("\n" %+% strrep("=", 80))
  message("Methylation Data Processing Pipeline")
  message(strrep("=", 80))
  
  # Step 1: Load methylation files
  message("\nStep 1: Loading methylation files...")
  methyl_list <- load_methylation_files(data_dir, pattern)
  
  # Step 2: Combine into matrix
  message("\nStep 2: Combining methylation files...")
  methyl_matrix <- combine_methylation_files(methyl_list)
  message(sprintf("Probe matrix: %d probes x %d samples",
                  nrow(methyl_matrix), ncol(methyl_matrix)))
  
  # Step 3: Load manifest
  message("\nStep 3: Loading HM450 annotation manifest...")
  manifest <- load_methylation_manifest(manifest_file)
  
  # Step 4: Annotate probes
  message("\nStep 4: Annotating probes to genes...")
  methyl_annotated <- annotate_probes_to_genes(methyl_matrix, manifest)
  
  # Step 5: Aggregate to gene level
  message("\nStep 5: Aggregating to gene level...")
  methyl_gene <- aggregate_to_gene_level(methyl_annotated, max_missing_rate)
  
  # Step 6: Standardize sample IDs
  message("\nStep 6: Standardizing sample IDs...")
  colnames(methyl_gene) <- standardize_sample_ids(colnames(methyl_gene))
  
  # Step 7: Add gene name prefix
  message("\nStep 7: Adding gene name prefixes...")
  rownames(methyl_gene) <- paste0(DEFAULT_GENE_PREFIX, rownames(methyl_gene))
  
  # Step 8: Save output
  message("\nStep 8: Saving results...")
  output_file <- sprintf("%s/COADREADaverage_methylation_beta_values_perGene%s.csv",
                        output_dir, OUTPUT_SUFFIX)
  write.csv(methyl_gene, file = output_file, row.names = TRUE)
  message(sprintf("Saved methylation data: %s", output_file))
  
  message("\n" %+% strrep("=", 80))
  message("Methylation processing complete!")
  message(strrep("=", 80) %+% "\n")
  
  return(methyl_gene)
}

# ---- Execute if script is run directly ----
if (!interactive()) {
  # Example usage
  methyl_data <- process_methylation_data(
    data_dir = "./data/methylation/",
    pattern = ".methylation_array.sesame.level3betas.txt$",
    manifest_file = "./data/HM450.hg38.manifest.gencode.v36.tsv",
    output_dir = "./results/"
  )
}
