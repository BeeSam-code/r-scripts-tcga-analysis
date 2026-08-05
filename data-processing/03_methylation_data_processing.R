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
#'   - tidyr (>= 1.1.0)
#'   - dplyr (>= 1.0.0)
#'   - readr (>= 2.0.0)
#'   - BiocManager, TCGAutils
#'
#' Input:
#'   - TCGA methylation array files (.sesame.level3betas.txt)
#'   - HM450 reference manifest (hg38 annotation)
#'
#' Output:
#'   - COADREADaverage_methylation_beta_values_perGene_Jan2024.csv
#'
#' ============================================================================

source('utils/load_libraries.R')

# ---- Configuration Constants ----
OUTPUT_SUFFIX <- "_Jan2024"
DEFAULT_MAX_MISSING_RATE <- 0.4
DEFAULT_GENE_PREFIX <- "methyl_"

# ---- Function: Load methylation files ----
load_methylation_files <- function(data_dir, pattern, recursive = TRUE) {
  file_names <- list.files(
    path = data_dir, pattern = pattern,
    full.names = TRUE, recursive = recursive
  )
  
  if (length(file_names) == 0) {
    stop(sprintf("No methylation files found: %s in %s", pattern, data_dir))
  }
  message(sprintf("Found %d methylation files", length(file_names)))
  
  methyl_list <- lapply(file_names, function(x) {
    tryCatch({
      sample_id <- sub('\\.*$', '', 
                       TCGAutils::UUIDtoBarcode(
                         basename(dirname(x)), from_type = "file_id"
                       )[[2]][1])
      
      df <- data.frame(
        readr::read_delim(x, delim = "\t",
          col_names = c("probeID", sample_id), skip = 0),
        row.names = 1, check.names = FALSE
      )
      return(df)
    }, error = function(e) {
      warning(sprintf("Error reading %s: %s", x, e$message))
      return(NULL)
    })
  })
  
  methyl_list <- Filter(Negate(is.null), methyl_list)
  return(methyl_list)
}

# ---- Function: Combine methylation files ----
combine_methylation_files <- function(methyl_list) {
  methyl_matrix <- do.call("cbindX", methyl_list)
  rownames(methyl_matrix) <- methyl_matrix$probeID
  methyl_matrix$probeID <- NULL
  return(methyl_matrix)
}

# ---- Function: Load manifest ----
load_methylation_manifest <- function(manifest_file) {
  manifest <- readr::read_tsv(manifest_file)
  manifest_cols <- manifest[, c("probeID", "CpG_chrm", "CpG_beg", 
                                "CpG_end", "genesUniq")]
  return(manifest_cols)
}

# ---- Function: Annotate probes to genes ----
annotate_probes_to_genes <- function(methyl_matrix, manifest) {
  methyl_matrix$probeID <- rownames(methyl_matrix)
  methyl_annotated <- merge(methyl_matrix, manifest, by = "probeID")
  message(sprintf("Annotated %d probes to genes", nrow(methyl_annotated)))
  return(methyl_annotated)
}

# ---- Function: Aggregate to gene level ----
aggregate_to_gene_level <- function(methyl_annotated, max_missing_rate = DEFAULT_MAX_MISSING_RATE) {
  message("Aggregating probe-level data to gene level...")
  
  methyl_values <- methyl_annotated[, !names(methyl_annotated) %in% 
                                      c("probeID", "CpG_chrm", "CpG_beg", 
                                        "CpG_end", "genesUniq")]
  
  methyl_long <- tidyr::separate_rows(methyl_annotated, genesUniq, sep=";")
  
  methyl_gene_level <- methyl_long %>%
    dplyr::select(-c(probeID, CpG_chrm, CpG_beg, CpG_end)) %>%
    dplyr::group_by(genesUniq) %>%
    dplyr::summarise_all(list(~mean(., na.rm = TRUE)), .groups = 'drop')
  
  methyl_gene_df <- as.data.frame(methyl_gene_level)
  rownames(methyl_gene_df) <- methyl_gene_df$genesUniq
  methyl_gene_df$genesUniq <- NULL
  
  completeness <- rowMeans(!is.na(methyl_gene_df))
  methyl_gene_filtered <- methyl_gene_df[completeness > (1 - max_missing_rate), ]
  
  message(sprintf("Genes after aggregation: %d (removed %d with >%s missing)",
                  nrow(methyl_gene_filtered),
                  nrow(methyl_gene_df) - nrow(methyl_gene_filtered),
                  max_missing_rate))
  return(methyl_gene_filtered)
}

# ---- Function: Standardize sample IDs ----
standardize_sample_ids <- function(sample_ids) {
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
  message("\n" %+% strrep("=", 80))
  message("Methylation Data Processing Pipeline")
  message(strrep("=", 80))
  
  message("\nStep 1: Loading methylation files...")
  methyl_list <- load_methylation_files(data_dir, pattern)
  
  message("\nStep 2: Combining methylation files...")
  methyl_matrix <- combine_methylation_files(methyl_list)
  message(sprintf("Probe matrix: %d probes x %d samples",
                  nrow(methyl_matrix), ncol(methyl_matrix)))
  
  message("\nStep 3: Loading HM450 annotation manifest...")
  manifest <- load_methylation_manifest(manifest_file)
  
  message("\nStep 4: Annotating probes to genes...")
  methyl_annotated <- annotate_probes_to_genes(methyl_matrix, manifest)
  
  message("\nStep 5: Aggregating to gene level...")
  methyl_gene <- aggregate_to_gene_level(methyl_annotated, max_missing_rate)
  
  message("\nStep 6: Standardizing sample IDs...")
  colnames(methyl_gene) <- standardize_sample_ids(colnames(methyl_gene))
  
  message("\nStep 7: Adding gene name prefixes...")
  rownames(methyl_gene) <- paste0(DEFAULT_GENE_PREFIX, rownames(methyl_gene))
  
  message("\nStep 8: Saving results...")
  output_file <- sprintf("%s/COADREADaverage_methylation_beta_values_perGene%s.csv",
                        output_dir, OUTPUT_SUFFIX)
  write.csv(methyl_gene, file = output_file, row.names = TRUE)
  message(sprintf("Saved: %s", output_file))
  
  message("\n" %+% strrep("=", 80))
  message("Methylation processing complete!\n")
  
  return(methyl_gene)
}
