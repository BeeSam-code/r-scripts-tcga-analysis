#' ============================================================================
#' Utility: data_utilities.R
#' Purpose: Common data manipulation and integration functions
#' Author: Sam Boudeau
#' Date: 2024-01-15
#'
#' Description:
#'   Provides utility functions for data manipulation, sample matching,
#'   and covariate matrix assembly across multiple data types.
#'
#' Functions:
#'   - intersect_samples(): Find common samples across datasets
#'   - standardize_sample_ids(): Normalize TCGA sample ID formats
#'   - create_covariate_matrix(): Assemble combined covariate matrix
#'   - impute_missing_values(): Handle missing data in covariates
#'
#' ============================================================================

library(dplyr)
library(data.table)

# ---- Function: Intersect samples across datasets ----
intersect_samples <- function(...) {
  #' Find samples common to all provided datasets
  #'
  #' @param ... Vector or list of sample IDs from different datasets.
  #'
  #' @return Character vector of common sample IDs.
  #'
  #' @examples
  #'   common_samples <- intersect_samples(
  #'     colnames(rnaseq_data),
  #'     colnames(cnv_data),
  #'     colnames(methyl_data)
  #'   )
  
  sample_lists <- list(...)
  common <- sample_lists[[1]]
  
  for (i in 2:length(sample_lists)) {
    common <- intersect(common, sample_lists[[i]])
  }
  
  return(common)
}

# ---- Function: Standardize TCGA sample IDs ----
standardize_sample_ids <- function(sample_ids, format = "tcga") {
  #' Convert sample IDs to standard TCGA format
  #'
  #' @param sample_ids Character vector of sample IDs.
  #' @param format Character. Target format ("tcga" = TCGA-XX-XXXX).
  #'
  #' @return Character vector with standardized IDs.
  #'
  #' @examples
  #'   std_ids <- standardize_sample_ids(c(
  #'     "TCGA.D8.A1JU",
  #'     "TCGA-D8-A1JU",
  #'     "TCGAD8A1JU"
  #'   ))
  
  # Replace dots with dashes
  standardized <- gsub("\\.", "-", sample_ids)
  
  # Keep only first 12 characters (TCGA-XX-XXXX)
  if (format == "tcga") {
    standardized <- gsub("^([^-]*-[^-]*-[^-]*).*", "\\1", standardized)
  }
  
  return(standardized)
}

# ---- Function: Create combined covariate matrix ----
create_covariate_matrix <- function(..., common_samples = NULL) {
  #' Combine multiple covariate data frames into single matrix
  #'
  #' @param ... Data frames with covariates (rows = variables, cols = samples).
  #' @param common_samples Character vector. Samples to retain (if NULL, uses all).
  #'
  #' @return Data frame with combined covariates.
  #'
  #' @examples
  #'   covariates <- create_covariate_matrix(
  #'     clinical_data,
  #'     pca_data,
  #'     methylation_data,
  #'     common_samples = common_samples
  #'   )
  
  cov_list <- list(...)
  
  # Combine all covariates
  combined <- do.call("rbind", cov_list)
  
  # Filter to common samples if provided
  if (!is.null(common_samples)) {
    combined <- combined[, colnames(combined) %in% common_samples]
  }
  
  return(combined)
}

# ---- Function: Impute missing values ----
impute_missing_values <- function(data, method = "mean", verbose = TRUE) {
  #' Handle missing values in covariate data
  #'
  #' @param data Data frame with potential missing values.
  #' @param method Character. Imputation method ("mean", "median", "zero").
  #' @param verbose Logical. Print imputation summary.
  #'
  #' @return Data frame with imputed values.
  #'
  #' @examples
  #'   covariates_imputed <- impute_missing_values(
  #'     covariate_data,
  #'     method = "mean"
  #'   )
  
  data_imputed <- data
  
  for (col in colnames(data_imputed)) {
    missing_count <- sum(is.na(data_imputed[[col]]))
    
    if (missing_count > 0) {
      if (method == "mean") {
        data_imputed[[col]][is.na(data_imputed[[col]])] <-
          mean(data_imputed[[col]], na.rm = TRUE)
      } else if (method == "median") {
        data_imputed[[col]][is.na(data_imputed[[col]])] <-
          median(data_imputed[[col]], na.rm = TRUE)
      } else if (method == "zero") {
        data_imputed[[col]][is.na(data_imputed[[col]])] <- 0
      }
      
      if (verbose) {
        message(sprintf("%s: imputed %d missing values using %s",
                       col, missing_count, method))
      }
    }
  }
  
  return(data_imputed)
}
