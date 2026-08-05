#' ============================================================================
#' Script: 02_eqtl_identification.R
#' Purpose: Identify expression quantitative trait loci (eQTLs)
#' Author: Sam Boudeau
#' Date: 2024-01-15
#' Last Updated: 2024-08-05
#'
#' Description:
#'   Performs eQTL analysis using MatrixEQTL through multiple steps:
#'   1. Formats genotype, expression, and covariate data
#'   2. Tests SNP-gene associations (cis and trans)
#'   3. Applies multiple testing corrections
#'   4. Annotates significant associations
#'   5. Generates results and visualization plots
#'
#' Dependencies:
#'   - MatrixEQTL (>= 2.1.0)
#'   - data.table (>= 1.13.0)
#'   - dplyr (>= 1.0.0)
#'   - readr (>= 2.0.0)
#'   - ggplot2 (>= 3.3.0)
#'
#' Input:
#'   - Genotype matrix (SNPs x samples, biallelic coding)
#'   - Expression matrix (genes x samples, normalized/INT-transformed)
#'   - Covariate matrix (age, gender, ancestry PCs, methylation, CNV)
#'   - SNP and gene position files
#'
#' Output:
#'   - Cis-eQTL results (SNPs within 1Mb of gene)
#'   - Trans-eQTL results (genome-wide associations)
#'   - Manhattan and Q-Q plots
#'   - Summary statistics
#'
#' Notes:
#'   - Original script from TCGA_analysis_inR refactored
#'   - Covariates: age, gender, genetic ancestry (PCs), methylation, CNV
#'   - Cis-distance threshold: 1Mb (1e6 bp)
#'   - Uses linear model (modelLINEAR) for association testing
#'
#' ============================================================================

source('utils/load_libraries.R')

# ---- Configuration Constants ----
CIS_DISTANCE <- 1e6  # 1 Megabase
PVALUE_THRESHOLD_CIS <- 1e-2
PVALUE_THRESHOLD_TRANS <- 1e-2
USE_MODEL <- "modelLINEAR"

# ---- Function: Load and prepare matrices ----
prepare_eqtl_matrices <- function(genotype_data, expression_data, 
                                  covariate_data) {
  #' Format data into SlicedData objects for MatrixEQTL
  #'
  #' @param genotype_data Data frame. Genotypes (SNPs x samples).
  #' @param expression_data Data frame. Expression (genes x samples).
  #' @param covariate_data Data frame. Covariates (covariates x samples).
  #'
  #' @return List with SlicedData objects.
  
  message("Preparing eQTL data matrices...")
  
  # Create SlicedData for SNPs
  snps_mat <- MatrixEQTL::SlicedData$new()
  snps_mat$CreateFromMatrix(as.matrix(genotype_data[, -1]))
  snps_mat$fileOmitCharacters <- "NA"
  snps_mat$fileSkipColumns <- 1
  snps_mat$fileSkipRows <- 1
  snps_mat$fileSliceSize <- 2000
  message(sprintf("SNPs matrix: %d SNPs x %d samples",
                  nrow(genotype_data), ncol(genotype_data) - 1))
  
  # Create SlicedData for genes
  genes_mat <- MatrixEQTL::SlicedData$new()
  genes_mat$CreateFromMatrix(as.matrix(expression_data[, -1]))
  genes_mat$fileOmitCharacters <- "NA"
  genes_mat$fileSkipColumns <- 1
  genes_mat$fileSkipRows <- 1
  genes_mat$fileSliceSize <- 2000
  message(sprintf("Genes matrix: %d genes x %d samples",
                  nrow(expression_data), ncol(expression_data) - 1))
  
  # Create SlicedData for covariates
  covar_mat <- MatrixEQTL::SlicedData$new()
  covar_mat$CreateFromMatrix(as.matrix(covariate_data[, -1]))
  covar_mat$fileOmitCharacters <- "NA"
  covar_mat$fileSkipColumns <- 1
  covar_mat$fileSkipRows <- 1
  covar_mat$fileSliceSize <- 2000
  message(sprintf("Covariates matrix: %d covariates x %d samples",
                  nrow(covariate_data), ncol(covariate_data) - 1))
  
  return(list(
    snps = snps_mat,
    genes = genes_mat,
    covariates = covar_mat
  ))
}

# ---- Function: Load position files ----
load_position_files <- function(snp_pos_file, gene_pos_file) {
  #' Load SNP and gene position information
  #'
  #' @param snp_pos_file Character. Path to SNP position file.
  #' @param gene_pos_file Character. Path to gene position file.
  #'
  #' @return List with SNP and gene position data frames.
  
  message("Loading position files...")
  
  snp_pos <- readr::read_tsv(snp_pos_file)
  colnames(snp_pos) <- c("snpid", "chr", "pos")
  message(sprintf("SNP positions: %d SNPs", nrow(snp_pos)))
  
  gene_pos <- readr::read_tsv(gene_pos_file)
  colnames(gene_pos) <- c("geneid", "chr", "left", "right")
  message(sprintf("Gene positions: %d genes", nrow(gene_pos)))
  
  return(list(snp_pos = as.data.frame(snp_pos),
              gene_pos = as.data.frame(gene_pos)))
}

# ---- Function: Run eQTL analysis ----
run_eqtl_analysis <- function(sliced_data, position_data, 
                              output_file_cis, output_file_trans) {
  #' Perform eQTL analysis using MatrixEQTL
  #'
  #' @param sliced_data List from prepare_eqtl_matrices().
  #' @param position_data List from load_position_files().
  #' @param output_file_cis Character. Output file for cis-eQTLs.
  #' @param output_file_trans Character. Output file for trans-eQTLs.
  #'
  #' @return MatrixEQTL results object.
  
  message("Running eQTL analysis...")
  
  # Run analysis
  me <- MatrixEQTL::Matrix_eQTL_main(
    snps = sliced_data$snps,
    gene = sliced_data$genes,
    cvrt = sliced_data$covariates,
    output_file_name = output_file_trans,
    output_file_name.cis = output_file_cis,
    pvOutputThreshold = PVALUE_THRESHOLD_TRANS,
    pvOutputThreshold.cis = PVALUE_THRESHOLD_CIS,
    useModel = USE_MODEL,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    snpspos = position_data$snp_pos,
    genepos = position_data$gene_pos,
    cisDist = CIS_DISTANCE
  )
  
  message(sprintf("Analysis complete.\ncis-eQTLs found: %d\ntrans-eQTLs found: %d",
                  nrow(me$cis$eqtls),
                  nrow(me$all$eqtls)))
  
  return(me)
}

# ---- Function: Summarize results ----
summarize_eqtl_results <- function(eqtl_results, output_file) {
  #' Create summary statistics of eQTL findings
  #'
  #' @param eqtl_results MatrixEQTL results object.
  #' @param output_file Character. Output file for summary.
  #'
  #' @return Invisible. Writes summary to file.
  
  message("Generating eQTL summary...")
  
  sink(output_file)
  
  cat("eQTL Analysis Summary\n")
  cat(strrep("=", 80), "\n\n")
  
  cat("CIS-EQTL RESULTS\n")
  cat("-", nrow(eqtl_results$cis$eqtls), "significant cis-eQTLs\n")
  cat("-", length(unique(eqtl_results$cis$eqtls$snps)), "unique SNPs\n")
  cat("-", length(unique(eqtl_results$cis$eqtls$gene)), "unique genes\n")
  cat("-", "Min p-value:", min(eqtl_results$cis$eqtls$pvalue), "\n")
  cat("-", "Max effect size:", max(abs(eqtl_results$cis$eqtls$beta)), "\n\n")
  
  cat("TRANS-EQTL RESULTS\n")
  cat("-", nrow(eqtl_results$all$eqtls), "significant trans-eQTLs\n")
  cat("-", "Min p-value:", min(eqtl_results$all$eqtls$pvalue), "\n")
  cat("-", "Max effect size:", max(abs(eqtl_results$all$eqtls$beta)), "\n\n")
  
  cat("FDR INFORMATION\n")
  cat("-", "cis-eQTL FDR:", eqtl_results$cis$ftdrcis, "\n")
  cat("-", "trans-eQTL FDR:", eqtl_results$all$fdr, "\n")
  
  sink()
  message(sprintf("Summary saved: %s", output_file))
}

# ---- Main processing pipeline ----
process_eqtl_analysis <- function(genotype_file,
                                  expression_file,
                                  covariate_file,
                                  snp_pos_file,
                                  gene_pos_file,
                                  output_dir = "./",
                                  output_prefix = "eqtl") {
  #' Complete eQTL analysis pipeline
  #'
  #' @param genotype_file Character. Path to genotype matrix.
  #' @param expression_file Character. Path to expression matrix.
  #' @param covariate_file Character. Path to covariate matrix.
  #' @param snp_pos_file Character. Path to SNP positions.
  #' @param gene_pos_file Character. Path to gene positions.
  #' @param output_dir Character. Directory for output files.
  #' @param output_prefix Character. Prefix for output files.
  #'
  #' @return List with eQTL results.
  
  message("\n" %+% strrep("=", 80))
  message("eQTL Analysis Pipeline")
  message(strrep("=", 80))
  
  # Step 1: Load data
  message("\nStep 1: Loading data matrices...")
  genotype <- readr::read_csv(genotype_file)
  expression <- readr::read_csv(expression_file)
  covariates <- readr::read_csv(covariate_file)
  
  # Step 2: Prepare matrices
  message("\nStep 2: Preparing eQTL data matrices...")
  sliced_data <- prepare_eqtl_matrices(
    as.data.frame(genotype),
    as.data.frame(expression),
    as.data.frame(covariates)
  )
  
  # Step 3: Load positions
  message("\nStep 3: Loading genomic position information...")
  position_data <- load_position_files(snp_pos_file, gene_pos_file)
  
  # Step 4: Run analysis
  message("\nStep 4: Running eQTL analysis...")
  cis_output <- sprintf("%s/%s_cis_eqtls_Jan2024.txt", output_dir, output_prefix)
  trans_output <- sprintf("%s/%s_trans_eqtls_Jan2024.txt", output_dir, output_prefix)
  
  eqtl_results <- run_eqtl_analysis(sliced_data, position_data,
                                     cis_output, trans_output)
  
  # Step 5: Summarize results
  message("\nStep 5: Generating summary statistics...")
  summary_file <- sprintf("%s/%s_summary_Jan2024.txt", output_dir, output_prefix)
  summarize_eqtl_results(eqtl_results, summary_file)
  
  message("\n" %+% strrep("=", 80))
  message("eQTL analysis complete!\n")
  
  return(eqtl_results)
}
