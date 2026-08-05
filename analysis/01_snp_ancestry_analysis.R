#' ============================================================================
#' Script: 01_snp_ancestry_analysis.R
#' Purpose: Perform SNP-based ancestry inference using PCA and admixture
#' Author: Sam Boudeau
#' Date: 2024-01-15
#' Last Updated: 2024-08-05
#'
#' Description:
#'   Performs SNP-based ancestry analysis through multiple steps:
#'   1. Converts VCF to efficient GDS format
#'   2. Performs PCA on genotype data
#'   3. Maps reference populations (CEU, YRI, JPT, CHB) to PCA space
#'   4. Estimates ancestry proportions for TCGA samples
#'   5. Classifies samples as European, African, or admixed
#'   6. Generates publication-quality PCA and admixture plots
#'
#' Dependencies:
#'   - SNPRelate (>= 1.20.0)
#'   - gdsfmt (>= 1.24.0)
#'   - dplyr (>= 1.0.0)
#'   - readr (>= 2.0.0)
#'   - ggplot2 (>= 3.3.0)
#'   - BiocManager
#'
#' Input:
#'   - VCF file with merged TCGA and 1000 Genomes samples
#'   - Clinical data with self-reported race/ethnicity
#'
#' Output:
#'   - SNPrelate_COADREAD_PCAJan2024.png: PCA plot
#'   - Ancestry classification table
#'   - Admixture proportion estimates
#'
#' Notes:
#'   - Original script from TCGA_analysis_inR refactored
#'   - Reference populations: CEU (European), YRI (African), JPT/CHB (Asian)
#'   - Ancestry threshold: 0.60 for classification
#'
#' ============================================================================

source('utils/load_libraries.R')

# ---- Configuration Constants ----
NUM_SNPS_SAMPLE <- 1000
ANCESTRY_THRESHOLD <- 0.60
PNG_WIDTH <- 1400
PNG_HEIGHT <- 960
PNG_RES <- 300

# ---- Function: Convert VCF to GDS ----
convert_vcf_to_gds <- function(vcf_file, gds_file) {
  #' Convert VCF format to efficient GDS format
  #'
  #' @param vcf_file Character. Path to VCF file.
  #' @param gds_file Character. Output GDS file path.
  #'
  #' @return Invisible. Creates GDS file.
  
  message(sprintf("Converting VCF to GDS: %s", vcf_file))
  SNPRelate::snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only")
  SNPRelate::snpgdsSummary(gds_file)
  message("VCF to GDS conversion complete")
}

# ---- Function: Perform PCA ----
perform_pca_analysis <- function(gds_file, num_snps = NUM_SNPS_SAMPLE) {
  #' Perform principal component analysis on genotype data
  #'
  #' @param gds_file Character. Path to GDS file.
  #' @param num_snps Integer. Number of SNPs to sample.
  #'
  #' @return List with PCA results.
  
  message("Opening GDS file...")
  genofile <- SNPRelate::snpgdsOpen(gds_file)
  
  # Sample SNPs for initial PCA
  message(sprintf("Sampling %d SNPs for PCA...", num_snps))
  snp_set <- sample(SNPRelate::read.gdsn(
    SNPRelate::index.gdsn(genofile, "snp.id")
  ), num_snps)
  
  # Perform PCA
  message("Performing PCA...")
  pca <- SNPRelate::snpgdsPCA(
    genofile,
    snp.id = snp_set,
    remove.monosnp = TRUE,
    algorithm = "exact",
    eigen.method = "DSPEVX",
    num.thread = 1L
  )
  
  SNPRelate::snpgdsClose(genofile)
  return(pca)
}

# ---- Function: Prepare metadata ----
prepare_ancestry_metadata <- function(pca_results, clinical_data) {
  #' Prepare population metadata for ancestry analysis
  #'
  #' @param pca_results List from perform_pca_analysis().
  #' @param clinical_data Data frame with clinical information.
  #'
  #' @return Data frame with sample IDs and population labels.
  
  data_ids <- pca_results$sample.id
  
  # Create population labels (self-reported race for TCGA)
  id_info <- clinical_data %>%
    dplyr::select(submitter_id.x, race) %>%
    dplyr::rename(Individual_ID = submitter_id.x, Population = race)
  
  # Define reference populations
  subpops <- c("YRI", "CEU", "JPT", "CHB")
  
  pops_metadata <- data.frame(
    Individual_ID = data_ids,
    Population = NA
  )
  
  # Map known populations
  for (i in 1:nrow(id_info)) {
    idx <- which(pops_metadata$Individual_ID == id_info$Individual_ID[i])
    if (length(idx) > 0) {
      pops_metadata$Population[idx] <- id_info$Population[i]
    }
  }
  
  return(pops_metadata)
}

# ---- Function: Estimate ancestry proportions ----
estimate_ancestry_admixture <- function(gds_file, pca_results, pops_metadata, 
                                        threshold = ANCESTRY_THRESHOLD) {
  #' Estimate ancestry proportions using SNPRelate
  #'
  #' @param gds_file Character. Path to GDS file.
  #' @param pca_results List from perform_pca_analysis().
  #' @param pops_metadata Data frame with population info.
  #' @param threshold Numeric. Threshold for ancestry classification.
  #'
  #' @return Data frame with ancestry estimates.
  
  message("Estimating ancestry proportions...")
  
  # Define reference populations
  subpops <- c("YRI", "CEU", "JPT", "CHB")
  groups <- list(
    CEU = pca_results$sample.id[pops_metadata$Population == "CEU"],
    YRI = pca_results$sample.id[pops_metadata$Population == "YRI"],
    ASN = pca_results$sample.id[pops_metadata$Population %in% c("CHB", "JPT")]
  )
  
  # Remove NA values
  groups_clean <- lapply(groups, function(x) x[!is.na(x)])
  
  # Estimate admixture proportions
  admix_prop <- SNPRelate::snpgdsAdmixProp(
    pca_results,
    groups = groups_clean,
    bound = TRUE
  )
  
  # Create results data frame
  admix_df <- as.data.frame(admix_prop)
  admix_df$Sample_ID <- rownames(admix_df)
  
  # Classify ancestry
  admix_df$Estimated_Ancestry <- NA
  admix_df$Estimated_Ancestry[admix_df$YRI >= threshold] <- "African"
  admix_df$Estimated_Ancestry[admix_df$CEU >= threshold] <- "European"
  admix_df$Estimated_Ancestry[is.na(admix_df$Estimated_Ancestry)] <- "Admixed"
  
  message(sprintf("Ancestry classification complete\nAfrican: %d, European: %d, Admixed: %d",
                  sum(admix_df$Estimated_Ancestry == "African"),
                  sum(admix_df$Estimated_Ancestry == "European"),
                  sum(admix_df$Estimated_Ancestry == "Admixed")))
  
  return(admix_df)
}

# ---- Function: Create PCA plot ----
plot_pca_ancestry <- function(pca_results, pops_metadata, output_file) {
  #' Create PCA scatter plot colored by population
  #'
  #' @param pca_results List from perform_pca_analysis().
  #' @param pops_metadata Data frame with population info.
  #' @param output_file Character. Path for output PNG file.
  #'
  #' @return Invisible. Creates PNG file.
  
  message("Creating PCA plot...")
  
  # Prepare data
  pc_percent <- pca_results$varprop * 100
  
  pca_df <- data.frame(
    Sample_ID = pca_results$sample.id,
    EV1 = signif(pca_results$eigenvect[, 1], digits = 2),
    EV2 = signif(pca_results$eigenvect[, 2], digits = 2),
    EV3 = signif(pca_results$eigenvect[, 3], digits = 2),
    Population = pops_metadata$Population[match(
      pca_results$sample.id,
      pops_metadata$Individual_ID
    )]
  )
  
  # Create plot
  png(output_file, pointsize = 10, width = PNG_WIDTH, height = PNG_HEIGHT, res = PNG_RES)
  
  plot(pca_df$EV2, pca_df$EV1,
       col = as.integer(as.factor(pca_df$Population)),
       xlab = sprintf("PC2 (%.2f%%)", pc_percent[2]),
       ylab = sprintf("PC1 (%.2f%%)", pc_percent[1]),
       cex = 1.5,
       pch = 19,
       main = "SNPRelate COAD-READ PCA",
       font.main = 2)
  
  legend("topright",
         legend = levels(as.factor(pca_df$Population)),
         pch = 19,
         col = 1:length(levels(as.factor(pca_df$Population))),
         cex = 0.7)
  
  dev.off()
  message(sprintf("PCA plot saved: %s", output_file))
}

# ---- Main processing pipeline ----
process_ancestry_analysis <- function(vcf_file,
                                      gds_file,
                                      clinical_file,
                                      output_dir = "./",
                                      output_prefix = "SNPrelate_COADREAD") {
  #' Complete ancestry analysis pipeline
  #'
  #' @param vcf_file Character. Path to merged VCF file.
  #' @param gds_file Character. Output GDS file path.
  #' @param clinical_file Character. Path to clinical data.
  #' @param output_dir Character. Directory for output files.
  #' @param output_prefix Character. Prefix for output files.
  #'
  #' @return List with PCA results and ancestry estimates.
  
  message("\n" %+% strrep("=", 80))
  message("SNP Ancestry Analysis Pipeline")
  message(strrep("=", 80))
  
  # Step 1: Convert VCF to GDS
  message("\nStep 1: Converting VCF to GDS...")
  convert_vcf_to_gds(vcf_file, gds_file)
  
  # Step 2: Perform PCA
  message("\nStep 2: Performing PCA analysis...")
  pca_results <- perform_pca_analysis(gds_file)
  
  # Step 3: Load clinical data
  message("\nStep 3: Loading clinical data...")
  clinical_data <- readr::read_csv(clinical_file)
  
  # Step 4: Prepare metadata
  message("\nStep 4: Preparing ancestry metadata...")
  pops_metadata <- prepare_ancestry_metadata(pca_results, clinical_data)
  
  # Step 5: Estimate ancestry
  message("\nStep 5: Estimating ancestry proportions...")
  admix_results <- estimate_ancestry_admixture(gds_file, pca_results, pops_metadata)
  
  # Step 6: Create visualization
  message("\nStep 6: Creating visualizations...")
  pca_plot_file <- sprintf("%s/%s_PCA_Jan2024.png", output_dir, output_prefix)
  plot_pca_ancestry(pca_results, pops_metadata, pca_plot_file)
  
  # Step 7: Save results
  message("\nStep 7: Saving results...")
  admix_file <- sprintf("%s/%s_ancestry_estimates_Jan2024.csv", output_dir, output_prefix)
  write.csv(admix_results, file = admix_file, row.names = FALSE)
  message(sprintf("Ancestry estimates saved: %s", admix_file))
  
  message("\n" %+% strrep("=", 80))
  message("Ancestry analysis complete!\n")
  
  return(list(
    pca_results = pca_results,
    pops_metadata = pops_metadata,
    admix_results = admix_results
  ))
}
