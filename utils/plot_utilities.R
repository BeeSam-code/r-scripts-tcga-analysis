#' ============================================================================
#' Utility: plot_utilities.R
#' Purpose: Visualization helper functions for TCGA analysis
#' Author: Sam Boudeau
#' Date: 2024-01-15
#'
#' Description:
#'   Provides functions for creating publication-quality plots including
#'   PCA scatter plots, admixture bar plots, Manhattan plots, and Q-Q plots.
#'
#' Functions:
#'   - plot_pca(): PCA scatter plot with population coloring
#'   - plot_admixture(): Stacked bar plot of ancestry proportions
#'   - plot_manhattan(): Manhattan plot for eQTL results
#'   - plot_qq(): Q-Q plot for p-value distributions
#'
#' ============================================================================

library(ggplot2)
library(dplyr)

# ---- Function: Create PCA plot ----
plot_pca <- function(pca_df, pc1 = "PC1", pc2 = "PC2",
                     color_var = "population", title = "Principal Component Analysis") {
  #' Create PCA scatter plot
  #'
  #' @param pca_df Data frame with PC columns and metadata.
  #' @param pc1 Character. Column name for first PC.
  #' @param pc2 Character. Column name for second PC.
  #' @param color_var Character. Column for coloring points.
  #' @param title Character. Plot title.
  #'
  #' @return ggplot object.
  #'
  #' @examples
  #'   p <- plot_pca(pca_results, color_var = "population")
  #'   print(p)
  
  p <- ggplot(pca_df, aes_string(x = pc2, y = pc1, color = color_var)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) +
    labs(title = title, color = "Population")
  
  return(p)
}

# ---- Function: Create admixture plot ----
plot_admixture <- function(admix_df, samples = NULL,
                           title = "Ancestry Admixture Proportions") {
  #' Create stacked bar plot of admixture proportions
  #'
  #' @param admix_df Data frame with ancestry proportions.
  #' @param samples Character vector. Samples to include (if NULL, uses all).
  #' @param title Character. Plot title.
  #'
  #' @return ggplot object.
  
  if (!is.null(samples)) {
    admix_df <- admix_df[admix_df$Sample_ID %in% samples, ]
  }
  
  # Prepare data for stacked bar plot
  admix_long <- admix_df %>%
    select(-Estimated_Ancestry) %>%
    pivot_longer(-Sample_ID, names_to = "Ancestry", values_to = "Proportion")
  
  p <- ggplot(admix_long, aes(x = Sample_ID, y = Proportion, fill = Ancestry)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = title, x = "Sample", y = "Ancestry Proportion")
  
  return(p)
}

# ---- Function: Create Manhattan plot ----
plot_manhattan <- function(eqtl_results, title = "Manhattan Plot",
                           pvalue_threshold = 1e-5) {
  #' Create Manhattan plot for eQTL results
  #'
  #' @param eqtl_results Data frame with eQTL results.
  #' @param title Character. Plot title.
  #' @param pvalue_threshold Numeric. Significance threshold line.
  #'
  #' @return ggplot object.
  
  # Add chromosome if not present
  if (!"chr" %in% colnames(eqtl_results)) {
    message("Warning: chromosome information not found in results")
    return(NULL)
  }
  
  p <- ggplot(eqtl_results, aes(x = pos, y = -log10(pvalue), color = chr)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "red") +
    facet_wrap(~chr, scales = "free_x", nrow = 2) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    labs(title = title, x = "Position", y = "-log10(p-value)")
  
  return(p)
}

# ---- Function: Create Q-Q plot ----
plot_qq <- function(pvalues, title = "Q-Q Plot") {
  #' Create Q-Q plot for p-value distribution
  #'
  #' @param pvalues Numeric vector of p-values.
  #' @param title Character. Plot title.
  #'
  #' @return ggplot object.
  
  # Calculate expected vs observed
  n <- length(pvalues)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(pvalues))
  
  qq_df <- data.frame(Expected = expected, Observed = observed)
  
  p <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = title,
         x = "Expected -log10(p-value)",
         y = "Observed -log10(p-value)")
  
  return(p)
}
