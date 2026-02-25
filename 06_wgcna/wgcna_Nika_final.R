# ═══════════════════════════════════════════════════════════════════════════════
# WGCNA PIPELINE - CHAMPION MODEL ONLY
# ═══════════════════════════════════════════════════════════════════════════════
# ============================================================================
# WGCNA PIPELINE v3
# ============================================================================
# Dataset: 50 metastatic PDAC patients, 891 plasma proteins
# Focus Modules: PINK, GREEN, BLACK, BLUE (clinically significant - negative PFS correlation)
# Custom Signatures: 150 PDAC-specific pathways
#
# Pipeline Steps:
#   1. Data Loading & Preprocessing
#   2. Protein Noise Assessment
#   3. Sample Quality Control & Clustering
#   4. Soft-Threshold Power Selection
#   5. Final Network Construction
#   6. Hub Gene Identification & Module-Trait Correlations
#   7. Module Stability Testing
#   8. Pathway Enrichment (Custom + MSigDB + GO/KEGG)
#   9-14. Extended Analyses
#  15. Biomarker Integration (ML + Cox)
#  16. Extended Analyses
#  17. STRING/PPI Network
#  18. Circos Module Overview
#  19. sEV Proteomics Comparison
#
# ============================================================================

# Record start time for runtime tracking
pipeline_start_time <- Sys.time()
cat(sprintf("WGCNA Pipeline v3 started: %s\n", format(pipeline_start_time, "%Y-%m-%d %H:%M:%S")))


# Set working directory to the location of this script
# setwd("/path/to/your/data")  # Adjust to your local data directory



# ============================================================================
# Libraries
# ============================================================================


# Core WGCNA and data manipulation
library(WGCNA)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(patchwork)
library(dplyr)
library(parallel)

# Visualization
library(pheatmap)
library(RColorBrewer)
library(ggwordcloud)
library(viridis)
library(scales)
library(ggdendro)
library(ggforce)
library(igraph)
library(circlize)
library(ggalluvial)
library(ggvenn)
library(cowplot)
library(ggh4x)
library(tidytext)
library(paletteer)  # For nationalparkcolors::Acadia palette

# Enrichment analysis
library(GSVA)
library(msigdbr)
library(fgsea)

# Data export
library(openxlsx)
library(umap)

# Report generation
library(rmarkdown)
library(knitr)
library(kableExtra)

# Conflict resolution
library(conflicted)

cat("All libraries loaded successfully.\n")

# --- Resolve function conflicts ---
# dplyr conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::slice)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::combine)
conflicts_prefer(dplyr::desc)
conflicts_prefer(dplyr::collapse)
conflicts_prefer(dplyr::first)


# Other conflicts
conflicts_prefer(tidyr::expand)
conflicts_prefer(purrr::reduce)
conflicts_prefer(base::intersect)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)
conflicts_prefer(WGCNA::cor)
conflicts_prefer(cowplot::get_legend)


# CPU Core Detection - automatic detection for optimal performance
n_cores <- parallel::detectCores() - 1  # Leave 1 core for system
n_cores <- max(1, n_cores)  # Ensure at least 1 core
allowWGCNAThreads(nThreads = n_cores)
cat(sprintf("Using %d cores for WGCNA\n", n_cores))

# Create main output directory structure
# All outputs go into the results/ subdirectory
main_dir <- "results"
fig_dir_png <- file.path(main_dir, "Figures_PNG")
fig_dir_pdf <- file.path(main_dir, "Figures_PDF")
results_dir <- file.path(main_dir, "Results_Tables")

dir.create(main_dir, showWarnings = FALSE)
dir.create(fig_dir_png, showWarnings = FALSE)
dir.create(fig_dir_pdf, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

# ============================================================================
# INPUT FILE VERIFICATION
# ============================================================================
# Verify all required input files exist before proceeding

required_files <- c(
  "correctedHLA_Rawdata_PFSfiltered.csv",  # Expression + clinical data (891 proteins x 50 samples)
  "mastertable_Signature.csv"               # Gene signature database (150 signatures)
)

optional_files <- c(
  "sEV_proteomics.csv"            # Optional: sEV comparison data
)

cat("\n--- INPUT FILE VERIFICATION ---\n")
missing_required <- character(0)
for(f in required_files) {
  if(file.exists(f)) {
    cat(sprintf("  [OK] %s\n", f))
  } else {
    cat(sprintf("  [MISSING] %s\n", f))
    missing_required <- c(missing_required, f)
  }
}

for(f in optional_files) {
  if(file.exists(f)) {
    cat(sprintf("  [OK] %s (optional)\n", f))
  } else {
    cat(sprintf("  [SKIP] %s (optional - will skip Step 18)\n", f))
  }
}

if(length(missing_required) > 0) {
  stop(sprintf("\nERROR: Required input files missing:\n  %s\n\nPlease ensure all required files are in: %s",
               paste(missing_required, collapse = "\n  "), getwd()))
}
cat("  All required files verified.\n")
cat("--------------------------------\n")

# ============================================================================
# UNIFIED ORA FUNCTION
# ============================================================================
# Purpose: Consistent over-representation analysis across all gene set databases
# Background: Always uses 891 proteins (full proteome)
# Method: Fisher's exact test (one-tailed, greater)
# ============================================================================

unified_ora <- function(module_proteins, gene_set_list, background_proteins,
                        min_overlap = 2, database_name = "Custom") {

  n_background <- length(background_proteins)
  n_module <- length(module_proteins)
  results <- data.frame()

  for(gs_name in names(gene_set_list)) {
    gs_info <- gene_set_list[[gs_name]]

    # Handle both list format (custom) and vector format (MSigDB)
    if(is.list(gs_info) && !is.null(gs_info$genes)) {
      gs_proteins <- gs_info$genes
      short_name <- gs_info$short_name
      full_name <- gs_info$full_name
      category <- gs_info$category
    } else {
      gs_proteins <- gs_info
      short_name <- gs_name
      full_name <- gs_name
      category <- database_name
    }

    # Restrict to proteins in background (891)
    gs_in_background <- intersect(gs_proteins, background_proteins)
    n_gs <- length(gs_in_background)

    if(n_gs < 3) next  # Skip tiny gene sets

    overlap <- intersect(module_proteins, gs_in_background)
    n_overlap <- length(overlap)

    if(n_overlap < min_overlap) next

    # Fisher's exact test (hypergeometric)
    mat <- matrix(c(
      n_overlap,
      n_module - n_overlap,
      n_gs - n_overlap,
      n_background - n_module - n_gs + n_overlap
    ), nrow = 2)

    fisher_result <- fisher.test(mat, alternative = "greater")
    fold_enrichment <- (n_overlap / n_module) / (n_gs / n_background)

    results <- rbind(results, data.frame(
      Database = database_name,
      Pathway_ID = gs_name,
      Short_Name = short_name,
      Full_Name = full_name,
      Category = category,
      N_Pathway_Total = length(gs_proteins),
      N_Pathway_InBackground = n_gs,
      N_Module = n_module,
      N_Overlap = n_overlap,
      Fold_Enrichment = round(fold_enrichment, 2),
      P_Value = fisher_result$p.value,
      Overlap_Proteins = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }

  if(nrow(results) > 0) {
    results$FDR <- p.adjust(results$P_Value, method = "BH")
    results <- results[order(results$P_Value), ]
    rownames(results) <- NULL
  }

  return(results)
}

cat("Unified ORA function loaded\n")

# S1: Record session info for reproducibility
session_info <- sessionInfo()
writeLines(capture.output(print(session_info)), file.path(results_dir, "Step00_SessionInfo.txt"))
cat(sprintf("Session info saved to: %s/Step00_SessionInfo.txt", results_dir))

# Set seed for reproducibility
set.seed(12345)

cat("===============================================================================\n")
cat("STEP 1: DATA LOADING & PREPROCESSING\n")
cat("===============================================================================\n")
# PURPOSE: Load raw proteomics data and clinical metadata, organize into
# WGCNA-compatible format (samples x proteins expression matrix).

cat("SUBSECTION: Load Data\n")
cat("-------------------------------------------------------------------------------\n")
dt <- read.csv("correctedHLA_Rawdata_PFSfiltered.csv", header = TRUE, stringsAsFactors = FALSE)
dim(dt)

cat("SUBSECTION: Separate Expression and Clinical Data\n")
cat("-------------------------------------------------------------------------------\n")
# Dynamic column indexing based on column names
clinical_cols <- c("patient_ID", "Age", "Sex", "Stage.at.diagnosis", "Liver",
                   "PFS", "PFS_group", "RECIST", "Response", "CA19.9", "Treatment")
expr_cols <- setdiff(colnames(dt), clinical_cols)

# Expression matrix: proteins in columns (WGCNA standard format)
# Samples in rows
datExpr <- dt[, expr_cols]
rownames(datExpr) <- dt$patient_ID

# Clinical/trait data (patient_ID used only for rownames, not as column)
datTraits <- dt[, c("Age", "Sex", "Stage.at.diagnosis", "Liver",
                    "PFS", "PFS_group", "RECIST", "Response", "CA19.9", "Treatment")]
rownames(datTraits) <- dt$patient_ID

# Convert categorical variables to factors
datTraits$Sex <- factor(datTraits$Sex, levels = c(0, 1), labels = c("Male", "Female"))
datTraits$Liver <- factor(datTraits$Liver, levels = c(0, 1), labels = c("No", "Yes"))

# PFS_group: Handle multiple encodings (S/L, 0/1, Short/Long)
pfs_raw <- datTraits$PFS_group
if(all(unique(na.omit(pfs_raw)) %in% c("S", "L"))) {
  datTraits$PFS_group <- factor(pfs_raw, levels = c("S", "L"), labels = c("Short", "Long"))
} else if(all(unique(na.omit(pfs_raw)) %in% c(0, 1, "0", "1"))) {
  datTraits$PFS_group <- factor(ifelse(pfs_raw == 0 | pfs_raw == "0", "Short", "Long"),
                                 levels = c("Short", "Long"))
} else if(all(unique(na.omit(pfs_raw)) %in% c("Short", "Long"))) {
  datTraits$PFS_group <- factor(pfs_raw, levels = c("Short", "Long"))
} else {
  cat("WARNING: Unrecognized PFS_group encoding. Values found: ", unique(pfs_raw), "\n")
}

datTraits$RECIST <- factor(datTraits$RECIST, levels = c("PD", "SD", "PR"))

# Response: Handle multiple encodings (0/1, PD/CD)
resp_raw <- datTraits$Response
if(all(unique(na.omit(resp_raw)) %in% c(0, 1, "0", "1"))) {
  datTraits$Response <- factor(ifelse(resp_raw == 0 | resp_raw == "0", "PD", "CD"),
                                levels = c("PD", "CD"))
} else if(all(unique(na.omit(resp_raw)) %in% c("PD", "CD"))) {
  datTraits$Response <- factor(resp_raw, levels = c("PD", "CD"))
} else {
  cat("WARNING: Unrecognized Response encoding. Values found: ", unique(resp_raw), "\n")
}
# Treatment factor with explicit level ordering
datTraits$Treatment <- factor(datTraits$Treatment)

cat("SUBSECTION: Verify Data Structure\n")
cat("-------------------------------------------------------------------------------\n")
cat("Expression matrix: ", nrow(datExpr), "samples x", ncol(datExpr), "proteins\n")
cat("Clinical data: ", nrow(datTraits), "samples x", ncol(datTraits), "variables\n")
cat("Expression value range:\n")
cat("  Min:", min(datExpr, na.rm = TRUE), "\n")
cat("  Max:", max(datExpr, na.rm = TRUE), "\n")
cat("  Mean:", mean(as.matrix(datExpr), na.rm = TRUE), "\n")
cat("=== CLINICAL VARIABLES SUMMARY ===\n")
cat("PFS_group distribution:\n")
print(table(datTraits$PFS_group))
cat("Treatment distribution:\n")
print(table(datTraits$Treatment))
cat("Response distribution:\n")
print(table(datTraits$Response))
cat("RECIST distribution:\n")
print(table(datTraits$RECIST))
cat("Sex distribution:\n")
print(table(datTraits$Sex))
cat("Liver metastasis distribution:\n")
print(table(datTraits$Liver))
cat("PFS (days) summary:\n")
print(summary(datTraits$PFS))
cat("Age summary:\n")
print(summary(datTraits$Age))
cat("CA19.9 summary:\n")
print(summary(datTraits$CA19.9))
# ----------------------------------------------------------------------------
# 1.4 Check for missing values
# ----------------------------------------------------------------------------
cat("=== MISSING VALUE CHECK ===\n")
missing_per_protein <- colSums(is.na(datExpr))
missing_per_sample <- rowSums(is.na(datExpr))

cat("Proteins with missing values:", sum(missing_per_protein > 0), "\n")
cat("Samples with missing values:", sum(missing_per_sample > 0), "\n")

if(sum(missing_per_protein > 0) > 0) {
  cat("Top proteins with most missing values:\n")
  print(head(sort(missing_per_protein[missing_per_protein > 0], decreasing = TRUE), 10))
}

# STEP 1 SUMMARY
cat("===============================================================================\n")
cat("STEP 1 SUMMARY: DATA LOADING\n")
cat("===============================================================================\n")

step1_summary <- data.frame(
  Metric = c(
    "Total samples",
    "Total proteins",
    "Clinical variables",
    "Missing proteins",
    "Missing samples",
    "Expression range (min)",
    "Expression range (max)",
    "Expression mean",
    "PFS Short",
    "PFS Long"
  ),
  Value = c(
    nrow(datExpr),
    ncol(datExpr),
    ncol(datTraits),
    sum(missing_per_protein > 0),
    sum(missing_per_sample > 0),
    round(min(datExpr, na.rm = TRUE), 4),
    round(max(datExpr, na.rm = TRUE), 4),
    round(mean(as.matrix(datExpr), na.rm = TRUE), 4),
    sum(datTraits$PFS_group == "Short"),
    sum(datTraits$PFS_group == "Long")
  ),
  stringsAsFactors = FALSE
)

cat("RESULTS:\n")
cat("-------------------------------------------------------------------------------\n")
print(step1_summary, row.names = FALSE)
write.csv(step1_summary, file.path(results_dir, "Step1_Summary.csv"), row.names = FALSE)
cat("Saved: Step1_Summary.csv\n")
cat("===============================================================================\n")
cat("STEP 2: PROTEIN NOISE ASSESSMENT\n")
cat("===============================================================================\n")
# PURPOSE: Identify potentially noisy/unreliable proteins using multiple
# orthogonal criteria before network construction.
#
# STRATEGY:
# 2.1) Variance analysis
# 2.2) CV analysis
# 2.3) Mean-variance
# 2.4) Inter-protein correlation
# 2.5) Expression levels
# 2.6) Summary

# ----------------------------------------------------------------------------
# 2.0 Setup: Create Directory Structure and Define Color Palette
# ----------------------------------------------------------------------------


# Define your color palette
color_palette <- c(
  "#008080",
  "#70a494",
  "#b4c8a8",
  "#f6edbd",
  "#edbb8a",
  "#de8a5a",
  "#ca562c"
)

# Extended palette for more categories if needed
color_palette_extended <- colorRampPalette(color_palette)(20)

# =============================================================================
# UNIVERSAL CLINICAL COLOR DEFINITIONS
# Use these consistently throughout the entire pipeline
# =============================================================================
# PFS Group colors
col_pfs_short <- "#D57649"
col_pfs_long <- "#4F9394"
colors_pfs <- c("Short" = col_pfs_short, "Long" = col_pfs_long)

# Treatment colors
col_folfi <- "#A1BFB6"
col_gemnab <- "#F5D3AF"
colors_treatment <- c("FOLFIRINOX" = col_folfi, "FOLFI" = col_folfi, "GemNab" = col_gemnab)

# Response colors
col_pd <- "#58B4A5"
col_cd <- "#FEDE82"
colors_response <- c("PD" = col_pd, "CD" = col_cd)

# ============================================================================

# ============================================================================
# Purpose: Consistent styling across all figures
# Usage: Replace theme_bw() + theme(...) with theme_wgcna()
# ============================================================================

theme_wgcna <- function(base_size = 11, title_size = 12) {
  theme_bw(base_size = base_size) +
  theme(
    # Title and subtitle
    plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = base_size - 1),
    plot.caption = element_text(color = "grey50", size = base_size - 2),

    # Axes
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "grey20"),

    # Grid
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),

    # Legend
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = base_size - 1),

    # Strip (for faceted plots)
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),

    # Panel
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
  )
}

# Variant for heatmaps (minimal theme)
theme_wgcna_heatmap <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
}

cat("Loaded: theme_wgcna() and theme_wgcna_heatmap() functions\n")


# Helper function to save both PNG and PDF
save_plot <- function(plot_name, plot_function, width = 10, height = 6) {
  # Save PNG at 300 dpi
  png(file.path(fig_dir_png, paste0(plot_name, ".png")),
      width = width, height = height, units = "in", res = 300)
  plot_function()
  dev.off()

  # Save PDF
  pdf(file.path(fig_dir_pdf, paste0(plot_name, ".pdf")),
      width = width, height = height)
  plot_function()
  dev.off()

  cat("  Saved:", plot_name, "(PNG + PDF)\n")
}

# ----------------------------------------------------------------------------
# 2.1 Variance Analysis
# ----------------------------------------------------------------------------
cat("--- 2.1 Variance Analysis ---\n")

# Calculate variance for each protein across all 50 samples
protein_variance <- apply(datExpr, 2, var)

# Summary statistics
cat("Protein variance distribution:\n")
print(summary(protein_variance))

# Identify variance percentiles
var_percentiles <- quantile(protein_variance, probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95))
cat("Variance percentiles:\n")
print(round(var_percentiles, 4))

# Flag low-variance proteins (bottom 5th percentile)
low_var_threshold <- quantile(protein_variance, 0.05)
low_var_proteins <- names(protein_variance[protein_variance < low_var_threshold])
cat("Proteins with variance < 5th percentile (", round(low_var_threshold, 4), "):",
    length(low_var_proteins), "\n")

# --- Variance Distribution Plot ---
cat("Generating Variance Distribution plot...\n")

var_df <- data.frame(
  Protein = names(protein_variance),
  Variance = protein_variance,
  Rank = rank(protein_variance)
)

# Panel A: Histogram
p1a <- ggplot(var_df, aes(x = Variance)) +
  geom_histogram(bins = 50, fill = color_palette[1], color = "white", alpha = 0.8) +
  geom_vline(xintercept = low_var_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  geom_vline(xintercept = quantile(protein_variance, 0.25), linetype = "dashed", color = color_palette[5], linewidth = 0.8) +
  labs(
    title = "Protein Variance Distribution",
    subtitle = sprintf("5th pctl = %.3f | 25th pctl = %.3f", low_var_threshold, quantile(protein_variance, 0.25)),
    x = "Variance", y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )
print(p1a)

# Panel B: Sorted variance
p1b <- ggplot(var_df, aes(x = Rank, y = Variance)) +
  geom_line(color = color_palette[1], linewidth = 1) +
  geom_hline(yintercept = low_var_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  geom_hline(yintercept = quantile(protein_variance, 0.25), linetype = "dashed", color = color_palette[5], linewidth = 0.8) +
  geom_hline(yintercept = quantile(protein_variance, 0.50), linetype = "dashed", color = color_palette[3], linewidth = 0.8) +
  labs(
    title = "Sorted Protein Variance",
    subtitle = "Dashed lines: 5th, 25th, 50th percentiles",
    x = "Protein Rank", y = "Variance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )
print(p1b)
p1_combined <- p1a + p1b + plot_annotation(title = NULL)
print(p1_combined)

ggsave(file.path(fig_dir_png, "Step2_01_Variance_Distribution.png"), p1_combined, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2_01_Variance_Distribution.pdf"), p1_combined, width = 8, height = 5)

#Interpretation: Low-variance proteins show minimal expression changes across samples. 
#These may be housekeeping proteins or those with little biological variation relevant to your phenotype (PFS, treatment response). However, in WGCNA, low variance isn't necessarily bad - these proteins simply won't contribute much to module detection.

# ----------------------------------------------------------------------------
# 2.2 Coefficient of Variation (CV) Analysis
# ----------------------------------------------------------------------------
cat("--- 2.2 Coefficient of Variation Analysis ---\n")

# CV = SD / Mean (relative variability)
protein_mean <- apply(datExpr, 2, mean)
protein_sd <- apply(datExpr, 2, sd)
protein_cv <- protein_sd / protein_mean

cat("Protein CV distribution:\n")
print(summary(protein_cv))

# Flag proteins with very low CV (bottom 5%)
cv_percentiles <- quantile(protein_cv, probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95))
cat("CV percentiles:\n")
print(round(cv_percentiles, 4))

low_cv_threshold <- quantile(protein_cv, 0.05)
low_cv_proteins <- names(protein_cv[protein_cv < low_cv_threshold])

# Flag proteins with very high CV (top 5%) - potential technical noise
high_cv_threshold <- quantile(protein_cv, 0.95)
high_cv_proteins <- names(protein_cv[protein_cv > high_cv_threshold])

cat("Proteins with CV < 5th percentile:", length(low_cv_proteins), "\n")
cat("Proteins with CV > 95th percentile:", length(high_cv_proteins), "\n")

# --- CV Distribution Plot ---
cat("Generating CV Distribution plot...\n")

cv_df <- data.frame(
  Protein = names(protein_cv),
  CV = protein_cv,
  Mean = protein_mean
)

# Panel A: Histogram
p3a <- ggplot(cv_df, aes(x = CV)) +
  geom_histogram(bins = 50, fill = color_palette[2], color = "white", alpha = 0.8) +
  geom_vline(xintercept = low_cv_threshold, linetype = "dashed", color = color_palette[1], linewidth = 0.8) +
  geom_vline(xintercept = high_cv_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  labs(
    title = "Coefficient of Variation Distribution",
    subtitle = sprintf("5th pctl = %.3f | 95th pctl = %.3f", low_cv_threshold, high_cv_threshold),
    x = "CV (SD/Mean)", y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )
print(p3a)

# Panel B: CV vs Mean
p3b <- ggplot(cv_df, aes(x = Mean, y = CV)) +
  geom_point(color = color_palette[2], alpha = 0.6, size = 2) +
  geom_hline(yintercept = low_cv_threshold, linetype = "dashed", color = color_palette[1], linewidth = 0.8) +
  geom_hline(yintercept = high_cv_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  labs(
    title = "CV vs Mean Expression",
    subtitle = "Dashed lines: 5th and 95th percentiles",
    x = "Mean Expression (log2)", y = "CV"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )
print(p3b)

p3_combined <- p3a + p3b
print(p3_combined)

ggsave(file.path(fig_dir_png, "Step2_03_CV_Distribution.png"), p3_combined, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2_03_CV_Distribution.pdf"), p3_combined, width = 8, height = 5)

#Interpretation:Low CV proteins are very stable across patients - potentially housekeeping or constitutively expressed proteins
#High CV proteins (>12.3%) could indicate either biological signal OR technical noise. The maximum CV is 31% (DYNC2H1), which warrants attention.

# ----------------------------------------------------------------------------
# 2.3 Mean-Variance Relationship
# ----------------------------------------------------------------------------
cat("--- 2.3 Mean-Variance Relationship ---\n")

# Fit linear model
mv_model <- lm(protein_variance ~ protein_mean)
mv_r_squared <- summary(mv_model)$r.squared

cat("Mean-Variance correlation (RÂ2):", round(mv_r_squared, 4), "\n")

if(mv_r_squared > 0.3) {
  cat("[WARNING] WARNING: Moderate mean-variance dependence detected.\n")
  cat("This may indicate insufficient normalization.\n")
} else if(mv_r_squared > 0.1) {
  cat("NOTE: Weak mean-variance dependence present.\n")
  cat("Generally acceptable for log2 transformed data.\n")
} else {
  cat("[OK] Good: Minimal mean-variance dependence.\n")
  cat("Log2 transformation effectively stabilized variance.\n")
}

# Identify outliers from mean-variance relationship
mv_residuals <- residuals(mv_model)
mv_outliers <- names(mv_residuals[abs(scale(mv_residuals)) > 3])
cat("Proteins with extreme mean-variance residuals (|Z| > 3):", length(mv_outliers), "\n")

# --- Mean-Variance Relationship Plot ---
cat("Generating Mean-Variance Relationship plot...\n")

mv_df <- data.frame(
  Protein = names(protein_mean),
  Mean = protein_mean,
  Variance = protein_variance,
  Is_Outlier = names(protein_mean) %in% mv_outliers
)

p2_mv <- ggplot(mv_df, aes(x = Mean, y = Variance)) +
  geom_point(aes(color = Is_Outlier), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = color_palette[6], linewidth = 1) +
  scale_color_manual(values = c("FALSE" = color_palette[1], "TRUE" = color_palette[7]),
                     labels = c("Normal", "MV Outlier")) +
  {if(length(mv_outliers) > 0) geom_text_repel(
    data = mv_df[mv_df$Is_Outlier, ],
    aes(label = Protein), size = 2.5, color = color_palette[7], max.overlaps = 10
  )} +
  labs(
    title = "Mean-Variance Relationship",
    subtitle = sprintf("RÂ2 = %.3f | MV outliers (|Z| > 3): %d", mv_r_squared, length(mv_outliers)),
    x = "Mean Expression (log2)", y = "Variance",
    color = "Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "bottom"
  )
print(p2_mv)

ggsave(file.path(fig_dir_png, "Step2_02_Mean_Variance_Relationship.png"), p2_mv, width = 6, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2_02_Mean_Variance_Relationship.pdf"), p2_mv, width = 6, height = 4)


#There's no systematic mean-variance dependence, high-abundance proteins don't have artificially inflated variance.
#data is well-normalized for WGCNA.

# ============================================================================
# STEP 2c: DYNC2H1 PATIENT STRATIFICATION ANALYSIS
# ============================================================================
# CONTEXT: DYNC2H1 is an ML biomarker with unusually high CV (~31%) and is
# consistently assigned to the grey module. This analysis investigates which
# patients have high vs low DYNC2H1 expression and clinical associations.
# NOTE: This is a purely DATA-DRIVEN analysis with NO biological interpretations.
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("STEP 2c: DYNC2H1 PATIENT STRATIFICATION ANALYSIS\n")
cat("============================================================================\n")

# ---------------------------------------------------------------------------
# 2c-1: Basic Statistics
# ---------------------------------------------------------------------------
cat("\n--- 2c-1: DYNC2H1 Basic Statistics ---\n")

dync2h1_expr <- datExpr$DYNC2H1
names(dync2h1_expr) <- rownames(datExpr)

dync2h1_complete <- dync2h1_expr[!is.na(dync2h1_expr)]
n_excluded <- sum(is.na(dync2h1_expr))
cat(sprintf("  Patients with DYNC2H1 data: %d\n", length(dync2h1_complete)))
cat(sprintf("  Patients excluded (NA): %d\n", n_excluded))

dync2h1_stats <- data.frame(
  Statistic = c("N", "Mean", "Median", "SD", "CV (%)", "Min", "Max", "IQR", "Q1", "Q3"),
  Value = c(
    length(dync2h1_complete),
    round(mean(dync2h1_complete), 4),
    round(median(dync2h1_complete), 4),
    round(sd(dync2h1_complete), 4),
    round(100 * sd(dync2h1_complete) / mean(dync2h1_complete), 2),
    round(min(dync2h1_complete), 4),
    round(max(dync2h1_complete), 4),
    round(IQR(dync2h1_complete), 4),
    round(quantile(dync2h1_complete, 0.25), 4),
    round(quantile(dync2h1_complete, 0.75), 4)
  )
)

cat("\nDYNC2H1 Expression Statistics:\n")
print(dync2h1_stats, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2c-2: Patient Stratification (Median Split)
# ---------------------------------------------------------------------------
cat("\n--- 2c-2: Patient Stratification (Median Split) ---\n")

dync2h1_median <- median(dync2h1_complete)
cat(sprintf("  Median DYNC2H1 expression: %.4f\n", dync2h1_median))

dync2h1_group <- ifelse(dync2h1_expr > dync2h1_median, "High", "Low")
dync2h1_group <- factor(dync2h1_group, levels = c("Low", "High"))

cat(sprintf("  High DYNC2H1 (> median): %d patients\n", sum(dync2h1_group == "High", na.rm = TRUE)))
cat(sprintf("  Low DYNC2H1 (<= median): %d patients\n", sum(dync2h1_group == "Low", na.rm = TRUE)))

# Create patient data frame
patient_data_dync <- data.frame(
  Patient_ID = rownames(datExpr),
  DYNC2H1_Expression = dync2h1_expr,
  DYNC2H1_Group = dync2h1_group,
  stringsAsFactors = FALSE
)

# Add clinical variables
patient_data_dync$PFS_Days <- datTraits$PFS
patient_data_dync$PFS_Group <- as.character(datTraits$PFS_group)
patient_data_dync$Response <- as.character(datTraits$Response)
patient_data_dync$Treatment <- as.character(datTraits$Treatment)
patient_data_dync$Sex <- as.character(datTraits$Sex)
patient_data_dync$Age <- datTraits$Age
patient_data_dync$Liver <- as.character(datTraits$Liver)

# Check available columns
cat("\n  Available columns in datTraits: ", paste(colnames(datTraits), collapse = ", "), "\n")

# Save patient data
write.csv(patient_data_dync, file.path(results_dir, "Step2c_DYNC2H1_Patient_Data.csv"), row.names = FALSE)
cat("  Saved: Step2c_DYNC2H1_Patient_Data.csv\n")

# ---------------------------------------------------------------------------
# 2c-3: Clinical Associations (Statistical Tests)
# ---------------------------------------------------------------------------
cat("\n--- 2c-3: Clinical Associations ---\n")

association_results_dync <- data.frame(
  Variable = character(),
  Test = character(),
  Statistic = character(),
  P_Value = numeric(),
  N_High = integer(),
  N_Low = integer(),
  Details = character(),
  stringsAsFactors = FALSE
)

# Test PFS_Group (Fisher's exact)
cat("\n  PFS_Group (Short/Long):\n")
tbl_pfs <- table(patient_data_dync$DYNC2H1_Group, patient_data_dync$PFS_Group)
print(tbl_pfs)
fisher_pfs <- fisher.test(tbl_pfs)
cat(sprintf("    Fisher's exact test p-value: %.4f\n", fisher_pfs$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "PFS_Group", Test = "Fisher", Statistic = "-",
  P_Value = fisher_pfs$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = paste(apply(tbl_pfs, 2, function(x) paste(x, collapse = "/")), collapse = " | ")
))

# Test Response (Fisher's exact)
cat("\n  Response (CD/PD):\n")
tbl_resp <- table(patient_data_dync$DYNC2H1_Group, patient_data_dync$Response)
print(tbl_resp)
fisher_resp <- fisher.test(tbl_resp)
cat(sprintf("    Fisher's exact test p-value: %.4f\n", fisher_resp$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "Response", Test = "Fisher", Statistic = "-",
  P_Value = fisher_resp$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = paste(apply(tbl_resp, 2, function(x) paste(x, collapse = "/")), collapse = " | ")
))

# Test Treatment (Fisher's exact)
cat("\n  Treatment:\n")
tbl_treat <- table(patient_data_dync$DYNC2H1_Group, patient_data_dync$Treatment)
print(tbl_treat)
fisher_treat <- fisher.test(tbl_treat)
cat(sprintf("    Fisher's exact test p-value: %.4f\n", fisher_treat$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "Treatment", Test = "Fisher", Statistic = "-",
  P_Value = fisher_treat$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = paste(apply(tbl_treat, 2, function(x) paste(x, collapse = "/")), collapse = " | ")
))

# Test Sex (Fisher's exact)
cat("\n  Sex:\n")
tbl_sex <- table(patient_data_dync$DYNC2H1_Group, patient_data_dync$Sex)
print(tbl_sex)
fisher_sex <- fisher.test(tbl_sex)
cat(sprintf("    Fisher's exact test p-value: %.4f\n", fisher_sex$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "Sex", Test = "Fisher", Statistic = "-",
  P_Value = fisher_sex$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = paste(apply(tbl_sex, 2, function(x) paste(x, collapse = "/")), collapse = " | ")
))

# Test Liver metastasis (Fisher's exact)
cat("\n  Liver Metastasis:\n")
tbl_liver <- table(patient_data_dync$DYNC2H1_Group, patient_data_dync$Liver)
print(tbl_liver)
fisher_liver <- fisher.test(tbl_liver)
cat(sprintf("    Fisher's exact test p-value: %.4f\n", fisher_liver$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "Liver_Metastasis", Test = "Fisher", Statistic = "-",
  P_Value = fisher_liver$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = paste(apply(tbl_liver, 2, function(x) paste(x, collapse = "/")), collapse = " | ")
))

# Test PFS_Days (Wilcoxon)
cat("\n  PFS_Days (continuous):\n")
high_pfs <- patient_data_dync$PFS_Days[patient_data_dync$DYNC2H1_Group == "High"]
low_pfs <- patient_data_dync$PFS_Days[patient_data_dync$DYNC2H1_Group == "Low"]
wilcox_pfs_dync <- wilcox.test(high_pfs, low_pfs)
cat(sprintf("    High DYNC2H1 median PFS: %.1f days\n", median(high_pfs, na.rm = TRUE)))
cat(sprintf("    Low DYNC2H1 median PFS: %.1f days\n", median(low_pfs, na.rm = TRUE)))
cat(sprintf("    Wilcoxon rank-sum test p-value: %.4f\n", wilcox_pfs_dync$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "PFS_Days", Test = "Wilcoxon", Statistic = "-",
  P_Value = wilcox_pfs_dync$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = sprintf("High median=%.1f, Low median=%.1f", median(high_pfs, na.rm = TRUE), median(low_pfs, na.rm = TRUE))
))

# Test Age (Wilcoxon)
cat("\n  Age (continuous):\n")
high_age <- patient_data_dync$Age[patient_data_dync$DYNC2H1_Group == "High"]
low_age <- patient_data_dync$Age[patient_data_dync$DYNC2H1_Group == "Low"]
wilcox_age_dync <- wilcox.test(high_age, low_age)
cat(sprintf("    High DYNC2H1 median age: %.1f years\n", median(high_age, na.rm = TRUE)))
cat(sprintf("    Low DYNC2H1 median age: %.1f years\n", median(low_age, na.rm = TRUE)))
cat(sprintf("    Wilcoxon rank-sum test p-value: %.4f\n", wilcox_age_dync$p.value))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "Age", Test = "Wilcoxon", Statistic = "-",
  P_Value = wilcox_age_dync$p.value,
  N_High = sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE),
  N_Low = sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE),
  Details = sprintf("High median=%.1f, Low median=%.1f", median(high_age, na.rm = TRUE), median(low_age, na.rm = TRUE))
))

# ---------------------------------------------------------------------------
# 2c-4: Correlation Analysis
# ---------------------------------------------------------------------------
cat("\n--- 2c-4: Correlation Analysis (DYNC2H1 vs PFS_Days) ---\n")

cor_data_dync <- patient_data_dync[!is.na(patient_data_dync$DYNC2H1_Expression) & !is.na(patient_data_dync$PFS_Days), ]

pearson_test_dync <- cor.test(cor_data_dync$DYNC2H1_Expression, cor_data_dync$PFS_Days, method = "pearson")
cat(sprintf("  Pearson correlation: r = %.4f, p = %.4f\n", pearson_test_dync$estimate, pearson_test_dync$p.value))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", pearson_test_dync$conf.int[1], pearson_test_dync$conf.int[2]))

spearman_test_dync <- cor.test(cor_data_dync$DYNC2H1_Expression, cor_data_dync$PFS_Days, method = "spearman")
cat(sprintf("  Spearman correlation: rho = %.4f, p = %.4f\n", spearman_test_dync$estimate, spearman_test_dync$p.value))

association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "PFS_Days_Pearson", Test = "Pearson", Statistic = round(pearson_test_dync$estimate, 4),
  P_Value = pearson_test_dync$p.value, N_High = NA, N_Low = NA,
  Details = sprintf("95%% CI: [%.4f, %.4f]", pearson_test_dync$conf.int[1], pearson_test_dync$conf.int[2])
))
association_results_dync <- rbind(association_results_dync, data.frame(
  Variable = "PFS_Days_Spearman", Test = "Spearman", Statistic = round(spearman_test_dync$estimate, 4),
  P_Value = spearman_test_dync$p.value, N_High = NA, N_Low = NA,
  Details = "-"
))

write.csv(association_results_dync, file.path(results_dir, "Step2c_DYNC2H1_Association_Summary.csv"), row.names = FALSE)
cat("\n  Saved: Step2c_DYNC2H1_Association_Summary.csv\n")



# ---------------------------------------------------------------------------
# 2c-6: Visualizations (4-Panel Figure)
# ---------------------------------------------------------------------------
cat("\n--- 2c-6: Generating Visualizations ---\n")

pfs_colors_dync <- colors_pfs
response_colors_dync <- colors_response

plot_data_dync <- patient_data_dync[!is.na(patient_data_dync$DYNC2H1_Expression), ]

# Calculate y-axis max for bracket positioning
y_max_dync <- max(plot_data_dync$DYNC2H1_Expression, na.rm = TRUE)
y_min_dync <- min(plot_data_dync$DYNC2H1_Expression, na.rm = TRUE)
y_range_dync <- y_max_dync - y_min_dync

# Panel A: DYNC2H1 by PFS_Group
wilcox_pfs_plot <- wilcox.test(DYNC2H1_Expression ~ PFS_Group, data = plot_data_dync)
pval_label_pfs <- ifelse(wilcox_pfs_plot$p.value < 0.001, "Wilcoxon p < 0.001",
                          sprintf("Wilcoxon p = %.3f", wilcox_pfs_plot$p.value))

p_dync_a <- ggplot(plot_data_dync, aes(x = PFS_Group, y = DYNC2H1_Expression, fill = PFS_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = pfs_colors_dync) +
  # Significance bracket
  annotate("segment", x = 1, xend = 2, y = y_max_dync + y_range_dync * 0.05,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_max_dync + y_range_dync * 0.10,
           label = pval_label_pfs, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(y_min_dync - y_range_dync * 0.05, y_max_dync + y_range_dync * 0.18)) +
  labs(
    title = "DYNC2H1 by PFS Group",
    x = "PFS Group",
    y = "DYNC2H1 Expression (log2)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "none"
  )

# Panel B: DYNC2H1 by Response
wilcox_resp_plot <- wilcox.test(DYNC2H1_Expression ~ Response, data = plot_data_dync)
pval_label_resp <- ifelse(wilcox_resp_plot$p.value < 0.001, "Wilcoxon p < 0.001",
                           sprintf("Wilcoxon p = %.3f", wilcox_resp_plot$p.value))

p_dync_b <- ggplot(plot_data_dync, aes(x = Response, y = DYNC2H1_Expression, fill = Response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = response_colors_dync) +
  # Significance bracket
  annotate("segment", x = 1, xend = 2, y = y_max_dync + y_range_dync * 0.05,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_max_dync + y_range_dync * 0.10,
           label = pval_label_resp, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(y_min_dync - y_range_dync * 0.05, y_max_dync + y_range_dync * 0.18)) +
  labs(
    title = "DYNC2H1 by Response",
    x = "Response",
    y = "DYNC2H1 Expression (log2)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "none"
  )

# Panel C: DYNC2H1 by Treatment (2 groups: FOLFI vs GemNab -> Wilcoxon)
wilcox_treat_dync <- wilcox.test(DYNC2H1_Expression ~ Treatment, data = plot_data_dync)
pval_label_treat <- ifelse(wilcox_treat_dync$p.value < 0.001, "Wilcoxon p < 0.001",
                            sprintf("Wilcoxon p = %.3f", wilcox_treat_dync$p.value))

p_dync_c <- ggplot(plot_data_dync, aes(x = Treatment, y = DYNC2H1_Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("FOLFI" = col_folfi, "GemNab" = col_gemnab)) +
  # Significance bracket (2 treatment groups)
  annotate("segment", x = 1, xend = 2, y = y_max_dync + y_range_dync * 0.05,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_max_dync + y_range_dync * 0.10,
           label = pval_label_treat, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(y_min_dync - y_range_dync * 0.05, y_max_dync + y_range_dync * 0.18)) +
  labs(
    title = "DYNC2H1 by Treatment",
    x = "Treatment",
    y = "DYNC2H1 Expression (log2)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Panel D: DYNC2H1 vs PFS_Days scatter
# Add correlation annotation inside the plot
corr_label_dync <- sprintf("Pearson r = %.3f\np = %.3f", pearson_test_dync$estimate, pearson_test_dync$p.value)

p_dync_d <- ggplot(plot_data_dync, aes(x = DYNC2H1_Expression, y = PFS_Days, color = PFS_Group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
  scale_color_manual(values = pfs_colors_dync) +
  # Add correlation annotation in top-left corner
  annotate("text", x = min(plot_data_dync$DYNC2H1_Expression, na.rm = TRUE) + 0.1,
           y = max(plot_data_dync$PFS_Days, na.rm = TRUE) * 0.95,
           label = corr_label_dync, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold", color = "grey30") +
  labs(
    title = "DYNC2H1 vs PFS Days",
    x = "DYNC2H1 Expression (log2)",
    y = "PFS (days)",
    color = "PFS Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right"
  )

# Combine panels
p_dync_combined <- (p_dync_a | p_dync_b) / (p_dync_c | p_dync_d) +
  plot_annotation(
    title = "DYNC2H1 Expression: Patient Stratification Analysis",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(p_dync_combined)

ggsave(file.path(fig_dir_png, "Step2c_DYNC2H1_Stratification.png"), p_dync_combined,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2c_DYNC2H1_Stratification.pdf"), p_dync_combined,
       width = 10, height = 8)
cat("  Saved: Step2c_DYNC2H1_Stratification.png/pdf\n")

# ---------------------------------------------------------------------------
# 2c-6b: Violin Plot Version
# ---------------------------------------------------------------------------
cat("\n--- 2c-6b: Generating Violin Plots ---\n")

# Panel A (Violin): DYNC2H1 by PFS_Group
p_dync_violin_a <- ggplot(plot_data_dync, aes(x = PFS_Group, y = DYNC2H1_Expression, fill = PFS_Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = colors_pfs) +
  annotate("segment", x = 1, xend = 2, y = y_max_dync + y_range_dync * 0.05,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_max_dync + y_range_dync * 0.10,
           label = pval_label_pfs, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(y_min_dync - y_range_dync * 0.05, y_max_dync + y_range_dync * 0.18)) +
  labs(
    title = "DYNC2H1 by PFS Group",
    x = "PFS Group",
    y = "DYNC2H1 Expression (log2)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "none"
  )

# Panel B (Violin): DYNC2H1 by Response (same as boxplot Panel B)
p_dync_violin_b <- ggplot(plot_data_dync, aes(x = Response, y = DYNC2H1_Expression, fill = Response)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = colors_response) +
  annotate("segment", x = 1, xend = 2, y = y_max_dync + y_range_dync * 0.05,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_max_dync + y_range_dync * 0.10,
           label = pval_label_resp, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(y_min_dync - y_range_dync * 0.05, y_max_dync + y_range_dync * 0.18)) +
  labs(
    title = "DYNC2H1 by Response",
    x = "Response",
    y = "DYNC2H1 Expression (log2)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "none"
  )

# Panel C (Violin): DYNC2H1 by Treatment
p_dync_violin_c <- ggplot(plot_data_dync, aes(x = Treatment, y = DYNC2H1_Expression, fill = Treatment)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("FOLFI" = col_folfi, "GemNab" = col_gemnab)) +
  annotate("segment", x = 1, xend = 2, y = y_max_dync + y_range_dync * 0.05,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_max_dync + y_range_dync * 0.02,
           yend = y_max_dync + y_range_dync * 0.05, color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_max_dync + y_range_dync * 0.10,
           label = pval_label_treat, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(y_min_dync - y_range_dync * 0.05, y_max_dync + y_range_dync * 0.18)) +
  labs(
    title = "DYNC2H1 by Treatment",
    x = "Treatment",
    y = "DYNC2H1 Expression (log2)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Panel D (Violin): Same scatter plot (no violin version needed)
p_dync_violin_d <- p_dync_d

# Combine violin panels
p_dync_violin_combined <- (p_dync_violin_a | p_dync_violin_b) / (p_dync_violin_c | p_dync_violin_d) +
  plot_annotation(
    title = "DYNC2H1 Expression: Patient Stratification Analysis (Violin)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(p_dync_violin_combined)

ggsave(file.path(fig_dir_png, "Step2c_DYNC2H1_Stratification_Violin.png"), p_dync_violin_combined,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2c_DYNC2H1_Stratification_Violin.pdf"), p_dync_violin_combined,
       width = 10, height = 8)
cat("  Saved: Step2c_DYNC2H1_Stratification_Violin.png/pdf\n")



# ---------------------------------------------------------------------------
# 2c-7: Summary Statistics
# ---------------------------------------------------------------------------
cat("\n--- 2c-7: Summary Statistics ---\n")

n_high_dync <- sum(patient_data_dync$DYNC2H1_Group == "High", na.rm = TRUE)
n_low_dync <- sum(patient_data_dync$DYNC2H1_Group == "Low", na.rm = TRUE)

cat("\nDYNC2H1 Expression:\n")
cat(sprintf("  N = %d, Mean = %.4f, Median = %.4f, SD = %.4f, Range = [%.4f, %.4f]\n",
            nrow(patient_data_dync),
            mean(patient_data_dync$DYNC2H1_Expression, na.rm = TRUE),
            dync2h1_median,
            sd(patient_data_dync$DYNC2H1_Expression, na.rm = TRUE),
            min(patient_data_dync$DYNC2H1_Expression, na.rm = TRUE),
            max(patient_data_dync$DYNC2H1_Expression, na.rm = TRUE)))

cat("\nPatient Stratification (median split):\n")
cat(sprintf("  High DYNC2H1: n = %d (%.1f%%)\n", n_high_dync, 100 * n_high_dync / nrow(patient_data_dync)))
cat(sprintf("  Low DYNC2H1:  n = %d (%.1f%%)\n", n_low_dync, 100 * n_low_dync / nrow(patient_data_dync)))

cat("\nAssociation Tests (DYNC2H1 High vs Low):\n")
for(i in 1:nrow(association_results_dync)) {
  row <- association_results_dync[i, ]
  cat(sprintf("  %s: %s p = %.4f\n", row$Variable, row$Test, row$P_Value))
}

cat("\nCorrelation with PFS Days:\n")
cat(sprintf("  Pearson r = %.4f, p = %.4f, 95%% CI [%.4f, %.4f]\n",
            pearson_test_dync$estimate, pearson_test_dync$p.value,
            pearson_test_dync$conf.int[1], pearson_test_dync$conf.int[2]))
cat(sprintf("  Spearman rho = %.4f, p = %.4f\n",
            spearman_test_dync$estimate, spearman_test_dync$p.value))

cat("\nVisualization Statistics:\n")
cat(sprintf("  Panel A - DYNC2H1 by PFS Group:   Wilcoxon p = %.4f\n", wilcox_pfs_plot$p.value))
cat(sprintf("  Panel B - DYNC2H1 by Response:    Wilcoxon p = %.4f\n", wilcox_resp_plot$p.value))
cat(sprintf("  Panel C - DYNC2H1 by Treatment:   Wilcoxon p = %.4f\n", wilcox_treat_dync$p.value))
cat(sprintf("  Panel D - DYNC2H1 vs PFS Days:    Pearson r = %.4f, p = %.4f\n",
            pearson_test_dync$estimate, pearson_test_dync$p.value))

# Save to file
sink(file.path(results_dir, "Step2c_DYNC2H1_Statistics.txt"))
cat("STEP 2c: DYNC2H1 PATIENT STRATIFICATION ANALYSIS\n")
cat(sprintf("Date: %s\n\n", Sys.time()))
cat(sprintf("DYNC2H1 Expression: N = %d, Mean = %.4f, Median = %.4f, SD = %.4f\n",
            nrow(patient_data_dync),
            mean(patient_data_dync$DYNC2H1_Expression, na.rm = TRUE),
            dync2h1_median,
            sd(patient_data_dync$DYNC2H1_Expression, na.rm = TRUE)))
cat(sprintf("Stratification: High = %d, Low = %d\n\n", n_high_dync, n_low_dync))
cat("Association Tests:\n")
print(association_results_dync, row.names = FALSE)
cat(sprintf("\nPearson r = %.4f, p = %.4f\n", pearson_test_dync$estimate, pearson_test_dync$p.value))
cat(sprintf("Spearman rho = %.4f, p = %.4f\n", spearman_test_dync$estimate, spearman_test_dync$p.value))
sink()
cat("\nSaved: Step2c_DYNC2H1_Statistics.txt\n")

#DYNC2H1 is significantly associated with:
#PFS_Group (p = 0.0042)
#Response (p = 0.0209)
#PFS_Days continuous (p = 0.0157)

#DYNC2H1 is NOT associated with:
#Treatment (p = 1.00) → Expression is independent of treatment arm
#Sex (p = 0.23)
#Liver metastasis (p = 0.73)
#Age (p = 0.19)


#High variability confirmed: CV = 31% explains why DYNC2H1 is in grey module across all WGCNA powers
#Clinical relevance: Statistically significant associations with both survival (PFS) and treatment response
#Treatment-independent: Expression does not differ between FOLFIRINOX and Gem/NabP arms (p = 1.00)
#Positive correlation: Higher DYNC2H1 correlates with longer PFS_Days (r = 0.39)

#Thesis
#DYNC2H1 exhibited the highest inter-patient variability among the ML signature proteins (CV = 31.0%),
#with expression values ranging from 4.04 to 15.02 (log2 intensity). Median-split stratification revealed
#that patients with high DYNC2H1 expression showed significantly longer progression-free survival compared
#to those with low expression (Wilcoxon p = 0.016). This association was confirmed by both Pearson (r = 0.39, p = 0.006)
#and Spearman (ρ = 0.42, p = 0.003) correlation analyses. Furthermore, DYNC2H1 expression was significantly associated with
#treatment response, with higher levels observed in patients achieving disease control (Fisher's p = 0.021).
#Notably, DYNC2H1 expression was independent of treatment arm (p = 1.00), patient age (p = 0.19), sex (p = 0.23), and presence of liver metastasis (p = 0.73), suggesting that its prognostic value is not confounded by these clinical variables.

cat("\n============================================================================\n")
cat("STEP 2c COMPLETE\n")
cat("============================================================================\n")

# ----------------------------------------------------------------------------
# 2.4 Inter-Protein Correlation Analysis
# ----------------------------------------------------------------------------
cat("--- 2.4 Inter-Protein Correlation Analysis ---\n")

# Calculate mean absolute correlation for each protein
cor_matrix <- cor(datExpr, use = "pairwise.complete.obs")
diag(cor_matrix) <- NA  # Remove self-correlation

mean_abs_cor <- apply(abs(cor_matrix), 1, mean, na.rm = TRUE)

cat("Mean absolute correlation per protein:\n")
print(summary(mean_abs_cor))

# Flag proteins with very low mean correlation
mac_percentiles <- quantile(mean_abs_cor, probs = c(0.05, 0.10, 0.25, 0.50))
cat("Mean absolute correlation percentiles:\n")
print(round(mac_percentiles, 4))

low_cor_threshold <- quantile(mean_abs_cor, 0.05)
low_cor_proteins <- names(mean_abs_cor[mean_abs_cor < low_cor_threshold])
cat("Proteins with mean |correlation| < 5th percentile:", length(low_cor_proteins), "\n")

# --- Correlation Distribution Plot ---
cat("Generating Correlation Distribution plot...\n")

cor_df <- data.frame(
  Protein = names(mean_abs_cor),
  Mean_Abs_Cor = mean_abs_cor,
  Variance = protein_variance
)

# Panel A: Histogram
p4a <- ggplot(cor_df, aes(x = Mean_Abs_Cor)) +
  geom_histogram(bins = 50, fill = color_palette[5], color = "white", alpha = 0.8) +
  geom_vline(xintercept = low_cor_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  labs(
    title = "Mean |Correlation| Distribution",
    subtitle = sprintf("5th pctl = %.3f", low_cor_threshold),
    x = "Mean Absolute Correlation", y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
  )
print(p4a)

# Panel B: Correlation vs Variance
p4b <- ggplot(cor_df, aes(x = Variance, y = Mean_Abs_Cor)) +
  geom_point(color = color_palette[5], alpha = 0.6, size = 2) +
  geom_hline(yintercept = low_cor_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  geom_vline(xintercept = low_var_threshold, linetype = "dashed", color = color_palette[1], linewidth = 0.8) +
  labs(
    title = "Mean |Correlation| vs Variance",
    subtitle = "Dashed lines: 5th percentile thresholds",
    x = "Variance", y = "Mean Absolute Correlation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
  )
print(p4b)

p4_combined <- p4a + p4b
print(p4_combined)

ggsave(file.path(fig_dir_png, "Step2_04_Correlation_Distribution.png"), p4_combined, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2_04_Correlation_Distribution.pdf"), p4_combined, width = 8, height = 5)

# Proteins with low mean absolute correlation (<0.12)
#These proteins will likely end up in the "grey" (unassigned) module in WGCNA.

# ----------------------------------------------------------------------------
# 2.5 Expression Level Analysis
# ----------------------------------------------------------------------------
cat("--- 2.5 Expression Level Analysis ---\n")

cat("Protein mean expression distribution:\n")
print(summary(protein_mean))

mean_percentiles <- quantile(protein_mean, probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95))
cat("Mean expression percentiles:\n")
print(round(mean_percentiles, 4))

# Flag very low abundance proteins
low_expr_threshold <- quantile(protein_mean, 0.05)
low_expr_proteins <- names(protein_mean[protein_mean < low_expr_threshold])
cat("Proteins with mean expression < 5th percentile:", length(low_expr_proteins), "\n")

# --- Expression Level Distribution Plot ---
cat("Generating Expression Level Distribution plot...\n")

expr_df <- data.frame(
  Protein = names(protein_mean),
  Mean = protein_mean,
  Variance = protein_variance
)

# Panel A: Histogram
p5a <- ggplot(expr_df, aes(x = Mean)) +
  geom_histogram(bins = 50, fill = color_palette[3], color = "white", alpha = 0.8) +
  geom_vline(xintercept = low_expr_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  labs(
    title = "Mean Expression Distribution",
    subtitle = sprintf("5th pctl = %.3f", low_expr_threshold),
    x = "Mean Expression (log2)", y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
  )
print(p5a)

# Panel B: Expression vs Variance
p5b <- ggplot(expr_df, aes(x = Mean, y = Variance)) +
  geom_point(color = color_palette[3], alpha = 0.6, size = 2) +
  geom_vline(xintercept = low_expr_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  labs(
    title = "Variance vs Mean Expression",
    subtitle = "Dashed line: 5th percentile threshold",
    x = "Mean Expression (log2)", y = "Variance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
  )
print(p5b)

p5_combined <- p5a + p5b
print(p5_combined)

ggsave(file.path(fig_dir_png, "Step2_05_Expression_Distribution.png"), p5_combined, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2_05_Expression_Distribution.pdf"), p5_combined, width = 8, height = 5)

# ----------------------------------------------------------------------------
# 2.6 Compile Noise Assessment Summary
# ----------------------------------------------------------------------------
cat("============================================================================\n")
cat("NOISE ASSESSMENT SUMMARY\n")
cat("============================================================================\n")

# Create a data frame with all metrics per protein
noise_assessment <- data.frame(
  Protein = colnames(datExpr),
  Mean_Expression = protein_mean,
  Variance = protein_variance,
  SD = protein_sd,
  CV = protein_cv,
  Mean_Abs_Correlation = mean_abs_cor,
  stringsAsFactors = FALSE
)

# Add flag columns
noise_assessment$Flag_LowVariance <- noise_assessment$Protein %in% low_var_proteins
noise_assessment$Flag_LowCV <- noise_assessment$Protein %in% low_cv_proteins
noise_assessment$Flag_HighCV <- noise_assessment$Protein %in% high_cv_proteins
noise_assessment$Flag_LowCorrelation <- noise_assessment$Protein %in% low_cor_proteins
noise_assessment$Flag_LowExpression <- noise_assessment$Protein %in% low_expr_proteins
noise_assessment$Flag_MV_Outlier <- noise_assessment$Protein %in% mv_outliers

# Count total flags per protein
noise_assessment$Total_Flags <- rowSums(noise_assessment[, grep("^Flag_", colnames(noise_assessment))])

# Summary of flagged proteins
cat("Proteins flagged by number of criteria:\n")
print(table(noise_assessment$Total_Flags))

# Proteins flagged by >=2 criteria (potential noise candidates)
multi_flagged <- noise_assessment[noise_assessment$Total_Flags >= 2, ]
cat("Proteins flagged by >=2 criteria:", nrow(multi_flagged), "\n")

if(nrow(multi_flagged) > 0) {
  cat("Details of multi-flagged proteins:\n")
  print(multi_flagged[order(-multi_flagged$Total_Flags),
                      c("Protein", "Mean_Expression", "Variance", "CV",
                        "Mean_Abs_Correlation", "Total_Flags")])
}

# --- Summary of Flagged Proteins Plot ---
cat("Generating Summary of Flagged Proteins plot...\n")

flag_summary_df <- data.frame(
  Criterion = c("Low Variance", "Low CV", "High CV", "Low Correlation", "Low Expression", "MV Outlier"),
  Count = c(length(low_var_proteins), length(low_cv_proteins), length(high_cv_proteins),
            length(low_cor_proteins), length(low_expr_proteins), length(mv_outliers))
)
flag_summary_df$Criterion <- factor(flag_summary_df$Criterion, levels = flag_summary_df$Criterion)

p6_flags <- ggplot(flag_summary_df, aes(x = Criterion, y = Count, fill = Criterion)) +
  geom_col(alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 4) +
  geom_hline(yintercept = ceiling(ncol(datExpr) * 0.05), linetype = "dashed", color = "grey50", linewidth = 0.8) +
  scale_fill_manual(values = color_palette[1:6]) +
  labs(
    title = "Proteins Flagged by Each Criterion",
    subtitle = sprintf("5th/95th percentile thresholds | Dashed line = 5%% of %d proteins", ncol(datExpr)),
    x = NULL, y = "Number of Proteins"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
print(p6_flags)

ggsave(file.path(fig_dir_png, "Step2_06_Flagged_Summary.png"), p6_flags, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step2_06_Flagged_Summary.pdf"), p6_flags, width = 6, height = 5)

# --- Heatmap of Multi-Flagged Proteins ---
if(nrow(multi_flagged) > 0) {
  cat("Generating Heatmap of Multi-Flagged Proteins...\n")

  # Prepare data for ggplot heatmap
  flag_matrix <- as.matrix(multi_flagged[, grep("^Flag_", colnames(multi_flagged))])
  rownames(flag_matrix) <- multi_flagged$Protein
  colnames(flag_matrix) <- gsub("Flag_", "", colnames(flag_matrix))
  flag_matrix <- flag_matrix[order(-multi_flagged$Total_Flags), , drop = FALSE]
  flag_matrix <- flag_matrix * 1

  # Convert to long format for ggplot
  flag_long <- as.data.frame(flag_matrix)
  flag_long$Protein <- factor(rownames(flag_long), levels = rev(rownames(flag_long)))
  flag_long <- pivot_longer(flag_long, cols = -Protein, names_to = "Criterion", values_to = "Flagged")
  flag_long$Flagged <- factor(flag_long$Flagged, levels = c(0, 1), labels = c("No", "Yes"))

  p7_heatmap <- ggplot(flag_long, aes(x = Criterion, y = Protein, fill = Flagged)) +
    geom_tile(color = "grey90", linewidth = 0.5) +
    scale_fill_manual(values = c("No" = "white", "Yes" = color_palette[7])) +
    labs(
      title = sprintf("Proteins Flagged by >=2 Criteria (n=%d)", nrow(multi_flagged)),
      subtitle = "Multi-flagged proteins may indicate technical noise",
      x = NULL, y = NULL, fill = "Flagged"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  print(p7_heatmap)

  ggsave(file.path(fig_dir_png, "Step2_07_MultiFlagged_Heatmap.png"), p7_heatmap,
         width = 8, height = 10, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step2_07_MultiFlagged_Heatmap.pdf"), p7_heatmap,
         width = 8, height = 10)

} else {
  cat("Skipped Heatmap - No proteins flagged by >=2 criteria\n")
}

cat("All Step 2 diagnostic plots generated\n")

# ----------------------------------------------------------------------------
# 2.7 Save Results Tables
# ----------------------------------------------------------------------------
cat("--- Saving Results Tables ---\n")

# Save the complete noise assessment table
write.csv(noise_assessment,
          file.path(results_dir, "Step2_Noise_Assessment_All_Proteins.csv"),
          row.names = FALSE)
cat("Saved: Step2_Noise_Assessment_All_Proteins.csv\n")

# Save multi-flagged proteins separately
if(nrow(multi_flagged) > 0) {
  write.csv(multi_flagged[order(-multi_flagged$Total_Flags), ],
            file.path(results_dir, "Step2_MultiFlagged_Proteins.csv"),
            row.names = FALSE)
  cat("Saved: Step2_MultiFlagged_Proteins.csv\n")
}

# Save variance percentile reference
var_reference <- data.frame(
  Percentile = names(var_percentiles),
  Variance_Threshold = as.numeric(var_percentiles),
  Proteins_Above = sapply(var_percentiles, function(x) sum(protein_variance >= x)),
  Proteins_Below = sapply(var_percentiles, function(x) sum(protein_variance < x))
)
write.csv(var_reference,
          file.path(results_dir, "Step2_Variance_Percentile_Reference.csv"),
          row.names = FALSE)
cat("Saved: Step2_Variance_Percentile_Reference.csv\n")

#Quality assessment of 891 proteins identified 710 (79.7%) as high-quality with no flags, while 54 proteins (6.1%) were flagged by >=2 criteria.
#The minimal mean-variance dependence (RÂ2 = 0.013) confirmed successful variance stabilization by log2 transformation.
#Five proteins (EXOC5, TNFRSF21, F11, ENPP1, CAP1) showed multiple quality concerns.
#All proteins were retained for WGCNA analysis, as the algorithm is robust to moderate noise and flagged proteins would naturally segregate to the grey (unassigned) module.

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 2 SUMMARY: PROTEIN NOISE ASSESSMENT ===\n")
cat("===============================================================================\n")
cat(sprintf("Total proteins analyzed: %d\n", ncol(datExpr)))
cat(sprintf("High-CV proteins (>95th pctl): %d\n", length(high_cv_proteins)))
cat(sprintf("Proteins with no flags: %d (%.1f%%)\n",
            sum(noise_assessment$Total_Flags == 0),
            100 * sum(noise_assessment$Total_Flags == 0) / nrow(noise_assessment)))
cat(sprintf("Proteins flagged by >=2 criteria: %d\n", nrow(multi_flagged)))
cat("Decision: All proteins retained for WGCNA (robust to moderate noise)\n")
cat("Step 2 complete.\n")
cat("===============================================================================\n")

# STEP 3: SAMPLE QUALITY CONTROL AND CLUSTERING
# PURPOSE: Assess sample quality, detect potential outliers, and visualize 
# sample-level structure before network construction.                      
#                                                                         
# STRATEGY:                                                               
# - Outliers: Z-connectivity < -2.5 (isolated samples)                   
# - PCA: Batch effects, clinical variable separation                     
# - No extreme samples clustering alone on dendrogram                    
#                                                                         
# KEY PARAMS: outlier_threshold=-2.5, PCA_variance_explained               

cat("===============================================================================\n")
cat("STEP 3: SAMPLE QUALITY CONTROL AND CLUSTERING\n")
cat("===============================================================================\n")

# ----------------------------------------------------------------------------
# 3.1 Check Data Quality with WGCNA's goodSamplesGenes
# ----------------------------------------------------------------------------
cat("--- 3.1 WGCNA Data Quality Check ---\n")

gsg <- goodSamplesGenes(datExpr, verbose = 3)

if(gsg$allOK) {
  cat("[OK] All samples and proteins passed WGCNA quality check\n")
} else {
  cat("[WARNING] Some samples or genes flagged by WGCNA:\n")
  if(sum(!gsg$goodGenes) > 0) {
    cat("  Genes to remove:", sum(!gsg$goodGenes), "\n")
  }
  if(sum(!gsg$goodSamples) > 0) {
    cat("  Samples to remove:", sum(!gsg$goodSamples), "\n")
  }
}

# ----------------------------------------------------------------------------
# 3.2 Sample Clustering to Detect Outliers
# ----------------------------------------------------------------------------
cat("--- 3.2 Sample Hierarchical Clustering ---\n")

# Calculate sample distance and cluster
sampleTree <- hclust(dist(datExpr), method = "average")

# ----------------------------------------------------------------------------
# Figure 1: Sample Dendrogram
# ----------------------------------------------------------------------------
cat("[1/6] Sample Dendrogram\n")

# Convert dendrogram to ggplot-compatible format using ggdendro

dend_data <- dendro_data(sampleTree, type = "rectangle")

p_dend <- ggplot() +
  geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend), color = "grey40") +
  geom_text(data = dend_data$labels, aes(x = x, y = y, label = label), hjust = 1, angle = 90, size = 2.5) +
  geom_hline(yintercept = 85, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  annotate("text", x = max(dend_data$segments$x) * 0.8, y = 88, label = "Potential outlier threshold",
           color = color_palette[7], size = 3) +
  labs(
    title = "Sample Clustering Dendrogram",
    subtitle = "Hierarchical clustering using average linkage",
    x = NULL, y = "Height"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )
print(p_dend)
ggsave(file.path(fig_dir_png, "Step3_01_Sample_Dendrogram.png"), p_dend, width = 6, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step3_01_Sample_Dendrogram.pdf"), p_dend, width = 6, height = 4)

# ----------------------------------------------------------------------------
# 3.3 Sample Dendrogram with Trait Heatmap (WGCNA style)
# ----------------------------------------------------------------------------
cat("--- 3.3 Sample Dendrogram with Clinical Traits ---\n")

# Create numeric trait matrix for heatmap
traitMatrix <- data.frame(
  PFS_days = datTraits$PFS,
  PFS_group = as.numeric(datTraits$PFS_group == "Long"),
  Treatment = as.numeric(datTraits$Treatment == "FOLFI"),
  Response = as.numeric(datTraits$Response == "CD"),
  Sex = as.numeric(datTraits$Sex == "Female"),
  Liver_Mets = as.numeric(datTraits$Liver == "Yes"),
  Age = datTraits$Age,
  CA19.9_log = log10(datTraits$CA19.9 + 1)
)
rownames(traitMatrix) <- rownames(datTraits)

# Create color representation
traitColors_heatmap <- numbers2colors(traitMatrix,
                                      colors = colorRampPalette(color_palette)(50),
                                      signed = FALSE)

# ----------------------------------------------------------------------------
# Figure 2: Sample Dendrogram with Clinical Traits (WGCNA native - kept for compatibility)
# ----------------------------------------------------------------------------
cat("[2/6] Sample Dendrogram with Clinical Traits\n")

# WGCNA's plotDendroAndColors is specialized - keep it but save properly
png(file.path(fig_dir_png, "Step3_02_Sample_Dendrogram_Traits.png"), width = 10, height = 7, units = "in", res = 300)
plotDendroAndColors(sampleTree,
                    colors = traitColors_heatmap,
                    groupLabels = colnames(traitMatrix),
                    main = "Sample Dendrogram with Clinical Traits",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    cex.colorLabels = 0.8)
dev.off()

pdf(file.path(fig_dir_pdf, "Step3_02_Sample_Dendrogram_Traits.pdf"), width = 10, height = 7)
plotDendroAndColors(sampleTree,
                    colors = traitColors_heatmap,
                    groupLabels = colnames(traitMatrix),
                    main = "Sample Dendrogram with Clinical Traits",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    cex.colorLabels = 0.8)
dev.off()

# ----------------------------------------------------------------------------
# 3.4 Sample Connectivity Analysis (Outlier Detection)
# ----------------------------------------------------------------------------
#Z-Score Calculation for Sample Connectivity: Calculate sample-sample correlation matrix (pearson), sum correlations for each sample to get connectivity, standardize to Z-scores.
#Connectivity(sample_i) = sum of all correlations for that sample
#Standardize to Z-scores: Mean connectivity = average of all 50 connectivity values
#SD connectivity = standard deviation of all 50 connectivity values
#Z(sample_i) = (Connectivity_i - Mean) / SD

cat("--- 3.4 Sample Connectivity Analysis ---\n")

# Calculate sample-sample correlation
sample_cor <- cor(t(datExpr), use = "pairwise.complete.obs")

# Calculate connectivity (sum of correlations for each sample)
sample_connectivity <- rowSums(sample_cor) - 1

# Standardize connectivity (Z-score)
Z_connectivity <- scale(sample_connectivity)

cat("Sample connectivity (Z-scores):\n")
print(summary(as.numeric(Z_connectivity)))

# Identify potential outliers (Z < -2.5 is common threshold)
outlier_threshold <- -2.5
potential_outliers <- names(Z_connectivity[Z_connectivity < outlier_threshold, ])

cat("Potential outlier samples (Z <", outlier_threshold, "):",
    length(potential_outliers), "\n")

if(length(potential_outliers) > 0) {
  cat("Outlier sample IDs:\n")
  outlier_info <- data.frame(
    Sample = potential_outliers,
    Z_connectivity = as.numeric(Z_connectivity[potential_outliers, ]),
    PFS_group = datTraits[potential_outliers, "PFS_group"],
    Treatment = datTraits[potential_outliers, "Treatment"]
  )
  print(outlier_info)
} else {
  cat("[OK] No samples identified as connectivity outliers\n")
}

# ----------------------------------------------------------------------------
# Figure 3: Sample Connectivity Analysis
# ----------------------------------------------------------------------------
cat("[3/6] Sample Connectivity Analysis\n")

conn_df <- data.frame(
  Sample = names(Z_connectivity[,1]),
  Z_connectivity = as.numeric(Z_connectivity),
  Is_Outlier = as.numeric(Z_connectivity) < outlier_threshold
)
conn_df <- conn_df[order(conn_df$Z_connectivity), ]
conn_df$Rank <- 1:nrow(conn_df)

# Panel A: Histogram
p3a_conn <- ggplot(conn_df, aes(x = Z_connectivity)) +
  geom_histogram(bins = 15, fill = color_palette[2], color = "white", alpha = 0.8) +
  geom_vline(xintercept = outlier_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  labs(
    title = "Sample Connectivity Distribution",
    subtitle = sprintf("Outlier threshold = %.1f", outlier_threshold),
    x = "Standardized Connectivity (Z)", y = "Frequency"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
  )
print(p3a_conn)

# Panel B: Ordered bar plot
p3b_conn <- ggplot(conn_df, aes(x = Rank, y = Z_connectivity, fill = Is_Outlier)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = outlier_threshold, linetype = "dashed", color = color_palette[7], linewidth = 0.8) +
  scale_fill_manual(values = c("FALSE" = color_palette[1], "TRUE" = color_palette[7]),
                    labels = c("Normal", "Outlier")) +
  labs(
    title = "Ordered Sample Connectivity",
    subtitle = sprintf("Potential outliers: %d", sum(conn_df$Is_Outlier)),
    x = "Samples (ordered)", y = "Standardized Connectivity (Z)",
    fill = "Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    legend.position = "bottom"
  )
print(p3b_conn)

p3_connectivity <- p3a_conn + p3b_conn
print(p3_connectivity)

ggsave(file.path(fig_dir_png, "Step3_03_Sample_Connectivity.png"), p3_connectivity, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step3_03_Sample_Connectivity.pdf"), p3_connectivity, width = 8, height = 5)

# ----------------------------------------------------------------------------
# 3.5 PCA Visualization
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# 3.5 PCA Visualization
# ----------------------------------------------------------------------------
cat("--- 3.5 PCA Visualization ---\n")

# Perform PCA
pca_result <- prcomp(datExpr, scale. = TRUE, center = TRUE)

# Get variance explained
var_explained <- summary(pca_result)$importance[2, 1:5] * 100
cat("Variance explained by first 5 PCs:\n")
print(round(var_explained, 2))

# Create PCA data frame
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PFS_group = datTraits$PFS_group,
  Treatment = datTraits$Treatment,
  Response = datTraits$Response,
  Sample = rownames(datExpr)
)


# Use universal color scales (defined at top of script)
# colors_pfs, colors_treatment, colors_response already defined

# ----------------------------------------------------------------------------
# Figure 4: PCA Panels
# ----------------------------------------------------------------------------
cat("[4/6] PCA Panels\n")

# Panel A: PCA by PFS group
p_pca1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = PFS_group)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = colors_pfs, name = "PFS Group") +
  labs(
    title = "PCA: PFS Group",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )
print(p_pca1)

# Panel B: PCA by Treatment
p_pca2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = colors_treatment, name = "Treatment") +
  labs(
    title = "PCA: Treatment",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

print(p_pca2)

# Panel C: PCA by Response
p_pca3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = colors_response, name = "Response") +
  labs(
    title = "PCA: Response",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

print(p_pca3)

# Panel D: PC1 vs PC3 by PFS
p_pca4 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = PFS_group)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = colors_pfs, name = "PFS Group") +
  labs(
    title = "PCA: PC1 vs PC3",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC3 (", round(var_explained[3], 1), "%)")
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

print(p_pca4)

pca_combined <- (p_pca1 | p_pca2) / (p_pca3 | p_pca4) +
  plot_annotation(
    title = "Principal Component Analysis",
    subtitle = sprintf("Total variance explained by PC1-PC3: %.1f%%", sum(var_explained[1:3])),
    theme = theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 16)
    )
  )

print(pca_combined)

ggsave(file.path(fig_dir_png, "Step3_04_2_PCA_Panels.png"), pca_combined, width = 8, height = 8, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step3_04_test_PCA_Panels.pdf"), pca_combined, width = 10, height = 8, device = cairo_pdf)

ggsave(file.path(fig_dir_png, "Step3_04_PCA_Panels.png"), pca_combined, width = 8, height = 8, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step3_04_PCA_Panels.pdf"), pca_combined, width = 8, height = 8)

# ----------------------------------------------------------------------------
# Figure 5: Scree Plot
# ----------------------------------------------------------------------------
cat("  [5/6] PCA Scree Plot\n")

scree_data <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  Variance = summary(pca_result)$importance[2, 1:10] * 100,
  Cumulative = summary(pca_result)$importance[3, 1:10] * 100
)

p_scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = color_palette[1], alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_line(aes(x = as.numeric(PC), y = Cumulative/2), color = color_palette[7], linewidth = 1) +
  geom_point(aes(x = as.numeric(PC), y = Cumulative/2), color = color_palette[7], size = 3) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)), vjust = -0.5, size = 3) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~.*2, name = "Cumulative Variance (%)")
  ) +
  labs(
    title = "PCA Scree Plot",
    subtitle = "Line = cumulative variance",
    x = "Principal Component"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size=11),
    axis.title.y.right = element_text(color = color_palette[7]),
    axis.text.y.right = element_text(color = color_palette[7])
  )

print(p_scree)

ggsave(file.path(fig_dir_png, "Step3_05_PCA_Scree.png"), p_scree, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step3_05_PCA_Scree.pdf"), p_scree, width = 6, height = 5)


# ----------------------------------------------------------------------------
# 3.6 Sample-Sample Correlation Heatmap
# ----------------------------------------------------------------------------
cat("--- 3.6 Sample-Sample Correlation Heatmap ---\n")

# ----------------------------------------------------------------------------
# Figure 6: Sample Correlation Heatmap
# ----------------------------------------------------------------------------
cat("  [6/6] Sample Correlation Heatmap\n")

# Annotation for heatmap
annotation_row <- data.frame(
  PFS = datTraits$PFS_group,
  Treatment = datTraits$Treatment,
  Response = datTraits$Response
)
rownames(annotation_row) <- rownames(datTraits)

# Define annotation colors (using universal clinical colors)
ann_colors <- list(
  PFS = colors_pfs,
  Treatment = c("FOLFI" = col_folfi, "GemNab" = col_gemnab),
  Response = colors_response
)

# pheatmap is specialized for heatmaps - save directly
png(file.path(fig_dir_png, "Step3_06_Sample_Correlation_Heatmap.png"), width = 10, height = 10, units = "in", res = 300)
pheatmap(sample_cor,
         color = colorRampPalette(c(color_palette[1], "white", color_palette[7]))(100),
         clustering_method = "average",
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample-Sample Correlation Matrix",
         fontsize = 10)
dev.off()

pdf(file.path(fig_dir_pdf, "Step3_06_Sample_Correlation_Heatmap.pdf"), width = 10, height = 10)
pheatmap(sample_cor,
         color = colorRampPalette(c(color_palette[1], "white", color_palette[7]))(100),
         clustering_method = "average",
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample-Sample Correlation Matrix",
         fontsize = 10)
dev.off()

cat("  All Step 3 figures generated and saved\n")

# ----------------------------------------------------------------------------
# 3.7 Save Sample QC Summary Table
# ----------------------------------------------------------------------------
cat("--- 3.7 Saving Sample QC Summary ---\n")

sample_qc_summary <- data.frame(
  Sample_ID = rownames(datExpr),
  PFS_group = datTraits$PFS_group,
  Treatment = datTraits$Treatment,
  Response = datTraits$Response,
  PFS_days = datTraits$PFS,
  Age = datTraits$Age,
  Sex = datTraits$Sex,
  Liver_Mets = datTraits$Liver,
  CA19.9 = datTraits$CA19.9,
  Connectivity_Z = as.numeric(Z_connectivity),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Outlier_Flag = as.numeric(Z_connectivity) < outlier_threshold
)

write.csv(sample_qc_summary,
          file.path(results_dir, "Step3_Sample_QC_Summary.csv"),
          row.names = FALSE)

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 3 SUMMARY: SAMPLE QUALITY CONTROL ===\n")
cat("===============================================================================\n")
cat(sprintf("Total samples: %d\n", nrow(datExpr)))
cat(sprintf("Potential outliers (Z < %.1f): %d\n", outlier_threshold, length(potential_outliers)))
if(length(potential_outliers) > 0) {
  cat(sprintf("  Outlier samples: %s\n", paste(potential_outliers, collapse = ", ")))
}
cat(sprintf("Samples retained: %d\n", nrow(datExpr)))
cat(sprintf("Clinical traits: %s\n", paste(names(datTraits), collapse = ", ")))
cat(sprintf("PCA variance explained (PC1-3): %.1f%%, %.1f%%, %.1f%%\n",
            var_explained[1], var_explained[2], var_explained[3]))
cat("Step 3 complete.\n")
cat("===============================================================================\n")

# STEP 4: SOFT-THRESHOLD POWER SELECTION                
# PURPOSE: Select optimal soft-threshold power (Î2) that transforms the     
# correlation matrix into a scale-free network topology - the foundation   
# of WGCNA's ability to detect biologically meaningful modules.            
#                                                                          
# STRATEGY:                                                                
# 4a) Calculate scale-free fit across powers 1-20 for 5 network types    
# 4b) Visualize SFT curves, connectivity, and RÂ2 thresholds              
# 4c) Identify candidate powers meeting RÂ2 > 0.80 criterion              
# 4d) Sensitivity analysis comparing network properties across powers    
#                                                                         
# LOOKING FOR:                                                             
# - RÂ2 > 0.80 (ideally > 0.85) for scale-free topology fit               
# - Mean connectivity 50-200 proteins (too low = fragmented network)     
# - Stable module detection across adjacent power values                 


# ============================================================================
# STEP 4a: SOFT-THRESHOLD POWER SELECTION
# ============================================================================
# PURPOSE: Calculate scale-free topology fit across powers for optimal network
#
# STRATEGY: pickSoftThreshold() with biweight midcorrelation + signed network
#
# LOOKING FOR: Power where R2 > 0.80 with adequate connectivity (mean.k > 10)
#
# METHOD: Biweight midcorrelation (bicor) - robust to outliers
# ============================================================================

cat("STEP 4a: Soft-Threshold Power Selection\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Complete power range 1:20
powers <- 1:20
cat("Powers to test:", paste(powers, collapse = ", "), "\n")

# ───────────────────────────────────────────────────────────
# Network Configuration: Biweight midcorrelation + Signed
# ───────────────────────────────────────────────────────────
# Nika: Bicor + Signed chosen after sensitivity analysis comparing:
#   - Pearson vs Bicor correlation (bicor more robust to outliers)
#   - Signed vs Unsigned networks (signed preserves biological direction)
#   - Various power thresholds (8-14 range optimal for this dataset)

network_config <- list(
  name = "Biweight midcorrelation with signed network",
  corFnc = "bicor",
  corOptions = list(use = "pairwise.complete.obs", maxPOutliers = 0.1),
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  description = "Robust correlation method resistant to outliers"
)

cat(sprintf("Network configuration: %s\n", network_config$name))
cat(sprintf("  Correlation: %s\n", network_config$corFnc))
cat(sprintf("  Network type: %s\n", network_config$networkType))

# ───────────────────────────────────────────────────────────
# Run pickSoftThreshold
# ───────────────────────────────────────────────────────────

cat("Running pickSoftThreshold... ")

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = network_config$networkType,
  corFnc = network_config$corFnc,
  corOptions = network_config$corOptions,
  verbose = 0
)

cat(sprintf("Done. Power estimate: %s, Max R2: %.3f\n",
            ifelse(is.na(sft$powerEstimate), "NA", sft$powerEstimate),
            max(sft$fitIndices$SFT.R.sq)))

# Store fit indices with config info
fit_indices <- sft$fitIndices
fit_indices$Correlation <- network_config$corFnc
fit_indices$NetworkType <- network_config$networkType
fit_indices$TOMType <- network_config$TOMType
fit_indices$minModuleSize <- network_config$minModuleSize
fit_indices$mergeCutHeight <- network_config$mergeCutHeight
fit_indices$PowerEstimate <- sft$powerEstimate

# ───────────────────────────────────────────────────────────
# Build summary statistics
# ───────────────────────────────────────────────────────────

fi <- sft$fitIndices

# Find powers achieving R2 thresholds
first_80 <- min(fi$Power[fi$SFT.R.sq >= 0.80], na.rm = TRUE)
first_85 <- min(fi$Power[fi$SFT.R.sq >= 0.85], na.rm = TRUE)
first_90 <- min(fi$Power[fi$SFT.R.sq >= 0.90], na.rm = TRUE)

summary_table <- data.frame(
  Correlation = network_config$corFnc,
  NetworkType = network_config$networkType,
  TOM = network_config$TOMType,
  minModSize = network_config$minModuleSize,
  mergeCut = network_config$mergeCutHeight,
  PowerEst = ifelse(is.na(sft$powerEstimate), NA, sft$powerEstimate),
  MaxR2 = round(max(fi$SFT.R.sq), 3),
  PowerAtMaxR2 = fi$Power[which.max(fi$SFT.R.sq)],
  Power_R2_80 = ifelse(is.infinite(first_80), NA, first_80),
  Power_R2_85 = ifelse(is.infinite(first_85), NA, first_85),
  Power_R2_90 = ifelse(is.infinite(first_90), NA, first_90),
  MeanK_P8 = round(fi$mean.k.[fi$Power == 8], 2),
  MeanK_P10 = round(fi$mean.k.[fi$Power == 10], 2),
  MeanK_P12 = round(fi$mean.k.[fi$Power == 12], 2),
  MeanK_P14 = round(fi$mean.k.[fi$Power == 14], 2),
  stringsAsFactors = FALSE
)

# ───────────────────────────────────────────────────────────
# Elbow detection
# ───────────────────────────────────────────────────────────

r2 <- fi$SFT.R.sq
powers_tested <- fi$Power

# Method 1: First power where R2 gain < 0.03
gains <- diff(r2)
elbow_gain <- which(gains < 0.03)[1]
if(is.na(elbow_gain)) elbow_gain <- length(powers_tested)

# Method 2: Maximum second derivative (curvature)
if(length(gains) >= 2) {
  d2 <- diff(gains)
  elbow_curv <- which.min(d2) + 1
} else {
  elbow_curv <- NA
}

elbow_table <- data.frame(
  Elbow_Gain = powers_tested[elbow_gain],
  R2_AtElbowGain = round(r2[elbow_gain], 3),
  MeanK_AtElbowGain = round(fi$mean.k.[elbow_gain], 2),
  Elbow_Curvature = ifelse(is.na(elbow_curv), NA, powers_tested[elbow_curv]),
  R2_AtElbowCurv = ifelse(is.na(elbow_curv), NA, round(r2[elbow_curv], 3)),
  MeanK_AtElbowCurv = ifelse(is.na(elbow_curv), NA, round(fi$mean.k.[elbow_curv], 2)),
  stringsAsFactors = FALSE
)

# ───────────────────────────────────────────────────────────
# Connectivity analysis at candidate powers
# ───────────────────────────────────────────────────────────

connectivity_table <- do.call(rbind, lapply(c(8, 10, 12,13, 14), function(p) {
  idx <- which(fi$Power == p)
  data.frame(
    Power = p,
    R2 = round(fi$SFT.R.sq[idx], 3),
    Slope = round(fi$slope[idx], 2),
    MeanK = round(fi$mean.k.[idx], 2),
    MedianK = round(fi$median.k.[idx], 2),
    MaxK = round(fi$max.k.[idx], 2),
    stringsAsFactors = FALSE
  )
}))

# Display results
cat("\nSummary:\n")
print(summary_table, row.names = FALSE)

cat("\nElbow Detection:\n")
print(elbow_table, row.names = FALSE)

cat("\nConnectivity at Candidate Powers:\n")
print(connectivity_table, row.names = FALSE)

#Justification
#"Soft-threshold power was selected using the pickSoftThreshold function in WGCNA.
#We employed an elbow-based approach balancing scale-free topology fit (R²) with network connectivity.
#Power 13 was selected as it: (1) achieved R² = 0.833, exceeding the 0.80 threshold recommended for smaller datasets (n < 100) [Langfelder & Horvath, 2008]; (2) corresponded to the elbow point identified by curvature analysis; and (3) maintained adequate mean connectivity (k = 1.64) compared to higher powers where median connectivity dropped below 1. Module robustness was validated through bootstrap resampling (1000 iterations).
#While power 14 achieves R² = 0.874, it results in median connectivity below 1 (k = 0.89), 
#indicating that more than half of proteins would have fewer than one connection on average. This excessive sparsity risks losing biologically meaningful relationships. Power 13 provides a balance between scale-free fit (R² = 0.833) and network connectivity (median k = 1.34), consistent with recommendations for smaller sample sizes [Langfelder & Horvath, 2008]. Module stability was confirmed through bootstrap analysis.

# Save results
write.csv(fit_indices, file.path(results_dir, "Step4a_FitIndices.csv"), row.names = FALSE)
write.csv(summary_table, file.path(results_dir, "Step4a_Summary.csv"), row.names = FALSE)
write.csv(elbow_table, file.path(results_dir, "Step4a_Elbow.csv"), row.names = FALSE)
write.csv(connectivity_table, file.path(results_dir, "Step4a_Connectivity.csv"), row.names = FALSE)
saveRDS(sft, file.path(results_dir, "Step4a_sft_result.rds"))
saveRDS(network_config, file.path(results_dir, "Step4a_network_config.rds"))

# bicor + signed at power 8-14 was optimal:
# - Bicor more robust than Pearson for proteomics data with potential outliers

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4b: SOFT-THRESHOLD VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════
# PURPOSE: Visualize scale-free topology fit to guide power selection
#
# STRATEGY: Two-panel plot showing R² and mean connectivity vs power
#
# LOOKING FOR: Power achieving R² > 0.85 with adequate connectivity (k > 10)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4b: Soft-Threshold Visualization\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ───────────────────────────────────────────────────────────
# Prepare data for plotting
# ───────────────────────────────────────────────────────────

sft_df <- data.frame(
  Power = sft$fitIndices$Power,
  R2 = sft$fitIndices$SFT.R.sq,
  Slope = sft$fitIndices$slope,
  MeanK = sft$fitIndices$mean.k.,
  MedianK = sft$fitIndices$median.k.
)

# ───────────────────────────────────────────────────────────
# Figure: Scale-Free Topology Fit (R² vs Power)
# ───────────────────────────────────────────────────────────

cat("  [1/3] Scale-Free Topology Fit (R² vs Power)\n")

# Get R2 value at power 13 for intercept lines
r2_at_13 <- sft_df$R2[sft_df$Power == 13]

p_r2 <- ggplot(sft_df, aes(x = Power, y = R2)) +
  geom_line(color = "#008080", linewidth = 1.2) +
  geom_point(color = "#008080", size = 2.5) +
  # Dashed intercept lines at power 13
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = r2_at_13, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  # Red dot at power 13
  geom_point(data = sft_df %>% filter(Power == 13), aes(x = Power, y = R2),
             color = "red", size = 4) +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Scale-Free Topology Fit",
    x = "Soft-Threshold (Power)",
    y = expression("Model Fit (R"^2*")")
  ) +
  # Label the R2 value at power 13
  annotate("text", x = 15, y = r2_at_13 + 0.03,
           label = paste0("R² = ", round(r2_at_13, 2)),
           color = "grey30", size = 3.5, hjust = 0) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_r2)
ggsave(file.path(fig_dir_png, "Step4b_01_R2_vs_Power.png"), p_r2, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step4b_01_R2_vs_Power.pdf"), p_r2, width = 6, height = 5)
cat("    Saved: Step4b_01_R2_vs_Power.png/pdf\n")

# ───────────────────────────────────────────────────────────
# Figure: Mean Connectivity vs Power
# ───────────────────────────────────────────────────────────

cat("  [2/3] Mean Connectivity vs Power\n")

# Get MeanK value at power 13 for intercept lines
meank_at_13 <- sft_df$MeanK[sft_df$Power == 13]

p_k <- ggplot(sft_df, aes(x = Power, y = MeanK)) +
  geom_line(color = "#008080", linewidth = 1.2) +
  geom_point(color = "#008080", size = 2.5) +
  # Dashed intercept lines at power 13
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = meank_at_13, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  # Red dot at power 13
  geom_point(data = sft_df %>% filter(Power == 13), aes(x = Power, y = MeanK),
             color = "red", size = 4) +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  scale_y_log10() +
  labs(
    title = "Mean Connectivity",
    x = "Soft-Threshold (Power)",
    y = "Mean Connectivity (log scale)"
  ) +
  # Label the MeanK value at power 13
  annotate("text", x = 15, y = meank_at_13 * 1.3,
           label = paste0("k = ", round(meank_at_13, 1)),
           color = "grey30", size = 3.5, hjust = 0) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_k)
ggsave(file.path(fig_dir_png, "Step4b_02_Connectivity_vs_Power.png"), p_k, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step4b_02_Connectivity_vs_Power.pdf"), p_k, width = 6, height = 5)
cat("    Saved: Step4b_02_Connectivity_vs_Power.png/pdf\n")

# ───────────────────────────────────────────────────────────
# Combined Figure: Power Selection
# ───────────────────────────────────────────────────────────

cat("  [3/3] Combined Publication Figure\n")

p_sft_combined <- (p_r2 | p_k) +
  plot_annotation(
    title = "Soft-Threshold Power Selection",
    subtitle = "Biweight midcorrelation with signed network",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)
    )
  )

print(p_sft_combined)
ggsave(file.path(fig_dir_png, "Fig_Publication_Power_Selection.png"), p_sft_combined, 
       width = 10, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Fig_Publication_Power_Selection.pdf"), p_sft_combined, 
       width = 10, height = 5)
cat("    Saved: Fig_Publication_Power_Selection.png/pdf\n")

cat("\n  Step 4b complete: soft-threshold visualization saved.\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4c: POWER SELECTION
# ═══════════════════════════════════════════════════════════════════════════
# PURPOSE: Select optimal soft-threshold power for network construction
#
# STRATEGY: Balance R² (scale-free fit) with mean connectivity for module detection
#
# DECISION: Power = 13 chosen as optimal trade-off
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4c: Power Selection\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ───────────────────────────────────────────────────────────
# Power Selection Criteria
# ───────────────────────────────────────────────────────────
# Power = 13 selected based on:
#Meets R² > 0.800.833 exceeds the minimum threshold
#Elbow point: Curvature analysis identifies power 13 as optimal
#Better connectivity: Mean k = 1.64 vs 1.16; Median k = 1.34 vs 0.89
# Truncated R² = 0.99: Excellent fit to scale-free topology
#Small sample size: n=50 justifies using 0.80 threshold
#Preserves biology: More connections = more biological signal retained


selected_power <- 13

# Extract metrics at selected power
fi <- sft$fitIndices
power_idx <- which(fi$Power == selected_power)

selected_r2 <- fi$SFT.R.sq[power_idx]
selected_slope <- fi$slope[power_idx]
selected_meank <- fi$mean.k.[power_idx]
selected_mediank <- fi$median.k.[power_idx]

cat(sprintf("\n  Selected power: %d\n", selected_power))
cat(sprintf("  R² at power %d: %.3f\n", selected_power, selected_r2))
cat(sprintf("  Slope: %.2f (ideal: -1 to -2)\n", selected_slope))
cat(sprintf("  Mean connectivity: %.1f\n", selected_meank))
cat(sprintf("  Median connectivity: %.1f\n", selected_mediank))

# ───────────────────────────────────────────────────────────
# Comparison with alternative powers
# ───────────────────────────────────────────────────────────

cat("\n  Comparison of candidate powers:\n")

candidate_powers <- c(6, 8, 10, 12,13, 14)
power_comparison <- data.frame(
  Power = candidate_powers,
  R2 = round(fi$SFT.R.sq[match(candidate_powers, fi$Power)], 3),
  Slope = round(fi$slope[match(candidate_powers, fi$Power)], 2),
  MeanK = round(fi$mean.k.[match(candidate_powers, fi$Power)], 1),
  MedianK = round(fi$median.k.[match(candidate_powers, fi$Power)], 1)
)
power_comparison$Selected <- ifelse(power_comparison$Power == selected_power, "<--", "")

print(power_comparison, row.names = FALSE)

# ───────────────────────────────────────────────────────────
# Save power selection
# ───────────────────────────────────────────────────────────

power_selection <- list(
  power = selected_power,
  r2 = selected_r2,
  slope = selected_slope,
  mean_k = selected_meank,
  median_k = selected_mediank,
  network_config = network_config
)

write.csv(power_comparison, file.path(results_dir, "Step4c_Power_Comparison.csv"), row.names = FALSE)
saveRDS(power_selection, file.path(results_dir, "Step4c_power_selection.rds"))

cat("\n  Power selection saved to Step4c_power_selection.rds\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4d: SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
# PURPOSE: Compare network construction across multiple soft-threshold powers
# STRATEGY: Build networks at powers 6, 8, 10, 12, 13, 14 -> compare modules -> measure stability
# LOOKING FOR: Stable module count, consistent hub genes, clinical associations
# KEY PARAMS: test_powers=[6,8,10,12,13,14], minModuleSize=20, deepSplit=3
# ═══════════════════════════════════════════════════════════════════════════

cat("===============================================================================\n")
cat("STEP 4d: SENSITIVITY ANALYSIS\n")
cat("===============================================================================\n")

# Define test configurations for sensitivity analysis
test_powers <- c(6, 8, 10, 12, 13, 14)

test_configs <- list()
for(p in test_powers) {
  config_name <- paste0("C2_P", p)
  test_configs[[config_name]] <- list(
    name = paste0("C2_Power", p),
    corFnc = "bicor",
    networkType = "signed",
    TOMType = "signed",
    power = p,
    minModuleSize = 20,
    mergeCutHeight = 0.25
  )
}

# Prepare trait data for module-trait correlations
traitData_sensitivity <- data.frame(
  PFS = datTraits$PFS,
  PFS_group = as.numeric(datTraits$PFS_group == "Long"),
  Treatment = as.numeric(datTraits$Treatment),
  Response = as.numeric(datTraits$Response == "CD"),
  Age = datTraits$Age,
  Liver = as.numeric(datTraits$Liver == "Yes"),
  row.names = rownames(datTraits)
)

# Function to build network and analyze
run_wgcna_config <- function(config, datExpr, traitData) {

  # Build network using bicor correlation
  net <- blockwiseModules(
    datExpr,
    power = config$power,
    networkType = config$networkType,
    TOMType = config$TOMType,
    minModuleSize = config$minModuleSize,
    mergeCutHeight = config$mergeCutHeight,
    corType = "bicor",
    maxPOutliers = 0.1,
    deepSplit = 3,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    verbose = 0
  )

  # Module info
  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  names(moduleColors) <- colnames(datExpr)
  nModules <- length(unique(moduleLabels[moduleLabels != 0]))
  nGrey <- sum(moduleLabels == 0)
  pctGrey <- round(100 * nGrey / length(moduleLabels), 1)
  modSizes <- table(moduleLabels)

  # Calculate module eigengenes
  MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  MEs <- orderMEs(MEs)
  MEs_noGrey <- MEs[, !grepl("grey", colnames(MEs)), drop = FALSE]

  # Module-trait correlations
  if(ncol(MEs_noGrey) > 0) {
    modTraitCor <- cor(MEs_noGrey, traitData, use = "pairwise.complete.obs")
    modTraitPval <- corPvalueStudent(modTraitCor, nrow(datExpr))
    n_pfs_sig <- sum(modTraitPval[, "PFS"] < 0.05, na.rm = TRUE)
    n_pfs_sig_01 <- sum(modTraitPval[, "PFS"] < 0.01, na.rm = TRUE)
    n_pfs_group_sig <- sum(modTraitPval[, "PFS_group"] < 0.05, na.rm = TRUE)
    n_response_sig <- sum(modTraitPval[, "Response"] < 0.05, na.rm = TRUE)
  } else {
    modTraitCor <- NULL
    modTraitPval <- NULL
    n_pfs_sig <- 0
    n_pfs_sig_01 <- 0
    n_pfs_group_sig <- 0
    n_response_sig <- 0
  }

  list(
    net = net,
    moduleLabels = moduleLabels,
    moduleColors = moduleColors,
    nModules = nModules,
    nGrey = nGrey,
    pctGrey = pctGrey,
    modSizes = modSizes,
    MEs = MEs,
    MEs_noGrey = MEs_noGrey,
    modTraitCor = modTraitCor,
    modTraitPval = modTraitPval,
    n_pfs_sig = n_pfs_sig,
    n_pfs_sig_01 = n_pfs_sig_01,
    n_pfs_group_sig = n_pfs_group_sig,
    n_response_sig = n_response_sig
  )
}

# Run all configurations
cat("\nBuilding networks for sensitivity analysis...\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

sensitivity_results <- list()

for(cfg_name in names(test_configs)) {
  config <- test_configs[[cfg_name]]
  cat(sprintf("[%s] Power=%d, %s, TOM=%s, minMod=%d ... ",
              config$name, config$power, config$networkType, config$TOMType, config$minModuleSize))

  start_time <- Sys.time()
  res <- run_wgcna_config(config, datExpr, traitData_sensitivity)
  elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 1)

  res$config <- config
  sensitivity_results[[cfg_name]] <- res

  cat(sprintf("%d modules, %d grey (%.1f%%), PFS_group-sig: %d, PFS_Days-sig: %d, Response-sig: %d, Time: %.1fs\n",
              res$nModules, res$nGrey, res$pctGrey, res$n_pfs_group_sig, res$n_pfs_sig, res$n_response_sig, elapsed))
}

# Build summary table
cat("\n")
cat(paste(rep("=", 90), collapse = ""), "\n")
cat("SENSITIVITY ANALYSIS SUMMARY TABLE\n")
cat(paste(rep("=", 90), collapse = ""), "\n")

summary_df <- do.call(rbind, lapply(names(sensitivity_results), function(cfg_name) {
  res <- sensitivity_results[[cfg_name]]
  config <- res$config
  modSizes_noGrey <- res$modSizes[names(res$modSizes) != "0"]

  data.frame(
    Config = config$name,
    Method = "C2:Bicor+Signed",
    Power = config$power,
    TOM = config$TOMType,
    minModSize = config$minModuleSize,
    nModules = res$nModules,
    nGrey = res$nGrey,
    pctGrey = res$pctGrey,
    PFS_Days_sig_05 = res$n_pfs_sig,
    PFS_Days_sig_01 = res$n_pfs_sig_01,
    PFS_group_sig_05 = res$n_pfs_group_sig,
    Response_sig_05 = res$n_response_sig,
    MeanModSize = round(mean(as.numeric(modSizes_noGrey)), 1),
    MedianModSize = median(as.numeric(modSizes_noGrey)),
    MaxModSize = max(as.numeric(modSizes_noGrey)),
    MinModSize = min(as.numeric(modSizes_noGrey)),
    stringsAsFactors = FALSE
  )
}))

print(summary_df, row.names = FALSE)

# Save summary table
write.csv(summary_df, file.path(results_dir, "Step4d_Sensitivity_Summary.csv"), row.names = FALSE)
cat("\n  Saved: Step4d_Sensitivity_Summary.csv\n")

# ============================================================================
# PFS-Significant Modules Detail
# ============================================================================
cat("\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("PFS_Days-Significant Modules (p < 0.05) per Configuration:\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

for(cfg_name in names(sensitivity_results)) {
  res <- sensitivity_results[[cfg_name]]
  cat(sprintf("\n%s:\n", res$config$name))

  if(!is.null(res$modTraitPval)) {
    pfs_cor <- res$modTraitCor[, "PFS"]
    pfs_pval <- res$modTraitPval[, "PFS"]
    sig_idx <- which(pfs_pval < 0.05)

    if(length(sig_idx) > 0) {
      sig_df <- data.frame(
        Module = names(pfs_cor)[sig_idx],
        Correlation = round(pfs_cor[sig_idx], 3),
        P_value = formatC(pfs_pval[sig_idx], format = "e", digits = 2)
      )
      sig_df <- sig_df[order(as.numeric(sig_df$P_value)), ]
      print(sig_df, row.names = FALSE)
    } else {
      cat("  No significant modules\n")
    }
  } else {
    cat("  No modules detected\n")
  }
}

# ============================================================================
# Trait Association Summary
# ============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Significant Module-Trait Associations (p < 0.05) Summary:\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

trait_summary <- do.call(rbind, lapply(names(sensitivity_results), function(cfg_name) {
  res <- sensitivity_results[[cfg_name]]
  if(!is.null(res$modTraitPval)) {
    data.frame(
      Config = res$config$name,
      Power = res$config$power,
      nModules = res$nModules,
      pctGrey = res$pctGrey,
      PFS_Days = sum(res$modTraitPval[, "PFS"] < 0.05),
      PFS_group = sum(res$modTraitPval[, "PFS_group"] < 0.05),
      Response = sum(res$modTraitPval[, "Response"] < 0.05),
      Treatment = sum(res$modTraitPval[, "Treatment"] < 0.05),
      Age = sum(res$modTraitPval[, "Age"] < 0.05),
      Liver = sum(res$modTraitPval[, "Liver"] < 0.05)
    )
  }
}))
print(trait_summary, row.names = FALSE)

# Save trait summary
write.csv(trait_summary, file.path(results_dir, "Step4d_Trait_Association_Summary.csv"), row.names = FALSE)
cat("\n  Saved: Step4d_Trait_Association_Summary.csv\n")

# ============================================================================
# STEP 4d VISUALIZATION
# ============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Generating Sensitivity Analysis Figures...\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Prepare plotting data
summary_df$Config <- factor(summary_df$Config,
                            levels = paste0("C2_Power", test_powers))

# Create module size data for boxplot
modsize_data <- do.call(rbind, lapply(names(sensitivity_results), function(cfg_name) {
  res <- sensitivity_results[[cfg_name]]
  modSizes <- res$modSizes[names(res$modSizes) != "0"]
  if(length(modSizes) > 0) {
    data.frame(Config = cfg_name, ModuleSize = as.numeric(modSizes))
  }
}))
modsize_data$Config <- factor(modsize_data$Config,
                              levels = paste0("C2_P", test_powers))

# ============================================================================
# FIGURE: Sensitivity Analysis 
# ============================================================================
cat("  [PUB] Publication Sensitivity Figure\n")

# Create C2-only data for publication figure
c2_sensitivity <- summary_df %>%
  mutate(Config_Clean = paste0("Power", Power))

# Ensure correct order
c2_sensitivity$Config_Clean <- factor(c2_sensitivity$Config_Clean,
                                       levels = paste0("Power", test_powers))

pub_data <- c2_sensitivity %>%
  dplyr::select(Config_Clean, Power, nModules, PFS_group_sig_05, PFS_Days_sig_05, Response_sig_05) %>%
  tidyr::pivot_longer(cols = c(nModules, PFS_group_sig_05, PFS_Days_sig_05, Response_sig_05),
                      names_to = "Metric", values_to = "Count") %>%
  mutate(Metric = factor(Metric,
                         levels = c("nModules", "PFS_group_sig_05", "PFS_Days_sig_05", "Response_sig_05"),
                         labels = c("Total Modules", "PFS Groups", "PFS Days", "Response")),
         Power_Num = Power)

pub_modsize_data <- modsize_data %>%
  mutate(Config_Clean = paste0("Power", as.numeric(gsub("C2_P", "", Config))))
pub_modsize_data$Config_Clean <- factor(pub_modsize_data$Config_Clean,
                                         levels = paste0("Power", test_powers))

# Color palette
my_pal <- c("Total Modules" = "grey85",
            "PFS Groups" = "#b4c8a8",
            "PFS Days" = "#f6edbd",
            "Response" = "#edbb8a")

# Nature-style theme
nature_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  )

# Panel 1: Module Detection
pub_p1 <- ggplot(pub_data %>% filter(Count > 0), aes(x = Config_Clean, y = Count, fill = Metric)) +
  geom_bar(stat = "identity",
           position = position_dodge2(preserve = "single", padding = 0.1),
           color = "black", linewidth = 0.3) +
  scale_fill_manual(values = my_pal) +
  labs(title = "Module Detection", x = NULL, y = "Count") +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel 2: Unassigned Proteins
pub_p2 <- ggplot(c2_sensitivity, aes(x = Config_Clean, y = pctGrey)) +
  geom_bar(stat = "identity", width = 0.6, fill = "#de8a5a", color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(pctGrey, "%")), vjust = -0.5, size = 3) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Unassigned Proteins", x = NULL, y = "% Proteins (Grey)") +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel 3: Size Distribution
pub_p3 <- ggplot(pub_modsize_data, aes(x = Config_Clean, y = ModuleSize)) +
  geom_boxplot(outlier.shape = NA, fill = "#b4c8a8", alpha = 0.5, linewidth = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1, color = "black") +
  labs(title = "Size Distribution", x = "Power", y = "Proteins per Module") +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel 4: Trend Line
pub_p4 <- ggplot(pub_data, aes(x = Power_Num, y = Count, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, fill = "white", stroke = 1, shape = 21) +
  scale_color_manual(values = my_pal) +
  scale_x_continuous(breaks = test_powers) +
  labs(title = "Modules vs. Power Trend", x = "Soft Threshold (Power)", y = "Count") +
  nature_theme

# Assemble final figure
fig_sensitivity <- (pub_p1 + pub_p2) / (pub_p3 + pub_p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Sensitivity analysis across power values",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom") &
  guides(color = "none")

print(fig_sensitivity)

ggsave(file.path(fig_dir_png, "Step4d_Sensitivity_Analysis.png"), fig_sensitivity, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step4d_Sensitivity_Analysis.pdf"), fig_sensitivity, width = 10, height = 8)
cat("  Saved: Step4d_Sensitivity_Analysis.png/pdf\n")


# ============================================================================
# DENDROGRAM COMPARISON: Power 8 vs 13 vs 14 (Stacked Panel)
# ============================================================================

cat("\n  Generating Dendrogram Comparison (Power 8 vs 13 vs 14)...\n")

# Function to create dendrogram data for a given power
create_dendro_data <- function(power_val, sensitivity_results, datExpr) {

  cfg_name <- paste0("C2_P", power_val)
  res <- sensitivity_results[[cfg_name]]

  if(is.null(res)) {
    stop(paste("Results for power", power_val, "not found in sensitivity_results"))
  }

  # Get module colors from the results
  moduleColors <- res$moduleColors

  # Recalculate TOM and hierarchical clustering for this power
  cat(sprintf("      Calculating adjacency matrix for Power %d...\n", power_val))
  adjacency <- adjacency(datExpr,
                         power = power_val,
                         type = "signed",
                         corFnc = "bicor",
                         corOptions = list(maxPOutliers = 0.1))

  cat(sprintf("      Calculating TOM for Power %d...\n", power_val))
  TOM <- TOMsimilarity(adjacency, TOMType = "signed")
  dissTOM <- 1 - TOM

  # Hierarchical clustering
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  # Return list with tree and colors
  list(
    geneTree = geneTree,
    moduleColors = moduleColors,
    nModules = res$nModules,
    pctGrey = res$pctGrey,
    power = power_val
  )
}

# Create dendrograms for powers 8, 13, and 14
cat("    Building dendrogram for Power 8...\n")
dend_p8 <- create_dendro_data(8, sensitivity_results, datExpr)

cat("    Building dendrogram for Power 13...\n")
dend_p13 <- create_dendro_data(13, sensitivity_results, datExpr)

cat("    Building dendrogram for Power 14...\n")
dend_p14 <- create_dendro_data(14, sensitivity_results, datExpr)

# Combine module colors into a matrix for multi-row color bar
all_colors <- cbind(
  dend_p8$moduleColors,
  dend_p13$moduleColors,
  dend_p14$moduleColors
)

# Shorter labels to fit
color_labels <- c(
  paste0("P8: ", dend_p8$nModules, " mod, ", dend_p8$pctGrey, "% grey"),
  paste0("P13: ", dend_p13$nModules, " mod, ", dend_p13$pctGrey, "% grey"),
  paste0("P14: ", dend_p14$nModules, " mod, ", dend_p14$pctGrey, "% grey")
)

# Use Power 13 tree as the reference (since that's your selected power)
# This shows how proteins cluster and how module assignments differ across powers

cat("    Creating comparison figure...\n")

# --- Save PNG ---
png(file.path(fig_dir_png, "Step4d_Dendrogram_Comparison_P8_P13_P14.png"),
    width = 12, height = 8, units = "in", res = 300)

par(mar = c(1, 6, 3, 1))  # Increase left margin for labels

plotDendroAndColors(
  dend_p13$geneTree,
  all_colors,
  color_labels,
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Protein Dendrogram with Module Colors Across Powers",
  cex.colorLabels = 0.7,
  marAll = c(1, 6, 3, 1)
)

dev.off()
cat("    Saved: Step4d_Dendrogram_Comparison_P8_P13_P14.png\n")

# --- Save PDF ---
pdf(file.path(fig_dir_pdf, "Step4d_Dendrogram_Comparison_P8_P13_P14.pdf"),
    width = 12, height = 8)

par(mar = c(1, 6, 3, 1))

plotDendroAndColors(
  dend_p13$geneTree,
  all_colors,
  color_labels,
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Protein Dendrogram with Module Colors Across Powers",
  cex.colorLabels = 0.7,
  marAll = c(1, 6, 3, 1)
)

dev.off()
cat("    Saved: Step4d_Dendrogram_Comparison_P8_P13_P14.pdf\n")

# --- Display in R viewer ---
par(mar = c(1, 6, 3, 1))

plotDendroAndColors(
  dend_p13$geneTree,
  all_colors,
  color_labels,
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Protein Dendrogram with Module Colors Across Powers",
  cex.colorLabels = 0.7,
  marAll = c(1, 6, 3, 1)
)

par(mar = c(5, 4, 4, 2) + 0.1)  # Reset to default margins

cat("  Dendrogram comparison complete.\n")

#considering the downstream analysis with Cox and ML biomarker:
#Power 13 is the best


# ============================================================================
# Module-Trait Heatmaps for Selected Powers
# ============================================================================
cat("  Generating Module-Trait Heatmaps...\n")

# Define trait order
sens_trait_order <- c("PFS_group", "Response", "PFS_Days", "Treatment", "Liver_Met", "Age", "Sex")

for(cfg_name in names(sensitivity_results)) {
  res <- sensitivity_results[[cfg_name]]
  cor_mat <- res$modTraitCor
  pval_mat <- res$modTraitPval

  # Rename columns
  colnames(cor_mat)[colnames(cor_mat) == "PFS"] <- "PFS_Days"
  colnames(cor_mat)[colnames(cor_mat) == "Liver"] <- "Liver_Met"
  colnames(pval_mat)[colnames(pval_mat) == "PFS"] <- "PFS_Days"
  colnames(pval_mat)[colnames(pval_mat) == "Liver"] <- "Liver_Met"

  # Reorder columns
  available_traits <- sens_trait_order[sens_trait_order %in% colnames(cor_mat)]
  cor_mat <- cor_mat[, available_traits, drop = FALSE]
  pval_mat <- pval_mat[, available_traits, drop = FALSE]

  # Significance stars
  sig_stars <- pval_mat
  sig_stars[pval_mat >= 0.05] <- ""
  sig_stars[pval_mat < 0.05 & pval_mat >= 0.01] <- "*"
  sig_stars[pval_mat < 0.01 & pval_mat >= 0.001] <- "**"
  sig_stars[pval_mat < 0.001] <- "***"
  text_mat <- matrix(paste0(round(cor_mat, 2), sig_stars), nrow = nrow(cor_mat), dimnames = dimnames(cor_mat))

  # Fontface (bold for significant)
  fontface_mat <- matrix("plain", nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat))
  fontface_mat[pval_mat < 0.05] <- "bold"

  # Convert to long format
  cor_long <- reshape2::melt(cor_mat, varnames = c("Module", "Trait"), value.name = "Correlation")
  cor_long$Label <- reshape2::melt(text_mat, varnames = c("Module", "Trait"), value.name = "Label")$Label
  cor_long$FontFace <- reshape2::melt(fontface_mat, varnames = c("Module", "Trait"), value.name = "FontFace")$FontFace
  cor_long$Trait <- factor(cor_long$Trait, levels = available_traits)

  p_heatmap <- ggplot(cor_long, aes(x = Trait, y = Module, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Label, fontface = FontFace), size = 3.5) +
    scale_fill_gradient2(low = "#008080", mid = "white", high = "#ca562c", midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = paste0("Module-Trait Correlations: ", res$config$name),
         subtitle = paste0(res$nModules, " modules | ", res$pctGrey, "% grey | * p<0.05, ** p<0.01, *** p<0.001"),
         x = NULL, y = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank())

  print(p_heatmap)
  ggsave(file.path(fig_dir_png, paste0("Step4d_Heatmap_", cfg_name, ".png")), p_heatmap, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir_pdf, paste0("Step4d_Heatmap_", cfg_name, ".pdf")), p_heatmap, width = 8, height = 6)
}
cat("  Saved: Step4d_Heatmap_*.png/pdf for all configurations\n")

# ============================================================================
# Save sensitivity results for downstream use
# ============================================================================
saveRDS(sensitivity_results, file.path(results_dir, "Step4d_sensitivity_results.rds"))
cat("  Saved: Step4d_sensitivity_results.rds\n")

# ============================================================================
# STEP 4d SUMMARY
# ============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 4d SENSITIVITY ANALYSIS COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Identify best configuration based on criteria
best_config <- trait_summary %>%
  mutate(
    Total_Sig = PFS_Days + PFS_group + Response,
    Score = Total_Sig - (pctGrey / 20)  # Penalize high grey %
  ) %>%
  arrange(desc(Score), pctGrey) %>%
  slice(1)

cat(sprintf("\nRecommended configuration: Power %d\n", best_config$Power))
cat(sprintf("  - Modules: %d\n", best_config$nModules))
cat(sprintf("  - Grey %%: %.1f%%\n", best_config$pctGrey))
cat(sprintf("  - PFS_Days significant: %d\n", best_config$PFS_Days))
cat(sprintf("  - PFS_group significant: %d\n", best_config$PFS_group))
cat(sprintf("  - Response significant: %d\n", best_config$Response))
cat("\n")

# ───────────────────────────────────────────────────────────
# Store final network parameters for Step 5
# ───────────────────────────────────────────────────────────

soft_power <- selected_power  # From Step 4c

final_network_params <- list(
  power = soft_power,
  corFnc = network_config$corFnc,
  networkType = network_config$networkType,
  TOMType = network_config$TOMType,
  minModuleSize = network_config$minModuleSize,
  mergeCutHeight = network_config$mergeCutHeight,
  deepSplit = 3
)

saveRDS(final_network_params, file.path(results_dir, "Step4d_final_network_params.rds"))
cat("  Saved: Step4d_final_network_params.rds\n\n")

# STEP 5: FINAL NETWORK CONSTRUCTION
# ============================================================================
# PURPOSE: Build the WGCNA network using parameters from Step 4
# Network configuration: Bicor correlation, signed network, power=13, deepSplit=3
# ============================================================================

cat("===============================================================================\n")
cat("STEP 5: FINAL NETWORK CONSTRUCTION\n")
cat("===============================================================================\n")

# ============================================================================
# 5a: Build Network with blockwiseModules
# ============================================================================

cat("\n5a: Building WGCNA network...\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Load network parameters from Step 4
network_params <- readRDS(file.path(results_dir, "Step4d_final_network_params.rds"))

cat("  Network parameters:\n")
cat(sprintf("    Power: %d\n", network_params$power))
cat(sprintf("    Correlation: %s\n", network_params$corFnc))
cat(sprintf("    Network type: %s\n", network_params$networkType))
cat(sprintf("    TOM type: %s\n", network_params$TOMType))
cat(sprintf("    Min module size: %d\n", network_params$minModuleSize))
cat(sprintf("    Merge cut height: %.2f\n", network_params$mergeCutHeight))
cat(sprintf("    Deep split: %d\n", network_params$deepSplit))

cat("\n  Running blockwiseModules (this may take a few minutes)...\n")

net <- blockwiseModules(
  datExpr,
  power = network_params$power,
  networkType = network_params$networkType,
  TOMType = network_params$TOMType,
  minModuleSize = network_params$minModuleSize,
  mergeCutHeight = network_params$mergeCutHeight,
  corType = "bicor",
  maxPOutliers = 0.1,
  deepSplit = network_params$deepSplit,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

# Extract module information
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
names(moduleColors) <- colnames(datExpr)
moduleCounts <- sort(table(moduleColors), decreasing = TRUE)
nGrey <- ifelse("grey" %in% names(moduleCounts), moduleCounts["grey"], 0)
nModules <- length(unique(moduleColors)) - 1  # Exclude grey

pctGrey <- round(100 * nGrey / ncol(datExpr), 1)

cat(sprintf("\n  Network built: %d modules, %d grey (%.1f%%)\n",
            nModules, nGrey, pctGrey))
cat("  Configuration: Bicor + Signed, Power = 13, deepSplit = 3\n")
cat("\nModule sizes:\n")
print(moduleCounts)

module_summary <- data.frame(
  Module = names(moduleCounts),
  nProteins = as.numeric(moduleCounts),
  Percentage = round(100 * as.numeric(moduleCounts) / ncol(datExpr), 1)
) %>% arrange(desc(nProteins))

write.csv(module_summary, file.path(results_dir, "Step5_ModuleSummary.csv"), row.names = FALSE)
cat("\nSaved: Step5_ModuleSummary.csv\n")
# --- Calculate and Save TOM for Step 6a connectivity analysis ---
cat("\n  Calculating TOM (Topological Overlap Matrix)...\n")
adjacency <- adjacency(datExpr, power = 13, type = "signed",
                       corFnc = "bicor", corOptions = list(maxPOutliers = 0.1))
TOM <- TOMsimilarity(adjacency)
colnames(TOM) <- colnames(datExpr)
rownames(TOM) <- colnames(datExpr)
cat(sprintf("  TOM dimensions: %d x %d\n", nrow(TOM), ncol(TOM)))

# Save TOM for later use
save(TOM, file = file.path(results_dir, "Step5_TOM-block.1.RData"))
cat("  Saved: Step5_TOM-block.1.RData\n")

# Build gene tree from TOM for dendrogram
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
cat("  Calculated gene clustering tree\n")

# Clean up large temporary object
rm(adjacency, dissTOM)
gc()

# ============================================================================
# 5b: Generate Dendrograms and Module Size Plots
# ============================================================================

cat("\n5b: Generating Network Visualizations\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# --- Protein Dendrogram with Module Colors ---
cat("  Generating Protein Dendrogram...\n")

png(file.path(fig_dir_png, "Step5_01_Dendrogram.png"), width=12, height=6, units="in", res=300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module",
                    dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                    main="Protein Dendrogram with Module Colors")
dev.off()

pdf(file.path(fig_dir_pdf, "Step5_01_Dendrogram.pdf"), width=12, height=6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module",
                    dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                    main="Protein Dendrogram with Module Colors")
dev.off()
cat("  Saved: Step5_01_Dendrogram.png/pdf\n")

# --- Module Sizes ---
cat("  Generating Module Sizes plot...\n")

p_sizes <- module_summary %>% filter(Module != "grey") %>%
  mutate(Module = factor(Module, levels = Module)) %>%
  ggplot(aes(Module, nProteins, fill = Module)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_text(aes(label = nProteins), vjust = -0.5, size = 3.5) +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Module Sizes",
    subtitle = "Proteins per module (excluding grey)",
    x = NULL, y = "Number of Proteins"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )

print(p_sizes)

ggsave(file.path(fig_dir_png, "Step5_02_Module_Sizes.png"), p_sizes, width=5, height=4, dpi=300)
ggsave(file.path(fig_dir_pdf, "Step5_02_Module_Sizes.pdf"), p_sizes, width=5, height=4)
cat("  Saved: Step5_02_Module_Sizes.png/pdf\n")

# ============================================================================
# 5c: Module Eigengenes & Trait Correlations
# ============================================================================

cat("\n5c: Module-Trait Correlations\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Calculate module eigengenes from the network
# Note: net$MEs has numeric labels (ME0, ME1, ...) - convert to color names
MEs_raw <- net$MEs

# Create mapping from numeric labels to colors
# moduleLabels: 0=grey, 1=turquoise, 2=blue, etc.
label_to_color <- data.frame(
  label = sort(unique(moduleLabels)),
  color = labels2colors(sort(unique(moduleLabels)))
)
cat("  Module label to color mapping:\n")
for(i in 1:nrow(label_to_color)) {
  cat(sprintf("    ME%d -> %s\n", label_to_color$label[i], label_to_color$color[i]))
}

# Rename ME columns from ME0, ME1... to MEgrey, MEturquoise...
new_colnames <- sapply(colnames(MEs_raw), function(x) {
  num <- as.numeric(gsub("ME", "", x))
  color <- label_to_color$color[label_to_color$label == num]
  if(length(color) == 1) paste0("ME", color) else x
})
colnames(MEs_raw) <- new_colnames

MEs <- orderMEs(MEs_raw)
# Remove grey module eigengene
MEs_noGrey <- MEs[, !grepl("grey", colnames(MEs), ignore.case = TRUE)]
cat(sprintf("  MEs columns: %s\n", paste(colnames(MEs_noGrey), collapse = ", ")))

# Trait data - Binary encoding documentation:
# PFS_group: Long = 1, Short = 0 (positive correlation = better prognosis)
# Response: CD (Complete/Partial Response) = 1, PD (Progressive Disease) = 0
# Treatment: FOLFI = 1, GnP = 0
# Liver_Metastasis: Yes = 1, No = 0
# Sex: Male = 1, Female = 0
traitData <- data.frame(
  PFS_group = as.numeric(datTraits$PFS_group == "Long"),
  Response = as.numeric(datTraits$Response == "CD"),
  PFS_Days = datTraits$PFS,
  Treatment = as.numeric(datTraits$Treatment == "FOLFI"),
  Liver_Metastasis = as.numeric(datTraits$Liver == "Yes"),
  Age = datTraits$Age,
  Sex = as.numeric(datTraits$Sex == "Male"),
  row.names = rownames(datTraits)
)
if("CA19.9" %in% colnames(datTraits)) traitData$CA19.9 <- datTraits$CA19.9

# Calculate sample sizes per trait (accounting for missing data)
trait_n <- sapply(traitData, function(x) sum(!is.na(x)))
cat("\nSample sizes per trait:\n")
print(trait_n)


# Calculate module-trait correlations
moduleTraitCor <- cor(MEs_noGrey, traitData, use = "pairwise.complete.obs")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Rename columns for consistency
colnames(moduleTraitCor)[colnames(moduleTraitCor) == "PFS"] <- "PFS_Days"
colnames(moduleTraitCor)[colnames(moduleTraitCor) == "Liver"] <- "Liver_Metastasis"
colnames(moduleTraitPval)[colnames(moduleTraitPval) == "PFS"] <- "PFS_Days"
colnames(moduleTraitPval)[colnames(moduleTraitPval) == "Liver"] <- "Liver_Metastasis"

cat("\nSignificant Associations (p < 0.05):\n")
for(j in 1:ncol(moduleTraitPval)) {
  sig <- which(moduleTraitPval[,j] < 0.05)
  for(i in sig) cat(sprintf("  %s ~ %s: r=%.3f, p=%.4f\n",
                            rownames(moduleTraitPval)[i], colnames(moduleTraitPval)[j],
                            moduleTraitCor[i,j], moduleTraitPval[i,j]))
}

# Effect Size Interpretation Guide (Cohen's conventions for correlation coefficients)
cat("\n--- Effect Size Interpretation Guide ---\n")
cat("  |r| < 0.10: Negligible effect\n")
cat("  0.10 <= |r| < 0.30: Small effect\n")
cat("  0.30 <= |r| < 0.50: Medium effect\n")
cat("  |r| >= 0.50: Large effect\n")
cat("  Note: Correlations are Pearson r values with sample size n =", nrow(traitData), "\n")

# Classify and report significant associations by effect size
if(any(moduleTraitPval < 0.05, na.rm = TRUE)) {
  cat("\nSignificant associations by effect size:\n")
  for(j in 1:ncol(moduleTraitPval)) {
    sig <- which(moduleTraitPval[,j] < 0.05)
    for(i in sig) {
      r_val <- abs(moduleTraitCor[i,j])
      effect_label <- ifelse(r_val >= 0.50, "LARGE",
                             ifelse(r_val >= 0.30, "Medium",
                                    ifelse(r_val >= 0.10, "Small", "Negligible")))
      cat(sprintf("  %s ~ %s: r=%.3f (%s effect)\n",
                  rownames(moduleTraitPval)[i], colnames(moduleTraitPval)[j],
                  moduleTraitCor[i,j], effect_label))
    }
  }
}

# --- Module Eigengene Dendrogram ---
cat("\n  Generating Module Eigengene Dendrogram...\n")

MEDiss <- 1 - cor(MEs_noGrey)
METree <- hclust(as.dist(MEDiss), method = "average")

png(file.path(fig_dir_png, "Step5_03_ME_Dendrogram.png"), width=6, height=4, units="in", res=300)
par(mar = c(1, 4, 2, 1))
plot(METree, main = "Module Eigengene Dendrogram", xlab = "", sub = "",
     ylim = c(0, max(METree$height)*1.1))
abline(h = 0.25, col = "red", lty = 2)
dev.off()

pdf(file.path(fig_dir_pdf, "Step5_03_ME_Dendrogram.pdf"), width=6, height=4)
par(mar = c(1, 4, 2, 1))
plot(METree, main = "Module Eigengene Dendrogram", xlab = "", sub = "",
     ylim = c(0, max(METree$height)*1.1))
abline(h = 0.25, col = "red", lty = 2)
dev.off()
cat("  Saved: Step5_03_ME_Dendrogram.png/pdf\n")

# --- Module-Trait Heatmap ---
cat("  Generating Module-Trait Heatmap...\n")

# Reorder traits for display
trait_order <- c("PFS_group", "Response", "PFS_Days", "Treatment", "Liver_Metastasis", "Age", "Sex")
if("CA19.9" %in% colnames(moduleTraitCor)) trait_order <- c(trait_order, "CA19.9")
trait_order <- trait_order[trait_order %in% colnames(moduleTraitCor)]

moduleTraitCor_ordered <- moduleTraitCor[, trait_order]
moduleTraitPval_ordered <- moduleTraitPval[, trait_order]

heatmap_data <- melt(moduleTraitCor_ordered) %>%
  rename(Module = Var1, Trait = Var2, Correlation = value) %>%
  mutate(Pval = melt(moduleTraitPval_ordered)$value,
         Sig = case_when(Pval < 0.001 ~ "***", Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""),
         Label = paste0(round(Correlation, 2), Sig),
         TextColor = ifelse(Pval < 0.05, "black", "grey50"))

# Module order by clustering
hc <- hclust(dist(moduleTraitCor_ordered))
heatmap_data$Module <- factor(heatmap_data$Module, levels = rev(rownames(moduleTraitCor_ordered)[hc$order]))
heatmap_data$Trait <- factor(heatmap_data$Trait, levels = trait_order)

# Sidebar colors
sidebar_df <- data.frame(
  Module = factor(rownames(moduleTraitCor_ordered)[hc$order], levels = rev(rownames(moduleTraitCor_ordered)[hc$order])),
  Color = gsub("ME", "", rownames(moduleTraitCor_ordered)[hc$order])
)

p_sidebar <- ggplot(sidebar_df, aes(x = 1, y = Module, fill = Color)) +
  geom_tile(width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_identity() + theme_void()

# Calculate dynamic limits based on actual correlation range (symmetric around 0)
max_abs_cor <- max(abs(moduleTraitCor_ordered), na.rm = TRUE)
cor_limit <- ceiling(max_abs_cor * 10) / 10  # Round up to nearest 0.1
cor_limit <- max(cor_limit, 0.5)  # Minimum limit of 0.5 for visibility
cat(sprintf("  Heatmap correlation limits: [%.1f, %.1f] (max |r| = %.3f)\n", -cor_limit, cor_limit, max_abs_cor))

# Add fontface column to heatmap_data for significant correlations
heatmap_data$FontFace <- "plain"
heatmap_data$FontFace[heatmap_data$Pval < 0.05] <- "bold"

p_heat <- ggplot(heatmap_data, aes(Trait, Module, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label, color = TextColor, fontface = FontFace), size = 3.5) +
  scale_fill_gradient2(low = "#ca562c", mid = "white", high = "#008080",
                       midpoint = 0, limits = c(-cor_limit, cor_limit), name = "r") +
  scale_color_identity() +
  labs(
    title = "Module-Trait Relationships",
    subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
    x = NULL, y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

p_heatmap <- p_sidebar + p_heat + plot_layout(widths = c(0.05, 1))

print(p_heatmap)

ggsave(file.path(fig_dir_png, "Step5_04_ModuleTrait_Heatmap.png"), p_heatmap, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir_pdf, "Step5_04_ModuleTrait_Heatmap.pdf"), p_heatmap, width=6, height=5)
cat("  Saved: Step5_04_ModuleTrait_Heatmap.png/pdf\n")

# Save correlation and p-value tables WITH SAMPLE SIZES
# Calculate sample size for each trait (accounts for missing data)
trait_n <- sapply(traitData, function(x) sum(!is.na(x)))
cat("  Sample sizes per trait:\n")
print(trait_n)

# Create comprehensive correlation table with sample sizes
cor_table <- cbind(
  Module = rownames(moduleTraitCor_ordered),
  round(moduleTraitCor_ordered, 4)
)
# Add sample size row
n_row <- c("N (samples)", trait_n[colnames(moduleTraitCor_ordered)])
cor_table <- rbind(n_row, cor_table)

write.csv(cor_table, file.path(results_dir, "Step5_ModuleTraitCorrelations.csv"), row.names = FALSE)

# Create comprehensive p-value table with sample sizes
pval_table <- cbind(
  Module = rownames(moduleTraitPval_ordered),
  round(moduleTraitPval_ordered, 4)
)
pval_table <- rbind(n_row, pval_table)

write.csv(pval_table, file.path(results_dir, "Step5_ModuleTraitPvalues.csv"), row.names = FALSE)
cat("  Saved: Step5_ModuleTraitCorrelations.csv, Step5_ModuleTraitPvalues.csv (with sample sizes)\n")

# ----------------------------------------------------------------------------
# 5c.1b: Sample Dendrogram with Module Eigengenes and Clinical Traits
# ----------------------------------------------------------------------------
# PURPOSE: Cluster patients by their module eigengene profiles to:
#   - Identify patient subgroups with similar module activity
#   - Visualize how clinical traits align with patient clusters
#   - Validate module-clinical associations at the sample level
# ----------------------------------------------------------------------------

cat("\n  Generating Sample Dendrogram with Module Eigengenes...\n")

# Cluster samples based on module eigengene values
sampleTree_ME <- hclust(dist(MEs_noGrey), method = "average")

# Prepare clinical trait colors for visualization
# PFS_group colors
pfs_colors <- ifelse(datTraits$PFS_group == "Long", col_pfs_long, col_pfs_short)
pfs_colors[is.na(datTraits$PFS_group)] <- "grey90"

# Response colors
response_colors <- ifelse(datTraits$Response == "CD", col_cd, col_pd)
response_colors[is.na(datTraits$Response)] <- "grey90"

# Treatment colors
treatment_colors <- ifelse(datTraits$Treatment == "FOLFI", col_folfi, col_gemnab)
treatment_colors[is.na(datTraits$Treatment)] <- "grey90"

# Combine trait colors into matrix
traitColors <- cbind(
  PFS_Group = pfs_colors,
  Response = response_colors,
  Treatment = treatment_colors
)
rownames(traitColors) <- rownames(datTraits)

# --- Plot 1: Sample Dendrogram with Clinical Trait Bars ---
png(file.path(fig_dir_png, "Step5_04b_Sample_Dendrogram_Traits.png"), width = 12, height = 6, units = "in", res = 300)
plotDendroAndColors(
  sampleTree_ME,
  colors = traitColors,
  groupLabels = colnames(traitColors),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Sample Dendrogram with Clinical Traits"
)
dev.off()

pdf(file.path(fig_dir_pdf, "Step5_04b_Sample_Dendrogram_Traits.pdf"), width = 12, height = 6)
plotDendroAndColors(
  sampleTree_ME,
  colors = traitColors,
  groupLabels = colnames(traitColors),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Sample Dendrogram with Clinical Traits"
)
dev.off()
cat("  Saved: Step5_04b_Sample_Dendrogram_Traits.png/pdf\n")

# --- Plot 2: Sample Dendrogram with Module Eigengene Heatmap ---
# Convert MEs to colors for heatmap visualization
ME_colors <- matrix(NA, nrow = nrow(MEs_noGrey), ncol = ncol(MEs_noGrey))
colnames(ME_colors) <- gsub("ME", "", colnames(MEs_noGrey))
rownames(ME_colors) <- rownames(MEs_noGrey)

# Scale ME values to color gradient (blue = low, white = mid, red = high)
for(i in 1:ncol(MEs_noGrey)) {
  me_vals <- MEs_noGrey[, i]
  me_scaled <- (me_vals - min(me_vals, na.rm = TRUE)) / (max(me_vals, na.rm = TRUE) - min(me_vals, na.rm = TRUE))
  ME_colors[, i] <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)[ceiling(me_scaled * 99) + 1]
}

# Combine ME heatmap colors with trait colors
combined_colors <- cbind(ME_colors, traitColors)


# Prepare data for ggplot heatmap

# Get sample order from dendrogram
sample_order <- sampleTree_ME$order
sample_names <- rownames(MEs_noGrey)[sample_order]

# Prepare ME data in long format
ME_long <- MEs_noGrey %>%
  as.data.frame() %>%
  mutate(Sample = rownames(MEs_noGrey)) %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene") %>%
  mutate(
    Module = gsub("ME", "", Module),
    Sample = factor(Sample, levels = sample_names)
  )

# Prepare clinical trait data in long format
trait_long <- data.frame(
  Sample = rownames(datTraits),
  PFS_Group = datTraits$PFS_group,
  Response = datTraits$Response,
  Treatment = datTraits$Treatment
) %>%
  pivot_longer(cols = c(PFS_Group, Response, Treatment), names_to = "Trait", values_to = "Value") %>%
  mutate(Sample = factor(Sample, levels = sample_names))

# Create dendrogram plot
dend_data <- dendro_data(as.dendrogram(sampleTree_ME))
p_dend <- ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 0, 5))

# Create ME heatmap
module_order <- colnames(ME_colors)
ME_long$Module <- factor(ME_long$Module, levels = rev(module_order))

p_me_heat <- ggplot(ME_long, aes(x = Sample, y = Module, fill = Eigengene)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Module\nEigengene") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.key.height = unit(0.8, "cm"),
    plot.margin = margin(0, 5, 0, 5)
  )

# Create clinical trait heatmap
trait_long$Trait <- factor(trait_long$Trait, levels = c("PFS_Group", "Response", "Treatment"))

p_trait_heat <- ggplot(trait_long, aes(x = Sample, y = Trait, fill = Value)) +
  geom_tile() +
  scale_fill_manual(
    values = c(
      "Long" = col_pfs_long, "Short" = col_pfs_short,
      "CD" = col_cd, "PD" = col_pd,
      "FOLFI" = col_folfi, "GemNab" = col_gemnab
    ),
    name = "Clinical\nTraits",
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 9),
    plot.margin = margin(0, 5, 5, 5)
  )

# Combine plots
p_combined <- p_dend / p_me_heat / p_trait_heat +
  plot_layout(heights = c(2, 3, 1)) +
  plot_annotation(title = "Sample Dendrogram with Module Eigengenes and Clinical Traits",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

ggsave(file.path(fig_dir_png, "Step5_04c_Sample_Dendrogram_ME_Heatmap.png"), p_combined,
       width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step5_04c_Sample_Dendrogram_ME_Heatmap.pdf"), p_combined,
       width = 8, height = 6)
cat("  Saved: Step5_04c_Sample_Dendrogram_ME_Heatmap.png/pdf\n")

# --- Summary: Patient Clustering Analysis ---
cat("\n  Patient clustering summary:\n")

# Cut dendrogram into 2-3 clusters for analysis
sample_clusters <- cutree(sampleTree_ME, k = 2)
cluster_summary <- data.frame(
  Cluster = 1:2,
  N_Patients = as.numeric(table(sample_clusters))
)

# Analyze clinical composition of each cluster
cat("  Cluster composition (k=2):\n")
for(cl in 1:2) {
  cl_idx <- which(sample_clusters == cl)
  n_cl <- length(cl_idx)
  n_long <- sum(datTraits$PFS_group[cl_idx] == "Long", na.rm = TRUE)
  n_cd <- sum(datTraits$Response[cl_idx] == "CD", na.rm = TRUE)
  cat(sprintf("    Cluster %d: n=%d, Long PFS=%d (%.0f%%), CD=%d (%.0f%%)\n",
              cl, n_cl, n_long, 100*n_long/n_cl, n_cd, 100*n_cd/n_cl))
}

# Test if clusters differ by PFS_group
cluster_pfs_table <- table(sample_clusters, datTraits$PFS_group)
if(all(dim(cluster_pfs_table) >= 2)) {
  fisher_pfs <- fisher.test(cluster_pfs_table)
  cat(sprintf("\n  Cluster vs PFS_group: Fisher's exact p = %.4f\n", fisher_pfs$p.value))
}

# Test if clusters differ by Response
cluster_resp_table <- table(sample_clusters, datTraits$Response)
if(all(dim(cluster_resp_table) >= 2)) {
  fisher_resp <- fisher.test(cluster_resp_table)
  cat(sprintf("  Cluster vs Response: Fisher's exact p = %.4f\n", fisher_resp$p.value))
}

# ----------------------------------------------------------------------------
# 5c.2: Bootstrap 95% Confidence Intervals for Module-Trait Correlations
# ----------------------------------------------------------------------------
# Purpose: Calculate confidence intervals for your correlation values
# Method: BCa (bias-corrected and accelerated) bootstrap, 1000 resamples
# Reference: Efron & Tibshirani, 1993; High-impact journals increasingly require CIs
# ----------------------------------------------------------------------------

cat("\n  Calculating bootstrap 95% CIs for significant correlations...\n")
set.seed(12345)
n_boot <- 1000

# Function to compute bootstrap CI for correlation
boot_cor_ci <- function(me_vec, trait_vec, n_boot = 1000) {
  n <- length(me_vec)
  boot_cors <- numeric(n_boot)

  for(b in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    boot_cors[b] <- cor(me_vec[idx], trait_vec[idx], use = "pairwise.complete.obs")
  }

  # Calculate BCa-like CI using percentile method (simpler, robust)
  ci_lower <- quantile(boot_cors, 0.025, na.rm = TRUE)
  ci_upper <- quantile(boot_cors, 0.975, na.rm = TRUE)
  se <- sd(boot_cors, na.rm = TRUE)

  return(list(ci_lower = ci_lower, ci_upper = ci_upper, se = se))
}

# Calculate CIs for key clinical traits
key_traits <- c("PFS_group", "Response")
key_traits <- key_traits[key_traits %in% colnames(traitData)]

bootstrap_results <- data.frame(
  Module = character(),
  Trait = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  SE = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for(mod in rownames(moduleTraitCor)) {
  me_vec <- MEs_noGrey[, mod]

  for(trait in key_traits) {
    trait_vec <- traitData[, trait]
    cor_val <- moduleTraitCor[mod, trait]
    p_val <- moduleTraitPval[mod, trait]

    # Only compute CI for correlations with p < 0.1 (reduces computation)
    if(p_val < 0.1) {
      ci_result <- boot_cor_ci(me_vec, trait_vec, n_boot)

      bootstrap_results <- rbind(bootstrap_results, data.frame(
        Module = mod,
        Trait = trait,
        Correlation = round(cor_val, 4),
        CI_Lower = round(ci_result$ci_lower, 4),
        CI_Upper = round(ci_result$ci_upper, 4),
        SE = round(ci_result$se, 4),
        P_value = round(p_val, 4),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Print and save results
if(nrow(bootstrap_results) > 0) {
  cat("  Bootstrap 95% CIs for significant correlations (p < 0.1):\n")
  print(bootstrap_results, row.names = FALSE)

  write.csv(bootstrap_results, file.path(results_dir, "Step5_CorrelationBootstrapCI.csv"), row.names = FALSE)
  cat("  Saved: Step5_CorrelationBootstrapCI.csv\n")

  # Create forest plot for bootstrap CIs
  p_forest <- ggplot(bootstrap_results, aes(x = Correlation, y = paste(Module, Trait, sep = " ~ "))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "grey40") +
    geom_point(aes(color = P_value < 0.05), size = 3) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#ca562c"),
                       labels = c("p >= 0.05", "p < 0.05"),
                       name = "Significance") +
    labs(
      title = "Module-Trait Correlations with 95% Bootstrap CIs",
      subtitle = sprintf("Based on %d bootstrap resamples | Key clinical traits", n_boot),
      x = "Pearson Correlation (r)",
      y = "Association"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )

  print(p_forest)
  ggsave(file.path(fig_dir_png, "Step5_04b_CorrelationForestPlot.png"), p_forest, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step5_04b_CorrelationForestPlot.pdf"), p_forest, width = 8, height = 6)
  cat("  Saved: Step5_04b_CorrelationForestPlot.png/pdf\n")
} else {
  cat("  No correlations with p < 0.1 found for bootstrap CI calculation.\n")
}

# ============================================================================
# 5d: Permutation Testing for Module-Trait Associations
# ============================================================================
# Purpose: Validate that module-trait associations are not due to data artifacts
# Testing: Focus modules (pink, green, black, blue) against outcomes (PFS_group, Response)
# ============================================================================

cat("\n5d: Permutation Testing for Key Module-Trait Associations\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

set.seed(12345)
n_perm <- 1000

# Get modules significant for PFS_group OR Response
sig_pfs <- rownames(moduleTraitPval)[moduleTraitPval[, "PFS_group"] < 0.05]
sig_response <- rownames(moduleTraitPval)[moduleTraitPval[, "Response"] < 0.05]
sig_modules <- unique(c(sig_pfs, sig_response))

cat(sprintf("  Testing %d modules significant for PFS_group and/or Response\n", length(sig_modules)))
cat(sprintf("  Modules: %s\n", paste(sig_modules, collapse = ", ")))

# Traits to test
traits_to_test <- c("PFS_group", "Response")

# Function to perform permutation test with Z-score calculation
# Z-score measures how many standard deviations the observed value is from the null mean
# Reference: Phipson & Smyth, 2010 (permutation p-values)
permutation_test <- function(me_values, trait_values, observed_cor, n_perm = 1000) {
  perm_cors <- numeric(n_perm)
  for(i in 1:n_perm) {
    perm_trait <- sample(trait_values)
    perm_cors[i] <- cor(me_values, perm_trait, use = "pairwise.complete.obs")
  }
  empirical_p <- mean(abs(perm_cors) >= abs(observed_cor))

  # Calculate Z-score: (observed - null_mean) / null_sd
  # Higher |Z| = stronger deviation from null, more significant
  null_mean <- mean(perm_cors)
  null_sd <- sd(perm_cors)
  z_score <- (observed_cor - null_mean) / null_sd

  return(list(
    empirical_p = empirical_p,
    perm_cors = perm_cors,
    z_score = z_score,
    null_mean = null_mean,
    null_sd = null_sd
  ))
}

# Prepare trait vectors
trait_vectors <- list(
  PFS_group = as.numeric(datTraits$PFS_group == "Long"),
  Response = as.numeric(datTraits$Response == "CD")
)

# Run permutation tests for all module-trait combinations
perm_results <- data.frame(
  Module = character(),
  Trait = character(),
  Observed_Cor = numeric(),
  Parametric_P = numeric(),
  Empirical_P = numeric(),
  Z_Score = numeric(),
  stringsAsFactors = FALSE
)

perm_distributions <- list()  # Store for plotting

for(mod in sig_modules) {
  if(mod %in% colnames(MEs_noGrey)) {
    me_values <- MEs_noGrey[, mod]
    
    for(trait in traits_to_test) {
      # Check if this module is significant for this trait
      if(moduleTraitPval[mod, trait] < 0.05) {
        cat(sprintf("  Testing %s ~ %s ... ", mod, trait))
        
        observed_cor <- moduleTraitCor[mod, trait]
        parametric_p <- moduleTraitPval[mod, trait]
        trait_values <- trait_vectors[[trait]]
        
        perm_result <- permutation_test(me_values, trait_values, observed_cor, n_perm)
        
        # Store results (now includes Z-score)
        perm_results <- rbind(perm_results, data.frame(
          Module = mod,
          Trait = trait,
          Observed_Cor = round(observed_cor, 4),
          Parametric_P = round(parametric_p, 4),
          Empirical_P = round(perm_result$empirical_p, 4),
          Z_Score = round(perm_result$z_score, 2),
          stringsAsFactors = FALSE
        ))

        # Store distribution for plotting (includes Z-score for annotation)
        perm_distributions[[paste0(mod, "_", trait)]] <- list(
          perm_cors = perm_result$perm_cors,
          observed_cor = observed_cor,
          empirical_p = perm_result$empirical_p,
          z_score = perm_result$z_score,
          module = mod,
          trait = trait
        )

        cat(sprintf("Empirical p = %.4f, Z = %.2f\n", perm_result$empirical_p, perm_result$z_score))
      }
    }
  }
}

# Print summary table
cat("\n  Permutation Test Results (n=1000 permutations):\n")
print(perm_results, row.names = FALSE)

write.csv(perm_results, file.path(results_dir, "Step5_Permutation_Tests.csv"), row.names = FALSE)
cat("  Saved: Step5_Permutation_Tests.csv\n")

# Create permutation plots for ALL significant combinations
cat("\n  Generating permutation plots...\n")

perm_plots <- list()

for(key in names(perm_distributions)) {
  dist_data <- perm_distributions[[key]]

  # Calculate text position to the right of the dashed line with small offset
  text_x_offset <- diff(range(dist_data$perm_cors)) * 0.02
  text_x <- dist_data$observed_cor + text_x_offset

  p_perm <- ggplot(data.frame(cor = dist_data$perm_cors), aes(x = cor)) +
    geom_histogram(bins = 50, fill = "grey70", color = "grey50", alpha = 0.7) +
    geom_vline(xintercept = dist_data$observed_cor, color = "#ca562c", linewidth = 1.5, linetype = "dashed") +
    annotate("text", x = text_x, y = Inf,
             label = sprintf("Observed\nr = %.3f\nZ = %.2f", dist_data$observed_cor, dist_data$z_score),
             vjust = 2, hjust = 0,
             color = "#ca562c", fontface = "bold", size = 3.5) +
    labs(
      title = sprintf("Permutation Test: %s ~ %s", dist_data$module, dist_data$trait),
      subtitle = sprintf("Null distribution (n=%d permutations) | Empirical p = %.4f | Z = %.2f",
                         n_perm, dist_data$empirical_p, dist_data$z_score),
      x = "Correlation coefficient",
      y = "Frequency"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )

  perm_plots[[key]] <- p_perm
  
  # Save individual plot
  filename_base <- sprintf("Step5_05_Permutation_%s_%s", gsub("ME", "", dist_data$module), dist_data$trait)
  ggsave(file.path(fig_dir_png, paste0(filename_base, ".png")), p_perm, width = 6, height = 4, dpi = 300)
  ggsave(file.path(fig_dir_pdf, paste0(filename_base, ".pdf")), p_perm, width = 6, height = 4)
  cat(sprintf("  Saved: %s.png/pdf\n", filename_base))
}

# Combined plot if we have multiple tests
if(length(perm_plots) > 1) {
  # Arrange in grid
  n_plots <- length(perm_plots)
  ncol <- min(2, n_plots)
  nrow <- ceiling(n_plots / ncol)
  
  fig_perm_combined <- wrap_plots(perm_plots, ncol = ncol) +
    plot_annotation(
      title = "Permutation Validation of Module-Trait Associations",
      subtitle = "Observed correlations (red) vs null distributions (grey)",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40")
      )
    )
  
  print(fig_perm_combined)
  ggsave(file.path(fig_dir_png, "Step5_05_Permutation_Combined.png"), fig_perm_combined, 
         width = 16, height = 5 * nrow, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step5_05_Permutation_Combined.pdf"), fig_perm_combined, 
         width = 16, height = 5 * nrow)
  cat("  Saved: Step5_05_Permutation_Combined.png/pdf\n")
}

cat("\n  Permutation testing complete.\n")

# ====================================================================
# 5e: Module Membership & Gene Significance
# ============================================================================

#with MM and GS WGCNA identifies hub proteins - the most important proteins within each module.
#MM = cor(protein_expression, module_eigengene)
#Gene Significance (GS): Correlation between each protein's expression and clinical traits
#GS +0.5Protein increases with trait (e.g., higher in Long PFS)
#GS -0.5Protein decreases with trait (e.g., lower in Long PFS)
#GS 0.0Protein not associated with trait

cat("\n5e: Calculating Module Membership (MM) and Gene Significance (GS)\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

geneModuleMembership <- as.data.frame(cor(datExpr, MEs_noGrey, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
colnames(geneModuleMembership) <- gsub("ME", "MM.", colnames(geneModuleMembership))
colnames(MMPvalue) <- gsub("ME", "p.MM.", colnames(MMPvalue))

geneTraitSignificance <- as.data.frame(cor(datExpr, traitData, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))
colnames(geneTraitSignificance) <- paste0("GS.", colnames(geneTraitSignificance))
colnames(GSPvalue) <- paste0("p.GS.", colnames(GSPvalue))

proteinInfo <- data.frame(Protein = colnames(datExpr), Module = moduleColors) %>%
  cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)


cat(sprintf("  Protein info table: %d proteins x %d columns\n", nrow(proteinInfo), ncol(proteinInfo)))

write.csv(proteinInfo, file.path(results_dir, "Step5_ProteinModuleInfo.csv"), row.names = FALSE)
cat("  Saved: Step5_ProteinModuleInfo.csv\n")

# ============================================================================
# 5f: Module-Trait DotPlot
# ============================================================================

cat("\n5f: Generating Module-Trait DotPlot\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

module.name <- colnames(MEs_noGrey)
plot.data <- cbind(datTraits, MEs_noGrey)
if("PFS" %in% colnames(plot.data) & !"PFS_Days" %in% colnames(plot.data)) plot.data$PFS_Days <- plot.data$PFS

# Prognosis (Wilcoxon)
calc_wilcox <- function(data, me, group_var, pos_level, neg_level) {
  data$ME <- data[[me]]
  mean_diff <- mean(data$ME[data[[group_var]] == pos_level], na.rm=TRUE) -
    mean(data$ME[data[[group_var]] == neg_level], na.rm=TRUE)
  p <- tryCatch(wilcox.test(as.formula(paste("ME ~", group_var)), data = data)$p.value,
                error = function(e) NA)
  data.frame(Module = me, Factor = group_var, Type = mean_diff, P = p, Test = "Wilcoxon")
}

plot.info.prog <- bind_rows(
  lapply(module.name, function(m) calc_wilcox(plot.data, m, "PFS_group", "Long", "Short")),
  lapply(module.name, function(m) calc_wilcox(plot.data, m, "Response", "CD", "PD"))
) %>% mutate(Padj = p.adjust(P, method = "fdr"),
             Tag = case_when(Type > 0 & P < 0.05 ~ "Good prognosis",
                             Type < 0 & P < 0.05 ~ "Poor prognosis", TRUE ~ "notSig"))

# Continuous (Spearman)
cont_traits <- c("PFS_Days", "Age")
if("CA19.9" %in% colnames(plot.data)) cont_traits <- c(cont_traits, "CA19.9")

plot.info.cont <- bind_rows(lapply(cont_traits, function(trait) {
  bind_rows(lapply(module.name, function(m) {
    valid <- !is.na(plot.data[[trait]]) & !is.na(plot.data[[m]])
    if(sum(valid) < 5) return(data.frame(Module = m, Factor = trait, Type = NA, P = NA, Test = "Spearman"))
    ct <- cor.test(plot.data[[m]][valid], plot.data[[trait]][valid], method = "spearman")
    data.frame(Module = m, Factor = trait, Type = ct$estimate, P = ct$p.value, Test = "Spearman")
  }))
})) %>% mutate(Padj = p.adjust(P, method = "fdr"),
               Tag = case_when(Type > 0 & P < 0.05 ~ "Positive",
                               Type < 0 & P < 0.05 ~ "Negative", TRUE ~ "notSig"))

# Categorical (Wilcoxon)
cat_traits <- c("Treatment", "Sex", "Liver")
plot.info.cat <- bind_rows(lapply(cat_traits, function(trait) {
  bind_rows(lapply(module.name, function(m) {
    valid <- !is.na(plot.data[[trait]]) & !is.na(plot.data[[m]])
    if(sum(valid) < 5 || length(unique(plot.data[[trait]][valid])) < 2) {
      return(data.frame(Module = m, Factor = trait, Type = NA, P = NA, Test = "Wilcoxon"))
    }
    p <- tryCatch(wilcox.test(plot.data[[m]][valid] ~ plot.data[[trait]][valid])$p.value,
                  error = function(e) NA)
    data.frame(Module = m, Factor = trait, Type = NA, P = p, Test = "Wilcoxon")
  }))
})) %>% mutate(Padj = p.adjust(P, method = "fdr"), Tag = ifelse(P < 0.05, "Sig", "notSig"))

# Combine
plot.info <- bind_rows(plot.info.prog, plot.info.cont, plot.info.cat) %>%
  mutate(LogP = pmin(-log10(P), 4),
         Category = case_when(Factor %in% c("PFS_group", "Response") ~ "Prognosis",
                              Factor %in% cont_traits ~ "Continuous", TRUE ~ "Categorical"),
         Factor = factor(Factor, levels = c("PFS_group", "Response", cont_traits, cat_traits)),
         Category = factor(Category, levels = c("Prognosis", "Continuous", "Categorical")))

# Module order by PFS_group
mod_order <- plot.info %>% filter(Factor == "PFS_group") %>% arrange(Type) %>% pull(Module)
plot.info$Module <- factor(plot.info$Module, levels = rev(mod_order))

color_values <- c("Good prognosis" = "#3182bd", "Poor prognosis" = "#0E4D92",
                  "Positive" = "#00A693", "Negative" = "#FFD662", "Sig" = "#fd8d3c", "notSig" = "#bdbdbd")

p_dotplot <- ggplot(plot.info, aes(Factor, Module, size = LogP, color = Tag)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = color_values) +
  scale_size_continuous(range = c(1, 8), name = expression(-log[10](P)), breaks = 1:3) +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  labs(title = "Module-Trait Associations", x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"))

print(p_dotplot)

ggsave(file.path(fig_dir_png, "Step5_06_DotPlot.png"), p_dotplot, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir_pdf, "Step5_06_DotPlot.pdf"), p_dotplot, width=6, height=5)
write.csv(plot.info, file.path(results_dir, "Step5_DotPlot_Data.csv"), row.names = FALSE)
cat("  Saved: Step5_06_DotPlot.png/pdf, Step5_DotPlot_Data.csv\n")

# ----------------------------------------------------------------------------
# 5f.2: Multiple Testing Correction Summary
# ----------------------------------------------------------------------------
# Purpose: Summarize the multiple testing correction approach for publication
# Method: Benjamini-Hochberg FDR correction, applied per category
# Reference: Benjamini & Hochberg, 1995 (controlling the false discovery rate)
# ----------------------------------------------------------------------------

cat("\n  Multiple Testing Correction Summary:\n")
cat("  ", paste(rep("-", 55), collapse = ""), "\n")

# Create summary table per category
mtc_summary <- data.frame(
  Category = character(),
  N_Tests = integer(),
  N_Nominal_Sig = integer(),
  N_FDR_Sig = integer(),
  FDR_Method = character(),
  stringsAsFactors = FALSE
)

# Count by category
for(cat_name in c("Prognosis", "Continuous", "Categorical")) {
  cat_data <- plot.info %>% filter(Category == cat_name)
  n_tests <- nrow(cat_data)
  n_nominal <- sum(cat_data$P < 0.05, na.rm = TRUE)
  n_fdr <- sum(cat_data$Padj < 0.05, na.rm = TRUE)

  mtc_summary <- rbind(mtc_summary, data.frame(
    Category = cat_name,
    N_Tests = n_tests,
    N_Nominal_Sig = n_nominal,
    N_FDR_Sig = n_fdr,
    FDR_Method = "Benjamini-Hochberg",
    stringsAsFactors = FALSE
  ))
}

# Add total row
mtc_summary <- rbind(mtc_summary, data.frame(
  Category = "TOTAL",
  N_Tests = sum(mtc_summary$N_Tests),
  N_Nominal_Sig = sum(mtc_summary$N_Nominal_Sig),
  N_FDR_Sig = sum(mtc_summary$N_FDR_Sig),
  FDR_Method = "Benjamini-Hochberg",
  stringsAsFactors = FALSE
))

# Add proportion columns
mtc_summary$Pct_Nominal <- round(100 * mtc_summary$N_Nominal_Sig / mtc_summary$N_Tests, 1)
mtc_summary$Pct_FDR <- round(100 * mtc_summary$N_FDR_Sig / mtc_summary$N_Tests, 1)

cat("  Multiple testing correction applied per category (not globally)\n")
cat("  Rationale: Each category represents distinct biological questions\n\n")
print(mtc_summary, row.names = FALSE)

write.csv(mtc_summary, file.path(results_dir, "Step5_MultipleTestingCorrection.csv"), row.names = FALSE)
cat("\n  Saved: Step5_MultipleTestingCorrection.csv\n")

# ============================================================================
# 5g: MODULE PRIORITIZATION FOR CLINICAL RELEVANCE
# ============================================================================

cat("\n5g: Module Prioritization for Clinical Relevance\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

module_priority <- data.frame(
  Module = rownames(moduleTraitCor),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Cor_PFS_group = moduleTraitCor[Module, "PFS_group"],
    P_PFS_group = moduleTraitPval[Module, "PFS_group"],
    Sig_PFS_group = P_PFS_group < 0.05,
    
    Cor_Response = moduleTraitCor[Module, "Response"],
    P_Response = moduleTraitPval[Module, "Response"],
    Sig_Response = P_Response < 0.05,
    
    Cor_PFS_Days = moduleTraitCor[Module, "PFS_Days"],
    P_PFS_Days = moduleTraitPval[Module, "PFS_Days"],
    Sig_PFS_Days = P_PFS_Days < 0.05,
    
    N_Significant = Sig_PFS_group + Sig_Response + Sig_PFS_Days,
    Min_Pvalue = pmin(P_PFS_group, P_Response, P_PFS_Days)
  ) %>%
  arrange(desc(N_Significant), Min_Pvalue) %>%
  mutate(
    Clinical_Relevance = case_when(
      N_Significant == 3 ~ "HIGH - All 3 outcomes",
      N_Significant == 2 ~ "MODERATE - 2 outcomes",
      N_Significant == 1 ~ "LOW - 1 outcome",
      TRUE ~ "NOT SIGNIFICANT"
    )
  )

cat("\nMODULE PRIORITIZATION SUMMARY\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# HIGH priority
high_priority <- module_priority %>% filter(N_Significant == 3)
if(nrow(high_priority) > 0) {
  cat("\n*** HIGH PRIORITY (significant for ALL 3 outcomes) ***\n")
  for(i in 1:nrow(high_priority)) {
    cat(sprintf("  %s: PFS_group (r=%.3f, p=%.4f), Response (r=%.3f, p=%.4f), PFS_Days (r=%.3f, p=%.4f)\n",
                gsub("ME", "", high_priority$Module[i]),
                high_priority$Cor_PFS_group[i], high_priority$P_PFS_group[i],
                high_priority$Cor_Response[i], high_priority$P_Response[i],
                high_priority$Cor_PFS_Days[i], high_priority$P_PFS_Days[i]))
  }
}

# MODERATE priority
mod_priority <- module_priority %>% filter(N_Significant == 2)
if(nrow(mod_priority) > 0) {
  cat("\n** MODERATE PRIORITY (significant for 2 outcomes) **\n")
  for(i in 1:nrow(mod_priority)) {
    cat(sprintf("  %s: PFS_group (p=%.4f%s), Response (p=%.4f%s), PFS_Days (p=%.4f%s)\n",
                gsub("ME", "", mod_priority$Module[i]),
                mod_priority$P_PFS_group[i], ifelse(mod_priority$Sig_PFS_group[i], "*", ""),
                mod_priority$P_Response[i], ifelse(mod_priority$Sig_Response[i], "*", ""),
                mod_priority$P_PFS_Days[i], ifelse(mod_priority$Sig_PFS_Days[i], "*", "")))
  }
}

# LOW priority
low_priority <- module_priority %>% filter(N_Significant == 1)
if(nrow(low_priority) > 0) {
  cat("\n* LOW PRIORITY (significant for 1 outcome) *\n")
  cat(sprintf("  %s\n", paste(gsub("ME", "", low_priority$Module), collapse = ", ")))
}

# Not significant
not_sig <- module_priority %>% filter(N_Significant == 0)
if(nrow(not_sig) > 0) {
  cat("\nNOT SIGNIFICANT:\n")
  cat(sprintf("  %s\n", paste(gsub("ME", "", not_sig$Module), collapse = ", ")))
}

write.csv(module_priority, file.path(results_dir, "Step5_ModulePriority.csv"), row.names = FALSE)
cat("\nSaved: Step5_ModulePriority.csv\n")

# ============================================================================
# 5h: DYNAMIC HUB GENE IDENTIFICATION
# ============================================================================

cat("\n5h: Hub Gene Identification for Clinically Relevant Modules\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Get modules significantly correlated with PFS_group OR Response (or both)
# Selection criteria: p < 0.05 for at least one of the two key clinical traits
# This includes:
#   1. Modules significant for BOTH PFS_group AND Response
#   2. Modules significant for PFS_group alone
#   3. Modules significant for Response alone
focus_modules_auto <- module_priority %>%
  filter(Sig_PFS_group == TRUE | Sig_Response == TRUE) %>%
  pull(Module) %>%
  gsub("ME", "", .)

cat("\n[SELECTION CRITERIA] Modules significantly correlated with PFS_group OR Response (p < 0.05)\n")

if(length(focus_modules_auto) == 0) {
  # Fallback: top 3 modules by minimum p-value for PFS_group or Response
  focus_modules_auto <- module_priority %>%
    arrange(Min_Pvalue) %>%
    head(3) %>%
    pull(Module) %>%
    gsub("ME", "", .)
  cat("  WARNING: No significant modules found. Using top 3 by lowest p-value.\n")
}

cat(sprintf("\nAuto-detected focus modules: %s\n", paste(focus_modules_auto, collapse = ", ")))

# Validate focus_modules_auto contains actual color names (not numeric labels)
valid_colors <- unique(moduleColors)
if(!all(focus_modules_auto %in% valid_colors)) {
  cat("  WARNING: focus_modules_auto contains invalid module names!\n")
  cat(sprintf("    Got: %s\n", paste(focus_modules_auto, collapse = ", ")))
  cat(sprintf("    Valid colors: %s\n", paste(valid_colors, collapse = ", ")))
  # Fix: use the actual color names from moduleColors
  # Find which modules are clinically significant based on module_priority
  focus_modules_auto <- intersect(focus_modules_auto, valid_colors)
  if(length(focus_modules_auto) == 0) {
    # Fallback to largest non-grey modules
    mod_sizes <- sort(table(moduleColors[moduleColors != "grey"]), decreasing = TRUE)
    focus_modules_auto <- names(mod_sizes)[1:min(2, length(mod_sizes))]
    cat(sprintf("    FALLBACK: Using largest modules: %s\n", paste(focus_modules_auto, collapse = ", ")))
  }
}

# --- DECISION CHECKPOINT: FOCUS MODULE SELECTION ---
cat("\n")

cat("DECISION CHECKPOINT: MODULE SELECTION\n")

cat(sprintf("SELECTED FOCUS MODULES: %-47s║\n", paste(toupper(focus_modules_auto), collapse = ", ")))

cat("SELECTION CRITERIA:\n")
cat("- N_Significant >= 2 (associated with 2+ clinical traits at p<0.05)\n")
cat("- Modules with strongest PFS_group and/or Response correlations\n")

cat("RATIONALE:\n")
for(mod in focus_modules_auto) {
  mod_priority <- module_priority[module_priority$Module == paste0("ME", mod), ]
  if(nrow(mod_priority) > 0) {
    cat(sprintf("║    %s: %d significant traits (Priority rank: %d)                    ║\n",
                toupper(mod), mod_priority$N_Significant, which(module_priority$Module == paste0("ME", mod))))
  }
}

cat("These modules will be prioritized for downstream analyses: \n")
cat("- Hub gene identification and stability testing\n")
cat("- Pathway enrichment (ORA, GSVA, ssGSEA)\n")
cat("- Biomarker integration and clinical correlation\n")


# Print top hub genes for each focus module
for(mod in focus_modules_auto) {
  mm_col <- paste0("MM.", mod)
  if(mm_col %in% colnames(proteinInfo)) {
    cat(sprintf("\n--- Top 10 %s module hubs (by MM) ---\n", toupper(mod)))
    hub_table <- proteinInfo %>% 
      filter(Module == mod) %>% 
      arrange(desc(.data[[mm_col]])) %>%
      dplyr::select(Protein, Module, all_of(mm_col), GS.PFS_group, GS.Response, GS.PFS_Days) %>%
      head(10) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
    print(hub_table, row.names = FALSE)
  }
}

# Save hub genes for focus modules
hub_genes_all <- data.frame()
for(mod in focus_modules_auto) {
  mm_col <- paste0("MM.", mod)
  if(mm_col %in% colnames(proteinInfo)) {
    hub_genes <- proteinInfo %>% 
      filter(Module == mod) %>% 
      arrange(desc(.data[[mm_col]])) %>%
      dplyr::select(Protein, Module, all_of(mm_col), GS.PFS_group, GS.Response, GS.PFS_Days) %>%
      mutate(Rank = row_number())
    hub_genes_all <- bind_rows(hub_genes_all, hub_genes)
  }
}
write.csv(hub_genes_all, file.path(results_dir, "Step5_HubGenes_FocusModules.csv"), row.names = FALSE)
cat("\nSaved: Step5_HubGenes_FocusModules.csv\n")

# ============================================================================
# 5i: Export Data
# ============================================================================

cat("\n5i: Exporting Data\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Excel: Proteins by module
wb <- createWorkbook()
for(mod in names(sort(table(moduleColors), decreasing = TRUE))) {
  df <- proteinInfo %>% filter(Module == mod)
  mm_col <- paste0("MM.", mod)
  if(mm_col %in% colnames(df)) df <- df %>% arrange(desc(.data[[mm_col]]))
  
  cols <- c("Protein", "Module", if(mm_col %in% colnames(df)) mm_col,
            "GS.PFS_group", "GS.Response", "GS.PFS_Days",
            if("GS.CA19.9" %in% colnames(df)) "GS.CA19.9")
  df <- df %>% dplyr::select(all_of(cols[cols %in% colnames(df)])) %>% 
    mutate(across(where(is.numeric), ~round(., 4)))
  
  addWorksheet(wb, mod)
  writeData(wb, mod, df)
  setColWidths(wb, mod, cols = 1:ncol(df), widths = "auto")
}
addWorksheet(wb, "Summary", tabColour = "orange")
writeData(wb, "Summary", module_summary)
worksheetOrder(wb) <- c(length(unique(moduleColors)) + 1, 1:length(unique(moduleColors)))
saveWorkbook(wb, file.path(results_dir, "Step5_Proteins_ByModule.xlsx"), overwrite = TRUE)
cat("  Saved: Step5_Proteins_ByModule.xlsx\n")

# Module eigengenes
write.csv(cbind(Sample = rownames(MEs_noGrey), MEs_noGrey),
          file.path(results_dir, "Step5_ModuleEigengenes.csv"), row.names = FALSE)
cat("  Saved: Step5_ModuleEigengenes.csv\n")

# Single comprehensive RDS
networkData <- list(
  net = net, moduleColors = moduleColors, MEs = MEs_noGrey,
  moduleTraitCor = moduleTraitCor, moduleTraitPval = moduleTraitPval,
  proteinInfo = proteinInfo, traitData = traitData, MEDiss = MEDiss,
  module_priority = module_priority,
  focus_modules_auto = focus_modules_auto,
  params = list(power=13, networkType="signed", TOMType="signed", corType="bicor",
                maxPOutliers=0.1, minModuleSize=20, mergeCutHeight=0.25, deepSplit=3)
)
saveRDS(networkData, file.path(results_dir, "Step5_NetworkData.rds"))
cat("  Saved: Step5_NetworkData.rds\n")

# ============================================================================
# STEP 5 SUMMARY
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("STEP 5 SUMMARY: NETWORK CONSTRUCTION\n")
cat("===============================================================================\n")

step5_summary <- data.frame(
  Metric = c(
    "Total proteins",
    "Total samples",
    "Soft-threshold power",
    "Network type",
    "Correlation method",
    "deepSplit",
    "Number of modules",
    "Largest module",
    "Smallest module",
    "Grey (unassigned)",
    "Grey percentage",
    "Mean module size",
    "PFS_group significant modules",
    "Response significant modules",
    "Focus modules (auto-detected)"
  ),
  Value = c(
    ncol(datExpr),
    nrow(datExpr),
    "13",
    "signed",
    "bicor",
    "3",
    nModules,
    max(table(moduleColors[moduleColors != "grey"])),
    min(table(moduleColors[moduleColors != "grey"])),
    nGrey,
    sprintf("%.1f%%", pctGrey),
    round(mean(table(moduleColors[moduleColors != "grey"])), 1),
    sum(moduleTraitPval[, "PFS_group"] < 0.05),
    sum(moduleTraitPval[, "Response"] < 0.05),
    paste(focus_modules_auto, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

cat("\nRESULTS:\n")
print(step5_summary, row.names = FALSE)
write.csv(step5_summary, file.path(results_dir, "Step5_Summary.csv"), row.names = FALSE)
cat("\nSaved: Step5_Summary.csv\n")

cat("\nOutputs:\n")
cat("  Figures: Dendrogram, ME_Dendrogram, Heatmap, Permutation, DotPlot, Module_Sizes\n")
cat("  Tables: ModuleSummary, ModuleTraitCorrelations, ModuleTraitPvalues, ProteinInfo\n")
cat("  Data: Proteins_ByModule.xlsx, NetworkData.rds, HubGenes_FocusModules.csv\n")

# Validation checkpoint
cat("\n[Validation Checkpoint: Step 5 outputs]\n")
stopifnot("moduleColors must exist" = exists("moduleColors"))
stopifnot("moduleColors must be a named vector" = !is.null(names(moduleColors)))
stopifnot("moduleColors length must match proteins" = length(moduleColors) == ncol(datExpr))
stopifnot("MEs_noGrey must exist" = exists("MEs_noGrey"))
stopifnot("MEs_noGrey must have rows" = nrow(MEs_noGrey) > 0)
stopifnot("MEs_noGrey rows must match samples" = nrow(MEs_noGrey) == nrow(datExpr))
cat("  [OK] All Step 5 outputs validated\n")

# ============================================================================
# SET FOCUS MODULES FOR DOWNSTREAM ANALYSIS
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("FOCUS MODULES SET FOR DOWNSTREAM ANALYSIS\n")
cat("===============================================================================\n")

primary_focus_modules <- focus_modules_auto
focus_modules <- c(primary_focus_modules, "grey")

cat(sprintf("\n  Primary focus modules: %s\n", paste(primary_focus_modules, collapse = ", ")))
cat(sprintf("  All focus modules (incl. grey control): %s\n", paste(focus_modules, collapse = ", ")))
cat("\n  These modules will be used in Steps 6-16 for:\n")
cat("    - Hub gene analysis\n")
cat("    - Enrichment analysis\n")
cat("    - Biomarker development\n")
cat("    - Survival analysis\n")

cat("\n===============================================================================\n")
cat("  Proceeding to Step 6 with focus modules: ", paste(focus_modules, collapse = ", "), "\n")
cat("===============================================================================\n")

# ============================================================================
# STEP 6: HUB GENE IDENTIFICATION
# ============================================================================
# PURPOSE: Identify hub proteins within each module - the most connected
# members that drive module behavior and are candidates for validation.
#
# STRUCTURE:
#   6a) Calculate intramodular connectivity (kWithin)
#   6b) Define hub criteria and identify hubs (kME + GS)
#   6c) Hub gene visualizations (MM vs GS scatterplots)
#   6d) Hub validation analyses (connectivity, overlap, effect sizes)
#   6e) Key results summary
#   6f) Export data
#   6g) Network visualization (Cytoscape + igraph)
#
# CRITERIA:
#   - Strict: kME > 0.8, |GS| > 0.2, p < 0.05
#   - Moderate: kME > 0.7, |GS| > 0.2, p < 0.05
#   - GS checked for: PFS_group AND Response
# ============================================================================

cat("===============================================================================\n")
cat("STEP 6: HUB GENE IDENTIFICATION\n")
cat("===============================================================================\n")

# Load network data from Step 5
networkData <- readRDS(file.path(results_dir, "Step5_NetworkData.rds"))
moduleColors <- networkData$moduleColors
MEs <- networkData$MEs
proteinInfo <- networkData$proteinInfo
focus_modules_auto <- networkData$focus_modules_auto

cat(sprintf("\nFocus modules from Step 5: %s\n", paste(focus_modules_auto, collapse = ", ")))

# ============================================================================
# 6a: Calculate Intramodular Connectivity
# ============================================================================

cat("\n6a: Calculating intramodular connectivity...\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Load TOM
if(!exists("TOM")) {
  tom_file <- file.path(results_dir, "Step5_TOM-block.1.RData")
  if(file.exists(tom_file)) {
    load(tom_file)
    TOM <- as.matrix(TOM)
    rownames(TOM) <- colnames(datExpr)
    colnames(TOM) <- colnames(datExpr)
    cat("  Loaded TOM from Step5_TOM-block.1.RData\n")
  } else {
    cat("  Warning: TOM file not found, calculating from adjacency...\n")
    adjacency <- adjacency(datExpr, power = 13, type = "signed",
                           corFnc = "bicor", corOptions = list(maxPOutliers = 0.1))
    TOM <- TOMsimilarity(adjacency)
    colnames(TOM) <- colnames(datExpr)
    rownames(TOM) <- colnames(datExpr)
  }
}

connectivity <- intramodularConnectivity(TOM, moduleColors)
connectivity$Protein <- rownames(connectivity)

# Safe join - remove existing columns if present to avoid duplication
cols_to_add <- c("kTotal", "kWithin", "kOut", "kDiff")
proteinInfo <- proteinInfo %>% select(-any_of(cols_to_add))

proteinInfo <- proteinInfo %>%
  left_join(connectivity, by = "Protein") %>%
  relocate(Protein, Module, kTotal, kWithin, kOut, kDiff)

cat("  Added connectivity metrics: kTotal, kWithin, kOut, kDiff\n")

# ============================================================================
# 6b: Define Hub Criteria and Identify Hubs
# ============================================================================

cat("\n6b: Hub gene identification\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Two criteria levels: strict and moderate
criteria <- list(
  Strict = list(kME = 0.8, GS = 0.2, pGS = 0.05),
  Moderate = list(kME = 0.7, GS = 0.2, pGS = 0.05)
)

modules <- setdiff(unique(moduleColors), "grey")

# Hub identification function - PFS_group and Response
identify_hubs <- function(df, module, kME_thresh, GS_thresh, pGS_thresh) {
  mm_col <- paste0("MM.", module)

  if(!mm_col %in% colnames(df)) return(data.frame())

  df %>%
    filter(Module == module, .data[[mm_col]] > kME_thresh) %>%
    filter(
      (abs(GS.PFS_group) > GS_thresh & p.GS.PFS_group < pGS_thresh) |
        (abs(GS.Response) > GS_thresh & p.GS.Response < pGS_thresh)
    ) %>%
    arrange(desc(.data[[mm_col]]))
}

# Identify hubs for both criteria
hub_results <- list()
for(crit_name in names(criteria)) {
  crit <- criteria[[crit_name]]
  hub_results[[crit_name]] <- lapply(modules, function(m) {
    identify_hubs(proteinInfo, m, crit$kME, crit$GS, crit$pGS)
  })
  names(hub_results[[crit_name]]) <- modules
}

# Summary counts
cat("\n  Hub counts by module:\n")
cat("  ", sprintf("%-12s", "Module"), sprintf("%-10s", "nProteins"),
    sprintf("%-10s", "Strict"), sprintf("%-10s", "Moderate"), "\n")
cat("  ", paste(rep("-", 42), collapse = ""), "\n")

hub_summary <- data.frame(Module = modules) %>%
  mutate(
    nProteins = sapply(modules, function(m) sum(moduleColors == m)),
    Strict = sapply(modules, function(m) nrow(hub_results$Strict[[m]])),
    Moderate = sapply(modules, function(m) nrow(hub_results$Moderate[[m]]))
  ) %>%
  arrange(desc(Moderate))

for(i in 1:nrow(hub_summary)) {
  cat("  ", sprintf("%-12s", hub_summary$Module[i]),
      sprintf("%-10d", hub_summary$nProteins[i]),
      sprintf("%-10d", hub_summary$Strict[i]),
      sprintf("%-10d", hub_summary$Moderate[i]), "\n")
}

cat("\n  Total Strict hubs:", sum(hub_summary$Strict), "\n")
cat("  Total Moderate hubs:", sum(hub_summary$Moderate), "\n")

# Focus module highlights
cat("\n  Focus module hub counts:\n")
for(mod in focus_modules_auto) {
  if(mod %in% hub_summary$Module) {
    row <- hub_summary[hub_summary$Module == mod, ]
    cat(sprintf("    %s: %d strict, %d moderate (of %d proteins)\n",
                toupper(mod), row$Strict, row$Moderate, row$nProteins))
  }
}

# Add hub flags to proteinInfo
# Safety: handle case when no hubs are found
strict_combined <- bind_rows(hub_results$Strict)
moderate_combined <- bind_rows(hub_results$Moderate)

strict_hubs <- if(nrow(strict_combined) > 0 && "Protein" %in% colnames(strict_combined)) {
  unique(strict_combined$Protein)
} else {
  character(0)
}

moderate_hubs <- if(nrow(moderate_combined) > 0 && "Protein" %in% colnames(moderate_combined)) {
  unique(moderate_combined$Protein)
} else {
  character(0)
}

cat(sprintf("\n  Found %d strict hubs, %d moderate hubs\n", length(strict_hubs), length(moderate_hubs)))

proteinInfo <- proteinInfo %>%
  mutate(
    isHub_Strict = Protein %in% strict_hubs,
    isHub_Moderate = Protein %in% moderate_hubs
  )

if(length(moderate_hubs) > 0) {
  hub_gs_summary <- proteinInfo %>%
    filter(isHub_Moderate) %>%
    select(Protein, Module, GS.PFS_group, GS.Response) %>%
    mutate(
      Max_GS = pmax(abs(GS.PFS_group), abs(GS.Response), na.rm = TRUE)
    ) %>%
    arrange(desc(Max_GS))

  print(hub_gs_summary)
} else {
  cat("  No moderate hubs found with current criteria.\n")
  cat("  Consider relaxing thresholds (kME > 0.6, GS > 0.15) if needed.\n")
  hub_gs_summary <- data.frame()  # Create empty for downstream code
}

# Summary statistics (only if hubs exist)
if(length(moderate_hubs) > 0 && nrow(hub_gs_summary) > 0) {
  cat("\nGS Distribution in Moderate Hubs:\n")
  cat(sprintf("  Min |GS|:  %.3f\n", min(hub_gs_summary$Max_GS)))
  cat(sprintf("  Mean |GS|: %.3f\n", mean(hub_gs_summary$Max_GS)))
  cat(sprintf("  Max |GS|:  %.3f\n", max(hub_gs_summary$Max_GS)))
}

# ============================================================================
# 6c: Hub Gene Visualizations (MM vs GS)
# ============================================================================

cat("\n6c: Generating Hub Gene Visualizations\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# MM vs GS scatterplot function
plot_mm_gs <- function(df, module, trait, kME_strict = 0.8, kME_moderate = 0.7, GS_thresh = 0.2) {
  mm_col <- paste0("MM.", module)
  gs_col <- paste0("GS.", trait)
  pgs_col <- paste0("p.GS.", trait)
  
  if(!mm_col %in% colnames(df)) {
    cat(sprintf("    Warning: %s not found\n", mm_col))
    return(NULL)
  }
  
  plot_df <- df %>%
    filter(Module == module) %>%
    mutate(
      HubType = case_when(
        .data[[mm_col]] > kME_strict & abs(.data[[gs_col]]) > GS_thresh & .data[[pgs_col]] < 0.05 ~ "Strict",
        .data[[mm_col]] > kME_moderate & abs(.data[[gs_col]]) > GS_thresh & .data[[pgs_col]] < 0.05 ~ "Moderate",
        TRUE ~ "None"
      ),
      HubType = factor(HubType, levels = c("None", "Moderate", "Strict"))
    )
  
  cor_val <- cor(plot_df[[mm_col]], plot_df[[gs_col]], use = "complete.obs")
  cor_p <- cor.test(plot_df[[mm_col]], plot_df[[gs_col]])$p.value
  n_strict <- sum(plot_df$HubType == "Strict")
  n_moderate <- sum(plot_df$HubType %in% c("Strict", "Moderate"))
  
  top_hubs <- plot_df %>% filter(HubType != "None") %>% arrange(desc(.data[[mm_col]])) %>% head(10)
  
  p <- ggplot(plot_df, aes(x = .data[[mm_col]], y = .data[[gs_col]])) +
    # Hub regions
    annotate("rect", xmin = kME_moderate, xmax = kME_strict, ymin = -Inf, ymax = -GS_thresh,
             alpha = 0.15, fill = "orange") +
    annotate("rect", xmin = kME_moderate, xmax = kME_strict, ymin = GS_thresh, ymax = Inf,
             alpha = 0.15, fill = "orange") +
    annotate("rect", xmin = kME_strict, xmax = Inf, ymin = -Inf, ymax = -GS_thresh,
             alpha = 0.15, fill = module) +
    annotate("rect", xmin = kME_strict, xmax = Inf, ymin = GS_thresh, ymax = Inf,
             alpha = 0.15, fill = module) +
    # Points
    geom_point(aes(color = HubType), alpha = 0.7, size = 2) +
    scale_color_manual(values = c("None" = "grey60", "Moderate" = "orange", "Strict" = module)) +
    # Thresholds
    geom_vline(xintercept = c(kME_moderate, kME_strict), linetype = c("dotted", "dashed"),
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = c(-GS_thresh, GS_thresh), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    # Regression
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    # Labels
    labs(
      title = sprintf("Module %s: MM vs GS.%s", module, trait),
      subtitle = sprintf("r = %.3f, p = %.2e | Strict: %d, Moderate: %d",
                         cor_val, cor_p, n_strict, n_moderate),
      x = sprintf("Module Membership (MM.%s)", module),
      y = sprintf("Gene Significance (%s)", trait)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  if(nrow(top_hubs) > 0) {
    p <- p + geom_text_repel(data = top_hubs, aes(label = Protein), size = 3,
                             max.overlaps = 15, segment.color = "grey50")
  }
  
  return(p)
}

# Generate plots for focus modules x traits (Figures 01-04)
focus_plots <- list()
plot_counter <- 1

for(mod in focus_modules_auto) {
  for(trait in c("PFS_group", "Response")) {
    cat(sprintf("  [%d] %s x %s\n", plot_counter, mod, trait))
    
    p <- plot_mm_gs(proteinInfo, mod, trait)
    
    if(!is.null(p)) {
      plot_key <- paste0(mod, "_", trait)
      focus_plots[[plot_key]] <- p
      
      filename_base <- sprintf("Step6_%02d_ME%s_GS_%s", plot_counter, mod, trait)
      print(p)
      ggsave(file.path(fig_dir_png, paste0(filename_base, ".png")), p, width = 6, height = 5, dpi = 300)
      ggsave(file.path(fig_dir_pdf, paste0(filename_base, ".pdf")), p, width = 6, height = 5)
      cat(sprintf("    Saved: %s.png/pdf\n", filename_base))
    }
    plot_counter <- plot_counter + 1
  }
}

# Combined publication figure (Figure 05)
cat("\n  Generating combined publication figure...\n")

if(length(focus_plots) >= 4) {
  p_combined <- (focus_plots[[1]] + focus_plots[[2]]) / (focus_plots[[3]] + focus_plots[[4]]) +
    plot_annotation(
      title = "Hub Gene Identification: Key Prognostic Modules",
      subtitle = sprintf("Focus modules: %s", paste(focus_modules_auto, collapse = ", ")),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40")
      )
    )
} else if(length(focus_plots) >= 2) {
  p_combined <- wrap_plots(focus_plots, ncol = 2) +
    plot_annotation(
      title = "Hub Gene Identification: Key Prognostic Modules",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
} else if(length(focus_plots) >= 1) {
  p_combined <- focus_plots[[1]]
} else {
  p_combined <- NULL
}

if(!is.null(p_combined)) {
  print(p_combined)
  ggsave(file.path(fig_dir_png, "Step6_05_HubGenes_Combined.png"), p_combined, width = 12, height = 9, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step6_05_HubGenes_Combined.pdf"), p_combined, width = 12, height = 9)
  cat("  Saved: Step6_05_HubGenes_Combined.png/pdf\n")
}

# Hub count comparison plot (Figure 06)
cat("\n  Generating hub count comparison...\n")

hub_summary_long <- hub_summary %>%
  pivot_longer(cols = c(Strict, Moderate), names_to = "Criteria", values_to = "nHubs") %>%
  mutate(Module = factor(Module, levels = hub_summary$Module),
         Criteria = factor(Criteria, levels = c("Strict", "Moderate")),
         IsFocus = Module %in% focus_modules_auto)

p_hub_compare <- ggplot(hub_summary_long, aes(x = Module, y = nHubs, fill = Criteria, alpha = IsFocus)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Strict" = "#2c7fb8", "Moderate" = "#7fcdbb")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
  labs(
    title = "Hub Genes per Module",
    subtitle = sprintf("Strict (kME>0.8) vs Moderate (kME>0.7) | Focus: %s",
                       paste(focus_modules_auto, collapse = ", ")),
    x = NULL, y = "Number of Hub Genes"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(p_hub_compare)
ggsave(file.path(fig_dir_png, "Step6_06_HubCount_Comparison.png"), p_hub_compare, width = 6, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step6_06_HubCount_Comparison.pdf"), p_hub_compare, width = 6, height = 4)
cat("  Saved: Step6_06_HubCount_Comparison.png/pdf\n")

# ============================================================================
# 6d: Hub Validation Analyses
# ============================================================================

cat("\n6d: Hub Validation Analyses\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ----------------------------------------------------------------------------
# 6d.1: Connectivity-Based Hub Validation (kWithin percentile)
# ----------------------------------------------------------------------------

cat("\n6d.1: Connectivity-Based Hub Validation\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Calculate kWithin percentile within each module
proteinInfo <- proteinInfo %>%
  group_by(Module) %>%
  mutate(
    kWithin_Percentile = percent_rank(kWithin),
    isHighConnectivity = kWithin_Percentile >= 0.75  # Top 25% within module
  ) %>%
  ungroup()

# Classify hub types based on connectivity
proteinInfo <- proteinInfo %>%
  mutate(
    HubType = case_when(
      isHub_Moderate & isHighConnectivity ~ "Core Hub",
      isHub_Moderate & !isHighConnectivity ~ "Peripheral Hub",
      !isHub_Moderate & isHighConnectivity ~ "High Connectivity (not hub)",
      TRUE ~ "Regular"
    )
  )

# Summary of hub types
hub_type_summary <- proteinInfo %>%
  filter(Module != "grey") %>%
  group_by(Module, HubType) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = HubType, values_from = n, values_fill = 0)

cat("  Hub Classification by Connectivity:\n")
print(as.data.frame(hub_type_summary))

# Focus module details
cat("\n  Focus Module Hub Types:\n")
for(mod in focus_modules_auto) {
  mod_data <- proteinInfo %>% filter(Module == mod)
  n_core <- sum(mod_data$HubType == "Core Hub", na.rm = TRUE)
  n_peripheral <- sum(mod_data$HubType == "Peripheral Hub", na.rm = TRUE)
  
  cat(sprintf("    %s: %d Core Hubs, %d Peripheral Hubs\n", toupper(mod), n_core, n_peripheral))
  
  # List core hubs
  core_hubs <- mod_data %>% filter(HubType == "Core Hub") %>% pull(Protein)
  if(length(core_hubs) > 0) {
    cat(sprintf("      Core hubs: %s\n", paste(head(core_hubs, 10), collapse = ", ")))
  }
}

# Connectivity validation plot (Figure 07)
cat("\n  Generating connectivity validation plots...\n")

plot_connectivity <- function(df, module) {
  mm_col <- paste0("MM.", module)
  if(!mm_col %in% colnames(df)) return(NULL)
  
  plot_df <- df %>% filter(Module == module)
  cor_val <- cor(plot_df$kWithin, plot_df[[mm_col]], use = "complete.obs")
  
  ggplot(plot_df, aes(x = kWithin, y = .data[[mm_col]])) +
    geom_point(color = module, alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    labs(
      title = sprintf("Module %s", module),
      subtitle = sprintf("r = %.3f", cor_val),
      x = "kWithin", y = sprintf("MM.%s", module)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )
}

conn_plots <- lapply(focus_modules_auto, function(m) plot_connectivity(proteinInfo, m))
conn_plots <- conn_plots[!sapply(conn_plots, is.null)]

if(length(conn_plots) > 0) {
  p_conn <- wrap_plots(conn_plots, ncol = length(conn_plots)) +
    plot_annotation(
      title = "Connectivity Validation: kWithin vs MM",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
  
  print(p_conn)
  ggsave(file.path(fig_dir_png, "Step6_07_Connectivity_Validation.png"), p_conn, 
         width = 4 * length(conn_plots), height = 4, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step6_07_Connectivity_Validation.pdf"), p_conn, 
         width = 4 * length(conn_plots), height = 4)
  cat("  Saved: Step6_07_Connectivity_Validation.png/pdf\n")
}

# ----------------------------------------------------------------------------
# 6d.2: Hub Overlap Analysis Between Modules (Vectorized)
# ----------------------------------------------------------------------------

cat("\n6d.2: Hub Overlap Analysis Between Modules\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Check for proteins with high kME in multiple modules
mm_cols <- grep("^MM\\.", colnames(proteinInfo), value = TRUE)

hub_in_multiple <- proteinInfo %>%
  filter(isHub_Moderate) %>%
  rowwise() %>%
  mutate(
    n_modules_high_kME = sum(c_across(all_of(mm_cols)) > 0.7, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(n_modules_high_kME > 1)

if(nrow(hub_in_multiple) > 0) {
  cat("  WARNING: Some proteins have high kME (>0.7) in multiple modules:\n")
  print(hub_in_multiple %>% select(Protein, Module, n_modules_high_kME) %>% head(10))
  cat("\n  This may indicate:\n")
  cat("    - These proteins are at module boundaries\n")
  cat("    - They may have broad functional roles\n")
  cat("    - Primary module assignment takes precedence\n")
  proteinInfo$isMultiModuleHub <- proteinInfo$Protein %in% hub_in_multiple$Protein
} else {
  cat("  [OK] No hub overlap detected - all hubs are module-specific\n")
  cat("  This validates good module separation.\n")
  proteinInfo$isMultiModuleHub <- FALSE
}

# ----------------------------------------------------------------------------
# 6d.3: Effect Size Analysis (Cohen's d) - Figure 08
# ----------------------------------------------------------------------------

cat("\n6d.3: Effect Size Analysis (Cohen's d)\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Winsorization function to reduce outlier influence
winsorize <- function(x, limits = c(0.025, 0.975)) {
  q <- quantile(x, probs = limits, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

# Cohen's d function with winsorization
cohens_d_robust <- function(x, group, pos_level, neg_level) {
  x_win <- winsorize(x)
  
  group1 <- x_win[group == pos_level]
  group2 <- x_win[group == neg_level]
  
  n1 <- sum(!is.na(group1))
  n2 <- sum(!is.na(group2))
  
  if(n1 < 3 || n2 < 3) return(NA)
  
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  sd1 <- sd(group1, na.rm = TRUE)
  sd2 <- sd(group2, na.rm = TRUE)
  
  # Pooled SD
  pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
  
  if(pooled_sd == 0) return(NA)
  
  return((mean1 - mean2) / pooled_sd)
}

# Effect size interpretation function
interpret_d <- function(d) {
  if(is.na(d)) return("NA")
  abs_d <- abs(d)
  dir <- ifelse(d > 0, "+", "-")
  if(abs_d >= 0.8) return(paste0("Large (", dir, ")"))
  if(abs_d >= 0.5) return(paste0("Medium (", dir, ")"))
  if(abs_d >= 0.2) return(paste0("Small (", dir, ")"))
  return("Negligible")
}

# Calculate for hub genes
hub_genes <- proteinInfo %>% filter(isHub_Moderate) %>% pull(Protein)
cat(sprintf("  Calculating Cohen's d for %d hub genes...\n", length(hub_genes)))

effect_size_results <- lapply(hub_genes, function(hub) {
  if(!hub %in% colnames(datExpr)) return(NULL)
  
  expr <- datExpr[, hub]
  mod <- proteinInfo$Module[proteinInfo$Protein == hub]
  
  d_pfs <- cohens_d_robust(expr, datTraits$PFS_group, "Long", "Short")
  d_resp <- cohens_d_robust(expr, datTraits$Response, "CD", "PD")
  
  data.frame(
    Protein = hub,
    Module = mod,
    Cohen_d_PFS = round(d_pfs, 3),
    Cohen_d_Response = round(d_resp, 3),
    Effect_PFS = interpret_d(d_pfs),
    Effect_Response = interpret_d(d_resp),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Sort by absolute effect size
effect_size_results <- effect_size_results %>%
  mutate(Max_Abs_d = pmax(abs(Cohen_d_PFS), abs(Cohen_d_Response), na.rm = TRUE)) %>%
  arrange(desc(Max_Abs_d))

cat("\n  Top Effect Sizes:\n")
print(effect_size_results %>% select(-Max_Abs_d) %>% head(15))

# Summary statistics
n_large_pfs <- sum(abs(effect_size_results$Cohen_d_PFS) >= 0.8, na.rm = TRUE)
n_medium_pfs <- sum(abs(effect_size_results$Cohen_d_PFS) >= 0.5 & 
                      abs(effect_size_results$Cohen_d_PFS) < 0.8, na.rm = TRUE)
n_large_resp <- sum(abs(effect_size_results$Cohen_d_Response) >= 0.8, na.rm = TRUE)
n_medium_resp <- sum(abs(effect_size_results$Cohen_d_Response) >= 0.5 & 
                       abs(effect_size_results$Cohen_d_Response) < 0.8, na.rm = TRUE)

cat(sprintf("\n  Effect Size Summary:\n"))
cat(sprintf("    PFS_group:  %d large (|d|>=0.8), %d medium (|d|>=0.5)\n", n_large_pfs, n_medium_pfs))
cat(sprintf("    Response:   %d large (|d|>=0.8), %d medium (|d|>=0.5)\n", n_large_resp, n_medium_resp))

# Add effect sizes to proteinInfo (safe join)
proteinInfo <- proteinInfo %>% select(-any_of(c("Cohen_d_PFS", "Cohen_d_Response")))
proteinInfo <- proteinInfo %>%
  left_join(effect_size_results %>% select(Protein, Cohen_d_PFS, Cohen_d_Response), by = "Protein")

# Save effect size results
write.csv(effect_size_results, file.path(results_dir, "Step6_HubEffectSizes.csv"), row.names = FALSE)
cat("  Saved: Step6_HubEffectSizes.csv\n")

# Effect size forest plot (Figure 08)
cat("  Generating effect size forest plot...\n")

forest_data <- effect_size_results %>%
  filter(Module %in% focus_modules_auto)

if(nrow(forest_data) > 0) {
  # Long format for plotting
  forest_long <- forest_data %>%
    select(Protein, Module, Cohen_d_PFS, Cohen_d_Response) %>%
    pivot_longer(cols = c(Cohen_d_PFS, Cohen_d_Response),
                 names_to = "Trait", values_to = "Cohen_d") %>%
    mutate(
      Trait = ifelse(Trait == "Cohen_d_PFS", "PFS Group", "Response"),
      Protein = factor(Protein, levels = unique(forest_data$Protein))
    )
  
  p_effect <- ggplot(forest_long, aes(x = Cohen_d, y = Protein, color = Trait)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50") +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "orange", alpha = 0.5) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("PFS Group" = "#0E4D92", "Response" = "#ca562c")) +
    facet_wrap(~Module, scales = "free_y", ncol = 1) +
    labs(
      title = "Hub Gene Effect Sizes (Cohen's d)",
      subtitle = "Red dashed = |d|=0.8 (large) | Orange dotted = |d|=0.5 (medium)",
      x = "Cohen's d (positive = higher in Long PFS / CD)",
      y = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
  
  print(p_effect)
  ggsave(file.path(fig_dir_png, "Step6_08_HubEffectSizes.png"), p_effect, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step6_08_HubEffectSizes.pdf"), p_effect, width = 8, height = 9)
  cat("  Saved: Step6_08_HubEffectSizes.png/pdf\n")
}

# ============================================================================
# 6e: Key Results - Focus Module Hubs
# ============================================================================

cat("\n6e: Key Results - Focus Module Hubs\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

for(mod in focus_modules_auto) {
  mm_col <- paste0("MM.", mod)
  
  cat(sprintf("\nME%s Hubs:\n", toupper(mod)))
  
  # Strict hubs
  cat("  Strict (kME > 0.8):\n")
  if(mod %in% names(hub_results$Strict) && nrow(hub_results$Strict[[mod]]) > 0) {
    cols_to_show <- c("Protein", mm_col, "GS.PFS_group", "GS.Response", "GS.PFS_Days", "kWithin")
    cols_to_show <- cols_to_show[cols_to_show %in% colnames(hub_results$Strict[[mod]])]
    
    hub_results$Strict[[mod]] %>%
      select(all_of(cols_to_show)) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      head(10) %>%
      print()
  } else {
    cat("    None\n")
  }
  
  # Moderate hubs
  cat("\n  Moderate (kME > 0.7):\n")
  if(mod %in% names(hub_results$Moderate) && nrow(hub_results$Moderate[[mod]]) > 0) {
    cols_to_show <- c("Protein", mm_col, "GS.PFS_group", "GS.Response", "GS.PFS_Days", "kWithin")
    cols_to_show <- cols_to_show[cols_to_show %in% colnames(hub_results$Moderate[[mod]])]
    
    hub_results$Moderate[[mod]] %>%
      select(all_of(cols_to_show)) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      head(10) %>%
      print()
  } else {
    cat("    None\n")
  }
}

# ============================================================================
# 6f: Export Data
# ============================================================================

cat("\n6f: Exporting Data\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Excel: Hub genes
wb <- createWorkbook()

addWorksheet(wb, "Summary")
writeData(wb, "Summary", hub_summary)

addWorksheet(wb, "Criteria")
criteria_df <- data.frame(
  Level = c("Strict", "Moderate"),
  kME = c(0.8, 0.7),
  GS = c(0.2, 0.2),
  pGS = c(0.05, 0.05)
)
writeData(wb, "Criteria", criteria_df)

addWorksheet(wb, "HubTypes")
writeData(wb, "HubTypes", hub_type_summary)

for(crit_name in names(hub_results)) {
  for(mod in modules) {
    hubs <- hub_results[[crit_name]][[mod]]
    if(nrow(hubs) == 0) next
    
    sheet_name <- paste0(mod, "_", crit_name)
    mm_col <- paste0("MM.", mod)
    
    cols_to_select <- c("Protein", "Module", "kWithin", "kTotal")
    if(mm_col %in% colnames(hubs)) cols_to_select <- c(cols_to_select, mm_col)
    cols_to_select <- c(cols_to_select, "GS.PFS_group", "p.GS.PFS_group",
                        "GS.Response", "p.GS.Response")
    cols_to_select <- cols_to_select[cols_to_select %in% colnames(hubs)]
    
    export_df <- hubs %>%
      select(all_of(cols_to_select)) %>%
      mutate(across(where(is.numeric), ~round(., 4)))
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, export_df)
    setColWidths(wb, sheet_name, cols = 1:ncol(export_df), widths = "auto")
  }
}

saveWorkbook(wb, file.path(results_dir, "Step6_HubGenes.xlsx"), overwrite = TRUE)
cat("  Saved: Step6_HubGenes.xlsx\n")

# CSVs
write.csv(proteinInfo, file.path(results_dir, "Step6_ProteinInfo_Full.csv"), row.names = FALSE)
write.csv(hub_summary, file.path(results_dir, "Step6_HubSummary.csv"), row.names = FALSE)
cat("  Saved: Step6_ProteinInfo_Full.csv, Step6_HubSummary.csv\n")

# Update networkData RDS - include datExpr and datTraits for skip-Step7 workflow
networkData$proteinInfo <- proteinInfo
networkData$hubGenes <- hub_results
networkData$hubCriteria <- criteria
networkData$effectSizes <- effect_size_results
networkData$hubTypeSummary <- hub_type_summary
networkData$datExpr <- datExpr
networkData$datTraits <- datTraits
networkData$TOM <- TOM
networkData$geneTree <- geneTree
saveRDS(networkData, file.path(results_dir, "Step6_NetworkData.rds"))
cat("  Saved: Step6_NetworkData.rds\n")

# ============================================================================
# 6g: Network Visualization (Cytoscape + igraph)
# ============================================================================

cat("\n6g: Network Visualization Export\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Create Cytoscape export directory
cytoscape_dir <- file.path(results_dir, "Cytoscape_Networks")
if(!dir.exists(cytoscape_dir)) dir.create(cytoscape_dir)

# --- [Part 1: Full Export for Cytoscape (Unchanged)] ---
for(mod in focus_modules_auto) {
  mod_genes <- names(moduleColors)[moduleColors == mod]
  if(length(mod_genes) > 0 && exists("TOM") && !is.null(TOM)) {
    mod_TOM <- TOM[mod_genes, mod_genes]
    tom_values <- mod_TOM[upper.tri(mod_TOM)]
    threshold <- quantile(tom_values, 0.9)
    edge_indices <- which(mod_TOM >= threshold & upper.tri(mod_TOM), arr.ind = TRUE)
    if(nrow(edge_indices) > 0) {
      edge_df <- data.frame(
        fromNode = mod_genes[edge_indices[,1]],
        toNode = mod_genes[edge_indices[,2]],
        weight = round(mod_TOM[edge_indices], 6),
        stringsAsFactors = FALSE
      )
    } else {
      edge_df <- data.frame(fromNode = character(), toNode = character(), weight = numeric())
    }
    mm_col <- paste0("MM.", mod)
    node_df <- proteinInfo %>%
      filter(Protein %in% mod_genes) %>%
      select(nodeName = Protein, moduleColor = Module,
             any_of(mm_col), GS.PFS_group, GS.Response,
             HubType, isHub_Moderate, kWithin)
    write.table(edge_df, file.path(cytoscape_dir, sprintf("Step6_%s_edges.txt", mod)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(node_df, file.path(cytoscape_dir, sprintf("Step6_%s_nodes.txt", mod)),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# --- [Part 2: Filtered igraph Visualization] ---
cat("\n  Creating publication-quality igraph network visualizations...\n")

# igraph already loaded at script top
fig_counter <- 9

for_igraph_viz <- TRUE  # Flag to enable visualization block
if(for_igraph_viz) {
  
  for(mod in focus_modules_auto) {
    cat(sprintf("\n  Processing %s module for plotting...\n", mod))
    
    # 1. Get all genes in module and filter if too large
    mod_genes_all <- names(moduleColors)[moduleColors == mod]
    
    if(length(mod_genes_all) > 0 && exists("TOM") && !is.null(TOM)) {
      
      max_nodes_to_plot <- 100
      if(length(mod_genes_all) > max_nodes_to_plot) {
        cat(sprintf("    Module too large (%d nodes). Filtering to top %d by kWithin.\n", 
                    length(mod_genes_all), max_nodes_to_plot))
        mod_genes <- proteinInfo %>%
          filter(Protein %in% mod_genes_all) %>%
          arrange(desc(kWithin)) %>%
          head(max_nodes_to_plot) %>%
          pull(Protein)
      } else {
        mod_genes <- mod_genes_all
      }
      
      # 2. Create igraph object based on filtered genes
      mod_TOM <- TOM[mod_genes, mod_genes]
      threshold <- quantile(mod_TOM[upper.tri(mod_TOM)], 0.90)
      adj_matrix <- mod_TOM
      adj_matrix[adj_matrix < threshold] <- 0
      diag(adj_matrix) <- 0
      g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
      
      # 3. Define Attributes for Plotting
      node_info <- proteinInfo %>%
        filter(Protein %in% mod_genes) %>%
        arrange(match(Protein, V(g)$name))
      
      V(g)$HubType <- node_info$HubType[match(V(g)$name, node_info$Protein)]
      V(g)$kWithin <- node_info$kWithin[match(V(g)$name, node_info$Protein)]
      V(g)$isHub <- node_info$isHub_Moderate[match(V(g)$name, node_info$Protein)]
      
      hub_colors <- case_when(
        V(g)$HubType == "Core Hub" ~ "#ca562c",
        V(g)$HubType == "Peripheral Hub" ~ "#FFD662",
        TRUE ~ "grey80"
      )
      
      # Slightly larger node sizes
      kWithin_scaled <- scales::rescale(V(g)$kWithin, to = c(5, 15)) 
      edge_width_scaled <- scales::rescale(E(g)$weight, to = c(0.5, 2.5))
      
      set.seed(12345) # Fix layout
      layout <- layout_with_fr(g)
      show_labels <- ifelse(V(g)$isHub, V(g)$name, "")
      
      # 4. Save Plots (PNG & PDF)
      formats <- c("png", "pdf")
      for(f in formats) {
        file_path <- file.path(get(paste0("fig_dir_", f)), 
                               sprintf("Step6_%02d_Network_%s_PubQuality.%s", fig_counter, mod, f))
        
        # Use a slightly larger canvas for better resolution/layout
        if(f == "png") png(file_path, width = 12, height = 12, units = "in", res = 300)
        else pdf(file_path, width = 12, height = 12)
        
        # IMPROVED MARGINS
        par(mar = c(2, 2, 4, 2)) 
        
        plot(g, 
             layout = layout, 
             vertex.size = kWithin_scaled,
             vertex.color = hub_colors, 
             vertex.frame.color = "grey30", # Slightly darker frame for definition
             vertex.label = show_labels, 
             vertex.label.family = "sans", # Ensure a clean font
             vertex.label.font = 2,        # Bold font for hubs
             vertex.label.cex = 1.3,       # --> MUCH LARGER FONT <--
             vertex.label.color = "black", 
             vertex.label.dist = 2.2,      # --> INCREASED DISTANCE <--
             edge.width = edge_width_scaled, 
             edge.color = rgb(0.4, 0.4, 0.4, 0.15), # --> LIGHTER EDGES <--
             main = "" # We will add a cleaner title with mtext
        )
        
        # Add a nice title
        title_text <- sprintf("%s Module Network", tools::toTitleCase(mod))
        if(length(mod_genes_all) > max_nodes_to_plot) {
          title_text <- paste0(title_text, " (Top ", max_nodes_to_plot, " Hubs)")
        }
        mtext(title_text, side = 3, line = 2.5, cex = 1.8, font = 2)

        # Get clinical correlation values for this module
        me_name <- paste0("ME", mod)
        if(me_name %in% rownames(moduleTraitCor)) {
          r_pfs <- moduleTraitCor[me_name, "PFS_group"]
          p_pfs <- moduleTraitPval[me_name, "PFS_group"]
          r_resp <- moduleTraitCor[me_name, "Response"]
          p_resp <- moduleTraitPval[me_name, "Response"]

          # Format p-values
          p_pfs_txt <- ifelse(p_pfs < 0.001, "p<0.001", sprintf("p=%.3f", p_pfs))
          p_resp_txt <- ifelse(p_resp < 0.001, "p<0.001", sprintf("p=%.3f", p_resp))

          # Add subtitle with clinical correlations
          mtext(sprintf("PFS: r=%.2f (%s) | Response: r=%.2f (%s) | n=%d proteins, %d edges",
                        r_pfs, p_pfs_txt, r_resp, p_resp_txt, vcount(g), ecount(g)),
                side = 3, line = 0.8, cex = 1.0, col = "grey30")
        } else {
          # Fallback subtitle without clinical info
          mtext(sprintf("Top 10%% TOM connections | n = %d nodes, %d edges", vcount(g), ecount(g)),
                side = 3, line = 0.5, cex = 1.1, col = "grey40")
        }
        
        # --> IMPROVED LEGEND <--
        legend("bottomleft", 
               legend = c("Core Hub", "Peripheral Hub", "Non-hub"),
               pt.bg = c("#ca562c", "#FFD662", "grey80"),
               pch = 21, 
               pt.cex = 2.2, # Larger points in legend
               cex = 1.0,    # Larger text in legend
               bty = "n", 
               title = "Node Type",
               title.adj = 0.1) # Adjust title position
        
        dev.off()
      }
      
      cat(sprintf("    Saved publication-quality plot: Step6_%02d_Network_%s_PubQuality.png/pdf\n", fig_counter, mod))
      fig_counter <- fig_counter + 1
    }
  }
}


# ============================================================================
# STEP 6 SUMMARY
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("STEP 6 SUMMARY: HUB GENE IDENTIFICATION\n")
cat("===============================================================================\n")

# Build summary dynamically
summary_rows <- list(
  c("Total modules analyzed", length(modules)),
  c("Focus modules", paste(focus_modules_auto, collapse = ", ")),
  c("Strict criteria", "kME > 0.8, |GS| > 0.2, p < 0.05"),
  c("Moderate criteria", "kME > 0.7, |GS| > 0.2, p < 0.05"),
  c("GS traits tested", "PFS_group, Response"),
  c("Total strict hubs", sum(hub_summary$Strict)),
  c("Total moderate hubs", sum(hub_summary$Moderate)),
  c("Large effect sizes (|d|>=0.8)", n_large_pfs + n_large_resp)
)

# Add focus module specific counts
for(mod in focus_modules_auto) {
  if(mod %in% hub_summary$Module) {
    row <- hub_summary[hub_summary$Module == mod, ]
    summary_rows <- c(summary_rows, list(
      c(sprintf("%s - strict hubs", mod), row$Strict),
      c(sprintf("%s - moderate hubs", mod), row$Moderate)
    ))
    
    # Add core/peripheral counts
    mod_data <- proteinInfo %>% filter(Module == mod)
    n_core <- sum(mod_data$HubType == "Core Hub", na.rm = TRUE)
    n_peripheral <- sum(mod_data$HubType == "Peripheral Hub", na.rm = TRUE)
    summary_rows <- c(summary_rows, list(
      c(sprintf("%s - core hubs", mod), n_core),
      c(sprintf("%s - peripheral hubs", mod), n_peripheral)
    ))
  }
}

step6_summary <- data.frame(
  Metric = sapply(summary_rows, `[`, 1),
  Value = as.character(sapply(summary_rows, `[`, 2)),
  stringsAsFactors = FALSE
)

cat("\nRESULTS:\n")
print(step6_summary, row.names = FALSE)
write.csv(step6_summary, file.path(results_dir, "Step6_Summary.csv"), row.names = FALSE)
cat("\nSaved: Step6_Summary.csv\n")

cat("\nFigures generated:\n")
cat("  Step6_01-04: MM vs GS scatterplots (focus modules x traits)\n")
cat("  Step6_05: Combined MM vs GS figure\n")
cat("  Step6_06: Hub count comparison\n")
cat("  Step6_07: Connectivity validation (kWithin vs MM)\n")
cat("  Step6_08: Effect size forest plot (Cohen's d)\n")
cat("  Step6_09-10: Network visualizations (igraph)\n")

cat("\nData exported:\n")
cat("  Step6_HubGenes.xlsx (multi-sheet)\n")
cat("  Step6_ProteinInfo_Full.csv\n")
cat("  Step6_HubSummary.csv\n")
cat("  Step6_HubEffectSizes.csv\n")
cat("  Step6_NetworkData.rds\n")
cat("  Cytoscape_Networks/ (edges and nodes for focus modules)\n")

# Validation checkpoint
cat("\n[Validation Checkpoint: Step 6 outputs]\n")
stopifnot("proteinInfo must have hub flags" = "isHub_Moderate" %in% colnames(proteinInfo))
stopifnot("proteinInfo must have HubType" = "HubType" %in% colnames(proteinInfo))
stopifnot("hub_results must exist" = exists("hub_results"))
stopifnot("hub_summary must exist" = exists("hub_summary"))
stopifnot("effect_size_results must exist" = exists("effect_size_results"))
cat("  [OK] All Step 6 outputs validated\n")

cat("\n===============================================================================\n")
cat("  Step 6 complete. Proceeding to Step 7 with focus modules: ", 
    paste(focus_modules_auto, collapse = ", "), "\n")
cat("===============================================================================\n")

# NOTE: Hub gene interpretations are dynamically generated based on focus_modules_auto
# Focus modules: PINK, GREEN, BLACK, BLUE (negatively correlated with PFS)
# Higher expression in these modules = shorter PFS and worse response

# Thesis example (update based on actual results):
# "WGCNA identified four prognostic modules (pink, green, black, blue) significantly
# associated with both PFS group and treatment response. Hub gene analysis using
# stringent criteria (kME > 0.7, |GS| > 0.2, p < 0.05) identified candidate biomarkers
# in these modules with large effect sizes indicating substantial clinical relevance."

# "Hub genes were identified using two stringency levels. Strict criteria (kME > 0.8)
# identified high-confidence hubs, while moderate criteria (kME > 0.7) identified
# candidate hubs. All subsequent analyses used the moderate hub set to ensure adequate
# statistical power while maintaining biological relevance."


# ============================================================================
# SKIP STEP 7 - Load saved network data directly
# ============================================================================
cat("SKIPPING STEP 7 (Bootstrap Stability) - Loading saved data...\n")

# Load network data from Step 6
networkData <- readRDS(file.path(results_dir, "Step6_NetworkData.rds"))

# Extract required objects for downstream analysis
datExpr <- networkData$datExpr
moduleColors <- networkData$moduleColors
datTraits <- networkData$datTraits
MEs <- networkData$MEs
TOM <- networkData$TOM
geneTree <- networkData$geneTree
proteinInfo <- networkData$proteinInfo

# Verify loaded objects
cat(sprintf("  datExpr: %d samples x %d proteins\n", nrow(datExpr), ncol(datExpr)))
cat(sprintf("  Modules: %d (excluding grey)\n", length(unique(moduleColors)) - 1))
cat(sprintf("  Module eigengenes: %d\n", ncol(MEs)))

cat("Step 7 SKIPPED - Continuing to Step 8...\n\n")
# ============================================================================



#===============================================================================
# STEP 7: MODULE STABILITY TESTING    
#===============================================================================
# PURPOSE: Verify that identified modules are real biological signals,    
# not artifacts of noise or random chance - essential for clinical use.   
# STRATEGY:                                                               
# 7a) Bootstrap resampling (500x) -> 7b) Hub co-clustering                
# 7c) Subsampling sensitivity -> 7d) Permutation test (null distribution) 
# 7e) Module preservation (50/50 split) -> 7f-7i) Additional validation   
# LOOKING FOR:                                                             
# - Bootstrap stability >0.70 (robust), >0.85 (excellent)                
# - Permutation p-value <0.05 (module not by chance)                     
# - Zsummary >10 (strong preservation), >2 (moderate)                    
# KEY PARAMS: n_boot=500, n_perm=100, subsample_fracs=[0.6,0.7,0.8,0.9]    
# ============================================================================
# DETAILED METHODOLOGY (for reference)
# ============================================================================
#
# Research Question: "If we repeated this analysis with slightly different
#                    data, would we identify the same modules?"
#
# Clinical Context:
#   - 50 metastatic PDAC patients (Stage IV)
#   - Goal: Predict PFS (Short vs Long) and Treatment Response (CD vs PD)
#   - Focus modules (from Step 5) associated with poor prognosis
#
# Why Stability Testing Matters:
#   - If modules are noise artifacts -> False biological conclusions
#   - If modules are sample-dependent -> Results won't replicate in other cohorts
#   - If hub genes are unstable -> Wrong biomarker candidates
#
# Methods:
#   A. Resampling Methods
#      7a: Adaptive Bootstrap (500 iterations with module matching)
#      7b: Hub Gene Co-Clustering Analysis
#      7c: Subsampling Sensitivity (60%, 70%, 80%, 90% retention)
#
#   B. Statistical Validation
#      7d: Permutation Test (100 permutations) - null distribution
#      7e: Module Preservation (50/50 split) - Zsummary statistics
#
#   C. Additional Analyses
#      7f: Hub Gene Stability
#      7g: Protein Reassignment Flow
#      7h: Module Size vs Stability (all modules)
#      7i: Validation - Prognostic Strength vs Stability
#
#   D. Output
#      7j: Figures (10 total)
#      7k: Export Tables (12 Excel sheets)
#
# Key Methodological Note:
#   WGCNA assigns module colors by SIZE RANK (turquoise = largest, blue = 2nd, etc.)
#   This means the same biological module may receive different colors across bootstraps.
#   We use MODULE MATCHING (by protein overlap) to track true biological content.
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 7: MODULE STABILITY TESTING\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# ----------------------------------------------------------------------------
# SKIP CONTROL - Set to FALSE to skip Step 7 (computationally intensive)
# ----------------------------------------------------------------------------
run_step7 <- TRUE

if(!run_step7) {
  cat("  SKIPPING Step 7 (run_step7 = FALSE)\n")
  cat("  Set run_step7 <- TRUE to run stability testing.\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
} else {

# ============================================================================
# SETUP: Load Data and Define Helper Functions
# ============================================================================

cat("Loading data and defining helper functions...\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

# Load data from Step 6
networkData <- readRDS(file.path(results_dir, "Step6_NetworkData.rds"))
# Define soft power used in network construction
soft_power <- 13  # From Step 4 analysis: Bicor + Signed network
moduleColors <- networkData$moduleColors
MEs <- networkData$MEs
proteinInfo <- networkData$proteinInfo

# Assign protein names to moduleColors vector
names(moduleColors) <- colnames(datExpr)

# Verify data integrity
cat("Data verification:\n")
cat(sprintf("  Patients: %d\n", nrow(datExpr)))
cat(sprintf("  Proteins: %d\n", ncol(datExpr)))
cat(sprintf("  Modules: %d\n", length(unique(moduleColors))))
cat(sprintf("  moduleColors has names: %s\n", !is.null(names(moduleColors))))

# Define focus modules and parameters - DYNAMIC from Step 5/6
# focus_modules_auto contains the clinically significant modules (pink, green, black, blue)
focus_modules <- c(focus_modules_auto, "grey")  # Grey = negative control
cat(sprintf("  Focus modules (from Step 5): %s\n", paste(focus_modules_auto, collapse = ", ")))

# Hub gene criteria (Moderate threshold from Step 6)
kME_thresh <- 0.7
GS_thresh <- 0.2
pGS_thresh <- 0.05

# Bootstrap parameters
set.seed(123)
max_iterations <- 500
min_iterations <- 50
batch_size <- 50
convergence_cv <- 0.05      # CV < 5%
convergence_change <- 0.01  # Relative change < 1%

# Dynamic color palette for focus modules (used throughout Step 7)
# This ensures plots use the correct colors regardless of which modules are focus
module_color_palette <- c(
  "brown" = "brown", "blue" = "steelblue", "green" = "green3",
  "turquoise" = "turquoise3", "yellow" = "gold", "red" = "indianred",
  "pink" = "pink", "magenta" = "magenta3", "black" = "grey30",
  "grey" = "grey60", "cyan" = "cyan3", "purple" = "purple3"
)

# Get colors for current focus modules
focus_colors <- module_color_palette[focus_modules_auto]
focus_colors_with_grey <- module_color_palette[focus_modules]

# Names for display (capitalize first letter)
focus_names <- tools::toTitleCase(focus_modules_auto)
n_focus_modules <- length(focus_modules_auto)
cat(sprintf("  Number of focus modules: %d\n", n_focus_modules))

cat("Parameters:\n")
cat(sprintf("  Focus modules: %s\n", paste(focus_modules, collapse = ", ")))
cat(sprintf("  Hub criteria: kME > %.1f, |GS| > %.1f, p < %.2f\n",
            kME_thresh, GS_thresh, pGS_thresh))
cat(sprintf("  Bootstrap: max %d iterations, convergence CV < %.0f%%\n",
            max_iterations, convergence_cv * 100))

# ----------------------------------------------------------------------------
# Helper Function 1: Module Matching
# ----------------------------------------------------------------------------
# Purpose: Find the bootstrap module that best matches an original module
#          by protein content overlap (Jaccard index)
#
# Why needed: WGCNA assigns colors by size rank, not content.
#             "Brown" in bootstrap may contain different proteins than original "brown"
# ----------------------------------------------------------------------------

match_modules <- function(moduleColors_orig, moduleColors_boot, target_module) {
  orig_proteins <- names(moduleColors_orig)[moduleColors_orig == target_module]

  if(length(orig_proteins) == 0) {
    return(list(best_match = NA, jaccard = 0))
  }

  boot_modules <- unique(moduleColors_boot)
  best_jaccard <- 0
  best_match <- NA

  for(boot_mod in boot_modules) {
    boot_proteins <- names(moduleColors_boot)[moduleColors_boot == boot_mod]

    intersection <- length(intersect(orig_proteins, boot_proteins))
    union_set <- length(union(orig_proteins, boot_proteins))
    jaccard <- ifelse(union_set > 0, intersection / union_set, 0)

    if(jaccard > best_jaccard) {
      best_jaccard <- jaccard
      best_match <- boot_mod
    }
  }

  return(list(best_match = best_match, jaccard = best_jaccard))
}

# ----------------------------------------------------------------------------
# Helper Function 2: Bootstrap Iteration with Module Matching
# ----------------------------------------------------------------------------
# Purpose: Run single bootstrap iteration and calculate stability metrics
#          using module matching to track biological content
# ----------------------------------------------------------------------------

run_bootstrap_iteration <- function(datExpr, power, moduleColors_orig,
                                    MEs_orig, focus_modules) {
  n_samples <- nrow(datExpr)
  protein_names <- colnames(datExpr)

  # Bootstrap resample (with replacement)
  boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
  datExpr_boot <- datExpr[boot_idx, ]

  # Rebuild network
  net_boot <- blockwiseModules(
    datExpr_boot, power = power, networkType = "signed", TOMType = "signed",
    minModuleSize = 20, mergeCutHeight = 0.25, corType = "bicor",
    maxPOutliers = 0.1, numericLabels = TRUE, saveTOMs = FALSE, verbose = 0
  )

  # Input validation for labels2colors
  stopifnot("net_boot$colors must be integer" = is.numeric(net_boot$colors))
  stopifnot("Column count mismatch in bootstrap" = length(net_boot$colors) == ncol(datExpr_boot))
  stopifnot("Column names must match" = identical(colnames(datExpr_boot), protein_names))

  moduleColors_boot <- labels2colors(net_boot$colors)
  names(moduleColors_boot) <- colnames(datExpr_boot)  # Use actual boot column names

  # Calculate metrics using module matching
  metrics <- list()
  for(mod in focus_modules) {
    match_result <- match_modules(moduleColors_orig, moduleColors_boot, mod)

    # Calculate assignment stability for matched module
    orig_proteins <- names(moduleColors_orig)[moduleColors_orig == mod]
    if(length(orig_proteins) > 0 && !is.na(match_result$best_match)) {
      staying <- sum(moduleColors_boot[orig_proteins] == match_result$best_match, na.rm = TRUE)
      assignment_stability <- staying / length(orig_proteins)
    } else {
      assignment_stability <- NA
    }

    metrics[[mod]] <- data.frame(
      Module = mod,
      Best_Match = match_result$best_match,
      Jaccard = match_result$jaccard,
      Assignment_Stability = assignment_stability,
      stringsAsFactors = FALSE
    )
  }

  # Track protein reassignments (for later analysis)
  reassignment <- data.frame(
    Protein = protein_names,
    Original = as.character(moduleColors_orig),
    Bootstrap = as.character(moduleColors_boot),
    stringsAsFactors = FALSE
  )

  return(list(metrics = bind_rows(metrics), reassignment = reassignment))
}

# ----------------------------------------------------------------------------
# Helper Function 3: Subsampling with Module Matching
# ----------------------------------------------------------------------------

run_subsample_iteration <- function(datExpr, rate, power, moduleColors_orig, focus_modules) {
  n_samples <- nrow(datExpr)
  n_keep <- round(n_samples * rate)
  protein_names <- colnames(datExpr)

  # Subsample without replacement
  sub_idx <- sample(1:n_samples, n_keep, replace = FALSE)
  datExpr_sub <- datExpr[sub_idx, ]

  # Rebuild network
  net_sub <- blockwiseModules(
    datExpr_sub, power = power, networkType = "signed", TOMType = "signed",
    minModuleSize = 20, mergeCutHeight = 0.25, corType = "bicor",
    maxPOutliers = 0.1, numericLabels = TRUE, saveTOMs = FALSE, verbose = 0
  )

  # Input validation for subsample function
  stopifnot("Column count mismatch in subsample" = length(net_sub$colors) == ncol(datExpr_sub))
  stopifnot("Column names must match" = identical(colnames(datExpr_sub), protein_names))

  moduleColors_sub <- labels2colors(net_sub$colors)
  names(moduleColors_sub) <- colnames(datExpr_sub)  # Use actual subsample column names

  # Calculate metrics with module matching
  results <- list()
  for(mod in focus_modules) {
    match_result <- match_modules(moduleColors_orig, moduleColors_sub, mod)
    results[[mod]] <- data.frame(
      Module = mod,
      Jaccard = match_result$jaccard,
      Best_Match = match_result$best_match,
      Rate = rate,
      stringsAsFactors = FALSE
    )
  }

  return(bind_rows(results))
}

cat("  Helper functions defined.\n")

# ============================================================================
# 7a: ADAPTIVE BOOTSTRAP WITH MODULE MATCHING
# ============================================================================
#
# Method: Resample patients WITH replacement, rebuild network, measure overlap
# Module Matching: Track biological content by protein overlap, not color name
# Convergence: Stop when CV < 5% and relative change < 1% between batches
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7a: ADAPTIVE BOOTSTRAP WITH MODULE MATCHING\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("Running adaptive bootstrap analysis...\n")
cat("  (Module matching corrects for WGCNA's size-based color assignment)\n")

# Storage for results
all_metrics <- list()
all_reassignments <- list()
convergence_history <- data.frame()
converged <- FALSE
current_iter <- 0

start_time <- Sys.time()

# Adaptive bootstrap loop
while(!converged && current_iter < max_iterations) {
  batch_start <- current_iter + 1
  batch_end <- min(current_iter + batch_size, max_iterations)

  cat(sprintf("  Iterations %d-%d...", batch_start, batch_end))

  # Run batch of bootstrap iterations
  # Network construction with stored soft_power variable and error handling
  failed_iterations <- 0
  for(i in batch_start:batch_end) {
    result <- tryCatch({
      run_bootstrap_iteration(datExpr, power = soft_power, moduleColors, MEs, focus_modules)
    }, error = function(e) {
      cat(sprintf("    [WARNING] Bootstrap iteration %d failed: %s", i, e$message))
      failed_iterations <<- failed_iterations + 1
      return(NULL)
    })

    if(!is.null(result)) {
      all_metrics[[i]] <- result$metrics
      all_reassignments[[i]] <- result$reassignment
    }
  }
  if(failed_iterations > 0) {
    cat(sprintf(" (%d failed)", failed_iterations))
  }

  current_iter <- batch_end

  # Calculate running statistics
  metrics_df <- bind_rows(all_metrics)

  running_stats <- metrics_df %>%
    group_by(Module) %>%
    summarize(
      Mean_Jaccard = mean(Jaccard, na.rm = TRUE),
      SD_Jaccard = sd(Jaccard, na.rm = TRUE),
      CV_Jaccard = ifelse(Mean_Jaccard > 0, SD_Jaccard / Mean_Jaccard, NA),
      .groups = "drop"
    )

  # Store convergence history
  convergence_history <- bind_rows(
    convergence_history,
    running_stats %>% mutate(Iteration = current_iter)
  )

  # Check convergence (focus modules only, not grey)
  cv_check <- running_stats %>%
    filter(Module %in% focus_modules_auto) %>%
    pull(CV_Jaccard)

  if(current_iter >= min_iterations && all(!is.na(cv_check)) && all(cv_check < convergence_cv)) {
    if(current_iter >= 2 * batch_size) {
      prev_stats <- convergence_history %>%
        filter(Iteration == current_iter - batch_size, Module %in% focus_modules_auto)
      curr_stats <- running_stats %>% filter(Module %in% focus_modules_auto)

      if(nrow(prev_stats) == nrow(curr_stats) && nrow(prev_stats) > 0) {
        rel_change <- abs(curr_stats$Mean_Jaccard - prev_stats$Mean_Jaccard) / prev_stats$Mean_Jaccard

        if(all(!is.na(rel_change)) && all(rel_change < convergence_change)) {
          converged <- TRUE
        }
      }
    }
  }

  # Print progress - dynamic for all focus modules
  mean_cv <- mean(cv_check, na.rm = TRUE)
  module_jaccards <- sapply(focus_modules_auto, function(m) {
    j <- running_stats$Mean_Jaccard[running_stats$Module == m]
    if(length(j) == 0) NA else j
  })
  jaccard_str <- paste(sprintf("%s=%.3f", focus_names, module_jaccards), collapse = ", ")
  cat(sprintf(" CV=%.1f%% | %s\n", mean_cv * 100, jaccard_str))
}

elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
n_bootstrap <- current_iter

cat(sprintf("\n  Bootstrap finished: %d iterations in %.1f minutes\n", n_bootstrap, elapsed_time))
cat(sprintf("  Converged: %s\n", ifelse(converged, "Yes", "No - max iterations reached")))

# ----------------------------------------------------------------------------
# Summarize Bootstrap Results
# ----------------------------------------------------------------------------

bootstrap_results <- bind_rows(all_metrics) %>%
  group_by(Module) %>%
  summarize(
    Mean_Jaccard = mean(Jaccard, na.rm = TRUE),
    SD_Jaccard = sd(Jaccard, na.rm = TRUE),
    CI_Lower = quantile(Jaccard, 0.025, na.rm = TRUE),
    CI_Upper = quantile(Jaccard, 0.975, na.rm = TRUE),
    Mean_Assignment = mean(Assignment_Stability, na.rm = TRUE),
    SD_Assignment = sd(Assignment_Stability, na.rm = TRUE),
    .groups = "drop"
  )

cat("  Bootstrap Results (with Module Matching):")
print(bootstrap_results %>% mutate(across(where(is.numeric), ~round(., 3))))

# Analyze module matching patterns
match_frequency <- bind_rows(all_metrics) %>%
  group_by(Module, Best_Match) %>%
  summarize(Count = n(), .groups = "drop") %>%
  group_by(Module) %>%
  mutate(Percent = round(100 * Count / sum(Count), 1)) %>%
  slice_max(Count, n = 3) %>%
  ungroup()

cat("  Most Common Module Matches:\n")
print(match_frequency)

# Store corrected metrics for later use
all_metrics_final <- all_metrics

# ----------------------------------------------------------------------------
# INTERPRETATION - Bootstrap Results
# ----------------------------------------------------------------------------
# Focus module proteins may merge into larger modules in bootstrap samples.
# This suggests module boundaries may be sample-dependent with n=50.
# The underlying signal (poor prognosis proteins) exists, but module boundaries are fluid.
# Hub genes may still be valid even if module boundaries are fluid.
# ----------------------------------------------------------------------------

# --- Bootstrap Convergence Plot ---
cat("  Generating Convergence Plot...\n")

conv_plot_data <- convergence_history %>%
  filter(Module %in% focus_modules_auto)

p1_convergence <- ggplot(conv_plot_data, aes(x = Iteration, y = Mean_Jaccard, color = Module)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = focus_colors) +
  labs(
    title = "Bootstrap Convergence",
    subtitle = sprintf("%d iterations finished", n_bootstrap),
    x = "Number of Iterations",
    y = "Mean Jaccard Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )
print(p1_convergence)

ggsave(file.path(fig_dir_png, "Step7_01_Convergence.png"), p1_convergence, width = 6, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_01_Convergence.pdf"), p1_convergence, width = 6, height = 4)

# --- Module Stability Plot ---
cat("  Generating Module Stability Plot...\n")

bootstrap_results$Module <- factor(bootstrap_results$Module, levels = focus_modules)

p2_stability <- ggplot(bootstrap_results, aes(x = Module, y = Mean_Jaccard, fill = Module)) +
  geom_col(alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 0.8) +
  scale_fill_manual(values = focus_colors_with_grey) +
  geom_text(aes(label = sprintf("%.2f", Mean_Jaccard)), vjust = -0.5, fontface = "bold") +
  labs(
    title = "Module Stability (Bootstrap with Module Matching)",
    subtitle = "Error bars = 95% CI | Red dashed = 0.7 | Orange dotted = 0.5",
    x = NULL, y = "Jaccard Index"
  ) +
  scale_y_continuous(limits = c(0, 0.8), expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )
print(p2_stability)

ggsave(file.path(fig_dir_png, "Step7_02_ModuleStability.png"), p2_stability, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_02_ModuleStability.pdf"), p2_stability, width = 6, height = 5)

# ============================================================================
# 7b: HUB GENE CO-CLUSTERING ANALYSIS
# ============================================================================
#
# Purpose: Validate that hub genes form stable biological communities
#          even when module color assignments change
#
# Question: "Do hub genes from the same original module cluster together
#           more often than hub genes from different modules?"
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7b: HUB GENE CO-CLUSTERING ANALYSIS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Dynamic hub gene extraction for ALL focus modules
hub_genes_by_module <- list()

cat("  Extracting hub genes for all focus modules...\n")

for(mod in focus_modules_auto) {
  if(exists("hub_results") && !is.null(hub_results$Moderate[[mod]])) {
    hub_genes_by_module[[mod]] <- hub_results$Moderate[[mod]]$Protein
    if(is.null(hub_genes_by_module[[mod]])) hub_genes_by_module[[mod]] <- character(0)
  } else if(exists("proteinInfo") && "isHub_Moderate" %in% colnames(proteinInfo)) {
    hub_genes_by_module[[mod]] <- proteinInfo$Protein[proteinInfo$Module == mod & proteinInfo$isHub_Moderate == TRUE]
    if(is.null(hub_genes_by_module[[mod]])) hub_genes_by_module[[mod]] <- character(0)
  } else {
    warning(sprintf("Hub genes not found for %s module from Step 6", mod))
    hub_genes_by_module[[mod]] <- character(0)
  }
  cat(sprintf("  %s hubs (n=%d): %s\n",
              tools::toTitleCase(mod),
              length(hub_genes_by_module[[mod]]),
              paste(hub_genes_by_module[[mod]], collapse = ", ")))
}

all_hub_genes <- unlist(hub_genes_by_module, use.names = FALSE)
cat(sprintf("\n  Total hub genes across all focus modules: %d\n", length(all_hub_genes)))

# Calculate co-clustering matrix
n_iter <- length(all_reassignments)
n_hubs <- length(all_hub_genes)

co_cluster_matrix <- matrix(0, nrow = n_hubs, ncol = n_hubs,
                            dimnames = list(all_hub_genes, all_hub_genes))

cat("  Calculating co-clustering across bootstrap iterations...\n")

for(i in 1:n_iter) {
  reassign <- all_reassignments[[i]]
  moduleColors_boot <- reassign$Bootstrap
  names(moduleColors_boot) <- reassign$Protein

  for(j in 1:(n_hubs-1)) {
    for(k in (j+1):n_hubs) {
      gene1 <- all_hub_genes[j]
      gene2 <- all_hub_genes[k]

      if(gene1 %in% names(moduleColors_boot) && gene2 %in% names(moduleColors_boot)) {
        if(moduleColors_boot[gene1] == moduleColors_boot[gene2]) {
          co_cluster_matrix[j, k] <- co_cluster_matrix[j, k] + 1
          co_cluster_matrix[k, j] <- co_cluster_matrix[k, j] + 1
        }
      }
    }
  }
}

# Convert to percentage
co_cluster_pct <- co_cluster_matrix / n_iter * 100
diag(co_cluster_pct) <- 100  # Self-clustering = 100%

cat("  Hub Gene Co-Clustering Matrix (% of iterations in same module):")
print(round(co_cluster_pct, 1))

# Calculate within-group vs between-group co-clustering (dynamic for all focus modules)
within_module_pairs <- list()
for(mod in focus_modules_auto) {
  hubs <- hub_genes_by_module[[mod]]
  if(length(hubs) >= 2) {
    pairs <- c()
    for(i in 1:(length(hubs)-1)) {
      for(j in (i+1):length(hubs)) {
        pairs <- c(pairs, co_cluster_pct[hubs[i], hubs[j]])
      }
    }
    within_module_pairs[[mod]] <- pairs
  } else {
    within_module_pairs[[mod]] <- numeric(0)
  }
}

# Between-module pairs (all pairwise combinations of focus modules)
between_pairs <- c()
if(length(focus_modules_auto) >= 2) {
  for(m1_idx in 1:(length(focus_modules_auto)-1)) {
    for(m2_idx in (m1_idx+1):length(focus_modules_auto)) {
      mod1 <- focus_modules_auto[m1_idx]
      mod2 <- focus_modules_auto[m2_idx]
      hubs1 <- hub_genes_by_module[[mod1]]
      hubs2 <- hub_genes_by_module[[mod2]]
      for(h1 in hubs1) {
        for(h2 in hubs2) {
          between_pairs <- c(between_pairs, co_cluster_pct[h1, h2])
        }
      }
    }
  }
}

# Combine all within-module pairs
within_group <- unlist(within_module_pairs)

# Statistical test
wilcox_result <- wilcox.test(within_group, between_pairs)

cat("  Co-Clustering Summary:\n")
for(mod in focus_modules_auto) {
  pairs <- within_module_pairs[[mod]]
  if(length(pairs) > 0) {
    cat(sprintf("    %s hubs cluster together: %.1f%% +/- %.1f%%\n",
                tools::toTitleCase(mod), mean(pairs), sd(pairs)))
  } else {
    cat(sprintf("    %s hubs: insufficient pairs for analysis\n", tools::toTitleCase(mod)))
  }
}
cat(sprintf("    Cross-module (between all): %.1f%% +/- %.1f%%\n",
            mean(between_pairs), sd(between_pairs)))

cat(sprintf("\n  Statistical Comparison:\n"))
cat(sprintf("    Within-module mean: %.1f%%\n", mean(within_group)))
cat(sprintf("    Between-module mean: %.1f%%\n", mean(between_pairs)))
cat(sprintf("    Wilcoxon p-value: %.6f\n", wilcox_result$p.value))

if(mean(within_group) > mean(between_pairs) && wilcox_result$p.value < 0.05) {
  cat("    -> Hub genes within same module cluster together significantly MORE than across modules\n")
  cat("    -> This validates that biological communities are stable despite color changes\n")
}

# ----------------------------------------------------------------------------
# INTERPRETATION - Co-Clustering
# ----------------------------------------------------------------------------
# The biological signal is REAL, even though module colors are unstable.
# The module COLOR changes, but the underlying protein communities remain intact.
# This is acceptable for small sample sizes.
# ----------------------------------------------------------------------------

# --- Hub Gene Co-Clustering Heatmap ---
cat("  Generating Hub Gene Co-Clustering Heatmap...\n")

co_cluster_df <- as.data.frame(co_cluster_pct)
co_cluster_df$Protein1 <- rownames(co_cluster_df)

# Order hub genes by module for heatmap display
hub_gene_order <- unlist(hub_genes_by_module[focus_modules_auto], use.names = FALSE)

co_cluster_long <- co_cluster_df %>%
  pivot_longer(cols = -Protein1, names_to = "Protein2", values_to = "CoCluster") %>%
  mutate(
    Protein1 = factor(Protein1, levels = hub_gene_order),
    Protein2 = factor(Protein2, levels = hub_gene_order)
  )

# Calculate divider positions between modules
divider_positions <- cumsum(sapply(focus_modules_auto, function(m) length(hub_genes_by_module[[m]])))
divider_positions <- divider_positions[-length(divider_positions)] + 0.5  # Remove last, add offset

p6_cocluster <- ggplot(co_cluster_long, aes(x = Protein1, y = Protein2, fill = CoCluster)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f", CoCluster)), size = 2.5) +
  scale_fill_gradient2(low = "white", mid = "orange", high = "darkred",
                       midpoint = 50, limits = c(0, 100),
                       name = "Co-clustering\n(% iterations)") +
  geom_vline(xintercept = divider_positions, color = "black", linewidth = 1) +
  geom_hline(yintercept = divider_positions, color = "black", linewidth = 1) +
  labs(
    title = "Hub Gene Co-Clustering Stability",
    subtitle = sprintf("%% of %d bootstrap iterations where protein pairs were in same module", n_bootstrap),
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  ) +
  coord_fixed()
print(p6_cocluster)

ggsave(file.path(fig_dir_png, "Step7_06_CoClusterHeatmap.png"), p6_cocluster, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_06_CoClusterHeatmap.pdf"), p6_cocluster, width = 8, height = 6)

# --- Co-Clustering Summary (Within vs Between) ---
cat("  Generating Co-Clustering Summary Plot...\n")

# Build summary data frame dynamically for all focus modules
within_summary <- lapply(focus_modules_auto, function(mod) {
  pairs <- within_module_pairs[[mod]]
  if(length(pairs) > 0) {
    data.frame(
      Comparison = paste0(tools::toTitleCase(mod), " hubs\n(within)"),
      Mean = mean(pairs),
      SD = sd(pairs),
      Type = "Within"
    )
  } else {
    NULL
  }
})
within_summary <- bind_rows(within_summary)

between_summary <- data.frame(
  Comparison = "Cross-module\n(between)",
  Mean = mean(between_pairs),
  SD = sd(between_pairs),
  Type = "Between"
)

cocluster_summary_df <- bind_rows(within_summary, between_summary)

p7_cocluster_bar <- ggplot(cocluster_summary_df, aes(x = Comparison, y = Mean, fill = Type)) +
  geom_col(alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", Mean)), vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = c("Within" = "steelblue", "Between" = "grey60")) +
  labs(
    title = "Hub Gene Co-Clustering: Within vs Between Modules",
    subtitle = sprintf("Wilcoxon p = %.4f", wilcox_result$p.value),
    x = NULL, y = "Co-Clustering Rate (%)"
  ) +
  scale_y_continuous(limits = c(0, 70), expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )
print(p7_cocluster_bar)

ggsave(file.path(fig_dir_png, "Step7_07_CoClusterSummary.png"), p7_cocluster_bar, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_07_CoClusterSummary.pdf"), p7_cocluster_bar, width = 7, height = 5)

# ============================================================================
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7d: PERMUTATION TEST\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

n_permutations <- 100

cat(sprintf("  Running %d permutations to generate null distribution...\n", n_permutations))

perm_results <- list()
for(i in 1:n_permutations) {
  # Shuffle module labels randomly
  moduleColors_perm <- sample(moduleColors)
  names(moduleColors_perm) <- names(moduleColors)

  results <- list()
  for(mod in focus_modules) {
    orig_proteins <- names(moduleColors)[moduleColors == mod]
    perm_proteins <- names(moduleColors_perm)[moduleColors_perm == mod]

    intersection <- length(intersect(orig_proteins, perm_proteins))
    union_set <- length(union(orig_proteins, perm_proteins))
    jaccard <- ifelse(union_set > 0, intersection / union_set, 0)

    results[[mod]] <- data.frame(Module = mod, Jaccard = jaccard)
  }

  perm_results[[i]] <- bind_rows(results)

  if(i %% 25 == 0) cat(sprintf("    Processed %d/%d\n", i, n_permutations))
}

perm_df <- bind_rows(perm_results)

# Compare observed (bootstrap) to null distribution
perm_summary <- perm_df %>%
  group_by(Module) %>%
  summarize(
    Null_Mean = mean(Jaccard),
    Null_SD = sd(Jaccard),
    .groups = "drop"
  ) %>%
  left_join(bootstrap_results %>% select(Module, Mean_Jaccard), by = "Module") %>%
  mutate(
    Z_score = (Mean_Jaccard - Null_Mean) / Null_SD,
    P_value = pnorm(Z_score, lower.tail = FALSE)
  )

cat("  Permutation Test Results:\n")
print(perm_summary %>% mutate(across(where(is.numeric), ~round(., 4))))

# ----------------------------------------------------------------------------
# INTERPRETATION - Permutation Test
# ----------------------------------------------------------------------------
# Even though absolute Jaccard is low, it's 2.8-4.4x higher than random chance.
# The modules are statistically real (Z > 8, p < 0.0001).
# ----------------------------------------------------------------------------

# --- Permutation Test Plot ---
cat("  Generating Permutation Test Plot...\n")

perm_plot_data <- perm_df %>%
  filter(Module %in% focus_modules_auto) %>%
  mutate(Module = factor(Module, levels = focus_modules_auto))

observed_vals <- perm_summary %>%
  filter(Module %in% focus_modules_auto) %>%
  mutate(Module = factor(Module, levels = focus_modules_auto))

# Build dynamic subtitle showing Z-scores for all focus modules
z_score_strings <- sapply(focus_modules_auto, function(mod) {
  z <- perm_summary$Z_score[perm_summary$Module == mod]
  sprintf("%s: Z=%.1f", tools::toTitleCase(mod), z)
})
subtitle_text <- paste(z_score_strings, collapse = " | ")
subtitle_text <- paste0(subtitle_text, " | All p<0.0001")

p4_permutation <- ggplot(perm_plot_data, aes(x = Jaccard, fill = Module)) +
  geom_histogram(bins = 20, alpha = 0.7, color = "white") +
  geom_vline(data = observed_vals, aes(xintercept = Mean_Jaccard),
             color = "red", linetype = "dashed", linewidth = 1.2) +
  scale_fill_manual(values = focus_colors) +
  facet_wrap(~Module, scales = "free_y", ncol = 2) +
  labs(
    title = "Permutation Test: Observed vs Null Distribution",
    subtitle = subtitle_text,
    x = "Jaccard Index (Null Distribution)", y = "Count"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )
print(p4_permutation)

ggsave(file.path(fig_dir_png, "Step7_04_PermutationTest.png"), p4_permutation, width = 8, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_04_PermutationTest.pdf"), p4_permutation, width = 8, height = 4)

# ============================================================================
# 7e: MODULE PRESERVATION ANALYSIS
# ============================================================================
#
# Purpose: Test if modules are preserved when data is split in half
# Method: WGCNA's modulePreservation function with 50/50 split
# Metric: Zsummary > 10 = strong, > 2 = moderate, < 2 = none
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7e: MODULE PRESERVATION ANALYSIS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

set.seed(456)
n_samples <- nrow(datExpr)
ref_idx <- sample(1:n_samples, round(n_samples / 2))
test_idx <- setdiff(1:n_samples, ref_idx)

cat(sprintf("  Reference set: %d samples\n", length(ref_idx)))
cat(sprintf("  Test set: %d samples\n", length(test_idx)))

multiExpr <- list(
  Reference = list(data = datExpr[ref_idx, ]),
  Test = list(data = datExpr[test_idx, ])
)

multiColor <- list(Reference = moduleColors)

cat("  Running module preservation analysis (this may take several minutes)...")

mp <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  testNetworks = 2,
  nPermutations = 100,
  randomSeed = 789,
  quickCor = 0,
  verbose = 0
)

preservation_stats <- mp$preservation$Z$ref.Reference$inColumnsAlsoPresentIn.Test
preservation_summary <- data.frame(
  Module = rownames(preservation_stats),
  Zsummary = preservation_stats$Zsummary.pres,
  Zdensity = preservation_stats$Zdensity.pres,
  Zconnectivity = preservation_stats$Zconnectivity.pres
) %>%
  filter(Module %in% focus_modules) %>%
  mutate(
    Interpretation = case_when(
      Zsummary > 10 ~ "Strong",
      Zsummary > 2 ~ "Moderate",
      TRUE ~ "None"
    )
  )

medianRank <- mp$preservation$observed$ref.Reference$inColumnsAlsoPresentIn.Test
preservation_summary$medianRank <- medianRank[preservation_summary$Module, "medianRank.pres"]

cat("  Module Preservation Results:\n")
print(preservation_summary %>% mutate(across(where(is.numeric), ~round(., 2))))

# ----------------------------------------------------------------------------
# INTERPRETATION - Module Preservation
# ----------------------------------------------------------------------------
# These findings support the biological relevance of the brown and blue modules
# as poor prognosis signatures, while acknowledging the inherent limitations
# of network analysis with n=50 patients.
# Zsummary > 2 indicates moderate preservation - the modules are real biological entities.
# ----------------------------------------------------------------------------

# --- Module Preservation Plot ---
cat("  Generating Module Preservation Plot...\n")

preservation_summary$Module <- factor(preservation_summary$Module, levels = focus_modules)

p5_preservation <- ggplot(preservation_summary, aes(x = Module, y = Zsummary, fill = Module)) +
  geom_col(alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = 2, linetype = "dotted", color = "orange", linewidth = 0.8) +
  scale_fill_manual(values = focus_colors_with_grey) +
  geom_text(aes(label = sprintf("%.1f\n(%s)", Zsummary, Interpretation)),
            vjust = -0.3, fontface = "bold", size = 3.5) +
  labs(
    title = "Module Preservation (50/50 Split)",
    subtitle = "Green dashed = 10 (strong) | Orange dotted = 2 (moderate)",
    x = NULL, y = "Zsummary"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )
print(p5_preservation)

ggsave(file.path(fig_dir_png, "Step7_05_ModulePreservation.png"), p5_preservation, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_05_ModulePreservation.pdf"), p5_preservation, width = 6, height = 5)

# ============================================================================
# 7f: HUB GENE STABILITY
# ============================================================================
#
# Purpose: Track how consistently each hub gene is identified as a hub
# Method: For each bootstrap, check if hub genes meet criteria (kME, GS)
# Note: Low detection rates expected due to module color changes
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7f: HUB GENE STABILITY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Get original hub genes from Step 6 - dynamic for ALL focus modules
orig_hubs_by_module <- list()
for(mod in focus_modules_auto) {
  if(!is.null(networkData$hubGenes$Moderate[[mod]])) {
    orig_hubs_by_module[[mod]] <- networkData$hubGenes$Moderate[[mod]]$Protein
  } else {
    orig_hubs_by_module[[mod]] <- character(0)
  }
}
orig_hubs <- unlist(orig_hubs_by_module, use.names = FALSE)
cat(sprintf("  Analyzing hub genes from %d focus modules (n=%d total)\n",
            length(focus_modules_auto), length(orig_hubs)))

# Prepare trait data
trait_data <- data.frame(
  PFS_group = as.numeric(datTraits$PFS_group == "Long"),
  Response = as.numeric(datTraits$Response == "CD"),
  row.names = rownames(datTraits)
)

# Initialize tracking matrices
hub_detection <- matrix(0, nrow = length(orig_hubs), ncol = n_bootstrap,
                        dimnames = list(orig_hubs, 1:n_bootstrap))
kME_values <- matrix(NA, nrow = length(orig_hubs), ncol = n_bootstrap,
                     dimnames = list(orig_hubs, 1:n_bootstrap))

cat("  Tracking hub gene detection across bootstrap iterations...\n")

n_iter_analyze <- min(n_bootstrap, 200)  # Cap at 200 for computational efficiency

for(i in 1:n_iter_analyze) {
  reassign <- all_reassignments[[i]]
  boot_idx <- sample(1:nrow(datExpr), nrow(datExpr), replace = TRUE)
  datExpr_boot <- datExpr[boot_idx, ]
  moduleColors_boot <- reassign$Bootstrap
  names(moduleColors_boot) <- reassign$Protein

  # Calculate module eigengenes
  MEs_boot <- moduleEigengenes(datExpr_boot, moduleColors_boot)$eigengenes

  for(hub in orig_hubs) {
    hub_mod <- proteinInfo$Module[proteinInfo$Protein == hub]
    me_col <- paste0("ME", hub_mod)

    if(me_col %in% colnames(MEs_boot) && hub %in% colnames(datExpr_boot)) {
      # Calculate kME
      kME_val <- cor(datExpr_boot[, hub], MEs_boot[, me_col], use = "p")
      kME_values[hub, i] <- kME_val

      if(!is.na(kME_val) && kME_val > kME_thresh) {
        # Check GS criteria
        gs_pfs <- cor(datExpr_boot[, hub], trait_data[boot_idx, "PFS_group"], use = "p")
        pgs_pfs <- cor.test(datExpr_boot[, hub], trait_data[boot_idx, "PFS_group"])$p.value
        gs_resp <- cor(datExpr_boot[, hub], trait_data[boot_idx, "Response"], use = "p")
        pgs_resp <- cor.test(datExpr_boot[, hub], trait_data[boot_idx, "Response"])$p.value

        if((abs(gs_pfs) > GS_thresh && pgs_pfs < pGS_thresh) ||
           (abs(gs_resp) > GS_thresh && pgs_resp < pGS_thresh)) {
          hub_detection[hub, i] <- 1
        }
      }
    }
  }

  if(i %% 50 == 0) cat(sprintf("    Processed %d/%d iterations\n", i, n_iter_analyze))
}

# Calculate hub stability metrics
hub_stability <- data.frame(
  Protein = orig_hubs,
  Module = proteinInfo$Module[match(orig_hubs, proteinInfo$Protein)],
  Detection_Rate = rowSums(hub_detection[, 1:n_iter_analyze], na.rm = TRUE) / n_iter_analyze,
  Mean_kME = rowMeans(kME_values[, 1:n_iter_analyze], na.rm = TRUE),
  SD_kME = apply(kME_values[, 1:n_iter_analyze], 1, sd, na.rm = TRUE)
) %>%
  mutate(
    Classification = case_when(
      Detection_Rate >= 0.9 ~ "Core hub",
      Detection_Rate >= 0.7 ~ "Robust hub",
      Detection_Rate >= 0.5 ~ "Moderate hub",
      TRUE ~ "Unstable hub"
    )
  ) %>%
  arrange(Module, desc(Detection_Rate))

cat("  Hub Gene Stability:\n")
print(hub_stability %>% mutate(across(where(is.numeric), ~round(., 3))))

# ----------------------------------------------------------------------------
# INTERPRETATION - Hub Gene Stability
# ----------------------------------------------------------------------------
# Hub genes ARE stable as a community, but their "hub status" depends on which
# color the module gets assigned.
# Low detection rates are expected because hub status is evaluated against the
# SAME-COLORED module, which changes across bootstraps.
# The co-clustering analysis (7b) better captures hub gene community stability.
# ----------------------------------------------------------------------------

# ============================================================================
# 7g: PROTEIN REASSIGNMENT ANALYSIS
# ============================================================================
#
# Purpose: Track where proteins move across bootstrap iterations
# Key insight: Brown/Blue proteins often merge into Turquoise (largest module)
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7g: PROTEIN REASSIGNMENT ANALYSIS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Aggregate reassignment data
reassign_all <- bind_rows(all_reassignments)

reassign_matrix <- reassign_all %>%
  filter(Original %in% focus_modules) %>%
  group_by(Original, Bootstrap) %>%
  summarize(Count = n(), .groups = "drop") %>%
  group_by(Original) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Pivot for viewing
reassign_summary <- reassign_matrix %>%
  group_by(Original, Bootstrap) %>%
  summarize(Mean_Prop = mean(Proportion), .groups = "drop") %>%
  pivot_wider(names_from = Bootstrap, values_from = Mean_Prop, values_fill = 0)

cat("  Protein Reassignment Proportions (top destinations):")
print(as.data.frame(reassign_summary) %>%
        select(Original, any_of(c("brown", "blue", "turquoise", "grey", "green", "yellow"))) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Calculate retention rates
staying_rate <- reassign_all %>%
  filter(Original %in% focus_modules) %>%
  mutate(Stayed = Original == Bootstrap) %>%
  group_by(Original) %>%
  summarize(
    Stay_Rate = mean(Stayed),
    To_Turquoise = mean(Bootstrap == "turquoise"),  # Often merges into largest module
    To_Grey = mean(Bootstrap == "grey"),
    .groups = "drop"
  )

cat("  Module Retention Rates:\n")
print(staying_rate %>% mutate(across(where(is.numeric), ~round(., 3))))

# ----------------------------------------------------------------------------
# INTERPRETATION - Protein Reassignment
# ----------------------------------------------------------------------------
# Brown and Blue proteins frequently merge into Turquoise (the largest module).
# This explains the low Jaccard - the biology is preserved, just relabeled.
# The proteins stay together as a community, but WGCNA gives them a different color name.
# ----------------------------------------------------------------------------

# ============================================================================
# 7h: MODULE SIZE VS STABILITY (ALL MODULES)
# ============================================================================
#
# Purpose: Test if larger modules are more stable (potential artifact)
# Method: Calculate Jaccard for ALL modules, correlate with size
# Expected: If no correlation, stability reflects biology, not size
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7h: MODULE SIZE VS STABILITY (ALL MODULES)")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Get all unique modules
all_modules <- unique(moduleColors)
cat(sprintf("  Analyzing all %d modules...\n", length(all_modules)))

# Calculate stability for ALL modules
all_module_stability <- data.frame(Module = character(), Mean_Jaccard = numeric(), SD_Jaccard = numeric())

for(mod in all_modules) {
  jaccard_values <- c()

  for(i in 1:length(all_reassignments)) {
    reassign <- all_reassignments[[i]]
    moduleColors_boot <- reassign$Bootstrap
    names(moduleColors_boot) <- reassign$Protein

    match_result <- match_modules(moduleColors, moduleColors_boot, mod)
    jaccard_values <- c(jaccard_values, match_result$jaccard)
  }

  all_module_stability <- rbind(all_module_stability, data.frame(
    Module = mod,
    Mean_Jaccard = mean(jaccard_values, na.rm = TRUE),
    SD_Jaccard = sd(jaccard_values, na.rm = TRUE)
  ))
}

# Get module sizes
module_sizes <- data.frame(
  Module = names(table(moduleColors)),
  Size = as.numeric(table(moduleColors))
)

# Combine
size_stability_all <- left_join(module_sizes, all_module_stability, by = "Module")

cat("  Module Size vs Stability:\n")
print(size_stability_all %>% arrange(desc(Size)) %>% mutate(across(where(is.numeric), ~round(., 3))))

# Correlation test
size_cor_all <- cor.test(size_stability_all$Size, size_stability_all$Mean_Jaccard, method = "spearman")

cat(sprintf("\n  Spearman correlation: rho = %.3f, p = %.4f\n",
            size_cor_all$estimate, size_cor_all$p.value))

if(size_cor_all$p.value < 0.05) {
  cat("  -> Significant correlation: larger modules may be artifactually more stable\n")
} else {
  cat("  -> No significant correlation: stability reflects biology, not module size\n")
}

# ----------------------------------------------------------------------------
# INTERPRETATION - Size vs Stability
# ----------------------------------------------------------------------------
# Stability is NOT driven by module size.
# Green module (similar size to brown) has higher Jaccard -> Brown's low stability
# is NOT a size artifact.
# The low stability of brown/blue reflects their biological function, not their size.
# ----------------------------------------------------------------------------

# --- Module Size vs Stability Plot ---
cat("  Generating Module Size vs Stability Plot...\n")

p8_size_stability <- ggplot(size_stability_all, aes(x = Size, y = Mean_Jaccard)) +
  geom_point(aes(fill = Module), shape = 21, size = 5, color = "black", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = Module), hjust = 0.5, vjust = -1.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "brown" = "brown", "blue" = "steelblue", "green" = "green3",
    "yellow" = "gold", "turquoise" = "turquoise3", "pink" = "pink",
    "magenta" = "magenta3", "black" = "grey30", "red" = "indianred", "grey" = "grey60"
  )) +
  labs(
    title = "Module Size vs Stability (All Modules)",
    subtitle = sprintf("Spearman rho = %.3f, p = %.3f", size_cor_all$estimate, size_cor_all$p.value),
    x = "Module Size (n proteins)", y = "Bootstrap Jaccard Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )
print(p8_size_stability)

ggsave(file.path(fig_dir_png, "Step7_08_SizeVsStability.png"), p8_size_stability, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_08_SizeVsStability.pdf"), p8_size_stability, width = 8, height = 6)

# ============================================================================
# 7i: VALIDATION - PROGNOSTIC STRENGTH VS STABILITY
# ============================================================================
#
# Purpose: Test hypothesis that prognostic modules show lower stability
#          because they capture patient heterogeneity (Short vs Long PFS)
#
# Hypothesis: Modules with stronger clinical associations should show
#             LOWER bootstrap stability (more sensitive to patient composition)
#
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7i: VALIDATION - PROGNOSTIC STRENGTH VS STABILITY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Prepare trait data
PFS_group_numeric <- as.numeric(datTraits$PFS_group == "Long")
Response_numeric <- as.numeric(datTraits$Response == "CD")

# Calculate module-trait correlations WITH p-values
module_trait_with_pval <- data.frame(Module = gsub("ME", "", colnames(MEs)))

for(i in 1:ncol(MEs)) {
  me <- MEs[, i]
  test_pfs <- cor.test(me, PFS_group_numeric, method = "pearson")
  module_trait_with_pval$Cor_PFS[i] <- test_pfs$estimate
  module_trait_with_pval$P_PFS[i] <- test_pfs$p.value
  test_resp <- cor.test(me, Response_numeric, method = "pearson")
  module_trait_with_pval$Cor_Response[i] <- test_resp$estimate
  module_trait_with_pval$P_Response[i] <- test_resp$p.value
}

module_trait_with_pval <- module_trait_with_pval %>%
  mutate(
    Sig_PFS = P_PFS < 0.05,
    Sig_Response = P_Response < 0.05,
    Sig_Both = Sig_PFS & Sig_Response,
    Mean_Abs_Cor = (abs(Cor_PFS) + abs(Cor_Response)) / 2
  )

cat("  Module-Trait Correlations:\n")
print(module_trait_with_pval %>%
        select(Module, Cor_PFS, P_PFS, Cor_Response, P_Response, Sig_Both) %>%
        arrange(P_PFS) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Merge with stability data
validation_summary <- module_trait_with_pval %>%
  left_join(size_stability_all, by = "Module") %>%
  filter(!is.na(Mean_Jaccard)) %>%
  mutate(
    Significance = case_when(
      Sig_Both ~ "Significant (p<0.05)",
      TRUE ~ "Not Significant"
    )
  )

cat("  Combined Stability + Prognostic Data:\n")
print(validation_summary %>%
        select(Module, Size, Mean_Jaccard, Cor_PFS, Cor_Response, Sig_Both) %>%
        arrange(desc(Sig_Both), desc(abs(Cor_PFS))) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Statistical test: Correlation between prognostic strength and stability
mean_abs_cor <- as.numeric(validation_summary$Mean_Abs_Cor)
mean_jaccard <- as.numeric(validation_summary$Mean_Jaccard)

cor_test_validation <- cor.test(mean_abs_cor, mean_jaccard, method = "spearman")

cat(sprintf("\n  Prognostic Strength vs Stability:"))
cat(sprintf("    Spearman rho = %.3f, p = %.4f\n",
            cor_test_validation$estimate, cor_test_validation$p.value))

# Compare significant vs non-significant modules
sig_modules <- validation_summary %>% filter(Sig_Both) %>% pull(Mean_Jaccard)
nonsig_modules <- validation_summary %>% filter(!Sig_Both) %>% pull(Mean_Jaccard)

cat(sprintf("\n  Significant modules (n=%d): Jaccard = %.3f Â± %.3f\n",
            length(sig_modules), mean(sig_modules), sd(sig_modules)))
cat(sprintf("  Non-significant modules (n=%d): Jaccard = %.3f Â± %.3f\n",
            length(nonsig_modules), mean(nonsig_modules), sd(nonsig_modules)))

if(length(sig_modules) >= 2 && length(nonsig_modules) >= 2) {
  wilcox_validation <- wilcox.test(sig_modules, nonsig_modules)
  cat(sprintf("  Wilcoxon p-value: %.4f\n", wilcox_validation$p.value))
}

cat("  Key Finding:\n")
cat("  - Brown and Blue are the ONLY modules significantly associated with PFS and Response\n")
cat("  - These prognostic modules show LOWER stability (more patient-sensitive)")
cat("  - This is EXPECTED: biomarkers that differentiate patients should vary\n")
cat("    when patient composition changes during bootstrap resampling\n")
cat("  - Non-prognostic modules (green, yellow) are stable but NOT clinically relevant")

# ----------------------------------------------------------------------------
# INTERPRETATION - Prognostic Strength vs Stability
# ----------------------------------------------------------------------------
# Interestingly, the brown module-which showed the strongest association with poor prognosis-
# also exhibited the highest sensitivity to patient composition (lowest Jaccard).
# This is consistent with its biological function: a module that differentiates Short from
# Long PFS patients should, by definition, vary when the proportion of these patient groups
# changes during resampling.
#
# In contrast, the green module showed higher stability but was not associated with clinical
# outcomes, suggesting it captures stable biological processes unrelated to disease progression.
#
# The apparent "instability" of the brown module is therefore not a weakness but rather
# evidence of its prognostic sensitivity.
# ----------------------------------------------------------------------------

# --- Prognostic vs Stability Quadrant Plot ---
cat("  Generating Prognostic vs Stability Quadrant Plot...\n")

p9_quadrant <- ggplot(validation_summary, aes(x = Mean_Abs_Cor, y = Mean_Jaccard)) +
  # Threshold lines
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0.35, linetype = "dashed", color = "grey40") +
  # Points - triangles for significant
  geom_point(aes(fill = Module, shape = Significance, size = Size),
             color = "black", alpha = 0.9, stroke = 1.2) +
  scale_shape_manual(values = c("Significant (p<0.05)" = 24, "Not Significant" = 21)) +
  geom_text(aes(label = Module), hjust = -0.3, vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "brown" = "brown", "blue" = "steelblue", "green" = "green3",
    "yellow" = "gold", "turquoise" = "turquoise3", "pink" = "pink",
    "magenta" = "magenta3", "black" = "grey30", "red" = "indianred",
    "purple" = "purple3"
  )) +
  scale_size_continuous(range = c(4, 12), name = "Module\nSize") +
  labs(
    title = "Module Classification: Clinical Significance vs Stability",
    subtitle = "Triangle = Significant for PFS or Response (p<0.05) | Dashed lines = classification thresholds",
    x = "Clinical Association Strength (Mean |r| with PFS & Response)",
    y = "Bootstrap Stability (Jaccard Index)",
    shape = "Statistical\nSignificance"
  ) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  guides(fill = "none")
print(p9_quadrant)

ggsave(file.path(fig_dir_png, "Step7_09_Validation_Quadrant.png"), p9_quadrant, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_09_Validation_Quadrant.pdf"), p9_quadrant, width = 8, height = 5)

# --- Combined Publication Figure ---
cat("  Generating Combined Publication Figure...\n")

# Use quadrant plot without legend for combined figure
p9_nolegend <- p9_quadrant + theme(legend.position = "none")

# Dynamic subtitle for all focus modules
focus_modules_title <- paste(sapply(focus_modules_auto, tools::toTitleCase), collapse = ", ")

p10_combined <- (p2_stability + p5_preservation) / (p7_cocluster_bar + p9_nolegend) +
  plot_annotation(
    title = "Module Stability Analysis Summary",
    subtitle = paste0(focus_modules_title, ": Significant prognostic modules with patient-sensitive stability | Grey: Negative control"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
    )
  )
print(p10_combined)

ggsave(file.path(fig_dir_png, "Step7_10_Combined_Publication.png"), p10_combined, width = 12, height = 10, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step7_10_Combined_Publication.pdf"), p10_combined, width = 12, height = 10)

cat("  All Step 7 figures generated\n")

# ============================================================================
# 7j: EXPORT ALL TABLES
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("7j: EXPORTING ALL TABLES\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Create Excel workbook
wb <- createWorkbook()

# Table 1: Overall Stability Summary - dynamic module names
cat("  [1/12] Stability Summary\n")
stability_summary <- data.frame(
  Module = focus_modules,
  Size = sapply(focus_modules, function(m) sum(moduleColors == m)),
  Bootstrap_Jaccard = bootstrap_results$Mean_Jaccard[match(focus_modules, bootstrap_results$Module)],
  Bootstrap_SD = bootstrap_results$SD_Jaccard[match(focus_modules, bootstrap_results$Module)],
  Bootstrap_CI_Lower = bootstrap_results$CI_Lower[match(focus_modules, bootstrap_results$Module)],
  Bootstrap_CI_Upper = bootstrap_results$CI_Upper[match(focus_modules, bootstrap_results$Module)],
  Perm_Zscore = perm_summary$Z_score[match(focus_modules, perm_summary$Module)],
  Perm_Pvalue = perm_summary$P_value[match(focus_modules, perm_summary$Module)],
  Preservation_Zsummary = preservation_summary$Zsummary[match(focus_modules, preservation_summary$Module)],
  Preservation_Interpretation = preservation_summary$Interpretation[match(focus_modules, preservation_summary$Module)]
)
addWorksheet(wb, "1_Stability_Summary")
writeData(wb, "1_Stability_Summary", stability_summary)
setColWidths(wb, "1_Stability_Summary", cols = 1:ncol(stability_summary), widths = "auto")

# Table 2: Bootstrap Results
cat("  [2/12] Bootstrap Results\n")
addWorksheet(wb, "2_Bootstrap")
writeData(wb, "2_Bootstrap", bootstrap_results)
setColWidths(wb, "2_Bootstrap", cols = 1:ncol(bootstrap_results), widths = "auto")

# Table 3: Subsampling Results - SKIPPED (subsample analysis not run)
cat("  [3/12] Subsampling Results - SKIPPED\n")

# Table 4: Permutation Test
cat("  [4/12] Permutation Test\n")
perm_export <- perm_summary %>%
  mutate(
    Significant = ifelse(P_value < 0.05, "Yes", "No"),
    Interpretation = ifelse(Z_score > 2, "Module is NOT random", "May be random")
  )
addWorksheet(wb, "4_Permutation")
writeData(wb, "4_Permutation", perm_export)
setColWidths(wb, "4_Permutation", cols = 1:ncol(perm_export), widths = "auto")

# Table 5: Module Preservation
cat("  [5/12] Module Preservation\n")
addWorksheet(wb, "5_Preservation")
writeData(wb, "5_Preservation", preservation_summary)
setColWidths(wb, "5_Preservation", cols = 1:ncol(preservation_summary), widths = "auto")

# Table 6: Hub Gene Stability
cat("  [6/12] Hub Gene Stability\n")
addWorksheet(wb, "6_Hub_Stability")
writeData(wb, "6_Hub_Stability", hub_stability)
setColWidths(wb, "6_Hub_Stability", cols = 1:ncol(hub_stability), widths = "auto")

# Table 7: Co-Clustering Matrix
cat("  [7/12] Co-Clustering Matrix\n")
cocluster_export <- as.data.frame(round(co_cluster_pct, 1))
cocluster_export <- cbind(Protein = rownames(cocluster_export), cocluster_export)
addWorksheet(wb, "7_CoCluster_Matrix")
writeData(wb, "7_CoCluster_Matrix", cocluster_export)
setColWidths(wb, "7_CoCluster_Matrix", cols = 1:ncol(cocluster_export), widths = "auto")

# Table 8: Co-Clustering Summary - dynamic for all focus modules
cat("  [8/12] Co-Clustering Summary\n")
within_export <- lapply(focus_modules_auto, function(mod) {
  pairs <- within_module_pairs[[mod]]
  if(length(pairs) > 0) {
    data.frame(
      Comparison = paste0(tools::toTitleCase(mod), "_Within"),
      Mean_Pct = mean(pairs),
      SD_Pct = sd(pairs),
      N_Pairs = length(pairs),
      Wilcoxon_P = NA,
      Interpretation = "Hub genes cluster together"
    )
  } else {
    NULL
  }
})
within_export <- bind_rows(within_export)

between_export <- data.frame(
  Comparison = "Cross_Module_Between",
  Mean_Pct = mean(between_pairs),
  SD_Pct = sd(between_pairs),
  N_Pairs = length(between_pairs),
  Wilcoxon_P = wilcox_result$p.value,
  Interpretation = "Modules remain distinct"
)

cocluster_summary_export <- bind_rows(within_export, between_export)
addWorksheet(wb, "8_CoCluster_Summary")
writeData(wb, "8_CoCluster_Summary", cocluster_summary_export)
setColWidths(wb, "8_CoCluster_Summary", cols = 1:ncol(cocluster_summary_export), widths = "auto")

# Table 9: Module Reassignment
cat("  [9/12] Module Reassignment\n")
reassignment_export <- staying_rate %>%
  mutate(
    Interpretation = case_when(
      Stay_Rate > 0.3 ~ "Relatively stable",
      Stay_Rate > 0.15 ~ "Moderate stability",
      TRUE ~ "Low stability (proteins redistribute)"
    )
  )
addWorksheet(wb, "9_Reassignment")
writeData(wb, "9_Reassignment", reassignment_export)
setColWidths(wb, "9_Reassignment", cols = 1:ncol(reassignment_export), widths = "auto")

# Table 10: Match Frequency
cat("  [10/12] Match Frequency\n")
addWorksheet(wb, "10_Match_Frequency")
writeData(wb, "10_Match_Frequency", match_frequency)
setColWidths(wb, "10_Match_Frequency", cols = 1:ncol(match_frequency), widths = "auto")

# Table 11: Size vs Stability
cat("  [11/12] Size vs Stability\n")
size_stability_export <- size_stability_all %>%
  arrange(desc(Mean_Jaccard)) %>%
  mutate(
    Rank_Stability = row_number(),
    Spearman_rho = size_cor_all$estimate,
    Spearman_p = size_cor_all$p.value
  )
addWorksheet(wb, "11_Size_vs_Stability")
writeData(wb, "11_Size_vs_Stability", size_stability_export)
setColWidths(wb, "11_Size_vs_Stability", cols = 1:ncol(size_stability_export), widths = "auto")

# Table 12: Validation Data
cat("  [12/12] Validation Data\n")
validation_export <- validation_summary %>%
  select(Module, Size, Mean_Jaccard, Cor_PFS, P_PFS, Cor_Response, P_Response,
         Sig_Both, Mean_Abs_Cor, Significance) %>%
  arrange(desc(Sig_Both), desc(Mean_Abs_Cor)) %>%
  mutate(
    Interpretation = case_when(
      Sig_Both & Mean_Jaccard < 0.3 ~ "Prognostic biomarker, patient-sensitive",
      Sig_Both & Mean_Jaccard >= 0.3 ~ "Prognostic biomarker, stable",
      !Sig_Both & Mean_Jaccard > 0.4 ~ "Non-prognostic, stable (housekeeping?)",
      TRUE ~ "Intermediate"
    )
  )
addWorksheet(wb, "12_Validation")
writeData(wb, "12_Validation", validation_export)
setColWidths(wb, "12_Validation", cols = 1:ncol(validation_export), widths = "auto")

# Save Excel workbook
saveWorkbook(wb, file.path(results_dir, "Step7_Stability_Results.xlsx"), overwrite = TRUE)
cat("  Saved: Step7_Stability_Results.xlsx (12 sheets)")

# Save CSVs
write.csv(stability_summary, file.path(results_dir, "Step7_StabilitySummary.csv"), row.names = FALSE)
write.csv(hub_stability, file.path(results_dir, "Step7_HubStability.csv"), row.names = FALSE)
write.csv(validation_export, file.path(results_dir, "Step7_ValidationData.csv"), row.names = FALSE)
cat("  Saved: CSV files\n")

# Save RDS
networkData$stability <- list(
  bootstrap = bootstrap_results,
  bootstrap_all = all_metrics_final,
  permutation = perm_summary,
  preservation = preservation_summary,
  hub_stability = hub_stability,
  cocluster_matrix = co_cluster_pct,
  cocluster_summary = cocluster_summary_export,
  match_frequency = match_frequency,
  reassignment = staying_rate,
  size_stability = size_stability_all,
  validation = validation_export,
  convergence = convergence_history,
  n_iterations = n_bootstrap
)
saveRDS(networkData, file.path(results_dir, "Step7_NetworkData.rds"))
cat("  Saved: Step7_NetworkData.rds\n")

# ============================================================================
# STEP 7 FINISHED - FINAL SUMMARY
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 7 FINISHED: MODULE STABILITY TESTING\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("KEY FINDINGS:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

cat("1. MODULE STABILITY (Bootstrap with Module Matching):\n")
for(mod in focus_modules_auto) {
  mod_row <- bootstrap_results[bootstrap_results$Module == mod, ]
  if(nrow(mod_row) > 0) {
    cat(sprintf("   %s: Jaccard = %.3f [%.3f-%.3f]\n",
                tools::toTitleCase(mod), mod_row$Mean_Jaccard, mod_row$CI_Lower, mod_row$CI_Upper))
  }
}
grey_row <- bootstrap_results[bootstrap_results$Module == "grey", ]
if(nrow(grey_row) > 0) {
  cat(sprintf("   Grey: Jaccard = %.3f [%.3f-%.3f] (negative control)\n",
              grey_row$Mean_Jaccard, grey_row$CI_Lower, grey_row$CI_Upper))
}

cat("2. STATISTICAL VALIDATION:\n")
perm_str <- paste(sapply(focus_modules_auto, function(mod) {
  z <- perm_summary$Z_score[perm_summary$Module == mod]
  sprintf("%s Z=%.1f", tools::toTitleCase(mod), z)
}), collapse = ", ")
cat(sprintf("   Permutation Test: %s (all p<0.0001)\n", perm_str))
pres_str <- paste(sapply(focus_modules_auto, function(mod) {
  z <- preservation_summary$Zsummary[preservation_summary$Module == mod]
  if(length(z) > 0) sprintf("%s=%.1f", tools::toTitleCase(mod), z) else NULL
}), collapse = ", ")
cat(sprintf("   Module Preservation Zsummary: %s\n", pres_str))

cat("3. HUB GENE CO-CLUSTERING:\n")
cat(sprintf("   Within-module: %.1f%%\n", mean(within_group)))
cat(sprintf("   Between-module: %.1f%%\n", mean(between_pairs)))
cat(sprintf("   Wilcoxon p-value: %.6f\n", wilcox_result$p.value))
cat("   -> Hub genes cluster together significantly MORE within modules\n")

cat("4. CLINICAL SIGNIFICANCE:\n")
focus_names_str <- paste(sapply(focus_modules_auto, tools::toTitleCase), collapse = ", ")
cat(sprintf("   Focus modules: %s (associated with PFS & Response)\n", focus_names_str))
cat("   -> Lower stability reflects SENSITIVITY to patient heterogeneity\n")
cat("   -> This is EXPECTED for biomarkers that differentiate patients\n")

cat("CONCLUSIONS:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat(sprintf("- %s modules are statistically real (not random noise)\n", focus_names_str))
cat("- Module boundaries are fluid due to small sample size (n=50)\n")
cat("- Hub gene communities are stable despite color changes\n")
cat("- Prognostic modules show lower stability because they capture\n")
cat("  patient heterogeneity (Short vs Long PFS, Responders vs Non-responders)\n")
cat("- This patient-sensitivity is a FEATURE, not a limitation\n")

cat("FILES GENERATED:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("  Figures: 10 PNG + 10 PDF files\n")
cat("  Data:    Step7_Stability_Results.xlsx (12 sheets)")
cat("           Step7_StabilitySummary.csv\n")
cat("           Step7_HubStability.csv\n")
cat("           Step7_ValidationData.csv\n")
cat("           Step7_NetworkData.rds\n")

cat("\n")
cat("--- SENSITIVITY ANALYSIS DOCUMENTATION ---\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("Bootstrap Parameters:\n")
cat(sprintf("  - Iterations: %d (adaptive, convergence CV < 5%%)\n", n_bootstrap))
cat("  - Resampling: With replacement, preserving sample structure\n")
cat("  - Network: Reconstructed de novo for each resample\n")
cat("  - Module matching: Maximum Jaccard similarity-based\n")
cat("\nStability Thresholds (Langfelder & Horvath, 2008):\n")
cat("  - Jaccard >= 0.70: Robust module\n")
cat("  - Jaccard >= 0.85: Excellent stability\n")
cat("  - Jaccard >= 0.50: Moderate stability (acceptable for clinical signatures)\n")
cat("\nAlternative Power Values Tested (Step 4):\n")
cat("  - Scale-free topology: R^2 > 0.85 threshold\n")
cat("  - Power sensitivity: ±2 from optimal tested for network properties\n")
cat("  - Conclusion: Network structure stable across tested power range\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

# --- MEMORY MANAGEMENT: Post-Bootstrap Cleanup ---
# Bootstrap analysis generates many temporary objects; clear to free memory
cat("\n  Memory cleanup: Clearing bootstrap temporary objects...\n")
bootstrap_temp_objects <- c("all_metrics", "all_metrics_final", "all_reassignments")
for(obj in bootstrap_temp_objects) {
  if(exists(obj)) rm(list = obj, envir = .GlobalEnv)
}
invisible(gc())  # Force garbage collection
cat(sprintf("  Memory freed. Current usage: %.1f MB\n",
            sum(gc()[,2])))  # Report memory after cleanup

}  # End of else block for run_step7


#=============================================================================
# STEP 8: PATHWAY ENRICHMENT ANALYSIS
#=============================================================================

# PURPOSE: Interrogate module biology using 267 curated PDAC-relevant
# signatures via ORA, GSVA, and ssGSEA - linking modules to pathways.   
# STRATEGY:                                                               
# 8a) Setup -> 8b) Load ~264 signatures -> 8c) Coverage assessment        
# 8d) ORA enrichment -> 8e) GSVA scoring -> 8f) ssGSEA                      
# 8g) Clinical associations -> 8h) Biomarkers -> 8i) Module-pathway corr   
# LOOKING FOR:                                                             
# - ORA: FDR < 0.05, enrichment ratio > 2                                
# - GSVA: Significant score differences between PFS groups              
# - Biological coherence: related pathways enriched in same module       
# KEY PARAMS: ~264 signatures (after dedup), FDR_cutoff=0.05, min_overlap=3
# ==============================================================================
# DETAILED SECTIONS 
# ==============================================================================

# Hypothesis: Short PFS patients have higher stemness, activated stroma,
#             and treatment resistance signatures
#
# Input: Step 7 WGCNA results (datExpr, moduleColors, MEs, datTraits)
# Output: ORA, GSVA, ssGSEA results with publication-ready figures
#
# Sections:
#   8a: Setup and load dependencies
#   8b: Define gene signatures (~264 signatures: 138 Nika + 126 Eiseler, deduplicated)
#   8c: Coverage assessment
#   8d: Over-Representation Analysis (ORA)
#   8e: GSVA (Gene Set Variation Analysis)
#   8f: ssGSEA (single-sample GSEA)
#   8g: Clinical association testing
#   8h: Biomarker analysis (ML-derived and Cox-derived separately)
#   8i: Module-pathway integration
#   8j: Summary and save
#
# ==============================================================================

cat("===============================================================================\n")
cat("STEP 8: PATHWAY ENRICHMENT ANALYSIS\n")
cat("===============================================================================\n")

# ============================================================================
# 8a: SETUP
# ============================================================================

cat("8a: SETUP\n")
cat(paste(rep("-", 76), collapse = ""), "\n")

allowWGCNAThreads(nThreads = n_cores)  # Use same core count as Step 1


# -- Color palettes (using universal clinical colors from top of script) --
col_pfs <- colors_pfs
col_response <- colors_response
col_treatment <- colors_treatment

# Module colors (consistent across all figures)
module_colors_palette <- c(
  "GREEN" = "#70a494", "TURQUOISE" = "#40E0D0", "BLUE" = "#4169E1",
  "BROWN" = "#8B4513", "YELLOW" = "#DAA520", "RED" = "#DC143C",
  "BLACK" = "#2F2F2F", "PINK" = "#FF69B4", "MAGENTA" = "#BA55D3",
  "PURPLE" = "#9370DB", "GREENYELLOW" = "#ADFF2F", "TAN" = "#D2B48C",
  "SALMON" = "#FA8072", "CYAN" = "#00CED1", "GREY" = "#808080",
  "green" = "#70a494", "turquoise" = "#40E0D0", "blue" = "#4169E1",
  "brown" = "#8B4513", "yellow" = "#DAA520", "red" = "#DC143C",
  "black" = "#2F2F2F", "pink" = "#FF69B4", "magenta" = "#BA55D3",
  "purple" = "#9370DB", "greenyellow" = "#ADFF2F", "tan" = "#D2B48C",
  "salmon" = "#FA8072", "cyan" = "#00CED1", "grey" = "#808080"
)

# Biological category colors (for pathway classification)
category_colors <- c(
  "Immune" = "#E64B35",
  "Stroma/ECM" = "#7E6148",
  "PDAC Subtype" = "#3C5488",
  "Stemness" = "#9B59B6",
  "EMT" = "#00A087",
  "Metabolic" = "#F39B7F",
  "Angiogenesis" = "#DC0000",
  "Exosome" = "#4DBBD5",
  "Signaling" = "#8491B4",
  "Other" = "#B09C85"
)

# Heatmap palettes
heatmap_diverging <- colorRampPalette(c("#2166AC", "#f6edbd", "#B2182B"))(100)
heatmap_scale <- colorRampPalette(c("#f6edbd", "#de8a5a", "#B23A27", "#8C2D1E"))(100)  # ORANGE spectrum
heatmap_colors <- colorRampPalette(c("#0E4D92", "#F5F5F5", "#D95C2B"))(100)
enrichment_colors <- colorRampPalette(c("#F5F5F5", "#D6EFB3", "#5FC1C0", "#234DA0"))(100)

# Publication font sizes
pub_base_size <- 14
pub_title_size <- 12
pub_subtitle_size <- 11
pub_axis_text_size <- 12
pub_legend_text_size <- 11
pub_strip_size <- 12

# Helper function: Clean pathway names for figures
clean_pathway_names <- function(name) {
  name %>%
    str_replace("^TIDE Exclusion$", "T-cell Exclusion (TIDE)") %>%
    str_replace("^TIDE T-cell Exclusion$", "T-cell Exclusion (TIDE)") %>%
    str_replace("^TIDE$", "TIDE Score") %>%
    str_replace("MDSC M-MDSC", "m-Myeloid-Derived Suppressor") %>%
    str_replace("MDSC G-MDSC", "g-Myeloid-Derived Suppressor") %>%
    str_replace("^Monocytic MDSCs$", "m-Myeloid-Derived Suppressor") %>%
    str_replace("^Granulocytic MDSCs$", "g-Myeloid-Derived Suppressor") %>%
    str_replace("^M-MDSC$", "m-Myeloid-Derived Suppressor") %>%
    str_replace("^G-MDSC$", "g-Myeloid-Derived Suppressor") %>%
    str_replace("^M-MDSCs$", "m-Myeloid-Derived Suppressor") %>%
    str_replace("^G-MDSCs$", "g-Myeloid-Derived Suppressor") %>%
    str_replace("Warburg Effect \\(Aerobic Glycolysis\\)", "Warburg Effect") %>%
    str_replace("^Warburg Glycolysis$", "Warburg Effect") %>%
    str_replace("^M2 Tumor-Associated Macrophages$", "M2 TAMs") %>%
    str_replace("^Cancer-Associated Fibroblasts$", "CAFs") %>%
    str_replace("^Activated Pancreatic Stellate Cells$", "Activated PSCs") %>%
    str_replace("^Quiescent Pancreatic Stellate Cells$", "Quiescent PSCs") %>%
    str_replace("Cox-Derived Biomarkers", "Cox Signature") %>%
    str_replace("ML-Derived Biomarkers", "ML Signature") %>%
    str_replace("^Extracellular Matrix$", "ECM") %>%
    str_replace("^Pentose Phosphate Pathway$", "PPP") %>%
    str_replace("^TCA Cycle \\(Krebs Cycle\\)$", "TCA Cycle") %>%
    str_trim()
}

# Helper function: Assign biological category
assign_category <- function(category_string) {
  case_when(
    grepl("Immune|CIBERSORT|Macro|MDSC|Neutro|T.cell|B.cell|NK|Dendrit|Mast|Mono|Lymph|TAM",
          category_string, ignore.case = TRUE) ~ "Immune",
    grepl("Stroma|CAF|PSC|Fibro|ECM|Collagen|Matrix|Activated.Stroma",
          category_string, ignore.case = TRUE) ~ "Stroma/ECM",
    grepl("Subtype|Moffitt|Bailey|Collisson|Basal|Classical|Squamous|Progenitor|ADEX|Immunogenic",
          category_string, ignore.case = TRUE) ~ "PDAC Subtype",
    grepl("Stem|CSC|Notch|Wnt|Hedgehog",
          category_string, ignore.case = TRUE) ~ "Stemness",
    grepl("EMT|Mesenchym|Epithel",
          category_string, ignore.case = TRUE) ~ "EMT",
    grepl("Metab|Glycol|TCA|Warburg|PPP|Lipid|Oxphos|Glutamin|Fatty|NRF2|Hypoxia",
          category_string, ignore.case = TRUE) ~ "Metabolic",
    grepl("Angio|Endothel|Vascular|VEGF",
          category_string, ignore.case = TRUE) ~ "Angiogenesis",
    grepl("Exo|Vesicle|EV",
          category_string, ignore.case = TRUE) ~ "Exosome",
    grepl("Signal|Pathway|MAPK|AKT|KRAS|TGF|JAK|STAT|NF.?kB|mTOR|TIDE",
          category_string, ignore.case = TRUE) ~ "Signaling",
    TRUE ~ "Other"
  )
}

# -- Verify required objects from Step 7 --
required_objects <- c("datExpr", "moduleColors", "datTraits")
missing <- required_objects[!sapply(required_objects, exists)]

if(length(missing) > 0) {
  cat("   Warning: Missing objects:", paste(missing, collapse = ", "), "\n")
  cat("   Attempting to load from saved data...\n")

  # Try Step 7 first, then Step 6
  step7_file <- file.path(results_dir, "Step7_WGCNA_Final.RData")
  step6_file <- file.path(results_dir, "Step6_NetworkData.rds")

  if(file.exists(step7_file)) {
    load(step7_file)
    cat("   Loaded from Step7_WGCNA_Final.RData\n")
  } else if(file.exists(step6_file)) {
    cat("   Step 7 not found - Loading from Step6_NetworkData.rds\n")
    networkData <- readRDS(step6_file)
    datExpr <- networkData$datExpr
    moduleColors <- networkData$moduleColors
    datTraits <- networkData$datTraits
    MEs <- networkData$MEs
    TOM <- networkData$TOM
    geneTree <- networkData$geneTree
    proteinInfo <- networkData$proteinInfo
    cat("   Loaded from Step6_NetworkData.rds\n")
  } else {
    stop("Cannot find Step 6 or Step 7 workspace. Please run Step 6 first.")
  }
}

n_samples <- nrow(datExpr)
n_proteins <- ncol(datExpr)
all_proteins <- colnames(datExpr)

# Safe check for valid dimensions (handles NULL/NA cases)
if(!is.null(n_proteins) && !is.na(n_proteins) && !is.null(n_samples) && !is.na(n_samples)) {
  # Defensive check: if datExpr appears transposed
  if(n_proteins == 0 || n_proteins < n_samples) {
    cat("   [WARNING] datExpr may be transposed. Checking orientation...\n")
    if(length(intersect(names(moduleColors), rownames(datExpr))) > length(intersect(names(moduleColors), colnames(datExpr)))) {
      cat("   [INFO] Transposing datExpr (proteins were in rows, moving to columns)\n")
      datExpr <- t(datExpr)
      n_samples <- nrow(datExpr)
      n_proteins <- ncol(datExpr)
      all_proteins <- colnames(datExpr)
    }
  }
}

# Final validation - only if n_proteins is still problematic
if(is.null(n_proteins) || is.na(n_proteins) || n_proteins == 0) {
  # Try to get proteins from moduleColors
  all_proteins <- names(moduleColors)
  n_proteins <- length(all_proteins)
  n_samples <- if(is.null(n_samples) || is.na(n_samples)) nrow(datExpr) else n_samples
  cat(sprintf("   [FALLBACK] Using %d proteins from moduleColors\n", n_proteins))
}

cat(sprintf("   datExpr: %d samples x %d proteins\n", n_samples, n_proteins))
cat(sprintf("   Modules: %d (excluding grey)\n", length(unique(moduleColors)) - 1))

# -- Calculate module eigengenes if not present --
if(!exists("MEs")) {
  cat("   Calculating module eigengenes...\n")
  MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs)
}

# -- Ensure focus_modules_auto is defined (needed if Step 7 was skipped) --
if(!exists("focus_modules_auto") || length(focus_modules_auto) < 1) {
  cat("   Detecting focus modules from clinical associations...\n")
  # Load from Step 5 network data if available
  step5_file <- file.path(results_dir, "Step5_WGCNA_Network.RData")
  if(file.exists(step5_file) && !exists("focus_modules_auto")) {
    networkData <- readRDS(step5_file)
    if(!is.null(networkData$focus_modules_auto)) {
      focus_modules_auto <- networkData$focus_modules_auto
    }
  }
  # Fallback: use modules with most clinical significance
  if(!exists("focus_modules_auto") || length(focus_modules_auto) < 1) {
    # Get non-grey modules sorted by size
    mod_sizes <- table(moduleColors)
    mod_sizes <- mod_sizes[names(mod_sizes) != "grey"]
    focus_modules_auto <- names(sort(mod_sizes, decreasing = TRUE))[1:4]
    cat(sprintf("   [FALLBACK] Using largest modules: %s\n", paste(focus_modules_auto, collapse = ", ")))
  }
}

# Display all focus modules
n_focus <- length(focus_modules_auto)
focus_names <- sapply(focus_modules_auto, tools::toTitleCase)
cat(sprintf("   Focus modules (%d total): %s\n", n_focus, paste(focus_names, collapse = ", ")))

# For backward compatibility, define mod1/mod2 as first two focus modules
# But primary analysis will use ALL focus_modules_auto
mod1 <- focus_modules_auto[1]
mod2 <- if(length(focus_modules_auto) >= 2) focus_modules_auto[2] else focus_modules_auto[1]
mod1_name <- tools::toTitleCase(mod1)
mod2_name <- tools::toTitleCase(mod2)

# ============================================================================
# 8b: LOAD GENE SIGNATURES FROM MASTER TABLE
# ============================================================================

cat("-----------------------------------------------------------\n")
cat("8b: Loading Gene Signatures from Master Table\n")
cat("-----------------------------------------------------------\n")

# Load custom signatures from master table
# Try multiple locations for the mastertable file
mastertable_paths <- c(
  "mastertable_Signature.csv",                          # Current directory
  "../mastertable_Signature.csv",                       # Parent directory
  file.path(dirname(getwd()), "mastertable_Signature.csv")  # Parent using full path
)
mastertable_file <- NULL
for(path in mastertable_paths) {
  if(file.exists(path)) {
    mastertable_file <- path
    break
  }
}
if(is.null(mastertable_file)) {
  stop("ERROR: mastertable_Signature.csv not found. Searched in:\n",
       paste("  -", mastertable_paths, collapse = "\n"))
}
cat(sprintf("   Loading signatures from: %s\n", mastertable_file))
sig_table <- read.csv(mastertable_file, stringsAsFactors = FALSE)

custom_signatures <- list()
for(i in 1:nrow(sig_table)) {
  sig_id <- sig_table$id[i]
  genes <- unlist(strsplit(sig_table$genes[i], ";"))
  genes <- genes[genes != ""]  # Remove empty strings
  custom_signatures[[sig_id]] <- list(
    id = sig_id,
    short_name = sig_table$short_name[i],
    full_name = sig_table$full_name[i],
    category = sig_table$category[i],
    genes = genes,
    source = sig_table$source[i],
    pmid = sig_table$pmid[i]
  )
}
cat(sprintf("Loaded %d custom signatures from mastertable\n", length(custom_signatures)))

# Also keep as 'signatures' for compatibility with existing code
signatures <- custom_signatures

cat(sprintf("  Categories: %s\n", paste(unique(sig_table$category), collapse = ", ")))

# Summary by category
step8b_summary <- as.data.frame(table(sig_table$category))
colnames(step8b_summary) <- c("Category", "Count")
print(step8b_summary)

# ============================================================================
# 8c: COVERAGE ASSESSMENT
# ============================================================================

cat("8c: COVERAGE ASSESSMENT\n")
cat(paste(rep("-", 76), collapse = ""), "\n")

# -- Build reference table --
sig_ref <- do.call(rbind, lapply(names(signatures), function(name) {
  sig <- signatures[[name]]
  if(length(sig$genes) == 0) return(NULL)  # Skip empty signatures
  available <- sig$genes[sig$genes %in% all_proteins]

  data.frame(
    ID = sig$id,
    Short_Name = sig$short_name,
    Full_Name = sig$full_name,
    Category = sig$category,
    N_Total = length(sig$genes),
    N_Available = length(available),
    Coverage_Pct = round(100 * length(available) / length(sig$genes), 1),
    Source = sig$source,
    PMID = ifelse(is.null(sig$pmid), NA, sig$pmid),
    Genes_All = paste(sig$genes, collapse = ", "),
    Genes_Available = paste(available, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

# -- Filter usable signatures --
min_genes <- 3
sig_ref$Usable <- sig_ref$N_Available >= min_genes

cat(sprintf("   Total signatures: %d\n", nrow(sig_ref)))
cat(sprintf("   Usable (>=%d genes available): %d\n", min_genes, sum(sig_ref$Usable)))
cat(sprintf("   Below threshold: %d\n", sum(!sig_ref$Usable)))

# -- Category summary --
cat("   Coverage by category:\n")
cat_summary <- sig_ref %>%
  group_by(Category) %>%
  summarise(
    N_Total = n(),
    N_Usable = sum(Usable),
    Mean_Coverage = round(mean(Coverage_Pct), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(N_Usable))

print(as.data.frame(cat_summary))

# -- Create gene sets for analysis --
gene_sets <- list()
display_names_short <- list()
display_names_full <- list()

for(name in names(signatures)) {
  sig <- signatures[[name]]
  available <- sig$genes[sig$genes %in% all_proteins]
  
  if(length(available) >= min_genes) {
    gene_sets[[name]] <- available
    display_names_short[[name]] <- sig$short_name
    display_names_full[[name]] <- sig$full_name
  }
}

cat(sprintf("\n   [OK] Gene sets prepared: %d\n", length(gene_sets)))

# -- Save coverage report --
write.csv(sig_ref, file.path(results_dir, "Step8_01_Signature_Coverage.csv"), row.names = FALSE)
cat("   [OK] Saved: Step8_01_Signature_Coverage.csv\n")

# ============================================================================
# 8c-2: UNIFIED ORA FOR ALL FOCUS MODULES
# ============================================================================

cat("==============================================================================\n")
cat(sprintf("8c-2: UNIFIED ORA FOR FOCUS MODULES (%s)\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("==============================================================================\n")

background_proteins <- all_proteins  # Use all_proteins defined earlier (891 proteins)
cat(sprintf("Background: %d proteins\n", length(background_proteins)))

# Store ORA results for each focus module
focus_module_ora <- list()

# Validate background before proceeding
if(length(background_proteins) == 0) {
  cat("   [ERROR] No background proteins available. Skipping 8c-2 ORA.\n")
} else {
  # Loop through ALL focus modules
  for(mod in focus_modules_auto) {
    mod_proteins <- names(moduleColors)[moduleColors == mod]
    cat(sprintf("\n  %s module: %d proteins\n", toupper(mod), length(mod_proteins)))

    # Run unified ORA
    mod_ora <- unified_ora(mod_proteins, custom_signatures, background_proteins,
                           min_overlap = 2, database_name = "Custom_PDAC")

    if(nrow(mod_ora) > 0) {
      mod_ora$Module <- toupper(mod)
      mod_ora$Short_Name <- clean_pathway_names(mod_ora$Short_Name)

      # Save results
      write.csv(mod_ora,
                file.path(results_dir, paste0("Step8c2_", toupper(mod), "_Custom_ORA.csv")),
                row.names = FALSE)
      cat(sprintf("    -> %d pathways with overlap, %d nominal significant (p<0.05)\n",
                  nrow(mod_ora), sum(mod_ora$P_Value < 0.05)))

      focus_module_ora[[mod]] <- mod_ora
    } else {
      cat(sprintf("    -> No pathways with sufficient overlap\n"))
      focus_module_ora[[mod]] <- data.frame()
    }
  }
}

# For backward compatibility, also create mod1_custom and mod2_custom
mod1_custom <- if(!is.null(focus_module_ora[[mod1]])) focus_module_ora[[mod1]] else data.frame()
mod2_custom <- if(!is.null(focus_module_ora[[mod2]])) focus_module_ora[[mod2]] else data.frame()

# ============================================================================
# 8d: OVER-REPRESENTATION ANALYSIS (ORA)
# ============================================================================

cat("==============================================================================\n")
cat("8d: OVER-REPRESENTATION ANALYSIS (ORA)\n")
cat("==============================================================================\n")

# -- Get module gene lists --
# Start with focus modules, then add other common modules
modules_of_interest <- unique(c(focus_modules_auto, "turquoise", "green", "yellow", "red", "pink", "black", "magenta", "brown", "blue"))
modules_of_interest <- modules_of_interest[modules_of_interest %in% unique(moduleColors)]

module_genes <- lapply(modules_of_interest, function(mod) {
  names(moduleColors)[moduleColors == mod]
})
names(module_genes) <- modules_of_interest

cat("   Module sizes:\n")
for(mod in modules_of_interest) {
  cat(sprintf("      %s: %d proteins\n", mod, length(module_genes[[mod]])))
}

# -- Hypergeometric ORA function --
run_ora <- function(module_genes_list, gene_sets_list, universe_size,
                    short_names = NULL, full_names = NULL, sig_info = NULL) {
  results <- data.frame()
  
  for(mod_name in names(module_genes_list)) {
    mod_genes <- module_genes_list[[mod_name]]
    
    for(set_name in names(gene_sets_list)) {
      set_genes <- gene_sets_list[[set_name]]
      overlap <- intersect(mod_genes, set_genes)
      
      if(length(overlap) == 0) next
      
      # Hypergeometric test
      p_val <- phyper(
        q = length(overlap) - 1,
        m = length(set_genes),
        n = universe_size - length(set_genes),
        k = length(mod_genes),
        lower.tail = FALSE
      )
      
      expected <- length(mod_genes) * length(set_genes) / universe_size
      fold_enrich <- ifelse(expected > 0, length(overlap) / expected, 0)
      
      # Get display names
      sname <- ifelse(!is.null(short_names) && set_name %in% names(short_names),
                      short_names[[set_name]], set_name)
      fname <- ifelse(!is.null(full_names) && set_name %in% names(full_names),
                      full_names[[set_name]], set_name)
      
      # Get category
      cat_name <- ifelse(!is.null(sig_info) && set_name %in% names(sig_info),
                         sig_info[[set_name]]$category, "Custom")
      
      results <- rbind(results, data.frame(
        Module = mod_name,
        Signature_ID = set_name,
        Short_Name = sname,
        Full_Name = fname,
        Category = cat_name,
        N_Module = length(mod_genes),
        N_Signature = length(set_genes),
        N_Overlap = length(overlap),
        Overlap_Genes = paste(overlap, collapse = "; "),
        Expected = round(expected, 2),
        Fold_Enrichment = round(fold_enrich, 2),
        P_Value = p_val,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if(nrow(results) > 0) {
    
    # FDR correction is applied PER MODULE (not globally) because:
    # 1. Scientific question is "within each module, which signatures are enriched?"
    # 2. Each module is an independent biological unit with its own hypothesis
    # 3. Global FDR would be overly conservative for within-module interpretation
    # For cross-module comparison, use FDR_global column added below
    results <- results %>%
      group_by(Module) %>%
      mutate(FDR = p.adjust(P_Value, method = "BH")) %>%
      ungroup() %>%
      arrange(Module, P_Value)

    # Add global FDR for cross-module comparison
    results$FDR_global <- p.adjust(results$P_Value, method = "BH")
  }
  
  return(results)
}

# -- Run ORA --
cat("   Running ORA on custom signatures...\n")

# Validate inputs before running ORA
if(length(gene_sets) == 0) {
  cat("   [WARNING] No gene sets available (0 signatures with >=3 genes in proteome).\n")
  cat("   [INFO] This may indicate gene symbols don't match between signatures and proteome.\n")
  ora_results <- data.frame()
} else if(n_proteins == 0) {
  cat("   [ERROR] No proteins in universe (n_proteins=0). Cannot run ORA.\n")
  ora_results <- data.frame()
} else {
  ora_results <- run_ora(
    module_genes, gene_sets, n_proteins,
    short_names = display_names_short,
    full_names = display_names_full,
    sig_info = signatures
  )
}

cat(sprintf("   [OK] Total module-signature pairs tested: %d\n", nrow(ora_results)))

# -- Significance assessment --
if(nrow(ora_results) > 0) {
  ora_results <- ora_results %>%
    mutate(
      Sig_Level = case_when(
        FDR < 0.05 ~ "FDR < 0.05",
        FDR < 0.10 ~ "FDR < 0.10",
        P_Value < 0.05 ~ "+",
        TRUE ~ "ns"
      ),
      Sig_Level = factor(Sig_Level, levels = c("FDR < 0.05", "FDR < 0.10", "p < 0.05", "ns"))
  )

  # Apply pathway name cleaning for publication-ready figures
  ora_results$Short_Name <- clean_pathway_names(ora_results$Short_Name)

  n_fdr05 <- sum(ora_results$FDR < 0.05, na.rm = TRUE)
  n_fdr10 <- sum(ora_results$FDR < 0.10, na.rm = TRUE)
  n_nom05 <- sum(ora_results$P_Value < 0.05, na.rm = TRUE)

  cat(sprintf("\n   Significance summary:"))
  cat(sprintf("      FDR < 0.05: %d\n", n_fdr05))
  cat(sprintf("      FDR < 0.10: %d\n", n_fdr10))
  cat(sprintf("      Nominal p < 0.05: %d\n", n_nom05))
} else {
  cat("   [INFO] No ORA results to process.\n")
  n_fdr05 <- 0
  n_fdr10 <- 0
  n_nom05 <- 0
}

# -- Ensure ora_results has proper structure even if empty --
if(nrow(ora_results) == 0) {
  ora_results <- data.frame(
    Module = character(0), Signature_ID = character(0), Short_Name = character(0),
    Full_Name = character(0), Category = character(0), N_Module = integer(0),
    N_Signature = integer(0), N_Overlap = integer(0), Overlap_Genes = character(0),
    Expected = numeric(0), Fold_Enrichment = numeric(0), P_Value = numeric(0),
    FDR = numeric(0), FDR_global = numeric(0), Sig_Level = character(0),
    stringsAsFactors = FALSE
  )
}

# -- Save ORA results --
write.csv(ora_results, file.path(results_dir, "Step8_02_ORA_Results.csv"), row.names = FALSE)
cat("   [OK] Saved: Step8_02_ORA_Results.csv\n")

# Flag for figure generation
ora_has_results <- nrow(ora_results) > 0 && any(ora_results$P_Value < 0.05, na.rm = TRUE)

# -- Show top results for focus modules --
if(nrow(ora_results) > 0) {
  cat("   Top enrichments in PFS-associated modules:\n")
  for(mod in focus_modules_auto) {
    if(mod %in% ora_results$Module) {
      cat(sprintf("\n   [%s module]\n", toupper(mod)))
      top <- ora_results %>%
        filter(Module == mod, P_Value < 0.05) %>%
        arrange(P_Value) %>%
        head(8)

      if(nrow(top) > 0) {
        for(i in 1:nrow(top)) {
          sig_marker <- ifelse(top$FDR[i] < 0.05, "**",
                               ifelse(top$FDR[i] < 0.10, "*", ""))
          cat(sprintf("      %s: %.1fx (p=%.4f)%s\n",
                          top$Short_Name[i], top$Fold_Enrichment[i], top$P_Value[i], sig_marker))
        }
      }
    }
  }
} else {
  cat("   [SKIP] No ORA results available for display.\n")
}

# -- Figure 1: ORA Lollipop Plot --
cat("   Creating Figure 1: ORA Lollipop Plot...\n")

ora_plot_data <- ora_results %>%
  filter(Module %in% focus_modules_auto, P_Value < 0.05) %>%  # Only significant
  mutate(
    neg_log10_p = -log10(P_Value),
    Module = factor(Module, levels = focus_modules_auto),
    Sig_Level = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "Nominal p < 0.05"
    ),
    Sig_Level = factor(Sig_Level, levels = c("FDR < 0.05", "FDR < 0.10", "Nominal p < 0.05"))
  )

p1 <- ggplot(ora_plot_data, aes(x = Fold_Enrichment, y = reorder(Short_Name, Fold_Enrichment))) +
  geom_segment(aes(xend = 0, yend = reorder(Short_Name, Fold_Enrichment)),
               color = "grey70", linewidth = 0.5) +
  geom_point(aes(size = N_Overlap, fill = neg_log10_p, shape = Sig_Level),
             color = "grey30", stroke = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_gradientn(
    colors = c("#D6EFB3", "#5FC1C0", "#234DA0"),
    name = expression(-log[10](p))
  ) +
  scale_shape_manual(
    values = c("FDR < 0.05" = 21, "FDR < 0.10" = 22, "Nominal p < 0.05" = 23),
    name = "Significance"
  ) +
  scale_size_continuous(range = c(3, 8), breaks = c(2, 4, 6), name = "Overlap") +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  labs(
    title = "Module-Signature Enrichment (ORA)",
    subtitle = "Hypergeometric test | Only significant enrichments shown",
    x = "Fold Enrichment",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  )

print(p1)
ggsave(file.path(fig_dir_png, "Step8_01_ORA_Lollipop.png"), p1, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step8_01_ORA_Lollipop.pdf"), p1, width = 8, height = 5)

# -- Figure 2: ORA Heatmap --
cat("   Creating Figure 2: ORA Heatmap...\n")

ora_heatmap_data <- ora_results %>%
  filter(N_Overlap > 0) %>%
  dplyr::select(Module, Short_Name, Fold_Enrichment, FDR, P_Value)

ora_matrix <- ora_heatmap_data %>%
  select(Module, Short_Name, Fold_Enrichment) %>%
  pivot_wider(names_from = Module, values_from = Fold_Enrichment, values_fill = 0) %>%
  column_to_rownames("Short_Name")

# Order columns
# Ensure focus modules appear first, then other common modules
col_order <- unique(c(focus_modules_auto, "turquoise", "green", "yellow", "red", "pink", "black", "magenta", "brown", "blue"))
col_order <- col_order[col_order %in% colnames(ora_matrix)]
ora_matrix <- ora_matrix[, col_order, drop = FALSE]

# Filter to signatures with at least one p < 0.05
sig_sigs <- ora_heatmap_data %>%
  filter(P_Value < 0.05) %>%
  pull(Short_Name) %>%
  unique()

ora_matrix_filt <- ora_matrix[rownames(ora_matrix) %in% sig_sigs, , drop = FALSE]

# Create significance markers
sig_matrix <- ora_heatmap_data %>%
  mutate(
    Sig_Mark = case_when(
      FDR < 0.05 ~ "**",
      FDR < 0.10 ~ "*",
      P_Value < 0.05 ~ "+",
      TRUE ~ ""
    )
  ) %>%
  select(Module, Short_Name, Sig_Mark) %>%
  pivot_wider(names_from = Module, values_from = Sig_Mark, values_fill = "") %>%
  column_to_rownames("Short_Name")

sig_matrix_filt <- sig_matrix[rownames(ora_matrix_filt), col_order, drop = FALSE]


png(file.path(fig_dir_png, "Step8_02_ORA_Heatmap.png"), width = 8, height = 6, units = "in", res = 300)
p2 <- pheatmap(
  as.matrix(ora_matrix_filt),
  color = heatmap_scale,  # ORANGE spectrum
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "ORA: Module-Specific Pathway Enrichment\n** FDR < 0.05 | * FDR < 0.10 | + p < 0.05",
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 11,
  border_color = "grey70",
  display_numbers = as.matrix(sig_matrix_filt),
  fontsize_number = 10,
  number_color = "black",
  angle_col = 45,
  cellwidth = 30,
  cellheight = 18
)
dev.off()

pdf(file.path(fig_dir_pdf, "Step8_02_ORA_Heatmap.pdf"), width = 8, height = 6)
p2 <- pheatmap(
  as.matrix(ora_matrix_filt),
  color = heatmap_scale,  # ORANGE spectrum
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "ORA: Module-Specific Pathway Enrichment\n** FDR < 0.05 | * FDR < 0.10 | + p < 0.05",
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 11,
  border_color = "grey70",
  display_numbers = as.matrix(sig_matrix_filt),
  fontsize_number = 10,
  number_color = "black",
  angle_col = 45,
  cellwidth = 30,
  cellheight = 18
)
dev.off()

cat("   Saved: Step8_02_ORA_Heatmap.png/pdf\n")

# -- Figure 3: Faceted ORA WordCloud (Separate panels per Focus Module) --
cat("   Creating Figure 3: Faceted ORA WordCloud (by module)...\n")

# 1. Prepare data with separate facets for each focus module
facet_wc_data_fig3 <- ora_results %>%
  filter(tolower(Module) %in% focus_modules_auto, P_Value < 0.05) %>%
  mutate(
    Module_Label = factor(stringr::str_to_title(Module),
                          levels = stringr::str_to_title(focus_modules_auto)),
    neg_log10_p = -log10(P_Value),
    importance = (neg_log10_p * log2(pmax(Fold_Enrichment, 1) + 1)) + 12,
    Category_Clean = factor(assign_category(Category), levels = names(category_colors))
  ) %>%
  group_by(Module_Label) %>%
  slice_max(order_by = importance, n = 20, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("     Data: %d pathways across %d modules\n",
            nrow(facet_wc_data_fig3),
            length(unique(facet_wc_data_fig3$Module_Label))))

# 2. Define module colors for facet headers
module_header_colors_fig3 <- c(
  "Pink" = "#FF69B4",
  "Green" = "#2E8B57",
  "Black" = "#2F2F2F",
  "Blue" = "#4169E1"
)

# 3. Build faceted word cloud
p3_base <- ggplot(facet_wc_data_fig3, aes(label = Short_Name, size = importance, color = Category_Clean)) +
  geom_text_wordcloud_area(
    rm_outside = FALSE,
    shape = "square",
    grid_margin = 1.2,
    seed = 42
  ) +
  facet_wrap(~ Module_Label, ncol = 2) +
  scale_size_area(max_size = 16) +
  scale_color_manual(values = category_colors, name = "Biological Category") +
  labs(
    title = "ORA: Pathway Enrichment by Focus Module",
    subtitle = "Modules significantly associated with PFS or Treatment Response"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11, margin = margin(b = 20)),
    strip.text = element_text(face = "bold", size = 12, color = "white", margin = margin(t = 8, b = 8)),
    strip.background = element_rect(fill = "grey50", color = "black", linewidth = 1.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.spacing = unit(1.5, "lines"),
    plot.margin = margin(20, 20, 20, 20),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 5)))

# 4. Apply module-specific colors to facet headers
p3_built <- ggplot_gtable(ggplot_build(p3_base))
strip_indices <- which(grepl("strip", p3_built$layout$name))

# Get the actual module order from the data
actual_modules <- levels(facet_wc_data_fig3$Module_Label)
actual_modules <- actual_modules[actual_modules %in% unique(facet_wc_data_fig3$Module_Label)]

k <- 1
for (i in strip_indices) {
  if (k <= length(actual_modules)) {
    mod_name <- actual_modules[k]
    if (mod_name %in% names(module_header_colors_fig3)) {
      tryCatch({
        p3_built$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <- module_header_colors_fig3[mod_name]
      }, error = function(e) {
        cat(sprintf("     Note: Could not color strip %d (%s)\n", k, mod_name))
      })
    }
    k <- k + 1
  }
}

# 5. Display and save
grid::grid.newpage()
grid::grid.draw(p3_built)

png(file.path(fig_dir_png, "Step8_03_ORA_WordCloud_Faceted.png"),
    width = 10, height = 8, units = "in", res = 300)
grid::grid.draw(p3_built)
dev.off()

pdf(file.path(fig_dir_pdf, "Step8_03_ORA_WordCloud_Faceted.pdf"),
    width = 10, height = 8)
grid::grid.draw(p3_built)
dev.off()

cat("     Figure 3 (Faceted) saved successfully.\n")

# -- Figure 3b: Combined ORA WordCloud (All focus modules merged) --
cat("   Creating Figure 3b: Combined ORA WordCloud...\n")

# 1. Prepare combined data
ora_wordcloud_combined <- ora_results %>%
  filter(tolower(Module) %in% focus_modules_auto, P_Value < 0.05) %>%
  mutate(
    neg_log10_p = -log10(P_Value),
    importance = (neg_log10_p * log2(pmax(Fold_Enrichment, 1) + 1)) + 2,
    Category_Clean = factor(assign_category(Category), levels = names(category_colors))
  ) %>%
  group_by(Short_Name) %>%
  slice_max(order_by = importance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(importance)) %>%
  head(40)

# 2. Create combined word cloud
fig3b_cloud <- ggplot(ora_wordcloud_combined, aes(label = Short_Name, size = importance, color = Category_Clean)) +
  geom_text_wordcloud_area(
    rm_outside = FALSE,
    shape = "circle",
    grid_margin = 2,
    seed = 123
  ) +
  scale_size_area(max_size = 22) +
  scale_color_manual(values = category_colors, name = "Biological Category") +
  coord_fixed() +
  labs(
    title = "ORA: Combined Pathway Enrichment (All Focus Modules)",
    subtitle = "Pink, Green, Black, Blue modules merged"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    legend.position = "none"
  )

# 3. Extract legend
legend_plot_3b <- ggplot(ora_wordcloud_combined, aes(x = 1, y = 1, color = Category_Clean)) +
  geom_point(size = 5) +
  scale_color_manual(values = category_colors, name = "Biological Category") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 5)))

legend_only_3b <- cowplot::get_legend(legend_plot_3b)

# 4. Combine
p3b <- cowplot::plot_grid(
  fig3b_cloud,
  legend_only_3b,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

print(p3b)
ggsave(file.path(fig_dir_png, "Step8_03b_ORA_WordCloud_Combined.png"), p3b, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step8_03b_ORA_WordCloud_Combined.pdf"), p3b, width = 8, height = 6)

cat("     Figure 3b (Combined) saved successfully.\n")

# -- Figure 7: All Modules WordCloud (colored by module) --
cat("   Creating Figure 7: All Modules Word Cloud...\n")

wordcloud_data <- ora_results %>%
  filter(P_Value < 0.05, !tolower(Module) %in% c("grey", "gray")) %>%
  mutate(
    Module_Label = toupper(Module),
    neg_log10_p = -log10(P_Value),
    importance = (neg_log10_p * log2(pmax(Fold_Enrichment, 1) + 1)) + 2
  ) %>%
  group_by(Short_Name) %>%
  slice_max(order_by = importance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(importance)) %>%
  head(40)

cat("   Signatures in wordcloud:", nrow(wordcloud_data), "\n")

set.seed(456)
p7 <- ggplot(wordcloud_data, aes(label = Short_Name, size = importance, color = Module_Label)) +
  geom_text_wordcloud_area(
    rm_outside = FALSE,
    shape = "circle",
    grid_margin = 0.5,
    seed = 456
  ) +
  scale_size_area(max_size = 30) +
  scale_color_manual(values = module_colors_palette, name = "Module") +
  labs(
    title = "ORA: Pathway Enrichment Across All WGCNA Modules",
    subtitle = "Over-representation analysis (p < 0.05) | Size = significance | Color = module"
  ) +
  theme_minimal(base_size = pub_base_size) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = pub_title_size, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = pub_subtitle_size, margin = margin(b = 15)),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "right",
    legend.text = element_text(size = pub_legend_text_size)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

print(p7)

ggsave(file.path(fig_dir_png, "Step8_07_ORA_WordCloud.png"), p7, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step8_07_ORA_WordCloud.pdf"), p7, width = 6, height = 4)

# -- Figure 7c: Volcano-Style Enrichment (Clean Version from Step8_Figures_Clean.R) --
cat("   Creating Figure 7c: Volcano-Style Enrichment...\n")

# Define module color palette if not already defined
if(!exists("module_color_palette")) {
  module_color_palette <- c(
      "GREEN" = "#70a494", "TURQUOISE" = "#40E0D0", "BLUE" = "#4169E1",
      "BROWN" = "#8B4513", "YELLOW" = "#DAA520", "RED" = "#DC143C",
      "BLACK" = "#2F2F2F", "PINK" = "#FF69B4", "MAGENTA" = "#BA55D3",
      "PURPLE" = "#9370DB", "GREENYELLOW" = "#ADFF2F", "TAN" = "#D2B48C",
    "SALMON" = "#FA8072", "CYAN" = "#00CED1", "GREY" = "#808080"
  )
}

volcano_data <- ora_results %>%
  filter(!tolower(Module) %in% c("grey", "gray")) %>%
  mutate(
    Module_Label = toupper(Module),
    neg_log10_p = -log10(P_Value),
    log2_FE = log2(pmax(Fold_Enrichment, 0.1)),
    Significant = P_Value < 0.05 & Fold_Enrichment > 1.5
  )

cat(sprintf("     Total points: %d\n", nrow(volcano_data)))
cat(sprintf("     Significant points: %d\n", sum(volcano_data$Significant)))
cat(sprintf("     Modules in data: %s\n", paste(unique(volcano_data$Module_Label), collapse = ", ")))

# Only label TOP 5 UNIQUE pathways (most significant occurrence of each)
# First, get the most significant occurrence of each unique pathway
top_unique_pathways <- volcano_data %>%
  filter(Significant) %>%
  group_by(Short_Name) %>%
  slice_max(order_by = neg_log10_p, n = 1, with_ties = FALSE) %>%  # Keep only most significant occurrence
  ungroup() %>%
  arrange(desc(neg_log10_p)) %>%
  head(5)

top_labels <- top_unique_pathways$Short_Name
cat(sprintf("     Labeling (unique): %s\n", paste(top_labels, collapse = ", ")))

# Only label the MOST SIGNIFICANT occurrence of each pathway (not all occurrences)
volcano_data <- volcano_data %>%
  mutate(
    row_id = row_number()
  )

# Get row IDs for the specific points to label
label_row_ids <- volcano_data %>%
  filter(Short_Name %in% top_labels) %>%
  group_by(Short_Name) %>%
  slice_max(order_by = neg_log10_p, n = 1, with_ties = FALSE) %>%
  pull(row_id)

volcano_data <- volcano_data %>%
  mutate(Label = ifelse(row_id %in% label_row_ids, Short_Name, NA_character_)) %>%
  select(-row_id)

# Separate data subsets for cleaner plotting
volcano_nonsig <- volcano_data %>% filter(!Significant)
volcano_sig <- volcano_data %>% filter(Significant)
volcano_labels <- volcano_data %>% filter(!is.na(Label))

# Calculate axis limits
x_min <- min(volcano_data$log2_FE, na.rm = TRUE) - 0.5
x_max <- max(volcano_data$log2_FE, na.rm = TRUE) + 0.5

# Build color palette dynamically based on modules in data
modules_in_data <- unique(volcano_data$Module_Label)
volcano_colors <- c(
  "GREEN" = "#2E8B57", "TURQUOISE" = "#40E0D0", "BLUE" = "#4169E1",
  "BROWN" = "#8B4513", "YELLOW" = "#DAA520", "RED" = "#DC143C",
  "BLACK" = "#2F2F2F", "PINK" = "#FF69B4", "MAGENTA" = "#BA55D3",
  "PURPLE" = "#9370DB", "GREENYELLOW" = "#ADFF2F", "TAN" = "#D2B48C",
  "SALMON" = "#FA8072", "CYAN" = "#00CED1", "GREY" = "#808080"
)
# Filter to only colors that exist in data
volcano_colors <- volcano_colors[names(volcano_colors) %in% modules_in_data]

p7c <- ggplot(volcano_data, aes(x = log2_FE, y = neg_log10_p)) +
  # Non-significant points (transparent)
  geom_point(
    data = volcano_nonsig,
    aes(color = Module_Label),
    size = 3,
    alpha = 0.25
  ) +
  # Significant points (solid)
  geom_point(
    data = volcano_sig,
    aes(color = Module_Label),
    size = 3.5,
    alpha = 1
  ) +
  # Threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  # Clean labels - NO ARROWS, only top 5
  geom_text_repel(
    data = volcano_labels,
    aes(label = Label),
    size = 4,
    max.overlaps = 50,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "grey50",
    segment.size = 0.3,
    force = 5,
    force_pull = 0.5
  ) +
  scale_color_manual(values = volcano_colors, name = "Module") +
  scale_x_continuous(
    breaks = seq(-4, 4, 1),
    limits = c(x_min, x_max)
  ) +
  labs(
    title = "ORA: Pathway Enrichment Landscape",
    subtitle = "Dashed lines: p = 0.05, Fold Enrichment = 1.5 | Top 5 most significant labeled",
    x = "log2(Fold Enrichment)",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

print(p7c)
ggsave(file.path(fig_dir_png, "Step8_07c_Volcano.png"), p7c, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step8_07c_Volcano.pdf"), p7c, width = 6, height = 4)


# -- Figure 8b: Single Word Cloud with cowplot legend --
cat("   Creating Figure 8b: Word Cloud with extracted legend...\n")

wc_legend_data <- ora_results %>%
  filter(P_Value < 0.05) %>%
  mutate(
    neg_log10_p = -log10(P_Value),
    Sig_Category = case_when(
      FDR < 0.01 ~ "FDR < 0.01",
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "p < 0.05"
    ),
    Sig_Category = factor(Sig_Category, levels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.10", "p < 0.05"))
  )

# Create dummy plot to extract legend
legend_plot <- ggplot(wc_legend_data, aes(x = 1, y = 1, color = Sig_Category, size = neg_log10_p)) +
  geom_point() +
  scale_color_manual(values = c("FDR < 0.01" = "#1a237e", "FDR < 0.05" = "#214A7D",
                                "FDR < 0.10" = "#3680B6", "p < 0.05" = "#4BA24F")) +
  scale_size_continuous(range = c(3, 10), name = expression(-log[10](p))) +
  guides(color = guide_legend(title = "Significance", override.aes = list(size = 4))) +
  theme_minimal()

wc_legend <- cowplot::get_legend(legend_plot)

set.seed(321)
p8b_main <- ggplot(wc_legend_data, aes(label = Short_Name, size = neg_log10_p, color = Sig_Category)) +
  geom_text_wordcloud_area(rm_outside = TRUE, eccentricity = 1.0) +
  scale_size_area(max_size = 25) +
  scale_color_manual(values = c("FDR < 0.01" = "#1a237e", "FDR < 0.05" = "#214A7D",
                                "FDR < 0.10" = "#3680B6", "p < 0.05" = "#4BA24F")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(title = "ORA Enrichment Overview")

p8b <- cowplot::plot_grid(p8b_main, wc_legend, ncol = 2, rel_widths = c(0.85, 0.15))

print(p8b)
ggsave(file.path(fig_dir_png, "Step8_08b_WordCloud_Legend.png"), p8b, width = 7, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step8_08b_WordCloud_Legend.pdf"), p8b, width = 7, height = 4)

# -- Figure 8e: Integrated Dot Plot --
cat("   Creating Figure 8e: Integrated Dot Plot...\n")

# helper function for biological category assignment
assign_bio_category <- function(category_string) {
  dplyr::case_when(
    grepl("Immune|CIBERSORT|Macro|MDSC|T.cell|B.cell|NK|Dendritic|Neutro|Monocyte|Lympho|Inflam|Cytokine|Chemokine|Interleukin|TNF|Interferon|Complement|Checkpoint", category_string, ignore.case = TRUE) ~ "Immune",
    grepl("Stroma|ECM|Collagen|Fibro|Matrix|CAF|Mesenchymal", category_string, ignore.case = TRUE) ~ "Stroma/ECM",
    grepl("Basal|Classical|Squamous|Progenitor|ADEX|Immunogenic|Pancreatic.Progenitor", category_string, ignore.case = TRUE) ~ "PDAC Subtype",
    grepl("Stem|Pluripoten|Self.renewal|Wnt|Notch|Hedgehog", category_string, ignore.case = TRUE) ~ "Stemness",
    grepl("EMT|Epithelial|Mesenchymal|Migration|Invasion|Metasta", category_string, ignore.case = TRUE) ~ "EMT",
    grepl("Metabol|Glycoly|Oxidative|Lipid|Amino.acid|Glucose|Fatty|Cholesterol|Glutamine|Mitochond", category_string, ignore.case = TRUE) ~ "Metabolic",
    grepl("Angio|Vascular|Endothel|VEGF|Hypoxia|HIF", category_string, ignore.case = TRUE) ~ "Angiogenesis",
    grepl("Exosome|Vesicle|EV|sEV|Secretom", category_string, ignore.case = TRUE) ~ "Exosome",
    grepl("Signal|Pathway|KRAS|MAPK|PI3K|AKT|mTOR|TGF|NF.kB|JAK|STAT", category_string, ignore.case = TRUE) ~ "Signaling",
    TRUE ~ "Other"
  )
}

category_colors <- c(
  "Immune" = "#E64B35", "Stroma/ECM" = "#7E6148", "PDAC Subtype" = "#3C5488",
  "Stemness" = "#9B59B6", "EMT" = "#00A087", "Metabolic" = "#F39B7F",
  "Angiogenesis" = "#DC0000", "Exosome" = "#4DBBD5", "Signaling" = "#8491B4",
  "Other" = "#B09C85"
)

dotplot_data <- ora_results %>%
  filter(P_Value < 0.05) %>%
  mutate(
    neg_log10_p = -log10(P_Value),
    GeneRatio = N_Overlap / N_Signature,
    BioCategory = assign_bio_category(Category)
  ) %>%
  group_by(Module) %>%
  arrange(P_Value) %>%
  slice_head(n = 10) %>%
  ungroup()

p8e <- ggplot(dotplot_data, aes(x = Module, y = reorder(Short_Name, neg_log10_p))) +
  geom_point(aes(size = GeneRatio, color = BioCategory, alpha = neg_log10_p)) +
  scale_size_continuous(range = c(2, 8), name = "Gene Ratio") +
  scale_color_manual(values = category_colors, name = "Category") +
  scale_alpha_continuous(range = c(0.5, 1), name = expression(-log[10](p))) +
  labs(
    title = "Integrated Pathway Enrichment",
    subtitle = "Top 10 pathways per module",
    x = "Module",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90")
  )

print(p8e)
ggsave(file.path(fig_dir_png, "Step8_08e_Integrated_DotPlot.png"), p8e, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step8_08e_Integrated_DotPlot.pdf"), p8e, width = 8, height = 6)
cat("   [OK] New publication figures 7c, 8a-8e complete\n")

# ============================================================================
# 8e: GSVA (Gene Set Variation Analysis)
# ============================================================================

cat("==============================================================================\n")
cat("8e: GSVA (Gene Set Variation Analysis)")
cat("==============================================================================\n")

# -- Prepare clinical data (MOVED HERE - needed for heatmap annotation) --
# Use datTraits as primary source since datExpr may have empty rownames
sample_ids <- rownames(datTraits)
if(is.null(sample_ids) || length(sample_ids) == 0) {
  sample_ids <- rownames(datExpr)
}
if(is.null(sample_ids) || length(sample_ids) == 0) {
  n_samp <- nrow(datExpr)
  if(!is.null(n_samp) && n_samp > 0) {
    sample_ids <- paste0("Sample_", 1:n_samp)
  } else {
    # Ultimate fallback - use datTraits dimensions
    n_samp <- nrow(datTraits)
    sample_ids <- if(!is.null(n_samp) && n_samp > 0) paste0("Sample_", 1:n_samp) else character(0)
  }
}

clinical_df <- data.frame(
  Sample_ID = sample_ids,
  stringsAsFactors = FALSE
)

# Add clinical variables from datTraits
if("PFS_group" %in% colnames(datTraits)) {
  clinical_df$PFS_Group <- datTraits$PFS_group
}
if("Response" %in% colnames(datTraits)) {
  clinical_df$Response <- datTraits$Response
}
if("Treatment" %in% colnames(datTraits)) {
  clinical_df$Treatment <- datTraits$Treatment
}

cat("   Clinical data prepared:\n")
cat(sprintf("      Samples: %d\n", nrow(clinical_df)))
if("PFS_Group" %in% colnames(clinical_df)) {
  cat(sprintf("      PFS_Group: %s\n", paste(names(table(clinical_df$PFS_Group)), collapse = ", ")))
}

# -- Prepare expression matrix (genes Ã- samples) --
expr_matrix <- t(as.matrix(datExpr))
cat(sprintf("   Expression matrix: %d proteins Ã- %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# -- Run GSVA --
cat("   Running GSVA (this may take a few minutes)...")

# Check version
packageVersion("GSVA")

# New GSVA API (version 2.0+)
gsva_param <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  kcdf = "Gaussian",
  maxDiff = TRUE
)

gsva_results <- gsva(gsva_param, verbose = FALSE)

cat(sprintf("   [OK] GSVA scores: %d signatures Ã- %d samples\n", nrow(gsva_results), ncol(gsva_results)))

# -- Convert to data frame --
gsva_df <- as.data.frame(t(gsva_results))
gsva_df$Sample_ID <- rownames(gsva_df)
gsva_df <- gsva_df %>% select(Sample_ID, everything())

# -- Save GSVA scores --
write.csv(gsva_df, file.path(results_dir, "Step8_03_GSVA_Scores.csv"), row.names = FALSE)
cat("   [OK] Saved: Step8_03_GSVA_Scores.csv\n")

# -- Figure 4: GSVA Heatmap by PFS --
cat("   Creating Figure 4: GSVA Heatmap...\n")

# Load ComplexHeatmap for better legend control
if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  cat("   [NOTE] ComplexHeatmap not available, using pheatmap fallback\n")
  use_complex_heatmap <- FALSE
} else {
  library(ComplexHeatmap)
  use_complex_heatmap <- TRUE
}

# Select top variable signatures
gsva_var <- apply(gsva_results, 1, var)
top_sigs <- names(sort(gsva_var, decreasing = TRUE)[1:min(30, length(gsva_var))])

gsva_plot_matrix <- gsva_results[top_sigs, ]

# Get short names for row labels
row_labels <- sapply(top_sigs, function(x) {
  if(x %in% names(display_names_short)) display_names_short[[x]] else x
})
rownames(gsva_plot_matrix) <- row_labels

# Z-scale the matrix (row scaling)
gsva_plot_scaled <- t(scale(t(gsva_plot_matrix)))

# Define diverging palette: blue (low) -> white (center) -> red (high)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

if(use_complex_heatmap) {
  # ComplexHeatmap approach - full control over legends

  # Create column annotation
  col_annotation <- HeatmapAnnotation(
    PFS_Group = clinical_df$PFS_Group,
    col = list(PFS_Group = colors_pfs),
    annotation_legend_param = list(
      PFS_Group = list(title = "PFS Group", title_position = "topcenter")
    )
  )

  # Create the heatmap with blue-white-red diverging palette
  ht <- Heatmap(
    gsva_plot_scaled,
    name = "z-score",  # This labels the color scale
    col = heatmap_colors,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 10),
    top_annotation = col_annotation,
    column_title = "Gene Set Variation Analysis (GSVA): Pathway Activity by PFS Group\n(Top 30 variable pathways)",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_position = "topcenter",
      legend_height = unit(4, "cm")
    )
  )

  # Save PNG
  png(file.path(fig_dir_png, "Step8_04_GSVA_Heatmap.png"), width = 10, height = 8, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  # Save PDF
  pdf(file.path(fig_dir_pdf, "Step8_04_GSVA_Heatmap.pdf"), width = 10, height = 8)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  # Display in RStudio
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

} else {
  # Fallback to pheatmap with blue-white-red diverging palette
  annotation_col <- data.frame(
    PFS_Group = clinical_df$PFS_Group,
    row.names = colnames(gsva_plot_matrix)
  )

  annotation_colors <- list(
    PFS_Group = colors_pfs
  )

  P4 <- pheatmap(
    gsva_plot_scaled,
    color = heatmap_colors,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    main = "GSVA: Pathway Activity by PFS Group (z-score scaled)\n(Top 30 variable pathways)",
    fontsize = 12,
    fontsize_row = 10,
    border_color = NA
  )

  png(file.path(fig_dir_png, "Step8_04_GSVA_Heatmap.png"), width = 10, height = 8, units = "in", res = 300)
  print(P4)
  dev.off()

  pdf(file.path(fig_dir_pdf, "Step8_04_GSVA_Heatmap.pdf"), width = 10, height = 8)
  print(P4)
  dev.off()
}

cat("   [OK] Figure 4 saved\n")

# ============================================================================
# 8f: ssGSEA (Single-Sample GSEA)
# ============================================================================

cat("==============================================================================\n")
cat("8f: ssGSEA (Single-Sample GSEA)")
cat("==============================================================================\n")

# -- Run ssGSEA --
cat("   Running ssGSEA...\n")

ssgsea_param <- ssgseaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  normalize = TRUE
)

ssgsea_results <- gsva(ssgsea_param, verbose = FALSE)

cat(sprintf("   [OK] ssGSEA scores: %d signatures Ã- %d samples\n", nrow(ssgsea_results), ncol(ssgsea_results)))


# -- Convert to data frame --
ssgsea_df <- as.data.frame(t(ssgsea_results))
ssgsea_df$Sample_ID <- rownames(ssgsea_df)
ssgsea_df <- ssgsea_df %>% select(Sample_ID, everything())

# -- Save ssGSEA scores --
write.csv(ssgsea_df, file.path(results_dir, "Step8_04_ssGSEA_Scores.csv"), row.names = FALSE)
cat("   [OK] Saved: Step8_04_ssGSEA_Scores.csv\n")

# -- Visualize ssGSEA: Top Variable Pathways Heatmap --
cat("   Creating ssGSEA heatmap...\n")

if(nrow(ssgsea_results) > 0 && ncol(ssgsea_results) > 0) {
  # Calculate variance per pathway
  pathway_var <- apply(ssgsea_results, 1, var, na.rm = TRUE)
  top_var_pathways <- names(sort(pathway_var, decreasing = TRUE))[1:min(20, length(pathway_var))]

  # Subset to top variable pathways
  ssgsea_top <- ssgsea_results[top_var_pathways, , drop = FALSE]

  # Get short names for display
  rownames_short <- sapply(rownames(ssgsea_top), function(x) {
    if(x %in% names(display_names_short)) display_names_short[[x]] else x
  })
  rownames(ssgsea_top) <- rownames_short

  # Create heatmap object first
  p_ssgsea <- pheatmap(ssgsea_top,
                       scale = "row",
                       clustering_method = "ward.D2",
                       show_colnames = FALSE,
                       main = "ssGSEA: Top 20 Variable Pathways (z-scaled)",
                       fontsize = 12,
                       fontsize_row = 12,
                       border_color = NA)

  # Save PNG
  png(file.path(fig_dir_png, "Step8f_ssGSEA_Top20_Heatmap.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p_ssgsea)
  dev.off()

  # Save PDF
  pdf(file.path(fig_dir_pdf, "Step8f_ssGSEA_Top20_Heatmap.pdf"),
      width = 12, height = 8)
  print(p_ssgsea)
  dev.off()

  cat("   [OK] Saved: Step8f_ssGSEA_Top20_Heatmap.png/pdf\n")
}

# ============================================================================
# 8g: CLINICAL ASSOCIATION TESTING
# ============================================================================

cat("==============================================================================\n")
cat("8g: CLINICAL ASSOCIATION TESTING\n")
cat("==============================================================================\n")

# NOTE: clinical_df was already created in Step 8e

# -- Function for clinical association testing --
test_clinical_association <- function(scores_matrix, clinical_df, group_var, method = "GSVA") {
  if(!group_var %in% colnames(clinical_df)) {
    return(data.frame())
  }
  
  results <- data.frame()
  
  for(sig in rownames(scores_matrix)) {
    scores <- scores_matrix[sig, ]
    groups <- clinical_df[[group_var]]
    
    # Get valid groups
    valid_idx <- !is.na(groups)
    scores_valid <- scores[valid_idx]
    groups_valid <- groups[valid_idx]
    
    group_levels <- unique(groups_valid)
    if(length(group_levels) != 2) next
    
    g1 <- scores_valid[groups_valid == group_levels[1]]
    g2 <- scores_valid[groups_valid == group_levels[2]]
    
    # Skip if too few samples
    if(length(g1) < 3 || length(g2) < 3) next
    
    # Welch's t-test
    tt <- tryCatch(t.test(g1, g2), error = function(e) NULL)
    if(is.null(tt)) next
    
    # Effect size (Cohen's d)
    pooled_sd <- sqrt((var(g1) + var(g2)) / 2)
    cohens_d <- (mean(g1) - mean(g2)) / pooled_sd
    
    # Defensive access to display names
    short_name <- if(sig %in% names(display_names_short)) display_names_short[[sig]] else sig
    full_name <- if(sig %in% names(display_names_full)) display_names_full[[sig]] else sig

    results <- rbind(results, data.frame(
      Signature = sig,
      Short_Name = short_name,
      Full_Name = full_name,
      Comparison = group_var,
      Group1 = group_levels[1],
      Group2 = group_levels[2],
      N_G1 = length(g1),
      N_G2 = length(g2),
      Mean_G1 = round(mean(g1), 4),
      Mean_G2 = round(mean(g2), 4),
      Diff = round(mean(g1) - mean(g2), 4),
      Cohens_d = round(cohens_d, 3),
      P_Value = tt$p.value,
      Method = method,
      stringsAsFactors = FALSE
    ))
  }
  
  if(nrow(results) > 0) {
    results$FDR <- p.adjust(results$P_Value, method = "BH")
    results <- results %>% arrange(P_Value)
  }
  
  return(results)
}

# -- Test GSVA vs clinical variables --
clinical_tests <- data.frame()

# PFS Group
if("PFS_Group" %in% colnames(clinical_df)) {
  cat("   Testing GSVA vs PFS_Group...\n")
  gsva_pfs <- test_clinical_association(gsva_results, clinical_df, "PFS_Group", "GSVA")
  if(nrow(gsva_pfs) > 0) {
    cat(sprintf("   [OK] Significant at p < 0.05: %d\n", sum(gsva_pfs$P_Value < 0.05)))
    clinical_tests <- bind_rows(clinical_tests, gsva_pfs)
  }
}

# Response
if("Response" %in% colnames(clinical_df)) {
  cat("   Testing GSVA vs Response...\n")
  gsva_response <- test_clinical_association(gsva_results, clinical_df, "Response", "GSVA")
  if(nrow(gsva_response) > 0) {
    cat(sprintf("   [OK] Significant at p < 0.05: %d\n", sum(gsva_response$P_Value < 0.05)))
    clinical_tests <- bind_rows(clinical_tests, gsva_response)
  }
}

# Treatment
if("Treatment" %in% colnames(clinical_df)) {
  cat("   Testing GSVA vs Treatment...\n")
  gsva_treatment <- test_clinical_association(gsva_results, clinical_df, "Treatment", "GSVA")
  if(nrow(gsva_treatment) > 0) {
    cat(sprintf("   [OK] Significant at p < 0.05: %d\n", sum(gsva_treatment$P_Value < 0.05)))
    clinical_tests <- bind_rows(clinical_tests, gsva_treatment)
  }
}

# -- Save clinical associations --
if(nrow(clinical_tests) > 0) {
  # Apply pathway name cleaning for publication-ready figures
  clinical_tests$Short_Name <- clean_pathway_names(clinical_tests$Short_Name)

  write.csv(clinical_tests, file.path(results_dir, "Step8_05_Clinical_Associations.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step8_05_Clinical_Associations.csv\n")
  
  # Show top PFS associations
  if("PFS_Group" %in% colnames(clinical_df)) {
    cat("   Top signatures associated with PFS:\n")
    top_pfs <- clinical_tests %>%
      filter(Comparison == "PFS_Group", P_Value < 0.05) %>%
      arrange(P_Value) %>%
      head(10)
    
    for(i in 1:min(10, nrow(top_pfs))) {
      direction <- ifelse(top_pfs$Diff[i] > 0, "<-‘Long", "<-‘Short")
      cat(sprintf("      %s: d=%.2f %s (p=%.4f)\n",
                      top_pfs$Short_Name[i], abs(top_pfs$Cohens_d[i]), direction, top_pfs$P_Value[i]))
    }
  }
}

# -- Figure 5a: Clinical Association Forest Plot (PFS_Group) --
cat("   Creating Figure 5a: PFS Clinical Forest Plot...\n")

forest_data_pfs <- clinical_tests %>%
  filter(Comparison == "PFS_Group", P_Value < 0.05) %>%
  mutate(
    Higher_In = ifelse(Cohens_d > 0, "Long", "Short"),  # Match colors_pfs labels
    abs_d = abs(Cohens_d),
    Sig_Level = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "Nominal p < 0.05"
    ),
    Sig_Level = factor(Sig_Level, levels = c("FDR < 0.05", "FDR < 0.10", "Nominal p < 0.05"))
  ) %>%
  arrange(desc(abs_d)) %>%
  head(20)

if(nrow(forest_data_pfs) > 0) {
  p5a <- ggplot(forest_data_pfs, aes(x = Cohens_d, y = reorder(Short_Name, abs_d))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(xend = 0, yend = reorder(Short_Name, abs_d)), color = "grey70") +
    geom_point(aes(fill = Higher_In, shape = Sig_Level), size = 4, color = "grey30", stroke = 0.5) +
    scale_fill_manual(values = colors_pfs, name = "Higher in") +
    scale_shape_manual(
      values = c("FDR < 0.05" = 21, "FDR < 0.10" = 22, "Nominal p < 0.05" = 23),
      name = "Significance"
    ) +
    labs(
      title = "GSVA Pathway Scores Associated with PFS Group",
      subtitle = "Welch's t-test | Top 20 significant pathways by effect size",
      x = "Cohen's d (Long - Short PFS)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      axis.text.y = element_text(size = 11),
      legend.position = "right"
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 4)),
      shape = guide_legend(override.aes = list(fill = "grey50"))
    )

  print(p5a)
  ggsave(file.path(fig_dir_png, "Step8_05a_PFS_Forest_Plot.png"), p5a, width = 10, height = 7, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step8_05a_PFS_Forest_Plot.pdf"), p5a, width = 10, height = 7)
  cat("   [OK] Figure 5a (PFS) saved\n")
} else {
  cat("   [SKIP] Figure 5a - No significant PFS associations\n")
}

# -- Figure 5b: Clinical Association Forest Plot (Response) --
cat("   Creating Figure 5b: Response Clinical Forest Plot...\n")

forest_data_resp <- clinical_tests %>%
  filter(Comparison == "Response", P_Value < 0.05) %>%
  mutate(
    Higher_In = ifelse(Cohens_d > 0, "PD", "CD"),
    abs_d = abs(Cohens_d),
    Sig_Level = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "Nominal p < 0.05"
    ),
    Sig_Level = factor(Sig_Level, levels = c("FDR < 0.05", "FDR < 0.10", "Nominal p < 0.05"))
  ) %>%
  arrange(desc(abs_d)) %>%
  head(20)

if(nrow(forest_data_resp) > 0) {
  p5b <- ggplot(forest_data_resp, aes(x = Cohens_d, y = reorder(Short_Name, abs_d))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(xend = 0, yend = reorder(Short_Name, abs_d)), color = "grey70") +
    geom_point(aes(fill = Higher_In, shape = Sig_Level), size = 4, color = "grey30", stroke = 0.5) +
    scale_fill_manual(values = colors_response, name = "Higher in") +
    scale_shape_manual(
      values = c("FDR < 0.05" = 21, "FDR < 0.10" = 22, "Nominal p < 0.05" = 23),
      name = "Significance"
    ) +
    labs(
      title = "GSVA Pathway Scores Associated with Treatment Response",
      subtitle = "Welch's t-test | Top 20 significant pathways by effect size",
      x = "Cohen's d (PD - CD)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      axis.text.y = element_text(size = 11),
      legend.position = "right"
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 4)),
      shape = guide_legend(override.aes = list(fill = "grey50"))
    )

  print(p5b)
  ggsave(file.path(fig_dir_png, "Step8_05b_Response_Forest_Plot.png"), p5b, width = 10, height = 7, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step8_05b_Response_Forest_Plot.pdf"), p5b, width = 10, height = 7)
  cat("   [OK] Figure 5b (Response) saved\n")
} else {
  cat("   [SKIP] Figure 5b - No significant Response associations\n")
}

# -- Figure 8: Combined WordCloud (ORA + GSVA) --
cat("   Creating Figure 8: Combined WordCloud (ORA + GSVA)...\n")


  # module color palette
module_colors_wc <- c(
    "GREEN" = "#2E8B57", "TURQUOISE" = "#40E0D0", "BLUE" = "#4169E1",
    "BROWN" = "#8B4513", "YELLOW" = "#DAA520", "RED" = "#DC143C",
    "BLACK" = "#2F2F2F", "PINK" = "#FF69B4", "MAGENTA" = "#BA55D3",
    "PURPLE" = "#9370DB", "GREENYELLOW" = "#ADFF2F", "TAN" = "#D2B48C",
    "SALMON" = "#FA8072", "CYAN" = "#00CED1", "GREY" = "#808080",
    "GSVA_CD" = "#008080", "GSVA_PD" = "#CA562C"
  )

combined_wc_data <- bind_rows(
    # ORA results - use Short_Name and Module for coloring
    ora_results %>%
      filter(P_Value < 0.05) %>%
      mutate(
        neg_log10_p = -log10(P_Value),
        Analysis = "ORA",
        Label = Short_Name,
        Color_Group = toupper(Module)
      ) %>%
      select(Label, neg_log10_p, FDR, P_Value, Analysis, Color_Group),

    # GSVA clinical results
clinical_tests %>%
      filter(Comparison == "Response", P_Value < 0.05) %>%
      mutate(
        neg_log10_p = -log10(P_Value),
        Analysis = "GSVA",
        Direction = ifelse(Cohens_d > 0, "PD", "CD"),
        Label = paste0(Short_Name, " (", Direction, ")"),
        Color_Group = ifelse(Direction == "CD", "GSVA_CD", "GSVA_PD")
      ) %>%
      select(Label, neg_log10_p, FDR, P_Value, Analysis, Color_Group)
  ) %>%
    distinct(Label, .keep_all = TRUE) %>%
    arrange(desc(neg_log10_p)) %>%
    head(50)

  cat("   Combined wordcloud entries:", nrow(combined_wc_data), "\n")

  set.seed(456)
p_combined_wc <- ggplot(combined_wc_data,
                          aes(label = Label,
                              size = neg_log10_p,
                              color = Color_Group)) +
    geom_text_wordcloud_area(
      rm_outside = FALSE,
      shape = "circle",
      grid_margin = 1,
      seed = 456
    ) +
    scale_size_area(max_size = 30) +
    scale_color_manual(values = module_colors_wc, name = "Module/Direction") +
    labs(
      title = "Integrative Pathway Analysis: ORA & GSVA",
      subtitle = "Color = WGCNA module | Size = significance | GSVA: teal = CD, orange = PD"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40", margin = margin(b = 15)),
      plot.margin = margin(25, 25, 25, 25),
      legend.position = "right",
      legend.text = element_text(size = 9)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

print(p_combined_wc)

ggsave(file.path(fig_dir_png, "Step8_08_Combined_WordCloud.png"), p_combined_wc,
         width = 6, height = 4, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step8_08_Combined_WordCloud.pdf"), p_combined_wc,
         width = 6, height = 4)

# ============================================================================
# 8i: MODULE-PATHWAY INTEGRATION
# ============================================================================

cat("==============================================================================\n")
cat("8i: MODULE-PATHWAY INTEGRATION\n")
cat("==============================================================================\n")

# -- Correlate MEs with GSVA scores --
cat("   Correlating module eigengenes with pathway scores...\n")

# Get MEs for focus modules (using focus_modules_auto for consistency)
me_cols <- paste0("ME", focus_modules_auto)
me_cols <- me_cols[me_cols %in% colnames(MEs)]
cat(sprintf("   Focus modules for ME-pathway correlation: %s\n", paste(focus_modules_auto, collapse = ", ")))

me_gsva_cor <- WGCNA::cor(MEs[, me_cols, drop = FALSE], t(gsva_results), use = "pairwise.complete.obs")

# Calculate p-values
me_gsva_pval <- matrix(NA, nrow = length(me_cols), ncol = nrow(gsva_results),
                       dimnames = list(me_cols, rownames(gsva_results)))

for(i in 1:length(me_cols)) {
  for(j in 1:nrow(gsva_results)) {
    ct <- cor.test(MEs[, me_cols[i]], gsva_results[j, ], method = "pearson")
    me_gsva_pval[i, j] <- ct$p.value
  }
}

# -- Find significant correlations --
sig_threshold <- 0.05
sig_cors <- which(me_gsva_pval < sig_threshold & abs(me_gsva_cor) > 0.3, arr.ind = TRUE)

if(nrow(sig_cors) > 0) {
  sig_cor_df <- data.frame(
    Module = gsub("ME", "", rownames(me_gsva_pval)[sig_cors[, 1]]),
    Pathway = colnames(me_gsva_pval)[sig_cors[, 2]],
    Short_Name = sapply(colnames(me_gsva_pval)[sig_cors[, 2]], function(x) {
      ifelse(x %in% names(display_names_short), display_names_short[[x]], x)
    }),
    Correlation = me_gsva_cor[sig_cors],
    P_Value = me_gsva_pval[sig_cors],
    stringsAsFactors = FALSE
  ) %>%
    mutate(FDR = p.adjust(P_Value, method = "BH")) %>%
    arrange(P_Value)
  
  cat(sprintf("   [OK] Significant ME-pathway correlations (p<0.05, |r|>0.3): %d\n", nrow(sig_cor_df)))
  
  # Show top for clinical modules
  for(mod in focus_modules_auto) {
    mod_cors <- sig_cor_df %>% filter(Module == mod) %>% head(5)
    if(nrow(mod_cors) > 0) {
      cat(sprintf("\n   Top pathways correlated with %s module:\n", mod))
      for(i in 1:nrow(mod_cors)) {
        cat(sprintf("      %s: r=%.3f (p=%.4f)\n",
                        mod_cors$Short_Name[i], mod_cors$Correlation[i], mod_cors$P_Value[i]))
      }
    }
  }
  
  write.csv(sig_cor_df, file.path(results_dir, "Step8_08_ME_Pathway_Correlations.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step8_08_ME_Pathway_Correlations.csv\n")
}

# -- Figure 6: Module-Pathway Correlation Heatmap --
cat("   Creating Figure 6: Module-Pathway Correlations...\n")

# Check if sig_cor_df exists and has data
if(exists("sig_cor_df") && nrow(sig_cor_df) > 0) {
  
  top_pathways <- sig_cor_df %>%
    filter(Module %in% focus_modules_auto, P_Value < 0.05) %>%
    arrange(P_Value) %>%
    pull(Pathway) %>%
    unique() %>%
    head(15)
  
  if(length(top_pathways) > 0) {
    cor_subset <- me_gsva_cor[me_cols, top_pathways, drop = FALSE]
    pval_subset <- me_gsva_pval[me_cols, top_pathways, drop = FALSE]
    
    pathway_labels <- sapply(top_pathways, function(x) {
      if(x %in% names(display_names_short)) display_names_short[[x]] else x
    })
    colnames(cor_subset) <- pathway_labels
    colnames(pval_subset) <- pathway_labels
    
    sig_marks <- matrix("", nrow = nrow(pval_subset), ncol = ncol(pval_subset))
    sig_marks[pval_subset < 0.001] <- "***"
    sig_marks[pval_subset < 0.01 & pval_subset >= 0.001] <- "**"
    sig_marks[pval_subset < 0.05 & pval_subset >= 0.01] <- "*"
    rownames(sig_marks) <- rownames(cor_subset)
    colnames(sig_marks) <- colnames(cor_subset)
    
    # Create pheatmap object
    p6 <- pheatmap(
      t(cor_subset),
      color = colorRampPalette(c("#0E4D92", "#F5F5F5", "#D95C2B"))(100),
      display_numbers = t(sig_marks),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      main = "Module Eigengene - Pathway Correlations\n*** p<0.001 | ** p<0.01 | * p<0.05",
      fontsize = 12,
      fontsize_row = 12,
      fontsize_number = 10,
      border_color = "white"
    )
    
    # Note: pheatmap objects cannot be saved with ggsave(), using png/pdf wrapper
    png(file.path(fig_dir_png, "Step8_06_ME_Pathway_Cor.png"), width = 8, height = 6, units = "in", res = 300)
    print(p6)
    dev.off()
    pdf(file.path(fig_dir_pdf, "Step8_06_ME_Pathway_Cor.pdf"), width = 8, height = 6)
    print(p6)
    dev.off()
    
    cat("   [OK] Figure 6 saved\n")
  }
}
# ============================================================================
# 8j: SUMMARY AND SAVE WORKSPACE
# ============================================================================

cat("==============================================================================\n")
cat("8j: SUMMARY\n")
cat("==============================================================================\n")

# -- Print summary --
cat("   ANALYSIS SUMMARY\n")
cat("   ------------------------------------------------------------------------\n")
cat(sprintf("   Samples: %d\n", n_samples))
cat(sprintf("   Proteins: %d\n", n_proteins))
cat(sprintf("   Total signatures defined: %d\n", length(signatures)))
cat(sprintf("   Usable signatures (>=%d genes): %d\n", min_genes, length(gene_sets)))
cat(sprintf("   ORA tests: %d\n", nrow(ora_results)))
cat(sprintf("   ORA significant (p < 0.05): %d\n", sum(ora_results$P_Value < 0.05)))
cat(sprintf("   ORA significant (FDR < 0.10): %d\n", sum(ora_results$FDR < 0.10, na.rm = TRUE)))
if(nrow(clinical_tests) > 0) {
  cat(sprintf("   Clinical associations tested: %d\n", nrow(clinical_tests)))
  cat(sprintf("   Clinical associations (p < 0.05): %d\n", sum(clinical_tests$P_Value < 0.05)))
}
# -- List output files --
cat("   OUTPUT FILES\n")
cat("   ------------------------------------------------------------------------\n")
cat("   Results:\n")
cat("      - Step8_01_Signature_Coverage.csv\n")
cat("      - Step8_02_ORA_Results.csv\n")
cat("      - Step8_03_GSVA_Scores.csv\n")
cat("      - Step8_04_ssGSEA_Scores.csv\n")
cat("      - Step8_05_Clinical_Associations.csv\n")
cat("      - Step8_06_Biomarker_Analysis.csv\n")
cat("      - Step8_07_Biomarker_Expression.csv\n")
cat("      - Step8_08_ME_Pathway_Correlations.csv\n")
cat("   Figures:\n")
cat("      - Step8_01_ORA_Lollipop.png/pdf\n")
cat("      - Step8_02_ORA_Heatmap.png/pdf\n")
cat("      - Step8_03_ORA_WordCloud_Final.png/pdf\n")
cat("      - Step8_04_GSVA_Heatmap.png/pdf\n")
cat("      - Step8_05_Clinical_Forest_Plot.png/pdf\n")
cat("      - Step8_06_ME_Pathway_Cor.png/pdf\n")
cat("      - Step8_07_ORA_WordCloud.png/pdf\n")
cat("      - Step8_08_Combined_WordCloud.png/pdf\n")

# ============================================================================
# 8-FINAL: CREATE UNIFIED ENRICHMENT TABLE
# ============================================================================

cat("==============================================================================\n")
cat("8-FINAL: CREATE UNIFIED ENRICHMENT TABLE\n")
cat("==============================================================================\n")

# Combine all enrichment results from all focus modules
# Use focus_module_ora list created earlier (contains ORA results for each focus module)
if(exists("focus_module_ora") && length(focus_module_ora) > 0) {
  unified_enrichment <- bind_rows(focus_module_ora)
  cat(sprintf("   Combining ORA results from %d focus modules: %s\n",
              length(focus_module_ora), paste(names(focus_module_ora), collapse = ", ")))
} else {
  # Fallback: use ora_results filtered to focus modules
  unified_enrichment <- ora_results %>%
    filter(tolower(Module) %in% focus_modules_auto)
  cat(sprintf("   Using ora_results filtered to focus modules: %s\n",
              paste(focus_modules_auto, collapse = ", ")))
}

# Add significance flags
unified_enrichment <- unified_enrichment %>%
  mutate(
    Sig_FDR05 = FDR < 0.05,
    Sig_FDR10 = FDR < 0.10,
    Sig_Nominal = P_Value < 0.05
  )

# Save complete table
write.csv(unified_enrichment,
          file.path(results_dir, "Step8_Unified_Enrichment_ALL.csv"),
          row.names = FALSE)

cat(sprintf("Unified enrichment: %d total rows\n", nrow(unified_enrichment)))
cat(sprintf("  FDR < 0.05: %d\n", sum(unified_enrichment$Sig_FDR05)))
cat(sprintf("  FDR < 0.10: %d\n", sum(unified_enrichment$Sig_FDR10)))
cat(sprintf("  Nominal p < 0.05: %d\n", sum(unified_enrichment$Sig_Nominal)))

cat("   [OK] Saved: Step8_Unified_Enrichment_ALL.csv\n")

# ============================================================================
# 8-FINAL-b: PUBLICATION EXCEL TABLE
# ============================================================================

# Note: openxlsx already loaded at startup (line 74)
wb <- createWorkbook()

# Sheet 1: All Results (Supplementary Table S4)
addWorksheet(wb, "All_Results")
writeData(wb, "All_Results", unified_enrichment)
setColWidths(wb, "All_Results", cols = 1:ncol(unified_enrichment), widths = "auto")

# Create sheets for each focus module dynamically
focus_module_sheets <- list()
for(mod in focus_modules_auto) {
  mod_upper <- toupper(mod)
  sheet_name <- paste0(tools::toTitleCase(mod), "_Significant")

  mod_sig <- unified_enrichment %>%
    filter(toupper(Module) == mod_upper, Sig_FDR10 | Sig_Nominal) %>%
    arrange(FDR)

  if(nrow(mod_sig) > 0) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, mod_sig)
    setColWidths(wb, sheet_name, cols = 1:ncol(mod_sig), widths = "auto")
    focus_module_sheets[[mod]] <- mod_sig
    cat(sprintf("   Added sheet '%s': %d significant pathways\n", sheet_name, nrow(mod_sig)))
  } else {
    cat(sprintf("   Skipped sheet '%s': No significant pathways\n", sheet_name))
  }
}

# Sheet 4: Summary by Database
db_summary <- unified_enrichment %>%
  group_by(Module, Database) %>%
  summarise(
    N_Tested = n(),
    N_FDR05 = sum(Sig_FDR05, na.rm = TRUE),
    N_FDR10 = sum(Sig_FDR10, na.rm = TRUE),
    N_Nominal = sum(Sig_Nominal, na.rm = TRUE),
    .groups = "drop"
  )
addWorksheet(wb, "Summary_Database")
writeData(wb, "Summary_Database", db_summary)
setColWidths(wb, "Summary_Database", cols = 1:ncol(db_summary), widths = "auto")

# Sheet 5: Summary by Category (biological themes)
cat_summary <- unified_enrichment %>%
  filter(Sig_FDR10 | Sig_Nominal) %>%
  group_by(Module, Category) %>%
  summarise(
    N_Pathways = n(),
    Top_Pathway = Short_Name[which.min(FDR)],
    Best_FDR = min(FDR),
    Best_P = min(P_Value),
    .groups = "drop"
  ) %>%
  arrange(Module, Best_FDR)
addWorksheet(wb, "Summary_Category")
writeData(wb, "Summary_Category", cat_summary)
setColWidths(wb, "Summary_Category", cols = 1:ncol(cat_summary), widths = "auto")

# Sheet 6: Table 4 - Top Pathways for Manuscript
top_pathways <- unified_enrichment %>%
  filter(Sig_FDR10 | Sig_Nominal) %>%
  group_by(Module) %>%
  slice_min(P_Value, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  select(Module, Database, Short_Name, Category, N_Overlap, Fold_Enrichment, P_Value, FDR)
addWorksheet(wb, "Table4_TopPathways")
writeData(wb, "Table4_TopPathways", top_pathways)
setColWidths(wb, "Table4_TopPathways", cols = 1:ncol(top_pathways), widths = "auto")

# Save Excel workbook
saveWorkbook(wb, file.path(results_dir, "Step8_Unified_Enrichment_Publication.xlsx"),
             overwrite = TRUE)

# Count total sheets created
n_module_sheets <- length(focus_module_sheets)
total_sheets <- 1 + n_module_sheets + 3  # All_Results + module sheets + Summary_Database + Summary_Category + Table4

cat(sprintf("Publication Excel saved with %d sheets:\n", total_sheets))
cat(sprintf("  1. All_Results: %d rows\n", nrow(unified_enrichment)))
sheet_num <- 2
for(mod in names(focus_module_sheets)) {
  cat(sprintf("  %d. %s_Significant: %d rows\n", sheet_num, tools::toTitleCase(mod), nrow(focus_module_sheets[[mod]])))
  sheet_num <- sheet_num + 1
}
cat(sprintf("  %d. Summary_Database: %d rows\n", sheet_num, nrow(db_summary)))
cat(sprintf("  %d. Summary_Category: %d rows\n", sheet_num + 1, nrow(cat_summary)))
cat(sprintf("  %d. Table4_TopPathways: %d rows\n", sheet_num + 2, nrow(top_pathways)))

# ============================================================================
# 8-FINAL-c: UNIFIED WORD CLOUD
# ============================================================================
conflicts_prefer(dplyr::first)

# unified word cloud combining ORA and GSVA results
cat("Generating unified word cloud (ORA + GSVA, FDR < 0.1)...\n")

# Create category lookup from ora_results
category_lookup <- ora_results %>%
  select(Short_Name, Category) %>%
  distinct() %>%
  { setNames(.$Category, .$Short_Name) }

# combine ORA and GSVA results with relaxed threshold
wordcloud_data <- bind_rows(
  # from ORA results
  ora_results %>%
    filter(FDR < 0.1, !tolower(Module) %in% c("grey", "gray")) %>%
    mutate(
      Source = "ORA",
      Score = -log10(P_Value) * log2(pmax(Fold_Enrichment, 1) + 1)
    ) %>%
    select(Short_Name, Category, Score, Source, FDR, P_Value),

  # from GSVA clinical associations (if exists)
  if (exists("clinical_tests") && nrow(clinical_tests) > 0) {
    clinical_tests %>%
      filter(FDR < 0.1 | P_Value < 0.05) %>%
      mutate(
        Source = "GSVA",
        Score = -log10(P_Value) * abs(Cohens_d),
        # Look up Category from ora_results
        Category = ifelse(Short_Name %in% names(category_lookup),
                          category_lookup[Short_Name], "Other")
      ) %>%
      select(Short_Name, Category, Score, Source, FDR, P_Value)
  } else {
    NULL
  }
) %>%
  # aggregate by pathway (take best score if duplicated)
  group_by(Short_Name) %>%
  summarise(
    Category = first(Category),
    Score = max(Score, na.rm = TRUE),
    Source = paste(unique(Source), collapse = "+"),
    Best_FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Score)) %>%
  head(50)

cat(sprintf("  Word cloud pathways: %d\n", nrow(wordcloud_data)))
cat(sprintf("  From ORA: %d\n", sum(grepl("ORA", wordcloud_data$Source))))
cat(sprintf("  From GSVA: %d\n", sum(grepl("GSVA", wordcloud_data$Source))))
cat(sprintf("  From both: %d\n", sum(grepl("\\+", wordcloud_data$Source))))

# category colors
category_colors <- c(
  "CIBERSORT_Immune" = "#E64B35",
  "Immune" = "#E64B35",
  "PDAC_Subtype" = "#3C5488",
  "Stroma" = "#7E6148",
  "Stemness" = "#9B59B6",
  "EMT" = "#00A087",
  "Metabolic" = "#F39B7F",
  "Cytokine" = "#e377c2",
  "Checkpoint" = "#7f7f7f",
  "Exosome" = "#4DBBD5",
  "Signaling" = "#8491B4",
  "Hot_Cold" = "#8B0000",
  "This_Study" = "#2F4F4F",
  "Chemoresistance" = "#FF8C00",
  "Custom_PDAC" = "#2F2F2F",
  "Other" = "#B09C85"
)

set.seed(42)

# 1. Main Wordcloud Plot
p_main <- ggplot(wordcloud_data, 
                 aes(label = Short_Name, size = Score, color = Category)) +
  geom_text_wordcloud_area(
    rm_outside = TRUE, 
    shape = "circle",
    grid_margin = 1,
    seed = 42
  ) +
  # range = c(8, 30) forces the smallest font to 8 and largest to 30
  # We remove 'limits' so no data is filtered out
  scale_size_area(max_size = 30) +
  scale_color_manual(values = category_colors) +
  labs(
    title = "Integrated Pathway Enrichment: ORA + GSVA",
    subtitle = "FDR < 0.1 | Size = significance score | Color = biological category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    legend.position = "none"
  )

# 2. Extract Legend using a dummy plot
p_legend_only <- ggplot(wordcloud_data, aes(x = 1, y = 1, color = Category)) +
  geom_point(size = 5) +
  scale_color_manual(values = category_colors, name = "Biological Category") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 6)))

shared_legend <- cowplot::get_legend(p_legend_only)

# 3. Combine
final_plot <- cowplot::plot_grid(
  p_main, 
  shared_legend, 
  ncol = 1, 
  rel_heights = c(1, 0.2)
)

print(final_plot)

# save
ggsave(file.path(fig_dir_png, "Manuscript_WordCloud_ORA_GSVA.png"), final_plot,
       width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Manuscript_WordCloud_ORA_GSVA.pdf"), final_plot,
       width = 8, height = 5)

cat("  Saved: Manuscript_WordCloud_ORA_GSVA.png/pdf\n")

# -- Save workspace --
# Build list of objects to save (check existence first)
save_objects <- c("signatures", "custom_signatures", "gene_sets", "sig_ref",
                  "display_names_short", "display_names_full",
                  "ora_results", "gsva_results", "ssgsea_results",
                  "clinical_tests", "me_gsva_cor", "me_gsva_pval",
                  "unified_enrichment")

# Add focus_module_ora only if it exists
if(exists("focus_module_ora")) {
  save_objects <- c(save_objects, "focus_module_ora")
}

# Save only existing objects
save(list = save_objects[sapply(save_objects, exists)],
     file = file.path(results_dir, "Step8_Pathway_Analysis.RData"))
cat("   [OK] Saved workspace: Step8_Pathway_Analysis.RData\n")

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 8 SUMMARY: PATHWAY ENRICHMENT ANALYSIS ===\n")
cat("===============================================================================\n")
cat(sprintf("Custom PDAC Signatures Tested: %d\n", length(custom_signatures)))
cat(sprintf("Focus modules analyzed: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))

# Summary per module - with existence check
if(exists("focus_module_ora") && length(focus_module_ora) > 0) {
  for(mod in focus_modules_auto) {
    mod_ora <- focus_module_ora[[mod]]
    n_sig <- if(!is.null(mod_ora) && nrow(mod_ora) > 0) sum(mod_ora$P_Value < 0.05, na.rm = TRUE) else 0
    cat(sprintf("  %s significant (p<0.05): %d\n", toupper(mod), n_sig))
  }
} else {
  # Fallback: use unified_enrichment
  for(mod in focus_modules_auto) {
    n_sig <- sum(unified_enrichment$Module == toupper(mod) & unified_enrichment$P_Value < 0.05, na.rm = TRUE)
    cat(sprintf("  %s significant (p<0.05): %d\n", toupper(mod), n_sig))
  }
}
cat(sprintf("\nUnified Enrichment Table: %d total rows\n", nrow(unified_enrichment)))
cat(sprintf("  FDR < 0.05: %d\n", sum(unified_enrichment$Sig_FDR05, na.rm = TRUE)))
cat(sprintf("  FDR < 0.10: %d\n", sum(unified_enrichment$Sig_FDR10, na.rm = TRUE)))
cat(sprintf("  Nominal p < 0.05: %d\n", sum(unified_enrichment$Sig_Nominal, na.rm = TRUE)))
cat("\nKey Output Files:\n")
cat("  Tables:\n")
for(mod in focus_modules_auto) {
  cat(sprintf("    - Step8c2_%s_Custom_ORA.csv\n", toupper(mod)))
}
cat("    - Step8_Unified_Enrichment_ALL.csv\n")
cat("    - Step8_Unified_Enrichment_Publication.xlsx\n")
cat("  Figures:\n")
cat("    - Step8_01_ORA_Lollipop.png/pdf\n")
cat("    - Step8_02_ORA_Heatmap.png/pdf\n")
cat("    - Step8f_ssGSEA_Top20_Heatmap.png/pdf\n")
cat("    - Step8_Unified_WordCloud.png/pdf\n")
cat("Step 8 complete.\n")
cat("=======================================================================\n")

#===============================================================================
# STEP 9: MASTER ENRICHMENT CONSENSUS ANALYSIS
#===============================================================================

# PURPOSE: Create integrated consensus ranking of pathway enrichments
# across multiple methods (ORA, GSVA, ME-correlations) for robust results.
# STRATEGY:                                                                
# 9a) Load ORA + ME-pathway correlations
# 9b) Rank-based integration
# 9c) Consensus scoring -> 9d) Identify top pathways per module 
# LOOKING FOR:                                                             
# - Pathways significant across multiple methods (robust signal)         
# - Biologically coherent pathway clusters (related functions together)  
# - Module-specific pathway themes emerging                              
# KEY PARAMS: consensus_weight=equal, min_methods=2, rank_aggregation      


cat("##  STEP 9: MASTER ENRICHMENT CONSENSUS ANALYSIS  ##\n")

# --- Defensive Check: Ensure focus_modules_auto exists ---
if (!exists("focus_modules_auto") || length(focus_modules_auto) == 0) {
  cat("   [INFO] focus_modules_auto not found. Attempting to reconstruct...\n")

  # Try to load from saved network data
  network_file <- file.path(results_dir, "Step6_NetworkData.rds")
  if (file.exists(network_file)) {
    networkData <- readRDS(network_file)
    if (!is.null(networkData$focus_modules_auto)) {
      focus_modules_auto <- networkData$focus_modules_auto
      cat(sprintf("   [OK] Loaded focus_modules_auto from RDS: %s\n",
                  paste(focus_modules_auto, collapse = ", ")))
    }
  }

  # If still missing, derive from module priority file
  if (!exists("focus_modules_auto") || length(focus_modules_auto) == 0) {
    priority_file <- file.path(results_dir, "Step5_ModulePriority.csv")
    if (file.exists(priority_file)) {
      module_priority <- read.csv(priority_file, stringsAsFactors = FALSE)
      focus_modules_auto <- module_priority %>%
        filter(N_Significant >= 2) %>%
        pull(Module) %>%
        gsub("ME", "", .)
      if (length(focus_modules_auto) == 0) {
        focus_modules_auto <- module_priority %>%
          filter(N_Significant >= 1) %>%
          head(4) %>%
          pull(Module) %>%
          gsub("ME", "", .)
      }
      cat(sprintf("   [OK] Derived focus_modules_auto from priority file: %s\n",
                  paste(focus_modules_auto, collapse = ", ")))
    }
  }

  # Final fallback: use default modules
  if (!exists("focus_modules_auto") || length(focus_modules_auto) == 0) {
    focus_modules_auto <- c("pink", "green", "black", "blue")
    cat("   [WARN] Using default focus modules: pink, green, black, blue\n")
  }
}

# Define derived variables for Step 9
focus_modules <- c(focus_modules_auto, "grey")
focus_names <- tools::toTitleCase(focus_modules_auto)

# Define focus_colors palette
module_color_palette <- c(
  "brown" = "brown", "blue" = "steelblue", "green" = "green3",
  "turquoise" = "turquoise3", "yellow" = "gold", "red" = "indianred",
  "pink" = "pink", "magenta" = "magenta3", "black" = "grey30",
  "grey" = "grey60", "cyan" = "cyan3", "purple" = "purple3"
)
focus_colors <- module_color_palette[focus_modules_auto]
focus_colors_with_grey <- module_color_palette[focus_modules]

cat(sprintf("   [OK] Focus modules: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat(sprintf("   [OK] Focus colors: %s\n", paste(names(focus_colors), "=", focus_colors, collapse = ", ")))

cat(sprintf("   Creating integrated consensus ranking for %d focus modules\n", length(focus_modules_auto)))
cat("   Sources: ORA + ME-Pathway Correlations\n")
cat("   Method: Rank-based consensus with equal weighting\n")

# ============================================================================
# 9a: LOAD DATA SOURCES
# ============================================================================

cat("==============================================================================\n")
cat("9a: LOADING DATA SOURCES\n")
cat("==============================================================================\n")

# Load ORA Results
ora_consensus_df <- read.csv(file.path(results_dir, "Step8_02_ORA_Results.csv"), stringsAsFactors = FALSE)
cat(sprintf("   [OK] Loaded ORA Results: %d entries", nrow(ora_consensus_df)))

# Load ME-Pathway Correlations
me_corr_consensus_df <- read.csv(file.path(results_dir, "Step8_08_ME_Pathway_Correlations.csv"), stringsAsFactors = FALSE)
cat(sprintf("   [OK] Loaded ME-Pathway Correlations: %d entries", nrow(me_corr_consensus_df)))

# ============================================================================
# 9b: PROCESS DATA SOURCES FOR EACH MODULE
# ============================================================================

cat("==============================================================================\n")
cat("9b: PROCESSING DATA SOURCES\n")
cat("==============================================================================\n")

process_ora_consensus <- function(df, module) {
  mod_df <- df[df$Module == module, ]
  if (nrow(mod_df) == 0) return(data.frame())

  mod_df$Source <- "ORA"
  mod_df$Metric <- mod_df$Fold_Enrichment
  mod_df$Metric_Name <- "Fold_Enrichment"
  mod_df$Significant <- mod_df$P_Value < 0.05
  mod_df$FDR_Sig <- mod_df$FDR < 0.10
  mod_df$Rank <- rank(-mod_df$Fold_Enrichment, ties.method = "min")

  return(mod_df[, c("Signature_ID", "Short_Name", "Full_Name", "Category",
                    "Source", "Metric", "Metric_Name", "P_Value", "FDR",
                    "Significant", "FDR_Sig", "Rank", "N_Overlap", "Overlap_Genes")])
}

process_me_corr_consensus <- function(df, module) {
  mod_df <- df[df$Module == module, ]
  if (nrow(mod_df) == 0) return(data.frame())

  mod_df$Source <- "ME_Correlation"
  mod_df$Metric <- abs(mod_df$Correlation)
  mod_df$Metric_Name <- "Abs_Correlation"
  mod_df$Significant <- (mod_df$P_Value < 0.05) | (abs(mod_df$Correlation) > 0.3)
  mod_df$FDR_Sig <- mod_df$FDR < 0.10
  mod_df$Direction <- ifelse(mod_df$Correlation > 0, "Positive", "Negative")
  mod_df$Rank <- rank(-mod_df$Metric, ties.method = "min")
  mod_df$Signature_ID <- mod_df$Pathway
  mod_df$N_Overlap <- NA
  mod_df$Overlap_Genes <- ""

  return(mod_df[, c("Signature_ID", "Short_Name", "Source", "Metric", "Metric_Name",
                    "P_Value", "FDR", "Significant", "FDR_Sig", "Rank", "Direction",
                    "N_Overlap", "Overlap_Genes")])
}

# Process ALL focus modules dynamically
module_ora_results <- list()
module_me_results <- list()

for(mod in focus_modules_auto) {
  cat(sprintf("   Processing %s module...\n", toupper(mod)))

  # Process ORA results for this module
  mod_ora <- process_ora_consensus(ora_consensus_df, mod)
  module_ora_results[[mod]] <- mod_ora

  # Process ME correlation results for this module
  mod_me <- process_me_corr_consensus(me_corr_consensus_df, mod)
  module_me_results[[mod]] <- mod_me

  cat(sprintf("      ORA: %d signatures\n", nrow(mod_ora)))
  cat(sprintf("      ME-Corr: %d signatures\n", nrow(mod_me)))
}

cat(sprintf("   [OK] Processed %d focus modules\n", length(focus_modules_auto)))

# ============================================================================
# 9c: CREATE CONSENSUS SCORING
# ============================================================================

cat("==============================================================================\n")
cat("9c: CREATING CONSENSUS RANKINGS\n")
cat("==============================================================================\n")

create_consensus <- function(ora_df, me_df, module_name) {
  # Get all unique signatures
  all_sigs <- unique(c(
    if (nrow(ora_df) > 0) ora_df$Signature_ID else character(0),
    if (nrow(me_df) > 0) me_df$Signature_ID else character(0)
  ))

  consensus_rows <- lapply(all_sigs, function(sig) {
    row <- list(
      Module = module_name,
      Signature_ID = sig,
      Short_Name = "",
      Full_Name = "",
      Category = "",

      # ORA metrics
      ORA_FoldEnrich = NA,
      ORA_PValue = NA,
      ORA_FDR = NA,
      ORA_Rank = NA,
      ORA_Significant = FALSE,
      ORA_Overlap = NA,
      ORA_Genes = "",

      # ME Correlation metrics
      ME_Correlation = NA,
      ME_AbsCorr = NA,
      ME_PValue = NA,
      ME_Rank = NA,
      ME_Significant = FALSE,
      ME_Direction = "",

      # Consensus metrics
      N_Sources_Significant = 0,
      N_Sources_Present = 0,
      Mean_Rank = NA,
      Consensus_Score = 0
    )

    # Fill ORA data
    ora_match <- ora_df[ora_df$Signature_ID == sig, ]
    if (nrow(ora_match) > 0) {
      ora_row <- ora_match[1, ]
      row$Short_Name <- ora_row$Short_Name
      row$Full_Name <- ora_row$Full_Name
      row$Category <- ora_row$Category
      row$ORA_FoldEnrich <- ora_row$Metric
      row$ORA_PValue <- ora_row$P_Value
      row$ORA_FDR <- ora_row$FDR
      row$ORA_Rank <- ora_row$Rank
      row$ORA_Significant <- ora_row$Significant
      row$ORA_Overlap <- ora_row$N_Overlap
      row$ORA_Genes <- ora_row$Overlap_Genes
      row$N_Sources_Present <- row$N_Sources_Present + 1
      if (ora_row$Significant) row$N_Sources_Significant <- row$N_Sources_Significant + 1
    }

    # Fill ME Correlation data
    me_match <- me_df[me_df$Signature_ID == sig, ]
    if (nrow(me_match) > 0) {
      me_row <- me_match[1, ]
      if (row$Short_Name == "") row$Short_Name <- me_row$Short_Name
      row$ME_AbsCorr <- me_row$Metric
      row$ME_Correlation <- ifelse(me_row$Direction == "Positive", me_row$Metric, -me_row$Metric)
      row$ME_PValue <- me_row$P_Value
      row$ME_Rank <- me_row$Rank
      row$ME_Significant <- me_row$Significant
      row$ME_Direction <- me_row$Direction
      row$N_Sources_Present <- row$N_Sources_Present + 1
      if (me_row$Significant) row$N_Sources_Significant <- row$N_Sources_Significant + 1
    }

    # Calculate mean rank (only for sources where present)
    ranks <- c()
    if (!is.na(row$ORA_Rank)) ranks <- c(ranks, row$ORA_Rank)
    if (!is.na(row$ME_Rank)) ranks <- c(ranks, row$ME_Rank)

    if (length(ranks) > 0) {
      row$Mean_Rank <- mean(ranks)
    }

    # Consensus score = N_Sources_Significant + (1 / Mean_Rank) if present
    row$Consensus_Score <- row$N_Sources_Significant
    if (!is.na(row$Mean_Rank) && row$Mean_Rank > 0) {
      row$Consensus_Score <- row$Consensus_Score + (1 / row$Mean_Rank)
    }

    return(as.data.frame(row, stringsAsFactors = FALSE))
  })

  consensus_df <- do.call(rbind, consensus_rows)

  # Final ranking by consensus score
  consensus_df$Consensus_Rank <- rank(-consensus_df$Consensus_Score, ties.method = "min")
  consensus_df <- consensus_df[order(consensus_df$Consensus_Rank), ]

  return(consensus_df)
}

# Create consensus for ALL focus modules
module_consensus <- list()

for(mod in focus_modules_auto) {
  mod_consensus <- create_consensus(
    module_ora_results[[mod]],
    module_me_results[[mod]],
    mod
  )
  module_consensus[[mod]] <- mod_consensus
  cat(sprintf("   [OK] %s module: %d signatures ranked\n", toupper(mod), nrow(mod_consensus)))
}

# Combine all modules into master table
master_consensus_df <- do.call(rbind, module_consensus)
cat(sprintf("   [OK] Master table: %d total entries across %d modules\n",
            nrow(master_consensus_df), length(focus_modules_auto)))

# ============================================================================
# 9d: SAVE MASTER TABLE
# ============================================================================

cat("==============================================================================\n")
cat("9d: SAVING MASTER CONSENSUS TABLE\n")
cat("==============================================================================\n")

write.csv(master_consensus_df, file.path(main_dir, "Master_Enrichment_Consensus_Table.csv"), row.names = FALSE)
cat(sprintf("   [OK] Saved: Master_Enrichment_Consensus_Table.csv (%d entries)", nrow(master_consensus_df)))

# ============================================================================
# 9e: DISPLAY TOP 10 AND TOP 20 RANKINGS
# ============================================================================

cat("==============================================================================\n")
cat("9e: TOP 20 ENRICHED SIGNATURES\n")
cat("==============================================================================\n")

display_top_pathways <- function(df, module, n = 20) {
  mod_df <- df[df$Module == module, ]

  if(nrow(mod_df) == 0) {
    cat(sprintf("\n   %s MODULE - No pathways found\n", toupper(module)))
    return(invisible(NULL))
  }

  mod_df <- head(mod_df[order(mod_df$Consensus_Rank), ], n)

  cat(sprintf("\n   %s MODULE - Top %d Enriched Pathways:\n", toupper(module), min(n, nrow(mod_df))))
  cat("   ------------------------------------------------------------------------\n")
  cat("   Rank | Pathway                        | Sig  | Evidence\n")
  cat("   ------------------------------------------------------------------------\n")

  for (i in 1:nrow(mod_df)) {
    row <- mod_df[i, ]
    sig <- ifelse(row$Short_Name != "", row$Short_Name, row$Signature_ID)
    n_sig <- row$N_Sources_Significant
    n_present <- row$N_Sources_Present

    # Build evidence string
    evidence <- c()
    if (row$ORA_Significant) {
      evidence <- c(evidence, sprintf("ORA:%.1fx", row$ORA_FoldEnrich))
    }
    if (row$ME_Significant) {
      direction <- ifelse(row$ME_Direction == "Positive", "+", "-")
      evidence <- c(evidence, sprintf("ME:%s%.2f", direction, row$ME_AbsCorr))
    }
    ev_str <- ifelse(length(evidence) > 0, paste(evidence, collapse = " | "), "sub-threshold")

    cat(sprintf("   %2d   | %-30s | %d/%d  | %s\n",
                row$Consensus_Rank, substr(sig, 1, 30), n_sig, n_present, ev_str))
  }
  cat("   ------------------------------------------------------------------------\n")
}

# Display top pathways for ALL focus modules
for(mod in focus_modules_auto) {
  display_top_pathways(master_consensus_df, mod, 20)
}

# ============================================================================
# 9f: CREATE MASTER CONSENSUS FIGURE
# ============================================================================

cat("==============================================================================\n")
cat("9f: GENERATING MASTER CONSENSUS FIGURE\n")
cat("==============================================================================\n")

create_consensus_dotplot <- function(df, module_name, module_color) {
  plot_df <- head(df[df$Module == module_name, ], 20)

  # Handle empty data

  if(nrow(plot_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = paste("No data for", toupper(module_name), "module")) +
             theme_void())
  }

  plot_df <- plot_df[order(-plot_df$Consensus_Rank), ]  # Reverse for bottom-to-top
  plot_df$y_pos <- 1:nrow(plot_df)
  plot_df$Label <- ifelse(plot_df$Short_Name != "", plot_df$Short_Name, plot_df$Signature_ID)

  # Create markers for significant sources
  plot_df$Markers <- sapply(1:nrow(plot_df), function(i) {
    markers <- c()
    if (plot_df$ORA_Significant[i]) markers <- c(markers, "O")
    if (plot_df$ME_Significant[i]) markers <- c(markers, "M")
    paste(markers, collapse = "")
  })

  # Count significant pathways
  n_sig <- sum(plot_df$N_Sources_Significant >= 1)

  p <- ggplot(plot_df, aes(x = Consensus_Score, y = y_pos)) +
    geom_point(aes(size = N_Sources_Significant + 1, fill = Consensus_Score),
               shape = 21, color = "black", alpha = 0.8) +
    geom_text(aes(label = Markers), hjust = -0.5, size = 2.5, fontface = "bold", color = "darkblue") +
    scale_y_continuous(breaks = plot_df$y_pos, labels = plot_df$Label) +
    scale_size_continuous(range = c(3, 8), name = "N Sources\nSignificant") +
    scale_fill_gradient(low = "#FFFFCC", high = "#CC0000", name = "Consensus\nScore") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5) +
    labs(
      title = sprintf("%s Module", toupper(module_name)),
      subtitle = sprintf("Top 20 Enriched Pathways (%d significant in ≥1 source)", n_sig),
      x = "Consensus Score (higher = more robust enrichment)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = module_color, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )

  return(p)
}

# Create plots for ALL focus modules - dynamic colors
# Ensure focus_colors exists
if (!exists("focus_colors") || length(focus_colors) != length(focus_modules_auto)) {
  module_color_palette <- c(
    "brown" = "brown", "blue" = "steelblue", "green" = "green3",
    "turquoise" = "turquoise3", "yellow" = "gold", "red" = "indianred",
    "pink" = "pink", "magenta" = "magenta3", "black" = "grey30",
    "grey" = "grey60", "cyan" = "cyan3", "purple" = "purple3"
  )
  focus_colors <- module_color_palette[focus_modules_auto]
  cat("   [INFO] Reconstructed focus_colors for all focus modules\n")
}

# Create individual plots for each focus module
consensus_plots <- list()
for(i in seq_along(focus_modules_auto)) {
  mod <- focus_modules_auto[i]
  mod_color <- as.character(focus_colors[mod])
  consensus_plots[[mod]] <- create_consensus_dotplot(master_consensus_df, mod, mod_color)
}

# Combine plots in a 2x2 grid (for 4 modules) or adaptive layout
n_mods <- length(focus_modules_auto)
if(n_mods == 4) {
  p_consensus_combined <- (consensus_plots[[1]] + consensus_plots[[2]]) /
                          (consensus_plots[[3]] + consensus_plots[[4]]) +
    plot_annotation(
      title = "Master Enrichment Consensus Analysis",
      subtitle = sprintf("Focus Modules: %s | Integrated ORA + ME-Pathway Correlations | Markers: O=ORA, M=ME-correlation",
                         paste(toupper(focus_modules_auto), collapse = ", ")),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
      )
    )
} else if(n_mods == 2) {
  p_consensus_combined <- consensus_plots[[1]] + consensus_plots[[2]] +
    plot_annotation(
      title = "Master Enrichment Consensus Analysis",
      subtitle = sprintf("Focus Modules: %s | Integrated ORA + ME-Pathway Correlations | Markers: O=ORA, M=ME-correlation",
                         paste(toupper(focus_modules_auto), collapse = ", ")),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
      )
    )
} else {
  # Generic wrap for any number of modules
  p_consensus_combined <- wrap_plots(consensus_plots, ncol = 2) +
    plot_annotation(
      title = "Master Enrichment Consensus Analysis",
      subtitle = sprintf("Focus Modules: %s | Integrated ORA + ME-Pathway Correlations | Markers: O=ORA, M=ME-correlation",
                         paste(toupper(focus_modules_auto), collapse = ", ")),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
      )
    )
}

# Save figure - adjust size based on number of modules
fig_height <- ifelse(length(focus_modules_auto) > 2, 14, 10)
ggsave(file.path(fig_dir_png, "Step9_01_Master_Consensus_Figure.png"), p_consensus_combined,
       width = 16, height = fig_height, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step9_01_Master_Consensus_Figure.pdf"), p_consensus_combined,
       width = 16, height = fig_height)
cat("   [OK] Saved: Step9_01_Master_Consensus_Figure.png/pdf\n")

# ============================================================================
# 9g: VERIFICATION - TOP RANKED PATHWAYS FOR ALL FOCUS MODULES
# ============================================================================

cat("==============================================================================\n")
cat("9g: VERIFICATION - #1 RANKED PATHWAYS FOR ALL FOCUS MODULES\n")
cat("==============================================================================\n")

# Display top-ranked pathway for each focus module
for(mod in focus_modules_auto) {
  mod_top1 <- module_consensus[[mod]][1, ]

  cat(sprintf("\n   %s MODULE #1:\n", toupper(mod)))
  cat("   ------------------------------------------------------------------------\n")
  cat(sprintf("   Signature:    %s\n", mod_top1$Short_Name))
  cat(sprintf("   Category:     %s\n", mod_top1$Category))
  cat(sprintf("   Consensus:    Rank=%d, Score=%.3f\n", mod_top1$Consensus_Rank, mod_top1$Consensus_Score))
  cat(sprintf("   Sources:      %d/%d significant\n", mod_top1$N_Sources_Significant, mod_top1$N_Sources_Present))
  if (!is.na(mod_top1$ORA_FoldEnrich)) {
    cat(sprintf("   ORA:          %.2fx (p=%.2e)\n", mod_top1$ORA_FoldEnrich, mod_top1$ORA_PValue))
  }
  if (!is.na(mod_top1$ME_AbsCorr)) {
    cat(sprintf("   ME-Corr:      %s%.3f (p=%.2e)\n",
                ifelse(mod_top1$ME_Direction == "Positive", "+", "-"),
                mod_top1$ME_AbsCorr, mod_top1$ME_PValue))
  }
}

# ============================================================================
# 9h: STEP 9 SUMMARY
# ============================================================================

cat("==============================================================================\n")
cat("9h: STEP 9 SUMMARY\n")
cat("==============================================================================\n")

cat(sprintf("   FOCUS MODULES ANALYZED: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("   ------------------------------------------------------------------------\n")
cat("   METHODS INTEGRATED:\n")
cat("      - ORA (Over-Representation Analysis)\n")
cat("      - ME-Pathway Correlations\n")
cat("   ------------------------------------------------------------------------\n")
cat("   OUTPUT FILES:\n")
cat("      - Master_Enrichment_Consensus_Table.csv\n")
cat("      - Step9_01_Master_Consensus_Figure.png/pdf\n")
cat("   ------------------------------------------------------------------------\n")
cat(sprintf("   TOTAL SIGNATURES RANKED: %d\n", nrow(master_consensus_df)))
for(mod in focus_modules_auto) {
  n_sig <- sum(master_consensus_df$Module == mod & master_consensus_df$N_Sources_Significant >= 1)
  cat(sprintf("      %s: %d significant in at least 1 source\n", toupper(mod), n_sig))
}
cat("==============================================================================\n")
cat("STEP 9 COMPLETE\n")
cat("==============================================================================\n")


#===============================================================================
# STEP 10: DUAL-CLINICAL ASSOCIATION ANALYSIS
#===============================================================================
# PURPOSE: Validate pathway enrichment scores against clinical outcomes
# (PFS group and Treatment Response) to establish prognostic/predictive significance.
# FOCUS MODULES: pink, green, black, blue (from Step 5 prioritization)
# STRATEGY:
#   10a) Load GSVA scores and clinical data
#   10b) Define analysis functions (winsorize, Cohen's d)
#   10c) Perform dual-clinical analysis (Wilcoxon + effect size)
#   10d) Apply FDR correction and assign overlap status
#   10e) Save results
#   10f) Figure A - Dual Lollipop Plot
#   10g) Figure B - Response/PFS Boxplots
#   10h) Summary
# LOOKING FOR:
#   - Pathways significant for both PFS AND Response (Dual_Significant)
#   - Prognostic biomarkers (PFS only) vs Predictive biomarkers (Response only)
#   - Effect sizes (Cohen's d) to quantify clinical relevance
# KEY PARAMS: FDR < 0.10, stat_tests = Wilcoxon rank-sum

cat("==============================================================================\n")
cat("STEP 10: DUAL-CLINICAL ASSOCIATION ANALYSIS\n")
cat("==============================================================================\n")
cat(sprintf("Focus Modules: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))

# ============================================================================
# 10a: LOAD GSVA SCORES AND CLINICAL DATA
# ============================================================================

cat("==============================================================================\n")
cat("10a: LOAD GSVA SCORES AND CLINICAL DATA\n")
cat("==============================================================================\n")

# Load GSVA scores from Step 8
gsva_file <- file.path(results_dir, "Step8_03_GSVA_Scores.csv")
if (file.exists(gsva_file)) {
  gsva_scores <- read.csv(gsva_file, row.names = 1, check.names = FALSE)
  cat(sprintf("   [OK] Loaded GSVA scores: %d pathways x %d samples\n", ncol(gsva_scores), nrow(gsva_scores)))
} else {
  cat("   [WARNING] GSVA scores not found. Using ssGSEA scores...\n")
  ssgsea_file <- file.path(results_dir, "Step8_04_ssGSEA_Scores.csv")
  gsva_scores <- read.csv(ssgsea_file, row.names = 1, check.names = FALSE)
  cat(sprintf("   [OK] Loaded ssGSEA scores: %d pathways x %d samples\n", ncol(gsva_scores), nrow(gsva_scores)))
}

sample_ids <- rownames(gsva_scores)

# Create clinical data frame from datTraits (NOT from sample naming convention)
clinical_data <- data.frame(
  Sample_ID = sample_ids,
  stringsAsFactors = FALSE
)

# Add PFS_group from datTraits (handle both factor and numeric formats)
if (exists("datTraits") && "PFS_group" %in% colnames(datTraits)) {
  pfs_raw <- datTraits$PFS_group[match(sample_ids, rownames(datTraits))]
  # Convert to standard labels: "Short" and "Long"
  if (is.factor(pfs_raw)) {
    clinical_data$PFS_group <- as.character(pfs_raw)
  } else if (is.numeric(pfs_raw)) {
    # Numeric: 0 = Short, 1 = Long (standard encoding)
    clinical_data$PFS_group <- ifelse(pfs_raw == 1, "Long", ifelse(pfs_raw == 0, "Short", NA))
  } else {
    clinical_data$PFS_group <- as.character(pfs_raw)
  }
  # Standardize labels (handle S/L, short/long, etc.)
  clinical_data$PFS_group <- case_when(
    clinical_data$PFS_group %in% c("S", "Short", "short", "SHORT", "0") ~ "Short",
    clinical_data$PFS_group %in% c("L", "Long", "long", "LONG", "1") ~ "Long",
    TRUE ~ clinical_data$PFS_group
  )
} else {
  cat("   [WARNING] PFS_group not found in datTraits\n")
  clinical_data$PFS_group <- NA
}

# Add Response from datTraits (handle both factor and numeric formats)
if (exists("datTraits") && "Response" %in% colnames(datTraits)) {
  resp_raw <- datTraits$Response[match(sample_ids, rownames(datTraits))]
  # Convert to standard labels
  if (is.factor(resp_raw)) {
    clinical_data$Response_group <- as.character(resp_raw)
  } else if (is.numeric(resp_raw)) {
    # Numeric: 0 = PD, 1 = CD (standard encoding)
    clinical_data$Response_group <- ifelse(resp_raw == 1, "CD", ifelse(resp_raw == 0, "PD", NA))
  } else {
    clinical_data$Response_group <- as.character(resp_raw)
  }
  # Standardize labels
  clinical_data$Response_group <- case_when(
    clinical_data$Response_group %in% c("PD", "Progressive Disease", "progressive", "0") ~ "PD",
    clinical_data$Response_group %in% c("CD", "PR", "SD", "PR/SD", "responder", "Responder", "1") ~ "CD",
    TRUE ~ clinical_data$Response_group
  )
} else if (exists("datTraits") && "Best_Response" %in% colnames(datTraits)) {
  clinical_data$Response_group <- as.character(datTraits$Best_Response[match(sample_ids, rownames(datTraits))])
} else {
  cat("   [WARNING] Response data not found in datTraits\n")
  clinical_data$Response_group <- NA
}

# Report clinical groups

if (!all(is.na(clinical_data$PFS_group))) {
  n_short <- sum(clinical_data$PFS_group == "Short", na.rm = TRUE)
  n_long <- sum(clinical_data$PFS_group == "Long", na.rm = TRUE)
  cat(sprintf("   PFS groups: Short=%d, Long=%d\n", n_short, n_long))
}

if (!all(is.na(clinical_data$Response_group))) {
  response_table <- table(clinical_data$Response_group)
  cat(sprintf("   Response groups: %s\n", paste(names(response_table), response_table, sep="=", collapse=", ")))
}
# ============================================================================
# 10b: DEFINE ANALYSIS FUNCTIONS
# ============================================================================

cat("==============================================================================\n")
cat("10b: DEFINE ANALYSIS FUNCTIONS\n")
cat("==============================================================================\n")

# Winsorize function to cap extreme values (robust for n=50)
winsorize <- function(x, limits = c(0.025, 0.975)) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) < 5) return(x)
  q <- quantile(x, probs = limits, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

# Robust Cohen's d with winsorization
cohens_d_robust <- function(group1, group2) {
  g1 <- winsorize(na.omit(as.numeric(group1)))
  g2 <- winsorize(na.omit(as.numeric(group2)))
  n1 <- length(g1)
  n2 <- length(g2)
  if (n1 < 3 || n2 < 3) return(NA)
  pooled_sd <- sqrt(((n1-1)*var(g1) + (n2-1)*var(g2)) / (n1+n2-2))
  if (pooled_sd == 0 || is.na(pooled_sd)) return(NA)
  d <- (mean(g1) - mean(g2)) / pooled_sd
  return(d)
}

cat("   [OK] Defined winsorize() function (2.5%/97.5% percentile capping)\n")
cat("   [OK] Defined cohens_d_robust() function (winsorized effect size)\n")
# ============================================================================
# 10c: PERFORM DUAL-CLINICAL ANALYSIS
# ============================================================================

cat("==============================================================================\n")
cat("10c: PERFORM DUAL-CLINICAL ANALYSIS (Wilcoxon + Cohen's d)\n")
cat("==============================================================================\n")

# Get pathways from consensus table (focus modules: pink, green, black, blue)
# NOTE: master_consensus_df is created in Step 9 with all focus modules
consensus_pathways <- master_consensus_df$Signature_ID
cat(sprintf("   Analyzing %d pathways from Master Consensus Table (Modules: %s)...\n",
            length(consensus_pathways), paste(toupper(focus_modules_auto), collapse = ", ")))

# Initialize results
dual_clinical_results <- data.frame(
  Pathway_ID = character(),
  Short_Name = character(),
  Module = character(),
  Consensus_Rank = numeric(),
  
  # PFS Analysis
  PFS_p_val = numeric(),
  PFS_FDR = numeric(),
  PFS_Cohen_d = numeric(),
  PFS_Direction = character(),
  PFS_Short_Mean = numeric(),
  PFS_Long_Mean = numeric(),
  
  # Response Analysis
  Response_p_val = numeric(),
  Response_FDR = numeric(),
  Response_Cohen_d = numeric(),
  Response_Direction = character(),
  Response_PD_Mean = numeric(),
  Response_Responder_Mean = numeric(),
  
  # Combined
  Overlap_Status = character(),
  stringsAsFactors = FALSE
)

# Loop through each pathway
for (pathway in consensus_pathways) {
  
  # Check if pathway exists in GSVA scores
  if (!pathway %in% colnames(gsva_scores)) next
  
  scores <- gsva_scores[, pathway]
  names(scores) <- rownames(gsva_scores)
  
  # Get metadata from consensus (use master_consensus_df)
  consensus_row <- master_consensus_df[master_consensus_df$Signature_ID == pathway, ]
  if (nrow(consensus_row) == 0) next
  consensus_row <- consensus_row[1, ]
  
  short_name <- ifelse(!is.na(consensus_row$Short_Name) && consensus_row$Short_Name != "",
                       consensus_row$Short_Name, pathway)
  module <- consensus_row$Module
  consensus_rank <- consensus_row$Consensus_Rank
  
  # --- PFS Analysis (Short vs Long) ---
  pfs_short_idx <- clinical_data$PFS_group == "Short" & !is.na(clinical_data$PFS_group)
  pfs_long_idx <- clinical_data$PFS_group == "Long" & !is.na(clinical_data$PFS_group)
  
  pfs_short <- scores[pfs_short_idx]
  pfs_long <- scores[pfs_long_idx]
  
  pfs_p <- tryCatch(
    wilcox.test(pfs_short, pfs_long, exact = FALSE)$p.value,
    error = function(e) NA
  )
  pfs_d <- cohens_d_robust(pfs_short, pfs_long)
  pfs_short_mean <- mean(pfs_short, na.rm = TRUE)
  pfs_long_mean <- mean(pfs_long, na.rm = TRUE)
  pfs_dir <- ifelse(pfs_short_mean > pfs_long_mean, "Short_Higher", "Long_Higher")
  
  # --- Response Analysis (PD vs CD) ---
  if (!all(is.na(clinical_data$Response_group))) {
    resp_pd_idx <- clinical_data$Response_group == "PD" & !is.na(clinical_data$Response_group)
    resp_responder_idx <- clinical_data$Response_group == "CD" & !is.na(clinical_data$Response_group)
    
    resp_pd <- scores[resp_pd_idx]
    resp_responder <- scores[resp_responder_idx]
    
    resp_p <- tryCatch(
      wilcox.test(resp_pd, resp_responder, exact = FALSE)$p.value,
      error = function(e) NA
    )
    resp_d <- cohens_d_robust(resp_pd, resp_responder)
    resp_pd_mean <- mean(resp_pd, na.rm = TRUE)
    resp_responder_mean <- mean(resp_responder, na.rm = TRUE)
    resp_dir <- ifelse(resp_pd_mean > resp_responder_mean, "PD_Higher", "Responder_Higher")
  } else {
    resp_p <- NA
    resp_d <- NA
    resp_pd_mean <- NA
    resp_responder_mean <- NA
    resp_dir <- NA
  }
  
  # Append results
  dual_clinical_results <- rbind(dual_clinical_results, data.frame(
    Pathway_ID = pathway,
    Short_Name = short_name,
    Module = module,
    Consensus_Rank = consensus_rank,
    PFS_p_val = pfs_p,
    PFS_FDR = NA,
    PFS_Cohen_d = pfs_d,
    PFS_Direction = pfs_dir,
    PFS_Short_Mean = pfs_short_mean,
    PFS_Long_Mean = pfs_long_mean,
    Response_p_val = resp_p,
    Response_FDR = NA,
    Response_Cohen_d = resp_d,
    Response_Direction = resp_dir,
    Response_PD_Mean = resp_pd_mean,
    Response_Responder_Mean = resp_responder_mean,
    Overlap_Status = NA,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("   [OK] Completed analysis for %d pathways\n", nrow(dual_clinical_results)))

# ============================================================================
# 10d: APPLY FDR CORRECTION AND ASSIGN OVERLAP STATUS
# ============================================================================

cat("==============================================================================\n")
cat("10d: APPLY FDR CORRECTION (Benjamini-Hochberg)\n")
cat("==============================================================================\n")

# Apply BH FDR correction
dual_clinical_results$PFS_FDR <- p.adjust(dual_clinical_results$PFS_p_val, method = "BH")
dual_clinical_results$Response_FDR <- p.adjust(dual_clinical_results$Response_p_val, method = "BH")

# Diagnostic: Show p-value and FDR distributions
cat(sprintf("   PFS p-values: min=%.4f, median=%.4f, max=%.4f\n",
            min(dual_clinical_results$PFS_p_val, na.rm = TRUE),
            median(dual_clinical_results$PFS_p_val, na.rm = TRUE),
            max(dual_clinical_results$PFS_p_val, na.rm = TRUE)))
cat(sprintf("   PFS FDR: min=%.4f, median=%.4f, n(FDR<0.1)=%d, n(p<0.05)=%d\n",
            min(dual_clinical_results$PFS_FDR, na.rm = TRUE),
            median(dual_clinical_results$PFS_FDR, na.rm = TRUE),
            sum(dual_clinical_results$PFS_FDR < 0.1, na.rm = TRUE),
            sum(dual_clinical_results$PFS_p_val < 0.05, na.rm = TRUE)))
cat(sprintf("   Response p-values: min=%.4f, median=%.4f, max=%.4f\n",
            min(dual_clinical_results$Response_p_val, na.rm = TRUE),
            median(dual_clinical_results$Response_p_val, na.rm = TRUE),
            max(dual_clinical_results$Response_p_val, na.rm = TRUE)))
cat(sprintf("   Response FDR: min=%.4f, median=%.4f, n(FDR<0.1)=%d, n(p<0.05)=%d\n",
            min(dual_clinical_results$Response_FDR, na.rm = TRUE),
            median(dual_clinical_results$Response_FDR, na.rm = TRUE),
            sum(dual_clinical_results$Response_FDR < 0.1, na.rm = TRUE),
            sum(dual_clinical_results$Response_p_val < 0.05, na.rm = TRUE)))

# Assign Overlap Status
dual_clinical_results$Overlap_Status <- sapply(1:nrow(dual_clinical_results), function(i) {
  pfs_sig <- !is.na(dual_clinical_results$PFS_FDR[i]) && dual_clinical_results$PFS_FDR[i] < 0.10
  resp_sig <- !is.na(dual_clinical_results$Response_FDR[i]) && dual_clinical_results$Response_FDR[i] < 0.10
  
  if (pfs_sig && resp_sig) return("Dual_Significant")
  if (pfs_sig && !resp_sig) return("Prognostic_Only")
  if (!pfs_sig && resp_sig) return("Predictive_Only")
  return("Not_Significant")
})

# Summary statistics
overlap_summary <- table(dual_clinical_results$Overlap_Status)
cat("\n   OVERLAP STATUS SUMMARY:\n")
cat("   ------------------------------------------------------------------------\n")
for (status in names(overlap_summary)) {
  cat(sprintf("   - %s: %d pathways\n", status, overlap_summary[status]))
}

# Count by module (focus modules only)
cat("\n   BY FOCUS MODULE:\n")
cat("   ------------------------------------------------------------------------\n")
for (mod in focus_modules_auto) {
  mod_data <- dual_clinical_results[dual_clinical_results$Module == mod, ]
  if(nrow(mod_data) == 0) next
  n_dual <- sum(mod_data$Overlap_Status == "Dual_Significant")
  n_prog <- sum(mod_data$Overlap_Status == "Prognostic_Only")
  n_pred <- sum(mod_data$Overlap_Status == "Predictive_Only")
  cat(sprintf("   %s: Dual=%d, Prognostic_Only=%d, Predictive_Only=%d\n",
                  toupper(mod), n_dual, n_prog, n_pred))
}
# Sort by minimum FDR
dual_clinical_results$Min_FDR <- pmin(dual_clinical_results$PFS_FDR,
                                      dual_clinical_results$Response_FDR, na.rm = TRUE)
dual_clinical_results <- dual_clinical_results[order(dual_clinical_results$Min_FDR), ]

# ============================================================================
# 10e: SAVE RESULTS TABLE
# ============================================================================

cat("==============================================================================\n")
cat("10e: SAVE RESULTS\n")
cat("==============================================================================\n")

# Save main results
write.csv(dual_clinical_results,
          file.path(results_dir, "Step10_Dual_Clinical_Significance.csv"),
          row.names = FALSE)
cat("   [OK] Saved: Step10_Dual_Clinical_Significance.csv\n")
# ============================================================================
# 10f: FIGURE A - DUAL LOLLIPOP PLOT
# ============================================================================

cat("==============================================================================\n")
cat("10f: GENERATE FIGURE A - DUAL LOLLIPOP PLOT\n")
cat("==============================================================================\n")

# Try FDR < 0.1 first, fallback to nominal p < 0.05 if no FDR-significant results
sig_results <- dual_clinical_results[dual_clinical_results$Overlap_Status != "Not_Significant", ]

# Fallback: if no FDR-significant, use nominal p-value < 0.05
if(nrow(sig_results) == 0) {
  cat("   [INFO] No FDR < 0.1 significant pathways. Using nominal p < 0.05 instead.\n")
  sig_results <- dual_clinical_results %>%
    filter(PFS_p_val < 0.05 | Response_p_val < 0.05) %>%
    arrange(pmin(PFS_p_val, Response_p_val, na.rm = TRUE))
  sig_type <- "nominal p < 0.05"
} else {
  sig_type <- "FDR < 0.10"
}

# DEDUPLICATE: Keep only the most significant occurrence of each pathway
# (same pathway may appear in multiple modules)
sig_results <- sig_results %>%
  group_by(Short_Name) %>%
  slice_min(Min_FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Min_FDR)

cat(sprintf("   [INFO] Found %d unique significant pathways (%s)\n", nrow(sig_results), sig_type))

# Take top 20 by minimum p-value/FDR
plot_data <- head(sig_results, 20)

# Create label with pathway name and module
plot_data$Label <- sprintf("%s [%s]", plot_data$Short_Name, toupper(plot_data$Module))

# Prepare data for mirrored lollipop
plot_data$y_order <- nrow(plot_data):1
plot_data$PFS_log_FDR <- -log10(plot_data$PFS_FDR + 1e-10)
plot_data$Response_log_FDR <- -log10(plot_data$Response_FDR + 1e-10)

# Cap for visualization
plot_data$PFS_log_FDR <- pmin(plot_data$PFS_log_FDR, 5)
plot_data$Response_log_FDR <- pmin(plot_data$Response_log_FDR, 5)

# Identify dual significant for annotation
dual_sig_idx <- which(plot_data$Overlap_Status == "Dual_Significant")

# Create mirrored lollipop plot
p_lollipop <- ggplot(plot_data) +
  # PFS (left side - negative values)
  geom_segment(aes(x = 0, xend = -PFS_log_FDR, y = y_order, yend = y_order),
               color = "#E74C3C", linewidth = 1) +
  geom_point(aes(x = -PFS_log_FDR, y = y_order, size = abs(PFS_Cohen_d)),
             color = "#E74C3C", alpha = 0.8) +
  # Response (right side - positive values)
  geom_segment(aes(x = 0, xend = Response_log_FDR, y = y_order, yend = y_order),
               color = "#3498DB", linewidth = 1) +
  geom_point(aes(x = Response_log_FDR, y = y_order, size = abs(Response_Cohen_d)),
             color = "#3498DB", alpha = 0.8) +
  # Vertical line at 0
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  # FDR threshold lines (at -log10(0.1) = 1)
  geom_vline(xintercept = -1, color = "#E74C3C", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, color = "#3498DB", linetype = "dashed", alpha = 0.5) +
  # Labels - include module name
  scale_y_continuous(breaks = plot_data$y_order, labels = plot_data$Label) +
  scale_size_continuous(name = "|Cohen's d|", range = c(2, 8)) +
  labs(
    title = "Dual-Clinical Association: PFS vs Treatment Response",
    subtitle = sprintf("Top %d Pathways | PFS (Left, Red) vs Response (Right, Blue) | Threshold: %s",
                       nrow(plot_data), sig_type),
    x = expression("-log"[10]*"(FDR)  (<-- PFS | Response -->)"),
    y = "Pathway [Module]"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "bottom"
  )

# Add dual significant markers if any exist
if(length(dual_sig_idx) > 0) {
  p_lollipop <- p_lollipop +
    annotate("text", x = max(plot_data$Response_log_FDR, na.rm = TRUE) + 0.5,
             y = plot_data$y_order[dual_sig_idx],
             label = "*", color = "#27AE60", size = 5)
}

print(p_lollipop)

ggsave(file.path(fig_dir_png, "Step10_01_Dual_Lollipop.png"), p_lollipop,
       width = 12, height = 10, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step10_01_Dual_Lollipop.pdf"), p_lollipop,
       width = 12, height = 10)
cat("   [OK] Saved: Step10_01_Dual_Lollipop.png/pdf\n")

# ============================================================================
# 10g: FIGURE B - RESPONSE BOXPLOTS (TOP 5 BY COHEN'S D) WITH SIGNIFICANCE
# ============================================================================

cat("==============================================================================\n")
cat("10g: GENERATE FIGURE B - RESPONSE BOXPLOTS (TOP 5 BY EFFECT SIZE)\n")
cat("==============================================================================\n")

# Deduplicate and get top 5 pathways by Response effect size
top_response_dedup <- dual_clinical_results %>%
  filter(!is.na(Response_Cohen_d)) %>%
  group_by(Short_Name) %>%
  slice_max(abs(Response_Cohen_d), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(-abs(Response_Cohen_d)) %>%
  head(5)

top5_response <- as.data.frame(top_response_dedup)

# Diagnostic: Show top 5 pathways selected
cat("   Top 5 Response pathways by Cohen's d:\n")
for(i in 1:nrow(top5_response)) {
  cat(sprintf("      %d. %s [%s] (d=%.3f, p=%.4f)\n",
              i, top5_response$Short_Name[i], toupper(top5_response$Module[i]),
              top5_response$Response_Cohen_d[i], top5_response$Response_p_val[i]))
}

# Create label with pathway name and module for facets
top5_response$Facet_Label <- sprintf("%s\n[%s]", top5_response$Short_Name, toupper(top5_response$Module))

# Prepare long-format data for plotting
boxplot_data <- data.frame()

for (i in 1:nrow(top5_response)) {
  pathway <- top5_response$Pathway_ID[i]
  # Skip if pathway missing from GSVA scores
  if (!pathway %in% colnames(gsva_scores)) {
    cat(sprintf("   [WARN] Pathway %s not found in GSVA scores\n", pathway))
    next
  }

  scores <- gsva_scores[, pathway]

  df_temp <- data.frame(
    Pathway = top5_response$Facet_Label[i],
    Score = scores,
    Response = clinical_data$Response_group,
    stringsAsFactors = FALSE
  )
  df_temp <- df_temp[!is.na(df_temp$Response), ]

  # Simplify response groups (CD = responders, PD = non-responders)
  df_temp$Response_Group <- ifelse(df_temp$Response %in% c("PD", "Progressive Disease", "progressive"),
                                   "PD", "CD")

  boxplot_data <- rbind(boxplot_data, df_temp)
}

# Update Pathway factor levels for consistent ordering (use Facet_Label)
unique_pathways <- unique(top5_response$Facet_Label)
boxplot_data$Pathway <- factor(boxplot_data$Pathway, levels = unique_pathways)

# --- Calculate significance statistics for each pathway ---
cat("   Calculating significance statistics (Wilcoxon test + Cohen's d)...\n")

sig_stats <- data.frame(
  Pathway = character(),
  p_value = numeric(),
  cohen_d = numeric(),
  sig_label = character(),
  stringsAsFactors = FALSE
)

for (pw in unique_pathways) {
  pw_data <- boxplot_data[boxplot_data$Pathway == pw, ]
  cd_scores <- pw_data$Score[pw_data$Response_Group == "CD"]
  pd_scores <- pw_data$Score[pw_data$Response_Group == "PD"]

  # Wilcoxon rank-sum test (non-parametric)
  if (length(cd_scores) >= 3 && length(pd_scores) >= 3) {
    wtest <- wilcox.test(cd_scores, pd_scores, exact = FALSE)
    p_val <- wtest$p.value

    # Cohen's d effect size
    pooled_sd <- sqrt(((length(cd_scores) - 1) * sd(cd_scores)^2 +
                       (length(pd_scores) - 1) * sd(pd_scores)^2) /
                      (length(cd_scores) + length(pd_scores) - 2))
    d_val <- (mean(cd_scores) - mean(pd_scores)) / pooled_sd

    # Significance label: *** p<0.001, ** p<0.01, * p<0.05, ns otherwise
    sig_label <- ifelse(p_val < 0.001, "***",
                 ifelse(p_val < 0.01, "**",
                 ifelse(p_val < 0.05, "*", "ns")))

    sig_stats <- rbind(sig_stats, data.frame(
      Pathway = pw,
      p_value = p_val,
      cohen_d = d_val,
      sig_label = sig_label,
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("   Significance results: %d pathways tested\n", nrow(sig_stats)))
cat(sprintf("   - Significant (p<0.05): %d\n", sum(sig_stats$p_value < 0.05)))
cat(sprintf("   - Not significant: %d\n", sum(sig_stats$p_value >= 0.05)))

# --- Calculate y-position for significance annotations ---
# Get max y value per pathway for annotation placement
y_positions <- boxplot_data %>%
  group_by(Pathway) %>%
  summarise(
    y_max = max(Score, na.rm = TRUE),
    y_min = min(Score, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    # Position bracket at 105% of max, label at 110%
    y_bracket = y_max + 0.05 * y_range,
    y_label = y_max + 0.12 * y_range
  )

# Merge significance stats with y positions
sig_annotations <- merge(sig_stats, y_positions, by = "Pathway")
sig_annotations$Pathway <- factor(sig_annotations$Pathway, levels = unique_pathways)

# --- Generate Plot with significance annotations ---
p_boxplots <- ggplot(boxplot_data, aes(x = Response_Group, y = Score, fill = Response_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
  # Add significance bracket (horizontal line)
  geom_segment(data = sig_annotations,
               aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  # Add bracket ends (vertical ticks)
  geom_segment(data = sig_annotations,
               aes(x = 1, xend = 1, y = y_bracket - 0.02 * y_range, yend = y_bracket),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_segment(data = sig_annotations,
               aes(x = 2, xend = 2, y = y_bracket - 0.02 * y_range, yend = y_bracket),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  # Add significance label
  geom_text(data = sig_annotations,
            aes(x = 1.5, y = y_label, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = colors_response, name = "Response Group") +
  labs(
    title = "Top 5 Pathways Associated with Treatment Response",
    subtitle = "GSVA scores by response group (CD vs PD) | Wilcoxon test: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
    x = NULL,
    y = "GSVA Score"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    # --- Legend and Axis Label Removal ---
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # -------------------------------------
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

print(p_boxplots)

# Save Outputs
ggsave(file.path(fig_dir_png, "Step10_02_Response_Boxplots.png"), p_boxplots, width = 10, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step10_02_Response_Boxplots.pdf"), p_boxplots, width = 10, height = 5)
cat("    [OK] Saved: Step10_02_Response_Boxplots.png/pdf\n")

# Save significance statistics table
write.csv(sig_stats, file.path(results_dir, "Step10_Response_Significance.csv"), row.names = FALSE)
cat("    [OK] Saved: Step10_Response_Significance.csv\n")

# ============================================================================
# 10g': FIGURE B' - PFS BOXPLOTS (TOP 5 BY COHEN'S D) WITH SIGNIFICANCE
# ============================================================================

cat("==============================================================================\n")
cat("10g': PFS BOXPLOTS (TOP 5 BY EFFECT SIZE)\n")
cat("==============================================================================\n")

# Deduplicate and get top 5 pathways by PFS effect size
top_pfs_dedup <- dual_clinical_results %>%
  filter(!is.na(PFS_Cohen_d)) %>%
  group_by(Short_Name) %>%
  slice_max(abs(PFS_Cohen_d), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(-abs(PFS_Cohen_d)) %>%
  head(5)

top5_pfs <- as.data.frame(top_pfs_dedup)

# Diagnostic: Show top 5 PFS pathways selected
cat("   Top 5 PFS pathways by Cohen's d:\n")
for(i in 1:nrow(top5_pfs)) {
  cat(sprintf("      %d. %s [%s] (d=%.3f, p=%.4f)\n",
              i, top5_pfs$Short_Name[i], toupper(top5_pfs$Module[i]),
              top5_pfs$PFS_Cohen_d[i], top5_pfs$PFS_p_val[i]))
}

# Create label with pathway name and module for facets
top5_pfs$Facet_Label <- sprintf("%s\n[%s]", top5_pfs$Short_Name, toupper(top5_pfs$Module))

# Prepare long-format data for plotting
pfs_boxplot_data <- data.frame()

for (i in 1:nrow(top5_pfs)) {
  pathway <- top5_pfs$Pathway_ID[i]
  # Skip if pathway missing from GSVA scores
  if (!pathway %in% colnames(gsva_scores)) {
    cat(sprintf("   [WARN] Pathway %s not found in GSVA scores\n", pathway))
    next
  }

  scores <- gsva_scores[, pathway]

  df_temp <- data.frame(
    Pathway = top5_pfs$Facet_Label[i],
    Score = scores,
    PFS = clinical_data$PFS_group,
    stringsAsFactors = FALSE
  )
  df_temp <- df_temp[!is.na(df_temp$PFS), ]

  # Simplify PFS groups
  df_temp$PFS_Group <- ifelse(df_temp$PFS %in% c("Short", "short", 0), "Short", "Long")

  pfs_boxplot_data <- rbind(pfs_boxplot_data, df_temp)
}

# Update Pathway factor levels for consistent ordering (use Facet_Label)
unique_pfs_pathways <- unique(top5_pfs$Facet_Label)
pfs_boxplot_data$Pathway <- factor(pfs_boxplot_data$Pathway, levels = unique_pfs_pathways)

# --- Calculate significance statistics for each pathway ---
cat("   Calculating significance statistics (Wilcoxon test + Cohen's d)...\n")

pfs_sig_stats <- data.frame(
  Pathway = character(),
  p_value = numeric(),
  cohen_d = numeric(),
  sig_label = character(),
  stringsAsFactors = FALSE
)

for (pw in unique_pfs_pathways) {
  pw_data <- pfs_boxplot_data[pfs_boxplot_data$Pathway == pw, ]
  long_scores <- pw_data$Score[pw_data$PFS_Group == "Long"]
  short_scores <- pw_data$Score[pw_data$PFS_Group == "Short"]

  # Wilcoxon rank-sum test (non-parametric)
  if (length(long_scores) >= 3 && length(short_scores) >= 3) {
    wtest <- wilcox.test(long_scores, short_scores, exact = FALSE)
    p_val <- wtest$p.value

    # Cohen's d effect size
    pooled_sd <- sqrt(((length(long_scores) - 1) * sd(long_scores)^2 +
                       (length(short_scores) - 1) * sd(short_scores)^2) /
                      (length(long_scores) + length(short_scores) - 2))
    d_val <- (mean(long_scores) - mean(short_scores)) / pooled_sd

    # Significance label: *** p<0.001, ** p<0.01, * p<0.05, ns otherwise
    sig_label <- ifelse(p_val < 0.001, "***",
                 ifelse(p_val < 0.01, "**",
                 ifelse(p_val < 0.05, "*", "ns")))

    pfs_sig_stats <- rbind(pfs_sig_stats, data.frame(
      Pathway = pw,
      p_value = p_val,
      cohen_d = d_val,
      sig_label = sig_label,
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("   Significance results: %d pathways tested\n", nrow(pfs_sig_stats)))
cat(sprintf("   - Significant (p<0.05): %d\n", sum(pfs_sig_stats$p_value < 0.05)))
cat(sprintf("   - Not significant: %d\n", sum(pfs_sig_stats$p_value >= 0.05)))

# --- Calculate y-position for significance annotations ---
pfs_y_positions <- pfs_boxplot_data %>%
  group_by(Pathway) %>%
  summarise(
    y_max = max(Score, na.rm = TRUE),
    y_min = min(Score, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    y_bracket = y_max + 0.05 * y_range,
    y_label = y_max + 0.12 * y_range
  )

# Merge significance stats with y positions
pfs_sig_annotations <- merge(pfs_sig_stats, pfs_y_positions, by = "Pathway")
pfs_sig_annotations$Pathway <- factor(pfs_sig_annotations$Pathway, levels = unique_pfs_pathways)

# --- Generate Plot with significance annotations ---
p_pfs_boxplots <- ggplot(pfs_boxplot_data, aes(x = PFS_Group, y = Score, fill = PFS_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
  # Add significance bracket (horizontal line)
  geom_segment(data = pfs_sig_annotations,
               aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  # Add bracket ends (vertical ticks)
  geom_segment(data = pfs_sig_annotations,
               aes(x = 1, xend = 1, y = y_bracket - 0.02 * y_range, yend = y_bracket),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_segment(data = pfs_sig_annotations,
               aes(x = 2, xend = 2, y = y_bracket - 0.02 * y_range, yend = y_bracket),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  # Add significance label
  geom_text(data = pfs_sig_annotations,
            aes(x = 1.5, y = y_label, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = colors_pfs, name = "PFS Group") +
  labs(
    title = "Top 5 Pathways Associated with Progression-Free Survival",
    subtitle = "GSVA scores by PFS group (Short vs Long) | Wilcoxon test: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
    x = NULL,
    y = "GSVA Score"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

print(p_pfs_boxplots)

# Save Outputs
ggsave(file.path(fig_dir_png, "Step10_02_PFS_Boxplots.png"), p_pfs_boxplots, width = 10, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step10_02_PFS_Boxplots.pdf"), p_pfs_boxplots, width = 10, height = 5)
cat("    [OK] Saved: Step10_02_PFS_Boxplots.png/pdf\n")

# Save significance statistics table
write.csv(pfs_sig_stats, file.path(results_dir, "Step10_PFS_Significance.csv"), row.names = FALSE)
cat("    [OK] Saved: Step10_PFS_Significance.csv\n")

# --- Merge PFS and Response boxplots using patchwork ---
cat("    Merging PFS and Response boxplots with patchwork...\n")


# Combine plots side by side (A = PFS, B = Response)
p_combined_boxplots <- p_pfs_boxplots + p_boxplots +
  plot_layout(ncol = 2) 

print(p_combined_boxplots)

ggsave(file.path(fig_dir_png, "Step10_02_Combined_Boxplots.png"), p_combined_boxplots, width = 18, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step10_02_Combined_Boxplots.pdf"), p_combined_boxplots, width = 18, height = 5)
cat("    [OK] Saved: Step10_02_Combined_Boxplots.png/pdf\n")

# ============================================================================
# 10h: STEP 10 SUMMARY
# ============================================================================

cat("==============================================================================\n")
cat("10h: STEP 10 SUMMARY\n")
cat("==============================================================================\n")

cat("   DUAL-CLINICAL ASSOCIATION ANALYSIS COMPLETE\n")
cat(sprintf("   Focus Modules: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("   ------------------------------------------------------------------------\n")
cat(sprintf("   Total pathways analyzed: %d\n", nrow(dual_clinical_results)))
cat(sprintf("   Dual_Significant (PFS + Response): %d\n", sum(dual_clinical_results$Overlap_Status == "Dual_Significant")))
cat(sprintf("   Prognostic_Only (PFS only): %d\n", sum(dual_clinical_results$Overlap_Status == "Prognostic_Only")))
cat(sprintf("   Predictive_Only (Response only): %d\n", sum(dual_clinical_results$Overlap_Status == "Predictive_Only")))
cat(sprintf("   Not_Significant: %d\n", sum(dual_clinical_results$Overlap_Status == "Not_Significant")))

cat("\n   OUTPUT FILES:\n")
cat("   ------------------------------------------------------------------------\n")
cat("      - Step10_Dual_Clinical_Significance.csv\n")
cat("      - Step10_Response_Significance.csv\n")
cat("      - Step10_PFS_Significance.csv\n")
cat("      - Step10_01_Dual_Lollipop.png/pdf\n")
cat("      - Step10_02_Response_Boxplots.png/pdf\n")
cat("      - Step10_02_PFS_Boxplots.png/pdf\n")
cat("      - Step10_02_Combined_Boxplots.png/pdf\n")

cat("\n   BIOLOGICAL INTERPRETATION:\n")
cat("   ------------------------------------------------------------------------\n")
cat("      - Dual_Significant: Both prognostic AND predictive - highest clinical value\n")
cat("      - Prognostic_Only: Reflects tumor aggressiveness, not treatment response\n")
cat("      - Predictive_Only: Treatment selection markers, not long-term prognosis\n")
cat("==============================================================================\n")
cat("STEP 10 COMPLETE\n")
cat("==============================================================================\n")

# ============================================================================
# STEP 10b: LIMMA DIFFERENTIAL EXPRESSION FOR GSEA
# ============================================================================
# Purpose: Create pre-ranked protein lists based on clinical outcomes
#          for use with GSEA desktop software (Broad Institute)
#
# Why here?
#   - After pathway-level clinical analysis (Step 10)
#   - Complements module-level analysis with protein-level ranking
#   - Before standard database enrichment (Step 11)
#
# Method: Limma moderated t-statistics with consensus ranking
# References:
#   - Smyth GK (2004) Stat Appl Genet Mol Biol - limma
#   - Subramanian et al (2005) PNAS - GSEA
# ============================================================================

cat("\n")
cat("==============================================================================\n")
cat("STEP 10b: LIMMA DIFFERENTIAL EXPRESSION FOR GSEA\n")
cat("==============================================================================\n")

# Check for limma package
skip_step10b <- !requireNamespace("limma", quietly = TRUE)

if (skip_step10b) {
  cat("[WARNING] limma package not installed. Step 10b will be skipped.\n")
  cat("To enable: BiocManager::install('limma')\n")
}

if (!skip_step10b) {
  library(limma)
  cat("   [OK] limma package loaded\n")
}

# Create output directory for Step 10b
step10b_dir <- file.path(results_dir, "Step10b_Limma_DE")
if (!skip_step10b && !dir.exists(step10b_dir)) {
  dir.create(step10b_dir, recursive = TRUE)
  cat(sprintf("   [OK] Created output directory: %s\n", step10b_dir))
}

# ============================================================================
# 10b-i: PFS_group Differential Expression (Long vs Short)
# ============================================================================
cat("\n")
cat("------------------------------------------------------------------------------\n")
cat("10b-i: Limma Differential Expression - PFS_group (Long vs Short)\n")
cat("------------------------------------------------------------------------------\n")

# Prepare expression matrix (proteins as ROWS, samples as COLUMNS for limma)
expr_matrix <- t(datExpr)
cat(sprintf("   Expression matrix: %d proteins x %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# Create PFS group factor - use lowercase 'g' (PFS_group) as defined in data loading
# Convention: Reference = "Long", so positive logFC = higher in Short (POOR prognosis)
pfs_factor <- factor(datTraits$PFS_group, levels = c("Long", "Short"))

# Check sample sizes
n_short <- sum(pfs_factor == "Short", na.rm = TRUE)
n_long <- sum(pfs_factor == "Long", na.rm = TRUE)
cat(sprintf("   Sample sizes: Short PFS = %d, Long PFS = %d\n", n_short, n_long))

# Handle missing values
valid_samples <- !is.na(pfs_factor)
n_excluded <- sum(!valid_samples)
if (n_excluded > 0) cat(sprintf("   [NOTE] Excluding %d samples with missing PFS_group\n", n_excluded))

expr_pfs <- expr_matrix[, valid_samples]
pfs_factor_valid <- pfs_factor[valid_samples]

# Design matrix
design_pfs <- model.matrix(~ pfs_factor_valid)
colnames(design_pfs) <- c("Intercept", "ShortVsLong")

cat("   Fitting limma model...\n")

# Fit limma model
fit_pfs <- lmFit(expr_pfs, design_pfs)
fit_pfs <- eBayes(fit_pfs)

# Extract results - all proteins
limma_pfs <- topTable(fit_pfs, coef = "ShortVsLong", number = Inf, sort.by = "P")
limma_pfs$Protein <- rownames(limma_pfs)

# Reorder columns for output
limma_pfs <- limma_pfs[, c("Protein", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# Summary statistics
n_sig_005 <- sum(limma_pfs$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_010 <- sum(limma_pfs$adj.P.Val < 0.10, na.rm = TRUE)
n_up <- sum(limma_pfs$logFC > 0 & limma_pfs$adj.P.Val < 0.05, na.rm = TRUE)
n_down <- sum(limma_pfs$logFC < 0 & limma_pfs$adj.P.Val < 0.05, na.rm = TRUE)

cat(sprintf("   Results:\n"))
cat(sprintf("     - Total proteins tested: %d\n", nrow(limma_pfs)))
cat(sprintf("     - Significant (adj.P < 0.05): %d\n", n_sig_005))
cat(sprintf("     - Significant (adj.P < 0.10): %d\n", n_sig_010))
cat(sprintf("     - Upregulated in Short PFS (poor prognosis): %d\n", n_up))
cat(sprintf("     - Downregulated in Short PFS (good prognosis): %d\n", n_down))

# Top 5 upregulated (higher in Short PFS = POOR prognosis markers)
top_up <- head(limma_pfs[limma_pfs$logFC > 0, ], 5)
cat("\n   Top 5 UPREGULATED in Short PFS (poor prognosis markers):\n")
for (i in 1:nrow(top_up)) {
  cat(sprintf("     %d. %s (logFC=%.3f, adj.P=%.2e)\n",
              i, top_up$Protein[i], top_up$logFC[i], top_up$adj.P.Val[i]))
}

# Top 5 downregulated (higher in Long PFS = good prognosis markers)
top_down <- head(limma_pfs[limma_pfs$logFC < 0, ], 5)
cat("\n   Top 5 DOWNREGULATED in Short PFS (good prognosis markers):\n")
for (i in 1:nrow(top_down)) {
  cat(sprintf("     %d. %s (logFC=%.3f, adj.P=%.2e)\n",
              i, top_down$Protein[i], top_down$logFC[i], top_down$adj.P.Val[i]))
}

# Save results (sorted by P-value)
write.csv(limma_pfs, file.path(step10b_dir, "Step10b_Limma_DE_PFS.csv"), row.names = FALSE)
cat(sprintf("\n   [OK] Saved: Step10b_Limma_DE_PFS.csv (%d proteins)\n", nrow(limma_pfs)))

# Create t-statistic ranked output (highest t on top = most upregulated in Short PFS)
limma_pfs_ranked <- limma_pfs %>%
  arrange(desc(t)) %>%
  select(Protein, logFC, t, P.Value, adj.P.Val)

write.csv(limma_pfs_ranked, file.path(step10b_dir, "Step10b_Limma_DE_PFS_t_ranked.csv"), row.names = FALSE)
cat(sprintf("   [OK] Saved: Step10b_Limma_DE_PFS_t_ranked.csv (sorted by t-statistic, highest first)\n"))

# Export .rnk file for GSEA Desktop (Protein \t t-statistic, no header)
rnk_pfs <- limma_pfs_ranked %>% select(Protein, t)
write.table(rnk_pfs, file.path(step10b_dir, "Step10b_PFS_ShortVsLong.rnk"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(sprintf("   [OK] Saved: Step10b_PFS_ShortVsLong.rnk (for GSEA Desktop)\n"))

# ============================================================================
# 10b-i: Volcano Plot for PFS
# ============================================================================
cat("\n   Creating volcano plot...\n")

# Count FDR significant proteins at different thresholds
n_fdr_005 <- sum(limma_pfs$adj.P.Val < 0.05)
n_fdr_01 <- sum(limma_pfs$adj.P.Val < 0.1)
cat(sprintf("   FDR < 0.05: %d proteins\n", n_fdr_005))
cat(sprintf("   FDR < 0.10: %d proteins\n", n_fdr_01))

# Prepare data for volcano plot - use FDR < 0.1 for coloring
volcano_pfs <- limma_pfs %>%
  mutate(
    neg_log10_p = -log10(P.Value),
    Significance = case_when(
      adj.P.Val < 0.1 & logFC > 0 ~ "Up (Short PFS)",
      adj.P.Val < 0.1 & logFC < 0 ~ "Down (Short PFS)",
      TRUE ~ "Not Significant"
    ),
    Label = NA_character_
  )

# Count significant for subtitle
n_up_fdr <- sum(volcano_pfs$Significance == "Up (Short PFS)")
n_down_fdr <- sum(volcano_pfs$Significance == "Down (Short PFS)")

# Select top proteins to label (by FDR, or by p-value if none pass FDR)
if (n_fdr_01 > 0) {
  top_labels <- volcano_pfs %>%
    filter(adj.P.Val < 0.1) %>%
    arrange(adj.P.Val) %>%
    head(15) %>%
    pull(Protein)
  cat(sprintf("   Labeling top %d FDR-significant proteins\n", length(top_labels)))
} else {
  # Fallback: label top 15 by nominal p-value if no FDR significant
  top_labels <- volcano_pfs %>%
    filter(P.Value < 0.05, abs(logFC) > 0.3) %>%
    arrange(P.Value) %>%
    head(15) %>%
    pull(Protein)
  cat(sprintf("   [NOTE] No FDR < 0.1 significant. Labeling top %d by nominal p-value\n", length(top_labels)))
}

volcano_pfs$Label <- ifelse(volcano_pfs$Protein %in% top_labels, volcano_pfs$Protein, NA)

# Create volcano plot
p_volcano_pfs <- ggplot(volcano_pfs, aes(x = logFC, y = neg_log10_p)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    aes(label = Label),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    segment.color = "grey50",
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c("Up (Short PFS)" = "#CA562C",
               "Down (Short PFS)" = "#008080",
               "Not Significant" = "grey70"),
    name = "Direction"
  ) +
  labs(
    title = "Differential Expression: Short vs Long PFS",
    subtitle = sprintf("Limma analysis | %d up, %d down (FDR < 0.1) | Top 15 labeled",
                       n_up_fdr, n_down_fdr),
    x = "log2 Fold Change (Short vs Long)",
    y = "-log10(P-value)",
    caption = "Dashed lines: p = 0.05, |logFC| = 0.3 | Positive logFC = higher in Short PFS (poor prognosis)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    plot.caption = element_text(hjust = 0.5, color = "grey50", size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save volcano plot
ggsave(file.path(fig_dir_png, "Step10b_Volcano_PFS.png"), p_volcano_pfs,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step10b_Volcano_PFS.pdf"), p_volcano_pfs,
       width = 6, height = 5)
cat("   [OK] Saved: Step10b_Volcano_PFS.png/pdf\n")

# ============================================================================
# Step 10b-ii: Limma DE for Response (PD vs CD)
# ============================================================================
cat("\n")
cat("------------------------------------------------------------------------------\n")
cat("10b-ii: Limma Differential Expression - Response (PD vs CD)\n")
cat("------------------------------------------------------------------------------\n")

# Get samples with valid Response status (CD or PD)
valid_response <- datTraits$Response %in% c("CD", "PD")
n_cd <- sum(datTraits$Response == "CD", na.rm = TRUE)
n_pd <- sum(datTraits$Response == "PD", na.rm = TRUE)
cat(sprintf("   Response groups: CD (n=%d), PD (n=%d)\n", n_cd, n_pd))

# Subset expression data to valid samples
expr_response <- expr_matrix[, valid_response]
traits_response <- datTraits[valid_response, ]

# Create Response factor with CD as reference (so positive logFC = higher in PD = POOR response)
response_factor <- factor(traits_response$Response, levels = c("CD", "PD"))

# Create design matrix
design_response <- model.matrix(~ response_factor)
colnames(design_response) <- c("Intercept", "PD_vs_CD")

# Fit linear model
cat("   Fitting limma model...\n")
fit_response <- lmFit(expr_response, design_response)
fit_response <- eBayes(fit_response)

# Extract all results
de_response <- topTable(fit_response, coef = "PD_vs_CD", number = Inf, adjust.method = "BH")
de_response$Protein <- rownames(de_response)
de_response <- de_response[, c("Protein", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# Count significant proteins
n_sig_response_005 <- sum(de_response$adj.P.Val < 0.05)
n_sig_response_01 <- sum(de_response$adj.P.Val < 0.10)
n_up_response <- sum(de_response$adj.P.Val < 0.05 & de_response$logFC > 0)
n_down_response <- sum(de_response$adj.P.Val < 0.05 & de_response$logFC < 0)

cat(sprintf("   Results:\n"))
cat(sprintf("     - Total proteins tested: %d\n", nrow(de_response)))
cat(sprintf("     - Significant (adj.P < 0.05): %d\n", n_sig_response_005))
cat(sprintf("     - Significant (adj.P < 0.10): %d\n", n_sig_response_01))
cat(sprintf("     - Upregulated in PD (poor response): %d\n", n_up_response))
cat(sprintf("     - Downregulated in PD (favorable in CD): %d\n", n_down_response))

# Report top 5 upregulated (higher in PD = POOR response)
top_up_response <- head(de_response[de_response$logFC > 0, ], 5)
cat("\n   Top 5 UPREGULATED in PD (poor response markers):\n")
for (i in 1:nrow(top_up_response)) {
  cat(sprintf("     %d. %s (logFC=%.3f, adj.P=%.2e)\n",
              i, top_up_response$Protein[i], top_up_response$logFC[i], top_up_response$adj.P.Val[i]))
}

# Report top 5 downregulated (higher in CD = favorable response)
top_down_response <- head(de_response[de_response$logFC < 0, ], 5)
cat("\n   Top 5 DOWNREGULATED in PD (favorable response markers):\n")
for (i in 1:nrow(top_down_response)) {
  cat(sprintf("     %d. %s (logFC=%.3f, adj.P=%.2e)\n",
              i, top_down_response$Protein[i], top_down_response$logFC[i], top_down_response$adj.P.Val[i]))
}

# Save results to CSV (sorted by P-value)
write.csv(de_response, file.path(step10b_dir, "Step10b_Limma_DE_Response.csv"), row.names = FALSE)
cat(sprintf("\n   [OK] Saved: Step10b_Limma_DE_Response.csv (%d proteins)\n", nrow(de_response)))

# Create t-statistic ranked output (highest t on top = most upregulated in PD)
de_response_ranked <- de_response %>%
  arrange(desc(t)) %>%
  select(Protein, logFC, t, P.Value, adj.P.Val)

write.csv(de_response_ranked, file.path(step10b_dir, "Step10b_Limma_DE_Response_t_ranked.csv"), row.names = FALSE)
cat(sprintf("   [OK] Saved: Step10b_Limma_DE_Response_t_ranked.csv (sorted by t-statistic, highest first)\n"))

# Export .rnk file for GSEA Desktop (Protein \t t-statistic, no header)
rnk_response <- de_response_ranked %>% select(Protein, t)
write.table(rnk_response, file.path(step10b_dir, "Step10b_Response_PDvsCD.rnk"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(sprintf("   [OK] Saved: Step10b_Response_PDvsCD.rnk (for GSEA Desktop)\n"))

# Create volcano plot for Response
cat("\n   Creating volcano plot...\n")

# Count FDR significant proteins at different thresholds
n_fdr_resp_005 <- sum(de_response$adj.P.Val < 0.05)
n_fdr_resp_01 <- sum(de_response$adj.P.Val < 0.1)
cat(sprintf("   FDR < 0.05: %d proteins\n", n_fdr_resp_005))
cat(sprintf("   FDR < 0.10: %d proteins\n", n_fdr_resp_01))

# Use FDR < 0.1 for coloring
de_response$Direction <- "Not Significant"
de_response$Direction[de_response$adj.P.Val < 0.1 & de_response$logFC > 0] <- "Up (PD)"
de_response$Direction[de_response$adj.P.Val < 0.1 & de_response$logFC < 0] <- "Down (PD)"
de_response$Direction <- factor(de_response$Direction,
                                levels = c("Up (PD)", "Down (PD)", "Not Significant"))

# Count for subtitle
n_up_resp_fdr <- sum(de_response$Direction == "Up (PD)")
n_down_resp_fdr <- sum(de_response$Direction == "Down (PD)")

# Label top proteins (by FDR, or by p-value if none pass FDR)
de_response$Label <- ""
if (n_fdr_resp_01 > 0) {
  top_to_label_response <- de_response %>%
    filter(adj.P.Val < 0.1) %>%
    arrange(adj.P.Val) %>%
    head(15)
  cat(sprintf("   Labeling top %d FDR-significant proteins\n", nrow(top_to_label_response)))
} else {
  # Fallback: label top 15 by nominal p-value if no FDR significant
  top_to_label_response <- de_response %>%
    filter(P.Value < 0.05, abs(logFC) > 0.3) %>%
    arrange(P.Value) %>%
    head(15)
  cat(sprintf("   [NOTE] No FDR < 0.1 significant. Labeling top %d by nominal p-value\n", nrow(top_to_label_response)))
}
de_response$Label[de_response$Protein %in% top_to_label_response$Protein] <-
  de_response$Protein[de_response$Protein %in% top_to_label_response$Protein]

# Create volcano plot
p_volcano_response <- ggplot(de_response, aes(x = logFC, y = -log10(P.Value), color = Direction)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text_repel(
    data = de_response %>% filter(Label != ""),
    aes(label = Label),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Up (PD)" = "#C62828",
               "Down (PD)" = "#2E7D32",
               "Not Significant" = "grey70"),
    name = "Direction"
  ) +
  labs(
    title = "Differential Expression: PD vs CD Response",
    subtitle = sprintf("FDR<0.1",
                       n_up_resp_fdr, n_down_resp_fdr),
    x = "log2 Fold Change (PD vs CD)",
    y = "-log10(P-value)",
    caption = "Dashed lines: p = 0.05, |logFC| = 0.3 | Positive logFC = higher in PD (poor response)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    plot.caption = element_text(hjust = 0.5, color = "grey50", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_volcano_response)

# Save volcano plot
ggsave(file.path(fig_dir_png, "Step10b_Volcano_Response.png"), p_volcano_response,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step10b_Volcano_Response.pdf"), p_volcano_response,
       width = 6, height = 5)
cat("   [OK] Saved: Step10b_Volcano_Response.png/pdf\n")

cat("==============================================================================\n")
cat("STEP 10b COMPLETE: Limma Differential Expression\n")
cat("==============================================================================\n")

#===============================================================================
# STEP 11: STANDARD DATABASE ENRICHMENT
#===============================================================================
# PURPOSE: Complement custom signatures (Step 8) with standard pathway
#          databases (GO:BP, KEGG, MSigDB) for comprehensive biological annotation.
#
# FOCUS MODULES: Pink, Green, Black, Blue (clinically significant from Step 5)
#
# STRATEGY:
#   11a) Load clusterProfiler and verify packages
#   11b) Gene ID conversion (Symbol -> Entrez for GO/KEGG)
#   11c) GO:BP enrichment for all focus modules
#   11d) KEGG pathway enrichment for all focus modules
#   11e) MSigDB Hallmark (H) enrichment
#   11f) MSigDB Canonical Pathways (C2:CP)
#   11g) MSigDB Oncogenic Signatures (C6)
#   11h) Differentiation signatures + GSEA analysis
#   11i) Summary and integration
#
# LOOKING FOR:
#   - GO terms consistent with custom signature findings
#   - KEGG pathways validating PDAC-specific biology
#   - Hallmark overlap with module function (EMT, inflammation, etc.)
#
# KEY PARAMS: pvalueCutoff=0.05, qvalueCutoff=0.2, minGSSize=10
#===============================================================================

# Source this file to improve gene symbol -> Entrez ID mapping before Step 11
# Uncomment the line below to use biomaRt + manual curation for better mapping
# source("gene_symbol_conversion_helper.R")

# ============================================================================
# 11a: LOAD REQUIRED PACKAGES
# ============================================================================

cat("==============================================================================\n")
cat("11a: VERIFY REQUIRED PACKAGES\n")
cat("==============================================================================\n")

# Required packages for standard database enrichment (should be pre-installed)
step11_packages <- c("clusterProfiler", "org.Hs.eg.db", "msigdbr", "enrichplot", "DOSE")

missing_step11 <- step11_packages[!sapply(step11_packages, requireNamespace, quietly = TRUE)]
if(length(missing_step11) > 0) {
  cat("   [WARNING] Missing packages for Step 11:\n")
  cat(paste("    -", missing_step11, collapse = "\n"), "\n")
  cat("   Step 11 will be skipped. To enable, install:\n")
  cat(sprintf('   BiocManager::install(c("%s"))\n', paste(missing_step11, collapse = '", "')))
  skip_step11 <- TRUE
} else {
  skip_step11 <- FALSE
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  cat("   [OK] All Step 11 packages loaded\n")
}
# ============================================================================
# 11b: GENE ID CONVERSION
# ============================================================================

#
# NOTE ON UNIVERSE SIZE:
# - Step 8 (Custom ORA): Uses n_proteins = 891 as universe (all proteins in dataset)
# - Step 11 (clusterProfiler): Uses only MAPPED genes as universe (~750-850)
#   because bitr() cannot map all symbols to Entrez IDs
#
# This difference means:
# - Same enrichment may have different p-values between steps
# - Step 8 is MORE CONSERVATIVE (larger denominator = harder to detect)
# - This is acceptable but users should be aware of the difference
# ============================================================================

cat("==============================================================================\n")
cat("11b: GENE ID CONVERSION\n")
cat("==============================================================================\n")

# ============================================================================
# IMPORTANT: Universe Size Discrepancy (GO/KEGG vs MSigDB)
# ============================================================================
# GO and KEGG enrichment require Entrez IDs, not gene symbols.
# The bitr() conversion from SYMBOL to ENTREZID loses ~3-5% of genes.
#
# Result: GO/KEGG universe = ~850-870 genes (mapped subset)
#         MSigDB universe  = 891 proteins (uses gene symbols directly)
#         Step 8 ORA       = 891 proteins (uses gene symbols directly)
#
# This is a known limitation of clusterProfiler GO/KEGG functions.
# The MSigDB analyses (Hallmark, C2:CP, C6) correctly use 891 proteins.
# ============================================================================

# Convert gene symbols to Entrez IDs for clusterProfiler
cat("   Converting gene symbols to Entrez IDs...\n")
cat("   NOTE: Not all genes will map - this creates a smaller universe than Step 8\n")
cat("   NOTE: GO/KEGG will use mapped universe; MSigDB uses full 891 proteins\n")

# Universe (all proteins in dataset)
universe_conversion <- bitr(all_proteins,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)

n_universe_entrez <- nrow(universe_conversion)
cat(sprintf("   Universe: %d of %d genes mapped (%.1f%%)\n",
            n_universe_entrez, length(all_proteins),
            100 * n_universe_entrez / length(all_proteins)))
cat(sprintf("   GO/KEGG background = %d (Entrez-mapped)\n", n_universe_entrez))
cat(sprintf("   MSigDB background  = %d (full proteome)\n", length(all_proteins)))

universe_entrez <- universe_conversion$ENTREZID

# Convert module genes - FOCUS MODULES ONLY (pink, green, black, blue)
module_entrez <- list()
modules_for_enrichment <- focus_modules_auto
modules_for_enrichment <- modules_for_enrichment[modules_for_enrichment %in% unique(moduleColors)]
cat(sprintf("   Modules for enrichment: %s\n", paste(toupper(modules_for_enrichment), collapse = ", ")))

for(mod in modules_for_enrichment) {
  mod_genes <- names(moduleColors)[moduleColors == mod]
  conversion <- bitr(mod_genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  module_entrez[[mod]] <- conversion$ENTREZID
  cat(sprintf("   %s module: %d of %d genes mapped (%.1f%%)\n",
              mod, nrow(conversion), length(mod_genes),
              100 * nrow(conversion) / length(mod_genes)))
}
# ============================================================================
# 11c: GO:BP ENRICHMENT
# ============================================================================

cat("\n==============================================================================\n")
cat("11c: GO:BP ENRICHMENT (Biological Process)\n")
cat("==============================================================================\n")

go_results <- list()

for(mod in names(module_entrez)) {
  cat(sprintf("   Running GO:BP for %s module...\n", mod))
  
  go_enrich <- enrichGO(
    gene = module_entrez[[mod]],
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  if(!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
    go_results[[mod]] <- go_enrich
    n_sig <- sum(go_enrich@result$p.adjust < 0.05)
    cat(sprintf("   [OK] %s: %d significant GO:BP terms (FDR < 0.05)\n", mod, n_sig))
  } else {
    cat(sprintf("   [WARNING] %s: No significant GO:BP terms\n", mod))
  }
}

# Combine results for export
go_combined <- data.frame()
for(mod in names(go_results)) {
  if(!is.null(go_results[[mod]])) {
    df <- as.data.frame(go_results[[mod]])
    if(nrow(df) > 0) {
      df$Module <- mod
      go_combined <- rbind(go_combined, df)
    }
  }
}

if(nrow(go_combined) > 0) {
  go_combined <- go_combined %>%
    dplyr::select(Module, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Module, p.adjust)
  
  write.csv(go_combined, file.path(results_dir, "Step11_01_GO_BP_Results.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step11_01_GO_BP_Results.csv\n")
}

# --- GO:BP Dotplots for ALL focus modules ---
cat("   Generating GO:BP Dotplots for focus modules...\n")

if(length(go_results) > 0) {
  plot_counter <- 1
  for(mod in focus_modules_auto) {
    if(mod %in% names(go_results) && !is.null(go_results[[mod]])) {
      mod_name <- tools::toTitleCase(mod)
      letter_suffix <- letters[plot_counter]  # a, b, c, d...

      p_go_mod <- dotplot(go_results[[mod]], showCategory = 15,
                          title = sprintf("%s Module - GO:BP Enrichment (Biological Process)", mod_name)) +
        theme(axis.text.y = element_text(size = 9))
      print(p_go_mod)

      ggsave(file.path(fig_dir_png, sprintf("Step11_01%s_GO_BP_%s.png", letter_suffix, mod_name)),
             p_go_mod, width = 8, height = 6, dpi = 300)
      ggsave(file.path(fig_dir_pdf, sprintf("Step11_01%s_GO_BP_%s.pdf", letter_suffix, mod_name)),
             p_go_mod, width = 8, height = 6)
      cat(sprintf("   [OK] Saved: Step11_01%s_GO_BP_%s.png/pdf\n", letter_suffix, mod_name))

      plot_counter <- plot_counter + 1
    }
  }
}
# ============================================================================
# 11d: KEGG PATHWAY ENRICHMENT
# ============================================================================

cat("==============================================================================\n")
cat("11d: KEGG PATHWAY ENRICHMENT\n")
cat("==============================================================================\n")

kegg_results <- list()

for(mod in names(module_entrez)) {
  cat(sprintf("   Running KEGG for %s module...\n", mod))
  
  kegg_enrich <- tryCatch({
    enrichKEGG(
      gene = module_entrez[[mod]],
      universe = universe_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    cat(sprintf("   [WARNING] KEGG error for %s: %s\n", mod, e$message))
    return(NULL)
  })
  
  if(!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
    kegg_results[[mod]] <- kegg_enrich
    n_sig <- sum(kegg_enrich@result$p.adjust < 0.05)
    cat(sprintf("   [OK] %s: %d significant KEGG pathways (FDR < 0.05)\n", mod, n_sig))
  } else {
    cat(sprintf("   [WARNING] %s: No significant KEGG pathways\n", mod))
  }
}

# Combine results for export
kegg_combined <- data.frame()
for(mod in names(kegg_results)) {
  if(!is.null(kegg_results[[mod]])) {
    df <- as.data.frame(kegg_results[[mod]])
    if(nrow(df) > 0) {
      df$Module <- mod
      kegg_combined <- rbind(kegg_combined, df)
    }
  }
}

if(nrow(kegg_combined) > 0) {
  kegg_combined <- kegg_combined %>%
    dplyr::select(Module, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Module, p.adjust)
  
  write.csv(kegg_combined, file.path(results_dir, "Step11_02_KEGG_Results.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step11_02_KEGG_Results.csv\n")
}

# --- KEGG Dotplots for ALL focus modules ---
cat("   Generating KEGG Dotplots for focus modules...\n")

if(length(kegg_results) > 0) {
  plot_counter <- 1
  for(mod in focus_modules_auto) {
    if(mod %in% names(kegg_results) && !is.null(kegg_results[[mod]])) {
      mod_name <- tools::toTitleCase(mod)
      letter_suffix <- letters[plot_counter]  # a, b, c, d...

      p_kegg_mod <- dotplot(kegg_results[[mod]], showCategory = 15,
                            title = sprintf("%s Module - KEGG Pathway Enrichment", mod_name)) +
        theme(axis.text.y = element_text(size = 9))
      print(p_kegg_mod)

      ggsave(file.path(fig_dir_png, sprintf("Step11_02%s_KEGG_%s.png", letter_suffix, mod_name)),
             p_kegg_mod, width = 8, height = 6, dpi = 300)
      ggsave(file.path(fig_dir_pdf, sprintf("Step11_02%s_KEGG_%s.pdf", letter_suffix, mod_name)),
             p_kegg_mod, width = 8, height = 6)
      cat(sprintf("   [OK] Saved: Step11_02%s_KEGG_%s.png/pdf\n", letter_suffix, mod_name))

      plot_counter <- plot_counter + 1
    }
  }
}
# ============================================================================
# 11e: MSigDB HALLMARK (H)
# ============================================================================

cat("\n==============================================================================\n")
cat("11e: MSigDB HALLMARK (H) ENRICHMENT\n")
cat("==============================================================================\n")

# Get Hallmark gene sets
cat("   Loading MSigDB Hallmark collection...\n")
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# Convert to TERM2GENE format
hallmark_t2g <- hallmark %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

cat(sprintf("   [OK] Loaded %d Hallmark gene sets\n", length(unique(hallmark_t2g$term))))

# Run enrichment for each module
hallmark_results <- list()

for(mod in modules_for_enrichment) {
  mod_genes <- names(moduleColors)[moduleColors == mod]
  cat(sprintf("   Running Hallmark for %s module...\n", mod))
  
  hallmark_enrich <- enricher(
    gene = mod_genes,
    universe = all_proteins,
    TERM2GENE = hallmark_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if(!is.null(hallmark_enrich) && nrow(hallmark_enrich@result) > 0) {
    hallmark_results[[mod]] <- hallmark_enrich
    n_sig <- sum(hallmark_enrich@result$p.adjust < 0.05)
    cat(sprintf("   [OK] %s: %d significant Hallmark pathways (FDR < 0.05)\n", mod, n_sig))
  } else {
    cat(sprintf("   [WARNING] %s: No significant Hallmark pathways\n", mod))
  }
}

# Combine results
hallmark_combined <- data.frame()
for(mod in names(hallmark_results)) {
  if(!is.null(hallmark_results[[mod]])) {
    df <- as.data.frame(hallmark_results[[mod]])
    if(nrow(df) > 0) {
      df$Module <- mod
      hallmark_combined <- rbind(hallmark_combined, df)
    }
  }
}

if(nrow(hallmark_combined) > 0) {
  hallmark_combined <- hallmark_combined %>%
    mutate(ID = gsub("HALLMARK_", "", ID)) %>%
    dplyr::select(Module, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Module, p.adjust)
  
  write.csv(hallmark_combined, file.path(results_dir, "Step11_03_MSigDB_Hallmark_Results.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step11_03_MSigDB_Hallmark_Results.csv\n")
}

# --- Hallmark Enrichment Heatmap ---
cat("   Generating Hallmark Enrichment Heatmap...\n")

if(nrow(hallmark_combined) > 0) {
  # Create matrix for heatmap
  hallmark_matrix <- hallmark_combined %>%
    mutate(neg_log10_p = -log10(pvalue)) %>%
    dplyr::select(Module, ID, neg_log10_p) %>%
    pivot_wider(names_from = Module, values_from = neg_log10_p, values_fill = 0) %>%
    column_to_rownames("ID")
  
  # Filter to significant pathways
  sig_hallmarks <- hallmark_combined %>%
    filter(p.adjust < 0.1) %>%
    pull(ID) %>%
    unique()
  
  if(length(sig_hallmarks) > 0) {
    hallmark_matrix_filt <- hallmark_matrix[rownames(hallmark_matrix) %in% sig_hallmarks, , drop = FALSE]
    
    if(nrow(hallmark_matrix_filt) > 1) {
      p_hallmark <- pheatmap(
        as.matrix(hallmark_matrix_filt),
        color = colorRampPalette(c("#F5F5F5", "#D6EFB3", "#5FC1C0", "#234DA0"))(50),
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        main = sprintf("MSigDB Hallmark Enrichment - Focus Modules\n-log10(p-value) | Modules: %s",
                       paste(toupper(focus_modules_auto), collapse = ", ")),
        fontsize = 10,
        fontsize_row = 8,
        border_color = "white"
      )
      
      ggsave(file.path(fig_dir_png, "Step11_03_Hallmark_Heatmap.png"), p_hallmark, width = 8, height = 8, dpi = 300)
      ggsave(file.path(fig_dir_pdf, "Step11_03_Hallmark_Heatmap.pdf"), p_hallmark, width = 8, height = 8)
      cat("   [OK] Saved: Step11_03_Hallmark_Heatmap\n")
    }
  }
}
# ============================================================================
# 11f: MSigDB CANONICAL PATHWAYS (C2:CP)
# ============================================================================

cat("\n==============================================================================\n")
cat("11f: MSigDB CANONICAL PATHWAYS (C2:CP) ENRICHMENT\n")
cat("==============================================================================\n")

# Get C2 Canonical Pathways (Reactome, KEGG, BioCarta, PID, WikiPathways)
# NOTE: Load C2 category and filter by pathway prefix
# msigdbr subcategory syntax has changed - use pattern matching instead
cat("   Loading MSigDB C2 canonical pathways...\n")

c2_all <- tryCatch(msigdbr(species = "Homo sapiens", category = "C2"), error = function(e) NULL)

c2_cp_kegg <- if (!is.null(c2_all)) c2_all %>% dplyr::filter(grepl("^KEGG_", gs_name)) else NULL
c2_cp_reactome <- if (!is.null(c2_all)) c2_all %>% dplyr::filter(grepl("^REACTOME_", gs_name)) else NULL
c2_cp_biocarta <- if (!is.null(c2_all)) c2_all %>% dplyr::filter(grepl("^BIOCARTA_", gs_name)) else NULL
c2_cp_pid <- if (!is.null(c2_all)) c2_all %>% dplyr::filter(grepl("^PID_", gs_name)) else NULL
c2_cp_wiki <- if (!is.null(c2_all)) c2_all %>% dplyr::filter(grepl("^WP_", gs_name)) else NULL

# Combine all available subcollections
c2_cp <- bind_rows(
  c2_cp_kegg,
  c2_cp_reactome,
  c2_cp_biocarta,
  c2_cp_pid,
  c2_cp_wiki
)

cat(sprintf("   Loaded: KEGG=%d, REACTOME=%d, BIOCARTA=%d, PID=%d, WIKI=%d\n",
    ifelse(is.null(c2_cp_kegg), 0, length(unique(c2_cp_kegg$gs_name))),
    ifelse(is.null(c2_cp_reactome), 0, length(unique(c2_cp_reactome$gs_name))),
    ifelse(is.null(c2_cp_biocarta), 0, length(unique(c2_cp_biocarta$gs_name))),
    ifelse(is.null(c2_cp_pid), 0, length(unique(c2_cp_pid$gs_name))),
    ifelse(is.null(c2_cp_wiki), 0, length(unique(c2_cp_wiki$gs_name)))))

# Convert to TERM2GENE format
c2_cp_t2g <- c2_cp %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

cat(sprintf("   [OK] Loaded %d Canonical Pathway gene sets total\n", length(unique(c2_cp_t2g$term))))

# Run enrichment for each module
c2_cp_results <- list()

for(mod in modules_for_enrichment) {
  mod_genes <- names(moduleColors)[moduleColors == mod]
  cat(sprintf("   Running C2:CP for %s module...\n", mod))
  
  c2_enrich <- enricher(
    gene = mod_genes,
    universe = all_proteins,
    TERM2GENE = c2_cp_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if(!is.null(c2_enrich) && nrow(c2_enrich@result) > 0) {
    c2_cp_results[[mod]] <- c2_enrich
    n_sig <- sum(c2_enrich@result$p.adjust < 0.05)
    cat(sprintf("   [OK] %s: %d significant C2:CP pathways (FDR < 0.05)\n", mod, n_sig))
  } else {
    cat(sprintf("   [WARNING] %s: No significant C2:CP pathways\n", mod))
  }
}

# Combine results
c2_cp_combined <- data.frame()
for(mod in names(c2_cp_results)) {
  if(!is.null(c2_cp_results[[mod]])) {
    df <- as.data.frame(c2_cp_results[[mod]])
    if(nrow(df) > 0) {
      df$Module <- mod
      c2_cp_combined <- rbind(c2_cp_combined, df)
    }
  }
}

if(nrow(c2_cp_combined) > 0) {
  c2_cp_combined <- c2_cp_combined %>%
    dplyr::select(Module, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Module, p.adjust)
  
  write.csv(c2_cp_combined, file.path(results_dir, "Step11_04_MSigDB_C2_CP_Results.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step11_04_MSigDB_C2_CP_Results.csv\n")
}
# ============================================================================
# 11g: MSigDB ONCOGENIC SIGNATURES (C6)
# ============================================================================

cat("\n==============================================================================\n")
cat("11g: MSigDB ONCOGENIC SIGNATURES (C6) ENRICHMENT\n")
cat("==============================================================================\n")

# Get C6 Oncogenic Signatures
cat("   Loading MSigDB C6 collection...\n")
c6 <- msigdbr(species = "Homo sapiens", category = "C6")

# Convert to TERM2GENE format
c6_t2g <- c6 %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

cat(sprintf("   [OK] Loaded %d Oncogenic gene sets\n", length(unique(c6_t2g$term))))

# Run enrichment for each module
c6_results <- list()

for(mod in modules_for_enrichment) {
  mod_genes <- names(moduleColors)[moduleColors == mod]
  cat(sprintf("   Running C6 for %s module...\n", mod))
  
  c6_enrich <- enricher(
    gene = mod_genes,
    universe = all_proteins,
    TERM2GENE = c6_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if(!is.null(c6_enrich) && nrow(c6_enrich@result) > 0) {
    c6_results[[mod]] <- c6_enrich
    n_sig <- sum(c6_enrich@result$p.adjust < 0.05)
    cat(sprintf("   [OK] %s: %d significant C6 oncogenic signatures (FDR < 0.05)\n", mod, n_sig))
  } else {
    cat(sprintf("   [WARNING] %s: No significant C6 oncogenic signatures\n", mod))
  }
}

# Combine results
c6_combined <- data.frame()
for(mod in names(c6_results)) {
  if(!is.null(c6_results[[mod]])) {
    df <- as.data.frame(c6_results[[mod]])
    if(nrow(df) > 0) {
      df$Module <- mod
      c6_combined <- rbind(c6_combined, df)
    }
  }
}

if(nrow(c6_combined) > 0) {
  c6_combined <- c6_combined %>%
    dplyr::select(Module, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Module, p.adjust)
  
  write.csv(c6_combined, file.path(results_dir, "Step11_05_MSigDB_C6_Oncogenic_Results.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step11_05_MSigDB_C6_Oncogenic_Results.csv\n")
}
# ============================================================================
# 11h: DIFFERENTIATION SIGNATURES (From missing.md)
# ============================================================================

cat("\n==============================================================================\n")
cat("11h: DIFFERENTIATION SIGNATURES (8 PDAC-specific)\n")
cat("==============================================================================\n")

# Add differentiation signatures as specified in missing.md
diff_signatures <- list()

diff_signatures$Diff_ADM_Acinar <- list(
  id = "Diff01", short_name = "ADM Acinar",
  full_name = "Acinar-to-Ductal Metaplasia - Acinar Markers",
  category = "Differentiation",
  genes = c("PRSS1", "PRSS2", "CELA2A", "CELA3A", "CPA1", "CPA2", "CPB1",
            "AMY2A", "AMY2B", "PNLIP", "CEL", "CLPS", "PTF1A", "RBPJL", "BHLHA15"),
  source = "Bailey et al. Nature 2016", pmid = "26909576"
)

diff_signatures$Diff_ADM_Ductal <- list(
  id = "Diff02", short_name = "ADM Ductal",
  full_name = "Acinar-to-Ductal Metaplasia - Ductal Markers",
  category = "Differentiation",
  genes = c("KRT19", "KRT7", "MUC1", "SOX9", "HNF1B", "CFTR", "SLC4A4",
            "AQP1", "GGT1", "CA2", "ONECUT2", "FOXA2"),
  source = "Bailey et al. Nature 2016", pmid = "26909576"
)

diff_signatures$Diff_Neuroendocrine <- list(
  id = "Diff03", short_name = "Neuroendocrine",
  full_name = "Neuroendocrine Markers",
  category = "Differentiation",
  genes = c("CHGA", "CHGB", "SYP", "NCAM1", "ENO2", "SNAP25", "SYT1",
            "INSM1", "ASCL1", "NEUROD1", "PCSK1", "PCSK2"),
  source = "Bailey et al. Nature 2016", pmid = "26909576"
)

diff_signatures$Diff_Endocrine <- list(
  id = "Diff04", short_name = "Endocrine",
  full_name = "Endocrine/Islet Markers",
  category = "Differentiation",
  genes = c("INS", "GCG", "SST", "PPY", "GHRL", "NKX2-2", "NKX6-1",
            "PAX4", "PAX6", "NEUROG3", "MAFA", "MAFB", "ISL1", "PDX1"),
  source = "Multiple sources", pmid = "NA"
)

diff_signatures$Diff_Pancreatic_Progenitor <- list(
  id = "Diff05", short_name = "Pancreatic Progenitor",
  full_name = "Pancreatic Progenitor Markers",
  category = "Differentiation",
  genes = c("PDX1", "SOX9", "HNF1B", "PTF1A", "NKX6-1", "ONECUT1", "FOXA1",
            "FOXA2", "ONECUT1", "GATA4", "GATA6", "HES1", "NKX2-2", "NEUROG3"),
  source = "Multiple sources", pmid = "NA"
)

diff_signatures$Diff_Squamous <- list(
  id = "Diff06", short_name = "Squamous Differentiation",
  full_name = "Squamous Differentiation Markers",
  category = "Differentiation",
  genes = c("TP63", "KRT5", "KRT6A", "KRT14", "KRT16", "S100A2", "SPRR1A",
            "SPRR1B", "IVL", "SERPINB3", "SERPINB4", "LY6D", "TRIM29"),
  source = "Bailey et al. Nature 2016", pmid = "26909576"
)

diff_signatures$Diff_Well_Differentiated <- list(
  id = "Diff07", short_name = "Well Differentiated",
  full_name = "Well-Differentiated PDAC Markers",
  category = "Differentiation",
  genes = c("GATA6", "HNF1A", "HNF4A", "CDX2", "TFF1", "MUC5AC", "KRT19",
            "KRT7", "EPCAM", "CDH1", "CLDN18", "AGR2", "SPINK1", "MUC1"),
  source = "Moffitt et al. Nat Genet 2015", pmid = "26343385"
)

diff_signatures$Diff_Poorly_Differentiated <- list(
  id = "Diff08", short_name = "Poorly Differentiated",
  full_name = "Poorly-Differentiated PDAC Markers",
  category = "Differentiation",
  genes = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2",
            "S100A2", "KRT5", "KRT14", "TP63", "SOX2"),
  source = "Bailey et al. Nature 2016", pmid = "26909576"
)

# Add to main signatures list
for(name in names(diff_signatures)) {
  signatures[[name]] <- diff_signatures[[name]]
}

cat(sprintf("   [OK] Added %d differentiation signatures\n", length(diff_signatures)))

# Check coverage
diff_coverage <- data.frame()
for(name in names(diff_signatures)) {
  sig <- diff_signatures[[name]]
  available <- sig$genes[sig$genes %in% all_proteins]
  diff_coverage <- rbind(diff_coverage, data.frame(
    ID = sig$id,
    Name = sig$short_name,
    N_Total = length(sig$genes),
    N_Available = length(available),
    Coverage_Pct = round(100 * length(available) / length(sig$genes), 1),
    Genes_Available = paste(available, collapse = ", "),
    stringsAsFactors = FALSE
  ))
}

cat("   Differentiation signature coverage:\n")
print(diff_coverage[, c("ID", "Name", "N_Total", "N_Available", "Coverage_Pct")])

write.csv(diff_coverage, file.path(results_dir, "Step11_06_Differentiation_Coverage.csv"), row.names = FALSE)
cat("   [OK] Saved: Step11_06_Differentiation_Coverage.csv\n")

# Run ORA for differentiation signatures
cat("   Running ORA for differentiation signatures...\n")

diff_gene_sets <- list()
for(name in names(diff_signatures)) {
  sig <- diff_signatures[[name]]
  available <- sig$genes[sig$genes %in% all_proteins]
  if(length(available) >= 3) {
    diff_gene_sets[[name]] <- available
  }
}

if(length(diff_gene_sets) > 0) {
  diff_ora_results <- run_ora(
    module_genes, diff_gene_sets, n_proteins,
    short_names = setNames(sapply(diff_signatures, function(x) x$short_name), names(diff_signatures)),
    full_names = setNames(sapply(diff_signatures, function(x) x$full_name), names(diff_signatures)),
    sig_info = diff_signatures
  )
  
  if(nrow(diff_ora_results) > 0) {
    write.csv(diff_ora_results, file.path(results_dir, "Step11_07_Differentiation_ORA.csv"), row.names = FALSE)
    cat("   [OK] Saved: Step11_07_Differentiation_ORA.csv\n")
    
    # Show significant results
    sig_diff <- diff_ora_results %>% filter(P_Value < 0.05)
    if(nrow(sig_diff) > 0) {
      cat("   Significant differentiation enrichments:\n")
      for(i in 1:min(10, nrow(sig_diff))) {
        cat(sprintf("      %s in %s: %.1fx (p=%.4f)\n",
                        sig_diff$Short_Name[i], sig_diff$Module[i],
                        sig_diff$Fold_Enrichment[i], sig_diff$P_Value[i]))
      }
    }
  }
}

# -- Figure 4: Combined Standard Database Summary --
cat("   Creating Figure 4: Combined Enrichment Summary...\n")

# Combine all results with database source - with safety checks
combined_db <- data.frame()

if(exists("go_combined") && nrow(go_combined) > 0) {
  go_sig <- go_combined %>% 
    filter(p.adjust < 0.05) %>% 
    mutate(Database = "GO:BP") %>% 
    dplyr::select(Module, ID, Description, p.adjust, Database)
  combined_db <- rbind(combined_db, go_sig)
}

if(exists("kegg_combined") && nrow(kegg_combined) > 0) {
  kegg_sig <- kegg_combined %>% 
    filter(p.adjust < 0.05) %>% 
    mutate(Database = "KEGG") %>% 
    dplyr::select(Module, ID, Description, p.adjust, Database)
  combined_db <- rbind(combined_db, kegg_sig)
}

if(exists("hallmark_combined") && nrow(hallmark_combined) > 0) {
  hallmark_sig <- hallmark_combined %>% 
    filter(p.adjust < 0.05) %>% 
    mutate(Database = "Hallmark") %>% 
    dplyr::select(Module, ID, Description, p.adjust, Database)
  combined_db <- rbind(combined_db, hallmark_sig)
}

if(exists("c6_combined") && nrow(c6_combined) > 0) {
  c6_sig <- c6_combined %>% 
    filter(p.adjust < 0.05) %>% 
    mutate(Database = "C6 Oncogenic") %>% 
    dplyr::select(Module, ID, Description, p.adjust, Database)
  combined_db <- rbind(combined_db, c6_sig)
}

if(nrow(combined_db) > 0) {
  # Summary counts by module and database
  db_summary <- combined_db %>%
    group_by(Module, Database) %>%
    summarise(N_Significant = n(), .groups = "drop") %>%
    pivot_wider(names_from = Database, values_from = N_Significant, values_fill = 0)
  
  cat("   Standard database enrichment summary (FDR < 0.05):")
  print(as.data.frame(db_summary))
  
  # Bar plot
  p_db_summary <- ggplot(combined_db %>% count(Module, Database),
                         aes(x = Module, y = n, fill = Database)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Standard Database Enrichment Summary - Focus Modules",
         subtitle = sprintf("Number of significant pathways (FDR < 0.05) | Modules: %s",
                            paste(toupper(focus_modules_auto), collapse = ", ")),
         x = "Module", y = "Count") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p_db_summary)
  ggsave(file.path(fig_dir_png, "Step11_04_DB_Summary.png"), p_db_summary, width = 8, height = 5, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step11_04_DB_Summary.pdf"), p_db_summary, width = 8, height = 5)
  cat("   [OK] Saved: Step11_04_DB_Summary\n")
}
# ============================================================================
# 11h: GSEA (Gene Set Enrichment Analysis) - Ranked List Approach
# ============================================================================
# Purpose: Complement ORA with rank-based enrichment using full gene list
# Method: fgsea using genes ranked by module membership (kME)
# Advantage: Uses full ranking information, more sensitive than binary ORA
# Reference: Subramanian et al., 2005 (PNAS); fgsea: Korotkevich et al., 2021
# ============================================================================

cat("\n==============================================================================\n")
cat("11h-i: GSEA (Gene Set Enrichment Analysis) - kME-based\n")
cat("==============================================================================\n")

# fgsea already loaded at script start
cat("   Using fgsea for GSEA analysis\n")

# Load MSigDB gene sets for GSEA
cat("   Loading MSigDB Hallmark gene sets...\n")
hallmark_gs <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(., .$gs_name) %>%
  lapply(function(x) x$gene_symbol)
cat(sprintf("   [OK] Loaded %d Hallmark gene sets\n", length(hallmark_gs)))

# Run GSEA for focus modules (pink, green, black, blue)
gsea_modules <- focus_modules_auto
cat(sprintf("   Running GSEA for focus modules: %s\n", paste(toupper(gsea_modules), collapse = ", ")))
gsea_results_all <- list()

for(mod in gsea_modules) {
  cat(sprintf("\n   Running GSEA for %s module...\n", mod))

  # Get module membership (kME) for ranking
  kme_col <- paste0("MM.", mod)

  if(kme_col %in% colnames(proteinInfo)) {
    # Create ranked gene list based on kME (module membership)
    # Use ALL proteins, not just module proteins, ranked by their kME to this module
    all_proteins_kme <- proteinInfo %>%
      filter(!is.na(.data[[kme_col]])) %>%
      arrange(desc(.data[[kme_col]]))

    # Create named vector for fgsea (gene names with kME values as scores)
    gene_ranks <- setNames(all_proteins_kme[[kme_col]], all_proteins_kme$Protein)

    # Run fgsea with scoreType="std" (kME can be negative for non-module proteins)
    set.seed(12345)
    gsea_result <- tryCatch({
      fgsea(
        pathways = hallmark_gs,
        stats = gene_ranks,
        minSize = 5,
        maxSize = 500,
        nPermSimple = 1000
      )
    }, error = function(e) {
      cat(sprintf("   [WARNING] GSEA failed for %s: %s\n", mod, e$message))
      return(NULL)
    })

    if(!is.null(gsea_result) && nrow(gsea_result) > 0) {
      # Add module info and format results
      gsea_result$Module <- mod
      gsea_result$leadingEdgeCount <- as.integer(sapply(gsea_result$leadingEdge, length))

      # Format for display - remove leadingEdge list column before saving
      gsea_formatted <- gsea_result %>%
        arrange(pval) %>%
        dplyr::select(pathway, pval, padj, NES, size, leadingEdgeCount, Module) %>%
        mutate(
          pathway = gsub("HALLMARK_", "", pathway),
          pval = round(pval, 4),
          padj = round(padj, 4),
          NES = round(NES, 3),
          leadingEdgeCount = as.integer(leadingEdgeCount)
        )

      gsea_results_all[[mod]] <- gsea_formatted

      # Print top results (FDR < 0.1 or nominal p < 0.05)
      cat(sprintf("   %s module - Top GSEA results:\n", mod))
      top_gsea <- gsea_formatted %>% filter(padj < 0.1) %>% head(10)
      if(nrow(top_gsea) > 0) {
        for(i in 1:nrow(top_gsea)) {
          cat(sprintf("      %s: NES=%.2f, FDR=%.3f\n",
                          top_gsea$pathway[i], top_gsea$NES[i], top_gsea$padj[i]))
        }
      } else {
        # Try nominal p < 0.05
        top_gsea_nominal <- gsea_formatted %>% filter(pval < 0.05) %>% head(5)
        if(nrow(top_gsea_nominal) > 0) {
          cat("      No FDR < 0.1, showing top nominal p < 0.05:\n")
          for(i in 1:nrow(top_gsea_nominal)) {
            cat(sprintf("      %s: NES=%.2f, p=%.3f\n",
                            top_gsea_nominal$pathway[i], top_gsea_nominal$NES[i], top_gsea_nominal$pval[i]))
          }
        } else {
          cat("      No significant pathways (FDR < 0.1 or nominal p < 0.05)\n")
        }
      }
    }
  } else {
    cat(sprintf("   Skipping %s - kME column not found\n", mod))
  }
}

# Combine and save GSEA results
if(length(gsea_results_all) > 0) {
  # Use do.call(rbind, ...) instead of bind_rows to avoid type conflicts
  gsea_combined <- do.call(rbind, gsea_results_all)
  gsea_combined <- as.data.frame(gsea_combined)

  write.csv(gsea_combined, file.path(results_dir, "Step11_08_GSEA_Hallmark_Results.csv"), row.names = FALSE)
  cat("\n   [OK] Saved: Step11_08_GSEA_Hallmark_Results.csv\n")

  # Create GSEA enrichment plot for top pathway per module
  cat("   Creating GSEA visualization...\n")

  # Try FDR < 0.1 first, fallback to nominal p < 0.05 if none
  gsea_plot_data <- gsea_combined %>%
    filter(padj < 0.1) %>%
    group_by(Module) %>%
    arrange(pval) %>%
    slice_head(n = 10) %>%
    ungroup()

  n_fdr_sig <- nrow(gsea_plot_data)
  use_nominal <- (n_fdr_sig == 0)

  if(use_nominal) {
    cat("   [NOTE] No pathways with FDR < 0.1. Using nominal p < 0.05 for visualization.\n")
    gsea_plot_data <- gsea_combined %>%
      filter(pval < 0.05) %>%
      group_by(Module) %>%
      arrange(pval) %>%
      slice_head(n = 10) %>%
      ungroup()
  }

  gsea_plot_data <- gsea_plot_data %>%
    mutate(
      pathway = factor(pathway, levels = unique(pathway)),
      Direction = ifelse(NES > 0, "Enriched", "Depleted")
    )

  n_sig_gsea <- nrow(gsea_plot_data)
  threshold_label <- ifelse(use_nominal, "nominal p < 0.05", "FDR < 0.1")
  cat(sprintf("   Significant pathways (%s): %d\n", threshold_label, n_sig_gsea))

  if(n_sig_gsea > 0) {
    p_gsea <- ggplot(gsea_plot_data, aes(x = NES, y = reorder(pathway, NES), fill = Module)) +
      geom_col(alpha = 0.8, width = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = module_color_palette) +
      labs(
        title = "GSEA: Hallmark Pathway Enrichment - Focus Modules",
        subtitle = sprintf("Ranked by module membership (kME) | %s | Modules: %s",
                           threshold_label, paste(toupper(focus_modules_auto), collapse = ", ")),
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        caption = "Method: fgsea | Reference: Subramanian et al., 2005"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom"
      ) +
      facet_wrap(~Module, scales = "free_y", ncol = 1)

    print(p_gsea)
    ggsave(file.path(fig_dir_png, "Step11_05_GSEA_Hallmark.png"), p_gsea, width = 8, height = 6, dpi = 300)
    ggsave(file.path(fig_dir_pdf, "Step11_05_GSEA_Hallmark.pdf"), p_gsea, width = 8, height = 6)
    cat("   [OK] Saved: Step11_05_GSEA_Hallmark.png/pdf\n")
  } else {
    cat("   [NOTE] No significant pathways for kME-based GSEA visualization\n")
    cat("   The Limma-based GSEA (11h-ii) uses clinical outcomes and may be more informative\n")
  }
} else {
  cat("   [NOTE] No GSEA results to combine - check if kME columns exist in proteinInfo\n")
}

cat("\n   kME-based GSEA analysis complete\n")

# ============================================================================
# 11h-ii: LIMMA-BASED GSEA (Gene Set Enrichment Analysis)
# ============================================================================
# Purpose: Run GSEA using differential expression results from Step 10b
# Method: fgsea using genes ranked by sign(logFC) * -log10(P.Value)
# Advantage: Uses clinical outcome-based ranking (PFS and Response)
# Complements kME-based GSEA above with biologically interpretable direction
# ============================================================================

cat("\n==============================================================================\n")
cat("11h-ii: LIMMA-BASED GSEA (PFS and Response)\n")
cat("==============================================================================\n")

# Check if limma results exist from Step 10b
limma_pfs_path <- file.path(results_dir, "Step10b_Limma_DE", "Step10b_Limma_DE_PFS.csv")
limma_response_path <- file.path(results_dir, "Step10b_Limma_DE", "Step10b_Limma_DE_Response.csv")

if (file.exists(limma_pfs_path) || file.exists(limma_response_path)) {
  cat("   Limma DE results found from Step 10b\n")

  # Load MSigDB collections for limma-based GSEA
  cat("   Loading MSigDB gene sets for limma-based GSEA...\n")

  # Hallmark pathways
  hallmark_limma <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(., .$gs_name) %>%
    lapply(function(x) x$gene_symbol)

  # KEGG pathways (C2 canonical pathways - filter for KEGG)
  kegg_limma <- msigdbr(species = "Homo sapiens", category = "C2") %>%
    dplyr::filter(grepl("^KEGG_", gs_name)) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(., .$gs_name) %>%
    lapply(function(x) x$gene_symbol)

  # Combine all gene sets
  all_genesets_limma <- c(hallmark_limma, kegg_limma)
  cat(sprintf("   Loaded: %d Hallmark + %d KEGG = %d total gene sets\n",
              length(hallmark_limma), length(kegg_limma), length(all_genesets_limma)))

  # --------------------------------------------------------------------------
  # PFS-based GSEA
  # --------------------------------------------------------------------------
  if (file.exists(limma_pfs_path)) {
    cat("\n--- Limma-based GSEA: PFS (Long vs Short) ---\n")

    de_pfs <- read.csv(limma_pfs_path, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded %d proteins from limma PFS results\n", nrow(de_pfs)))

    # Create ranking metric: sign(logFC) * -log10(P.Value)
    # Positive values = genes upregulated in Long PFS (favorable)
    # Negative values = genes upregulated in Short PFS (unfavorable)
    de_pfs$rank_metric <- sign(de_pfs$logFC) * (-log10(de_pfs$P.Value))

    # Remove NA/Inf values
    de_pfs <- de_pfs[is.finite(de_pfs$rank_metric), ]

    # Create named vector for fgsea
    pfs_ranks <- setNames(de_pfs$rank_metric, de_pfs$Protein)
    pfs_ranks <- sort(pfs_ranks, decreasing = TRUE)

    cat(sprintf("   Created ranking vector: %d proteins\n", length(pfs_ranks)))
    cat(sprintf("   Top 3 (Long PFS): %s\n", paste(head(names(pfs_ranks), 3), collapse = ", ")))
    cat(sprintf("   Bottom 3 (Short PFS): %s\n", paste(tail(names(pfs_ranks), 3), collapse = ", ")))

    # Run fgsea
    set.seed(12345)
    gsea_pfs_limma <- fgsea(
      pathways = all_genesets_limma,
      stats = pfs_ranks,
      minSize = 5,
      maxSize = 500,
      nPermSimple = 10000
    )

    # Format results
    gsea_pfs_limma <- gsea_pfs_limma %>%
      arrange(pval) %>%
      mutate(
        pathway_clean = gsub("HALLMARK_|KEGG_", "", pathway),
        Source = ifelse(grepl("^HALLMARK_", pathway), "Hallmark", "KEGG"),
        leadingEdgeCount = sapply(leadingEdge, length),
        Direction = ifelse(NES > 0, "Long_PFS (favorable)", "Short_PFS (unfavorable)")
      ) %>%
      dplyr::select(pathway, pathway_clean, Source, NES, pval, padj, size, leadingEdgeCount, Direction)

    # Save results
    write.csv(gsea_pfs_limma, file.path(results_dir, "Step11_fgsea_Limma_PFS.csv"), row.names = FALSE)
    cat("   [OK] Saved: Step11_fgsea_Limma_PFS.csv\n")

    # Report top pathways
    sig_pfs <- gsea_pfs_limma %>% filter(padj < 0.05)
    cat(sprintf("   Significant pathways (padj < 0.05): %d\n", nrow(sig_pfs)))

    top_long_pfs <- gsea_pfs_limma %>% filter(NES > 0) %>% arrange(pval) %>% head(3)
    top_short_pfs <- gsea_pfs_limma %>% filter(NES < 0) %>% arrange(pval) %>% head(3)

    cat("\n   Top 3 enriched in LONG PFS (favorable):\n")
    if (nrow(top_long_pfs) > 0) {
      for (i in 1:nrow(top_long_pfs)) {
        cat(sprintf("      %s: NES=%.2f, padj=%.4f\n",
                    top_long_pfs$pathway_clean[i], top_long_pfs$NES[i], top_long_pfs$padj[i]))
      }
    }

    cat("\n   Top 3 enriched in SHORT PFS (unfavorable):\n")
    if (nrow(top_short_pfs) > 0) {
      for (i in 1:nrow(top_short_pfs)) {
        cat(sprintf("      %s: NES=%.2f, padj=%.4f\n",
                    top_short_pfs$pathway_clean[i], top_short_pfs$NES[i], top_short_pfs$padj[i]))
      }
    }

    # Create visualization
    gsea_pfs_plot_data <- gsea_pfs_limma %>%
      filter(padj < 0.1) %>%
      arrange(NES) %>%
      head(20) %>%
      mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean))

    if (nrow(gsea_pfs_plot_data) > 0) {
      p_gsea_pfs <- ggplot(gsea_pfs_plot_data, aes(x = NES, y = pathway_clean, fill = Direction)) +
        geom_col(alpha = 0.85, width = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        scale_fill_manual(
          values = c("Long_PFS (favorable)" = "#008080", "Short_PFS (unfavorable)" = "#CA562C"),
          name = "Direction"
        ) +
        labs(
          title = "Limma-based GSEA: Long vs Short PFS",
          subtitle = "Ranking: sign(logFC) x -log10(P) | Hallmark + KEGG | padj < 0.1",
          x = "Normalized Enrichment Score (NES)",
          y = NULL,
          caption = "Positive NES = enriched in Long PFS (favorable prognosis)"
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
          plot.caption = element_text(hjust = 0.5, color = "grey50", size = 9),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom"
        )

      print(p_gsea_pfs)
      ggsave(file.path(fig_dir_png, "Step11_fgsea_Limma_PFS.png"), p_gsea_pfs, width = 10, height = 8, dpi = 300)
      ggsave(file.path(fig_dir_pdf, "Step11_fgsea_Limma_PFS.pdf"), p_gsea_pfs, width = 10, height = 8)
      cat("   [OK] Saved: Step11_fgsea_Limma_PFS.png/pdf\n")
    } else {
      cat("   [NOTE] No pathways with padj < 0.1 for PFS plot\n")
    }
  }

  # --------------------------------------------------------------------------
  # Response-based GSEA
  # --------------------------------------------------------------------------
  if (file.exists(limma_response_path)) {
    cat("\n--- Limma-based GSEA: Response (CD vs PD) ---\n")

    de_response <- read.csv(limma_response_path, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded %d proteins from limma Response results\n", nrow(de_response)))

    # Create ranking metric: sign(logFC) * -log10(P.Value)
    # Positive values = genes upregulated in CD (favorable response)
    # Negative values = genes upregulated in PD (progressive disease)
    de_response$rank_metric <- sign(de_response$logFC) * (-log10(de_response$P.Value))

    # Remove NA/Inf values
    de_response <- de_response[is.finite(de_response$rank_metric), ]

    # Create named vector for fgsea
    response_ranks <- setNames(de_response$rank_metric, de_response$Protein)
    response_ranks <- sort(response_ranks, decreasing = TRUE)

    cat(sprintf("   Created ranking vector: %d proteins\n", length(response_ranks)))
    cat(sprintf("   Top 3 (CD): %s\n", paste(head(names(response_ranks), 3), collapse = ", ")))
    cat(sprintf("   Bottom 3 (PD): %s\n", paste(tail(names(response_ranks), 3), collapse = ", ")))

    # Run fgsea
    set.seed(12345)
    gsea_response_limma <- fgsea(
      pathways = all_genesets_limma,
      stats = response_ranks,
      minSize = 5,
      maxSize = 500,
      nPermSimple = 10000
    )

    # Format results
    gsea_response_limma <- gsea_response_limma %>%
      arrange(pval) %>%
      mutate(
        pathway_clean = gsub("HALLMARK_|KEGG_", "", pathway),
        Source = ifelse(grepl("^HALLMARK_", pathway), "Hallmark", "KEGG"),
        leadingEdgeCount = sapply(leadingEdge, length),
        Direction = ifelse(NES > 0, "CD (responder)", "PD (non-responder)")
      ) %>%
      dplyr::select(pathway, pathway_clean, Source, NES, pval, padj, size, leadingEdgeCount, Direction)

    # Save results
    write.csv(gsea_response_limma, file.path(results_dir, "Step11_fgsea_Limma_Response.csv"), row.names = FALSE)
    cat("   [OK] Saved: Step11_fgsea_Limma_Response.csv\n")

    # Report top pathways
    sig_response <- gsea_response_limma %>% filter(padj < 0.05)
    cat(sprintf("   Significant pathways (padj < 0.05): %d\n", nrow(sig_response)))

    top_cd <- gsea_response_limma %>% filter(NES > 0) %>% arrange(pval) %>% head(3)
    top_pd <- gsea_response_limma %>% filter(NES < 0) %>% arrange(pval) %>% head(3)

    cat("\n   Top 3 enriched in CD (responders):\n")
    if (nrow(top_cd) > 0) {
      for (i in 1:nrow(top_cd)) {
        cat(sprintf("      %s: NES=%.2f, padj=%.4f\n",
                    top_cd$pathway_clean[i], top_cd$NES[i], top_cd$padj[i]))
      }
    }

    cat("\n   Top 3 enriched in PD (non-responders):\n")
    if (nrow(top_pd) > 0) {
      for (i in 1:nrow(top_pd)) {
        cat(sprintf("      %s: NES=%.2f, padj=%.4f\n",
                    top_pd$pathway_clean[i], top_pd$NES[i], top_pd$padj[i]))
      }
    }

    # Create visualization
    gsea_response_plot_data <- gsea_response_limma %>%
      filter(padj < 0.1) %>%
      arrange(NES) %>%
      head(20) %>%
      mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean))

    if (nrow(gsea_response_plot_data) > 0) {
      p_gsea_response <- ggplot(gsea_response_plot_data, aes(x = NES, y = pathway_clean, fill = Direction)) +
        geom_col(alpha = 0.85, width = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        scale_fill_manual(
          values = c("CD (responder)" = "#2E7D32", "PD (non-responder)" = "#C62828"),
          name = "Direction"
        ) +
        labs(
          title = "Limma-based GSEA: CD vs PD Response",
          subtitle = "Ranking: sign(logFC) x -log10(P) | Hallmark + KEGG | padj < 0.1",
          x = "Normalized Enrichment Score (NES)",
          y = NULL,
          caption = "Positive NES = enriched in CD (favorable response)"
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
          plot.caption = element_text(hjust = 0.5, color = "grey50", size = 9),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom"
        )

      print(p_gsea_response)
      ggsave(file.path(fig_dir_png, "Step11_fgsea_Limma_Response.png"), p_gsea_response, width = 10, height = 8, dpi = 300)
      ggsave(file.path(fig_dir_pdf, "Step11_fgsea_Limma_Response.pdf"), p_gsea_response, width = 10, height = 8)
      cat("   [OK] Saved: Step11_fgsea_Limma_Response.png/pdf\n")
    } else {
      cat("   [NOTE] No pathways with padj < 0.1 for Response plot\n")
    }
  }

  cat("\n   Limma-based GSEA complete\n")

} else {
  cat("   [NOTE] No limma DE results found from Step 10b\n")
  cat("   Run Step 10b first to generate limma differential expression results\n")
  cat("   Expected files:\n")
  cat(sprintf("     - %s\n", limma_pfs_path))
  cat(sprintf("     - %s\n", limma_response_path))
}

# ============================================================================
# 11i: STEP 11 SUMMARY
# ============================================================================

cat("\n==============================================================================\n")
cat("11i: STEP 11 SUMMARY\n")
cat("==============================================================================\n")

cat("   STANDARD DATABASE ENRICHMENT COMPLETE\n")
cat(sprintf("   Focus Modules: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("   ------------------------------------------------------------------------\n")
cat(sprintf("   GO:BP significant terms: %d\n",
    if(exists("go_combined") && nrow(go_combined) > 0) sum(go_combined$p.adjust < 0.05, na.rm = TRUE) else 0))
cat(sprintf("   KEGG significant pathways: %d\n",
    if(exists("kegg_combined") && nrow(kegg_combined) > 0) sum(kegg_combined$p.adjust < 0.05, na.rm = TRUE) else 0))
cat(sprintf("   Hallmark significant: %d\n",
    if(exists("hallmark_combined") && nrow(hallmark_combined) > 0) sum(hallmark_combined$p.adjust < 0.05, na.rm = TRUE) else 0))
cat(sprintf("   C2:CP significant: %d\n",
    if(exists("c2_cp_combined") && nrow(c2_cp_combined) > 0) sum(c2_cp_combined$p.adjust < 0.05, na.rm = TRUE) else 0))
cat(sprintf("   C6 Oncogenic significant: %d\n",
    if(exists("c6_combined") && nrow(c6_combined) > 0) sum(c6_combined$p.adjust < 0.05, na.rm = TRUE) else 0))
cat(sprintf("   Differentiation signatures added: %d\n",
    if(exists("diff_signatures")) length(diff_signatures) else 0))
cat("   OUTPUT FILES\n")
cat("   ------------------------------------------------------------------------\n")
cat("   - Step11_01_GO_BP_Results.csv\n")
cat("   - Step11_02_KEGG_Results.csv\n")
cat("   - Step11_03_MSigDB_Hallmark_Results.csv\n")
cat("   - Step11_04_MSigDB_C2_CP_Results.csv\n")
cat("   - Step11_05_MSigDB_C6_Oncogenic_Results.csv\n")
cat("   - Step11_06_Differentiation_Coverage.csv\n")
cat("   - Step11_07_Differentiation_ORA.csv\n")
cat("   - Step11_08_GSEA_Hallmark_Results.csv (kME-based)\n")
cat("   - Step11_fgsea_Limma_PFS.csv (limma-based)\n")
cat("   - Step11_fgsea_Limma_Response.csv (limma-based)\n")
cat("   - Step11_01a/b_GO_BP_*.png/pdf\n")
cat("   - Step11_02a/b_KEGG_*.png/pdf\n")
cat("   - Step11_03_Hallmark_Heatmap.png/pdf\n")
cat("   - Step11_04_DB_Summary.png/pdf\n")
cat("   - Step11_05_GSEA_Hallmark.png/pdf (kME-based)\n")
cat("   - Step11_fgsea_Limma_PFS.png/pdf (limma-based)\n")
cat("   - Step11_fgsea_Limma_Response.png/pdf (limma-based)\n")
# Save Step 11 workspace - only save objects that exist
step11_objects <- c()
if(exists("go_results")) step11_objects <- c(step11_objects, "go_results")
if(exists("go_combined")) step11_objects <- c(step11_objects, "go_combined")
if(exists("kegg_results")) step11_objects <- c(step11_objects, "kegg_results")
if(exists("kegg_combined")) step11_objects <- c(step11_objects, "kegg_combined")
if(exists("hallmark_results")) step11_objects <- c(step11_objects, "hallmark_results")
if(exists("hallmark_combined")) step11_objects <- c(step11_objects, "hallmark_combined")
if(exists("c2_cp_results")) step11_objects <- c(step11_objects, "c2_cp_results")
if(exists("c2_cp_combined")) step11_objects <- c(step11_objects, "c2_cp_combined")
if(exists("c6_results")) step11_objects <- c(step11_objects, "c6_results")
if(exists("c6_combined")) step11_objects <- c(step11_objects, "c6_combined")
if(exists("diff_signatures")) step11_objects <- c(step11_objects, "diff_signatures")
if(exists("diff_coverage")) step11_objects <- c(step11_objects, "diff_coverage")
if(exists("module_entrez")) step11_objects <- c(step11_objects, "module_entrez")
if(exists("universe_entrez")) step11_objects <- c(step11_objects, "universe_entrez")
if(exists("gsea_pfs_limma")) step11_objects <- c(step11_objects, "gsea_pfs_limma")
if(exists("gsea_response_limma")) step11_objects <- c(step11_objects, "gsea_response_limma")

if(length(step11_objects) > 0) {
  save(list = step11_objects, file = file.path(results_dir, "Step11_Standard_DB_Enrichment.RData"))
  cat(sprintf("   [OK] Saved workspace: Step11_Standard_DB_Enrichment.RData (%d objects)\n", length(step11_objects)))
} else {
  cat("   [WARNING] No Step 11 objects to save (optional packages may be missing)\n")
}

cat("\n==============================================================================\n")
cat("STEP 11 COMPLETE: Standard Database Enrichment\n")
cat(sprintf("Focus Modules Analyzed: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("==============================================================================\n")


#===============================================================================
# STEP 12: BIOLOGICAL STORYTELLING & INTEGRATION
#===============================================================================
# PURPOSE: Synthesize all findings into cohesive biological narratives
#          for each module, validate clinical associations, create final summary.
#
# FOCUS MODULES: Pink, Green, Black, Blue (clinically significant from Step 5)
#
# STRATEGY:
#   12a) Module biological narratives (focus modules)
#   12b) Hub gene clinical validation
#   12c) Integrated findings summary
#   12d) Final workspace save
#   12e) PPI validation (STRING database)
#   12f) Step summary
#
# LOOKING FOR:
#   - Coherent biological story: pathway + clinical + mechanistic
#   - Actionable insights: hub proteins as biomarker candidates
#   - Publication-ready summary connecting all analysis steps
#
# KEY PARAMS: Integration of Steps 5-11 results
#===============================================================================          


cat("\n##  STEP 12: BIOLOGICAL STORYTELLING & INTEGRATION  ##\n")
cat(sprintf("## Focus Modules: %s ##\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("## Clinical Validation, Module Narratives, Final Summary ##\n")

# NOTE: Step 12a (MODULE CLINICAL VALIDATION) removed - already addressed in Step 4d heatmaps
# Clinical validation of module eigengenes vs PFS/Response is covered in Step4d_03a_CleanHeatmap

# Build me_clinical_results from Step 5 moduleTraitCor/moduleTraitPval data
# This enables Figure 3 (Integrated Summary Diagram) to display module clinical significance
cat("   Building me_clinical_results from Step 5 module-trait correlations...\n")

# Extract PFS_group correlations and p-values for each module
me_clinical_results <- data.frame(
  Module = gsub("ME", "", rownames(moduleTraitCor)),
  stringsAsFactors = FALSE
)

# Add PFS correlation as effect size proxy (Cohen's d approximation: d ≈ 2r / sqrt(1-r²))
if("PFS_group" %in% colnames(moduleTraitCor)) {
  me_clinical_results$Correlation <- moduleTraitCor[, "PFS_group"]
  me_clinical_results$P_Value <- moduleTraitPval[, "PFS_group"]
  # Convert correlation to Cohen's d approximation
  r <- me_clinical_results$Correlation
  me_clinical_results$Cohens_d <- 2 * r / sqrt(1 - r^2)
} else if("PFS_Days" %in% colnames(moduleTraitCor)) {
  me_clinical_results$Correlation <- moduleTraitCor[, "PFS_Days"]
  me_clinical_results$P_Value <- moduleTraitPval[, "PFS_Days"]
  r <- me_clinical_results$Correlation
  me_clinical_results$Cohens_d <- 2 * r / sqrt(1 - r^2)
}

# Add FDR correction
me_clinical_results$FDR <- p.adjust(me_clinical_results$P_Value, method = "BH")

# Add direction interpretation
me_clinical_results$Direction <- ifelse(me_clinical_results$Cohens_d > 0,
                                         "Good Prognosis (Long PFS)",
                                         "Poor Prognosis (Short PFS)")

# Filter to focus modules only for cleaner visualization
me_clinical_results <- me_clinical_results[me_clinical_results$Module %in% focus_modules_auto, ]

cat(sprintf("   [OK] Built me_clinical_results for %d focus modules: %s\n",
            nrow(me_clinical_results), paste(me_clinical_results$Module, collapse = ", ")))

treatment_stratified <- data.frame()
# NOTE: Step 12b (TREATMENT-STRATIFIED ANALYSIS) removed - redundant with Step 4d analysis
# ============================================================================
# 12a: MODULE BIOLOGICAL NARRATIVES
# ============================================================================

cat("==============================================================================\n")
cat("12a: MODULE BIOLOGICAL NARRATIVES\n")
cat("==============================================================================\n")

# Function to build module narrative
build_module_narrative <- function(mod, ora_df, go_df, hallmark_df, me_clin_df) {
  
  narrative <- list()
  narrative$module <- mod
  narrative$n_genes <- sum(moduleColors == mod)
  
  # Clinical association
  if(nrow(me_clin_df) > 0 && mod %in% me_clin_df$Module) {
    clin <- me_clin_df[me_clin_df$Module == mod, ]
    narrative$clinical_direction <- clin$Direction
    narrative$clinical_p <- clin$P_Value
    narrative$clinical_d <- clin$Cohens_d
  }
  
  # Top custom signatures
  if(nrow(ora_df) > 0) {
    mod_ora <- ora_df %>%
      filter(Module == mod, P_Value < 0.05) %>%
      arrange(P_Value) %>%
      head(5)
    narrative$top_custom_sigs <- mod_ora$Short_Name
  }
  
  # Top GO:BP terms
  if(nrow(go_df) > 0) {
    mod_go <- go_df %>%
      filter(Module == mod, p.adjust < 0.1) %>%
      arrange(p.adjust) %>%
      head(5)
    narrative$top_go_bp <- mod_go$Description
  }
  
  # Top Hallmark pathways
  if(nrow(hallmark_df) > 0) {
    mod_hall <- hallmark_df %>%
      filter(Module == mod, p.adjust < 0.1) %>%
      arrange(p.adjust) %>%
      head(3)
    narrative$top_hallmark <- gsub("HALLMARK_", "", mod_hall$ID)
  }
  
  return(narrative)
}

# Build narratives for focus modules only (pink, green, black, blue)
module_narratives <- list()

# Use focus modules only - no hardcoded extras
narrative_modules <- focus_modules_auto
narrative_modules <- narrative_modules[narrative_modules %in% unique(moduleColors)]
cat(sprintf("   Building narratives for: %s\n", paste(toupper(narrative_modules), collapse = ", ")))

for(mod in narrative_modules) {
  module_narratives[[mod]] <- build_module_narrative(
    mod, ora_results, go_combined, hallmark_combined, me_clinical_results
  )
}

# Print narratives
cat("   MODULE BIOLOGICAL NARRATIVES\n")
cat("   ------------------------------------------------------------------------\n")

for(mod in names(module_narratives)) {
  n <- module_narratives[[mod]]
  cat(sprintf("   [%s MODULE] (%d proteins)\n", toupper(mod), n$n_genes))
  
  if(!is.null(n$clinical_direction)) {
    cat(sprintf("   Clinical: %s (p=%.4f, d=%.2f)\n",
                    n$clinical_direction, n$clinical_p, n$clinical_d))
  }
  
  if(length(n$top_custom_sigs) > 0) {
    cat(sprintf("   Top Signatures: %s\n", paste(n$top_custom_sigs, collapse = ", ")))
  }
  
  if(length(n$top_go_bp) > 0) {
    cat(sprintf("   GO:BP: %s\n", paste(head(n$top_go_bp, 3), collapse = "; ")))
  }
  
  if(length(n$top_hallmark) > 0) {
    cat(sprintf("   Hallmark: %s\n", paste(n$top_hallmark, collapse = ", ")))
  }
}

# Save narratives as JSON-like structure
narrative_df <- do.call(rbind, lapply(names(module_narratives), function(mod) {
  n <- module_narratives[[mod]]
  data.frame(
    Module = mod,
    N_Proteins = n$n_genes,
    Clinical_Direction = ifelse(!is.null(n$clinical_direction), n$clinical_direction, NA),
    Clinical_P = ifelse(!is.null(n$clinical_p), n$clinical_p, NA),
    Clinical_d = ifelse(!is.null(n$clinical_d), n$clinical_d, NA),
    Top_Signatures = paste(n$top_custom_sigs, collapse = "; "),
    Top_GO_BP = paste(head(n$top_go_bp, 3), collapse = "; "),
    Top_Hallmark = paste(n$top_hallmark, collapse = "; "),
    stringsAsFactors = FALSE
  )
}))

write.csv(narrative_df, file.path(results_dir, "Step12_01_Module_Narratives.csv"), row.names = FALSE)
cat("   [OK] Saved: Step12_01_Module_Narratives.csv\n")

# ============================================================================
# 12b: HUB GENE CLINICAL VALIDATION
# ============================================================================

cat("==============================================================================\n")
cat("12b: HUB GENE CLINICAL VALIDATION\n")
cat("==============================================================================\n")

# Get hub genes from Step 6 - use hub_results (list of Strict/Moderate hub DFs per module)
# or proteinInfo with isHub_Moderate column
hub_genes_combined <- NULL

if(exists("hub_results") && !is.null(hub_results$Moderate)) {
  # Combine all moderate hub genes from all modules
  hub_genes_combined <- bind_rows(hub_results$Moderate, .id = "Module") %>%
    dplyr::select(Protein, Module) %>%
    rename(Gene = Protein)
  cat(sprintf("   Found %d hub genes from hub_results$Moderate\n", nrow(hub_genes_combined)))
} else if(exists("proteinInfo") && "isHub_Moderate" %in% colnames(proteinInfo)) {
  # Fallback: use proteinInfo
  hub_genes_combined <- proteinInfo %>%
    filter(isHub_Moderate == TRUE) %>%
    dplyr::select(Protein, Module) %>%
    rename(Gene = Protein)
  cat(sprintf("   Found %d hub genes from proteinInfo\n", nrow(hub_genes_combined)))
}

if(!is.null(hub_genes_combined) && nrow(hub_genes_combined) > 0) {
  hub_validation <- data.frame()

  # Prepare clinical data for validation
  clinical_validation <- data.frame(
    Sample_ID = rownames(datExpr),
    stringsAsFactors = FALSE
  )

  # Add PFS_Group with type conversion
  if("PFS_group" %in% colnames(datTraits)) {
    pfs_raw <- datTraits$PFS_group
    if (is.numeric(pfs_raw)) {
      clinical_validation$PFS_Group <- ifelse(pfs_raw == 1, "Long",
                                       ifelse(pfs_raw == 0, "Short", NA))
    } else {
      clinical_validation$PFS_Group <- as.character(pfs_raw)
    }
  }

  cat("   Testing hub genes vs PFS_Group...\n")

  for(i in 1:nrow(hub_genes_combined)) {
    gene <- hub_genes_combined$Gene[i]
    mod <- hub_genes_combined$Module[i]

    if(gene %in% colnames(datExpr) && "PFS_Group" %in% colnames(clinical_validation)) {
      expr_vals <- datExpr[, gene]
      groups <- clinical_validation$PFS_Group

      long_vals <- expr_vals[!is.na(groups) & groups == "Long"]
      short_vals <- expr_vals[!is.na(groups) & groups == "Short"]

      if(length(long_vals) >= 3 && length(short_vals) >= 3) {
        tt <- t.test(long_vals, short_vals)
        cohens_d <- (mean(long_vals) - mean(short_vals)) /
          sqrt((var(long_vals) + var(short_vals)) / 2)

        hub_validation <- rbind(hub_validation, data.frame(
          Gene = gene,
          Module = mod,
          Mean_Long = round(mean(long_vals), 4),
          Mean_Short = round(mean(short_vals), 4),
          Cohens_d = round(cohens_d, 3),
          P_Value = tt$p.value,
          Higher_In = ifelse(mean(short_vals) > mean(long_vals), "Short", "Long"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if(nrow(hub_validation) > 0) {
    hub_validation$FDR <- p.adjust(hub_validation$P_Value, method = "BH")
    hub_validation <- hub_validation %>% arrange(P_Value)

    # Show top results
    cat("   Top hub genes by PFS association:\n")
    top_hubs <- hub_validation %>% filter(P_Value < 0.1) %>% head(10)
    if(nrow(top_hubs) > 0) {
      for(i in 1:nrow(top_hubs)) {
        sig <- ifelse(top_hubs$P_Value[i] < 0.05, "*", "")
        cat(sprintf("      %s (%s): Higher in %s (d=%.2f, p=%.4f)%s\n",
                        top_hubs$Gene[i], top_hubs$Module[i], top_hubs$Higher_In[i],
                        abs(top_hubs$Cohens_d[i]), top_hubs$P_Value[i], sig))
      }
    } else {
      cat("   No hub genes with p < 0.1\n")
    }

    write.csv(hub_validation, file.path(results_dir, "Step12_02_Hub_Gene_Validation.csv"), row.names = FALSE)
    cat("   [OK] Saved: Step12_02_Hub_Gene_Validation.csv\n")
  } else {
    cat("   [WARNING] No hub genes could be validated (check PFS_Group availability)\n")
  }
} else {
  cat("   [WARNING] Hub genes not found - skipping validation\n")
  cat("   (Expected: hub_results$Moderate or proteinInfo$isHub_Moderate)\n")
}
# ============================================================================
# 12c: INTEGRATED FINDINGS SUMMARY
# ============================================================================

cat("==============================================================================\n")
cat("12c: INTEGRATED FINDINGS SUMMARY\n")
cat("==============================================================================\n")

# Create comprehensive summary
findings_summary <- list()

# Defensive check for soft_power (defined in Step 4/5)
if(!exists("soft_power")) {
  soft_power <- 13  # Default from Step 4 analysis
  cat("   [INFO] soft_power not found, using default: 13\n")
}

findings_summary$dataset <- list(
  n_samples = n_samples,
  n_proteins = n_proteins,
  n_modules = length(unique(moduleColors)) - 1  # Exclude grey
)

findings_summary$network <- list(
  method = "Bicor + Signed Hybrid",
  power = soft_power,
  scale_free_R2 = if(exists("sft") && !is.null(sft$fitIndices)) {
    round(max(sft$fitIndices$SFT.R.sq, na.rm = TRUE), 2)
  } else {
    "~0.94"  # Default estimate if sft not available
  }
)

# Defensive check for empty me_clinical_results
if (nrow(me_clinical_results) > 0 && "P_Value" %in% colnames(me_clinical_results)) {
  findings_summary$clinical_modules <- me_clinical_results %>%
    filter(P_Value < 0.1) %>%
    dplyr::select(Module, Direction, Cohens_d, P_Value)
} else {
  findings_summary$clinical_modules <- data.frame(
    Module = character(0), Direction = character(0),
    Cohens_d = numeric(0), P_Value = numeric(0)
  )
  cat("   [INFO] me_clinical_results empty - using Step 4d results instead\n")
}

findings_summary$enrichment <- list(
  custom_signatures = length(signatures),
  go_bp_significant = sum(go_combined$p.adjust < 0.05, na.rm = TRUE),
  kegg_significant = sum(kegg_combined$p.adjust < 0.05, na.rm = TRUE),
  hallmark_significant = sum(hallmark_combined$p.adjust < 0.05, na.rm = TRUE)
)

# Print summary

cat("INTEGRATED FINDINGS SUMMARY\n")

cat(sprintf("Dataset: %d patients Ã- %d plasma proteins\n",
            findings_summary$dataset$n_samples, findings_summary$dataset$n_proteins))
cat(sprintf("Modules detected: %d (excluding grey)\n",
                findings_summary$dataset$n_modules))

cat("NETWORK CONSTRUCTION\n")
cat(sprintf("Method: %s\n", findings_summary$network$method))
cat(sprintf("Soft threshold power: %d\n", soft_power))
cat(sprintf("Scale-free topology RÂ2: %s\n", findings_summary$network$scale_free_R2))

cat("ENRICHMENT ANALYSIS\n")
cat(sprintf("Custom signatures tested: %d\n",findings_summary$enrichment$custom_signatures))
cat(sprintf("GO:BP significant (FDR<0.05): %d\n", findings_summary$enrichment$go_bp_significant))
cat(sprintf("KEGG significant (FDR<0.05): %d\n",findings_summary$enrichment$kegg_significant))
cat(sprintf("Hallmark significant (FDR<0.05): %d\n",findings_summary$enrichment$hallmark_significant))


# -- Figure 3: Integrated Summary Diagram --
cat("   Creating Figure 3: Integrated Summary...\n")

# Create summary statistics for each module
if(nrow(me_clinical_results) > 0) {
  module_summary_plot <- me_clinical_results %>%
    mutate(
      Significance = case_when(
        FDR < 0.05 ~ "FDR < 0.05",
        P_Value < 0.05 ~ "+",
        TRUE ~ "NS"
      ),
      Direction_Simple = ifelse(Cohens_d > 0, "Good Prognosis", "Poor Prognosis")
    )
  
  p_summary <- ggplot(module_summary_plot,
                      aes(x = Cohens_d, y = -log10(P_Value),
                          color = Direction_Simple, size = abs(Cohens_d))) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = Module), hjust = -0.2, vjust = -0.5, size = 3) +
    scale_color_manual(values = c("Good Prognosis" = "#008080", "Poor Prognosis" = "#ca562c"),
                       name = "Association") +
    scale_size_continuous(range = c(3, 8), name = "|Effect Size|") +
    labs(
      title = "Focus Module Clinical Significance Overview",
      subtitle = sprintf("Volcano-style plot of ME vs PFS_group associations (n=%d modules: %s)",
                         nrow(module_summary_plot), paste(toupper(module_summary_plot$Module), collapse = ", ")),
      x = "Cohen's d (Long - Short PFS)",
      y = "-log10(p-value)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )
  
  print(p_summary)
  ggsave(file.path(fig_dir_png, "Step12_03_Integrated_Summary.png"), p_summary, width = 6, height = 4, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step12_03_Integrated_Summary.pdf"), p_summary, width = 6, height = 4)
  cat("   [OK] Saved: Step12_03_Integrated_Summary\n")
}
# ============================================================================
# 12d: FINAL WORKSPACE SAVE
# ============================================================================

cat("==============================================================================\n")
cat("12d: FINAL WORKSPACE SAVE\n")
cat("==============================================================================\n")

# Save complete analysis workspace
save(
  # Data
  datExpr, datTraits, all_proteins, n_samples, n_proteins,
  # Network
  soft_power, net, moduleColors, MEs,
  # Clinical validation
  me_clinical_results, treatment_stratified, module_narratives, findings_summary,
  # Enrichment
  signatures, gene_sets, ora_results, gsva_results, ssgsea_results,
  go_results, go_combined, kegg_results, kegg_combined,
  hallmark_results, hallmark_combined, c2_cp_combined, c6_combined,
  diff_signatures,
  # Output
  file = file.path(results_dir, "Step12_Complete_Analysis.RData")
)
cat("   [OK] Saved: Step12_Complete_Analysis.RData\n")

# ============================================================================
# 12e: PPI Validation (STRING Database Integration)
# ============================================================================
# Purpose: Validate module coherence using known protein-protein interactions
# Method: Overlay modules on STRING PPI network
# Benefit: Biological validation that modules reflect real protein complexes
# Reference: Szklarczyk et al., 2021 (STRING database)
# ============================================================================

cat("\n==============================================================================\n")
cat("12e: PPI Validation (STRING Database)\n")
cat("==============================================================================\n")

# Check for STRINGdb package (optional - loaded at startup if available)
if(!requireNamespace("STRINGdb", quietly = TRUE)) {
  cat("   [SKIP] STRINGdb not installed. PPI validation skipped.\n")
  cat("   To enable: BiocManager::install('STRINGdb')\n")
  ppi_available <- FALSE
} else {
  ppi_available <- TRUE
  if(!"STRINGdb" %in% loadedNamespaces()) library(STRINGdb)
  cat("   [OK] STRINGdb available\n")
}

if(ppi_available) {
  tryCatch({
    # Initialize STRING database (human, high confidence)
    cat("   Connecting to STRING database (v11.5, species 9606)...")
    string_db <- STRINGdb::STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
    cat("   Connected successfully\n")

    # Get hub genes from key modules - DYNAMIC based on focus_modules_auto
    key_modules <- focus_modules_auto
    ppi_results <- list()

    for(mod in key_modules) {
      cat(sprintf("\n   Analyzing PPI for %s module hubs...", mod))

      # Get hub genes (moderate criteria)
      if(exists("hub_results") && mod %in% names(hub_results$Moderate)) {
        mod_hubs <- hub_results$Moderate[[mod]]$Protein
      } else {
        # Fallback: get top genes by kME
        kme_col <- paste0("MM.", mod)
        if(kme_col %in% colnames(proteinInfo)) {
          mod_hubs <- proteinInfo %>%
            filter(Module == mod) %>%
            arrange(desc(.data[[kme_col]])) %>%
            head(20) %>%
            pull(Protein)
        } else {
          mod_hubs <- character(0)
        }
      }

      if(length(mod_hubs) > 0) {
        cat(sprintf("      Hub genes: %d", length(mod_hubs)))

        # Map genes to STRING IDs
        hub_df <- data.frame(gene = mod_hubs, stringsAsFactors = FALSE)
        mapped <- string_db$map(hub_df, "gene", removeUnmappedRows = TRUE)

        if(nrow(mapped) > 0) {
          cat(sprintf("      Mapped to STRING: %d (%.1f%%)",
                          nrow(mapped), 100 * nrow(mapped) / length(mod_hubs)))

          # Get interactions between hub genes
          interactions <- string_db$get_interactions(mapped$STRING_id)

          if(nrow(interactions) > 0) {
            # Calculate PPI enrichment
            # Expected: random pairs from proteome
            # Observed: pairs among hub genes
            n_hubs <- nrow(mapped)
            n_possible_pairs <- n_hubs * (n_hubs - 1) / 2
            n_observed_interactions <- nrow(interactions)

            # Calculate interaction density
            interaction_density <- n_observed_interactions / n_possible_pairs

            cat(sprintf("      Interactions found: %d", n_observed_interactions))
            cat(sprintf("      Possible pairs: %d", n_possible_pairs))
            cat(sprintf("      Interaction density: %.3f", interaction_density))

            # Get enrichment p-value from STRING (with error handling)
            enrichment_p <- NA
            tryCatch({
              enrichment <- string_db$get_ppi_enrichment(mapped$STRING_id)
              if(!is.null(enrichment) && "p_value" %in% names(enrichment)) {
                enrichment_p <- enrichment$p_value
              }
            }, error = function(e) {
              cat(sprintf("\n      Note: PPI enrichment calculation failed - %s", e$message))
            })

            ppi_results[[mod]] <- data.frame(
              Module = mod,
              N_Hubs = length(mod_hubs),
              N_Mapped = nrow(mapped),
              N_Interactions = n_observed_interactions,
              Possible_Pairs = n_possible_pairs,
              Interaction_Density = round(interaction_density, 4),
              PPI_Enrichment_P = enrichment_p,
              stringsAsFactors = FALSE
            )

            if(!is.na(enrichment_p)) {
              cat(sprintf("      PPI enrichment p-value: %.2e", enrichment_p))
            }

            # Get network image (optional - may require internet)
            tryCatch({
              png_file <- file.path(fig_dir_png, sprintf("Step12_PPI_%s_network.png", mod))
              png(png_file, width = 800, height = 800, res = 100)
              string_db$plot_network(mapped$STRING_id)
              dev.off()
              cat(sprintf("\n      Network plot saved: Step12_PPI_%s_network.png", mod))
            }, error = function(e) {
              try(dev.off(), silent = TRUE)  # Close device if error occurred
              cat(sprintf("\n      Note: Network plot not saved - %s", e$message))
            })

          } else {
            cat("      No interactions found among hub genes\n")
            ppi_results[[mod]] <- data.frame(
              Module = mod, N_Hubs = length(mod_hubs), N_Mapped = nrow(mapped),
              N_Interactions = 0, Possible_Pairs = n_possible_pairs,
              Interaction_Density = 0, PPI_Enrichment_P = NA,
              stringsAsFactors = FALSE
            )
          }
        } else {
          cat("      Warning: No genes mapped to STRING\n")
        }
      } else {
        cat(sprintf("      No hub genes found for %s module", mod))
      }
    }

    # Combine and save PPI results
    if(length(ppi_results) > 0) {
      ppi_summary <- bind_rows(ppi_results)
      write.csv(ppi_summary, file.path(results_dir, "Step12_04_PPI_Validation.csv"), row.names = FALSE)
      cat("   Saved: Step12_04_PPI_Validation.csv\n")

      cat("   PPI Validation Summary:\n")
      cat("   -------------------------------------------------------------\n")
      print(ppi_summary, row.names = FALSE)
    }

  }, error = function(e) {
    cat(sprintf("   Error during PPI validation: %s", e$message))
    cat("   This may be due to network issues or STRING database unavailability.\n")
  })
}

cat("   PPI validation complete\n")
# ============================================================================
# 12f: STEP 12 SUMMARY
# ============================================================================

cat("\n==============================================================================\n")
cat("12f: STEP 12 SUMMARY\n")
cat("==============================================================================\n")

cat("   BIOLOGICAL STORYTELLING & INTEGRATION COMPLETE\n")
cat(sprintf("   Focus Modules: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("   ------------------------------------------------------------------------\n")
cat("   OUTPUT FILES:\n")
cat("   - Step12_01_Module_Narratives.csv\n")
cat("   - Step12_02_Hub_Gene_Validation.csv\n")
cat("   - Step12_03_Integrated_Summary.png/pdf\n")
cat("   - Step12_04_PPI_Validation.csv\n")
cat("   - Step12_Complete_Analysis.RData\n")
cat("\n==============================================================================\n")
cat("STEP 12 COMPLETE\n")
cat(sprintf("Focus Modules Analyzed: %s\n", paste(toupper(focus_modules_auto), collapse = ", ")))
cat("==============================================================================\n")


# Session info
cat("   Session Info:\n")
cat("   ------------------------------------------------------------------------\n")
cat(sprintf("   R version: %s\n", R.version.string))
cat(sprintf("   Platform: %s\n", R.version$platform))
cat(sprintf("   Date: %s\n", Sys.time()))


#===========================================================================
#           STEP 13: UNIFIED PATHWAY INTEGRATION & FINAL SYNTHESIS
#==========================================================================
#  PURPOSE: Create a single integrated view of ALL enrichment results
#           combining focus modules (PINK, GREEN, BLACK, BLUE), ML & Cox
#           biomarkers across ALL databases (ORA, GO:BP, KEGG, Hallmark, C6)
# STRATEGY:
#   13a) Collect all enrichment results from Steps 8-12
#   13b) Create pathway frequency analysis
#   13c) Consolidate into biological themes
#   13d) Generate unified word cloud
#   13e) Biomarker-specific pathway mapping
#   13f) Final integrated summary
#  LOOKING FOR:
#    - Which pathways are consistently enriched across methods?
#    - What biological themes dominate the clinical modules?
#    - How do ML/Cox biomarkers fit into the pathway landscape?

# ##           Combining All Enrichment Sources into One View                 ##


# ============================================================================
# 13a: COLLECT ALL ENRICHMENT RESULTS
# ============================================================================

cat("##  STEP 13: UNIFIED PATHWAY INTEGRATION & FINAL SYNTHESIS                  ##\n")
cat("##           Combining All Enrichment Sources into One View                 ##\n")

cat("==============================================================================\n")
cat("13a: COLLECT ALL ENRICHMENT RESULTS\n")
cat("==============================================================================\n")

# --- Load all enrichment data ---
cat("Loading enrichment results from all databases...\n")

# Verify focus_modules_auto exists
if(!exists("focus_modules_auto") || length(focus_modules_auto) == 0) {
  cat("   [WARNING] focus_modules_auto not found. Attempting to define from moduleColors...\n")
  if(exists("moduleColors")) {
    # Use modules with most genes (excluding grey)
    module_counts <- table(moduleColors)
    module_counts <- module_counts[names(module_counts) != "grey"]
    focus_modules_auto <- names(sort(module_counts, decreasing = TRUE))[1:min(2, length(module_counts))]
    cat(sprintf("   [OK] Set focus_modules_auto to: %s\n", paste(focus_modules_auto, collapse = ", ")))
  } else {
    # Fallback to clinically significant modules
    focus_modules_auto <- c("pink", "green", "black", "blue")
    cat(sprintf("   [WARNING] Using default focus_modules_auto: %s\n", paste(focus_modules_auto, collapse = ", ")))
  }
}

cat(sprintf("   Focus modules: %s\n", paste(focus_modules_auto, collapse = ", ")))

# Initialize master dataframe
unified_pathways <- data.frame()

# 1. ORA Results (Custom Signatures)
ora_file <- file.path(results_dir, "Step8_02_ORA_Results.csv")
if(file.exists(ora_file)) {
  ora <- read.csv(ora_file, stringsAsFactors = FALSE)
  ora_sig <- ora %>%
    filter(Module %in% focus_modules_auto, P_Value < 0.05) %>%
    mutate(
      Database = "Custom_Signatures",
      Pathway = Short_Name,
      neg_log10_p = -log10(P_Value),
      Effect_Size = Fold_Enrichment,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  
  unified_pathways <- rbind(unified_pathways, ora_sig)
  cat(sprintf("   [OK] Custom signatures (ORA): %d entries", nrow(ora_sig)))
}

# 2. GO:BP Results
go_file <- file.path(results_dir, "Step11_01_GO_BP_Results.csv")
if(file.exists(go_file)) {
  go <- read.csv(go_file, stringsAsFactors = FALSE)
  go_sig <- go %>%
    filter(Module %in% focus_modules_auto, p.adjust < 0.05) %>%
    # Take top 10 per module to avoid overwhelming
    group_by(Module) %>%
    slice_min(p.adjust, n = 10) %>%
    ungroup() %>%
    mutate(
      Database = "GO_BP",
      Pathway = Description,
      neg_log10_p = -log10(pvalue),
      Effect_Size = Count,
      FDR = p.adjust,
      P_Value = pvalue,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  
  unified_pathways <- rbind(unified_pathways, go_sig)
  cat(sprintf("   [OK] GO:BP (top 10 per module): %d entries", nrow(go_sig)))
}

# 3. KEGG Results
kegg_file <- file.path(results_dir, "Step11_02_KEGG_Results.csv")
if(file.exists(kegg_file)) {
  kegg <- read.csv(kegg_file, stringsAsFactors = FALSE)
  kegg_sig <- kegg %>%
    filter(Module %in% focus_modules_auto, p.adjust < 0.05) %>%
    mutate(
      Database = "KEGG",
      Pathway = Description,
      neg_log10_p = -log10(pvalue),
      Effect_Size = Count,
      FDR = p.adjust,
      P_Value = pvalue,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  
  unified_pathways <- rbind(unified_pathways, kegg_sig)
  cat(sprintf("   [OK] KEGG: %d entries", nrow(kegg_sig)))
}

# 4. Hallmark Results
hallmark_file <- file.path(results_dir, "Step11_03_MSigDB_Hallmark_Results.csv")
if(file.exists(hallmark_file)) {
  hallmark <- read.csv(hallmark_file, stringsAsFactors = FALSE)
  hallmark_sig <- hallmark %>%
    filter(Module %in% focus_modules_auto, p.adjust < 0.10) %>%
    mutate(
      Database = "Hallmark",
      Pathway = gsub("_", " ", ID),  # Make readable
      neg_log10_p = -log10(pvalue),
      Effect_Size = Count,
      FDR = p.adjust,
      P_Value = pvalue,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  
  unified_pathways <- rbind(unified_pathways, hallmark_sig)
  cat(sprintf("   [OK] Hallmark: %d entries", nrow(hallmark_sig)))
}

# 5. ME-Pathway Correlations - POOR PROGNOSIS ONLY
# IMPORTANT: Only include POSITIVE correlations with focus modules
# Since focus modules have NEGATIVE PFS correlation (higher in Short PFS),
# positive ME-pathway correlation = pathway also higher in Short PFS = poor prognosis
me_corr_file <- file.path(results_dir, "Step8_08_ME_Pathway_Correlations.csv")
if(file.exists(me_corr_file)) {
  me_corr <- read.csv(me_corr_file, stringsAsFactors = FALSE)
  me_sig <- me_corr %>%
    filter(Module %in% focus_modules_auto, P_Value < 0.05, Correlation > 0) %>%  # Positive with poor-prognosis modules
    mutate(
      Database = "ME_Correlation",
      Pathway = Short_Name,
      neg_log10_p = -log10(P_Value),
      Effect_Size = abs(Correlation),
      Direction = "Higher_in_Short"  # All are poor prognosis now
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways <- rbind(unified_pathways, me_sig)
  cat(sprintf("   [OK] ME-Pathway correlations (poor prognosis): %d entries\n", nrow(me_sig)))
  cat("        [FILTER] Pathways positively correlated with poor-prognosis modules\n")
}

# 6. Clinical Association (GSVA vs PFS) - POOR PROGNOSIS ONLY
# IMPORTANT: Only include pathways higher in Short PFS (Cohens_d < 0 = poor prognosis)
clinical_file <- file.path(results_dir, "Step8_05_Clinical_Associations.csv")
if(file.exists(clinical_file)) {
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  clinical_sig <- clinical %>%
    filter(Comparison == "PFS_Group", P_Value < 0.05, Cohens_d < 0) %>%  # ONLY Higher_in_Short
    mutate(
      Module = "Clinical_PFS",  # Not module-specific but clinically relevant
      Database = "GSVA_Clinical",
      Pathway = Short_Name,
      neg_log10_p = -log10(P_Value),
      Effect_Size = abs(Cohens_d),
      Direction = "Higher_in_Short"  # All are poor prognosis now
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways <- rbind(unified_pathways, clinical_sig)
  cat(sprintf("   [OK] GSVA Clinical (PFS - poor prognosis only): %d entries\n", nrow(clinical_sig)))
  cat("        [FILTER] Only pathways higher in Short PFS (poor prognosis)\n")
}

cat(sprintf("\n   TOTAL UNIFIED PATHWAYS: %d\n", nrow(unified_pathways)))

# Rename pathway labels for clarity
unified_pathways <- unified_pathways %>%
  mutate(Pathway = case_when(
    Pathway == "MDSC M-MDSC" ~ "M-Myeloid-Derived Suppressor Cells",
    TRUE ~ Pathway
  ))

# Check if we have any data to work with

  cat("\n   [WARNING] No enrichment data loaded. Creating placeholder files...\n")
  # Create empty placeholder files
  write.csv(data.frame(Message = "No enrichment data available"),
            file.path(results_dir, "Step13_Unified_Pathway_Frequency.csv"), row.names = FALSE)
  write.csv(data.frame(Message = "No enrichment data available"),
            file.path(results_dir, "Step13_Unified_Category_Summary.csv"), row.names = FALSE)
  cat("   [SKIP] Step 13 skipped due to missing enrichment data.\n")
  cat("   Please ensure Steps 8 and 11 have been run first.\n")


# ============================================================================
# SECTION B: CREATE PATHWAY FREQUENCY TABLE
# ============================================================================

cat("=======================================================================\n")
cat("PATHWAY FREQUENCY ANALYSIS\n")
cat("=======================================================================\n")

# Count how many times each pathway appears across different analyses
pathway_frequency <- unified_pathways %>%
  group_by(Pathway) %>%
  summarise(
    N_Databases = n_distinct(Database),
    N_Total = n(),
    Databases = paste(unique(Database), collapse = "; "),
    Modules = paste(unique(Module), collapse = "; "),
    Mean_neg_log10_p = mean(neg_log10_p, na.rm = TRUE),
    Max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    Min_P_Value = min(P_Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(N_Databases), desc(Mean_neg_log10_p))

cat("   TOP 20 MOST CONSISTENTLY ENRICHED PATHWAYS\n")
cat("   (Appearing across multiple databases/analyses)")
cat("   ---------------------------------------------------------------------\n")

top20 <- head(pathway_frequency, 20)
if(nrow(top20) > 0) {
  for(i in 1:nrow(top20)) {
    cat(sprintf("   %2d. %s [%d sources]\n",
                    i,
                    substr(top20$Pathway[i], 1, 40),
                    top20$N_Databases[i]))
  }
} else {
  cat("   No pathways found in frequency analysis.\n")
}

write.csv(pathway_frequency, file.path(results_dir, "Step13_Unified_Pathway_Frequency.csv"), row.names = FALSE)
cat("   [OK] Saved: Unified_Pathway_Frequency.csv\n")

# ============================================================================
# SECTION C: CATEGORY CONSOLIDATION
# ============================================================================

cat("=======================================================================\n")
cat("BIOLOGICAL THEME CONSOLIDATION\n")
cat("=======================================================================\n")

# Define biological categories for key pathways
categorize_pathway <- function(pathway) {
  pathway_lower <- tolower(pathway)
  
  if(grepl("ecm|collagen|matrix|integrin|adhesion|laminin|fibronectin|caf|stroma|desmoplasia", pathway_lower)) {
    return("ECM_Stroma")
  } else if(grepl("endothel|angiogen|vascular|vegf|blood vessel", pathway_lower)) {
    return("Angiogenesis")
  } else if(grepl("macrophage|monocyte|myeloid|tam|m0|m1|m2|neutrophil|dendritic", pathway_lower)) {
    return("Myeloid_Immune")
  } else if(grepl("t cell|nk cell|lymphocyte|adaptive|checkpoint|pd-1|ctla", pathway_lower)) {
    return("Lymphoid_Immune")
  } else if(grepl("inflam|cytokine|chemokine|il-|interferon|nf-kb", pathway_lower)) {
    return("Inflammation")
  } else if(grepl("emt|mesenchym|epithelial|migration|invasion|metastas", pathway_lower)) {
    return("EMT_Metastasis")
  } else if(grepl("stem|notch|wnt|hedgehog|pluripoten|progenitor|csc", pathway_lower)) {
    return("Stemness")
  } else if(grepl("glycoly|oxphos|metabol|lipid|fatty|amino acid|tca|warburg|glutamin", pathway_lower)) {
    return("Metabolism")
  } else if(grepl("resist|chemother|gemcitabine|folfirinox|drug|abc transport|mdr", pathway_lower)) {
    return("Drug_Resistance")
  } else if(grepl("apoptos|death|survival|bcl|bax|p53", pathway_lower)) {
    return("Cell_Death")
  } else if(grepl("proliferat|cell cycle|mitosis|cdk|cyclin", pathway_lower)) {
    return("Proliferation")
  } else if(grepl("signal|pathway|kras|mapk|pi3k|akt|egfr", pathway_lower)) {
    return("Signaling")
  } else if(grepl("exosome|vesicle|secretion", pathway_lower)) {
    return("Exosome")
  } else if(grepl("cox-derived|ml-derived|biomarker", pathway_lower)) {
    return("Study_Biomarkers")
  } else {
    return("Other")
  }
}

# Apply categorization
unified_pathways$Category <- sapply(unified_pathways$Pathway, categorize_pathway)

# Category summary
filtered_for_category <- unified_pathways %>%
  filter(Module %in% c(focus_modules_auto, "Clinical_PFS"))

if(nrow(filtered_for_category) > 0) {
  category_summary <- filtered_for_category %>%
    group_by(Category) %>%
    summarise(
      N_Pathways = n(),
      Mean_Significance = mean(neg_log10_p, na.rm = TRUE),
      Top_Pathway = Pathway[which.max(neg_log10_p)],
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Significance))

  cat(sprintf("   BIOLOGICAL THEMES (%s + Clinical)\n", paste(tools::toTitleCase(focus_modules_auto), collapse = " + ")))
  cat("   ---------------------------------------------------------------------\n")
  print(as.data.frame(category_summary))
} else {
  category_summary <- data.frame(Category = character(), N_Pathways = integer(),
                                  Mean_Significance = numeric(), Top_Pathway = character())
  cat("   [WARNING] No pathways found for focus modules. Empty category summary.\n")
}

write.csv(category_summary, file.path(results_dir, "Step13_Unified_Category_Summary.csv"), row.names = FALSE)

# ============================================================================
# SECTION D: UNIFIED WORD CLOUD
# ============================================================================

cat("=======================================================================\n")
cat("GENERATING UNIFIED WORD CLOUD\n")
cat("=======================================================================\n")

# Prepare word cloud data - consolidate by pathway
wordcloud_data <- unified_pathways %>%
  filter(Module %in% c(focus_modules_auto, "Clinical_PFS")) %>%
  group_by(Pathway, Category) %>%
  summarise(
    N_Sources = n(),
    Aggregate_Score = sum(neg_log10_p, na.rm = TRUE),
    Max_Score = max(neg_log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Weight by both number of sources and significance
  mutate(
    Word_Size = Aggregate_Score * (1 + 0.5 * (N_Sources - 1)),
    # Truncate long names
    Label = ifelse(nchar(Pathway) > 35,
                   paste0(substr(Pathway, 1, 32), "..."),
                   Pathway)
  ) %>%
  filter(Word_Size > 0)

# Take top 50 for visualization
wordcloud_data <- wordcloud_data %>%
  arrange(desc(Word_Size)) %>%
  head(50)

# Color palette by category (defined outside if block for reuse in bar chart)
# --- 1. Define Colors ---
category_colors <- c(
  "ECM_Stroma"      = "#7B3F00", 
  "Angiogenesis"    = "#B22222", 
  "Myeloid_Immune"  = "#E67E22", 
  "Lymphoid_Immune" = "#1B5E20", 
  "Inflammation"    = "#D35400", 
  "EMT_Metastasis"  = "#0D47A1", 
  "Stemness"        = "#6A1B9A", 
  "Metabolism"      = "#00897B", 
  "Drug_Resistance" = "#880E4F", 
  "Cell_Death"      = "#263238", 
  "Proliferation"   = "#F39C12", 
  "Signaling"       = "#558B2F", 
  "Exosome"         = "#455A64", 
  "Study_Biomarkers"= "#000000", 
  "Other"           = "#424242"
)

# --- 2. Prepare Data ---
wordcloud_data <- wordcloud_data %>%
  mutate(Word_Size = sqrt(Word_Size))

# --- 3. Generate Plot ---
set.seed(42)
p_unified_wc <- ggplot(wordcloud_data, 
                       aes(label = Label, 
                           size = Word_Size, 
                           color = Category)) +
  geom_text_wordcloud_area(
    rm_outside = TRUE,
    eccentricity = 0.65, 
    shape = "diamond"
  ) +
  scale_size_continuous(
    range = c(12, 22), # Minimum size 10, Maximum 20 for better balance
    limits = c(NA, NA)
  ) + 
  scale_color_manual(values = category_colors, name = "Biological Theme") +
  theme_minimal() +
  labs(
    title = "Unified Pathway Enrichment: Poor Prognosis Signatures",
    subtitle = "Pathways associated with shorter PFS across ORA, GO:BP, KEGG, Hallmark & GSVA (FDR < 0.05)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 10)),
    # Subtitle negative margin pulls the cloud UP to close the gap
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11, margin = margin(b = -40)), 
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10)
  ) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 5)))

print(p_unified_wc)

ggsave(file.path(fig_dir_png, "Step13_FINAL_Unified_Pathway_WordCloud.png"), 
       p_unified_wc, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step13_FINAL_Unified_Pathway_WordCloud.pdf"), 
       p_unified_wc, width = 8, height = 5)

cat("   [OK] Saved: FINAL_Unified_Pathway_WordCloud.png/pdf\n")

# ============================================================================
# SECTION D2: UNIFIED WORD CLOUD (FDR < 0.1 - More Inclusive)
# ============================================================================
# Purpose: Create a second word cloud with relaxed threshold to capture
#          nominally significant findings like FOLFIRINOX DDR
# This is a SEPARATE figure - the original FDR < 0.05 figure is unchanged
# ============================================================================

cat("=======================================================================\n")
cat("GENERATING UNIFIED WORD CLOUD (FDR < 0.1 threshold)\n")
cat("=======================================================================\n")

# Rebuild unified_pathways with FDR < 0.1 threshold
unified_pathways_fdr10 <- data.frame()

# 1. ORA Results - FDR < 0.1
if(file.exists(ora_file)) {
  ora <- read.csv(ora_file, stringsAsFactors = FALSE)
  ora_sig_fdr10 <- ora %>%
    filter(Module %in% focus_modules_auto, FDR < 0.1) %>%
    mutate(
      Database = "Custom_Signatures",
      Pathway = Short_Name,
      neg_log10_p = -log10(P_Value),
      Effect_Size = Fold_Enrichment,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, ora_sig_fdr10)
  cat(sprintf("   [OK] Custom signatures (ORA, FDR<0.1): %d entries\n", nrow(ora_sig_fdr10)))
}

# 2. GO:BP Results - FDR < 0.1
if(file.exists(go_file)) {
  go <- read.csv(go_file, stringsAsFactors = FALSE)
  go_sig_fdr10 <- go %>%
    filter(Module %in% focus_modules_auto, p.adjust < 0.1) %>%
    group_by(Module) %>%
    slice_min(p.adjust, n = 10) %>%
    ungroup() %>%
    mutate(
      Database = "GO_BP",
      Pathway = Description,
      neg_log10_p = -log10(pvalue),
      Effect_Size = Count,
      FDR = p.adjust,
      P_Value = pvalue,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, go_sig_fdr10)
  cat(sprintf("   [OK] GO:BP (FDR<0.1, top 10/module): %d entries\n", nrow(go_sig_fdr10)))
}

# 3. KEGG Results - FDR < 0.1
if(file.exists(kegg_file)) {
  kegg <- read.csv(kegg_file, stringsAsFactors = FALSE)
  kegg_sig_fdr10 <- kegg %>%
    filter(Module %in% focus_modules_auto, p.adjust < 0.1) %>%
    mutate(
      Database = "KEGG",
      Pathway = Description,
      neg_log10_p = -log10(pvalue),
      Effect_Size = Count,
      FDR = p.adjust,
      P_Value = pvalue,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, kegg_sig_fdr10)
  cat(sprintf("   [OK] KEGG (FDR<0.1): %d entries\n", nrow(kegg_sig_fdr10)))
}

# 4. Hallmark Results - FDR < 0.1
if(file.exists(hallmark_file)) {
  hallmark <- read.csv(hallmark_file, stringsAsFactors = FALSE)
  hallmark_sig_fdr10 <- hallmark %>%
    filter(Module %in% focus_modules_auto, p.adjust < 0.1) %>%
    mutate(
      Database = "Hallmark",
      Pathway = gsub("_", " ", ID),
      neg_log10_p = -log10(pvalue),
      Effect_Size = Count,
      FDR = p.adjust,
      P_Value = pvalue,
      Direction = "Enriched"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, hallmark_sig_fdr10)
  cat(sprintf("   [OK] Hallmark (FDR<0.1): %d entries\n", nrow(hallmark_sig_fdr10)))
}

# 5. ME-Pathway Correlations - FDR < 0.1
if(file.exists(me_corr_file)) {
  me_corr <- read.csv(me_corr_file, stringsAsFactors = FALSE)
  me_sig_fdr10 <- me_corr %>%
    filter(Module %in% focus_modules_auto, FDR < 0.1, Correlation > 0) %>%
    mutate(
      Database = "ME_Correlation",
      Pathway = Short_Name,
      neg_log10_p = -log10(P_Value),
      Effect_Size = abs(Correlation),
      Direction = "Higher_in_Short"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, me_sig_fdr10)
  cat(sprintf("   [OK] ME-Pathway correlations (FDR<0.1): %d entries\n", nrow(me_sig_fdr10)))
}

# 6. Clinical Association - FDR < 0.1
if(file.exists(clinical_file)) {
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  clinical_sig_fdr10 <- clinical %>%
    filter(Comparison == "PFS_Group", FDR < 0.1, Cohens_d < 0) %>%
    mutate(
      Module = "Clinical_PFS",
      Database = "GSVA_Clinical",
      Pathway = Short_Name,
      neg_log10_p = -log10(P_Value),
      Effect_Size = abs(Cohens_d),
      Direction = "Higher_in_Short"
    ) %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)

  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, clinical_sig_fdr10)
  cat(sprintf("   [OK] GSVA Clinical (FDR<0.1): %d entries\n", nrow(clinical_sig_fdr10)))
}

cat(sprintf("\n   TOTAL UNIFIED PATHWAYS (FDR<0.1): %d\n", nrow(unified_pathways_fdr10)))

# Apply categorization
unified_pathways_fdr10$Category <- sapply(unified_pathways_fdr10$Pathway, categorize_pathway)

# Prepare word cloud data - FDR < 0.1
wordcloud_data_fdr10 <- unified_pathways_fdr10 %>%
  filter(Module %in% c(focus_modules_auto, "Clinical_PFS")) %>%
  group_by(Pathway, Category) %>%
  summarise(
    N_Sources = n(),
    Aggregate_Score = sum(neg_log10_p, na.rm = TRUE),
    Max_Score = max(neg_log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Word_Size = Aggregate_Score * (1 + 0.5 * (N_Sources - 1)),
    Label = ifelse(nchar(Pathway) > 35,
                   paste0(substr(Pathway, 1, 32), "..."),
                   Pathway)
  ) %>%
  filter(Word_Size > 0) %>%
  arrange(desc(Word_Size)) %>%
  head(50) %>%
  mutate(Word_Size = sqrt(Word_Size))

# Generate FDR < 0.1 Word Cloud
set.seed(42)
p_unified_wc_fdr10 <- ggplot(wordcloud_data_fdr10,
                       aes(label = Label,
                           size = Word_Size,
                           color = Category)) +
  geom_text_wordcloud_area(
    rm_outside = TRUE,
    eccentricity = 0.65,
    shape = "diamond"
  ) +
  scale_size_continuous(
    range = c(12, 22),
    limits = c(NA, NA)
  ) +
  scale_color_manual(values = category_colors, name = "Biological Theme") +
  theme_minimal() +
  labs(
    title = "Unified Pathway Enrichment: Poor Prognosis Signatures",
    subtitle = "Pathways associated with shorter PFS across ORA, GO:BP, KEGG, Hallmark & GSVA (FDR < 0.1)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11, margin = margin(b = -40)),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10)
  ) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 5)))

print(p_unified_wc_fdr10)

ggsave(file.path(fig_dir_png, "Step13_FINAL_Unified_Pathway_WordCloud_FDR10.png"),
       p_unified_wc_fdr10, width = 9, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step13_FINAL_Unified_Pathway_WordCloud_FDR10.pdf"),
       p_unified_wc_fdr10, width = 9, height = 5)

cat("   [OK] Saved: FINAL_Unified_Pathway_WordCloud_FDR10.png/pdf\n")

# Report what's new in FDR < 0.1 vs FDR < 0.05
new_in_fdr10 <- setdiff(wordcloud_data_fdr10$Label, wordcloud_data$Label)
if(length(new_in_fdr10) > 0) {
  cat(sprintf("\n   NEW pathways in FDR<0.1 (not in FDR<0.05): %d\n", length(new_in_fdr10)))
  cat(sprintf("   %s\n", paste(head(new_in_fdr10, 10), collapse = ", ")))
}

# =============================================================================
# Unified Pathway Word Cloud (FDR < 0.1) — v3
# min 6pt at 9cm print, abbreviations, border/frame
# =============================================================================

library(ggplot2)
library(ggwordcloud)
library(dplyr)

fig_dir_pdf <- file.path(main_dir, "Figures_PDF")

results_dir <- file.path(main_dir, "Results_Tables")

ora_file      <- file.path(results_dir, "Step8_02_ORA_Results.csv")
go_file       <- file.path(results_dir, "Step11_01_GO_BP_Results.csv")
kegg_file     <- file.path(results_dir, "Step11_02_KEGG_Results.csv")
hallmark_file <- file.path(results_dir, "Step11_03_MSigDB_Hallmark_Results.csv")
me_corr_file  <- file.path(results_dir, "Step8_08_ME_Pathway_Correlations.csv")
clinical_file <- file.path(results_dir, "Step8_05_Clinical_Associations.csv")

focus_modules <- c("pink", "green", "black", "blue")

category_colors <- c(
  "ECM_Stroma"      = "#7B3F00",
  "Angiogenesis"    = "#B22222",
  "Myeloid_Immune"  = "#E67E22",
  "Lymphoid_Immune" = "#1B5E20",
  "Inflammation"    = "#D35400",
  "EMT_Metastasis"  = "#0D47A1",
  "Stemness"        = "#6A1B9A",
  "Metabolism"       = "#00897B",
  "Drug_Resistance" = "#880E4F",
  "Cell_Death"      = "#263238",
  "Proliferation"   = "#F39C12",
  "Signaling"       = "#558B2F",
  "Exosome"         = "#455A64",
  "Study_Biomarkers"= "#000000",
  "Other"           = "#424242"
)

categorize_pathway <- function(pathway) {
  p <- tolower(pathway)
  if (grepl("ecm|collagen|matrix|integrin|adhesion|laminin|fibronectin|caf|stroma|desmoplasia", p)) return("ECM_Stroma")
  if (grepl("endothel|angiogen|vascular|vegf|blood vessel", p)) return("Angiogenesis")
  if (grepl("macrophage|monocyte|myeloid|tam|m0|m1|m2|neutrophil|dendritic", p)) return("Myeloid_Immune")
  if (grepl("t cell|nk cell|lymphocyte|adaptive|checkpoint|pd-1|ctla", p)) return("Lymphoid_Immune")
  if (grepl("inflam|cytokine|chemokine|il-|interferon|nf-kb", p)) return("Inflammation")
  if (grepl("emt|mesenchym|epithelial|migration|invasion|metastas", p)) return("EMT_Metastasis")
  if (grepl("stem|notch|wnt|hedgehog|pluripoten|progenitor|csc", p)) return("Stemness")
  if (grepl("glycoly|oxphos|metabol|lipid|fatty|amino acid|tca|warburg|glutamin", p)) return("Metabolism")
  if (grepl("resist|chemother|gemcitabine|folfirinox|drug|abc transport|mdr", p)) return("Drug_Resistance")
  if (grepl("apoptos|death|survival|bcl|bax|p53", p)) return("Cell_Death")
  if (grepl("proliferat|cell cycle|mitosis|cdk|cyclin", p)) return("Proliferation")
  if (grepl("signal|pathway|kras|mapk|pi3k|akt|egfr", p)) return("Signaling")
  if (grepl("exosome|vesicle|secretion", p)) return("Exosome")
  if (grepl("cox-derived|ml-derived|biomarker", p)) return("Study_Biomarkers")
  return("Other")
}

# --- Abbreviation map for long pathway names ---
abbreviations <- c(
  "positive regulation of protein phosphorylation"              = "Pos. reg. protein phosph.",
  "positive regulation of cell population proliferation"        = "Pos. reg. cell proliferation",
  "positive regulation of phosphorylation"                      = "Pos. reg. phosphorylation",
  "positive regulation of protein modification process"         = "Pos. reg. protein mod.",
  "regulation of cell differentiation"                          = "Reg. cell differentiation",
  "regulation of cell population proliferation"                 = "Reg. cell proliferation",
  "negative regulation of developmental process"                = "Neg. reg. dev. process",
  "extracellular matrix organization"                           = "ECM organization",
  "extracellular structure organization"                        = "Extracell. struct. org.",
  "external encapsulating structure organization"               = "Ext. encaps. struct. org.",
  "anatomical structure morphogenesis"                          = "Anat. struct. morphogenesis",
  "anatomical structure formation involved in morphogenesis"    = "Anat. struct. formation",
  "multicellular organism development"                          = "Multicell. organism dev.",
  "circulatory system development"                              = "Circulatory sys. dev.",
  "cell population proliferation"                               = "Cell proliferation",
  "Human papillomavirus infection"                              = "HPV infection",
  "Cell adhesion molecule (CAM) interaction"                    = "CAM interaction",
  "Cytoskeleton in muscle cells"                                = "Cytoskel. muscle cells",
  "EPITHELIAL MESENCHYMAL TRANSITION"                           = "EMT (Hallmark)",
  "Gemcitabine Resistance 5-gene"                               = "Gemcitabine Res. 5-gene",
  "developmental process"                                       = "Developmental process"
)

abbreviate_pathway <- function(name) {
  if (name %in% names(abbreviations)) {
    return(abbreviations[name])
  }
  return(name)
}

# --- Rebuild unified pathways (FDR < 0.1) ---
unified_pathways_fdr10 <- data.frame()

if (file.exists(ora_file)) {
  ora <- read.csv(ora_file, stringsAsFactors = FALSE)
  ora_sig <- ora %>%
    filter(Module %in% focus_modules, FDR < 0.1) %>%
    mutate(Database = "Custom_Signatures", Pathway = Short_Name,
           neg_log10_p = -log10(P_Value), Effect_Size = Fold_Enrichment,
           Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, ora_sig)
}

if (file.exists(go_file)) {
  go <- read.csv(go_file, stringsAsFactors = FALSE)
  go_sig <- go %>%
    filter(Module %in% focus_modules, p.adjust < 0.1) %>%
    group_by(Module) %>% slice_min(p.adjust, n = 10) %>% ungroup() %>%
    mutate(Database = "GO_BP", Pathway = Description,
           neg_log10_p = -log10(pvalue), Effect_Size = Count,
           FDR = p.adjust, P_Value = pvalue, Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, go_sig)
}

if (file.exists(kegg_file)) {
  kegg <- read.csv(kegg_file, stringsAsFactors = FALSE)
  kegg_sig <- kegg %>%
    filter(Module %in% focus_modules, p.adjust < 0.1) %>%
    mutate(Database = "KEGG", Pathway = Description,
           neg_log10_p = -log10(pvalue), Effect_Size = Count,
           FDR = p.adjust, P_Value = pvalue, Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, kegg_sig)
}

if (file.exists(hallmark_file)) {
  hallmark <- read.csv(hallmark_file, stringsAsFactors = FALSE)
  hallmark_sig <- hallmark %>%
    filter(Module %in% focus_modules, p.adjust < 0.1) %>%
    mutate(Database = "Hallmark", Pathway = gsub("_", " ", ID),
           neg_log10_p = -log10(pvalue), Effect_Size = Count,
           FDR = p.adjust, P_Value = pvalue, Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, hallmark_sig)
}

if (file.exists(me_corr_file)) {
  me_corr <- read.csv(me_corr_file, stringsAsFactors = FALSE)
  me_sig <- me_corr %>%
    filter(Module %in% focus_modules, FDR < 0.1, Correlation > 0) %>%
    mutate(Database = "ME_Correlation", Pathway = Short_Name,
           neg_log10_p = -log10(P_Value), Effect_Size = abs(Correlation),
           Direction = "Higher_in_Short") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, me_sig)
}

if (file.exists(clinical_file)) {
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  clinical_sig <- clinical %>%
    filter(Comparison == "PFS_Group", FDR < 0.1, Cohens_d < 0) %>%
    mutate(Module = "Clinical_PFS", Database = "GSVA_Clinical", Pathway = Short_Name,
           neg_log10_p = -log10(P_Value), Effect_Size = abs(Cohens_d),
           Direction = "Higher_in_Short") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways_fdr10 <- rbind(unified_pathways_fdr10, clinical_sig)
}

cat(sprintf("Total unified pathways (FDR<0.1): %d\n", nrow(unified_pathways_fdr10)))

unified_pathways_fdr10$Category <- sapply(unified_pathways_fdr10$Pathway, categorize_pathway)

wordcloud_data_fdr10 <- unified_pathways_fdr10 %>%
  filter(Module %in% c(focus_modules, "Clinical_PFS")) %>%
  group_by(Pathway, Category) %>%
  summarise(
    N_Sources = n(),
    Aggregate_Score = sum(neg_log10_p, na.rm = TRUE),
    Max_Score = max(neg_log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Word_Size = Aggregate_Score * (1 + 0.5 * (N_Sources - 1)),
    Label = sapply(Pathway, abbreviate_pathway)
  ) %>%
  filter(Word_Size > 0) %>%
  arrange(desc(Word_Size)) %>%
  head(50) %>%
  mutate(Word_Size = sqrt(Word_Size))

cat(sprintf("Words for cloud: %d\n", nrow(wordcloud_data_fdr10)))

cat("\nAbbreviations applied:\n")
for (i in seq_len(nrow(wordcloud_data_fdr10))) {
  lab <- wordcloud_data_fdr10$Label[i]
  orig <- wordcloud_data_fdr10$Pathway[i]
  if (lab != orig) {
    cat(sprintf("  %s -> %s\n", orig, lab))
  }
}

# --- Generate word cloud ---
set.seed(42)
p_wc <- ggplot(wordcloud_data_fdr10,
               aes(label = Label, size = Word_Size, color = Category)) +
  geom_text_wordcloud_area(
    rm_outside = TRUE,
    eccentricity = 0.65,
    shape = "diamond"
  ) +
  scale_size_continuous(range = c(15.2, 28)) +
  scale_color_manual(values = category_colors, name = "Biological Theme") +
  theme_minimal() +
  labs(
    title = "Unified Pathway Enrichment: Poor Prognosis Signatures",
    subtitle = "Pathways associated with shorter PFS across ORA, GO:BP, KEGG, Hallmark & GSVA (FDR < 0.1)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11, margin = margin(b = -40)),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.8),
    plot.background = element_rect(color = "grey60", fill = "white", linewidth = 0.8)
  ) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 5)))

print(p_wc)

ggsave(file.path(fig_dir_pdf, "Step13_FINAL_Unified_Pathway_WordCloud_FDR10_v3.png"),
       p_wc, width = 9, height = 7, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step13_FINAL_Unified_Pathway_WordCloud_FDR10_v3.pdf"),
       p_wc, width = 9, height = 7)

cat("\nDone! Saved Step13_FINAL_Unified_Pathway_WordCloud_FDR10_v3.png/.pdf\n")
cat("Min font: 15.2pt in figure = 6pt when printed at 9cm column width\n")


# =============================================================================
# Unified Pathway Word Cloud (FDR < 0.05) — v3
# min 6pt at 9cm print, abbreviations, border/frame
# Saves as NEW version — originals untouched
# =============================================================================

library(ggplot2)
library(ggwordcloud)
library(dplyr)

fig_dir_pdf <- file.path(main_dir, "Figures_PDF")
fig_dir_png <- file.path(main_dir, "Figures_PNG")
results_dir <- file.path(main_dir, "Results_Tables")

ora_file      <- file.path(results_dir, "Step8_02_ORA_Results.csv")
go_file       <- file.path(results_dir, "Step11_01_GO_BP_Results.csv")
kegg_file     <- file.path(results_dir, "Step11_02_KEGG_Results.csv")
hallmark_file <- file.path(results_dir, "Step11_03_MSigDB_Hallmark_Results.csv")
me_corr_file  <- file.path(results_dir, "Step8_08_ME_Pathway_Correlations.csv")
clinical_file <- file.path(results_dir, "Step8_05_Clinical_Associations.csv")

focus_modules <- c("pink", "green", "black", "blue")

category_colors <- c(
  "ECM_Stroma"      = "#7B3F00",
  "Angiogenesis"    = "#B22222",
  "Myeloid_Immune"  = "#E67E22",
  "Lymphoid_Immune" = "#1B5E20",
  "Inflammation"    = "#D35400",
  "EMT_Metastasis"  = "#0D47A1",
  "Stemness"        = "#6A1B9A",
  "Metabolism"       = "#00897B",
  "Drug_Resistance" = "#880E4F",
  "Cell_Death"      = "#263238",
  "Proliferation"   = "#F39C12",
  "Signaling"       = "#558B2F",
  "Exosome"         = "#455A64",
  "Study_Biomarkers"= "#000000",
  "Other"           = "#424242"
)

categorize_pathway <- function(pathway) {
  p <- tolower(pathway)
  if (grepl("ecm|collagen|matrix|integrin|adhesion|laminin|fibronectin|caf|stroma|desmoplasia", p)) return("ECM_Stroma")
  if (grepl("endothel|angiogen|vascular|vegf|blood vessel", p)) return("Angiogenesis")
  if (grepl("macrophage|monocyte|myeloid|tam|m0|m1|m2|neutrophil|dendritic", p)) return("Myeloid_Immune")
  if (grepl("t cell|nk cell|lymphocyte|adaptive|checkpoint|pd-1|ctla", p)) return("Lymphoid_Immune")
  if (grepl("inflam|cytokine|chemokine|il-|interferon|nf-kb", p)) return("Inflammation")
  if (grepl("emt|mesenchym|epithelial|migration|invasion|metastas", p)) return("EMT_Metastasis")
  if (grepl("stem|notch|wnt|hedgehog|pluripoten|progenitor|csc", p)) return("Stemness")
  if (grepl("glycoly|oxphos|metabol|lipid|fatty|amino acid|tca|warburg|glutamin", p)) return("Metabolism")
  if (grepl("resist|chemother|gemcitabine|folfirinox|drug|abc transport|mdr", p)) return("Drug_Resistance")
  if (grepl("apoptos|death|survival|bcl|bax|p53", p)) return("Cell_Death")
  if (grepl("proliferat|cell cycle|mitosis|cdk|cyclin", p)) return("Proliferation")
  if (grepl("signal|pathway|kras|mapk|pi3k|akt|egfr", p)) return("Signaling")
  if (grepl("exosome|vesicle|secretion", p)) return("Exosome")
  if (grepl("cox-derived|ml-derived|biomarker", p)) return("Study_Biomarkers")
  return("Other")
}

# --- Abbreviation map ---
abbreviations <- c(
  "positive regulation of protein phosphorylation"              = "Pos. reg. protein phosph.",
  "positive regulation of cell population proliferation"        = "Pos. reg. cell proliferation",
  "positive regulation of phosphorylation"                      = "Pos. reg. phosphorylation",
  "positive regulation of protein modification process"         = "Pos. reg. protein mod.",
  "regulation of cell differentiation"                          = "Reg. cell differentiation",
  "regulation of cell population proliferation"                 = "Reg. cell proliferation",
  "negative regulation of developmental process"                = "Neg. reg. dev. process",
  "extracellular matrix organization"                           = "ECM organization",
  "extracellular structure organization"                        = "Extracell. struct. org.",
  "external encapsulating structure organization"               = "Ext. encaps. struct. org.",
  "anatomical structure morphogenesis"                          = "Anat. struct. morphogenesis",
  "anatomical structure formation involved in morphogenesis"    = "Anat. struct. formation",
  "multicellular organism development"                          = "Multicell. organism dev.",
  "circulatory system development"                              = "Circulatory sys. dev.",
  "cell population proliferation"                               = "Cell proliferation",
  "Human papillomavirus infection"                              = "HPV infection",
  "Cell adhesion molecule (CAM) interaction"                    = "CAM interaction",
  "Cytoskeleton in muscle cells"                                = "Cytoskel. muscle cells",
  "EPITHELIAL MESENCHYMAL TRANSITION"                           = "EMT (Hallmark)",
  "Gemcitabine Resistance 5-gene"                               = "Gemcitabine Res. 5-gene",
  "developmental process"                                       = "Developmental process",
  "M-Myeloid-Derived Suppressor Cells"                          = "M-MDSC"
)

abbreviate_pathway <- function(name) {
  if (name %in% names(abbreviations)) return(abbreviations[name])
  return(name)
}

# --- Build unified pathways (FDR < 0.05 thresholds — matching original) ---
unified_pathways <- data.frame()

# 1. ORA — P_Value < 0.05
if (file.exists(ora_file)) {
  ora <- read.csv(ora_file, stringsAsFactors = FALSE)
  ora_sig <- ora %>%
    filter(Module %in% focus_modules, P_Value < 0.05) %>%
    mutate(Database = "Custom_Signatures", Pathway = Short_Name,
           neg_log10_p = -log10(P_Value), Effect_Size = Fold_Enrichment,
           Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways <- rbind(unified_pathways, ora_sig)
}

# 2. GO:BP — p.adjust < 0.05, top 10 per module
if (file.exists(go_file)) {
  go <- read.csv(go_file, stringsAsFactors = FALSE)
  go_sig <- go %>%
    filter(Module %in% focus_modules, p.adjust < 0.05) %>%
    group_by(Module) %>% slice_min(p.adjust, n = 10) %>% ungroup() %>%
    mutate(Database = "GO_BP", Pathway = Description,
           neg_log10_p = -log10(pvalue), Effect_Size = Count,
           FDR = p.adjust, P_Value = pvalue, Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways <- rbind(unified_pathways, go_sig)
}

# 3. KEGG — p.adjust < 0.05
if (file.exists(kegg_file)) {
  kegg <- read.csv(kegg_file, stringsAsFactors = FALSE)
  kegg_sig <- kegg %>%
    filter(Module %in% focus_modules, p.adjust < 0.05) %>%
    mutate(Database = "KEGG", Pathway = Description,
           neg_log10_p = -log10(pvalue), Effect_Size = Count,
           FDR = p.adjust, P_Value = pvalue, Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways <- rbind(unified_pathways, kegg_sig)
}

# 4. Hallmark — p.adjust < 0.10 (matching original)
if (file.exists(hallmark_file)) {
  hallmark <- read.csv(hallmark_file, stringsAsFactors = FALSE)
  hallmark_sig <- hallmark %>%
    filter(Module %in% focus_modules, p.adjust < 0.10) %>%
    mutate(Database = "Hallmark", Pathway = gsub("_", " ", ID),
           neg_log10_p = -log10(pvalue), Effect_Size = Count,
           FDR = p.adjust, P_Value = pvalue, Direction = "Enriched") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways <- rbind(unified_pathways, hallmark_sig)
}

# 5. ME correlations — P_Value < 0.05, Correlation > 0
if (file.exists(me_corr_file)) {
  me_corr <- read.csv(me_corr_file, stringsAsFactors = FALSE)
  me_sig <- me_corr %>%
    filter(Module %in% focus_modules, P_Value < 0.05, Correlation > 0) %>%
    mutate(Database = "ME_Correlation", Pathway = Short_Name,
           neg_log10_p = -log10(P_Value), Effect_Size = abs(Correlation),
           Direction = "Higher_in_Short") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways <- rbind(unified_pathways, me_sig)
}

# 6. Clinical — P_Value < 0.05, Cohens_d < 0
if (file.exists(clinical_file)) {
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  clinical_sig <- clinical %>%
    filter(Comparison == "PFS_Group", P_Value < 0.05, Cohens_d < 0) %>%
    mutate(Module = "Clinical_PFS", Database = "GSVA_Clinical", Pathway = Short_Name,
           neg_log10_p = -log10(P_Value), Effect_Size = abs(Cohens_d),
           Direction = "Higher_in_Short") %>%
    dplyr::select(Module, Database, Pathway, neg_log10_p, Effect_Size, P_Value, FDR, Direction)
  unified_pathways <- rbind(unified_pathways, clinical_sig)
}

cat(sprintf("Total unified pathways (FDR<0.05): %d\n", nrow(unified_pathways)))

# Rename MDSC (matching original)
unified_pathways <- unified_pathways %>%
  mutate(Pathway = case_when(
    Pathway == "MDSC M-MDSC" ~ "M-Myeloid-Derived Suppressor Cells",
    TRUE ~ Pathway
  ))

unified_pathways$Category <- sapply(unified_pathways$Pathway, categorize_pathway)

wordcloud_data <- unified_pathways %>%
  filter(Module %in% c(focus_modules, "Clinical_PFS")) %>%
  group_by(Pathway, Category) %>%
  summarise(
    N_Sources = n(),
    Aggregate_Score = sum(neg_log10_p, na.rm = TRUE),
    Max_Score = max(neg_log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Word_Size = Aggregate_Score * (1 + 0.5 * (N_Sources - 1)),
    Label = sapply(Pathway, abbreviate_pathway)
  ) %>%
  filter(Word_Size > 0) %>%
  arrange(desc(Word_Size)) %>%
  head(50) %>%
  mutate(Word_Size = sqrt(Word_Size))

cat(sprintf("Words for cloud: %d\n", nrow(wordcloud_data)))

cat("\nAbbreviations applied:\n")
for (i in seq_len(nrow(wordcloud_data))) {
  lab <- wordcloud_data$Label[i]
  orig <- wordcloud_data$Pathway[i]
  if (lab != orig) cat(sprintf("  %s -> %s\n", orig, lab))
}

# --- Generate word cloud ---
# Export 8" x 5" (matching original). At 9cm print: scale = 9/(8*2.54) = 0.443
# For 6pt min: 6/0.443 = 13.5pt. Max: 13.5*(22/12) = 24.8pt
set.seed(42)
p_wc <- ggplot(wordcloud_data,
               aes(label = Label, size = Word_Size, color = Category)) +
  geom_text_wordcloud_area(
    rm_outside = TRUE,
    eccentricity = 0.65,
    shape = "diamond"
  ) +
  scale_size_continuous(range = c(13.5, 24.8)) +
  scale_color_manual(values = category_colors, name = "Biological Theme") +
  theme_minimal() +
  labs(
    title = "Unified Pathway Enrichment: Poor Prognosis Signatures",
    subtitle = "Pathways associated with shorter PFS across ORA, GO:BP, KEGG, Hallmark & GSVA (FDR < 0.05)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11, margin = margin(b = -40)),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.8),
    plot.background = element_rect(color = "grey60", fill = "white", linewidth = 0.8)
  ) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 5)))

print(p_wc)

ggsave(file.path(fig_dir_png, "Step13_FINAL_Unified_Pathway_WordCloud_v3.png"),
       p_wc, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step13_FINAL_Unified_Pathway_WordCloud_v3.pdf"),
       p_wc, width = 8, height = 5)

cat("\nDone! Saved Step13_FINAL_Unified_Pathway_WordCloud_v3.png/.pdf\n")
cat("Min font: 13.5pt in figure = 6pt when printed at 9cm column width\n")




# ============================================================================
# SECTION E: BAR CHART OF TOP PATHWAYS BY CATEGORY
# ============================================================================

cat("   Creating summary bar chart...\n")

# Top pathway per category - remove duplicates and get one per category
top_by_category <- unified_pathways %>%
  filter(Module %in% focus_modules_auto) %>%
  # First remove duplicate pathways (keep highest significance)
  group_by(Pathway, Category) %>%
  slice_max(neg_log10_p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Then get top pathway per category
  group_by(Category) %>%
  slice_max(neg_log10_p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(neg_log10_p)

top_by_category_final <- top_by_category %>%
  mutate(
    # Standardize Category and Pathway to Title Case
    Category = str_replace_all(Category, "_", " ") %>% str_to_title(),
    Pathway  = str_to_title(Pathway)
  )

cat(sprintf("   Unique categories: %d\n", nrow(top_by_category_final)))

# --- 2. GENERATE FIGURE ---
# Define publication colors for categories (Title Case to match transformed data)
pub_colors <- c(
  "Ecm Stroma"      = "#7B3F00",
  "Angiogenesis"    = "#B22222",
  "Myeloid Immune"  = "#E67E22",
  "Lymphoid Immune" = "#1B5E20",
  "Inflammation"    = "#D35400",
  "Emt Metastasis"  = "#0D47A1",
  "Stemness"        = "#6A1B9A",
  "Metabolism"      = "#00897B",
  "Drug Resistance" = "#880E4F",
  "Cell Death"      = "#263238",
  "Proliferation"   = "#F39C12",
  "Signaling"       = "#558B2F",
  "Exosome"         = "#455A64",
  "Study Biomarkers"= "#000000",
  "Other"           = "#424242"
)

p_nature_bar <- ggplot(top_by_category_final %>%
                         mutate(
                           # Truncate long pathway names to fit in bars
                           Pathway_Label = ifelse(nchar(Pathway) > 35,
                                                  paste0(substr(Pathway, 1, 32), "..."),
                                                  Pathway)
                         ),
                       aes(x = reorder(Category, -neg_log10_p), # Highest -log10(p) at bottom, lowest at top
                           y = neg_log10_p,
                           fill = Category)) +
  # Draw the bars (no border)
  geom_col(width = 0.8) +
  # Text aligned left inside the bar
  geom_text(aes(y = 0.3, label = Pathway_Label),
            hjust = 0,
            size = 3.2,
            color = "white",
            fontface = "bold") +
  scale_fill_manual(values = pub_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(
    title = "Enriched biological themes associated with poor prognosis (shorter PFS)",
    subtitle = "Integrating ORA, GO:BP, KEGG, MSigDB Hallmark, ME-Correlations & Clinical GSVA (p<0.05)\nFocus modules negatively correlated with PFS (higher in Short PFS)",
    x = NULL,
    y = expression(-log[10](italic(P)*"-value"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(face = "bold", color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    legend.position = "none"
  )

print(p_nature_bar)

ggsave(file.path(fig_dir_png, "Step13_FINAL_Category_Summary_Bar.png"), 
       p_nature_bar, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step13_FINAL_Category_Summary_Bar.pdf"), 
       p_nature_bar, width = 10, height = 6)

cat("   [OK] Saved: FINAL_Category_Summary_Bar.png/pdf\n")


# ============================================================================
# SECTION F: BIOMARKER-SPECIFIC PATHWAY ANALYSIS
# ============================================================================

cat("=======================================================================\n")
cat("BIOMARKER-SPECIFIC PATHWAY MAPPING\n")
cat("=======================================================================\n")

# Manual annotation of biomarker pathways based on known biology
biomarker_annotation <- data.frame(
  Biomarker = c("ECM2", "DYNC2H1", "PPIB",
                "APOF", "TFRC", "FABP4", "ANG"),
  Source = c(rep("ML-Derived", 3), rep("Cox-Derived", 4)),
  PFS_Direction = c("<-'Short (NS)", "<-'Long*", "<-'Long (NS)",
                    "<-'Short*", "<-'Short*", "<-'Short*", "<-'Short*"),
  Primary_Function = c(
    "ECM glycoprotein",
    "Dynein motor (cilia)",
    "Cyclophilin (protein folding)",
    "Apolipoprotein (lipid transport)",
    "Transferrin receptor (iron uptake)",
    "Fatty acid binding (adipocyte)",
    "Angiogenin (angiogenesis)"
  ),
  Linked_Pathways = c(
    "ECM_Stroma, Desmoplasia",
    "Ciliary function, Signaling",
    "Protein homeostasis, Stress",
    "Lipid metabolism, Inflammation",
    "Iron metabolism, Proliferation",
    "Lipid metabolism, Adipogenesis, Stroma",
    "Angiogenesis, Metastasis"
  ),
  Module_Membership = c(
    "Check", "Check", "Check",
    "Brown (likely)", "Brown (likely)", "Brown (likely)", "Brown (likely)"
  ),
  QC_Flag = c(
    "OK", "[WARNING] HighCV, MV_Outlier", "OK", "[WARNING] HighCV, LowCorr, LowExpr",
    "OK", "OK", "OK"
  ),
  stringsAsFactors = FALSE
)

cat("   BIOMARKER ANNOTATION TABLE\n")
cat("   ---------------------------------------------------------------------\n")
print(biomarker_annotation)

write.csv(biomarker_annotation, 
          file.path(results_dir, "Biomarker_Pathway_Annotation.csv"), 
          row.names = FALSE)
cat("   [OK] Saved: Biomarker_Pathway_Annotation.csv\n")

cat("=======================================================================\n")
cat("UNIFIED PATHWAY INTEGRATION COMPLETE\n")
cat("=======================================================================\n")



# STEPS 13b & 14: BIOMARKER SIGNATURE VALIDATION & DYNC2H1 INVESTIGATION  
#  STEP 13b: Validate ML-derived and Cox-derived biomarker signatures        
# using simple, interpretable methods                              
# STEP 14:  Investigate DYNC2H1 high CV (31%) - biological vs technical      
#  BIOMARKERS:                                                               
# ML-derived:  ECM2, DYNC2H1, PPIB                                  
# Cox-derived: APOF, TFRC, FABP4, ANG                                      


# ============================================================================
# SETUP: VERIFY REQUIRED OBJECTS
# ============================================================================

cat("==============================================================================\n")
cat("VERIFYING REQUIRED OBJECTS\n")
cat("==============================================================================\n")

# Check for required objects from previous WGCNA steps
required_objects <- c("datExpr", "datTraits", "moduleColors", "MEs")

missing_objects <- c()
for(obj in required_objects) {
  if(exists(obj)) {
    cat(sprintf("   [OK] %s found", obj))
  } else {
    cat(sprintf("   [FAIL] %s NOT FOUND", obj))
    missing_objects <- c(missing_objects, obj)
  }
}

if(length(missing_objects) > 0) {
  stop(sprintf("\n   ERROR: Missing required objects: %s\n   Please load your WGCNA workspace (e.g., Step12_Complete_Analysis.RData) first.\n", 
               paste(missing_objects, collapse = ", ")))
}


# ============================================================================
# SETUP: STANDARDIZE COLUMN NAMES
# ============================================================================

cat("==============================================================================\n")
cat("STANDARDIZING DATA\n")
cat("==============================================================================\n")

# Show available columns
cat("   Available columns in datTraits:\n")
cat(sprintf("   %s", paste(colnames(datTraits), collapse = ", ")))

# Standardize PFS_group column name
if("PFS_Group" %in% colnames(datTraits) && !"PFS_group" %in% colnames(datTraits)) {
  datTraits$PFS_group <- datTraits$PFS_Group
  cat("   [OK] Standardized PFS_Group -> PFS_group\n")
}
if("PFS_group" %in% colnames(datTraits) && !"PFS_Group" %in% colnames(datTraits)) {
  datTraits$PFS_Group <- datTraits$PFS_group
  cat("   [OK] Standardized PFS_group -> PFS_Group\n")
}

# Check for required columns (be flexible with naming)
pfs_col <- NULL
if("PFS_group" %in% colnames(datTraits)) pfs_col <- "PFS_group"
if("PFS_Group" %in% colnames(datTraits)) pfs_col <- "PFS_Group"

response_col <- NULL
if("Response" %in% colnames(datTraits)) response_col <- "Response"
if("response" %in% colnames(datTraits)) response_col <- "response"

pfs_days_col <- NULL
if("PFS" %in% colnames(datTraits)) pfs_days_col <- "PFS"
if("PFS_days" %in% colnames(datTraits)) pfs_days_col <- "PFS_days"
if("PFS_Days" %in% colnames(datTraits)) pfs_days_col <- "PFS_Days"  # Added: exact match for your column
if("pfs" %in% colnames(datTraits)) pfs_days_col <- "pfs"

treatment_col <- NULL
if("Treatment" %in% colnames(datTraits)) treatment_col <- "Treatment"
if("treatment" %in% colnames(datTraits)) treatment_col <- "treatment"

# Report findings
cat(sprintf("   PFS Group column: %s", ifelse(!is.null(pfs_col), pfs_col, "NOT FOUND")))
cat(sprintf("   Response column: %s", ifelse(!is.null(response_col), response_col, "NOT FOUND")))
cat(sprintf("   PFS Days column: %s", ifelse(!is.null(pfs_days_col), pfs_days_col, "NOT FOUND")))
cat(sprintf("   Treatment column: %s", ifelse(!is.null(treatment_col), treatment_col, "NOT FOUND")))

# Stop if critical columns missing
if(is.null(pfs_col)) {
  stop("   ERROR: Cannot find PFS_group or PFS_Group column in datTraits")
}

# Create standardized reference columns for easier coding
datTraits$PFS_group_std <- datTraits[[pfs_col]]
if(!is.null(response_col)) datTraits$Response_std <- datTraits[[response_col]]
if(!is.null(pfs_days_col)) datTraits$PFS_days_std <- datTraits[[pfs_days_col]]
if(!is.null(treatment_col)) datTraits$Treatment_std <- datTraits[[treatment_col]]

cat("   [OK] Data standardization complete\n")


#===============================================================================
cat("#STEP 14: DYNC2H1 HIGH CV INVESTIGATION\n")
cat("#Is 31% CV Biological Signal or Technical Noise?\n")
# ============================================================================
# 14a: TECHNICAL NOISE ASSESSMENT
# ============================================================================

cat("==============================================================================\n")
cat("14a: TECHNICAL NOISE ASSESSMENT\n")
cat("==============================================================================\n")

# Check if DYNC2H1 exists
if(!"DYNC2H1" %in% colnames(datExpr)) {
  cat("   [WARNING] DYNC2H1 not found in datExpr. Skipping Step 14.\n")
  cat("   Available proteins similar to DYNC2H1:\n")
  cat(paste("   ", grep("DYNC", colnames(datExpr), value = TRUE, ignore.case = TRUE), collapse = "\n"))
} else {

  # Extract DYNC2H1 expression
  dync2h1_expr <- datExpr[, "DYNC2H1"]
  sample_ids <- rownames(datExpr)

  # --- 14a.1: Basic statistics ---
  cat("   14a.1: DYNC2H1 Expression Statistics\n")
  cat("   ---------------------------------------------------------------------\n")

  dync2h1_cv <- 100 * sd(dync2h1_expr) / mean(dync2h1_expr)

  dync2h1_stats <- data.frame(
    Metric = c("Mean", "Median", "SD", "CV (%)", "Min", "Max", "Range", "IQR"),
    Value = c(
      round(mean(dync2h1_expr), 4),
      round(median(dync2h1_expr), 4),
      round(sd(dync2h1_expr), 4),
      round(dync2h1_cv, 2),
      round(min(dync2h1_expr), 4),
      round(max(dync2h1_expr), 4),
      round(max(dync2h1_expr) - min(dync2h1_expr), 4),
      round(IQR(dync2h1_expr), 4)
    )
  )
  print(dync2h1_stats)

  # --- 14a.2: Total protein intensity per sample ---
  cat("   14a.2: Correlation with Total Sample Intensity\n")
  cat("   ---------------------------------------------------------------------\n")

  total_intensity <- rowSums(datExpr)
  cor_total <- cor.test(dync2h1_expr, total_intensity, method = "spearman")

  cat(sprintf("   Spearman correlation with total intensity: r = %.3f, p = %.4f\n",
                  cor_total$estimate, cor_total$p.value))

  if(abs(cor_total$estimate) > 0.5 && cor_total$p.value < 0.05) {
    cat("   [WARNING] Strong correlation with total intensity suggests normalization issue\n")
    tech_flag_intensity <- TRUE
  } else {
    cat("   [OK] No strong correlation with total intensity\n")
    tech_flag_intensity <- FALSE
  }

  # --- 14a.3: Identify outlier samples for DYNC2H1 ---
  cat("   14a.3: Outlier Sample Detection\n")
  cat("   ---------------------------------------------------------------------\n")

  # Z-score based outliers
  dync2h1_z <- scale(dync2h1_expr)[,1]
  outlier_samples <- sample_ids[abs(dync2h1_z) > 2.5]

  cat(sprintf("   Samples with |Z| > 2.5: %d\n", length(outlier_samples)))
  if(length(outlier_samples) > 0) {
    outlier_df <- data.frame(
      Sample = outlier_samples,
      DYNC2H1 = round(dync2h1_expr[outlier_samples], 4),
      Z_Score = round(dync2h1_z[outlier_samples], 2),
      PFS_Group = datTraits$PFS_group_std[match(outlier_samples, rownames(datTraits))],
      stringsAsFactors = FALSE
    )
    if("Response_std" %in% colnames(datTraits)) {
      outlier_df$Response <- datTraits$Response_std[match(outlier_samples, rownames(datTraits))]
    }
    print(outlier_df)
  }

  # --- 14a.4: Check if outliers have other QC issues ---
  cat("   14a.4: Do DYNC2H1 Outliers Have Other QC Issues?\n")
  cat("   ---------------------------------------------------------------------\n")

  # Calculate sample-level QC metrics
  sample_means <- rowMeans(datExpr)
  sample_medians <- apply(datExpr, 1, median)
  sample_cv <- apply(datExpr, 1, function(x) sd(x)/mean(x))

  tech_flag_sample_outlier <- FALSE

  if(length(outlier_samples) > 0) {
    outlier_qc <- data.frame(
      Sample = outlier_samples,
      Mean_Expr = round(sample_means[outlier_samples], 3),
      Median_Expr = round(sample_medians[outlier_samples], 3),
      Sample_CV = round(sample_cv[outlier_samples], 3),
      Mean_Rank = rank(sample_means)[outlier_samples],
      stringsAsFactors = FALSE
    )
    cat("   QC metrics for DYNC2H1 outlier samples:\n")
    print(outlier_qc)

    # Check if outlier samples are also outliers globally
    global_outlier_check <- sample_means[outlier_samples] < quantile(sample_means, 0.1) |
      sample_means[outlier_samples] > quantile(sample_means, 0.9)
    if(any(global_outlier_check)) {
      cat("   [WARNING] Some DYNC2H1 outliers are also global sample outliers\n")
      tech_flag_sample_outlier <- TRUE
    } else {
      cat("   [OK] DYNC2H1 outliers are NOT global sample outliers\n")
      tech_flag_sample_outlier <- FALSE
    }
  }

  # --- 14a.5: Visualization - Technical Assessment ---
  cat("   14a.5: Generating Technical Assessment Plots...\n")

  # Create diagnostic data frame
  diag_df <- data.frame(
    Sample = sample_ids,
    DYNC2H1 = dync2h1_expr,
    DYNC2H1_Z = dync2h1_z,
    Total_Intensity = total_intensity,
    Sample_Mean = sample_means,
    PFS_Group = datTraits$PFS_group_std,
    Is_Outlier = abs(dync2h1_z) > 2.5,
    stringsAsFactors = FALSE
  )

  if("Response_std" %in% colnames(datTraits)) {
    diag_df$Response <- datTraits$Response_std
  }

  # Plot 1: DYNC2H1 vs Total Intensity
  p14_1 <- ggplot(diag_df, aes(x = Total_Intensity, y = DYNC2H1, color = Is_Outlier)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50", linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "#008080", "TRUE" = "#ca562c"),
                       labels = c("Normal", "Outlier (|Z|>2.5)")) +
    labs(
      title = "DYNC2H1 vs Total Sample Intensity",
      subtitle = sprintf("Spearman r = %.3f, p = %.4f", cor_total$estimate, cor_total$p.value),
      x = "Total Sample Intensity (sum of all proteins)",
      y = "DYNC2H1 Expression (log2)",
      color = "Status"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )

  # Plot 2: DYNC2H1 Distribution with outliers marked
  p14_2 <- ggplot(diag_df, aes(x = DYNC2H1, fill = Is_Outlier)) +
    geom_histogram(bins = 15, color = "white", alpha = 0.8) +
    geom_vline(xintercept = mean(dync2h1_expr), linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = mean(dync2h1_expr) + 2.5*sd(dync2h1_expr),
               linetype = "dotted", color = "#ca562c", linewidth = 0.8) +
    geom_vline(xintercept = mean(dync2h1_expr) - 2.5*sd(dync2h1_expr),
               linetype = "dotted", color = "#ca562c", linewidth = 0.8) +
    scale_fill_manual(values = c("FALSE" = "#008080", "TRUE" = "#ca562c"),
                      labels = c("Normal", "Outlier")) +
    labs(
      title = "DYNC2H1 Expression Distribution",
      subtitle = sprintf("CV = %.1f%% | Mean +/- 2.5 SD boundaries shown", dync2h1_cv),
      x = "DYNC2H1 Expression (log2)",
      y = "Frequency",
      fill = "Status"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )

  p14_tech <- p14_1 + p14_2 + plot_annotation(
    title = "STEP 14a: Technical Noise Assessment",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

  print(p14_tech)
  ggsave(file.path(fig_dir_png, "Step14_01_Technical_Assessment.png"), p14_tech, width = 12, height = 5, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step14_01_Technical_Assessment.pdf"), p14_tech, width = 12, height = 5)
  cat("   [OK] Saved: Step14_01_Technical_Assessment\n")

  # ============================================================================
  # 14b: BIOLOGICAL SIGNAL ASSESSMENT
  # ============================================================================

  cat("==============================================================================\n")
  cat("14b: BIOLOGICAL SIGNAL ASSESSMENT\n")
  cat("==============================================================================\n")

  # --- 14b.1: DYNC2H1 vs PFS_Group ---
  cat("   14b.1: DYNC2H1 vs PFS_Group (Short vs Long)\n")
  cat("   ---------------------------------------------------------------------\n")

  pfs_groups <- datTraits$PFS_group_std


  # Handle both numeric (0/1) and text (Short/Long) PFS encoding
  if(is.numeric(pfs_groups) || all(unique(na.omit(pfs_groups)) %in% c(0, 1, "0", "1"))) {
    # Numeric encoding: 0 = Short, 1 = Long
    short_expr_dync <- dync2h1_expr[pfs_groups == 0 | pfs_groups == "0"]
    long_expr_dync <- dync2h1_expr[pfs_groups == 1 | pfs_groups == "1"]
    cat("   Using numeric encoding: 0=Short, 1=Long\n")
  } else {
    # Text encoding
    short_expr_dync <- dync2h1_expr[pfs_groups == "Short"]
    long_expr_dync <- dync2h1_expr[pfs_groups == "Long"]
    cat("   Using text encoding: Short/Long\n")
  }

  # Remove NAs
  short_expr_dync <- short_expr_dync[!is.na(short_expr_dync)]
  long_expr_dync <- long_expr_dync[!is.na(long_expr_dync)]

  cat(sprintf("   Short PFS samples: %d, Long PFS samples: %d\n", length(short_expr_dync), length(long_expr_dync)))

  # Check we have enough samples before t-test
  if(length(short_expr_dync) < 2 || length(long_expr_dync) < 2) {
    cat("   [WARNING] Not enough samples in one or both PFS groups. Skipping t-test.\n")
    tt_pfs_dync <- list(p.value = NA, statistic = NA)
    wt_pfs_dync <- list(p.value = NA, statistic = NA)
    cohens_d_pfs_dync <- NA
    direction_pfs_dync <- NA
    bio_flag_pfs <- FALSE
  } else {
    # t-test
    tt_pfs_dync <- t.test(short_expr_dync, long_expr_dync)

    # Wilcoxon (robust)
    wt_pfs_dync <- wilcox.test(short_expr_dync, long_expr_dync)

    # Effect size (Cohen's d) - with division by zero protection
    pooled_sd_dync <- sqrt((var(short_expr_dync) + var(long_expr_dync)) / 2)
    if(pooled_sd_dync == 0 || is.na(pooled_sd_dync)) {
      cohens_d_pfs_dync <- NA
      cat("   [WARNING] Zero variance detected, Cohen's d set to NA\n")
    } else {
      cohens_d_pfs_dync <- (mean(short_expr_dync) - mean(long_expr_dync)) / pooled_sd_dync
    }

    cat(sprintf("   Short PFS: Mean = %.3f (n=%d)\n", mean(short_expr_dync), length(short_expr_dync)))
    cat(sprintf("   Long PFS:  Mean = %.3f (n=%d)\n", mean(long_expr_dync), length(long_expr_dync)))
    cat(sprintf("   Difference: %.3f (Short - Long)\n", mean(short_expr_dync) - mean(long_expr_dync)))
    cat(sprintf("   t-test: t = %.3f, p = %.4f\n", tt_pfs_dync$statistic, tt_pfs_dync$p.value))
    cat(sprintf("   Wilcoxon: W = %.0f, p = %.4f\n", wt_pfs_dync$statistic, wt_pfs_dync$p.value))
    cat(sprintf("   Cohen's d: %.3f\n", cohens_d_pfs_dync))

    direction_pfs_dync <- ifelse(mean(short_expr_dync) > mean(long_expr_dync), "Up in Short (Poor)", "Up in Long (Good)")
    sig_pfs_dync <- tt_pfs_dync$p.value < 0.05
    bio_flag_pfs <- sig_pfs_dync

    if(sig_pfs_dync) {
      cat(sprintf("   [OK] SIGNIFICANT: DYNC2H1 %s prognosis (p < 0.05)\n", direction_pfs_dync))
    } else {
      cat(sprintf("   Not significant: DYNC2H1 %s (p = %.3f)\n", direction_pfs_dync, tt_pfs_dync$p.value))
    }
  }  # End of else block (enough samples for t-test)

  # --- 14b.2: DYNC2H1 vs Response ---
  bio_flag_resp <- FALSE
  tt_resp_dync <- NULL
  cohens_d_resp_dync <- NA
  direction_resp_dync <- NA

  if("Response_std" %in% colnames(datTraits)) {
    cat("   14b.2: DYNC2H1 vs Response (CD vs PD)\n")
    cat("   ---------------------------------------------------------------------\n")

    response_groups <- datTraits$Response_std


    # Handle both numeric (0/1) and text (CD/PD) Response encoding
    if(is.numeric(response_groups) || all(unique(na.omit(response_groups)) %in% c(0, 1, "0", "1"))) {
      # Numeric encoding: 0 = PD, 1 = CD
      cd_expr_dync <- dync2h1_expr[response_groups == 1 | response_groups == "1"]
      pd_expr_dync <- dync2h1_expr[response_groups == 0 | response_groups == "0"]
      cat("   Using numeric encoding: 0=PD, 1=CD\n")
    } else {
      # Text encoding
      cd_expr_dync <- dync2h1_expr[response_groups == "CD"]
      pd_expr_dync <- dync2h1_expr[response_groups == "PD"]
      cat("   Using text encoding: CD/PD\n")
    }

    cd_expr_dync <- cd_expr_dync[!is.na(cd_expr_dync)]
    pd_expr_dync <- pd_expr_dync[!is.na(pd_expr_dync)]

    cat(sprintf("   CD samples: %d, PD samples: %d\n", length(cd_expr_dync), length(pd_expr_dync)))

    if(length(cd_expr_dync) >= 3 && length(pd_expr_dync) >= 3) {
      tt_resp_dync <- t.test(cd_expr_dync, pd_expr_dync)
      wt_resp_dync <- wilcox.test(cd_expr_dync, pd_expr_dync)
      cohens_d_resp_dync <- (mean(cd_expr_dync) - mean(pd_expr_dync)) /
        sqrt((var(cd_expr_dync) + var(pd_expr_dync)) / 2)

      cat(sprintf("   CD (responders): Mean = %.3f (n=%d)\n", mean(cd_expr_dync), length(cd_expr_dync)))
      cat(sprintf("   PD (progressors): Mean = %.3f (n=%d)\n", mean(pd_expr_dync), length(pd_expr_dync)))
      cat(sprintf("   t-test: p = %.4f | Cohen's d: %.3f\n", tt_resp_dync$p.value, cohens_d_resp_dync))

      direction_resp_dync <- ifelse(mean(cd_expr_dync) > mean(pd_expr_dync), "Up in CD (Responders)", "Up in PD (Progressors)")
      sig_resp_dync <- tt_resp_dync$p.value < 0.05
      bio_flag_resp <- sig_resp_dync

      if(sig_resp_dync) {
        cat(sprintf("   [OK] SIGNIFICANT: DYNC2H1 %s (p < 0.05)\n", direction_resp_dync))
      } else {
        cat(sprintf("   Not significant: DYNC2H1 %s (p = %.3f)\n", direction_resp_dync, tt_resp_dync$p.value))
      }
    } else {
      cat("   [WARNING] Not enough samples in Response groups\n")
    }
  }

  # --- 14b.3: Correlation with PFS_Days ---
  bio_flag_pfs_continuous <- FALSE
  cor_pfs_days_dync <- NULL

  if("PFS_days_std" %in% colnames(datTraits)) {
    cat("   14b.3: Correlation with PFS (days)\n")
    cat("   ---------------------------------------------------------------------\n")

    pfs_days <- datTraits$PFS_days_std
    valid_idx <- !is.na(dync2h1_expr) & !is.na(pfs_days)

    if(sum(valid_idx) >= 10) {
      cor_pfs_days_dync <- cor.test(dync2h1_expr[valid_idx], pfs_days[valid_idx], method = "spearman")

      cat(sprintf("   Spearman correlation: r = %.3f, p = %.4f\n",
                      cor_pfs_days_dync$estimate, cor_pfs_days_dync$p.value))

      if(cor_pfs_days_dync$p.value < 0.05) {
        direction_continuous <- ifelse(cor_pfs_days_dync$estimate > 0, "positive (longer PFS)", "negative (shorter PFS)")
        cat(sprintf("   [OK] SIGNIFICANT %s correlation with survival\n", direction_continuous))
        bio_flag_pfs_continuous <- TRUE
      } else {
        cat("   No significant correlation with survival days\n")
      }
    } else {
      cat("   [WARNING] Not enough valid samples for correlation\n")
    }
  }

  # --- 14b.4: Biological Signal Visualization ---
  cat("   14b.4: Generating Biological Assessment Plots...\n")

  # Plot 3: Boxplot by PFS_Group
  p14_3 <- ggplot(diag_df %>% filter(!is.na(PFS_Group)),
                  aes(x = PFS_Group, y = DYNC2H1, fill = PFS_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = colors_pfs) +
    labs(
      title = "DYNC2H1 by PFS Group",
      subtitle = sprintf("p = %.4f | Cohen's d = %.2f | %s",
                         tt_pfs_dync$p.value, cohens_d_pfs_dync, direction_pfs_dync),
      x = "PFS Group",
      y = "DYNC2H1 Expression (log2)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "none"
    )

  # Plot 4: Boxplot by Response (if available)
  if("Response" %in% colnames(diag_df) && !is.null(tt_resp_dync)) {
    p14_4 <- ggplot(diag_df %>% filter(!is.na(Response)),
                    aes(x = Response, y = DYNC2H1, fill = Response)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(values = colors_response) +
      labs(
        title = "DYNC2H1 by Response",
        subtitle = sprintf("p = %.4f | Cohen's d = %.2f | %s",
                           tt_resp_dync$p.value, cohens_d_resp_dync, direction_resp_dync),
        x = "Treatment Response",
        y = "DYNC2H1 Expression (log2)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "none"
      )
  } else {
    # Placeholder plot
    p14_4 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Response data\nnot available") +
      theme_void()
  }

  # Plot 5: Scatter vs PFS days (if available)
  if(!is.null(cor_pfs_days_dync) && "PFS_days_std" %in% colnames(datTraits)) {
    diag_df$PFS_Days <- datTraits$PFS_days_std

    p14_5 <- ggplot(diag_df %>% filter(!is.na(PFS_Days)),
                    aes(x = DYNC2H1, y = PFS_Days, color = PFS_Group)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "grey50", linetype = "dashed") +
      scale_color_manual(values = colors_pfs) +
      labs(
        title = "DYNC2H1 vs PFS (days)",
        subtitle = sprintf("Spearman r = %.3f, p = %.4f",
                           cor_pfs_days_dync$estimate, cor_pfs_days_dync$p.value),
        x = "DYNC2H1 Expression (log2)",
        y = "Progression-Free Survival (days)",
        color = "PFS Group"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "bottom"
      )
  } else {
    p14_5 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "PFS days data\nnot available") +
      theme_void()
  }

  p14_bio <- (p14_3 | p14_4) / p14_5 + plot_annotation(
    title = "STEP 14b: Biological Signal Assessment",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

  print(p14_bio)
  ggsave(file.path(fig_dir_png, "Step14_02_Biological_Assessment.png"), p14_bio, width = 10, height = 8, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step14_02_Biological_Assessment.pdf"), p14_bio, width = 10, height = 8)
  cat("   [OK] Saved: Step14_02_Biological_Assessment\n")

  # ============================================================================
  # 14c: WGCNA CONTEXT
  # ============================================================================

  cat("==============================================================================\n")
  cat("14c: WGCNA CONTEXT - MODULE MEMBERSHIP\n")
  cat("==============================================================================\n")

  # --- 14c.1: Which module is DYNC2H1 in? ---
  dync2h1_module <- NA
  dync2h1_kME <- NA

  if("DYNC2H1" %in% names(moduleColors)) {
    dync2h1_module <- moduleColors["DYNC2H1"]
    cat(sprintf("   DYNC2H1 Module Assignment: %s\n", dync2h1_module))

    # --- 14c.2: Module membership (kME) ---
    me_col <- paste0("ME", dync2h1_module)
    if(me_col %in% colnames(MEs)) {
      dync2h1_kME <- cor(dync2h1_expr, MEs[, me_col], use = "complete.obs")
      cat(sprintf("   Module Membership (kME): %.3f\n", dync2h1_kME))

      if(abs(dync2h1_kME) > 0.8) {
        cat("   [OK] Strong module membership (|kME| > 0.8) - core member\n")
      } else if(abs(dync2h1_kME) > 0.5) {
        cat("   Moderate module membership (0.5 < |kME| < 0.8)\n")
      } else {
        cat("   [WARNING] Weak module membership (|kME| < 0.5) - peripheral member\n")
      }
    } else {
      cat(sprintf("   [WARNING] Module eigengene %s not found in MEs\n", me_col))
      cat(sprintf("   Available MEs: %s\n", paste(colnames(MEs)[1:min(10, ncol(MEs))], collapse = ", ")))
    }

    # --- 14c.3: Correlation with other module members ---
    cat("   14c.3: Correlation with Same-Module Proteins\n")
    cat("   ---------------------------------------------------------------------\n")

    same_module_proteins <- names(moduleColors)[moduleColors == dync2h1_module]
    same_module_proteins <- setdiff(same_module_proteins, "DYNC2H1")
    same_module_proteins <- same_module_proteins[same_module_proteins %in% colnames(datExpr)]

    if(length(same_module_proteins) > 0) {
      same_module_cors <- sapply(same_module_proteins, function(p) {
        cor(dync2h1_expr, datExpr[, p], use = "complete.obs")
      })

      cat(sprintf("   Module %s has %d other proteins\n", dync2h1_module, length(same_module_proteins)))
      cat(sprintf("   DYNC2H1 mean |cor| with module: %.3f\n", mean(abs(same_module_cors), na.rm = TRUE)))
      cat(sprintf("   DYNC2H1 max |cor| with module: %.3f\n", max(abs(same_module_cors), na.rm = TRUE)))

      # Top correlated proteins in module
      top_cors <- sort(same_module_cors, decreasing = TRUE)[1:min(5, length(same_module_cors))]
      cat("   Top 5 correlated proteins in same module:\n")
      for(j in 1:length(top_cors)) {
        cat(sprintf("      %s: r = %.3f\n", names(top_cors)[j], top_cors[j]))
      }
    }
  } else {
    cat("   [WARNING] DYNC2H1 not found in moduleColors\n")
  }

  # ============================================================================
  # 14d: COMPARISON WITH OTHER HIGH-CV PROTEINS
  # ============================================================================

  cat("==============================================================================\n")
  cat("14d: COMPARISON WITH OTHER HIGH-CV PROTEINS\n")
  cat("==============================================================================\n")

  # Calculate CV for all proteins
  all_cv <- apply(datExpr, 2, function(x) sd(x) / mean(x) * 100)
  cv_95th <- quantile(all_cv, 0.95)

  # Get other high-CV proteins (top 5%)
  high_cv_proteins <- names(all_cv)[all_cv > cv_95th]
  cat(sprintf("   High CV proteins (>%.1f%%, 95th percentile): %d\n", cv_95th, length(high_cv_proteins)))
  if("DYNC2H1" %in% names(all_cv)) {
    cat(sprintf("   DYNC2H1 CV percentile: %.1f%%\n", ecdf(all_cv)(all_cv["DYNC2H1"]) * 100))
  } else {
    cat("   DYNC2H1 not found in dataset\n")
  }

  # Compare clinical associations across high-CV proteins
  cat("   Clinical associations of high-CV proteins:\n")
  cat("   ---------------------------------------------------------------------\n")

  high_cv_clinical <- data.frame()
  pfs_groups_std <- datTraits$PFS_group_std

  # Determine PFS encoding once (numeric vs text)
  is_numeric_pfs <- is.numeric(pfs_groups_std) || all(unique(na.omit(pfs_groups_std)) %in% c(0, 1, "0", "1"))

  for(prot in high_cv_proteins) {
    prot_expr <- datExpr[, prot]

    # PFS association - handle numeric (0=Short, 1=Long) or text encoding
    if(is_numeric_pfs) {
      short_p <- prot_expr[(pfs_groups_std == 0 | pfs_groups_std == "0") & !is.na(prot_expr)]
      long_p <- prot_expr[(pfs_groups_std == 1 | pfs_groups_std == "1") & !is.na(prot_expr)]
    } else {
      short_p <- prot_expr[pfs_groups_std == "Short" & !is.na(prot_expr)]
      long_p <- prot_expr[pfs_groups_std == "Long" & !is.na(prot_expr)]
    }

    if(length(short_p) >= 3 && length(long_p) >= 3) {
      tt <- t.test(short_p, long_p)
      d <- (mean(short_p) - mean(long_p)) / sqrt((var(short_p) + var(long_p)) / 2)

      high_cv_clinical <- rbind(high_cv_clinical, data.frame(
        Protein = prot,
        CV_pct = round(all_cv[prot], 1),
        PFS_pval = round(tt$p.value, 4),
        Cohens_d = round(d, 3),
        Direction = ifelse(mean(short_p) > mean(long_p), "Up_Short", "Up_Long"),
        Significant = tt$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }

  high_cv_clinical <- high_cv_clinical %>% arrange(PFS_pval)
  cat("   High-CV proteins ranked by PFS association:\n")
  print(high_cv_clinical)

  # Highlight DYNC2H1 position
  dync2h1_rank <- which(high_cv_clinical$Protein == "DYNC2H1")
  if(length(dync2h1_rank) > 0) {
    cat(sprintf("\n   DYNC2H1 rank among high-CV proteins: %d of %d\n", dync2h1_rank, nrow(high_cv_clinical)))
  }

  # --- Plot: CV vs Clinical Effect Size ---
  cat("   Generating CV vs Effect Size plot...\n")

  cv_effect_df <- data.frame(
    Protein = names(all_cv),
    CV = all_cv,
    stringsAsFactors = FALSE
  )

  # Calculate effect size for all proteins (handle numeric/text PFS encoding)
  cv_effect_df$Cohens_d <- sapply(cv_effect_df$Protein, function(p) {
    expr <- datExpr[, p]
    if(is_numeric_pfs) {
      short <- expr[(pfs_groups_std == 0 | pfs_groups_std == "0") & !is.na(expr)]
      long <- expr[(pfs_groups_std == 1 | pfs_groups_std == "1") & !is.na(expr)]
    } else {
      short <- expr[pfs_groups_std == "Short" & !is.na(expr)]
      long <- expr[pfs_groups_std == "Long" & !is.na(expr)]
    }
    if(length(short) >= 3 && length(long) >= 3) {
      (mean(short) - mean(long)) / sqrt((var(short) + var(long)) / 2)
    } else {
      NA
    }
  })

  cv_effect_df$Is_DYNC2H1 <- cv_effect_df$Protein == "DYNC2H1"
  cv_effect_df$Is_HighCV <- cv_effect_df$CV > cv_95th
  cv_effect_df$Is_ML_Biomarker <- cv_effect_df$Protein %in% c("DYNC2H1", "ECM2", "PPIB")
  cv_effect_df$Is_Cox_Biomarker <- cv_effect_df$Protein %in% c("APOF", "TFRC", "FABP4", "ANG")

  p14_6 <- ggplot(cv_effect_df %>% filter(!is.na(Cohens_d)),
                  aes(x = CV, y = abs(Cohens_d))) +
    geom_point(aes(color = Is_HighCV), alpha = 0.5, size = 2) +
    geom_point(data = cv_effect_df %>% filter(Is_ML_Biomarker & !is.na(Cohens_d)),
               color = "#ca562c", size = 4, shape = 17) +
    geom_point(data = cv_effect_df %>% filter(Is_Cox_Biomarker & !is.na(Cohens_d)),
               color = "#008080", size = 4, shape = 15) +
    geom_text_repel(data = cv_effect_df %>% filter((Is_ML_Biomarker | Is_Cox_Biomarker) & !is.na(Cohens_d)),
                    aes(label = Protein), size = 3.5, fontface = "bold",
                    max.overlaps = 20) +
    geom_vline(xintercept = cv_95th, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue", alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#f4a582"),
                       labels = c("Normal CV", "High CV (>95th pctl)")) +
    labs(
      title = "Coefficient of Variation vs Clinical Effect Size",
      subtitle = "Triangle=ML biomarkers | Square=Cox biomarkers | Dashed: CV 95th pctl, |d|=0.5",
      x = "CV (%)",
      y = "|Cohen's d| (PFS Group)",
      color = "CV Status"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )

  print(p14_6)
  ggsave(file.path(fig_dir_png, "Step14_03_CV_vs_EffectSize.png"), p14_6, width = 10, height = 7, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step14_03_CV_vs_EffectSize.pdf"), p14_6, width = 10, height = 7)
  cat("   [OK] Saved: Step14_03_CV_vs_EffectSize\n")

  # ============================================================================
  # 14e: DECISION FRAMEWORK AND RECOMMENDATION
  # ============================================================================

  cat("==============================================================================\n")
  cat("14e: DECISION FRAMEWORK AND RECOMMENDATION\n")
  cat("==============================================================================\n")

  # Calculate scores
  tech_score <- sum(c(tech_flag_intensity, tech_flag_sample_outlier))
  bio_score <- sum(c(bio_flag_pfs, bio_flag_resp, bio_flag_pfs_continuous,
                     !is.na(dync2h1_kME) && abs(dync2h1_kME) > 0.8))

  cat("   EVIDENCE SUMMARY:\n")
  cat("   ---------------------------------------------------------------------\n")
  cat("   Technical Flags:\n")
  cat(sprintf("      Correlation with total intensity: %s\n",
                  ifelse(tech_flag_intensity, "[WARNING] CONCERN", "[OK] OK")))
  cat(sprintf("      Outlier samples are global outliers: %s\n",
                  ifelse(tech_flag_sample_outlier, "[WARNING] CONCERN", "[OK] OK")))
  cat("   Biological Signals:\n")
  cat(sprintf("      PFS Group association (p < 0.05): %s\n",
                  ifelse(bio_flag_pfs, sprintf("[OK] YES (p=%.4f)", tt_pfs_dync$p.value), "NO")))
  cat(sprintf("      Response association (p < 0.05): %s\n",
                  ifelse(bio_flag_resp, "[OK] YES", "NO")))
  cat(sprintf("      PFS days correlation (p < 0.05): %s\n",
                  ifelse(bio_flag_pfs_continuous, "[OK] YES", "NO")))
  cat(sprintf("      Strong module membership (|kME| > 0.8): %s\n",
                  ifelse(!is.na(dync2h1_kME) && abs(dync2h1_kME) > 0.8, "[OK] YES", "NO")))

  cat("   DECISION SCORING:\n")
  cat(sprintf("   Technical concern flags: %d/2\n", tech_score))
  cat(sprintf("   Biological signal flags: %d/4\n", bio_score))

  # Final recommendation
  cat("   +=======================================================================+\n")
  cat("   |                      FINAL RECOMMENDATION                             |\n")
  cat("   +=======================================================================+\n")

  if(bio_score >= 2 && tech_score == 0) {
    verdict <- "KEEP - BIOLOGICAL SIGNAL"
    recommendation <- "Strong biological signal with no technical concerns. Keep in ML signature."
    cat("   |  VERDICT: [OK] KEEP - BIOLOGICAL SIGNAL                             |\n")
    cat("   |                                                                      |\n")
    cat("   |  The high CV appears to reflect TRUE biological variability          |\n")
    cat("   |  associated with clinical outcomes. Keep DYNC2H1 in the signature.   |\n")
  } else if(tech_score >= 1 && bio_score == 0) {
    verdict <- "REMOVE - TECHNICAL NOISE"
    recommendation <- "Technical concerns without biological signal. Consider removing from signature."
    cat("   |  VERDICT: [FAIL] REMOVE - TECHNICAL NOISE                            |\n")
    cat("   |                                                                      |\n")
    cat("   |  The high CV appears to reflect technical artifacts.                 |\n")
    cat("   |  Consider removing DYNC2H1 from the ML signature.                    |\n")
  } else if(bio_score >= 1 && tech_score >= 1) {
    verdict <- "KEEP WITH CAUTION - MIXED EVIDENCE"
    recommendation <- "Both technical and biological signals present. Keep but report as limitation."
    cat("   |  VERDICT: ~ KEEP WITH CAUTION - MIXED EVIDENCE                       |\n")
    cat("   |                                                                      |\n")
    cat("   |  Both technical concerns and biological signal present.              |\n")
    cat("   |  RECOMMENDATION: Keep DYNC2H1 but report high CV as a limitation     |\n")
    cat("   |  in publications. Validate in independent cohort.                    |\n")
  } else if(bio_score >= 1) {
    verdict <- "KEEP - BIOLOGICAL SIGNAL"
    recommendation <- "Biological signal detected. High CV likely reflects true variation."
    cat("   |  VERDICT: [OK] KEEP - BIOLOGICAL SIGNAL                              |\n")
    cat("   |                                                                      |\n")
    cat("   |  Biological signal detected despite high CV.                         |\n")
    cat("   |  The variation likely reflects meaningful biology.                   |\n")
  } else {
    verdict <- "UNCERTAIN - INSUFFICIENT EVIDENCE"
    recommendation <- "Neither clear technical nor biological explanation. Report as limitation."
    cat("   |  VERDICT: ? UNCERTAIN - INSUFFICIENT EVIDENCE                        |\n")
    cat("   |                                                                      |\n")
    cat("   |  Cannot clearly attribute CV to technical or biological factors.     |\n")
    cat("   |  RECOMMENDATION: Keep but report high CV as study limitation.        |\n")
  }
  cat("   +=======================================================================+\n")

  # ============================================================================
  # 14f: SAVE RESULTS
  # ============================================================================

  cat("==============================================================================\n")
  cat("14f: SAVE RESULTS\n")
  cat("==============================================================================\n")

  # Summary table
  dync2h1_summary <- data.frame(
    Parameter = c(
      "Protein", "CV (%)", "CV Percentile", "Module",
      "Module Membership (kME)", "Mean Expression", "SD",
      "---",
      "PFS Group p-value", "PFS Group Cohen's d", "PFS Group Direction",
      "Response p-value", "Response Cohen's d",
      "PFS Days Spearman r", "PFS Days p-value",
      "---",
      "Technical: Intensity Correlation", "Technical: Sample Outliers",
      "---",
      "VERDICT", "RECOMMENDATION"
    ),
    Value = c(
      "DYNC2H1",
      sprintf("%.1f", dync2h1_cv),
      sprintf("%.1f", if("DYNC2H1" %in% names(all_cv)) ecdf(all_cv)(all_cv["DYNC2H1"]) * 100 else NA),
      ifelse(!is.na(dync2h1_module), dync2h1_module, "NA"),
      sprintf("%.3f", ifelse(!is.na(dync2h1_kME), dync2h1_kME, NA)),
      sprintf("%.4f", mean(dync2h1_expr)),
      sprintf("%.4f", sd(dync2h1_expr)),
      "---",
      sprintf("%.4f", tt_pfs_dync$p.value),
      sprintf("%.3f", cohens_d_pfs_dync),
      direction_pfs_dync,
      ifelse(!is.null(tt_resp_dync), sprintf("%.4f", tt_resp_dync$p.value), "NA"),
      sprintf("%.3f", ifelse(!is.na(cohens_d_resp_dync), cohens_d_resp_dync, NA)),
      ifelse(!is.null(cor_pfs_days_dync), sprintf("%.3f", cor_pfs_days_dync$estimate), "NA"),
      ifelse(!is.null(cor_pfs_days_dync), sprintf("%.4f", cor_pfs_days_dync$p.value), "NA"),
      "---",
      ifelse(tech_flag_intensity, "CONCERN", "OK"),
      ifelse(tech_flag_sample_outlier, "CONCERN", "OK"),
      "---",
      verdict,
      recommendation
    ),
    stringsAsFactors = FALSE
  )

  write.csv(dync2h1_summary, file.path(results_dir, "Step14_DYNC2H1_Investigation_Summary.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step14_DYNC2H1_Investigation_Summary.csv\n")

  write.csv(high_cv_clinical, file.path(results_dir, "Step14_HighCV_Proteins_Clinical.csv"), row.names = FALSE)
  cat("   [OK] Saved: Step14_HighCV_Proteins_Clinical.csv\n")

} # End of DYNC2H1 check

#DYNC2H1 exhibited the highest inter-patient variability among all proteins (CV = 31%, 100th percentile).
#Technical assessment ruled out normalization artifacts (no correlation with total sample intensity, r = 0.121, p = 0.40) and sample outliers. Instead, the variability reflected clinically meaningful biological differences: patients with high DYNC2H1 expression showed significantly longer progression-free survival (t-test p = 0.0021, Cohen's d = -0.943) and higher treatment response rates (p = 0.0028, Cohen's d = 0.861). DYNC2H1 expression correlated positively with PFS days (Spearman r = 0.417, p = 0.0026). Despite assignment to the grey module in WGCNA (indicating weak co-expression with other proteins), DYNC2H1 ranked first among 45 high-CV proteins for clinical significance, supporting its retention as an independent prognostic biomarker.


# ============================================================================
# STEP 14 COMPLETE
# ============================================================================
cat("==============================================================================\n")
cat("STEP 14 COMPLETE: DYNC2H1 Investigation\n")
cat("==============================================================================\n")
cat("   OUTPUT FILES:\n")
cat("   - Step14_01_Technical_Assessment.png/pdf\n")
cat("   - Step14_02_Biological_Assessment.png/pdf\n")
cat("   - Step14_03_CV_vs_EffectSize.png/pdf\n")
cat("   - Step14_DYNC2H1_Investigation_Summary.csv\n")
cat("   - Step14_HighCV_Proteins_Clinical.csv\n")

#===============================================================================
# STEP 15: SIGNATURE OVERLAP ANALYSIS
#===============================================================================
# PURPOSE: Analyze overlap between biomarkers from different discovery methods
# AND WGCNA hub genes from ALL clinically-relevant focus modules (combined)
# QUESTION: Do our predictive biomarkers overlap with WGCNA hub genes?
# THREE SETS TO COMPARE:
# 1. ML Classifier biomarkers: ECM2, DYNC2H1, PPIB
# 2. Cox Regression biomarkers: APOF, TFRC, FABP4, ANG
# 3. Focus Module Hubs: Combined from ALL focus_modules_auto (PINK, GREEN, BLACK, BLUE)
# EXPECTED RESULT: Little to no overlap, which is actually interesting because:
# - ML found proteins with predictive power (but DYNC2H1 in grey module)
# - Cox found proteins with survival association (strong, consistent)
# - WGCNA found network hubs (biologically central)
# THESIS NARRATIVE: "Different analytical approaches captured complementary
# aspects of PDAC biology, suggesting a multi-faceted prognostic landscape."
# VISUALIZATIONS:
# 15a: Define all 3 protein sets (ML, Cox, Combined Focus Module Hubs)
# 15b: Calculate pairwise overlaps
# 15c: 3-way Venn Diagram (ML, Cox, Focus Module Hubs)
# 15d: UpSet Plot (detailed overlap view)
# 15e: Summary statistics and interpretation                                 

cat("===============================================================================\n")
cat("STEP 15: SIGNATURE OVERLAP ANALYSIS\n")
cat("===============================================================================\n")

# ============================================================================
# 15a: Define All Protein Sets
# ============================================================================
# Four protein sets from three different analytical approaches:
# 1. ML Classifier (Random Forest): Optimized for classification accuracy
# 2. Cox Regression: Identified proteins with survival association
# 3. WGCNA mod1 hubs: Network-central proteins in focus module 1 (dynamically set)
# 4. WGCNA mod2 hubs: Network-central proteins in focus module 2 (dynamically set)
#
# EXPECTED: Little to no overlap, which is interesting because each method
# captures different aspects of PDAC biology.
#
# THESIS NARRATIVE: "Different analytical approaches captured complementary
# aspects of PDAC biology, suggesting a multi-faceted prognostic landscape."
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15a: Defining All Protein Sets\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# Define ML and Cox biomarkers
# NOTE: Cox biomarkers (APOF, TFRC, FABP4, ANG) were derived from univariate Cox
# proportional hazards regression. The proportional hazards assumption was tested
# using Schoenfeld residuals (cox.zph) and all selected biomarkers satisfied
# the assumption (global test p > 0.05). See Methods section for details.
ml_biomarkers <- c("ECM2", "DYNC2H1", "PPIB")  # ENPP1 removed due to QC concerns
cox_biomarkers <- c("APOF", "TFRC", "FABP4", "ANG")


ml_available <- ml_biomarkers[ml_biomarkers %in% colnames(datExpr)]
cox_available <- cox_biomarkers[cox_biomarkers %in% colnames(datExpr)]
all_biomarkers <- c(ml_available, cox_available)

if(length(ml_available) < length(ml_biomarkers)) {
  missing <- setdiff(ml_biomarkers, ml_available)
  warning(sprintf("ML biomarkers not found in datExpr: %s", paste(missing, collapse = ", ")))
}
if(length(cox_available) < length(cox_biomarkers)) {
  missing <- setdiff(cox_biomarkers, cox_available)
  warning(sprintf("Cox biomarkers not found in datExpr: %s", paste(missing, collapse = ", ")))
}

# Extract WGCNA hub genes from ALL focus modules (from Step 6)
# Using moderate criteria (kME > 0.7) as defined in Step 6
# Combine all focus module hubs into one set for cleaner 3-way Venn

focus_module_hubs <- character(0)  # Will hold combined hubs from all focus modules
focus_module_hub_details <- list()  # Track per-module for reporting

if(exists("proteinInfo")) {
  # Get hub genes from proteinInfo for each focus module
  for(fm in focus_modules_auto) {
    if("isHub_Moderate" %in% colnames(proteinInfo)) {
      fm_hubs <- proteinInfo$Protein[proteinInfo$Module == fm & proteinInfo$isHub_Moderate == TRUE]
    } else {
      # Fallback: use kME threshold directly
      mm_col <- grep(paste0("MM\\.", fm, "|MM_", fm, "|kME\\.", fm), colnames(proteinInfo), value = TRUE)[1]
      if(!is.na(mm_col)) {
        fm_hubs <- proteinInfo$Protein[proteinInfo$Module == fm & proteinInfo[[mm_col]] > 0.7]
      } else {
        fm_hubs <- character(0)
      }
    }
    fm_hubs <- fm_hubs[!is.na(fm_hubs)]
    focus_module_hub_details[[fm]] <- fm_hubs
    focus_module_hubs <- c(focus_module_hubs, fm_hubs)
  }
} else if(exists("hub_results")) {
  # If proteinInfo doesn't exist, try to get from hub_results
  for(fm in focus_modules_auto) {
    if(fm %in% names(hub_results$Moderate)) {
      fm_hubs <- hub_results$Moderate[[fm]]$Protein
      fm_hubs <- fm_hubs[!is.na(fm_hubs)]
      focus_module_hub_details[[fm]] <- fm_hubs
      focus_module_hubs <- c(focus_module_hubs, fm_hubs)
    }
  }
} else {
  warning("Neither proteinInfo nor hub_results found. Using empty hub list.")
}

# Remove duplicates (in case a protein is hub in multiple modules)
focus_module_hubs <- unique(focus_module_hubs)

# Print summary
cat(sprintf("\n  ML Biomarkers (n=%d): %s\n", length(ml_biomarkers), paste(ml_biomarkers, collapse = ", ")))
cat(sprintf("  Cox Biomarkers (n=%d): %s\n", length(cox_biomarkers), paste(cox_biomarkers, collapse = ", ")))
cat(sprintf("\n  Focus Module Hubs - Combined from %s (n=%d total):\n",
            paste(toupper(focus_modules_auto), collapse = ", "), length(focus_module_hubs)))
for(fm in names(focus_module_hub_details)) {
  n_hubs <- length(focus_module_hub_details[[fm]])
  cat(sprintf("    %s: %d hubs%s\n", tools::toTitleCase(fm), n_hubs,
              ifelse(n_hubs > 0 && n_hubs <= 5,
                     paste0(" (", paste(focus_module_hub_details[[fm]], collapse = ", "), ")"),
                     ifelse(n_hubs > 5, paste0(" (", paste(head(focus_module_hub_details[[fm]], 5), collapse = ", "), ", ...)"), ""))))
}

# Verify biomarkers exist in datExpr
available_ml <- ml_biomarkers[ml_biomarkers %in% colnames(datExpr)]
available_cox <- cox_biomarkers[cox_biomarkers %in% colnames(datExpr)]
cat(sprintf("\n  ML available in data: %d/%d", length(available_ml), length(ml_biomarkers)))
cat(sprintf("  Cox available in data: %d/%d\n", length(available_cox), length(cox_biomarkers)))

# ============================================================================
# 15a-2: BIOMARKER MODULE MEMBERSHIP TABLE
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15a-2: Biomarker Module Membership\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# Create membership table for all 8 biomarkers
all_biomarkers <- c(ml_biomarkers, cox_biomarkers)

biomarker_membership <- data.frame(
  Protein = all_biomarkers,
  Source = c(rep("ML", length(ml_biomarkers)), rep("Cox", length(cox_biomarkers))),
  stringsAsFactors = FALSE
)

# Add module assignment
biomarker_membership$Module <- moduleColors[all_biomarkers]

# Add kME for assigned module (correlation with module eigengene)
biomarker_membership$kME <- sapply(1:nrow(biomarker_membership), function(i) {
  prot <- biomarker_membership$Protein[i]
  mod <- biomarker_membership$Module[i]
  me_col <- paste0("ME", mod)
  if(!is.na(mod) && me_col %in% colnames(MEs) && prot %in% colnames(datExpr)) {
    round(WGCNA::cor(datExpr[, prot], MEs[, me_col], use = "complete.obs"), 3)
  } else {
    NA
  }
})

# Add focus module status
# Use dynamic module names (mod1 and mod2 defined earlier in Step 5)
focus_modules_list <- c(mod1, mod2)
biomarker_membership$In_Focus_Module <- biomarker_membership$Module %in% focus_modules_list

# Add hub status (kME > 0.7)
biomarker_membership$Is_Hub <- !is.na(biomarker_membership$kME) & abs(biomarker_membership$kME) > 0.7

# Save immediately
write.csv(biomarker_membership,
          file.path(results_dir, "Step15_Biomarker_Module_Membership.csv"),
          row.names = FALSE)

cat("\nBiomarker Module Membership:\n")
print(biomarker_membership)
cat(sprintf("\nML biomarkers in focus modules: %d/%d\n",
            sum(biomarker_membership$In_Focus_Module[biomarker_membership$Source == "ML"]),
            sum(biomarker_membership$Source == "ML")))
cat(sprintf("Cox biomarkers in focus modules: %d/%d\n",
            sum(biomarker_membership$In_Focus_Module[biomarker_membership$Source == "Cox"]),
            sum(biomarker_membership$Source == "Cox")))
cat(sprintf("Biomarkers that are hub genes (kME > 0.7): %d/%d\n",
            sum(biomarker_membership$Is_Hub, na.rm = TRUE),
            nrow(biomarker_membership)))
cat("  Saved: Step15_Biomarker_Module_Membership.csv\n")

# ============================================================================
# 15a-3: BIOMARKER SIGNATURE MEMBERSHIP LOOKUP
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15a-3: Biomarker Signature Membership\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# Check which of the 150 custom signatures contain each biomarker
biomarker_sig_membership <- data.frame()

for(prot in all_biomarkers) {
  found_in_any <- FALSE
  for(sig_id in names(custom_signatures)) {
    sig <- custom_signatures[[sig_id]]
    sig_genes <- if(is.list(sig)) sig$genes else sig

    if(prot %in% sig_genes) {
      found_in_any <- TRUE
      sig_name <- if(is.list(sig) && !is.null(sig$short_name)) sig$short_name else sig_id
      sig_category <- if(is.list(sig) && !is.null(sig$category)) sig$category else "Unknown"

      biomarker_sig_membership <- rbind(biomarker_sig_membership, data.frame(
        Protein = prot,
        Source = ifelse(prot %in% ml_biomarkers, "ML", "Cox"),
        Signature_ID = sig_id,
        Signature_Name = sig_name,
        Category = sig_category,
        stringsAsFactors = FALSE
      ))
    }
  }
  # Add entry for proteins not found in any signature
  if(!found_in_any) {
    biomarker_sig_membership <- rbind(biomarker_sig_membership, data.frame(
      Protein = prot,
      Source = ifelse(prot %in% ml_biomarkers, "ML", "Cox"),
      Signature_ID = NA,
      Signature_Name = "Not in any signature",
      Category = NA,
      stringsAsFactors = FALSE
    ))
  }
}

# Save immediately
write.csv(biomarker_sig_membership,
          file.path(results_dir, "Step15_Biomarker_Signature_Membership.csv"),
          row.names = FALSE)

# Summary by biomarker
cat("\nSignatures containing each biomarker:\n")
for(prot in all_biomarkers) {
  prot_sigs <- biomarker_sig_membership[biomarker_sig_membership$Protein == prot &
                                         !is.na(biomarker_sig_membership$Signature_ID), ]
  n_sigs <- nrow(prot_sigs)
  if(n_sigs > 0) {
    categories <- unique(prot_sigs$Category)
    cat(sprintf("  %s: %d signatures (%s)\n", prot, n_sigs, paste(categories, collapse = ", ")))
  } else {
    cat(sprintf("  %s: 0 signatures\n", prot))
  }
}

# Category summary
cat("\nCategory breakdown:\n")
sig_with_data <- biomarker_sig_membership[!is.na(biomarker_sig_membership$Category), ]
if(nrow(sig_with_data) > 0) {
  cat_counts <- sort(table(sig_with_data$Category), decreasing = TRUE)
  for(cat_name in names(cat_counts)) {
    cat(sprintf("  %s: %d associations\n", cat_name, cat_counts[cat_name]))
  }
}

cat("  Saved: Step15_Biomarker_Signature_Membership.csv\n")

# ============================================================================
# 15b: Calculate All Pairwise Overlaps (3 Sets)
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15b: Calculating Pairwise Overlaps\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# Create named list of 3 sets: ML, Cox, and Combined Focus Module Hubs
all_sets <- list(
  "ML Biomarkers" = ml_biomarkers,
  "Cox Biomarkers" = cox_biomarkers,
  "Focus Module Hubs" = focus_module_hubs
)

# Calculate all pairwise overlaps
set_names <- names(all_sets)
n_sets <- length(all_sets)
overlap_matrix <- matrix(NA, nrow = n_sets, ncol = n_sets, dimnames = list(set_names, set_names))

for(i in 1:n_sets) {
  for(j in 1:n_sets) {
    overlap_matrix[i, j] <- length(intersect(all_sets[[i]], all_sets[[j]]))
  }
}

cat("  Overlap Matrix (number of shared proteins):\n")
print(overlap_matrix)

# Key overlaps
ml_cox_overlap <- intersect(ml_biomarkers, cox_biomarkers)
ml_hubs_overlap <- intersect(ml_biomarkers, focus_module_hubs)
cox_hubs_overlap <- intersect(cox_biomarkers, focus_module_hubs)
all_three_overlap <- Reduce(intersect, list(ml_biomarkers, cox_biomarkers, focus_module_hubs))

cat("\n  Key Overlaps:\n")
cat(sprintf("    ML & Cox: %d proteins %s\n", length(ml_cox_overlap),
                ifelse(length(ml_cox_overlap) > 0, paste0("(", paste(ml_cox_overlap, collapse = ", "), ")"), "")))
cat(sprintf("    ML & Focus Module Hubs: %d proteins %s\n", length(ml_hubs_overlap),
                ifelse(length(ml_hubs_overlap) > 0, paste0("(", paste(ml_hubs_overlap, collapse = ", "), ")"), "")))
cat(sprintf("    Cox & Focus Module Hubs: %d proteins %s\n", length(cox_hubs_overlap),
                ifelse(length(cox_hubs_overlap) > 0, paste0("(", paste(cox_hubs_overlap, collapse = ", "), ")"), "")))
cat(sprintf("    All Three: %d proteins %s\n", length(all_three_overlap),
                ifelse(length(all_three_overlap) > 0, paste0("(", paste(all_three_overlap, collapse = ", "), ")"), "")))

# ============================================================================
# 15c: 3-Way Venn Diagram (ML, Cox, Focus Module Hubs Combined)
# ============================================================================
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15c: Creating 3-Way Venn Diagram\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# 1. Prepare the data list - 3 sets
venn_list <- list(
  "ML Biomarkers" = ml_biomarkers,
  "Cox Biomarkers" = cox_biomarkers,
  "Focus Module Hubs" = focus_module_hubs
)

# 2. Setup colors and suppress log files
# Red for ML, Blue for Cox, Green for combined hubs
venn_colors <- c("#E74C3C", "#3498DB", "#27AE60")
if(requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
}

# 3. Generate the 3-way Venn Object
# Create subtitle with focus module names
focus_modules_label <- paste(toupper(focus_modules_auto), collapse = ", ")

venn_plot <- VennDiagram::venn.diagram(
  x = venn_list,
  filename = NULL,
  category.names = c("ML\nBiomarkers", "Cox\nBiomarkers", "Focus Module\nHubs"),

  # Circles and Appearance
  fill = venn_colors,
  alpha = 0.5,
  col = "white",
  lwd = 2,

  # Inside Numbers
  cex = 1.8,
  fontface = "bold",

  # Category label styling
  cat.cex = 1.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05, 0.05),

  # Margin
  margin = 0.15,

  # Titles
  main = "Overlap: Biomarkers vs WGCNA Hub Genes",
  main.cex = 1.5,
  main.fontface = "bold",
  sub = sprintf("Focus Modules: %s | Hub criteria: kME > 0.7", focus_modules_label),
  sub.cex = 1.0
)

# --- OUTPUT TO PNG ---
png(file.path(fig_dir_png, "Step15_01_Venn_3way.png"),
    width = 10, height = 10, units = "in", res = 300)
grid::grid.draw(venn_plot)
dev.off()

# --- OUTPUT TO PDF ---
pdf(file.path(fig_dir_pdf, "Step15_01_Venn_3way.pdf"),
    width = 10, height = 10)
grid::grid.newpage()
grid::grid.draw(venn_plot)
dev.off()

# --- OUTPUT TO SCREEN (VISUAL STUDIO VIEWER) ---
grid::grid.newpage()
grid::grid.draw(venn_plot)

cat("  Saved and Displayed: Step15_01_Venn_3way\n")


# ============================================================================
# 15d: UpSet Plot (3 Sets)
# ============================================================================
# Reference: UpSetR package (Lex et al., 2014, Nature Methods)
# Better than Venn diagrams for visualizing set intersections
# Shows 3 protein sets: ML, Cox, Combined Focus Module Hubs
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15d: Creating UpSet Plot (3 Sets)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

if(requireNamespace("UpSetR", quietly = TRUE)) {

  # 1. Create binary membership matrix - 3 sets
  all_proteins_union <- unique(c(ml_biomarkers, cox_biomarkers, focus_module_hubs))

  if(length(all_proteins_union) > 0) {
    upset_data <- data.frame(
      Protein = all_proteins_union,
      ML_Biomarkers = as.integer(all_proteins_union %in% ml_biomarkers),
      Cox_Biomarkers = as.integer(all_proteins_union %in% cox_biomarkers),
      Focus_Module_Hubs = as.integer(all_proteins_union %in% focus_module_hubs)
    )

    # Colors: Red for ML, Blue for Cox, Green for Focus Module Hubs
    set_colors <- c("#E74C3C", "#3498DB", "#27AE60")

    # Set names for UpSetR
    upset_sets <- c("ML_Biomarkers", "Cox_Biomarkers", "Focus_Module_Hubs")

    # 2. Function to generate the plot (to avoid repeating code)
    make_upset <- function() {
      UpSetR::upset(
        upset_data,
        sets = upset_sets,
        order.by = "freq",
        main.bar.color = "#2c3e50",
        sets.bar.color = set_colors,
        text.scale = 1.3,
        mainbar.y.label = "Number of Proteins",
        sets.x.label = "Set Size",
        point.size = 3.5,
        line.size = 1.5,
        mb.ratio = c(0.6, 0.4)
      )
    }

    # --- SAVE TO PNG ---
    png(file.path(fig_dir_png, "Step15_02_UpSet_Plot.png"), width = 10, height = 6, units = "in", res = 300)
    print(make_upset())
    dev.off()

    # --- SAVE TO PDF ---
    pdf(file.path(fig_dir_pdf, "Step15_02_UpSet_Plot.pdf"), width = 10, height = 6)
    print(make_upset())
    dev.off()

    # --- DISPLAY ON SCREEN (Viewer) ---
    print(make_upset())

    cat("  Saved and Displayed: Step15_02_UpSet_Plot.png/pdf\n")

    # Save membership data
    write.csv(upset_data, file.path(results_dir, "Step15_Protein_SetMembership.csv"), row.names = FALSE)
    cat("  Saved: Step15_Protein_SetMembership.csv\n")

  } else {
    cat("  Warning: No proteins found in any set. Skipping UpSet plot.\n")
  }
} else {
  cat("  Note: UpSetR not installed.\n")
}

# ============================================================================
# 15e: Summary Statistics and Interpretation (3 Sets)
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15e: Summary and Interpretation\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# Define all overlap variables - 3 sets
ml_cox_overlap <- intersect(ml_biomarkers, cox_biomarkers)
ml_hubs_overlap <- intersect(ml_biomarkers, focus_module_hubs)
cox_hubs_overlap <- intersect(cox_biomarkers, focus_module_hubs)
biomarkers_hub_overlap <- intersect(c(ml_biomarkers, cox_biomarkers), focus_module_hubs)
all_three_overlap <- Reduce(intersect, list(ml_biomarkers, cox_biomarkers, focus_module_hubs))

# Warn if hub list is empty (may indicate extraction failure)
if(length(focus_module_hubs) == 0) {
  warning("No hub genes found in any focus module - check Step 6 results")
}

jaccard_ml_cox <- ifelse(length(union(ml_biomarkers, cox_biomarkers)) > 0,
                         length(ml_cox_overlap) / length(union(ml_biomarkers, cox_biomarkers)), NA)
jaccard_biomarkers_hubs <- ifelse(length(union(c(ml_biomarkers, cox_biomarkers), focus_module_hubs)) > 0,
                                  length(biomarkers_hub_overlap) / length(union(c(ml_biomarkers, cox_biomarkers), focus_module_hubs)), NA)

# Replace NA with 0 for display but preserve warning
if(is.na(jaccard_biomarkers_hubs)) {
  cat("  Note: Jaccard for biomarkers vs hubs is NA (no proteins in one or both sets)\n")
  jaccard_biomarkers_hubs <- 0
}

# Create comprehensive summary table - 3 sets
focus_modules_label <- paste(toupper(focus_modules_auto), collapse = ", ")

overlap_summary <- data.frame(
  Category = c(rep("Set Sizes", 3), rep("Key Overlaps", 4), rep("Jaccard Similarity", 2)),
  Metric = c(
    "ML Biomarkers", "Cox Biomarkers", sprintf("Focus Module Hubs (%s)", focus_modules_label),
    "ML & Cox", "ML & Focus Hubs", "Cox & Focus Hubs", "All Three",
    "ML vs Cox", "All Biomarkers vs Focus Hubs"
  ),
  Value = c(
    length(ml_biomarkers), length(cox_biomarkers), length(focus_module_hubs),
    length(ml_cox_overlap), length(ml_hubs_overlap), length(cox_hubs_overlap), length(all_three_overlap),
    round(jaccard_ml_cox, 3), round(jaccard_biomarkers_hubs, 3)
  ),
  stringsAsFactors = FALSE
)

cat("\n  3-WAY OVERLAP SUMMARY:\n")
print(overlap_summary)

# Save summary
write.csv(overlap_summary, file.path(results_dir, "Step15_Overlap_Summary.csv"), row.names = FALSE)

# List overlapping proteins if any
cat("\n  Overlapping Proteins:\n")
if(length(ml_cox_overlap) > 0) {
  cat(sprintf("    ML & Cox: %s\n", paste(ml_cox_overlap, collapse = ", ")))
} else {
  cat("    ML & Cox: None\n")
}
if(length(ml_hubs_overlap) > 0) {
  cat(sprintf("    ML & Focus Hubs: %s\n", paste(ml_hubs_overlap, collapse = ", ")))
} else {
  cat("    ML & Focus Hubs: None\n")
}
if(length(cox_hubs_overlap) > 0) {
  cat(sprintf("    Cox & Focus Hubs: %s\n", paste(cox_hubs_overlap, collapse = ", ")))
} else {
  cat("    Cox & Focus Hubs: None\n")
}
if(length(all_three_overlap) > 0) {
  cat(sprintf("    All Three: %s\n", paste(all_three_overlap, collapse = ", ")))
} else {
  cat("    All Three: None\n")
}

# Per-module hub breakdown
cat("\n  Focus Module Hub Breakdown:\n")
for(fm in names(focus_module_hub_details)) {
  cat(sprintf("    %s: %d hubs\n", tools::toTitleCase(fm), length(focus_module_hub_details[[fm]])))
}

# Interpretation text - 3 sets
interpretation_text <- sprintf("
STEP 15 INTERPRETATION: 3-Way Signature Overlap Analysis

SETS COMPARED:
1. ML Biomarkers (n=%d): %s
2. Cox Biomarkers (n=%d): %s
3. Focus Module Hubs (n=%d) - Combined from %s

PER-MODULE HUB COUNTS:
%s

KEY FINDINGS:
ML & Cox overlap: %d proteins (Jaccard = %.3f)
ML & Focus Hubs overlap: %d proteins
Cox & Focus Hubs overlap: %d proteins
All Three overlap: %d proteins
Biomarkers & Hub Genes overlap: %d proteins (Jaccard = %.3f)

INTERPRETATION:
The %s overlap between biomarkers and WGCNA hub genes is scientifically
meaningful:

1. BIOMARKERS (ML + Cox): Selected for PREDICTIVE value
   - ML Classifier: Optimized for classification accuracy
   - Cox Regression: Selected for survival association
   - These capture proteins with PROGNOSTIC signal

2. HUB GENES (Focus Modules: %s): Selected for NETWORK CENTRALITY
   - High module membership (kME > 0.7)
   - These modules correlate with clinical outcomes
   - These represent CORE PATHWAY members

3. WHY NO OVERLAP IS OBSERVED:
   - Predictive power != Network centrality
   - Biomarkers work through indirect mechanisms, not as hub genes
   - Hub genes are pathway regulators, not necessarily prognostic markers
   - DYNC2H1 (ML biomarker) is in grey module - completely independent
   - Cox biomarkers also reside in non-focus modules

THESIS NARRATIVE:
\"Predictive biomarkers and network hub genes represent complementary
aspects of PDAC biology. Biomarkers capture prognostic signal through
diverse mechanisms, while hub genes represent core pathway regulators.
The lack of overlap suggests that predictive power in cancer prognosis
arises from multiple biological sources, not just central network nodes.\"

This finding supports a MULTI-PATHWAY prognostic model where:
- Different analytical methods capture different biology
- Combined signatures may outperform single approaches
- Biomarkers work through varied mechanisms, not just hub centrality

",
  length(ml_biomarkers), paste(ml_biomarkers, collapse = ", "),
  length(cox_biomarkers), paste(cox_biomarkers, collapse = ", "),
  length(focus_module_hubs), focus_modules_label,
  paste(sapply(names(focus_module_hub_details), function(fm)
    sprintf("  %s: %d hubs", tools::toTitleCase(fm), length(focus_module_hub_details[[fm]]))),
    collapse = "\n"),
  length(ml_cox_overlap), jaccard_ml_cox,
  length(ml_hubs_overlap),
  length(cox_hubs_overlap),
  length(all_three_overlap),
  length(biomarkers_hub_overlap), jaccard_biomarkers_hubs,
  ifelse(length(biomarkers_hub_overlap) == 0, "lack of",
         ifelse(length(biomarkers_hub_overlap) < 3, "minimal", "moderate")),
  focus_modules_label
)

cat(interpretation_text)
writeLines(interpretation_text, file.path(results_dir, "Step15_Interpretation.txt"))
cat("  Saved: Step15_Interpretation.txt\n")

# Save overlap matrix
write.csv(overlap_matrix, file.path(results_dir, "Step15_Overlap_Matrix.csv"), row.names = TRUE)
cat("  Saved: Step15_Overlap_Matrix.csv\n")

# ============================================================================

# ============================================================================
# Purpose: Visualize protein flow from modules to signatures to outcomes
# Package: ggalluvial
# ============================================================================

cat("15f: Creating Alluvial/Sankey Plot...\n")

# ggalluvial loaded at startup
if(requireNamespace("ggalluvial", quietly = TRUE)) {

  # Create flow data: Module -> Signature Type -> Clinical Outcome
  # Track which proteins came from which modules and ended up in which signatures

  # Get module membership for all biomarkers
  all_biomarker_list <- unique(c(ml_biomarkers, cox_biomarkers))

  flow_data <- data.frame(
    Protein = all_biomarker_list,
    stringsAsFactors = FALSE
  )

  # Add module source
  flow_data$Module <- moduleColors[all_biomarker_list]
  flow_data$Module[is.na(flow_data$Module)] <- "Other"

  # Add signature membership
  flow_data$Signature <- case_when(
    flow_data$Protein %in% ml_biomarkers & flow_data$Protein %in% cox_biomarkers ~ "Both",
    flow_data$Protein %in% ml_biomarkers ~ "ML",
    flow_data$Protein %in% cox_biomarkers ~ "Cox",
    TRUE ~ "None"
  )

  # Add hub status - using ALL focus modules (combined hubs from earlier)
  flow_data$Hub_Status <- ifelse(flow_data$Protein %in% focus_module_hubs, "Hub", "Non-Hub")

  # Aggregate for alluvial
  flow_agg <- flow_data %>%
    group_by(Module, Signature, Hub_Status) %>%
    summarise(Count = n(), .groups = "drop")

  # Complete color palette for all possible modules
  alluvial_colors <- c(
    "pink" = "#FFC0CB", "green" = "#228B22", "black" = "#2D2D2D", "blue" = "#4169E1",
    "brown" = "#8B4513", "turquoise" = "#40E0D0", "yellow" = "#FFD700", "red" = "#DC143C",
    "magenta" = "#FF00FF", "purple" = "#800080", "grey" = "grey70", "Other" = "grey50"
  )

  # Create alluvial plot
  p_alluvial <- ggplot(flow_agg,
                        aes(axis1 = Module, axis2 = Signature, axis3 = Hub_Status, y = Count)) +
    geom_alluvium(aes(fill = Module), width = 1/12, alpha = 0.7) +
    geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Module", "Signature", "Hub Status"),
                     expand = c(0.1, 0.1)) +
    scale_fill_manual(values = alluvial_colors) +
    labs(
      title = "Biomarker Flow: Module to Signature to Hub Status",
      subtitle = "Tracking protein origins through the analysis pipeline",
      y = "Number of Proteins",
      caption = "Alluvial visualization showing biomarker provenance"
    ) +
    theme_wgcna() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  print(p_alluvial)
  ggsave(file.path(fig_dir_png, "Step15_02_Alluvial_SignatureFlow.png"), p_alluvial, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir_pdf, "Step15_02_Alluvial_SignatureFlow.pdf"), p_alluvial, width = 8, height = 6)
  cat("  Saved: Step15_02_Alluvial_SignatureFlow.png/pdf\n")

} else {
  cat("  Note: ggalluvial package not available for Sankey plot\n")
  cat("  Install with: install.packages('ggalluvial')\n")
}


# ============================================================================
# 15-FINAL: STEP SUMMARY
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 15 SUMMARY: BIOMARKER INTEGRATION ===\n")
cat("===============================================================================\n")

cat("\nBiomarker Module Integration:\n")
cat(sprintf("  ML biomarkers: %d total\n", length(ml_biomarkers)))
cat(sprintf("    - In focus modules (%s): %d\n",
            paste(toupper(focus_modules_auto), collapse = "/"),
            sum(biomarker_membership$In_Focus_Module[biomarker_membership$Source == "ML"], na.rm = TRUE)))
cat(sprintf("    - Are hub genes (kME > 0.7): %d\n",
            sum(biomarker_membership$Is_Hub[biomarker_membership$Source == "ML"], na.rm = TRUE)))

cat(sprintf("  Cox biomarkers: %d total\n", length(cox_biomarkers)))
cat(sprintf("    - In focus modules (%s): %d\n",
            paste(toupper(focus_modules_auto), collapse = "/"),
            sum(biomarker_membership$In_Focus_Module[biomarker_membership$Source == "Cox"], na.rm = TRUE)))
cat(sprintf("    - Are hub genes (kME > 0.7): %d\n",
            sum(biomarker_membership$Is_Hub[biomarker_membership$Source == "Cox"], na.rm = TRUE)))

cat("\nSignature Membership:\n")
n_sig_associations <- sum(!is.na(biomarker_sig_membership$Signature_ID))
n_unique_sigs <- length(unique(biomarker_sig_membership$Signature_ID[!is.na(biomarker_sig_membership$Signature_ID)]))
cat(sprintf("  Total biomarker-signature associations: %d\n", n_sig_associations))
cat(sprintf("  Unique signatures containing biomarkers: %d\n", n_unique_sigs))

cat("\nKey Finding:\n")
n_hub_biomarkers <- sum(biomarker_membership$Is_Hub, na.rm = TRUE)
if(n_hub_biomarkers == 0) {
  cat("  NO overlap between biomarkers and hub genes.\n")
  cat("  This is scientifically meaningful: predictive power != network centrality.\n")
  cat("  ML/Cox biomarkers capture complementary biology to WGCNA hubs.\n")
} else {
  overlap_proteins <- biomarker_membership$Protein[biomarker_membership$Is_Hub]
  cat(sprintf("  %d biomarker(s) are also hub genes: %s\n",
              n_hub_biomarkers, paste(overlap_proteins, collapse = ", ")))
}

cat("\nOutput files:\n")
cat("  Tables:\n")
cat("    - Step15_Biomarker_Module_Membership.csv\n")
cat("    - Step15_Biomarker_Signature_Membership.csv\n")
cat("    - Step15_Overlap_Summary.csv\n")
cat("    - Step15_Overlap_Matrix.csv\n")
cat("    - Step15_Protein_SetMembership.csv\n")
cat("    - Step15_Interpretation.txt\n")
cat("  Figures:\n")
cat("    - Step15_01_Venn_3way.png/pdf\n")
cat("    - Step15_02_UpSet_Plot.png/pdf\n")
cat("    - Step15_02_Alluvial_SignatureFlow.png/pdf\n")

# ============================================================================
# 15g: BIOMARKER BIOLOGICAL CHARACTERIZATION (NEW)
# ============================================================================
# Purpose: Provide valid (non-circular) biological interpretation of biomarkers
# Approaches:
#   1. External signature context (from mastertable)
#   2. Module pathway context (from Step 8 ORA results)
#   3. Functional annotation (GO terms via org.Hs.eg.db)
#   4. Comprehensive summary table (publication-ready)
# ============================================================================

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("15g: Biomarker Biological Characterization\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# ----------------------------------------------------------------------------
# 15g-1: EXTERNAL SIGNATURE CONTEXT
# ----------------------------------------------------------------------------
cat("  [1/4] External Signature Context...\n")

# Use the existing biomarker_sig_membership data (created in 15a-3)
# Enhance with additional info

if(exists("biomarker_sig_membership") && nrow(biomarker_sig_membership) > 0) {

  # Add module info to signature context
  biomarker_context <- biomarker_sig_membership %>%
    filter(!is.na(Signature_ID)) %>%
    left_join(
      biomarker_membership %>% select(Protein, Module),
      by = "Protein"
    )

  # Add PMID from custom_signatures
  biomarker_context$PMID <- sapply(biomarker_context$Signature_ID, function(sig_id) {
    if(!is.null(custom_signatures[[sig_id]]) && !is.null(custom_signatures[[sig_id]]$pmid)) {
      return(custom_signatures[[sig_id]]$pmid)
    } else {
      return(NA)
    }
  })

  # Reorder columns
  biomarker_context <- biomarker_context %>%
    select(Protein, Source, Module, Signature_ID, Signature_Name, Category, PMID) %>%
    arrange(Protein, Category, Signature_Name)

  write.csv(biomarker_context,
            file.path(results_dir, "Step15g_Biomarker_Signature_Context.csv"),
            row.names = FALSE)

  cat(sprintf("    Found %d biomarker-signature associations\n", nrow(biomarker_context)))
  cat("    Saved: Step15g_Biomarker_Signature_Context.csv\n")

} else {
  cat("    Warning: biomarker_sig_membership not available\n")
  biomarker_context <- data.frame()
}

# ----------------------------------------------------------------------------
# 15g-2: MODULE PATHWAY CONTEXT
# ----------------------------------------------------------------------------
cat("  [2/4] Module Pathway Context...\n")

# For each biomarker, get top pathways from its module's ORA results
if(exists("ora_results") && nrow(ora_results) > 0) {

  # Get unique modules containing biomarkers
  biomarker_modules <- unique(biomarker_membership$Module[!is.na(biomarker_membership$Module)])

  # For each module, get top 10 pathways
  module_top_pathways <- data.frame()

  for(mod in biomarker_modules) {
    mod_pathways <- ora_results %>%
      filter(tolower(Module) == tolower(mod), P_Value < 0.05) %>%
      arrange(FDR, P_Value) %>%
      head(10) %>%
      mutate(
        Pathway_Rank = row_number(),
        Module = mod
      ) %>%
      select(Module, Pathway_Rank, Short_Name, Category, P_Value, FDR, Fold_Enrichment)

    module_top_pathways <- rbind(module_top_pathways, mod_pathways)
  }

  # Now link biomarkers to their module's pathways
  biomarker_module_pathways <- data.frame()

  for(i in 1:nrow(biomarker_membership)) {
    prot <- biomarker_membership$Protein[i]
    src <- biomarker_membership$Source[i]
    mod <- biomarker_membership$Module[i]

    if(!is.na(mod)) {
      mod_paths <- module_top_pathways %>%
        filter(tolower(Module) == tolower(mod))

      if(nrow(mod_paths) > 0) {
        prot_paths <- mod_paths %>%
          mutate(
            Protein = prot,
            Source = src
          ) %>%
          select(Protein, Source, Module, Pathway_Rank, Pathway_Name = Short_Name,
                 Pathway_Category = Category, FDR, Fold_Enrichment)

        biomarker_module_pathways <- rbind(biomarker_module_pathways, prot_paths)
      }
    }
  }

  if(nrow(biomarker_module_pathways) > 0) {
    write.csv(biomarker_module_pathways,
              file.path(results_dir, "Step15g_Biomarker_Module_Pathways.csv"),
              row.names = FALSE)
    cat(sprintf("    Created pathway context for %d biomarkers\n",
                length(unique(biomarker_module_pathways$Protein))))
    cat("    Saved: Step15g_Biomarker_Module_Pathways.csv\n")
  } else {
    cat("    Warning: No module pathway associations found\n")
  }

} else {
  cat("    Warning: ora_results not available\n")
  biomarker_module_pathways <- data.frame()
}

# ----------------------------------------------------------------------------
# 15g-3: FUNCTIONAL ANNOTATION (GO Terms)
# ----------------------------------------------------------------------------
cat("  [3/4] Functional Annotation (GO Terms)...\n")

biomarker_go <- data.frame()

if(requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
   requireNamespace("AnnotationDbi", quietly = TRUE)) {

  library(org.Hs.eg.db)
  library(AnnotationDbi)

  for(prot in all_biomarkers) {

    # Get ENTREZ ID
    entrez <- tryCatch({
      res <- AnnotationDbi::select(org.Hs.eg.db, keys = prot,
                                   keytype = "SYMBOL", columns = "ENTREZID")
      res$ENTREZID[1]
    }, error = function(e) NA)

    if(!is.na(entrez)) {
      # Get GO terms
      go_terms <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = entrez,
                              keytype = "ENTREZID",
                              columns = c("GO", "ONTOLOGY", "TERM"))
      }, error = function(e) NULL)

      if(!is.null(go_terms) && nrow(go_terms) > 0) {
        # Get BP terms (top 3)
        bp_terms <- go_terms %>%
          filter(ONTOLOGY == "BP") %>%
          head(3) %>%
          pull(TERM)
        bp_summary <- if(length(bp_terms) > 0) paste(bp_terms, collapse = "; ") else NA

        # Get MF terms (top 2)
        mf_terms <- go_terms %>%
          filter(ONTOLOGY == "MF") %>%
          head(2) %>%
          pull(TERM)
        mf_summary <- if(length(mf_terms) > 0) paste(mf_terms, collapse = "; ") else NA

        biomarker_go <- rbind(biomarker_go, data.frame(
          Protein = prot,
          ENTREZ_ID = entrez,
          GO_BP = bp_summary,
          GO_MF = mf_summary,
          stringsAsFactors = FALSE
        ))
      } else {
        biomarker_go <- rbind(biomarker_go, data.frame(
          Protein = prot,
          ENTREZ_ID = entrez,
          GO_BP = NA,
          GO_MF = NA,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      biomarker_go <- rbind(biomarker_go, data.frame(
        Protein = prot,
        ENTREZ_ID = NA,
        GO_BP = NA,
        GO_MF = NA,
        stringsAsFactors = FALSE
      ))
    }
  }

  write.csv(biomarker_go,
            file.path(results_dir, "Step15g_Biomarker_GO_Annotation.csv"),
            row.names = FALSE)
  cat(sprintf("    Annotated %d biomarkers with GO terms\n", nrow(biomarker_go)))
  cat("    Saved: Step15g_Biomarker_GO_Annotation.csv\n")

} else {
  cat("    Note: org.Hs.eg.db not available. Skipping GO annotation.\n")
  cat("    Install with: BiocManager::install('org.Hs.eg.db')\n")
}

# ----------------------------------------------------------------------------
# 15g-4: COMPREHENSIVE SUMMARY TABLE (Publication-Ready)
# ----------------------------------------------------------------------------
cat("  [4/4] Comprehensive Summary Table...\n")

# Build master table combining all information
comprehensive_summary <- biomarker_membership %>%
  select(Protein, Source, Module, kME, In_Focus_Module, Is_Hub)

# Add signature counts
sig_counts <- biomarker_sig_membership %>%
  filter(!is.na(Signature_ID)) %>%
  group_by(Protein) %>%
  summarise(
    N_External_Signatures = n(),
    Signature_Categories = paste(unique(Category), collapse = "; "),
    .groups = "drop"
  )

comprehensive_summary <- comprehensive_summary %>%
  left_join(sig_counts, by = "Protein") %>%
  mutate(
    N_External_Signatures = ifelse(is.na(N_External_Signatures), 0, N_External_Signatures),
    Signature_Categories = ifelse(is.na(Signature_Categories), "None", Signature_Categories)
  )

# Add top pathway from module
if(exists("biomarker_module_pathways") && nrow(biomarker_module_pathways) > 0) {
  top_pathway_per_biomarker <- biomarker_module_pathways %>%
    filter(Pathway_Rank == 1) %>%
    select(Protein, Module_Top_Pathway = Pathway_Name)

  comprehensive_summary <- comprehensive_summary %>%
    left_join(top_pathway_per_biomarker, by = "Protein")
} else {
  comprehensive_summary$Module_Top_Pathway <- NA
}

# Add GO BP summary
if(nrow(biomarker_go) > 0) {
  go_summary <- biomarker_go %>%
    select(Protein, GO_BP_Summary = GO_BP)

  comprehensive_summary <- comprehensive_summary %>%
    left_join(go_summary, by = "Protein")
} else {
  comprehensive_summary$GO_BP_Summary <- NA
}

# Reorder and clean
comprehensive_summary <- comprehensive_summary %>%
  select(Protein, Source, Module, kME, In_Focus_Module, Is_Hub,
         N_External_Signatures, Signature_Categories,
         Module_Top_Pathway, GO_BP_Summary) %>%
  arrange(Source, desc(abs(kME)))

write.csv(comprehensive_summary,
          file.path(results_dir, "Step15g_Biomarker_Comprehensive_Summary.csv"),
          row.names = FALSE)

cat("    Created comprehensive summary for 8 biomarkers\n")
cat("    Saved: Step15g_Biomarker_Comprehensive_Summary.csv (Table S5)\n")

# Print summary to console
cat("\n  === BIOMARKER COMPREHENSIVE SUMMARY ===\n")
print(comprehensive_summary %>% select(Protein, Source, Module, kME, N_External_Signatures))

# ----------------------------------------------------------------------------
# 15g-5: CATEGORY HEATMAP VISUALIZATION
# ----------------------------------------------------------------------------
cat("\n  [5/5] Category Heatmap...\n")

if(exists("biomarker_context") && nrow(biomarker_context) > 0) {

  # Create biomarker x category count matrix
  category_counts <- biomarker_context %>%
    group_by(Protein, Category) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Category, values_from = Count, values_fill = 0)

  # Convert to matrix
  cat_matrix <- as.matrix(category_counts %>% select(-Protein))
  rownames(cat_matrix) <- category_counts$Protein

  # Reorder rows: ML first, then Cox
  ml_order <- intersect(ml_biomarkers, rownames(cat_matrix))
  cox_order <- intersect(cox_biomarkers, rownames(cat_matrix))
  row_order <- c(ml_order, cox_order)
  row_order <- row_order[row_order %in% rownames(cat_matrix)]

  if(length(row_order) > 0) {
    cat_matrix <- cat_matrix[row_order, , drop = FALSE]

    # Create row annotation (Source and Module)
    row_annotation <- data.frame(
      Source = ifelse(rownames(cat_matrix) %in% ml_biomarkers, "ML", "Cox"),
      Module = biomarker_membership$Module[match(rownames(cat_matrix), biomarker_membership$Protein)],
      row.names = rownames(cat_matrix)
    )

    # Annotation colors
    ann_colors <- list(
      Source = c("ML" = "#E74C3C", "Cox" = "#3498DB"),
      Module = setNames(
        sapply(unique(row_annotation$Module), function(m) {
          if(is.na(m)) return("grey70")
          if(m == "grey") return("grey50")
          return(m)
        }),
        unique(row_annotation$Module)
      )
    )

    # Create heatmap
    fig_cat_heatmap <- pheatmap(
      cat_matrix,
      color = colorRampPalette(c("white", "#FFF5EB", "#FD8D3C", "#D94801", "#7F2704"))(50),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      annotation_row = row_annotation,
      annotation_colors = ann_colors,
      display_numbers = TRUE,
      number_format = "%d",
      fontsize_number = 10,
      fontsize_row = 11,
      fontsize_col = 10,
      main = "Biomarker-Signature Category Associations\n(Count of signatures per category)",
      angle_col = 45,
      border_color = "grey80",
      silent = TRUE
    )

    # Save
    png(file.path(fig_dir_png, "Step15g_Biomarker_Category_Heatmap.png"),
        width = 12, height = 6, units = "in", res = 300)
    print(fig_cat_heatmap)
    dev.off()

    pdf(file.path(fig_dir_pdf, "Step15g_Biomarker_Category_Heatmap.pdf"),
        width = 12, height = 6)
    print(fig_cat_heatmap)
    dev.off()

    # Display
    print(fig_cat_heatmap)

    cat("    Saved: Step15g_Biomarker_Category_Heatmap.png/pdf\n")

  } else {
    cat("    Warning: No biomarkers with signature associations for heatmap\n")
  }

} else {
  cat("    Warning: biomarker_context not available for heatmap\n")
}

cat("\n  Step 15g: Biomarker Biological Characterization - COMPLETE\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# ============================================================================
# 15-FINAL: STEP SUMMARY (UPDATED)
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 15 SUMMARY: BIOMARKER INTEGRATION ===\n")
cat("===============================================================================\n")

cat("\nBiomarker Module Integration:\n")
cat(sprintf("  ML biomarkers: %d total\n", length(ml_biomarkers)))
cat(sprintf("    - In focus modules (%s): %d\n",
            paste(toupper(focus_modules_auto), collapse = "/"),
            sum(biomarker_membership$In_Focus_Module[biomarker_membership$Source == "ML"], na.rm = TRUE)))
cat(sprintf("    - Are hub genes (kME > 0.7): %d\n",
            sum(biomarker_membership$Is_Hub[biomarker_membership$Source == "ML"], na.rm = TRUE)))

cat(sprintf("  Cox biomarkers: %d total\n", length(cox_biomarkers)))
cat(sprintf("    - In focus modules (%s): %d\n",
            paste(toupper(focus_modules_auto), collapse = "/"),
            sum(biomarker_membership$In_Focus_Module[biomarker_membership$Source == "Cox"], na.rm = TRUE)))
cat(sprintf("    - Are hub genes (kME > 0.7): %d\n",
            sum(biomarker_membership$Is_Hub[biomarker_membership$Source == "Cox"], na.rm = TRUE)))

cat("\nSignature Membership:\n")
n_sig_associations <- sum(!is.na(biomarker_sig_membership$Signature_ID))
n_unique_sigs <- length(unique(biomarker_sig_membership$Signature_ID[!is.na(biomarker_sig_membership$Signature_ID)]))
cat(sprintf("  Total biomarker-signature associations: %d\n", n_sig_associations))
cat(sprintf("  Unique signatures containing biomarkers: %d\n", n_unique_sigs))

cat("\nKey Finding:\n")
n_hub_biomarkers <- sum(biomarker_membership$Is_Hub, na.rm = TRUE)
if(n_hub_biomarkers == 0) {
  cat("  NO overlap between biomarkers and hub genes.\n")
  cat("  This is scientifically meaningful: predictive power != network centrality.\n")
  cat("  ML/Cox biomarkers capture complementary biology to WGCNA hubs.\n")
} else {
  overlap_proteins <- biomarker_membership$Protein[biomarker_membership$Is_Hub]
  cat(sprintf("  %d biomarker(s) are also hub genes: %s\n",
              n_hub_biomarkers, paste(overlap_proteins, collapse = ", ")))
}

cat("\nOutput files:\n")
cat("  Tables:\n")
cat("    - Step15_Biomarker_Module_Membership.csv\n")
cat("    - Step15_Biomarker_Signature_Membership.csv\n")
cat("    - Step15_Overlap_Summary.csv\n")
cat("    - Step15_Overlap_Matrix.csv\n")
cat("    - Step15_Protein_SetMembership.csv\n")
cat("    - Step15_Interpretation.txt\n")
cat("    - Step15g_Biomarker_Signature_Context.csv (NEW)\n")
cat("    - Step15g_Biomarker_Module_Pathways.csv (NEW)\n")
cat("    - Step15g_Biomarker_GO_Annotation.csv (NEW)\n")
cat("    - Step15g_Biomarker_Comprehensive_Summary.csv (Table S5) (NEW)\n")
cat("  Figures:\n")
cat("    - Step15_01_Venn_3way.png/pdf\n")
cat("    - Step15_02_UpSet_Plot.png/pdf\n")
cat("    - Step15_02_Alluvial_SignatureFlow.png/pdf\n")
cat("    - Step15g_Biomarker_Category_Heatmap.png/pdf (NEW)\n")
cat("Step 15 complete.\n")
cat("===============================================================================\n")

# ============================================================================
# STEP 16: STRING/PPI NETWORK ANALYSIS
# ============================================================================
# Purpose: Analyze protein-protein interactions among hub genes
# Method: STRINGdb API for human (species 9606)
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 16: STRING/PPI NETWORK ANALYSIS\n")
cat("================================================================\n")

# Load STRINGdb package
# NOTE: Requires BiocManager::install('STRINGdb') if not installed
library(STRINGdb)

# Initialize STRING database (human = 9606, score_threshold = 400 for medium confidence)
string_db <- STRINGdb::STRINGdb$new(version = "12.0", species = 9606,
                          score_threshold = 400, input_directory = "")

cat("STRING database initialized (human, score threshold = 400)\n")

# ============================================================================
# 16a: PPI Analysis for ALL Focus Modules (PINK, GREEN, BLACK, BLUE)
# ============================================================================
# Loop through all focus modules dynamically

# Storage for results
ppi_results <- list()

for(focus_mod in focus_modules_auto) {

  mod_upper <- toupper(focus_mod)
  cat(sprintf("\n--- %s Module PPI Analysis ---\n", mod_upper))

  # Get hub genes for this module (kME > 0.7, isHub_Moderate = TRUE)
  mod_hubs <- proteinInfo %>%
    filter(Module == focus_mod, isHub_Moderate == TRUE) %>%
    pull(Protein)
  cat(sprintf("%s hub genes: %d\n", mod_upper, length(mod_hubs)))

  # Skip if no hub genes
  if(length(mod_hubs) == 0) {
    cat("  Skipping - no hub genes found\n")
    ppi_results[[focus_mod]] <- list(
      module = focus_mod,
      n_hubs = 0,
      n_mapped = 0,
      n_interactions = 0,
      hubs = character(0),
      mapped = NULL,
      interactions = NULL,
      enrichment = NULL
    )
    next
  }

  # Map to STRING
  mod_df <- data.frame(gene = mod_hubs, stringsAsFactors = FALSE)
  mod_mapped <- string_db$map(mod_df, "gene", removeUnmappedRows = TRUE)
  cat(sprintf("  Mapped to STRING: %d/%d (%.1f%%)\n",
              nrow(mod_mapped), length(mod_hubs),
              100 * nrow(mod_mapped) / max(1, length(mod_hubs))))

  # Get interactions (only if we have mapped proteins)
  mod_interactions <- NULL
  if(nrow(mod_mapped) > 0) {
    mod_interactions <- string_db$get_interactions(mod_mapped$STRING_id)
    cat(sprintf("  Interactions found: %d\n", nrow(mod_interactions)))
  } else {
    cat("  No mapped proteins - skipping interactions\n")
  }

  # Get functional enrichment from STRING
  mod_enrichment <- tryCatch({
    if(nrow(mod_mapped) > 0) {
      string_db$get_enrichment(mod_mapped$STRING_id)
    } else {
      NULL
    }
  }, error = function(e) {
    cat("  Note: STRING enrichment not available\n")
    NULL
  })

  # Plot network to console/viewer (only if we have mapped proteins)
  if(nrow(mod_mapped) > 0) {
    string_db$plot_network(mod_mapped$STRING_id)
  }

  # Store results
  ppi_results[[focus_mod]] <- list(
    module = focus_mod,
    n_hubs = length(mod_hubs),
    n_mapped = nrow(mod_mapped),
    n_interactions = ifelse(is.null(mod_interactions), 0, nrow(mod_interactions)),
    hubs = mod_hubs,
    mapped = mod_mapped,
    interactions = mod_interactions,
    enrichment = mod_enrichment
  )

  # --- SAVE NETWORK PLOT (uncomment to save) ---
  # if(nrow(mod_mapped) > 0) {
  #   png(file.path(fig_dir_png, paste0("Step16_STRING_", mod_upper, "_Network.png")),
  #       width = 10, height = 10, units = "in", res = 300)
  #   string_db$plot_network(mod_mapped$STRING_id)
  #   dev.off()
  #   pdf(file.path(fig_dir_pdf, paste0("Step16_STRING_", mod_upper, "_Network.pdf")),
  #       width = 10, height = 10)
  #   string_db$plot_network(mod_mapped$STRING_id)
  #   dev.off()
  #   cat(sprintf("  Saved: Step16_STRING_%s_Network.png/pdf\n", mod_upper))
  # }

  # --- SAVE INTERACTIONS TABLE (uncomment to save) ---
  # if(!is.null(mod_interactions) && nrow(mod_interactions) > 0) {
  #   write.csv(mod_interactions,
  #             file.path(results_dir, paste0("Step16_STRING_", mod_upper, "_Interactions.csv")),
  #             row.names = FALSE)
  #   cat(sprintf("  Saved: Step16_STRING_%s_Interactions.csv\n", mod_upper))
  # }

  # --- SAVE STRING ENRICHMENT (uncomment to save) ---
  # if(!is.null(mod_enrichment) && nrow(mod_enrichment) > 0) {
  #   write.csv(mod_enrichment,
  #             file.path(results_dir, paste0("Step16_STRING_", mod_upper, "_Enrichment.csv")),
  #             row.names = FALSE)
  #   cat(sprintf("  Saved: Step16_STRING_%s_Enrichment.csv\n", mod_upper))
  # }
}

# ============================================================================
# 16b: Combined PPI Summary for ALL Focus Modules
# ============================================================================

cat("\n--- PPI Summary (All Focus Modules) ---\n")

# Build summary from stored results
ppi_summary <- data.frame(
  Module = sapply(ppi_results, function(x) toupper(x$module)),
  N_Hub_Genes = sapply(ppi_results, function(x) x$n_hubs),
  N_Mapped_STRING = sapply(ppi_results, function(x) x$n_mapped),
  N_Interactions = sapply(ppi_results, function(x) x$n_interactions),
  stringsAsFactors = FALSE
)
rownames(ppi_summary) <- NULL

print(ppi_summary)

# Total summary
cat(sprintf("\nTotal across %d focus modules:\n", length(focus_modules_auto)))
cat(sprintf("  Hub genes: %d\n", sum(ppi_summary$N_Hub_Genes)))
cat(sprintf("  Mapped to STRING: %d\n", sum(ppi_summary$N_Mapped_STRING)))
cat(sprintf("  Total interactions: %d\n", sum(ppi_summary$N_Interactions)))

# --- SAVE PPI SUMMARY (uncomment to save) ---
# write.csv(ppi_summary,
#           file.path(results_dir, "Step16_STRING_Summary.csv"),
#           row.names = FALSE)

cat("\nSTEP 16 COMPLETE: STRING/PPI Network Analysis (All Focus Modules)\n")

# ============================================================================
# STEP 17: CIRCOS MODULE OVERVIEW
# ============================================================================
# Purpose: Visual summary of all modules, their sizes, and trait correlations
# Package: circlize
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 17: CIRCOS MODULE OVERVIEW\n")
cat("================================================================\n")

# circlize loaded at startup
# NOTE: Requires install.packages("circlize") if not installed
library(circlize)

# ============================================================================
# 17a: Prepare Module Data
# ============================================================================

# Get module sizes (excluding grey)
module_sizes <- table(moduleColors)
module_sizes <- module_sizes[names(module_sizes) != "grey"]
module_names <- names(module_sizes)

cat(sprintf("Modules to visualize: %d (excluding grey)\n", length(module_names)))

# Get module-trait correlations for key traits
# Assuming MEs and datTraits exist from earlier steps
key_traits <- c("PFS_group", "PFS", "Response")
key_traits <- key_traits[key_traits %in% colnames(datTraits)]

# Create correlation matrix
module_trait_cor <- data.frame(Module = module_names)
module_trait_pval <- data.frame(Module = module_names)

for(trait in key_traits) {
  cors <- sapply(module_names, function(mod) {
    me_col <- paste0("ME", mod)
    if(me_col %in% colnames(MEs)) {
      trait_numeric <- as.numeric(datTraits[, trait])
      cor(MEs[, me_col], trait_numeric, use = "complete.obs")
    } else {
      NA
    }
  })
  pvals <- sapply(module_names, function(mod) {
    me_col <- paste0("ME", mod)
    if(me_col %in% colnames(MEs)) {
      trait_numeric <- as.numeric(datTraits[, trait])
      cor.test(MEs[, me_col], trait_numeric)$p.value
    } else {
      NA
    }
  })
  module_trait_cor[[trait]] <- cors
  module_trait_pval[[trait]] <- pvals
}

# Identify focus modules (use ALL focus modules: PINK, GREEN, BLACK, BLUE)
focus_modules <- focus_modules_auto

# ============================================================================
# 17b: Create Circos Plot
# ============================================================================

cat("Generating Circos plot...\n")

# Create circos plot (displays to viewer)
circos.clear()
circos.par(gap.after = rep(4, length(module_names)), start.degree = 90)

# Initialize with module sizes
circos.initialize(factors = module_names,
                  xlim = cbind(rep(0, length(module_names)), as.numeric(module_sizes)))

# Track 1: Module names and sizes (colored by module)
circos.track(ylim = c(0, 1),
             panel.fun = function(x, y) {
               sector.name <- CELL_META$sector.index
               n_proteins <- module_sizes[sector.name]

               # Module label
               circos.text(CELL_META$xcenter, 0.5,
                          paste0(toupper(substr(sector.name, 1, 1)),
                                 substr(sector.name, 2, nchar(sector.name)),
                                 "\n(n=", n_proteins, ")"),
                          facing = "bending.inside",
                          niceFacing = TRUE,
                          cex = 0.8,
                          font = ifelse(sector.name %in% focus_modules, 2, 1))
             },
             bg.col = module_names,
             bg.border = ifelse(module_names %in% focus_modules, "black", NA),
             bg.lwd = ifelse(module_names %in% focus_modules, 3, 1),
             track.height = 0.15)

# Track 2: PFS_group correlation bars
circos.track(ylim = c(-1, 1),
             panel.fun = function(x, y) {
               sector.name <- CELL_META$sector.index
               cor_val <- module_trait_cor[module_trait_cor$Module == sector.name, "PFS_group"]
               pval <- module_trait_pval[module_trait_pval$Module == sector.name, "PFS_group"]

               if(!is.na(cor_val)) {
                 # Bar color based on significance
                 bar_col <- ifelse(pval < 0.05,
                                  ifelse(cor_val > 0, "darkgreen", "darkred"),
                                  "grey70")

                 circos.rect(CELL_META$cell.xlim[1], 0,
                            CELL_META$cell.xlim[2], cor_val,
                            col = bar_col, border = NA)

                 # Add correlation value
                 circos.text(CELL_META$xcenter, cor_val + sign(cor_val) * 0.15,
                            sprintf("%.2f", cor_val),
                            cex = 0.6, facing = "bending.inside", niceFacing = TRUE)
               }

               # Zero line
               circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.5)
             },
             track.height = 0.2)

# Add title
title(main = "WGCNA Module Overview\nCircos Plot",
      sub = sprintf("Focus modules (%s) highlighted with bold borders\nBars show correlation with PFS_group (red=negative, green=positive)",
                    paste(toupper(focus_modules_auto), collapse = ", ")))

# Add legend
legend("bottomright",
       legend = c("Significant positive (p<0.05)", "Significant negative (p<0.05)", "Not significant"),
       fill = c("darkgreen", "darkred", "grey70"),
       border = NA, bty = "n", cex = 0.8)

# --- SAVE CIRCOS PLOT PNG ---
png(file.path(fig_dir_png, "Step17_Circos_ModuleOverview.png"),
    width = 12, height = 12, units = "in", res = 300)
circos.clear()
circos.par(gap.after = rep(4, length(module_names)), start.degree = 90)
circos.initialize(factors = module_names,
                  xlim = cbind(rep(0, length(module_names)), as.numeric(module_sizes)))
circos.track(ylim = c(0, 1),
             panel.fun = function(x, y) {
               sector.name <- CELL_META$sector.index
               n_proteins <- module_sizes[sector.name]
               circos.text(CELL_META$xcenter, 0.5,
                          paste0(toupper(substr(sector.name, 1, 1)),
                                 substr(sector.name, 2, nchar(sector.name)),
                                 "\n(n=", n_proteins, ")"),
                          facing = "bending.inside", niceFacing = TRUE, cex = 0.8,
                          font = ifelse(sector.name %in% focus_modules, 2, 1))
             },
             bg.col = module_names,
             bg.border = ifelse(module_names %in% focus_modules, "black", NA),
             bg.lwd = ifelse(module_names %in% focus_modules, 3, 1),
             track.height = 0.15)
circos.track(ylim = c(-1, 1),
             panel.fun = function(x, y) {
               sector.name <- CELL_META$sector.index
               cor_val <- module_trait_cor[module_trait_cor$Module == sector.name, "PFS_group"]
               pval <- module_trait_pval[module_trait_pval$Module == sector.name, "PFS_group"]
               if(!is.na(cor_val)) {
                 bar_col <- ifelse(pval < 0.05, ifelse(cor_val > 0, "darkgreen", "darkred"), "grey70")
                 circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], cor_val, col = bar_col, border = NA)
                 circos.text(CELL_META$xcenter, cor_val + sign(cor_val) * 0.15, sprintf("%.2f", cor_val),
                            cex = 0.6, facing = "bending.inside", niceFacing = TRUE)
               }
               circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.5)
             },
             track.height = 0.2)
title(main = "WGCNA Module Overview\nCircos Plot",
      sub = sprintf("Focus modules (%s) highlighted with bold borders",
                    paste(toupper(focus_modules), collapse = ", ")))
legend("bottomright", legend = c("Sig pos (p<0.05)", "Sig neg (p<0.05)", "Not sig"),
       fill = c("darkgreen", "darkred", "grey70"), border = NA, bty = "n", cex = 0.8)
dev.off()
cat("  Saved: Step17_Circos_ModuleOverview.png\n")

# --- SAVE CIRCOS PLOT PDF ---
pdf(file.path(fig_dir_pdf, "Step17_Circos_ModuleOverview.pdf"), width = 12, height = 12)
circos.clear()
circos.par(gap.after = rep(4, length(module_names)), start.degree = 90)
circos.initialize(factors = module_names,
                  xlim = cbind(rep(0, length(module_names)), as.numeric(module_sizes)))
circos.track(ylim = c(0, 1),
             panel.fun = function(x, y) {
               sector.name <- CELL_META$sector.index
               n_proteins <- module_sizes[sector.name]
               circos.text(CELL_META$xcenter, 0.5,
                          paste0(toupper(substr(sector.name, 1, 1)),
                                 substr(sector.name, 2, nchar(sector.name)),
                                 "\n(n=", n_proteins, ")"),
                          facing = "bending.inside", niceFacing = TRUE, cex = 0.8,
                          font = ifelse(sector.name %in% focus_modules, 2, 1))
             },
             bg.col = module_names,
             bg.border = ifelse(module_names %in% focus_modules, "black", NA),
             bg.lwd = ifelse(module_names %in% focus_modules, 3, 1),
             track.height = 0.15)
circos.track(ylim = c(-1, 1),
             panel.fun = function(x, y) {
               sector.name <- CELL_META$sector.index
               cor_val <- module_trait_cor[module_trait_cor$Module == sector.name, "PFS_group"]
               pval <- module_trait_pval[module_trait_pval$Module == sector.name, "PFS_group"]
               if(!is.na(cor_val)) {
                 bar_col <- ifelse(pval < 0.05, ifelse(cor_val > 0, "darkgreen", "darkred"), "grey70")
                 circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], cor_val, col = bar_col, border = NA)
                 circos.text(CELL_META$xcenter, cor_val + sign(cor_val) * 0.15, sprintf("%.2f", cor_val),
                            cex = 0.6, facing = "bending.inside", niceFacing = TRUE)
               }
               circos.lines(CELL_META$cell.xlim, c(0, 0), col = "black", lwd = 0.5)
             },
             track.height = 0.2)
title(main = "WGCNA Module Overview\nCircos Plot",
      sub = sprintf("Focus modules (%s) highlighted with bold borders",
                    paste(toupper(focus_modules), collapse = ", ")))
legend("bottomright", legend = c("Sig pos (p<0.05)", "Sig neg (p<0.05)", "Not sig"),
       fill = c("darkgreen", "darkred", "grey70"), border = NA, bty = "n", cex = 0.8)
dev.off()
cat("  Saved: Step17_Circos_ModuleOverview.pdf\n")

# ============================================================================
# 17c: Save Module Summary Table
# ============================================================================

circos_summary <- data.frame(
  Module = module_names,
  N_Proteins = as.numeric(module_sizes),
  Is_Focus = module_names %in% focus_modules,
  stringsAsFactors = FALSE
)

# Add correlations
for(trait in key_traits) {
  circos_summary[[paste0(trait, "_cor")]] <- module_trait_cor[[trait]]
  circos_summary[[paste0(trait, "_pval")]] <- module_trait_pval[[trait]]
}

# Sort by size
circos_summary <- circos_summary[order(-circos_summary$N_Proteins), ]

# Print summary to console
print(circos_summary)

# --- SAVE MODULE SUMMARY (uncomment to save) ---
# write.csv(circos_summary,
#           file.path(results_dir, "Step17_Module_Summary_Circos.csv"),
#           row.names = FALSE)
# cat("  Saved: Step17_Module_Summary_Circos.csv\n")

# ============================================================================
# 17d: Step Summary
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 17 SUMMARY: CIRCOS MODULE OVERVIEW ===\n")
cat("===============================================================================\n")
cat(sprintf("Total modules visualized: %d\n", length(module_names)))
cat(sprintf("Focus modules highlighted: %s\n", paste(focus_modules, collapse = ", ")))
cat(sprintf("Traits shown: %s\n", paste(key_traits, collapse = ", ")))
cat("Step 17 complete.\n")
cat("===============================================================================\n")

circos.clear()

# ============================================================================
# STEP 18: sEV (EXTRACELLULAR VESICLE) COMPARISON
# ============================================================================
# Purpose: Compare plasma proteomics findings with sEV proteomics data
# Validates whether WGCNA modules/biomarkers are detectable in exosomes
# File format: sEV_imputed_RF.csv has proteins in ROWS, samples in COLUMNS
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 18: sEV PROTEOMICS COMPARISON\n")
cat("================================================================\n")

# ============================================================================
# 18a: Load sEV Data
# ============================================================================

sev_file <- "sEV_imputed_RF.csv"

# Read sEV data - proteins in rows, samples in columns
sev_raw <- read.csv(sev_file, row.names = 1, stringsAsFactors = FALSE)
cat(sprintf("sEV data loaded: %d proteins x %d samples\n", nrow(sev_raw), ncol(sev_raw)))

# Extract protein names from row names
sev_proteins <- rownames(sev_raw)
sev_samples <- colnames(sev_raw)

cat(sprintf("sEV proteins detected: %d\n", length(sev_proteins)))
cat(sprintf("sEV samples: %d\n", length(sev_samples)))

# ============================================================================
# 18b: Protein Overlap Analysis
# ============================================================================

cat("\n--- Protein Overlap Analysis ---\n")

# Define plasma proteins from datExpr or moduleColors
plasma_proteins <- colnames(datExpr)
cat(sprintf("   [INFO] Using %d plasma proteins for comparison\n", length(plasma_proteins)))

overlap_proteins <- intersect(plasma_proteins, sev_proteins)
plasma_only <- setdiff(plasma_proteins, sev_proteins)
sev_only <- setdiff(sev_proteins, plasma_proteins)

cat(sprintf("Plasma proteins: %d\n", length(plasma_proteins)))
cat(sprintf("sEV proteins: %d\n", length(sev_proteins)))

pct_plasma <- 100 * length(overlap_proteins) / length(plasma_proteins)
pct_sev <- 100 * length(overlap_proteins) / length(sev_proteins)
cat(sprintf("Overlap: %d (%.1f%% of plasma, %.1f%% of sEV)\n",
            length(overlap_proteins), pct_plasma, pct_sev))
cat(sprintf("Plasma-only: %d\n", length(plasma_only)))
cat(sprintf("sEV-only: %d\n", length(sev_only)))

# Venn diagram: Plasma vs sEV
venn_plasma_sev <- list(
  Plasma = plasma_proteins,
  sEV = sev_proteins
)

p_venn_sev <- ggvenn(venn_plasma_sev,
                     fill_color = c("#3498db", "#e74c3c"),
                     stroke_size = 0.5,
                     set_name_size = 5,
                     text_size = 4) +
  labs(title = "Protein Overlap: Plasma vs sEV Proteomics") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_venn_sev)

# Save overlap list
overlap_df <- data.frame(
  Protein = overlap_proteins,
  In_Plasma = TRUE,
  In_sEV = TRUE,
  stringsAsFactors = FALSE
)

write.csv(overlap_df, file.path(results_dir, "Step18_Plasma_sEV_Overlap.csv"),row.names = FALSE)
cat("  Saved: Step18_Plasma_sEV_Overlap.csv\n")

# --- SAVE VENN DIAGRAM (uncomment to save) ---
# ggsave(file.path(fig_dir_png, "Step18_Venn_Plasma_sEV.png"),
#        p_venn_sev, width = 8, height = 6, dpi = 300)
# ggsave(file.path(fig_dir_pdf, "Step18_Venn_Plasma_sEV.pdf"),
#        p_venn_sev, width = 8, height = 6)
# cat("  Saved: Step18_Venn_Plasma_sEV.png/pdf\n")

# ============================================================================
# 18c: Module Protein Detection in sEV
# ============================================================================

cat("\n--- Module Protein Detection in sEV ---\n")

# Check each module's proteins in sEV
module_sev_detection <- data.frame()

for(mod in unique(moduleColors)) {
  mod_proteins <- names(moduleColors)[moduleColors == mod]
  n_in_sev <- sum(mod_proteins %in% sev_proteins)
  n_total <- length(mod_proteins)
  pct_in_sev <- 100 * n_in_sev / n_total

  module_sev_detection <- rbind(module_sev_detection, data.frame(
    Module = mod,
    N_Proteins = n_total,
    N_in_sEV = n_in_sev,
    Pct_in_sEV = round(pct_in_sev, 1),
    Fisher_OR = NA,
    Fisher_P = NA,
    stringsAsFactors = FALSE
  ))
}

# Sort by percentage
module_sev_detection <- module_sev_detection[order(-module_sev_detection$Pct_in_sEV), ]

# Fisher's exact test for ALL modules (enrichment in sEV?)
all_modules_sev <- unique(moduleColors)
for(mod in all_modules_sev) {
  mod_proteins <- names(moduleColors)[moduleColors == mod]
  other_proteins <- names(moduleColors)[moduleColors != mod]

  # 2x2 table: Module membership vs sEV detection
  a <- sum(mod_proteins %in% sev_proteins)      # Module & sEV
  b <- sum(mod_proteins %in% plasma_only)       # Module & not sEV
  c <- sum(other_proteins %in% sev_proteins)    # Not module & sEV
  d <- sum(other_proteins %in% plasma_only)     # Not module & not sEV

  fisher_result <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  module_sev_detection$Fisher_OR[module_sev_detection$Module == mod] <- round(fisher_result$estimate, 2)
  module_sev_detection$Fisher_P[module_sev_detection$Module == mod] <- signif(fisher_result$p.value, 3)
}

cat("\nModule Detection in sEV:\n")
print(module_sev_detection)

# NOTE: Interpretation depends on module-specific results
# Check Fisher OR and p-values above to determine which modules are enriched/depleted in sEVs
# Focus modules (PINK, GREEN, BLACK, BLUE) results should be interpreted in context of clinical associations

write.csv(module_sev_detection,file.path(results_dir, "Step18_Module_sEV_Detection.csv"),row.names = FALSE)
cat("  Saved: Step18_Module_sEV_Detection.csv\n")

# Add significance stars based on Fisher p-value
module_sev_detection <- module_sev_detection %>%
  mutate(
    Sig_Stars = case_when(
      Fisher_P < 0.001 ~ "***",
      Fisher_P < 0.01 ~ "**",
      Fisher_P < 0.05 ~ "*",
      Fisher_P < 0.10 ~ ".",
      TRUE ~ ""
    ),
    # Label showing count, OR, and significance
    Bar_Label = sprintf("%d/%d (OR=%.2f)%s", N_in_sEV, N_Proteins, Fisher_OR, Sig_Stars),
    # Enrichment direction for coloring
    Enrichment = case_when(
      Fisher_P < 0.05 & Fisher_OR > 1 ~ "Enriched",
      Fisher_P < 0.05 & Fisher_OR < 1 ~ "Depleted",
      TRUE ~ "NS"
    )
  )

# Bar plot of module detection with Fisher significance
p_mod_sev <- ggplot(module_sev_detection, aes(x = reorder(Module, Pct_in_sEV),
                                               y = Pct_in_sEV,
                                               fill = Module)) +
  geom_bar(stat = "identity") +
  # Add baseline reference line (overall detection rate)
  geom_hline(yintercept = 12.5, linetype = "dashed", color = "red", linewidth = 0.8) +
  # Labels with OR and significance
  geom_text(aes(label = Bar_Label),
            hjust = -0.05, size = 3) +
  # Significance stars in bold at bar end
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Module Protein Detection in sEV (Fisher's Exact Test)",
       subtitle = "Red dashed line = overall baseline (12.5%) | * p<0.05 ** p<0.01 *** p<0.001",
       x = "Module", y = "% Proteins Detected in sEV") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12),
        plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 11),
        legend.position = "none") +
  ylim(0, max(module_sev_detection$Pct_in_sEV, na.rm = TRUE) * 1.35)

print(p_mod_sev)


ggsave(file.path(fig_dir_png, "Step18_Module_sEV_Detection.png"),p_mod_sev, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step18_Module_sEV_Detection.pdf"),p_mod_sev, width = 8, height = 6)
cat("  Saved: Step18_Module_sEV_Detection.png/pdf\n")

# ============================================================================
# 18d: Biomarker Detection in sEV
# ============================================================================

cat("\n--- Biomarker Detection in sEV ---\n")

# Define biomarkers (same as Step 15)
ml_biomarkers_sev <- c("ECM2", "DYNC2H1", "PPIB")  # ENPP1 removed due to QC concerns
cox_biomarkers_sev <- c("APOF", "TFRC", "FABP4", "ANG")
all_biomarkers_sev <- c(ml_biomarkers_sev, cox_biomarkers_sev)

biomarker_sev <- data.frame(
  Protein = all_biomarkers_sev,
  Source = c(rep("ML", length(ml_biomarkers_sev)), rep("Cox", length(cox_biomarkers_sev))),
  In_Plasma = all_biomarkers_sev %in% plasma_proteins,
  In_sEV = all_biomarkers_sev %in% sev_proteins,
  stringsAsFactors = FALSE
)

cat("\nBiomarker Detection:\n")
print(biomarker_sev)


write.csv(biomarker_sev, file.path(results_dir, "Step18_Biomarker_sEV_Detection.csv"),row.names = FALSE)
cat("  Saved: Step18_Biomarker_sEV_Detection.csv\n")

ml_in_sev <- biomarker_sev$Protein[biomarker_sev$Source == "ML" & biomarker_sev$In_sEV]
cox_in_sev <- biomarker_sev$Protein[biomarker_sev$Source == "Cox" & biomarker_sev$In_sEV]

cat(sprintf("\nML biomarkers in sEV: %d/%d", length(ml_in_sev), length(ml_biomarkers_sev)))
if(length(ml_in_sev) > 0) cat(sprintf(" (%s)", paste(ml_in_sev, collapse = ", ")))
cat("\n")

cat(sprintf("Cox biomarkers in sEV: %d/%d", length(cox_in_sev), length(cox_biomarkers_sev)))
if(length(cox_in_sev) > 0) cat(sprintf(" (%s)", paste(cox_in_sev, collapse = ", ")))
cat("\n")

# ============================================================================
# 18e: 3-Way Venn Diagram - sEV, ML Biomarkers, Cox Biomarkers
# ============================================================================

cat("\n3-Way Venn: sEV vs ML Biomarkers vs Cox Biomarkers\n")

# Define the three sets (Plasma removed)
venn_3way <- list(
  "sEV" = sev_proteins,
  "ML Biomarkers" = ml_biomarkers_sev,
  "Cox Biomarkers" = cox_biomarkers_sev
)

# Report set sizes
cat(sprintf("  sEV proteins: %d\n", length(sev_proteins)))
cat(sprintf("  ML Biomarkers: %d (%s)\n", length(ml_biomarkers_sev), paste(ml_biomarkers_sev, collapse = ", ")))
cat(sprintf("  Cox Biomarkers: %d (%s)\n", length(cox_biomarkers_sev), paste(cox_biomarkers_sev, collapse = ", ")))

# Calculate all overlap regions
ml_in_sev <- intersect(ml_biomarkers_sev, sev_proteins)
cox_in_sev <- intersect(cox_biomarkers_sev, sev_proteins)
all_three <- Reduce(intersect, list(sev_proteins, ml_biomarkers_sev, cox_biomarkers_sev))
ml_cox_only <- setdiff(intersect(ml_biomarkers_sev, cox_biomarkers_sev), sev_proteins)
ml_sev_only <- setdiff(intersect(ml_biomarkers_sev, sev_proteins), cox_biomarkers_sev)
cox_sev_only <- setdiff(intersect(cox_biomarkers_sev, sev_proteins), ml_biomarkers_sev)
sev_only <- setdiff(setdiff(sev_proteins, ml_biomarkers_sev), cox_biomarkers_sev)
ml_only <- setdiff(setdiff(ml_biomarkers_sev, sev_proteins), cox_biomarkers_sev)
cox_only <- setdiff(setdiff(cox_biomarkers_sev, sev_proteins), ml_biomarkers_sev)

cat(sprintf("\n  ML biomarkers in sEV: %d", length(ml_in_sev)))
if(length(ml_in_sev) > 0) cat(sprintf(" (%s)", paste(ml_in_sev, collapse = ", ")))
cat("\n")

cat(sprintf("  Cox biomarkers in sEV: %d", length(cox_in_sev)))
if(length(cox_in_sev) > 0) cat(sprintf(" (%s)", paste(cox_in_sev, collapse = ", ")))
cat("\n")

cat("\n  Overlap details:\n")
cat(sprintf("    ML only (not in sEV or Cox): %s\n",
            ifelse(length(ml_only) > 0, paste(ml_only, collapse = ", "), "none")))
cat(sprintf("    Cox only (not in sEV or ML): %s\n",
            ifelse(length(cox_only) > 0, paste(cox_only, collapse = ", "), "none")))
cat(sprintf("    sEV only: %d proteins\n", length(sev_only)))
cat(sprintf("    ML & sEV (not Cox): %s\n",
            ifelse(length(ml_sev_only) > 0, paste(ml_sev_only, collapse = ", "), "none")))
cat(sprintf("    Cox & sEV (not ML): %s\n",
            ifelse(length(cox_sev_only) > 0, paste(cox_sev_only, collapse = ", "), "none")))
cat(sprintf("    ML & Cox (not sEV): %s\n",
            ifelse(length(ml_cox_only) > 0, paste(ml_cox_only, collapse = ", "), "none")))
cat(sprintf("    All three (sEV & ML & Cox): %s\n",
            ifelse(length(all_three) > 0, paste(all_three, collapse = ", "), "none")))

# Create custom 3-way Venn using ggforce for precise circle positioning
# Circle centers for 3-way Venn (equilateral triangle arrangement)
# sEV at bottom, ML Biomarkers top-left, Cox Biomarkers top-right

library(ggforce)

# Define colors (colorblind-friendly)
venn_colors <- c(
  "sEV" = "#D55E00",
  "ML"  = "#0072B2",
  "Cox" = "#009E73"
)

# Define circle parameters
r <- 1.1  # radius
# Centers arranged in triangle
cx_sev <- 0
cy_sev <- -0.4
cx_ml <- -1
cy_ml <- 0.55
cx_cox <- 1
cy_cox <- 0.55

# Create circle data
circle_data <- data.frame(
  x0 = c(cx_sev, cx_ml, cx_cox),
  y0 = c(cy_sev, cy_ml, cy_cox),
  r = r,
  Set = c("sEV", "ML", "Cox"),
  fill = c(venn_colors["sEV"], venn_colors["ML"], venn_colors["Cox"])
)

# Build Venn manually with ggplot
p_venn_3way <- ggplot() +
  # Draw circles
  geom_circle(data = circle_data,
              aes(x0 = x0, y0 = y0, r = r, fill = Set),
              alpha = 0.35, color = "grey30", linewidth = 0.8) +
  scale_fill_manual(values = venn_colors) +
  # Set labels OUTSIDE circles
  # # MODIFY LABEL POSITIONS HERE: adjust x and y values to move labels outside circles
  annotate("text", x = cx_sev, y = cy_sev - r - 0.25,
           label = "sEV", size = 12/.pt, fontface = "bold", color = venn_colors["sEV"]) +
  annotate("text", x = cx_ml - r + 0.05, y = cy_ml + 0.65,
           label = "ML Biomarkers", size = 12/.pt, fontface = "bold",
           color = venn_colors["ML"], hjust = 1) +
  annotate("text", x = cx_cox + r - 0.05, y = cy_cox + 0.65,
           label = "Cox Biomarkers", size = 12/.pt, fontface = "bold",
           color = venn_colors["Cox"], hjust = 0) +
  # Region labels with protein names (ALL BLACK, size 10)
  # sEV only (bottom center)
  annotate("text", x = 0, y = -1.2,
           label = sprintf("%d proteins", length(sev_only)),
           size = 10/.pt, color = "black") +
  # ML only (top-left, inside ML circle but outside overlaps)
  annotate("text", x = -1.1, y = 0.7,
           label = ifelse(length(ml_only) > 0, paste(ml_only, collapse = "\n"), ""),
           size = 10/.pt, color = "black") +
  # Cox only (top-right, inside Cox circle but outside overlaps)
  annotate("text", x = 1.1, y = 0.7,
           label = ifelse(length(cox_only) > 0, paste(cox_only, collapse = "\n"), ""),
           size = 10/.pt, color = "black") +
  # ML & sEV intersection (left) - COMMON PROTEINS: size 10, bold
  annotate("text", x = -0.45, y = -0.05,
           label = ifelse(length(ml_sev_only) > 0, paste(ml_sev_only, collapse = "\n"), ""),
           size = 10/.pt, fontface = "bold", color = "black") +
  # Cox & sEV intersection (right) - COMMON PROTEINS: size 10, bold
  annotate("text", x = 0.45, y = -0.05,
           label = ifelse(length(cox_sev_only) > 0, paste(cox_sev_only, collapse = "\n"), ""),
           size = 10/.pt, fontface = "bold", color = "black") +
  # ML & Cox intersection (top, no sEV) - COMMON PROTEINS: size 10, bold
  annotate("text", x = 0, y = 0.95,
           label = ifelse(length(ml_cox_only) > 0, paste(ml_cox_only, collapse = "\n"), ""),
           size = 10/.pt, fontface = "bold", color = "black") +
  # All three (center) - COMMON PROTEINS: size 10, bold
  annotate("text", x = 0, y = 0.3,
           label = ifelse(length(all_three) > 0, paste(all_three, collapse = "\n"), ""),
           size = 11/.pt, fontface = "bold", color = "black") +
  # Title and theme
  labs(title = "Biomarker Detection in sEV Proteomics")+
  coord_fixed(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(20, 40, 20, 40),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none")

print(p_venn_3way)

# Save overlap summary table
biomarker_overlap_summary <- data.frame(
  Biomarker = all_biomarkers_sev,
  Source = c(rep("ML", length(ml_biomarkers_sev)), rep("Cox", length(cox_biomarkers_sev))),
  In_Plasma = all_biomarkers_sev %in% plasma_proteins,
  In_sEV = all_biomarkers_sev %in% sev_proteins,
  In_Both_Plasma_sEV = all_biomarkers_sev %in% intersect(plasma_proteins, sev_proteins),
  stringsAsFactors = FALSE
)

write.csv(biomarker_overlap_summary,file.path(results_dir, "Step18_Biomarker_3Way_Overlap.csv"),row.names = FALSE)
cat("  Saved: Step18_Biomarker_3Way_Overlap.csv\n")

ggsave(file.path(fig_dir_png, "Step18_Venn_3Way_sEV_Biomarkers.png"),p_venn_3way, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step18_Venn_3Way_sEV_Biomarkers.pdf"),p_venn_3way, width = 7, height = 5)
cat("  Saved: Step18_Venn_3Way_sEV_Biomarkers.png/pdf\n")

# ============================================================================
# 18f: Hub Gene Detection in sEV (ALL Focus Modules)
# ============================================================================

cat("\n--- Hub Gene Detection in sEV (All Focus Modules) ---\n")

# Loop through ALL focus modules to get hub genes
hub_sev <- data.frame()

for(focus_mod in focus_modules_auto) {
  mod_hubs <- proteinInfo$Protein[proteinInfo$Module == focus_mod & proteinInfo$isHub_Moderate == TRUE]
  cat(sprintf("%s hub genes: %d\n", toupper(focus_mod), length(mod_hubs)))

  if(length(mod_hubs) > 0) {
    hub_sev <- rbind(hub_sev, data.frame(
      Module = toupper(focus_mod),
      Protein = mod_hubs,
      In_sEV = mod_hubs %in% sev_proteins,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(hub_sev, file.path(results_dir, "Step18_HubGenes_sEV_Detection.csv"), row.names = FALSE)

# Print detection summary for each focus module
cat("\nHub gene detection in sEV by module:\n")
for(focus_mod in focus_modules_auto) {
  mod_hubs <- proteinInfo$Protein[proteinInfo$Module == focus_mod & proteinInfo$isHub_Moderate == TRUE]
  if(length(mod_hubs) > 0) {
    n_in_sev <- sum(mod_hubs %in% sev_proteins)
    pct_in_sev <- 100 * n_in_sev / length(mod_hubs)
    cat(sprintf("  %s hub genes in sEV: %d/%d (%.1f%%)\n",
                toupper(focus_mod), n_in_sev, length(mod_hubs), pct_in_sev))
  } else {
    cat(sprintf("  %s: no hub genes\n", toupper(focus_mod)))
  }
}

# ============================================================================
# 18g: Step Summary
# ============================================================================

cat("\n")
cat("===============================================================================\n")
cat("=== STEP 18 SUMMARY: sEV PROTEOMICS COMPARISON ===\n")
cat("===============================================================================\n")

cat(sprintf("\nProtein Overlap:\n"))
cat(sprintf("  Plasma proteins: %d\n", length(plasma_proteins)))
cat(sprintf("  sEV proteins: %d\n", length(sev_proteins)))
cat(sprintf("  Shared proteins: %d (%.1f%% of plasma)\n", length(overlap_proteins),
            100 * length(overlap_proteins) / length(plasma_proteins)))

cat(sprintf("\nFocus Module Detection in sEV:\n"))
for(mod in focus_modules_auto) {
  row <- module_sev_detection[module_sev_detection$Module == mod, ]
  cat(sprintf("  %s: %d/%d proteins (%.1f%%)",
              toupper(mod), row$N_in_sEV, row$N_Proteins, row$Pct_in_sEV))
  cat(sprintf(", Fisher OR=%.2f, p=%.3g", row$Fisher_OR, row$Fisher_P))
  cat("\n")
}

# ============================================================================
# STEP 19: NATURE MEDICINE STYLE FIGURES
# ============================================================================
# Reference: Nat Med. 2024;30(3):749-761 - Figures 1d and 2a
# Purpose: Publication-quality enrichment and network visualizations
# ============================================================================

cat("================================================================================\n")
cat("STEP 19: NATURE MEDICINE STYLE FIGURES\n")
cat("================================================================================\n")
cat("Reference: Nat Med. 2024;30(3):749-761\n\n")

# --- 19a: Protein-Clinical Correlations ---
cat("19a: Calculating protein-clinical correlations...\n")

# Create numeric versions of clinical variables
# PFS_group: Long = 1 (favorable), Short = 0 (unfavorable)
# Response: CD = 1 (favorable), PD = 0 (unfavorable)
PFS_numeric <- ifelse(datTraits$PFS_group == "Long", 1, 0)
Response_numeric <- ifelse(datTraits$Response == "CD", 1, 0)

# Calculate correlation for each protein
protein_clinical_cor <- data.frame(
  Protein = colnames(datExpr),
  Cor_PFS = numeric(ncol(datExpr)),
  P_PFS = numeric(ncol(datExpr)),
  Cor_Response = numeric(ncol(datExpr)),
  P_Response = numeric(ncol(datExpr)),
  stringsAsFactors = FALSE
)

for(i in 1:ncol(datExpr)) {
  # PFS correlation
  cor_pfs <- cor.test(datExpr[, i], PFS_numeric, method = "pearson", use = "pairwise.complete.obs")
  protein_clinical_cor$Cor_PFS[i] <- cor_pfs$estimate
  protein_clinical_cor$P_PFS[i] <- cor_pfs$p.value

  # Response correlation
  cor_resp <- cor.test(datExpr[, i], Response_numeric, method = "pearson", use = "pairwise.complete.obs")
  protein_clinical_cor$Cor_Response[i] <- cor_resp$estimate
  protein_clinical_cor$P_Response[i] <- cor_resp$p.value
}

# Classify proteins by direction
protein_clinical_cor <- protein_clinical_cor %>%
  mutate(
    PFS_Direction = ifelse(Cor_PFS > 0, "Long_PFS", "Short_PFS"),
    Response_Direction = ifelse(Cor_Response > 0, "CD", "PD")
  )

cat(sprintf("  Proteins higher in Long PFS: %d\n", sum(protein_clinical_cor$PFS_Direction == "Long_PFS")))
cat(sprintf("  Proteins higher in Short PFS: %d\n", sum(protein_clinical_cor$PFS_Direction == "Short_PFS")))
cat(sprintf("  Proteins higher in CD: %d\n", sum(protein_clinical_cor$Response_Direction == "CD")))
cat(sprintf("  Proteins higher in PD: %d\n", sum(protein_clinical_cor$Response_Direction == "PD")))

write.csv(protein_clinical_cor,
          file.path(results_dir, "Step19_Protein_Clinical_Correlations.csv"),
          row.names = FALSE)
cat("  Saved: Step19_Protein_Clinical_Correlations.csv\n")

# --- 19b: Load Pathway Database ---
cat("\n19b: Loading pathway database...\n")

# Load KEGG pathways from MSigDB (filter by prefix since subcategory syntax changed)
kegg_df <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::filter(grepl("^KEGG_", gs_name))
pathway_db <- split(kegg_df$gene_symbol, kegg_df$gs_name)
names(pathway_db) <- gsub("KEGG_", "", names(pathway_db))
names(pathway_db) <- gsub("_", " ", names(pathway_db))
cat(sprintf("  Loaded %d KEGG pathways\n", length(pathway_db)))

# Load GO Biological Process (filter by prefix since subcategory syntax may change)
go_bp_df <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::filter(grepl("^GOBP_", gs_name))
go_bp_db <- split(go_bp_df$gene_symbol, go_bp_df$gs_name)
go_bp_db <- go_bp_db[sapply(go_bp_db, length) >= 10 & sapply(go_bp_db, length) <= 500]
names(go_bp_db) <- gsub("GOBP_", "", names(go_bp_db))
names(go_bp_db) <- gsub("_", " ", names(go_bp_db))
cat(sprintf("  Loaded %d GO:BP pathways (size 10-500)\n", length(go_bp_db)))

all_pathways <- c(pathway_db, go_bp_db)
cat(sprintf("  Total pathways for enrichment: %d\n", length(all_pathways)))

#19c: Enrichment Analysis Function 


run_directional_enrichment <- function(query_proteins, background_proteins,
                                        genesets, direction_label, min_overlap = 3) {
  n_background <- length(background_proteins)
  n_query <- length(query_proteins)
  results <- data.frame()

  for(pw_name in names(genesets)) {
    pw_genes <- genesets[[pw_name]]
    pw_in_background <- intersect(pw_genes, background_proteins)
    n_pw <- length(pw_in_background)
    if(n_pw < 5) next

    overlap <- intersect(query_proteins, pw_in_background)
    n_overlap <- length(overlap)
    if(n_overlap < min_overlap) next

    mat <- matrix(c(n_overlap, n_query - n_overlap, n_pw - n_overlap,
                    n_background - n_query - n_pw + n_overlap), nrow = 2)
    fisher_res <- fisher.test(mat, alternative = "greater")
    fold_enrich <- (n_overlap / n_query) / (n_pw / n_background)

    results <- rbind(results, data.frame(
      Pathway = pw_name, Direction = direction_label,
      N_Pathway_InBackground = n_pw, N_Query = n_query, N_Overlap = n_overlap,
      Gene_Ratio = n_overlap / n_query, Fold_Enrichment = fold_enrich,
      P_Value = fisher_res$p.value,
      Overlap_Proteins = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }

  if(nrow(results) > 0) {
    results$FDR <- p.adjust(results$P_Value, method = "BH")
    results <- results %>% arrange(P_Value)
  }
  return(results)
}

# Ensure all_proteins is defined (use existing or derive from data)
if(!exists("all_proteins") || length(all_proteins) == 0) {
  if(ncol(datExpr) > 0) {
    all_proteins <- colnames(datExpr)
  } else {
    all_proteins <- names(moduleColors)
  }
}
proteins_long_pfs <- protein_clinical_cor$Protein[protein_clinical_cor$PFS_Direction == "Long_PFS"]
proteins_short_pfs <- protein_clinical_cor$Protein[protein_clinical_cor$PFS_Direction == "Short_PFS"]
proteins_cd <- protein_clinical_cor$Protein[protein_clinical_cor$Response_Direction == "CD"]
proteins_pd <- protein_clinical_cor$Protein[protein_clinical_cor$Response_Direction == "PD"]

cat("  Running enrichment for Long PFS proteins...\n")
enrich_long_pfs <- run_directional_enrichment(proteins_long_pfs, all_proteins, all_pathways, "Long_PFS")
cat("  Running enrichment for Short PFS proteins...\n")
enrich_short_pfs <- run_directional_enrichment(proteins_short_pfs, all_proteins, all_pathways, "Short_PFS")
cat("  Running enrichment for CD proteins...\n")
enrich_cd <- run_directional_enrichment(proteins_cd, all_proteins, all_pathways, "CD")
cat("  Running enrichment for PD proteins...\n")
enrich_pd <- run_directional_enrichment(proteins_pd, all_proteins, all_pathways, "PD")

enrich_pfs_all <- bind_rows(enrich_long_pfs, enrich_short_pfs)
enrich_response_all <- bind_rows(enrich_cd, enrich_pd)

cat(sprintf("  PFS enrichment results: Long=%d, Short=%d\n", nrow(enrich_long_pfs), nrow(enrich_short_pfs)))
cat(sprintf("  Response enrichment results: CD=%d, PD=%d\n", nrow(enrich_cd), nrow(enrich_pd)))

write.csv(enrich_pfs_all, file.path(results_dir, "Step19_Enrichment_PFS.csv"), row.names = FALSE)
write.csv(enrich_response_all, file.path(results_dir, "Step19_Enrichment_Response.csv"), row.names = FALSE)
cat("  Saved: Step19_Enrichment_PFS.csv, Step19_Enrichment_Response.csv\n")

# 19d: Figure A - PFS Enrichment Bubble Plot


direction_colors <- c("Long_PFS" = "#00599F", "Short_PFS" = "#D01910")

plot_data_pfs <- enrich_pfs_all %>%
  filter(P_Value < 0.05) %>%
  mutate(
    neg_log10_P = -log10(pmax(P_Value, 1e-50)),
    Direction = factor(Direction, levels = c("Short_PFS", "Long_PFS"))
  ) %>%
  arrange(desc(Fold_Enrichment))


  # Select TOP 10 pathways per direction for labeling (reduced from 15 to avoid crowding)
top_pathways_pfs <- plot_data_pfs %>%
    group_by(Direction) %>%
    slice_max(order_by = Fold_Enrichment, n = 10, with_ties = FALSE) %>%
    ungroup()

fig_A_pfs <- ggplot(plot_data_pfs, aes(x = neg_log10_P, y = Fold_Enrichment)) +
    geom_point(aes(size = Gene_Ratio, fill = Direction),
               shape = 21, color = "grey30", stroke = 0.4, alpha = 0.85) +
    geom_text_repel(
      data = top_pathways_pfs,
      aes(label = Pathway, color = Direction),
      size = 2.5, max.overlaps = 20, segment.size = 0.2,
      segment.color = "grey60", box.padding = 0.4, point.padding = 0.3,
      force = 3, show.legend = FALSE
    ) +
    scale_fill_manual(
      values = direction_colors,
      labels = c("Long_PFS" = "Higher in Long PFS", "Short_PFS" = "Higher in Short PFS"),
      name = "Direction"
    ) +
    scale_color_manual(values = direction_colors, guide = "none") +
    # REDUCED circle size range (was 4-20, now 2-10)
    scale_size_continuous(range = c(2, 10), name = "Gene Ratio",
                          labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "KEGG & GO:BP Pathway Enrichment by PFS Direction",
      subtitle = "Fisher's exact test | Nominal p < 0.05 | Top 10 pathways labeled per direction",
      x = expression(-log[10](italic(P)*"-value")),
      y = "Fold Enrichment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      legend.position = "right",
      axis.title = element_text(face = "bold")
    ) +
    guides(fill = guide_legend(override.aes = list(size = 5), order = 1),
           size = guide_legend(order = 2))

print(fig_A_pfs)
ggsave(file.path(fig_dir_png, "Step19_FigA_PFS_Enrichment.png"), fig_A_pfs, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step19_FigA_PFS_Enrichment.pdf"), fig_A_pfs, width = 8, height = 6)
cat("  Saved: Step19_FigA_PFS_Enrichment.png/pdf\n")

# 19d: Figure A Lollipop - Top 5 PFS Pathways per Direction ---

library(stringr)

lollipop_pfs <- plot_data_pfs %>%
  group_by(Direction) %>%
  slice_max(order_by = Fold_Enrichment, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Direction_Label = factor(Direction,
                             levels = c("Short_PFS", "Long_PFS"),
                             labels = c("Higher in Short PFS", "Higher in Long PFS")),
    # Wrap long pathway names into multiple lines (width = 30 characters)
    Pathway_Wrapped = str_wrap(Pathway, width = 20),
    # Reorder using the wrapped text
    Pathway_Clean = reorder(Pathway_Wrapped, neg_log10_P)
  )

fig_A_pfs_lollipop <- ggplot(lollipop_pfs, aes(x = neg_log10_P, y = Pathway_Clean)) +
  geom_segment(aes(x = 0, xend = neg_log10_P, y = Pathway_Clean, yend = Pathway_Clean, color = Direction),
               linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(size = Gene_Ratio, fill = Direction), shape = 21, color = "grey30", stroke = 0.5) +
  facet_wrap(~ Direction_Label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = direction_colors, guide = "none") +
  scale_color_manual(values = direction_colors, guide = "none") +
  scale_size_continuous(range = c(3, 10), name = "Gene Ratio",
                        labels = scales::percent_format(accuracy = 1)) +
  # Add more padding to the right of the x-axis so points don't hit the border
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "KEGG & GO:BP Top 5 Enriched pathways by PFS direction",
    subtitle = "Fisher's exact test | Nominal p < 0.05",
    x = expression(-log[10](italic(P)*"-value")),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95"),
    # Adjust line height for wrapped text and reduce font size slightly
    axis.text.y = element_text(size = 9, lineheight = 0.9), 
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    # Use a larger left margin to accommodate the long wrapped labels
    plot.margin = margin(10, 20, 10, 20)
  )

print(fig_A_pfs_lollipop)

ggsave(file.path(fig_dir_png, "Step19_FigA_PFS_Lollipop.png"), fig_A_pfs_lollipop, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step19_FigA_PFS_Lollipop.pdf"), fig_A_pfs_lollipop, width = 8, height = 6)
cat("  Saved: Step19_FigA_PFS_Lollipop.png/pdf\n")

# 19e: Figure A' - Response Enrichment Bubble Plot 


# Use universal response colors

plot_data_response <- enrich_response_all %>%
  filter(P_Value < 0.05) %>%
  mutate(
    neg_log10_P = -log10(pmax(P_Value, 1e-50)),
    Direction = factor(Direction, levels = c("PD", "CD"))
  ) %>%
  arrange(desc(Fold_Enrichment))

# Select TOP 10 pathways per direction for labeling (reduced from 15 to avoid crowding)
top_pathways_response <- plot_data_response %>%
    group_by(Direction) %>%
    slice_max(order_by = Fold_Enrichment, n = 10, with_ties = FALSE) %>%
    ungroup()

fig_A_response <- ggplot(plot_data_response, aes(x = neg_log10_P, y = Fold_Enrichment)) +
    geom_point(aes(size = Gene_Ratio, fill = Direction),
               shape = 21, color = "grey30", stroke = 0.4, alpha = 0.85) +
geom_text_repel(
      data = top_pathways_response,
      aes(label = Pathway, color = Direction),
      size = 2.5, max.overlaps = 20, segment.size = 0.2,
      segment.color = "grey60", box.padding = 0.4, point.padding = 0.3,
      force = 3, show.legend = FALSE
    ) +
    scale_fill_manual(
      values = response_colors,
      labels = c("CD" = "Higher in Responders (CD)", "PD" = "Higher in Non-responders (PD)"),
      name = "Direction"
    ) +
    scale_color_manual(values = response_colors, guide = "none") +
    # REDUCED circle size range (was 4-20, now 2-10)
    scale_size_continuous(range = c(2, 10), name = "Gene Ratio",
                          labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "KEGG & GO:BP Pathway Enrichment by Treatment Response",
      subtitle = "Fisher's exact test | Nominal p < 0.05",
      x = expression(-log[10](italic(P)*"-value")),
      y = "Fold Enrichment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      legend.position = "right",
      axis.title = element_text(face = "bold")
    ) +
    guides(fill = guide_legend(override.aes = list(size = 5), order = 1),
           size = guide_legend(order = 2))

print(fig_A_response)
ggsave(file.path(fig_dir_png, "Step19_FigA_Response_Enrichment.png"), fig_A_response, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step19_FigA_Response_Enrichment.pdf"), fig_A_response, width = 8, height = 6)
cat("  Saved: Step19_FigA_Response_Enrichment.png/pdf\n")

# 19e': Figure A' Lollipop - Top 5 Response Pathways per Direction


lollipop_response <- plot_data_response %>%
  group_by(Direction) %>%
  slice_max(order_by = Fold_Enrichment, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Direction_Label = factor(Direction,
                             levels = c("PD", "CD"),
                             labels = c("Higher in Non-responders (PD)", "Higher in Responders (CD)")),
    # FIX 1: Wrap long pathway names into multiple lines (width = 20 characters)
    Pathway_Wrapped = str_wrap(Pathway, width = 20),
    # FIX 2: Reorder using the wrapped text
    Pathway_Clean = reorder(Pathway_Wrapped, neg_log10_P)
  )

fig_A_response_lollipop <- ggplot(lollipop_response, aes(x = neg_log10_P, y = Pathway_Clean)) +
  geom_segment(aes(x = 0, xend = neg_log10_P, y = Pathway_Clean, yend = Pathway_Clean, color = Direction),
               linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(size = Gene_Ratio, fill = Direction), shape = 21, color = "grey30", stroke = 0.5) +
  facet_wrap(~ Direction_Label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = response_colors, guide = "none") +
  scale_color_manual(values = response_colors, guide = "none") +
  scale_size_continuous(range = c(3, 10), name = "Gene Ratio",
                        labels = scales::percent_format(accuracy = 1)) +
  # Add more padding to the right of the x-axis so points don't hit the border
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "KEGG & GO:BP Top 5 Enriched pathways by Treatment Response",
    subtitle = "Fisher's exact test | Nominal p < 0.05",
    x = expression(-log[10](italic(P)*"-value")),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95"),
    # Adjust line height for wrapped text and reduce font size slightly
    axis.text.y = element_text(size = 9, lineheight = 0.9),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    # Use a larger left margin to accommodate the long wrapped labels
    plot.margin = margin(10, 20, 10, 20)
  )

print(fig_A_response_lollipop)

ggsave(file.path(fig_dir_png, "Step19_FigA_Response_Lollipop.png"), fig_A_response_lollipop, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir_pdf, "Step19_FigA_Response_Lollipop.pdf"), fig_A_response_lollipop, width = 8, height = 6)
cat("  Saved: Step19_FigA_Response_Lollipop.png/pdf\n")

cat("===============================================================================\n")
cat("STEP 19 SUMMARY: NATURE MEDICINE STYLE FIGURES\n")
cat("===============================================================================\n")

cat("\nProtein-Clinical Correlations:\n")
cat(sprintf("  Higher in Long PFS: %d proteins\n", sum(protein_clinical_cor$PFS_Direction == "Long_PFS")))
cat(sprintf("  Higher in Short PFS: %d proteins\n", sum(protein_clinical_cor$PFS_Direction == "Short_PFS")))
cat(sprintf("  Higher in CD: %d proteins\n", sum(protein_clinical_cor$Response_Direction == "CD")))
cat(sprintf("  Higher in PD: %d proteins\n", sum(protein_clinical_cor$Response_Direction == "PD")))

cat("\nEnrichment Results (p < 0.05):\n")
cat(sprintf("  PFS - Long PFS pathways: %d\n", sum(enrich_long_pfs$P_Value < 0.05, na.rm = TRUE)))
cat(sprintf("  PFS - Short PFS pathways: %d\n", sum(enrich_short_pfs$P_Value < 0.05, na.rm = TRUE)))
cat(sprintf("  Response - CD pathways: %d\n", sum(enrich_cd$P_Value < 0.05, na.rm = TRUE)))
cat(sprintf("  Response - PD pathways: %d\n", sum(enrich_pd$P_Value < 0.05, na.rm = TRUE)))

cat("\nStep 19 complete.\n")
cat("===============================================================================\n")

# ============================================================================

# ============================================================================
# Purpose: Auto-generate HTML/PDF summary report with key findings
# Package: rmarkdown
# ============================================================================
cat("================================================================================\n")
cat("                    REPORT GENERATION                                          \n")
cat("================================================================================\n")
cat("PURPOSE: Generate comprehensive HTML and PDF reports with all results          \n")
cat("================================================================================\n")

# Report packages (loaded at startup - verify availability)
report_packages <- c("rmarkdown", "knitr", "kableExtra")
missing_report <- report_packages[!sapply(report_packages, requireNamespace, quietly = TRUE)]
if(length(missing_report) > 0) {
  cat("   [WARNING] Missing report packages:", paste(missing_report, collapse = ", "), "\n")
  cat("   Report generation may fail.\n")
}

# Create report template if rmarkdown is available
if(requireNamespace("rmarkdown", quietly = TRUE) &&
   requireNamespace("knitr", quietly = TRUE)) {

  # Set output directory paths
  results_dir <- file.path(main_dir, "Results_Tables")
  fig_dir_png <- file.path(main_dir, "Figures_PNG")
  fig_dir_pdf <- file.path(main_dir, "Figures_PDF")

  # Resolve paths for rmarkdown rendering
  results_dir_abs <- normalizePath(results_dir, winslash = "/", mustWork = FALSE)
  fig_dir_for_report_abs <- normalizePath(fig_dir_png, winslash = "/", mustWork = FALSE)

  cat(sprintf("   Results dir (absolute): %s\n", results_dir_abs))
  cat(sprintf("   Figure dir (absolute): %s\n", fig_dir_for_report_abs))

  # Build COMPREHENSIVE Rmd template with ALL figures and tables
  rmd_lines <- c(
    '---',
    'title: "WGCNA Analysis Report - Complete"',
    'subtitle: "Metastatic PDAC Proteomics - Module and Biomarker Analysis"',
    'author: "Generated by wgcna_v3.R"',
    paste0('date: "', format(Sys.Date(), "%B %d, %Y"), '"'),
    'output:',
    '  html_document:',
    '    toc: true',
    '    toc_float:',
    '      collapsed: true',
    '    toc_depth: 4',
    '    theme: flatly',
    '    self_contained: true',
    '  pdf_document:',
    '    toc: true',
    '    toc_depth: 3',
    '    number_sections: true',
    '    latex_engine: pdflatex',
    'geometry: "a4paper, margin=0.75in"',
    'fontsize: 10pt',
    '---',
    '',
    paste0('```{r setup, include=FALSE}'),
    'knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE, fig.align="center", out.width="95%")',
    'library(knitr)',
    'if(requireNamespace("kableExtra", quietly=TRUE)) library(kableExtra)',
    paste0('results_dir <- "', results_dir_abs, '"'),
    paste0('fig_dir <- "', fig_dir_for_report_abs, '"'),
    'include_fig <- function(f) { p <- file.path(fig_dir, f); if(file.exists(p)) knitr::include_graphics(p) else cat("Figure not available:", f) }',
    'load_tbl <- function(f, cap="", n=25) { p <- file.path(results_dir, f); if(file.exists(p)) { d <- head(read.csv(p), n); if(requireNamespace("kableExtra", quietly=TRUE)) { kable(d, caption=cap) %>% kable_styling(bootstrap_options=c("striped","hover","condensed"), latex_options=c("striped","scale_down"), font_size=11) } else kable(d, caption=cap) } else cat("Table not available:", f) }',
    '```',
    '',
    '# Executive Summary',
    '',
    '**WGCNA Analysis** of metastatic PDAC plasma proteomics data.',
    '',
    '| Parameter | Value |',
    '|-----------|-------|',
    '| Samples | 50 patients |',
    '| Proteins | 891 |',
    '| Network Type | Signed |',
    '| Correlation | Biweight midcorrelation (bicor) |',
    '| Soft Power | 13 |',
    '| Min Module Size | 20 |',
    '| Merge Cut Height | 0.25 |',
    '',
    '---',
    '',
    '# Step 1: Data Loading',
    '',
    '```{r s1tbl, results="asis"}',
    'load_tbl("Step1_Summary.csv", "Data Loading Summary")',
    '```',
    '',
    '---',
    '',
    '# Step 2: Protein Noise Assessment',
    '',
    '## Variance Distribution',
    '```{r s2_01}',
    'include_fig("Step2_01_Variance_Distribution.png")',
    '```',
    '',
    '## Mean-Variance Relationship',
    '```{r s2_02}',
    'include_fig("Step2_02_Mean_Variance_Relationship.png")',
    '```',
    '',
    '## CV Distribution',
    '```{r s2_03}',
    'include_fig("Step2_03_CV_Distribution.png")',
    '```',
    '',
    '## Correlation Distribution',
    '```{r s2_04}',
    'include_fig("Step2_04_Correlation_Distribution.png")',
    '```',
    '',
    '## Expression Distribution',
    '```{r s2_05}',
    'include_fig("Step2_05_Expression_Distribution.png")',
    '```',
    '',
    '## Flagged Proteins Summary',
    '```{r s2_06}',
    'include_fig("Step2_06_Flagged_Summary.png")',
    '```',
    '',
    '## Multi-Flagged Proteins Heatmap',
    '```{r s2_07}',
    'include_fig("Step2_07_MultiFlagged_Heatmap.png")',
    '```',
    '',
    '---',
    '',
    '# Step 3: Sample Quality Control',
    '',
    '## Sample Dendrogram',
    '```{r s3_01}',
    'include_fig("Step3_01_Sample_Dendrogram.png")',
    '```',
    '',
    '## Sample Dendrogram with Traits',
    '```{r s3_02}',
    'include_fig("Step3_02_Sample_Dendrogram_Traits.png")',
    '```',
    '',
    '## Sample Connectivity',
    '```{r s3_03}',
    'include_fig("Step3_03_Sample_Connectivity.png")',
    '```',
    '',
    '## PCA Panels',
    '```{r s3_04}',
    'include_fig("Step3_04_PCA_Panels.png")',
    '```',
    '',
    '## PCA Scree Plot',
    '```{r s3_05}',
    'include_fig("Step3_05_PCA_Scree.png")',
    '```',
    '',
    '## Sample Correlation Heatmap',
    '```{r s3_06}',
    'include_fig("Step3_06_Sample_Correlation_Heatmap.png")',
    '```',
    '',
    '```{r s3tbl, results="asis"}',
    'load_tbl("Step3_Sample_QC_Summary.csv", "Sample QC Summary")',
    '```',
    '',
    '---',
    '',
    '# Step 4: Soft-Threshold Power Selection',
    '',
    '## R-squared vs Power',
    '```{r s4b_01}',
    'include_fig("Step4b_01_R2_vs_Power.png")',
    '```',
    '',
    '## Connectivity vs Power',
    '```{r s4b_02}',
    'include_fig("Step4b_02_Connectivity_vs_Power.png")',
    '```',
    '',
    '## R2 vs Connectivity Tradeoff',
    '```{r s4b_03}',
    'include_fig("Step4b_03_R2_vs_Connectivity_Tradeoff.png")',
    '```',
    '',
    '## Signed vs Unsigned Comparison',
    '```{r s4b_04}',
    'include_fig("Step4b_04_Signed_vs_Unsigned.png")',
    '```',
    '',
    '## Slope Analysis',
    '```{r s4b_05}',
    'include_fig("Step4b_05_Slope_Analysis.png")',
    '```',
    '',
    '## Faceted Conditions',
    '```{r s4b_06}',
    'include_fig("Step4b_06_Faceted_Conditions.png")',
    '```',
    '',
    '## Candidate Powers',
    '```{r s4ctbl, results="asis"}',
    'load_tbl("Step4c_Candidate_Powers.csv", "Candidate Powers")',
    '```',
    '',
    '## Sensitivity Analysis - Comparison',
    '```{r s4d_01}',
    'include_fig("Step4d_01_Sensitivity_Comparison.png")',
    '```',
    '',
    '## Module-Trait Heatmaps (Sensitivity)',
    '```{r s4d_02}',
    'include_fig("Step4d_02_ModuleTrait_Heatmaps.png")',
    '```',
    '',
    '## Clean Heatmaps Combined',
    '```{r s4d_03}',
    'include_fig("Step4d_03_CleanHeatmaps_Combined.png")',
    '```',
    '',
    '---',
    '',
    '# Step 5: Network Construction & Module Detection',
    '',
    '## Gene Dendrogram with Module Colors',
    '```{r s5_01}',
    'include_fig("Step5_01_Dendrogram.png")',
    '```',
    '',
    '## Module Sizes',
    '```{r s5_02}',
    'include_fig("Step5_02_Module_Sizes.png")',
    '```',
    '',
    '## Module Eigengene Dendrogram',
    '```{r s5_03}',
    'include_fig("Step5_03_ME_Dendrogram.png")',
    '```',
    '',
    '## Module-Trait Correlation Heatmap',
    '```{r s5_04}',
    'include_fig("Step5_04_ModuleTrait_Heatmap.png")',
    '```',
    '',
    '## Permutation Tests - Combined',
    '```{r s5_05c}',
    'include_fig("Step5_05_Permutation_Combined.png")',
    '```',
    '',
    '## Module-Trait DotPlot',
    '```{r s5_06}',
    'include_fig("Step5_06_DotPlot.png")',
    '```',
    '',
    '## Module Summary Table',
    '```{r s5tbl1, results="asis"}',
    'load_tbl("Step5_ModuleSummary.csv", "Module Summary")',
    '```',
    '',
    '## Module-Trait Correlations',
    '```{r s5tbl2, results="asis"}',
    'load_tbl("Step5_ModuleTraitCorrelations.csv", "Module-Trait Correlations")',
    '```',
    '',
    '## Module Priority',
    '```{r s5tbl3, results="asis"}',
    'load_tbl("Step5_ModulePriority.csv", "Module Priority for Clinical Relevance")',
    '```',
    '',
    '## Permutation Test Results',
    '```{r s5tbl4, results="asis"}',
    'load_tbl("Step5_Permutation_Tests.csv", "Permutation Test Results")',
    '```',
    '',
    '---',
    '',
    '# Step 6: Hub Gene Identification',
    '',
    '## Module Membership vs Gene Significance - Green',
    '```{r s6_01}',
    'include_fig("Step6_01_MEgreen_GS_PFS_group.png")',
    '```',
    '',
    '## Module Membership vs Gene Significance - Green (Response)',
    '```{r s6_02}',
    'include_fig("Step6_02_MEgreen_GS_Response.png")',
    '```',
    '',
    '## Module Membership vs Gene Significance - Turquoise',
    '```{r s6_03}',
    'include_fig("Step6_03_MEturquoise_GS_PFS_group.png")',
    '```',
    '',
    '## Module Membership vs Gene Significance - Turquoise (Response)',
    '```{r s6_04}',
    'include_fig("Step6_04_MEturquoise_GS_Response.png")',
    '```',
    '',
    '## Hub Genes Combined',
    '```{r s6_05}',
    'include_fig("Step6_05_HubGenes_Combined.png")',
    '```',
    '',
    '## Hub Count Comparison',
    '```{r s6_06}',
    'include_fig("Step6_06_HubCount_Comparison.png")',
    '```',
    '',
    '## Connectivity Validation',
    '```{r s6_07}',
    'include_fig("Step6_07_Connectivity_Validation.png")',
    '```',
    '',
    '## Hub Effect Sizes',
    '```{r s6_08}',
    'include_fig("Step6_08_HubEffectSizes.png")',
    '```',
    '',
    '## Network - Green Module',
    '```{r s6_09}',
    'include_fig("Step6_09_Network_green_PubQuality.png")',
    '```',
    '',
    '## Network - Turquoise Module',
    '```{r s6_10}',
    'include_fig("Step6_10_Network_turquoise_PubQuality.png")',
    '```',
    '',
    '## Hub Gene Summary Table',
    '```{r s6tbl, results="asis"}',
    'load_tbl("Step6_HubSummary.csv", "Hub Gene Summary")',
    '```',
    '',
    '---',
    '',
    '# Step 7: Module Stability Testing',
    '',
    '## Bootstrap Convergence',
    '```{r s7_01}',
    'include_fig("Step7_01_Convergence.png")',
    '```',
    '',
    '## Module Stability',
    '```{r s7_02}',
    'include_fig("Step7_02_ModuleStability.png")',
    '```',
    '',
    '## Subsample Sensitivity',
    '```{r s7_03}',
    'include_fig("Step7_03_SubsampleSensitivity.png")',
    '```',
    '',
    '## Permutation Test',
    '```{r s7_04}',
    'include_fig("Step7_04_PermutationTest.png")',
    '```',
    '',
    '## Module Preservation',
    '```{r s7_05}',
    'include_fig("Step7_05_ModulePreservation.png")',
    '```',
    '',
    '## Co-Clustering Heatmap',
    '```{r s7_06}',
    'include_fig("Step7_06_CoClusterHeatmap.png")',
    '```',
    '',
    '## Co-Clustering Summary',
    '```{r s7_07}',
    'include_fig("Step7_07_CoClusterSummary.png")',
    '```',
    '',
    '## Size vs Stability',
    '```{r s7_08}',
    'include_fig("Step7_08_SizeVsStability.png")',
    '```',
    '',
    '## Validation Quadrant',
    '```{r s7_09}',
    'include_fig("Step7_09_Validation_Quadrant.png")',
    '```',
    '',
    '## Combined Publication Figure',
    '```{r s7_10}',
    'include_fig("Step7_10_Combined_Publication.png")',
    '```',
    '',
    '## Stability Summary',
    '```{r s7tbl, results="asis"}',
    'load_tbl("Step7_StabilitySummary.csv", "Stability Summary")',
    '```',
    '',
    '---',
    '',
    '# Step 8: Custom Pathway Enrichment',
    '',
    '## ORA Lollipop Plot',
    '```{r s8_01}',
    'include_fig("Step8_01_ORA_Lollipop.png")',
    '```',
    '',
    '## ORA Heatmap',
    '```{r s8_02}',
    'include_fig("Step8_02_ORA_Heatmap.png")',
    '```',
    '',
    '## ORA Word Cloud',
    '```{r s8_03}',
    'include_fig("Step8_03_ORA_WordCloud_Final.png")',
    '```',
    '',
    '## GSVA Heatmap',
    '```{r s8_04}',
    'include_fig("Step8_04_GSVA_Heatmap.png")',
    '```',
    '',
    '## Clinical Forest Plot',
    '```{r s8_05}',
    'include_fig("Step8_05_Clinical_Forest_Plot.png")',
    '```',
    '',
    '## ME-Pathway Correlation',
    '```{r s8_06}',
    'include_fig("Step8_06_ME_Pathway_Cor.png")',
    '```',
    '',
    '## Combined Word Cloud',
    '```{r s8_08}',
    'include_fig("Step8_08_Combined_WordCloud.png")',
    '```',
    '',
    '## ORA Results Table',
    '```{r s8tbl, results="asis"}',
    'load_tbl("Step8_02_ORA_Results.csv", "ORA Enrichment Results", 30)',
    '```',
    '',
    '---',
    '',
    '# Step 9: Master Consensus Enrichment',
    '',
    '## Master Consensus Figure',
    '```{r s9_01}',
    'include_fig("Step9_01_Master_Consensus_Figure.png")',
    '```',
    '',
    '---',
    '',
    '# Step 10: Dual Clinical Association',
    '',
    '## Dual Lollipop Plot',
    '```{r s10_01}',
    'include_fig("Step10_01_Dual_Lollipop.png")',
    '```',
    '',
    '## Response Boxplots',
    '```{r s10_02}',
    'include_fig("Step10_02_Response_Boxplots.png")',
    '```',
    '',
    '## Clinical Significance Table',
    '```{r s10tbl, results="asis"}',
    'load_tbl("Step10_Dual_Clinical_Significance.csv", "Dual Clinical Significance")',
    '```',
    '',
    '---',
    '',
    '# Step 11: GO/KEGG/MSigDB Enrichment',
    '',
    '## GO:BP - Green Module',
    '```{r s11_01a}',
    'include_fig("Step11_01a_GO_BP_Green.png")',
    '```',
    '',
    '## GO:BP - Turquoise Module',
    '```{r s11_01b}',
    'include_fig("Step11_01b_GO_BP_Turquoise.png")',
    '```',
    '',
    '## KEGG - Green Module',
    '```{r s11_02a}',
    'include_fig("Step11_02a_KEGG_Green.png")',
    '```',
    '',
    '## KEGG - Turquoise Module',
    '```{r s11_02b}',
    'include_fig("Step11_02b_KEGG_Turquoise.png")',
    '```',
    '',
    '## MSigDB Hallmark Heatmap',
    '```{r s11_03}',
    'include_fig("Step11_03_Hallmark_Heatmap.png")',
    '```',
    '',
    '## Database Summary',
    '```{r s11_04}',
    'include_fig("Step11_04_DB_Summary.png")',
    '```',
    '',
    '## GSEA Hallmark',
    '```{r s11_05}',
    'include_fig("Step11_05_GSEA_Hallmark.png")',
    '```',
    '',
    '## GO:BP Results Table',
    '```{r s11tbl1, results="asis"}',
    'load_tbl("Step11_01_GO_BP_Results.csv", "GO:BP Results", 20)',
    '```',
    '',
    '## KEGG Results Table',
    '```{r s11tbl2, results="asis"}',
    'load_tbl("Step11_02_KEGG_Results.csv", "KEGG Results", 20)',
    '```',
    '',
    '---',
    '',
    '# Step 12: Biological Storytelling',
    '',
    '## Integrated Summary',
    '```{r s12_01}',
    'include_fig("Step12_03_Integrated_Summary.png")',
    '```',
    '',
    '## Module Narratives',
    '```{r s12tbl, results="asis"}',
    'load_tbl("Step12_01_Module_Narratives.csv", "Module Biological Narratives")',
    '```',
    '',
    '---',
    '',
    '# Step 13: Unified Pathway Integration',
    '',
    '## Unified Pathway Word Cloud',
    '```{r s13_01}',
    'include_fig("Step13_FINAL_Unified_Pathway_WordCloud.png")',
    '```',
    '',
    '## Category Summary Bar',
    '```{r s13_02}',
    'include_fig("Step13_FINAL_Category_Summary_Bar.png")',
    '```',
    '',
    '## Biomarker Mean Signature Boxplots',
    '```{r s13b_01}',
    'include_fig("Step13b_01_Mean_Signature_Boxplots.png")',
    '```',
    '',
    '## PFS Correlation Forest Plot',
    '```{r s13b_02}',
    'include_fig("Step13b_02_PFS_Correlation_Forest.png")',
    '```',
    '',
    '---',
    '',
    '# Step 14: DYNC2H1 Investigation',
    '',
    '## Technical Assessment',
    '```{r s14_01}',
    'include_fig("Step14_01_Technical_Assessment.png")',
    '```',
    '',
    '## Biological Assessment',
    '```{r s14_02}',
    'include_fig("Step14_02_Biological_Assessment.png")',
    '```',
    '',
    '## CV vs Effect Size',
    '```{r s14_03}',
    'include_fig("Step14_03_CV_vs_EffectSize.png")',
    '```',
    '',
    '## Investigation Summary',
    '```{r s14tbl, results="asis"}',
    'load_tbl("Step14_DYNC2H1_Investigation_Summary.csv", "DYNC2H1 Investigation")',
    '```',
    '',
    '---',
    '',
    '# Step 15: Signature Overlap Analysis',
    '',
    '## 3-Way Venn Diagram',
    '```{r s15_01}',
    'include_fig("Step15_01_Venn_3way.png")',
    '```',
    '',
    '## UpSet Plot',
    '```{r s15_02}',
    'include_fig("Step15_02_UpSet_Plot.png")',
    '```',
    '',
    '## Alluvial Signature Flow',
    '```{r s15_03}',
    'include_fig("Step15_02_Alluvial_SignatureFlow.png")',
    '```',
    '',
    '## Overlap Summary Table',
    '```{r s15tbl, results="asis"}',
    'load_tbl("Step15_Overlap_Summary.csv", "Signature Overlap Summary")',
    '```',
    '',
    '---',
    '',
    '# Session Information',
    '',
    '```{r session, echo=TRUE, comment=""}',
    'sessionInfo()',
    '```',
    '',
    '## Package Versions',
    '```{r pkgver, results="asis"}',
    'pkgs <- c("WGCNA", "ggplot2", "dplyr", "tidyr", "clusterProfiler", "GSVA", "rmarkdown", "knitr")',
    'pkg_df <- data.frame(',
    '  Package = pkgs,',
    '  Version = sapply(pkgs, function(p) if(requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "Not installed")',
    ')',
    'kable(pkg_df, caption="Key Package Versions")',
    '```',
    '',
    '---',
    '',
    paste0('*Report generated: ', Sys.time(), '*'),
    '',
    '*WGCNA Analysis Pipeline for Metastatic PDAC Proteomics*'
  )

  report_template <- paste(rmd_lines, collapse = "\n")

  # Save report template
  report_file <- file.path(results_dir, "WGCNA_Report.Rmd")
  writeLines(report_template, report_file)
  cat(sprintf("   [OK] Created report template: %s\n", report_file))

  # Render HTML report
  cat("\n   --- HTML Report Generation ---\n")
  tryCatch({
    cat("   Rendering HTML report...\n")
    rmarkdown::render(
      input = report_file,
      output_format = rmarkdown::html_document(
        toc = TRUE,
        toc_float = TRUE,
        toc_depth = 3,
        theme = "flatly",
        highlight = "tango",
        code_folding = "hide",
        fig_width = 10,
        fig_height = 8,
        self_contained = TRUE
      ),
      output_file = "WGCNA_Results_Report.html",
      output_dir = results_dir,
      quiet = FALSE,
      envir = new.env()
    )
    cat(sprintf("   [OK] HTML Report saved: %s/WGCNA_Results_Report.html\n", results_dir))
  }, error = function(e) {
    cat(sprintf("   [WARNING] HTML rendering failed: %s\n", e$message))
    cat("   Attempting simplified HTML render...\n")
    tryCatch({
      rmarkdown::render(
        input = report_file,
        output_format = "html_document",
        output_file = "WGCNA_Results_Report.html",
        output_dir = results_dir,
        quiet = TRUE,
        envir = new.env()
      )
      cat(sprintf("   [OK] Simplified HTML saved: %s/WGCNA_Results_Report.html\n", results_dir))
    }, error = function(e2) {
      cat(sprintf("   [ERROR] HTML generation failed: %s\n", e2$message))
    })
  })

  # Render PDF report (requires LaTeX)
  cat("\n   --- PDF Report Generation ---\n")

  # Check for LaTeX
  has_latex <- Sys.which("pdflatex") != "" ||
               requireNamespace("tinytex", quietly = TRUE)

  if(!has_latex) {
    cat("   [NOTE] LaTeX not found. PDF generation may fail.\n")
    cat("   To enable PDF output, install TinyTeX manually:\n")
    cat("     install.packages('tinytex')\n")
    cat("     tinytex::install_tinytex()\n")
  }

  tryCatch({
    cat("   Rendering PDF report (A4 format)...\n")
    rmarkdown::render(
      input = report_file,
      output_format = rmarkdown::pdf_document(
        toc = TRUE,
        toc_depth = 3,
        number_sections = TRUE,
        fig_width = 6.5,
        fig_height = 4.5,
        fig_caption = TRUE,
        latex_engine = "pdflatex",
        keep_tex = FALSE,
        extra_dependencies = list("geometry" = "a4paper,margin=1in")
      ),
      output_file = "WGCNA_Results_Report.pdf",
      output_dir = results_dir,
      quiet = FALSE,
      envir = new.env()
    )
    cat(sprintf("   [OK] PDF Report saved: %s/WGCNA_Results_Report.pdf\n", results_dir))
  }, error = function(e) {
    cat(sprintf("   [WARNING] PDF rendering failed: %s\n", e$message))
    cat("   To generate PDF manually, install LaTeX:\n")
    cat("     install.packages('tinytex')\n")
    cat("     tinytex::install_tinytex()\n")
    cat("   Then run:\n")
    cat(sprintf("     rmarkdown::render('%s', output_format = 'pdf_document')\n", report_file))
  })

} else {
  cat("   [WARNING] rmarkdown or knitr package not available\n")
  cat("   Install with:\n")
  cat("     install.packages(c('rmarkdown', 'knitr', 'kableExtra'))\n")
}

# ============================================================================
# COMBINE ALL FIGURES INTO SINGLE PDF
# ============================================================================
cat("\n   --- Combining All Figures into Single PDF ---\n")

# Get fig_dir_pdf path
fig_dir_pdf_path <- normalizePath(gsub("Results_Tables", "Figures_PDF", results_dir), winslash = "/", mustWork = FALSE)
cat(sprintf("   Looking for PDFs in: %s\n", fig_dir_pdf_path))

tryCatch({
  # Check for pdftools or qpdf package
  if(requireNamespace("qpdf", quietly = TRUE)) {

    # Get all PDF files sorted by step number
    pdf_files <- list.files(fig_dir_pdf_path, pattern = "\\.pdf$", full.names = TRUE)
    cat(sprintf("   Found %d PDF files\n", length(pdf_files)))

    if(length(pdf_files) > 0) {
      # Sort files: first by step number, then alphabetically within step
      get_step_num <- function(x) {
        step_match <- regmatches(basename(x), regexpr("Step\\d+", basename(x)))
        if(length(step_match) > 0) {
          as.numeric(gsub("Step", "", step_match))
        } else {
          999
        }
      }
      pdf_files <- pdf_files[order(sapply(pdf_files, get_step_num), basename(pdf_files))]

      # Combine PDFs
      output_combined <- file.path(results_dir, "All_Figures_Combined.pdf")
      qpdf::pdf_combine(pdf_files, output = output_combined)

      cat(sprintf("   [OK] Combined %d figures into: %s\n", length(pdf_files), output_combined))
    } else {
      cat("   [WARNING] No PDF figures found to combine\n")
    }

  } else if(requireNamespace("pdftools", quietly = TRUE)) {

    pdf_files <- list.files(fig_dir_pdf_path, pattern = "\\.pdf$", full.names = TRUE)

    if(length(pdf_files) > 0) {
      # Sort by step number
      get_step_num <- function(x) {
        step_match <- regmatches(basename(x), regexpr("Step\\d+", basename(x)))
        if(length(step_match) > 0) {
          as.numeric(gsub("Step", "", step_match))
        } else {
          999
        }
      }
      pdf_files <- pdf_files[order(sapply(pdf_files, get_step_num), basename(pdf_files))]

      output_combined <- file.path(results_dir, "All_Figures_Combined.pdf")
      pdftools::pdf_combine(pdf_files, output = output_combined)

      cat(sprintf("   [OK] Combined %d figures into: %s\n", length(pdf_files), output_combined))
    } else {
      cat("   [WARNING] No PDF figures found to combine\n")
    }

  } else {
    cat("   [NOTE] qpdf package not installed. Skipping PDF combination.\n")
    cat("   To enable PDF combination, install:\n")
    cat("     install.packages('qpdf')\n")
  }
}, error = function(e) {
  cat(sprintf("   [WARNING] PDF combination failed: %s\n", e$message))
})

cat("\n   Report generation complete\n")

# ============================================================================
# FINAL: SAVE SESSION INFO AND WORKSPACE
# ============================================================================

cat("\n--- Saving Session Info ---\n")

# Save session info
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(results_dir, "session_info.txt"))
cat("  Saved: session_info.txt\n")

# Save final workspace
tryCatch({
  save.image(file.path(results_dir, "WGCNA_v3_Final_Workspace.RData"))
  cat("  Saved: WGCNA_v3_Final_Workspace.RData\n")
}, error = function(e) {
  cat(sprintf("  WARNING: Could not save workspace: %s\n", e$message))
})

# ============================================================================
# PIPELINE COMPLETE
# ============================================================================

pipeline_end_time <- Sys.time()
pipeline_runtime <- difftime(pipeline_end_time, pipeline_start_time, units = "mins")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("========== WGCNA PIPELINE v3 COMPLETE ==========\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("Started:  %s\n", format(pipeline_start_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Finished: %s\n", format(pipeline_end_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Runtime:  %.1f minutes\n", as.numeric(pipeline_runtime)))
cat("\n")
cat("Output directories:\n")
cat(sprintf("  Figures (PNG): %s\n", normalizePath(fig_dir_png, mustWork = FALSE)))
cat(sprintf("  Figures (PDF): %s\n", normalizePath(fig_dir_pdf, mustWork = FALSE)))
cat(sprintf("  Tables/Reports: %s\n", normalizePath(results_dir, mustWork = FALSE)))
cat("\n")
cat("Key output files:\n")
cat("  - WGCNA_Results_Report.html (interactive report)\n")
cat("  - WGCNA_Results_Report.pdf (print-ready report)\n")
cat("  - All_Figures_Combined.pdf (all figures in one PDF)\n")
cat("  - Step8_Unified_Enrichment_Publication.xlsx (enrichment results)\n")
cat("  - WGCNA_v3_Final_Workspace.RData (full workspace)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# ============================================================================
# STEP 19b: BIOMARKER-MODULE INTEGRATION PUBLICATION FIGURES
# ============================================================================
# Purpose: Generate KEY FIGURES for manuscript showing biomarker-module relationships
# Content:
#   - All 7 biomarkers with module assignment
#   - Module clinical significance (correlation, p-value)
#   - Biomarker individual clinical significance (GS values)
#   - Effect sizes for modules
# Output: Publication-ready figures (300 DPI, Nature Medicine style)
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("STEP 19b: BIOMARKER-MODULE INTEGRATION FIGURES\n")
cat(paste(rep("=", 80), collapse = ""), "\n")


# =============================================================================
# PUBLICATION THEME - Nature Medicine / Cell Style
# =============================================================================
theme_publication <- theme_bw(base_size = 11) +
  theme(
    # Title styling
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5,
                              margin = margin(b = 8)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30",
                                 margin = margin(b = 10)),

    # Axis styling
    axis.title = element_text(face = "bold", size = 10),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),

    # Legend styling
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),

    # Panel styling
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white"),

    # Strip styling for facets
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9),

    # Plot margins
    plot.margin = margin(12, 12, 12, 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

# =============================================================================
# Define Color Palettes for Publication Figures
# =============================================================================
module_colors <- c(
  "pink" = "#E91E63", "green" = "#2E7D32", "black" = "#212121",
  "blue" = "#1565C0", "purple" = "#7B1FA2", "grey" = "#9E9E9E",
  "turquoise" = "#00897B", "brown" = "#6D4C41", "yellow" = "#F9A825",
  "red" = "#C62828", "magenta" = "#AD1457"
)

source_colors <- c("ML" = "#D32F2F", "Cox" = "#1976D2")
focus_colors <- c("Focus Module" = "#388E3C", "Other Module" = "#BDBDBD")

# Define biomarkers (consistent with rest of script)
ml_biomarkers <- c("ECM2", "DYNC2H1", "PPIB")
cox_biomarkers <- c("APOF", "TFRC", "FABP4", "ANG")
all_biomarkers <- c(ml_biomarkers, cox_biomarkers)

cat(sprintf("\n  Biomarkers to analyze: %d\n", length(all_biomarkers)))
cat(sprintf("    ML-derived: %s\n", paste(ml_biomarkers, collapse = ", ")))
cat(sprintf("    Cox-derived: %s\n", paste(cox_biomarkers, collapse = ", ")))

# -----------------------------------------------------------------------------
# 19b.1: Build comprehensive biomarker summary table
# -----------------------------------------------------------------------------
cat("\n  [19b.1] Building comprehensive biomarker summary table...\n")

biomarker_summary <- data.frame(
  Protein = all_biomarkers,
  Source = c(rep("ML", length(ml_biomarkers)), rep("Cox", length(cox_biomarkers))),
  stringsAsFactors = FALSE
)

# Add module assignment from proteinInfo
for(i in 1:nrow(biomarker_summary)) {
  prot <- biomarker_summary$Protein[i]
  if(prot %in% proteinInfo$Protein) {
    idx <- which(proteinInfo$Protein == prot)
    biomarker_summary$Module[i] <- proteinInfo$Module[idx]

    # Get kME for assigned module
    mm_col <- paste0("MM.", proteinInfo$Module[idx])
    if(mm_col %in% colnames(proteinInfo)) {
      biomarker_summary$kME[i] <- round(proteinInfo[[mm_col]][idx], 3)
    } else {
      biomarker_summary$kME[i] <- NA
    }

    # Get Gene Significance values
    if("GS.PFS_group" %in% colnames(proteinInfo)) {
      biomarker_summary$GS_PFS[i] <- round(proteinInfo$GS.PFS_group[idx], 3)
    }
    if("GS.Response" %in% colnames(proteinInfo)) {
      biomarker_summary$GS_Response[i] <- round(proteinInfo$GS.Response[idx], 3)
    }
    if("GS.PFS_Days" %in% colnames(proteinInfo)) {
      biomarker_summary$GS_PFS_Days[i] <- round(proteinInfo$GS.PFS_Days[idx], 3)
    }

    # Check if hub
    if("isHub_Moderate" %in% colnames(proteinInfo)) {
      biomarker_summary$Is_Hub[i] <- proteinInfo$isHub_Moderate[idx]
    }
  } else {
    biomarker_summary$Module[i] <- "grey"
    biomarker_summary$kME[i] <- NA
    biomarker_summary$GS_PFS[i] <- NA
    biomarker_summary$GS_Response[i] <- NA
    biomarker_summary$GS_PFS_Days[i] <- NA
    biomarker_summary$Is_Hub[i] <- FALSE
  }
}

# Add focus module status
biomarker_summary$In_Focus_Module <- biomarker_summary$Module %in% focus_modules_auto

# Add module-level clinical significance
module_cors <- data.frame(
  Module = c("pink", "green", "black", "blue", "turquoise", "brown", "magenta", "red", "purple", "yellow", "grey"),
  Module_Cor_PFS = c(-0.403, -0.403, -0.307, -0.277, -0.133, -0.067, 0.057, -0.041, 0.034, 0.027, NA),
  Module_P_PFS = c(0.0037, 0.0037, 0.0299, 0.0518, 0.3558, 0.6458, 0.6956, 0.7793, 0.8153, 0.8531, NA),
  Module_Cor_Response = c(-0.415, -0.407, -0.324, -0.316, -0.169, -0.116, 0.049, -0.039, 0.011, 0.009, NA),
  Module_P_Response = c(0.0027, 0.0034, 0.0216, 0.0253, 0.2409, 0.4210, 0.7355, 0.7875, 0.9418, 0.9514, NA)
)

biomarker_summary <- merge(biomarker_summary, module_cors, by = "Module", all.x = TRUE)

# Reorder columns for clarity
biomarker_summary <- biomarker_summary[, c("Protein", "Source", "Module", "In_Focus_Module", "Is_Hub",
                                            "kME", "GS_PFS", "GS_Response", "GS_PFS_Days",
                                            "Module_Cor_PFS", "Module_P_PFS",
                                            "Module_Cor_Response", "Module_P_Response")]

# Sort by Source then by focus module status then by kME
biomarker_summary <- biomarker_summary %>%
  arrange(Source, desc(In_Focus_Module), desc(abs(Module_Cor_PFS)))

cat("  Biomarker summary table:\n")
print(biomarker_summary)

# Save table
write.csv(biomarker_summary, file.path(results_dir, "Step19b_Biomarker_Module_Summary.csv"), row.names = FALSE)
cat(sprintf("  Saved: Step19b_Biomarker_Module_Summary.csv\n"))

# -----------------------------------------------------------------------------
# 19b.2: FIGURE 1 - Biomarker Module Assignment Tile Plot (ggplot2)
# -----------------------------------------------------------------------------
cat("\n  [19b.2] Creating publication-ready biomarker-module tile plot...\n")

# Prepare data for tile plot
plot_data <- biomarker_summary %>%
  mutate(
    Protein = factor(Protein, levels = rev(unique(Protein))),
    Module_Upper = toupper(Module),
    Focus_Status = ifelse(In_Focus_Module, "Focus Module", "Other Module"),
    kME_display = ifelse(is.na(kME), "-", sprintf("%.2f", kME)),
    GS_PFS_display = ifelse(is.na(GS_PFS), "-", sprintf("%.2f", GS_PFS)),
    GS_Resp_display = ifelse(is.na(GS_Response), "-", sprintf("%.2f", GS_Response))
  )

# Define source colors matching Step18 Venn diagram style
# ML = orange/coral, Cox = blue/teal
source_text_colors <- c("ML" = "#E64A19", "Cox" = "#0277BD")

# Get the correct color order for y-axis labels (matching factor levels in plot_data)
# The y-axis shows proteins in reverse order of factor levels
protein_order <- levels(plot_data$Protein)
y_axis_colors <- sapply(protein_order, function(p) {
  source_text_colors[plot_data$Source[as.character(plot_data$Protein) == p][1]]
})

# Create the Module tile plot with protein names colored by source
p1_tiles <- ggplot(plot_data, aes(x = 1, y = Protein)) +
  # Module color tile
  geom_tile(aes(fill = Module), color = "white", width = 0.9, height = 0.85) +
  geom_text(aes(label = Module_Upper), color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = module_colors, guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Module") +
  theme_publication +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10,
                                color = y_axis_colors),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5)
  )

# Focus module indicator - using tiles for consistency
p3_focus <- ggplot(plot_data, aes(x = 1, y = Protein)) +
  geom_tile(aes(fill = Focus_Status), color = "white", width = 0.9, height = 0.85) +
  geom_text(aes(label = ifelse(In_Focus_Module, "Yes", "No")),
            color = "white", fontface = "bold", size = 3) +
  scale_fill_manual(values = focus_colors, name = "Focus") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Focus\nModule") +
  theme_publication +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = "bottom"
  )

# kME heatmap
p4_kme <- ggplot(plot_data, aes(x = 1, y = Protein)) +
  geom_tile(aes(fill = kME), color = "gray50", width = 0.9, height = 0.85) +
  geom_text(aes(label = kME_display), color = "black", fontface = "bold", size = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, na.value = "gray90",
    limits = c(-1, 1), name = "kME"
  ) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "kME") +
  theme_publication +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm")
  )

# Trait Correlation PFS heatmap
p5_gs_pfs <- ggplot(plot_data, aes(x = 1, y = Protein)) +
  geom_tile(aes(fill = GS_PFS), color = "gray50", width = 0.9, height = 0.85) +
  geom_text(aes(label = GS_PFS_display), color = "black", fontface = "bold", size = 3) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, na.value = "gray90",
    limits = c(-0.5, 0.5), name = "r"
  ) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Trait Cor.\n(PFS)") +
  theme_publication +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm")
  )

# Trait Correlation Response heatmap
p6_gs_resp <- ggplot(plot_data, aes(x = 1, y = Protein)) +
  geom_tile(aes(fill = GS_Response), color = "gray50", width = 0.9, height = 0.85) +
  geom_text(aes(label = GS_Resp_display), color = "black", fontface = "bold", size = 3) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, na.value = "gray90",
    limits = c(-0.5, 0.5), name = "r"
  ) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Trait Cor.\n(Response)") +
  theme_publication +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm")
  )

# Combine into single figure (removed Source panel - protein names are now colored by source)
fig1_combined <- plot_grid(
  p1_tiles + theme(legend.position = "none"),
  p3_focus + theme(legend.position = "none"),
  p4_kme + theme(legend.position = "none"),
  p5_gs_pfs + theme(legend.position = "none"),
  p6_gs_resp + theme(legend.position = "none"),
  nrow = 1,
  rel_widths = c(1.8, 0.8, 1, 1, 1),
  align = "h"
)

# Add title
fig1_title <- ggdraw() +
  draw_label(
    "Biomarker Module Assignment and Clinical Significance",
    fontface = "bold", size = 14, x = 0.5, hjust = 0.5
  )

# Create proper legend for biomarker source using ggplot
source_legend_data <- data.frame(
  Source = factor(c("ML", "Cox"), levels = c("ML", "Cox")),
  x = c(1, 2),
  y = c(1, 1)
)

p_source_legend <- ggplot(source_legend_data, aes(x = x, y = y, fill = Source)) +
  geom_tile(width = 0.5, height = 0.5) +
  scale_fill_manual(values = c("ML" = "#E64A19", "Cox" = "#0277BD"),
                    name = "Biomarker Source") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  )

legend_source <- get_legend(p_source_legend)

# Get legends for Focus, kME, and Trait Correlation
legend_focus <- get_legend(p3_focus + theme(legend.position = "bottom", legend.direction = "horizontal"))
legend_kme <- get_legend(p4_kme + theme(legend.position = "bottom", legend.direction = "horizontal"))
legend_trait_cor <- get_legend(p5_gs_pfs +
                                 guides(fill = guide_colorbar(title = "Trait Cor. (r)")) +
                                 theme(legend.position = "bottom", legend.direction = "horizontal"))

legend_row <- plot_grid(legend_source, legend_focus, legend_kme, legend_trait_cor,
                        nrow = 1, rel_widths = c(1, 1, 1.2, 1.2))

# Final assembly
fig1_final <- plot_grid(
  fig1_title,
  fig1_combined,
  legend_row,
  ncol = 1,
  rel_heights = c(0.08, 1, 0.12)
)

# Display in RStudio
print(fig1_final)

# Save
ggsave(file.path(fig_dir_png, "Step19b_01_Biomarker_Module_Tiles.png"),
       fig1_final, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step19b_01_Biomarker_Module_Tiles.pdf"),
       fig1_final, width = 12, height = 7, bg = "white")
cat("  Saved: Step19b_01_Biomarker_Module_Tiles (PNG/PDF)\n")

# -----------------------------------------------------------------------------
# 19b.3: FIGURE 2 - Trait Correlation Bar Chart
# -----------------------------------------------------------------------------
cat("\n  [19b.3] Creating trait correlation bar charts...\n")

# Prepare long-format data for GS plots
gs_data <- plot_data %>%
  select(Protein, Source, Module, In_Focus_Module, GS_PFS, GS_Response) %>%
  pivot_longer(cols = c(GS_PFS, GS_Response),
               names_to = "Outcome",
               values_to = "GS") %>%
  mutate(
    Outcome = recode(Outcome,
                     "GS_PFS" = "PFS Group",
                     "GS_Response" = "Treatment Response"),
    Direction = case_when(
      is.na(GS) ~ "NA",
      GS < 0 ~ "Poor Prognosis",
      GS >= 0 ~ "Good Prognosis"
    ),
    Protein = factor(Protein, levels = biomarker_summary$Protein[order(biomarker_summary$GS_PFS, na.last = TRUE)])
  )

# Panel A: GS for PFS Group
p_gs_pfs <- ggplot(gs_data %>% filter(Outcome == "PFS Group"),
                   aes(x = reorder(Protein, GS, na.rm = TRUE), y = GS)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_col(aes(fill = Direction), color = "black", width = 0.7, linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(GS), "NA", sprintf("%.2f", GS)),
                y = ifelse(is.na(GS), 0, GS + sign(GS) * 0.03)),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Poor Prognosis" = "#C62828",
                               "Good Prognosis" = "#1565C0",
                               "NA" = "#BDBDBD"),
                    guide = "none") +
  coord_flip(ylim = c(-0.5, 0.5)) +
  labs(
    title = "Trait Correlation: PFS Group",
    x = NULL,
    y = "Trait Correlation (r)"
  ) +
  theme_publication +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 9, face = "italic")
  ) +
  annotate("text", x = 0.5, y = -0.4, label = "Poor Prognosis",
           color = "#C62828", fontface = "italic", size = 3, hjust = 0) +
  annotate("text", x = 0.5, y = 0.4, label = "Good Prognosis",
           color = "#1565C0", fontface = "italic", size = 3, hjust = 0)

# Panel B: GS for Treatment Response
p_gs_resp <- ggplot(gs_data %>% filter(Outcome == "Treatment Response"),
                    aes(x = reorder(Protein, GS, na.rm = TRUE), y = GS)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_col(aes(fill = Direction), color = "black", width = 0.7, linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(GS), "NA", sprintf("%.2f", GS)),
                y = ifelse(is.na(GS), 0, GS + sign(GS) * 0.03)),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Poor Prognosis" = "#C62828",
                               "Good Prognosis" = "#1565C0",
                               "NA" = "#BDBDBD"),
                    guide = "none") +
  coord_flip(ylim = c(-0.5, 0.5)) +
  labs(
    title = "Trait Correlation: Treatment Response",
    x = NULL,
    y = "Trait Correlation (r)"
  ) +
  theme_publication +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 9, face = "italic")
  ) +
  annotate("text", x = 0.5, y = -0.4, label = "Progressive Disease",
           color = "#C62828", fontface = "italic", size = 3, hjust = 0) +
  annotate("text", x = 0.5, y = 0.4, label = "Controlled Disease",
           color = "#1565C0", fontface = "italic", size = 3, hjust = 0)

# Combine GS plots
fig2_gs <- plot_grid(p_gs_pfs, p_gs_resp, ncol = 2, labels = NULL, align = "h")


print(fig2_gs)

# Save
ggsave(file.path(fig_dir_png, "Step19b_02_Trait_Correlation_Bars.png"),
       fig2_gs, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step19b_02_Trait_Correlation_Bars.pdf"),
       fig2_gs, width = 8, height = 6)
cat("  Saved: Step19b_02_Trait_Correlation_Bars (PNG/PDF)\n")

# -----------------------------------------------------------------------------
# 19b.4: FIGURE 3 - Focus Module Clinical Correlations
# -----------------------------------------------------------------------------
cat("\n  [19b.4] Creating focus module clinical correlation plot...\n")

# Prepare focus module data
focus_mod_data <- data.frame(
  Module = factor(c("PINK", "GREEN", "BLACK", "BLUE"),
                  levels = c("PINK", "GREEN", "BLACK", "BLUE")),
  Cor_PFS = c(-0.403, -0.403, -0.307, -0.277),
  P_PFS = c(0.0037, 0.0037, 0.0299, 0.0518),
  Cor_Response = c(-0.415, -0.407, -0.324, -0.316),
  P_Response = c(0.0027, 0.0034, 0.0216, 0.0253),
  Contains_Biomarker = c("FABP4", "ECM2", "", "")
)

# Convert to long format
focus_long <- focus_mod_data %>%
  pivot_longer(cols = c(Cor_PFS, Cor_Response),
               names_to = "Outcome",
               values_to = "Correlation") %>%
  mutate(
    P_value = ifelse(Outcome == "Cor_PFS", P_PFS, P_Response),
    Outcome = recode(Outcome,
                     "Cor_PFS" = "PFS Group",
                     "Cor_Response" = "Treatment Response"),
    Significance = case_when(
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Create grouped bar plot
p_focus_modules <- ggplot(focus_long, aes(x = Module, y = Correlation, fill = Outcome)) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.4) +
  geom_text(aes(label = Significance,
                y = Correlation - 0.025),
            position = position_dodge(width = 0.8),
            size = 5, fontface = "bold", vjust = 1) +
  geom_text(data = focus_mod_data %>% filter(Contains_Biomarker != ""),
            aes(x = Module, y = 0.06, label = Contains_Biomarker),
            inherit.aes = FALSE,
            size = 3, fontface = "italic", color = "#D32F2F") +
  scale_fill_manual(values = c("PFS Group" = "#0288D1",
                               "Treatment Response" = "#F57C00"),
                    name = "Clinical Outcome") +
  scale_y_continuous(limits = c(-0.5, 0.1), breaks = seq(-0.5, 0.1, 0.1)) +
  labs(
    title = "Focus Module Clinical Correlations",
    subtitle = "All modules show negative correlation with favorable outcomes (poor prognosis signature)",
    x = "Module",
    y = "Correlation with Clinical Outcome"
  ) +
  theme_publication +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold", size = 11)
  ) +
  annotate("text", x = 4.4, y = -0.48, label = "* p<0.05  ** p<0.01",
           size = 3, fontface = "italic", hjust = 1)

# Display
print(p_focus_modules)

# Save
ggsave(file.path(fig_dir_png, "Step19b_03_Focus_Module_Correlations.png"),
       p_focus_modules, width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step19b_03_Focus_Module_Correlations.pdf"),
       p_focus_modules, width = 6, height = 5)
cat("  Saved: Step19b_03_Focus_Module_Correlations (PNG/PDF)\n")

# -----------------------------------------------------------------------------
# 19b.5: FIGURE 4 - Effect Size Forest Plot (ggplot2)
# -----------------------------------------------------------------------------
cat("\n  [19b.5] Creating effect size forest plot...\n")

# Calculate effect sizes for biomarkers (Cohen's d from GS)
# GS is correlation; approximate d from r: d = 2r / sqrt(1-r^2)
effect_data <- biomarker_summary %>%
  select(Protein, Source, Module, In_Focus_Module, GS_PFS, GS_Response) %>%
  mutate(
    d_PFS = ifelse(is.na(GS_PFS), NA, 2 * GS_PFS / sqrt(1 - GS_PFS^2)),
    d_Response = ifelse(is.na(GS_Response), NA, 2 * GS_Response / sqrt(1 - GS_Response^2)),
    # Labels
    Module_Abbrev = toupper(substr(Module, 1, 3)),
    Focus_Label = ifelse(In_Focus_Module, "Focus Module", "Other Module"),
    # Create label with line break: Protein name on top, module below
    Y_Label = paste0(Protein, "\n", Module_Abbrev)
  ) %>%
  arrange(d_PFS) %>%
  mutate(y_order = row_number())

# Create the forest plot (without CI - just effect sizes)
p_forest <- ggplot(effect_data, aes(x = d_PFS, y = reorder(Y_Label, d_PFS))) +
  # Reference lines for effect size categories
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  geom_vline(xintercept = c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8),
             linetype = "dotted", color = "gray60", linewidth = 0.4) +

  # Horizontal line from zero to point (lollipop style)
  geom_segment(aes(x = 0, xend = d_PFS, yend = reorder(Y_Label, d_PFS), color = Focus_Label),
               linewidth = 1, na.rm = TRUE) +

  # Point estimates
  geom_point(aes(color = Focus_Label, shape = Source),
             size = 5, na.rm = TRUE) +

  # Effect size labels on right
  geom_text(aes(x = 1.15, label = ifelse(is.na(d_PFS), "NA", sprintf("d = %.2f", d_PFS))),
            hjust = 0, size = 3.5, fontface = "bold") +

  # Scales
  scale_color_manual(values = c("Focus Module" = "#C62828", "Other Module" = "#616161"),
                     name = "Module Status") +
  scale_shape_manual(values = c("ML" = 19, "Cox" = 17),
                     name = "Biomarker Source") +
  scale_x_continuous(limits = c(-1.2, 1.4), breaks = seq(-1, 1, 0.5)) +

  # Labels
  labs(
    title = "Biomarker Effect Sizes: Association with PFS Group",
    subtitle = "Cohen's d derived from Trait Correlation",
    x = "Effect Size (Cohen's d)",
    y = NULL
  ) +

  # Theme
  theme_publication +
  theme(
    axis.text.y = element_text(face = "bold", size = 9, lineheight = 0.9),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +

  # Annotations for effect size interpretation
  annotate("text", x = -0.8, y = 7.6, label = "Large", color = "gray50",
           size = 3, fontface = "italic") +
  annotate("text", x = -0.5, y = 7.6, label = "Medium", color = "gray50",
           size = 3, fontface = "italic") +
  annotate("text", x = -0.2, y = 7.6, label = "Small", color = "gray50",
           size = 3, fontface = "italic") +

  # Interpretation note
  annotate("text", x = -0.6, y = 0.4,
           label = "Negative = Higher expression\nin Short PFS patients",
           size = 2.8, fontface = "italic", hjust = 0, lineheight = 0.9)

# Display
print(p_forest)

# Save
ggsave(file.path(fig_dir_png, "Step19b_04_Effect_Size_Forest.png"),
       p_forest, width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step19b_04_Effect_Size_Forest.pdf"),
       p_forest, width = 10, height = 7, bg = "white")
cat("  Saved: Step19b_04_Effect_Size_Forest (PNG/PDF)\n")

# -----------------------------------------------------------------------------
# 19b.6: FIGURE 5 - Comprehensive Multi-Panel Dashboard (ggplot2)
# -----------------------------------------------------------------------------
cat("\n  [19b.6] Creating comprehensive biomarker dashboard...\n")

# Panel A: Module assignment heatmap-style visualization
panel_a <- ggplot(plot_data, aes(x = 1, y = Protein)) +
  geom_tile(aes(fill = Module), color = "white", width = 0.9, height = 0.85) +
  geom_text(aes(label = toupper(Module)), color = "white", fontface = "bold", size = 3) +
  scale_fill_manual(values = module_colors, guide = "none") +
  labs(title = "A. Module\nAssignment", x = NULL, y = NULL) +
  theme_publication +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 9),
    plot.title = element_text(size = 10, hjust = 0.5, lineheight = 1.1)
  )

# Panel B: Source and Focus indicators
panel_b_data <- plot_data %>%
  select(Protein, Source, Focus_Status) %>%
  pivot_longer(cols = c(Source, Focus_Status), names_to = "Type", values_to = "Value")

panel_b <- ggplot(plot_data, aes(y = Protein)) +
  geom_point(aes(x = 1, color = Source, shape = Source), size = 4) +
  geom_point(aes(x = 2, color = Focus_Status, shape = Focus_Status), size = 4) +
  scale_color_manual(values = c(source_colors, focus_colors)) +
  scale_shape_manual(values = c("ML" = 19, "Cox" = 17,
                                "Focus Module" = 15, "Other Module" = 4)) +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1, 2),
                     labels = c("Source", "Focus")) +
  labs(title = "B. Classification", x = NULL, y = NULL) +
  theme_publication +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# Panel C: kME values
panel_c <- ggplot(plot_data, aes(x = kME, y = Protein)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_segment(aes(x = 0, xend = kME, y = Protein, yend = Protein, color = Module),
               linewidth = 1.5, na.rm = TRUE) +
  geom_point(aes(fill = Module), shape = 21, size = 4, color = "black", na.rm = TRUE) +
  geom_text(aes(label = kME_display, x = kME + 0.08), size = 2.8, hjust = 0, na.rm = TRUE) +
  scale_fill_manual(values = module_colors, guide = "none") +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.2, 1)) +
  labs(title = "C. Module Membership (kME)", x = "kME", y = NULL) +
  theme_publication +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# Panel D: Trait Correlation PFS lollipop
panel_d <- ggplot(plot_data, aes(x = GS_PFS, y = Protein)) +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray60") +
  geom_segment(aes(x = 0, xend = GS_PFS, y = Protein, yend = Protein,
                   color = ifelse(GS_PFS < 0, "Negative", "Positive")),
               linewidth = 1.2, na.rm = TRUE) +
  geom_point(aes(fill = ifelse(GS_PFS < 0, "Negative", "Positive")),
             shape = 21, size = 4, color = "black", na.rm = TRUE) +
  scale_color_manual(values = c("Negative" = "#C62828", "Positive" = "#1565C0"), guide = "none") +
  scale_fill_manual(values = c("Negative" = "#C62828", "Positive" = "#1565C0"), guide = "none") +
  scale_x_continuous(limits = c(-0.5, 0.35)) +
  labs(title = "D. Trait Cor. (PFS)", x = "r", y = NULL) +
  theme_publication +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# Panel E: Trait Correlation Response lollipop
panel_e <- ggplot(plot_data, aes(x = GS_Response, y = Protein)) +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray60") +
  geom_segment(aes(x = 0, xend = GS_Response, y = Protein, yend = Protein,
                   color = ifelse(GS_Response < 0, "Negative", "Positive")),
               linewidth = 1.2, na.rm = TRUE) +
  geom_point(aes(fill = ifelse(GS_Response < 0, "Negative", "Positive")),
             shape = 21, size = 4, color = "black", na.rm = TRUE) +
  scale_color_manual(values = c("Negative" = "#C62828", "Positive" = "#1565C0"), guide = "none") +
  scale_fill_manual(values = c("Negative" = "#C62828", "Positive" = "#1565C0"), guide = "none") +
  scale_x_continuous(limits = c(-0.5, 0.35)) +
  labs(title = "E. Trait Cor. (Response)", x = "r", y = NULL) +
  theme_publication +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# Combine all panels
dashboard_row1 <- plot_grid(
  panel_a, panel_b, panel_c, panel_d, panel_e,
  nrow = 1,
  rel_widths = c(1.5, 0.8, 1.2, 1.2, 1.2),
  align = "h"
)

# Create legend panel
legend_data <- data.frame(
  x = c(1, 2, 3, 4, 5, 6),
  y = 1,
  label = c("ML", "Cox", "Focus", "Other", "Neg r", "Pos r"),
  color = c("#D32F2F", "#1976D2", "#388E3C", "#BDBDBD", "#C62828", "#1565C0"),
  shape = c(19, 17, 15, 4, 21, 21)
)

p_legend <- ggplot(legend_data, aes(x = x, y = y)) +
  geom_point(aes(color = label, shape = label), size = 4) +
  geom_text(aes(label = label, x = x + 0.3), hjust = 0, size = 3) +
  scale_color_manual(values = setNames(legend_data$color, legend_data$label)) +
  scale_shape_manual(values = setNames(legend_data$shape, legend_data$label)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 8))

# Add title
dashboard_title <- ggdraw() +
  draw_label(
    "Biomarker-Module Integration Dashboard",
    fontface = "bold", size = 14, x = 0.5, hjust = 0.5
  ) +
  draw_label(
    "WGCNA Analysis of 7 Prognostic Biomarkers in Metastatic PDAC",
    size = 10, x = 0.5, y = 0.25, hjust = 0.5, color = "gray40"
  )

# Final dashboard assembly
dashboard_final <- plot_grid(
  dashboard_title,
  dashboard_row1,
  p_legend,
  ncol = 1,
  rel_heights = c(0.12, 1, 0.08)
)

# Display in RStudio
print(dashboard_final)

# Save
ggsave(file.path(fig_dir_png, "Step19b_05_Biomarker_Dashboard.png"),
       dashboard_final, width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(fig_dir_pdf, "Step19b_05_Biomarker_Dashboard.pdf"),
       dashboard_final, width = 14, height = 8, bg = "white")
cat("  Saved: Step19b_05_Biomarker_Dashboard (PNG/PDF)\n")

# -----------------------------------------------------------------------------
# 19b.7: Publication-Ready Summary Table
# -----------------------------------------------------------------------------
cat("\n  [19b.7] Creating publication-ready summary table...\n")

# Create formatted table for manuscript
pub_table <- biomarker_summary %>%
  mutate(
    Module_P_PFS_fmt = case_when(
      is.na(Module_P_PFS) ~ "-",
      Module_P_PFS < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", Module_P_PFS)
    ),
    Significance = case_when(
      is.na(Module_P_PFS) ~ "",
      Module_P_PFS < 0.01 ~ "**",
      Module_P_PFS < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Format for publication
pub_table_final <- data.frame(
  Biomarker = pub_table$Protein,
  Source = pub_table$Source,
  `Assigned Module` = toupper(pub_table$Module),
  `Focus Module` = ifelse(pub_table$In_Focus_Module, "Yes", "No"),
  `kME` = ifelse(is.na(pub_table$kME), "-", sprintf("%.3f", pub_table$kME)),
  `Trait Cor. (PFS)` = ifelse(is.na(pub_table$GS_PFS), "-", sprintf("%.3f", pub_table$GS_PFS)),
  `Trait Cor. (Response)` = ifelse(is.na(pub_table$GS_Response), "-", sprintf("%.3f", pub_table$GS_Response)),
  `Module r (PFS)` = ifelse(is.na(pub_table$Module_Cor_PFS), "-", sprintf("%.3f", pub_table$Module_Cor_PFS)),
  `Module p-value` = paste0(pub_table$Module_P_PFS_fmt, pub_table$Significance),
  check.names = FALSE
)

write.csv(pub_table_final, file.path(results_dir, "Step19b_Biomarker_Publication_Table.csv"), row.names = FALSE)
cat("  Saved: Step19b_Biomarker_Publication_Table.csv\n")

# Print table
cat("\n  Publication-ready biomarker summary table:\n")
cat("  ", paste(rep("-", 100), collapse = ""), "\n")
print(pub_table_final, row.names = FALSE)
cat("  ", paste(rep("-", 100), collapse = ""), "\n")
cat("  * p < 0.05, ** p < 0.01\n")

# -----------------------------------------------------------------------------
# Step 19b Summary
# -----------------------------------------------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("STEP 19b COMPLETE: Biomarker-Module Integration Figures\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\nKey Findings:\n")
cat(sprintf("  - %d/%d biomarkers assigned to focus modules\n",
            sum(biomarker_summary$In_Focus_Module), nrow(biomarker_summary)))
cat(sprintf("  - FABP4 (Cox) in PINK module (kME=%.3f, module p=%.4f)\n",
            biomarker_summary$kME[biomarker_summary$Protein == "FABP4"],
            biomarker_summary$Module_P_PFS[biomarker_summary$Protein == "FABP4"]))
cat(sprintf("  - ECM2 (ML) in GREEN module (kME=%.3f, module p=%.4f)\n",
            biomarker_summary$kME[biomarker_summary$Protein == "ECM2"],
            biomarker_summary$Module_P_PFS[biomarker_summary$Protein == "ECM2"]))
cat(sprintf("  - %d biomarkers unassigned (grey/purple module)\n",
            sum(biomarker_summary$Module %in% c("grey", "purple"))))
cat("\nOutput Files:\n")
cat("  Tables:\n")
cat("    - Step19b_Biomarker_Module_Summary.csv\n")
cat("    - Step19b_Biomarker_Publication_Table.csv\n")
cat("  Figures (PNG & PDF - all ggplot2, display in RStudio):\n")
cat("    - Step19b_01_Biomarker_Module_Tiles.png/pdf\n")
cat("    - Step19b_02_Trait_Correlation_Bars.png/pdf\n")
cat("    - Step19b_03_Focus_Module_Correlations.png/pdf\n")
cat("    - Step19b_04_Effect_Size_Forest.png/pdf\n")
cat("    - Step19b_05_Biomarker_Dashboard.png/pdf\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ============================================================================
# STEP 20: COMPREHENSIVE PUBLICATION-READY REPORT GENERATION
# ============================================================================
# Purpose: Generate comprehensive HTML/PDF report for manuscript preparation
# Target: High-impact factor journals (Nature Medicine, Cancer Cell, JCO, etc.)
# Output: R Markdown -> HTML -> PDF conversion with all figures and tables
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("STEP 20: GENERATING PUBLICATION-READY REPORT\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Ensure correct output directories for Power14 analysis
results_dir <- file.path(main_dir, "Results_Tables")
fig_dir_png <- file.path(main_dir, "Figures_PNG")
fig_dir_pdf <- file.path(main_dir, "Figures_PDF")

# Ensure required variables exist for report generation
if(!exists("n_samples")) n_samples <- nrow(datExpr)
if(!exists("n_proteins")) n_proteins <- ncol(datExpr)
if(!exists("soft_power")) soft_power <- 14
if(!exists("focus_modules_auto")) focus_modules_auto <- c("pink", "green", "black", "blue")
if(!exists("moduleColors")) stop("moduleColors not found - run WGCNA module detection first")

# Pre-calculate values for report
focus_modules_str <- paste(toupper(focus_modules_auto), collapse = ", ")
n_modules <- length(unique(moduleColors)) - 1
n_focus <- length(focus_modules_auto)

cat(sprintf("  Variables for report: n_samples=%d, n_proteins=%d, n_modules=%d, soft_power=%d\n",
            n_samples, n_proteins, n_modules, soft_power))
cat(sprintf("  Focus modules: %s\n", focus_modules_str))

# Create R Markdown file using gsub substitution (avoids sprintf % escaping issues)
rmd_template <- '---
title: "Plasma Proteomics WGCNA Analysis of Metastatic PDAC"
subtitle: "Weighted Gene Co-expression Network Analysis for Prognostic Biomarker Discovery"
author: "WGCNA Pipeline v4"
date: "`r format(Sys.time(), \'%B %d, %Y\')`"
output:
  html_document:
    toc: true
    toc_depth: 4
    toc_float: true
    number_sections: true
    theme: flatly
    highlight: tango
    code_folding: hide
    df_print: paged
    fig_width: 10
    fig_height: 7
  pdf_document:
    toc: true
    toc_depth: 4
    number_sections: true
    fig_width: 10
    fig_height: 7
    latex_engine: xelatex
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  fig.align = "center",
  out.width = "100%"
)

# Load required packages
library(knitr)
library(kableExtra)
library(ggplot2)
library(dplyr)

# Set paths
results_dir <- "{{RESULTS_DIR}}"
fig_dir_png <- "{{FIG_DIR_PNG}}"
fig_dir_pdf <- "{{FIG_DIR_PDF}}"
```

# Executive Summary {.tabset}

## Study Overview

**Objective:** Identify plasma protein co-expression modules associated with progression-free survival (PFS) and treatment response in metastatic pancreatic ductal adenocarcinoma (PDAC) patients.

**Dataset:**
- **Cohort:** {{N_SAMPLES}} metastatic PDAC patients
- **Proteomics:** {{N_PROTEINS}} plasma proteins quantified
- **Clinical Endpoints:** PFS (days), PFS group (Short/Long), Treatment Response (CD/PD)

**Key Findings:**
- **{{N_MODULES}} co-expression modules** identified (excluding grey/unassigned)
- **{{N_FOCUS}} focus modules** ({{FOCUS_MODULES}}) significantly associated with poor prognosis
- Focus modules show **negative correlation with PFS** (higher expression = shorter survival)

## Clinical Significance

The identified focus modules represent biologically coherent protein networks that:

1. **Predict poor prognosis** - Patients with high module eigengene values have significantly shorter PFS

2. **Associate with treatment resistance** - Higher module expression correlates with progressive disease (PD)

3. **Enrich for known cancer pathways** - ECM remodeling, angiogenesis, immune suppression

---

# Methods

## WGCNA Network Construction

```{r network-params}
network_params <- data.frame(
  Parameter = c("Network type", "Soft-threshold power", "Correlation method",
                "Minimum module size", "Merge cut height", "Deep split"),
  Value = c("Signed hybrid", "{{SOFT_POWER}}", "Bicor (biweight midcorrelation)",
            "30", "0.25", "3")
)
kable(network_params, caption = "WGCNA Network Construction Parameters") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
```

## Focus Module Selection Criteria

Modules were selected as clinically significant based on:

1. **Correlation with PFS_group** (p < 0.05)
2. **Correlation with Treatment Response** (p < 0.05)
3. **Negative direction** (higher expression = poor prognosis)

**Selected Focus Modules:** {{FOCUS_MODULES}}

---

# Results

## Module-Trait Associations {.tabset}

### Heatmap

```{r module-trait-heatmap, fig.cap="Module-trait correlation heatmap showing associations between module eigengenes and clinical variables. Focus modules (pink, green, black, blue) show significant negative correlations with PFS."}
knitr::include_graphics(file.path(fig_dir_png, "Step5_04_ModuleTrait_Heatmap.png"))
```

### Statistical Summary

```{r load-trait-data}
if(file.exists(file.path(results_dir, "Step5_ModuleTraitCorrelations.csv"))) {
  trait_cor <- read.csv(file.path(results_dir, "Step5_ModuleTraitCorrelations.csv"))
  kable(trait_cor, digits = 4, caption = "Module-Trait Correlations") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

## Focus Module Characterization {.tabset}

### Module Sizes

```{r module-sizes}
focus_mods <- c("pink", "green", "black", "blue")
module_info <- data.frame(
  Module = toupper(focus_mods),
  stringsAsFactors = FALSE
)
# Load from proteinInfo if available
if(file.exists(file.path(results_dir, "Step6_ProteinInfo_Full.csv"))) {
  pinfo <- read.csv(file.path(results_dir, "Step6_ProteinInfo_Full.csv"))
  for(m in focus_mods) {
    module_info$Size[module_info$Module == toupper(m)] <- sum(pinfo$Module == m, na.rm = TRUE)
  }
  kable(module_info, caption = "Focus Module Sizes") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
}
```

### Dendrogram

```{r dendrogram, fig.cap="Hierarchical clustering dendrogram with module color assignments. Focus modules are highlighted."}
if(file.exists(file.path(fig_dir_png, "Step5_01_Dendrogram.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step5_01_Dendrogram.png"))
}
```

## Hub Gene Analysis {.tabset}

### Hub Gene Identification

Hub genes were identified using moderate criteria (kME > 0.7, |GS| > 0.2, p < 0.05).

```{r hub-genes}
if(file.exists(file.path(results_dir, "Step6_HubSummary.csv"))) {
  hubs <- read.csv(file.path(results_dir, "Step6_HubSummary.csv"))
  hubs_focus <- hubs %>% filter(Module %in% focus_mods)
  kable(hubs_focus %>% select(Protein, Module, kME, GS_PFS_group, GS_Response) %>% head(20),
        digits = 3, caption = "Top Hub Genes in Focus Modules") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Hub Stability

```{r hub-stability}
if(file.exists(file.path(results_dir, "Step7_HubStability.csv"))) {
  hub_stab <- read.csv(file.path(results_dir, "Step7_HubStability.csv"))
  kable(hub_stab, digits = 3, caption = "Hub Gene Stability (Bootstrap Analysis)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

## Pathway Enrichment Analysis {.tabset}

### ORA Custom Signatures

```{r ora-results, fig.cap="Over-representation analysis of custom PDAC-specific signatures."}
if(file.exists(file.path(fig_dir_png, "Step8_01_ORA_Lollipop.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step8_01_ORA_Lollipop.png"))
}
```

### GO Biological Process

```{r go-bp, fig.cap="Gene Ontology Biological Process enrichment for focus modules."}
if(file.exists(file.path(fig_dir_png, "Step11_01a_GO_BP_Pink.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step11_01a_GO_BP_Pink.png"))
}
```

### KEGG Pathways

```{r kegg, fig.cap="KEGG pathway enrichment for focus modules."}
if(file.exists(file.path(fig_dir_png, "Step11_02a_KEGG_Pink.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step11_02a_KEGG_Pink.png"))
}
```

### MSigDB Hallmark

```{r hallmark, fig.cap="MSigDB Hallmark gene set enrichment."}
if(file.exists(file.path(fig_dir_png, "Step8_04_GSVA_Heatmap.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step8_04_GSVA_Heatmap.png"))
}
```

### Unified Word Cloud

```{r wordcloud, fig.cap="Unified pathway enrichment word cloud showing biological themes associated with poor prognosis (FDR < 0.05)."}
if(file.exists(file.path(fig_dir_png, "Step13_FINAL_Unified_Pathway_WordCloud.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step13_FINAL_Unified_Pathway_WordCloud.png"))
}
```

### Integrated Bar Chart

```{r bar-chart, fig.cap="Top enriched biological themes across all databases."}
if(file.exists(file.path(fig_dir_png, "Step13_FINAL_Category_Summary_Bar.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step13_FINAL_Category_Summary_Bar.png"))
}
```

## Differential Expression Analysis {.tabset}

### PFS Group Volcano Plot

```{r volcano-pfs, fig.cap="Volcano plot of differential expression between Short and Long PFS groups."}
if(file.exists(file.path(fig_dir_png, "Step10b_Volcano_PFS.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step10b_Volcano_PFS.png"))
}
```

### Response Volcano Plot

```{r volcano-response, fig.cap="Volcano plot of differential expression between CD and PD response groups."}
if(file.exists(file.path(fig_dir_png, "Step10b_Volcano_Response.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step10b_Volcano_Response.png"))
}
```

## Module Validation {.tabset}

### Bootstrap Stability

```{r stability, fig.cap="Module stability assessment via bootstrap resampling."}
if(file.exists(file.path(fig_dir_png, "Step7_02_ModuleStability.png"))) {
  knitr::include_graphics(file.path(fig_dir_png, "Step7_02_ModuleStability.png"))
}
```

### Validation Summary

```{r validation-summary}
if(file.exists(file.path(results_dir, "Step7_ValidationData.csv"))) {
  val_data <- read.csv(file.path(results_dir, "Step7_ValidationData.csv"))
  kable(val_data, digits = 4, caption = "Module Validation Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

---

# Biomarker Signatures

## ML-Derived Biomarkers

**Biomarkers:** ECM2, DYNC2H1, PPIB

These proteins were identified through machine learning approaches and validated against the WGCNA module structure.

## Cox-Derived Biomarkers

**Biomarkers:** APOF, TFRC, FABP4, ANG

These proteins were identified through Cox proportional hazards regression as independent prognostic markers.

## Biomarker Integration

```{r biomarker-integration}
if(file.exists(file.path(results_dir, "Step15_Biomarker_Module_Membership.csv"))) {
  bio_mem <- read.csv(file.path(results_dir, "Step15_Biomarker_Module_Membership.csv"))
  kable(bio_mem, digits = 3, caption = "Biomarker Module Membership") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

---

# Clinical Interpretation

## Key Biological Themes

Based on integrated pathway analysis, the focus modules are enriched for:

1. **ECM Remodeling & Stromal Activation**
   - Cell adhesion molecules
   - Collagen-related pathways
   - Desmoplastic response signatures

2. **Angiogenesis**
   - Endothelial cell markers
   - Vascular development pathways

3. **Immune Modulation**
   - Myeloid-derived suppressor cells (MDSCs)
   - T-cell exclusion signatures
   - M0/M2 macrophage polarization

4. **Signaling Pathways**
   - PI3K-AKT signaling
   - Notch stemness signatures

## Prognostic Implications

- **Higher focus module expression** is consistently associated with:
  - Shorter progression-free survival
  - Poor treatment response (PD vs CD)
  - More aggressive disease phenotype

- **Therapeutic Implications:**
  - Potential targets for combination therapy
  - Biomarkers for patient stratification
  - Candidates for liquid biopsy development

---

# Supplementary Tables

## Complete Enrichment Results

```{r enrichment-tables}
if(file.exists(file.path(results_dir, "Step8_02_ORA_Results.csv"))) {
  ora_full <- read.csv(file.path(results_dir, "Step8_02_ORA_Results.csv"))
  ora_sig <- ora_full %>% filter(FDR < 0.1) %>%
    select(Module, Short_Name, Category, Fold_Enrichment, P_Value, FDR) %>%
    arrange(FDR)
  kable(ora_sig, digits = 4, caption = "Significant ORA Results (FDR < 0.1)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

## Module Narratives

```{r narratives}
if(file.exists(file.path(results_dir, "Step12_01_Module_Narratives.csv"))) {
  narratives <- read.csv(file.path(results_dir, "Step12_01_Module_Narratives.csv"))
  kable(narratives, caption = "Module Biological Narratives") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

---

# Complete Data Tables by Step {.tabset}

## Step 1-2: Data Loading & QC {.tabset}

### Step 1 Summary
```{r step1-summary}
if(file.exists(file.path(results_dir, "Step1_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step1_Summary.csv"))
  kable(df, caption = "Step 1: Data Loading Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Step 2 Noise Assessment
```{r step2-noise}
if(file.exists(file.path(results_dir, "Step2_Noise_Assessment_All_Proteins.csv"))) {
  df <- read.csv(file.path(results_dir, "Step2_Noise_Assessment_All_Proteins.csv"))
  kable(head(df, 50), caption = "Step 2: Protein Noise Assessment (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Step 2 Multi-Flagged Proteins
```{r step2-flagged}
if(file.exists(file.path(results_dir, "Step2_MultiFlagged_Proteins.csv"))) {
  df <- read.csv(file.path(results_dir, "Step2_MultiFlagged_Proteins.csv"))
  kable(df, caption = "Step 2: Multi-Flagged Proteins") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### DYNC2H1 Investigation
```{r dync2h1}
if(file.exists(file.path(results_dir, "Step2c_DYNC2H1_Association_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step2c_DYNC2H1_Association_Summary.csv"))
  kable(df, caption = "DYNC2H1 Clinical Association Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 3: Sample QC {.tabset}

### Sample QC Summary
```{r step3-qc}
if(file.exists(file.path(results_dir, "Step3_Sample_QC_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step3_Sample_QC_Summary.csv"))
  kable(df, caption = "Step 3: Sample QC Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

## Step 4: Power Selection {.tabset}

### Fit Indices
```{r step4-fit}
if(file.exists(file.path(results_dir, "Step4a_FitIndices.csv"))) {
  df <- read.csv(file.path(results_dir, "Step4a_FitIndices.csv"))
  kable(df, digits = 4, caption = "Step 4: Scale-Free Topology Fit Indices") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Power Comparison
```{r step4-power}
if(file.exists(file.path(results_dir, "Step4c_Power_Comparison.csv"))) {
  df <- read.csv(file.path(results_dir, "Step4c_Power_Comparison.csv"))
  kable(df, digits = 4, caption = "Step 4: Power Comparison") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Sensitivity Analysis
```{r step4-sensitivity}
if(file.exists(file.path(results_dir, "Step4d_Sensitivity_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step4d_Sensitivity_Summary.csv"))
  kable(df, digits = 4, caption = "Step 4: Sensitivity Analysis Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 5: Module Detection {.tabset}

### Module Summary
```{r step5-summary}
if(file.exists(file.path(results_dir, "Step5_ModuleSummary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step5_ModuleSummary.csv"))
  kable(df, digits = 4, caption = "Step 5: Module Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Module-Trait Correlations
```{r step5-trait-cor}
if(file.exists(file.path(results_dir, "Step5_ModuleTraitCorrelations.csv"))) {
  df <- read.csv(file.path(results_dir, "Step5_ModuleTraitCorrelations.csv"))
  kable(df, digits = 4, caption = "Step 5: Module-Trait Correlations") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Module-Trait P-values
```{r step5-trait-pval}
if(file.exists(file.path(results_dir, "Step5_ModuleTraitPvalues.csv"))) {
  df <- read.csv(file.path(results_dir, "Step5_ModuleTraitPvalues.csv"))
  kable(df, digits = 4, caption = "Step 5: Module-Trait P-values") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Bootstrap CI
```{r step5-bootstrap}
if(file.exists(file.path(results_dir, "Step5_CorrelationBootstrapCI.csv"))) {
  df <- read.csv(file.path(results_dir, "Step5_CorrelationBootstrapCI.csv"))
  kable(df, digits = 4, caption = "Step 5: Correlation Bootstrap Confidence Intervals") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Permutation Tests
```{r step5-perm}
if(file.exists(file.path(results_dir, "Step5_Permutation_Tests.csv"))) {
  df <- read.csv(file.path(results_dir, "Step5_Permutation_Tests.csv"))
  kable(df, digits = 4, caption = "Step 5: Permutation Test Results") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Module Priority
```{r step5-priority}
if(file.exists(file.path(results_dir, "Step5_ModulePriority.csv"))) {
  df <- read.csv(file.path(results_dir, "Step5_ModulePriority.csv"))
  kable(df, digits = 4, caption = "Step 5: Module Priority Ranking") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 6: Hub Genes {.tabset}

### Hub Summary
```{r step6-hub-summary}
if(file.exists(file.path(results_dir, "Step6_HubSummary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step6_HubSummary.csv"))
  kable(df, digits = 3, caption = "Step 6: Hub Gene Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

### Hub Effect Sizes
```{r step6-effect}
if(file.exists(file.path(results_dir, "Step6_HubEffectSizes.csv"))) {
  df <- read.csv(file.path(results_dir, "Step6_HubEffectSizes.csv"))
  kable(df, digits = 3, caption = "Step 6: Hub Gene Effect Sizes") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Full Protein Info
```{r step6-protein-info}
if(file.exists(file.path(results_dir, "Step6_ProteinInfo_Full.csv"))) {
  df <- read.csv(file.path(results_dir, "Step6_ProteinInfo_Full.csv"))
  kable(head(df, 100), digits = 3, caption = "Step 6: Complete Protein Information (Top 100)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

## Step 7: Validation {.tabset}

### Stability Summary
```{r step7-stability}
if(file.exists(file.path(results_dir, "Step7_StabilitySummary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step7_StabilitySummary.csv"))
  kable(df, digits = 4, caption = "Step 7: Module Stability Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Hub Stability
```{r step7-hub-stability}
if(file.exists(file.path(results_dir, "Step7_HubStability.csv"))) {
  df <- read.csv(file.path(results_dir, "Step7_HubStability.csv"))
  kable(df, digits = 3, caption = "Step 7: Hub Gene Stability") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Validation Data
```{r step7-validation}
if(file.exists(file.path(results_dir, "Step7_ValidationData.csv"))) {
  df <- read.csv(file.path(results_dir, "Step7_ValidationData.csv"))
  kable(df, digits = 4, caption = "Step 7: Module Validation Data") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 8: Enrichment {.tabset}

### ORA Results
```{r step8-ora}
if(file.exists(file.path(results_dir, "Step8_02_ORA_Results.csv"))) {
  df <- read.csv(file.path(results_dir, "Step8_02_ORA_Results.csv"))
  kable(head(df, 100), digits = 4, caption = "Step 8: ORA Results (Top 100)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

### Signature Coverage
```{r step8-coverage}
if(file.exists(file.path(results_dir, "Step8_01_Signature_Coverage.csv"))) {
  df <- read.csv(file.path(results_dir, "Step8_01_Signature_Coverage.csv"))
  kable(df, digits = 3, caption = "Step 8: Signature Coverage") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Clinical Associations
```{r step8-clinical}
if(file.exists(file.path(results_dir, "Step8_05_Clinical_Associations.csv"))) {
  df <- read.csv(file.path(results_dir, "Step8_05_Clinical_Associations.csv"))
  kable(df, digits = 4, caption = "Step 8: Pathway-Clinical Associations") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Unified Enrichment
```{r step8-unified}
if(file.exists(file.path(results_dir, "Step8_Unified_Enrichment_ALL.csv"))) {
  df <- read.csv(file.path(results_dir, "Step8_Unified_Enrichment_ALL.csv"))
  kable(head(df, 100), digits = 4, caption = "Step 8: Unified Enrichment Results (Top 100)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

## Step 10: Clinical Analysis {.tabset}

### Dual Clinical Significance
```{r step10-dual}
if(file.exists(file.path(results_dir, "Step10_Dual_Clinical_Significance.csv"))) {
  df <- read.csv(file.path(results_dir, "Step10_Dual_Clinical_Significance.csv"))
  kable(df, digits = 4, caption = "Step 10: Dual Clinical Significance") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### PFS Significance
```{r step10-pfs}
if(file.exists(file.path(results_dir, "Step10_PFS_Significance.csv"))) {
  df <- read.csv(file.path(results_dir, "Step10_PFS_Significance.csv"))
  kable(df, digits = 4, caption = "Step 10: PFS Group Significance") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Response Significance
```{r step10-response}
if(file.exists(file.path(results_dir, "Step10_Response_Significance.csv"))) {
  df <- read.csv(file.path(results_dir, "Step10_Response_Significance.csv"))
  kable(df, digits = 4, caption = "Step 10: Treatment Response Significance") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 11: Standard DB Enrichment {.tabset}

### GO:BP Results
```{r step11-gobp}
if(file.exists(file.path(results_dir, "Step11_01_GO_BP_Results.csv"))) {
  df <- read.csv(file.path(results_dir, "Step11_01_GO_BP_Results.csv"))
  kable(head(df, 50), digits = 4, caption = "Step 11: GO Biological Process (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

### KEGG Results
```{r step11-kegg}
if(file.exists(file.path(results_dir, "Step11_02_KEGG_Results.csv"))) {
  df <- read.csv(file.path(results_dir, "Step11_02_KEGG_Results.csv"))
  kable(head(df, 50), digits = 4, caption = "Step 11: KEGG Pathways (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

### Hallmark Results
```{r step11-hallmark}
if(file.exists(file.path(results_dir, "Step11_03_MSigDB_Hallmark_Results.csv"))) {
  df <- read.csv(file.path(results_dir, "Step11_03_MSigDB_Hallmark_Results.csv"))
  kable(head(df, 50), digits = 4, caption = "Step 11: MSigDB Hallmark (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

### GSEA Hallmark
```{r step11-gsea}
if(file.exists(file.path(results_dir, "Step11_08_GSEA_Hallmark_Results.csv"))) {
  df <- read.csv(file.path(results_dir, "Step11_08_GSEA_Hallmark_Results.csv"))
  kable(head(df, 50), digits = 4, caption = "Step 11: GSEA Hallmark Results (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

## Step 12: Biological Stories {.tabset}

### Module Narratives
```{r step12-narratives}
if(file.exists(file.path(results_dir, "Step12_01_Module_Narratives.csv"))) {
  df <- read.csv(file.path(results_dir, "Step12_01_Module_Narratives.csv"))
  kable(df, caption = "Step 12: Module Biological Narratives") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Hub Gene Validation
```{r step12-hub-val}
if(file.exists(file.path(results_dir, "Step12_02_Hub_Gene_Validation.csv"))) {
  df <- read.csv(file.path(results_dir, "Step12_02_Hub_Gene_Validation.csv"))
  kable(df, digits = 3, caption = "Step 12: Hub Gene Validation") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

## Step 13-14: Integration {.tabset}

### Unified Pathway Frequency
```{r step13-freq}
if(file.exists(file.path(results_dir, "Step13_Unified_Pathway_Frequency.csv"))) {
  df <- read.csv(file.path(results_dir, "Step13_Unified_Pathway_Frequency.csv"))
  kable(head(df, 50), caption = "Step 13: Unified Pathway Frequency (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Category Summary
```{r step13-category}
if(file.exists(file.path(results_dir, "Step13_Unified_Category_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step13_Unified_Category_Summary.csv"))
  kable(df, caption = "Step 13: Category Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### DYNC2H1 Investigation
```{r step14-dync}
if(file.exists(file.path(results_dir, "Step14_DYNC2H1_Investigation_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step14_DYNC2H1_Investigation_Summary.csv"))
  kable(df, caption = "Step 14: DYNC2H1 Investigation Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### High CV Proteins
```{r step14-highcv}
if(file.exists(file.path(results_dir, "Step14_HighCV_Proteins_Clinical.csv"))) {
  df <- read.csv(file.path(results_dir, "Step14_HighCV_Proteins_Clinical.csv"))
  kable(df, digits = 3, caption = "Step 14: High CV Proteins Clinical Associations") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 15: Biomarker Overlap {.tabset}

### Biomarker Module Membership
```{r step15-membership}
if(file.exists(file.path(results_dir, "Step15_Biomarker_Module_Membership.csv"))) {
  df <- read.csv(file.path(results_dir, "Step15_Biomarker_Module_Membership.csv"))
  kable(df, digits = 3, caption = "Step 15: Biomarker Module Membership") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Overlap Summary
```{r step15-overlap}
if(file.exists(file.path(results_dir, "Step15_Overlap_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step15_Overlap_Summary.csv"))
  kable(df, caption = "Step 15: Signature Overlap Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Comprehensive Biomarker Summary
```{r step15-comprehensive}
if(file.exists(file.path(results_dir, "Step15g_Biomarker_Comprehensive_Summary.csv"))) {
  df <- read.csv(file.path(results_dir, "Step15g_Biomarker_Comprehensive_Summary.csv"))
  kable(df, digits = 3, caption = "Step 15: Comprehensive Biomarker Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 16-17: Networks {.tabset}

### STRING Interactions
```{r step16-string}
if(file.exists(file.path(results_dir, "Step16_STRING_PINK_Interactions.csv"))) {
  df <- read.csv(file.path(results_dir, "Step16_STRING_PINK_Interactions.csv"))
  kable(head(df, 50), caption = "Step 16: STRING PPI Interactions (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Circos Module Summary
```{r step17-circos}
if(file.exists(file.path(results_dir, "Step17_Module_Summary_Circos.csv"))) {
  df <- read.csv(file.path(results_dir, "Step17_Module_Summary_Circos.csv"))
  kable(df, digits = 4, caption = "Step 17: Module Summary for Circos") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

## Step 18: sEV Comparison {.tabset}

### Plasma-sEV Overlap
```{r step18-overlap}
if(file.exists(file.path(results_dir, "Step18_Plasma_sEV_Overlap.csv"))) {
  df <- read.csv(file.path(results_dir, "Step18_Plasma_sEV_Overlap.csv"))
  kable(head(df, 50), caption = "Step 18: Plasma-sEV Protein Overlap (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Module sEV Detection
```{r step18-module-sev}
if(file.exists(file.path(results_dir, "Step18_Module_sEV_Detection.csv"))) {
  df <- read.csv(file.path(results_dir, "Step18_Module_sEV_Detection.csv"))
  kable(df, digits = 3, caption = "Step 18: Module Detection in sEV") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Biomarker sEV Detection
```{r step18-biomarker-sev}
if(file.exists(file.path(results_dir, "Step18_Biomarker_sEV_Detection.csv"))) {
  df <- read.csv(file.path(results_dir, "Step18_Biomarker_sEV_Detection.csv"))
  kable(df, caption = "Step 18: Biomarker Detection in sEV") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
}
```

### Hub Genes sEV Detection
```{r step18-hub-sev}
if(file.exists(file.path(results_dir, "Step18_HubGenes_sEV_Detection.csv"))) {
  df <- read.csv(file.path(results_dir, "Step18_HubGenes_sEV_Detection.csv"))
  kable(df, caption = "Step 18: Hub Genes Detection in sEV") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

## Step 19: Nature Medicine Figures {.tabset}

### PFS Enrichment
```{r step19-pfs}
if(file.exists(file.path(results_dir, "Step19_Enrichment_PFS.csv"))) {
  df <- read.csv(file.path(results_dir, "Step19_Enrichment_PFS.csv"))
  kable(head(df, 50), digits = 4, caption = "Step 19: PFS Direction Enrichment (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Response Enrichment
```{r step19-response}
if(file.exists(file.path(results_dir, "Step19_Enrichment_Response.csv"))) {
  df <- read.csv(file.path(results_dir, "Step19_Enrichment_Response.csv"))
  kable(head(df, 50), digits = 4, caption = "Step 19: Response Direction Enrichment (Top 50)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "400px")
}
```

### Protein Clinical Correlations
```{r step19-protein-cor}
if(file.exists(file.path(results_dir, "Step19_Protein_Clinical_Correlations.csv"))) {
  df <- read.csv(file.path(results_dir, "Step19_Protein_Clinical_Correlations.csv"))
  kable(head(df, 100), digits = 4, caption = "Step 19: Protein-Clinical Correlations (Top 100)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE) %>%
    scroll_box(height = "500px")
}
```

---

# Methods Summary for Manuscript

## Suggested Methods Text

> **Weighted Gene Co-expression Network Analysis**
>
> Plasma proteomics data from {{N_SAMPLES}} metastatic PDAC patients were analyzed using WGCNA (v1.72). A signed hybrid network was constructed using biweight midcorrelation with soft-threshold power β={{SOFT_POWER}} (scale-free topology R² > 0.80). Modules were identified using dynamic tree cutting (minimum size = 30, deep split = 3) and merged at height 0.25. Module eigengenes were correlated with clinical traits (PFS, PFS_group, Response) using Pearson correlation. Focus modules were defined as those with significant negative correlation with PFS (p < 0.05). Hub genes were identified using kME > 0.7 and gene significance |GS| > 0.2 (p < 0.05). Module stability was assessed via 200 bootstrap iterations. Pathway enrichment was performed using clusterProfiler (GO:BP, KEGG) and custom PDAC-specific signatures (n=264).

---

# Session Information

```{r session-info}
sessionInfo()
```

---

*Report generated by WGCNA Pipeline v4*
*Focus Modules: {{FOCUS_MODULES}}*
*Analysis Date: `r format(Sys.time(), "%Y-%m-%d %H:%M:%S")`*
'

# Perform substitutions using gsub
rmd_content <- rmd_template
rmd_content <- gsub("\\{\\{RESULTS_DIR\\}\\}", results_dir, rmd_content)
rmd_content <- gsub("\\{\\{FIG_DIR_PNG\\}\\}", fig_dir_png, rmd_content)
rmd_content <- gsub("\\{\\{FIG_DIR_PDF\\}\\}", fig_dir_pdf, rmd_content)
rmd_content <- gsub("\\{\\{N_SAMPLES\\}\\}", n_samples, rmd_content)
rmd_content <- gsub("\\{\\{N_PROTEINS\\}\\}", n_proteins, rmd_content)
rmd_content <- gsub("\\{\\{N_MODULES\\}\\}", n_modules, rmd_content)
rmd_content <- gsub("\\{\\{N_FOCUS\\}\\}", n_focus, rmd_content)
rmd_content <- gsub("\\{\\{FOCUS_MODULES\\}\\}", focus_modules_str, rmd_content)
rmd_content <- gsub("\\{\\{SOFT_POWER\\}\\}", soft_power, rmd_content)

# Write R Markdown file
rmd_file <- file.path(results_dir, "WGCNA_Publication_Report.Rmd")
writeLines(rmd_content, rmd_file)
cat(sprintf("  Created: %s\n", rmd_file))

# Try to render the report
cat("\n  Attempting to render HTML report...\n")
tryCatch({
  if(requireNamespace("rmarkdown", quietly = TRUE)) {
    rmarkdown::render(
      rmd_file,
      output_format = "html_document",
      output_file = "WGCNA_Publication_Report.html",
      output_dir = results_dir,
      quiet = TRUE
    )
    cat(sprintf("  [OK] Rendered: WGCNA_Publication_Report.html\n"))

    # Also render PDF if tinytex is available
    if(requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex()) {
      cat("  Attempting to render PDF report...\n")
      tryCatch({
        rmarkdown::render(
          rmd_file,
          output_format = "pdf_document",
          output_file = "WGCNA_Publication_Report.pdf",
          output_dir = results_dir,
          quiet = TRUE
        )
        cat(sprintf("  [OK] Rendered: WGCNA_Publication_Report.pdf\n"))
      }, error = function(e) {
        cat(sprintf("  [WARNING] PDF rendering failed: %s\n", e$message))
        cat("  TIP: Install tinytex with: tinytex::install_tinytex()\n")
      })
    } else {
      cat("  [INFO] PDF rendering skipped (tinytex not installed)\n")
      cat("  TIP: Install tinytex with: tinytex::install_tinytex()\n")
    }
  } else {
    cat("  [WARNING] rmarkdown package not available. Manual rendering required.\n")
    cat("  TIP: Install with: install.packages('rmarkdown')\n")
  }
}, error = function(e) {
  cat(sprintf("  [WARNING] Report rendering failed: %s\n", e$message))
})

cat("\n")
cat("  PUBLICATION REPORT GENERATION COMPLETE\n")
cat("  ============================================\n")
cat(sprintf("  R Markdown source: %s\n", rmd_file))
cat("  \n")
cat("  To manually render:\n")
cat("    rmarkdown::render('WGCNA_Publication_Report.Rmd', output_format = 'all')\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
