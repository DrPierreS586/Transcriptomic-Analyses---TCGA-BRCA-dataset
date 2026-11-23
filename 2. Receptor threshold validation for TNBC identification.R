################Threshold Validation (ER/PR/HER2 cutoffs for TCGA-BRCA PAM50 Basal/Normal like Patients) - Pierre SARR - March 18, 2025################

rm(list = ls())

# ============================================================
# PACKAGES
# ============================================================
install.packages(c("readxl", "ggplot2", "gridExtra", "dplyr"))

library(readxl)
library(ggplot2)
library(gridExtra)
library(dplyr)

# ============================================================
# UPLOAD DATA
# ============================================================

path <- "/Users/drpierred.sarr/Desktop/PAM50 Bas-Norm like for TNBC prediction.xlsx"

df_Basal  <- read_excel(path, sheet = "Basal")
df_Normal <- read_excel(path, sheet = "Normal")

cat("Dimensions Basal  :", dim(df_Basal),  "\n")
cat("Dimensions Normal :", dim(df_Normal), "\n")
cat("Premières colonnes :", colnames(df_Basal)[1:5], "\n")
cat("Premiers gènes     :", head(df_Basal$Gene), "\n")

# ============================================================
# EXTRAIRE ESR1, PGR, ERBB2
# ============================================================

# Gènes en lignes → extraire valeurs numériques par gène
get_expr <- function(df, gene) {
  row <- df %>% dplyr::filter(Gene == gene)
  if (nrow(row) == 0) {
    cat("ATTENTION : gène", gene, "introuvable\n")
    return(NULL)
  }
  as.numeric(row[, -which(colnames(df) == "Gene")])
}

ESR1_Basal   <- get_expr(df_Basal,  "ESR1")
PGR_Basal    <- get_expr(df_Basal,  "PGR")
ERBB2_Basal  <- get_expr(df_Basal,  "ERBB2")

ESR1_Normal  <- get_expr(df_Normal, "ESR1")
PGR_Normal   <- get_expr(df_Normal, "PGR")
ERBB2_Normal <- get_expr(df_Normal, "ERBB2")

# Combiner Basal + Normal
ESR1_all  <- c(ESR1_Basal,  ESR1_Normal)
PGR_all   <- c(PGR_Basal,   PGR_Normal)
ERBB2_all <- c(ERBB2_Basal, ERBB2_Normal)

# Vérifications
cat("\n--- Statistiques descriptives ---\n")
cat("ESR1  — min:", round(min(ESR1_all,  na.rm = TRUE), 2),
            "max:", round(max(ESR1_all,  na.rm = TRUE), 2),
           "mean:", round(mean(ESR1_all, na.rm = TRUE), 2), "\n")
cat("PGR   — min:", round(min(PGR_all,   na.rm = TRUE), 2),
            "max:", round(max(PGR_all,   na.rm = TRUE), 2),
           "mean:", round(mean(PGR_all,  na.rm = TRUE), 2), "\n")
cat("ERBB2 — min:", round(min(ERBB2_all, na.rm = TRUE), 2),
            "max:", round(max(ERBB2_all, na.rm = TRUE), 2),
           "mean:", round(mean(ERBB2_all,na.rm = TRUE), 2), "\n")

# ============================================================
# BREAKDOWN PAR CUTOFF
# ============================================================

cat("\n--- Breakdown par cutoff ---\n")
cat("ESR1  < 4  :", sum(ESR1_all  < 4,  na.rm = TRUE), "samples |",
    round(sum(ESR1_all  < 4,  na.rm = TRUE) / length(ESR1_all)  * 100, 1), "%\n")
cat("ESR1  >= 4 :", sum(ESR1_all  >= 4, na.rm = TRUE), "samples |",
    round(sum(ESR1_all  >= 4, na.rm = TRUE) / length(ESR1_all)  * 100, 1), "%\n")
cat("PGR   < 3  :", sum(PGR_all   < 3,  na.rm = TRUE), "samples |",
    round(sum(PGR_all   < 3,  na.rm = TRUE) / length(PGR_all)   * 100, 1), "%\n")
cat("PGR   >= 3 :", sum(PGR_all   >= 3, na.rm = TRUE), "samples |",
    round(sum(PGR_all   >= 3, na.rm = TRUE) / length(PGR_all)   * 100, 1), "%\n")
cat("ERBB2 < 10 :", sum(ERBB2_all < 10, na.rm = TRUE), "samples |",
    round(sum(ERBB2_all < 10, na.rm = TRUE) / length(ERBB2_all) * 100, 1), "%\n")
cat("ERBB2 >= 10:", sum(ERBB2_all >= 10, na.rm = TRUE), "samples |",
    round(sum(ERBB2_all >= 10, na.rm = TRUE) / length(ERBB2_all) * 100, 1), "%\n")

# ============================================================
# HISTOGRAMMES — GGPLOT2 (version publication)
# ============================================================

df_plot <- data.frame(
  ESR1  = ESR1_all,
  PGR   = PGR_all,
  ERBB2 = ERBB2_all
)

genes   <- c("ESR1", "PGR", "ERBB2")
cutoffs <- c(4, 3, 10)
colors  <- c("#AED6F1", "#A9DFBF", "#FAD7A0")
labels  <- c("ER (ESR1)", "PR (PGR)", "HER2 (ERBB2)")

library(ggplot2)
library(gridExtra)

plots <- lapply(seq_along(genes), function(i) {
  gene   <- genes[i]
  cutoff <- cutoffs[i]
  color  <- colors[i]
  label  <- labels[i]

  # Proportion sous le cutoff
  n_below <- sum(df_plot[[gene]] < cutoff, na.rm = TRUE)
  n_total <- sum(!is.na(df_plot[[gene]]))
  pct     <- round(n_below / n_total * 100, 1)

  ggplot(df_plot, aes(x = .data[[gene]])) +
    geom_histogram(bins = 50, fill = color, color = "white") +
    geom_vline(xintercept = cutoff,
               color      = "red",
               linewidth  = 0.8,
               linetype   = "dashed") +
    annotate("text",
             x      = cutoff + diff(range(df_plot[[gene]], na.rm = TRUE)) * 0.05,
             y      = Inf,
             label  = paste0("Cutoff = ", cutoff,
                             "\n(", pct, "% negative)"),
             vjust  = 1.5,
             hjust  = 0,
             color  = "red",
             size   = 3.2,
             fontface = "bold") +
    labs(
      title = paste(label, "— Basal-like & Normal-like"),
      x     = "log2 Expression (TPM)",
      y     = "Number of Samples"
    ) +
    theme_classic() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text   = element_text(size = 10),
      axis.title  = element_text(size = 11)
    )
})

# Afficher
grid.arrange(grobs = plots, ncol = 3)

# ============================================================
# EXPORT PDF
# ============================================================

output_dir <- dirname(path)

pdf(paste0(output_dir, "/Threshold_Validation_ESR1_PGR_ERBB2.pdf"),
    width = 14, height = 5)
grid.arrange(grobs = plots, ncol = 3,
             top = "TNBC Empirical Thresholds Validation")
dev.off()
cat("PDF sauvegardé : Threshold_Validation_ESR1_PGR_ERBB2.pdf\n")

# ============================================================
# CORRECTIF — Labels complets + marges élargies
# ============================================================

plots <- lapply(seq_along(genes), function(i) {
  gene   <- genes[i]
  cutoff <- cutoffs[i]
  color  <- colors[i]
  label  <- labels[i]

  n_below <- sum(df_plot[[gene]] < cutoff, na.rm = TRUE)
  n_total <- sum(!is.na(df_plot[[gene]]))
  pct     <- round(n_below / n_total * 100, 1)

  ggplot(df_plot, aes(x = .data[[gene]])) +
    geom_histogram(bins = 50, fill = color, color = "white") +
    geom_vline(xintercept = cutoff,
               color      = "red",
               linewidth  = 0.8,
               linetype   = "dashed") +
    annotate("text",
             x        = cutoff,
             y        = Inf,
             label    = paste0("Cutoff = ", cutoff, "\n(", pct, "% negative)"),
             vjust    = 1.5,
             hjust    = -0.1,        # légèrement à droite du cutoff
             color    = "red",
             size     = 3.2,
             fontface = "bold") +
    labs(
      title = paste(label, " "),
      x     = "log2 Expression (TPM)",
      y     = "Number of Samples"
    ) +
    theme_classic() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text    = element_text(size = 10),
      axis.title   = element_text(size = 11),
      plot.margin  = margin(10, 30, 10, 10)  # marge droite élargie
    ) +
    coord_cartesian(clip = "off")  # ← clé : empêche la coupure des labels
})

# Afficher
grid.arrange(grobs = plots, ncol = 3)

# Export PDF — largeur augmentée
pdf(paste0(dirname(path), "/Threshold_Validation_ESR1_PGR_ERBB2.pdf"),
    width  = 16,   # augmenté de 14 à 16
    height = 5)
grid.arrange(grobs = plots, ncol = 3,
             top = "TNBC Empirical Thresholds Validation")
dev.off()
cat("PDF sauvegardé\n")