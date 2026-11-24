
####################################################################################################
# EFFECT SIZES (VOLCANO PLOTS) AND DIFFERENTIAL EXPRESSION ANALYSIS
#				A. Pairwise comparisons & Volcano plots
#   1. BAA-TNBC vs ASN-TNBC
#   2. BAA-TNBC vs HLW-TNBC
#   3. BAA-TNBC vs NHW-TNBC
#   4. BAA-TNBC vs non-BAA TNBC
#
# 				B. DEG identification
#.  1. BAA-TNBC vs non-BAA TNBC only
####################################################################################################


#---------------------------------------------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
#---------------------------------------------------------------------------------------------------

library(readxl)
library(limma)
library(ggplot2)
library(ggrepel)


#---------------------------------------------------------------------------------------------------
# 2. IMPORT AND PREPARE DATASETS
#---------------------------------------------------------------------------------------------------

# Set working directory
setwd("/Users/drpierred.sarr/Desktop/TCGA-BRCA")

# Import dataset (genes in rows, samples in columns)
BAA <- read_excel("BAA.xlsx", sheet = "BAA") |> as.data.frame()
nonBAA <- read_excel("nonBAA.xlsx", sheet = "nonBAA") |> as.data.frame()
ASN <- read_excel("ASN.xlsx", sheet = "ASN") |> as.data.frame()
HLW <- read_excel("HLW.xlsx", sheet = "HLW") |> as.data.frame()
NHW <- read_excel("NHW.xlsx", sheet = "NHW") |> as.data.frame()

# Extract gene identifiers and sample names, then transpose expression matrices
prepare_data <- function(df, group_name) {

  genes <- make.unique(as.character(df[[1]]))
  samples <- colnames(df)[-1]

  expression <- data.frame(t(as.matrix(df[, -1, drop = FALSE])))
  colnames(expression) <- genes
  rownames(expression) <- samples

  expression$group <- group_name

  return(expression)
}

# Prepare all datasets
BAA <- prepare_data(BAA, "BAA")
nonBAA <- prepare_data(nonBAA, "nonBAA")
ASN <- prepare_data(ASN, "ASN")
HLW <- prepare_data(HLW, "HLW")
NHW <- prepare_data(NHW, "NHW")


#---------------------------------------------------------------------------------------------------
# 3. COMBINE DATASETS AND PREPARE EXPRESSION MATRIX
#---------------------------------------------------------------------------------------------------

# Retain genes shared across all five datasets
common_genes <- Reduce(
  intersect,
  list(
    setdiff(colnames(BAA), "group"),
    setdiff(colnames(nonBAA), "group"),
    setdiff(colnames(ASN), "group"),
    setdiff(colnames(HLW), "group"),
    setdiff(colnames(NHW), "group")
  )
)

# Combine all samples
data <- rbind(
  BAA[, c(common_genes, "group")],
  nonBAA[, c(common_genes, "group")],
  ASN[, c(common_genes, "group")],
  HLW[, c(common_genes, "group")],
  NHW[, c(common_genes, "group")]
)

# Convert expression data to numeric matrix (genes in rows, samples in columns)
expression_matrix <- t(as.matrix(data[, common_genes]))

# Define group factor
group <- factor(
  data$group,
  levels = c("BAA", "nonBAA", "ASN", "HLW", "NHW")
)


#---------------------------------------------------------------------------------------------------
# 4. DIFFERENTIAL EXPRESSION ANALYSIS USING LIMMA
#---------------------------------------------------------------------------------------------------

# Create design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Define pairwise contrasts (BAA is the reference group in all comparisons)
contrasts <- makeContrasts(
  BAA_vs_ASN = BAA - ASN,
  BAA_vs_HLW = BAA - HLW,
  BAA_vs_NHW = BAA - NHW,
  BAA_vs_nonBAA = BAA - nonBAA,
  levels = design
)

# Fit linear model and apply empirical Bayes moderation
fit <- lmFit(expression_matrix, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)


#---------------------------------------------------------------------------------------------------
# 5. EXTRACT RESULTS FOR EACH PAIRWISE COMPARISON
#---------------------------------------------------------------------------------------------------

results_BAA_ASN <- topTable(
  fit, coef = "BAA_vs_ASN",
  number = Inf, adjust.method = "BH", sort.by = "P"
)

results_BAA_HLW <- topTable(
  fit, coef = "BAA_vs_HLW",
  number = Inf, adjust.method = "BH", sort.by = "P"
)

results_BAA_NHW <- topTable(
  fit, coef = "BAA_vs_NHW",
  number = Inf, adjust.method = "BH", sort.by = "P"
)

results_BAA_nonBAA <- topTable(
  fit, coef = "BAA_vs_nonBAA",
  number = Inf, adjust.method = "BH", sort.by = "P"
)


#---------------------------------------------------------------------------------------------------
# 6. GENERATE VOLCANO PLOTS
#---------------------------------------------------------------------------------------------------

# Create output directory
dir.create("results", showWarnings = FALSE)

# Volcano plot function
plot_volcano <- function(results, comparison, y_metric = "FDR") {

  results$Gene <- rownames(results)

  # Define y-axis using either adjusted or unadjusted p-values
  if (y_metric == "FDR") {
    results$y_value <- -log10(pmax(results$adj.P.Val, .Machine$double.xmin))
    y_label <- expression(-log[10](FDR))
  } else {
    results$y_value <- -log10(pmax(results$P.Value, .Machine$double.xmin))
    y_label <- expression(-log[10](p-value))
  }

  # Classify genes according to statistical significance and effect size
  results$Regulation <- "Not significant"

  results$Regulation[
    results$y_value > -log10(0.05) & results$logFC > 1
  ] <- "Upregulated in BAA"

  results$Regulation[
    results$y_value > -log10(0.05) & results$logFC < -1
  ] <- "Downregulated in BAA"

  results$Regulation <- factor(
    results$Regulation,
    levels = c(
      "Upregulated in BAA",
      "Downregulated in BAA",
      "Not significant"
    )
  )

  # Generate volcano plot
  volcano <- ggplot(
    results,
    aes(x = logFC, y = y_value, color = Regulation)
  ) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_color_manual(
      values = c(
        "Upregulated in BAA" = "#D73027",
        "Downregulated in BAA" = "#4575B4",
        "Not significant" = "grey75"
      )
    ) +
    labs(
      title = comparison,
      x = expression(log[2]~Fold~Change),
      y = y_label,
      color = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom"
    )

  print(volcano)

  # Save volcano plot
  ggsave(
    filename = paste0(
      "results/Volcano_",
      gsub("[^A-Za-z0-9]+", "_", comparison),
      ".png"
    ),
    plot = volcano,
    width = 7,
    height = 6,
    dpi = 300
  )
}


# Generate the four volcano plots
plot_volcano(
  results_BAA_ASN,
  "BAA-TNBC vs ASN-TNBC",
  y_metric = "FDR"
)

plot_volcano(
  results_BAA_HLW,
  "BAA-TNBC vs HLW-TNBC",
  y_metric = "p-value"
)

plot_volcano(
  results_BAA_NHW,
  "BAA-TNBC vs NHW-TNBC",
  y_metric = "FDR"
)

plot_volcano(
  results_BAA_nonBAA,
  "BAA-TNBC vs non-BAA TNBC",
  y_metric = "FDR"
)


#---------------------------------------------------------------------------------------------------
#  IDENTIFICATION AND EXPORT OF DEGs: BAA-TNBC vs NON-BAA TNBC
#---------------------------------------------------------------------------------------------------

# Identify significantly differentially expressed genes
# Criteria: FDR < 0.05 and |log2FC| > 1

DEGs <- results_BAA_nonBAA[
  results_BAA_nonBAA$adj.P.Val < 0.05 &
    abs(results_BAA_nonBAA$logFC) > 1,
]

# Separate upregulated and downregulated genes in BAA-TNBC
DEGs_up <- DEGs[
  DEGs$logFC > 1,
  c("logFC", "adj.P.Val", "P.Value"),
  drop = FALSE
]

DEGs_down <- DEGs[
  DEGs$logFC < -1,
  c("logFC", "adj.P.Val", "P.Value"),
  drop = FALSE
]

# Add gene identifiers as a column
DEGs_up$Gene <- rownames(DEGs_up)
DEGs_down$Gene <- rownames(DEGs_down)

# Arrange output columns
DEGs_up <- DEGs_up[, c("Gene", "logFC", "adj.P.Val", "P.Value")]
DEGs_down <- DEGs_down[, c("Gene", "logFC", "adj.P.Val", "P.Value")]

# Rename columns for clarity
colnames(DEGs_up) <- c("Gene", "Log2FC", "FDR", "P_value")
colnames(DEGs_down) <- c("Gene", "Log2FC", "FDR", "P_value")

# Export upregulated and downregulated gene lists as tab-delimited text files
write.table(
  DEGs_up,
  file = "results/DEGs_BAA_vs_nonBAA_upregulated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  DEGs_down,
  file = "results/DEGs_BAA_vs_nonBAA_downregulated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Display DEG counts
cat("Upregulated DEGs in BAA-TNBC:", nrow(DEGs_up), "\n")
cat("Downregulated DEGs in BAA-TNBC:", nrow(DEGs_down), "\n")