
################R Script TCGA-BRCA - Predicted TNBC - Pierre SARR, March 08, 2026################
							++++Addition Volcano plot nonBAA
													################

source("http://www.bioconductor.org/biocLite.R")
biocLite(c("limma"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)


rm(list=ls())
library(readxl)
library(openxlsx)
library(dplyr)
								######### I. DEGs & EFFECT SIZES ####################
#Dataset 1: ASN
df_ASN <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/ASN.xlsx", sheet = "ASN")
df_ASN_transposed <- data.frame(t(df_ASN))
rows<-colnames(df_ASN)[2:10]
cols <- df_ASN_transposed[1,]
df_ASN_transposed<- lapply(df_ASN_transposed[-1,], as.numeric)
df_ASN_transposed <- data.frame(df_ASN_transposed)
colnames(df_ASN_transposed) <- df_ASN$samples
rownames(df_ASN_transposed) <- rows
df_ASN_transposed$group <- "ASN"

#Dataset 2: BAA
df_BAA <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/BAA.xlsx", sheet = "BAA")
df_BAA_transposed <- data.frame(t(df_BAA))
rows<-colnames(df_BAA)[2:64]
cols <- df_BAA_transposed[1,]
df_BAA_transposed<- lapply(df_BAA_transposed[-1,], as.numeric)
df_BAA_transposed <- data.frame(df_BAA_transposed)
colnames(df_BAA_transposed) <- df_BAA$samples
rownames(df_BAA_transposed) <- rows
df_BAA_transposed$group <- "BAA"

#Dataset 3: HLW
df_HLW <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/HLW.xlsx", sheet = "HLW")
df_HLW_transposed <- data.frame(t(df_HLW))
rows<-colnames(df_HLW)[2:6]
cols <- df_HLW_transposed[1,]
df_HLW_transposed<- lapply(df_HLW_transposed[-1,], as.numeric)
df_HLW_transposed <- data.frame(df_HLW_transposed)
colnames(df_HLW_transposed) <- df_HLW $samples
rownames(df_HLW_transposed) <- rows
df_HLW_transposed$group <- "HLW"

#Dataset 4: NHW
df_NHW <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/NHW.xlsx", sheet = "NHW")
df_NHW_transposed <- data.frame(t(df_NHW))
rows<-colnames(df_NHW)[2:97]
cols <- df_NHW_transposed[1,]
df_NHW_transposed<- lapply(df_NHW_transposed[-1,], as.numeric)
df_NHW_transposed <- data.frame(df_NHW_transposed)
colnames(df_NHW_transposed) <- df_NHW $samples
rownames(df_NHW_transposed) <- rows
df_NHW_transposed$group <- "NHW"

#Dataset 5: nonBAA
df_nonBAA <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/nonBAA.xlsx", sheet = "nonBAA")
df_nonBAA_transposed <- data.frame(t(df_nonBAA))
rows<-colnames(df_nonBAA)[2:111]
cols <- df_nonBAA_transposed[1,]
df_nonBAA_transposed<- lapply(df_nonBAA_transposed[-1,], as.numeric)
df_nonBAA_transposed <- data.frame(df_nonBAA_transposed)
colnames(df_nonBAA_transposed) <- df_nonBAA $samples
rownames(df_nonBAA_transposed) <- rows
df_nonBAA_transposed$group <- "nonBAA"

#Pairwise dataset combination
#1. BAA vs ASN
df_BAA_ASN <- rbind(df_BAA_transposed,df_ASN_transposed)
df_BAA_ASN$condition <- factor(df_BAA_ASN$group, levels = c("BAA", "ASN"))

#2. BAA vs HLW
df_BAA_HLW <- rbind(df_BAA_transposed,df_HLW_transposed)
df_BAA_HLW$condition <- factor(df_BAA_HLW$group, levels = c("BAA", "HLW"))

#3. BAA vs NHW
df_BAA_NHW <- rbind(df_BAA_transposed,df_NHW_transposed)
df_BAA_NHW $condition <- factor(df_BAA_NHW $group, levels = c("BAA", "NHW"))

#4. BAA vs nonBAA
df_BAA_nonBAA <- rbind(df_BAA_transposed,df_nonBAA_transposed)
df_BAA_nonBAA $condition <- factor(df_BAA_nonBAA $group, levels = c("BAA", "nonBAA"))

#Identification of numeric column
#1. BAA vs ASN
numeric_cols <- sapply(df_BAA_ASN, is.numeric)
design <- model.matrix(~df_BAA_ASN$condition)

#2. BAA vs HLW
numeric_cols <- sapply(df_BAA_HLW, is.numeric)
design <- model.matrix(~df_BAA_HLW$condition)

#3. BAA vs NHW
numeric_cols <- sapply(df_BAA_NHW, is.numeric)
design <- model.matrix(~df_BAA_NHW$condition)

#4. BAA vs nonBAA
numeric_cols <- sapply(df_BAA_nonBAA, is.numeric)
design <- model.matrix(~df_BAA_nonBAA$condition)
					
#Voom (Données déjà log2-transformées)
#1. BAA vs ASN

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Data preparation ---
condition <- df_BAA_ASN$condition
num_df_BAA_ASN <- df_BAA_ASN[, !colnames(df_BAA_ASN) %in% c("condition", "group")]

dim(num_df_BAA_ASN)
length(condition)

# --- Design ---
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "ASN"

nrow(design)
ncol(t(num_df_BAA_ASN))

# --- Limma ---
variances   <- apply(t(num_df_BAA_ASN), 1, var)
mat_filtree <- t(num_df_BAA_ASN)[variances > 0, ]

contrast_matrix <- makeContrasts(BAA - ASN, levels = design)

fit  <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters -Test 1 Volcano--
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Parameters -Test 2 Volcano--
fdr_threshold  <- 0.25
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
) %>%
  mutate(categorie = case_when(
    FDR >= fdr_threshold & abs(logFC) <= fc_threshold ~ "Not significant",
    FDR <  fdr_threshold & abs(logFC) <= fc_threshold ~ "Significant (weak effect)",
    FDR <  fdr_threshold & abs(logFC) >  fc_threshold ~ "Significant (strong effect)",
    FDR >= fdr_threshold & abs(logFC) >  fc_threshold ~ "Not significant (strong effect)",
    TRUE ~ "Other"
  ))

# --- Volcano plot ---
top_genes <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(FDR) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = categorie), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Not significant"                 = "black",
      "Significant (weak effect)"       = "blue",
      "Significant (strong effect)"     = "red",
      "Not significant (strong effect)" = "green4"
    ),
    name = "Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    title = "BAA vs ASN"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(FDR)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_ASN.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_ASN.csv\n")
  
#2. BAA vs HLW
# Séparer la colonne condition des données numériques
condition <- df_BAA_HLW$condition

# Matrice numérique uniquement (sans la colonne condition/group)
num_df_BAA_HLW <- df_BAA_HLW[, !colnames(df_BAA_HLW) %in% c("condition", "group")]

# Vérifier
dim(num_df_BAA_HLW)
length(condition)

design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "HLW"

# Vérifier l'alignement
nrow(design)
ncol(t(num_df_BAA_HLW))

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Filtrer gènes à variance nulle
variances <- apply(t(num_df_BAA_HLW), 1, var)
mat_filtree <- t(num_df_BAA_HLW)[variances > 0, ]

# Contraste et modèle limma
contrast_matrix <- makeContrasts(BAA - HLW, levels = design)

fit <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters -Test 1 Volcano--
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Parameters -Test 2 Volcano--
fdr_threshold  <- 0.5
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
) %>%
  mutate(categorie = case_when(
    FDR >= fdr_threshold & abs(logFC) <= fc_threshold ~ "Not significant",
    FDR <  fdr_threshold & abs(logFC) <= fc_threshold ~ "Significant (weak effect)",
    FDR <  fdr_threshold & abs(logFC) >  fc_threshold ~ "Significant (strong effect)",
    FDR >= fdr_threshold & abs(logFC) >  fc_threshold ~ "Not significant (strong effect)",
    TRUE ~ "Other"
  ))

top_genes <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(FDR) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = categorie), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Not significant"                 = "black",
      "Significant (weak effect)"       = "blue",
      "Significant (strong effect)"     = "red",
      "Not significant (strong effect)" = "green4"
    ),
    name = "Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    title = "BAA vs HLW"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(FDR)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_HLW.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_HLW.csv\n")

#3. BAA vs NHW
condition <- df_BAA_NHW$condition
num_df_BAA_NHW <- df_BAA_NHW[, !colnames(df_BAA_NHW) %in% c("condition", "group")]

dim(num_df_BAA_NHW)
length(condition)

design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "NHW"

nrow(design)
ncol(t(num_df_BAA_NHW))

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Filtrer gènes à variance nulle
variances <- apply(t(num_df_BAA_NHW), 1, var)
mat_filtree <- t(num_df_BAA_NHW)[variances > 0, ]

# Contraste et modèle limma
contrast_matrix <- makeContrasts(BAA - NHW, levels = design)

fit <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters ---
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
) %>%
  mutate(categorie = case_when(
    FDR >= fdr_threshold & abs(logFC) <= fc_threshold ~ "Not significant",
    FDR <  fdr_threshold & abs(logFC) <= fc_threshold ~ "Significant (weak effect)",
    FDR <  fdr_threshold & abs(logFC) >  fc_threshold ~ "Significant (strong effect)",
    FDR >= fdr_threshold & abs(logFC) >  fc_threshold ~ "Not significant (strong effect)",
    TRUE ~ "Other"
  ))

# --- Volcano plot ---
top_genes <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(FDR) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = categorie), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Not significant"                 = "black",
      "Significant (weak effect)"       = "blue",
      "Significant (strong effect)"     = "red",
      "Not significant (strong effect)" = "green4"
    ),
    name = "Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    title = "BAA vs NHW"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(FDR)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_NHW.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_NHW.csv\n")

# 4. BAA vs nonBAA
condition <- df_BAA_nonBAA$condition
num_df_BAA_nonBAA <- df_BAA_nonBAA[, !colnames(df_BAA_nonBAA) %in% c("condition", "group")]

dim(num_df_BAA_nonBAA)
length(condition)

design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "nonBAA"

nrow(design)
ncol(t(num_df_BAA_nonBAA))

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Filtrer gènes à variance nulle
variances <- apply(t(num_df_BAA_nonBAA), 1, var)
mat_filtree <- t(num_df_BAA_nonBAA)[variances > 0, ]

# Contraste et modèle limma
contrast_matrix <- makeContrasts(BAA - nonBAA, levels = design)

fit <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters ---
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
) %>%
  mutate(categorie = case_when(
    FDR >= fdr_threshold & abs(logFC) <= fc_threshold ~ "Not significant",
    FDR <  fdr_threshold & abs(logFC) <= fc_threshold ~ "Significant (weak effect)",
    FDR <  fdr_threshold & abs(logFC) >  fc_threshold ~ "Significant (strong effect)",
    FDR >= fdr_threshold & abs(logFC) >  fc_threshold ~ "Not significant (strong effect)",
    TRUE ~ "Other"
  ))

# --- Volcano plot ---
top_genes <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(FDR) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = categorie), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Not significant"                 = "black",
      "Significant (weak effect)"       = "blue",
      "Significant (strong effect)"     = "red",
      "Not significant (strong effect)" = "green4"
    ),
    name = "Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    title = "BAA vs non-BAA"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(FDR < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(FDR)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_Others.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_Others.csv\n")


# ============================================================
# BAA vs nonBAA - Script readapté avec conditions & threshold Yann Ilboudo (Y pvalue et non FDR)
# ============================================================

source("http://www.bioconductor.org/biocLite.R")
biocLite(c("limma"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)


rm(list=ls())
library(readxl)
library(openxlsx)
library(dplyr)

#Dataset 1: ASN
df_ASN <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/ASN.xlsx", sheet = "ASN")
df_ASN_transposed <- data.frame(t(df_ASN))
rows<-colnames(df_ASN)[2:10]
cols <- df_ASN_transposed[1,]
df_ASN_transposed<- lapply(df_ASN_transposed[-1,], as.numeric)
df_ASN_transposed <- data.frame(df_ASN_transposed)
colnames(df_ASN_transposed) <- df_ASN$samples
rownames(df_ASN_transposed) <- rows
df_ASN_transposed$group <- "ASN"

#Dataset 2: BAA
df_BAA <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/BAA.xlsx", sheet = "BAA")
df_BAA_transposed <- data.frame(t(df_BAA))
rows<-colnames(df_BAA)[2:64]
cols <- df_BAA_transposed[1,]
df_BAA_transposed<- lapply(df_BAA_transposed[-1,], as.numeric)
df_BAA_transposed <- data.frame(df_BAA_transposed)
colnames(df_BAA_transposed) <- df_BAA$samples
rownames(df_BAA_transposed) <- rows
df_BAA_transposed$group <- "BAA"

#Dataset 3: HLW
df_HLW <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/HLW.xlsx", sheet = "HLW")
df_HLW_transposed <- data.frame(t(df_HLW))
rows<-colnames(df_HLW)[2:6]
cols <- df_HLW_transposed[1,]
df_HLW_transposed<- lapply(df_HLW_transposed[-1,], as.numeric)
df_HLW_transposed <- data.frame(df_HLW_transposed)
colnames(df_HLW_transposed) <- df_HLW $samples
rownames(df_HLW_transposed) <- rows
df_HLW_transposed$group <- "HLW"

#Dataset 4: NHW
df_NHW <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/NHW.xlsx", sheet = "NHW")
df_NHW_transposed <- data.frame(t(df_NHW))
rows<-colnames(df_NHW)[2:97]
cols <- df_NHW_transposed[1,]
df_NHW_transposed<- lapply(df_NHW_transposed[-1,], as.numeric)
df_NHW_transposed <- data.frame(df_NHW_transposed)
colnames(df_NHW_transposed) <- df_NHW $samples
rownames(df_NHW_transposed) <- rows
df_NHW_transposed$group <- "NHW"

#Dataset 5: nonBAA
df_nonBAA <- read_excel ("/Users/drpierred.sarr/Desktop/TCGA-BRCA/nonBAA.xlsx", sheet = "nonBAA")
df_nonBAA_transposed <- data.frame(t(df_nonBAA))
rows<-colnames(df_nonBAA)[2:111]
cols <- df_nonBAA_transposed[1,]
df_nonBAA_transposed<- lapply(df_nonBAA_transposed[-1,], as.numeric)
df_nonBAA_transposed <- data.frame(df_nonBAA_transposed)
colnames(df_nonBAA_transposed) <- df_nonBAA $samples
rownames(df_nonBAA_transposed) <- rows
df_nonBAA_transposed$group <- "nonBAA"

#Pairwise dataset combination
#1. BAA vs ASN
df_BAA_ASN <- rbind(df_BAA_transposed,df_ASN_transposed)
df_BAA_ASN$condition <- factor(df_BAA_ASN$group, levels = c("BAA", "ASN"))

#2. BAA vs HLW
df_BAA_HLW <- rbind(df_BAA_transposed,df_HLW_transposed)
df_BAA_HLW$condition <- factor(df_BAA_HLW$group, levels = c("BAA", "HLW"))

#3. BAA vs NHW
df_BAA_NHW <- rbind(df_BAA_transposed,df_NHW_transposed)
df_BAA_NHW $condition <- factor(df_BAA_NHW $group, levels = c("BAA", "NHW"))

#4. BAA vs nonBAA
df_BAA_nonBAA <- rbind(df_BAA_transposed,df_nonBAA_transposed)
df_BAA_nonBAA $condition <- factor(df_BAA_nonBAA $group, levels = c("BAA", "nonBAA"))

#Identification of numeric column
#1. BAA vs ASN
numeric_cols <- sapply(df_BAA_ASN, is.numeric)
design <- model.matrix(~df_BAA_ASN$condition)

#2. BAA vs HLW
numeric_cols <- sapply(df_BAA_HLW, is.numeric)
design <- model.matrix(~df_BAA_HLW$condition)

#3. BAA vs NHW
numeric_cols <- sapply(df_BAA_NHW, is.numeric)
design <- model.matrix(~df_BAA_NHW$condition)

#4. BAA vs nonBAA
numeric_cols <- sapply(df_BAA_nonBAA, is.numeric)
design <- model.matrix(~df_BAA_nonBAA$condition)


#########################################################################################
# ============================================================
# BAA vs ASN
# ============================================================

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Data preparation ---
condition <- df_BAA_ASN$condition
num_df_BAA_ASN <- df_BAA_ASN[, !colnames(df_BAA_ASN) %in% c("condition", "group")]

dim(num_df_BAA_ASN)
length(condition)

# --- Design ---
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "ASN"

# --- Limma ---
variances   <- apply(t(num_df_BAA_ASN), 1, var)
mat_filtree <- t(num_df_BAA_ASN)[variances > 0, ]

contrast_matrix <- makeContrasts(BAA - ASN, levels = design)

fit  <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters ---
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
)

# --- Top genes ---
top_genes <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(pvalue) %>%
  head(10)

# --- Volcano plot ---
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = case_when(
    -log10(pvalue) < 1.25 & logFC >= -1 & logFC <= 1      ~ "Black",
    -log10(pvalue) >= 1.25 & logFC >= -1 & logFC <= 1     ~ "Blue",
    -log10(pvalue) >= 1.25 & (logFC < -1 | logFC > 1)     ~ "Red",
    -log10(pvalue) < 1.25  & (logFC < -1 | logFC > 1)     ~ "Green",
    TRUE                                                    ~ "Fallback"
  )),
  alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Black"    = "black",
      "Blue"     = "blue",
      "Red"      = "red",
      "Green"    = "green4",
      "Fallback" = "gray60"
    ),
    labels = c(
      "Black"    = "Not significant",
      "Blue"     = "Significant & weak effect",
      "Red"      = "Significant & strong effect",
      "Green"    = "Not significant & strong effect"
    ),
    name = "Significance Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(p-value)",
    title = "BAA vs ASN"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(pvalue)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_ASN.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_ASN.csv\n")

# ============================================================
# BAA vs HLW
# ============================================================

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Data preparation ---
condition <- df_BAA_HLW$condition
num_df_BAA_HLW <- df_BAA_HLW[, !colnames(df_BAA_HLW) %in% c("condition", "group")]

dim(num_df_BAA_HLW)
length(condition)

# --- Design ---
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "HLW"

# --- Limma ---
variances   <- apply(t(num_df_BAA_HLW), 1, var)
mat_filtree <- t(num_df_BAA_HLW)[variances > 0, ]

contrast_matrix <- makeContrasts(BAA - HLW, levels = design)

fit  <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters ---
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
)

# --- Top genes ---
top_genes <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(pvalue) %>%
  head(10)

# --- Volcano plot ---
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = case_when(
    -log10(pvalue) < 1.25 & logFC >= -1 & logFC <= 1      ~ "Black",
    -log10(pvalue) >= 1.25 & logFC >= -1 & logFC <= 1     ~ "Blue",
    -log10(pvalue) >= 1.25 & (logFC < -1 | logFC > 1)     ~ "Red",
    -log10(pvalue) < 1.25  & (logFC < -1 | logFC > 1)     ~ "Green",
    TRUE                                                    ~ "Fallback"
  )),
  alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Black"    = "black",
      "Blue"     = "blue",
      "Red"      = "red",
      "Green"    = "green4",
      "Fallback" = "gray60"
    ),
    labels = c(
      "Black"    = "Not significant",
      "Blue"     = "Significant & weak effect",
      "Red"      = "Significant & strong effect",
      "Green"    = "Not significant & strong effect"
    ),
    name = "Significance Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(p-value)",
    title = "BAA vs HLW"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(pvalue)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_HLW.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_HLW.csv\n")

# ============================================================
# BAA vs NHW
# ============================================================

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Data preparation ---
condition <- df_BAA_NHW$condition
num_df_BAA_NHW <- df_BAA_NHW[, !colnames(df_BAA_NHW) %in% c("condition", "group")]

dim(num_df_BAA_NHW)
length(condition)

# --- Design ---
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "NHW"

# --- Limma ---
variances   <- apply(t(num_df_BAA_NHW), 1, var)
mat_filtree <- t(num_df_BAA_NHW)[variances > 0, ]

contrast_matrix <- makeContrasts(BAA - NHW, levels = design)

fit  <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters ---
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
)

# --- Top genes ---
top_genes <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(pvalue) %>%
  head(10)

# --- Volcano plot ---
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = case_when(
    -log10(pvalue) < 1.25 & logFC >= -1 & logFC <= 1      ~ "Black",
    -log10(pvalue) >= 1.25 & logFC >= -1 & logFC <= 1     ~ "Blue",
    -log10(pvalue) >= 1.25 & (logFC < -1 | logFC > 1)     ~ "Red",
    -log10(pvalue) < 1.25  & (logFC < -1 | logFC > 1)     ~ "Green",
    TRUE                                                    ~ "Fallback"
  )),
  alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Black"    = "black",
      "Blue"     = "blue",
      "Red"      = "red",
      "Green"    = "green4",
      "Fallback" = "gray60"
    ),
    labels = c(
      "Black"    = "Not significant",
      "Blue"     = "Significant & weak effect",
      "Red"      = "Significant & strong effect",
      "Green"    = "Not significant & strong effect"
    ),
    name = "Significance Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(p-value)",
    title = "BAA vs NHW"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(pvalue)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_NHW.csv", row.names = FALSE)
cat("Fichic sauvegarde : DEGs_BAA_vs_NHW.csv\n")

# ============================================================
# BAA vs nonBAA
# ============================================================

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Data preparation ---
condition <- df_BAA_nonBAA$condition
num_df_BAA_nonBAA <- df_BAA_nonBAA[, !colnames(df_BAA_nonBAA) %in% c("condition", "group")]

dim(num_df_BAA_nonBAA)
length(condition)

# --- Design ---
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)  # "BAA", "nonBAA"

# --- Limma ---
variances   <- apply(t(num_df_BAA_nonBAA), 1, var)
mat_filtree <- t(num_df_BAA_nonBAA)[variances > 0, ]

contrast_matrix <- makeContrasts(BAA - nonBAA, levels = design)

fit  <- lmFit(mat_filtree, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

resultats <- topTable(fit2, adjust = "BH", number = Inf)

# --- Parameters ---
fdr_threshold  <- 0.05
fc_threshold   <- 1
pval_threshold <- -log10(fdr_threshold)

# --- Volcano data ---
volcano_data <- data.frame(
  gene   = rownames(resultats),
  logFC  = resultats$logFC,
  pvalue = resultats$P.Value,
  FDR    = resultats$adj.P.Val,
  stringsAsFactors = FALSE
)

# --- Top genes ---
top_genes <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::arrange(pvalue) %>%
  head(10)

# --- Volcano plot ---
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = case_when(
    -log10(pvalue) < 1.25 & logFC >= -1 & logFC <= 1      ~ "Black",
    -log10(pvalue) >= 1.25 & logFC >= -1 & logFC <= 1     ~ "Blue",
    -log10(pvalue) >= 1.25 & (logFC < -1 | logFC > 1)     ~ "Red",
    -log10(pvalue) < 1.25  & (logFC < -1 | logFC > 1)     ~ "Green",
    TRUE                                                    ~ "Fallback"
  )),
  alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Black"    = "black",
      "Blue"     = "blue",
      "Red"      = "red",
      "Green"    = "green4",
      "Fallback" = "gray60"
    ),
    labels = c(
      "Black"    = "Not significant",
      "Blue"     = "Significant & weak effect",
      "Red"      = "Significant & strong effect",
      "Green"    = "Not significant & strong effect"
    ),
    name = "Significance Categories"
  ) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data         = top_genes,
    aes(label    = gene),
    size         = 3.5,
    box.padding  = 0.5,
    max.overlaps = Inf,
    color        = "black"
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(p-value)",
    title = "BAA vs nonBAA"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(volcano_plot)

# --- DEGs list ---
DEGs <- volcano_data %>%
  dplyr::filter(pvalue < fdr_threshold & abs(logFC) > fc_threshold) %>%
  dplyr::mutate(
    status = case_when(
      logFC >  fc_threshold ~ "Upregulated",
      logFC < -fc_threshold ~ "Downregulated"
    )
  ) %>%
  dplyr::select(gene, status, logFC, pvalue, FDR) %>%
  dplyr::arrange(pvalue)

print(DEGs, n = Inf)
cat("\n--- Summary ---\n")
cat("Total DEGs    :", nrow(DEGs), "\n")
cat("Upregulated   :", sum(DEGs$status == "Upregulated"), "\n")
cat("Downregulated :", sum(DEGs$status == "Downregulated"), "\n")

write.csv(DEGs, "DEGs_BAA_vs_nonBAA.csv", row.names = FALSE)
cat("Fichier sauvegarde : DEGs_BAA_vs_nonBAA.csv\n")

								######### II. GENE ANNOTATION####################

library(clusterProfiler)
library(org.Hs.eg.db)

condition1_up_down_UP <- condition1_up_down$gene[condition1_up_down$logFC > 0]
condition1_up_down_DOWN <- condition1_up_down$gene[condition1_up_down$logFC < 0]

gene_df_up <- bitr(condition1_up_down_UP, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_list_up <- gene_df_up$ENTREZID

gene_df_down <- bitr(condition1_up_down_DOWN, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_list_down <- gene_df_down$ENTREZID


###########################################
condition2_up_down_UP <- condition2_up_down$gene[condition2_up_down$logFC > 0]
condition2_up_down_DOWN <- condition2_up_down$gene[condition2_up_down$logFC < 0]

###########################################
condition3_up_down_UP <- condition3_up_down$gene[condition3_up_down$logFC > 0]
condition3_up_down_DOWN <- condition3_up_down$gene[condition3_up_down$logFC < 0]

###########################################

condition123_up_down_UP <- c(condition3_up_down_UP, condition1_up_down_UP, condition2_up_down_UP)
condition123_up_down_DOWN <- c(condition3_up_down_DOWN, condition1_up_down_DOWN, condition2_up_down_DOWN)

gene_df_up <- bitr(condition123_up_down_UP, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_list_up <- gene_df_up$ENTREZID

gene_df_down <- bitr(condition123_up_down_DOWN, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_list_down <- gene_df_down$ENTREZID

ego <- enrichGO(
  gene          = gene_list_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert back to gene symbols
)

ego_down <- enrichGO(
  gene          = gene_list_down,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert back to gene symbols
)
#####
ego <- enrichGO(
  gene          = gene_list_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert back to gene symbols
)

ego_down <- enrichGO(
  gene          = gene_list_down,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert back to gene symbols
)

# View results
head(ego)
dotplot(ego)

head(ego_down)
dotplot(ego_down)
#####
ego <- enrichGO(
  gene          = gene_list_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert back to gene symbols
)

ego_down <- enrichGO(
  gene          = gene_list_down,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # convert back to gene symbols
)


###

# View results
head(ego)
dotplot(ego)

head(ego_down)
dotplot(ego_down)

#####
ekegg <- enrichKEGG(
  gene          = gene_list_up,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2              # convert back to gene symbols
)


ekegg_down <- enrichKEGG(
  gene          = gene_list_down,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2              # convert back to gene symbols
)


head(ekegg)
dotplot(ekegg)

head(ekegg_down)
dotplot(ekegg_down)


####################################################
genelist_up <- rbind(condition3_up_down[condition3_up_down$logFC > 0,], condition2_up_down[condition2_up_down$logFC > 0,], condition1_up_down[condition1_up_down$logFC > 0,])
genelist_down <- rbind(condition3_up_down[condition3_up_down$logFC < 0,], condition2_up_down[condition2_up_down$logFC < 0,], condition1_up_down[condition1_up_down$logFC < 0,])

#write.table( condition1,"DEGs_BAA_ASN.txt",  quote = FALSE, row.names = FALSE, sep = '\t')

#write.table( condition2, "DEGs_BAA_HLW.txt", quote = FALSE, row.names = FALSE, sep = '\t')


genelist_up_new <- genelist_up[!duplicated(genelist_up$gene),]
genelist_down_new <- genelist_down[!duplicated(genelist_down$gene),]

write.table( genelist_up_new,"DEGs_BAA_others_up.txt",  quote = FALSE, row.names = FALSE, sep = '\t')

write.table( genelist_down_new, "DEGs_BAA_others_down.txt", quote = FALSE, row.names = FALSE, sep = '\t')

###################################################################################################
#DEGs
# Filtrer par p-value ajustée (FDR) < 0.05
Signif_DEGs <- results[results$adj.P.Val < 0.05, ]

# Voir combien de gènes sont identifiés
nrow(Signif_DEGs)

# Gènes UP (plus exprimés dans BPA) : logFC > 1 et p-adj < 0.05
DEGs_up <- results[results$adj.P.Val < 0.05 & results$logFC > 1, ]

# Gènes DOWN (moins exprimés dans BPA) : logFC < -1 et p-adj < 0.05
DEGs_down <- results[results$adj.P.Val < 0.05 & results$logFC < -1, ]

results$diff_status <- "Non-significatif"
results$diff_status[results$adj.P.Val < 0.05 & results$logFC > 1] <- "UP"
results$diff_status[results$adj.P.Val < 0.05 & results$logFC < -1] <- "DOWN"

# Compter le nombre de gènes dans chaque catégorie
table(results$diff_status)

# Export
write.csv(Signif_DEGs, "List_DEGs_BPA_vs_CTRL.csv")
head(List_DEGs_BPA_vs_CTRL.csv)
