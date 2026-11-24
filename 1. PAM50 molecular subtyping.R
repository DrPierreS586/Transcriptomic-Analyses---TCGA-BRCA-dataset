# #######################################################################################################################
# 				TCGA-BRCA MOLECULAR SUBTYPING USING THE PAM50 MODEL (Script for BAA Patients only)
# #######################################################################################################################

# ============================================================
# INSTALL AND LOAD REQUIRED PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Install devtools for GitHub package installation
install.packages("devtools")
library(devtools)

# Install genefu from GitHub
install_github("bhklab/genefu")

BiocManager::install("openxlsx")

# Load required libraries
library(openxlsx)
library(genefu)
library(dplyr)
library(data.table)
library(pheatmap)

# ============================================================
# LOAD PAM50 REFERENCE MODELS
# ============================================================
data(pam50)
data(pam50.robust)

# ============================================================
# LOAD GENE EXPRESSION DATA
# ============================================================

file_path <- "/Users/drpierred.sarr/Desktop/PAM50-BAA.xlsx"

data <- read.xlsx(file_path, sheet = 1)
head(data)

# ============================================================
# DATA PREPROCESSING AND STRUCTURING
# ============================================================
# Convert the gene expression data into a matrix
# Rows represent samples and columns represent genes

expression_data=as.matrix(t(data[,-1]))

# Assign gene names as column names
colnames(expression_data)<-data[[1]]

# Extract gene annotation information
gene_info <- select(data, c("Gene"))

# ============================================================
# PAM50 MOLECULAR SUBTYPE PREDICTION
# ============================================================
# Perform molecular subtyping using the PAM50 classifier

pam50_subtype_predictions <- molecular.subtyping(
  sbt.model = "pam50",
  data = expression_data,
  annot = gene_info,
  do.mapping = FALSE)

# ============================================================
# EXTRACT AND DISPLAY MOLECULAR SUBTYPE PREDICTIONS
# ============================================================
# Display the predicted molecular subtypes

as.data.frame(pam50_subtype_predictions$subtype)

# ============================================================
# EXTRACT AND DISPLAY SUBTYPE PROBABILITIES
# ============================================================

# Display the probability scores assigned to each PAM50 subtype

as.data.frame(pam50_subtype_predictions$subtype.proba)


# ============================================================
# EXTRACT AND DISPLAY FINAL SUBTYPE CLASSIFICATIONS
# ============================================================
# Display the final categorical subtype predictions

as.data.frame(pam50_subtype_predictions$subtype.crisp)

# ============================================================
# VISUALIZE PAM50 SUBTYPE PREDICTION PROBABILITIES
# ============================================================
# Convert subtype probabilities into a data frame

m=as.data.frame(pam50_subtype_predictions$subtype.proba)

# Generate a heatmap of subtype prediction probabilities
# Row-wise scaling is applied to facilitate visualization
pheatmap(
  m,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "#FF1493"))(75)
)

