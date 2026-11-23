#R-script
#Cancer subtyping using PAM50 model
#Resource: https://www.bioconductor.org/packages/release/bioc/manuals/genefu/man/genefu.pdf
#Dr.NaziaPathan
#Dec15th_2025

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

#Install the package Devtools
install.packages("devtools")
library(devtools)

# Install genefu from GitHub
install_github("bhklab/genefu")

BiocManager::install("openxlsx")

library(openxlsx)
library(genefu)
library(dplyr)
library(data.table)
library(pheatmap)

data(pam50)
data(pam50.robust)

#load data

file_path <- "/Users/drpierred.sarr/Desktop/PAM50-BAA.xlsx"

data <- read.xlsx(file_path, sheet = 1)

head(data)

#data structuring
expression_data=as.matrix(t(data[,-1]))
colnames(expression_data)<-data[[1]]
gene_info <- select(data, c("Gene"))

#main subtyping function

pam50_subtype_predictions <- molecular.subtyping(
  sbt.model = "pam50",
  data = expression_data,
  annot = gene_info,
  do.mapping = FALSE)

# Display the subtypes
as.data.frame(pam50_subtype_predictions$subtype)

# Display the subtype probabilities
as.data.frame(pam50_subtype_predictions$subtype.proba)


# Display the subtypes predictions
as.data.frame(pam50_subtype_predictions$subtype.crisp)

# Display the subtypes
m=as.data.frame(pam50_subtype_predictions$subtype.proba)
pheatmap(m, scale = "row", colorRampPalette(c("navy", "white", "#FF1493"))(75))



