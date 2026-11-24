#######################################################################################
#			A. DEG IDENTIFICATION - BAA vs Others
## ######################################################################################

						## MISSING SECTION !!!! YANN'S TO COMPLETE

#########################################################################################################
#			B. BAA vs Others DEG ANNOTATIONS - GO & KEGG - Significant DEG Filtered by FDR
## ########################################################################################################

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

path <- "/Users/drpierred.sarr/Desktop/TCGA-BRCA/DEGs Yann/"

# ============================================================
# Read and filter input files
# ============================================================

all_up <- read.table(
  paste0(path, "DEGs_BAA_others_Upregulated.txt"),
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

all_down <- read.table(
  paste0(path, "DEGs_BAA_others_Downregulated.txt"),
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

# Check the columns (the uploaded files contain all upregulated/downregulated DEGs; annotation is performed using significant DEGs only)
cat("Colonnes UP   :", colnames(all_up), "\n")
cat("Colonnes DOWN :", colnames(all_down), "\n")
cat("Gènes totaux UP   :", nrow(all_up), "\n")
cat("Gènes totaux DOWN :", nrow(all_down), "\n")

#  --- Filter significant DEGs ---
DEGs_up_df <- all_up %>%
  dplyr::filter(FDR < 0.05 & logFC > 1)

DEGs_down_df <- all_down %>%
  dplyr::filter(FDR < 0.05 & logFC < -1)

cat("DEGs UP filtrés   :", nrow(DEGs_up_df), "\n")
cat("DEGs DOWN filtrés :", nrow(DEGs_down_df), "\n")

# --- Extract gene symbols ---
DEGs_up   <- DEGs_up_df$gene
DEGs_down <- DEGs_down_df$gene

# ============================================================
# Convert gene SYMBOLs to ENTREZID
# ============================================================

genes_up_entrez <- bitr(DEGs_up,
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org.Hs.eg.db)$ENTREZID

genes_down_entrez <- bitr(DEGs_down,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org.Hs.eg.db)$ENTREZID

cat("Up entrez   :", length(genes_up_entrez), "\n")
cat("Down entrez :", length(genes_down_entrez), "\n")

# ============================================================
# GO-BP
# ============================================================

go_bp_up <- enrichGO(gene          = genes_up_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)

go_bp_down <- enrichGO(gene          = genes_down_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2,
                       readable      = TRUE)

cat("GO-BP UP   :", nrow(go_bp_up), "\n")
cat("GO-BP DOWN :", nrow(go_bp_down), "\n")

dotplot(go_bp_up,   showCategory = 10, title = "GO-BP | Upregulated")
dotplot(go_bp_down, showCategory = 10, title = "GO-BP | Downregulated")

# ============================================================
# GO-CC
# ============================================================

go_cc_up <- enrichGO(gene          = genes_up_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)

go_cc_down <- enrichGO(gene          = genes_down_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2,
                       readable      = TRUE)

cat("GO-CC UP   :", nrow(go_cc_up), "\n")
cat("GO-CC DOWN :", nrow(go_cc_down), "\n")

dotplot(go_cc_up,   showCategory = 10, title = "GO-CC | Upregulated")
dotplot(go_cc_down, showCategory = 10, title = "GO-CC | Downregulated")

# ============================================================
# GO-MF
# ============================================================

go_mf_up <- enrichGO(gene          = genes_up_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)

go_mf_down <- enrichGO(gene          = genes_down_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2,
                       readable      = TRUE)

cat("GO-MF UP   :", nrow(go_mf_up), "\n")
cat("GO-MF DOWN :", nrow(go_mf_down), "\n")

dotplot(go_mf_up,   showCategory = 10, title = "GO-MF | Upregulated")
dotplot(go_mf_down, showCategory = 10, title = "GO-MF | Downregulated")
# ============================================================
# KEGG
# ============================================================

kegg_up <- enrichKEGG(gene          = genes_up_entrez,
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

kegg_down <- enrichKEGG(gene          = genes_down_entrez,
                        organism      = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2)

cat("KEGG UP   :", nrow(kegg_up), "\n")
cat("KEGG DOWN :", nrow(kegg_down), "\n")

dotplot(kegg_up,   showCategory = 10, title = "KEGG | Upregulated")
dotplot(kegg_down, showCategory = 10, title = "KEGG | Downregulated")

# ============================================================
# Export CSV
# ============================================================

if (!is.null(go_bp_up)   && nrow(go_bp_up) > 0)   write.csv(as.data.frame(go_bp_up),   paste0(path, "GO_BP_UP_BAA_vs_Others.csv"),   row.names = FALSE)
if (!is.null(go_bp_down) && nrow(go_bp_down) > 0)  write.csv(as.data.frame(go_bp_down), paste0(path, "GO_BP_DOWN_BAA_vs_Others.csv"), row.names = FALSE)
if (!is.null(go_cc_up)   && nrow(go_cc_up) > 0)   write.csv(as.data.frame(go_cc_up),   paste0(path, "GO_CC_UP_BAA_vs_Others.csv"),   row.names = FALSE)
if (!is.null(go_cc_down) && nrow(go_cc_down) > 0)  write.csv(as.data.frame(go_cc_down), paste0(path, "GO_CC_DOWN_BAA_vs_Others.csv"), row.names = FALSE)
if (!is.null(go_mf_up)   && nrow(go_mf_up) > 0)   write.csv(as.data.frame(go_mf_up),   paste0(path, "GO_MF_UP_BAA_vs_Others.csv"),   row.names = FALSE)
if (!is.null(go_mf_down) && nrow(go_mf_down) > 0)  write.csv(as.data.frame(go_mf_down), paste0(path, "GO_MF_DOWN_BAA_vs_Others.csv"), row.names = FALSE)
if (!is.null(kegg_up)    && nrow(kegg_up) > 0)    write.csv(as.data.frame(kegg_up),    paste0(path, "KEGG_UP_BAA_vs_Others.csv"),    row.names = FALSE)
if (!is.null(kegg_down)  && nrow(kegg_down) > 0)   write.csv(as.data.frame(kegg_down),  paste0(path, "KEGG_DOWN_BAA_vs_Others.csv"),  row.names = FALSE)

cat("\n--- Summary of generated files ---\n")
cat("PDF 1 : GO_BP_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 2 : GO_BP_Downregulated_BAA_vs_Others.pdf\n")
cat("PDF 3 : GO_CC_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 4 : GO_CC_Downregulated_BAA_vs_Others.pdf\n")
cat("PDF 5 : GO_MF_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 6 : GO_MF_Downregulated_BAA_vs_Others.pdf\n")
cat("PDF 7 : KEGG_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 8 : KEGG_Downregulated_BAA_vs_Others.pdf\n")

# NOTE: If the TXT files already contain only significant DEGs, remove the dplyr::filter() steps
# and directly assign the gene symbols:
# DEGs_up <- all_up$gene
# DEGs_down <- all_down$gene


#################################################################################################
#			C. BAA vs Others DEG ANNOTATIONS - GO & KEGG - Significant DEG Filtered by Nominal P-value
## ###############################################################################################
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

path <- "/Users/drpierred.sarr/Desktop/TCGA-BRCA/DEGs Yann/"

# ============================================================
# Read input files
# ============================================================

all_up <- read.table(
  paste0(path, "DEGs_BAA_others_Upregulated.txt"),
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

all_down <- read.table(
  paste0(path, "DEGs_BAA_others_Downregulated.txt"),
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

cat("Colonnes UP   :", colnames(all_up), "\n")
cat("Colonnes DOWN :", colnames(all_down), "\n")
cat("Gènes totaux UP   :", nrow(all_up), "\n")
cat("Gènes totaux DOWN :", nrow(all_down), "\n")

# ============================================================
# Filter significant DEGs based on p-values
# ============================================================

DEGs_up_df <- all_up %>%
  dplyr::filter(pvalue < 0.05 & logFC > 1)

DEGs_down_df <- all_down %>%
  dplyr::filter(pvalue < 0.05 & logFC < -1)

cat("DEGs UP filtrés   :", nrow(DEGs_up_df), "\n")
cat("DEGs DOWN filtrés :", nrow(DEGs_down_df), "\n")

DEGs_up   <- DEGs_up_df$gene
DEGs_down <- DEGs_down_df$gene

# ============================================================
# Convert gene SYMBOLs to ENTREZID
# ============================================================

genes_up_entrez <- bitr(DEGs_up,
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org.Hs.eg.db)$ENTREZID

genes_down_entrez <- bitr(DEGs_down,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org.Hs.eg.db)$ENTREZID

cat("Up entrez   :", length(genes_up_entrez), "\n")
cat("Down entrez :", length(genes_down_entrez), "\n")

# ============================================================
# GO-BP
# ============================================================

go_bp_up <- enrichGO(gene          = genes_up_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.20,
                     readable      = TRUE)

go_bp_down <- enrichGO(gene          = genes_down_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.20,
                       readable      = TRUE)

cat("GO-BP UP   :", nrow(go_bp_up), "\n")
cat("GO-BP DOWN :", nrow(go_bp_down), "\n")

dotplot(go_bp_up,   showCategory = 10, title = "GO-BP | Upregulated")
dotplot(go_bp_down, showCategory = 10, title = "GO-BP | Downregulated")
# ============================================================
# GO-CC
# ============================================================

go_cc_up <- enrichGO(gene          = genes_up_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.20,
                     readable      = TRUE)

go_cc_down <- enrichGO(gene          = genes_down_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.20,
                       readable      = TRUE)

cat("GO-CC UP   :", nrow(go_cc_up), "\n")
cat("GO-CC DOWN :", nrow(go_cc_down), "\n")

dotplot(go_cc_up,   showCategory = 10, title = "GO-CC | Upregulated")
dotplot(go_cc_down, showCategory = 10, title = "GO-CC | Downregulated")

# ============================================================
# GO-MF
# ============================================================

go_mf_up <- enrichGO(gene          = genes_up_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.20,
                     readable      = TRUE)

go_mf_down <- enrichGO(gene          = genes_down_entrez,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.20,
                       readable      = TRUE)

cat("GO-MF UP   :", nrow(go_mf_up), "\n")
cat("GO-MF DOWN :", nrow(go_mf_down), "\n")

dotplot(go_mf_up,   showCategory = 10, title = "GO-MF | Upregulated")
dotplot(go_mf_down, showCategory = 10, title = "GO-MF | Downregulated")
# ============================================================
# KEGG
# ============================================================

kegg_up <- enrichKEGG(gene          = genes_up_entrez,
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.20)

kegg_down <- enrichKEGG(gene          = genes_down_entrez,
                        organism      = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.20)

cat("KEGG UP   :", nrow(kegg_up), "\n")
cat("KEGG DOWN :", nrow(kegg_down), "\n")

dotplot(kegg_up,   showCategory = 10, title = "KEGG | Upregulated")
dotplot(kegg_down, showCategory = 10, title = "KEGG | Downregulated")

# ============================================================
# Export CSV
# ============================================================

if (!is.null(go_bp_up)   && nrow(go_bp_up) > 0)   write.csv(as.data.frame(go_bp_up),   paste0(path, "GO_BP_UP_BAA_vs_Others.csv"),   row.names = FALSE)
if (!is.null(go_bp_down) && nrow(go_bp_down) > 0)  write.csv(as.data.frame(go_bp_down), paste0(path, "GO_BP_DOWN_BAA_vs_Others.csv"), row.names = FALSE)
if (!is.null(go_cc_up)   && nrow(go_cc_up) > 0)   write.csv(as.data.frame(go_cc_up),   paste0(path, "GO_CC_UP_BAA_vs_Others.csv"),   row.names = FALSE)
if (!is.null(go_cc_down) && nrow(go_cc_down) > 0)  write.csv(as.data.frame(go_cc_down), paste0(path, "GO_CC_DOWN_BAA_vs_Others.csv"), row.names = FALSE)
if (!is.null(go_mf_up)   && nrow(go_mf_up) > 0)   write.csv(as.data.frame(go_mf_up),   paste0(path, "GO_MF_UP_BAA_vs_Others.csv"),   row.names = FALSE)
if (!is.null(go_mf_down) && nrow(go_mf_down) > 0)  write.csv(as.data.frame(go_mf_down), paste0(path, "GO_MF_DOWN_BAA_vs_Others.csv"), row.names = FALSE)
if (!is.null(kegg_up)    && nrow(kegg_up) > 0)    write.csv(as.data.frame(kegg_up),    paste0(path, "KEGG_UP_BAA_vs_Others.csv"),    row.names = FALSE)
if (!is.null(kegg_down)  && nrow(kegg_down) > 0)   write.csv(as.data.frame(kegg_down),  paste0(path, "KEGG_DOWN_BAA_vs_Others.csv"),  row.names = FALSE)

cat("\n--- Summary of generated files ---\n")
cat("PDF 1 : GO_BP_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 2 : GO_BP_Downregulated_BAA_vs_Others.pdf\n")
cat("PDF 3 : GO_CC_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 4 : GO_CC_Downregulated_BAA_vs_Others.pdf\n")
cat("PDF 5 : GO_MF_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 6 : GO_MF_Downregulated_BAA_vs_Others.pdf\n")
cat("PDF 7 : KEGG_Upregulated_BAA_vs_Others.pdf\n")
cat("PDF 8 : KEGG_Downregulated_BAA_vs_Others.pdf\n")















