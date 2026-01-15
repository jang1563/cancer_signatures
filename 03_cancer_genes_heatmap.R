# =============================================================================
# Cancer Gene Expression Heatmaps
# =============================================================================
#
# Purpose: Generate heatmaps showing differential expression of oncogenes and
#          tumor suppressor genes across different cell types at various
#          post-spaceflight timepoints
#
# Author: [Your Name]
# Date: April 2025
# Repository: https://github.com/[username]/cancer_signatures
#
# Input:
#   - Differential expression results from scRNA-seq analysis
#   - Gene lists: oncogene/TSG curated list, OncoKB database
#
# Output:
#   - Heatmaps (PDF/TIFF) showing log2FC across cell types
#   - Intermediate CSV files for downstream analysis
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(tidyr)
library(ggsignif)
library(circlize)
library(ComplexHeatmap)

# -----------------------------------------------------------------------------
# 2. Create Output Directories
# -----------------------------------------------------------------------------
dir.create('output/2025_01_03', recursive = TRUE, showWarnings = FALSE)
dir.create('intermediate/2025_04_24', recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3. Define Cancer Gene Lists
# -----------------------------------------------------------------------------

# Oncogenes and Tumor Suppressor Genes (curated list)
oncogene.tsg <- c(
  "ALOX15B", "CYP27B1", "HMGCS2", "PGD", "CNTN2", "TNFAIP3", "SERPINA1",
  "ALDH3A1", "ACLY", "ADIPOQ", "AGT", "AICDA", "NLRP3", "FOLR1", "AREG",
  "NPPA", "CYR61", "ANG", "ANGPTL4", "VIP", "IL10", "COL2A1", "TFAP2A",
  "GDNF", "BCAN", "ADRB1", "GREM1", "GRP", "IBSP", "TFF1", "MYB", "C5",
  "C6", "PSMA7", "MUC16", "S100A8", "CASP14", "CAV3", "CCL2", "CCL21",
  "NOV", "ACKR2", "CCR7", "CD276", "CD300A", "CD38", "CD80", "TNFRSF10D",
  "CD97", "CDKN1A", "CDX2", "CEACAM1", "IL12B", "CHI3L1", "CHIT1", "CKS2",
  "S100A9", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLU", "COL18A1", "CSTA",
  "CTSG", "CTSK", "CTSL", "CXCL10", "CXCL12", "CXCL5", "CXCL6", "DAB2",
  "POSTN", "DIO2", "DKK1", "DNER", "EGF", "EMP1", "ENPEP", "EPHA2",
  "TNFRSF10B", "EFNA1", "FAS", "LTF", "FCGR3A", "TNFRSF6B", "FGF18",
  "FGF7", "FGF9", "FIGF", "VEGFD", "FLT3LG", "FN1", "FOLH1", "FPR2",
  "GAL", "GAS6", "GDF15", "GJB2", "MARCO", "HAS2", "HDGF", "GPC3",
  "HIF1A", "HMGA1", "HOXB7", "HOXB9", "SDC1", "IFI27", "CXCL9", "IFITM3",
  "IGF1", "IGF2", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6",
  "IGFBP7", "IL17A", "IL17B", "IL18", "IL1A", "IL1B", "IL1RN", "IL27",
  "IL32", "IL4", "IL6", "INSL4", "CXCL11", "IRF1", "ISG15", "ITGB4",
  "JAG1", "PXN", "KIT", "KLK10", "KLK6", "GPNMB", "INHBA", "LAMA5",
  "LAMB3", "LAMC2", "LGALS3", "LGALS7", "LGALS3BP", "LIF", "LTA", "LTB",
  "LTBP2", "ADRB2", "LY6E", "MIF", "CXCL1", "CXCL2", "CXCL3", "CCL7",
  "MMP1", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3",
  "MMP7", "MMP9", "MSMB", "MT2A", "GNAS", "MUC1", "MX1", "NAMPT",
  "CCL13", "CCL11", "NCR1", "NGFR", "NID1", "NME1", "NOS2", "NPY",
  "OAS1", "CCL18", "OLFM4", "OPN3", "OXTR", "P4HA1", "P4HA2", "PAM",
  "PAPPA", "PCDH7", "SERPINE1", "SERPINE2", "PECAM1", "CFD", "PF4",
  "PLAUR", "PLA2G2A", "PLAT", "PLAU", "CXCL4", "PMP22", "PPBP", "STC2",
  "CDKN2A", "LGALS1", "LBP", "PRSS23", "PRTN3", "PTGS1", "PTGS2",
  "PTH", "PTMA", "PTN", "PTX3", "RAC2", "CCL5", "REG1A", "REN", "CCL3",
  "S100A12", "S100A4", "S100A6", "S100A7", "S100B", "S100P", "SAA1",
  "SCG2", "MIF4GD", "SCGB2A1", "SCGB2A2", "SCTR", "SDCBP", "SEL1L3",
  "SELE", "SELL", "SELP", "SEMA4D", "SERPINA3", "SERPINB13", "SERPINB1",
  "SERPINB2", "SERPINB5", "SERPINB9", "SERPIND1", "CD320", "SLCO4C1",
  "SPINK1", "SPOCK1", "SPP1", "RARRES2", "CCL4", "TACSTD2", "TF",
  "TFPI", "TFPI2", "TGM2", "THBD", "THBS1", "THBS2", "TIMP1", "TIMP2",
  "TIMP3", "TIMP4", "TINAGL1", "TNC", "TNF", "TNFAIP2", "TNFRSF11B",
  "TNFRSF1B", "TNFSF10", "TNFSF13B", "TNFSF4", "TNFSF9", "TP73",
  "TPSAB1", "TPSB2", "TSC22D1", "TSPAN8", "ANXA1", "TYROBP", "UCN",
  "VCAN", "VCAM1", "VEGFA", "VEGFC", "VGF", "VTN", "CCL23", "WIF1",
  "CCL24", "CCL26", "CXCL13"
)

# OncoKB gene list (abbreviated for space - full list in actual analysis)
oncokb <- c(
  "ABL1", "ACVR1", "AGO2", "AKT1", "AKT2", "AKT3", "ALK", "ALOX12B",
  "AMER1", "APC", "AR", "ARAF", "ARID1A", "ARID1B", "ARID2", "ARID5B",
  "ASXL1", "ASXL2", "ATM", "ATR", "ATRX", "AURKA", "AURKB", "AXIN1",
  "AXIN2", "AXL", "B2M", "BAP1", "BARD1", "BBC3", "BCL10", "BCL2",
  "BCL2L1", "BCL2L11", "BCL2L12", "BCL6", "BCOR", "BIRC3", "BLM",
  "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD4", "BRIP1", "BTG1", "BTK",
  "CALR", "CARD11", "CARM1", "CASP8", "CBFB", "CBL", "CCND1", "CCND2",
  "CCND3", "CCNE1", "CD274", "CD276", "CD79A", "CD79B", "CDC73", "CDH1",
  "CDK12", "CDK4", "CDK6", "CDK8", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B",
  "CDKN2C", "CEBPA", "CENPA", "CHEK1", "CHEK2", "CIC", "CREBBP", "CRKL",
  "CRLF2", "CSDE1", "CSF1R", "CSF3R", "CTCF", "CTLA4", "CTNNB1", "CUL3",
  "CXCR4", "CYLD", "DAXX", "DDR1", "DDR2", "DDX3X", "DDX41", "DICER1",
  "DNMT1", "DNMT3A", "DNMT3B", "E2F3", "EED", "EGFR", "EIF1AX", "ELF3",
  "EP300", "EPCAM", "EPHA2", "ERBB2", "ERBB3", "ERBB4", "ERG", "ESR1",
  "ETS1", "ETV1", "ETV4", "ETV5", "ETV6", "EZH1", "EZH2", "FAM175A",
  "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL",
  "FANCM", "FAS", "FAT1", "FBXW7", "FGF10", "FGF19", "FGF3", "FGF4",
  "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FH", "FLCN", "FLI1", "FLT1",
  "FLT3", "FLT4", "FOXA1", "FOXL2", "FOXO1", "FOXP1", "GATA1", "GATA2",
  "GATA3", "GLI1", "GNA11", "GNA13", "GNAQ", "GNAS", "GREM1", "GRIN2A",
  "GRM3", "GSK3B", "H3C14", "H3C2", "H3C3", "HGF", "HLA-A", "HLA-B",
  "HNF1A", "HOXB13", "HRAS", "IDH1", "IDH2", "IGF1R", "IKBKE", "IKZF1",
  "IL10", "IL7R", "INHA", "INHBA", "INPP4B", "INSR", "IRF4", "IRS1",
  "IRS2", "JAK1", "JAK2", "JAK3", "JUN", "KDM5A", "KDM5C", "KDM6A",
  "KDR", "KEAP1", "KIF5B", "KIT", "KLF4", "KMT2A", "KMT2B", "KMT2C",
  "KMT2D", "KRAS", "LATS1", "LATS2", "LYN", "LZTR1", "MAF", "MALT1",
  "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAP3K14", "MAPK1",
  "MAPK3", "MAX", "MCL1", "MDM2", "MDM4", "MED12", "MEF2B", "MEN1",
  "MERTK", "MET", "MGA", "MITF", "MLH1", "MLH3", "MPL", "MRE11A", "MSH2",
  "MSH3", "MSH6", "MTOR", "MUTYH", "MYC", "MYCL", "MYCN", "MYD88",
  "MYOD1", "NBN", "NCOR1", "NF1", "NF2", "NFE2L2", "NFKBIA", "NKX2-1",
  "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NPM1", "NRAS", "NRG1", "NSD1",
  "NSD2", "NSD3", "NTRK1", "NTRK2", "NTRK3", "PAK1", "PAK3", "PALB2",
  "PARP1", "PAX3", "PAX5", "PAX8", "PBRM1", "PDCD1", "PDGFRA", "PDGFRB",
  "PDPK1", "PGR", "PHF6", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG",
  "PIK3R1", "PIK3R2", "PIM1", "PLCG2", "PMAIP1", "PMS1", "PMS2", "POLD1",
  "POLE", "POT1", "PPARG", "PPM1D", "PPP2R1A", "PRDM1", "PREX2", "PRKAR1A",
  "PTCH1", "PTEN", "PTPN11", "PTPRD", "PTPRS", "PTPRT", "RAC1", "RAD21",
  "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAF1", "RARA", "RASA1",
  "RB1", "RBM10", "RECQL4", "REL", "RET", "RHEB", "RHOA", "RICTOR",
  "RIT1", "RNF43", "ROS1", "RPTOR", "RUNX1", "RUNX1T1", "RXRA", "SDHA",
  "SDHAF2", "SDHB", "SDHC", "SDHD", "SETD2", "SF3B1", "SGK1", "SH2B3",
  "SH2D1A", "SLX4", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1",
  "SMC1A", "SMC3", "SMO", "SOCS1", "SOX10", "SOX17", "SOX2", "SOX9",
  "SPEN", "SPOP", "SPRED1", "SRC", "SRSF2", "SS18", "STAG2", "STAT3",
  "STAT5A", "STAT5B", "STAT6", "STK11", "SUZ12", "SYK", "TBL1XR1", "TBX3",
  "TCEB1", "TCF3", "TCF7L2", "TEK", "TERT", "TET1", "TET2", "TGFBR1",
  "TGFBR2", "TMEM127", "TNFAIP3", "TNFRSF14", "TOP1", "TOP2A", "TP53",
  "TP53BP1", "TP63", "TRAF2", "TRAF3", "TRAF7", "TSC1", "TSC2", "TSHR",
  "U2AF1", "VEGFA", "VHL", "WT1", "WWTR1", "XBP1", "XIAP", "XPO1",
  "YAP1", "YES1", "ZFHX3", "ZNF217", "ZNF703"
)

# -----------------------------------------------------------------------------
# 4. Define Cell Types
# -----------------------------------------------------------------------------
cell_list <- c("pbmc", "CD4_T", "CD8_T", "other_T", "B",
               "NK", "CD14_Mono", "CD16_Mono", "DC", "other")
all_list <- cell_list
gene.list <- c('oncogene.tsg', 'oncokb')

# -----------------------------------------------------------------------------
# 5. Load Differential Expression Data
# -----------------------------------------------------------------------------
for (i in cell_list) {
  # R+1 (immediately post-flight)
  x1 <- read.csv(paste0(
    "../../../Differential_expression/with_filter/celltype/pval/update/all.list/update/DEGs_",
    i, "_R+1_preflight.csv"
  ))

  # R+45&R+82 (long-term post-flight)
  x2 <- read.csv(paste0(
    "../../../Differential_expression/with_filter/celltype/pval/update/all.list/update/DEGs_",
    i, "_R+45&R+82_preflight.csv"
  ))

  assign(paste0(i, '.R1'), x1)
  assign(paste0(i, '.R45&R82'), x2)
}

# -----------------------------------------------------------------------------
# 6. Filter for Cancer Genes
# -----------------------------------------------------------------------------
for (i in all_list) {
  x1 <- get(paste0(i, '.R1'))
  x2 <- get(paste0(i, '.R45&R82'))

  for (j in gene.list) {
    A <- get(j)

    y1 <- x1 %>% filter(X %in% A)
    y2 <- x2 %>% filter(X %in% A)

    # Sort by gene name for consistent ordering
    y1 <- y1[order(y1$X, decreasing = TRUE), ]
    y2 <- y2[order(y2$X, decreasing = TRUE), ]

    assign(paste0(i, ".", j, ".R1"), y1)
    assign(paste0(i, ".", j, ".R45_R82"), y2)
  }
}

# -----------------------------------------------------------------------------
# 7. Create Heatmap Data Matrices
# -----------------------------------------------------------------------------

# Oncogene/TSG - R+1
celltype.oncogene.tsg.R1 <- data.frame(
  CD4_T = CD4_T.oncogene.tsg.R1$avg_log2FC,
  CD8_T = CD8_T.oncogene.tsg.R1$avg_log2FC,
  other_T = other_T.oncogene.tsg.R1$avg_log2FC,
  B = B.oncogene.tsg.R1$avg_log2FC,
  NK = NK.oncogene.tsg.R1$avg_log2FC,
  CD14_Mono = CD14_Mono.oncogene.tsg.R1$avg_log2FC,
  CD16_Mono = CD16_Mono.oncogene.tsg.R1$avg_log2FC,
  DC = DC.oncogene.tsg.R1$avg_log2FC,
  other = other.oncogene.tsg.R1$avg_log2FC
)

# Oncogene/TSG - R+45&R+82
celltype.oncogene.tsg.R45_R82 <- data.frame(
  CD4_T = CD4_T.oncogene.tsg.R45_R82$avg_log2FC,
  CD8_T = CD8_T.oncogene.tsg.R45_R82$avg_log2FC,
  other_T = other_T.oncogene.tsg.R45_R82$avg_log2FC,
  B = B.oncogene.tsg.R45_R82$avg_log2FC,
  NK = NK.oncogene.tsg.R45_R82$avg_log2FC,
  CD14_Mono = CD14_Mono.oncogene.tsg.R45_R82$avg_log2FC,
  CD16_Mono = CD16_Mono.oncogene.tsg.R45_R82$avg_log2FC,
  DC = DC.oncogene.tsg.R45_R82$avg_log2FC,
  other = other.oncogene.tsg.R45_R82$avg_log2FC
)

# OncoKB - R+1
celltype.oncokb.R1 <- data.frame(
  CD4_T = CD4_T.oncokb.R1$avg_log2FC,
  CD8_T = CD8_T.oncokb.R1$avg_log2FC,
  other_T = other_T.oncokb.R1$avg_log2FC,
  B = B.oncokb.R1$avg_log2FC,
  NK = NK.oncokb.R1$avg_log2FC,
  CD14_Mono = CD14_Mono.oncokb.R1$avg_log2FC,
  CD16_Mono = CD16_Mono.oncokb.R1$avg_log2FC,
  DC = DC.oncokb.R1$avg_log2FC,
  other = other.oncokb.R1$avg_log2FC
)

# OncoKB - R+45&R+82
celltype.oncokb.R45_R82 <- data.frame(
  CD4_T = CD4_T.oncokb.R45_R82$avg_log2FC,
  CD8_T = CD8_T.oncokb.R45_R82$avg_log2FC,
  other_T = other_T.oncokb.R45_R82$avg_log2FC,
  B = B.oncokb.R45_R82$avg_log2FC,
  NK = NK.oncokb.R45_R82$avg_log2FC,
  CD14_Mono = CD14_Mono.oncokb.R45_R82$avg_log2FC,
  CD16_Mono = CD16_Mono.oncokb.R45_R82$avg_log2FC,
  DC = DC.oncokb.R45_R82$avg_log2FC,
  other = other.oncokb.R45_R82$avg_log2FC
)

# Set row names
row.names(celltype.oncogene.tsg.R1) <- CD4_T.oncogene.tsg.R1$X
row.names(celltype.oncogene.tsg.R45_R82) <- CD4_T.oncogene.tsg.R45_R82$X
row.names(celltype.oncokb.R1) <- CD4_T.oncokb.R1$X
row.names(celltype.oncokb.R45_R82) <- CD4_T.oncokb.R45_R82$X

# -----------------------------------------------------------------------------
# 8. Save Intermediate Results
# -----------------------------------------------------------------------------
write.csv(celltype.oncokb.R1, 'intermediate/2025_04_24/I4.oncokb.R1.csv')
write.csv(celltype.oncokb.R45_R82, 'intermediate/2025_04_24/I4.oncokb.R45_R82.csv')

# -----------------------------------------------------------------------------
# 9. Define Color Palette for Heatmaps
# -----------------------------------------------------------------------------
paletteLength <- 200
low <- colorRampPalette(c("black", "navy", "white"))(100)
high <- colorRampPalette(c("white", "firebrick1", "firebrick2", "firebrick3", "firebrick4"))(100)
combined <- c(low, high)
myColor2 <- combined

# Calculate color breaks for each dataset
myBreaks.oncogene.tsg.R1 <- c(
  seq(min(celltype.oncogene.tsg.R1), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(celltype.oncogene.tsg.R1) / paletteLength, max(celltype.oncogene.tsg.R1),
      length.out = floor(paletteLength / 2))
)

myBreaks.oncogene.tsg.R45_R82 <- c(
  seq(min(celltype.oncogene.tsg.R45_R82), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(celltype.oncogene.tsg.R45_R82) / paletteLength, max(celltype.oncogene.tsg.R45_R82),
      length.out = floor(paletteLength / 2))
)

myBreaks.oncokb.R1 <- c(
  seq(min(celltype.oncokb.R1), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(celltype.oncokb.R1) / paletteLength, max(celltype.oncokb.R1),
      length.out = floor(paletteLength / 2))
)

myBreaks.oncokb.R45_R82 <- c(
  seq(min(celltype.oncokb.R45_R82), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(celltype.oncokb.R45_R82) / paletteLength, max(celltype.oncokb.R45_R82),
      length.out = floor(paletteLength / 2))
)

# -----------------------------------------------------------------------------
# 10. Generate Heatmaps - Oncogene/TSG R+1
# -----------------------------------------------------------------------------
pdf("output/2025_01_03/I4.oncogene.tsg.R1.ver2.pdf", width = 6, height = 20)

Heatmap(
  as.matrix(celltype.oncogene.tsg.R1),
  column_labels = colnames(celltype.oncogene.tsg.R1),
  border = TRUE,
  use_raster = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  name = "Log2FC",
  col = colorRamp2(c(min(celltype.oncogene.tsg.R1), 0, max(celltype.oncogene.tsg.R1)),
                   c("navy", "white", "firebrick3")),
  column_title = "Oncogene/TSG Expression - R+1 vs Preflight",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

dev.off()

# TIFF version
tiff("output/2025_01_03/I4.oncogene.tsg.R1.ver2.tiff", width = 6, height = 20,
     units = 'in', res = 300, compression = 'lzw')

Heatmap(
  as.matrix(celltype.oncogene.tsg.R1),
  column_labels = colnames(celltype.oncogene.tsg.R1),
  border = TRUE,
  use_raster = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  name = "Log2FC",
  col = colorRamp2(c(min(celltype.oncogene.tsg.R1), 0, max(celltype.oncogene.tsg.R1)),
                   c("navy", "white", "firebrick3")),
  column_title = "Oncogene/TSG Expression - R+1 vs Preflight",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

dev.off()

# -----------------------------------------------------------------------------
# 11. Generate Heatmaps - Oncogene/TSG R+45&R+82
# -----------------------------------------------------------------------------
pdf("output/2025_01_03/I4.oncogene.tsg.R45_R82.ver2.pdf", width = 6, height = 20)

Heatmap(
  as.matrix(celltype.oncogene.tsg.R45_R82),
  column_labels = colnames(celltype.oncogene.tsg.R45_R82),
  border = TRUE,
  use_raster = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  name = "Log2FC",
  col = colorRamp2(c(min(celltype.oncogene.tsg.R45_R82), 0, max(celltype.oncogene.tsg.R45_R82)),
                   c("navy", "white", "firebrick3")),
  column_title = "Oncogene/TSG Expression - R+45&R+82 vs Preflight",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

dev.off()

# -----------------------------------------------------------------------------
# 12. Generate Heatmaps - OncoKB R+1
# -----------------------------------------------------------------------------
pdf("output/2025_01_03/I4.oncokb.R1.ver2.pdf", width = 6, height = 28)

Heatmap(
  as.matrix(celltype.oncokb.R1),
  column_labels = colnames(celltype.oncokb.R1),
  border = TRUE,
  use_raster = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  name = "Log2FC",
  col = colorRamp2(c(min(celltype.oncokb.R1), 0, max(celltype.oncokb.R1)),
                   c("navy", "white", "firebrick3")),
  column_title = "OncoKB Gene Expression - R+1 vs Preflight",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

dev.off()

# -----------------------------------------------------------------------------
# 13. Generate Heatmaps - OncoKB R+45&R+82
# -----------------------------------------------------------------------------
pdf("output/2025_01_03/I4.oncokb.R45_R82.ver2.pdf", width = 6, height = 28)

Heatmap(
  as.matrix(celltype.oncokb.R45_R82),
  column_labels = colnames(celltype.oncokb.R45_R82),
  border = TRUE,
  use_raster = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  name = "Log2FC",
  col = colorRamp2(c(min(celltype.oncokb.R45_R82), 0, max(celltype.oncokb.R45_R82)),
                   c("navy", "white", "firebrick3")),
  column_title = "OncoKB Gene Expression - R+45&R+82 vs Preflight",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

dev.off()

# -----------------------------------------------------------------------------
# 14. Session Info
# -----------------------------------------------------------------------------
sessionInfo()
