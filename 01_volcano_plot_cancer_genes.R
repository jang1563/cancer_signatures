# =============================================================================
# Volcano Plot Analysis for Cancer-Related Gene Expression
# =============================================================================
#
# Purpose: Generate volcano plots for differential gene expression analysis
#          focusing on oncogenes and tumor suppressor genes in various cell types
#
# Author: [Your Name]
# Date: January 2025
# Repository: https://github.com/[username]/cancer_signatures
#
# Input:
#   - Differential expression results from scRNA-seq analysis
#   - Gene lists: oncogene/TSG list, OncoKB database
#
# Output:
#   - Volcano plots (PDF/TIFF) for each cell type and gene set
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
library(fgsea)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(ensembldb)
library(tidyr)
library(Seurat)
library(GenomicRanges)
library(ggsignif)
library(stringr)
library(ggpubr)
library(patchwork)

# Set global theme
theme_set(theme_classic(base_size = 20))

# -----------------------------------------------------------------------------
# 2. Create Output Directories
# -----------------------------------------------------------------------------
dir.create('output/2025_01_03/', recursive = TRUE, showWarnings = FALSE)

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

# OncoKB gene list
oncokb <- c(
  "ABL1", "ACVR1", "AGO2", "AKT1", "AKT2", "AKT3", "ALK", "ALOX12B",
  "AMER1", "ANKRD11", "APC", "AR", "ARAF", "ARID1A", "ARID1B", "ARID2",
  "ARID5B", "ASXL1", "ASXL2", "ATM", "ATR", "ATRX", "AURKA", "AURKB",
  "AXIN1", "AXIN2", "AXL", "B2M", "BAP1", "BARD1", "BBC3", "BCL10",
  "BCL2", "BCL2L1", "BCL2L11", "BCL2L12", "BCL6", "BCOR", "BIRC3",
  "BLM", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD4", "BRIP1", "BTG1",
  "BTK", "CALR", "CARD11", "CARM1", "CASP8", "CBFB", "CBL", "CCND1",
  "CCND2", "CCND3", "CCNE1", "CD274", "CD276", "CD79A", "CD79B",
  "CDC73", "CDH1", "CDK12", "CDK4", "CDK6", "CDK8", "CDKN1A", "CDKN1B",
  "CDKN2A", "CDKN2B", "CDKN2C", "CEBPA", "CENPA", "CHD4", "CHEK1",
  "CHEK2", "CIC", "CMPK1", "CREBBP", "CRKL", "CRLF2", "CSDE1", "CSF1R",
  "CSF3R", "CTCF", "CTLA4", "CTNNB1", "CUL3", "CXCR4", "CYLD", "CYP17A1",
  "DAXX", "DCUN1D1", "DDR1", "DDR2", "DDX3X", "DDX41", "DICER1",
  "DIS3", "DNAJB1", "DNMT1", "DNMT3A", "DNMT3B", "E2F3", "EED",
  "EGFL7", "EGFR", "EIF1AX", "EIF4A2", "EIF4E", "ELF3", "ELOC", "EP300",
  "EPCAM", "EPHA2", "EPHA3", "EPHA5", "EPHA7", "EPHB1", "EPHB4",
  "ERBB2", "ERBB3", "ERBB4", "ERCC2", "ERCC3", "ERCC4", "ERCC5",
  "ERF", "ERG", "ERRFI1", "ESR1", "ETS1", "ETV1", "ETV4", "ETV5",
  "ETV6", "EZH1", "EZH2", "FAM175A", "FAM46C", "FAM58A", "FANCA",
  "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL",
  "FANCM", "FAS", "FAT1", "FBXO11", "FBXW7", "FES", "FGF10", "FGF14",
  "FGF19", "FGF23", "FGF3", "FGF4", "FGF6", "FGFR1", "FGFR2", "FGFR3",
  "FGFR4", "FH", "FHIT", "FLCN", "FLI1", "FLT1", "FLT3", "FLT4",
  "FOXA1", "FOXL2", "FOXO1", "FOXP1", "FRS2", "FUBP1", "FYN", "GATA1",
  "GATA2", "GATA3", "GID4", "GLI1", "GNA11", "GNA13", "GNAQ", "GNAS",
  "GPS2", "GREM1", "GRIN2A", "GRM3", "GSK3B", "H3C14", "H3C2", "H3C3",
  "HDAC1", "HDAC4", "HDAC7", "HGF", "HIST1H1C", "HIST1H2BD", "HIST1H3A",
  "HIST1H3B", "HIST1H3C", "HIST1H3D", "HIST1H3E", "HIST1H3F",
  "HIST1H3G", "HIST1H3H", "HIST1H3I", "HIST1H3J", "HIST2H3C",
  "HIST2H3D", "HIST3H3", "HLA-A", "HLA-B", "HNF1A", "HOXB13", "HRAS",
  "HSD3B1", "HSP90AA1", "ICK", "ICOSLG", "ID3", "IDH1", "IDH2",
  "IFITM1", "IGF1", "IGF1R", "IGF2", "IKBKE", "IKZF1", "IL10", "IL7R",
  "INHA", "INHBA", "INPP4A", "INPP4B", "INPPL1", "INSR", "IRF2",
  "IRF4", "IRS1", "IRS2", "JAK1", "JAK2", "JAK3", "JUN", "KAT6A",
  "KDM5A", "KDM5C", "KDM6A", "KDR", "KEAP1", "KEL", "KIF5B", "KIT",
  "KLF4", "KLHL6", "KMT2A", "KMT2B", "KMT2C", "KMT2D", "KNSTRN",
  "KRAS", "LATS1", "LATS2", "LMO1", "LYN", "LZTR1", "MACC1", "MAF",
  "MALT1", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAP3K14",
  "MAPK1", "MAPK3", "MAX", "MCL1", "MDM2", "MDM4", "MED12", "MEF2B",
  "MEN1", "MERTK", "MET", "MGA", "MICA", "MICB", "MITF", "MLH1",
  "MLH3", "MLL", "MLL2", "MLL3", "MPL", "MRE11", "MRE11A", "MSH2",
  "MSH3", "MSH6", "MSI1", "MSI2", "MST1", "MST1R", "MTAP", "MTOR",
  "MUTYH", "MYC", "MYCL", "MYCL1", "MYCN", "MYD88", "MYOD1", "NBN",
  "NCOR1", "NCOR2", "NEGR1", "NF1", "NF2", "NFE2L2", "NFKBIA", "NKX2-1",
  "NKX3-1", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NPM1", "NR2F2",
  "NRAS", "NRG1", "NSD1", "NSD2", "NSD3", "NT5C2", "NTHL1", "NTRK1",
  "NTRK2", "NTRK3", "NUF2", "NUP93", "NUP98", "P2RY8", "PAK1", "PAK3",
  "PAK5", "PAK7", "PALB2", "PARK2", "PARP1", "PARP2", "PARP3", "PAX3",
  "PAX5", "PAX7", "PAX8", "PBRM1", "PDCD1", "PDCD1LG2", "PDGFRA",
  "PDGFRB", "PDPK1", "PGR", "PHF6", "PHOX2B", "PIK3C2B", "PIK3C2G",
  "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2",
  "PIK3R3", "PIM1", "PLCG2", "PLK2", "PMAIP1", "PMS1", "PMS2", "PNRC1",
  "POLD1", "POLE", "POT1", "PPARG", "PPM1D", "PPP2R1A", "PPP2R2A",
  "PPP4R2", "PPP6C", "PRDM1", "PRDM14", "PRDM16", "PREX2", "PRKAR1A",
  "PRKCI", "PRKD1", "PRKDC", "PRSS8", "PTCH1", "PTEN", "PTPN11",
  "PTPN2", "PTPN6", "PTPRD", "PTPRS", "PTPRT", "QKI", "RAB35",
  "RAC1", "RAC2", "RAD21", "RAD50", "RAD51", "RAD51B", "RAD51C",
  "RAD51D", "RAD52", "RAD54L", "RAF1", "RANBP2", "RARA", "RASA1",
  "RASGEF1A", "RB1", "RBM10", "RECQL", "RECQL4", "REL", "REST",
  "RET", "RFWD2", "RHEB", "RHOA", "RICTOR", "RIT1", "RNF43", "ROS1",
  "RPS6KA4", "RPS6KB2", "RPTOR", "RRAGC", "RRAS", "RRAS2", "RSPO2",
  "RUNX1", "RUNX1T1", "RXRA", "RYBP", "SDHA", "SDHAF2", "SDHB",
  "SDHC", "SDHD", "SETBP1", "SETD2", "SETD8", "SF3B1", "SGK1",
  "SH2B3", "SH2D1A", "SHQ1", "SLIT2", "SLX4", "SMAD2", "SMAD3",
  "SMAD4", "SMARCA4", "SMARCB1", "SMARCD1", "SMC1A", "SMC3", "SMO",
  "SMPD1", "SMYD3", "SNCAIP", "SOCS1", "SOX10", "SOX17", "SOX2",
  "SOX9", "SPEN", "SPOP", "SPRED1", "SRC", "SRSF2", "SS18", "STAG2",
  "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "STK11", "STK19",
  "STK40", "SUFU", "SUZ12", "SYK", "TAP1", "TAP2", "TBL1XR1", "TBX3",
  "TCEB1", "TCF3", "TCF7L2", "TDGF1", "TEK", "TERC", "TERT", "TET1",
  "TET2", "TGFBR1", "TGFBR2", "TIPARP", "TMEM127", "TMPRSS2", "TNFAIP3",
  "TNFRSF14", "TNK2", "TOP1", "TOP2A", "TP53", "TP53BP1", "TP63",
  "TRAF2", "TRAF3", "TRAF5", "TRAF7", "TSC1", "TSC2", "TSHR", "U2AF1",
  "U2AF2", "UPF1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VHL", "VTCN1",
  "WHSC1", "WHSC1L1", "WT1", "WWTR1", "XBP1", "XIAP", "XPA", "XPC",
  "XPO1", "XRCC2", "YAP1", "YES1", "ZBTB2", "ZFHX3", "ZMYM3", "ZNF217",
  "ZNF703"
)

# -----------------------------------------------------------------------------
# 4. Define Cell Types
# -----------------------------------------------------------------------------
cell_list <- c('pbmc', 'CD4_T', 'CD8_T', 'other_T',
               'B', 'NK', 'CD14_Mono', 'CD16_Mono',
               'DC', 'other')

# -----------------------------------------------------------------------------
# 5. Load Differential Expression Data
# -----------------------------------------------------------------------------
for (i in cell_list) {
  # R+1 vs preflight (immediately post-flight)
  x1 <- read.csv(paste0(
    "../../../Differential_expression/with_filter/celltype/pval/update/all.list/update/DEGs_",
    i, "_R+1_preflight.csv"
  ))

  # R+45&R+82 vs preflight (long-term post-flight)
  x2 <- read.csv(paste0(
    "../../../Differential_expression/with_filter/celltype/pval/update/all.list/update/DEGs_",
    i, "_R+45&R+82_preflight.csv"
  ))

  # Assign to named variables
  assign(paste0(i, '.FP1'), x1)
  assign(paste0(i, '.LP3'), x2)
}

# Uppercase PBMC for consistency
PBMC.FP1 <- pbmc.FP1
PBMC.LP3 <- pbmc.LP3

# -----------------------------------------------------------------------------
# 6. Filter Data for Cancer Genes
# -----------------------------------------------------------------------------
gene.list <- c('oncogene.tsg', 'oncokb')
cell_list <- c('PBMC', 'CD4_T', 'CD8_T', 'other_T',
               'B', 'NK', 'CD14_Mono', 'CD16_Mono',
               'DC', 'other')

for (i in cell_list) {
  x <- get(paste0(i, '.FP1'))
  y <- get(paste0(i, '.LP3'))

  for (j in gene.list) {
    A <- get(j)

    x1 <- x %>% dplyr::filter(X %in% A)
    y1 <- y %>% dplyr::filter(X %in% A)

    # Assign direction based on log2FC threshold
    x1$DIRECTION <- ifelse(x1$avg_log2FC > 0.25, 'UP',
                           ifelse(x1$avg_log2FC < -0.25, 'DOWN', 'NONE'))
    y1$DIRECTION <- ifelse(y1$avg_log2FC > 0.25, 'UP',
                           ifelse(y1$avg_log2FC < -0.25, 'DOWN', 'NONE'))

    assign(paste0(i, ".", j, ".FP1"), x1)
    assign(paste0(i, ".", j, ".LP3"), y1)
  }
}

# -----------------------------------------------------------------------------
# 7. Generate Volcano Plots
# -----------------------------------------------------------------------------
for (i in cell_list) {
  for (j in gene.list) {
    x1 <- get(paste0(i, ".", j, ".FP1"))
    y1 <- get(paste0(i, ".", j, ".LP3"))

    # Handle zero p-values (set to minimum representable value)
    x1$p_val_adj[x1$p_val_adj == 0] <- 1e-320
    y1$p_val_adj[y1$p_val_adj == 0] <- 1e-320

    options(repr.plot.width = 8, repr.plot.height = 8)

    # ----- Immediately Post-flight (R+1) -----
    g1 <- ggplot(x1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DIRECTION)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NONE" = "grey60")) +
      geom_text_repel(
        data = x1 %>% filter(DIRECTION != 'NONE'),
        aes(label = X),
        size = 4,
        max.overlaps = 20
      ) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      labs(
        title = paste0(i, " - ", j, " (R+1 vs Preflight)"),
        x = "Log2 Fold Change",
        y = "-Log10(Adjusted P-value)"
      ) +
      theme_classic(base_size = 14) +
      theme(legend.position = "bottom")

    # Save plot
    ggsave(
      filename = paste0("output/2025_01_03/volcano_", i, ".", j, ".FP1.pdf"),
      plot = g1,
      width = 8,
      height = 8
    )

    ggsave(
      filename = paste0("output/2025_01_03/volcano_", i, ".", j, ".FP1.tiff"),
      plot = g1,
      width = 8,
      height = 8,
      dpi = 300,
      compression = 'lzw'
    )

    # ----- Long-term Post-flight (R+45&R+82) -----
    g2 <- ggplot(y1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DIRECTION)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NONE" = "grey60")) +
      geom_text_repel(
        data = y1 %>% filter(DIRECTION != 'NONE'),
        aes(label = X),
        size = 4,
        max.overlaps = 20
      ) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      labs(
        title = paste0(i, " - ", j, " (R+45&R+82 vs Preflight)"),
        x = "Log2 Fold Change",
        y = "-Log10(Adjusted P-value)"
      ) +
      theme_classic(base_size = 14) +
      theme(legend.position = "bottom")

    # Save plot
    ggsave(
      filename = paste0("output/2025_01_03/volcano_", i, ".", j, ".LP3.pdf"),
      plot = g2,
      width = 8,
      height = 8
    )

    ggsave(
      filename = paste0("output/2025_01_03/volcano_", i, ".", j, ".LP3.tiff"),
      plot = g2,
      width = 8,
      height = 8,
      dpi = 300,
      compression = 'lzw'
    )
  }
}

# -----------------------------------------------------------------------------
# 8. Session Info
# -----------------------------------------------------------------------------
sessionInfo()
