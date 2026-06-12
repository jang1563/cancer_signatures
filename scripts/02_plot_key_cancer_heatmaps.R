#!/usr/bin/env Rscript

# Regenerate manuscript-style key cancer gene heatmaps from processed tables.

required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "circlize", "ComplexHeatmap", "grid")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Install required R packages before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE))
  }
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (basename(wd) == "scripts") {
    return(normalizePath(file.path(wd, ".."), mustWork = TRUE))
  }
  wd
}

root_dir <- repo_root()
data_dir <- file.path(root_dir, "data", "processed")
fig_dir <- file.path(root_dir, "outputs", "figures")
table_dir <- file.path(root_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(data_dir, "i4_key_cancer_genes_log2fc.csv")
if (!file.exists(input_path)) stop("Missing processed table: ", input_path, call. = FALSE)

df <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "Gene", "Cell_Type",
  "log2FC_R1_vs_Preflight", "log2FC_R45_R82_vs_Preflight",
  "In_Key_Cancer"
)
missing <- setdiff(required, names(df))
if (length(missing)) {
  stop("Missing columns in key cancer table: ", paste(missing, collapse = ", "), call. = FALSE)
}

cell_keep <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")
cell_labels <- c(
  CD4_T = "CD4 T",
  CD8_T = "CD8 T",
  other_T = "other T",
  B = "B",
  NK = "NK",
  CD14_Mono = "CD14 Mono",
  CD16_Mono = "CD16 Mono",
  DC = "DC",
  other = "other"
)

plot_df <- df %>%
  filter(Cell_Type %in% cell_keep, In_Key_Cancer %in% c(TRUE, "TRUE", "True", "true", 1, "1")) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = cell_keep, labels = unname(cell_labels[cell_keep]))) %>%
  pivot_longer(
    cols = c("log2FC_R1_vs_Preflight", "log2FC_R45_R82_vs_Preflight"),
    names_to = "comparison",
    values_to = "log2FC"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      log2FC_R1_vs_Preflight = "R+1 vs Pre-flight",
      log2FC_R45_R82_vs_Preflight = "R+45/R+82 vs Pre-flight"
    ),
    comparison = factor(comparison, levels = c("R+1 vs Pre-flight", "R+45/R+82 vs Pre-flight")),
    Gene = factor(Gene, levels = rev(unique(Gene)))
  )

write.csv(
  as.data.frame(plot_df),
  file.path(table_dir, "fig3b_key_cancer_log2fc_long_table.csv"),
  row.names = FALSE
)

fig3b <- ggplot(plot_df, aes(x = Cell_Type, y = Gene, fill = log2FC)) +
  geom_tile(color = "grey82", linewidth = 0.25) +
  facet_grid(. ~ comparison) +
  scale_fill_gradient2(
    low = "#1F3B8B",
    mid = "white",
    high = "#B8232A",
    midpoint = 0,
    limits = c(-0.35, 0.35),
    oob = scales::squish,
    name = "log2FC"
  ) +
  labs(x = "Cell type", y = "Gene") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(face = "italic", size = 9),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10, color = "white"),
    strip.background = element_rect(fill = "#243B6B", color = NA),
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(fig_dir, "fig3b_i4_key_cancer_log2fc_heatmap.pdf"),
  fig3b,
  width = 10.5,
  height = 6.5,
  useDingbats = FALSE,
  bg = "white"
)
ggsave(
  file.path(fig_dir, "fig3b_i4_key_cancer_log2fc_heatmap.png"),
  fig3b,
  width = 10.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

wide <- df %>%
  filter(Cell_Type %in% cell_keep, In_Key_Cancer %in% c(TRUE, "TRUE", "True", "true", 1, "1")) %>%
  select(Gene, Cell_Type, log2FC_R45_R82_vs_Preflight) %>%
  pivot_wider(names_from = Cell_Type, values_from = log2FC_R45_R82_vs_Preflight)

mat <- as.matrix(wide[, cell_keep])
rownames(mat) <- wide$Gene
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

write.csv(
  mat,
  file.path(table_dir, "fig3c_source_matrix_log2FC_R45R82_vs_Pre.csv")
)

col_fun <- circlize::colorRamp2(
  c(-0.1, 0, 0.3),
  c("#1F3B8B", "white", "#B8232A")
)

ht <- Heatmap(
  mat,
  name = "Log2 Fold\nChange",
  col = col_fun,
  border = TRUE,
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_labels = unname(cell_labels[colnames(mat)]),
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_title = "Key Cancer Genes\nR+45/R+82 vs Pre-flight (log2FC)",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Log2 Fold\nChange",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    at = c(-0.1, 0, 0.1, 0.2, 0.3),
    legend_height = unit(4, "cm")
  ),
  use_raster = FALSE
)

draw_heatmap <- function() {
  grid.newpage()
  draw(
    ht,
    newpage = FALSE,
    heatmap_legend_side = "right",
    merge_legends = TRUE,
    padding = unit(c(12, 24, 22, 20), "mm")
  )
}

pdf(
  file.path(fig_dir, "fig3c_i4_key_cancer_log2fc_r45_r82_heatmap.pdf"),
  width = 8,
  height = 9.5,
  useDingbats = FALSE
)
draw_heatmap()
dev.off()

png(
  file.path(fig_dir, "fig3c_i4_key_cancer_log2fc_r45_r82_heatmap.png"),
  width = 8,
  height = 9.5,
  units = "in",
  res = 300,
  bg = "white",
  type = "cairo"
)
draw_heatmap()
dev.off()

message("Figures written to: ", fig_dir)
message("Source tables written to: ", table_dir)
