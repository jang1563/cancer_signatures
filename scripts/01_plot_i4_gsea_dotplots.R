#!/usr/bin/env Rscript

# Regenerate manuscript-style Inspiration4 fGSEA dot plots from processed tables.

required_packages <- c("ggplot2", "dplyr", "forcats", "scales", "grid")
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
  library(forcats)
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

cell_levels <- c(
  "PBMC", "CD4 T", "CD8 T", "other T", "B",
  "NK", "CD14 Mono", "CD16 Mono", "DC", "other"
)
time_levels <- c("Immediately Post-flight", "Long-term Post-flight")

read_fgsea <- function(filename, collection) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) stop("Missing processed table: ", path, call. = FALSE)

  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("pathway", "padj", "NES", "timepoint", "celltype")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("Missing columns in ", filename, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }

  df$info <- collection
  df$padj <- as.numeric(df$padj)
  df$NES <- as.numeric(df$NES)
  df$timepoint <- factor(df$timepoint, levels = time_levels)
  df$celltype <- factor(df$celltype, levels = cell_levels)
  df$neg_log10_padj <- pmin(-log10(pmax(df$padj, .Machine$double.xmin)), 50)
  df
}

clean_pathway_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x
}

drop_large_gene_column <- function(df) {
  df[, setdiff(names(df), "leadingEdge"), drop = FALSE]
}

write_table <- function(df, filename) {
  write.csv(drop_large_gene_column(as.data.frame(df)), file.path(table_dir, filename), row.names = FALSE)
}

plot_dot <- function(plot_df, title, subtitle, output_stub, width = 12.5, height = 10.5) {
  if (!nrow(plot_df)) stop("No rows available for ", output_stub, call. = FALSE)

  nes_limit <- ceiling(max(abs(plot_df$NES), na.rm = TRUE) * 10) / 10
  p <- ggplot(
    plot_df,
    aes(x = celltype, y = pathway, size = neg_log10_padj, fill = NES)
  ) +
    geom_point(shape = 21, alpha = 0.92, color = "grey15", stroke = 0.28) +
    facet_grid(. ~ timepoint) +
    scale_y_discrete(labels = clean_pathway_label) +
    scale_fill_gradient2(
      low = "#2C5AA0",
      mid = "#F8F8F8",
      high = "#B8322A",
      midpoint = 0,
      limits = c(-nes_limit, nes_limit),
      oob = scales::squish,
      name = "NES"
    ) +
    scale_size_continuous(
      range = c(2.4, 8.2),
      breaks = c(10, 20, 30, 40),
      name = "-log10(adj. p)"
    ) +
    labs(title = title, subtitle = subtitle, x = "Cell type", y = "Pathway") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey25", margin = margin(b = 7)),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5),
      axis.text.y = element_text(size = 7.4),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25, linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(8, "mm"),
      strip.text = element_text(size = 10.5, face = "bold", color = "white"),
      strip.background = element_rect(fill = "#243B6B", color = NA),
      legend.title = element_text(size = 9.5, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.box = "vertical",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 10, 8, 8)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(36, "mm")),
      size = guide_legend(order = 2, override.aes = list(fill = "white"))
    )

  ggsave(file.path(fig_dir, paste0(output_stub, ".pdf")), p, width = width, height = height, useDingbats = FALSE, bg = "white")
  ggsave(file.path(fig_dir, paste0(output_stub, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  invisible(p)
}

build_top_plot_df <- function(df, padj_cutoff, positive_n, negative_n) {
  sig <- df %>% filter(padj < padj_cutoff)
  best_per_pathway <- sig %>%
    group_by(pathway) %>%
    arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(direction = if_else(NES >= 0, "Positive NES", "Negative NES"))

  selected <- bind_rows(
    best_per_pathway %>%
      filter(direction == "Positive NES") %>%
      arrange(padj, desc(abs(NES)), pathway) %>%
      slice_head(n = positive_n),
    best_per_pathway %>%
      filter(direction == "Negative NES") %>%
      arrange(padj, desc(abs(NES)), pathway) %>%
      slice_head(n = negative_n)
  ) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(desc(NES), padj, pathway) %>%
    mutate(display_order = row_number())

  plot_df <- sig %>%
    filter(pathway %in% selected$pathway) %>%
    mutate(pathway = factor(pathway, levels = rev(selected$pathway)))

  list(plot_df = plot_df, selected = selected, significant = sig)
}

hallmark <- read_fgsea("i4_fgsea_hallmark.csv", "Hallmark")
c4 <- read_fgsea("i4_fgsea_c4.csv", "C4")
c6 <- read_fgsea("i4_fgsea_c6.csv", "C6")

hallmark_sig <- hallmark %>%
  filter(padj < 0.05) %>%
  mutate(pathway = fct_reorder(pathway, NES))
write_table(hallmark_sig, "fig4c_hallmark_significant_padj_lt_0.05.csv")
plot_dot(
  hallmark_sig,
  "Hallmark pathway activity after spaceflight",
  "fGSEA results for Inspiration4 immune cell subsets; points show adjusted p < 0.05",
  "fig4c_i4_hallmark_dotplot",
  width = 12.5,
  height = 7.5
)

c6_top <- build_top_plot_df(c6, padj_cutoff = 0.05, positive_n = 20, negative_n = 20)
write_table(c6_top$significant, "fig4a_c6_significant_padj_lt_0.05.csv")
write_table(c6_top$selected, "fig4a_c6_selected_pathways.csv")
plot_dot(
  c6_top$plot_df,
  "C6 oncogenic signatures after spaceflight",
  "Top positive and negative C6 pathways selected by minimum adjusted p-value",
  "fig4a_i4_c6_top_pathways_dotplot",
  width = 12.5,
  height = 10.5
)

c4_top <- build_top_plot_df(c4, padj_cutoff = 0.001, positive_n = 20, negative_n = 20)
write_table(c4_top$significant, "fig4b_c4_significant_padj_lt_0.001.csv")
write_table(c4_top$selected, "fig4b_c4_selected_pathways.csv")
plot_dot(
  c4_top$plot_df,
  "C4 immune pathway shifts after spaceflight",
  "Top positive and negative C4 pathways selected by minimum adjusted p-value; points show adjusted p < 0.001",
  "fig4b_i4_c4_top_pathways_dotplot",
  width = 12.5,
  height = 10.5
)

message("Figures written to: ", fig_dir)
message("Source tables written to: ", table_dir)
