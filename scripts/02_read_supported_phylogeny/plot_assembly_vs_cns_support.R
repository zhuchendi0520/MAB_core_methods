#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}

trailing <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(trailing) >= 1) trailing[1] else script_dir
output_dir <- if (length(trailing) >= 2) trailing[2] else input_dir

whole_file <- file.path(input_dir, "assembly_vs_mapping_summary.csv")
core_file <- file.path(input_dir, "core_assembly_vs_mapping_summary.csv")

whole <- read.csv(whole_file, check.names = FALSE)
core <- read.csv(core_file, check.names = FALSE)

whole$total <- whole$assembly_mutations_vs_reference
core$total <- core$core_gene_assembly_mutations_vs_reference

whole$Dataset <- "Whole genome"
core$Dataset <- "Core genes"

common_cols <- c(
  "sample",
  "Dataset",
  "total",
  "cns_supported",
  "cns_inconsistent",
  "cns_missing"
)

df <- rbind(
  whole[, common_cols],
  core[, common_cols]
)

df <- df[df$total > 0, ]

long <- rbind(
  data.frame(
    sample = df$sample,
    Dataset = df$Dataset,
    Category = "Supported",
    Proportion = df$cns_supported / df$total
  ),
  data.frame(
    sample = df$sample,
    Dataset = df$Dataset,
    Category = "Inconsistent",
    Proportion = df$cns_inconsistent / df$total
  ),
  data.frame(
    sample = df$sample,
    Dataset = df$Dataset,
    Category = "Missing / low coverage",
    Proportion = df$cns_missing / df$total
  )
)

long$Dataset <- factor(
  long$Dataset,
  levels = c("Whole genome", "Core genes")
)

long$Category <- factor(
  long$Category,
  levels = c("Supported", "Inconsistent", "Missing / low coverage")
)

category_cols <- c(
  "Supported" = "#6BAED6",
  "Inconsistent" = "#E57373",
  "Missing / low coverage" = "#BDBDBD"
)

category_labels <- c(
  "Supported" = "Supported",
  "Inconsistent" = "Inconsistent",
  "Missing / low coverage" = "Missing"
)

summary_df <- aggregate(
  Proportion ~ Dataset + Category,
  data = long,
  FUN = function(x) median(x, na.rm = TRUE)
)

summary_df$label <- sprintf("%.1f%%", summary_df$Proportion * 100)

make_plot <- function(dataset_name, title_text) {
  plot_df <- long[long$Dataset == dataset_name, ]

  ggplot(
    plot_df,
    aes(
      x = Category,
      y = Proportion,
      fill = Category,
      color = Category
    )
  ) +
    geom_violin(
      width = 0.90,
      alpha = 0.28,
      linewidth = 0,
      trim = FALSE,
      scale = "width",
      adjust = 1.10
    ) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      alpha = 0.92,
      linewidth = 0.95,
      color = "black"
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 3.9,
      stroke = 1.0,
      fill = "white",
      color = "black"
    ) +
    scale_fill_manual(values = category_cols) +
    scale_color_manual(values = category_cols) +
    scale_x_discrete(labels = category_labels) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = function(x) paste0(x * 100, "%"),
      expand = c(0.01, 0)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = title_text,
      x = NULL,
      y = "Proportion (%)"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 22,
        face = "bold",
        color = "black",
        hjust = 0.5,
        margin = margin(0, 0, 10, 0)
      ),
      axis.text.x = element_text(
        size = 16,
        color = "black",
        lineheight = 1
      ),
      axis.text.y = element_text(
        size = 15,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.9,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.75,
        color = "black"
      ),
      axis.ticks.length = unit(4, "pt"),
      plot.margin = margin(12, 16, 12, 14)
    )
}

p_whole <- make_plot(
  "Whole genome",
  "Whole-genome assembly SNPs"
)

p_core <- make_plot(
  "Core genes",
  "Core-gene assembly SNPs"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "assembly_vs_cns_support_whole_genome_boxplot.pdf"),
  p_whole,
  width = 6.4,
  height = 4.8
)

ggsave(
  file.path(output_dir, "assembly_vs_cns_support_whole_genome_boxplot.png"),
  p_whole,
  width = 6.4,
  height = 4.8,
  dpi = 300
)

ggsave(
  file.path(output_dir, "assembly_vs_cns_support_core_genes_boxplot.pdf"),
  p_core,
  width = 6.4,
  height = 4.8
)

ggsave(
  file.path(output_dir, "assembly_vs_cns_support_core_genes_boxplot.png"),
  p_core,
  width = 6.4,
  height = 4.8,
  dpi = 300
)

write.csv(
  summary_df,
  file.path(output_dir, "assembly_vs_cns_support_median_summary.csv"),
  row.names = FALSE
)
