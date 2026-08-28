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
input_csv <- if (length(trailing) >= 1) {
  trailing[1]
} else {
  file.path(script_dir, "core_problematic_sites.gene_summary.csv")
}
output_dir <- if (length(trailing) >= 2) trailing[2] else dirname(input_csv)

df <- read.csv(input_csv, check.names = FALSE)
df <- df[df$assembly_mutations_vs_reference > 0, ]

extract_gene_number <- function(x) {
  as.numeric(sub(".*MAB_([0-9]+).*", "\\1", x))
}

df$gene_number <- extract_gene_number(df$locus_tag)
df <- df[order(df$gene_number, df$locus_tag), ]
df$gene_rank <- seq_len(nrow(df))
df$problematic_percent <- df$problematic_rate * 100
df$band <- factor(floor((df$gene_rank - 1) / 120) %% 2)

median_problematic <- median(df$problematic_percent, na.rm = TRUE)

p <- ggplot(
  df,
  aes(
    x = gene_rank,
    y = problematic_percent,
    color = band
  )
) +
  geom_point(
    size = 2,
    alpha = 0.72,
    stroke = 0
  ) +
  geom_hline(
    yintercept = median_problematic,
    linetype = "dashed",
    linewidth = 0.55,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0, ceiling(max(df$gene_rank) / 1000) * 1000, 1000),
    expand = expansion(mult = c(0.003, 0.01))
  ) +
  scale_color_manual(
    values = c("0" = "#1F1F1F", "1" = "#8A8A8A")
  ) +
  scale_y_continuous(
    limits = c(0, max(75, ceiling(max(df$problematic_percent, na.rm = TRUE) / 10) * 10)),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = "Core genes ordered by locus tag",
    y = "Problematic assembly-derived SNPs"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10.5, color = "black"),
    axis.title = element_text(size = 12.5, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(3, "pt"),
    plot.margin = margin(6, 10, 6, 6)
  )

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "core_problematic_rate_manhattan.pdf"),
  p,
  width = 6.8,
  height = 5.1
)

ggsave(
  file.path(output_dir, "core_problematic_rate_manhattan.png"),
  p,
  width = 6.8,
  height = 5.1,
  dpi = 300
)
