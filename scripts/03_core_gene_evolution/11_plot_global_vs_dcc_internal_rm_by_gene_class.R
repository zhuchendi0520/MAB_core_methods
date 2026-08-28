#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript 11_plot_global_vs_dcc_internal_rm_by_gene_class.R ",
    "<global_rm.csv> <dcc_internal_summary.csv> <output_directory>"
  )
}

global_file <- args[[1]]
dcc_file <- args[[2]]
out_dir <- args[[3]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

class_order <- c("ES", "GA", "GD", "NE")
class_labels <- c(
  ES = "Essential",
  GA = "Growth\nadvantage",
  GD = "Growth\ndefect",
  NE = "Non-essential"
)

status_order <- c("Global", "DCC internal")
status_cols <- c(
  "Global" = "#4FA3D1",
  "DCC internal" = "#C77DBB"
)

global <- read_csv(global_file, show_col_types = FALSE) %>%
  transmute(
    gene = trimws(gene),
    HMM_class = trimws(as.character(HMM_class)),
    Global = rm_m0_as1
  ) %>%
  filter(HMM_class %in% class_order)

dcc <- read_csv(dcc_file, show_col_types = FALSE) %>%
  transmute(
    gene = trimws(gene),
    rSNP_post = as.numeric(rSNP_post),
    mSNP_post = as.numeric(mSNP_post),
    `DCC internal` = rSNP_post / if_else(mSNP_post == 0 | is.na(mSNP_post), 1, mSNP_post)
  )

paired <- global %>%
  inner_join(dcc, by = "gene") %>%
  filter(!is.na(Global), !is.na(`DCC internal`)) %>%
  mutate(
    HMM_class = factor(HMM_class, levels = class_order)
  )

plot_df <- paired %>%
  pivot_longer(
    cols = all_of(status_order),
    names_to = "status",
    values_to = "rm"
  ) %>%
  mutate(
    status = factor(status, levels = status_order),
    HMM_class = factor(HMM_class, levels = class_order),
    x = as.numeric(HMM_class),
    x_density = x - 0.20,
    x_box = if_else(status == "Global", x + 0.06, x + 0.18)
  )

summary_df <- plot_df %>%
  group_by(HMM_class, status) %>%
  summarise(
    n_genes = n(),
    median_rm = median(rm, na.rm = TRUE),
    IQR_rm = IQR(rm, na.rm = TRUE),
    .groups = "drop"
  )

stats_df <- paired %>%
  group_by(HMM_class) %>%
  summarise(
    n_genes = n(),
    median_global = median(Global, na.rm = TRUE),
    median_DCC_internal = median(`DCC internal`, na.rm = TRUE),
    median_difference_DCC_minus_global = median(`DCC internal` - Global, na.rm = TRUE),
    wilcoxon_V = unname(wilcox.test(`DCC internal`, Global, paired = TRUE, exact = FALSE)$statistic),
    p_value = wilcox.test(`DCC internal`, Global, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj_BH < 0.0001 ~ "****",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x = as.numeric(factor(HMM_class, levels = class_order)),
    xstart = x + 0.03,
    xend = x + 0.21,
    x_label = (xstart + xend) / 2,
    y = c(31, 31, 31, 31)
  )

write_csv(paired, file.path(out_dir, "core_gene_global_vs_DCC_internal_rm_by_HMM_paired.csv"))
write_csv(summary_df, file.path(out_dir, "core_gene_global_vs_DCC_internal_rm_by_HMM_summary.csv"))
write_csv(stats_df, file.path(out_dir, "core_gene_global_vs_DCC_internal_rm_by_HMM_paired_wilcoxon.csv"))

rm_trans <- scales::pseudo_log_trans(sigma = 0.03, base = 10)
y_break_values <- c(0, 0.1, 1, 10)
y_breaks <- rm_trans$transform(y_break_values)

plot_df <- plot_df %>%
  mutate(y_plot = rm_trans$transform(rm))

stats_df <- stats_df %>%
  mutate(y_plot = rm_trans$transform(y))

make_half_violin <- function(df, max_width = 0.23) {
  out <- list()
  split_df <- split(df, list(df$HMM_class, df$status), drop = TRUE)
  for (nm in names(split_df)) {
    d <- split_df[[nm]]
    vals <- d$y_plot[is.finite(d$y_plot)]
    if (length(unique(vals)) < 2) next
    den <- density(vals, n = 256, adjust = 1.05, from = min(vals), to = max(vals))
    x_base <- unique(d$x_density)[1]
    direction <- -1
    widths <- den$y / max(den$y) * max_width
    edge <- tibble(
      HMM_class = unique(d$HMM_class)[1],
      status = unique(d$status)[1],
      x = x_base + direction * widths,
      y_plot = den$x,
      part = "edge"
    )
    base <- tibble(
      HMM_class = unique(d$HMM_class)[1],
      status = unique(d$status)[1],
      x = x_base,
      y_plot = rev(den$x),
      part = "base"
    )
    out[[nm]] <- bind_rows(edge, base) %>%
      mutate(poly_id = nm)
  }
  bind_rows(out)
}

violin_df <- make_half_violin(plot_df)

p <- ggplot(plot_df, aes(y = y_plot, fill = status, color = status)) +
  geom_polygon(
    data = violin_df,
    aes(x = x, y = y_plot, group = poly_id, fill = status),
    inherit.aes = FALSE,
    alpha = 0.50,
    color = "#2C2C2C",
    linewidth = 0.34,
    show.legend = TRUE
  ) +
  geom_point(
    aes(x = x_box, y = y_plot),
    position = position_jitter(width = 0.026, height = 0, seed = 20260811),
    size = 0.68,
    alpha = 0.90,
    stroke = 0,
    show.legend = FALSE
  ) +
  geom_boxplot(
    aes(x = x_box, y = y_plot, group = interaction(HMM_class, status), fill = status, color = status),
    width = 0.085,
    outlier.shape = NA,
    linewidth = 0.36,
    alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = rm_trans$transform(1),
    linetype = "dashed",
    linewidth = 0.42,
    color = "#666666"
  ) +
  geom_segment(
    data = stats_df,
    aes(x = xstart, xend = xend, y = y_plot, yend = y_plot),
    inherit.aes = FALSE,
    linewidth = 0.38,
    color = "black"
  ) +
  geom_text(
    data = stats_df,
    aes(x = x_label, y = y_plot, label = significance),
    inherit.aes = FALSE,
    vjust = -0.30,
    size = 6.0,
    fontface = "bold",
    color = "black"
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = c("0", "0.1", "1", "10")
  ) +
  scale_x_continuous(
    breaks = seq_along(class_order),
    labels = unname(class_labels[class_order]),
    limits = c(0.55, 4.55)
  ) +
  scale_fill_manual(values = status_cols, breaks = status_order, name = "r/m source") +
  scale_color_manual(values = status_cols, breaks = status_order, name = "r/m source") +
  coord_cartesian(ylim = rm_trans$transform(c(0, 44)), clip = "off") +
  labs(
    x = NULL,
    y = "Recombination-to-mutation ratio"
  ) +
  theme_classic(base_size = 19, base_family = "Helvetica") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 13, color = "black", face = "plain"),
    legend.text = element_text(size = 13, color = "black"),
    legend.key.size = unit(0.27, "cm"),
    legend.spacing.x = unit(0.10, "cm"),
    axis.title.y = element_text(face = "bold", size = 18, margin = margin(r = 8)),
    axis.text.x = element_text(face = "bold", size = 18, color = "black", lineheight = 0.88),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 0.95, color = "black"),
    axis.ticks = element_line(linewidth = 0.82, color = "black"),
    axis.ticks.length = unit(0.12, "cm"),
    plot.margin = margin(6, 12, 8, 10)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 0.80, color = "#2C2C2C", linewidth = 0.30),
      nrow = 1,
      title.position = "left"
    ),
    color = "none"
  )

ggsave(
  file.path(out_dir, "core_gene_global_vs_DCC_internal_rm_HMM_split_violin.png"),
  p,
  width = 8.45,
  height = 5.0,
  dpi = 450,
  bg = "white"
)
ggsave(
  file.path(out_dir, "core_gene_global_vs_DCC_internal_rm_HMM_split_violin.pdf"),
  p,
  width = 8.45,
  height = 5.0,
  useDingbats = FALSE,
  bg = "white"
)

message("Saved outputs to: ", out_dir)
message("Genes compared: ", nrow(paired))
