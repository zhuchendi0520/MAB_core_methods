#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 01_plot_treemmer_sampling_schematic.R OUTPUT_DIR")
}
outdir <- args[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cluster_sizes <- c(8, 6, 9, 5, 7, 6, 8, 5, 7, 6, 8, 5)
cluster_subsp <- c("ABS", "ABS", "ABS", "MAS", "MAS", "ABS", "BOL", "MAS", "ABS", "MAS", "BOL", "ABS")
subsp_cols <- c("ABS" = "#F28B63", "MAS" = "#256DB3", "BOL" = "#7C2C8E")
branch_col <- "#66BD73"

segments_full <- tibble()
tips_full <- tibble()
cluster_info <- tibble()

add_seg <- function(tbl, panel, x, y, xend, yend, cluster) {
  bind_rows(tbl, tibble(panel = panel, x = x, y = y, xend = xend, yend = yend, cluster = cluster))
}

y0 <- 0
for (i in seq_along(cluster_sizes)) {
  size <- cluster_sizes[[i]]
  ys <- y0:(y0 + size - 1)
  root_y <- mean(ys)
  stem_x <- 0.28 + 0.045 * ((i - 1) %% 4)
  tip_x <- 1.00
  segments_full <- add_seg(segments_full, "Full tree", stem_x, min(ys), stem_x, max(ys), i)
  segments_full <- add_seg(segments_full, "Full tree", 0.14, root_y, stem_x, root_y, i)
  close_cluster <- i %in% c(2, 3, 6, 9, 11)
  keep_idx <- c(ceiling(size / 2))
  if (size >= 8) keep_idx <- c(keep_idx, 2)
  for (j in seq_along(ys)) {
    branch_x <- 0.62 + 0.035 * ((j + i) %% 3)
    segments_full <- add_seg(segments_full, "Full tree", stem_x, ys[[j]], branch_x, ys[[j]], i)
    segments_full <- add_seg(segments_full, "Full tree", branch_x, ys[[j]], tip_x, ys[[j]], i)
    tips_full <- bind_rows(
      tips_full,
      tibble(
        panel = "Full tree",
        x = tip_x,
        y = ys[[j]],
        cluster = i,
        subsp = cluster_subsp[[i]],
        selected = j %in% keep_idx,
        close_cluster = close_cluster
      )
    )
  }
  cluster_info <- bind_rows(
    cluster_info,
    tibble(cluster = i, y_min = min(ys) - 0.45, y_max = max(ys) + 0.45, y_mid = root_y, close_cluster = close_cluster)
  )
  y0 <- y0 + size + 1
}

root_y <- mean(cluster_info$y_mid)
segments_full <- add_seg(segments_full, "Full tree", 0.04, min(cluster_info$y_mid), 0.04, max(cluster_info$y_mid), 0)
segments_full <- add_seg(segments_full, "Full tree", 0.00, root_y, 0.04, root_y, 0)
for (yy in cluster_info$y_mid) {
  segments_full <- add_seg(segments_full, "Full tree", 0.04, yy, 0.14, yy, 0)
}

selected_full <- tips_full %>% filter(selected)
rep_y <- selected_full %>%
  arrange(y) %>%
  mutate(rep_y = (row_number() - 1) / (n() - 1) * (max(tips_full$y) * 0.86) + max(tips_full$y) * 0.06)

rep_segments <- tibble()
rep_tips <- tibble()
rep_clusters <- rep_y %>%
  group_by(cluster, subsp) %>%
  summarise(y_mid = mean(rep_y), y_min = min(rep_y), y_max = max(rep_y), .groups = "drop")

for (i in seq_len(nrow(rep_clusters))) {
  row <- rep_clusters[i, ]
  stem_x <- 1.70 + 0.04 * ((row$cluster - 1) %% 3)
  rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", stem_x, row$y_min, stem_x, row$y_max, row$cluster)
  rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", 1.50, row$y_mid, stem_x, row$y_mid, row$cluster)
}
for (i in seq_len(nrow(rep_y))) {
  row <- rep_y[i, ]
  branch_x <- 1.98 + 0.03 * (i %% 2)
  rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", 1.76, row$rep_y, branch_x, row$rep_y, row$cluster)
  rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", branch_x, row$rep_y, 2.26, row$rep_y, row$cluster)
  rep_tips <- bind_rows(rep_tips, tibble(panel = "Treemmer-pruned tree", x = 2.26, y = row$rep_y, subsp = row$subsp))
}
rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", 1.42, min(rep_clusters$y_mid), 1.42, max(rep_clusters$y_mid), 0)
rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", 1.36, mean(rep_clusters$y_mid), 1.42, mean(rep_clusters$y_mid), 0)
for (yy in rep_clusters$y_mid) {
  rep_segments <- add_seg(rep_segments, "Treemmer-pruned tree", 1.42, yy, 1.50, yy, 0)
}

highlight_df <- cluster_info %>%
  filter(close_cluster) %>%
  left_join(tibble(cluster = seq_along(cluster_subsp), subsp = cluster_subsp), by = "cluster") %>%
  mutate(
    xmin = 0.54,
    xmax = 1.07,
    fill_col = recode(subsp, "ABS" = "#F28B63", "MAS" = "#40C9D0", "BOL" = "#C77CFF")
  )

legend_df <- tibble(subsp = names(subsp_cols), x = c(0.10, 0.40, 0.70), y = -5.2)

p <- ggplot() +
  geom_rect(
    data = highlight_df,
    aes(xmin = xmin, xmax = xmax, ymin = y_min, ymax = y_max),
    fill = highlight_df$fill_col,
    color = highlight_df$fill_col,
    linewidth = 0.28,
    alpha = 0.22
  ) +
  geom_segment(
    data = segments_full,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = branch_col,
    linewidth = 0.75,
    lineend = "round"
  ) +
  geom_point(
    data = tips_full,
    aes(x = x, y = y, fill = subsp),
    shape = 21,
    color = "white",
    stroke = 0.18,
    size = 2.05
  ) +
  geom_point(
    data = selected_full,
    aes(x = x, y = y, fill = subsp),
    shape = 21,
    color = "white",
    stroke = 0.25,
    size = 2.7
  ) +
  geom_segment(
    data = rep_segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = branch_col,
    linewidth = 0.85,
    lineend = "round"
  ) +
  geom_point(
    data = rep_tips,
    aes(x = x, y = y, fill = subsp),
    shape = 21,
    color = "white",
    stroke = 0.25,
    size = 3.0
  ) +
  annotate("text", x = 0.02, y = max(tips_full$y) + 5.3, label = "Full tree", hjust = 0, size = 3.0, fontface = "bold") +
  annotate("text", x = 0.02, y = max(tips_full$y) + 2.8, label = "6,733 genomes", hjust = 0, size = 2.35, color = "#555555") +
  annotate("text", x = 0.58, y = max(tips_full$y) + 2.8, label = "closely related\nredundant tips", hjust = 0, size = 2.2, color = "#A33A34", lineheight = 0.92) +
  annotate("text", x = 1.35, y = max(tips_full$y) + 5.3, label = "After Treemmer pruning", hjust = 0, size = 3.0, fontface = "bold") +
  annotate("text", x = 1.35, y = max(tips_full$y) + 2.8, label = "1,130 representatives; 90% diversity retained", hjust = 0, size = 2.35, color = "#555555") +
  annotate(
    "segment",
    x = 1.10, xend = 1.29, y = max(tips_full$y) * 0.55, yend = max(tips_full$y) * 0.55,
    arrow = arrow(length = unit(0.13, "inches"), type = "closed"),
    color = "#333333",
    linewidth = 0.55
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, fill = subsp),
    shape = 21,
    color = "white",
    stroke = 0.25,
    size = 2.3
  ) +
  geom_text(
    data = legend_df,
    aes(x = x + 0.04, y = y, label = subsp),
    hjust = 0,
    vjust = 0.5,
    size = 2.1,
    color = "#333333"
  ) +
  scale_fill_manual(values = subsp_cols) +
  coord_cartesian(xlim = c(-0.02, 2.42), ylim = c(-7.0, max(tips_full$y) + 8), clip = "off") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  )

png <- file.path(outdir, "Fig1A_treemmer_two_panel_schematic.png")
pdf <- file.path(outdir, "Fig1A_treemmer_two_panel_schematic.pdf")
svg <- file.path(outdir, "Fig1A_treemmer_two_panel_schematic.svg")

ggsave(png, p, width = 7.4, height = 2.8, dpi = 600, bg = "white")
ggsave(pdf, p, width = 7.4, height = 2.8, useDingbats = FALSE, bg = "white")
ggsave(svg, p, width = 7.4, height = 2.8, bg = "white")

message(png)
message(pdf)
message(svg)
