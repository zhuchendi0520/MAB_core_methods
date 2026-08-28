library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_dcc_pns_switch_convergence_heatmap_fullpost_horizontal.R INPUT_CSV OUTPUT_DIR")
}
infile <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dcc_order <- c("DCC1", "DCC2", "DCC4", "DCC5", "DCC3", "DCC6", "DCC7")

df <- read.csv(infile, check.names = FALSE)
df$gene <- as.character(df$gene)
df$Preferred_name <- as.character(df$Preferred_name)
df$display_gene <- ifelse(
  is.na(df$Preferred_name) | df$Preferred_name == "" | df$Preferred_name == "nan",
  df$gene,
  df$Preferred_name
)
# Correct the displayed symbol for MAB_1915 without altering its gene ID.
df$display_gene[df$gene == "MAB_1915"] <- "fadD"
df$display_gene <- make.unique(df$display_gene)

keep_genes <- names(which(tapply(df$switch_call != "unchanged", df$gene, any)))
plot_df <- df[df$gene %in% keep_genes, ]
plot_df$DCC <- factor(plot_df$DCC, levels = rev(dcc_order))
plot_df$pNS_pre <- as.numeric(plot_df$pNS_pre)
plot_df$pNS_post_full <- as.numeric(plot_df$pNS_post_full)
# Manual review: ubiA in DCC6 is classified as purifying-to-positive.
plot_df$switch_call[plot_df$gene == "MAB_0173" & plot_df$DCC == "DCC6"] <- "negative->positive"

# Present genes by the direction observed in the pooled analysis: genes with an
# overall purifying-to-positive switch first, followed by positive-to-purifying.
# Within each block, prioritize genes showing that direction in more DCCs.
gene_summary <- aggregate(
  cbind(
    orange_n = plot_df$switch_call == "negative->positive",
    blue_n = plot_df$switch_call == "positive->negative"
  ),
  by = list(gene = plot_df$gene),
  FUN = sum
)
direction_map <- plot_df[!duplicated(plot_df$gene), c("gene", "total_switch_direction")]
gene_summary <- merge(gene_summary, direction_map, by = "gene", all.x = TRUE, sort = FALSE)
gene_summary$direction_rank <- ifelse(
  gene_summary$total_switch_direction == "purifying_to_positive", 1, 2
)
gene_summary$primary_n <- ifelse(
  gene_summary$direction_rank == 1, gene_summary$orange_n, gene_summary$blue_n
)
gene_summary$secondary_n <- ifelse(
  gene_summary$direction_rank == 1, gene_summary$blue_n, gene_summary$orange_n
)
gene_summary <- gene_summary[
  order(
    gene_summary$direction_rank,
    -gene_summary$primary_n,
    gene_summary$secondary_n,
    gene_summary$gene
  ),
]
gene_order <- gene_summary$gene
gene_label_map <- plot_df[match(gene_order, plot_df$gene), c("gene", "display_gene")]

plot_df$gene <- factor(plot_df$gene, levels = gene_order)
plot_df$display_gene_ordered <- factor(
  plot_df$gene,
  levels = gene_order,
  labels = gene_label_map$display_gene
)
plot_df$state_raw <- ifelse(
  plot_df$switch_call == "negative->positive",
  "Purifying \u2192 Positive selection",
  ifelse(
    plot_df$switch_call == "positive->negative",
    "Positive \u2192 Purifying selection",
    "Unchanged"
  )
)
plot_df$state <- factor(
  plot_df$state_raw,
  levels = c(
    "Purifying \u2192 Positive selection",
    "Positive \u2192 Purifying selection",
    "Unchanged"
  )
)

state_cols <- c(
  "Purifying \u2192 Positive selection" = "#E59673",
  "Positive \u2192 Purifying selection" = "#8FA7DB",
  "Unchanged" = "#FFFFFF"
)

p <- ggplot(plot_df, aes(x = display_gene_ordered, y = DCC, fill = state)) +
  geom_tile(color = "#D8D8D8", linewidth = 0.28, width = 0.95, height = 0.95) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  scale_x_discrete(position = "top") +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(color = "black", size = 8.5, angle = 90, hjust = 0, vjust = 0.5, face = "italic"),
    axis.text.y = element_text(color = "black", size = 8.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 8.5, color = "black"),
    panel.border = element_blank(),
    plot.margin = margin(8, 8, 5, 6)
  ) +
  guides(fill = guide_legend(override.aes = list(color = "#D8D8D8", linewidth = 0.28), nrow = 1))

width <- max(8.5, 0.18 * length(unique(plot_df$gene)) + 1.1)
ggsave(file.path(out_dir, "DCC_pNS_switch_convergence_heatmap_fullpost_simple_horizontal.png"), p, width = width, height = 3.2, dpi = 600)
ggsave(file.path(out_dir, "DCC_pNS_switch_convergence_heatmap_fullpost_simple_horizontal.pdf"), p, width = width, height = 3.2)
