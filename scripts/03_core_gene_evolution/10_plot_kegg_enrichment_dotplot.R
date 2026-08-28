library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(
  length(args) >= 1,
  args[1],
  "kegg_enrichment_plot_data.csv"
)

out_dir <- ifelse(length(args) >= 2, args[2], dirname(input_file))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, check.names = FALSE)

df <- df %>%
  mutate(
    direction = factor(
      direction,
      levels = c("increased", "decreased"),
      labels = c("r/m increased", "r/m decreased")
    ),
    term_type = factor(term_type, levels = c("Pathway", "Module")),
    label = str_replace_all(label, " - ", "-"),
    label = forcats::fct_reorder(label, neg_log10_fdr, .desc = FALSE),
    gene_ratio_percent = gene_ratio * 100,
    fdr_for_plot = pmax(fdr_BH, 1e-6),
    significance = if_else(fdr_BH <= 0.05, "FDR <= 0.05", "FDR > 0.05")
  )

plot_df <- df %>%
  group_by(direction) %>%
  arrange(fdr_BH, p_value, .by_group = TRUE) %>%
  mutate(label = factor(label, levels = rev(unique(label)))) %>%
  ungroup()

global_max_x <- max(plot_df$gene_ratio_percent, na.rm = TRUE)

make_plot <- function(sub_df, plot_title, x_max) {
  sub_df <- sub_df %>%
    arrange(fdr_BH, p_value) %>%
    mutate(label = factor(label, levels = rev(unique(label))))

  ggplot(
    sub_df,
    aes(
      x = gene_ratio_percent,
      y = label
    )
  ) +
    geom_point(
      aes(
        size = target_genes_in_term,
        fill = -log10(fdr_for_plot),
        shape = term_type
      ),
      color = "black",
      stroke = 0.28,
      alpha = 0.92
    ) +
    scale_x_continuous(
      limits = c(0, x_max * 1.16),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_shape_manual(
      values = c("Pathway" = 21, "Module" = 22),
      name = "KEGG class"
    ) +
    scale_fill_gradient(
      low = "#F7E6A3",
      high = "#D95F59",
      name = expression(-log[10]("FDR"))
    ) +
    scale_size_continuous(
      range = c(2.2, 6.4),
      name = "Gene count"
    ) +
    labs(
      title = plot_title,
      x = "Gene ratio (%)",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 16.6, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 12.3, color = "black"),
      axis.text.y = element_text(size = 11.3, color = "black"),
      axis.title.x = element_text(size = 13.4, face = "bold"),
      axis.line = element_line(linewidth = 0.52, color = "black"),
      axis.ticks = element_line(linewidth = 0.42, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 11.3, face = "bold"),
      legend.text = element_text(size = 11.0, color = "black"),
      legend.key.size = unit(0.43, "cm"),
      legend.spacing.y = unit(0.06, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    guides(
      fill = guide_colorbar(
        order = 1,
        barheight = unit(2.56, "cm"),
        barwidth = unit(0.43, "cm")
      ),
      size = guide_legend(order = 2, override.aes = list(alpha = 1)),
      shape = guide_legend(order = 3, override.aes = list(size = 4.1))
    )
}

p_increased <- make_plot(
  plot_df %>% filter(direction == "r/m increased"),
  "r/m increased",
  global_max_x
)

p_decreased <- make_plot(
  plot_df %>% filter(direction == "r/m decreased"),
  "r/m decreased",
  max((plot_df %>% filter(direction == "r/m decreased"))$gene_ratio_percent, na.rm = TRUE)
)

png_file_increased <- file.path(out_dir, "rm_direction_KEGG_enrichment_dotplot_increased.png")
pdf_file_increased <- file.path(out_dir, "rm_direction_KEGG_enrichment_dotplot_increased.pdf")
png_file_decreased <- file.path(out_dir, "rm_direction_KEGG_enrichment_dotplot_decreased.png")
pdf_file_decreased <- file.path(out_dir, "rm_direction_KEGG_enrichment_dotplot_decreased.pdf")

ggsave(png_file_increased, p_increased, width = 7.8, height = 4.6, dpi = 350)
ggsave(pdf_file_increased, p_increased, width = 7.8, height = 4.6)
ggsave(png_file_decreased, p_decreased, width = 8.6, height = 6.4, dpi = 350)
ggsave(pdf_file_decreased, p_decreased, width = 8.6, height = 6.4)

print(p_increased)
print(p_decreased)
message("Saved: ", png_file_increased)
message("Saved: ", pdf_file_increased)
message("Saved: ", png_file_decreased)
message("Saved: ", pdf_file_decreased)
