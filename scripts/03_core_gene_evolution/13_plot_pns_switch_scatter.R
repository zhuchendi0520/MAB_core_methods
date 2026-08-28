suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: Rscript 13_plot_pns_switch_scatter.R ",
    "<pns_switch_posterior.csv> <output_directory>"
  )
}

input_csv <- args[[1]]
out_dir <- args[[2]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

orange <- "#E59673"
blue <- "#8FA7DB"

df <- read.csv(input_csv, check.names = FALSE)

df$pNS_pre_plot <- ifelse(is.na(df$pNS_pre) | df$pNS_pre <= 0, 0.01, df$pNS_pre)
df$pNS_post_plot <- ifelse(is.na(df$pNS_post) | df$pNS_post <= 0, 0.01, df$pNS_post)
df$log10_pre <- log10(df$pNS_pre_plot)
df$log10_post <- log10(df$pNS_post_plot)

df <- df[is.finite(df$log10_pre) & is.finite(df$log10_post), ]

df$Switch <- NA_character_
df$Switch[
  !is.na(df$posterior_switch_probability) &
    df$posterior_switch_probability > 0.9 &
    df$posterior_switch_direction == "purifying_to_positive"
] <- "negative->positive"
df$Switch[
  !is.na(df$posterior_switch_probability) &
    df$posterior_switch_probability > 0.9 &
    df$posterior_switch_direction == "positive_to_purifying"
] <- "positive->negative"

df$Switch <- factor(df$Switch, levels = c("negative->positive", "positive->negative"))

bg_df <- df[is.na(df$Switch), ]
switch_df <- df[!is.na(df$Switch), ]
switch_df <- switch_df[order(switch_df$Switch), ]

bg_df$confidence_for_size <- bg_df$posterior_switch_probability
bg_df$confidence_for_size[is.na(bg_df$confidence_for_size)] <- 0
bg_df$size_for_plot <- 0.65 + pmin(bg_df$confidence_for_size, 0.9) / 0.9 * (2.2 - 0.65)

size_legend <- data.frame(
  probability = c("0.9", "0.5", "0.1"),
  x = c(1.16, 1.16, 1.16),
  y = c(-1.20, -1.46, -1.70),
  point_size = 0.65 + c(0.9, 0.5, 0.1) / 0.9 * (2.2 - 0.65),
  label_x = c(1.31, 1.31, 1.31)
)

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, color = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, color = "#777777") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.45, color = "#999999") +
  geom_point(
    data = bg_df,
    aes(x = log10_post, y = log10_pre, size = size_for_plot),
    shape = 21,
    stroke = 0.28,
    color = "#CFCFCF",
    fill = "#E3E3E3",
    alpha = 0.30
  ) +
  geom_point(
    data = switch_df,
    aes(x = log10_post, y = log10_pre, color = Switch, fill = Switch),
    shape = 21,
    size = 2.9,
    stroke = 0.55,
    alpha = 0.95
  ) +
  annotate(
    "rect",
    xmin = 0.88, xmax = 1.78,
    ymin = -1.92, ymax = -0.98,
    fill = "white",
    color = NA,
    alpha = 0.86
  ) +
  annotate(
    "text",
    x = 1.33, y = -1.04,
    label = "Posterior switch\nprobability",
    size = 2.55,
    fontface = "bold",
    lineheight = 0.9,
    hjust = 0.5
  ) +
  geom_point(
    data = size_legend,
    aes(x = x, y = y, size = point_size),
    shape = 21,
    fill = "#E3E3E3",
    color = "#CFCFCF",
    alpha = 0.65,
    stroke = 0.28,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = size_legend,
    aes(x = label_x, y = y, label = probability),
    size = 2.65,
    fontface = "bold",
    color = "black",
    hjust = 0,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c("negative->positive" = orange, "positive->negative" = blue),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = c("negative->positive" = "#F3B18E", "positive->negative" = "#AFC1E9"),
    drop = TRUE
  ) +
  scale_size_identity(guide = "none") +
  coord_equal(xlim = c(-2.2, 2.0), ylim = c(-2.2, 2.0), clip = "off") +
  labs(
    x = expression(log[10]("post-expansion pN/pS")),
    y = expression(log[10]("pre-expansion pN/pS")),
    color = NULL,
    fill = NULL
  ) +
  guides(
    color = "none",
    fill = guide_legend(override.aes = list(size = 3.0, alpha = 1))
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 9.5),
    axis.title = element_text(face = "bold", size = 11),
    axis.line = element_line(linewidth = 0.75),
    axis.ticks = element_line(linewidth = 0.55),
    legend.position = c(0.28, 0.88),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.key.size = unit(0.42, "cm"),
    legend.text = element_text(size = 8.8),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(file.path(out_dir, "Fig4B_pNS_switch_scatter_y_pre.png"), p, width = 5.0, height = 5.0, dpi = 600)
ggsave(file.path(out_dir, "Fig4B_pNS_switch_scatter_y_pre.pdf"), p, width = 5.0, height = 5.0)

write.csv(
  switch_df[, c(
    "DCC", "locus", "gene", "product", "pNS_pre", "pNS_post",
    "posterior_switch_probability", "posterior_switch_direction", "Switch"
  )],
  file.path(out_dir, "Fig4B_pNS_switch_scatter_y_pre_genes.csv"),
  row.names = FALSE
)
