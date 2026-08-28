suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript 12_plot_rm_pre_post_by_gene_class.R ",
    "<rm_statistics.csv> <output_directory>"
  )
}

input_file <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(
  input_file,
  check.names = FALSE
)

hmm_candidates <- names(df)[startsWith(names(df), "HMM call")]
if (length(hmm_candidates) == 0) {
  stop("No column beginning with 'HMM call' was found")
}
hmm_col <- hmm_candidates[[1]]

df_plot <- df %>%
  filter(
    .data[[hmm_col]] %in% c("ES", "NE"),
    !is.na(`r/m_raw_pre`),
    !is.na(`r/m_raw_post`),
    !is.na(`r/m_post-pre`),
    `r/m_post-pre` != 0
  ) %>%
  mutate(
    HMM = factor(
      .data[[hmm_col]],
      levels = c("ES", "NE")
    ),
    Direction = if_else(
      `r/m_post-pre` > 0,
      "Increased r/m",
      "Decreased r/m"
    ),
    Direction = factor(
      Direction,
      levels = c("Decreased r/m", "Increased r/m")
    ),
    Group = interaction(HMM, Direction, sep = " | ", lex.order = TRUE),
    GENE = ATCC19977
  )

df_long <- df_plot %>%
  select(
    GENE,
    HMM,
    Direction,
    Group,
    `r/m_raw_pre`,
    `r/m_raw_post`
  ) %>%
  pivot_longer(
    cols = c(`r/m_raw_pre`, `r/m_raw_post`),
    names_to = "Phase",
    values_to = "RM"
  ) %>%
  mutate(
    Phase = factor(
      Phase,
      levels = c("r/m_raw_pre", "r/m_raw_post"),
      labels = c("Pre-expansion", "Post-expansion")
    )
  )

group_levels <- c(
  "ES | Decreased r/m",
  "ES | Increased r/m",
  "NE | Decreased r/m",
  "NE | Increased r/m"
)

df_long$Group <- factor(df_long$Group, levels = group_levels)

group_fill <- c(
  "ES | Decreased r/m" = "#F2B8A0",
  "ES | Increased r/m" = "#F2B8A0",
  "NE | Decreased r/m" = "#BEC8E6",
  "NE | Increased r/m" = "#BEC8E6"
)

group_line <- c(
  "ES | Decreased r/m" = "#DD8F70",
  "ES | Increased r/m" = "#DD8F70",
  "NE | Decreased r/m" = "#8FA3D6",
  "NE | Increased r/m" = "#8FA3D6"
)

panel_title <- c(
  "ES | Decreased r/m" = "ES genes: decreased r/m",
  "ES | Increased r/m" = "ES genes: increased r/m",
  "NE | Decreased r/m" = "NE genes: decreased r/m",
  "NE | Increased r/m" = "NE genes: increased r/m"
)

plot_list <- list()

for (grp in group_levels) {
  sub <- df_long %>%
    filter(Group == grp)

  n_gene <- n_distinct(sub$GENE)

  p <- ggplot(
    sub,
    aes(
      x = Phase,
      y = RM
    )
  ) +
    geom_line(
      aes(group = GENE),
      color = "grey78",
      alpha = 0.35,
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_boxplot(
      width = 0.24,
      outlier.shape = NA,
      linewidth = 0.75,
      fill = group_fill[grp],
      color = group_line[grp],
      alpha = 0.95
    ) +
    geom_point(
      shape = 21,
      fill = "white",
      color = group_line[grp],
      stroke = 0.45,
      size = 1.6,
      alpha = 0.9,
      position = position_jitter(
        width = 0.035,
        height = 0
      )
    ) +
    labs(
      title = paste0(panel_title[grp], " (n=", n_gene, ")"),
      x = NULL,
      y = "r/m"
    ) +
    theme_classic(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0
      ),
      axis.text.x = element_text(
        size = 10.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 10.5,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 12
      ),
      axis.line = element_line(
        linewidth = 0.6,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.5,
        color = "black"
      ),
      plot.margin = margin(6, 8, 6, 6)
    )

  plot_list[[grp]] <- p
}

final_plot <- (
  plot_list[["ES | Decreased r/m"]] +
    plot_list[["ES | Increased r/m"]] +
    plot_list[["NE | Decreased r/m"]] +
    plot_list[["NE | Increased r/m"]]
) +
  plot_layout(ncol = 2)

pdf_file <- file.path(
  out_dir,
  "RM_ES_NE_pre_post_direction_boxplot.pdf"
)

png_file <- file.path(
  out_dir,
  "RM_ES_NE_pre_post_direction_boxplot.png"
)

ggsave(
  pdf_file,
  final_plot,
  width = 7.2,
  height = 6.2
)

ggsave(
  png_file,
  final_plot,
  width = 7.2,
  height = 6.2,
  dpi = 300
)

print(final_plot)
message("Saved: ", pdf_file)
message("Saved: ", png_file)
