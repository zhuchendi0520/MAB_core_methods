#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_plot_DCC_pangenome_composition_polar.R INPUT_CSV OUTPUT_DIR")
}
input_csv <- args[1]
out_dir <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dcc_order <- paste0("DCC", 1:7)
component_order <- c("Core", "Soft Core", "Shell", "Cloud")

palette <- c(
  "Core" = "#717CB5",
  "Soft Core" = "#D5E2B5",
  "Shell" = "#F1DDA2",
  "Cloud" = "#DEA187"
)

wide <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("DCC", component_order)
missing <- setdiff(required, names(wide))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

long <- wide %>%
  select(all_of(required)) %>%
  pivot_longer(
    cols = all_of(component_order),
    names_to = "Component",
    values_to = "Count"
  ) %>%
  mutate(
    DCC = factor(DCC, levels = dcc_order),
    Component = factor(Component, levels = component_order),
    Count = as.numeric(Count)
  ) %>%
  group_by(DCC) %>%
  mutate(
    Total = sum(Count, na.rm = TRUE),
    Percentage = Count / Total * 100
  ) %>%
  arrange(DCC, Component) %>%
  mutate(
    xmax = cumsum(Percentage),
    xmin = lag(xmax, default = 0)
  ) %>%
  ungroup() %>%
  mutate(
    ring = 8 - match(as.character(DCC), dcc_order),
    ymin = ring - 0.36,
    ymax = ring + 0.36
  )

write.csv(
  long %>%
    select(DCC, Component, Count, Total, Percentage),
  file.path(out_dir, "DCC_pangenome_composition_percent.csv"),
  row.names = FALSE
)

dcc_labels <- data.frame(
  DCC = dcc_order,
  ring = 7:1,
  x = -7.2
)

tick_labels <- data.frame(
  x = c(0, 25, 50, 75, 100),
  y = 7.78,
  label = paste0(c(0, 25, 50, 75, 100), "%")
)

outer_arc <- data.frame(
  x = seq(0, 100, length.out = 401),
  y = 7.58
)

p <- ggplot() +
  geom_rect(
    data = long,
    aes(
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax,
      fill = Component
    ),
    color = "#262626",
    linewidth = 0.62
  ) +
  geom_path(
    data = outer_arc,
    aes(x = x, y = y),
    linewidth = 0.55,
    color = "#262626"
  ) +
  geom_text(
    data = dcc_labels,
    aes(x = x, y = ring, label = DCC),
    family = "Helvetica",
    fontface = "bold",
    size = 4.4,
    hjust = 1.15
  ) +
  geom_text(
    data = tick_labels,
    aes(x = x, y = y, label = label),
    family = "Helvetica",
    size = 3.6
  ) +
  scale_fill_manual(
    values = palette,
    breaks = component_order,
    drop = FALSE
  ) +
  scale_x_continuous(limits = c(-14, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3.0, 8.15), expand = c(0, 0)) +
  coord_polar(
    theta = "x",
    start = -14 / 114 * 2 * pi,
    direction = 1,
    clip = "off"
  ) +
  labs(
    title = "Pangenome composition across DCCs",
    fill = NULL
  ) +
  theme_void(base_family = "Helvetica", base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold", size = 18, hjust = 0.5,
      margin = margin(b = 8)
    ),
    legend.position = c(0.18, 0.78),
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(0.62, "cm"),
    legend.spacing.y = grid::unit(0.08, "cm"),
    plot.margin = margin(16, 18, 16, 18)
  )

ggsave(
  file.path(out_dir, "DCC_pangenome_composition_percent_polar.png"),
  p,
  width = 8.5,
  height = 8.5,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(out_dir, "DCC_pangenome_composition_percent_polar.pdf"),
  p,
  width = 8.5,
  height = 8.5,
  device = "pdf",
  bg = "white"
)

message("Saved outputs to: ", out_dir)
