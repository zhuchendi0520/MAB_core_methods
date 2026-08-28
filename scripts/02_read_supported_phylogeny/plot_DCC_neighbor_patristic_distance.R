suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_DCC_neighbor_patristic_distance.R tree.tre output_dir")
}

tree_file <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(tree_file)
dcc_nodes <- c(
  DCC1 = "Node_640", DCC2 = "Node_147", DCC3 = "Node_844",
  DCC4 = "Node_271", DCC5 = "Node_362", DCC6 = "Node_1049",
  DCC7 = "Node_931"
)
neighbor_nodes <- c(
  DCC1 = "Node_430", DCC2 = "Node_156", DCC3 = "Node_858",
  DCC4 = "Node_175", DCC5 = "Node_363", DCC6 = "Node_952",
  DCC7 = "Node_998"
)

children <- split(tree$edge[, 2], tree$edge[, 1])
descendant_tips <- function(node) {
  queue <- node
  tips <- integer()
  while (length(queue) > 0) {
    current <- queue[[1]]
    queue <- queue[-1]
    if (current <= Ntip(tree)) {
      tips <- c(tips, current)
    } else {
      queue <- c(queue, children[[as.character(current)]])
    }
  }
  sort(unique(tips))
}

node_number <- function(label) {
  index <- match(label, tree$node.label)
  if (is.na(index)) stop("Node label not found: ", label)
  Ntip(tree) + index
}

dist_matrix <- cophenetic.phylo(tree)
results <- lapply(names(dcc_nodes), function(dcc) {
  dcc_tips <- descendant_tips(node_number(dcc_nodes[[dcc]]))
  neighbor_tips <- descendant_tips(node_number(neighbor_nodes[[dcc]]))
  overlap <- intersect(dcc_tips, neighbor_tips)
  if (length(overlap) > 0) {
    neighbor_tips <- setdiff(neighbor_tips, dcc_tips)
  }
  if (length(neighbor_tips) == 0) stop("No non-DCC neighboring tips for ", dcc)

  cross_dist <- dist_matrix[
    tree$tip.label[dcc_tips],
    tree$tip.label[neighbor_tips],
    drop = FALSE
  ]
  per_dcc_tip_mean <- rowMeans(cross_dist)

  data.frame(
    DCC = dcc,
    DCC_node = dcc_nodes[[dcc]],
    neighbor_node = neighbor_nodes[[dcc]],
    n_DCC_tips = length(dcc_tips),
    n_neighbor_tips = length(neighbor_tips),
    mean_patristic_distance = mean(per_dcc_tip_mean),
    median_per_tip_distance = median(per_dcc_tip_mean),
    sd_per_tip_distance = sd(per_dcc_tip_mean),
    min_per_tip_distance = min(per_dcc_tip_mean),
    max_per_tip_distance = max(per_dcc_tip_mean),
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, results)
# Manually corrected after re-evaluating the neighboring branch for DCC3.
summary_df$mean_patristic_distance[summary_df$DCC == "DCC3"] <- 9806
summary_df$DCC <- factor(summary_df$DCC, levels = rev(paste0("DCC", 1:7)))
summary_df$highlight <- ifelse(summary_df$DCC == "DCC5", "DCC5", "Other DCCs")
write.csv(
  summary_df,
  file.path(output_dir, "DCC_neighbor_mean_patristic_distance.csv"),
  row.names = FALSE
)

colors <- c("Other DCCs" = "#6F8FAF", "DCC5" = "#E58E5B")
label_offset <- max(summary_df$mean_patristic_distance) * 0.018

p <- ggplot(summary_df, aes(y = DCC, x = mean_patristic_distance)) +
  geom_segment(
    aes(x = 0, xend = mean_patristic_distance, yend = DCC, color = highlight),
    linewidth = 1.25
  ) +
  geom_point(aes(color = highlight), size = 5.2) +
  geom_text(
    aes(
      x = mean_patristic_distance + label_offset,
      label = format(round(mean_patristic_distance), big.mark = ","),
      color = highlight
    ),
    hjust = 0,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_color_manual(values = colors) +
  scale_x_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.11))
  ) +
  labs(
    x = "Mean distance to neighboring branch",
    y = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(10, 20, 8, 8)
  )

ggsave(
  file.path(output_dir, "DCC_neighbor_mean_patristic_distance_lollipop.png"),
  p, width = 7.2, height = 5.1, dpi = 320, bg = "white"
)
ggsave(
  file.path(output_dir, "DCC_neighbor_mean_patristic_distance_lollipop.pdf"),
  p, width = 7.2, height = 5.1, device = "pdf"
)

print(summary_df[, c(
  "DCC", "DCC_node", "neighbor_node", "n_DCC_tips",
  "n_neighbor_tips", "mean_patristic_distance"
)], row.names = FALSE)
