suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_DCC_node_to_tip_and_TBL.R tree.tre target_nodes.txt output_dir")
}

tree_file <- args[[1]]
target_file <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(tree_file)
targets <- read.table(target_file, header = FALSE, stringsAsFactors = FALSE)
colnames(targets) <- c("DCC", "target_node")
targets$target_node[targets$DCC == "DCC1"] <- "Node_640"
targets$target_node[targets$DCC == "DCC3"] <- "Node_844"

dcc_levels <- paste0("DCC", 1:7)
targets <- targets[match(dcc_levels, targets$DCC), , drop = FALSE]

dcc_colors <- c(
  DCC1 = "#5B86C4", DCC2 = "#D889AE", DCC3 = "#E8A15A",
  DCC4 = "#D9CD45", DCC5 = "#72B7A1", DCC6 = "#9B82C9",
  DCC7 = "#CF7F72"
)

children <- split(tree$edge[, 2], tree$edge[, 1])
descendant_tips <- function(node) {
  current <- node
  tips <- integer()
  while (length(current) > 0) {
    x <- current[[1]]
    current <- current[-1]
    if (x <= Ntip(tree)) {
      tips <- c(tips, x)
    } else {
      current <- c(current, children[[as.character(x)]])
    }
  }
  sort(unique(tips))
}

depth <- node.depth.edgelength(tree)
terminal_edge <- tree$edge[, 2] <= Ntip(tree)
tbl_lookup <- setNames(tree$edge.length[terminal_edge], tree$edge[terminal_edge, 2])

records <- list()
for (i in seq_len(nrow(targets))) {
  dcc <- targets$DCC[[i]]
  label <- targets$target_node[[i]]
  label_index <- match(label, tree$node.label)
  if (is.na(label_index)) stop("Node label not found: ", label)
  node_number <- Ntip(tree) + label_index
  tips <- descendant_tips(node_number)

  records[[i]] <- data.frame(
    DCC = dcc,
    target_node = label,
    tip = tree$tip.label[tips],
    node_to_tip_branch_length = depth[tips] - depth[node_number],
    TBL = unname(tbl_lookup[as.character(tips)]),
    stringsAsFactors = FALSE
  )
}

dat <- do.call(rbind, records)
dat$DCC <- factor(dat$DCC, levels = dcc_levels)
write.csv(dat, file.path(output_dir, "DCC_node_to_tip_and_TBL_by_strain.csv"), row.names = FALSE)

summary_fun <- function(x) {
  c(
    n = length(x), mean = mean(x), median = median(x), sd = sd(x),
    Q1 = unname(quantile(x, 0.25)), Q3 = unname(quantile(x, 0.75)),
    min = min(x), max = max(x)
  )
}
summary_rows <- lapply(split(dat, dat$DCC), function(x) {
  a <- summary_fun(x$node_to_tip_branch_length)
  b <- summary_fun(x$TBL)
  data.frame(
    DCC = as.character(x$DCC[[1]]), target_node = x$target_node[[1]],
    n_tips = a[["n"]],
    node_to_tip_mean = a[["mean"]], node_to_tip_median = a[["median"]],
    node_to_tip_sd = a[["sd"]], node_to_tip_Q1 = a[["Q1"]],
    node_to_tip_Q3 = a[["Q3"]], node_to_tip_min = a[["min"]],
    node_to_tip_max = a[["max"]],
    TBL_mean = b[["mean"]], TBL_median = b[["median"]], TBL_sd = b[["sd"]],
    TBL_Q1 = b[["Q1"]], TBL_Q3 = b[["Q3"]], TBL_min = b[["min"]],
    TBL_max = b[["max"]]
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[match(dcc_levels, summary_df$DCC), ]
write.csv(summary_df, file.path(output_dir, "DCC_node_to_tip_and_TBL_summary.csv"), row.names = FALSE)

base_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(8, 12, 8, 8)
  )

make_boxplot <- function(plot_dat, y, ylab, upper_limit) {
  ggplot(plot_dat, aes(x = DCC, y = .data[[y]], fill = DCC, color = DCC)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.58, linewidth = 0.75) +
    geom_jitter(width = 0.14, height = 0, size = 1.65, alpha = 0.72, stroke = 0) +
    scale_fill_manual(values = dcc_colors) +
    scale_color_manual(values = dcc_colors) +
    scale_y_continuous(limits = c(0, upper_limit), expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = ylab) +
    base_theme
}

node_tip_plot_dat <- dat[!is.na(dat$node_to_tip_branch_length) &
                           dat$node_to_tip_branch_length <= 500, , drop = FALSE]
tbl_plot_dat <- dat[!is.na(dat$TBL) & dat$TBL <= 200, , drop = FALSE]

p_node_tip <- make_boxplot(
  node_tip_plot_dat,
  "node_to_tip_branch_length",
  "Branch length from DCC ancestral node to tip",
  500
)
p_tbl <- make_boxplot(tbl_plot_dat, "TBL", "Terminal branch length (TBL)", 200)

ggsave(file.path(output_dir, "DCC_node_to_tip_branch_length_boxplot.png"), p_node_tip,
       width = 5.3, height = 6.3, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "DCC_node_to_tip_branch_length_boxplot.pdf"), p_node_tip,
       width = 5.3, height = 6.3, device = "pdf")
ggsave(file.path(output_dir, "DCC_TBL_boxplot.png"), p_tbl,
       width = 5.3, height = 6.3, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "DCC_TBL_boxplot.pdf"), p_tbl,
       width = 5.3, height = 6.3, device = "pdf")

cat("DCC tip counts:\n")
print(summary_df[, c("DCC", "target_node", "n_tips")], row.names = FALSE)
