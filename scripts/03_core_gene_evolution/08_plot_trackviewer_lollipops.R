suppressPackageStartupMessages({
  library(trackViewer)
  library(GenomicRanges)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste(
    "Usage: Rscript plot_selected_genes_trackviewer_lollipop.R",
    "mutations.csv uniprot_features.tsv gene_metadata.csv output_dir"
  ))
}

mutation_file <- args[[1]]
uniprot_file <- args[[2]]
metadata_file <- args[[3]]
output_dir <- args[[4]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "individual_png"), showWarnings = FALSE)
dir.create(file.path(output_dir, "individual_pdf"), showWarnings = FALSE)

target_genes <- c(
  "MAB_3029", "MAB_0674", "MAB_1638", "MAB_0173", "MAB_4596c",
  "MAB_3404c", "MAB_1678c", "MAB_1915", "MAB_1499", "MAB_0401",
  "MAB_4595c", "MAB_1060", "MAB_0673", "MAB_2928", "MAB_4470c",
  "MAB_0101", "MAB_4294", "MAB_1076", "MAB_2819", "MAB_4075",
  "MAB_4559", "MAB_0114", "MAB_1218", "MAB_0019", "MAB_3045c",
  "MAB_1635", "MAB_1140", "MAB_3675", "MAB_2689"
)

mutation_colors <- c(
  "Nonsynonymous" = "#E9A184",
  "Synonymous" = "#89BDD2"
)

feature_colors <- c(
  "Domain" = "#A9C8D3",
  "Region" = "#CFDDB2",
  "Transmembrane" = "#F1D18A",
  "Coiled coil" = "#D8B8D1",
  "Repeat" = "#B9C3E1",
  "Motif" = "#E7B6A5"
)

block_palette <- c(
  "#A9C8D3", "#CFDDB2", "#F1D18A", "#D8B8D1", "#B9C3E1",
  "#E7B6A5", "#9FC8B8", "#E6C6A8", "#B7B3D7", "#C7D8E8"
)

mut <- read.csv(mutation_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
uni <- read.delim(
  uniprot_file, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)

extract_gene <- function(x) {
  hit <- regmatches(x, regexpr("MAB_[0-9]+c?", x))
  ifelse(length(hit) == 0 || hit == "", NA_character_, hit)
}
uni$gene <- vapply(uni$`Gene Names`, extract_gene, character(1))

clean_aa <- function(consequence, aa_pos) {
  parts <- strsplit(consequence, "-", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(paste0("aa", aa_pos))
  from <- sub(".*/", "", parts[[2]])
  to <- sub(".*/", "", parts[[3]])
  paste0(from, aa_pos, to)
}

mutation_class <- ifelse(mut$SNP_TYPE == "NSY", "Nonsynonymous", "Synonymous")
mut$mutation_class <- mutation_class
mut$aa_position <- as.integer(mut$AA_OR_CDS_POSITION)
mut$event_ref <- mut$EVENT_REF_ATCC_ORIENTATION
mut$event_alt <- mut$EVENT_ALT_ATCC_ORIENTATION

key_cols <- c(
  "ATCC19977_GENE", "aa_position", "ATCC19977_LOC", "event_ref",
  "event_alt", "mutation_class", "CONSEQUENCE"
)
split_key <- interaction(mut[key_cols], drop = TRUE, lex.order = TRUE)
agg_list <- lapply(split(mut, split_key), function(x) {
  data.frame(
    ATCC19977_GENE = x$ATCC19977_GENE[[1]],
    aa_position = x$aa_position[[1]],
    ATCC19977_LOC = x$ATCC19977_LOC[[1]],
    event_ref = x$event_ref[[1]],
    event_alt = x$event_alt[[1]],
    mutation_class = x$mutation_class[[1]],
    CONSEQUENCE = x$CONSEQUENCE[[1]],
    frequency = nrow(x),
    dcc_count = length(unique(x$DCC)),
    dccs = paste(sort(unique(x$DCC)), collapse = ";"),
    stringsAsFactors = FALSE
  )
})
agg <- do.call(rbind, agg_list)
agg$mutation_label <- mapply(clean_aa, agg$CONSEQUENCE, agg$aa_position)
agg$nucleotide_change <- paste0(agg$event_ref, ">", agg$event_alt)
agg <- agg[agg$ATCC19977_GENE %in% target_genes, ]
agg <- agg[order(match(agg$ATCC19977_GENE, target_genes), agg$aa_position), ]

# Collapse alternative nucleotide changes at the same amino-acid position for
# plotting. The detailed event-level table above is still written unchanged.
plot_groups <- split(
  agg,
  interaction(agg$ATCC19977_GENE, agg$aa_position, drop = TRUE, lex.order = TRUE)
)
plot_agg <- do.call(rbind, lapply(plot_groups, function(x) {
  cls <- if (any(x$mutation_class == "Nonsynonymous")) {
    "Nonsynonymous"
  } else {
    "Synonymous"
  }
  aa_labels <- unique(x$mutation_label[x$mutation_class == cls])
  aa_labels <- aa_labels[!is.na(aa_labels) & aa_labels != ""]
  # Label every observed amino-acid site, including singletons.
  display_label <- paste(head(aa_labels, 2L), collapse = "/")
  data.frame(
    ATCC19977_GENE = x$ATCC19977_GENE[[1]],
    aa_position = x$aa_position[[1]],
    mutation_class = cls,
    frequency = sum(x$frequency),
    dcc_count = length(unique(unlist(strsplit(x$dccs, ";", fixed = TRUE)))),
    display_label = display_label,
    stringsAsFactors = FALSE
  )
}))
plot_agg <- plot_agg[
  order(match(plot_agg$ATCC19977_GENE, target_genes), plot_agg$aa_position),
]

parse_uniprot_column <- function(gene, column, type_label) {
  row <- uni[uni$gene == gene, , drop = FALSE]
  if (nrow(row) == 0 || !(column %in% colnames(row))) return(NULL)
  text <- row[[column]][[1]]
  if (is.na(text) || text == "") return(NULL)

  prefix <- switch(type_label,
    "Domain" = "DOMAIN",
    "Region" = "REGION",
    "Transmembrane" = "TRANSMEM",
    "Coiled coil" = "COILED",
    "Repeat" = "REPEAT",
    "Motif" = "MOTIF"
  )
  starts <- gregexpr(paste0(prefix, " [0-9]+(?:\\.\\.[0-9]+)?"), text, perl = TRUE)[[1]]
  if (starts[[1]] == -1) return(NULL)
  tokens <- regmatches(text, gregexpr(paste0(prefix, " [0-9]+(?:\\.\\.[0-9]+)?"), text, perl = TRUE))[[1]]
  note_matches <- regmatches(text, gregexpr('/note="[^"]+"', text, perl = TRUE))[[1]]

  out <- vector("list", length(tokens))
  for (i in seq_along(tokens)) {
    coord <- sub(paste0("^", prefix, " "), "", tokens[[i]])
    parts <- as.integer(strsplit(coord, "..", fixed = TRUE)[[1]])
    label <- if (length(note_matches) >= i && note_matches[[1]] != "") {
      gsub('^/note="|"$', "", note_matches[[i]])
    } else {
      type_label
    }
    out[[i]] <- data.frame(
      gene = gene,
      feature_type = type_label,
      start = parts[[1]],
      end = ifelse(length(parts) == 2, parts[[2]], parts[[1]]),
      label = label,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

parsed_features <- do.call(rbind, Filter(Negate(is.null), unlist(lapply(target_genes, function(gene) {
  list(
    parse_uniprot_column(gene, "Domain [FT]", "Domain"),
    parse_uniprot_column(gene, "Region", "Region"),
    parse_uniprot_column(gene, "Transmembrane", "Transmembrane"),
    parse_uniprot_column(gene, "Coiled coil", "Coiled coil"),
    parse_uniprot_column(gene, "Repeat", "Repeat"),
    parse_uniprot_column(gene, "Motif", "Motif")
  )
}), recursive = FALSE)))

if (is.null(parsed_features)) {
  parsed_features <- data.frame(
    gene = character(), feature_type = character(), start = integer(),
    end = integer(), label = character()
  )
}
write.csv(
  parsed_features,
  file.path(output_dir, "uniprot_selected_genes_parsed_blocks.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  agg,
  file.path(output_dir, "selected_genes_lollipop_aggregated_mutations.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  plot_agg,
  file.path(output_dir, "selected_genes_lollipop_amino_acid_sites.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

make_track_objects <- function(gene) {
  dat <- plot_agg[plot_agg$ATCC19977_GENE == gene, , drop = FALSE]
  # Keep every observed site. Singleton mutations are biologically relevant
  # on the pre-expansion branches and must not be removed by recurrence filters.
  gene_meta <- meta[meta$gene == gene, , drop = FALSE]
  protein_length <- if (nrow(gene_meta) > 0) gene_meta$protein_length_aa[[1]] else max(dat$aa_position)
  product <- if (nrow(gene_meta) > 0) gene_meta$product[[1]] else ""

  labels <- dat$display_label
  snp <- GRanges(
    "protein",
    IRanges(dat$aa_position, width = 1, names = labels)
  )
  snp$score <- dat$frequency
  snp$color <- unname(mutation_colors[dat$mutation_class])
  snp$border <- "#4A504E"
  snp$alpha <- 0.92
  snp$SNPsideID <- "top"
  snp$label.parameter.rot <- 58

  block_df <- parsed_features[parsed_features$gene == gene, , drop = FALSE]
  block_labels <- character()
  block_colors <- character()
  if (nrow(block_df) > 0) {
    block_labels <- unique(block_df$label)
    block_colors <- setNames(
      rep(block_palette, length.out = length(block_labels)),
      block_labels
    )
    blocks <- GRanges(
      "protein",
      IRanges(
        block_df$start,
        end = block_df$end,
        names = rep("", nrow(block_df))
      )
    )
    blocks$fill <- unname(block_colors[block_df$label])
    blocks$color <- "#596864"
    blocks$height <- 0.028
    blocks$featureLayerID <- "UniProt block"
    features <- blocks
  } else {
    features <- GRanges(
      "protein",
      IRanges(1, width = protein_length, names = "")
    )
    features$fill <- "transparent"
    features$color <- "#596864"
    features$height <- 0
    features$featureLayerID <- "protein line"
  }

  list(
    snp = snp,
    features = features,
    protein_length = protein_length,
    product = product,
    mutation_classes = unique(dat$mutation_class),
    block_labels = block_labels,
    block_colors = block_colors
  )
}

draw_gene <- function(gene, newpage = TRUE) {
  obj <- make_track_objects(gene)
  x_ticks <- pretty(c(1, obj$protein_length), n = 6)
  x_ticks <- x_ticks[x_ticks >= 1 & x_ticks <= obj$protein_length]

  legend_labels <- c(obj$mutation_classes, obj$block_labels)
  legend_fill <- c(
    unname(mutation_colors[obj$mutation_classes]),
    unname(obj$block_colors[obj$block_labels])
  )
  legend_info <- list(labels = legend_labels, fill = legend_fill)

  lolliplot(
    obj$snp,
    obj$features,
    ranges = GRanges("protein", IRanges(1, obj$protein_length)),
    xaxis = x_ticks,
    cex = 0.82,
    legend = legend_info,
    legendPosition = "top",
    ylab = "Mutation frequency",
    ylab.gp = gpar(fontsize = 11, fontface = "bold"),
    xaxis.gp = gpar(fontsize = 9.5),
    yaxis.gp = gpar(fontsize = 9.5),
    jitter = c("node", "label"),
    label_on_feature = TRUE,
    lollipop_style_switch_limit = 100,
    dashline.col = "#D5D9D7",
    newpage = newpage
  )
  grid.text(
    gene, x = 0.04, y = 0.965, just = c("left", "top"),
    gp = gpar(fontsize = 15, fontface = "bold")
  )
  grid.text(
    "Protein position (amino acid)", x = 0.53, y = 0.035,
    gp = gpar(fontsize = 10.5, fontface = "bold")
  )
}

for (gene in target_genes) {
  png(
    file.path(output_dir, "individual_png", paste0(gene, "_lollipop.png")),
    width = 8.5, height = 5.6, units = "in", res = 320, bg = "white"
  )
  draw_gene(gene)
  dev.off()

  pdf(
    file.path(output_dir, "individual_pdf", paste0(gene, "_lollipop.pdf")),
    width = 8.5, height = 5.6, bg = "white"
  )
  draw_gene(gene)
  dev.off()
}

combined_pdf <- file.path(output_dir, "pre_29_genes_all_mutations_trackViewer_lollipop.pdf")
pdf(combined_pdf, width = 8.5, height = 5.6, onefile = TRUE, useDingbats = FALSE)
for (gene in target_genes) draw_gene(gene)
dev.off()

cat("Genes plotted:", length(target_genes), "\n")
cat("Aggregated mutation sites:", nrow(agg), "\n")
cat("UniProt blocks:", nrow(parsed_features), "\n")
cat("Combined PDF:", combined_pdf, "\n")
