library(grid)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 04_plot_pangenome_filtering_workflow.R OUTPUT_DIR")
}
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

box_grob <- function(x, y, w, h, title, subtitle = "", fill = "#F7F7F7",
                     col = "#333333", title_size = 11, subtitle_size = 8.6) {
  grobTree(
    roundrectGrob(
      x = unit(x, "npc"), y = unit(y, "npc"),
      width = unit(w, "npc"), height = unit(h, "npc"),
      r = unit(0.018, "npc"),
      gp = gpar(fill = fill, col = col, lwd = 1.1)
    ),
    textGrob(
      title,
      x = unit(x, "npc"), y = unit(y + h * 0.16, "npc"),
      gp = gpar(fontsize = title_size, fontface = "bold", col = "#111111")
    ),
    textGrob(
      subtitle,
      x = unit(x, "npc"), y = unit(y - h * 0.15, "npc"),
      gp = gpar(fontsize = subtitle_size, col = "#222222", lineheight = 0.95)
    )
  )
}

arrow_grob <- function(x1, y1, x2, y2, col = "#555555", lwd = 1.5) {
  segmentsGrob(
    x0 = unit(x1, "npc"), y0 = unit(y1, "npc"),
    x1 = unit(x2, "npc"), y1 = unit(y2, "npc"),
    arrow = arrow(type = "closed", length = unit(0.12, "inches")),
    gp = gpar(col = col, lwd = lwd, lineend = "round")
  )
}

draw_workflow <- function() {
  grid.newpage()
  pushViewport(viewport(width = 0.98, height = 0.96))

  green <- "#DDEBC8"
  blue <- "#DDEAF8"
  orange <- "#F8D9C8"
  yellow <- "#F5E5B5"
  grey <- "#ECECEC"
  core_col <- "#7888C4"
  acc_col <- "#E6C16F"

  grid.text(
    "Pangenome construction and stable core-genome definition",
    x = unit(0.04, "npc"), y = unit(0.95, "npc"),
    just = "left",
    gp = gpar(fontsize = 14, fontface = "bold")
  )

  # Top row
  grid.draw(box_grob(0.13, 0.82, 0.18, 0.11, "DCC-stratified trees",
                     "reference-mapping trees\nfor each DCC", grey))
  grid.draw(box_grob(0.36, 0.82, 0.19, 0.11, "Treemmer sampling",
                     "retain 90% phylogenetic\ndiversity per DCC", green))
  grid.draw(box_grob(0.59, 0.82, 0.18, 0.11, "Representative set",
                     "1,130 genomes", green))
  grid.draw(box_grob(0.82, 0.82, 0.20, 0.11, "Panaroo pangenome",
                     "34,571 gene clusters", blue))

  grid.draw(arrow_grob(0.22, 0.82, 0.265, 0.82))
  grid.draw(arrow_grob(0.455, 0.82, 0.50, 0.82))
  grid.draw(arrow_grob(0.68, 0.82, 0.72, 0.82))

  grid.text(
    "Panaroo gene-cluster classification",
    x = unit(0.55, "npc"), y = unit(0.665, "npc"),
    gp = gpar(fontsize = 10, col = "#333333")
  )

  # Split branches
  grid.draw(box_grob(0.25, 0.62, 0.20, 0.11, "Initial core genome",
                     "3,436 genes", "#E7EAF7", core_col))
  grid.draw(box_grob(0.76, 0.62, 0.22, 0.11, "Accessory genome",
                     "31,135 gene clusters", "#F8EDCC", acc_col))
  grid.draw(arrow_grob(0.77, 0.765, 0.34, 0.675, core_col, 1.5))
  grid.draw(arrow_grob(0.82, 0.765, 0.76, 0.675, acc_col, 1.5))

  # Core path
  grid.draw(box_grob(0.25, 0.44, 0.22, 0.12, "Core alignment QC",
                     "remove low-quality genes\nand genes with gap >=10%",
                     orange, "#333333", 10.8, 8.3))
  grid.draw(box_grob(0.50, 0.44, 0.20, 0.11, "QC-passed core",
                     "3,195 genes", orange))
  grid.draw(box_grob(0.74, 0.44, 0.22, 0.12, "Presence validation",
                     "recheck in remaining\n6,733 genomes", yellow))
  grid.draw(box_grob(0.50, 0.23, 0.28, 0.13, "Stable lower-bound core genome",
                     "present in >99% genomes\n3,001 genes",
                     "#DEEAD6", "#5F8F5F", 11, 8.8))

  grid.draw(arrow_grob(0.25, 0.565, 0.25, 0.50, core_col))
  grid.draw(arrow_grob(0.36, 0.44, 0.40, 0.44, core_col))
  grid.draw(arrow_grob(0.60, 0.44, 0.63, 0.44, core_col))
  grid.draw(arrow_grob(0.70, 0.38, 0.55, 0.295, core_col))

  # Accessory path
  grid.draw(box_grob(0.76, 0.22, 0.22, 0.15, "Accessory analysis",
                     "DCC-specific/shared genes,\ngain/loss and functional\nannotation",
                     "#F8EDCC", acc_col, 10.8, 8.2))
  grid.draw(arrow_grob(0.76, 0.565, 0.76, 0.305, acc_col))

  grid.text(
    "Downstream analyses use the filtered 3,001-gene core as a conservative, stable core-genome set.",
    x = unit(0.50, "npc"), y = unit(0.08, "npc"),
    gp = gpar(fontsize = 9.3, col = "#333333")
  )

  popViewport()
}

png(file.path(out_dir, "pangenome_filtering_workflow.png"), width = 6600, height = 3600, res = 600)
draw_workflow()
dev.off()

pdf(file.path(out_dir, "pangenome_filtering_workflow.pdf"), width = 11, height = 6)
draw_workflow()
dev.off()

svg(file.path(out_dir, "pangenome_filtering_workflow.svg"), width = 11, height = 6)
draw_workflow()
dev.off()
