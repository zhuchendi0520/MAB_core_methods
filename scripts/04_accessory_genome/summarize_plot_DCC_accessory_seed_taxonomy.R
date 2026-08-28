#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "DCC_accessory_seed_species_summary.csv"
outdir <- if (length(args) >= 2) args[[2]] else dirname(normalizePath(input))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dcc_levels <- c("DCC1", "DCC2", "DCC4", "DCC5", "DCC3", "DCC6", "DCC7")
top_n_genera <- 15
mycobacterial_genera <- c(
  "Mycobacterium", "Mycobacteroides", "Mycolicibacterium",
  "Mycolicibacter", "Mycolicibacillus"
)

clean_species <- function(x) {
  words <- str_split(str_squish(x), "\\s+")[[1]]
  if (length(words) < 2) return(words[[1]])
  if (words[[2]] == "sp.") return(paste(words[1:2], collapse = " "))
  paste(words[1:2], collapse = " ")
}

raw <- read_csv(input, show_col_types = FALSE) %>%
  mutate(
    DCC = factor(DCC, levels = dcc_levels),
    seed_scientific_name = na_if(str_squish(seed_scientific_name), ""),
    genus = case_when(
      is.na(seed_scientific_name) ~ "Unclassified",
      str_detect(seed_scientific_name, regex("^(uncultured|unclassified|environmental)", TRUE)) ~ "Unclassified",
      TRUE ~ word(seed_scientific_name, 1)
    ),
    species = if_else(
      is.na(seed_scientific_name),
      "Unclassified",
      vapply(seed_scientific_name, clean_species, character(1))
    ),
    genus_broad = if_else(
      genus %in% mycobacterial_genera,
      "Mycobacterium sensu lato",
      genus
    )
  )

# Genus-level composition among all eggNOG seed-annotated accessory genes.
genus_all <- raw %>%
  group_by(DCC, genus = genus_broad) %>%
  summarise(gene_count = sum(gene_count, na.rm = TRUE), .groups = "drop")

top_genera <- genus_all %>%
  group_by(genus) %>%
  summarise(total_gene_count = sum(gene_count), .groups = "drop") %>%
  arrange(desc(total_gene_count), genus) %>%
  slice_head(n = top_n_genera) %>%
  pull(genus)

genus_summary <- genus_all %>%
  mutate(genus_group = if_else(genus %in% top_genera, genus, "Other genera")) %>%
  group_by(DCC, genus = genus_group) %>%
  summarise(gene_count = sum(gene_count), .groups = "drop") %>%
  group_by(DCC) %>%
  mutate(
    seed_annotated_total = sum(gene_count),
    percent_within_DCC = 100 * gene_count / seed_annotated_total
  ) %>%
  ungroup()

genus_order <- genus_summary %>%
  group_by(genus) %>%
  summarise(total_gene_count = sum(gene_count), .groups = "drop") %>%
  arrange(total_gene_count) %>%
  pull(genus)

genus_plot_df <- genus_summary %>%
  complete(DCC = factor(dcc_levels, levels = dcc_levels), genus = genus_order,
           fill = list(gene_count = 0, percent_within_DCC = 0)) %>%
  mutate(genus = factor(genus, levels = genus_order))

write_csv(genus_summary %>% arrange(DCC, desc(percent_within_DCC)),
          file.path(outdir, "DCC_accessory_seed_genus_summary.csv"))

heat_scale <- scale_fill_gradientn(
  colours = c("#FFFFFF", "#DCEEF4", "#86BED1", "#397D9C", "#173F5F"),
  values = scales::rescale(c(0, 2, 10, 35, 100)),
  limits = c(0, 100),
  name = "Seed-assigned\ngenes (%)"
)

p_genus <- ggplot(genus_plot_df, aes(DCC, genus, fill = percent_within_DCC)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  heat_scale +
  coord_fixed(ratio = 0.72) +
  labs(x = NULL, y = NULL, title = "Accessory-gene eggNOG seed assignments by genus") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "italic"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(file.path(outdir, "DCC_accessory_seed_genus_heatmap.pdf"), p_genus,
       width = 7.2, height = 7.0, device = "pdf")
ggsave(file.path(outdir, "DCC_accessory_seed_genus_heatmap.png"), p_genus,
       width = 7.2, height = 7.0, dpi = 400, bg = "white")

# Species composition within Mycobacterium sensu lato seed assignments.
myco_summary <- raw %>%
  filter(genus %in% mycobacterial_genera) %>%
  group_by(DCC, genus, species) %>%
  summarise(gene_count = sum(gene_count, na.rm = TRUE), .groups = "drop") %>%
  group_by(DCC) %>%
  mutate(
    mycobacterial_seed_total = sum(gene_count),
    percent_within_mycobacterial_seeds = 100 * gene_count / mycobacterial_seed_total
  ) %>%
  ungroup()

species_order <- myco_summary %>%
  group_by(species) %>%
  summarise(total_gene_count = sum(gene_count), .groups = "drop") %>%
  arrange(total_gene_count) %>%
  pull(species)

myco_plot_df <- myco_summary %>%
  complete(DCC = factor(dcc_levels, levels = dcc_levels), species = species_order,
           fill = list(gene_count = 0, percent_within_mycobacterial_seeds = 0)) %>%
  mutate(species = factor(species, levels = species_order))

write_csv(myco_summary %>% arrange(DCC, desc(percent_within_mycobacterial_seeds)),
          file.path(outdir, "DCC_accessory_seed_mycobacterial_species_summary.csv"))

p_species <- ggplot(myco_plot_df, aes(DCC, species, fill = percent_within_mycobacterial_seeds)) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  heat_scale +
  coord_fixed(ratio = 0.56) +
  labs(
    x = NULL, y = NULL,
    title = "Mycobacterial eggNOG seed assignments by species",
    subtitle = "Percentages are normalized within mycobacterial seed assignments for each DCC"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "italic", size = 8.5),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 12, 8, 8)
  )

species_height <- max(7.2, 0.27 * length(species_order) + 2.0)
ggsave(file.path(outdir, "DCC_accessory_seed_mycobacterial_species_heatmap.pdf"), p_species,
       width = 7.5, height = species_height, device = "pdf")
ggsave(file.path(outdir, "DCC_accessory_seed_mycobacterial_species_heatmap.png"), p_species,
       width = 7.5, height = species_height, dpi = 400, bg = "white")

message("Wrote genus and mycobacterial-species summaries and heatmaps to: ", outdir)
