#!/usr/bin/env Rscript ✳️

# all 78 scaffolds, split into three rows by
# telomere-based classification. Genes/Exons/TEs/Cov/GC content

library(tidyverse)
library(patchwork)
library(ragg)
library(svglite)

# ══════════════════════════════════════════════════════════════════════════════
# COLOURS & PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

COLOURS <- list(
  coverage = "grey55",
  gc       = "#2166ac",
  gene     = "#d6604d",
  exon     = "#f4a343",
  te       = "#4dac26",
  telo     = "black",
  chrom    = "grey75"
)

GC_MIN    <- 0.36
GC_MAX    <- 0.48
BASE_SIZE <- 11

# ══════════════════════════════════════════════════════════════════════════════
# SCAFFOLD CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════

class_df <- read.table("scaffold_classification.tsv", header = TRUE, sep = "\t")

karyo_all <- read.table("scaffold_lengths.txt", header = FALSE, col.names = c("scaffold", "length"))

karyo_all <- karyo_all |>
  left_join(class_df |> select(scaffold, class), by = "scaffold") |>
  mutate(class = replace_na(class, "No_telomere"))

order_complete <- karyo_all |> filter(class == "Complete") |>
  arrange(desc(length)) |> pull(scaffold)
order_single   <- karyo_all |> filter(class == "Single_telomere") |>
  arrange(desc(length)) |> pull(scaffold)
order_none     <- karyo_all |> filter(class %in% c("No_telomere", "Internal_only")) |>
  arrange(desc(length)) |> pull(scaffold)

message("Complete (T2T): ", length(order_complete))
message("Single telomere: ", length(order_single))
message("No telomere / other: ", length(order_none))
message("Total: ", length(order_complete) + length(order_single) + length(order_none))

# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

read_track <- function(file, track_name, scaffold_order) {
  read.table(file, header = FALSE,
             col.names = c("scaffold", "start", "end", "value")) |>
    filter(scaffold %in% scaffold_order) |>
    mutate(
      mid      = (start + end) / 2,
      property = track_name,
      scaffold = factor(scaffold, levels = scaffold_order)
    )
}

mb_breaks <- function(limits) {
  hi <- limits[2]
  if (hi <= 1e6) return(c(0, hi))
  seq(0, floor(hi / 1e6) * 1e6, by = 1e6)
}

mb_labels <- function(breaks) {
  mb <- breaks / 1e6
  paste0(ifelse(mb == floor(mb), as.integer(mb), round(mb, 2)), " Mb")
}

x_scale_mb <- scale_x_continuous(
  breaks = mb_breaks,
  labels = mb_labels,
  expand = expansion(0, 0)
)

base_theme <- theme_bw(base_size = BASE_SIZE) + theme(
  panel.grid.minor    = element_blank(),
  panel.grid.major.x  = element_blank(),
  panel.spacing.x     = unit(1, "pt"),
  panel.spacing.y     = unit(1, "pt"),
  strip.background    = element_rect(fill = "grey95", colour = NA),
  strip.text.y        = element_text(size = BASE_SIZE, angle = 0, hjust = 0),
  axis.text           = element_text(size = BASE_SIZE - 1),
  axis.title          = element_text(size = BASE_SIZE),
  axis.text.x         = element_blank(),
  axis.ticks.x        = element_blank(),
  axis.title.x        = element_blank()
)

scaffold_labels <- theme(
  strip.text.x = element_text(size = BASE_SIZE - 3, angle = 90,
                              hjust = 0, margin = margin(b = 2))
)
no_scaffold_labels <- theme(strip.text.x = element_blank())

make_facet <- function() facet_grid(. ~ scaffold, scales = "free_x", space = "free_x")

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED HEATMAP+BAR FACTORY (genes / exons / TEs)
# Legend stays attached to each panel here — we strip/collect it later
# ══════════════════════════════════════════════════════════════════════════════

make_combined_plot <- function(df, track_colour, track_label) {
  ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
              fill = track_colour, alpha = 0.06, colour = NA) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value, fill = value),
              colour = NA) +
    scale_fill_gradient(
      low = adjustcolor(track_colour, alpha.f = 0.15), high = track_colour,
      limits = c(0, 1), name = track_label,
      guide = guide_colorbar(barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm"),
                             title.position = "top", title.hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "50%", "100%"), expand = expansion(0, 0)) +
    x_scale_mb + make_facet() +
    labs(y = track_label) +
    base_theme + no_scaffold_labels +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = BASE_SIZE - 1),
          legend.text  = element_text(size = BASE_SIZE - 2))
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILD ONE FULL ROW
# legend.position = "none" on EVERY panel inside the row — legends
# only exist on a separate standalone set of plots built once at the end
# ══════════════════════════════════════════════════════════════════════════════

build_row <- function(scaffold_order, row_title) {

  if (length(scaffold_order) == 0) return(NULL)

  karyo <- karyo_all |>
    filter(scaffold %in% scaffold_order) |>
    mutate(scaffold = factor(scaffold, levels = scaffold_order))

  tel <- read.table("telomeres_highlight.txt", header = FALSE,
                    col.names = c("scaffold", "start", "end")) |>
    filter(scaffold %in% scaffold_order) |>
    mutate(scaffold = factor(scaffold, levels = scaffold_order))

  coverage <- read_track("coverage_100kb.txt",       "Coverage", scaffold_order)
  gc       <- read_track("GC_100kb.txt",              "GC %",    scaffold_order)
  gene     <- read_track("gene_fraction_100kb.txt",   "Genes",   scaffold_order)
  exon     <- read_track("exon_bp_100kb.bed",         "Exons",   scaffold_order)
  te       <- read_track("TE_fraction_100kb.txt",     "TEs",     scaffold_order)

  # ---- Coverage ----
  cov_max    <- ceiling(max(coverage$value, na.rm = TRUE) / 50) * 50
  cov_median <- median(coverage$value, na.rm = TRUE)
  cov_breaks <- pretty(c(0, cov_max), n = 3)

  p_cov <- ggplot(coverage) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value),
              fill = COLOURS$coverage, colour = NA, alpha = 0.85) +
    geom_hline(yintercept = cov_median, linewidth = 0.25,
               linetype = "dashed", colour = "grey30") +
    scale_y_continuous(limits = c(0, cov_max), breaks = cov_breaks,
                       labels = paste0(cov_breaks, "x"), expand = expansion(0, 0)) +
    x_scale_mb + make_facet() +
    labs(y = "Coverage", title = row_title) +
    base_theme + scaffold_labels +
    theme(plot.title = element_text(size = BASE_SIZE + 2, face = "bold", hjust = 0))

  # ---- GC ----
  p_gc <- ggplot(gc) +
    geom_hline(yintercept = mean(c(GC_MIN, GC_MAX)), linewidth = 0.25,
               colour = "grey70", linetype = "dashed") +
    geom_rect(aes(xmin = start, xmax = end,
                 ymin = mean(c(GC_MIN, GC_MAX)),
                 ymax = pmin(pmax(value, GC_MIN), GC_MAX)),
              fill = COLOURS$gc, colour = NA, alpha = 0.75) +
    scale_y_continuous(limits = c(GC_MIN, GC_MAX),
                       breaks = c(GC_MIN, mean(c(GC_MIN, GC_MAX)), GC_MAX),
                       labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(0, 0)) +
    x_scale_mb + make_facet() +
    labs(y = "GC %") +
    base_theme + no_scaffold_labels

  p_gene <- make_combined_plot(gene, COLOURS$gene, "Genes")
  p_exon <- make_combined_plot(exon, COLOURS$exon, "Exons")
  p_te   <- make_combined_plot(te,   COLOURS$te,   "TEs")

  # ---- Ideogram ----
  ideo <- ggplot() +
    geom_segment(data = karyo, aes(x = 0, xend = length, y = 0.5, yend = 0.5),
                 linewidth = 2, colour = COLOURS$chrom, lineend = "round") +
    geom_segment(data = tel, aes(x = start, xend = end, y = 0.5, yend = 0.5),
                 linewidth = 3, colour = COLOURS$telo, lineend = "round") +
    geom_text(data = karyo,
              aes(x = length / 2, y = 0.08,
                  label = paste0(round(length / 1e6, 2), " Mb")),
              size = (BASE_SIZE - 4) / .pt, colour = "grey30", vjust = 1) +
    x_scale_mb +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    facet_grid(. ~ scaffold, scales = "free_x", space = "free_x") +
    labs(x = "Position", y = NULL) +
    base_theme + no_scaffold_labels +
    theme(
      axis.text.x  = element_text(size = BASE_SIZE - 3, angle = 35, hjust = 1),
      axis.ticks.x = element_line(linewidth = 0.3),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      panel.grid   = element_blank()
    )

  # Strip legends from EVERY panel in this row — they'll be added back
  # once, externally, at the bottom of the whole figure
  (p_cov / p_gc / p_gene / p_exon / p_te / ideo) +
    plot_layout(heights = c(0.35, 0.35, 0.9, 0.9, 0.9, 0.3)) &
    theme(plot.margin = margin(b = 25, t = 5),
          legend.position = "none")
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILD ALL THREE ROWS (no legends — stripped inside build_row via `&`)
# ══════════════════════════════════════════════════════════════════════════════

row1 <- build_row(order_complete,
                  paste0("A. Complete chromosomes (T2T)  —  n = ", length(order_complete)))
row2 <- build_row(order_single,
                  paste0("B. Single telomere  —  n = ", length(order_single)))
row3 <- build_row(order_none,
                  paste0("C. No telomere / fragments  —  n = ", length(order_none)))

# ══════════════════════════════════════════════════════════════════════════════
# BUILD ONE STANDALONE LEGEND STRIP
# Use representative data (any scaffold set works — colour scale is fixed
# whose ONLY purpose is to produce a clean horizontal legend each.
# ══════════════════════════════════════════════════════════════════════════════

legend_gene <- make_combined_plot(
  read_track("gene_fraction_100kb.txt", "Genes", order_complete),
  COLOURS$gene, "Genes"
) + theme(legend.position = "bottom")

legend_exon <- make_combined_plot(
  read_track("Spa255.hifiasm.v25.ragtag.pecat.braker.dedup_exon_bp_100kb.bed",
            "Exons", order_complete),
  COLOURS$exon, "Exons"
) + theme(legend.position = "bottom")

legend_te <- make_combined_plot(
  read_track("TE_fraction_100kb.txt", "TEs", order_complete),
  COLOURS$te, "TEs"
) + theme(legend.position = "bottom")

# Extract ONLY the legend grob from each (cowplot-free approach using
# patchwork's own guide-extraction via a 1-row dummy layout)
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  legend_idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[legend_idx]]
}

legend_gene_grob <- extract_legend(legend_gene)
legend_exon_grob <- extract_legend(legend_exon)
legend_te_grob   <- extract_legend(legend_te)

legend_strip <- wrap_elements(full = legend_gene_grob) +
                wrap_elements(full = legend_exon_grob) +
                wrap_elements(full = legend_te_grob) +
  plot_layout(nrow = 1)

# ══════════════════════════════════════════════════════════════════════════════
# STACK: 3 rows + 1 shared legend strip at the bottom
# ══════════════════════════════════════════════════════════════════════════════

final_plot <- wrap_elements(full = row1) /
              wrap_elements(full = row2) /
              wrap_elements(full = row3) /
              legend_strip +
  plot_layout(heights = c(1, 1, 1, 0.06))

# ══════════════════════════════════════════════════════════════════════════════
# EXPORT
# ══════════════════════════════════════════════════════════════════════════════

total_n    <- length(order_complete) + length(order_single) + length(order_none)
fig_width  <- max(30, total_n * 0.35)
fig_height <- 31   # +1" to accommodate the legend strip

svglite::svglite("genome_browser_supplement_all78.svg",
                 width = fig_width, height = fig_height)
print(final_plot)
dev.off()

ggsave("genome_browser_supplement_all78.png",
       final_plot,
       device = ragg::agg_png,
       width  = fig_width,
       height = fig_height,
       dpi    = 300,
       limitsize = FALSE)

message("Done → genome_browser_supplement_all78.svg / .png")
message("Figure size: ", fig_width, " x ", fig_height, " inches")

