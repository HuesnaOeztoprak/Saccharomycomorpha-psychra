# =============================================================================
# Transcriptome-Wide Expression Analysis ✳️
#
# Sections:
#   0.  Load libraries
#   1.  Load data
#   2.  Clean IDs
#   3.  Label functional categories (transcript level)
#   4.  Collapse to gene level (max-TPM isoform per gene)
#   5.  Log-transform TPM
#   6.  Summary statistics
#   7.  Proportion of total transcriptional output
#   8.  Pairwise Wilcoxon + effect sizes
#   9.  Fisher enrichment (top-10% expressed genes)
#  10.  Colour palette
#  11.  PANEL A — Annotation-based overview figure (Venn + bar)
#  12.  PANEL B — Global expression violin/boxplot
#  13.  CAZyme family annotation
#  14.  PANEL C — Top CAZyme families (expression)
#  15.  PFAM-rescue of unknown TC families
#  16.  TC class harmonisation (consistent subclass level, remove PFAM suffix)
#  17.  Transporter family analysis
#  18.  PANEL D — Top transporter families (expression)
#  19.  Supplementary: per-family CAZyme violin
#  20.  Supplementary: expression breadth stacked bar
#  21.  Supplementary tables
#  22.  Length-bias check
#  23.  Expression breadth stats
#  24.  Combined main figure (patchwork)
#
# =============================================================================
# 0. LOAD LIBRARIES

library(ggplot2)      # plotting engine
library(dplyr)        # data manipulation
library(rstatix)      # Wilcoxon tests + effect sizes
library(ggpubr)       # stat_compare_means for significance brackets
library(patchwork)    # multi-panel figures
library(pheatmap)     # heatmaps
library(ggvenn)       # Venn diagrams

# =============================================================================
# 1. LOAD DATA

setwd("~/Desktop/Postdoc/Saccaromycomorpha/expression")
quant    <- read.table("quant.sf",                      header = TRUE, sep = "\t")
secreted <- read.table("final_secreted_Protein_IDs.txt", header = TRUE, stringsAsFactors = FALSE)
cazy     <- read.table("CAZymes_high_confidence.tsv",   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
trans    <- read.table("transporter_final.tsv",         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(secreted)[1] <- "Name"
colnames(cazy)[1]     <- "Name"
colnames(trans)[1]    <- "Name"

# Genome-level annotation counts (ground truth from annotation pipeline).
# These are used in Panel A (annotation overview). 
GENOME_total_genes        <- 18680
GENOME_secreted           <- 1015
GENOME_cazymes            <- 109
GENOME_secreted_cazymes   <- 64
GENOME_transporters       <- 724

# =============================================================================
# 2. CLEAN IDs

quant$Name    <- trimws(quant$Name)
secreted$Name <- trimws(secreted$Name)

# Gene-level IDs: strip isoform suffix (g79.t2 → g79).
# Gene-level matching for secreted automatically recovers cases where secretion
# was predicted on .t1 but only .t2 is present in the transcriptome.
quant$Gene    <- sub("\\.t[0-9]+$", "", quant$Name)
secreted$Gene <- sub("\\.t[0-9]+$", "", secreted$Name)

# Diagnostic
n_sec_total   <- length(unique(secreted$Gene))
n_sec_inQuant <- length(intersect(unique(secreted$Gene), unique(quant$Gene)))
cat(sprintf("Secreted predicted: %d | In transcriptome: %d | Absent: %d (%.1f%%)\n\n",
    n_sec_total, n_sec_inQuant, n_sec_total - n_sec_inQuant,
    100 * (n_sec_total - n_sec_inQuant) / n_sec_total))


# =============================================================================
# 3. LABEL FUNCTIONAL CATEGORIES (transcript level)

quant$Secreted    <- quant$Gene %in% secreted$Gene
quant$CAZyme      <- quant$Name %in% cazy$Name
quant$Transporter <- quant$Name %in% trans$Name

# Hierarchical priority: Secreted CAZyme > CAZyme > Transporter > Secreted > background.
# Note: 2 genes with Secreted=TRUE are overwritten to 'Transporter' — they retain
# Secreted=TRUE boolean flag and are counted in overlap diagnostics.
quant$Category <- "background"
quant$Category[quant$Secreted]                 <- "Secreted"
quant$Category[quant$CAZyme]                   <- "CAZyme"
quant$Category[quant$Secreted & quant$CAZyme]  <- "Secreted CAZyme"
quant$Category[quant$Transporter]              <- "Transporter"

cat("Transcript counts per category (pre-collapse):\n")
print(table(quant$Category))
cat("\n")


# =============================================================================
# 4. COLLAPSE TO GENE LEVEL — MAX-TPM ISOFORM

# Note: 712 genes have >1 isoform (792 extra transcripts). Retain max-TPM isoform
# per gene to avoid double-counting in statistical comparisons.

quant_gene <- quant %>%
  group_by(Gene) %>%
  slice_max(TPM, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Gene-level counts (post-collapse):\n")
print(table(quant_gene$Category))
cat("Total unique genes:", nrow(quant_gene), "\n\n")


# =============================================================================
# 5. LOG-TRANSFORM

# log10(TPM+1): +1 pseudocount keeps zero-expressed genes finite.
# Use this column directly on plot axes — do NOT use scale_y_log10()
# which silently drops TPM=0 rows and generates warnings.

quant_gene$logTPM <- log10(quant_gene$TPM + 1)


# =============================================================================
# 6. SUMMARY STATISTICS

summary_stats <- quant_gene %>%
  group_by(Category) %>%
  summarise(
    n             = n(),
    n_zero        = sum(TPM == 0),
    pct_zero      = round(100 * mean(TPM == 0), 1),
    mean_logTPM   = round(mean(logTPM), 3),
    median_logTPM = round(median(logTPM), 3),
    mean_TPM      = round(mean(TPM), 1),
    median_TPM    = round(median(TPM), 1)
  )

cat("Expression summary by category:\n")
print(as.data.frame(summary_stats))
write.table(summary_stats, "summary_stats_by_category.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")


# =============================================================================
# 7. PROPORTION OF TOTAL TRANSCRIPTIONAL OUTPUT

total_TPM <- sum(quant_gene$TPM)

tpo <- quant_gene %>%
  group_by(Category) %>%
  summarise(
    n_genes   = n(),
    sum_TPM   = round(sum(TPM), 1),
    pct_genes = round(100 * n() / nrow(quant_gene), 2),
    pct_TPM   = round(100 * sum(TPM) / total_TPM, 2)
  ) %>%
  arrange(desc(pct_TPM))

cat("Proportion of total transcriptional output:\n")
print(as.data.frame(tpo))
write.table(tpo, "transcriptional_output_proportions.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")


# =============================================================================
# 8. PAIRWISE WILCOXON + EFFECT SIZES

# Non-parametric test chosen for: right-skewed TPM, zero-inflation (~18%),
# unequal variances. BH correction for multiple comparisons.
# Effect size r (rank-biserial): |r|~0.1 small, 0.3 medium, 0.5 large.

pairwise_stats <- quant_gene %>%
  pairwise_wilcox_test(logTPM ~ Category, p.adjust.method = "BH") %>%
  add_significance("p.adj")

pairwise_effects <- quant_gene %>%
  wilcox_effsize(logTPM ~ Category, paired = FALSE)

pairwise_full <- pairwise_stats %>%
  left_join(pairwise_effects %>% select(group1, group2, effsize, magnitude),
            by = c("group1", "group2"))

cat("Pairwise Wilcoxon (BH-corrected) with effect sizes:\n")
print(as.data.frame(
  pairwise_full %>% select(group1, group2, n1, n2, p, p.adj, p.adj.signif, effsize, magnitude)
))
write.table(pairwise_full, "pairwise_wilcox_effectsizes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")


# =============================================================================
# 9. FISHER ENRICHMENT — TOP-10% EXPRESSED GENES

# Gene-level 90th percentile TPM threshold, propagated to all enrichment tests
# and to the CAZyme family bar chart.

threshold_90 <- quantile(quant_gene$TPM, 0.90)
cat(sprintf("Gene-level TPM 90th percentile: %.2f\n\n", threshold_90))

quant_gene$HighExpr <- quant_gene$TPM > threshold_90

for (flag in c("CAZyme", "Secreted", "Transporter")) {
  ft <- fisher.test(table(quant_gene[[flag]], quant_gene$HighExpr))
  cat(sprintf("Fisher — %s in top 10%%: OR=%.3f [%.3f-%.3f], p=%.3e\n",
      flag, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
}
cat("\n")


# =============================================================================
# 10. COLOUR PALETTE AND FACTOR ORDER

category_colors <- c(
  "background"      = "#6F6F6F",
  "Secreted"        = "#BDBDBD",
  "CAZyme"          = "#6A0DAD",
  "Secreted CAZyme" = "#C8A2E8",
  "Transporter"     = "#5FA8A8"
)

category_order <- c("background", "Secreted", "CAZyme", "Secreted CAZyme", "Transporter")
quant_gene$Category <- factor(quant_gene$Category, levels = category_order)


# =============================================================================
# 11. PANEL A — ANNOTATION-BASED OVERVIEW (Venn + bar)

# This panel uses GENOME ANNOTATION COUNTS (not expression data).

# Build sets from genome annotation counts using integer vectors as proxies
# (ggvenn works on named lists of sets — we simulate from known counts)
secreted_ids  <- 1:GENOME_secreted                   # 1–1015
cazyme_ids    <- (GENOME_secreted - GENOME_secreted_cazymes + 1):
                  (GENOME_secreted - GENOME_secreted_cazymes + GENOME_cazymes)
# cazyme_ids spans: non-secreted CAZymes + secreted CAZymes (overlap region)
# Overlap: last 64 IDs of secreted_ids = first 64 IDs of cazyme_ids
overlap_start <- GENOME_secreted - GENOME_secreted_cazymes + 1
cazyme_ids    <- c(overlap_start:(overlap_start + GENOME_cazymes - 1))
# Ensure exactly 64 overlap
cat(sprintf("Venn check — overlap: %d (should be %d)\n",
    length(intersect(secreted_ids, cazyme_ids)), GENOME_secreted_cazymes))

venn_annot <- list(
  "Secreted proteins\n(n = 1,015)" = secreted_ids,
  "CAZymes\n(n = 109)"             = cazyme_ids
)

pVenn_annot <- ggvenn(
  venn_annot,
  fill_color    = c("#BDBDBD", "#6A0DAD"),
  fill_alpha    = 0.55,
  stroke_size   = 0.5,
  set_name_size = 3.8,
  text_size     = 3.8
) +
  # Annotation: add transporter count as text below the circles
  annotate("text", x = 0, y = -2.5, size = 4, fontface = "bold",
           label = sprintf("Transporters: %d (shown separately)", GENOME_transporters)) +
  annotate("text", x = 0, y = -3.2, size = 3.2, color = "grey40",
           label = sprintf("Total genes: %d", GENOME_total_genes)) +
  ggtitle("Predicted Proteome Composition") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.margin = margin(5, 5, 25, 5))

# Companion bar: gene counts per functional group (annotation-level)
annot_bar_df <- data.frame(
  Category = factor(
    c("background", "Secreted", "CAZyme", "Secreted\nCAZyme", "Transporter"),
    levels = c("background", "Secreted", "CAZyme", "Secreted\nCAZyme", "Transporter")
  ),
  n = c(
    GENOME_total_genes - GENOME_secreted - (GENOME_cazymes - GENOME_secreted_cazymes) - GENOME_transporters,
    GENOME_secreted - GENOME_secreted_cazymes,   # secreted only
    GENOME_cazymes - GENOME_secreted_cazymes,    # CAZyme only
    GENOME_secreted_cazymes,
    GENOME_transporters
  ),
  fill_col = c("#6F6F6F", "#BDBDBD", "#6A0DAD", "#C8A2E8", "#5FA8A8")
)

pAnnotBar <- ggplot(annot_bar_df, aes(x = Category, y = n, fill = Category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = setNames(annot_bar_df$fill_col, annot_bar_df$Category)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic(base_size = 12) +
  ylab("Number of genes") +
  xlab("") +
  ggtitle("Functional annotation (all predicted genes)") +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# Combine Venn + bar side by side as Panel A
pA_overview <- pVenn_annot | pAnnotBar
ggsave("PanelA_annotation_overview.pdf", pA_overview, width = 11, height = 5)


# =============================================================================
# 12. PANEL B — GLOBAL EXPRESSION VIOLIN + BOXPLOT

# Uses quant.sf numbers
# Significance brackets: each category vs. background only.

comparisons_vs_bg <- list(
  c("background", "Secreted"),
  c("background", "CAZyme"),
  c("background", "Secreted CAZyme"),
  c("background", "Transporter")
)

pB <- ggplot(quant_gene, aes(x = Category, y = logTPM, fill = Category)) +
  geom_violin(trim = FALSE, alpha = 0.65) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = category_colors) +
  stat_compare_means(comparisons = comparisons_vs_bg, method = "wilcox.test",
                     label = "p.signif", step.increase = 0.07) +
  theme_classic(base_size = 13) +
  ylab("log10(TPM + 1)") +
  xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("PanelB_global_expression.pdf", pB, width = 7, height = 5.5)


# =============================================================================
# 13. CAZYME FAMILY ANNOTATION
# =============================================================================
# Extract family from dbCAN2 output. Priority: HMMER > DIAMOND > Hotpep.

cazy_raw        <- read.table("CAZymes_high_confidence.tsv", header = TRUE,
                               sep = "\t", stringsAsFactors = FALSE)
cazy_raw$gene_id <- sub("\\.t[0-9]+$", "", cazy_raw$Gene.ID)

fam_hmmer   <- ifelse(cazy_raw$HMMER   != "N", sub("\\(.*", "", cazy_raw$HMMER),   NA)
fam_diamond <- ifelse(cazy_raw$DIAMOND != "N", sub(".*\\+",  "", cazy_raw$DIAMOND), NA)
fam_hotpep  <- ifelse(cazy_raw$Hotpep  != "N", sub("\\(.*",  "", cazy_raw$Hotpep),  NA)

cazy_raw$family <- dplyr::coalesce(fam_hmmer, fam_diamond, fam_hotpep)
cazy_clean <- cazy_raw %>% filter(!is.na(family)) %>% distinct(gene_id, family)

cat(sprintf("CAZyme genes with family: %d | Unique families: %d\n\n",
    nrow(cazy_clean), length(unique(cazy_clean$family))))

quant_gene$gene_id <- sub("\\.t[0-9]+$", "", quant_gene$Name)
expr_cazy <- quant_gene %>% inner_join(cazy_clean, by = "gene_id")

# Mark high-expression using the gene-level 90th percentile (consistent with Fisher test)
expr_cazy$HighExpr_gene <- expr_cazy$TPM > threshold_90

family_expression <- expr_cazy %>%
  group_by(family) %>%
  summarise(
    n_genes    = n(),
    total_TPM  = sum(TPM),
    mean_TPM   = mean(TPM),
    median_TPM = median(TPM),
    HighFamily = any(HighExpr_gene)     # TRUE if any member gene in top 10%
  ) %>%
  arrange(desc(total_TPM))

write.table(family_expression, "CAZyme_family_expression.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# =============================================================================
# 14. PANEL C — TOP CAZYME FAMILIES (expression)

# Show all families (47 total)
# Colour: purple = has at least one gene in top-10% expressed (gene-level threshold)

n_families_show <- min(nrow(family_expression), 47)
top_cazy_plot   <- family_expression %>% slice(1:n_families_show)

pC <- ggplot(top_cazy_plot,
             aes(x = reorder(family, total_TPM), y = total_TPM, fill = HighFamily)) +
  geom_col() +
  scale_fill_manual(
    values = c("FALSE" = "#C8A2E8", "TRUE" = "#6A0DAD"),
    labels = c("FALSE" = "No gene in top 10%", "TRUE" = "Has gene in top 10%"),
    name = ""
  ) +
  coord_flip() +
  theme_classic(base_size = 11) +
  ylab("Total TPM") +
  xlab("CAZyme family") +
  theme(legend.position = "bottom")

ggsave("PanelC_CAZyme_families_expression.pdf", pC, width = 7, height = 8)


# =============================================================================
# 15. PFAM-RESCUE OF UNKNOWN TC FAMILIES

# Curated Pfam → TC subclass lookup (TCDB.org; Saier et al. 2021).
# ~28% of eukaryotic transporters remain unresolvable after PFAM rescue - this is expected in divergent lineages and not an annotation quality issue.

# Curated Pfam → TC subclass lookup
pfam_to_tc <- c(
  "MFS_1"           = "2.A.1",
  "MFS_2"           = "2.A.1",
  "Sugar_tr"        = "2.A.1",
  "Lactate_perm"    = "2.A.1",
  "MFS_3"           = "2.A.1",
  "Aa_trans"        = "2.A.3",
  "Trp_Tyr_perm"    = "2.A.3",
  "Mito_carr"       = "2.A.29",
  "Nuc_H_symport"   = "2.A.39",
  "ABC_tran"        = "3.A.1",
  "ABC_membrane"    = "3.A.1",
  "ABC_membrane_2"  = "3.A.1",
  "TOBE_2"          = "3.A.1",
  "Ammonium_transp" = "1.A.11",
  "Form_Nir_trans"  = "1.A.16",
  "MIP"             = "1.A.8",
  "Ion_trans"       = "1.A.1",
  "E1-E2_ATPase"    = "3.A.3",
  "Cation_ATPase"   = "3.A.3",
  "HAD"             = "3.A.3"
)

rescue_tc_from_pfam <- function(pfam_string, lookup = pfam_to_tc) {
  if (is.na(pfam_string) || pfam_string %in% c("-", "")) return(NA_character_)
  domains <- trimws(unlist(strsplit(pfam_string, ",")))
  for (d in domains) {
    if (d %in% names(lookup)) return(lookup[[d]])
  }
  return(NA_character_)
}

trans$TC_rescued <- sapply(trans$PFAM, rescue_tc_from_pfam)

# Combine eggNOG + PFAM
trans$TC_resolved <- dplyr::case_when(
  trans$TC_family != "unknown" ~ trans$TC_family,
  !is.na(trans$TC_rescued)     ~ trans$TC_rescued,
  TRUE                         ~ "unknown"
)

# Track annotation source
trans$TC_source <- dplyr::case_when(
  trans$TC_family != "unknown" ~ "eggNOG",
  !is.na(trans$TC_rescued)     ~ "PFAM",
  TRUE                         ~ "unresolved"
)

cat("TC source summary:\n")
print(table(trans$TC_source))
cat("\n")


# =============================================================================
# 16. TC CLASS HARMONISATION

# Goal: consistent subclass level throughout (X.A.XX format where known).
# Extract TC levels
get_tc_levels <- function(tc) {
  if (is.na(tc) || tc == "unknown") {
    return(list(class="unknown", subclass="unknown", family="unknown", resolution="unknown"))
  }
  
  parts <- unlist(strsplit(tc, "\\."))
  n <- length(parts)
  
  class <- ifelse(n >= 2, paste(parts[1:2], collapse="."), tc)
  subclass <- ifelse(n >= 3, paste(parts[1:3], collapse="."), class)
  family <- ifelse(n >= 4, paste(parts[1:4], collapse="."), subclass)
  
  resolution <- ifelse(n >= 4, "family",
                  ifelse(n == 3, "subclass",
                  ifelse(n == 2, "class", "other")))
  
  return(list(class=class, subclass=subclass, family=family, resolution=resolution))
}

tc_parsed <- lapply(trans$TC_resolved, get_tc_levels)

trans$TC_class     <- sapply(tc_parsed, `[[`, "class")
trans$TC_subclass  <- sapply(tc_parsed, `[[`, "subclass")
trans$TC_family_lv <- sapply(tc_parsed, `[[`, "family")
trans$TC_resolution<- sapply(tc_parsed, `[[`, "resolution")

cat("TC resolution levels:\n")
print(table(trans$TC_resolution))
cat("\n")


# =============================================================================
# 17. TRANSPORTER FAMILY ANALYSIS

expr_trans <- quant_gene %>%
  filter(Name %in% trans$Name) %>%
  left_join(trans, by = "Name")

# =========================
# MAIN TABLE (subclass level)

trans_subclass <- expr_trans %>%
  group_by(TC_subclass) %>%
  summarise(
    n_genes    = n(),
    total_TPM  = sum(TPM),
    mean_TPM   = mean(TPM),
    pct_eggnog = round(100 * mean(TC_source == "eggNOG"), 1),
    pct_highres= round(100 * mean(TC_resolution %in% c("subclass","family")), 1)
  ) %>%
  arrange(desc(total_TPM))

write.table(trans_subclass, "Transporter_subclass_expression.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)


# TOP 20 subclass for plotting
top20_trans <- trans_subclass %>% slice(1:20)

# Stats: transporter vs background
wilcox_trans  <- wilcox.test(logTPM ~ Transporter, data = quant_gene)
effsize_trans <- quant_gene %>% wilcox_effsize(logTPM ~ Transporter)

cat(sprintf("Transporter vs background: W=%.0f, p=%.2e, r=%.3f (%s)\n\n",
    wilcox_trans$statistic,
    wilcox_trans$p.value,
    effsize_trans$effsize,
    effsize_trans$magnitude))
# =============================================================================
# 18. PANEL D — TOP TRANSPORTER FAMILIES (expression)

# Labels are harmonised TC subclass names (no PFAM-rescue suffix in the figure).

pD <- ggplot(top20_trans,
             aes(x = reorder(TC_subclass, total_TPM), y = total_TPM)) +
  geom_col(fill = "#5FA8A8") +
  coord_flip() +
  theme_classic(base_size = 11) +
  ylab("Total TPM") +
  xlab("Transporter subclass (TC)") +
  ggtitle("Top transporter subclasses (expression)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("PanelD_Transporter_subclasses.pdf", pD, width = 7, height = 6)

ggplot(trans_subclass,
       aes(x = reorder(TC_subclass, mean_TPM), y = mean_TPM)) +
  geom_col(fill = "#8888AA") +
  coord_flip() +
  theme_classic() +
  ylab("Mean TPM per gene")
  
#panel D  
# =============================================================================
# PANEL D — Transporter subclasses (bubble plot)

# ---------------------------
# Ordering (TC hierarchy)
# ---------------------------
trans_subclass <- trans_subclass %>%
  mutate(
    TC_class_num = case_when(
      TC_subclass == "unknown" ~ 99,
      TRUE ~ suppressWarnings(
        as.numeric(sub("^([0-9]+).*", "\\1", TC_subclass))
      )
    )
  ) %>%
  arrange(TC_class_num, TC_subclass) %>%
  mutate(
    TC_subclass = factor(TC_subclass, levels = unique(TC_subclass))
  )

# ---------------------------
# Plot
D <- ggplot(trans_subclass,
       aes(x = TC_subclass,
           y = mean_TPM,
           size = n_genes,
           color = total_TPM)) +
  geom_point(alpha = 0.95) +

  # Size = gene count
  scale_size_continuous(name = "Gene count") +

  # Saturated teal gradient
  scale_color_gradient(
    low  = "#5FA8A8",
    high = "#0F4C5C",
    name = "Total TPM"
  ) +

  coord_flip() +
  theme_classic(base_size = 12) +

  ylab("Mean TPM per gene") +
  xlab("Transporter subclass (TC)") +

  theme(
    legend.position = "right"
  )

# ---------------------------
# Save
ggsave("PanelD_Transporter_subclasses_BLOB.pdf", D, width = 7, height = 6)


## log transformerd version
  ggplot(trans_subclass,
       aes(x = reorder(TC_subclass, mean_TPM),
           y = mean_TPM,
           size = n_genes,
           color = total_TPM)) +
  geom_point(alpha = 0.9) +

  # Size = gene count
  scale_size_continuous(name = "Gene count") +

  # Color = single-hue gradient (light → dark teal)
  scale_color_gradient(
    low  = "#D6ECEC",   # very light teal
    high = "#5FA8A8",   # your transporter color
    name = "Total TPM"
  ) +

  # Optional: log scale
  scale_y_log10() +

  coord_flip() +
  theme_classic(base_size = 12) +

  ylab("Mean TPM per gene (log scale)") +
  xlab("Transporter subclass (TC)") +

  theme(
    legend.position = "right"
  )
  ##
  
# =============================================================================
# 19. SUPPLEMENTARY FIGURE S1 — PER-FAMILY CAZYME VIOLIN

pSup1 <- ggplot(expr_cazy,
                aes(x = reorder(family, TPM, FUN = median), y = logTPM, fill = family)) +
  geom_violin(trim = FALSE, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
               alpha = 0.6, show.legend = FALSE) +
  ylab("log10(TPM + 1)") + xlab("CAZyme family") +
  coord_flip() +
  theme_classic(base_size = 10) +
  ggtitle("Expression distribution per CAZyme family (all families)")

ggsave("SupFigS1_CAZyme_family_violins.pdf", pSup1, width = 8, height = 12)


# =============================================================================
# 20. SUPPLEMENTARY FIGURE S2 — EXPRESSION BREADTH STACKED BAR

# Quantifies the bimodal pattern in secreted proteins: many lowly/not expressed,
# a subset highly expressed (the osmotrophy "workhorses").
# Threshold for 'high' = gene-level 90th percentile (consistent throughout).

breadth_df <- quant_gene %>%
  mutate(ExprClass = case_when(
    TPM == 0                      ~ "Not expressed",
    TPM > 0 & TPM <= 1            ~ "Low (0 < TPM <= 1)",
    TPM > 1 & TPM <= threshold_90 ~ "Moderate",
    TPM > threshold_90            ~ "High (top 10%)"
  )) %>%
  mutate(ExprClass = factor(ExprClass,
    levels = c("Not expressed", "Low (0 < TPM <= 1)", "Moderate", "High (top 10%)")))

breadth_summary <- breadth_df %>%
  group_by(Category, ExprClass) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(pct = 100 * n / sum(n))

pSup2 <- ggplot(breadth_summary, aes(x = Category, y = pct, fill = ExprClass)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(
    "Not expressed"       = "#D9D9D9",
    "Low (0 < TPM <= 1)"  = "#A8D8B9",
    "Moderate"            = "#5FA8A8",
    "High (top 10%)"      = "#6A0DAD"
  ), name = "Expression class") +
  theme_classic(base_size = 12) +
  ylab("Percentage of genes (%)") + xlab("") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("SupFigS2_expression_breadth.pdf", pSup2, width = 7, height = 5)


# =============================================================================
# 21. SUPPLEMENTARY TABLES

# S1: Top 20 expressed CAZyme genes
write.table(
  expr_cazy %>% arrange(desc(TPM)) %>%
    select(Name, Gene, family, Secreted, Category, TPM, logTPM) %>% head(20),
  "SupTable_S1_Top20_CAZyme_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# S2: Top 20 expressed transporter genes
write.table(
  expr_trans %>% arrange(desc(TPM)) %>%
    select(Name, TC_harmonised, TC_source, confidence, TPM, logTPM) %>% head(20),
  "SupTable_S2_Top20_Transporter_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# S3: Top 20 expressed secreted proteins (all)
write.table(
  quant_gene %>% filter(Secreted) %>% arrange(desc(TPM)) %>%
    select(Name, Gene, Category, TPM, logTPM) %>% head(20),
  "SupTable_S3_Top20_Secreted_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Top 10 expressed CAZyme genes:\n")
print(as.data.frame(
  expr_cazy %>% arrange(desc(TPM)) %>%
    select(Name, family, Secreted, Category, TPM) %>% head(10)
))
cat("\nTop 10 expressed secreted proteins:\n")
print(as.data.frame(
  quant_gene %>% filter(Secreted) %>% arrange(desc(TPM)) %>%
    select(Name, Category, TPM) %>% head(10)
))
cat("\n")


# =============================================================================
# 22. LENGTH-BIAS CHECK

# TPM corrects for length — Spearman rho should be negligible (<0.1).

len_cor <- cor.test(quant_gene$Length, quant_gene$logTPM, method = "spearman")
cat(sprintf("Length-bias check: rho=%.4f, p=%.3e\n", len_cor$estimate, len_cor$p.value))
cat(ifelse(abs(len_cor$estimate) < 0.1,
    "Interpretation: negligible — no systematic length bias.\n\n",
    "WARNING: non-negligible correlation. Review TPM normalisation.\n\n"))


# =============================================================================
# 23. EXPRESSION BREADTH STATS (for results text)


sec_genes   <- quant_gene %>% filter(Secreted)
n_sec_total <- nrow(sec_genes)

cat("Expression breadth of secreted proteins:\n")
cat(sprintf("  Not expressed (TPM=0):         %d (%.1f%%)\n",
    sum(sec_genes$TPM == 0),         100*mean(sec_genes$TPM == 0)))
cat(sprintf("  Low (0 < TPM <= 1):            %d (%.1f%%)\n",
    sum(sec_genes$TPM > 0 & sec_genes$TPM <= 1), 100*mean(sec_genes$TPM > 0 & sec_genes$TPM <= 1)))
cat(sprintf("  Moderate (1 to 90th pctile):   %d (%.1f%%)\n",
    sum(sec_genes$TPM > 1 & sec_genes$TPM <= threshold_90),
    100*mean(sec_genes$TPM > 1 & sec_genes$TPM <= threshold_90)))
cat(sprintf("  Highly expressed (top 10%%):   %d (%.1f%%)\n\n",
    sum(sec_genes$TPM > threshold_90), 100*mean(sec_genes$TPM > threshold_90)))


# =============================================================================
# 24. COMBINED MAIN FIGURE

# Panel A: Annotation overview (Venn + gene count bar)
# Panel B: Expression violin/boxplot
# Panel C: CAZyme family expression bar
# Panel D: Transporter family expression bar

main_figure <- pA_overview / pB / (pC | pD) +
  plot_layout(heights = c(0.9, 1.2, 1.4)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 15))
  )

ggsave("Figure1_combined.pdf", main_figure, width = 12, height = 20)
ggsave("Figure1_combined.png", main_figure, width = 12, height = 20, dpi = 300)


# =============================================================================
# FILE MANIFEST, uff finally done :)

cat("\n=== OUTPUT FILES ===\n")
cat("Main figure:    Figure1_combined.pdf/png\n")
cat("Panels:         PanelA_annotation_overview.pdf\n")
cat("                PanelB_global_expression.pdf\n")
cat("                PanelC_CAZyme_families_expression.pdf\n")
cat("                PanelD_Transporter_top20.pdf\n")
cat("Supp figures:   SupFigS1_CAZyme_family_violins.pdf\n")
cat("                SupFigS2_expression_breadth.pdf\n")
cat("Stats:          pairwise_wilcox_effectsizes.tsv\n")
cat("                summary_stats_by_category.tsv\n")
cat("                transcriptional_output_proportions.tsv\n")
cat("                CAZyme_family_expression.tsv\n")
cat("                Transporter_family_expression.tsv\n")
cat("Supp tables:    SupTable_S1_Top20_CAZyme_genes.tsv\n")
cat("                SupTable_S2_Top20_Transporter_genes.tsv\n")
cat("                SupTable_S3_Top20_Secreted_genes.tsv\n")


