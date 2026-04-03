# ============================================
# Crohn's Disease Metagenomics Analysis
# Carlos Simon Foundation Technical Assignment
# Author: Larisa Atanasiu
# Date: April 2026
# Data: HMP2 study (Lloyd-Price et al., Nature 2019)
# ============================================

library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(Maaslin2)

# ============================================
# GLOBAL SETTINGS
# ============================================

# Color palette — consistent across all figures
COLORS <- c("Control" = "#2166AC", "Crohn" = "#D6604D")

# Publication quality theme
theme_publication <- function() {
  theme_bw() +
    theme(
      text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

# Output directory
dir.create("~/crohns_mgx/figures", showWarnings = FALSE)

# ============================================
# LOAD DATA
# ============================================

abund <- read.table("~/crohns_mgx/taxonomy/species_abundance.txt",
                    header = FALSE, sep = "\t", quote = "")

colnames(abund) <- c("clade_name", "NCBI_tax_id",
                     "SRR5983484", "SRR5983464", "SRR5983451",
                     "SRR5983438", "SRR5983306", "SRR5983287")

# Extract clean species names
abund$species <- gsub(".*s__", "", abund$clade_name)
abund$species <- gsub("_", " ", abund$species)
rownames(abund) <- abund$species

# Keep only abundance columns
abund_mat <- abund[, c("SRR5983484", "SRR5983464", "SRR5983451",
                       "SRR5983438", "SRR5983306", "SRR5983287")]

cat("=== Data Summary ===\n")
cat("Species detected:", nrow(abund_mat), "\n")
cat("Samples:", ncol(abund_mat), "\n\n")

# ============================================
# METADATA
# Healthy controls: SRR5983438, SRR5983451, SRR5983484
# Crohn's disease:  SRR5983306, SRR5983464, SRR5983287
# Source: Assignment specification (HMP2 study)
# QC flags: SRR5983484 = low sequencing depth (0.5M reads)
#           SRR5983464 = high host contamination (35% human reads)
# ============================================

metadata <- data.frame(
  sample = c(
    # Healthy controls
    "SRR5983438", "SRR5983451", "SRR5983484",
    # Crohn's disease
    "SRR5983306", "SRR5983464", "SRR5983287"
  ),
  condition = c(
    # Healthy controls
    "Control", "Control", "Control",
    # Crohn's disease
    "Crohn", "Crohn", "Crohn"
  ),
  note = c(
    # Healthy controls flags
    "", "", "LOW_DEPTH",
    # Crohn's disease flags
    "", "HIGH_HOST", ""
  )
)
rownames(metadata) <- metadata$sample

# Create short sample labels for plots
metadata$label <- paste0(
  ifelse(metadata$condition == "Control", "HC", "CD"),
  c("1", "2", "3*", "1", "2*", "3")
)

cat("=== Sample Metadata ===\n")
print(metadata)
cat("\n* = flagged sample\n\n")

# ============================================
# ALPHA DIVERSITY
# Higher diversity generally associated with gut health
# Crohn's disease typically shows reduced diversity
# ============================================

shannon <- diversity(t(abund_mat), index = "shannon")
richness <- specnumber(t(abund_mat))

alpha_df <- data.frame(
  sample = names(shannon),
  shannon = shannon,
  richness = richness
)
alpha_df <- merge(alpha_df, metadata, by = "sample")

cat("=== Alpha Diversity ===\n")
print(alpha_df[, c("sample", "condition", "shannon", "richness", "note")])
cat("\n")

# Shannon diversity plot
p_shannon <- ggplot(alpha_df, aes(x = condition, y = shannon,
                                  fill = condition, color = condition)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 4, alpha = 0.8) +
  geom_text(aes(label = label), nudge_x = 0.28, size = 3.5, color = "black") +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  theme_publication() +
  labs(title = "Alpha Diversity",
       subtitle = "Shannon Diversity Index",
       x = NULL,
       y = "Shannon Index") +
  stat_compare_means(method = "wilcox.test", size = 4) +
  theme(legend.position = "none")

# Species richness plot
p_richness <- ggplot(alpha_df, aes(x = condition, y = richness,
                                   fill = condition, color = condition)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 4, alpha = 0.8) +
  geom_text(aes(label = label), nudge_x = 0.28, size = 3.5, color = "black") +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  theme_publication() +
  labs(title = "Alpha Diversity",
       subtitle = "Species Richness",
       x = NULL,
       y = "Number of Species") +
  stat_compare_means(method = "wilcox.test", size = 4) +
  theme(legend.position = "none")

# Combine into one figure
p_alpha_combined <- ggarrange(p_shannon, p_richness,
                              ncol = 2, nrow = 1,
                              common.legend = TRUE)

print(p_alpha_combined)
ggsave("~/crohns_mgx/figures/alpha_diversity.pdf",
       p_alpha_combined, width = 10, height = 5)
cat("Alpha diversity saved\n\n")

# ============================================
# BETA DIVERSITY
# Bray-Curtis dissimilarity captures differences in
# both species presence and relative abundance
# ============================================

bc_dist <- vegdist(t(abund_mat), method = "bray")

pcoa <- cmdscale(bc_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  sample = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)
pcoa_df <- merge(pcoa_df, metadata, by = "sample")

var_explained <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

perm <- adonis2(bc_dist ~ condition,
                data = metadata[rownames(as.matrix(bc_dist)), ],
                permutations = 999)

cat("=== PERMANOVA Results ===\n")
print(perm)
cat("\n")

p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2,
                              color = condition,
                              fill = condition,
                              shape = condition,
                              label = label)) +
  geom_point(size = 6, alpha = 0.8) +
  geom_text(nudge_y = 0.025, size = 3.5, color = "black") +
  stat_ellipse(aes(fill = condition), geom = "polygon",
               alpha = 0.1, level = 0.95, linetype = "dashed") +
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  scale_shape_manual(values = c("Control" = 21, "Crohn" = 24)) +
  theme_publication() +
  labs(title = "Beta Diversity",
       subtitle = paste0("Bray-Curtis PCoA | PERMANOVA: R² = ",
                         round(perm$R2[1], 3),
                         ", p = ", perm$`Pr(>F)`[1]),
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)"),
       color = "Condition", fill = "Condition", shape = "Condition")

print(p_beta)
ggsave("~/crohns_mgx/figures/beta_pcoa.pdf", p_beta, width = 8, height = 6)
cat("Beta diversity saved\n\n")

# ============================================
# TAXONOMIC COMPOSITION - TOP TAXA BARPLOT
# ============================================

mean_abund <- rowMeans(abund_mat)
top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:15])

plot_data <- as.data.frame(t(abund_mat[top_taxa, ]))
plot_data$sample <- rownames(plot_data)

plot_long <- melt(plot_data, id.vars = "sample",
                  variable.name = "taxon",
                  value.name = "abundance")
plot_long <- merge(plot_long, metadata, by = "sample")

# Use a nice color palette for taxa
taxa_colors <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Set2"))

p_taxa <- ggplot(plot_long, aes(x = label, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~condition, scales = "free_x") +
  scale_fill_manual(values = taxa_colors) +
  theme_publication() +
  labs(title = "Microbial Community Composition",
       subtitle = "Top 15 species by mean relative abundance",
       x = "Sample",
       y = "Relative Abundance (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 8, face = "italic")) +
  guides(fill = guide_legend(ncol = 1))

print(p_taxa)
ggsave("~/crohns_mgx/figures/top_taxa.pdf", p_taxa, width = 14, height = 7)
cat("Top taxa saved\n\n")

# ============================================
# FAECALIBACTERIUM PRAUSNITZII
# Key butyrate-producing anti-inflammatory commensal
# Consistently depleted in IBD (Sokol et al. 2008)
# Produces butyrate via SCFA pathway
# Has direct anti-inflammatory effects on epithelium
# ============================================

faec_rows <- grep("Faecalibacterium", rownames(abund_mat), ignore.case = TRUE)
cat("=== Faecalibacterium species detected ===\n")
print(rownames(abund_mat)[faec_rows])

faec_abundance <- colSums(abund_mat[faec_rows, , drop = FALSE])

faec_df <- data.frame(
  sample = names(faec_abundance),
  abundance = faec_abundance
)
faec_df <- merge(faec_df, metadata, by = "sample")

cat("\n=== Faecalibacterium abundance per sample ===\n")
print(faec_df[order(faec_df$condition), c("sample", "label", "condition", "abundance", "note")])
cat("\n")

p_faec <- ggplot(faec_df, aes(x = condition, y = abundance,
                              fill = condition, color = condition)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 4, alpha = 0.8) +
  geom_text(aes(label = label), nudge_x = 0.28, size = 3.5, color = "black") +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  theme_publication() +
  labs(title = expression(italic("Faecalibacterium prausnitzii")*" Abundance"),
       subtitle = "Anti-inflammatory butyrate producer | Expected depletion in IBD",
       x = NULL,
       y = "Relative Abundance (%)") +
  stat_compare_means(method = "wilcox.test", size = 4) +
  theme(legend.position = "none")

print(p_faec)
ggsave("~/crohns_mgx/figures/faecalibacterium.pdf", p_faec, width = 6, height = 5)
cat("Faecalibacterium saved\n\n")

# ============================================
# DIFFERENTIAL ABUNDANCE - MaAsLin2
# Linear model testing condition effect on each species
# Note: With n=3 per group, statistical power is very limited
# Results shown as exploratory trends only
# ============================================

dir.create("~/crohns_mgx/maaslin2_relaxed", showWarnings = FALSE, recursive = TRUE)

fit <- Maaslin2(
  input_data = t(abund_mat),
  input_metadata = metadata,
  output = "~/crohns_mgx/maaslin2_relaxed",
  fixed_effects = "condition",
  reference = "condition,Control",
  min_prevalence = 0.1,
  min_abundance = 0.001,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.3
)

all_results <- read.table("~/crohns_mgx/maaslin2_relaxed/all_results.tsv",
                          header = TRUE, sep = "\t")
all_results <- all_results[order(all_results$pval), ]

cat("=== Top 20 Differentially Abundant Species ===\n")
print(head(all_results[, c("feature", "coef", "pval", "qval")], 20))

cat("\n=== Faecalibacterium prausnitzii ===\n")
print(all_results[grep("Faecalibacterium", all_results$feature),
                  c("feature", "coef", "pval", "qval")])
cat("\n")

# Clean species names for plot
top20 <- head(all_results, 20)
top20$feature_clean <- gsub("\\.", " ", top20$feature)
top20$feature_clean <- gsub("  ", " ", top20$feature_clean)

p_diff <- ggplot(top20, aes(x = reorder(feature_clean, coef),
                            y = coef,
                            fill = coef > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr),
                width = 0.3, color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "#2166AC"),
                    labels = c("TRUE" = "Enriched in Crohn",
                               "FALSE" = "Depleted in Crohn")) +
  theme_publication() +
  labs(title = "Differential Abundance Analysis",
       subtitle = "Top 20 species ranked by p-value | No FDR correction (n=3 per group)",
       x = NULL,
       y = "Coefficient (Crohn vs Control)",
       fill = "Direction") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "italic", size = 9))

print(p_diff)
ggsave("~/crohns_mgx/figures/differential_abundance.pdf",
       p_diff, width = 10, height = 8)
cat("Differential abundance saved\n\n")

# ============================================
# SESSION INFO - for reproducibility
# ============================================

cat("=== Session Info ===\n")
print(sessionInfo())

cat("\n=== All analyses complete! ===\n")
cat("Figures saved to: ~/crohns_mgx/figures/\n")