# ============================================================
# assignment2.R
# Salmon (quant.sf) -> edgeR DE -> plots
# + GO ORA (BP) for functional analysis
#
# Saves:
#   figs/*.png  (for GitHub README, high-res)
#
# Also saves:
#   results/DE_*.csv
#   results/GO_BP_ORA_*.csv
# ============================================================

# ---------------------------
# 0) Packages
# ---------------------------
# If missing, run once:
# install.packages(c("tximport","ggplot2","pheatmap"))
# install.packages("BiocManager")
# BiocManager::install(c("edgeR","clusterProfiler","org.Sc.sgd.db","enrichplot"))

suppressPackageStartupMessages({
  library(tximport)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(enrichplot)
})

set.seed(1)

# ---------------------------
# 1) Output folders
# ---------------------------
dir.create("figs", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

PNG_DPI <- 300

# ggsave helper: PNG
save_gg_dual <- function(plot, filename_base, w=5, h=4) {
  png_file <- file.path("figs", paste0(filename_base, ".png"))
  ggsave(filename = png_file, plot = plot, width = w, height = h, dpi = PNG_DPI)
}

# pheatmap helper: PNG
save_pheatmap_dual <- function(p, filename_base, w=6, h=5) {
  png_file <- file.path("figs", paste0(filename_base, ".png"))
  png(png_file, width = w, height = h, units = "in", res = PNG_DPI)
  grid::grid.newpage(); grid::grid.draw(p$gtable)
  dev.off()
}

# ---------------------------
# 2) colors + theme
# ---------------------------
pal_group <- c(
  Early  = "#4E79A7",
  Thin   = "#59A14F", 
  Mature = "#F28E2B"   
)

pal_de <- c(
  "Not significant" = "#B8B8B8",
  "Up"   = "#D64F4F",
  "Down" = "#3B6FB6"
)

theme_nature <- function(base_size = 10, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      plot.title = element_text(face = "plain", size = base_size + 1),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# ---------------------------
# 3) Sample metadata
# ---------------------------
# Early:  65 64 63
# Thin:   62 61 60
# Mature: 59 58 57
samples <- data.frame(
  sample = c("SRR10551665","SRR10551664","SRR10551663",
             "SRR10551662","SRR10551661","SRR10551660",
             "SRR10551659","SRR10551658","SRR10551657"),
  group  = c("Early","Early","Early",
             "Thin","Thin","Thin",
             "Mature","Mature","Mature"),
  stringsAsFactors = FALSE
)
samples$group <- factor(samples$group, levels = c("Early","Thin","Mature"))

# ---------------------------
# 4) Import Salmon quant.sf
# ---------------------------
files <- file.path(paste0(samples$sample, "_quant"), "quant.sf")
names(files) <- samples$sample

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

# txOut=TRUE keeps transcript-level features
txi <- tximport(files, type = "salmon", txOut = TRUE)

# ---------------------------
# 5) edgeR workflow
# ---------------------------
y <- DGEList(counts = txi$counts, group = samples$group)

design <- model.matrix(~ group, data = samples)

keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)       # TMM
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

message("Design columns: ", paste(colnames(design), collapse = ", "))

# ---------------------------
# 6) Differential expression
# ---------------------------
lrt_E_T <- glmLRT(fit, coef = "groupThin")       # Thin vs Early
lrt_E_M <- glmLRT(fit, coef = "groupMature")     # Mature vs Early
lrt_T_M <- glmLRT(fit, contrast = c(0, -1, 1))   # Mature - Thin

tab_E_T <- topTags(lrt_E_T, n = Inf)$table
tab_E_M <- topTags(lrt_E_M, n = Inf)$table
tab_T_M <- topTags(lrt_T_M, n = Inf)$table

write.csv(tab_E_T, file.path("results", "DE_Early_vs_Thin_edgeR.csv"))
write.csv(tab_E_M, file.path("results", "DE_Early_vs_Mature_edgeR.csv"))
write.csv(tab_T_M, file.path("results", "DE_Thin_vs_Mature_edgeR.csv"))

cat("\nSignificant genes (FDR < 0.05):\n")
cat("Early vs Thin  :", sum(tab_E_T$FDR < 0.05), "\n")
cat("Early vs Mature:", sum(tab_E_M$FDR < 0.05), "\n")
cat("Thin vs Mature :", sum(tab_T_M$FDR < 0.05), "\n")

# ---------------------------
# 7) logCPM for plots
# ---------------------------
logCPM <- cpm(y, log = TRUE, prior.count = 2)

anno_col <- data.frame(Group = samples$group)
rownames(anno_col) <- samples$sample
anno_col <- anno_col[colnames(logCPM), , drop = FALSE]
ann_colors <- list(Group = pal_group)

# ---------------------------
# 8) PCA plot
# ---------------------------
pca <- prcomp(t(logCPM), scale. = FALSE)
pve <- (pca$sdev^2) / sum(pca$sdev^2) * 100

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  group = anno_col$Group
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3.0) +
  scale_color_manual(values = pal_group) +
  labs(
    title = "PCA",
    x = sprintf("PC1 (%.1f%%)", pve[1]),
    y = sprintf("PC2 (%.1f%%)", pve[2]),
    color = NULL
  ) +
  theme_nature(base_size = 10) +
  theme(legend.position = "right")

save_gg_dual(p_pca, "Fig_PCA", w=4.6, h=3.8)

# ---------------------------
# 9) Volcano plots
# ---------------------------
make_volcano_plot <- function(tab, title_text) {
  df <- tab
  df$Direction <- "Not significant"
  df$Direction[df$FDR < 0.05 & df$logFC > 0] <- "Up"
  df$Direction[df$FDR < 0.05 & df$logFC < 0] <- "Down"
  df$Direction <- factor(df$Direction, levels = c("Not significant","Down","Up"))
  
  ggplot(df, aes(x = logFC, y = -log10(FDR), color = Direction)) +
    geom_point(alpha = 0.75, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = pal_de) +
    labs(title = title_text, x = "log2 fold change", y = "-log10(FDR)", color = NULL) +
    theme_nature(base_size = 10) +
    theme(legend.position = "right")
}

save_gg_dual(make_volcano_plot(tab_E_T, "Volcano: Early vs Thin"),
             "Fig_Volcano_Early_vs_Thin", 4.6, 3.8)
save_gg_dual(make_volcano_plot(tab_E_M, "Volcano: Early vs Mature"),
             "Fig_Volcano_Early_vs_Mature", 4.6, 3.8)
save_gg_dual(make_volcano_plot(tab_T_M, "Volcano: Thin vs Mature"),
             "Fig_Volcano_Thin_vs_Mature", 4.6, 3.8)

# ---------------------------
# 10) MA plots
# ---------------------------
make_ma_plot <- function(tab, title_text) {
  df <- tab
  df$Direction <- "Not significant"
  df$Direction[df$FDR < 0.05 & df$logFC > 0] <- "Up"
  df$Direction[df$FDR < 0.05 & df$logFC < 0] <- "Down"
  df$Direction <- factor(df$Direction, levels = c("Not significant","Down","Up"))
  
  ggplot(df, aes(x = logCPM, y = logFC, color = Direction)) +
    geom_point(alpha = 0.75, size = 1.2) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = pal_de) +
    labs(title = title_text, x = "Average logCPM", y = "log2 fold change", color = NULL) +
    theme_nature(base_size = 10) +
    theme(legend.position = "right")
}

save_gg_dual(make_ma_plot(tab_E_T, "MA: Early vs Thin"),
             "Fig_MA_Early_vs_Thin", 4.6, 3.8)
save_gg_dual(make_ma_plot(tab_E_M, "MA: Early vs Mature"),
             "Fig_MA_Early_vs_Mature", 4.6, 3.8)
save_gg_dual(make_ma_plot(tab_T_M, "MA: Thin vs Mature"),
             "Fig_MA_Thin_vs_Mature", 4.6, 3.8)

# ---------------------------
# 11) Heatmap (Top 50 DEGs, Early vs Thin)
# ---------------------------
top_n <- 50
top_genes <- rownames(tab_E_T)[1:top_n]
heatmap_data <- logCPM[top_genes, , drop = FALSE]

hm_cols <- colorRampPalette(c("#2C7FB8", "white", "#D95F0E"))(101)

p_hm <- pheatmap(
  heatmap_data,
  annotation_col = anno_col,
  annotation_colors = ann_colors,
  color = hm_cols,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 8,
  fontsize_col = 8,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  silent = TRUE
)

save_pheatmap_dual(p_hm, "Fig_Heatmap_Top50_Early_vs_Thin", w=5.6, h=5.2)

# ---------------------------
# 12) Functional analysis (ORA): GO BP enrichment
# ---------------------------
# org.Sc.sgd.db is required for yeast GO mapping
if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
  stop("Missing package: org.Sc.sgd.db\nInstall with:\nBiocManager::install('org.Sc.sgd.db')")
}
suppressPackageStartupMessages(library(org.Sc.sgd.db))

run_go_bp_ora <- function(tab, contrast_name) {
  sig <- subset(tab, FDR < 0.05)
  genes_orf <- sub("_mRNA$", "", rownames(sig))  # YEL069C_mRNA -> YEL069C
  
  ego_bp <- enrichGO(
    gene          = genes_orf,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  write.csv(as.data.frame(ego_bp),
            file.path("results", paste0("GO_BP_ORA_", contrast_name, ".csv")),
            row.names = FALSE)
  
  p_go <- dotplot(ego_bp, showCategory = 12) +
    labs(title = paste0("GO BP ORA: ", gsub("_", " ", contrast_name)),
         x = "Gene ratio", y = NULL) +
    theme_nature(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))
  
  save_gg_dual(p_go, paste0("Fig_GO_BP_ORA_", contrast_name), w=6.2, h=4.2)
  
  ego_bp
}

ego_E_T <- run_go_bp_ora(tab_E_T, "Early_vs_Thin")
ego_E_M <- run_go_bp_ora(tab_E_M, "Early_vs_Mature")
ego_T_M <- run_go_bp_ora(tab_T_M, "Thin_vs_Mature")

# ---------------------------
# 13) Done
# ---------------------------
cat("\nDONE ✅\n")
cat("Figures saved to: figs/  (PNG for README + PDF for submission)\n")
cat("Tables saved to:  results/\n")
cat("\nREADME usage example:\n")
cat("![PCA](figs/Fig_PCA.png)\n")
cat("![Heatmap](figs/Fig_Heatmap_Top50_Early_vs_Thin.png)\n")