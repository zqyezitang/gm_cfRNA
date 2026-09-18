library(DESeq2)
library(readr)
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(cowplot)


counts_path <- "Tissue_samples_COUNTS.csv"
tpm_path <- "Tissue_samples_TPM.csv"
#meta_path <- "Sample_tissue_info.xlsx"
meta_path <- "Tissue_paired_info.xlsx"

#
outdir <- "DEG_results_paired"
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "volcano"), showWarnings = FALSE)
dir.create(file.path(outdir, "heatmap"), showWarnings = FALSE)
dir.create(file.path(outdir, "sig_gene_lists"), showWarnings = FALSE)
dir.create(file.path(outdir, "combined_plots"), showWarnings = FALSE)

#
outlier_samples <- c("2488-CQCU034A", "2546-CQCU092A", "MPU-CQCH1499A",
                     "2518-CQCU064A", "MPU-CQCU780A", "MPU-CQCU941A",
                     "MPU-CQCH425A", "MPU-CQCH543A")

#
count_data <- read_csv(counts_path, show_col_types = FALSE)
gene_ids <- count_data[[1]]
count_mat <- as.data.frame(count_data[,-1])
rownames(count_mat) <- gene_ids
colnames(count_mat) <- colnames(count_data)[-1]
count_mat <- count_mat[, !colnames(count_mat) %in% outlier_samples]

tpm_data <- read_csv(tpm_path, show_col_types = FALSE)
gene_ids_tpm <- tpm_data[[1]]
expr_data <- as.data.frame(tpm_data[,-1])
rownames(expr_data) <- gene_ids_tpm
colnames(expr_data) <- colnames(tpm_data)[-1]
expr_data <- expr_data[, !colnames(expr_data) %in% outlier_samples]

# 
meta <- read_excel(meta_path)
meta <- meta %>% filter(sample %in% colnames(count_mat))
count_mat <- count_mat[, meta$sample]
expr_data <- expr_data[, meta$sample]

#
get_pca_plot <- function(df_expr, df_meta, group_col, shape_col = NULL, title) {
  expr_t <- as.data.frame(t(df_expr))
  expr_t$sample <- rownames(expr_t)
  df <- left_join(expr_t, df_meta, by = "sample")

  rownames(df) <- df$sample
  expr_mat <- df %>% select(where(is.numeric)) %>% as.matrix()
  expr_log <- log2(expr_mat + 1)
  expr_log <- expr_log[, apply(expr_log, 2, function(x) sd(x) > 0)]

  pca <- prcomp(expr_log, center = TRUE, scale. = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2]
  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$sample <- rownames(pca_df)
  pca_df <- left_join(pca_df, df_meta, by = "sample")

  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
    geom_point(aes(shape = if (!is.null(shape_col)) .data[[shape_col]] else NULL),
               size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("red", "blue")) +
    labs(title = title,
         x = paste0("PC1 (", round(var_exp[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2]*100, 1), "%)")) +
    theme_minimal(base_size = 11)
}

# 
run_analysis <- function(count_sub, expr_sub, meta_sub,name_col, group_col, contrast_pair, prefix) {
  dds <- DESeqDataSetFromMatrix(countData = round(count_sub),
                                colData = meta_sub,
                                design = as.formula(paste0("~ ",name_col," + ", group_col)))
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(group_col, contrast_pair[1], contrast_pair[2]))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  meta_sub[[group_col]] <- factor(meta_sub[[group_col]])

  write.csv(res_df, file = file.path(outdir, paste0("pvalue_DEG_", prefix, ".csv")), row.names = FALSE)
  res_sig <- res_df %>% filter(padj < 0.001 & abs(log2FoldChange) > 1)
  write.csv(res_sig, file = file.path(outdir, paste0("pvalue_DEG_", prefix, "_sig.csv")), row.names = FALSE)


  gene_map <- tibble(
    ensembl_id = rownames(count_mat),
    gene_id = gene_ids
  )

  res_up <- res_sig %>% filter(log2FoldChange > 1)
  up_df <- gene_map %>%
    filter(ensembl_id %in% res_up$gene) %>%
    left_join(res_up %>% select(gene, log2FoldChange), by = c("ensembl_id" = "gene")) %>%
    select(gene_id, log2FoldChange)

  res_down <- res_sig %>% filter(log2FoldChange < -1)
  down_df <- gene_map %>%
    filter(ensembl_id %in% res_down$gene) %>%
    left_join(res_down %>% select(gene, log2FoldChange), by = c("ensembl_id" = "gene")) %>%
    select(gene_id, log2FoldChange)

  write.table(up_df, file = file.path(outdir, "sig_gene_lists", paste0("pvalue_sig_genes_up_", prefix, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  write.table(down_df, file = file.path(outdir, "sig_gene_lists", paste0("pvalue_sig_genes_down_", prefix, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


  res_df <- res_df %>%
    mutate(sig = ifelse(padj < 0.001 & abs(log2FoldChange) > 1,
                        ifelse(log2FoldChange > 1, "Up", "Down"), "NS"))
 # top_up <- res_df %>% filter(sig == "Up") %>% arrange(pvalue) %>% slice(1:5)
 # top_down <- res_df %>% filter(sig == "Down") %>% arrange(pvalue) %>% slice(1:5)

  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.001), linetype = "dashed") +
  #  geom_label_repel(data = rbind(top_up, top_down), aes(label = gene),
  #                   size = 1.5, fontface = "italic", fill = "white", color = "black",
   #                  box.padding = 0.2, label.size = 0.1, segment.color = "black") +
    theme_minimal(base_size = 10) +
    labs(title = paste0("Volcano: ", prefix),
         x = "log2 Fold Change", y = "-log10 adjusted p-value")

  ggsave(file.path(outdir, "volcano", paste0("padj_volcano_", prefix, ".pdf")),
         volcano, width = 8, height = 5)

  if (nrow(res_sig) >= 2) {
    mat <- assay(vst(dds, blind = TRUE))
    heat_genes <- intersect(res_sig$gene, rownames(mat))
    mat <- mat[heat_genes[1:min(50, length(heat_genes))], ]
    mat <- t(scale(t(mat)))

   # pdf(file.path(outdir, "heatmap", paste0("padj_Heatmap_", prefix, ".pdf")), width = 8, height = 10)
    anno_col <- meta_sub %>% select(all_of(group_col)) %>% as.data.frame()
    rownames(anno_col) <- meta_sub$sample

    if (ncol(anno_col) > 0 && all(!is.na(anno_col[[1]]))) {
      pheatmap(mat,
               annotation_col = anno_col,
               color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
               show_rownames = TRUE, show_colnames = FALSE,
               filename = paste0(outdir,'/heatmap/Tissue_heatmap.pdf'),width = 8, height = 10
             #  main = paste0("Top DEGs Heatmap: ", prefix)
               )
    } else {
      pheatmap(mat,
               color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
               show_rownames = TRUE, show_colnames = FALSE,
               main = paste0("Top DEGs Heatmap: ", prefix))
    }

   # dev.off()
  }

  expr_sub <- expr_sub[, meta_sub$sample]
  group_values <- unique(meta_sub[[group_col]])
  pca_title <- paste0("PCA: ", prefix)
  pca_plot <- get_pca_plot(expr_sub, meta_sub, group_col,
                           shape_col = "treat_time",
                           title = pca_title)

  combined_plot <- plot_grid(volcano, pca_plot, ncol = 2, rel_widths = c(1, 1))
  ggsave(file.path(outdir, "combined_plots", paste0("pvalue_volcano_PCA_", prefix, ".pdf")),
         combined_plot, width = 12, height = 5)

}
# 🎯 Cancer vs Control
cancer_types <- unique(meta$group)
cancer_types <- setdiff(cancer_types, "Cervical_adjacent")

for (ct in cancer_types) {
  meta_sub <- meta %>% filter(group %in% c("Cervical_adjacent", ct))
  meta_sub$name <- as.factor(meta_sub$name)
  meta_sub$group <- as.factor(meta_sub$group)
  run_analysis(count_mat[, meta_sub$sample], expr_data[, meta_sub$sample],
               meta_sub,"name", "group", c(ct, "Cervical_adjacent"), paste0(ct, "_vs_adjacent"))
}



###################comparing with cfRNA
library(confintr)

# ============================================================
res1 <- read.csv("pvalue_DEG_Cervical_tumor_vs_adjacent_sig.csv", stringsAsFactors = FALSE)
res2 <- read.csv("CC_vs_Control_ALL_results.csv", stringsAsFactors = FALSE)

common_genes <- intersect(res1$Gene, res2$Gene)
sub1 <- res1[res1$Gene %in% common_genes, c("Gene", "log2FoldChange")]
sub2 <- res2[res2$Gene %in% common_genes, c("Gene", "log2FoldChange")]
merged <- merge(sub1, sub2, by = "Gene", suffixes = c("_1", "_2"))


#
lfc1 <- merged$log2FoldChange_1
lfc2 <- merged$log2FoldChange_2

# ============================================================
library(confintr)

# 3.1 Spearman
spearman_test <- cor.test(lfc1, lfc2, method = "spearman")
cat("Spearman rho =", round(spearman_test$estimate, 4),
    ", P =", format.pval(spearman_test$p.value, digits = 3),
    ", n =", length(lfc1), "\n")

# 3.2 Pearson
pearson_test <- cor.test(lfc1, lfc2, method = "pearson")
cat("Pearson r =", round(pearson_test$estimate, 4),
    ", P =", format.pval(pearson_test$p.value, digits = 3),
    ", n =", length(lfc1), "\n")

# ============================================================
#Permutation Test
# ============================================================
library(CarletonStats)

set.seed(123)
perm_spearman <- permTestCor(lfc1, lfc2, method = "spearman",
                             B = 9999, alternative = "two.sided")
print(perm_spearman)

# 
set.seed(123)
perm_pearson <- permTestCor(lfc1, lfc2, method = "pearson",
                            B = 9999, alternative = "two.sided")
print(perm_pearson)

# 5.2 置换分布直方图（手动实现，便于自定义绘图）
set.seed(123)
B <- 9999
perm_rho <- numeric(B)
for (i in 1:B) {
  perm_rho[i] <- cor(lfc1, sample(lfc2), method = "spearman")
}
observed_rho <- cor(lfc1, lfc2, method = "spearman")

pdf("permutation_distribution.pdf", width = 8, height = 6)
hist(perm_rho, breaks = 50, col = "lightblue",
     #xlim=range(-0.02,0.08),
     main = "Permutation distribution of Spearman's rho",
     xlab = "Spearman's rho (permuted)")
abline(v = observed_rho, col = "red", lwd = 2)
abline(v = -observed_rho, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = c(paste0("Observed rho = ", round(observed_rho, 3)),
                  "Null distribution"),
       col = c("red", "lightblue"), lwd = c(2, NA), pch = c(NA, 15))
dev.off()

perm_p <- (sum(abs(perm_rho) >= abs(observed_rho)) + 1) / (B + 1)
cat("\n P :", perm_p, "\n")

# ============================================================
# 
# ============================================================
summary_table <- data.frame(
  Method = c("Spearman", "Pearson"),
  Coefficient = c(round(spearman_test$estimate, 3),
                  round(pearson_test$estimate, 3)),
  CI_lower = c(round(spearman_ci$interval[1], 3),
               round(pearson_ci$interval[1], 3)),
  CI_upper = c(round(spearman_ci$interval[2], 3),
               round(pearson_ci$interval[2], 3)),
  P_value = c(format.pval(spearman_test$p.value, digits = 3),
              format.pval(pearson_test$p.value, digits = 3)),
  N_genes = c(length(lfc1), length(lfc1))
)
print(summary_table)

write.csv(summary_table, "correlation_summary.csv", row.names = FALSE)

############# Scatter plot of fold-change
data <- data.frame(
  x = lfc1,
  y = lfc2,
  gene =merged$Gene
  #  condition = rep(c("Treatment", "Control"), each = 50)
)

data$quadrant <- case_when(
  data$x >= 0 & data$y >= 0 ~ "Q1 (+,+)",
  data$x < 0 & data$y >= 0 ~ "Q2 (-,+)",
  data$x < 0 & data$y < 0 ~ "Q3 (-,-)",
  data$x >= 0 & data$y < 0 ~ "Q4 (+,-)"
)
valid_genes <- c("STAT1","IFI30", "CLEC7A","IGLV2-14","IFI44L") #top 5
top_genes <- data %>% 
  filter(gene %in% valid_genes)

pdf('scatter_FC.pdf')
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(alpha = 0.7, size = 2,color="lightblue") +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  #  fontface = "italic",
                  size = 3.5,
                  fontface = "bold.italic",
                  color = "red",
                  box.padding = 0.6,
                  point.padding = 0.5,
                  arrow = arrow(length = unit(0.01, "npc"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", linewidth = 1) +
  labs(title = "",
       x = "Log2(fold-change) in tissue", y = "Log2(fold-change) in plasma") +
  theme_classic() +
  coord_cartesian(ylim = c(-1.4, 1.4))
p
dev.off()
