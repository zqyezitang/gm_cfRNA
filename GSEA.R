# ============================================================
# Gynecologic cancer GSEA analysis (TPM, all genes)
# Purpose: use the same sample grouping as DEG analysis and perform
#          GO Biological Process GSEA with the COMPLETE ranked gene list
# Groups: Cervical_cancer, CIN, Endometrial_cancer, Ovarian_cancer vs Control
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

set.seed(123)

# -----------------------------
# 1. Paths and parameters
# -----------------------------
WORK_DIR <- "C:/Users/wusihan/Desktop/Gyne_GSEA_Revision_AllGenes"
setwd(WORK_DIR)

GROUP_FILE <- "group_info.csv"

# Prefer the cleaned TPM file generated during QC; otherwise use the original file.
TPM_FILE <- if (file.exists("Samples_combined_TPM_clean_166.csv")) {
  "Samples_combined_TPM_clean_166.csv"
} else {
  "Samples_combined_TPM.csv"
}

OUT_DIR <- "Results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CONTROL_LABEL <- "Control"
GSEA_ONT <- "BP"
MIN_GS_SIZE <- 10
MAX_GS_SIZE <- 500
PADJ_CUTOFF <- 0.05
TOP_N_PLOT <- 8

cat("TPM file:", TPM_FILE, "\n")
cat("Group file:", GROUP_FILE, "\n")

# -----------------------------
# 2. Read sample grouping
# -----------------------------
group_info <- fread(GROUP_FILE)

stopifnot(all(c("Sample", "Group") %in% colnames(group_info)))

group_info[, Sample := as.character(Sample)]
group_info[, Group  := as.character(Group)]
group_info <- unique(group_info[!is.na(Sample) & Sample != "" &
                                !is.na(Group)  & Group  != ""])

if (!(CONTROL_LABEL %in% group_info$Group)) {
  stop("Control group not found in group_info.csv")
}

cat("\nSample numbers by group:\n")
print(table(group_info$Group))
cat("Total samples:", nrow(group_info), "\n")

# -----------------------------
# 3. Read TPM matrix (only required samples)
# -----------------------------
header <- names(fread(TPM_FILE, nrows = 0))
gene_col <- header[1]

missing_samples <- setdiff(group_info$Sample, header)
if (length(missing_samples) > 0) {
  stop("Samples in group_info.csv not found in TPM matrix: ",
       paste(missing_samples, collapse = ", "))
}

sample_cols <- group_info$Sample
use_cols <- c(gene_col, sample_cols)
tpm <- fread(TPM_FILE, select = use_cols,
             na.strings = c("", "NA", "N/A", "NULL", "/"))

setnames(tpm, gene_col, "GeneID")
tpm[, GeneID := as.character(GeneID)]

# Convert TPM columns to numeric.
for (j in sample_cols) {
  set(tpm, j = j,
      value = suppressWarnings(as.numeric(as.character(tpm[[j]]))))
}

expr_raw <- as.matrix(tpm[, ..sample_cols])
storage.mode(expr_raw) <- "double"

# Remove rows that cannot be analyzed (invalid gene ID / NA / Inf / all zero).
bad_row <- is.na(tpm$GeneID) | tpm$GeneID == "" |
           !apply(is.finite(expr_raw), 1, all) |
           rowSums(expr_raw, na.rm = TRUE) == 0

cat("\nGenes before QC:", nrow(tpm), "\n")
cat("Rows removed by QC:", sum(bad_row), "\n")

tpm <- tpm[!bad_row]
expr_raw <- expr_raw[!bad_row, , drop = FALSE]
gene_ids <- tpm$GeneID

cat("Genes after QC:", nrow(expr_raw), "\n")

# -----------------------------
# 4. log2(TPM + 1) and limma design
# -----------------------------
expr <- log2(expr_raw + 1)
rownames(expr) <- gene_ids

# Keep expression columns in exactly the same order as group_info.
expr <- expr[, sample_cols, drop = FALSE]

group_factor <- factor(group_info$Group)
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)

fit <- lmFit(expr, design)

case_groups <- setdiff(levels(group_factor), CONTROL_LABEL)
if (length(case_groups) == 0) stop("No disease groups found.")

contrast_strings <- paste0("`", case_groups, "`-`", CONTROL_LABEL, "`")
contrast_matrix <- makeContrasts(contrasts = contrast_strings, levels = design)
colnames(contrast_matrix) <- paste0(case_groups, "_vs_", CONTROL_LABEL)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 5. Gene ID handling
# -----------------------------
# Detect Ensembl IDs automatically. If Ensembl, map them to SYMBOL.
is_ensembl <- mean(grepl("^ENSG", gene_ids)) > 0.5

if (is_ensembl) {
  gene_ids_clean <- sub("\\..*$", "", gene_ids)
  symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(gene_ids_clean),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  gene_symbol <- unname(symbol_map[gene_ids_clean])
} else {
  gene_symbol <- gene_ids
}

# Helper: build a complete ranked gene list for GSEA.
# Duplicate symbols are reduced by keeping the row with the largest absolute t statistic.
make_ranked_list <- function(stat_vec, symbols) {
  d <- data.frame(
    SYMBOL = as.character(symbols),
    stat = as.numeric(stat_vec),
    stringsAsFactors = FALSE
  )

  d <- d[!is.na(d$SYMBOL) & d$SYMBOL != "" & is.finite(d$stat), ]
  d <- d[order(abs(d$stat), decreasing = TRUE), ]
  d <- d[!duplicated(d$SYMBOL), ]
  d <- d[order(d$stat, decreasing = TRUE), ]

  gene_list <- d$stat
  names(gene_list) <- d$SYMBOL
  gene_list
}

# -----------------------------
# 6. GSEA for every disease group vs Control
# -----------------------------
all_gsea_results <- list()
gsea_objects <- list()
ranked_lists <- list()

for (comparison in colnames(contrast_matrix)) {
  cat("\n============================================================\n")
  cat("Running:", comparison, "\n")
  cat("============================================================\n")

  # limma moderated t statistic for ALL usable genes
  stat_vec <- fit2$t[, comparison]
  gene_list <- make_ranked_list(stat_vec, gene_symbol)
  ranked_lists[[comparison]] <- gene_list

  cat("Genes used in ranked list:", length(gene_list), "\n")

  # Save complete ranked list (proof that no top-1000 truncation was used)
  fwrite(
    data.table(SYMBOL = names(gene_list), t_statistic = as.numeric(gene_list)),
    file.path(OUT_DIR, paste0(comparison, "_ALL_genes_ranked.csv"))
  )

  gsea <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = GSEA_ONT,
    minGSSize = MIN_GS_SIZE,
    maxGSSize = MAX_GS_SIZE,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE
  )

  gsea_objects[[comparison]] <- gsea
  saveRDS(gsea, file.path(OUT_DIR, paste0(comparison, "_gsea_object.rds")))

  res <- as.data.frame(gsea)
  if (nrow(res) > 0) {
    res$Comparison <- comparison
    all_gsea_results[[comparison]] <- res

    write.csv(
      res,
      file.path(OUT_DIR, paste0(comparison, "_GSEA_GO_BP.csv")),
      row.names = FALSE
    )

    # Top 3 enrichment curves by adjusted P value
    top_ids <- head(res$ID[order(res$p.adjust)], 3)
    if (length(top_ids) > 0) {
      p_curve <- gseaplot2(
        gsea,
        geneSetID = top_ids,
        title = paste0(comparison, " - Top 3 GO-BP pathways"),
        pvalue_table = TRUE
      )
      ggsave(
        file.path(OUT_DIR, paste0(comparison, "_top3_GSEA_curves.png")),
        p_curve, width = 12, height = 8, dpi = 300
      )
    }
  }
}

# -----------------------------
# 7. Combine and save result tables
# -----------------------------
if (length(all_gsea_results) == 0) {
  stop("No GSEA result was returned.")
}

all_results <- bind_rows(all_gsea_results)
all_results <- all_results %>%
  select(Comparison, everything())

write.csv(
  all_results,
  file.path(OUT_DIR, "01_ALL_GSEA_results.csv"),
  row.names = FALSE
)

sig_results <- all_results %>%
  filter(!is.na(p.adjust), p.adjust < PADJ_CUTOFF) %>%
  arrange(Comparison, p.adjust)

write.csv(
  sig_results,
  file.path(OUT_DIR, "02_Significant_GSEA_padj_lt_0.05.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Bubble plot across groups
# -----------------------------
plot_df <- all_results %>%
  filter(!is.na(p.adjust), is.finite(NES), p.adjust > 0) %>%
  group_by(Comparison) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = TOP_N_PLOT) %>%
  ungroup() %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    Comparison = factor(Comparison, levels = colnames(contrast_matrix))
  )

# Order pathway labels for readability.
pathway_order <- plot_df %>%
  arrange(Comparison, p.adjust) %>%
  pull(Description) %>%
  unique()
plot_df$Description <- factor(plot_df$Description, levels = rev(pathway_order))

p_bubble <- ggplot(
  plot_df,
  aes(x = Comparison, y = Description,
      size = neglog10_padj, color = NES)
) +
  geom_point(alpha = 0.95) +
  scale_color_gradient2(
    low = "#2448A6", mid = "white", high = "#D7191C", midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(name = "-log10(adj.P)") +
  labs(
    title = "Top GSEA pathways across groups",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUT_DIR, "04_GSEA_bubbleplot.png"),
  p_bubble, width = 13, height = 10, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "04_GSEA_bubbleplot.pdf"),
  p_bubble, width = 13, height = 10
)

# -----------------------------
# 9. NES heatmap
# -----------------------------
heat_df <- plot_df %>%
  select(Comparison, Description, NES) %>%
  distinct() %>%
  pivot_wider(names_from = Comparison, values_from = NES)

if (nrow(heat_df) > 1) {
  heat_mat <- as.data.frame(heat_df)
  rownames(heat_mat) <- heat_mat$Description
  heat_mat$Description <- NULL
  heat_mat <- as.matrix(heat_mat)

  png(file.path(OUT_DIR, "03_NES_heatmap.png"),
      width = 2400, height = 2400, res = 300)
  pheatmap(
    heat_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Top GSEA pathways: NES",
    fontsize_row = 8,
    fontsize_col = 10
  )
  dev.off()
}

# -----------------------------
# 10. Simple analysis report
# -----------------------------
report_file <- file.path(OUT_DIR, "00_READ_ME_FIRST.txt")
sink(report_file)
cat("GSEA analysis summary\n")
cat("=====================\n")
cat("TPM file:", TPM_FILE, "\n")
cat("Total samples:", nrow(group_info), "\n")
cat("Groups:\n")
print(table(group_info$Group))
cat("\nGenes after expression QC:", nrow(expr), "\n")
cat("Ranking metric: limma moderated t statistic\n")
cat("GSEA input: complete ranked gene list (no top-1000 truncation)\n")
cat("Ontology: GO Biological Process (BP)\n")
cat("Adjusted P-value method: BH\n")
cat("Significant cutoff: p.adjust <", PADJ_CUTOFF, "\n\n")

cat("Significant pathways by comparison:\n")
if (nrow(sig_results) > 0) {
  print(table(sig_results$Comparison))
} else {
  cat("No pathways passed the adjusted P-value cutoff.\n")
}

cat("\nTop pathways by adjusted P value:\n")
print(
  all_results %>%
    arrange(p.adjust) %>%
    select(Comparison, Description, NES, p.adjust) %>%
    head(20)
)
sink()

# Save R session/package versions for reproducibility.
sink(file.path(OUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nAnalysis completed. Results saved in:", OUT_DIR, "\n")
