# =====================================================================
# Gynecologic tumor RNA-seq differential expression analysis
# Revision version: covariates in DESeq2 design + scaled age/BMI
#
# Comparisons:
#   CIN vs Control
#   CC  (Cervical_cancer) vs Control
#   OC  (Ovarian_cancer) vs Control
#   EC  (Endometrial_cancer) vs Control
#


# =====================================================================
# 0. USER-ADJUSTABLE PARAMETERS
# =====================================================================
# If the teacher changes the cutoff later, usually only change these values.
PADJ_CUTOFF <- 0.10
FC_CUTOFF   <- 1.20

# Input file names
CLINICAL_FILE <- file.path(
  "input",
  "TableS1. Clinical information of samples.xlsx"
)

COUNT_FILE <- file.path(
  "input",
  "Paper_sample_Count.txt"
)

# Clinical Excel sheet
CLINICAL_SHEET <- "Samples_batch"

# Output directory
RESULT_DIR <- "result"

# =====================================================================
# 1. PACKAGE CHECK
# =====================================================================

required_cran <- c("readxl")
required_bioc <- c("DESeq2")

# Install CRAN packages if missing
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install BiocManager and Bioconductor packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
})


# =====================================================================
# 2. CHECK PATHS
# =====================================================================

dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(CLINICAL_FILE)) {
  stop(
    "\nCannot find clinical file:\n",
    CLINICAL_FILE,
    "\n\nPlease check the working directory with getwd()."
  )
}

if (!file.exists(COUNT_FILE)) {
  stop(
    "\nCannot find count file:\n",
    COUNT_FILE,
    "\n\nPlease check the working directory with getwd()."
  )
}

cat("\n============================================================\n")
cat("PROJECT CHECK\n")
cat("============================================================\n")
cat("Working directory:\n", getwd(), "\n")
cat("Clinical file: OK\n")
cat("Count file:    OK\n\n")


# =====================================================================
# 3. READ CLINICAL INFORMATION
# =====================================================================

meta <- read_excel(
  CLINICAL_FILE,
  sheet = CLINICAL_SHEET
)

meta <- as.data.frame(meta, stringsAsFactors = FALSE)

required_meta_cols <- c(
  "sample",
  "cancer",
  "batch",
  "age",
  "height",
  "weight"
)

missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing required clinical columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

# Convert sample and grouping columns to character
meta$sample  <- trimws(as.character(meta$sample))
meta$cancer  <- trimws(as.character(meta$cancer))
meta$batch   <- trimws(as.character(meta$batch))

if (anyDuplicated(meta$sample)) {
  duplicated_samples <- unique(meta$sample[duplicated(meta$sample)])
  stop(
    "Duplicated sample IDs were found in the clinical table: ",
    paste(duplicated_samples, collapse = ", ")
  )
}

cat("Original clinical samples:", nrow(meta), "\n")
cat("\nCancer groups in the clinical table:\n")
print(table(meta$cancer, useNA = "ifany"))


# =====================================================================
# 4. READ RAW COUNT MATRIX
# =====================================================================

count_df <- read.delim(
  COUNT_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (!"Gene" %in% colnames(count_df)) {
  stop(
    "Column 'Gene' was not found in the count matrix. ",
    "Please confirm that the first column is named Gene."
  )
}

if (anyDuplicated(count_df$Gene)) {
  duplicated_genes <- unique(count_df$Gene[duplicated(count_df$Gene)])
  stop(
    "Duplicated gene IDs/names were found in the count matrix. ",
    "Examples: ",
    paste(head(duplicated_genes, 20), collapse = ", ")
  )
}

rownames(count_df) <- count_df$Gene
count_df$Gene <- NULL

count_mat <- as.matrix(count_df)

# Validate count matrix
suppressWarnings(storage.mode(count_mat) <- "numeric")

if (anyNA(count_mat)) {
  stop("NA values were found in the count matrix.")
}

if (any(!is.finite(count_mat))) {
  stop("Non-finite values were found in the count matrix.")
}

if (any(count_mat < 0)) {
  stop("Negative values were found in the count matrix.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop(
    "The expression matrix contains non-integer values. ",
    "DESeq2 requires raw integer counts."
  )
}

storage.mode(count_mat) <- "integer"

cat("\nOriginal count samples:", ncol(count_mat), "\n")
cat("Original genes:", nrow(count_mat), "\n")


# =====================================================================
# 5. REMOVE THE 8 ABNORMAL SAMPLES
# =====================================================================

abnormal_samples <- c(
  "2488-CQCU034A",
  "2546-CQCU092A",
  "MPU-CQCH1499A",
  "2518-CQCU064A",
  "MPU-CQCU780A",
  "MPU-CQCU941A",
  "MPU-CQCH425A",
  "MPU-CQCH543A"
)

abnormal_check <- data.frame(
  sample = abnormal_samples,
  found_in_clinical_before_removal =
    abnormal_samples %in% meta$sample,
  found_in_count_before_removal =
    abnormal_samples %in% colnames(count_mat),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("ABNORMAL SAMPLE CHECK\n")
cat("============================================================\n")
print(abnormal_check)

write.csv(
  abnormal_check,
  file.path(RESULT_DIR, "01_abnormal_sample_check.csv"),
  row.names = FALSE
)

# Remove abnormal samples from both data sources if present
meta <- meta[
  !meta$sample %in% abnormal_samples,
  ,
  drop = FALSE
]

count_mat <- count_mat[
  ,
  !colnames(count_mat) %in% abnormal_samples,
  drop = FALSE
]

# Safety check
if (any(abnormal_samples %in% meta$sample)) {
  stop("Some abnormal samples were not successfully removed from metadata.")
}

if (any(abnormal_samples %in% colnames(count_mat))) {
  stop("Some abnormal samples were not successfully removed from count matrix.")
}


# =====================================================================
# 6. MATCH CLINICAL SAMPLES WITH COUNT-MATRIX SAMPLES
# =====================================================================

clinical_only <- setdiff(meta$sample, colnames(count_mat))
count_only    <- setdiff(colnames(count_mat), meta$sample)

cat("\n============================================================\n")
cat("SAMPLE MATCHING CHECK\n")
cat("============================================================\n")

cat("Clinical-only samples:", length(clinical_only), "\n")
if (length(clinical_only) > 0) {
  print(clinical_only)
}

cat("\nCount-only samples:", length(count_only), "\n")
if (length(count_only) > 0) {
  print(count_only)
}

# Save mismatch lists
max_len <- max(length(clinical_only), length(count_only), 1)

mismatch_table <- data.frame(
  clinical_only = c(
    clinical_only,
    rep(NA_character_, max_len - length(clinical_only))
  ),
  count_only = c(
    count_only,
    rep(NA_character_, max_len - length(count_only))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mismatch_table,
  file.path(RESULT_DIR, "02_sample_mismatch.csv"),
  row.names = FALSE
)

# Samples must have BOTH expression and clinical information
common_samples <- intersect(meta$sample, colnames(count_mat))

if (length(common_samples) == 0) {
  stop("No matched samples remain after sample matching.")
}

# Preserve the sample order from the clinical table
meta_use <- meta[
  meta$sample %in% common_samples,
  ,
  drop = FALSE
]

count_mat <- count_mat[
  ,
  meta_use$sample,
  drop = FALSE
]

stopifnot(
  identical(
    meta_use$sample,
    colnames(count_mat)
  )
)

cat("\nMatched samples before covariate QC:", nrow(meta_use), "\n")


# =====================================================================
# 7. CALCULATE BMI
# =====================================================================

meta_use$age    <- suppressWarnings(as.numeric(meta_use$age))
meta_use$height <- suppressWarnings(as.numeric(meta_use$height))
meta_use$weight <- suppressWarnings(as.numeric(meta_use$weight))

# Height handling:
# - Values <= 3 are interpreted as meters (e.g. 1.62)
# - Values > 3 are interpreted as centimeters (e.g. 162)
meta_use$height_m <- ifelse(
  !is.na(meta_use$height) & meta_use$height > 3,
  meta_use$height / 100,
  meta_use$height
)

meta_use$BMI <- meta_use$weight / (meta_use$height_m ^ 2)

# Save BMI calculation table for record keeping
bmi_check <- meta_use[
  ,
  c(
    "sample",
    "cancer",
    "batch",
    "age",
    "height",
    "height_m",
    "weight",
    "BMI"
  ),
  drop = FALSE
]

write.csv(
  bmi_check,
  file.path(RESULT_DIR, "03_BMI_calculation_check.csv"),
  row.names = FALSE
)


# =====================================================================
# 8. CHECK REQUIRED COVARIATES
# =====================================================================

required_analysis_cols <- c(
  "sample",
  "cancer",
  "batch",
  "age",
  "BMI"
)

bad_covariate_rows <- !complete.cases(
  meta_use[, required_analysis_cols, drop = FALSE]
)

if (any(bad_covariate_rows)) {
  bad_samples <- meta_use$sample[bad_covariate_rows]

  write.csv(
    meta_use[
      bad_covariate_rows,
      required_analysis_cols,
      drop = FALSE
    ],
    file.path(
      RESULT_DIR,
      "ERROR_samples_with_missing_covariates.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Missing cancer/batch/age/BMI values were found. ",
    "Affected samples: ",
    paste(bad_samples, collapse = ", "),
    "\nSee result/ERROR_samples_with_missing_covariates.csv"
  )
}

# Basic plausibility warnings only.
# These DO NOT automatically remove samples.
if (any(meta_use$height_m < 1.2 | meta_use$height_m > 2.2)) {
  warning(
    "Some converted heights are outside 1.2-2.2 m. ",
    "Please inspect 03_BMI_calculation_check.csv."
  )
}

if (any(meta_use$BMI < 10 | meta_use$BMI > 60)) {
  warning(
    "Some BMI values are outside 10-60 kg/m^2. ",
    "Please inspect 03_BMI_calculation_check.csv."
  )
}


# =====================================================================
# 9. DEFINE GROUPS AND COVARIATES
# =====================================================================

expected_groups <- c(
  "Control",
  "CIN",
  "Cervical_cancer",
  "Ovarian_cancer",
  "Endometrial_cancer"
)

unexpected_groups <- setdiff(
  unique(meta_use$cancer),
  expected_groups
)

if (length(unexpected_groups) > 0) {
  stop(
    "Unexpected cancer group(s) found: ",
    paste(unexpected_groups, collapse = ", ")
  )
}

missing_groups <- setdiff(
  expected_groups,
  unique(meta_use$cancer)
)

if (length(missing_groups) > 0) {
  stop(
    "Required cancer group(s) missing: ",
    paste(missing_groups, collapse = ", ")
  )
}

# Categorical variables
meta_use$cancer <- factor(
  meta_use$cancer,
  levels = expected_groups
)

meta_use$batch <- factor(meta_use$batch)

# Optional batch reference level
if ("batch1" %in% levels(meta_use$batch)) {
  meta_use$batch <- relevel(
    meta_use$batch,
    ref = "batch1"
  )
}

# Continuous covariates:
# Scale = center to mean 0 AND divide by standard deviation.
#
# IMPORTANT:
# We DO NOT pre-correct the expression matrix.
# Age and BMI are only transformed as model covariates.
if (sd(meta_use$age) == 0) {
  stop("Age has zero standard deviation and cannot be scaled.")
}

if (sd(meta_use$BMI) == 0) {
  stop("BMI has zero standard deviation and cannot be scaled.")
}

meta_use$age_scaled <- as.numeric(
  scale(meta_use$age)
)

meta_use$BMI_scaled <- as.numeric(
  scale(meta_use$BMI)
)

# Confirm scaling
scaling_check <- data.frame(
  variable = c("age_scaled", "BMI_scaled"),
  mean = c(
    mean(meta_use$age_scaled),
    mean(meta_use$BMI_scaled)
  ),
  sd = c(
    sd(meta_use$age_scaled),
    sd(meta_use$BMI_scaled)
  )
)

cat("\n============================================================\n")
cat("SCALE CHECK\n")
cat("============================================================\n")
print(scaling_check)

write.csv(
  scaling_check,
  file.path(RESULT_DIR, "04_scaling_check.csv"),
  row.names = FALSE
)


# =====================================================================
# 10. SAVE FINAL SAMPLE INFORMATION USED IN DESEQ2
# =====================================================================

sample_info_used <- meta_use[
  ,
  c(
    "sample",
    "cancer",
    "batch",
    "age",
    "age_scaled",
    "height",
    "height_m",
    "weight",
    "BMI",
    "BMI_scaled"
  ),
  drop = FALSE
]

write.csv(
  sample_info_used,
  file.path(RESULT_DIR, "05_samples_used_for_DESeq2.csv"),
  row.names = FALSE
)

group_number <- as.data.frame(
  table(meta_use$cancer),
  stringsAsFactors = FALSE
)

colnames(group_number) <- c(
  "cancer",
  "n"
)

cat("\n============================================================\n")
cat("FINAL SAMPLE NUMBERS USED\n")
cat("============================================================\n")
print(group_number)

write.csv(
  group_number,
  file.path(RESULT_DIR, "06_sample_number_by_group.csv"),
  row.names = FALSE
)

# Batch distribution
batch_table <- as.data.frame.matrix(
  table(
    meta_use$cancer,
    meta_use$batch
  )
)

batch_table$cancer <- rownames(batch_table)
rownames(batch_table) <- NULL

batch_table <- batch_table[
  ,
  c(
    "cancer",
    setdiff(colnames(batch_table), "cancer")
  ),
  drop = FALSE
]

write.csv(
  batch_table,
  file.path(RESULT_DIR, "07_batch_distribution_by_group.csv"),
  row.names = FALSE
)

cat("\nBatch distribution by disease group:\n")
print(table(meta_use$cancer, meta_use$batch))


# =====================================================================
# 11. CHECK THE DESIGN MATRIX
# =====================================================================

# Teacher's requested concept:
# Do not separately correct expression values for covariates.
# Put batch + scaled age + scaled BMI directly into the differential
# expression design together with cancer group.

design_formula <- ~ batch + age_scaled + BMI_scaled + cancer

design_matrix <- model.matrix(
  design_formula,
  data = meta_use
)

design_rank <- qr(design_matrix)$rank

cat("\n============================================================\n")
cat("DESIGN MATRIX CHECK\n")
cat("============================================================\n")
cat(
  "Design formula: ~ batch + age_scaled + BMI_scaled + cancer\n"
)
cat("Design matrix columns:", ncol(design_matrix), "\n")
cat("Design matrix rank:", design_rank, "\n")

if (design_rank < ncol(design_matrix)) {
  stop(
    "The design matrix is not full rank. ",
    "There may be confounding among batch/cancer/covariates."
  )
}


# =====================================================================
# 12. GENE FILTER
# =====================================================================

# To avoid adding an arbitrary expression filter not requested by the teacher,
# only genes with zero counts in ALL retained samples are removed here.
keep_gene <- rowSums(count_mat) > 0

count_mat_use <- count_mat[
  keep_gene,
  ,
  drop = FALSE
]

cat("\nGenes before all-zero filtering:", nrow(count_mat), "\n")
cat("Genes entering DESeq2:", nrow(count_mat_use), "\n")


# =====================================================================
# 13. BUILD AND RUN DESEQ2
# =====================================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_mat_use,
  colData = meta_use,
  design = design_formula
)

cat("\n============================================================\n")
cat("RUNNING DESEQ2\n")
cat("============================================================\n")
cat(
  "Model: counts ~ batch + age_scaled + BMI_scaled + cancer\n"
)
cat(
  "No pre-correction of expression matrix is performed.\n\n"
)

dds <- DESeq(dds)

cat("\nDESeq2 fitting completed successfully.\n")


# =====================================================================
# 14. DEFINE FOUR DISEASE-VS-CONTROL COMPARISONS
# =====================================================================

comparisons <- c(
  CIN = "CIN",
  CC  = "Cervical_cancer",
  OC  = "Ovarian_cancer",
  EC  = "Endometrial_cancer"
)

LOG2FC_CUTOFF <- log2(FC_CUTOFF)

cat("\n============================================================\n")
cat("DEG THRESHOLD\n")
cat("============================================================\n")
cat("Adjusted P-value <", PADJ_CUTOFF, "\n")
cat("Absolute fold-change >=", FC_CUTOFF, "\n")
cat(
  "Equivalent absolute log2FoldChange >=",
  round(LOG2FC_CUTOFF, 4),
  "\n"
)


# =====================================================================
# 15. RUN CONTRASTS AND SAVE RESULTS
# =====================================================================

summary_list <- list()

for (short_name in names(comparisons)) {

  case_group <- comparisons[[short_name]]
  comparison_name <- paste0(
    short_name,
    "_vs_Control"
  )

  cat("\n------------------------------------------------------------\n")
  cat("Analyzing:", comparison_name, "\n")
  cat("------------------------------------------------------------\n")

  # Disease group divided by Control:
  # positive log2FoldChange = higher in disease group
  # negative log2FoldChange = lower in disease group
  res <- results(
    dds,
    contrast = c(
      "cancer",
      case_group,
      "Control"
    ),
    alpha = PADJ_CUTOFF
  )

  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)

  # Convert DESeq2 log2FC into ordinary case/control fold change.
  res_df$FC_case_over_control <- 2 ^ res_df$log2FoldChange

  # DEG definition
  significant <- (
    !is.na(res_df$padj) &
    res_df$padj < PADJ_CUTOFF &
    abs(res_df$log2FoldChange) >= LOG2FC_CUTOFF
  )

  up <- (
    significant &
    res_df$log2FoldChange >= LOG2FC_CUTOFF
  )

  down <- (
    significant &
    res_df$log2FoldChange <= -LOG2FC_CUTOFF
  )

  res_df$DEG_status <- "Not_significant"
  res_df$DEG_status[up]   <- "Up"
  res_df$DEG_status[down] <- "Down"

  # Move Gene to the first column
  res_df <- res_df[
    ,
    c(
      "Gene",
      setdiff(colnames(res_df), "Gene")
    ),
    drop = FALSE
  ]

  # Sort by adjusted P-value
  res_df <- res_df[
    order(
      res_df$padj,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  # Save all DESeq2 results
  all_result_file <- file.path(
    RESULT_DIR,
    paste0(
      comparison_name,
      "_ALL_results.csv"
    )
  )

  write.csv(
    res_df,
    all_result_file,
    row.names = FALSE
  )

  # Save significant DEG list
  deg_df <- res_df[
    res_df$DEG_status %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]

  deg_file <- file.path(
    RESULT_DIR,
    paste0(
      comparison_name,
      "_DEG_padj",
      PADJ_CUTOFF,
      "_FC",
      FC_CUTOFF,
      ".csv"
    )
  )

  write.csv(
    deg_df,
    deg_file,
    row.names = FALSE
  )

  # Summary
  summary_list[[short_name]] <- data.frame(
    comparison = comparison_name,
    disease_group_in_clinical_table = case_group,
    n_case = sum(meta_use$cancer == case_group),
    n_control = sum(meta_use$cancer == "Control"),
    padj_cutoff = PADJ_CUTOFF,
    fold_change_cutoff = FC_CUTOFF,
    abs_log2FC_cutoff = LOG2FC_CUTOFF,
    DEG_total = sum(significant),
    DEG_up = sum(up),
    DEG_down = sum(down),
    stringsAsFactors = FALSE
  )

  cat(
    "Total DEG:", sum(significant),
    " | Up:", sum(up),
    " | Down:", sum(down),
    "\n"
  )
}


# =====================================================================
# 16. FINAL DEG-NUMBER SUMMARY
# =====================================================================

deg_summary <- do.call(
  rbind,
  summary_list
)

rownames(deg_summary) <- NULL

summary_file <- file.path(
  RESULT_DIR,
  "FINAL_DEG_number_summary_padj0.1_FC1.5.csv"
)

write.csv(
  deg_summary,
  summary_file,
  row.names = FALSE
)

cat("\n\n")
cat("============================================================\n")
cat("FINAL DEG COUNTS\n")
cat("padj < ", PADJ_CUTOFF,
    " AND |FC| >= ", FC_CUTOFF, "\n", sep = "")
cat("============================================================\n")

print(
  deg_summary[
    ,
    c(
      "comparison",
      "n_case",
      "n_control",
      "DEG_total",
      "DEG_up",
      "DEG_down"
    )
  ]
)

cat("\nMain summary file:\n")
cat(summary_file, "\n")


# =====================================================================
# 17. SAVE SESSION INFORMATION FOR REPRODUCIBILITY
# =====================================================================

session_file <- file.path(
  RESULT_DIR,
  "99_sessionInfo.txt"
)

sink(session_file)
print(sessionInfo())
sink()

cat("\nSession information saved to:\n")
cat(session_file, "\n")


# =====================================================================
# 18. FINAL MESSAGE
# =====================================================================

cat("\n============================================================\n")
cat("ANALYSIS FINISHED\n")
cat("============================================================\n")
cat(
  "Please STOP here and first report the four DEG numbers to the teacher.\n"
)
cat(
  "Do not proceed to downstream analysis until the DEG cutoff is confirmed.\n"
)
cat("\nMost important file:\n")
cat(
  "result/FINAL_DEG_number_summary_padj0.1_FC1.5.csv\n"
)
cat("============================================================\n")
