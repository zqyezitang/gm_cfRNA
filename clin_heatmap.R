# 加载包
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)

clinical <- data.frame(
  sample <- meta$sample,
  Age <- meta$age,
  BMI <- as.numeric(meta$weight)/as.numeric(meta$height)^2,
  group <- meta$cancer,
  Stage <- meta$stage_raw,
  Tumor_marker <- meta$tumor_marker_status,
  HPV <- meta$HPV,
  Differentiation <- meta$differentiation,
  Molecular_subtype <- meta$Molecular_subtype,
  BRCA_mutation <- meta$BRCA_status,
  stringsAsFactors = FALSE
)
colnames(clinical) <- c('sample','Age','BMI','group','Stage','Tumor_marker','HPV','Differentiation','Molecular_subtype','BRCA_mutation')

# ==========  ==========
continuous_vars <- c("Age", "BMI")
categorical_vars <- c("Stage", "Tumor_marker", "HPV", 
                      "Differentiation", "Molecular_subtype", "BRCA_mutation")

# ==========  ==========
cont_matrix <- t(clinical[, continuous_vars, drop = FALSE])
colnames(cont_matrix) <- clinical$sample

cat_matrix_raw <- t(clinical[, categorical_vars, drop = FALSE])
colnames(cat_matrix_raw) <- clinical$sample

# ==========  ==========
category_color_maps <- list()

stage_levels <- c("I", "II", "III", "IV")
stage_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D")
names(stage_colors) <- stage_levels
category_color_maps[["Stage"]] <- stage_colors

tm_levels <- c("no", "yes")
tm_colors <- c("white", "#3182BD")
names(tm_colors) <- tm_levels
category_color_maps[["Tumor_marker"]] <- tm_colors

hpv_levels <- c("no", "yes")
hpv_colors <- c("white", "#31A354")
names(hpv_colors) <- hpv_levels
category_color_maps[["HPV"]] <- hpv_colors

diff_levels <- c("low", "mid", "high")
diff_colors <- c("#FEE391", "#FEC44F", "#D95F0E")
names(diff_colors) <- diff_levels
category_color_maps[["Differentiation"]] <- diff_colors

brca_levels <- c("no", "yes")
brca_colors <- c("white", "#54278F")
names(brca_colors) <- brca_levels
category_color_maps[["BRCA_mutation"]] <- brca_colors

subtype_levels <- c("POLEmut", "p53abn", "NSMP", "MMRd","pMMR")
#subtype_colors <- c("#FBC4C4", "#F48CB6", "#C77CBF", "#7B3294","#54278F")
subtype_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
names(subtype_colors) <- subtype_levels
category_color_maps[["Molecular_subtype"]] <- subtype_colors

group_colors <- c("Control" = "#E41A1C", "CIN" = "#377EB8","Cervical_cancer"="#4DAF4A",
                  "Endometrial_cancer"="#984EA3","Ovarian_cancer"="#FFFF33")

#
col_annotation <- HeatmapAnnotation(
  Group = clinical$group,
  col = list(Group = group_colors),
  annotation_name_side = "left",  
  annotation_legend_param = list(
    Group = list(title = "Sample Group", at = c("Control", "CIN","Cervical_cancer","Endometrial_cancer","Ovarian_cancer"))
  ),
  show_annotation_name = TRUE,
  gp = gpar(col = "black"),
  border = TRUE
)

# 
heatmap_list <- list()
row_height <- unit(0.4, "cm")  

# 
for (var in continuous_vars) {
  values <- as.numeric(cont_matrix[var, ])
  min_val <- min(values, na.rm = TRUE)
  max_val <- max(values, na.rm = TRUE)
  
  mat <- matrix(values, nrow = 1, ncol = length(values))
  rownames(mat) <- var
  colnames(mat) <- colnames(cont_matrix)
  
  if (min_val == max_val) {
    col_fun <- colorRamp2(c(min_val, max_val), c("white", "gray"))
  } else {
    col_fun <- colorRamp2(
      c(min_val, (min_val + max_val)/2, max_val),
      c("#4575B4", "#FFFFBF", "#D73027")
    )
  }
  
  hm <- Heatmap(
    mat,
    name = var,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = (var == continuous_vars[1]),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(title = var, direction = "horizontal"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    na_col = "gray",
    height = row_height,
    top_annotation = if (var == continuous_vars[1]) col_annotation else NULL  # 只在第一个热图添加注释
  )
  heatmap_list[[var]] <- hm
}

# 
for (var in categorical_vars) {
  raw_values <- cat_matrix_raw[var, ]
  level_order <- names(category_color_maps[[var]])
  
  values_numeric <- rep(NA, length(raw_values))
  for (i in seq_along(level_order)) {
    values_numeric[raw_values == level_order[i]] <- i
  }
  
  mat <- matrix(values_numeric, nrow = 1, ncol = length(values_numeric))
  rownames(mat) <- var
  colnames(mat) <- colnames(cat_matrix_raw)
  
  col_fun <- colorRamp2(seq_along(level_order), category_color_maps[[var]])
  
  hm <- Heatmap(
    mat,
    name = var,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      title = var,
      at = seq_along(level_order),
      labels = level_order,
      direction = "horizontal"
    ),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    na_col = "gray90",
    height = row_height
  )
  heatmap_list[[var]] <- hm
}

# ==========  ==========
combined_heatmap <- NULL
for (hm in heatmap_list) {
  if (is.null(combined_heatmap)) {
    combined_heatmap <- hm
  } else {
    combined_heatmap <- combined_heatmap %v% hm
  }
}

# ==========  ==========
pdf('clinical_heatmap.pdf',width = 9,height = 4)
draw(combined_heatmap, 
     merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()
