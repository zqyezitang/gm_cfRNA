#################Code for ploting figures

######################### Figure 1
###Fig1b. details in clin_heatmap.R

###Fig1C. PCA plot of all samples
tpm_path <- "D:/项目/澳门理工大学宫颈癌cfRNA/分析结果/监测文章/Paper_sample_TPM.csv"
meta_path <- "D:/项目/澳门理工大学宫颈癌cfRNA/分析结果/监测文章/TableS1. Clinical information of samples.xlsx"

# 
outlier_samples <- c("2488-CQCU034A", "2546-CQCU092A", "MPU-CQCH1499A",
                     "2518-CQCU064A", "MPU-CQCU780A", "MPU-CQCU941A",
                     "MPU-CQCH425A", "MPU-CQCH543A")

#
tpm_data <- read_csv(tpm_path, show_col_types = FALSE)
gene_ids <- tpm_data[[1]]
expr_data <- as.data.frame(tpm_data[, -1])
rownames(expr_data) <- gene_ids
colnames(expr_data) <- colnames(tpm_data)[-1]
expr_data <- expr_data[, !colnames(expr_data) %in% outlier_samples]

expr_t <- as.data.frame(t(expr_data))
expr_t$sample <- rownames(expr_t)

#
meta <- read_excel(meta_path)
meta <- meta %>% filter(sample %in% expr_t$sample)
merged <- left_join(expr_t, meta, by = "sample")
merged

dir.create("PCA_results", showWarnings = FALSE)

run_pca_plot <- function(df, group_col, shape_col = NULL, title, filename) {
  rownames(df) <- df$sample
  expr_mat <- df %>% select(where(is.numeric)) %>% as.matrix()
  expr_mat_log <- log2(expr_mat + 1)
  expr_mat_log <- expr_mat_log[, apply(expr_mat_log, 2, function(x) sd(x, na.rm = TRUE) > 0)]
  
  pca_res <- prcomp(expr_mat_log, center = TRUE, scale. = TRUE)
  var_exp <- summary(pca_res)$importance[2, 1:2]
  
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$sample <- rownames(pca_df)
  
  meta_cols <- df[, setdiff(colnames(df), colnames(expr_mat_log)), drop = FALSE]
  meta_cols$sample <- rownames(meta_cols)
  pca_df <- left_join(pca_df, meta_cols, by = "sample")
  write.table(pca_df,file='PCA.txt',sep='\t')
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
    geom_point(aes(shape = if (!is.null(shape_col)) .data[[shape_col]] else NULL),
               size = 2.5, alpha = 0.85) +
    scale_color_manual(values = c("#4DAF4A", "#377EB8","#E41A1C","#984EA3","#FFFF33")) +
    labs(title = title,
         x = paste0("PC1 (", round(var_exp[1] * 100, 1), "% variance)"),
         y = paste0("PC2 (", round(var_exp[2] * 100, 1), "% variance)")) +
    theme_minimal(base_size = 12)
  
  ggsave(filename = file.path("PCA_results", filename),
         plot = p, width = 8, height = 5)
}

run_pca_plot(df=merged, group_col = 'cancer', shape_col = "treat_time", title='', file='PCA_total.pdf')


###Figure1d. details in GSEA.R


############################### Figure 2
#fig2C legend for GO
library('RColorBrewer')
colors <- c(brewer.pal(9,'Set1'),brewer.pal(12,'Set3'))
name <- c("cluster1",'cluster2','cluster3','cluster4','cluster5')
pdf('Fig2c.cluster_legend.pdf')
plot.new()
legend('topleft',pch=15,col=colors,legend = name,cex=0.8,bty='n')
dev.off()

#########################Figure 4I-H
tissue_exp <- data.frame(read.table('D:/项目/澳门理工大学宫颈癌cfRNA/分析结果/监测文章/HPA/tissue_sum_exp.txt',header = T,row.names = 1,sep="\t",check.names = F))
tissue_exp <- tissue_exp[meta$sample,]  #
tissue_exp$cancer=meta$cancer
tissue_exp$treattime=meta$treat_time


#diff between groups
#### DE of scores between groups
sum_control <- tissue_exp[tissue_exp$cancer=="Control",]
sum_CC <- tissue_exp[(tissue_exp$cancer=="Cervical_cancer" & tissue_exp$treattime=="pretherapy"),]
sum_OV <- tissue_exp[(tissue_exp$cancer=="Ovarian_cancer" & tissue_exp$treattime=="pretherapy"),]
sum_EM <- tissue_exp[(tissue_exp$cancer=="Endometrial_cancer" & tissue_exp$treattime=="pretherapy"),]
sum_CIN <- tissue_exp[(tissue_exp$cancer=="CIN"),]
p_values <- c()
labels <- c()

##CC
#for (i in seq(3,33))
for (i in seq(1,12))  
{
  p <- wilcox.test(as.numeric(sum_CC[,i]),as.numeric(sum_control[,i]))$p.value
  print(colnames(sum_CC)[i])
  print(p)
  p_values <- c(p_values,p)
  labels <- c(labels,paste0("CC_",colnames(sum_CC)[i]))
}
##CIN
#for (i in seq(1,12))  
#{
#  p <- wilcox.test(as.numeric(sum_CIN[,i]),as.numeric(sum_control[,i]))$p.value
#  print(colnames(sum_CIN)[i])
#  print(p)
#}
##OV
for (i in seq(1,12))  
{
  p <- wilcox.test(as.numeric(sum_OV[,i]),as.numeric(sum_control[,i]))$p.value
  print(colnames(sum_OV)[i])
  print(p)
  p_values <- c(p_values,p)
  labels <- c(labels,paste0("OC_",colnames(sum_CC)[i]))
}
##EM
for (i in seq(1,12))  
{
  p <- wilcox.test(as.numeric(sum_EM[,i]),as.numeric(sum_control[,i]))$p.value
  print(colnames(sum_EM)[i])
  print(p,i)
  p_values <- c(p_values,p)
  labels <- c(labels,paste0("EC_",colnames(sum_CC)[i]))
}

padj <- p.adjust(p_values, method = "BH")


#plots of diff tissue, p<0.01
Cervix_CC<-data.frame(score=c(sum_control[,3],sum_CC[,3],sum_EM[,3],sum_OV[,3]),
                      cancer=as.factor(c(rep('Control',27),rep('CC',40),rep('EC',28),rep('OC',36))))
Ovary_OV <- data.frame(score=c(sum_control[,10],sum_CC[,10],sum_EM[,10],sum_OV[,10]),
                       cancer=as.factor(c(rep('Control',27),rep('CC',40),rep('EC',28),rep('OC',36))))
bone_OV <- data.frame(score=c(sum_control[,2],sum_CC[,2],sum_EM[,2],sum_OV[,2]),
                      cancer=as.factor(c(rep('Control',27),rep('CC',40),rep('EC',28),rep('OC',36))))
Endo_EM <- data.frame(score=c(sum_control[,4],sum_CC[,4],sum_EM[,4],sum_OV[,4]),
                      cancer=as.factor(c(rep('Control',27),rep('CC',40),rep('EC',28),rep('OC',36))))

pdf('HPA_tissue_boxplot.pdf')
y_max <- max(Ovary_OV$score)
y_pos <- y_max * 1.05
ggplot(Ovary_OV, aes(x = factor(cancer,level=c('Control','CC','EC','OC')), y = score)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(fill = factor(cancer,level=c('Control','CC','EC','OC')))) + # 先画小提琴图展示密度
  geom_boxplot(width = 0.1, outlier.shape = NA) +                   
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +               
  labs(x = "", y = "Score of tissue",
       title = "Ovary",fill=NULL) +
  theme_classic() +
  #  geom_segment(aes(x = 1, xend = 4, y = y_pos* 1.08, yend = y_pos* 1.08)) +
  # 
  annotate("text", x = 4, y = y_pos * 1.1, label = "*", size = 6,col='red')
ggplot(Endo_EM, aes(x = factor(cancer,level=c('Control','CC','EC','OC')), y = score)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(fill = factor(cancer,level=c('Control','CC','EC','OC')))) + # 
  geom_boxplot(width = 0.1, outlier.shape = NA) +                   # 
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +               # 
  labs(x = "", y = "Score of tissue",
       title = "Endometrium",fill=NULL) +
  theme_classic() 
ggplot(Cervix_CC, aes(x = factor(cancer,level=c('Control','CC','EC','OC')), y = score)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(fill = factor(cancer,level=c('Control','CC','EC','OC')))) + # 
  geom_boxplot(width = 0.1, outlier.shape = NA) +                   # 
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +               # 
  labs(x = "", y = "Score of tissue",
       title = "Cervix",fill=NULL) +
  theme_classic() 
ggplot(bone_OV, aes(x = factor(cancer,level=c('Control','CC','EC','OC')), y = score)) +
  geom_violin(trim = FALSE, alpha = 0.4, aes(fill = factor(cancer,level=c('Control','CC','EC','OC')))) + # 
  geom_boxplot(width = 0.1, outlier.shape = NA) +                   # 
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +               # 
  labs(x = "", y = "Score of tissue",
       title = "Bone marrow",fill=NULL) +
  theme_classic() +
  annotate("text", x = 4, y = max(bone_OV$score) * 1.1, label = "*", size = 6,col='red')
dev.off()
