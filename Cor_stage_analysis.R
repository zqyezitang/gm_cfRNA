#########################COR WITH LM22 SCORE
library(ggplot2)
library(reshape2)
library(corrplot)
library(ggpubr)
library(pheatmap)

outlier_samples <- c("2488-CQCU034A", "2546-CQCU092A", "MPU-CQCH1499A", "2518-CQCU064A", "MPU-CQCU780A", "MPU-CQCU941A", "MPU-CQCH425A", "MPU-CQCH543A")
setwd('D:/项目/澳门理工大学宫颈癌cfRNA/分析结果/监测文章/Origin_LM22')
origin_path <- 'CIBERSORTx_Job7_Results.csv'
origin_data <- read_csv(origin_path, show_col_types = FALSE)
#count_data <- read_csv(tpm_path, show_col_types = FALSE)
sample <- origin_data[[1]]
origin_mat <- data.frame(origin_data[, -1])
#count_mat <- data.frame(log2(count_data[, -1]+1))
rownames(origin_mat) <- sample
colnames(origin_mat) <- colnames(origin_data)[-1]
origin_mat <- origin_mat[, !(colnames(origin_mat) %in% outlier_samples)]

meta_path <- "D:/项目/澳门理工大学宫颈癌cfRNA/分析结果/监测文章/TableS1. Clinical information of samples.xlsx"
meta <- read_excel(meta_path)
meta <- meta %>% filter(sample %in% rownames(origin_mat))
meta$sample <- as.character(meta$sample)
origin_mat <- origin_mat[meta$sample,]

cancer <- origin_mat[,1:22]
cancer$Stage <- meta$stage_raw
cancer$cancer <- meta$cancer 
cancer$age <- as.numeric(meta$age)
cancer$BMI <- as.numeric(meta$weight)/as.numeric(meta$height)^2
cancer$ANC <- as.numeric(meta$ANC)
cancer$ALC <- as.numeric(meta$ALC)
cancer$HPV <- meta$HPV
cancer$tumor_marker_status <- meta$tumor_marker_status
cancer$differentiation <- as.factor(meta$differentiation)
cancer$BRCA_status <- meta$BRCA_status
cancer$Molecular_subtype <- meta$Molecular_subtype

CIN_data <- cancer[cancer$cancer=="CIN",]
CC_data <- cancer[cancer$cancer=="Cervical_cancer",]
OC_data <- cancer[cancer$cancer=="Ovarian_cancer",]
EC_data <- cancer[cancer$cancer=="Endometrial_cancer",]


# 
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggpubr")) install.packages("ggpubr")
library(tidyverse)
library(ggpubr)

# 
analyze_stage_correlation <- function(data, stage_col = "Stage", cancer,
                                      method = "spearman") {
  
    stage_numeric <- data %>%
    mutate(Stage_num = case_when(
      Stage == "I" ~ 1,
      Stage == "II" ~ 2,
      Stage == "III" ~ 3,
      Stage == "IV" ~ 4,
      TRUE ~ NA_real_
    ))
  
 
  score_cols <- data %>%
    select(-all_of(stage_col)) %>%
    select(where(is.numeric)) %>%
    names()
  
  correlation_results <- map_df(score_cols, function(score) {
    test <- cor.test(stage_numeric[[score]], 
                     stage_numeric$Stage_num, 
                     method = method,
                     use = "complete.obs")
    
    tibble(
      Score = score,
      Correlation = test$estimate,
      P_value = test$p.value,
      Significance = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01 ~ "**",
        test$p.value < 0.05 ~ "*"
     #   TRUE ~ "ns"
      )
    )
  })
  
 p <- ggplot(correlation_results, 
              aes(x = reorder(Score, Correlation), 
                  y = Correlation,
                  fill = ifelse(Correlation > 0, "Positive", "Negative"))) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = Significance, y = Correlation + sign(Correlation) * 0.05),
              size = 4, vjust = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Positive" = "#E64B35", "Negative" = "#4DBBD5"),
                      name = "") +
    labs(title = "",
         subtitle = cancer,
         x = "Immune signature",
         y = paste0("Correlation : (", toupper(method), ")")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right") +
    coord_flip()  
  
  # 
  print(correlation_results)
  
  # 
  return(list(results = correlation_results, plot = p))
}

# 
 result <- analyze_stage_correlation(CC_data, stage_col = "Stage",cancer='CC')
 pdf('CC_LM22_stage_cor_bar.pdf')
 print(result$plot)  # 
 dev.off()
 result <- analyze_stage_correlation(OC_data, stage_col = "Stage",cancer='OC')
 pdf('OC_LM22_stage_cor_bar.pdf')
 print(result$plot)  
 dev.off()
 result <- analyze_stage_correlation(EC_data, stage_col = "Stage",cancer='EC')
 pdf('EC_LM22_stage_cor_bar.pdf')
 print(result$plot)  
 dev.off()
  
 
 
 #####################boxplot of significant
 boxplot(`Macrophages M1`~Stage,data=CC_data)
 boxplot(`Dendritic cells resting`~Stage,data=OC_data)
 boxplot(`Macrophages M1`~Stage,data=EC_data)
 
 #################### plots for manuscript
 library(tidyverse)
 library(ggpubr)   # 
 library(rstatix)  # 
 
 # 
plot_data <- CC_data %>%
   filter(Stage %in% c("I", "II", "III", "IV")) %>%
   mutate(Stage = factor(Stage, levels = c("I", "II", "III", "IV")))
score <- "Macrophages M1"   
 
plot_data <- OC_data %>%
  filter(Stage %in% c("I", "II", "III", "IV")) %>%
  mutate(Stage = factor(Stage, levels = c("I", "II", "III", "IV")))
score <- "Dendritic cells resting"   

plot_data <- EC_data %>%
  filter(Stage %in% c("I", "II", "III", "IV")) %>%
  mutate(Stage = factor(Stage, levels = c("I", "II", "III", "IV")))
score <- "Macrophages M1"   


 # 2.  Wilcoxon  + BH 
stat_adj <- plot_data %>%
   pairwise_wilcox_test(
     as.formula(paste0("`", score, "` ~ Stage")),
     comparisons = list(c("I","II"), c("II","III"), c("III","IV")),
     p.adjust.method = "BH"
   ) %>%
   add_significance("p.adj")     # 
 
print(stat_adj)
 
 # 3.plots
pdf('CC_sig_cell_boxplot_by_stage.pdf')
pdf('OC_sig_cell_boxplot_by_stage.pdf')
pdf('EC_sig_cell_boxplot_by_stage.pdf')
ggplot(plot_data, aes(x = Stage, y = .data[[score]], fill = Stage)) +
   geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
   geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
   stat_pvalue_manual(
     stat_adj,
     label    = "p.adj.signif",   # 
     y.position = max(plot_data[[score]], na.rm = TRUE) *
       c(1.02, 1.08, 1.14),   # 
     tip.length = 0.01,
     size = 4
   ) +
   scale_fill_manual(values = c("#4DBBD5","#00A087","#E64B35","#3C5488")) +
   labs(title = score, x = "Stage", y = score) +
   theme_classic(base_size = 13) +
   theme(legend.position = "none",
         plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off() 
 
 
 
 ###################### cor with clinical parameter
 calculate_clinical_correlations <- function(data, 
                                             clinical_vars,  # 
                                             score_vars,     # 
                                             method = "spearman") {
   
   # 
   results <- data.frame()
   
   # 
   for(clinical in clinical_vars) {
     for(score in score_vars) {
       # 
       if(sum(!is.na(data[[clinical]])) < 3 | sum(!is.na(data[[score]])) < 3) next
       
       # 
       test <- cor.test(data[[clinical]], data[[score]], 
                        method = method, 
                        use = "complete.obs")
       
       # 
       results <- rbind(results, data.frame(
         Clinical_Var = clinical,
         Score_Var = score,
         Correlation = test$estimate,
         P_value = test$p.value,
         Significance = case_when(
           test$p.value < 0.001 ~ "***",
           test$p.value < 0.01 ~ "**",
           test$p.value < 0.05 ~ "*",
           TRUE ~ "ns"
         ),
         Method = method
       ))
     }
   }
   
   # 
   results$FDR <- p.adjust(results$P_value, method = "fdr")
   results$Significance_FDR <- case_when(
     results$FDR < 0.001 ~ "***",
     results$FDR < 0.01 ~ "**",
     results$FDR < 0.05 ~ "*",
     TRUE ~ "ns"
   )
   
   return(results)
 }
 
 ###
 plot_correlation_heatmap <- function(correlation_results, 
                                      title = "") {
   
   # 
   heatmap_data <- correlation_results %>%
     select(Clinical_Var, Score_Var, Correlation) %>%
     pivot_wider(names_from = Score_Var, values_from = Correlation) %>%
     column_to_rownames("Clinical_Var")
   
   # 
   sig_matrix <- correlation_results %>%
     select(Clinical_Var, Score_Var, Significance) %>%
     pivot_wider(names_from = Score_Var, values_from = Significance) %>%
     column_to_rownames("Clinical_Var")
   
   # 
   corrplot(as.matrix(heatmap_data), 
            method = "color",
            type = "full",
            order = "hclust",
            tl.col = "black",
            tl.cex = 0.8,
            cl.cex = 0.8,
            addCoef.col = "black",
            number.cex = 0.7,
            p.mat = as.matrix(correlation_results %>%
                                select(Clinical_Var, Score_Var, P_value) %>%
                                pivot_wider(names_from = Score_Var, values_from = P_value) %>%
                                column_to_rownames("Clinical_Var")),
            sig.level = 0.05,
            insig = "label_sig",
            title = title,
            mar = c(0,0,2,0))
 }
 
 # 3.2 
 plot_correlation_bubble <- function(correlation_results, 
                                     p_threshold = 0.05) {
   
   # 
   sig_results <- correlation_results %>%
     filter(P_value < p_threshold)
   
   if(nrow(sig_results) == 0) {
     stop("no significant result")
   }
   
   # 
   p <- ggplot(sig_results, aes(x = Score_Var, y = Clinical_Var)) +
     geom_point(aes(size = abs(Correlation), 
                    color = Correlation,
                    alpha = 0.7)) +
     scale_color_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                           midpoint = 0, name = "Correlation") +
     scale_size_continuous(range = c(3, 12), name = "|Correlation|") +
     geom_text(aes(label = Significance), vjust = 0.8, hjust = 1.2, size = 5) +
     labs(title = "",
          x = "",
          y = "") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           panel.grid.minor = element_blank(),
           legend.position = "right")
   
   return(p)
 }
 
 # 
 plot_significant_scatters <- function(data, 
                                       correlation_results,
                                       p_threshold = 0.05) {
   
   # 
   sig_results <- correlation_results %>%
     filter(P_value < p_threshold)
   
   if(nrow(sig_results) == 0) {
     return(NULL)
   }
   
   # 
   plots <- list()
   
   for(i in 1:nrow(sig_results)) {
     clinical <- sig_results$Clinical_Var[i]
     score <- sig_results$Score_Var[i]
     corr <- round(sig_results$Correlation[i], 3)
     p_val <- sig_results$P_value[i]
     
     # 
     plot_data <- data[, c(clinical, score)] %>% drop_na()
     
     # 
     if(is.numeric(plot_data[[clinical]])) {
       # 
       p <- ggplot(plot_data, aes(x = .data[[clinical]], y = .data[[score]])) +
         geom_point(alpha = 0.6, size = 3, color = "#4DBBD5") +
         geom_smooth(method = "lm", se = TRUE, color = "#E64B35", alpha = 0.2) +
         stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
         labs(title = paste0(score, " vs ", clinical),
              subtitle = paste0("r = ", corr, ", p = ", format(p_val, scientific = TRUE, digits = 3)),
              x = clinical,
              y = score) +
         theme_bw() +
         theme(plot.title = element_text(hjust = 0.5, face = "bold"))
       
     } else {
       # 
       p <- ggplot(plot_data, aes(x = .data[[clinical]], y = .data[[score]], 
                                  fill = .data[[clinical]])) +
         geom_boxplot(alpha = 0.7) +
         geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
         stat_compare_means(method = "wilcox.test", label = "p.signif") +
         labs(title = '',
              subtitle = paste0("p = ", format(p_val, scientific = TRUE, digits = 3)),
              x = clinical,
              y = score) +
         theme_bw() +
         theme(plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.position = "none")
     }
     
     plots[[i]] <- p
   }
   
   return(plots)
 }
 

#################
analyze_mixed_clinical_data <- function(data, 
                                        continuous_vars,   # 
                                        categorical_vars,  # 
                                        score_vars) {      # 
  
  results_list <- list()
  
  # 1. 
  if(length(continuous_vars) > 0) {
    continuous_results <- calculate_clinical_correlations(
      data, continuous_vars, score_vars, method = "spearman"
    )
    results_list$continuous <- continuous_results
  }
  
  # 2. 
  if(length(categorical_vars) > 0) {
    categorical_results <- data.frame()
    
    for(cat_var in categorical_vars) {
      for(score in score_vars) {
        # 
        plot_data <- data[, c(cat_var, score)] %>% drop_na()
        
        # 
        n_groups <- length(unique(plot_data[[cat_var]]))
        
        if(n_groups == 2) {
          # 
          test <- wilcox.test(plot_data[[score]] ~ plot_data[[cat_var]])
          effect_size <- abs(qnorm(test$p.value/2)) / sqrt(nrow(plot_data))
          test_name <- "Wilcoxon"
       #   test <- wilcox.test(as.formula(paste(score, "~", cat_var)), 
      #                        data = plot_data)
      #    effect_size <- abs(qnorm(test$p.value/2)) / sqrt(nrow(plot_data))
      #    test_name <- "Wilcoxon"
        } else {
          # 
          test <- kruskal.test(plot_data[[score]] ~ plot_data[[cat_var]])
         # test <- kruskal.test(as.formula(paste(score, "~", cat_var)), 
         #                      data = plot_data)
          # 
          effect_size <- test$statistic / (nrow(plot_data) - 1)
          test_name <- "Kruskal-Wallis"
        }
        
        #
        categorical_results <- rbind(categorical_results, data.frame(
          Clinical_Var = cat_var,
          Score_Var = score,
          Test = test_name,
          Statistic = test$statistic,
          P_value = test$p.value,
          Effect_Size = effect_size,
          Significance = case_when(
            test$p.value < 0.001 ~ "***",
            test$p.value < 0.01 ~ "**",
            test$p.value < 0.05 ~ "*",
            TRUE ~ "ns"
          )
        ))
      }
    }
    
    # 
    categorical_results$FDR <- p.adjust(categorical_results$P_value, method = "fdr")
    results_list$categorical <- categorical_results
  }
  
  return(results_list)
}

#

 
clinical_categorical <- c("Stage", "Gender", "Smoking")


clinical_continuous <- c("age", "BMI","ANC","ALC")
clinical_categorical <- c("Stage", "HPV", "tumor_marker_status")
score_vars <- colnames(CC_data)[1:22]

CC_corr <- calculate_clinical_correlations(
  CC_data, 
  clinical_continuous, 
  score_vars,
  method = "spearman"
)
CC_corr_mix <- analyze_mixed_clinical_data(
  CC_data, 
  clinical_continuous,
  clinical_categorical,
  score_vars
)
EC_corr <- calculate_clinical_correlations(
  EC_data, 
  clinical_vars, 
  score_vars,
  method = "spearman"
)
OC_corr <- calculate_clinical_correlations(
  OC_data, 
  clinical_vars, 
  score_vars,
  method = "spearman"
)
pdf('CC_clin_cor_LM22.pdf')
bubble_plot <- plot_correlation_bubble(CC_corr, p_threshold = 0.05)
print(bubble_plot)

# 
scatter_plots <- plot_significant_scatters(CC_data, CC_corr, p_threshold = 0.05)
if(!is.null(scatter_plots)) {
  # 
  ggarrange(plotlist = scatter_plots[1:min(4, length(scatter_plots))], 
            ncol = 2, nrow = 2)
}
dev.off() 
pdf('OC_clin_cor_LM22.pdf')
bubble_plot <- plot_correlation_bubble(OC_corr, p_threshold = 0.05)
print(bubble_plot)

# 
scatter_plots <- plot_significant_scatters(OC_data, OC_corr, p_threshold = 0.05)
if(!is.null(scatter_plots)) {
  # 
  ggarrange(plotlist = scatter_plots[1:min(8, length(scatter_plots))], 
            ncol = 2, nrow = 2)
}
dev.off() 
pdf('EC_clin_cor_LM22.pdf')
bubble_plot <- plot_correlation_bubble(EC_corr, p_threshold = 0.05)
print(bubble_plot)

# 
scatter_plots <- plot_significant_scatters(EC_data, EC_corr, p_threshold = 0.05)
if(!is.null(scatter_plots)) {
  # 
  ggarrange(plotlist = scatter_plots[1:min(8, length(scatter_plots))], 
            ncol = 2, nrow = 2)
}
dev.off() 

mixed_results <- analyze_mixed_clinical_data(
  example_data,
  clinical_continuous,
  clinical_categorical,
  score_vars
)

print(CC_corr_mix$categorical %>% filter(P_value < 0.05))


write.csv(continuous_corr, "clinical_score_correlations.csv", row.names = FALSE)
write.csv(mixed_results$categorical, "categorical_analysis_results.csv", row.names = FALSE)




plot_group_comparison_heatmap <- function(data, clinical_vars, score_vars) {
  
  effect_matrix <- matrix(NA, nrow = length(clinical_vars), ncol = length(score_vars))
  p_matrix <- matrix(NA, nrow = length(clinical_vars), ncol = length(score_vars))
  
  for(i in seq_along(clinical_vars)) {
    for(j in seq_along(score_vars)) {
      clinical <- clinical_vars[i]
      score <- score_vars[j]
      
      plot_data <- data[, c(clinical, score)] %>% drop_na()
      
      if(length(unique(plot_data[[clinical]])) == 2) {
        test <- wilcox.test(plot_data[[score]] ~ plot_data[[clinical]])
        # 
        effect_matrix[i, j] <- abs(test$statistic / (nrow(plot_data)^2) - 0.5) * 2
        p_matrix[i, j] <- test$p.value
      } else if(length(unique(plot_data[[clinical]])) > 2) {
        # 
        test <- kruskal.test(plot_data[[score]] ~ plot_data[[clinical]])
        effect_matrix[i, j] <- test$statistic / (nrow(plot_data) - 1)
        p_matrix[i, j] <- test$p.value
      }
    }
  }
  
  rownames(effect_matrix) <- clinical_vars
  colnames(effect_matrix) <- score_vars
  
  # 
  sig_matrix <- matrix("", nrow = length(clinical_vars), ncol = length(score_vars))
  sig_matrix[p_matrix < 0.05] <- "*"
  sig_matrix[p_matrix < 0.01] <- "**"
  sig_matrix[p_matrix < 0.001] <- "***"
  
  # 
  library(reshape2)
  library(scales)
  
  plot_data <- melt(effect_matrix) %>%
    rename(Clinical = Var1, Score = Var2, Effect_Size = value) %>%
    mutate(Significance = as.vector(sig_matrix))
  
  ggplot(plot_data, aes(x = Score, y = Clinical, fill = Effect_Size)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = paste0(round(Effect_Size, 2), Significance)), 
              size = 1, color = "black") +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                         midpoint = 0.1, name = "") +
    labs(title = "",
         x = "",
         y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# 
 heatmap_plot <- plot_group_comparison_heatmap(CC_data, 
                                               c("Stage", "HPV", "tumor_marker_status"), 
                                               score_vars)
pdf('Heatmap_CC_cor.pdf')
print(heatmap_plot)
dev.off()
heatmap_plot <- plot_group_comparison_heatmap(EC_data, 
                                              c("Stage", "tumor_marker_status","Molecular_subtype","differentiation"), 
                                              score_vars)
pdf('Heatmap_EC_cor.pdf')
print(heatmap_plot)
dev.off()
heatmap_plot <- plot_group_comparison_heatmap(OC_data, 
                                              c("Stage", "tumor_marker_status","BRCA_status","differentiation"), 
                                              score_vars)
pdf('Heatmap_OC_cor.pdf')
print(heatmap_plot)
dev.off()


#### box for each group
plot_lollipop_summary <- function(data, clinical_var, score_vars) {
  
  # 
  summary_data <- data %>%
    select(all_of(c(clinical_var, score_vars))) %>%
    group_by(.data[[clinical_var]]) %>%
    summarise(across(all_of(score_vars), 
                     list(median = ~median(., na.rm = TRUE),
                          q1 = ~quantile(., 0.25, na.rm = TRUE),
                          q3 = ~quantile(., 0.75, na.rm = TRUE)),
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(cols = starts_with(score_vars),
                 names_to = c("Score", ".value"),
                 names_pattern = "(.*)_(median|q1|q3)")

  # 
  ggplot(summary_data, aes(x = .data[[clinical_var]], 
                           y = median, 
                           color = .data[[clinical_var]])) +
    geom_linerange(aes(ymin = q1, ymax = q3), size = 1.5, alpha = 0.5) +
    geom_point(size = 4) +
    facet_wrap(~Score, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = paste0(clinical_var, " "),
         x = clinical_var,
         y = "") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
}


pdf('CC_lollipop_sig_score_group.pdf')
lollipop_plot <- plot_lollipop_summary(CC_data, "tumor_marker_status", c("B cells memory"))
print(lollipop_plot)
lollipop_plot <- plot_lollipop_summary(CC_data, "Stage", c("Neutrophils"))
print(lollipop_plot)
dev.off() 
pdf('EC_lollipop_sig_score_group.pdf')
lollipop_plot <- plot_lollipop_summary(EC_data, "differentiation", c("T cells CD4 memory activated"))
 print(lollipop_plot)
lollipop_plot <- plot_lollipop_summary(EC_data, "Stage", c("T cells gamma delta"))
print(lollipop_plot)
dev.off() 
pdf('OC_lollipop_sig_score_group.pdf')
lollipop_plot <- plot_lollipop_summary(OC_data, "differentiation", c("Monocytes"))
print(lollipop_plot)
lollipop_plot <- plot_lollipop_summary(OC_data, "BRCA_status", c("Plasma cells"))
print(lollipop_plot)
lollipop_plot <- plot_lollipop_summary(OC_data, "Stage", c("B cells memory","T cells CD8","Dendritic cells resting"))
print(lollipop_plot)
dev.off() 

plot_lollipop_summary <- function(data, clinical_var, score_vars) {
  
  # 
  summary_data <- data %>%
    select(all_of(c(clinical_var, score_vars))) %>%
    group_by(.data[[clinical_var]]) %>%
    summarise(across(all_of(score_vars), 
                     list(median = ~median(., na.rm = TRUE),
                          q1 = ~quantile(., 0.25, na.rm = TRUE),
                          q3 = ~quantile(., 0.75, na.rm = TRUE)),
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(cols = starts_with(score_vars),
                 names_to = c("Score", ".value"),
                 names_pattern = "(.*)_(median|q1|q3)")
  summary_data[[clinical_var]] <- factor(
    summary_data[[clinical_var]],
    levels = c("low", "mid", "high","NA")   # 
  )
  
  # 
  ggplot(summary_data, aes(x = .data[[clinical_var]], 
                           y = median, 
                           color = .data[[clinical_var]])) +
    geom_linerange(aes(ymin = q1, ymax = q3), size = 1.5, alpha = 0.5) +
    geom_point(size = 4) +
    facet_wrap(~Score, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = paste0(clinical_var, " "),
         x = clinical_var,
         y = "") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
}
pdf('EC_differentiation_lollipop_sig_score_group.pdf')
lollipop_plot <- plot_lollipop_summary(EC_data, "differentiation", c("T cells CD4 memory activated"))
print(lollipop_plot)
dev.off()
pdf('OC_differentiation_lollipop_sig_score_group.pdf')
lollipop_plot <- plot_lollipop_summary(OC_data, "differentiation", c("Monocytes"))
print(lollipop_plot)
dev.off()
