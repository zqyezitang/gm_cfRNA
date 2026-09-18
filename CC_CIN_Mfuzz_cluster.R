##################Analysis code of cfRNA
library(readr)
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Mfuzz)
library(tidyverse)
library(clinfun)  
library(limma)    

mat = expr_data[,1:102]
groups = as.factor(meta$treat_time[1:102])
#groups = as.factor(meta$cancer[1:102])
group_score = as.numeric(groups)  # 1,2,3,4


# 1. Non-parametric monotonic trend per gene
jt_pval <- apply(mat, 1, function(y) {
  jonckheere.test(y, group_score, alternative = "two.sided")$p.value
})

# 2. FDR correction
jt_qval <- p.adjust(jt_pval, method = "fdr")

# 3. Effect size (Kendall's tau) for direction
tau <- apply(mat, 1, function(y) cor(y, group_score, method = "kendall"))
jt <- cbind(jt_pval,tau)

# 4. Select significant monotonic genes
sig_genes <- rownames(mat)[jt_pval < 0.05 & abs(tau) > 0.2]
write.csv(jt, "jt_results.csv", row.names = TRUE)
sig_jt <- jt[sig_genes,,drop=FALSE]

write.csv(sig_jt, "jt_sig_results.csv", row.names = TRUE)


# 5. Now run Mfuzz ONLY on sig_genes for shape visualization
expr_mfuzz <- expr_data[sig_genes, , drop = FALSE]

# ========== Mfuzz ==========
#mfuzz for total genes
control_means <- rowMeans(expr_data[, 1:27], na.rm = TRUE)
LSIL_means <- rowMeans(expr_data[, 28:37], na.rm = TRUE)
HSIL_means <- rowMeans(expr_data[, 38:62], na.rm = TRUE)
CIN_means <- rowMeans(expr_data[, 28:62], na.rm = TRUE)
CC_means <- rowMeans(expr_data[, 63:102], na.rm = TRUE)

#mfuzz for sig genes
control_means <- rowMeans(expr_mfuzz[, 1:27], na.rm = TRUE)
LSIL_means <- rowMeans(expr_mfuzz[, 28:37], na.rm = TRUE)
HSIL_means <- rowMeans(expr_mfuzz[, 38:62], na.rm = TRUE)
CIN_means <- rowMeans(expr_mfuzz[, 28:62], na.rm = TRUE)
CC_means <- rowMeans(expr_mfuzz[, 63:102], na.rm = TRUE)


library(Mfuzz)
library(tidyverse)

exp_data <- as.matrix(cbind(control_means,LSIL_means,HSIL_means,CC_means))
time_points <- c("Control","LSIL","HSIL","CC")

exp_data <- as.matrix(cbind(control_means,CIN_means,CC_means))
time_points <- c("Control","CIN","Cervical")
exp_data <- log2(exp_data+1)

#
eset <- new("ExpressionSet", exprs = as.matrix(exp_data))
eset <- filter.NA(eset, thres = 0.25)      
eset <- fill.NA(eset, mode = "knn")        
eset <- filter.std(eset, min.std = 0)      
eset <- standardise(eset)                  

m <- mestimate(eset)

c_range <- 4:12
pc_values <- numeric(length(c_range))
withinss_values <- numeric(length(c_range))
for(i in seq_along(c_range)) {
  set.seed(123)
  cl <- mfuzz(eset, c = c_range[i], m = m)
  centers <- cl$centers
  membership <- cl$membership
  expression <- exprs(eset)
  total_withinss <- 0
  for(j in 1:c_range[i]) {
    cluster_genes <- membership[, j] > 0.5
    if(sum(cluster_genes) > 0) {
      distances <- apply(expression[cluster_genes, ], 1, 
                         function(x) sum((x - centers[j, ])^2))
      total_withinss <- total_withinss + sum(distances * membership[cluster_genes, j])
    }
  }
  withinss_values[i] <- total_withinss
  cat("c =", c_range[i], ", WithinSS =", round(withinss_values[i], 2), "\n")
}
plot(c_range, withinss_values, type = "b", 
     xlab = "Number of clusters (c)", 
     ylab = "Within-cluster Sum of Squares",
     main = "Elbow Method for Optimal c",
     pch = 19, col = "darkgreen", lwd = 2)
grid()

diff_withinss <- diff(withinss_values)
elbow_point <- c_range[which.min(diff_withinss) + 1]
abline(v = elbow_point, col = "red", lty = 2, lwd = 2)
text(elbow_point, withinss_values[elbow_point - min(c_range) + 1], 
     paste("Elbow at c =", elbow_point), pos = 4, col = "red")
cat("\nc:", elbow_point, "\n")


set.seed(123)
c=5
cl <- mfuzz(eset, c = 5, m = m)

pdf("mfuzz_clusters.pdf", width = 8, height = 5)
mfuzz.plot2(eset, cl = cl, mfrow = c(2, 3),
            time.labels = time_points,
            min.mem = 0.5,  
            x11 = FALSE)
dev.off()

gene_cluster <- data.frame(
  gene = rownames(exprs(eset)),
  cluster = cl$cluster,
  membership_max = apply(cl$membership, 1, max)
)
write.csv(gene_cluster, "mfuzz_results.csv", row.names = FALSE)
