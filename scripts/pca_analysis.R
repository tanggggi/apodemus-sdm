# PCA分析脚本 - 对最终10个环境变量进行PCA分析
# 作者: 影子 (OpenClaw助手)
# 日期: 2026-04-04

# 1. 加载必要的包
library(tidyverse)
library(FactoMineR) # 用于PCA分析
library(factoextra) # 用于可视化
library(ggplot2)
library(patchwork) # 用于组合图表

# 2. 设置工作目录和输出目录
output_dir <- "/Users/hjhj/Desktop/毕设文件存放/提交版/10.（补）PCA分析"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 3. 准备环境数据
# 读取完整数据集
data_path <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/完整数据集_10变量/完整数据集_10变量.csv"
df_env_raw <- read.csv(data_path)

cat("数据基本信息:\n")
cat("行数:", nrow(df_env_raw), "\n")
cat("列数:", ncol(df_env_raw), "\n")
cat("列名:", names(df_env_raw), "\n")

# 选择10个环境变量进行PCA分析
env_vars <- c("elevation", "soil_type", "Contrast", "cv", "Entropy", 
              "std", "bio15", "bio19", "bio1", "bio12")

df_env <- df_env_raw[, env_vars]

# 检查数据摘要
cat("\n环境变量数据摘要:\n")
print(summary(df_env))

# 检查缺失值
missing_values <- sum(is.na(df_env))
cat("\n缺失值总数:", missing_values, "\n")
if (missing_values > 0) {
  cat("警告: 数据中存在缺失值，PCA分析前将删除含有缺失值的行\n")
  df_env <- na.omit(df_env)
  cat("删除缺失值后行数:", nrow(df_env), "\n")
}

# 4. 执行PCA分析 (自动对数据进行标准化)
cat("\n执行PCA分析...\n")
pca_result <- PCA(df_env, scale.unit = TRUE, graph = FALSE)

# 5. PCA结果摘要
cat("\n=== PCA分析结果摘要 ===\n")

# 特征值
eigenvalues <- pca_result$eig
cat("特征值:\n")
print(eigenvalues[, 1:3])

# 方差解释比例
var_explained <- eigenvalues[, 2]
cum_var_explained <- eigenvalues[, 3]

cat("\n前5个主成分解释的方差比例:\n")
for (i in 1:min(5, nrow(eigenvalues))) {
  cat(sprintf("PC%d: %.2f%% (累积: %.2f%%)\n", 
              i, var_explained[i], cum_var_explained[i]))
}

# 确定保留的主成分数 (特征值>1准则)
n_components <- sum(eigenvalues[, 1] > 1)
cat(sprintf("\n根据特征值>1准则，保留 %d 个主成分\n", n_components))

# 累积方差解释达到80%的主成分数
n_80 <- which(cum_var_explained >= 80)[1]
cat(sprintf("累积方差解释达到80%%需要前 %d 个主成分\n", n_80))

# 6. 绘制顶刊风格图表

## 图1: 特征值与方差解释图（碎石图）
cat("\n生成图1: 特征值与方差解释图...\n")
p_eigen <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50), 
                   barfill = "#2E5A88", barcolor = "#2E5A88",
                   xlab = "主成分", ylab = "方差解释比例 (%)",
                   main = "PCA碎石图: 特征值与方差解释比例") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 保存图1
ggsave(file.path(output_dir, "fig_pca_eigenvalue.png"), 
       p_eigen, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_pca_eigenvalue.pdf"), 
       p_eigen, width = 6, height = 4)

## 图2: 相关性圈图 (Correlation Circle)
cat("生成图2: 相关性圈图...\n")
p_var <- fviz_pca_var(pca_result, 
                     col.var = "contrib", # 根据贡献着色
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # 顶刊常用配色
                     repel = TRUE, # 避免标签重叠
                     title = "环境变量相关性圈图") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(color = "贡献度")

# 保存图2
ggsave(file.path(output_dir, "fig_pca_correlation_circle.png"), 
       p_var, width = 7, height = 6, dpi = 300)
ggsave(file.path(output_dir, "fig_pca_correlation_circle.pdf"), 
       p_var, width = 7, height = 6)

## 图3: 双标图 (Biplot)
# 为了创建有意义的双标图，我们需要对样本进行分组
# 这里使用k-means聚类创建3个分组
cat("生成图3: 双标图...\n")
set.seed(123)
# 使用前2个主成分进行聚类
pca_scores <- pca_result$ind$coord[, 1:2]
kmeans_result <- kmeans(pca_scores, centers = 3, nstart = 25)
df_env_raw$cluster <- as.factor(kmeans_result$cluster)

# 创建双标图
p_biplot <- fviz_pca_biplot(pca_result,
                           col.ind = df_env_raw$cluster, # 按聚类分组着色样本点
                           palette = c("#F8766D", "#00BFC4", "#7CAE00"), # 顶刊常用颜色
                           addEllipses = TRUE, # 添加置信椭圆
                           ellipse.level = 0.68, # 1个标准差椭圆
                           ellipse.type = "confidence",
                           col.var = "black", # 变量向量为黑色
                           repel = TRUE,
                           title = "PCA双标图: 样本分布与变量关系") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(color = "样本聚类", fill = "样本聚类")

# 保存图3
ggsave(file.path(output_dir, "fig_pca_biplot.png"), 
       p_biplot, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "fig_pca_biplot.pdf"), 
       p_biplot, width = 8, height = 6)

# 7. 创建组合图
cat("生成组合图...\n")
p_combined <- (p_eigen + labs(title = "(a) 碎石图")) / 
              (p_var + labs(title = "(b) 相关性圈图")) /
              (p_biplot + labs(title = "(c) 双标图"))

ggsave(file.path(output_dir, "fig_pca_combined.png"), 
       p_combined, width = 10, height = 12, dpi = 300)

# 8. 生成详细的PCA结果表格
cat("生成PCA结果表格...\n")

# 变量贡献度表格
var_contrib <- pca_result$var$contrib
var_cos2 <- pca_result$var$cos2

# 创建贡献度表格
contrib_table <- data.frame(
  变量 = rownames(var_contrib),
  PC1_贡献度 = sprintf("%.2f%%", var_contrib[, 1]),
  PC2_贡献度 = sprintf("%.2f%%", var_contrib[, 2]),
  PC3_贡献度 = sprintf("%.2f%%", var_contrib[, 3]),
  PC1_代表质量 = sprintf("%.3f", var_cos2[, 1]),
  PC2_代表质量 = sprintf("%.3f", var_cos2[, 2]),
  PC3_代表质量 = sprintf("%.3f", var_cos2[, 3])
)

# 保存贡献度表格
write.csv(contrib_table, 
          file.path(output_dir, "pca_variable_contributions.csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")

# 特征值表格
eigen_table <- data.frame(
  主成分 = paste0("PC", 1:nrow(eigenvalues)),
  特征值 = sprintf("%.3f", eigenvalues[, 1]),
  方差解释比例 = sprintf("%.2f%%", eigenvalues[, 2]),
  累积方差解释比例 = sprintf("%.2f%%", eigenvalues[, 3])
)

write.csv(eigen_table, 
          file.path(output_dir, "pca_eigenvalues.csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")

# 9. 生成分析报告
cat("生成分析报告...\n")
report_text <- paste0(
  "# PCA分析报告\n\n",
  "## 数据基本信息\n",
  "- 数据文件: ", data_path, "\n",
  "- 样本数量: ", nrow(df_env_raw), "\n",
  "- 环境变量: ", paste(env_vars, collapse = ", "), "\n",
  "- 分析时间: ", Sys.time(), "\n\n",
  
  "## PCA分析结果\n",
  "### 主成分选择\n",
  "- 特征值>1的主成分数: ", n_components, "\n",
  "- 累积方差解释达到80%所需主成分数: ", n_80, "\n\n",
  
  "### 前5个主成分的方差解释\n",
  paste(sapply(1:min(5, nrow(eigenvalues)), function(i) {
    sprintf("- PC%d: %.2f%% (累积: %.2f%%)", i, var_explained[i], cum_var_explained[i])
  }), collapse = "\n"), "\n\n",
  
  "## 生成图表\n",
  "1. `fig_pca_eigenvalue.png` - 特征值碎石图\n",
  "2. `fig_pca_correlation_circle.png` - 相关性圈图\n",
  "3. `fig_pca_biplot.png` - 双标图\n",
  "4. `fig_pca_combined.png` - 组合图\n\n",
  
  "## 结果文件\n",
  "1. `pca_variable_contributions.csv` - 变量贡献度表格\n",
  "2. `pca_eigenvalues.csv` - 特征值表格\n\n",
  
  "## 生态学解释建议\n",
  "### 相关性圈图解读\n",
  "- **变量向量角度**: 锐角表示正相关，钝角表示负相关\n",
  "- **向量长度**: 越长表示该变量在前两个主成分平面上被解释得越好\n",
  "- **聚类模式**: 聚集在一起的变量可能具有相似的生态效应\n\n",
  
  "### 双标图解读\n",
  "- **样本分布**: 观察样本点在主成分空间中的分布模式\n",
  "- **变量驱动**: 样本点朝向某个变量向量方向，表示该变量对该组样本影响较大\n",
  "- **生态梯度**: 主成分轴可能代表重要的生态梯度（如温度梯度、水分梯度等）"
)

writeLines(report_text, file.path(output_dir, "pca_analysis_report.md"))

# 10. 完成
cat("\n=== PCA分析完成 ===\n")
cat("输出目录:", output_dir, "\n")
cat("生成的文件:\n")
cat("- fig_pca_eigenvalue.png/pdf\n")
cat("- fig_pca_correlation_circle.png/pdf\n")
cat("- fig_pca_biplot.png/pdf\n")
cat("- fig_pca_combined.png\n")
cat("- pca_variable_contributions.csv\n")
cat("- pca_eigenvalues.csv\n")
cat("- pca_analysis_report.md\n")

cat("\n分析完成！\n")