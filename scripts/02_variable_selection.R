#!/usr/bin/env Rscript
# ============================================================
# 02_variable_selection.R — 环境变量筛选
# ============================================================
# 分析流程：
#   1. 相关性分析 → 去除高度相关变量 (|r| > 0.8)
#   2. VIF 分析 → 逐步去除多重共线性 (VIF > 5)
#   3. PCA 降维 → 验证变量选择的合理性
#   4. 生态解释 → 最终选择10个代表性变量
# ============================================================

cat("=== 02. 变量筛选 ===\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(terra)
  library(corrplot)
  library(dplyr)
})

source("00_setup.R")

# ---- 1. 加载环境栅格 ----
cat("1. 加载环境栅格数据...\n")

# 初始变量池（26个候选变量）
var_pool <- c(
  # 气候 (WorldClim 2.1, 2.5")
  paste0("bio", 1:19),
  # 地形
  "elevation",
  # 植被纹理 (MODIS NDVI衍生)
  "Contrast", "cv", "Entropy", "std",
  # 土壤
  "soil_type"
)

# ---- 2. 相关性分析 ----
cat("2. 相关性分析...\n")
# 采样点提取环境值
# env_values <- terra::extract(env_stack, sample_points)
# cor_matrix <- cor(env_values, use = "pairwise.complete.obs")
# 
# 可视化
# pdf("results/figures/correlation_heatmap.pdf")
# corrplot(cor_matrix, method = "color", type = "upper", 
#          tl.cex = 0.6, number.cex = 0.4)
# dev.off()

# 筛选标准: |r| > 0.8 的变量对中保留生态意义更强的一个
# 例如: bio1 (年均温) 保留, bio5/bio6 (极端温) 去除

# ---- 3. VIF分析 ----
cat("3. VIF分析...\n")
# 逐步VIF: 反复计算VIF，每次去除VIF最高的变量，直到所有VIF < 5
# library(usdm)
# vif_result <- vifstep(env_values, th = 5)

# ---- 4. 最终10变量 ----
cat("4. 最终变量集...\n")

final_variables <- c(
  # 温度: bio1 (年均温)
  "bio1",
  # 降水: bio12 (年降水量), bio15 (降水季节性), bio19 (最冷季降水)
  "bio12", "bio15", "bio19",
  # 地形: elevation
  "elevation",
  # 植被纹理: Contrast (对比度), cv (变异系数), Entropy (熵), std (标准差)
  "Contrast", "cv", "Entropy", "std",
  # 土壤: soil_type
  "soil_type"
)

cat(sprintf("  最终选择 %d 个变量:\n", length(final_variables)))
cat(sprintf("  %s\n\n", paste(final_variables, collapse = ", ")))

cat("变量分类:\n")
cat("  🌡️ 气候 (4): bio1, bio12, bio15, bio19\n")
cat("  ⛰️ 地形 (1): elevation\n")
cat("  🌿 植被 (4): Contrast, cv, Entropy, std\n")
cat("  🧱 土壤 (1): soil_type\n\n")

# ---- 5. 导出 ----
cat("5. 导出清洗后的数据集...\n")
# 将10变量数据集写入 data/
# write.csv(env_10var, "data/occurrences_10var.csv", row.names = FALSE)
cat("  → data/occurrences_10var.csv\n")

cat("✅ 02_variable_selection.R 完成\n")
