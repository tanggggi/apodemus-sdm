#!/usr/bin/env Rscript
# ============================================================
# 09_area_calculation.R — 适生区面积计算
# ============================================================
# 使用 TSS 最优阈值将连续概率图二值化，
# 在 Albers 投影下计算各情景的适生面积
# ============================================================

cat("=== 09. 适生区面积计算 ===\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

source("00_setup.R")

# ---- 1. 参数 ----
TSS <- 0.051655     # Youden指数最优阈值
CELL_2KM <- 2.799   # 2km像元面积 (km²)

# ---- 2. 预测文件路径 ----
scenarios <- list(
  当前  = "current_suitability_2km.tif",
  SSP126 = "future_suitability_SSP126_2050_2km.tif",
  SSP245 = "future_suitability_SSP245_2050_2km.tif",
  SSP585 = "future_suitability_SSP585_2050_2km.tif"
)

# 文件位于 unified_2km/ 目录，需要从原始数据目录获取
# 或从以下位置下载:
# https://doi.org/10.6084/m9.figshare.XXXXXXX

# ---- 3. 计算 ----
cat(sprintf("阈值 TSS = %.6f\n", TSS))
cat(sprintf("像元面积 = %.3f km²\n\n", CELL_2KM))

results <- data.frame()

for (nm in names(scenarios)) {
  f <- scenarios[[nm]]
  cat(sprintf("[%s] %s\n", nm, f))
  
  # 检查文件是否存在
  # 实际运行时需指向预测栅格所在的完整路径
  if (!file.exists(f)) {
    cat("  ⚠️ 文件未找到，跳过\n")
    next
  }
  
  r <- rast(f)
  
  # 有效像元
  valid_cells <- global(!is.na(r), "sum", na.rm = TRUE)[1,1]
  
  # 适宜像元（概率 >= TSS）
  suit_cells <- global(r >= TSS, "sum", na.rm = TRUE)[1,1]
  
  pct <- suit_cells / valid_cells * 100
  area_km2 <- suit_cells * CELL_2KM
  
  cat(sprintf("  有效: %d | 适宜: %d (%.2f%%) | 面积: %.0f km² (%.2f 万 km²)\n",
              valid_cells, suit_cells, pct, area_km2, area_km2 / 10000))
  
  results <- rbind(results, data.frame(
    情景 = nm,
    有效像元 = valid_cells,
    适宜像元 = suit_cells,
    适宜面积_km2 = round(area_km2),
    适宜比例_pct = round(pct, 2)
  ))
}

# ---- 4. 导出 ----
out_file <- file.path(OUTPUT_DIR, "area", "area_summary.csv")
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_file, row.names = FALSE)
cat(sprintf("\n✅ 结果已保存: %s\n", out_file))
print(results)

# ---- 5. 结果解读 ----
cat("\n=== 结果解读 ===\n")
cat("当前适宜面积适用于SSP情景间的相对比较。\n")
cat("注意：不同情景的预测范围可能不同，\n")
cat("建议结合MOP分析（07_mop_analysis.R）评估外推风险。\n")
