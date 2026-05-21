#!/usr/bin/env Rscript
# ============================================================
# 00_setup.R — 环境配置
# 运行前修改 BASE_DIR 为你的本地数据路径
# ============================================================
# 使用方式：
#   source("scripts/00_setup.R")
#   然后在脚本中使用 DATA_DIR, MODEL_DIR, OUTPUT_DIR 等变量
# ============================================================

# 仓库根目录（自动检测）
REPO_DIR <- normalizePath(file.path(getwd(), ".."))

# ============================================================
# ⚠️ 修改这里：指向你的数据目录
# 栅格数据较大（>10GB），未纳入版本控制
# 从以下来源下载原始数据：
#   - WorldClim: https://worldclim.org
#   - GBIF: https://gbif.org
#   - FAO HWSD: https://fao.org/soils-portal
# ============================================================
DATA_DIR <- file.path(dirname(REPO_DIR), "毕设文件存放", "变量筛选")

# 子目录
dirs <- list(
  raw_data    = file.path(DATA_DIR, "原始数据"),
  climate     = file.path(DATA_DIR, "气候变化预测", "climate"),
  soil        = file.path(DATA_DIR, "MOP分析_cleaned"),
  env_10var   = file.path(DATA_DIR, "完整数据集_10变量"),
  models      = file.path(DATA_DIR, "brms贝叶斯分析"),
  unified_2km = file.path(DATA_DIR, "未来情景_2km_unified"),
  mop_output  = file.path(DATA_DIR, "MOP分析_final"),
  china_shp   = file.path(REPO_DIR, "data", "shapefiles")
)

# ============================================================
# 依赖检查
# ============================================================
required_pkgs <- c("terra", "sf", "maxnet", "mgcv", "pROC", 
                   "ggplot2", "dplyr", "corrplot", "tidyr", "stringr")

check_packages <- function() {
  missing <- character()
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  if (length(missing) > 0) {
    cat("⚠️ 缺少以下R包，请运行 install.packages(c(",
        paste(shQuote(missing), collapse = ", "), "))\n")
  } else {
    cat("✅ 所有R包已安装\n")
  }
}

# ============================================================
# 输出目录
# ============================================================
OUTPUT_DIR <- file.path(REPO_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("仓库目录: %s\n", REPO_DIR))
cat(sprintf("数据目录: %s\n", DATA_DIR))
cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
