#!/usr/bin/env Rscript
# ============================================================
# 01_data_preparation.R — GBIF数据获取与清洗
# ============================================================
# 从GBIF下载大林姬鼠(Apodemus peninsulae)出现记录，
# 清洗无效坐标、重复记录、环境变量提取
# ============================================================

cat("=== 01. 数据准备 ===\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(stringr)
})

# ---- Config ----
source("00_setup.R")

# ---- 1. 从GBIF下载数据 ----
cat("1. 下载GBIF数据...\n")
# 可以使用 rgbif 包自动获取，也可以手动下载CSV
# install.packages("rgbif")
# 
# library(rgbif)
# occ_download_get("KEY", path = dirs$raw_data)
# 
# 手动下载: https://doi.org/10.15468/dl.xxxxxx

gbif_file <- file.path(dirs$raw_data, "occurrence.txt")
if (!file.exists(gbif_file)) {
  cat("  ⚠️ 未找到GBIF文件，请手动下载后放在:", gbif_file, "\n")
  cat("  下载链接: https://www.gbif.org/species/2437763\n")
  cat("  筛选条件: 中国(CN), 坐标存在, 1960-2024年\n")
}

# ---- 2. 清洗 ----
cat("2. 清洗数据...\n")
clean_occurrences <- function(df) {
  df %>%
    # 去除无效坐标
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
    filter(decimalLongitude >= 70, decimalLongitude <= 140) %>%
    filter(decimalLatitude >= 15, decimalLatitude <= 55) %>%
    # 去除坐标精度问题
    filter(coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters)) %>%
    # 去除重复
    distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
    # 去除疑似坐标倒置 / 海面点
    filter(!is.na(countryCode), countryCode == "CN")
}

# ---- 3. 空间稀疏化 ----
cat("3. 空间稀疏化（~5km）...\n")
spatial_thin <- function(df, min_dist_km = 5) {
  # 使用 spThin 包或手动实现
  # install.packages("spThin")
  df  # 占位，实际使用时解注释
}

# ---- 4. 环境变量提取 ----
cat("4. 提取环境变量值...\n")
extract_env <- function(df, env_stack) {
  coords <- df[, c("decimalLongitude", "decimalLatitude")]
  colnames(coords) <- c("longitude", "latitude")
  values <- terra::extract(env_stack, coords, ID = FALSE)
  cbind(df, values) %>% filter(complete.cases(.))
}

cat("\n✅ 01_data_preparation.R 完成\n")
cat("建议工作流:\n")
cat("  1. 从GBIF下载原始CSV → data/raw/\n")
cat("  2. 运行清洗函数 → data/occurrences_10var.csv\n")
cat("  3. 提取环境变量 → 为02_variable_selection.R做准备\n")
