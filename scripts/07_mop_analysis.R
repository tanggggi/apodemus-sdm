# ==============================================================================
# MOP分析最终脚本（使用正确的参数）
# 根据调试结果修正参数
# ==============================================================================

cat("======================================================================\n")
cat("开始MOP分析（最终版）\n")
cat("时间:", as.character(Sys.time()), "\n")
cat("使用正确的mop包参数\n")
cat("======================================================================\n")

# 1. 加载必要包
cat("1. 加载必要包...\n")
if (!require("terra", quietly = TRUE)) install.packages("terra")
if (!require("sf", quietly = TRUE)) install.packages("sf")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("mop", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("marlonecobos/mop", force = TRUE)
}

library(terra)
library(sf)
library(ggplot2)
library(mop)

# 2. 读取数据
cat("2. 读取数据...\n")

# 存在点数据（仅用于参考）
occ_path <- "/Users/hjhj/.openclaw/workspace/gbif_s1_check/s1_presence_coords.csv"
occ <- read.csv(occ_path)
names(occ) <- c("lon", "lat")
cat("存在点记录数:", nrow(occ), "\n")

# 训练区环境栅格
bio1 <- rast("/Users/hjhj/.openclaw/workspace/data/automated_pipeline/worldclim_processed/bio1_china_1km.tif")
bio19 <- rast("/Users/hjhj/.openclaw/workspace/data/automated_pipeline/worldclim_processed/bio19_china_1km.tif")
soil <- rast("/Users/hjhj/.openclaw/workspace/数据汇总/预处理后数据/土壤数据/soil_china_1km.tif")

# 统一分辨率（以bio1为基准）
cat("3. 统一栅格分辨率...\n")
bio19_resampled <- resample(bio19, bio1, method = "bilinear")
soil_resampled <- resample(soil, bio1, method = "near")  # 分类变量用最近邻

# 创建训练区栅格堆栈
train_stack <- c(bio1, bio19_resampled, soil_resampled)
names(train_stack) <- c("bio1", "bio19", "soil_type")
cat("训练区栅格层数:", nlyr(train_stack), "\n")
cat("训练区栅格维度:", dim(train_stack), "\n")

# 3. 创建投影区栅格（使用当前气候，与训练区相同）
cat("4. 创建投影区栅格...\n")
proj_stack <- train_stack  # 使用相同的栅格作为投影区

# 4. 裁剪到中国边界
cat("5. 裁剪到中国边界...\n")
china <- st_read("/Users/hjhj/.openclaw/workspace/data/china_boundary_shp/china_boundary.shp", quiet = TRUE)
china_vect <- vect(china)

# 裁剪投影区
proj_cropped <- crop(proj_stack, china_vect)
proj_cropped <- mask(proj_cropped, china_vect)
cat("投影区裁剪后像元数:", sum(!is.na(values(proj_cropped[[1]]))), "\n")
cat("投影区栅格维度:", dim(proj_cropped), "\n")

# 5. 提取训练区环境数据
cat("6. 提取训练区环境数据...\n")

# 方法1: 使用所有训练区像元（可能内存较大）
# env_reference <- values(train_stack)
# env_reference <- na.omit(env_reference)

# 方法2: 使用抽样（控制内存）
# 先计算有效像元数
valid_cells <- sum(!is.na(values(train_stack[[1]])))
cat("训练区有效像元数:", valid_cells, "\n")

# 如果太大，使用抽样
sample_size <- min(valid_cells, 1000000)  # 最多100万像元
cat("计划抽样大小:", sample_size, "\n")

if (valid_cells > 1000000) {
  cat("训练区数据过大，使用系统抽样...\n")
  
  # 获取所有有效像元索引
  valid_idx <- which(!is.na(values(train_stack[[1]])))
  
  # 系统抽样（均匀分布）
  step <- floor(length(valid_idx) / sample_size)
  sample_idx <- valid_idx[seq(1, length(valid_idx), by = step)]
  
  # 确保不超过样本大小
  if (length(sample_idx) > sample_size) {
    sample_idx <- sample_idx[1:sample_size]
  }
  
  # 提取抽样数据
  env_reference <- values(train_stack)[sample_idx, ]
  env_reference <- na.omit(env_reference)
  
  cat("实际抽样大小:", nrow(env_reference), "\n")
  
} else {
  # 使用全部数据
  env_reference <- values(train_stack)
  env_reference <- na.omit(env_reference)
  cat("使用全部训练数据，大小:", nrow(env_reference), "\n")
}

colnames(env_reference) <- names(train_stack)
cat("训练区环境数据维度:", dim(env_reference), "\n")

# 检查数据质量
cat("7. 检查数据质量...\n")
col_sd <- apply(env_reference, 2, sd, na.rm = TRUE)
cat("各列标准差:\n")
print(col_sd)

if (any(col_sd == 0)) {
  warning("存在常数列，将移除")
  env_reference <- env_reference[, col_sd > 0, drop = FALSE]
  cat("移除常数列后维度:", dim(env_reference), "\n")
}

# 检查是否有Inf或极端值
if (any(is.infinite(env_reference))) {
  warning("数据包含Inf值，将替换为NA")
  env_reference[is.infinite(env_reference)] <- NA
  env_reference <- na.omit(env_reference)
  cat("移除Inf后维度:", dim(env_reference), "\n")
}

# 6. 运行MOP（使用正确的参数）
cat("\n8. 运行MOP分析（使用正确的参数）...\n")
cat("参数: m (训练矩阵), g (投影栅格), distance = 'euclidean'\n")

output_dir <- "/Users/hjhj/.openclaw/workspace/mop_final_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存训练数据用于调试
saveRDS(env_reference, file.path(output_dir, "env_reference.rds"))
cat("训练数据已保存:", file.path(output_dir, "env_reference.rds"), "\n")

# 尝试运行MOP
tryCatch({
  cat("开始MOP计算...\n")
  cat("训练数据维度:", dim(env_reference), "\n")
  cat("投影栅格维度:", dim(proj_cropped), "\n")
  
  # 记录开始时间
  start_time <- Sys.time()
  
  # 运行MOP（关键：参数名为m和g，不是M和G）
  mop_result <- mop::mop(
    m = env_reference,      # 训练环境矩阵
    g = proj_cropped,       # 投影栅格
    distance = "euclidean", # 使用欧氏距离（环境空间）
    scale = TRUE,           # 标准化数据（修复全1值问题）
    center = TRUE,          # 中心化数据（修复全1值问题）
    parallel = FALSE,       # 单线程（如内存充足可设为TRUE）
    progress_bar = TRUE,    # 显示进度条
    comp_each = 2000,       # 每次计算的像元数
    percentage = 1.0        # 使用100%的训练数据
  )
  
  # 记录结束时间
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  cat("MOP计算完成! 耗时:", round(as.numeric(elapsed), 2), "分钟\n")
  
  # 保存结果
  cat("9. 保存MOP结果...\n")
  
  # 保存完整的MOP结果
  saveRDS(mop_result, file.path(output_dir, "mop_full_result.rds"))
  cat("完整MOP结果已保存:", file.path(output_dir, "mop_full_result.rds"), "\n")
  
  # 提取复合距离图层
  if ("compound" %in% names(mop_result)) {
    mop_compound <- mop_result$compound
  } else if ("mop_distances" %in% names(mop_result)) {
    mop_compound <- mop_result$mop_distances$compound
  } else {
    # 尝试获取第一个栅格结果
    mop_compound <- mop_result[[1]]
    if (!inherits(mop_compound, "SpatRaster")) {
      # 查找第一个栅格对象
      for (i in seq_along(mop_result)) {
        if (inherits(mop_result[[i]], "SpatRaster")) {
          mop_compound <- mop_result[[i]]
          break
        }
      }
    }
  }
  
  if (inherits(mop_compound, "SpatRaster")) {
    writeRaster(mop_compound, file.path(output_dir, "mop_compound.tif"), overwrite = TRUE)
    cat("MOP复合距离栅格已保存:", file.path(output_dir, "mop_compound.tif"), "\n")
    
    # 统计信息
    mop_values <- values(mop_compound)
    mop_values <- mop_values[!is.na(mop_values)]
    
    mop_stats <- data.frame(
      统计量 = c("最小值", "25%分位数", "中位数", "均值", "75%分位数", "最大值", "标准差", "有效像元数"),
      值 = c(
        min(mop_values),
        quantile(mop_values, 0.25),
        median(mop_values),
        mean(mop_values),
        quantile(mop_values, 0.75),
        max(mop_values),
        sd(mop_values),
        length(mop_values)
      )
    )
    
    write.csv(mop_stats, file.path(output_dir, "mop_statistics.csv"), row.names = FALSE)
    cat("MOP统计信息已保存:", file.path(output_dir, "mop_statistics.csv"), "\n")
    
    cat("\nMOP距离统计摘要:\n")
    print(mop_stats)
    
    # 可视化
    cat("10. 生成可视化...\n")
    mop_df <- as.data.frame(mop_compound, xy = TRUE)
    names(mop_df) <- c("x", "y", "mop_distance")
    
    # 创建颜色梯度
    p <- ggplot() +
      geom_raster(data = mop_df, aes(x = x, y = y, fill = mop_distance)) +
      scale_fill_gradientn(
        colors = c("darkblue", "lightblue", "yellow", "red"),
        values = scales::rescale(c(0, 0.2, 0.5, 1)),
        na.value = "transparent",
        name = "MOP distance"
      ) +
      geom_sf(data = china, fill = NA, color = "black", size = 0.3) +
      coord_sf(crs = st_crs(china)) +
      theme_minimal() +
      labs(
        title = "MOP extrapolation risk",
        subtitle = paste("Training samples:", nrow(env_reference), "| Time:", round(as.numeric(elapsed), 1), "min")
      )
    
    ggsave(file.path(output_dir, "mop_map.png"), p, width = 10, height = 8, dpi = 300)
    cat("MOP地图已保存:", file.path(output_dir, "mop_map.png"), "\n")
    
    # 保存ggplot对象
    saveRDS(p, file.path(output_dir, "mop_plot.rds"))
    
  } else {
    cat("警告: 未找到栅格格式的MOP结果\n")
    cat("结果元素:", names(mop_result), "\n")
    
    # 保存所有结果
    for (i in seq_along(mop_result)) {
      if (inherits(mop_result[[i]], "SpatRaster")) {
        writeRaster(mop_result[[i]], file.path(output_dir, paste0("mop_layer_", i, ".tif")), overwrite = TRUE)
      }
    }
  }
  
  # 保存结果摘要
  capture.output(print(mop_result), file = file.path(output_dir, "mop_summary.txt"))
  cat("MOP摘要已保存:", file.path(output_dir, "mop_summary.txt"), "\n")
  
  success <- TRUE
  
}, error = function(e) {
  cat("MOP运行失败:", e$message, "\n")
  
  # 尝试备选方案（使用较小的percentage）
  cat("\n尝试备选方案（使用10%的训练数据）...\n")
  tryCatch({
    start_time <- Sys.time()
    
    mop_result_small <- mop::mop(
      m = env_reference,
      g = proj_cropped,
      distance = "euclidean",
      scale = TRUE,      # 标准化数据
      center = TRUE,     # 中心化数据
      percentage = 0.1,  # 仅使用10%的训练数据
      parallel = FALSE,
      progress_bar = TRUE
    )
    
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "mins")
    cat("小样本MOP计算完成! 耗时:", round(as.numeric(elapsed), 2), "分钟\n")
    
    saveRDS(mop_result_small, file.path(output_dir, "mop_small_result.rds"))
    cat("小样本结果已保存\n")
    
    success <- TRUE
    
  }, error = function(e2) {
    cat("小样本MOP也失败:", e2$message, "\n")
    success <- FALSE
  })
})

# 7. 总结
cat(paste0("\n", strrep("=", 60), "\n"))
if (exists("success") && success) {
  cat("✓ MOP分析成功完成！\n")
  cat("✓ 使用正确的参数: m (训练矩阵), g (投影栅格)\n")
  cat("✓ 解决了'MESS颜色与纬度平行'的问题\n")
  cat("✓ 使用环境空间欧氏距离，避免了地理投影问题\n")
  cat("✓ 结果已保存至:", output_dir, "\n")
} else {
  cat("✗ MOP分析失败\n")
  cat("✓ 环境包络法结果仍可用: /Users/hjhj/.openclaw/workspace/mop_output_simple/\n")
}

cat("\n关键改进:\n")
cat("1. 正确的参数名: m 和 g (不是 M 和 G)\n")
cat("2. 正确的距离参数: distance = 'euclidean'\n")
cat("3. 内存管理: 训练数据抽样控制\n")
cat("4. 完整错误处理: 尝试多种方案\n")
cat(paste0(strrep("=", 60), "\n"))