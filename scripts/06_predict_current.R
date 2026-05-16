#!/usr/bin/env Rscript
# 修复版：当前分布预测（确保栅格保存正确）
# 基于用户提供的修复方案

cat("=== 修复版当前分布预测（移除土壤类型11128） ===\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("模型: 移除准完全分离的土壤类型11128后的稳健贝叶斯模型\n")
cat("策略: 分批处理，确保栅格保存正确\n")
cat("========================================\n\n")

# 加载必要的包
cat("1. 加载R包...\n")
if (!requireNamespace("terra", quietly = TRUE)) {
  install.packages("terra", repos = "https://cloud.r-project.org/")
}
if (!requireNamespace("brms", quietly = TRUE)) {
  install.packages("brms", repos = "https://cloud.r-project.org/")
}
library(terra)
library(brms)

# 设置路径
cat("2. 设置文件路径...\n")
params_path <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/brms贝叶斯分析/scaling_parameters_without_11128.rds"
model_path <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/brms贝叶斯分析/final_model_without_11128.rds"
bio1_path <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/中国生境适宜度图_1km_fast/bio1_1km.tif"
bio19_path <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/中国生境适宜度图_1km_fast/bio19_1km.tif"
soil_path <- "/Users/hjhj/.openclaw/workspace/数据汇总/预处理后数据/土壤数据/soil_china_1km.tif"
output_dir <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/brms贝叶斯分析/current_prediction_fixed"

# 检查文件是否存在
for (path in c(params_path, model_path, bio1_path, bio19_path, soil_path)) {
  if (!file.exists(path)) {
    stop("文件不存在: ", path)
  }
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("   输出目录:", output_dir, "\n")

# 加载参数和模型
cat("\n3. 加载参数和模型...\n")
scaling <- readRDS(params_path)
model <- readRDS(model_path)

cat("   ✅ 加载成功\n")
cat("   bio1_mean =", scaling$bio1_mean, "\n")
cat("   bio1_sd =", scaling$bio1_sd, "\n")
cat("   bio19_mean =", scaling$bio19_mean, "\n")
cat("   bio19_sd =", scaling$bio19_sd, "\n")
cat("   土壤类型水平:", paste(scaling$soil_levels, collapse=", "), "\n")

# 读取栅格
cat("\n4. 读取栅格数据...\n")
bio1 <- rast(bio1_path)
bio19 <- rast(bio19_path)
soil <- rast(soil_path)

cat("   原始栅格信息:\n")
cat("     bio1: 范围", paste(round(ext(bio1), 2), collapse=", "), 
    "，像元数", ncell(bio1), "\n")
cat("     bio19: 范围", paste(round(ext(bio19), 2), collapse=", "), 
    "，像元数", ncell(bio19), "\n")
cat("     soil: 范围", paste(round(ext(soil), 2), collapse=", "), 
    "，像元数", ncell(soil), "\n")

# 确保对齐
cat("\n5. 对齐栅格...\n")
ref <- bio1
bio19_aligned <- resample(bio19, ref, method = "bilinear")
soil_aligned <- resample(soil, ref, method = "near")

cat("   ✅ 栅格对齐完成\n")
cat("     对齐后范围一致:", compareGeom(ref, bio19_aligned, soil_aligned, stopOnError = FALSE), "\n")

# 创建结果栅格
result <- ref
values(result) <- NA

# 分批处理
cat("\n6. 设置分批处理参数...\n")
total_cells <- ncell(ref)
batch_size <- 50000  # 适当增大批次
n_batches <- ceiling(total_cells / batch_size)

cat("   总像元数:", total_cells, "\n")
cat("   批次大小:", batch_size, "\n")
cat("   批次数:", n_batches, "\n")

# 创建进度文件
progress_file <- file.path(output_dir, "progress_fixed.txt")
writeLines(paste("开始时间:", Sys.time()), progress_file)

# 主循环
cat("\n7. 开始分批处理...\n")
cat("   ========================================\n")

total_start_time <- Sys.time()

for (i in 1:n_batches) {
  batch_start_time <- Sys.time()
  
  start <- (i-1) * batch_size + 1
  end <- min(i * batch_size, total_cells)
  
  cat("\n   --- 批次", i, "/", n_batches, " ---\n")
  cat("   处理像元:", start, "-", end, " (", end - start + 1, "个)\n")
  
  # 提取当前批次的值
  bio1_vals <- as.vector(values(ref)[start:end])
  bio19_vals <- as.vector(values(bio19_aligned)[start:end])
  soil_vals <- as.vector(values(soil_aligned)[start:end])
  
  # 去除 NA
  valid <- !is.na(bio1_vals) & !is.na(bio19_vals) & !is.na(soil_vals)
  n_valid <- sum(valid)
  
  if (n_valid == 0) {
    cat("   所有像元都是NA，跳过\n")
    
    # 更新进度文件
    progress_line <- sprintf("批次 %04d/%04d: 0个有效像元，跳过", i, n_batches)
    write(progress_line, progress_file, append = TRUE)
    
    next
  }
  
  cat("   有效像元数:", n_valid, "\n")
  
  # 创建数据框
  df <- data.frame(
    bio1 = bio1_vals[valid],
    bio19 = bio19_vals[valid],
    soil = soil_vals[valid]
  )
  
  # 标准化连续变量
  df$bio1_scaled <- (df$bio1 - scaling$bio1_mean) / scaling$bio1_sd
  df$bio19_scaled <- (df$bio19 - scaling$bio19_mean) / scaling$bio19_sd
  
  # 土壤类型映射（与训练数据一致）
  soil_levels <- scaling$soil_levels
  soil_levels_numeric <- as.numeric(soil_levels[soil_levels != "Other"])
  
  # 四舍五入土壤值
  df$soil_int <- round(df$soil)
  
  # 映射到训练数据的类别
  df$soil_type_merged <- ifelse(
    df$soil_int %in% soil_levels_numeric,
    as.character(df$soil_int),
    "Other"
  )
  
  # 特别处理11128
  df$soil_type_merged[df$soil_int == 11128] <- "Other"
  
  # 转换为因子
  df$soil_type_merged <- factor(df$soil_type_merged, levels = soil_levels)
  
  # 预测
  cat("   运行模型预测...\n")
  pred <- fitted(model, newdata = df, scale = "response")[, "Estimate"]
  
  # 验证预测结果
  cat("   预测值统计: 均值=", round(mean(pred), 6), 
      "，范围=[", round(min(pred), 6), ",", round(max(pred), 6), "]\n")
  
  # 赋值回栅格
  idx <- start:end
  result[idx[valid]] <- pred
  
  # 验证赋值
  assigned_vals <- values(result)[idx[valid]]
  n_assigned <- sum(!is.na(assigned_vals))
  
  if (n_assigned == length(pred)) {
    cat("   ✅ 预测值成功赋值到栅格\n")
  } else {
    cat("   ⚠️ 警告: 赋值数量不匹配 (", n_assigned, "/", length(pred), ")\n")
  }
  
  # 每批都保存当前批次结果（用于调试）
  batch_file <- file.path(output_dir, sprintf("batch_%04d_debug.csv", i))
  batch_data <- data.frame(
    x = xyFromCell(ref, idx[valid])[, 1],
    y = xyFromCell(ref, idx[valid])[, 2],
    pred_prob = pred
  )
  write.csv(batch_data, batch_file, row.names = FALSE)
  
  # 清理内存
  rm(df, pred, batch_data)
  gc()
  
  # 每100批保存一次中间栅格
  if (i %% 100 == 0) {
    temp_raster_path <- file.path(output_dir, "temp_result_fixed.tif")
    writeRaster(result, temp_raster_path, overwrite = TRUE)
    cat("   💾 中间栅格已保存:", temp_raster_path, "\n")
  }
  
  # 计算批次用时
  batch_end_time <- Sys.time()
  batch_time <- difftime(batch_end_time, batch_start_time, units = "secs")
  
  # 更新进度文件
  progress_line <- sprintf("批次 %04d/%04d: 处理 %d 像元，用时 %.1f 秒", 
                          i, n_batches, n_valid, batch_time)
  write(progress_line, progress_file, append = TRUE)
  
  cat("   ✅ 批次", i, "完成，用时:", round(batch_time, 2), "秒\n")
  cat("      速度:", round(n_valid / as.numeric(batch_time), 0), "像元/秒\n")
}

# 最终保存
cat("\n8. 保存最终结果...\n")
final_raster_path <- file.path(output_dir, "current_suitability_1km_fixed.tif")
writeRaster(result, final_raster_path, overwrite = TRUE)
cat("   ✅ 最终栅格已保存到:", final_raster_path, "\n")

# 计算统计信息
cat("\n9. 计算统计信息...\n")
vals <- values(result)
vals_valid <- vals[!is.na(vals)]
n_valid_total <- length(vals_valid)

if (n_valid_total > 0) {
  stats <- list(
    n_cells = n_valid_total,
    mean_prob = mean(vals_valid),
    sd_prob = sd(vals_valid),
    min_prob = min(vals_valid),
    max_prob = max(vals_valid),
    median_prob = median(vals_valid),
    q01 = quantile(vals_valid, 0.01),
    q05 = quantile(vals_valid, 0.05),
    q25 = quantile(vals_valid, 0.25),
    q75 = quantile(vals_valid, 0.75),
    q95 = quantile(vals_valid, 0.95),
    q99 = quantile(vals_valid, 0.99)
  )
  
  cat("   预测概率统计:\n")
  cat("     有效像元数:", stats$n_cells, "/", total_cells, 
      "(", round(stats$n_cells/total_cells*100, 1), "%)\n")
  cat("     均值:", format(stats$mean_prob, digits = 6), "\n")
  cat("     中位数:", format(stats$median_prob, digits = 6), "\n")
  cat("     范围: [", format(stats$min_prob, digits = 6), 
      ", ", format(stats$max_prob, digits = 6), "]\n", sep = "")
  cat("     标准差:", format(stats$sd_prob, digits = 6), "\n")
  cat("     重要分位数:\n")
  cat("       1%:", format(stats$q01, digits = 6), "\n")
  cat("       5%:", format(stats$q05, digits = 6), "\n")
  cat("       25%:", format(stats$q25, digits = 6), "\n")
  cat("       75%:", format(stats$q75, digits = 6), "\n")
  cat("       95%:", format(stats$q95, digits = 6), "\n")
  cat("       99%:", format(stats$q99, digits = 6), "\n")
  
  # 保存统计信息
  stats_file <- file.path(output_dir, "prediction_statistics_fixed.txt")
  sink(stats_file)
  cat("=== 修复版预测统计 ===\n")
  cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("模型: 移除土壤类型11128后的稳健贝叶斯模型\n")
  cat("总有效像元数:", stats$n_cells, "/", total_cells, "\n")
  cat("预测概率统计:\n")
  cat("  均值:", stats$mean_prob, "\n")
  cat("  中位数:", stats$median_prob, "\n")
  cat("  最小值:", stats$min_prob, "\n")
  cat("  最大值:", stats$max_prob, "\n")
  cat("  标准差:", stats$sd_prob, "\n")
  cat("  1%分位数:", stats$q01, "\n")
  cat("  5%分位数:", stats$q05, "\n")
  cat("  25%分位数:", stats$q25, "\n")
  cat("  75%分位数:", stats$q75, "\n")
  cat("  95%分位数:", stats$q95, "\n")
  cat("  99%分位数:", stats$q99, "\n")
  sink()
  cat("   ✅ 统计信息保存到:", stats_file, "\n")
} else {
  cat("   ⚠️ 警告: 没有有效预测值\n")
}

# 生成可视化
cat("\n10. 生成可视化...\n")

if (n_valid_total > 0) {
  # 图1: 固定范围 0-0.05 (5%)
  plot_path1 <- file.path(output_dir, "suitability_fixed_range_0_5.png")
  png(plot_path1, width = 1200, height = 800, res = 150)
  plot(result, 
       range = c(0, 0.05),
       col = hcl.colors(100, "YlOrRd", rev = TRUE),
       main = "大林姬鼠当前生境适宜度（修复版，0-5%范围）",
       sub = paste("均值:", format(stats$mean_prob, digits = 4)))
  dev.off()
  cat("   ✅ 固定范围图保存到:", plot_path1, "\n")
  
  # 图2: 使用99%分位数作为上限
  plot_path2 <- file.path(output_dir, "suitability_99percentile.png")
  png(plot_path2, width = 1200, height = 800, res = 150)
  plot(result, 
       range = c(0, stats$q99),
       col = hcl.colors(100, "Viridis"),
       main = "大林姬鼠当前生境适宜度（修复版，99%分位数范围）",
       sub = paste("范围: 0 - ", format(stats$q99, digits = 4)))
  dev.off()
  cat("   ✅ 99%分位数图保存到:", plot_path2, "\n")
  
  # 图3: 分位数着色
  plot_path3 <- file.path(output_dir, "suitability_quantile.png")
  png(plot_path3, width = 1200, height = 800, res = 150)
  
  # 计算分位数断点
  n_classes <- 100
  quants <- quantile(vals_valid, probs = seq(0, 1, length.out = n_classes + 1), na.rm = TRUE)
  
  plot(result, 
       breaks = quants,
       col = hcl.colors(n_classes, "Plasma"),
       main = "大林姬鼠当前生境适宜度（修复版，分位数着色）",
       sub = "突出相对差异")
  dev.off()
  cat("   ✅ 分位数图保存到:", plot_path3, "\n")
} else {
  cat("   ⚠️ 没有有效数据，跳过可视化\n")
}

# 检查模型校准
cat("\n11. 检查模型校准...\n")

training_data_path <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/brms贝叶斯分析/data_for_brms_without_11128.csv"
if (file.exists(training_data_path)) {
  training_data <- read.csv(training_data_path)
  
  if ("response" %in% colnames(training_data)) {
    presence_rate <- mean(training_data$response, na.rm = TRUE)
    n_total_train <- nrow(training_data)
    n_presence <- sum(training_data$response == 1, na.rm = TRUE)
    n_absence <- sum(training_data$response == 0, na.rm = TRUE)
    
    cat("   训练数据统计:\n")
    cat("     总样本数:", n_total_train, "\n")
    cat("     presence数:", n_presence, "\n")
    cat("     absence数:", n_absence, "\n")
    cat("     presence比例:", format(presence_rate, digits = 4), 
        "(", round(presence_rate*100, 2), "%)\n")
    
    if (n_valid_total > 0) {
      pred_mean <- stats$mean_prob
      calibration_ratio <- pred_mean / presence_rate
      
      cat("   预测统计:\n")
      cat("     预测概率均值:", format(pred_mean, digits = 4), 
          "(", round(pred_mean*100, 2), "%)\n")
      cat("     校准比例（预测/训练）:", format(calibration_ratio, digits = 3), "\n")
      
      if (abs(calibration_ratio - 1) < 0.3) {
        cat("     ✅ 校准合理（预测均值与训练presence比例接近）\n")
      } else if (calibration_ratio < 0.5) {
        cat("     ⚠️ 可能低估：预测均值远低于训练presence比例\n")
        cat("       可能原因: 背景点比例过高，模型整体预测概率偏低\n")
      } else if (calibration_ratio > 2) {
        cat("     ⚠️ 可能高估：预测均值远高于训练presence比例\n")
      } else {
        cat("     ⚠️ 校准有偏差，但在可接受范围内\n")
      }
      
      # 保存校准信息
      cal_file <- file.path(output_dir, "calibration_info.txt")
      sink(cal_file)
      cat("=== 模型校准信息 ===\n")
      cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      cat("训练数据:", training_data_path, "\n")
      cat("总样本数:", n_total_train, "\n")
      cat("presence数:", n_presence, "\n")
      cat("absence数:", n_absence, "\n")
      cat("训练presence比例:", presence_rate, "\n")
      cat("预测概率均值:", pred_mean, "\n")
      cat("校准比例（预测/训练）:", calibration_ratio, "\n")
      cat("背景点比例:", n_absence/n_total_train, "\n")
      sink()
      cat("   ✅ 校准信息保存到:", cal_file, "\n")
    }
  }
}

# 计算总用时
total_end_time <- Sys.time()
total_time <- difftime(total_end_time, total_start_time, units = "mins")

cat("\n12. 完成总结:\n")
cat("   ✅ 数据处理: 完成\n")
cat("   ✅ 模型预测: 完成\n")
cat("   ✅ 结果保存: 栅格文件 + 统计数据 + 可视化\n")
cat("   ✅ 总用时:", round(total_time, 2), "分钟\n")

cat("\n输出文件清单:\n")
cat("   预测栅格:", final_raster_path, "\n")
cat("   进度日志:", progress_file, "\n")
if (n_valid_total > 0) {
  cat("   统计信息:", stats_file, "\n")
  cat("   可视化图:", plot_path1, "\n")
  cat("              ", plot_path2, "\n")
  cat("              ", plot_path3, "\n")
}
cat("   输出目录:", output_dir, "\n")

cat("\n=== 修复版当前分布预测完成 ===\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")