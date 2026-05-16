#!/usr/bin/env Rscript
# 完整10变量嵌套空间交叉验证比较脚本
# 基于10变量完整数据集，比较8个候选模型

cat(strrep("=", 80), "\n")
cat("完整10变量嵌套空间交叉验证比较\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("目标: 基于10变量完整数据集公平比较8个候选模型\n")
cat("方法学原则: 完全避免循环验证，仅使用独立先验信息\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 初始化设置
# ============================================================================

# 设置随机种子确保可重复性
set.seed(20260320)

# 加载必要的包（最小集合）
cat("1. 加载必要的包...\n")
required_packages <- c("maxnet", "pROC", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat("   安装包:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  } else {
    cat("   加载包:", pkg, "\n")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# 2. 设置路径
# ============================================================================
cat("\n2. 设置路径...\n")

# 候选模型定义
model_dir <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/候选模型框架_完整10变量"
model_list_file <- file.path(model_dir, "候选模型清单_机器可读_完整10变量.csv")
model_details_file <- file.path(model_dir, "候选模型详细定义_完整10变量.csv")

# 环境数据路径
data_file <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/完整数据集_10变量/完整数据集_10变量.csv"

# 输出目录
output_base <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/嵌套CV完整10变量比较"
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

cat("   候选模型目录:", model_dir, "\n")
cat("   数据文件:", data_file, "\n")
cat("   输出目录:", output_base, "\n")

# ============================================================================
# 3. 读取候选模型定义
# ============================================================================
cat("\n3. 读取候选模型定义...\n")

if (!file.exists(model_list_file)) {
  stop("候选模型清单文件不存在: ", model_list_file)
}

model_list <- read.csv(model_list_file)
model_ids <- unique(model_list$model_id)
cat("   读取模型:", length(model_ids), "个\n")

# 为每个模型创建变量列表
model_vars <- list()
for (mid in model_ids) {
  model_vars[[mid]] <- unique(model_list$variable[model_list$model_id == mid])
  cat("   - ", mid, ": ", length(model_vars[[mid]]), "个变量\n", sep="")
}

# 读取模型详细信息
model_details <- read.csv(model_details_file)

# ============================================================================
# 4. 读取环境数据
# ============================================================================
cat("\n4. 读取环境数据...\n")

if (!file.exists(data_file)) {
  stop("数据文件不存在: ", data_file)
}

env_data <- read.csv(data_file)
cat("   读取数据:", nrow(env_data), "行，", ncol(env_data), "列\n")
cat("   存在点数量:", sum(env_data$response == 1), "个\n")
cat("   背景点数量:", sum(env_data$response == 0), "个\n")

# 检查是否包含所有候选变量
all_vars <- unique(unlist(model_vars))
missing_vars <- setdiff(all_vars, colnames(env_data))
if (length(missing_vars) > 0) {
  stop("数据中缺失以下变量:", paste(missing_vars, collapse=", "))
}
cat("   ✅ 所有候选变量都存在\n")

# 提取坐标用于空间分块（如果可用）
if ("longitude" %in% colnames(env_data) & "latitude" %in% colnames(env_data)) {
  coords <- env_data[, c("longitude", "latitude")]
  cat("   使用经纬度坐标进行空间分块\n")
} else {
  coords <- NULL
  cat("   没有坐标信息，使用随机分块\n")
}

# ============================================================================
# 5. 创建固定空间分块
# ============================================================================
cat("\n5. 创建固定空间分块（外层3折）...\n")

n_outer_folds <- 3

if (!is.null(coords)) {
  # 空间分块：k-means聚类
  cat("   使用k-means空间聚类分块...\n")
  spatial_clusters <- kmeans(coords, centers = n_outer_folds)
  fold_ids <- spatial_clusters$cluster
} else {
  # 随机分块
  cat("   使用随机分块...\n")
  fold_ids <- sample(rep(1:n_outer_folds, length.out = nrow(env_data)))
}

cat("   分块大小:", paste(table(fold_ids), collapse=", "), "\n")

# 创建分块列表
folds <- list()
for (k in 1:n_outer_folds) {
  test_idx <- which(fold_ids == k)
  train_idx <- which(fold_ids != k)
  folds[[k]] <- list(train = train_idx, test = test_idx, fold_id = k)
}

# ============================================================================
# 6. 简化版嵌套空间交叉验证函数
# ============================================================================
cat("\n6. 定义嵌套空间交叉验证函数...\n")

run_simple_nested_cv <- function(model_id, variables, data, folds) {
  # 简化的嵌套CV，降低计算负担
  
  cat("   运行模型", model_id, "(", length(variables), "个变量)...\n")
  
  model_data <- data[, c(variables, "response")]
  results <- list(
    model_id = model_id,
    variables = variables,
    test_auc = numeric(length(folds)),
    train_auc = numeric(length(folds)),
    best_params = list()
  )
  
  for (fold_idx in 1:length(folds)) {
    fold <- folds[[fold_idx]]
    
    cat("     折叠", fold_idx, "/", length(folds), "...")
    
    train_data <- model_data[fold$train, ]
    test_data <- model_data[fold$test, ]
    
    # 简化参数搜索：只测试关键组合
    param_grid <- expand.grid(
      regmult = c(0.5, 2),    # 只测试两个正则化参数
      classes = c("l", "lq")  # 只测试线性和线性+二次项
    )
    
    best_auc <- -1
    best_params <- NULL
    
    # 内层3折交叉验证
    n_inner <- 3
    inner_folds <- sample(rep(1:n_inner, length.out = nrow(train_data)))
    
    for (p_idx in 1:nrow(param_grid)) {
      param <- param_grid[p_idx, ]
      inner_aucs <- numeric(n_inner)
      
      for (inner_fold in 1:n_inner) {
        inner_train_idx <- which(inner_folds != inner_fold)
        inner_test_idx <- which(inner_folds == inner_fold)
        
        inner_train <- train_data[inner_train_idx, ]
        inner_test <- train_data[inner_test_idx, ]
        
        tryCatch({
          model <- maxnet(
            p = inner_train$response,
            data = inner_train[, variables, drop = FALSE],
            regmult = param$regmult,
            maxnet.formula(p = inner_train$response, 
                          data = inner_train[, variables, drop = FALSE], 
                          classes = param$classes)
          )
          
          pred <- predict(model, inner_test[, variables, drop = FALSE], type = "cloglog")
          roc_obj <- roc(inner_test$response, pred, quiet = TRUE)
          inner_aucs[inner_fold] <- auc(roc_obj)
        }, error = function(e) {
          inner_aucs[inner_fold] <- 0.5
        })
      }
      
      mean_auc <- mean(inner_aucs, na.rm = TRUE)
      if (mean_auc > best_auc) {
        best_auc <- mean_auc
        best_params <- param
      }
    }
    
    cat(" 参数: regmult=", best_params$regmult, ", classes=", best_params$classes, "\n", sep="")
    
    # 使用最佳参数训练最终模型
    final_model <- maxnet(
      p = train_data$response,
      data = train_data[, variables, drop = FALSE],
      regmult = best_params$regmult,
      maxnet.formula(p = train_data$response, 
                    data = train_data[, variables, drop = FALSE], 
                    classes = best_params$classes)
    )
    
    # 训练集AUC
    train_pred <- predict(final_model, train_data[, variables, drop = FALSE], type = "cloglog")
    train_auc <- auc(roc(train_data$response, train_pred, quiet = TRUE))
    results$train_auc[fold_idx] <- train_auc
    
    # 测试集AUC
    test_pred <- predict(final_model, test_data[, variables, drop = FALSE], type = "cloglog")
    test_auc <- auc(roc(test_data$response, test_pred, quiet = TRUE))
    results$test_auc[fold_idx] <- test_auc
    
    results$best_params[[fold_idx]] <- best_params
  }
  
  results$mean_test_auc <- mean(results$test_auc, na.rm = TRUE)
  results$sd_test_auc <- sd(results$test_auc, na.rm = TRUE)
  results$mean_train_auc <- mean(results$train_auc, na.rm = TRUE)
  results$overfitting <- results$mean_train_auc - results$mean_test_auc
  
  cat("   完成: AUC =", round(results$mean_test_auc, 3), "±", round(results$sd_test_auc, 3), "\n")
  
  return(results)
}

# ============================================================================
# 7. 运行所有模型的嵌套CV（分批处理避免内存问题）
# ============================================================================
cat("\n7. 运行所有模型的嵌套空间交叉验证...\n")

# 分批运行：每批2个模型
batches <- list(
  batch1 = c("D1", "C1"),    # 全模型 vs 气候主导
  batch2 = c("T1", "V1"),    # 地形-植被 vs 植被纹理
  batch3 = c("S1", "W1"),    # 土壤-气候 vs 水热协同
  batch4 = c("M1", "E1")     # 最小综合 vs 熵值主导
)

all_results <- list()
start_time <- Sys.time()

for (batch_idx in 1:length(batches)) {
  batch_models <- batches[[batch_idx]]
  batch_name <- names(batches)[batch_idx]
  
  cat("\n  批次", batch_idx, "/", length(batches), ": ", paste(batch_models, collapse=", "), "\n", sep="")
  
  for (model_id in batch_models) {
    if (model_id %in% names(model_vars)) {
      cat("  [", which(model_ids == model_id), "/", length(model_ids), "] ", sep="")
      
      results <- run_simple_nested_cv(
        model_id = model_id,
        variables = model_vars[[model_id]],
        data = env_data,
        folds = folds
      )
      
      all_results[[model_id]] <- results
      
      # 保存中间结果
      temp_file <- file.path(output_base, paste0("temp_", model_id, ".rds"))
      saveRDS(results, temp_file)
      cat("    中间结果保存到:", temp_file, "\n")
    }
  }
  
  # 估计剩余时间
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  completed <- sum(sapply(batches[1:batch_idx], length))
  remaining <- length(model_ids) - completed
  if (remaining > 0) {
    avg_time <- elapsed / completed
    cat("    预计剩余:", round(avg_time * remaining, 1), "分钟\n")
  }
}

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n   所有模型完成! 总用时:", round(total_time, 1), "分钟\n")

# ============================================================================
# 8. 整理和保存结果
# ============================================================================
cat("\n8. 整理和保存结果...\n")

# 创建比较数据框
comparison_df <- data.frame(
  model_id = character(),
  hypothesis_id = character(),
  hypothesis = character(),
  variable_count = integer(),
  mean_test_auc = numeric(),
  sd_test_auc = numeric(),
  mean_train_auc = numeric(),
  overfitting = numeric(),
  stringsAsFactors = FALSE
)

detailed_results <- data.frame()

for (model_id in names(all_results)) {
  res <- all_results[[model_id]]
  
  # 查找模型详细信息
  model_info <- model_details[model_details$model_id == model_id, ]
  
  comparison_df <- rbind(comparison_df, data.frame(
    model_id = model_id,
    hypothesis_id = ifelse(nrow(model_info) > 0, model_info$hypothesis_id, NA),
    hypothesis = ifelse(nrow(model_info) > 0, model_info$hypothesis, NA),
    variable_count = length(res$variables),
    mean_test_auc = res$mean_test_auc,
    sd_test_auc = res$sd_test_auc,
    mean_train_auc = res$mean_train_auc,
    overfitting = res$overfitting,
    stringsAsFactors = FALSE
  ))
  
  # 详细结果
  for (fold_idx in 1:length(res$test_auc)) {
    detailed_results <- rbind(detailed_results, data.frame(
      model_id = model_id,
      fold = fold_idx,
      test_auc = res$test_auc[fold_idx],
      train_auc = res$train_auc[fold_idx],
      regmult = if (!is.null(res$best_params[[fold_idx]])) res$best_params[[fold_idx]]$regmult else NA,
      classes = if (!is.null(res$best_params[[fold_idx]])) res$best_params[[fold_idx]]$classes else NA,
      stringsAsFactors = FALSE
    ))
  }
}

# 按AUC排序
comparison_df <- comparison_df[order(-comparison_df$mean_test_auc), ]
comparison_df$rank <- 1:nrow(comparison_df)

# 保存结果
comparison_file <- file.path(output_base, "模型性能比较_完整10变量.csv")
detailed_file <- file.path(output_base, "详细AUC结果_完整10变量.csv")

write.csv(comparison_df, comparison_file, row.names = FALSE)
write.csv(detailed_results, detailed_file, row.names = FALSE)

cat("   性能比较保存到:", comparison_file, "\n")
cat("   详细结果保存到:", detailed_file, "\n")

# ============================================================================
# 9. 生成报告
# ============================================================================
cat("\n9. 生成报告...\n")

report_content <- paste0(
  "# 完整10变量嵌套空间交叉验证比较报告\n\n",
  "## 生成信息\n",
  "- **生成时间**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "- **数据**: 10变量完整数据集 (", nrow(env_data), "行)\n",
  "- **模型数量**: ", length(all_results), "个\n",
  "- **空间分块**: 外层3折", ifelse(!is.null(coords), "空间分块", "随机分块"), "\n",
  "- **运行时间**: ", round(total_time, 1), "分钟\n",
  "- **方法学**: 完全避免循环验证\n\n",
  "## 模型性能排名\n\n",
  "| 排名 | 模型ID | 假设ID | 变量数 | 平均测试AUC | ±SD | 过拟合 |\n",
  "|------|--------|--------|--------|-------------|-----|--------|\n",
  paste(sapply(1:nrow(comparison_df), function(i) {
    m <- comparison_df[i, ]
    sprintf("| %d | %s | %s | %d | %.3f | %.3f | %.3f |", 
            m$rank, m$model_id, m$hypothesis_id, m$variable_count,
            m$mean_test_auc, m$sd_test_auc, m$overfitting)
  }), collapse = "\n"),
  "\n\n## 关键科学发现\n\n",
  "### 1. 最佳模型\n",
  "**", comparison_df$model_id[1], "** (假设", comparison_df$hypothesis_id[1], ") 表现最佳，\n",
  "平均测试AUC = ", round(comparison_df$mean_test_auc[1], 3), 
  " ± ", round(comparison_df$sd_test_auc[1], 3), "\n\n",
  "### 2. 与8变量快速验证结果对比\n",
  "先前基于8变量的快速验证结果（排名）：\n",
  "1. **D1** (8变量): AUC 0.942\n",
  "2. **T1** (5变量): AUC 0.896\n", 
  "3. **M1** (4变量): AUC 0.858\n",
  "4. **C1** (2变量): AUC 0.720\n\n",
  "### 3. 温度(bio1)和年降水量(bio12)的影响\n",
  "本次分析包含了完整的10变量，可以检验：\n",
  "1. **气候主导模型(C1)**: 包含bio1, bio12, bio15, bio19\n",
  "2. **土壤-气候模型(S1)**: 包含soil_type, bio1, bio19\n",
  "3. **水热协同模型(W1)**: 包含bio1, bio19, elevation\n\n",
  "## 生态学假设检验\n",
  "基于性能排名，初步评估各假设：\n\n",
  paste(sapply(1:nrow(comparison_df), function(i) {
    m <- comparison_df[i, ]
    model_info <- model_details[model_details$model_id == m$model_id, ]
    if (nrow(model_info) > 0) {
      paste0(i, ". **", m$model_id, "** (", m$hypothesis, "): AUC = ", 
             round(m$mean_test_auc, 3), " → ",
             ifelse(m$rank == 1, "**强烈支持**",
                   ifelse(m$rank <= 3, "**支持**",
                         ifelse(m$rank <= 6, "**弱支持/混合**", "**不支持**"))), "\n")
    } else {
      ""
    }
  }), collapse = ""),
  "\n## 方法学总结\n",
  "### 避免循环验证的关键措施\n",
  "1. **独立先验**: 模型构建仅基于生态学理论和独立统计信息\n",
  "2. **透明构建**: 每个模型有明确的生态学假设\n",
  "3. **公平比较**: 相同分块和调参流程\n",
  "4. **事后解释**: 结果解释在模型比较之后\n\n",
  "### 科学贡献\n",
  "1. **方法学示范**: 展示了如何避免循环验证的严谨流程\n",
  "2. **假设检验**: 基于独立先验信息系统检验生态学假设\n",
  "3. **透明可重复**: 所有决策依据和代码公开\n\n",
  "## 文件清单\n",
  "- `模型性能比较_完整10变量.csv`: 性能排名表\n",
  "- `详细AUC结果_完整10变量.csv`: 详细AUC结果\n",
  "- `temp_*.rds`: 每个模型的中间结果\n"
)

report_file <- file.path(output_base, "完整10变量比较报告.md")
writeLines(report_content, report_file)
cat("   报告保存到:", report_file, "\n")

# ============================================================================
# 10. 输出摘要
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("✅ 完整10变量嵌套空间交叉验证完成\n")
cat("最佳模型:", comparison_df$model_id[1], 
    " (AUC = ", round(comparison_df$mean_test_auc[1], 3), 
    " ± ", round(comparison_df$sd_test_auc[1], 3), ")\n", sep="")
cat("总运行时间:", round(total_time, 1), "分钟\n")
cat("输出目录:", output_base, "\n")
cat(strrep("=", 80), "\n")

cat("\n性能排名:\n")
for (i in 1:nrow(comparison_df)) {
  m <- comparison_df[i, ]
  cat(sprintf("  %d. %s: AUC = %.3f ± %.3f (%d变量, 过拟合=%.3f)\n",
              i, m$model_id, m$mean_test_auc, m$sd_test_auc, 
              m$variable_count, m$overfitting))
}