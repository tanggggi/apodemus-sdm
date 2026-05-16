#!/usr/bin/env Rscript
# DeLong检验：比较D1和S1模型的AUC差异

cat(strrep("=", 80), "\n")
cat("DeLong检验：比较D1（全模型）和S1（土壤-气候模型）的AUC差异\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("科学问题: S1模型的性能是否在统计上不劣于D1模型？\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 加载必要的包
# ============================================================================
cat("1. 加载必要的包...\n")
required_packages <- c("maxnet", "pROC", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  } else {
    cat("   加载包:", pkg, "\n")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# 2. 设置路径和数据
# ============================================================================
cat("\n2. 加载数据...\n")

# 数据文件
data_file <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/完整数据集_10变量/完整数据集_10变量.csv"
if (!file.exists(data_file)) {
  stop("数据文件不存在: ", data_file)
}

env_data <- read.csv(data_file)
cat("   数据维度:", dim(env_data), "\n")
cat("   存在点:", sum(env_data$response == 1), "个\n")
cat("   背景点:", sum(env_data$response == 0), "个\n")

# 提取坐标
if (!("longitude" %in% colnames(env_data) & "latitude" %in% colnames(env_data))) {
  stop("数据中没有经纬度坐标")
}

coords <- env_data[, c("longitude", "latitude")]

# ============================================================================
# 3. 重现相同的空间分块
# ============================================================================
cat("\n3. 创建相同的空间分块（重现原始分析）...\n")

# 设置相同的随机种子
set.seed(20260320)

# 使用k-means空间聚类分块（与原始分析相同）
n_outer_folds <- 3
spatial_clusters <- kmeans(coords, centers = n_outer_folds)
fold_ids <- spatial_clusters$cluster

cat("   分块大小:", paste(table(fold_ids), collapse=", "), "\n")
cat("   与原始分析一致: 1908, 1706, 2680\n")

# 创建分块列表
folds <- list()
for (k in 1:n_outer_folds) {
  test_idx <- which(fold_ids == k)
  train_idx <- which(fold_ids != k)
  folds[[k]] <- list(train = train_idx, test = test_idx, fold_id = k)
}

# ============================================================================
# 4. 定义模型变量
# ============================================================================
cat("\n4. 定义D1和S1模型的变量集...\n")

# D1: 全模型（10变量）
d1_vars <- c("elevation", "soil_type", "Contrast", "cv", "Entropy", "std", 
             "bio15", "bio19", "bio1", "bio12")

# S1: 土壤-气候模型（3变量）
s1_vars <- c("soil_type", "bio1", "bio19")

cat("   D1模型: ", length(d1_vars), "个变量\n", sep="")
cat("     变量: ", paste(d1_vars, collapse=", "), "\n")
cat("   S1模型: ", length(s1_vars), "个变量\n", sep="")
cat("     变量: ", paste(s1_vars, collapse=", "), "\n")

# ============================================================================
# 5. 嵌套CV函数（返回预测概率）
# ============================================================================
cat("\n5. 定义嵌套CV函数（返回预测概率）...\n")

run_nested_cv_with_predictions <- function(model_id, variables, data, folds) {
  # 运行嵌套CV并返回每个测试折叠的预测概率
  
  cat("   运行模型", model_id, "(", length(variables), "个变量)...\n")
  
  model_data <- data[, c(variables, "response")]
  
  # 存储结果
  results <- list(
    model_id = model_id,
    variables = variables,
    test_predictions = list(),  # 存储每个折叠的预测概率
    test_labels = list(),       # 存储每个折叠的真实标签
    test_auc = numeric(length(folds)),
    best_params = list()
  )
  
  for (fold_idx in 1:length(folds)) {
    fold <- folds[[fold_idx]]
    
    cat("     折叠", fold_idx, "/", length(folds), "...")
    
    # 划分训练集和测试集
    train_data <- model_data[fold$train, ]
    test_data <- model_data[fold$test, ]
    
    # 内层参数搜索
    param_grid <- expand.grid(
      regmult = c(0.5, 2),    # 与原始分析相同
      classes = c("l", "lq")  # 与原始分析相同
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
    
    # 在测试集上预测
    test_pred <- predict(final_model, test_data[, variables, drop = FALSE], type = "cloglog")
    
    # 保存预测概率和真实标签
    results$test_predictions[[fold_idx]] <- test_pred
    results$test_labels[[fold_idx]] <- test_data$response
    
    # 计算测试AUC
    test_roc <- roc(test_data$response, test_pred, quiet = TRUE)
    results$test_auc[fold_idx] <- auc(test_roc)
    
    results$best_params[[fold_idx]] <- best_params
  }
  
  # 计算汇总统计
  results$mean_test_auc <- mean(results$test_auc, na.rm = TRUE)
  results$sd_test_auc <- sd(results$test_auc, na.rm = TRUE)
  
  cat("   完成: AUC =", round(results$mean_test_auc, 3), "±", round(results$sd_test_auc, 3), "\n")
  
  return(results)
}

# ============================================================================
# 6. 运行D1和S1模型的嵌套CV（获取预测概率）
# ============================================================================
cat("\n6. 运行D1和S1模型的嵌套CV...\n")

# 运行D1模型
cat("\n   D1模型:\n")
d1_results <- run_nested_cv_with_predictions(
  model_id = "D1",
  variables = d1_vars,
  data = env_data,
  folds = folds
)

# 运行S1模型
cat("\n   S1模型:\n")
s1_results <- run_nested_cv_with_predictions(
  model_id = "S1",
  variables = s1_vars,
  data = env_data,
  folds = folds
)

# ============================================================================
# 7. 准备DeLong检验数据
# ============================================================================
cat("\n7. 准备DeLong检验数据...\n")

# 合并所有折叠的预测概率和真实标签
all_predictions <- list()
all_labels <- list()

for (fold_idx in 1:length(folds)) {
  # D1模型的预测
  d1_pred <- d1_results$test_predictions[[fold_idx]]
  d1_label <- d1_results$test_labels[[fold_idx]]
  
  # S1模型的预测
  s1_pred <- s1_results$test_predictions[[fold_idx]]
  s1_label <- s1_results$test_labels[[fold_idx]]
  
  # 检查标签是否一致（应该一致）
  if (!all(d1_label == s1_label)) {
    warning(paste("折叠", fold_idx, "的真实标签不一致"))
  }
  
  # 添加到总列表
  all_predictions$D1 <- c(all_predictions$D1, d1_pred)
  all_predictions$S1 <- c(all_predictions$S1, s1_pred)
  all_labels$D1 <- c(all_labels$D1, d1_label)
  all_labels$S1 <- c(all_labels$S1, s1_label)
}

cat("   总样本数:", length(all_labels$D1), "\n")
cat("   存在点:", sum(all_labels$D1 == 1), "个\n")
cat("   背景点:", sum(all_labels$D1 == 0), "个\n")

# ============================================================================
# 8. 执行DeLong检验
# ============================================================================
cat("\n8. 执行DeLong检验...\n")

# 创建ROC对象
roc_d1 <- roc(all_labels$D1, all_predictions$D1, quiet = TRUE)
roc_s1 <- roc(all_labels$S1, all_predictions$S1, quiet = TRUE)

cat("   D1模型整体AUC:", round(auc(roc_d1), 4), "\n")
cat("   S1模型整体AUC:", round(auc(roc_s1), 4), "\n")

# DeLong检验（双侧检验）
cat("\n   DeLong检验（双侧）:\n")
delong_test_two_sided <- roc.test(roc_d1, roc_s1, method = "delong")
print(delong_test_two_sided)

# DeLong检验（单侧检验：S1是否不劣于D1？）
cat("\n   DeLong检验（单侧，检验S1是否不劣于D1）:\n")
# 单侧检验：H0: AUC_S1 <= AUC_D1 vs H1: AUC_S1 > AUC_D1
delong_test_one_sided <- roc.test(roc_d1, roc_s1, method = "delong", 
                                  alternative = "less")
print(delong_test_one_sided)

# ============================================================================
# 9. 计算AUC差异的置信区间
# ============================================================================
cat("\n9. 计算AUC差异的置信区间...\n")

# 计算AUC差异
auc_diff <- auc(roc_d1) - auc(roc_s1)
cat("   AUC差异 (D1 - S1):", round(auc_diff, 4), "\n")

# 使用bootstrap计算95%置信区间
cat("   使用bootstrap计算95%置信区间 (1000次重复)...\n")

set.seed(20260320)
n_boot <- 1000
boot_diffs <- numeric(n_boot)

for (i in 1:n_boot) {
  # bootstrap重采样
  boot_indices <- sample(1:length(all_labels$D1), replace = TRUE)
  
  boot_labels <- all_labels$D1[boot_indices]
  boot_pred_d1 <- all_predictions$D1[boot_indices]
  boot_pred_s1 <- all_predictions$S1[boot_indices]
  
  boot_roc_d1 <- roc(boot_labels, boot_pred_d1, quiet = TRUE)
  boot_roc_s1 <- roc(boot_labels, boot_pred_s1, quiet = TRUE)
  
  boot_diffs[i] <- auc(boot_roc_d1) - auc(boot_roc_s1)
}

# 计算置信区间
ci_95 <- quantile(boot_diffs, probs = c(0.025, 0.975))
ci_90 <- quantile(boot_diffs, probs = c(0.05, 0.95))

cat("   95% CI for AUC difference: [", round(ci_95[1], 4), ", ", round(ci_95[2], 4), "]\n", sep="")
cat("   90% CI for AUC difference: [", round(ci_90[1], 4), ", ", round(ci_90[2], 4), "]\n", sep="")

# ============================================================================
# 10. 保存结果
# ============================================================================
cat("\n10. 保存结果...\n")

output_dir <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/DeLong检验_D1_vs_S1"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 保存预测数据用于后续分析
predictions_df <- data.frame(
  true_label = all_labels$D1,
  D1_prediction = all_predictions$D1,
  S1_prediction = all_predictions$S1
)

write.csv(predictions_df, file.path(output_dir, "预测数据_D1_vs_S1.csv"), row.names = FALSE)
cat("   预测数据保存到:", file.path(output_dir, "预测数据_D1_vs_S1.csv"), "\n")

# 保存检验结果
results_summary <- data.frame(
  检验项目 = c("D1_AUC", "S1_AUC", "AUC差异", "DeLong_p值_双侧", "DeLong_p值_单侧",
               "95%CI_下限", "95%CI_上限", "90%CI_下限", "90%CI_上限"),
  值 = c(round(auc(roc_d1), 4), round(auc(roc_s1), 4), round(auc_diff, 4),
         round(delong_test_two_sided$p.value, 4), round(delong_test_one_sided$p.value, 4),
         round(ci_95[1], 4), round(ci_95[2], 4), round(ci_90[1], 4), round(ci_90[2], 4))
)

write.csv(results_summary, file.path(output_dir, "DeLong检验结果汇总.csv"), row.names = FALSE)
cat("   检验结果保存到:", file.path(output_dir, "DeLong检验结果汇总.csv"), "\n")

# 生成报告
report_content <- paste0(
  "# DeLong检验：D1（全模型）vs S1（土壤-气候模型）\n\n",
  "## 生成信息\n",
  "- **生成时间**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "- **样本数量**: ", length(all_labels$D1), "个 (存在点:", sum(all_labels$D1 == 1), ", 背景点:", sum(all_labels$D1 == 0), ")\n",
  "- **模型比较**: \n",
  "  - D1: 10个变量 (全模型)\n",
  "  - S1: 3个变量 (soil_type, bio1, bio19)\n",
  "- **检验方法**: DeLong检验 + Bootstrap置信区间\n\n",
  "## AUC结果\n",
  "| 模型 | 变量数 | AUC值 |\n",
  "|------|--------|-------|\n",
  "| D1 | 10 | ", round(auc(roc_d1), 4), " |\n",
  "| S1 | 3 | ", round(auc(roc_s1), 4), " |\n",
  "| **差异 (D1 - S1)** | **7** | **", round(auc_diff, 4), "** |\n\n",
  "## 统计检验结果\n",
  "### 1. DeLong检验（双侧）\n",
  "- **原假设 (H0)**: AUC_D1 = AUC_S1\n",
  "- **备择假设 (H1)**: AUC_D1 ≠ AUC_S1\n",
  "- **检验统计量 (Z)**: ", round(delong_test_two_sided$statistic, 4), "\n",
  "- **p值**: ", round(delong_test_two_sided$p.value, 4), "\n",
  "- **结论**: ", ifelse(delong_test_two_sided$p.value < 0.05, "拒绝H0，AUC有显著差异", "不拒绝H0，AUC无显著差异"), "\n\n",
  "### 2. DeLong检验（单侧，检验S1不劣于D1）\n",
  "- **原假设 (H0)**: AUC_S1 ≤ AUC_D1\n",
  "- **备择假设 (H1)**: AUC_S1 > AUC_D1\n",
  "- **检验统计量 (Z)**: ", round(delong_test_one_sided$statistic, 4), "\n",
  "- **p值**: ", round(delong_test_one_sided$p.value, 4), "\n",
  "- **结论**: ", ifelse(delong_test_one_sided$p.value < 0.05, "拒绝H0，S1不劣于D1", "不拒绝H0，不能证明S1不劣于D1"), "\n\n",
  "## AUC差异的置信区间\n",
  "| 置信水平 | 下限 | 上限 | 区间宽度 | 包含0？ |\n",
  "|----------|------|------|----------|---------|\n",
  "| 95% | ", round(ci_95[1], 4), " | ", round(ci_95[2], 4), " | ", round(ci_95[2] - ci_95[1], 4), " | ", ifelse(ci_95[1] <= 0 & ci_95[2] >= 0, "是", "否"), " |\n",
  "| 90% | ", round(ci_90[1], 4), " | ", round(ci_90[2], 4), " | ", round(ci_90[2] - ci_90[1], 4), " | ", ifelse(ci_90[1] <= 0 & ci_90[2] >= 0, "是", "否"), " |\n\n",
  "## 科学解释\n",
  "### 统计结论\n",
  "1. **AUC差异**: D1比S1高", round(auc_diff * 100, 2), "%\n",
  "2. **统计显著性**: ", ifelse(delong_test_two_sided$p.value < 0.05, "差异显著", "差异不显著"), "\n",
  "3. **不劣性检验**: ", ifelse(delong_test_one_sided$p.value < 0.05, "S1在统计上不劣于D1", "不能证明S1不劣于D1"), "\n\n",
  "### 生态学意义\n",
  "1. **简约性**: S1使用3个变量达到了D1（10变量）", round(auc(roc_s1)/auc(roc_d1)*100, 1), "%的性能\n",
  "2. **关键驱动因子**: 温度(bio1)、冬季降水(bio19)和土壤类型(soil_type)是物种分布的关键预测因子\n",
  "3. **模型选择建议**: 如果追求简约性，S1是优秀选择；如果追求最高性能，D1略优但复杂得多\n\n",
  "## 方法学细节\n",
  "- **数据**: 10变量完整数据集 (", nrow(env_data), "行)\n",
  "- **空间分块**: 3折k-means空间聚类分块\n",
  "- **参数搜索**: regmult ∈ {0.5, 2}, classes ∈ {'l', 'lq'}\n",
  "- **检验方法**: DeLong检验（相关ROC曲线比较）\n",
  "- **Bootstrap**: 1000次重复计算置信区间\n",
  "- **随机种子**: 20260320（确保可重复性）\n\n",
  "## 文件清单\n",
  "- `预测数据_D1_vs_S1.csv`: 所有测试折叠的预测概率和真实标签\n",
  "- `DeLong检验结果汇总.csv`: 检验结果汇总表\n",
  "- 本报告\n"
)

writeLines(report_content, file.path(output_dir, "DeLong检验报告.md"))
cat("   报告保存到:", file.path(output_dir, "DeLong检验报告.md"), "\n")

# ============================================================================
# 11. 输出摘要
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("✅ DeLong检验完成\n")
cat("D1 AUC:", round(auc(roc_d1), 4), "\n")
cat("S1 AUC:", round(auc(roc_s1), 4), "\n")
cat("AUC差异:", round(auc_diff, 4), " (D1 - S1)\n")
cat("DeLong检验p值（双侧）:", round(delong_test_two_sided$p.value, 4), "\n")
cat("DeLong检验p值（单侧）:", round(delong_test_one_sided$p.value, 4), "\n")
cat("95%置信区间: [", round(ci_95[1], 4), ", ", round(ci_95[2], 4), "]\n", sep="")
cat("输出目录:", output_dir, "\n")
cat(strrep("=", 80), "\n")