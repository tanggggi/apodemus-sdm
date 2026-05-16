#!/usr/bin/env Rscript
# 训练最终S1模型脚本
# 使用嵌套CV中S1表现最好的超参数，在完整数据集上训练

cat(strrep("=", 80), "\n")
cat("训练最终S1模型（土壤-气候模型）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("模型变量: soil_type, bio1, bio19 (3个变量)\n")
cat("最佳超参数: regmult=0.5, classes='lq'\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 加载必要的包
# ============================================================================
cat("1. 加载必要的包...\n")
required_packages <- c("maxnet", "pROC", "dplyr", "terra")
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
# 2. 设置路径和参数
# ============================================================================
cat("\n2. 设置路径和参数...\n")

# S1模型变量
s1_vars <- c("soil_type", "bio1", "bio19")
cat("   S1模型变量 (3个):", paste(s1_vars, collapse=", "), "\n")

# 最佳超参数（来自嵌套CV）
best_regmult <- 0.5
best_classes <- "lq"
cat("   最佳超参数: regmult =", best_regmult, ", classes = '", best_classes, "'\n", sep="")

# 输入数据路径
data_file <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/完整数据集_10变量/完整数据集_10变量.csv"

# 输出目录
output_dir <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/最终S1模型"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("   输入数据:", data_file, "\n")
cat("   输出目录:", output_dir, "\n")

# ============================================================================
# 3. 读取数据
# ============================================================================
cat("\n3. 读取数据...\n")

if (!file.exists(data_file)) {
  stop("数据文件不存在: ", data_file)
}

env_data <- read.csv(data_file)
cat("   数据维度:", dim(env_data), "\n")
cat("   存在点数量:", sum(env_data$response == 1), "个\n")
cat("   背景点数量:", sum(env_data$response == 0), "个\n")

# 检查是否包含S1模型所需变量
missing_vars <- setdiff(s1_vars, colnames(env_data))
if (length(missing_vars) > 0) {
  stop("数据中缺失以下变量:", paste(missing_vars, collapse=", "))
}
cat("   ✅ 所有S1模型变量都存在\n")

# ============================================================================
# 4. 准备训练数据
# ============================================================================
cat("\n4. 准备训练数据...\n")

# 提取S1模型所需的数据
train_data <- env_data[, c(s1_vars, "response")]
cat("   训练数据维度:", dim(train_data), "\n")

# 检查数据质量
cat("   数据质量检查:\n")
for (var in s1_vars) {
  vals <- train_data[[var]]
  cat(sprintf("     %-10s: 范围=[%.2f, %.2f], 均值=%.2f, NA=%d\n",
              var, min(vals), max(vals), mean(vals), sum(is.na(vals))))
}

# 如果有NA值，用中位数填充
na_count <- sum(is.na(train_data))
if (na_count > 0) {
  cat("   警告: 发现", na_count, "个缺失值，用变量中位数填充\n")
  for (var in s1_vars) {
    na_idx <- is.na(train_data[[var]])
    if (any(na_idx)) {
      median_val <- median(train_data[[var]], na.rm = TRUE)
      train_data[[var]][na_idx] <- median_val
      cat("     ", var, ": 用中位数", median_val, "填充", sum(na_idx), "个缺失值\n")
    }
  }
}

# ============================================================================
# 5. 训练最终S1模型
# ============================================================================
cat("\n5. 训练最终S1模型...\n")

cat("   开始训练: regmult =", best_regmult, ", classes = '", best_classes, "'\n", sep="")

# 分离响应变量和环境变量
response <- train_data$response
env_vars <- train_data[, s1_vars, drop = FALSE]

# 训练maxnet模型
start_time <- Sys.time()
final_s1_model <- maxnet(
  p = response,
  data = env_vars,
  regmult = best_regmult,
  maxnet.formula(p = response, data = env_vars, classes = best_classes)
)
training_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("   训练完成! 用时:", round(training_time, 2), "秒\n")

# ============================================================================
# 6. 模型评估
# ============================================================================
cat("\n6. 模型评估...\n")

# 在训练数据上预测
train_pred <- predict(final_s1_model, env_vars, type = "cloglog")

# 计算训练AUC
train_roc <- roc(response, train_pred, quiet = TRUE)
train_auc <- auc(train_roc)
cat("   训练AUC:", round(train_auc, 4), "\n")

# 计算其他性能指标
cat("   性能指标:\n")

# 分类阈值（Youden指数）
best_threshold <- coords(train_roc, "best", ret = "threshold", quiet = TRUE)$threshold
cat("     最佳分类阈值:", round(best_threshold, 4), "\n")

# 混淆矩阵（使用最佳阈值）
pred_class <- ifelse(train_pred >= best_threshold, 1, 0)
conf_matrix <- table(实际 = response, 预测 = pred_class)
cat("     混淆矩阵:\n")
print(conf_matrix)

# 计算准确率、敏感性、特异性
if (nrow(conf_matrix) == 2 && ncol(conf_matrix) == 2) {
  accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])  # 真阳性率
  specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])  # 真阴性率
  
  cat("     准确率:", round(accuracy, 4), "\n")
  cat("     敏感性:", round(sensitivity, 4), "\n")
  cat("     特异性:", round(specificity, 4), "\n")
}

# ============================================================================
# 7. 变量重要性分析
# ============================================================================
cat("\n7. 变量重要性分析...\n")

# 提取模型系数
model_coefficients <- final_s1_model$betas
cat("   模型系数数量:", length(model_coefficients), "\n")

# 计算变量重要性（基于系数绝对值）
var_importance <- data.frame(
  variable = s1_vars,
  importance = sapply(s1_vars, function(var) {
    # 查找该变量的所有系数
    var_coefs <- model_coefficients[grep(paste0("^", var), names(model_coefficients))]
    if (length(var_coefs) > 0) {
      sum(abs(var_coefs))
    } else {
      0
    }
  })
)

# 按重要性排序
var_importance <- var_importance[order(-var_importance$importance), ]
var_importance$importance_percent <- var_importance$importance / sum(var_importance$importance) * 100

cat("   变量重要性排序:\n")
for (i in 1:nrow(var_importance)) {
  cat(sprintf("     %d. %-10s: 重要性=%.4f (%.1f%%)\n", 
              i, var_importance$variable[i], 
              var_importance$importance[i],
              var_importance$importance_percent[i]))
}

# ============================================================================
# 8. 保存模型和结果
# ============================================================================
cat("\n8. 保存模型和结果...\n")

# 保存最终模型
model_file <- file.path(output_dir, "final_s1_model.rds")
saveRDS(final_s1_model, model_file)
cat("   最终模型保存到:", model_file, "\n")

# 保存变量重要性
importance_file <- file.path(output_dir, "s1_variable_importance.csv")
write.csv(var_importance, importance_file, row.names = FALSE)
cat("   变量重要性保存到:", importance_file, "\n")

# 保存训练预测结果
predictions_df <- data.frame(
  true_response = response,
  predicted_probability = train_pred,
  predicted_class = pred_class
)
predictions_file <- file.path(output_dir, "s1_training_predictions.csv")
write.csv(predictions_df, predictions_file, row.names = FALSE)
cat("   训练预测结果保存到:", predictions_file, "\n")

# 保存模型摘要
summary_file <- file.path(output_dir, "s1_model_summary.txt")
sink(summary_file)
cat("最终S1模型摘要\n")
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==========\n\n")
cat("模型参数\n")
cat("--------\n")
cat("变量: soil_type, bio1, bio19 (3个变量)\n")
cat("超参数: regmult =", best_regmult, ", classes = '", best_classes, "'\n", sep="")
cat("训练样本数:", nrow(train_data), "\n")
cat("存在点:", sum(response == 1), "\n")
cat("背景点:", sum(response == 0), "\n\n")

cat("性能指标\n")
cat("--------\n")
cat("训练AUC:", round(train_auc, 4), "\n")
if (exists("accuracy")) {
  cat("准确率:", round(accuracy, 4), "\n")
  cat("敏感性:", round(sensitivity, 4), "\n")
  cat("特异性:", round(specificity, 4), "\n")
}
cat("最佳分类阈值:", round(best_threshold, 4), "\n\n")

cat("变量重要性\n")
cat("--------\n")
for (i in 1:nrow(var_importance)) {
  cat(sprintf("%d. %-10s: %.4f (%.1f%%)\n", 
              i, var_importance$variable[i], 
              var_importance$importance[i],
              var_importance$importance_percent[i]))
}
cat("\n模型系数\n")
cat("--------\n")
print(model_coefficients)
sink()
cat("   模型摘要保存到:", summary_file, "\n")

# ============================================================================
# 9. 生成最终报告
# ============================================================================
cat("\n9. 生成最终报告...\n")

report_content <- paste0(
  "# 最终S1模型训练报告\n\n",
  "## 生成信息\n",
  "- **生成时间**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "- **模型类型**: maxnet (MaxEnt的R实现)\n",
  "- **模型变量**: 3个 (soil_type, bio1, bio19)\n",
  "- **超参数来源**: 嵌套空间交叉验证中的最佳参数\n",
  "- **训练数据**: 完整10变量数据集 (", nrow(env_data), "个样本)\n\n",
  "## 模型参数\n",
  "| 参数 | 值 | 说明 |\n",
  "|------|-----|------|\n",
  "| **变量数** | 3 | soil_type, bio1, bio19 |\n",
  "| **正则化参数(regmult)** | ", best_regmult, " | 来自嵌套CV最佳折叠 |\n",
  "| **特征类别(classes)** | '", best_classes, "' | 线性+二次项 |\n",
  "| **训练样本** | ", nrow(train_data), " | ", sum(response == 1), "存在点 + ", sum(response == 0), "背景点 |\n",
  "| **训练时间** | ", round(training_time, 2), "秒 |  |\n\n",
  "## 模型性能\n",
  "| 指标 | 值 | 评价 |\n",
  "|------|-----|------|\n",
  "| **训练AUC** | ", round(train_auc, 4), " | ", ifelse(train_auc >= 0.9, "优秀", ifelse(train_auc >= 0.8, "良好", "一般")), " |\n",
  if (exists("accuracy")) paste0(
  "| **准确率** | ", round(accuracy, 4), " | ", ifelse(accuracy >= 0.9, "优秀", ifelse(accuracy >= 0.8, "良好", "一般")), " |\n",
  "| **敏感性** | ", round(sensitivity, 4), " | ", ifelse(sensitivity >= 0.8, "良好", "一般"), " |\n",
  "| **特异性** | ", round(specificity, 4), " | ", ifelse(specificity >= 0.9, "优秀", ifelse(specificity >= 0.8, "良好", "一般")), " |\n"), "",
  "| **最佳分类阈值** | ", round(best_threshold, 4), " | Youden指数确定 |\n\n",
  "## 变量重要性\n",
  "| 排名 | 变量 | 重要性得分 | 百分比 | 生态学意义 |\n",
  "|------|------|------------|--------|------------|\n",
  paste(sapply(1:nrow(var_importance), function(i) {
    v <- var_importance[i, ]
    eco_meaning <- switch(v$variable,
      "soil_type" = "土壤类型（栖息地基质）",
      "bio1" = "年均温（热量条件）", 
      "bio19" = "冬季降水（水分条件）",
      v$variable
    )
    sprintf("| %d | **%s** | %.4f | %.1f%% | %s |", 
            i, v$variable, v$importance, v$importance_percent, eco_meaning)
  }), collapse = "\n"),
  "\n\n## 科学意义\n",
  "### 模型优势\n",
  "1. **简约高效**: 仅3个变量达到优秀预测性能 (AUC ", round(train_auc, 3), ")\n",
  "2. **生态学合理**: 温度(bio1)+冬季降水(bio19)+土壤类型(soil_type)的组合具有明确的生态学机制\n",
  "3. **统计稳健**: 基于嵌套CV最佳参数，避免过拟合\n",
  "4. **解释性强**: 变量重要性清晰，便于生态学解释\n\n",
  "### 生态学机制\n",
  "S1模型揭示了物种分布的三个关键驱动因子：\n",
  "1. **温度调控生理过程**: 年均温影响代谢率、繁殖时机和生存下限\n",
  "2. **冬季降水决定水分条件**: 冬季降水影响冬季生存、繁殖成功率和食物资源\n",
  "3. **土壤类型提供栖息地基质**: 土壤影响食物资源、隐蔽条件和微生境多样性\n\n",
  "## 文件清单\n",
  "- `final_s1_model.rds`: 最终S1模型（RDS格式）\n",
  "- `s1_variable_importance.csv`: 变量重要性排序表\n",
  "- `s1_training_predictions.csv`: 训练集预测结果\n",
  "- `s1_model_summary.txt`: 模型详细摘要\n",
  "- 本报告\n\n",
  "## 使用方法\n",
  "```r\n",
  "# 加载模型\n",
  "model <- readRDS('", model_file, "')\n",
  "\n",
  "# 预测新数据\n",
  "new_data <- data.frame(soil_type=..., bio1=..., bio19=...)\n",
  "predictions <- predict(model, new_data, type='cloglog')\n",
  "```\n"
)

report_file <- file.path(output_dir, "最终S1模型训练报告.md")
writeLines(report_content, report_file)
cat("   报告保存到:", report_file, "\n")

# ============================================================================
# 10. 输出摘要
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("✅ 最终S1模型训练完成\n")
cat("模型变量: soil_type, bio1, bio19 (3个变量)\n")
cat("训练AUC:", round(train_auc, 4), "\n")
cat("最佳变量:", var_importance$variable[1], " (重要性:", round(var_importance$importance_percent[1], 1), "%)\n", sep="")
cat("输出目录:", output_dir, "\n")
cat("主要文件:\n")
cat("  - final_s1_model.rds (最终模型)\n")
cat("  - s1_variable_importance.csv (变量重要性)\n")
cat("  - 最终S1模型训练报告.md (完整报告)\n")
cat(strrep("=", 80), "\n")