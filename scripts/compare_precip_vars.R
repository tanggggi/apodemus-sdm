#!/usr/bin/env Rscript
# 降水变量对照测试：bio19 vs bio12 vs bio15
# 固定 bio1 + soil_type，替换降水变量
# 方法论与原8模型嵌套CV一致

cat(strrep("=", 70), "\n")
cat("降水变量对照测试：bio19 vs bio12 vs bio15\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# 1. 加载包
# ============================================================================
cat("1. 加载包...\n")
required_pkgs <- c("maxnet", "pROC")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}
cat("   完成\n\n")

# ============================================================================
# 2. 读取数据
# ============================================================================
cat("2. 读取数据...\n")
data_file <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/完整数据集_10变量/完整数据集_10变量.csv"
env_data <- read.csv(data_file)
cat("   行数:", nrow(env_data), "\n")
cat("   存在点:", sum(env_data$response == 1), "\n")
cat("   背景点:", sum(env_data$response == 0), "\n\n")

# ============================================================================
# 3. 定义三个模型
# ============================================================================
cat("3. 定义对照模型...\n")
models <- list(
  S1      = list(name = "S1 (bio1+bio19+soil)",     vars = c("bio1", "bio19", "soil_type")),
  T_bio12 = list(name = "T_bio12 (bio1+bio12+soil)", vars = c("bio1", "bio12", "soil_type")),
  T_bio15 = list(name = "T_bio15 (bio1+bio15+soil)", vars = c("bio1", "bio15", "soil_type"))
)

for (m in names(models)) {
  cat("   - ", models[[m]]$name, ": ", 
      paste(models[[m]]$vars, collapse = ", "), "\n", sep = "")
}
cat("\n")

# ============================================================================
# 4. 空间分块（与原8模型一致：k-means 3折）
# ============================================================================
cat("4. 空间分块（k-means 3折）...\n")
set.seed(20260320)
coords <- env_data[, c("longitude", "latitude")]
n_folds <- 3
spatial_clusters <- kmeans(coords, centers = n_folds)
fold_ids <- spatial_clusters$cluster
cat("   分块大小:", paste(table(fold_ids), collapse = ", "), "\n")

folds <- list()
for (k in 1:n_folds) {
  folds[[k]] <- list(
    train = which(fold_ids != k),
    test  = which(fold_ids == k)
  )
}
cat("   完成\n\n")

# ============================================================================
# 5. 嵌套CV函数（与原脚本一致）
# ============================================================================
cat("5. 运行嵌套空间交叉验证...\n")

run_nested_cv <- function(model_name, variables, data, folds) {
  
  model_data <- data[, c(variables, "response")]
  
  test_aucs  <- numeric(length(folds))
  train_aucs <- numeric(length(folds))
  
  for (k in 1:length(folds)) {
    train_data <- model_data[folds[[k]]$train, ]
    test_data  <- model_data[folds[[k]]$test, ]
    
    # 内层：超参数搜索（与原脚本一致）
    param_grid <- expand.grid(
      regmult = c(0.5, 2),
      classes  = c("l", "lq")
    )
    
    best_auc <- -1
    best_params <- NULL
    
    for (p in 1:nrow(param_grid)) {
      inner_folds <- 3
      # 使用随机分块（与原8模型脚本一致）
      inner_fold_ids <- sample(rep(1:inner_folds, length.out = nrow(train_data)))
      
      inner_aucs <- numeric(inner_folds)
      
      for (ifold in 1:inner_folds) {
        inner_train <- train_data[inner_fold_ids != ifold, ]
        inner_test  <- train_data[inner_fold_ids == ifold, ]
        
        tryCatch({
          m <- maxnet(
            p = inner_train$response,
            data = inner_train[, variables, drop = FALSE],
            regmult = param_grid$regmult[p],
            classes = as.character(param_grid$classes[p])
          )
          pred <- predict(m, inner_test[, variables, drop = FALSE], type = "cloglog")
          inner_aucs[ifold] <- auc(roc(inner_test$response, pred, quiet = TRUE))
        }, error = function(e) {
          inner_aucs[ifold] <- NA
        })
      }
      
      mean_auc <- mean(inner_aucs, na.rm = TRUE)
      if (!is.na(mean_auc) && mean_auc > best_auc) {
        best_auc <- mean_auc
        best_params <- list(regmult = param_grid$regmult[p],
                           classes  = as.character(param_grid$classes[p]))
      }
    }
    
    # 外层：用最优参数训练
    final_model <- maxnet(
      p = train_data$response,
      data = train_data[, variables, drop = FALSE],
      regmult = best_params$regmult,
      classes = best_params$classes
    )
    
    train_pred <- predict(final_model, train_data[, variables, drop = FALSE], type = "cloglog")
    test_pred  <- predict(final_model, test_data[, variables, drop = FALSE], type = "cloglog")
    
    train_aucs[k] <- auc(roc(train_data$response, train_pred, quiet = TRUE))
    test_aucs[k]  <- auc(roc(test_data$response, test_pred, quiet = TRUE))
    
    cat("     折叠", k, ": 测试AUC =", round(test_aucs[k], 4), 
        "| 超参数 regmult=", best_params$regmult, "classes=", best_params$classes, "\n")
  }
  
  list(
    name       = model_name,
    variables  = variables,
    test_auc   = mean(test_aucs),
    test_sd    = sd(test_aucs),
    train_auc  = mean(train_aucs),
    overfit    = mean(train_aucs) - mean(test_aucs),
    test_aucs  = test_aucs,
    train_aucs = train_aucs
  )
}

# ============================================================================
# 6. 运行所有模型
# ============================================================================
all_results <- list()

for (m in names(models)) {
  cat("\n  ▶", models[[m]]$name, "\n")
  all_results[[m]] <- run_nested_cv(
    model_name = models[[m]]$name,
    variables   = models[[m]]$vars,
    data        = env_data,
    folds       = folds
  )
}

# ============================================================================
# 7. 汇总结果
# ============================================================================
cat("\n\n", strrep("=", 70), "\n")
cat("结果汇总\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-30s | %-20s | AUC(mean±sd)      | overfit\n", "模型", "降水变量"))
cat(strrep("-", 85), "\n")

for (m in names(models)) {
  r <- all_results[[m]]
  cat(sprintf("%-30s | %-20s | %.4f ± %.4f | %.4f\n",
              r$name,
              setdiff(r$variables, c("bio1", "soil_type")),
              r$test_auc, r$test_sd,
              r$overfit))
}

cat("\n\n--- 折叠级结果 ---\n")
cat(sprintf("%-30s | %s\n", "模型", paste0("折叠", 1:n_folds, " AUC", collapse = " | ")))
cat(strrep("-", 85), "\n")
for (m in names(models)) {
  cat(sprintf("%-30s | %s\n", models[[m]]$name,
              paste(sprintf("%.4f", all_results[[m]]$test_aucs), collapse = " | ")))
}

# ============================================================================
# 8. 保存结果
# ============================================================================
output_dir <- "/Users/hjhj/Desktop/毕设文件存放/变量筛选/降水变量对照测试"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

results_df <- data.frame(
  model     = sapply(all_results, `[[`, "name"),
  variables = sapply(all_results, function(x) paste(x$variables, collapse = ", ")),
  precip_var = sapply(all_results, function(x) {
    setdiff(x$variables, c("bio1", "soil_type"))
  }),
  auc_mean  = sapply(all_results, `[[`, "test_auc"),
  auc_sd    = sapply(all_results, `[[`, "test_sd"),
  overfit   = sapply(all_results, `[[`, "overfit"),
  row.names = NULL
)

write.csv(results_df, file.path(output_dir, "降水变量对照测试_结果汇总.csv"), row.names = FALSE)
cat("\n\n结果已保存至:", output_dir, "\n")

# 打印关键结论
cat("\n\n📋 关键结论:\n")
cat(strrep("-", 50), "\n")
s1_auc <- all_results[["S1"]]$test_auc
for (m in setdiff(names(models), "S1")) {
  r <- all_results[[m]]
  delta <- r$test_auc - s1_auc
  cat(sprintf("  %s vs S1: ΔAUC = %+.4f\n", r$name, delta))
}

cat("\n完成!\n")
