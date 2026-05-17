# Apodemus peninsulae 分布域内的生态因子分析

*大林姬鼠（*Apodemus peninsulae*）分布域内的生态因子分析*

## 项目概述

本仓库包含山东大学本科毕业论文《大林姬鼠分布域内的生态因子分析》的全部分析代码与数据。

**核心目标**：利用物种分布模型（SDM）分析大林姬鼠在当前及未来（SSP126、SSP245、SSP585三个排放情景，2050年代）的生境适宜性变化，为鼠传疾病风险预警提供基础数据。

## 论文信息

- **题目**：大林姬鼠分布域内的生态因子分析
- **作者**：汤泳思
- **学校**：山东大学
- **年份**：2026
- **关键词**：大林姬鼠；物种分布模型；MaxEnt；泊松点过程；气候变化；MOP

## 研究流程

```
GBIF物种数据 + WorldClim气候数据 + FAO土壤数据 + MODIS纹理数据
        │
        ▼
    数据清洗与预处理
        │
        ▼
    变量筛选（相关性分析 → VIF → PCA → 10变量）
        │
        ▼
    8候选模型嵌套交叉验证（MaxNet）
        │
        ▼
    S1模型（最优模型：气候+土壤+纹理）
        │
        ├──→ 当前生境适宜性预测
        │
        ├──→ 未来3情景预测（SSP126/245/585）
        │
        └──→ PPP泊松点过程模型（soil_type交互项）
                │
                └──→ MOP外推风险评估
```

## 目录结构

```
├── README.md                          # 本文件
├── DESCRIPTION                        # R 依赖描述
├── data/                              
│   ├── occurrences_raw.csv            # 89条大林姬鼠存在记录（GBIF）
│   ├── occurrences_10var.csv          # 10变量完整建模数据集（存在+背景点）
│   └── env_variables.txt              # 环境变量来源与处理说明
├── scripts/                          
│   ├── 03_candidate_models.R          # 核心：8候选模型 + 嵌套交叉验证
│   ├── 04_final_model.R              # 最优S1模型训练
│   ├── 05_ppp_model.R                # PPP泊松点过程模型（组分交互）
│   ├── 06_predict_current.R          # 当前生境适宜性预测
│   ├── 07_mop_analysis.R             # MOP外推风险评估
│   ├── 08_figures.R                  # 论文图片生成
│   ├── delong_test.R                 # DeLong检验（模型间AUC差异）
│   ├── pca_analysis.R                # PCA主成分分析
│   └── compare_precip_vars.R         # 降水变量对比验证
├── results/
│   ├── model_performance.csv         # 各候选模型性能指标
│   └── candidate_models.csv          # 候选模型变量组合定义
└── paper/
    └── (论文PDF将在此处)
```

## 依赖环境

- **R >= 4.2.0**
- 关键 R 包：

```r
install.packages(c(
  "terra",      # 栅格数据处理
  "maxnet",     # MaxEnt模型（纯R实现）
  "mgcv",       # PPP模型的GAM框架
  "pROC",       # ROC曲线与AUC
  "ggplot2",    # 数据可视化
  "corrplot",   # 相关性热图
  "dplyr",      # 数据操作
  "glmnet",     # LASSO/弹性网（M1参考模型）
  "randomForest", # 随机森林（T1参考模型）
  "sf",         # 空间矢量数据
  "rnaturalearth" # 中国边界数据
))
```

## 复现步骤

1. **克隆仓库**
   ```bash
   git clone https://github.com/汤泳思用户名/apodemus-sdm.git
   cd apodemus-sdm
   ```

2. **安装R依赖**（见上）

3. **下载环境数据**
   从 [WorldClim](https://worldclim.org/data/cmip6_2.5m.html) 下载当前和未来气候数据，详见 `data/env_variables.txt`

4. **运行分析**
   ```r
   source("scripts/03_candidate_models.R")    # 候选模型比较
   source("scripts/04_final_model.R")         # 训练S1
   source("scripts/05_ppp_model.R")           # PPP建模
   source("scripts/06_predict_current.R")     # 当前预测
   source("scripts/07_mop_analysis.R")        # MOP分析
   source("scripts/08_figures.R")             # 出图
   ```

## 关键发现

| 模型 | AUC | 特征 | 说明 |
|:-----|:----|:-----|:-----|
| **S1** | **0.917** | bio1 + bio19 + soil_type + 纹理 | **最优模型** |
| C1 | 0.813 | 仅气候变量 | 纯气候表现一般 |
| V1 | 0.757 | 仅纹理变量 | 纹理单独不够 |

- **最关键的预测变量**：年均温（bio1），置换重要性 >62%
- **未来趋势**：三个情景下适生区均呈扩张趋势，SSP585 扩张最显著
- **冬季条件**：最冷季降水（bio19）× 土壤类型交互项捕获了冬季种子埋藏-气味挥发的生态机制

## 引用

```bibtex
@thesis{tang2026apodemus,
  author  = {汤泳思},
  title   = {大林姬鼠分布域内的生态因子分析},
  school  = {山东大学},
  year    = {2026}
}
```

## 许可证

MIT License
