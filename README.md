# Apodemus peninsulae 分布域内的生态因子分析

*大林姬鼠（*Apodemus peninsulae*）分布域内的生态因子分析*

## 项目概述

本仓库包含山东大学本科毕业论文《大林姬鼠分布域内的生态因子分析》的全部分析代码与结果。

**核心目标**：利用物种分布模型（SDM）分析大林姬鼠在当前（1970-2000）及未来（2050s，SSP126/SSP245/SSP585三个排放情景）的生境适宜性变化，为鼠传疾病风险预警提供基础数据。

## 论文信息

- **题目**：大林姬鼠分布域内的生态因子分析
- **作者**：汤泳思
- **学校**：山东大学
- **年份**：2026
- **关键词**：大林姬鼠；物种分布模型；MaxEnt；泊松点过程模型；气候变化；MOP

## 研究流程

```
GBIF 83条存在记录 + WorldClim气候 + FAO土壤 + MODIS纹理 + SRTM地形
        │
        ▼
    数据清洗与稀疏化
        │
        ▼
    变量筛选（26→10变量：相关性分析 → VIF → PCA）
        │
        ▼
    8候选模型嵌套交叉验证（maxnet）
        │
        ▼
    ┌─ S1 MaxEnt模型（最优变量组合）──────────────┐
    │   - 当前生境适宜性预测（AUC = 0.889）       │
    │   - 未来3情景（SSP126/245/585）预测         │
    │                                             │
    └─ PPP泊松点过程模型（组分布交互）──────────┐
        - mgcv::gam(binomial(cloglog))          │
        - bio1 + bio19 × soil_type_merged       │
        - MOP外推风险评估                        │
        └──→ 面积变化分析 ──→ 论文图表 ──→ 结论
```

## 目录结构

```
├── README.md
├── DESCRIPTION                   # R 依赖描述
├── data/
│   ├── README.md                 # 📌 数据获取与处理说明
│   ├── occurrences_raw.csv       # GBIF 83条存在记录
│   ├── occurrences_10var.csv     # 含10环境变量的建模数据
│   └── env_variables.txt         # 变量来源
├── scripts/
│   ├── 00_setup.R                # 📌 路径配置（运行前编辑）
│   ├── 01_data_preparation.R     # 📌 GBIF下载与清洗
│   ├── 02_variable_selection.R   # 📌 变量筛选流程
│   ├── 03_candidate_models.R     # 8候选模型嵌套CV
│   ├── 04_final_model.R          # 最优S1 MaxEnt模型
│   ├── 05_ppp_model.R            # ⭐ PPP泊松点过程（核心模型）
│   ├── 06_predict_current.R      # 当前生境预测
│   ├── 07_mop_analysis.R         # MOP外推风险评估
│   ├── 08_figures.R              # 论文图表生成
│   ├── 09_area_calculation.R     # 📌 面积计算
│   ├── compare_precip_vars.R     # 降水变量对照
│   ├── delong_test.R             # DeLong检验
│   └── pca_analysis.R            # PCA分析
├── results/
│   ├── model_performance.csv     # 候选模型性能
│   ├── candidate_models.csv      # 模型定义
│   ├── figures/                  # 📊 最终论文图（16张PDF）
│   └── area/                     # 📊 面积汇总结果
├── models/
│   └── README.md                 # 模型文件说明
└── paper/                        # 论文源文件
```

## 依赖环境

- **R >= 4.2.0** (推荐 4.3+)
- **Python >= 3.8**（部分可视化脚本）

### R包安装

```r
install.packages(c(
  "terra", "sf", "maxnet", "mgcv", "pROC", "ggplot2",
  "corrplot", "dplyr", "tidyr", "stringr", "rnaturalearth",
  "glmnet", "randomForest"
))
```

## 🚀 复现步骤

### 1. 克隆仓库

```bash
git clone https://github.com/tanggggi/apodemus-sdm.git
cd apodemus-sdm
```

### 2. 配置数据路径

⚠️ **栅格数据（~10GB+）因体积过大未纳入 Git 仓库。**

编辑 `scripts/00_setup.R`，将 `DATA_DIR` 指向你的本地数据目录：

```r
DATA_DIR <- "/你的路径/毕设文件存放/变量筛选"
```

所有分析脚本中的硬编码路径均需根据你的本地环境修改。详见 `data/README.md`。

### 3. 下载环境数据

| 数据 | 来源 | 用途 |
|------|------|------|
| WorldClim 2.1 当前 | https://worldclim.org | bio1/bio12/bio15/bio19 |
| CMIP6 未来 (CNRM-CM6-1) | https://worldclim.org/data/cmip6_2.5m.html | SSP126/SSP245/SSP585 |
| FAO HWSD 土壤 | https://fao.org/soils-portal | 土壤类型分类 |
| MODIS NDVI | https://lpdaac.usgs.gov | 植被纹理指标 |
| SRTM DEM | https://www.gscloud.cn | 海拔数据 |

### 4. 运行分析

推荐按编号顺序执行：

```r
source("scripts/00_setup.R")           # 配置（先改路径）
source("scripts/01_data_preparation.R") # 数据清洗（准备原始GBIF后）
source("scripts/02_variable_selection.R") # 变量筛选
source("scripts/03_candidate_models.R")  # 候选模型比较
source("scripts/04_final_model.R")       # 训练S1模型
source("scripts/05_ppp_model.R")         # ⭐ PPP模型
source("scripts/06_predict_current.R")   # 当前预测
source("scripts/07_mop_analysis.R")      # MOP分析
source("scripts/08_figures.R")           # 出图
source("scripts/09_area_calculation.R")  # 面积计算
```

## 关键发现

### 模型性能

| 模型 | AUC | 变量 | 说明 |
|:-----|:----|:-----|:-----|
| **S1 (MaxEnt)** | **0.889** | 10个环境变量 | 最优MaxEnt |
| **PPP (泊松点过程)** | **交互项ΔAIC=-3.9** | bio1 + bio19 × soil_grp | 最优生态模型 |
| C1 (气候主导) | 0.813 | 仅气候 | 纯气候表现一般 |
| V1 (纹理综合) | 0.757 | 仅纹理 | 纹理单独不够 |

### 未来适生区变化

| 情景 | 当前面积(km²) | 未来面积(km²) | 净变化 |
|:-----|:------------:|:-------------:|:-----:|
| SSP126 | 1,718,060 | 2,277,874 | **+559,814** |
| SSP245 | 1,718,060 | 2,298,466 | **+580,406** |
| SSP585 | 1,718,060 | 3,057,488 | **+1,339,428** |

### 主要结论

1. **所有情景下大林姬鼠适宜区均呈北扩趋势**，SSP585 扩张最显著
2. **年均温（bio1）是最重要的单变量**，贡献 >62%
3. **冬季降水与土壤类型有显著交互作用**：体现了种子埋藏-气味挥发生态机制
4. **MOP分析**表明SSP585存在较高外推风险

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
