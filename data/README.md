# 数据目录说明

本目录包含分析使用的物种出现数据和环境变量描述。

## 文件列表

| 文件 | 说明 | 来源 |
|------|------|------|
| `occurrences_raw.csv` | 原始GBIF出现记录 | [GBIF](https://www.gbif.org/species/2437763) |
| `occurrences_10var.csv` | 清洗后含10个环境变量的出现记录 | 由 `01_data_preparation.R` + `02_variable_selection.R` 生成 |
| `env_variables.txt` | 环境变量来源说明 | 本仓库 |

## 栅格数据（未包含在仓库中）

所有 `.tif` 栅格数据因体积过大（>10GB）未纳入版本控制。
请从以下来源下载并预处理：

### 1. 气候数据
- **来源**: [WorldClim 2.1](https://worldclim.org)
- **分辨率**: 2.5 arc-minutes (~5km)
- **时期**: 当前 (1970-2000) + 未来 (2041-2060)
- **情景**: SSP126, SSP245, SSP585
- **模型**: CNRM-CM6-1
- **变量**: bio1 (年均温), bio12 (年降水量), bio15 (降水季节性), bio19 (最冷季降水)

### 2. 土壤数据
- **来源**: FAO Harmonized World Soil Database (HWSD) v1.2
- **处理**: 栅格化 → 提取土壤类型 → 合并低频类别（<2个存在点）→ "other"组

### 3. 植被纹理指标
- **来源**: MODIS NDVI 时间序列
- **变量**: std (标准差), Contrast (对比度), cv (变异系数), Entropy (熵值)

### 4. 数字高程模型
- **来源**: SRTM 90m (地理空间数据云)
- **处理**: 重采样至2.5 arc-minutes

### 5. 中国边界
- **来源**: [GADM](https://gadm.org) 或自然资源部标准地图

## 统一分辨率

所有栅格最终统一至:
- **分辨率**: 0.0166667° (~2km)
- **范围**: 73.6-134.75°E, 18.22-53°N
- **坐标系**: WGS84 (EPSG:4326)
- **投影**: Albers等面积 (EPSG: China Albers, 用于面积计算)

## 运行脚本前

编辑 `scripts/00_setup.R` 中的 `DATA_DIR` 路径指向你的本地数据目录。
