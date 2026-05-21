# 模型文件

由于体积较大（.rds 文件通常 50-200MB），模型文件未纳入版本控制。

## 最终模型

| 模型 | 文件 | 算法 | AUC |
|------|------|------|:---:|
| S1 (最终MaxEnt) | `final_s1_model.rds` | maxnet | 0.889 |
| PPP (泊松点过程) | `ppp_interaction_model.rds` | mgcv::gam(binomial(cloglog)) | 0.864 |

## 如何复现

运行以下脚本即可重新训练：

1. **S1 MaxEnt**: `scripts/04_final_model.R`
2. **PPP**: `scripts/05_ppp_model.R`

训练后的模型会自动保存到本地输出目录。
