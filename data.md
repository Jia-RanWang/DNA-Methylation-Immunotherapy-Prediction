# 原始数据说明
## 1. 数据来源
本研究所使用的数据来源于美国国家癌症研究所（National Cancer Institute，NCI）建立的癌症基因组图谱数据库（The Cancer Genome Atlas，TCGA），全部数据通过 Genomic Data Commons（GDC）官方数据平台获取。
- GDC 数据门户：<https://portal.gdc.cancer.gov/>
- TCGA 项目页面：<https://portal.gdc.cancer.gov/projects/TCGA-XXX>
- R语言下载工具：`TCGAbiolinks`
- TCGA项目编号：`TCGA-XXX`
- 数据下载日期：按需自行填写

由于TCGA/GDC原始文件数量庞大、单文件体积巨大，同时部分数据受NCI访问规范约束，**本仓库不存储、不上传任何原始测序/组学数据文件**。
原始数据可自行前往GDC平台公开下载，搭配本仓库完整代码完成下载、质控、筛选与全套生物信息分析流程。

## 2. 原始数据类型
### 2.1 基因表达数据
用途：
- 构建免疫基因表达特征GEP；
- 计算TGF-β通路活性评分；
- 差异表达分析；
- 机器学习特征构建；
- 功能富集通路解析。
标准化处理、注释版本、过滤逻辑均以仓库R代码内参数为准。

### 2.2 体细胞突变数据
采用MAF（Mutation Annotation Format）标准突变注释文件，核心用于**TMB肿瘤突变负荷**计算。
仅统计体细胞非同义突变用于指标换算，突变筛选规则、TMB计算公式完整写在项目代码中。

### 2.3 临床及样本注释数据
用途：
- 区分原发肿瘤/正常/转移样本；
- 多组学样本ID匹配；
- 构建分组标签、临床分层；
- 样本质控与异常剔除。

### 2.4 通路与基因集注释数据
内置固定免疫基因集（GEP、TGF-β）、GO/KEGG/Reactome注释库，用于评分计算与富集分析，注释版本、基因列表统一固化在脚本内。

## 3. 原始数据不上传仓库说明
本仓库完全不包含以下大型原始文件：
- TCGA原始甲基化、表达矩阵；
- MAF突变注释文件；
- 临床样本原始表格；
- FASTQ/BAM原始测序reads；
- IDAT甲基化芯片原始文件；
- GDC批量下载的大型归档文件。

### 不上传原始数据核心原因
1. 原始组学数据单文件可达数十GB，GitHub仓库存在存储容量限制；
2. TCGA/GDC数据受美国NCI数据使用规范约束，不适合开源仓库直接分发；
3. 所有原始数据均可通过官方渠道免费公开获取，无需随代码分发；
4. 仓库提供完整下载+预处理代码，可一键复现数据整理流程；
5. 海量原始文件会大幅拉低代码拉取、克隆、浏览速度。

> 补充声明：本研究数据集不含任何患者个人可识别隐私信息；文中“原始数据”仅指代公共数据库下载的标准化组学输入文件，不含测序仪原始输出文件。

## 4. 完整数据下载&分析执行流程
项目按模块化拆分分析脚本，Linux/macOS终端可直接批量运行：
```bash
Rscript code/01_download_data.R
Rscript code/02_preprocess_and_labels.R
Rscript code/03_limma_diff.R
Rscript code/04_feature_selection.R
Rscript code/05_ML_model.R
Rscript code/06_SHAP_analysis.R
Rscript code/07_enrichment.R
Rscript code/08_plot_figures.R
```

脚本功能对照表：
| 脚本 | 核心功能 |
|---|---|
| `01_download_data.R` | GDC批量下载甲基化/表达/突变/临床数据 |
| `02_preprocess_and_labels.R` | 样本筛选、多组学匹配、TCGA虚拟响应标签构建 |
| `03_limma_diff.R` | 甲基化差异位点DMP筛选、火山图绘制 |
| `04_feature_selection.R` | LASSO/RF-RFE/XGBoost等多算法特征筛选 |
| `05_ML_model.R` | 8类机器学习训练、交叉验证、模型评估 |
| `06_SHAP_analysis.R` | Kernel SHAP模型可解释性分析、关键CpG挖掘 |
| `07_enrichment.R` | GO/KEGG/Reactome/DO多数据库富集 |
| `08_plot_figures.R` | 论文全套可视化图表输出 |

> 依赖关系：`01_download_data.R` 为前置步骤，其余脚本依赖下载完成后的raw原始数据。

## 5. 项目标准数据目录结构
```text
data/
├── raw/                # GDC下载原始数据，不上传GitHub
│   ├── expression/
│   ├── mutation/
│   ├── clinical/
│   └── sample_annotation/
└── processed/         # 清洗、校正后的中间特征矩阵
    ├── expression/
    ├── methylation/
    ├── mutation/
    ├── clinical/
    └── feature_matrix/
```

目录说明：
| 目录 | 说明 |
|---|---|
| `data/raw/` | 原始未处理GDC下载文件，被.gitignore屏蔽 |
| `data/processed/` | 标准化、批次校正、插补后可直接建模矩阵 |
| `data/raw/methylation/` | TCGA 450K原始甲基化β矩阵 |
| `data/raw/mutation/` | MAF体细胞突变文件 |
| `data/processed/feature_matrix/` | 筛选完成的CpG特征集（建模输入） |

## 6. 数据可复现性说明
GDC数据库会持续更新样本、注释与工作流版本，**不同时间下载数据会造成样本量、指标轻微浮动**。
如需100%复现论文结果，建议完整记录以下信息留存备查：
1. TCGA官方项目编号
2. GDC数据精确下载日期
3. GDC下载manifest清单文件ID
4. 芯片注释库版本（hg19 450K/850K）
5. R主版本、Bioconductor版本
6. 全部生信/机器学习R包版本
7. 样本剔除筛选标准
8. TMB/GEP/TGF-β分位数阈值

在完全一致的数据集、软件环境、筛选规则下，可完整复现论文全部分析流程；若数据库/软件版本发生迭代，模型AUC、富集条目会出现小幅数值偏差。

## 7. 环境与下载元数据保存规范
推荐运行代码后自动导出环境记录，存放路径：
```text
results/
└── metadata/
    ├── raw_data_metadata.json  # 数据来源、下载时间记录
    ├── gdc_manifest.tsv        # GDC下载文件清单
    ├── package_versions.txt    # 全部依赖包版本
    └── sessionInfo.txt         # R完整会话环境
```

### 导出R会话信息代码
```r
# 保存完整R环境信息
writeLines(
  capture.output(sessionInfo()),
  "results/metadata/sessionInfo.txt"
)
```

### 批量记录核心包版本代码
```r
packages <- c(
  "TCGAbiolinks", "dplyr", "data.table", "stringr", "limma",
  "glmnet", "randomForest", "xgboost", "kernelshap", "clusterProfiler"
)

package_versions <- data.frame(
  package = packages,
  version = vapply(
    packages,
    function(x) if (requireNamespace(x, quietly = TRUE)) as.character(packageVersion(x)) else NA,
    character(1)
  ),
  row.names = NULL
)

write.table(
  package_versions,
  "results/metadata/package_versions.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)
```

## 8. .gitignore 原始数据屏蔽配置
复制以下内容到项目根目录`.gitignore`，防止误上传超大原始数据：
```gitignore
# TCGA/GDC 原始组学数据（禁止上传）
data/raw/*
!data/raw/.gitkeep

# 大型中间特征矩阵
data/processed/*.rds
data/processed/*.RData
data/processed/*.csv
data/processed/*.tsv

# 全部输出结果文件夹
results/

# R运行缓存文件
*.RData
*.Rhistory
renv/
.Rprofile

# 系统隐藏文件
.DS_Store
Thumbs.db
```

### 保留空目录命令（终端执行）
```bash
mkdir -p data/raw data/processed results/metadata
touch data/raw/.gitkeep
```
`.gitkeep`文件会被Git识别，仅保留文件夹结构，不存储内部原始数据。

## 9. 数据获取最终声明
本项目所有TCGA组学数据均可通过NCI GDC官网免费公开获取，仓库不提供任何原始数据分发服务。
仓库完整配套**数据下载、质控、特征筛选、建模、可解释性、富集绘图**全套可运行R代码，仅需自行下载原始数据即可完整复现论文全部分析。
论文所有图表、指标均基于本研究下载时段的GDC数据版本；若后续数据库、注释包、R包更新，重复分析结果会存在小幅数值差异。

