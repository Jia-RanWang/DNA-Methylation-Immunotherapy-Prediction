# DNA-Methylation-Immunotherapy-Prediction
Reproducible R workflow for cancer immunotherapy response prediction based on DNA methylation, combining 9 GEO EPIC 850K clinical cohorts and TCGA pan-cancer methylation dataset, covering feature selection, machine learning modeling and epigenetic immune escape functional enrichment.


# 基于DNA甲基化异常的癌症免疫治疗反应预测
[![R >= 4.2.0](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-latest-green)](https://bioconductor.org/)
[![License MIT](https://img.shields.io/badge/license-MIT-orange)](LICENSE)
[![Paper HIT 2026](https://img.shields.io/badge/Thesis-Harbin%20Institute%20of%20Technology-red)](https://hit.edu.cn/)

## 项目简介
本仓库为哈尔滨工业大学生物工程本科毕业论文《基于DNA甲基化异常的癌症免疫治疗反应预测》完整可复现生物信息分析代码，包含全部标准化分析流程、机器学习建模、模型可解释性(SHAP)、多维度功能富集全套R脚本。

### 研究核心创新点
1. **跨平台泛癌数据集整合**：融合Illumina 450K(TCGA)与EPIC 850K(GEO)两套甲基化芯片数据，使用ComBat、BMIQ消除平台/批次偏差；
2. **虚拟响应标签扩充策略**：基于TMB、GEP、TGF-β三大免疫标志物构建TCGA免疫治疗虚拟标签，弥补临床真实响应样本稀缺问题；
3. **多层级特征筛选框架**：limma差异甲基化初筛 → LASSO正则精筛 → 6种机器学习特征筛选横向对比，验证**RF-RFE**为最优CpG筛选算法；
4. **多模型泛化评估**：搭建8类经典机器学习分类模型，采用三分区（训练集/验证集/独立外部测试集）严格规避数据泄露；
5. **表观免疫逃逸机制解析**：基于Kernel SHAP挖掘cg09115473、cg16313807等核心CpG，阐明**HLA家族高甲基化抑制抗原呈递**是免疫治疗耐药核心表观机制。

## 目录
- [1 数据集与数据来源](#1-数据集与数据来源)
- [2 系统软硬件环境](#2-系统软硬件环境)
- [3 依赖R包完整清单&安装脚本](#3-依赖r包完整清单安装脚本)
- [4 完整标准化分析流程](#4-完整标准化分析流程)
- [5 代码文件说明](#5-代码文件说明)
- [6 输入文件规范](#6-输入文件规范)
- [7 输出结果目录结构](#7-输出结果目录结构)
- [8 复现完整运行步骤](#8-复现完整运行步骤)
- [9 模型核心指标结果摘要](#9-模型核心指标结果摘要)
- [10 项目限制与注意事项](#10-项目限制与注意事项)
- [11 后续拓展研究方向](#11-后续拓展研究方向)
- [12 论文引用规范](#12-论文引用规范)
- [13 开源许可](#13-开源许可)

## 1 数据集与数据来源
### 1.1 GEO数据集（含临床真实免疫治疗响应标签，Anti-PD-1/±CTLA4）
共9套Illumina EPIC 850K甲基化数据集：
1. GSE119144：非小细胞肺癌(NSCLC) Anti-PD-1
2. GSE126043：非小细胞肺癌(NSCLC) Anti-PD-1
3. GSE264158：黑色素瘤(SKCM) Anti-PD-1
4. GSE175699：黑色素瘤(SKCM) Anti-PD-1
5. GSE181781：黑色素瘤(SKCM) Anti-PD-1 ± CTLA-4
6. GSE235122：黑色素瘤(SKCM) Anti-PD-1
7. GSE172468：肉瘤(SARC) Anti-PD-1
8. GSE277573：头颈部鳞癌(HNSC) Anti-PD-1
9. GSE305240：胃腺癌(STAD) Anti-PD-1

> 独立外部测试集拆分：GSE119144、GSE277573完全脱离训练流程，仅用于最终泛化性能盲验证；剩余7套GEO参与训练/验证集分层抽样。

### 1.2 TCGA泛癌数据集（Illumina 450K，构建虚拟免疫响应标签）
数据源：[GDC TCGA Portal](https://portal.gdc.cancer.gov/)
包含三类配套组学数据：
1. DNA甲基化450K β值矩阵；
2. MAF体细胞突变文件（用于TMB计算）；
3. HTSeq基因表达矩阵（用于GEP、TGF-β通路评分）；
筛选后有效样本：234例原发性肿瘤，虚拟响应(R)106例、虚拟无响应(NR)128例。

### 1.3 芯片注释数据库
- IlluminaHumanMethylation450kanno.ilmn12.hg19
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19 / ilm10b2.hg19

### 1.4 公共通路注释库
GO、KEGG、Reactome、DO疾病本体、UCSC基因组、ENCODE组蛋白修饰数据库

## 2 系统软硬件环境
### 操作系统
macOS / Linux（推荐，支持并行大规模矩阵运算）；Windows兼容但大矩阵步骤内存压力更高。
### 基础配置要求
- R版本 ≥ 4.2.0
- 内存 ≥ 32GB（BMIQ归一化、ComBat批次校正、LASSO、XGBoost/LightGBM建模内存占用极高）
- CPU：多核心处理器，脚本内置`parallel`/`doParallel`并行加速
- 磁盘空间：≥50GB（原始IDAT/甲基化矩阵+中间结果+绘图文件）
### 工作目录
代码内置`setwd()`固定路径，运行前必须修改为本地数据存储路径，避免文件读取失败。

## 3 依赖R包完整清单&安装脚本
### 3.1 CRAN常规包一键安装
```r
install.packages(c(
  "dplyr", "data.table", "stringr", "ggplot2", "ggpubr", "tidyverse",
  "caret", "impute", "FactoMineR", "car", "multtest", "parallel",
  "foreach", "doParallel", "limma", "tibble", "reshape2", "glmnet",
  "randomForest", "xgboost", "Boruta", "lightgbm", "pROC", "knitr",
  "kableExtra", "tidyr", "fmsb", "rmda", "gridExtra", "boot", "RColorBrewer",
  "mlr3", "mlr3learners", "mlr3tuning", "mlr3extralearners", "kknn",
  "patchwork", "stats", "rms", "dcurves", "scales", "grid", "PMCMRplus",
  "kernelshap", "shapviz", "pheatmap"
))
```

### 3.2 Bioconductor生信专用包安装
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "TCGAbiolinks", "GenomicRanges", "biomaRt", "SummarizedExperiment",
  "ChAMP", "sva", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "AnnotationDbi",
  "clusterProfiler", "org.Hs.eg.db", "pathview", "minfi", "enrichplot",
  "ReactomePA", "DOSE"
))
```

### 3.3 编译类包额外依赖说明
`xgboost`、`lightgbm`、`ranger`需系统预装C++编译环境：
- Linux：`gcc g++ libomp-dev`
- macOS：Xcode命令行工具 + OpenMP
- Windows：Rtools42

## 4 完整标准化分析流程
本项目严格遵循论文技术路线，共11个连续分析阶段：
1. **多源数据下载与整合**
   批量读取GEO IDAT/甲基化矩阵、TCGA多组学数据，过滤缺失、异常、非原发肿瘤样本；
2. **三大免疫标志物定量计算**
   TMB突变负荷、GEP免疫效应评分、TGF-β抑制通路评分；
3. **TCGA虚拟免疫响应标签构建**
   30%/70%分位数分层，`TMB高+GEP高+TGF-β低=响应(R)`，反之无响应(NR)；
4. **甲基化数据全局预处理**
   探针/样本缺失过滤 → KNN插补 → BMIQ信号归一化 → ComBat批次校正消除平台偏差；
5. **差异甲基化位点(DMP)筛选**
   CV低变异过滤 → logit转换M值 → limma经验贝叶斯差异分析(FDR<0.05, |Δβ|>0.2) → 高相关冗余探针去重；
6. **多层级特征筛选**
   LASSO回归初筛 → RF/XGBoost/ElasticNet/RF-RFE/Boruta/LightGBM 6种算法并行筛选；
7. **特征筛选方法横向评估**
   统一逻辑回归模型，验证集AUC/F1/DCA对比，确定RF-RFE 32个CpG为最优特征集；
8. **8类机器学习模型训练与超参优化**
   LASSO、ElasticNet、随机森林、XGBoost、LightGBM、SVM、KNN、神经网络；5折网格搜索超参数；
9. **多维度模型性能评估**
   独立外部测试集计算AUC、准确率、灵敏度、特异度、MCC；绘制ROC/校准曲线/DCA/混淆矩阵；DeLong、Friedman统计学检验；
10. **Kernel SHAP模型可解释分析**
    全局特征重要性、SHAP蜂群图/力图/瀑布图/依赖图，锁定cg09115473、cg16313807核心位点；
11. **靶基因注释与多数据库功能富集**
    CpG-基因映射 → GO/KEGG/Reactome/DO富集 → 基因组分布分析 → HLA免疫逃逸机制绘图。

## 5 代码文件说明
```
./
├── 毕业论文代码.R           # 唯一完整主分析脚本，按论文章节分段，带标准化中文注释
├── data/                    # 原始输入数据目录（需自行下载GEO/TCGA放入）
├── output/                  # 全部分析输出结果（运行自动生成子文件夹）
│   ├── 01_data_preprocess/ # 预处理后甲基化矩阵、样本分组表
│   ├── 02_signature_label/  # TMB/GEP/TGF-β评分、TCGA虚拟标签文件
│   ├── 03_diff_methyl/     # limma差异分析表、火山图
│   ├── 04_feature_select/  # 7种算法筛选CpG列表、特征重要性图表
│   ├── 05_ml_model/        # 8类机器学习模型权重、预测结果、评估指标
│   ├── 06_shap_analysis/   # SHAP数值文件、全套SHAP可视化图表
│   ├── 07_enrichment/      # GO/KEGG/Reactome/DO富集表格与气泡图/条形图
│   └── figures_paper/      # 论文全部标准绘图（技术路线、ROC、热图、机制图等）
├── LICENSE                  # MIT开源协议
└── README.md                # 本项目说明文档
```

## 6 输入文件规范
所有文件命名、存储路径与代码内`readRDS()`/`read.csv()`严格匹配，缺失会直接报错：
1. GEO原始数据：各GSE ID IDAT文件/β值矩阵csv；
2. TCGA多组学：
   - 甲基化450K RDS矩阵；
   - MAF突变注释文件；
   - HTSeq基因表达count矩阵；
   - TCGA临床样本信息表；
3. 芯片注释文件：450K/850K官方探针注释R包（代码自动加载，无需手动下载）；
4. 样本分组标签：GEO真实R/NR标签、TCGA虚拟标签（脚本自动生成，无需提前准备）。

## 7 输出结果目录结构
运行脚本后自动递归创建output下分层文件夹，每一步结果独立归档，无文件覆盖冲突：
- 数值结果：csv/RDS格式，包含全部指标、CpG列表、富集表格，可直接导入GraphPad/AI绘图；
- 可视化文件：PDF/PNG双格式，高分辨率300DPI，适配毕业论文插图要求；
- 统计日志：每阶段输出日志txt，记录样本数量、探针数量、AUC等关键指标，便于核对论文数值。

## 8 复现完整运行步骤
1. **数据准备**
   从GEO数据库下载9套850K数据集、GDC下载TCGA泛癌450K/突变/表达数据，放入`./data/`目录，保持原文件名；
2. **环境配置**
   打开R/RStudio，依次执行[3 依赖R包](#3-依赖r包完整清单安装脚本)两段安装代码，确认无报错；
3. **路径修改**
   打开`毕业论文代码.R`，修改首行`setwd("你的本地项目根目录")`；
4. **分段运行脚本**
   按代码内置章节顺序依次执行，大矩阵预处理、建模步骤建议单独分段运行，中途保存R会话；
5. **结果校验**
   全部代码运行完成后，核对`output/05_ml_model/test_performance.csv`中随机森林测试集AUC≈0.8295，与论文结果一致即复现成功。

## 9 模型核心指标结果摘要
### 最优特征筛选：RF-RFE
32个核心CpG位点，验证集AUC=0.8804，综合指标优于LASSO/XGBoost等6种方法。
### 最优预测模型：随机森林（独立外部测试集）
- AUC = 0.8295
- 准确率 = 0.7931
- 灵敏度 = 0.5500
- 特异度 = 0.9211
- MCC = 0.5232
### 核心表观机制
cg09115473(TRIO增强子区)、cg16313807为TOP2预测CpG；HLA-I/II家族基因高甲基化沉默抗原呈递通路，是肿瘤免疫逃逸关键表观驱动。

## 10 项目限制与注意事项
1. **数据获取限制**
   GEO/TCGA原始测序数据受公共数据库访问规范约束，本仓库仅提供分析代码，不存储原始测序大文件，需研究者自行下载；
2. **环境一致性**
   R版本、R包版本差异会造成AUC微小浮动，建议使用renv锁定包版本实现100%完全复现；
3. **内存限制**
   32GB以下内存会导致ComBat、SHAP分析崩溃，建议使用服务器或高性能工作站运行；
4. **样本分区规范**
   严格预先拆分外部测试集，全程无测试集数据泄露，禁止用测试集调参；
5. **虚拟标签局限性**
   TCGA无真实免疫治疗临床数据，TMB/GEP/TGF-β构建虚拟标签仅用于模型训练，临床验证仍需前瞻性真实队列。

## 11 后续拓展研究方向
1. 整合RNA-seq、蛋白质组、单细胞多组学构建多模态融合预测模型；
2. 前瞻性多中心临床队列验证32个CpG甲基化特征的临床预测效能；
3. 体外细胞实验验证cg09115473、TRIO基因调控T细胞活化的表观机制；
4. 联合去甲基化药物与ICIs协同治疗的体外/动物模型验证；
5. 开发网页在线预测工具，输入患者甲基化矩阵自动输出免疫治疗响应风险。

## 12 论文引用规范
### 本科毕业论文引用格式
王嘉然. 基于DNA甲基化异常的癌症免疫治疗反应预测[D]. 哈尔滨: 哈尔滨工业大学, 2026.

### GitHub项目引用
若本代码用于科研产出，请同时引用本仓库与毕业论文：
```
Wang J. Code for Prediction of Cancer Immunotherapy Response Based on Abnormal DNA Methylation[EB/OL]. GitHub, 2026. https://github.com/Jia-RanWang/DNA-Methylation-Immunotherapy-Prediction
```

## 13 开源许可
本项目采用 **MIT License**，可自由商用/学术使用、修改分发，但请保留原作者与论文引用声明。详见仓库内`LICENSE`文件。
