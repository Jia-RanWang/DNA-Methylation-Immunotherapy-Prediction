options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_ROOT <- Sys.getenv(
  "THESIS_ROOT",
  "/Users/mac/Desktop/TCGA_Thesis"
)

TCGA_INPUT_DIR <- Sys.getenv(
  "TCGA_INPUT_DIR",
  "/Users/mac/Desktop/TCGA_DATA"
)

GEO_INPUT_DIR <- Sys.getenv(
  "GEO_INPUT_DIR",
  "/Users/mac/Desktop/GEO_DATA"
)

RESULT_DIR <- file.path(PROJECT_ROOT, "results")
TABLE_DIR <- file.path(RESULT_DIR, "tables")
RDATA_DIR <- file.path(RESULT_DIR, "rdata")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDATA_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "S4Vectors",
  "ChAMP",
  "dplyr",
  "tidyr",
  "stringr",
  "data.table",
  "caret",
  "impute",
  "sva",
  "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包：",
    paste(missing_packages, collapse = ", ")
  )
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(S4Vectors)
library(ChAMP)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(caret)
library(impute)
library(sva)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ==================================================================
# 一、计算TCGA TMB
# ==================================================================

mutation_dir <- file.path(TCGA_INPUT_DIR, "TCGA_Mutations_All")

read_rdata_object <- function(file, candidates = character()) {
  env <- new.env(parent = emptyenv())
  object_names <- load(file, envir = env)

  for (object_name in candidates) {
    if (exists(object_name, envir = env, inherits = FALSE)) {
      return(get(object_name, envir = env, inherits = FALSE))
    }
  }

  for (object_name in object_names) {
    object <- get(object_name, envir = env, inherits = FALSE)

    if (is.data.frame(object) &&
        "Tumor_Sample_Barcode" %in% colnames(object)) {
      return(object)
    }
  }

  NULL
}

mutation_files <- list.files(
  mutation_dir,
  pattern = "_mutations\\.rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

non_synonymous_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Ins",
  "Frame_Shift_Del",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins"
)

tmb_result_list <- list()

for (file in mutation_files) {
  cancer_type <- sub(
    "TCGA-([A-Z]+)_.*",
    "\\1",
    basename(file)
  )

  message("计算TMB：", cancer_type)

  mut_data <- read_rdata_object(
    file,
    candidates = c("mut_data", "brca_mut_data", "data")
  )

  if (is.null(mut_data)) {
    warning("未找到可用突变对象：", file)
    next
  }

  required_cols <- c(
    "Tumor_Sample_Barcode",
    "Variant_Classification"
  )

  if (!all(required_cols %in% colnames(mut_data))) {
    warning("突变文件缺少必要列：", file)
    next
  }

  sample_type_code <- substr(
    mut_data$Tumor_Sample_Barcode,
    14,
    15
  )

  tumor_data <- mut_data[
    sample_type_code %in% sprintf("%02d", 1:9),
    ,
    drop = FALSE
  ]

  if (nrow(tumor_data) == 0) {
    next
  }

  tumor_non_synonymous <- tumor_data[
    tumor_data$Variant_Classification %in%
      non_synonymous_classes,
    ,
    drop = FALSE
  ]

  if (nrow(tumor_non_synonymous) == 0) {
    next
  }

  mutation_counts <- table(
    tumor_non_synonymous$Tumor_Sample_Barcode
  )

  tmb_result_list[[cancer_type]] <- data.frame(
    Sample = names(mutation_counts),
    TMB = as.numeric(mutation_counts) / 38,
    CancerType = cancer_type,
    SampleCount = as.numeric(mutation_counts),
    stringsAsFactors = FALSE
  )
}

if (length(tmb_result_list) == 0) {
  stop("没有计算得到TMB结果")
}

tmb_results <- do.call(rbind, tmb_result_list)

write.csv(
  tmb_results,
  file.path(TABLE_DIR, "TCGA_All_Cancer_TMB_Results.csv"),
  row.names = FALSE
)

save(
  tmb_results,
  file = file.path(RDATA_DIR, "TCGA_All_Cancer_TMB_Results.RData")
)

# ==================================================================
# 二、计算GEP和TGF-β评分
# ==================================================================

expression_dir <- file.path(TCGA_INPUT_DIR, "TCGA_Expression_All")

immune_genes <- c(
  "CD3D", "CD3E", "CD8A", "CD8B",
  "GZMA", "GZMB", "PRF1", "IFNG",
  "TBX21", "STAT1", "CXCL9", "CXCL10"
)

tgfb_genes <- c(
  "TGFB1", "TGFB2", "TGFB3",
  "TGFBR1", "TGFBR2",
  "SMAD2", "SMAD3", "SMAD4", "SMAD7"
)

clean_gene_names <- function(x) {
  x <- sub("\\..*$", "", x)
  toupper(x)
}

extract_expression_matrix <- function(object) {
  if (inherits(object, "SummarizedExperiment")) {
    assay_names <- SummarizedExperiment::assayNames(object)

    if ("unstranded" %in% assay_names) {
      expression_matrix <- SummarizedExperiment::assay(
        object,
        "unstranded"
      )
    } else if ("counts" %in% assay_names) {
      expression_matrix <- SummarizedExperiment::assay(
        object,
        "counts"
      )
    } else {
      expression_matrix <- SummarizedExperiment::assay(
        object,
        1
      )
    }

    row_data <- as.data.frame(
      SummarizedExperiment::rowData(object)
    )

    gene_symbols <- NULL

    for (candidate in c(
      "gene_name",
      "gene_symbol",
      "external_gene_name"
    )) {
      if (candidate %in% colnames(row_data)) {
        gene_symbols <- as.character(row_data[[candidate]])
        break
      }
    }

    if (is.null(gene_symbols)) {
      gene_symbols <- rownames(expression_matrix)
    }

    return(list(
      matrix = as.matrix(expression_matrix),
      symbols = clean_gene_names(gene_symbols)
    ))
  }

  if (is.matrix(object) || is.data.frame(object)) {
    object <- as.data.frame(object)

    if ("unstranded" %in% colnames(object)) {
      expression_matrix <- as.matrix(
        object[, setdiff(
          colnames(object),
          c("gene_id", "gene_name", "gene_type")
        )]
      )

      gene_symbols <- if ("gene_name" %in% colnames(object)) {
        object$gene_name
      } else {
        object$gene_id
      }

      rownames(expression_matrix) <- object$gene_id

      return(list(
        matrix = expression_matrix,
        symbols = clean_gene_names(gene_symbols)
      ))
    }

    return(list(
      matrix = as.matrix(object),
      symbols = clean_gene_names(rownames(object))
    ))
  }

  stop("不支持的表达数据对象类型")
}

convert_ensembl_to_symbol <- function(ensembl_ids) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    return(clean_gene_names(ensembl_ids))
  }

  ensembl_ids <- sub("\\..*$", "", ensembl_ids)
  unique_ids <- unique(ensembl_ids)

  mapping <- tryCatch({
    mart <- biomaRt::useMart(
      "ensembl",
      dataset = "hsapiens_gene_ensembl"
    )

    biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol"
      ),
      filters = "ensembl_gene_id",
      values = unique_ids,
      mart = mart
    )
  }, error = function(e) {
    NULL
  })

  if (is.null(mapping) || nrow(mapping) == 0) {
    return(clean_gene_names(ensembl_ids))
  }

  mapping <- mapping[
    !duplicated(mapping$ensembl_gene_id),
    ,
    drop = FALSE
  ]

  result <- mapping$hgnc_symbol[
    match(ensembl_ids, mapping$ensembl_gene_id)
  ]

  result[is.na(result) | result == ""] <- "UNKNOWN"
  clean_gene_names(result)
}

expression_files <- list.files(
  expression_dir,
  pattern = "_expression\\.rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

expression_result_list <- list()

for (file in expression_files) {
  cancer_type <- sub(
    "TCGA-([A-Z]+)_.*",
    "\\1",
    basename(file)
  )

  message("计算表达指标：", cancer_type)

  exp_data <- read_rdata_object(
    file,
    candidates = c(
      "exp_data",
      "brca_expr_data",
      "laml_se",
      "data"
    )
  )

  if (is.null(exp_data)) {
    warning("未找到表达对象：", file)
    next
  }

  expression_info <- tryCatch(
    extract_expression_matrix(exp_data),
    error = function(e) NULL
  )

  if (is.null(expression_info)) {
    next
  }

  expression_matrix <- expression_info$matrix
  gene_symbols <- expression_info$symbols

  if (length(gene_symbols) != nrow(expression_matrix)) {
    gene_symbols <- convert_ensembl_to_symbol(
      rownames(expression_matrix)
    )
  }

  sample_names <- colnames(expression_matrix)

  sample_type_code <- substr(sample_names, 14, 15)
  tumor_index <- sample_type_code %in% sprintf("%02d", 1:9)

  if (!any(tumor_index)) {
    next
  }

  expression_matrix <- expression_matrix[
    ,
    tumor_index,
    drop = FALSE
  ]

  sample_names <- sample_names[tumor_index]

  available_immune <- intersect(
    immune_genes,
    gene_symbols
  )

  available_tgfb <- intersect(
    tgfb_genes,
    gene_symbols
  )

  if (length(available_immune) >= 3) {
    immune_index <- which(
      gene_symbols %in% available_immune
    )

    gep_scores <- colMeans(
      expression_matrix[immune_index, , drop = FALSE],
      na.rm = TRUE
    )
  } else {
    gep_scores <- rep(NA_real_, length(sample_names))
  }

  if (length(available_tgfb) >= 3) {
    tgfb_index <- which(
      gene_symbols %in% available_tgfb
    )

    tgfb_scores <- colMeans(
      expression_matrix[tgfb_index, , drop = FALSE],
      na.rm = TRUE
    )
  } else {
    tgfb_scores <- rep(NA_real_, length(sample_names))
  }

  expression_result_list[[cancer_type]] <- data.frame(
    Sample = sample_names,
    GEP = as.numeric(gep_scores),
    TGFB = as.numeric(tgfb_scores),
    CancerType = cancer_type,
    GEP_Gene_Count = length(available_immune),
    TGFB_Gene_Count = length(available_tgfb),
    stringsAsFactors = FALSE
  )
}

if (length(expression_result_list) == 0) {
  stop("没有得到GEP/TGF-β结果")
}

expression_results <- do.call(
  rbind,
  expression_result_list
)

write.csv(
  expression_results,
  file.path(TABLE_DIR, "TCGA_All_Cancer_GEP_TGFB_Results.csv"),
  row.names = FALSE
)

save(
  expression_results,
  file = file.path(
    RDATA_DIR,
    "TCGA_All_Cancer_GEP_TGFB_Results.RData"
  )
)

# ==================================================================
# 三、构建虚拟响应标签
# ==================================================================

tmb_df <- read.csv(
  file.path(TABLE_DIR, "TCGA_All_Cancer_TMB_Results.csv"),
  stringsAsFactors = FALSE
)

expr_df <- read.csv(
  file.path(
    TABLE_DIR,
    "TCGA_All_Cancer_GEP_TGFB_Results.csv"
  ),
  stringsAsFactors = FALSE
)

tmb_df$Sample_Core <- substr(tmb_df$Sample, 1, 15)
expr_df$Sample_Core <- substr(expr_df$Sample, 1, 15)

combined_df <- merge(
  tmb_df[, c("Sample_Core", "TMB", "CancerType")],
  expr_df[, c("Sample_Core", "GEP", "TGFB")],
  by = "Sample_Core",
  all = FALSE
)

combined_df <- combined_df[
  !duplicated(combined_df$Sample_Core),
  ,
  drop = FALSE
]

tmb_30 <- quantile(
  combined_df$TMB,
  0.3,
  na.rm = TRUE
)

tmb_70 <- quantile(
  combined_df$TMB,
  0.7,
  na.rm = TRUE
)

gep_30 <- quantile(
  combined_df$GEP,
  0.3,
  na.rm = TRUE
)

gep_70 <- quantile(
  combined_df$GEP,
  0.7,
  na.rm = TRUE
)

tgfb_30 <- quantile(
  combined_df$TGFB,
  0.3,
  na.rm = TRUE
)

tgfb_70 <- quantile(
  combined_df$TGFB,
  0.7,
  na.rm = TRUE
)

threshold_table <- data.frame(
  Marker = c("TMB", "GEP", "TGFB"),
  Low_Threshold = c(
    tmb_30,
    gep_30,
    tgfb_30
  ),
  High_Threshold = c(
    tmb_70,
    gep_70,
    tgfb_70
  )
)

write.csv(
  threshold_table,
  file.path(TABLE_DIR, "Virtual_Label_Thresholds.csv"),
  row.names = FALSE
)

combined_df <- combined_df %>%
  mutate(
    TMB_Status = case_when(
      TMB >= tmb_70 ~ "High",
      TMB <= tmb_30 ~ "Low",
      TRUE ~ "Middle"
    ),
    GEP_Status = case_when(
      GEP >= gep_70 ~ "High",
      GEP <= gep_30 ~ "Low",
      TRUE ~ "Middle"
    ),
    TGFB_Status = case_when(
      TGFB <= tgfb_30 ~ "Low",
      TGFB >= tgfb_70 ~ "High",
      TRUE ~ "Middle"
    ),
    Response_Label = case_when(
      TMB_Status == "High" &
        GEP_Status == "High" &
        TGFB_Status == "Low" ~ "R",
      TMB_Status == "Low" &
        GEP_Status == "Low" &
        TGFB_Status == "High" ~ "NR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Response_Label))

write.csv(
  combined_df[
    ,
    c("Sample_Core", "CancerType", "Response_Label")
  ],
  file.path(TABLE_DIR, "TCGA_Virtual_Response_Labels.csv"),
  row.names = FALSE
)

save(
  combined_df,
  file = file.path(RDATA_DIR, "TCGA_Virtual_Response_Labels.RData")
)

# ==================================================================
# 四、读取GEO和TCGA甲基化原始β矩阵
# ==================================================================

geo_ids <- c(
  "GSE126043",
  "GSE264158",
  "GSE175699",
  "GSE181781",
  "GSE172468",
  "GSE235122",
  "GSE305240"
)

geo_raw_list <- list()
geo_pd_list <- list()

for (gse_id in geo_ids) {
  gse_dir <- file.path(
    GEO_INPUT_DIR,
    paste0(gse_id, "_RAW")
  )

  if (!dir.exists(gse_dir)) {
    warning("GEO目录不存在：", gse_dir)
    next
  }

  message("读取GEO：", gse_id)

  geo_object <- ChAMP::champ.load(
    directory = gse_dir,
    arraytype = "EPIC",
    filterBeads = FALSE,
    filterDetP = FALSE
  )

  geo_raw_list[[gse_id]] <- geo_object$beta
  geo_pd_list[[gse_id]] <- geo_object$pd %>%
    mutate(GSE_ID = gse_id)
}

label_df <- read.csv(
  file.path(TABLE_DIR, "TCGA_Virtual_Response_Labels.csv"),
  stringsAsFactors = FALSE
)

tcga_beta_files <- list.files(
  TCGA_INPUT_DIR,
  pattern = "TCGA-.*_methy_beta.*\\.rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

extract_beta_matrix <- function(object) {
  if (inherits(object, "SummarizedExperiment")) {
    return(as.matrix(SummarizedExperiment::assay(object)))
  }

  if (is.matrix(object)) {
    return(object)
  }

  if (is.data.frame(object)) {
    numeric_columns <- vapply(
      object,
      is.numeric,
      logical(1)
    )

    matrix_data <- as.matrix(
      object[, numeric_columns, drop = FALSE]
    )

    if (!is.null(rownames(object))) {
      rownames(matrix_data) <- rownames(object)
    }

    return(matrix_data)
  }

  stop("不能提取β矩阵")
}

tcga_raw_beta <- list()
tcga_raw_pd <- list()

for (file in tcga_beta_files) {
  cancer_type <- stringr::str_match(
    basename(file),
    "TCGA-(.*?)_methy_beta"
  )[, 2]

  env <- new.env(parent = emptyenv())
  object_names <- load(file, envir = env)

  beta_object <- NULL

  for (object_name in c("beta_expr", "data", "beta")) {
    if (exists(object_name, envir = env, inherits = FALSE)) {
      beta_object <- get(
        object_name,
        envir = env,
        inherits = FALSE
      )
      break
    }
  }

  if (is.null(beta_object)) {
    next
  }

  beta_matrix <- tryCatch(
    extract_beta_matrix(beta_object),
    error = function(e) NULL
  )

  if (is.null(beta_matrix)) {
    next
  }

  sample_names <- colnames(beta_matrix)
  sample_core <- substr(sample_names, 1, 15)

  match_index <- match(
    sample_core,
    label_df$Sample_Core
  )

  sample_group <- label_df$Response_Label[match_index]

  pd <- data.frame(
    Sample_Name = sample_names,
    Sample_Group = sample_group,
    GSE_ID = "TCGA",
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(Sample_Group))

  if (nrow(pd) == 0) {
    next
  }

  pd <- pd %>%
    arrange(match(
      Sample_Name,
      colnames(beta_matrix)
    ))

  tcga_raw_beta[[cancer_type]] <- beta_matrix[
    ,
    pd$Sample_Name,
    drop = FALSE
  ]

  tcga_raw_pd[[cancer_type]] <- pd
}

save(
  geo_raw_list,
  file = file.path(RDATA_DIR, "geo_raw_list.RData")
)

save(
  geo_pd_list,
  file = file.path(RDATA_DIR, "geo_pd_list.RData")
)

save(
  tcga_raw_beta,
  file = file.path(RDATA_DIR, "tcga_raw_beta.RData")
)

save(
  tcga_raw_pd,
  file = file.path(RDATA_DIR, "tcga_raw_pd.RData")
)

# ==================================================================
# 五、合并数据并划分训练集、验证集
# ==================================================================

all_geo_pd <- bind_rows(
  lapply(geo_pd_list, function(x) {
    x[, c("Sample_Name", "Sample_Group", "GSE_ID")]
  })
)

all_tcga_pd <- bind_rows(
  lapply(tcga_raw_pd, function(x) {
    x[, c("Sample_Name", "Sample_Group", "GSE_ID")]
  })
)

all_pd <- bind_rows(
  all_geo_pd,
  all_tcga_pd
)

extract_beta_by_samples <- function(pd_df, geo_list, tcga_list) {
  beta_list <- vector("list", nrow(pd_df))
  names(beta_list) <- pd_df$Sample_Name

  for (i in seq_len(nrow(pd_df))) {
    sample_name <- pd_df$Sample_Name[i]
    gse_id <- pd_df$GSE_ID[i]

    if (gse_id == "TCGA") {
      for (cancer_type in names(tcga_list)) {
        matrix_data <- tcga_list[[cancer_type]]

        if (sample_name %in% colnames(matrix_data)) {
          beta_list[[i]] <- matrix_data[
            ,
            sample_name,
            drop = FALSE
          ]
          break
        }
      }
    } else if (gse_id %in% names(geo_list)) {
      matrix_data <- geo_list[[gse_id]]

      if (sample_name %in% colnames(matrix_data)) {
        beta_list[[i]] <- matrix_data[
          ,
          sample_name,
          drop = FALSE
        ]
      }
    }
  }

  beta_list <- beta_list[
    !vapply(beta_list, is.null, logical(1))
  ]

  if (length(beta_list) == 0) {
    stop("没有找到可合并的β矩阵")
  }

  common_probes <- Reduce(
    intersect,
    lapply(beta_list, rownames)
  )

  beta_combined <- do.call(
    cbind,
    lapply(beta_list, function(x) {
      x[common_probes, , drop = FALSE]
    })
  )

  beta_combined
}

all_beta_raw <- extract_beta_by_samples(
  all_pd,
  geo_raw_list,
  tcga_raw_beta
)

all_pd <- all_pd[
  match(colnames(all_beta_raw), all_pd$Sample_Name),
  ,
  drop = FALSE
]

save(
  all_beta_raw,
  all_pd,
  file = file.path(RDATA_DIR, "all_combined_raw.RData")
)

train_idx <- createDataPartition(
  y = paste(all_pd$Sample_Group),
  p = 0.7,
  list = FALSE
)

train_beta_raw <- all_beta_raw[, train_idx, drop = FALSE]
train_pd_raw <- all_pd[train_idx, , drop = FALSE]

val_beta_raw <- all_beta_raw[, -train_idx, drop = FALSE]
val_pd_raw <- all_pd[-train_idx, , drop = FALSE]

save(
  train_beta_raw,
  train_pd_raw,
  file = file.path(RDATA_DIR, "train_beta_raw.RData")
)

save(
  val_beta_raw,
  val_pd_raw,
  file = file.path(RDATA_DIR, "val_beta_raw.RData")
)

# ==================================================================
# 六、训练集和验证集预处理
# ==================================================================

epic_annotation <- getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b2.hg19
)

hm450_annotation <- getAnnotation(
  IlluminaHumanMethylation450kanno.ilmn12.hg19
)

ann_epic <- as.data.frame(epic_annotation) %>%
  tibble::rownames_to_column("Name") %>%
  mutate(Array = "EPIC")

ann_450k <- as.data.frame(hm450_annotation) %>%
  tibble::rownames_to_column("Name") %>%
  mutate(Array = "450K")

ann_combined <- bind_rows(
  ann_epic,
  ann_450k
) %>%
  distinct(Name, .keep_all = TRUE)

ann_filt <- ann_combined %>%
  filter(
    !str_detect(Name, "^ch\\."),
    !str_detect(chr, "X|Y"),
    is.na(Probe_rs) | Probe_rs == ""
  )

probe_na_rate_train <- rowMeans(is.na(train_beta_raw))

probes_to_keep_train_step1 <- names(
  probe_na_rate_train[
    probe_na_rate_train <= 0.05
  ]
)

train_beta_filtered_probe <- train_beta_raw[
  probes_to_keep_train_step1,
  ,
  drop = FALSE
]

sample_na_rate_train <- colMeans(
  is.na(train_beta_filtered_probe)
)

samples_to_keep_train <- names(
  sample_na_rate_train[
    sample_na_rate_train <= 0.1
  ]
)

train_beta_filtered_sample <- train_beta_filtered_probe[
  ,
  samples_to_keep_train,
  drop = FALSE
]

train_pd_filtered <- train_pd_raw %>%
  filter(Sample_Name %in% samples_to_keep_train) %>%
  arrange(match(
    Sample_Name,
    colnames(train_beta_filtered_sample)
  ))

probes_to_keep_train <- intersect(
  rownames(train_beta_filtered_sample),
  ann_filt$Name
)

train_beta_filtered_final <- train_beta_filtered_sample[
  probes_to_keep_train,
  ,
  drop = FALSE
]

train_beta_imputed <- t(
  impute.knn(
    t(train_beta_filtered_final),
    k = 10
  )$data
)

train_beta_norm <- champ.norm(
  beta = train_beta_imputed,
  arraytype = "EPIC",
  method = "BMIQ",
  cores = 2
)

train_beta_norm <- as.matrix(train_beta_norm)

batch_train <- train_pd_filtered$GSE_ID
mod_matrix_train <- model.matrix(
  ~ Sample_Group,
  data = train_pd_filtered
)

train_beta_combat <- ComBat(
  dat = train_beta_norm,
  batch = batch_train,
  mod = mod_matrix_train,
  par.prior = TRUE,
  prior.plots = FALSE
)

train_beta_combat_df <- as.data.frame(
  t(train_beta_combat)
) %>%
  tibble::rownames_to_column("Sample_Name")

train_data_final <- train_pd_filtered %>%
  inner_join(
    train_beta_combat_df,
    by = "Sample_Name"
  )

save(
  train_data_final,
  file = file.path(RDATA_DIR, "train_data.RData")
)

# -------------------- 验证集 --------------------

val_beta_filtered_probe <- val_beta_raw[
  probes_to_keep_train_step1,
  ,
  drop = FALSE
]

sample_na_rate_val <- colMeans(
  is.na(val_beta_filtered_probe)
)

samples_to_keep_val <- names(
  sample_na_rate_val[
    sample_na_rate_val <= 0.1
  ]
)

val_beta_filtered_sample <- val_beta_filtered_probe[
  ,
  samples_to_keep_val,
  drop = FALSE
]

val_pd_filtered <- val_pd_raw %>%
  filter(Sample_Name %in% samples_to_keep_val) %>%
  arrange(match(
    Sample_Name,
    colnames(val_beta_filtered_sample)
  ))

val_beta_filtered_final <- val_beta_filtered_sample[
  probes_to_keep_train,
  ,
  drop = FALSE
]

val_beta_imputed <- t(
  impute.knn(
    t(val_beta_filtered_final),
    k = 10
  )$data
)

val_beta_norm <- champ.norm(
  beta = val_beta_imputed,
  arraytype = "EPIC",
  method = "BMIQ",
  cores = 2
)

val_beta_norm <- as.matrix(val_beta_norm)

batch_val <- val_pd_filtered$GSE_ID
mod_matrix_val <- model.matrix(
  ~ Sample_Group,
  data = val_pd_filtered
)

val_beta_combat <- ComBat(
  dat = val_beta_norm,
  batch = batch_val,
  mod = mod_matrix_val,
  par.prior = TRUE,
  prior.plots = FALSE
)

val_beta_combat_df <- as.data.frame(
  t(val_beta_combat)
) %>%
  tibble::rownames_to_column("Sample_Name")

val_data_final <- val_pd_filtered %>%
  inner_join(
    val_beta_combat_df,
    by = "Sample_Name"
  )

save(
  val_data_final,
  file = file.path(RDATA_DIR, "val_data.RData")
)

# ==================================================================
# 七、外部测试集预处理
# ==================================================================

test_geo_ids <- c(
  "GSE119144",
  "GSE277573"
)

for (gse_id in test_geo_ids) {
  test_dir <- file.path(
    GEO_INPUT_DIR,
    paste0(gse_id, "_RAW")
  )

  if (!dir.exists(test_dir)) {
    warning("测试集目录不存在：", test_dir)
    next
  }

  message("处理外部测试集：", gse_id)

  test_object <- ChAMP::champ.load(
    directory = test_dir,
    arraytype = "EPIC"
  )

  test_beta_norm <- champ.norm(
    beta = test_object$beta,
    arraytype = "EPIC",
    method = "BMIQ",
    cores = 8
  )

  norm_df <- as.data.frame(
    t(test_beta_norm)
  ) %>%
    tibble::rownames_to_column("Sample_Name")

  final_df <- test_object$pd %>%
    select(Sample_Name, Sample_Group) %>%
    inner_join(norm_df, by = "Sample_Name") %>%
    relocate(Sample_Name, Sample_Group)

  save(
    final_df,
    file = file.path(
      RDATA_DIR,
      paste0(gse_id, "_test_data.RData")
    )
  )

  write.csv(
    final_df,
    file.path(
      TABLE_DIR,
      paste0(gse_id, "_test_data.csv")
    ),
    row.names = FALSE
  )
}

save(
  probes_to_keep_train_step1,
  probes_to_keep_train,
  samples_to_keep_train,
  samples_to_keep_val,
  file = file.path(RDATA_DIR, "preprocess_rules.RData")
)

message("02_preprocess_and_labels.R运行完成。")