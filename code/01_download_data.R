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

dir.create(PROJECT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "results", "rdata"), recursive = TRUE, showWarnings = FALSE)

TCGA_MUTATION_DIR <- file.path(TCGA_INPUT_DIR, "TCGA_Mutations_All")
TCGA_EXPRESSION_DIR <- file.path(TCGA_INPUT_DIR, "TCGA_Expression_All")

dir.create(TCGA_MUTATION_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TCGA_EXPRESSION_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "S4Vectors", "data.table")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包，请先安装：",
    paste(missing_packages, collapse = ", ")
  )
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(S4Vectors)
library(data.table)

cancer_types <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC",
  "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
  "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD",
  "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD",
  "TCGA-PCPG", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM",
  "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

query_has_data <- function(query) {
  if (is.null(query$results)) {
    return(FALSE)
  }

  result_tables <- query$results[
    vapply(query$results, is.data.frame, logical(1))
  ]

  if (length(result_tables) == 0) {
    return(FALSE)
  }

  any(vapply(result_tables, nrow, integer(1)) > 0)
}

download_mutation_data <- function(project_id) {
  output_file <- file.path(
    TCGA_MUTATION_DIR,
    paste0(project_id, "_mutations.rdata")
  )

  if (file.exists(output_file)) {
    message("已存在，跳过突变数据：", project_id)
    return(data.frame(
      Project = project_id,
      Type = "Mutation",
      Status = "Existing",
      Message = "",
      stringsAsFactors = FALSE
    ))
  }

  result <- tryCatch({
    options(timeout = 1800)

    query <- GDCquery(
      project = project_id,
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      access = "open"
    )

    if (!query_has_data(query)) {
      return(data.frame(
        Project = project_id,
        Type = "Mutation",
        Status = "NoData",
        Message = "",
        stringsAsFactors = FALSE
      ))
    }

    GDCdownload(
      query,
      directory = TCGA_MUTATION_DIR
    )

    mut_data <- GDCprepare(
      query,
      directory = TCGA_MUTATION_DIR
    )

    save(
      mut_data,
      file = output_file
    )

    data.frame(
      Project = project_id,
      Type = "Mutation",
      Status = "Downloaded",
      Message = paste0("n=", nrow(mut_data)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Project = project_id,
      Type = "Mutation",
      Status = "Failed",
      Message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  options(timeout = 60)
  result
}

download_expression_data <- function(project_id) {
  output_file <- file.path(
    TCGA_EXPRESSION_DIR,
    paste0(project_id, "_expression.rdata")
  )

  if (file.exists(output_file)) {
    message("已存在，跳过表达数据：", project_id)
    return(data.frame(
      Project = project_id,
      Type = "Expression",
      Status = "Existing",
      Message = "",
      stringsAsFactors = FALSE
    ))
  }

  workflows <- c("STAR - Counts", "HTSeq - Counts")

  for (workflow in workflows) {
    result <- tryCatch({
      options(timeout = 1800)

      query <- GDCquery(
        project = project_id,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = workflow,
        access = "open"
      )

      if (!query_has_data(query)) {
        stop("当前workflow没有可用数据")
      }

      GDCdownload(
        query,
        directory = TCGA_EXPRESSION_DIR
      )

      exp_data <- GDCprepare(
        query,
        directory = TCGA_EXPRESSION_DIR
      )

      save(
        exp_data,
        file = output_file
      )

      options(timeout = 60)

      return(data.frame(
        Project = project_id,
        Type = "Expression",
        Status = "Downloaded",
        Message = workflow,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      NULL
    })

    if (!is.null(result)) {
      return(result)
    }
  }

  options(timeout = 60)

  data.frame(
    Project = project_id,
    Type = "Expression",
    Status = "Failed",
    Message = "STAR - Counts和HTSeq - Counts均不可用",
    stringsAsFactors = FALSE
  )
}

download_results <- list()

for (project_id in cancer_types) {
  message("处理突变数据：", project_id)
  download_results[[length(download_results) + 1]] <-
    download_mutation_data(project_id)

  Sys.sleep(3)

  message("处理表达数据：", project_id)
  download_results[[length(download_results) + 1]] <-
    download_expression_data(project_id)

  Sys.sleep(3)
}

download_manifest <- do.call(rbind, download_results)

write.csv(
  download_manifest,
  file.path(
    PROJECT_ROOT,
    "results",
    "tables",
    "TCGA_Download_Manifest.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------
# LAML本地RNA-seq计数文件整合
# ------------------------------------------------------------------

LAML_RAW_DIR <- Sys.getenv(
  "LAML_RAW_DIR",
  "/Users/ranmac/Downloads/gdc_download_20251017_125311.407386"
)

laml_output <- file.path(
  TCGA_EXPRESSION_DIR,
  "TCGA-LAML_expression.rdata"
)

assemble_laml_expression <- function(input_dir, output_file) {
  if (!dir.exists(input_dir)) {
    message("LAML本地下载目录不存在，跳过：", input_dir)
    return(FALSE)
  }

  sample_folders <- list.dirs(
    input_dir,
    full.names = TRUE,
    recursive = FALSE
  )

  sample_tables <- list()

  for (folder in sample_folders) {
    count_file <- list.files(
      folder,
      pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
      full.names = TRUE
    )

    if (length(count_file) == 0) {
      next
    }

    count_table <- tryCatch(
      read.delim(
        count_file[1],
        skip = 1,
        header = TRUE,
        sep = "\t",
        check.names = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(count_table)) {
      next
    }

    if (!all(c("gene_id", "unstranded") %in% colnames(count_table))) {
      next
    }

    count_table <- count_table[, c("gene_id", "unstranded")]
    sample_id <- sub(
      "^.*(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$",
      "\\1",
      basename(folder)
    )

    colnames(count_table)[2] <- sample_id
    sample_tables[[sample_id]] <- count_table
  }

  if (length(sample_tables) == 0) {
    message("没有找到LAML计数文件")
    return(FALSE)
  }

  expr_matrix <- Reduce(
    function(x, y) merge(x, y, by = "gene_id"),
    sample_tables
  )

  annotation <- tryCatch({
    query <- GDCquery(
      project = "TCGA-LAML",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      access = "open"
    )

    GDCprepare(
      query,
      directory = TCGA_EXPRESSION_DIR,
      summarizedExperiment = FALSE
    )
  }, error = function(e) {
    NULL
  })

  if (!is.null(annotation) &&
      all(c("gene_id", "gene_name", "gene_type") %in%
          colnames(annotation))) {
    annotation <- annotation[
      !duplicated(annotation$gene_id),
      c("gene_id", "gene_name", "gene_type")
    ]

    expr_matrix <- merge(
      annotation,
      expr_matrix,
      by = "gene_id"
    )

    row_data <- expr_matrix[, c("gene_name", "gene_type")]
    assay_data <- as.matrix(
      expr_matrix[, setdiff(
        colnames(expr_matrix),
        c("gene_id", "gene_name", "gene_type")
      )]
    )

    rownames(assay_data) <- expr_matrix$gene_id
  } else {
    assay_data <- as.matrix(expr_matrix[, -1])
    rownames(assay_data) <- expr_matrix$gene_id

    row_data <- DataFrame(
      gene_name = rownames(assay_data),
      gene_type = NA_character_
    )
  }

  col_data <- DataFrame(
    sample_id = colnames(assay_data)
  )

  rownames(col_data) <- colnames(assay_data)

  laml_se <- SummarizedExperiment(
    assays = list(counts = assay_data),
    rowData = row_data,
    colData = col_data
  )

  exp_data <- laml_se

  save(
    exp_data,
    laml_se,
    file = output_file
  )

  message("LAML表达矩阵已保存：", output_file)
  TRUE
}

if (!file.exists(laml_output)) {
  assemble_laml_expression(
    input_dir = LAML_RAW_DIR,
    output_file = laml_output
  )
}

message("01_download_data.R运行完成。")