options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_ROOT <- Sys.getenv(
  "THESIS_ROOT",
  "/Users/mac/Desktop/TCGA_Thesis"
)

RESULT_DIR <- file.path(PROJECT_ROOT, "results")
TABLE_DIR <- file.path(RESULT_DIR, "tables")
FIGURE_DIR <- file.path(RESULT_DIR, "figures")
RDATA_DIR <- file.path(RESULT_DIR, "rdata")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDATA_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "dplyr",
  "limma",
  "tibble",
  "ggplot2",
  "reshape2"
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

library(dplyr)
library(limma)
library(tibble)
library(ggplot2)
library(reshape2)

load(
  file.path(RDATA_DIR, "train_data.RData")
)

train_data <- train_data_final
train_data$Sample_ID <- train_data$GSE_ID

sample_info <- train_data %>%
  select(
    Sample_Name,
    Sample_ID,
    Sample_Group
  )

probe_matrix <- train_data %>%
  select(
    -Sample_Name,
    -Sample_ID,
    -Sample_Group,
    -GSE_ID
  )

# ==================================================================
# 一、低变异探针过滤
# ==================================================================

probe_cv <- apply(
  probe_matrix,
  2,
  function(x) {
    sd(x, na.rm = TRUE) /
      mean(x, na.rm = TRUE)
  }
)

high_cv_probes <- names(
  probe_cv[probe_cv > 0.1]
)

probe_matrix_filtered <- probe_matrix[
  ,
  high_cv_probes,
  drop = FALSE
]

write.csv(
  data.frame(
    Probe_ID = high_cv_probes,
    CV = probe_cv[high_cv_probes]
  ),
  file.path(TABLE_DIR, "High_CV_Probes.csv"),
  row.names = FALSE
)

# ==================================================================
# 二、Beta值转换为M值
# ==================================================================

epsilon <- 1e-6

m_matrix <- log2(
  (as.matrix(probe_matrix_filtered) + epsilon) /
    (1 - as.matrix(probe_matrix_filtered) + epsilon)
)

group <- factor(train_data$Sample_Group)

if (length(levels(group)) != 2) {
  stop("Sample_Group必须包含两个分组")
}

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- levels(group)

contrast_formula <- paste(
  levels(group)[2],
  levels(group)[1],
  sep = "-"
)

contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)

fit <- lmFit(
  t(m_matrix),
  design
)

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(fit2)

results_limma <- topTable(
  fit2,
  number = Inf,
  sort.by = "P"
)

# ==================================================================
# 三、Delta Beta筛选
# ==================================================================

mean_group1_beta <- apply(
  probe_matrix_filtered[
    train_data$Sample_Group == levels(group)[1],
    ,
    drop = FALSE
  ],
  2,
  mean,
  na.rm = TRUE
)

mean_group2_beta <- apply(
  probe_matrix_filtered[
    train_data$Sample_Group == levels(group)[2],
    ,
    drop = FALSE
  ],
  2,
  mean,
  na.rm = TRUE
)

results_limma <- results_limma %>%
  rownames_to_column("Probe_ID") %>%
  mutate(
    Mean_Beta_Group1 = mean_group1_beta[Probe_ID],
    Mean_Beta_Group2 = mean_group2_beta[Probe_ID],
    Delta_Beta = Mean_Beta_Group2 - Mean_Beta_Group1
  )

fdr_threshold <- 0.05
delta_beta_threshold <- 0.2

sig_probes_df <- results_limma %>%
  filter(
    adj.P.Val < fdr_threshold,
    abs(Delta_Beta) > delta_beta_threshold
  )

sig_probes <- sig_probes_df$Probe_ID

write.csv(
  results_limma,
  file.path(TABLE_DIR, "limma_deltabeta_full_results.csv"),
  row.names = FALSE
)

write.csv(
  sig_probes_df,
  file.path(TABLE_DIR, "limma_deltabeta_sig_results.csv"),
  row.names = FALSE
)

saveRDS(
  results_limma,
  file.path(RDATA_DIR, "limma_deltabeta_full_results.rds")
)

saveRDS(
  sig_probes,
  file.path(RDATA_DIR, "limma_deltabeta_sig_probes.rds")
)

# ==================================================================
# 四、相关性过滤
# ==================================================================

if (length(sig_probes) > 1) {
  sig_probe_matrix <- probe_matrix_filtered[
    ,
    sig_probes,
    drop = FALSE
  ]

  cor_matrix <- cor(
    sig_probe_matrix,
    use = "pairwise.complete.obs"
  )

  cor_matrix[
    lower.tri(cor_matrix, diag = TRUE)
  ] <- NA

  highly_correlated_pairs <- reshape2::melt(
    cor_matrix,
    na.rm = TRUE
  ) %>%
    filter(abs(value) > 0.8)

  if (nrow(highly_correlated_pairs) > 0) {
    fdr_lookup <- setNames(
      sig_probes_df$adj.P.Val,
      sig_probes_df$Probe_ID
    )

    probes_to_remove <- character(0)

    for (i in seq_len(nrow(highly_correlated_pairs))) {
      p1 <- as.character(
        highly_correlated_pairs$Var1[i]
      )

      p2 <- as.character(
        highly_correlated_pairs$Var2[i]
      )

      if (p1 %in% probes_to_remove ||
          p2 %in% probes_to_remove) {
        next
      }

      if (fdr_lookup[p1] > fdr_lookup[p2]) {
        probes_to_remove <- c(
          probes_to_remove,
          p1
        )
      } else {
        probes_to_remove <- c(
          probes_to_remove,
          p2
        )
      }
    }

    final_sig_probes <- setdiff(
      sig_probes,
      unique(probes_to_remove)
    )
  } else {
    final_sig_probes <- sig_probes
  }
} else {
  final_sig_probes <- sig_probes
}

saveRDS(
  final_sig_probes,
  file.path(RDATA_DIR, "limma_deltabeta_sig_probes_cor_filtered.rds")
)

# 后续分析沿用相关性过滤后的探针
sig_probes <- final_sig_probes

# ==================================================================
# 五、生成训练集和验证集差异探针矩阵
# ==================================================================

train_limma_filtered <- cbind(
  sample_info,
  probe_matrix_filtered[
    ,
    sig_probes,
    drop = FALSE
  ]
)

saveRDS(
  train_limma_filtered,
  file.path(RDATA_DIR, "train_limma_filtered.rds")
)

load(
  file.path(RDATA_DIR, "val_data.RData")
)

val_data <- val_data_final

val_sample_info <- val_data %>%
  select(
    Sample_Name,
    GSE_ID,
    Sample_Group
  ) %>%
  rename(Sample_ID = GSE_ID)

available_val_probes <- intersect(
  sig_probes,
  colnames(val_data)
)

val_limma_filtered <- cbind(
  val_sample_info,
  val_data[
    ,
    available_val_probes,
    drop = FALSE
  ]
)

saveRDS(
  val_limma_filtered,
  file.path(RDATA_DIR, "val_limma_filtered.rds")
)

# ==================================================================
# 六、火山图
# ==================================================================

volcano_data <- results_limma %>%
  mutate(
    log_p = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    significance = case_when(
      adj.P.Val < fdr_threshold &
        Delta_Beta > delta_beta_threshold ~
        "Hypermethylated",
      adj.P.Val < fdr_threshold &
        Delta_Beta < -delta_beta_threshold ~
        "Hypomethylated",
      TRUE ~ "Not Significant"
    )
  )

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = Delta_Beta,
    y = log_p,
    color = significance
  )
) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_vline(
    xintercept = c(
      -delta_beta_threshold,
      delta_beta_threshold
    ),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_manual(
    values = c(
      Hypermethylated = "red",
      Hypomethylated = "blue",
      `Not Significant` = "grey"
    )
  ) +
  labs(
    title = "Volcano Plot based on Delta Beta",
    subtitle = paste(
      "Comparison:",
      contrast_formula
    ),
    x = "Delta Beta",
    y = "-log10(FDR)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(FIGURE_DIR, "limma_DeltaBeta_Volcano_Plot.pdf"),
  volcano_plot,
  width = 10,
  height = 8
)

ggsave(
  file.path(FIGURE_DIR, "limma_DeltaBeta_Volcano_Plot.png"),
  volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)

save(
  results_limma,
  fit,
  fit2,
  probe_matrix_filtered,
  sig_probes,
  file = file.path(RDATA_DIR, "limma_analysis_objects.RData")
)

message("03_limma_diff.R运行完成。")