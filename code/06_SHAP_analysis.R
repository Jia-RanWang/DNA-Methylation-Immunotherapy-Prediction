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
  "kernelshap",
  "shapviz",
  "mlr3",
  "mlr3learners",
  "dplyr",
  "data.table",
  "pheatmap"
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

library(kernelshap)
library(shapviz)
library(mlr3)
library(mlr3learners)
library(dplyr)
library(data.table)
library(pheatmap)

load(
  file.path(
    RDATA_DIR,
    "CV_Tuning_Results.RData"
  )
)

train_raw <- as.data.frame(
  readRDS(
    file.path(
      RDATA_DIR,
      "train_lasso_filtered.rds"
    )
  )
)

load(
  file.path(
    RDATA_DIR,
    "GSE119144_test_data.RData"
  )
)

test_raw <- final_df

top_features <- as.character(
  readRDS(
    file.path(
      RDATA_DIR,
      "RFE_Optimal_Probes.rds"
    )
  )
)

process_shap_data <- function(data, probes) {
  common_probes <- intersect(
    probes,
    colnames(data)
  )

  data %>%
    select(
      Sample_Name,
      Sample_Group,
      all_of(common_probes)
    ) %>%
    mutate(
      Sample_Group = factor(
        ifelse(
          Sample_Group == "R",
          1,
          0
        ),
        levels = c(0, 1)
      )
    ) %>%
    select(-Sample_Name)
}

train_data <- process_shap_data(
  train_raw,
  top_features
)

test_data <- process_shap_data(
  test_raw,
  top_features
)

best_model_id <- "RandomForest"
rf_best_params <- best_hyperparams[[best_model_id]]

task_train <- TaskClassif$new(
  id = "train_for_shap",
  backend = train_data,
  target = "Sample_Group",
  positive = "1"
)

final_rf_learner <- lrn(
  "classif.ranger",
  id = best_model_id,
  predict_type = "prob",
  importance = "impurity"
)

final_rf_learner$param_set$values <- rf_best_params
final_rf_learner$train(task_train)

x_train <- train_data %>%
  select(-Sample_Group)

x_test <- test_data %>%
  select(-Sample_Group)

set.seed(123)

background_index <- sample(
  seq_len(nrow(x_train)),
  size = min(100, nrow(x_train))
)

background_x <- x_train[
  background_index,
  ,
  drop = FALSE
]

prediction_function <- function(object, newdata) {
  newdata <- as.data.frame(newdata)

  prediction <- object$predict_newdata(
    newdata = newdata
  )

  as.numeric(
    prediction$prob[, "1"]
  )
}

explain_kernel <- kernelshap(
  object = final_rf_learner,
  X = x_test,
  pred_fun = prediction_function,
  bg_X = background_x
)

shap_value <- shapviz(explain_kernel)
shap_matrix <- shap_value$S

shap_importance_df <- data.frame(
  ProbeID = colnames(shap_matrix),
  Importance = colMeans(
    abs(shap_matrix),
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Importance))

write.csv(
  shap_importance_df,
  file.path(TABLE_DIR, "shap_importance.csv"),
  row.names = FALSE
)

saveRDS(
  shap_importance_df,
  file.path(RDATA_DIR, "shap_importance.rds")
)

saveRDS(
  explain_kernel,
  file.path(RDATA_DIR, "SHAP_explain_kernel.rds")
)

# ==================================================================
# SHAP图形
# ==================================================================

pdf(
  file.path(
    FIGURE_DIR,
    "SHAP_RandomForest_sv_force.pdf"
  ),
  width = 7,
  height = 5
)

print(
  sv_force(
    shap_value,
    row_id = min(2, nrow(shap_matrix)),
    size = 9
  ) +
    ggtitle("Force Plot for RandomForest") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
)

dev.off()

pdf(
  file.path(
    FIGURE_DIR,
    "SHAP_RandomForest_importance_beeswarm.pdf"
  ),
  width = 7,
  height = 5
)

print(
  sv_importance(
    shap_value,
    kind = "beeswarm",
    show_numbers = FALSE
  ) +
    ggtitle("Feature Importance (RandomForest)") +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
)

dev.off()

pdf(
  file.path(
    FIGURE_DIR,
    "SHAP_RandomForest_importance_bar.pdf"
  ),
  width = 7,
  height = 5
)

print(
  sv_importance(
    shap_value,
    kind = "bar",
    show_numbers = FALSE,
    fill = "#fca50a"
  ) +
    ggtitle("Feature Importance (RandomForest)") +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
)

dev.off()

top_probes <- head(
  shap_importance_df$ProbeID,
  10
)

dir.create(
  file.path(FIGURE_DIR, "SHAP_Dependence_Plots"),
  recursive = TRUE,
  showWarnings = FALSE
)

for (probe in top_probes) {
  dependence_plot <- sv_dependence(
    shap_value,
    v = probe
  ) +
    theme_bw() +
    ggtitle(
      paste(
        "SHAP Dependence for",
        probe
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      )
    )

  ggsave(
    file.path(
      FIGURE_DIR,
      "SHAP_Dependence_Plots",
      paste0(
        "SHAP_dependence_",
        probe,
        ".pdf"
      )
    ),
    dependence_plot,
    width = 6,
    height = 5
  )
}

pdf(
  file.path(
    FIGURE_DIR,
    "SHAP_RandomForest_waterfall.pdf"
  ),
  width = 6,
  height = 6
)

print(
  sv_waterfall(
    shap_value,
    row_id = min(2, nrow(shap_matrix)),
    fill_colors = c("#f7d13d", "#a52c60")
  ) +
    theme_bw() +
    ggtitle("Waterfall Plot for RandomForest") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
)

dev.off()

top20_probes <- head(
  shap_importance_df$ProbeID,
  min(20, ncol(shap_matrix))
)

shap_top20 <- shap_matrix[
  ,
  top20_probes,
  drop = FALSE
]

if (is.null(rownames(shap_top20))) {
  rownames(shap_top20) <- paste0(
    "Sample_",
    seq_len(nrow(shap_top20))
  )
}

set.seed(123)

sample_ids <- sample(
  rownames(shap_top20),
  size = min(20, nrow(shap_top20))
)

pdf(
  file.path(
    FIGURE_DIR,
    "SHAP_RandomForest_heatmap.pdf"
  ),
  width = 10,
  height = 8
)

pheatmap(
  shap_top20[sample_ids, , drop = FALSE],
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "SHAP Value Heatmap",
  angle_col = 45,
  fontsize_col = 9,
  fontsize_row = 10
)

dev.off()

message("06_SHAP_analysis.R运行完成。")