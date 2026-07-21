options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_ROOT <- Sys.getenv(
  "THESIS_ROOT",
  "/Users/mac/Desktop/TCGA_Thesis"
)

RESULT_DIR <- file.path(PROJECT_ROOT, "results")
TABLE_DIR <- file.path(RESULT_DIR, "tables")
RDATA_DIR <- file.path(RESULT_DIR, "rdata")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDATA_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "mlr3",
  "mlr3learners",
  "mlr3tuning",
  "mlr3extralearners",
  "dplyr",
  "data.table",
  "pROC",
  "xgboost",
  "lightgbm",
  "kknn"
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

library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3extralearners)
library(dplyr)
library(data.table)
library(pROC)
library(xgboost)
library(lightgbm)
library(kknn)

train_lasso_filtered <- as.data.frame(
  readRDS(
    file.path(
      RDATA_DIR,
      "train_lasso_filtered.rds"
    )
  )
)

val_data_raw <- as.data.frame(
  readRDS(
    file.path(
      RDATA_DIR,
      "val_limma_filtered.rds"
    )
  )
)

test_rdata <- file.path(
  RDATA_DIR,
  "GSE119144_test_data.RData"
)

if (!file.exists(test_rdata)) {
  stop("未找到GSE119144测试集")
}

load(test_rdata)
test_data_raw <- final_df

top_probes <- as.character(
  readRDS(
    file.path(
      RDATA_DIR,
      "RFE_Optimal_Probes.rds"
    )
  )
)

process_data <- function(data, probes) {
  common_features <- intersect(
    probes,
    colnames(data)
  )

  result <- data %>%
    select(
      Sample_Name,
      Sample_Group,
      all_of(common_features)
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

  result
}

train_data <- process_data(
  train_lasso_filtered,
  top_probes
)

val_data <- process_data(
  val_data_raw,
  top_probes
)

test_data <- process_data(
  test_data_raw,
  top_probes
)

feature_columns <- setdiff(
  colnames(train_data),
  "Sample_Group"
)

if (!identical(
  feature_columns,
  setdiff(colnames(val_data), "Sample_Group")
)) {
  stop("训练集和验证集特征不一致")
}

if (!identical(
  feature_columns,
  setdiff(colnames(test_data), "Sample_Group")
)) {
  stop("训练集和测试集特征不一致")
}

# ==================================================================
# 一、模型和参数空间
# ==================================================================

task_cv <- TaskClassif$new(
  id = "immunotherapy_response_cv",
  backend = train_data,
  target = "Sample_Group",
  positive = "1"
)

task_cv$col_roles$stratum <- "Sample_Group"

learners <- list(
  LASSO = lrn(
    "classif.glmnet",
    predict_type = "prob",
    alpha = 1
  ),
  ElasticNet = lrn(
    "classif.glmnet",
    predict_type = "prob"
  ),
  RandomForest = lrn(
    "classif.ranger",
    predict_type = "prob",
    importance = "impurity"
  ),
  Xgboost = lrn(
    "classif.xgboost",
    predict_type = "prob"
  ),
  LightGBM = lrn(
    "classif.lightgbm",
    predict_type = "prob"
  ),
  SVM = lrn(
    "classif.svm",
    predict_type = "prob",
    kernel = "radial",
    type = "C-classification"
  ),
  NeuralNetwork = lrn(
    "classif.nnet",
    predict_type = "prob",
    maxit = 100,
    MaxNWts = 5000
  ),
  KNN = lrn(
    "classif.kknn",
    predict_type = "prob"
  )
)

mtry_upper <- max(
  2,
  as.integer(
    sqrt(ncol(train_data) - 1)
  )
)

param_grids <- list(
  LASSO = ps(
    lambda = p_dbl(
      1e-3,
      10,
      logscale = TRUE
    )
  ),
  ElasticNet = ps(
    lambda = p_dbl(
      1e-3,
      10,
      logscale = TRUE
    ),
    alpha = p_dbl(
      0.1,
      0.9
    )
  ),
  RandomForest = ps(
    num.trees = p_int(
      30,
      150
    ),
    mtry = p_int(
      2,
      mtry_upper
    ),
    min.node.size = p_int(
      15,
      25
    )
  ),
  Xgboost = ps(
    nrounds = p_int(
      20,
      80
    ),
    eta = p_dbl(
      0.005,
      0.05
    ),
    max_depth = p_int(
      1,
      2
    )
  ),
  LightGBM = ps(
    num_iterations = p_int(
      20,
      80
    ),
    learning_rate = p_dbl(
      0.005,
      0.05
    ),
    max_depth = p_int(
      1,
      2
    )
  ),
  SVM = ps(
    cost = p_dbl(
      1e-4,
      0.1,
      logscale = TRUE
    ),
    gamma = p_dbl(
      1e-7,
      1e-4,
      logscale = TRUE
    )
  ),
  NeuralNetwork = ps(
    size = p_int(
      1,
      4
    ),
    decay = p_dbl(
      0.05,
      1,
      logscale = TRUE
    )
  ),
  KNN = ps(
    k = p_int(
      7,
      25
    ),
    distance = p_int(
      1,
      2
    )
  )
)

resampling_cv <- rsmp(
  "cv",
  folds = 5
)

resampling_cv$instantiate(task_cv)

best_hyperparams <- list()
cv_results_list <- list()
all_cv_pred_list <- list()
cv_roc_list <- list()

# ==================================================================
# 二、交叉验证和超参数调优
# ==================================================================

for (model_id in names(learners)) {
  message("调优模型：", model_id)

  learner <- learners[[model_id]]
  search_space <- param_grids[[model_id]]

  auto_tuner <- AutoTuner$new(
    learner = learner,
    resampling = resampling_cv,
    measure = msr("classif.auc"),
    tuner = tnr(
      "grid_search",
      resolution = 5
    ),
    terminator = trm("none"),
    search_space = search_space,
    store_models = TRUE
  )

  auto_tuner$train(task_cv)

  best_learner <- auto_tuner$learner$clone()

  best_hyperparams[[model_id]] <-
    best_learner$param_set$values

  cv_result <- resample(
    task = task_cv,
    learner = best_learner,
    resampling = resampling_cv,
    store_models = TRUE
  )

  pred_list <- lapply(
    seq_len(resampling_cv$iters),
    function(fold) {
      test_task <- task_cv$clone()
      test_task$filter(
        resampling_cv$test_set(fold)
      )

      pred <- cv_result$learners[[fold]]$predict(
        test_task
      )

      data.table(
        truth = as.character(pred$truth),
        prob.1 = pred$prob[, task_cv$positive],
        response = as.character(pred$response),
        fold = fold
      )
    }
  )

  cv_pred_dt <- rbindlist(
    pred_list,
    fill = TRUE
  )

  all_cv_pred_list[[model_id]] <- cv_pred_dt

  roc_obj <- pROC::roc(
    cv_pred_dt$truth,
    cv_pred_dt$prob.1,
    levels = c("0", "1"),
    direction = "<",
    quiet = TRUE
  )

  cv_confusion <- table(
    Truth = factor(
      cv_pred_dt$truth,
      levels = c("0", "1")
    ),
    Prediction = factor(
      cv_pred_dt$response,
      levels = c("0", "1")
    )
  )

  tn <- cv_confusion["0", "0"]
  fp <- cv_confusion["0", "1"]
  fn <- cv_confusion["1", "0"]
  tp <- cv_confusion["1", "1"]

  accuracy <- (tp + tn) / sum(cv_confusion)

  sensitivity <- if ((tp + fn) > 0) {
    tp / (tp + fn)
  } else {
    0
  }

  specificity <- if ((tn + fp) > 0) {
    tn / (tn + fp)
  } else {
    0
  }

  precision <- if ((tp + fp) > 0) {
    tp / (tp + fp)
  } else {
    0
  }

  f1 <- if ((precision + sensitivity) > 0) {
    2 * precision * sensitivity /
      (precision + sensitivity)
  } else {
    0
  }

  mcc_den <- sqrt(
    (tp + fp) *
      (tp + fn) *
      (tn + fp) *
      (tn + fn)
  )

  mcc <- if (mcc_den > 0) {
    (tp * tn - fp * fn) / mcc_den
  } else {
    0
  }

  cv_results_list[[model_id]] <- data.frame(
    Model = model_id,
    CV_AUC = round(as.numeric(roc_obj$auc), 4),
    CV_Accuracy = round(accuracy, 4),
    CV_Sensitivity = round(sensitivity, 4),
    CV_Specificity = round(specificity, 4),
    CV_Precision = round(precision, 4),
    CV_F1 = round(f1, 4),
    CV_MCC = round(mcc, 4),
    stringsAsFactors = FALSE
  )

  cv_roc_list[[model_id]] <- data.frame(
    model = model_id,
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    auc = as.numeric(roc_obj$auc)
  )
}

cv_results_df <- bind_rows(
  cv_results_list
) %>%
  arrange(desc(CV_AUC))

write.csv(
  cv_results_df,
  file.path(TABLE_DIR, "All_Models_CV_Performance.csv"),
  row.names = FALSE
)

save(
  best_hyperparams,
  all_cv_pred_list,
  cv_roc_list,
  cv_results_df,
  file = file.path(
    RDATA_DIR,
    "CV_Tuning_Results.RData"
  )
)

# ==================================================================
# 三、全量训练、验证集阈值优化和外部测试集评估
# ==================================================================

task_full <- TaskClassif$new(
  id = "train_full",
  backend = train_data,
  target = "Sample_Group",
  positive = "1"
)

val_task <- TaskClassif$new(
  id = "validation",
  backend = val_data,
  target = "Sample_Group",
  positive = "1"
)

test_task <- TaskClassif$new(
  id = "test",
  backend = test_data,
  target = "Sample_Group",
  positive = "1"
)

val_results_list <- list()
val_roc_list <- list()
val_pred_list <- list()

test_results_list <- list()
test_roc_list <- list()
test_pred_list <- list()

calculate_metrics <- function(truth, predicted) {
  confusion <- table(
    Truth = factor(
      truth,
      levels = c("0", "1")
    ),
    Prediction = factor(
      predicted,
      levels = c("0", "1")
    )
  )

  tn <- confusion["0", "0"]
  fp <- confusion["0", "1"]
  fn <- confusion["1", "0"]
  tp <- confusion["1", "1"]

  accuracy <- (tp + tn) / sum(confusion)

  sensitivity <- if ((tp + fn) > 0) {
    tp / (tp + fn)
  } else {
    0
  }

  specificity <- if ((tn + fp) > 0) {
    tn / (tn + fp)
  } else {
    0
  }

  precision <- if ((tp + fp) > 0) {
    tp / (tp + fp)
  } else {
    0
  }

  f1 <- if ((precision + sensitivity) > 0) {
    2 * precision * sensitivity /
      (precision + sensitivity)
  } else {
    0
  }

  mcc_den <- sqrt(
    (tp + fp) *
      (tp + fn) *
      (tn + fp) *
      (tn + fn)
  )

  mcc <- if (mcc_den > 0) {
    (tp * tn - fp * fn) / mcc_den
  } else {
    0
  }

  data.frame(
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    MCC = mcc
  )
}

for (model_id in names(learners)) {
  message("训练模型：", model_id)

  learner <- learners[[model_id]]$clone()
  learner$param_set$values <-
    best_hyperparams[[model_id]]

  learner$train(task_full)

  val_pred <- learner$predict(val_task)

  val_pred_dt <- data.table(
    truth = as.character(val_pred$truth),
    prob.1 = val_pred$prob[, "1"],
    response = as.character(val_pred$response)
  )

  val_pred_list[[model_id]] <- val_pred_dt

  val_roc_obj <- pROC::roc(
    val_pred_dt$truth,
    val_pred_dt$prob.1,
    levels = c("0", "1"),
    direction = "<",
    quiet = TRUE
  )

  val_roc_list[[model_id]] <- data.frame(
    model = model_id,
    fpr = 1 - val_roc_obj$specificities,
    tpr = val_roc_obj$sensitivities,
    auc = as.numeric(val_roc_obj$auc)
  )

  thresholds <- seq(
    0.01,
    0.99,
    by = 0.01
  )

  f1_scores <- vapply(
    thresholds,
    function(threshold) {
      pred_label <- ifelse(
        val_pred_dt$prob.1 >= threshold,
        "1",
        "0"
      )

      metric <- calculate_metrics(
        val_pred_dt$truth,
        pred_label
      )

      metric$F1
    },
    numeric(1)
  )

  optimal_threshold <- thresholds[
    which.max(f1_scores)
  ]

  test_pred <- learner$predict(test_task)

  test_pred_dt <- data.table(
    truth = as.character(test_pred$truth),
    prob.1 = test_pred$prob[, "1"]
  )

  test_pred_dt$response <- ifelse(
    test_pred_dt$prob.1 >= optimal_threshold,
    "1",
    "0"
  )

  test_pred_list[[model_id]] <- test_pred_dt

  test_roc_obj <- pROC::roc(
    test_pred_dt$truth,
    test_pred_dt$prob.1,
    levels = c("0", "1"),
    direction = "<",
    quiet = TRUE
  )

  test_metric <- calculate_metrics(
    test_pred_dt$truth,
    test_pred_dt$response
  )

  test_roc_list[[model_id]] <- data.frame(
    model = model_id,
    fpr = 1 - test_roc_obj$specificities,
    tpr = test_roc_obj$sensitivities,
    auc = as.numeric(test_roc_obj$auc)
  )

  test_results_list[[model_id]] <- data.frame(
    Model = model_id,
    Test_AUC = round(
      as.numeric(test_roc_obj$auc),
      4
    ),
    Test_Accuracy = round(
      test_metric$Accuracy,
      4
    ),
    Test_Sensitivity = round(
      test_metric$Sensitivity,
      4
    ),
    Test_Specificity = round(
      test_metric$Specificity,
      4
    ),
    Test_Precision = round(
      test_metric$Precision,
      4
    ),
    Test_F1 = round(
      test_metric$F1,
      4
    ),
    Test_MCC = round(
      test_metric$MCC,
      4
    ),
    Optimal_Threshold = round(
      optimal_threshold,
      2
    ),
    stringsAsFactors = FALSE
  )

  val_metric <- calculate_metrics(
    val_pred_dt$truth,
    val_pred_dt$response
  )

  val_results_list[[model_id]] <- data.frame(
    Model = model_id,
    Val_AUC = round(
      as.numeric(val_roc_obj$auc),
      4
    ),
    Val_Accuracy = round(
      val_metric$Accuracy,
      4
    ),
    Val_Sensitivity = round(
      val_metric$Sensitivity,
      4
    ),
    Val_Specificity = round(
      val_metric$Specificity,
      4
    ),
    Val_Precision = round(
      val_metric$Precision,
      4
    ),
    Val_F1 = round(
      val_metric$F1,
      4
    ),
    Val_MCC = round(
      val_metric$MCC,
      4
    ),
    stringsAsFactors = FALSE
  )
}

val_results_df <- bind_rows(
  val_results_list
) %>%
  arrange(desc(Val_AUC))

test_results_df <- bind_rows(
  test_results_list
) %>%
  arrange(desc(Test_AUC))

write.csv(
  val_results_df,
  file.path(
    TABLE_DIR,
    "All_Models_Validation_Performance.csv"
  ),
  row.names = FALSE
)

write.csv(
  test_results_df,
  file.path(
    TABLE_DIR,
    "All_Models_Test_Performance.csv"
  ),
  row.names = FALSE
)

save(
  val_results_df,
  val_roc_list,
  val_pred_list,
  file = file.path(
    RDATA_DIR,
    "Validation_Evaluation_Results.RData"
  )
)

save(
  test_results_df,
  test_roc_list,
  test_pred_list,
  file = file.path(
    RDATA_DIR,
    "Test_Evaluation_Results.RData"
  )
)

message("05_ML_model.R运行完成。")