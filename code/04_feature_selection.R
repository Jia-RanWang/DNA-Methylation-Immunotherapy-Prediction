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
  "dplyr",
  "tibble",
  "limma",
  "glmnet",
  "randomForest",
  "xgboost",
  "caret",
  "Boruta",
  "lightgbm",
  "reshape2",
  "rmda"
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
library(tibble)
library(limma)
library(glmnet)
library(randomForest)
library(xgboost)
library(caret)
library(Boruta)
library(lightgbm)
library(reshape2)
library(rmda)

train_limma_filtered <- readRDS(
  file.path(RDATA_DIR, "train_limma_filtered.rds")
)

val_limma_filtered <- readRDS(
  file.path(RDATA_DIR, "val_limma_filtered.rds")
)

train_data <- train_limma_filtered

probe_matrix <- train_data %>%
  select(
    -Sample_Name,
    -Sample_ID,
    -Sample_Group
  )

x <- as.matrix(probe_matrix)
y <- factor(train_data$Sample_Group)

if (!all(levels(y) %in% c("NR", "R"))) {
  y <- factor(y, levels = c("NR", "R"))
}

y_bin <- ifelse(y == "R", 1, 0)

save_feature <- function(object, filename) {
  saveRDS(
    object,
    file.path(RDATA_DIR, filename)
  )

  saveRDS(
    object,
    file.path(PROJECT_ROOT, filename)
  )
}

# ==================================================================
# 一、LASSO
# ==================================================================

set.seed(123)

lasso_cv <- cv.glmnet(
  x = x,
  y = y_bin,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "auc"
)

final_lasso <- glmnet(
  x = x,
  y = y_bin,
  family = "binomial",
  alpha = 1,
  lambda = lasso_cv$lambda.min
)

lasso_coef <- as.matrix(
  coef(final_lasso)
)

lasso_coef_df <- data.frame(
  Probe_ID = rownames(lasso_coef),
  LASSO_Coefficient = as.numeric(lasso_coef[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(
    Probe_ID != "(Intercept)",
    LASSO_Coefficient != 0
  ) %>%
  arrange(
    desc(abs(LASSO_Coefficient))
  )

core_probes <- head(
  lasso_coef_df$Probe_ID,
  min(30, nrow(lasso_coef_df))
)

save_feature(
  core_probes,
  "LASSO_Optimal_Probes.rds"
)

save_feature(
  core_probes,
  "LASSO_Optimal_Probes2.rds"
)

save_feature(
  lasso_coef_df,
  "LASSO_Feature_Importance.rds"
)

train_lasso_filtered <- train_data %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(core_probes)
  )

save_feature(
  train_lasso_filtered,
  "train_lasso_filtered.rds"
)

save(
  lasso_cv,
  final_lasso,
  file = file.path(
    RDATA_DIR,
    "LASSO_CV_Objects.RData"
  )
)

# ==================================================================
# 二、随机森林
# ==================================================================

set.seed(123)

rf_model <- randomForest(
  x = as.data.frame(
    train_lasso_filtered %>%
      select(-Sample_Name, -Sample_Group)
  ),
  y = factor(train_lasso_filtered$Sample_Group),
  ntree = 500,
  mtry = sqrt(ncol(core_probes)),
  importance = TRUE,
  proximity = FALSE,
  do.trace = 100,
  cv.fold = 10
)

rf_importance <- as.data.frame(
  importance(rf_model)
) %>%
  rownames_to_column("Probe_ID") %>%
  arrange(
    desc(MeanDecreaseGini)
  )

rf_importance_mean <- mean(
  rf_importance$MeanDecreaseGini,
  na.rm = TRUE
)

rf_optimal_probes <- rf_importance %>%
  filter(
    MeanDecreaseGini > rf_importance_mean
  ) %>%
  slice_head(n = 30) %>%
  pull(Probe_ID)

save_feature(
  rf_importance,
  "RF_Feature_Importance.rds"
)

save_feature(
  rf_optimal_probes,
  "RF_Optimal_Probes.rds"
)

train_rf_filtered <- train_lasso_filtered %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(rf_optimal_probes)
  )

save_feature(
  train_rf_filtered,
  "train_rf_filtered.rds"
)

# ==================================================================
# 三、XGBoost
# ==================================================================

x_lasso <- train_lasso_filtered %>%
  select(
    -Sample_Name,
    -Sample_Group
  ) %>%
  as.matrix()

y_lasso <- ifelse(
  train_lasso_filtered$Sample_Group == "R",
  1,
  0
)

dtrain <- xgb.DMatrix(
  data = x_lasso,
  label = y_lasso
)

xgb_params <- list(
  booster = "gbtree",
  max_depth = 3,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  objective = "binary:logistic",
  eval_metric = "auc"
)

set.seed(123)

xgb_cv <- xgb.cv(
  params = xgb_params,
  data = dtrain,
  nrounds = 500,
  nfold = 10,
  early_stopping_rounds = 50,
  verbose = 100
)

best_iter <- if (!is.null(xgb_cv$best_iteration)) {
  xgb_cv$best_iteration
} else if (!is.null(xgb_cv$early_stop$best_iteration)) {
  xgb_cv$early_stop$best_iteration
} else {
  100
}

xgb_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = best_iter,
  evals = list(train = dtrain),
  verbose = 100
)

xgb_importance <- xgb.importance(
  model = xgb_model
) %>%
  rename(
    Probe_ID = Feature,
    Weight = Cover
  ) %>%
  arrange(
    desc(Gain)
  )

xgb_gain_threshold <- quantile(
  xgb_importance$Gain,
  0.7,
  na.rm = TRUE
)

xgb_optimal_probes <- xgb_importance %>%
  filter(Gain > xgb_gain_threshold) %>%
  slice_head(n = 30) %>%
  pull(Probe_ID)

save_feature(
  xgb_importance,
  "XGBoost_Feature_Importance.rds"
)

save_feature(
  xgb_optimal_probes,
  "XGBoost_Optimal_Probes.rds"
)

train_xgb_filtered <- train_lasso_filtered %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(xgb_optimal_probes)
  )

save_feature(
  train_xgb_filtered,
  "train_xgb_filtered.rds"
)

save(
  xgb_cv,
  xgb_model,
  file = file.path(
    RDATA_DIR,
    "XGBoost_Objects.RData"
  )
)

# ==================================================================
# 四、Elastic Net
# ==================================================================

set.seed(123)

cv_enet <- cv.glmnet(
  x = x_lasso,
  y = y_lasso,
  family = "binomial",
  alpha = 0.5,
  nfolds = 10,
  type.measure = "auc"
)

enet_coef <- as.matrix(
  coef(
    cv_enet$glmnet.fit,
    s = cv_enet$lambda.1se
  )
)

enet_importance <- data.frame(
  Probe_ID = rownames(enet_coef),
  Coefficient = as.numeric(enet_coef[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(
    Probe_ID != "(Intercept)",
    Coefficient != 0
  ) %>%
  arrange(
    desc(abs(Coefficient))
  )

enet_optimal_probes <- head(
  enet_importance$Probe_ID,
  min(30, nrow(enet_importance))
)

save_feature(
  enet_importance,
  "ElasticNet_Feature_Importance.rds"
)

save_feature(
  enet_optimal_probes,
  "ElasticNet_Optimal_Probes.rds"
)

train_enet_filtered <- train_lasso_filtered %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(enet_optimal_probes)
  )

save_feature(
  train_enet_filtered,
  "train_enet_filtered.rds"
)

save(
  cv_enet,
  file = file.path(
    RDATA_DIR,
    "ElasticNet_CV_Objects.RData"
  )
)

# ==================================================================
# 五、RF-RFE
# ==================================================================

rfe_input <- train_lasso_filtered %>%
  select(
    -Sample_Name,
    -Sample_Group
  )

rfe_y <- factor(
  train_lasso_filtered$Sample_Group,
  levels = c("NR", "R")
)

rf_funcs_auc <- rfFuncs
rf_funcs_auc$summary <- twoClassSummary

set.seed(123)

rfe_control <- rfeControl(
  functions = rf_funcs_auc,
  method = "cv",
  number = 10,
  returnResamp = "final",
  verbose = FALSE
)

rfe_sizes <- c(27, 28, 29, 30, 31, 32, 33)

rfe_result <- tryCatch({
  rfe(
    x = rfe_input,
    y = rfe_y,
    sizes = rfe_sizes,
    rfeControl = rfe_control,
    metric = "ROC"
  )
}, error = function(e) {
  warning("RFE运行失败，使用RF重要性排序替代：", conditionMessage(e))
  NULL
})

if (!is.null(rfe_result)) {
  best_features <- predictors(rfe_result)

  set.seed(123)

  rf_best <- randomForest(
    x = rfe_input[, best_features, drop = FALSE],
    y = rfe_y,
    importance = TRUE,
    ntree = 500,
    proximity = FALSE
  )

  rfe_importance <- as.data.frame(
    importance(rf_best)
  ) %>%
    rownames_to_column("Probe_ID") %>%
    arrange(desc(MeanDecreaseGini))
} else {
  rfe_importance <- rf_importance
}

rfe_top_probes <- head(
  rfe_importance$Probe_ID,
  min(32, nrow(rfe_importance))
)

save_feature(
  rfe_result,
  "RFE_Result.rds"
)

save_feature(
  rfe_importance,
  "RFE_Feature_Importance.rds"
)

save_feature(
  rfe_top_probes,
  "RFE_Optimal_Probes.rds"
)

train_rfe_filtered <- train_lasso_filtered %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(rfe_top_probes)
  )

save_feature(
  train_rfe_filtered,
  "train_rfe_filtered.rds"
)

# ==================================================================
# 六、Boruta
# ==================================================================

boruta_x <- train_lasso_filtered %>%
  select(
    -Sample_Name,
    -Sample_Group
  )

boruta_y <- factor(
  train_lasso_filtered$Sample_Group,
  levels = c("NR", "R")
)

set.seed(123)

boruta_result <- Boruta(
  x = boruta_x,
  y = boruta_y,
  doTrace = 1,
  ntree = 500,
  maxRuns = 100,
  pValue = 0.01
)

boruta_initial_probes <- getSelectedAttributes(
  boruta_result,
  withTentative = FALSE
)

boruta_importance <- attStats(
  boruta_result
) %>%
  rownames_to_column("Probe_ID") %>%
  rename(
    MeanDecreaseAccuracy = meanImp,
    Decision = decision
  ) %>%
  arrange(
    desc(MeanDecreaseAccuracy)
  )

boruta_optimal_probes <- boruta_importance %>%
  filter(
    Probe_ID %in% boruta_initial_probes
  ) %>%
  slice_head(n = 30) %>%
  pull(Probe_ID)

save_feature(
  boruta_result,
  "Boruta_Result.rds"
)

save_feature(
  boruta_importance,
  "Boruta_Feature_Importance.rds"
)

save_feature(
  boruta_optimal_probes,
  "Boruta_Optimal_Probes.rds"
)

train_boruta_filtered <- train_lasso_filtered %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(boruta_optimal_probes)
  )

save_feature(
  train_boruta_filtered,
  "train_boruta_filtered.rds"
)

# ==================================================================
# 七、LightGBM
# ==================================================================

lgb_x <- train_lasso_filtered %>%
  select(
    -Sample_Name,
    -Sample_Group
  ) %>%
  as.matrix()

lgb_y <- ifelse(
  train_lasso_filtered$Sample_Group == "R",
  1,
  0
)

lgb_train <- lgb.Dataset(
  data = lgb_x,
  label = lgb_y,
  free_raw_data = FALSE
)

lgb_params <- list(
  objective = "binary",
  metric = "auc",
  boosting_type = "gbdt",
  learning_rate = 0.05,
  num_leaves = 31,
  max_depth = -1,
  min_child_samples = 20,
  subsample = 0.8,
  colsample_bytree = 0.8,
  reg_alpha = 0.1,
  reg_lambda = 0.1,
  verbosity = -1,
  num_threads = max(1, parallel::detectCores() - 1)
)

set.seed(123)

lgb_cv_result <- lgb.cv(
  params = lgb_params,
  data = lgb_train,
  nfold = 10,
  nrounds = 500,
  early_stopping_rounds = 50,
  stratified = TRUE,
  eval_freq = 10
)

best_lgb_iter <- if (!is.null(lgb_cv_result$best_iter)) {
  lgb_cv_result$best_iter
} else {
  100
}

lgb_final_model <- lgb.train(
  params = lgb_params,
  data = lgb_train,
  nrounds = best_lgb_iter,
  verbose = -1
)

lgb_importance <- lgb.importance(
  model = lgb_final_model,
  percentage = TRUE
) %>%
  rename(
    Probe_ID = Feature,
    Gain_Importance = Gain
  ) %>%
  arrange(
    desc(Gain_Importance)
  )

lgb_importance$Cumulative_Gain <- cumsum(
  lgb_importance$Gain_Importance
)

lgb_optimal_probes <- lgb_importance %>%
  filter(Cumulative_Gain <= 80) %>%
  slice_head(n = 30) %>%
  pull(Probe_ID)

if (length(lgb_optimal_probes) < 5) {
  lgb_optimal_probes <- lgb_importance %>%
    slice_head(n = 10) %>%
    pull(Probe_ID)
}

save_feature(
  lgb_importance,
  "LightGBM_Feature_Importance.rds"
)

save_feature(
  lgb_optimal_probes,
  "LightGBM_Optimal_Probes.rds"
)

train_lgb_filtered <- train_lasso_filtered %>%
  select(
    Sample_Name,
    Sample_Group,
    all_of(lgb_optimal_probes)
  )

save_feature(
  train_lgb_filtered,
  "train_lgb_filtered.rds"
)

save(
  lgb_cv_result,
  lgb_final_model,
  file = file.path(
    RDATA_DIR,
    "LightGBM_Objects.RData"
  )
)

# ==================================================================
# 八、特征筛选方法验证
# ==================================================================

method_names <- c(
  "LASSO",
  "RF",
  "XGBoost",
  "ElasticNet",
  "RFE",
  "Boruta",
  "LightGBM"
)

method_file_names <- c(
  LASSO = "LASSO_Optimal_Probes.rds",
  RF = "RF_Optimal_Probes.rds",
  XGBoost = "XGBoost_Optimal_Probes.rds",
  ElasticNet = "ElasticNet_Optimal_Probes.rds",
  RFE = "RFE_Optimal_Probes.rds",
  Boruta = "Boruta_Optimal_Probes.rds",
  LightGBM = "LightGBM_Optimal_Probes.rds"
)

evaluation_list <- list()

for (method in method_names) {
  message("验证特征筛选方法：", method)

  probe_file <- file.path(
    RDATA_DIR,
    method_file_names[[method]]
  )

  if (!file.exists(probe_file)) {
    next
  }

  opt_probes <- readRDS(probe_file)

  common_probes <- Reduce(
    intersect,
    list(
      opt_probes,
      colnames(train_limma_filtered),
      colnames(val_limma_filtered)
    )
  )

  if (length(common_probes) == 0) {
    next
  }

  train_features <- train_limma_filtered %>%
    select(
      all_of(common_probes),
      Sample_Group
    ) %>%
    mutate(
      Sample_Group = ifelse(
        Sample_Group == "R",
        1,
        0
      )
    )

  val_features <- val_limma_filtered %>%
    select(
      all_of(common_probes),
      Sample_Group
    ) %>%
    mutate(
      Sample_Group = ifelse(
        Sample_Group == "R",
        1,
        0
      )
    )

  x_train <- as.matrix(
    train_features %>% select(-Sample_Group)
  )

  y_train <- train_features$Sample_Group

  x_val <- as.matrix(
    val_features %>% select(-Sample_Group)
  )

  y_val <- val_features$Sample_Group

  set.seed(123)

  start_time <- Sys.time()

  cv_lr <- cv.glmnet(
    x_train,
    y_train,
    family = "binomial",
    alpha = 1
  )

  lr_model <- glmnet(
    x_train,
    y_train,
    family = "binomial",
    alpha = 1,
    lambda = cv_lr$lambda.min,
    standardize = TRUE
  )

  val_time <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  pred_prob <- as.numeric(
    predict(
      lr_model,
      newx = x_val,
      type = "response"
    )
  )

  pred_label <- ifelse(
    pred_prob > 0.5,
    1,
    0
  )

  roc_object <- pROC::roc(
    y_val,
    pred_prob,
    quiet = TRUE
  )

  confusion <- table(
    Truth = factor(y_val, levels = c(0, 1)),
    Prediction = factor(pred_label, levels = c(0, 1))
  )

  tn <- confusion["0", "0"]
  fp <- confusion["0", "1"]
  fn <- confusion["1", "0"]
  tp <- confusion["1", "1"]

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

  evaluation_list[[method]] <- data.frame(
    Method = method,
    Feature_Count = length(common_probes),
    Val_AUC = round(as.numeric(roc_object$auc), 4),
    Val_Accuracy = round(
      mean(pred_label == y_val),
      4
    ),
    Val_Sensitivity = round(sensitivity, 4),
    Val_Specificity = round(specificity, 4),
    Val_Precision = round(precision, 4),
    Val_F1 = round(f1, 4),
    Train_Time = round(val_time, 2),
    stringsAsFactors = FALSE
  )
}

feature_evaluation_results <- bind_rows(
  evaluation_list
) %>%
  arrange(desc(Val_AUC))

write.csv(
  feature_evaluation_results,
  file.path(
    TABLE_DIR,
    "Feature_Methods_Evaluation_Results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  feature_evaluation_results,
  file.path(
    RDATA_DIR,
    "Feature_Methods_Evaluation_Results.rds"
  )
)

# ==================================================================
# 九、特征筛选方法DCA
# ==================================================================

dca_list <- list()

for (method in method_names) {
  probe_file <- file.path(
    RDATA_DIR,
    method_file_names[[method]]
  )

  if (!file.exists(probe_file)) {
    next
  }

  opt_probes <- readRDS(probe_file)

  common_probes <- Reduce(
    intersect,
    list(
      opt_probes,
      colnames(train_limma_filtered),
      colnames(val_limma_filtered)
    )
  )

  if (length(common_probes) < 2) {
    next
  }

  train_dca <- train_limma_filtered %>%
    select(
      all_of(common_probes),
      Sample_Group
    ) %>%
    mutate(
      Result = ifelse(
        Sample_Group == "R",
        1,
        0
      )
    )

  val_dca <- val_limma_filtered %>%
    select(
      all_of(common_probes),
      Sample_Group
    ) %>%
    mutate(
      Result = ifelse(
        Sample_Group == "R",
        1,
        0
      )
    )

  x_train <- as.matrix(
    train_dca %>% select(all_of(common_probes))
  )

  y_train <- train_dca$Result

  x_val <- as.matrix(
    val_dca %>% select(all_of(common_probes))
  )

  lr_model <- glmnet(
    x_train,
    y_train,
    family = "binomial",
    alpha = 0,
    lambda = cv.glmnet(
      x_train,
      y_train,
      family = "binomial",
      alpha = 0
    )$lambda.min
  )

  val_dca$pred_prob <- as.numeric(
    predict(
      lr_model,
      newx = x_val,
      type = "response"
    )
  )

  dca_data <- val_dca %>%
    select(Result, pred_prob)

  dca_list[[method]] <- decision_curve(
    Result ~ pred_prob,
    data = dca_data,
    study.design = "cohort",
    bootstraps = 500,
    confidence.intervals = "none",
    thresholds = seq(0, 1, by = 0.01
    )
  )
}

save(
  dca_list,
  file = file.path(
    RDATA_DIR,
    "Feature_Selection_DCA_Objects.RData"
  )
)

message("04_feature_selection.R运行完成。")