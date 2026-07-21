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

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "data.table",
  "pROC",
  "patchwork",
  "RColorBrewer",
  "fmsb",
  "scales",
  "caret",
  "rmda",
  "PMCMRplus",
  "reshape2",
  "grid",
  "gridExtra",
  "viridis"
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

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pROC)
library(patchwork)
library(RColorBrewer)
library(fmsb)
library(scales)
library(caret)
library(rmda)
library(PMCMRplus)
library(reshape2)
library(grid)
library(gridExtra)
library(viridis)

# ==================================================================
# 一、TCGA关键指标图
# ==================================================================

tmb_data <- read.csv(
  file.path(
    TABLE_DIR,
    "TCGA_All_Cancer_TMB_Results.csv"
  )
)

expression_data <- read.csv(
  file.path(
    TABLE_DIR,
    "TCGA_All_Cancer_GEP_TGFB_Results.csv"
  )
)

plot_tcga_score <- function(data,
                            score_column,
                            score_label,
                            title,
                            output_prefix) {
  data <- data %>%
    filter(
      !is.na(.data[[score_column]])
    )

  stats <- data %>%
    group_by(CancerType) %>%
    summarise(
      median_score = round(
        median(
          .data[[score_column]],
          na.rm = TRUE
        ),
        2
      ),
      sample_n = sum(
        !is.na(.data[[score_column]])
      ),
      .groups = "drop"
    ) %>%
    arrange(median_score) %>%
    mutate(
      CancerLabel = paste0(
        CancerType,
        " (n=",
        sample_n,
        ")"
      )
    )

  data$CancerLabel <- factor(
    paste0(
      data$CancerType,
      " (n=",
      ave(
        !is.na(data[[score_column]]),
        data$CancerType,
        FUN = sum
      ),
      ")"
    ),
    levels = stats$CancerLabel
  )

  label_data <- stats

  plot <- ggplot(
    data,
    aes(
      x = CancerLabel,
      y = .data[[score_column]],
      fill = CancerLabel
    )
  ) +
    geom_jitter(
      alpha = 0.4,
      size = 1.2,
      color = "gray70",
      width = 0.2
    ) +
    geom_boxplot(
      width = 0.6,
      color = "black",
      linewidth = 0.8,
      outlier.shape = NA,
      alpha = 0.8
    ) +
    geom_text(
      data = label_data,
      aes(
        x = CancerLabel,
        y = median_score,
        label = paste0(
          "Med: ",
          median_score
        )
      ),
      inherit.aes = FALSE,
      color = "darkred",
      size = 2.8,
      fontface = "bold"
    ) +
    scale_fill_viridis_d(
      option = "plasma",
      begin = 0.1,
      end = 0.9
    ) +
    labs(
      x = "Cancer Type (TCGA) + Sample Size",
      y = score_label,
      title = title
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      ),
      axis.title = element_text(
        size = 14,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        size = 9
      ),
      axis.text.y = element_text(
        size = 11
      ),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(
        10,
        10,
        30,
        20
      )
    )

  ggsave(
    file.path(
      FIGURE_DIR,
      paste0(output_prefix, ".png")
    ),
    plot,
    width = 16,
    height = 9,
    dpi = 300,
    bg = "white"
  )

  write.csv(
    stats,
    file.path(
      TABLE_DIR,
      paste0(output_prefix, "_Stats.csv")
    ),
    row.names = FALSE
  )

  plot
}

tmb_plot <- plot_tcga_score(
  tmb_data,
  "TMB",
  "Tumor Mutational Burden",
  "TMB Distribution Across TCGA Cancers",
  "TCGA_TMB_Distribution_Boxplot"
)

tmb_plot <- tmb_plot +
  geom_hline(
    yintercept = 10,
    color = "firebrick",
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 1,
    y = 15,
    label = "TMB-H Threshold (≥10 mutations/Mb)",
    color = "firebrick",
    size = 2.8,
    fontface = "bold",
    hjust = 0
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "TCGA_TMB_Distribution_Boxplot.png"
  ),
  tmb_plot,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

gep_plot <- plot_tcga_score(
  expression_data,
  "GEP",
  "Immune Gene Expression Profile (GEP) Score",
  "GEP Distribution Across TCGA Cancers",
  "TCGA_GEP_Distribution_Boxplot"
)

tgfb_plot <- plot_tcga_score(
  expression_data,
  "TGFB",
  "TGF-β Signaling Pathway Expression Score",
  "TGF-β Distribution Across TCGA Cancers",
  "TCGA_TGFB_Distribution_Boxplot"
)

# ==================================================================
# 二、Limma火山图
# ==================================================================

limma_results <- readRDS(
  file.path(
    RDATA_DIR,
    "limma_deltabeta_full_results.rds"
  )
)

limma_data <- limma_results %>%
  mutate(
    log_p = -log10(
      pmax(
        adj.P.Val,
        .Machine$double.xmin
      )
    ),
    significance = case_when(
      adj.P.Val < 0.05 &
        Delta_Beta > 0.2 ~ "Hypermethylated",
      adj.P.Val < 0.05 &
        Delta_Beta < -0.2 ~ "Hypomethylated",
      TRUE ~ "Not Significant"
    )
  )

volcano_plot <- ggplot(
  limma_data,
  aes(
    x = Delta_Beta,
    y = log_p,
    color = significance
  )
) +
  geom_point(
    alpha = 0.5,
    size = 1.2
  ) +
  geom_vline(
    xintercept = c(-0.2, 0.2),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
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
    legend.position = "bottom"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "limma_DeltaBeta_Volcano_Plot.pdf"
  ),
  volcano_plot,
  width = 10,
  height = 8
)

ggsave(
  file.path(
    FIGURE_DIR,
    "limma_DeltaBeta_Volcano_Plot.png"
  ),
  volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ==================================================================
# 三、特征筛选方法图
# ==================================================================

plot_importance <- function(data,
                            id_column,
                            value_column,
                            title,
                            x_label,
                            output_file,
                            fill_color) {
  data <- data %>%
    arrange(
      desc(abs(.data[[value_column]]))
    ) %>%
    slice_head(n = 15)

  plot <- ggplot(
    data,
    aes(
      x = reorder(
        .data[[id_column]],
        abs(.data[[value_column]])
      ),
      y = abs(.data[[value_column]])
    )
  ) +
    geom_col(
      fill = fill_color
    ) +
    coord_flip() +
    labs(
      title = title,
      x = x_label,
      y = "Importance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )

  ggsave(
    file.path(
      FIGURE_DIR,
      output_file
    ),
    plot,
    width = 10,
    height = 8
  )
}

importance_specs <- list(
  list(
    file = "LASSO_Feature_Importance.rds",
    id = "Probe_ID",
    value = "LASSO_Coefficient",
    title = "Top 15 Probes by LASSO Importance",
    output = "LASSO_Feature_Importance_Top15.pdf",
    color = "#1f77b4"
  ),
  list(
    file = "RF_Feature_Importance.rds",
    id = "Probe_ID",
    value = "MeanDecreaseGini",
    title = "Top 15 Probes by Random Forest Importance",
    output = "RF_Feature_Importance_Top15.pdf",
    color = "#2ca02c"
  ),
  list(
    file = "XGBoost_Feature_Importance.rds",
    id = "Probe_ID",
    value = "Gain",
    title = "Top 15 Probes by XGBoost Importance",
    output = "XGBoost_Feature_Importance_Top15.pdf",
    color = "#ff7f0e"
  ),
  list(
    file = "ElasticNet_Feature_Importance.rds",
    id = "Probe_ID",
    value = "Coefficient",
    title = "Top 15 Probes by Elastic Net Importance",
    output = "ElasticNet_Feature_Importance_Top15.pdf",
    color = "#FF6347"
  ),
  list(
    file = "RFE_Feature_Importance.rds",
    id = "Probe_ID",
    value = "MeanDecreaseGini",
    title = "Top 15 Probes by RFE Importance",
    output = "RFE_Feature_Importance_Top15.pdf",
    color = "#9b59b6"
  ),
  list(
    file = "Boruta_Feature_Importance.rds",
    id = "Probe_ID",
    value = "MeanDecreaseAccuracy",
    title = "Top 15 Probes by Boruta Importance",
    output = "Boruta_Feature_Importance_Top15.pdf",
    color = "#F0E442"
  ),
  list(
    file = "LightGBM_Feature_Importance.rds",
    id = "Probe_ID",
    value = "Gain_Importance",
    title = "Top 15 Probes by LightGBM Importance",
    output = "LightGBM_Feature_Importance_Top15.pdf",
    color = "#00CED1"
  )
)

for (spec in importance_specs) {
  importance_file <- file.path(
    RDATA_DIR,
    spec$file
  )

  if (!file.exists(importance_file)) {
    next
  }

  importance_data <- readRDS(importance_file)

  if (!all(
    c(spec$id, spec$value) %in%
      colnames(importance_data)
  )) {
    next
  }

  plot_importance(
    importance_data,
    spec$id,
    spec$value,
    spec$title,
    "Probe ID",
    spec$output,
    spec$color
  )
}

if (file.exists(
  file.path(
    RDATA_DIR,
    "LASSO_CV_Objects.RData"
  )
)) {
  load(
    file.path(
      RDATA_DIR,
      "LASSO_CV_Objects.RData"
    )
  )

  pdf(
    file.path(
      FIGURE_DIR,
      "LASSO_CV_Plot.pdf"
    ),
    width = 8,
    height = 6
  )

  plot(lasso_cv)
  dev.off()
}

if (file.exists(
  file.path(
    RDATA_DIR,
    "RFE_Result.rds"
  )
)) {
  rfe_result <- readRDS(
    file.path(
      RDATA_DIR,
      "RFE_Result.rds"
    )
  )

  if (!is.null(rfe_result)) {
    rfe_performance <- rfe_result$results

    rfe_plot <- ggplot(
      rfe_performance,
      aes(
        x = Variables,
        y = ROC
      )
    ) +
      geom_point(
        size = 3,
        color = "#1f77b4"
      ) +
      geom_line(
        color = "#1f77b4",
        linewidth = 0.8
      ) +
      labs(
        title = "RFE Feature Selection Performance",
        x = "Number of Variables",
        y = "ROC"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )

    ggsave(
      file.path(
        FIGURE_DIR,
        "RFE_Performance_Plot.pdf"
      ),
      rfe_plot,
      width = 10,
      height = 6
    )
  }
}

# ==================================================================
# 四、特征筛选方法性能比较
# ==================================================================

feature_evaluation_file <- file.path(
  TABLE_DIR,
  "Feature_Methods_Evaluation_Results.csv"
)

if (file.exists(feature_evaluation_file)) {
  feature_eval <- read.csv(
    feature_evaluation_file,
    stringsAsFactors = FALSE
  )

  feature_eval_long <- feature_eval %>%
    pivot_longer(
      cols = c(
        Val_AUC,
        Val_Accuracy,
        Val_Sensitivity,
        Val_Specificity,
        Val_Precision,
        Val_F1
      ),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = sub(
        "^Val_",
        "",
        Metric
      )
    )

  auc_order <- feature_eval %>%
    arrange(desc(Val_AUC)) %>%
    pull(Method)

  feature_eval_long$Method <- factor(
    feature_eval_long$Method,
    levels = auc_order
  )

  performance_plot <- ggplot(
    feature_eval_long,
    aes(
      x = Method,
      y = Value,
      fill = Metric
    )
  ) +
    geom_col(
      position = "dodge",
      alpha = 0.8
    ) +
    labs(
      title = "Performance Comparison of Feature Selection Methods",
      x = "Feature Selection Method",
      y = "Performance",
      fill = "Metric"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )

  ggsave(
    file.path(
      FIGURE_DIR,
      "Feature_Methods_Performance_Plot.pdf"
    ),
    performance_plot,
    width = 10,
    height = 6
  )

  radar_metric <- feature_eval %>%
    select(
      Method,
      Val_AUC,
      Val_Accuracy,
      Val_Sensitivity,
      Val_Specificity,
      Val_Precision,
      Val_F1
    ) %>%
    mutate(
      across(
        -Method,
        ~ scales::rescale(.x, to = c(0, 1))
      )
    )

  radar_matrix <- radar_metric %>%
    tibble::column_to_rownames("Method")

  radar_data <- rbind(
    rep(1, ncol(radar_matrix)),
    rep(0, ncol(radar_matrix)),
    radar_matrix
  )

  radar_colors <- brewer.pal(
    max(3, nrow(radar_matrix)),
    "Paired"
  )[seq_len(nrow(radar_matrix))]

  pdf(
    file.path(
      FIGURE_DIR,
      "Feature_Selection_Radar_Plot.pdf"
    ),
    width = 9,
    height = 7
  )

  par(
    mar = c(1, 2, 3, 5)
  )

  radarchart(
    radar_data,
    axistype = 1,
    maxmin = TRUE,
    pcol = radar_colors,
    plwd = 2,
    plty = 1,
    cglcol = "gray70",
    cglty = 2,
    vlcex = 1,
    vlabels = c(
      "AUC",
      "Accuracy",
      "Sensitivity",
      "Specificity",
      "Precision",
      "F1"
    ),
    title = "Feature Selection Performance"
  )

  legend(
    "bottomright",
    legend = rownames(radar_matrix),
    col = radar_colors,
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.8
  )

  dev.off()
}

# ==================================================================
# 五、机器学习模型ROC曲线
# ==================================================================

load(
  file.path(
    RDATA_DIR,
    "Validation_Evaluation_Results.RData"
  )
)

load(
  file.path(
    RDATA_DIR,
    "Test_Evaluation_Results.RData"
  )
)

val_pred_all <- rbindlist(
  val_pred_list,
  idcol = "model"
) %>%
  mutate(
    dataset = "Validation Set"
  )

test_pred_all <- rbindlist(
  test_pred_list,
  idcol = "model"
) %>%
  mutate(
    dataset = "Test Set"
  )

model_order <- c(
  "ElasticNet",
  "KNN",
  "LASSO",
  "LightGBM",
  "NeuralNetwork",
  "RandomForest",
  "SVM",
  "Xgboost"
)

model_order <- model_order[
  model_order %in% unique(val_pred_all$model) &
    model_order %in% unique(test_pred_all$model)
]

make_roc_list <- function(pred_data, models) {
  roc_list <- lapply(
    models,
    function(model_name) {
      model_data <- pred_data[
        pred_data$model == model_name,
        ,
        drop = FALSE
      ]

      pROC::roc(
        model_data$truth,
        model_data$prob.1,
        levels = c("0", "1"),
        direction = "<",
        name = model_name,
        quiet = TRUE
      )
    }
  )

  names(roc_list) <- models
  roc_list
}

val_roc_objects <- make_roc_list(
  val_pred_all,
  model_order
)

test_roc_objects <- make_roc_list(
  test_pred_all,
  model_order
)

paired_colors <- brewer.pal(
  length(model_order),
  "Paired"
)

names(paired_colors) <- model_order

roc_labels <- function(roc_objects) {
  vapply(
    names(roc_objects),
    function(model_name) {
      roc_object <- roc_objects[[model_name]]
      auc_value <- as.numeric(auc(roc_object))
      auc_ci <- ci.auc(roc_object)

      sprintf(
        "%s (AUC=%.3f, 95%%CI: %.3f-%.3f)",
        model_name,
        auc_value,
        auc_ci[1],
        auc_ci[3]
      )
    },
    character(1)
  )
}

val_roc_plot <- ggroc(
  val_roc_objects,
  legacy.axes = TRUE,
  linewidth = 1
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_manual(
    values = paired_colors,
    labels = roc_labels(val_roc_objects)
  ) +
  labs(
    title = "ROC Curve (Validation Set)",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_blank()
  )

test_roc_plot <- ggroc(
  test_roc_objects,
  legacy.axes = TRUE,
  linewidth = 1
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_manual(
    values = paired_colors,
    labels = roc_labels(test_roc_objects)
  ) +
  labs(
    title = "ROC Curve (Test Set)",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_blank()
  )

roc_combined <- val_roc_plot +
  test_roc_plot +
  plot_annotation(
    title = "ROC Curves of All Models"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "All_Models_ROC_Combined.pdf"
  ),
  roc_combined,
  width = 14,
  height = 6
)

ggsave(
  file.path(
    FIGURE_DIR,
    "All_Models_ROC_Combined.png"
  ),
  roc_combined,
  width = 14,
  height = 6,
  dpi = 300
)

# ==================================================================
# 六、校准曲线
# ==================================================================

calibration_data <- function(pred_data) {
  pred_data %>%
    rename(
      Result = truth,
      prob = prob.1
    ) %>%
    mutate(
      Result = factor(
        Result,
        levels = c("0", "1"),
        labels = c("No", "Yes")
      )
    ) %>%
    select(
      model,
      Result,
      prob,
      dataset
    )
}

calib_val <- calibration_data(val_pred_all)
calib_test <- calibration_data(test_pred_all)

plot_calibration <- function(data, title) {
  output_list <- list()

  for (model_name in unique(data$model)) {
    model_data <- data[
      data$model == model_name,
      ,
      drop = FALSE
    ]

    calibration_object <- caret::calibration(
      Result ~ prob,
      data = model_data,
      class = "Yes",
      cuts = 4
    )

    calibration_df <- as.data.frame(
      calibration_object$data
    ) %>%
      na.omit() %>%
      mutate(
        model = model_name
      )

    output_list[[model_name]] <- calibration_df
  }

  plot_data <- bind_rows(output_list)

  ggplot(
    plot_data,
    aes(
      x = midpoint,
      y = Percent,
      color = model,
      group = model
    )
  ) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotdash"
    ) +
    labs(
      title = title,
      x = "Predicted Probability",
      y = "Observed Event Percentage",
      color = "Model"
    ) +
    scale_color_brewer(
      palette = "Paired"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
}

val_calibration_plot <- plot_calibration(
  calib_val,
  "Calibration Curve (Validation Set)"
)

test_calibration_plot <- plot_calibration(
  calib_test,
  "Calibration Curve (Test Set)"
)

calibration_combined <- val_calibration_plot +
  test_calibration_plot +
  plot_annotation(
    title = "Calibration Curves of All Models"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "All_Models_Calibration_Val_Test_Combined.pdf"
  ),
  calibration_combined,
  width = 14,
  height = 6
)

ggsave(
  file.path(
    FIGURE_DIR,
    "All_Models_Calibration_Val_Test_Combined.png"
  ),
  calibration_combined,
  width = 14,
  height = 6,
  dpi = 300
)

# ==================================================================
# 七、Platt校准和Brier Score
# ==================================================================

calibrate_platt <- function(val_data, test_data) {
  result_list <- list()

  for (model_name in unique(val_data$model)) {
    val_model <- val_data %>%
      filter(model == model_name) %>%
      mutate(
        y = ifelse(Result == "Yes", 1, 0),
        prob_clipped = pmin(
          pmax(prob, 1e-6),
          1 - 1e-6
        )
      )

    test_model <- test_data %>%
      filter(model == model_name) %>%
      mutate(
        y = ifelse(Result == "Yes", 1, 0),
        prob_clipped = pmin(
          pmax(prob, 1e-6),
          1 - 1e-6
        )
      )

    platt_model <- glm(
      y ~ prob_clipped,
      data = val_model,
      family = binomial()
    )

    test_model$prob <- predict(
      platt_model,
      newdata = test_model,
      type = "response"
    )

    result_list[[model_name]] <- test_model %>%
      select(
        model,
        Result,
        prob,
        dataset
      )
  }

  bind_rows(result_list)
}

calib_test_platt <- calibrate_platt(
  calib_val,
  calib_test
)

test_calibration_platt <- plot_calibration(
  calib_test_platt,
  "Test Set Calibration after Platt Scaling"
)

platt_comparison <- test_calibration_plot +
  test_calibration_platt +
  plot_annotation(
    title = "Test Set Calibration: Before vs After Platt Scaling"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "Test_Calibration_Before_After_Platt.pdf"
  ),
  platt_comparison,
  width = 14,
  height = 6
)

calc_brier <- function(data) {
  data %>%
    mutate(
      y = ifelse(Result == "Yes", 1, 0)
    ) %>%
    group_by(
      model,
      dataset
    ) %>%
    summarise(
      Brier_Score = mean(
        (prob - y)^2,
        na.rm = TRUE
      ),
      n = n(),
      .groups = "drop"
    )
}

brier_before <- calc_brier(
  bind_rows(calib_val, calib_test)
)

brier_after <- calc_brier(
  calib_test_platt
)

write.csv(
  brier_before,
  file.path(
    TABLE_DIR,
    "Brier_Score_Before_Calibration.csv"
  ),
  row.names = FALSE
)

write.csv(
  brier_after,
  file.path(
    TABLE_DIR,
    "Brier_Score_After_Platt.csv"
  ),
  row.names = FALSE
)

brier_comparison <- brier_before %>%
  filter(dataset == "Test Set") %>%
  select(
    model,
    Brier_Before = Brier_Score
  ) %>%
  inner_join(
    brier_after %>%
      select(
        model,
        Brier_After = Brier_Score
      ),
    by = "model"
  ) %>%
  mutate(
    Improvement = Brier_Before - Brier_After
  )

write.csv(
  brier_comparison,
  file.path(
    TABLE_DIR,
    "Brier_Score_Comparison.csv"
  ),
  row.names = FALSE
)

# ==================================================================
# 八、机器学习DCA和预测概率分布
# ==================================================================

make_dca_data <- function(pred_list, dataset_name) {
  lapply(
    names(pred_list),
    function(model_name) {
      pred_list[[model_name]] %>%
        mutate(
          Result = as.integer(
            as.character(truth)
          ),
          dataset = dataset_name,
          model = model_name
        ) %>%
        select(
          model,
          Result,
          prob.1,
          dataset
        )
    }
  ) %>%
    rbindlist()
}

val_dca_data <- make_dca_data(
  val_pred_list,
  "Validation Set"
)

test_dca_data <- make_dca_data(
  test_pred_list,
  "Test Set"
)

plot_dca <- function(data, dataset_name) {
  models <- model_order[
    model_order %in% unique(data$model)
  ]

  dca_objects <- list()

  for (model_name in models) {
    model_data <- data[
      data$model == model_name,
      ,
      drop = FALSE
    ]

    dca_objects[[model_name]] <- decision_curve(
      Result ~ prob.1,
      data = model_data,
      study.design = "cohort",
      bootstraps = 500,
      confidence.intervals = "none",
      thresholds = seq(0, 1, by = 0.01)
    )
  }

  plot_decision_curve(
    dca_objects,
    curve.names = models,
    cost.benefit.axis = FALSE,
    lwd = 2,
    col = brewer.pal(
      length(models),
      "Paired"
    ),
    legend.position = "topright",
    main = paste(
      "Decision Curve Analysis",
      dataset_name
    )
  )
}

pdf(
  file.path(
    FIGURE_DIR,
    "All_Models_DCA_Combined.pdf"
  ),
  width = 10,
  height = 14
)

layout(
  matrix(
    c(1, 2),
    nrow = 2
  )
)

plot_dca(
  val_dca_data,
  "Validation Set"
)

plot_dca(
  test_dca_data,
  "Test Set"
)

dev.off()
