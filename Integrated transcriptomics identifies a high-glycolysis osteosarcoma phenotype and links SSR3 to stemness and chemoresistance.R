# ================================================================
# 完整分析脚本：GGPS 分类器构建与嵌套验证
# 对应手稿：Integrated transcriptomics identifies a high-glycolysis osteosarcoma phenotype...
# 数据来源：codeocean_9535428, GSE152048, GSE162454
# ================================================================

# ================================================================
# 第一部分：环境配置与数据加载
# ================================================================

library(Seurat)
library(randomForest)
library(pROC)
library(glmnet)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(PRROC)
library(reshape2)
library(tidyr)

set.seed(123)  # 全局随机种子

# 加载数据
data <- readRDS("scRNAmalfinal.rds")
cat("原始Seurat对象加载成功，细胞数:", ncol(data), "\n")

# ================================================================
# 第二部分：META-Glycolysis基因集定义（341个基因）
# ================================================================

meta_genes_341 <- c(
  "ACSS2", "GCK", "PGK2", "PGK1", "PDHB", "PDHA1", "PDHA2", "PGM2", "TPI1",
  "ACSS1", "FBP1", "ADH1B", "HK2", "ADH1C", "HK1", "HK3", "ADH4", "PGAM2",
  "ADH5", "PGAM1", "ADH1A", "ALDOC", "ALDH7A1", "LDHAL6B", "PKLR", "LDHAL6A",
  "ENO1", "PKM", "PFKP", "BPGM", "PCK2", "PCK1", "ALDH1B1", "ALDH2", "ALDH3A1",
  "AKR1A1", "FBP2", "PFKM", "PFKL", "LDHC", "GAPDH", "ENO3", "ENO2", "PGAM4",
  "ADH7", "ADH6", "LDHB", "ALDH1A3", "ALDH3B1", "ALDH3B2", "ALDH9A1", "ALDH3A2",
  "GALM", "ALDOA", "DLD", "DLAT", "ALDOB", "G6PC2", "LDHA", "G6PC1", "PGM1",
  "GPI", "ACTN3", "ADPGK", "APP", "ARL2", "ARNT", "CBFA2T3", "DDIT4", "DHTKD1",
  "EIF6", "ENO4", "EP300", "FLCN", "FOXK1", "FOXK2", "GALK1", "GAPDHS", "GIT1",
  "GPD1", "HDAC4", "HIF1A", "HKDC1", "HTR2A", "IER3", "IFNG", "IGF1", "INS",
  "INSR", "JMJD8", "KAT2B", "LIPA", "MLST8", "MLXIPL", "MTOR", "NCOR1", "NUPR1",
  "OGDH", "OGDHL", "OGT", "P2RX7", "PFKFB1", "PFKFB2", "PPARA", "PPP2CA",
  "PRKAA1", "PRKAA2", "PRKACA", "PRKAG1", "PRKAG2", "PRKAG3", "PRXL2C", "PSEN1",
  "RPTOR", "SIRT6", "SLC2A6", "SLC4A1", "SLC4A4", "STAT3", "TIGAR", "TREX1",
  "UCP2", "ZBTB20", "ZBTB7A", "AAAS", "GCKR", "GNPDA1", "GNPDA2", "NDC1",
  "NUP107", "NUP133", "NUP153", "NUP155", "NUP160", "NUP188", "NUP205", "NUP210",
  "NUP214", "NUP35", "NUP37", "NUP42", "NUP43", "NUP50", "NUP54", "NUP58",
  "NUP62", "NUP85", "NUP88", "NUP93", "NUP98", "PFKFB3", "PFKFB4", "PGM2L1",
  "PGP", "POM121", "POM121C", "PPP2CB", "PPP2R1A", "PPP2R1B", "PPP2R5D",
  "PRKACB", "PRKACG", "RAE1", "RANBP2", "SEC13", "SEH1L", "TPR", "ABCB6",
  "ADORA2B", "AGL", "AGRN", "AK3", "AK4", "ALG1", "ANG", "ANGPTL4", "ANKZF1",
  "ARPP19", "ARTN", "AURKA", "B3GALT6", "B3GAT1", "B3GAT3", "B3GNT3", "B4GALT1",
  "B4GALT2", "B4GALT4", "B4GALT7", "BIK", "BPNT1", "CACNA1H", "CAPN5", "CASP6",
  "CD44", "CDK1", "CENPA", "CHPF", "CHPF2", "CHST1", "CHST12", "CHST2", "CHST4",
  "CHST6", "CITED2", "CLDN3", "CLDN9", "CLN6", "COG2", "COL5A1", "COPB2",
  "CTH", "CXCR4", "CYB5A", "DCN", "DEPDC1", "DPYSL4", "DSC2", "ECD", "EFNA3",
  "EGFR", "EGLN3", "ELF3", "ERO1A", "EXT1", "EXT2", "FAM162A", "FKBP4", "FUT8",
  "G6PD", "GAL3ST1", "GALE", "GALK2", "GCLC", "GFPT1", "GFUS", "GLCE", "GLRX",
  "GMPPA", "GMPPB", "GNE", "GOT1", "GOT2", "GPC1", "GPC3", "GPC4", "GPR87",
  "GUSB", "GYS1", "GYS2", "HAX1", "HDLBP", "HMMR", "HOMER1", "HS2ST1", "HS6ST2",
  "HSPA5", "IDH1", "IDUA", "IGFBP3", "IL13RA1", "IRS2", "ISG20", "KDELR3",
  "KIF20A", "KIF2A", "LCT", "LHPP", "LHX9", "MDH1", "MDH2", "ME1", "ME2",
  "MED24", "MERTK", "MET", "MIF", "MIOX", "MPI", "MXI1", "NANP", "NASP",
  "NDST3", "NDUFV3", "NOL3", "NSDHL", "NT5E", "P4HA1", "P4HA2", "PAM",
  "PAXIP1", "PC", "PDK3", "PGLS", "PHKA2", "PKP2", "PLOD1", "PLOD2", "PMM2",
  "POLR3K", "PPFIA4", "PPIA", "PRPS1", "PSMC4", "PYGB", "PYGL", "QSOX1",
  "RARS1", "RBCK1", "RPE", "RRAGD", "SAP30", "SDC1", "SDC2", "SDC3", "SDHC",
  "SLC16A3", "SLC25A10", "SLC25A13", "SLC35A3", "SLC37A4", "SOD1", "SOX9",
  "SPAG4", "SRD5A3", "STC1", "STC2", "STMN1", "TALDO1", "TFF3", "TGFA",
  "TGFBI", "TKTL1", "TPBG", "TPST1", "TXN", "UGP2", "VCAN", "VEGFA", "VLDLR",
  "XYLT2", "ZNF292"
)

# ================================================================
# 第三部分：构建表达矩阵和基础数据框
# ================================================================

expr_matrix <- GetAssayData(data, layer = "data")
meta_genes_present <- intersect(meta_genes_341, rownames(expr_matrix))
cat("实际存在的META基因数:", length(meta_genes_present), "\n")

expr_sub <- expr_matrix[meta_genes_present, ]
expr_df <- as.data.frame(t(as.matrix(expr_sub)))
expr_df$group <- data@meta.data$group
expr_df$dataset <- data@meta.data$dataset

# 只保留 H 和 L
expr_df <- expr_df %>% filter(group %in% c("HGlyOS", "LGlyOS"))
expr_df$group <- factor(expr_df$group, levels = c("LGlyOS", "HGlyOS"))
cat("总细胞数（H+L）:", nrow(expr_df), "\n")

# ================================================================
# 第四部分：CGAS阈值比较（训练集 vs 全局）
# ================================================================

# 检查 CGAS 列
if ("CGAS" %in% colnames(data@meta.data)) {
  expr_df$CGAS <- data@meta.data[rownames(expr_df), "CGAS"]
} else if ("Scoring" %in% colnames(data@meta.data)) {
  expr_df$CGAS <- data@meta.data[rownames(expr_df), "Scoring"]
} else {
  stop("未找到 CGAS 或 Scoring 列")
}

datasets <- c("codeocean_9535428", "GSE152048", "GSE162454")
thresholds_global <- quantile(expr_df$CGAS, c(0.1, 0.9), na.rm = TRUE)

cat("\n========== CGAS阈值对比 ==========\n")
cat("全局阈值: [", round(thresholds_global[1], 4), ", ", round(thresholds_global[2], 4), "]\n")

for (test_ds in datasets) {
  train_ds <- setdiff(datasets, test_ds)
  train_data <- expr_df %>% filter(dataset %in% train_ds)
  test_data  <- expr_df %>% filter(dataset == test_ds)
  
  thresholds_train <- quantile(train_data$CGAS, c(0.1, 0.9), na.rm = TRUE)
  
  cat("\n测试集:", test_ds, "\n")
  cat("  训练集阈值: [", round(thresholds_train[1], 4), ", ", round(thresholds_train[2], 4), "]\n")
  cat("  差值: [", round(thresholds_train[1] - thresholds_global[1], 4), ", ", 
      round(thresholds_train[2] - thresholds_global[2], 4), "]\n")
}

# ================================================================
# 第五部分：嵌套特征选择（核心分析）
# ================================================================

# 参数设置
n_folds <- 10           # LASSO内部CV折数
ntree_rf <- 500        # 随机森林树的数量
classification_threshold <- 0.5  # 二分类阈值

nested_true_list <- list()
nested_pred_list <- list()
nested_patient_list <- list()
selected_genes_list <- list()
results_auc <- list()

for (test_ds in datasets) {
  train_ds <- setdiff(datasets, test_ds)
  
  train_data <- expr_df %>% filter(dataset %in% train_ds)
  test_data  <- expr_df %>% filter(dataset == test_ds)
  
  cat("\n========== 测试集:", test_ds, "==========\n")
  cat("训练集细胞数:", nrow(train_data), " (H:", sum(train_data$group == "HGlyOS"), ", L:", sum(train_data$group == "LGlyOS"), ")\n")
  cat("测试集细胞数:", nrow(test_data), " (H:", sum(test_data$group == "HGlyOS"), ", L:", sum(test_data$group == "LGlyOS"), ")\n")
  
  # 获取测试集患者ID
  cell_ids <- rownames(test_data)
  patient_ids <- data@meta.data[cell_ids, "orig.ident"]
  cat("测试集患者数:", length(unique(patient_ids)), "\n")
  
  # ---------- 步骤1：训练集内部差异表达分析 ----------
  deg_results <- data.frame(gene = meta_genes_present, p_val = NA)
  for (gene in meta_genes_present) {
    expr_h <- train_data[[gene]][train_data$group == "HGlyOS"]
    expr_l <- train_data[[gene]][train_data$group == "LGlyOS"]
    if (length(expr_h) > 0 && length(expr_l) > 0) {
      deg_results$p_val[deg_results$gene == gene] <- wilcox.test(expr_h, expr_l)$p.value
    }
  }
  deg_filtered <- deg_results %>% filter(!is.na(p_val)) %>% arrange(p_val) %>% head(100)
  
  # ---------- 步骤2：LASSO筛选（lambda.min规则） ----------
  if (nrow(deg_filtered) > 10) {
    x_train <- as.matrix(train_data[, deg_filtered$gene])
    y_train <- as.numeric(train_data$group == "HGlyOS")
    set.seed(123)
    cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = n_folds)
    # lambda.min 选择规则（最小交叉验证误差）
    coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
    lasso_genes <- rownames(coefs)[coefs[, 1] != 0][-1]
    if (length(lasso_genes) < 5) lasso_genes <- deg_filtered$gene[1:20]
    fold_genes <- lasso_genes
  } else {
    fold_genes <- deg_filtered$gene[1:20]
  }
  
  selected_genes_list[[test_ds]] <- fold_genes
  cat("本折筛选基因数:", length(fold_genes), "\n")
  
  # ---------- 步骤3：随机森林训练（参数详见上方） ----------
  set.seed(123)
  rf <- randomForest(
    x = train_data[, fold_genes],
    y = train_data$group,
    ntree = ntree_rf,
    mtry = floor(sqrt(length(fold_genes))),  # RF默认参数
    importance = TRUE
  )
  
  # ---------- 步骤4：测试集预测 ----------
  pred_prob <- predict(rf, newdata = test_data[, fold_genes], type = "prob")[, "HGlyOS"]
  # 分类阈值 = 0.5（默认）
  pred_class <- ifelse(pred_prob >= classification_threshold, "HGlyOS", "LGlyOS")
  true_label <- test_data$group
  
  # 保存结果
  nested_true_list[[test_ds]] <- true_label
  nested_pred_list[[test_ds]] <- pred_prob
  nested_patient_list[[test_ds]] <- patient_ids
  
  auc_val <- roc(true_label, pred_prob, levels = c("LGlyOS", "HGlyOS"), direction = "<")$auc
  results_auc[[test_ds]] <- auc_val
  cat("AUC =", round(auc_val, 4), "\n")
}

# ================================================================
# 第六部分：共识基因与性能指标
# ================================================================

common_genes <- Reduce(intersect, selected_genes_list)
cat("\n共识基因数:", length(common_genes), "\n")
cat("共识基因:", paste(common_genes, collapse=", "), "\n")

# 计算完整性能指标
calc_metrics <- function(true_label, pred_prob, pos_class = "HGlyOS", threshold = 0.5) {
  roc_obj <- roc(true_label, pred_prob, levels = c("LGlyOS", "HGlyOS"), direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  y_binary <- as.numeric(true_label == pos_class)
  pr_obj <- pr.curve(scores.class0 = pred_prob[y_binary == 1], 
                     scores.class1 = pred_prob[y_binary == 0], 
                     curve = FALSE)
  auprc_val <- pr_obj$auc.integral
  
  pred_class <- ifelse(pred_prob >= threshold, "HGlyOS", "LGlyOS")
  cm <- table(Predicted = pred_class, Actual = true_label)
  
  tp <- ifelse("HGlyOS" %in% rownames(cm) & "HGlyOS" %in% colnames(cm), cm["HGlyOS", "HGlyOS"], 0)
  fp <- ifelse("HGlyOS" %in% rownames(cm) & "LGlyOS" %in% colnames(cm), cm["HGlyOS", "LGlyOS"], 0)
  fn <- ifelse("LGlyOS" %in% rownames(cm) & "HGlyOS" %in% colnames(cm), cm["LGlyOS", "HGlyOS"], 0)
  tn <- ifelse("LGlyOS" %in% rownames(cm) & "LGlyOS" %in% colnames(cm), cm["LGlyOS", "LGlyOS"], 0)
  
  tp <- as.numeric(tp); fp <- as.numeric(fp); fn <- as.numeric(fn); tn <- as.numeric(tn)
  
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)
  f1 <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
  
  sensitivity <- recall
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), 0)
  bal_acc <- (sensitivity + specificity) / 2
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(denom > 0, (tp * tn - fp * fn) / denom, 0)
  
  return(data.frame(
    AUC = round(auc_val, 4),
    AUPRC = round(auprc_val, 4),
    F1 = round(f1, 4),
    Balanced_Accuracy = round(bal_acc, 4),
    Accuracy = round(accuracy, 4),
    MCC = round(mcc, 4),
    Sensitivity = round(sensitivity, 4),
    Specificity = round(specificity, 4),
    TP = tp, TN = tn, FP = fp, FN = fn
  ))
}

metrics_list <- list()
for (ds in names(nested_true_list)) {
  metrics_list[[ds]] <- calc_metrics(nested_true_list[[ds]], nested_pred_list[[ds]])
}
metrics_df <- do.call(rbind, metrics_list)
metrics_df$Dataset <- names(nested_true_list)
metrics_df <- metrics_df[, c("Dataset", setdiff(colnames(metrics_df), "Dataset"))]
cat("\n========== 完整性能指标 ==========\n")
print(metrics_df, digits = 4)

# ================================================================
# 第七部分：患者级AUC计算
# ================================================================

per_patient_auc <- function(true_label, pred_prob, patient_ids) {
  patients <- unique(patient_ids)
  results <- list()
  for (p in patients) {
    idx <- patient_ids == p
    if (length(unique(true_label[idx])) < 2) next
    roc_obj <- roc(true_label[idx], pred_prob[idx], 
                   levels = c("LGlyOS", "HGlyOS"), 
                   direction = "<", 
                   quiet = TRUE)
    results[[p]] <- as.numeric(auc(roc_obj))
  }
  return(data.frame(Patient = names(results), AUC = unlist(results)))
}

all_patient_results <- list()
for (ds in names(nested_true_list)) {
  df <- per_patient_auc(nested_true_list[[ds]], nested_pred_list[[ds]], nested_patient_list[[ds]])
  all_patient_results[[ds]] <- df
}
combined_df <- do.call(rbind, all_patient_results)
combined_df$Dataset <- rep(names(all_patient_results), 
                           times = sapply(all_patient_results, nrow))

cat("\n========== 患者级AUC汇总 ==========\n")
print(combined_df)

# ================================================================
# 第八部分：患者聚类Bootstrap置信区间
# ================================================================

n_boot <- 2000
set.seed(123)
boot_auc_results <- list()

for (ds in names(nested_true_list)) {
  true_label <- nested_true_list[[ds]]
  pred_prob <- nested_pred_list[[ds]]
  patient_ids <- nested_patient_list[[ds]]
  unique_patients <- unique(patient_ids)
  
  boot_auc <- replicate(n_boot, {
    boot_patients <- sample(unique_patients, replace = TRUE)
    boot_cells <- patient_ids %in% boot_patients
    if (sum(boot_cells) < 5 || length(unique(true_label[boot_cells])) < 2) return(NA)
    roc_obj <- tryCatch(
      roc(true_label[boot_cells], pred_prob[boot_cells], 
          levels = c("LGlyOS", "HGlyOS"), direction = "<", quiet = TRUE),
      error = function(e) return(NULL)
    )
    if (is.null(roc_obj)) return(NA)
    as.numeric(auc(roc_obj))
  })
  
  boot_auc <- boot_auc[!is.na(boot_auc)]
  if (length(boot_auc) > 50) {
    ci <- quantile(boot_auc, c(0.025, 0.975))
    roc_orig <- roc(true_label, pred_prob, levels = c("LGlyOS", "HGlyOS"), direction = "<", quiet = TRUE)
    auc_orig <- as.numeric(auc(roc_orig))
    boot_auc_results[[ds]] <- data.frame(
      Dataset = ds,
      AUC = round(auc_orig, 4),
      CI_lower = round(ci[1], 4),
      CI_upper = round(ci[2], 4),
      n_patients = length(unique_patients),
      n_cells = length(true_label),
      n_boot_valid = length(boot_auc)
    )
  }
}

boot_ci_df <- do.call(rbind, boot_auc_results)
cat("\n========== Bootstrap置信区间 ==========\n")
print(boot_ci_df)

# ================================================================
# 第九部分：保存结果
# ================================================================

saveRDS(nested_true_list, "nested_true_list.rds")
saveRDS(nested_pred_list, "nested_pred_list.rds")
saveRDS(nested_patient_list, "nested_patient_list.rds")
saveRDS(selected_genes_list, "selected_genes_list.rds")
saveRDS(results_auc, "results_auc.rds")
saveRDS(metrics_df, "metrics_df.rds")
saveRDS(boot_ci_df, "boot_ci_df.rds")
saveRDS(combined_df, "patient_auc_df.rds")

cat("\n所有结果已保存完成！\n")