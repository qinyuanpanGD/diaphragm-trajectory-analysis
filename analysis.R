#!/usr/bin/env Rscript
# Consensus k-means of DE, DTee and DTF trajectories on days 1, 4 and 7.

packages <- c("dplyr", "ggplot2", "readr", "tidyr", "pheatmap", "patchwork")
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analysis.R INPUT.csv OUTPUT_DIR")

input_file <- args[1L]

output_dir <- args[2L]

if (!file.exists(input_file)) {
  stop("未找到输入文件：", input_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

heatmap_dir <- file.path(output_dir, "consensus_heatmaps")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

original_used_csv          <- file.path(output_dir, "original_data_used_long.csv")
wide_raw_csv               <- file.path(output_dir, "patient_wide_features_raw.csv")
wide_scaled_csv            <- file.path(output_dir, "patient_wide_features_scaled.csv")
cluster_assignment_csv     <- file.path(output_dir, "cluster_assignments_all_k.csv")
trajectory_summary_csv     <- file.path(output_dir, "trajectory_summary_selectedK.csv")
diagnostics_csv            <- file.path(output_dir, "consensus_diagnostics.csv")

fit_statistics_csv         <- file.path(output_dir, "model_fit_statistics_AIC_BIC_ICL_entropy.csv")
fit_statistics_pdf         <- file.path(output_dir, "model_fit_statistics_AIC_BIC_ICL.pdf")
fit_statistics_png         <- file.path(output_dir, "model_fit_statistics_AIC_BIC_ICL.png")

selection_score_csv        <- file.path(output_dir, "auto_k_selection_score.csv")

rdata_file                 <- file.path(output_dir, "clustering_result.RData")
summary_txt                <- file.path(output_dir, "clustering_summary.txt")

combined_heatmap_pdf       <- file.path(output_dir, "consensus_heatmaps_all_k.pdf")

delta_area_pdf             <- file.path(output_dir, "delta_area.pdf")
delta_area_png             <- file.path(output_dir, "delta_area.png")

cdf_pdf                    <- file.path(output_dir, "consensus_CDF.pdf")
cdf_png                    <- file.path(output_dir, "consensus_CDF.png")

cluster_consensus_pdf      <- file.path(output_dir, "cluster_consensus_barplot.pdf")
cluster_consensus_png      <- file.path(output_dir, "cluster_consensus_barplot.png")

model_selection_pdf        <- file.path(output_dir, "model_selection_combined.pdf")
model_selection_png        <- file.path(output_dir, "model_selection_combined.png")

trajectory_pdf             <- file.path(output_dir, "trajectory_plots_selectedK.pdf")
trajectory_png             <- file.path(output_dir, "trajectory_plots_selectedK.png")

set.seed(12345)

time_points <- c(1, 4, 7)

target_vars <- c("DE", "DTee", "DTF")

max_k  <- 6
reps   <- 1000
p_item <- 0.8

k_sel <- NA_integer_

require_complete_timepoints <- TRUE

kmeans_nstart   <- 50
kmeans_iter_max <- 100

cluster_colors_fixed <- c(
  "1" = "#c30078",
  "2" = "#4370b4",
  "3" = "#549f9a"
)

make_cluster_palette <- function(cluster_levels) {
  cluster_levels <- as.character(cluster_levels)
  cluster_levels <- cluster_levels[!is.na(cluster_levels)]
  cluster_levels <- unique(cluster_levels)

  suppressWarnings({
    numeric_levels <- as.numeric(cluster_levels)
  })

  if (all(is.finite(numeric_levels))) {
    cluster_levels <- cluster_levels[order(numeric_levels)]
  } else {
    cluster_levels <- sort(cluster_levels)
  }

  pal <- rep(NA_character_, length(cluster_levels))
  names(pal) <- cluster_levels

  fixed_levels <- intersect(names(cluster_colors_fixed), cluster_levels)
  pal[fixed_levels] <- cluster_colors_fixed[fixed_levels]

  missing_levels <- names(pal)[is.na(pal)]

  if (length(missing_levels) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_levels), palette = "Dark 3")
    names(extra_cols) <- missing_levels
    pal[missing_levels] <- extra_cols
  }

  pal
}

theme_bw_nogrid <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

cat("=== 1. 读取数据 ===\n")

raw_data <- readr::read_csv(input_file, show_col_types = FALSE)

cat("输入文件：", input_file, "\n")
cat("原始维度：", nrow(raw_data), "x", ncol(raw_data), "\n")

id_col <- if ("patient_id" %in% names(raw_data)) {
  "patient_id"
} else if ("ID" %in% names(raw_data)) {
  "ID"
} else if ("id" %in% names(raw_data)) {
  "id"
} else {
  NA_character_
}

time_col <- if ("time" %in% names(raw_data)) {
  "time"
} else if ("Time" %in% names(raw_data)) {
  "Time"
} else if ("day" %in% names(raw_data)) {
  "day"
} else if ("Day" %in% names(raw_data)) {
  "Day"
} else {
  NA_character_
}

if (is.na(id_col)) {
  stop("数据中未找到患者 ID 列，请确认是否存在 patient_id、ID 或 id。")
}

if (is.na(time_col)) {
  stop("数据中未找到时间列，请确认是否存在 time、Time、day 或 Day。")
}

var_aliases <- list(
  "DE"   = c("DE", "de"),
  "DTee" = c("DTee", "DTEE", "DT_EE", "DT-EE", "dtee"),
  "DTF"  = c("DTF", "dtf")
)

resolve_var <- function(label, data_names) {
  candidates <- var_aliases[[label]]
  found <- candidates[candidates %in% data_names]

  if (length(found) == 0) {
    return(NA_character_)
  } else {
    return(found[1])
  }
}

source_vars <- sapply(target_vars, resolve_var, data_names = names(raw_data))
names(source_vars) <- target_vars

missing_vars <- names(source_vars)[is.na(source_vars)]

if (length(missing_vars) > 0) {
  stop(
    "数据中缺少以下变量：",
    paste(missing_vars, collapse = ", "),
    "\n请检查 CSV 列名，或在 var_aliases 中添加对应别名。"
  )
}

cat("变量匹配结果：\n")
print(data.frame(
  target_label = target_vars,
  csv_column   = unname(source_vars),
  stringsAsFactors = FALSE
))

safe_vars <- make.names(target_vars, unique = TRUE)

indicator_map <- data.frame(
  indicator_var   = safe_vars,
  indicator_label = target_vars,
  source_column   = unname(source_vars),
  stringsAsFactors = FALSE
)

cat("\n=== 2. 长表数据预处理 ===\n")

dat <- data.frame(
  patient_id = as.character(raw_data[[id_col]]),
  time_num   = readr::parse_number(as.character(raw_data[[time_col]])),
  stringsAsFactors = FALSE
)

for (i in seq_along(target_vars)) {
  dat[[safe_vars[i]]] <- readr::parse_number(as.character(raw_data[[source_vars[i]]]))
}

dat <- dat %>%
  filter(!is.na(patient_id), !is.na(time_num)) %>%
  filter(time_num %in% time_points) %>%
  tidyr::drop_na(all_of(safe_vars)) %>%
  arrange(patient_id, time_num)

dat <- dat %>%
  group_by(patient_id, time_num) %>%
  summarise(
    across(all_of(safe_vars), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

if (nrow(dat) == 0) {
  stop("预处理后无可用数据。")
}

readr::write_excel_csv(dat, original_used_csv)

cat("长表数据行数：", nrow(dat), "\n")
cat("患者数：", n_distinct(dat$patient_id), "\n")
cat("时间点：", paste(sort(unique(dat$time_num)), collapse = ", "), "\n")

cat("\n=== 3. 构建患者级宽表轨迹特征矩阵 ===\n")

dat_wide <- dat %>%
  mutate(time_lab = paste0("D", time_num)) %>%
  tidyr::pivot_wider(
    id_cols = patient_id,
    names_from = time_lab,
    values_from = all_of(safe_vars),
    names_glue = "{.value}_{time_lab}"
  )

feature_vars <- unlist(
  lapply(safe_vars, function(v) {
    paste0(v, "_D", time_points)
  })
)

for (v in setdiff(feature_vars, names(dat_wide))) {
  dat_wide[[v]] <- NA_real_
}

dat_wide <- dat_wide %>%
  select(patient_id, all_of(feature_vars)) %>%
  arrange(patient_id)

missing_count_before <- sum(is.na(dat_wide[feature_vars]))

if (require_complete_timepoints) {
  complete_idx <- complete.cases(dat_wide[feature_vars])
  dropped_n <- sum(!complete_idx)
  dat_wide <- dat_wide[complete_idx, , drop = FALSE]

  cat("要求完整时间点：是\n")
  cat("因缺失任一轨迹特征而删除患者数：", dropped_n, "\n")
} else {
  cat("要求完整时间点：否\n")
  cat("宽表缺失值数量：", missing_count_before, "\n")
  cat("使用列中位数填补宽表缺失值。\n")

  for (v in feature_vars) {
    if (any(is.na(dat_wide[[v]]))) {
      med_v <- median(dat_wide[[v]], na.rm = TRUE)

      if (!is.finite(med_v)) {
        med_v <- 0
      }

      dat_wide[[v]][is.na(dat_wide[[v]])] <- med_v
    }
  }
}

if (nrow(dat_wide) < max_k) {
  stop("宽表纳入患者数少于 max_k，请降低 max_k 或检查数据完整性。")
}

patient_ids <- dat_wide$patient_id
patient_n   <- length(patient_ids)

feature_raw <- as.data.frame(dat_wide[feature_vars])

feature_sd <- sapply(feature_raw, sd, na.rm = TRUE)
zero_feature_vars <- names(feature_sd)[!is.finite(feature_sd) | feature_sd == 0]

if (length(zero_feature_vars) > 0) {
  cat("以下宽表特征方差为 0，将从聚类特征中删除：\n")
  print(zero_feature_vars)
}

keep_feature_vars <- setdiff(feature_vars, zero_feature_vars)

if (length(keep_feature_vars) < 2) {
  stop("可用于聚类的非零方差特征少于 2 个，无法聚类。")
}

feature_scaled_mat <- scale(feature_raw[keep_feature_vars])
feature_scaled_mat <- as.matrix(feature_scaled_mat)
rownames(feature_scaled_mat) <- patient_ids

if (any(!is.finite(feature_scaled_mat))) {
  stop("标准化后的宽表特征存在非有限值，请检查数据。")
}

dat_wide_scaled <- data.frame(
  patient_id = patient_ids,
  feature_scaled_mat,
  check.names = FALSE
)

readr::write_excel_csv(dat_wide, wide_raw_csv)
readr::write_excel_csv(dat_wide_scaled, wide_scaled_csv)

cat("最终纳入患者数：", patient_n, "\n")
cat("原始宽表特征数：", length(feature_vars), "\n")
cat("用于聚类的特征数：", length(keep_feature_vars), "\n")

cat("\n=== 4. 开始 consensus k-means 聚类 ===\n")

k_values <- 2:max_k
k_names  <- paste0("k", k_values)

cat("k 范围：", paste(k_values, collapse = ", "), "\n")
cat("重抽样次数：", reps, "\n")
cat("每次抽样比例：", p_item, "\n")
cat("k-means nstart：", kmeans_nstart, "\n")
cat("k-means iter.max：", kmeans_iter_max, "\n")

make_zero_mat <- function() {
  matrix(
    0,
    nrow = patient_n,
    ncol = patient_n,
    dimnames = list(patient_ids, patient_ids)
  )
}

same_counts <- lapply(k_names, function(x) make_zero_mat())
pair_counts <- lapply(k_names, function(x) make_zero_mat())

names(same_counts) <- k_names
names(pair_counts) <- k_names

success_counts <- setNames(integer(length(k_values)), k_names)
failure_counts <- setNames(integer(length(k_values)), k_names)

start_time <- Sys.time()

pb <- txtProgressBar(min = 0, max = reps, style = 3)

for (r in seq_len(reps)) {
  set.seed(12345 + r)

  samp_n <- floor(patient_n * p_item)
  samp_n <- max(samp_n, max_k + 1)
  samp_n <- min(samp_n, patient_n)

  samp_ids <- sample(patient_ids, samp_n, replace = FALSE)
  idx <- match(samp_ids, patient_ids)

  x_sub <- feature_scaled_mat[samp_ids, , drop = FALSE]

  for (kk in k_values) {
    k_name <- paste0("k", kk)

    fit <- tryCatch(
      {
        stats::kmeans(
          x = x_sub,
          centers = kk,
          nstart = kmeans_nstart,
          iter.max = kmeans_iter_max
        )
      },
      error = function(e) NULL
    )

    if (is.null(fit)) {
      failure_counts[k_name] <- failure_counts[k_name] + 1
      next
    }

    clu <- as.character(fit$cluster)
    names(clu) <- samp_ids

    same_sub <- outer(clu, clu, FUN = "==") * 1

    pair_counts[[k_name]][idx, idx] <- pair_counts[[k_name]][idx, idx] + 1
    same_counts[[k_name]][idx, idx] <- same_counts[[k_name]][idx, idx] + same_sub

    success_counts[k_name] <- success_counts[k_name] + 1
  }

  setTxtProgressBar(pb, r)
}

close(pb)

end_time <- Sys.time()

cat("\n共识聚类完成，耗时：",
    round(as.numeric(difftime(end_time, start_time, units = "mins")), 2),
    "分钟\n")

cat("各 k 成功次数：\n")
print(success_counts)

cat("各 k 失败次数：\n")
print(failure_counts)

cat("\n=== 5. 生成共识矩阵 ===\n")

consensus_mats <- lapply(k_names, function(k_name) {
  mat <- same_counts[[k_name]] / pair_counts[[k_name]]

  mat[!is.finite(mat)] <- NA
  mat[is.na(mat)] <- 0

  diag(mat) <- 1
  mat <- (mat + t(mat)) / 2

  mat
})

names(consensus_mats) <- k_names

cat("\n=== 6. 根据共识矩阵生成最终类别 ===\n")

final_clusters <- list()
hclust_list    <- list()

for (kk in k_values) {
  k_name <- paste0("k", kk)

  mat <- consensus_mats[[k_name]]

  dissim <- 1 - mat
  dissim <- (dissim + t(dissim)) / 2
  diag(dissim) <- 0

  hc <- hclust(as.dist(dissim), method = "ward.D2")

  hclust_list[[k_name]] <- hc

  clu <- cutree(hc, k = kk)
  clu <- clu[patient_ids]

  final_clusters[[k_name]] <- clu
}

cluster_assignments <- data.frame(
  patient_id = patient_ids,
  stringsAsFactors = FALSE
)

for (kk in k_values) {
  k_name <- paste0("k", kk)

  cluster_assignments[[paste0("cluster_", k_name)]] <-
    as.integer(final_clusters[[k_name]][patient_ids])
}

readr::write_excel_csv(cluster_assignments, cluster_assignment_csv)

cat("\n=== 7. 共识矩阵诊断 ===\n")

diagnostics_df <- do.call(
  rbind,
  lapply(k_values, function(kk) {
    k_name <- paste0("k", kk)

    mat <- consensus_mats[[k_name]]
    vals <- mat[upper.tri(mat)]
    vals <- vals[is.finite(vals)]

    data.frame(
      k = kk,
      success = success_counts[k_name],
      failure = failure_counts[k_name],
      offdiag_min = min(vals, na.rm = TRUE),
      offdiag_q25 = quantile(vals, 0.25, na.rm = TRUE),
      offdiag_median = median(vals, na.rm = TRUE),
      offdiag_mean = mean(vals, na.rm = TRUE),
      offdiag_q75 = quantile(vals, 0.75, na.rm = TRUE),
      offdiag_max = max(vals, na.rm = TRUE),
      offdiag_nonzero_prop = mean(vals > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

readr::write_excel_csv(diagnostics_df, diagnostics_csv)

print(diagnostics_df)

cat("\n=== 8. 计算 AIC/BIC/ICL/Entropy 模型统计量 ===\n")

calc_soft_membership_from_consensus <- function(consensus_mat, clu) {
  ids <- names(clu)
  cluster_levels <- sort(unique(as.integer(clu)))

  membership <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(cluster_levels),
    dimnames = list(ids, paste0("Class_", cluster_levels))
  )

  for (g in cluster_levels) {
    ids_g <- ids[as.integer(clu) == g]
    col_g <- paste0("Class_", g)

    for (i in ids) {
      use_ids <- ids_g

      if (as.integer(clu[i]) == g && length(ids_g) > 1) {
        use_ids <- setdiff(ids_g, i)
      }

      if (length(use_ids) == 0) {
        membership[i, col_g] <- 1
      } else {
        membership[i, col_g] <- mean(consensus_mat[i, use_ids], na.rm = TRUE)
      }
    }
  }

  row_sums <- rowSums(membership, na.rm = TRUE)

  for (i in seq_len(nrow(membership))) {
    if (!is.finite(row_sums[i]) || row_sums[i] <= 0) {
      membership[i, ] <- 0
      own_col <- paste0("Class_", as.integer(clu[rownames(membership)[i]]))
      membership[i, own_col] <- 1
    } else {
      membership[i, ] <- membership[i, ] / row_sums[i]
    }
  }

  membership
}

calc_kmeans_fit_statistics <- function(kk) {
  k_name <- paste0("k", kk)
  clu <- final_clusters[[k_name]][patient_ids]

  n <- nrow(feature_scaled_mat)
  p <- ncol(feature_scaled_mat)

  cluster_factor <- factor(as.integer(clu), levels = seq_len(kk))
  cluster_sizes <- as.integer(table(cluster_factor))

  centers <- rowsum(feature_scaled_mat[patient_ids, , drop = FALSE], group = cluster_factor)
  centers <- centers / cluster_sizes

  fitted_center <- centers[as.character(cluster_factor), , drop = FALSE]
  wss <- sum((feature_scaled_mat[patient_ids, , drop = FALSE] - fitted_center)^2)

  sigma2 <- wss / (n * p)
  sigma2 <- max(sigma2, .Machine$double.eps)

  loglik <- -0.5 * n * p * (log(2 * pi * sigma2) + 1)

  npar <- kk * p + 1 + (kk - 1)

  AIC <- -2 * loglik + 2 * npar
  BIC <- -2 * loglik + log(n) * npar

  membership <- calc_soft_membership_from_consensus(consensus_mats[[k_name]], clu)
  eps <- 1e-12
  classification_entropy <- -sum(membership * log(pmax(membership, eps)), na.rm = TRUE)
  entropy_std <- 1 - classification_entropy / (n * log(kk))
  entropy_std <- pmin(pmax(entropy_std, 0), 1)

  ICL <- BIC + 2 * classification_entropy

  class_props <- as.numeric(cluster_sizes / n)
  class_cols <- as.list(rep("/", max_k))
  names(class_cols) <- paste0("Class_", seq_len(max_k))

  for (j in seq_len(kk)) {
    class_cols[[paste0("Class_", j)]] <- round(class_props[j], 3)
  }

  out <- data.frame(
    Model = paste0("ConsensusKmeans-", kk),
    Number_of_Class = kk,
    Loglik = round(loglik, 1),
    AIC = round(AIC, 1),
    BIC = round(BIC, 1),
    ICL = round(ICL, 1),
    Entropy = round(entropy_std, 3),
    WSS = round(wss, 3),
    Sigma2 = signif(sigma2, 4),
    Npar = npar,
    stringsAsFactors = FALSE
  )

  cbind(out, as.data.frame(class_cols, stringsAsFactors = FALSE))
}

fit_statistics_df <- do.call(
  rbind,
  lapply(k_values, calc_kmeans_fit_statistics)
)

fit_statistics_df <- fit_statistics_df %>%
  mutate(
    Delta_AIC = AIC - min(AIC, na.rm = TRUE),
    Delta_BIC = BIC - min(BIC, na.rm = TRUE)
  )

cat("\n=== 9. 计算 cluster-consensus 并自动选择最佳 K ===\n")

calc_cluster_consensus <- function(kk) {
  k_name <- paste0("k", kk)

  mat <- consensus_mats[[k_name]]
  clu <- final_clusters[[k_name]]

  out <- lapply(sort(unique(clu)), function(cc) {
    ids <- names(clu)[clu == cc]

    if (length(ids) <= 1) {
      cc_value <- NA_real_
    } else {
      submat <- mat[ids, ids, drop = FALSE]
      cc_value <- mean(submat[upper.tri(submat)], na.rm = TRUE)
    }

    data.frame(
      k = kk,
      cluster = cc,
      n = length(ids),
      cluster_consensus = cc_value,
      singleton = length(ids) <= 1,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

cluster_consensus_df <- do.call(
  rbind,
  lapply(k_values, calc_cluster_consensus)
)

cluster_summary_by_k <- cluster_consensus_df %>%
  group_by(k) %>%
  summarise(
    mean_cluster_consensus = mean(cluster_consensus, na.rm = TRUE),
    min_cluster_consensus  = min(cluster_consensus, na.rm = TRUE),
    min_cluster_n          = min(n, na.rm = TRUE),
    singleton_n            = sum(singleton, na.rm = TRUE),
    .groups = "drop"
  )

min_cluster_n_required <- ifelse(
  patient_n < 40,
  2,
  max(5, ceiling(0.05 * patient_n))
)

selection_df <- fit_statistics_df %>%
  left_join(cluster_summary_by_k, by = c("Number_of_Class" = "k")) %>%
  mutate(
    rank_BIC      = rank(BIC, ties.method = "average"),
    rank_ICL      = rank(ICL, ties.method = "average"),
    rank_Entropy  = rank(-Entropy, ties.method = "average"),
    rank_CC       = rank(-mean_cluster_consensus, ties.method = "average"),
    rank_MinCC    = rank(-min_cluster_consensus, ties.method = "average"),
    small_cluster_penalty = ifelse(min_cluster_n < min_cluster_n_required, 10, 0),
    singleton_penalty     = singleton_n * 10,
    auto_score = rank_BIC + rank_ICL + rank_Entropy + rank_CC + rank_MinCC +
      small_cluster_penalty + singleton_penalty
  ) %>%
  arrange(auto_score, Number_of_Class)

k_sel <- selection_df$Number_of_Class[1]

fit_statistics_df <- fit_statistics_df %>%
  mutate(
    Auto_Selected = ifelse(Number_of_Class == k_sel, "Yes", "No")
  )

selection_df <- selection_df %>%
  mutate(
    Auto_Selected = ifelse(Number_of_Class == k_sel, "Yes", "No")
  )

readr::write_excel_csv(fit_statistics_df, fit_statistics_csv)
readr::write_excel_csv(selection_df, selection_score_csv)

cat("自动选择的最终 K =", k_sel, "\n")
cat("最小类别样本量要求：", min_cluster_n_required, "\n")
cat("自动选择评分表：\n")
print(selection_df)

cat("模型统计量：\n")
print(fit_statistics_df)

selected_dir <- file.path(output_dir, paste0("K", k_sel, "_auto_selected_results"))
dir.create(selected_dir, recursive = TRUE, showWarnings = FALSE)

selected_cluster_assignment_csv        <- file.path(selected_dir, paste0("K", k_sel, "_cluster_assignments.csv"))
selected_patient_membership_csv        <- file.path(selected_dir, paste0("K", k_sel, "_patient_consensus_membership.csv"))
selected_cluster_size_csv              <- file.path(selected_dir, paste0("K", k_sel, "_cluster_size.csv"))
selected_consensus_matrix_csv          <- file.path(selected_dir, paste0("K", k_sel, "_consensus_matrix.csv"))
selected_original_long_csv             <- file.path(selected_dir, paste0("K", k_sel, "_original_data_used_long_with_cluster.csv"))
selected_wide_raw_with_cluster_csv     <- file.path(selected_dir, paste0("K", k_sel, "_patient_wide_features_raw_with_cluster.csv"))
selected_wide_scaled_with_cluster_csv  <- file.path(selected_dir, paste0("K", k_sel, "_patient_wide_features_scaled_with_cluster.csv"))
selected_feature_profile_raw_csv       <- file.path(selected_dir, paste0("K", k_sel, "_feature_profile_raw.csv"))
selected_feature_profile_scaled_csv    <- file.path(selected_dir, paste0("K", k_sel, "_feature_profile_scaled.csv"))
selected_model_fit_statistics_csv      <- file.path(selected_dir, paste0("K", k_sel, "_model_fit_statistics.csv"))
selected_trajectory_summary_csv        <- file.path(selected_dir, paste0("K", k_sel, "_trajectory_summary.csv"))
selected_heatmap_pdf                   <- file.path(selected_dir, paste0("K", k_sel, "_consensus_matrix_heatmap.pdf"))
selected_heatmap_png                   <- file.path(selected_dir, paste0("K", k_sel, "_consensus_matrix_heatmap.png"))
selected_trajectory_pdf                <- file.path(selected_dir, paste0("K", k_sel, "_trajectory_plots.pdf"))
selected_trajectory_png                <- file.path(selected_dir, paste0("K", k_sel, "_trajectory_plots.png"))

cat("\n=== 10. 单独导出自动选择 K 的结果 ===\n")

selected_name <- paste0("k", k_sel)
selected_cluster_labels <- final_clusters[[selected_name]][patient_ids]

selected_patient_clusters <- data.frame(
  patient_id = patient_ids,
  cluster = as.integer(selected_cluster_labels),
  stringsAsFactors = FALSE
)

selected_cluster_size <- selected_patient_clusters %>%
  count(cluster, name = "n") %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(cluster)

selected_membership_mat <- calc_soft_membership_from_consensus(
  consensus_mats[[selected_name]],
  selected_cluster_labels
)

selected_patient_membership <- data.frame(
  patient_id = rownames(selected_membership_mat),
  cluster = as.integer(selected_cluster_labels[rownames(selected_membership_mat)]),
  selected_membership_mat,
  check.names = FALSE
) %>%
  arrange(patient_id)

selected_consensus_df <- data.frame(
  patient_id = rownames(consensus_mats[[selected_name]]),
  consensus_mats[[selected_name]],
  check.names = FALSE
)

selected_wide_raw_with_cluster <- dat_wide %>%
  left_join(selected_patient_clusters, by = "patient_id") %>%
  select(patient_id, cluster, everything())

selected_wide_scaled_with_cluster <- dat_wide_scaled %>%
  left_join(selected_patient_clusters, by = "patient_id") %>%
  select(patient_id, cluster, everything())

selected_original_long_with_cluster <- dat %>%
  left_join(selected_patient_clusters, by = "patient_id") %>%
  arrange(cluster, patient_id, time_num)

make_feature_profile <- function(feature_df, value_vars, cluster_df) {
  feature_df %>%
    left_join(cluster_df, by = "patient_id") %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      across(
        all_of(value_vars),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

selected_feature_profile_raw <- make_feature_profile(
  dat_wide,
  feature_vars,
  selected_patient_clusters
)

selected_feature_profile_scaled <- make_feature_profile(
  dat_wide_scaled,
  keep_feature_vars,
  selected_patient_clusters
)

selected_model_fit_statistics <- fit_statistics_df %>%
  filter(Number_of_Class == k_sel)

readr::write_excel_csv(selected_patient_clusters, selected_cluster_assignment_csv)
readr::write_excel_csv(selected_patient_membership, selected_patient_membership_csv)
readr::write_excel_csv(selected_cluster_size, selected_cluster_size_csv)
readr::write_excel_csv(selected_consensus_df, selected_consensus_matrix_csv)

readr::write_excel_csv(selected_original_long_with_cluster, selected_original_long_csv)
readr::write_excel_csv(selected_wide_raw_with_cluster, selected_wide_raw_with_cluster_csv)
readr::write_excel_csv(selected_wide_scaled_with_cluster, selected_wide_scaled_with_cluster_csv)
readr::write_excel_csv(selected_feature_profile_raw, selected_feature_profile_raw_csv)
readr::write_excel_csv(selected_feature_profile_scaled, selected_feature_profile_scaled_csv)
readr::write_excel_csv(selected_model_fit_statistics, selected_model_fit_statistics_csv)

cat("自动选择 K 的类别分布：\n")
print(selected_cluster_size)

cat("自动选择 K 的模型统计量：\n")
print(selected_model_fit_statistics)

cat("\n=== 11. 绘制各 k 值共识矩阵热图 ===\n")

plot_consensus_heatmap <- function(kk) {
  k_name <- paste0("k", kk)

  mat <- consensus_mats[[k_name]]
  hc  <- hclust_list[[k_name]]

  clu <- final_clusters[[k_name]][rownames(mat)]

  ann <- data.frame(
    cluster = factor(clu)
  )

  rownames(ann) <- rownames(mat)

  heat_cols <- colorRampPalette(
    c("white", "#d7e8f4", "#1f78b4", "#0b3dff")
  )(100)

  ann_colors <- list(
    cluster = make_cluster_palette(levels(ann$cluster))
  )

  ph <- pheatmap::pheatmap(
    mat,
    color             = heat_cols,
    breaks            = seq(0, 1, length.out = 101),
    cluster_rows      = hc,
    cluster_cols      = hc,
    annotation_row    = ann,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    main              = paste0("consensus matrix k=", kk),
    border_color      = NA,
    silent            = TRUE
  )

  ph
}

pdf(combined_heatmap_pdf, width = 10, height = 7)

for (kk in k_values) {
  ph <- plot_consensus_heatmap(kk)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
}

dev.off()

for (kk in k_values) {
  ph <- plot_consensus_heatmap(kk)

  pdf_file <- file.path(heatmap_dir, paste0("consensus_matrix_k", kk, ".pdf"))
  png_file <- file.path(heatmap_dir, paste0("consensus_matrix_k", kk, ".png"))

  pdf(pdf_file, width = 10, height = 7)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  png(png_file, width = 3000, height = 2100, res = 300)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  if (kk == k_sel) {
    pdf(selected_heatmap_pdf, width = 10, height = 7)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()

    png(selected_heatmap_png, width = 3000, height = 2100, res = 300)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
  }
}

cat("\n=== 12. 绘制 AIC/BIC/ICL 图 ===\n")

fit_plot_df <- fit_statistics_df %>%
  select(Number_of_Class, AIC, BIC, ICL) %>%
  pivot_longer(cols = c(AIC, BIC, ICL), names_to = "criterion", values_to = "value")

p_fit_statistics <- ggplot(
  fit_plot_df,
  aes(x = Number_of_Class, y = value, color = criterion, group = criterion)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = k_sel, linetype = "dashed") +
  scale_x_continuous(breaks = k_values) +
  labs(
    title = "Approximate AIC/BIC/ICL for consensus k-means",
    subtitle = paste0("Dashed line indicates auto-selected k = ", k_sel),
    x = "Number of classes, k",
    y = "Criterion value",
    color = "Criterion"
  ) +
  theme_bw_nogrid(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(fit_statistics_pdf, p_fit_statistics, width = 10, height = 6)
ggsave(fit_statistics_png, p_fit_statistics, width = 10, height = 6, dpi = 300)

cat("\n=== 13. 绘制 consensus CDF 与 Delta area ===\n")

get_consensus_values <- function(kk) {
  k_name <- paste0("k", kk)

  mat <- consensus_mats[[k_name]]
  vals <- mat[upper.tri(mat)]
  vals <- vals[is.finite(vals)]
  vals <- pmin(pmax(vals, 0), 1)

  vals
}

cdf_grid <- seq(0, 1, length.out = 201)

cdf_df <- do.call(
  rbind,
  lapply(k_values, function(kk) {
    vals <- get_consensus_values(kk)
    ec <- stats::ecdf(vals)

    data.frame(
      k = factor(kk, levels = k_values),
      consensus_index = cdf_grid,
      CDF = ec(cdf_grid)
    )
  })
)

p_cdf <- ggplot(cdf_df, aes(x = consensus_index, y = CDF, color = k)) +
  geom_line(linewidth = 1) +
  labs(
    title = "consensus CDF",
    x = "consensus index",
    y = "CDF",
    color = "k"
  ) +
  theme_bw_nogrid(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(cdf_pdf, p_cdf, width = 10, height = 6)
ggsave(cdf_png, p_cdf, width = 10, height = 6, dpi = 300)

cdf_area <- sapply(k_values, function(kk) {
  vals <- get_consensus_values(kk)
  mean(1 - vals, na.rm = TRUE)
})

delta_area <- rep(NA_real_, length(k_values))
delta_area[1] <- cdf_area[1]

for (i in 2:length(k_values)) {
  if (cdf_area[i - 1] == 0) {
    delta_area[i] <- NA_real_
  } else {
    delta_area[i] <- (cdf_area[i] - cdf_area[i - 1]) / cdf_area[i - 1]
  }
}

df_delta <- data.frame(
  k = k_values,
  cdf_area = as.numeric(cdf_area),
  delta_area = as.numeric(delta_area)
)

p_delta <- ggplot(df_delta, aes(x = k, y = delta_area)) +
  geom_line() +
  geom_point(size = 3, shape = 21, fill = "white") +
  geom_vline(xintercept = k_sel, linetype = "dashed") +
  scale_x_continuous(breaks = k_values) +
  labs(
    title = "Delta area",
    subtitle = paste0("Dashed line indicates auto-selected k = ", k_sel),
    x = "k",
    y = "relative change in area under CDF curve"
  ) +
  theme_bw_nogrid(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(delta_area_pdf, p_delta, width = 10, height = 6)
ggsave(delta_area_png, p_delta, width = 10, height = 6, dpi = 300)

cat("\n=== 14. 绘制 cluster-consensus 柱状图 ===\n")

cluster_consensus_plot_df <- cluster_consensus_df %>%
  mutate(
    k = factor(k, levels = k_values),
    cluster = factor(cluster),
    plot_value = ifelse(is.na(cluster_consensus), 0, cluster_consensus),
    label_text = ifelse(
      singleton,
      paste0("n=", n, "\nsingle"),
      paste0("n=", n)
    )
  )

cluster_consensus_palette <- make_cluster_palette(levels(cluster_consensus_plot_df$cluster))

p_cluster_consensus <- ggplot(
  cluster_consensus_plot_df,
  aes(x = cluster, y = plot_value, fill = cluster)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = label_text),
    vjust = -0.2,
    size = 3
  ) +
  facet_grid(. ~ k, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cluster_consensus_palette) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(
    title = "cluster-consensus",
    x = "cluster",
    y = "cluster consensus"
  ) +
  theme_bw_nogrid(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(cluster_consensus_pdf, p_cluster_consensus, width = 12, height = 5)
ggsave(cluster_consensus_png, p_cluster_consensus, width = 12, height = 5, dpi = 300)

cat("\n=== 15. 输出合并模型选择图 ===\n")

p_model_selection <- p_fit_statistics / p_delta / p_cdf / p_cluster_consensus +
  patchwork::plot_layout(heights = c(1, 1, 1, 0.8))

ggsave(model_selection_pdf, p_model_selection, width = 12, height = 18)
ggsave(model_selection_png, p_model_selection, width = 12, height = 18, dpi = 300)

cat("\n=== 16. 绘制自动选择 K 的三指标轨迹图 ===\n")

cluster_labels <- final_clusters[[paste0("k", k_sel)]]

patient_clusters <- data.frame(
  patient_id = names(cluster_labels),
  cluster = factor(cluster_labels),
  stringsAsFactors = FALSE
)

traj_data <- dat %>%
  left_join(patient_clusters, by = "patient_id")

traj_long <- traj_data %>%
  pivot_longer(
    cols = all_of(safe_vars),
    names_to = "indicator_var",
    values_to = "value"
  ) %>%
  left_join(indicator_map, by = "indicator_var")

traj_summary <- traj_long %>%
  group_by(indicator_label, time_num, cluster) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sd = ifelse(is.na(sd), 0, sd),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

readr::write_excel_csv(traj_summary, trajectory_summary_csv)
readr::write_excel_csv(traj_summary, selected_trajectory_summary_csv)

selected_cluster_palette <- make_cluster_palette(levels(traj_summary$cluster))

p_traj <- ggplot(
  traj_summary,
  aes(x = time_num, y = mean, color = cluster, group = cluster)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = cluster),
    alpha = 0.2,
    linetype = 0
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~ indicator_label,
    scales = "free_y",
    ncol = 2,
    strip.position = "left"
  ) +
  scale_x_continuous(breaks = time_points) +
  scale_color_manual(values = selected_cluster_palette) +
  scale_fill_manual(values = selected_cluster_palette) +
  labs(
    title = paste0("Trajectories by Cluster, auto-selected k = ", k_sel),
    x = "Day",
    y = NULL,
    color = "Cluster",
    fill = "Cluster"
  ) +
  theme_bw_nogrid(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold")
  )

ggsave(trajectory_pdf, p_traj, width = 11, height = 8)
ggsave(trajectory_png, p_traj, width = 11, height = 8, dpi = 300)
ggsave(selected_trajectory_pdf, p_traj, width = 11, height = 8)
ggsave(selected_trajectory_png, p_traj, width = 11, height = 8, dpi = 300)

p_traj_noCI <- ggplot(
  traj_summary,
  aes(x = time_num, y = mean, color = cluster, group = cluster)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~ indicator_label,
    scales = "free_y",
    ncol = 2,
    strip.position = "left"
  ) +
  scale_x_continuous(breaks = time_points) +
  scale_color_manual(values = selected_cluster_palette) +
  labs(
    title = paste0("Trajectories by Cluster, auto-selected k = ", k_sel, " (no CI)"),
    x = "Day",
    y = "Mean",
    color = "Cluster"
  ) +
  theme_bw_nogrid(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold")
  )

trajectory_noCI_pdf <- file.path(selected_dir, paste0("K", k_sel, "_trajectory_plots_noCI.pdf"))
trajectory_noCI_png <- file.path(selected_dir, paste0("K", k_sel, "_trajectory_plots_noCI.png"))

ggsave(trajectory_noCI_pdf, p_traj_noCI, width = 11, height = 8)
ggsave(trajectory_noCI_png, p_traj_noCI, width = 11, height = 8, dpi = 300)

traj_long_std <- traj_long %>%
  group_by(indicator_label) %>%
  mutate(z_value = as.numeric(scale(value))) %>%
  ungroup()

traj_std_summary <- traj_long_std %>%
  group_by(cluster, indicator_label, time_num) %>%
  summarise(
    mean = mean(z_value, na.rm = TRUE),
    sd = sd(z_value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sd = ifelse(is.na(sd), 0, sd),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

for (clu_id in unique(traj_std_summary$cluster)) {
  clu_label <- as.character(clu_id)

  df_raw <- traj_summary %>% filter(cluster == clu_id)

  p_raw_CI <- ggplot(
    df_raw,
    aes(x = time_num, y = mean, color = indicator_label, group = indicator_label)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper, fill = indicator_label),
      alpha = 0.2,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = time_points) +
    labs(
      title = paste0("Cluster ", clu_label, " trajectories (raw, CI)"),
      x = "Day",
      y = "Value",
      color = "Indicator",
      fill = "Indicator"
    ) +
    theme_bw_nogrid(base_size = 12)

  ggsave(
    file.path(selected_dir, paste0("Cluster", clu_label, "_raw_CI.pdf")),
    p_raw_CI,
    width = 10,
    height = 6
  )

  p_raw_noCI <- ggplot(
    df_raw,
    aes(x = time_num, y = mean, color = indicator_label, group = indicator_label)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = time_points) +
    labs(
      title = paste0("Cluster ", clu_label, " trajectories (raw, no CI)"),
      x = "Day",
      y = "Value",
      color = "Indicator"
    ) +
    theme_bw_nogrid(base_size = 12)

  ggsave(
    file.path(selected_dir, paste0("Cluster", clu_label, "_raw_noCI.pdf")),
    p_raw_noCI,
    width = 10,
    height = 6
  )

  df_std <- traj_std_summary %>% filter(cluster == clu_id)

  p_std_CI <- ggplot(
    df_std,
    aes(x = time_num, y = mean, color = indicator_label, group = indicator_label)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper, fill = indicator_label),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = time_points) +
    labs(
      title = paste0("Cluster ", clu_label, " trajectories (standardized, CI)"),
      x = "Day",
      y = "Standardized value (Z)",
      color = "Indicator",
      fill = "Indicator"
    ) +
    theme_bw_nogrid(base_size = 12)

  ggsave(
    file.path(selected_dir, paste0("Cluster", clu_label, "_standardized_CI.pdf")),
    p_std_CI,
    width = 10,
    height = 6
  )

  p_std_noCI <- ggplot(
    df_std,
    aes(x = time_num, y = mean, color = indicator_label, group = indicator_label)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = time_points) +
    labs(
      title = paste0("Cluster ", clu_label, " trajectories (standardized, no CI)"),
      x = "Day",
      y = "Standardized value (Z)",
      color = "Indicator"
    ) +
    theme_bw_nogrid(base_size = 12)

  ggsave(
    file.path(selected_dir, paste0("Cluster", clu_label, "_standardized_noCI.pdf")),
    p_std_noCI,
    width = 10,
    height = 6
  )
}

p_std_all_CI <- ggplot(
  traj_std_summary,
  aes(x = time_num, y = mean, color = indicator_label, group = indicator_label)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = indicator_label),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ cluster) +
  scale_x_continuous(breaks = time_points) +
  labs(
    title = paste0("All Clusters Standardized Trajectories, auto-selected k = ", k_sel, " (CI)"),
    x = "Day",
    y = "Z-score",
    color = "Indicator",
    fill = "Indicator"
  ) +
  theme_bw_nogrid(base_size = 12)

p_std_all_noCI <- ggplot(
  traj_std_summary,
  aes(x = time_num, y = mean, color = indicator_label, group = indicator_label)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ cluster) +
  scale_x_continuous(breaks = time_points) +
  labs(
    title = paste0("All Clusters Standardized Trajectories, auto-selected k = ", k_sel, " (no CI)"),
    x = "Day",
    y = "Z-score",
    color = "Indicator"
  ) +
  theme_bw_nogrid(base_size = 12)

ggsave(
  file.path(selected_dir, paste0("K", k_sel, "_AllClusters_standardized_CI.pdf")),
  p_std_all_CI,
  width = 13,
  height = 8
)

ggsave(
  file.path(selected_dir, paste0("K", k_sel, "_AllClusters_standardized_noCI.pdf")),
  p_std_all_noCI,
  width = 13,
  height = 8
)

cat("\n=== 17. 保存 RData 与 summary ===\n")

save(
  raw_data,
  dat,
  dat_wide,
  dat_wide_scaled,
  feature_scaled_mat,
  target_vars,
  source_vars,
  indicator_map,
  consensus_mats,
  final_clusters,
  hclust_list,
  diagnostics_df,
  fit_statistics_df,
  cluster_consensus_df,
  selection_df,
  k_sel,
  selected_patient_clusters,
  selected_cluster_size,
  traj_summary,
  traj_std_summary,
  file = rdata_file
)

summary_lines <- c(
  "三指标多变量轨迹 consensus k-means 聚类分析完成",
  paste0("输入文件：", input_file),
  paste0("输出目录：", output_dir),
  paste0("聚类指标：", paste(target_vars, collapse = ", ")),
  paste0("时间点：", paste(time_points, collapse = ", ")),
  paste0("最终纳入患者数：", patient_n),
  paste0("原始宽表特征数：", length(feature_vars)),
  paste0("用于聚类的非零方差特征数：", length(keep_feature_vars)),
  paste0("k 搜索范围：", paste(k_values, collapse = ", ")),
  paste0("重抽样次数：", reps),
  paste0("每次抽样比例：", p_item),
  paste0("k-means nstart：", kmeans_nstart),
  paste0("k-means iter.max：", kmeans_iter_max),
  paste0("自动选择的最终 K：", k_sel),
  paste0("最小类别样本量要求：", min_cluster_n_required),
  "",
  "固定 Cluster 颜色：",
  "Cluster 1：#c30078",
  "Cluster 2：#4370b4",
  "Cluster 3：#549f9a",
  "",
  "图形更新：",
  "ggplot 图已去掉网格线",
  "三指标轨迹图的指标分面标签已放在图左侧",
  "原始单位带 95% CI 轨迹图已去掉 Mean (95% CI) 字样",
  "",
  "自动选择 K 的类别分布：",
  paste(capture.output(print(selected_cluster_size)), collapse = "\n"),
  "",
  "自动选择 K 的模型统计量：",
  paste(capture.output(print(selected_model_fit_statistics)), collapse = "\n")
)

writeLines(summary_lines, summary_txt)

cat("\n===== 全部完成！ =====\n")
cat("输出目录：", output_dir, "\n")
cat("自动选择的最终 K：", k_sel, "\n")
cat("自动选择 K 的结果目录：", selected_dir, "\n")
cat("RData 已保存：", rdata_file, "\n")
cat("Summary 已保存：", summary_txt, "\n")
