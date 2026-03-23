# -*- coding: UTF-8 -*-
# ============================================================
# 单细胞 RNA-seq 完整 QC Pipeline
# 包含：空细胞检测 -> 环境RNA去除 -> QC指标(mt+Hb) ->
#       保存原始 -> QC过滤 -> 保存过滤后 -> DoubletFinder -> 合并
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(clustree)
library(DoubletFinder)
library(igraph)
library(cluster)
library(bluster)
library(SeuratWrappers)
library(DropletUtils)
library(SoupX)
library(celda)
library(scater)
library(SingleCellExperiment)

# ============================================================
# 全局参数设置
# ============================================================
set.seed(42)

projectPath <- "/mnt2/wanggd_group/zjj/ZengMin20260310"
setwd(projectPath)

# 确保输出目录存在
dir.create(file.path(projectPath, "Output"),               showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "EmptyDrops"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "AmbientRNA"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "DoubletFinder"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "CellCycle"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Data"),                  showWarnings = FALSE, recursive = TRUE)

# ---- 动态日期标签（格式：YYMMDD） ----
date_tag <- format(Sys.Date(), "%y%m%d")
message(sprintf("Date tag: %s", date_tag))

# ---- QC 阈值（集中管理） ----
MIN_FEATURES   <- 200
MAX_FEATURES   <- 5000
MAX_MT_PERCENT <- 20
MAX_HB_PERCENT <- 5
MIN_COUNTS     <- 500
PC_DIMS        <- 1:20
MIN_NOVELTY    <- 0.8  # 作为示例阈值

# ---- 物种相关 Pattern ----
# 小鼠: mt -> "^mt-",  Hb -> "^Hb[ab]-"
MT_PATTERN <- "^MT-"     # 如犬/其他需自行调整
HB_PATTERN <- "^HBA|^HBB"

# ---- EmptyDrops 参数 ----
EMPTY_DROPS_FDR   <- 0.01
EMPTY_DROPS_LOWER <- 100

# ---- DecontX 参数 ----
DECONTX_CONTAMINATION_THRESHOLD <- 0.5

# ---- 动态文件名生成函数 ----
make_filename <- function(n_samples, suffix, prefix = "", ext = "Rdata") {
  fname <- paste0(prefix, date_tag, "_", n_samples, "Samples_", suffix, ".", ext)
  return(fname)
}

# ============================================================
# 读取 Cell Ranger 的初始细胞数（若有）
# 这里实现一个健壮函数，尽量从 outs/stats.txt 或 metrics.json 获取
# 返回整数或 NA
get_cellranger_cells <- function(sample_dir) {
  # 假设 sample_dir 指向 Data/QuantifyRawdata/样本名
  # 常见的 stats 文件路径：
  # outs/stats.txt 或 outs/metrics.json
  stats_path1 <- file.path(sample_dir, "outs", "stats.txt")
  metrics_path <- file.path(sample_dir, "outs", "metrics.json")
  n_cells <- NA

  if (file.exists(stats_path1)) {
    # 尝试读取 stats.txt, 找到包含 "number of reads" 或 "reads" 的行，需要看实际文件
    stat_lines <- readLines(stats_path1, warn = FALSE)
    # 常见字段：" Peanut": 不同版本可能不同，请尝试匹配 "nCells"、"cell barcodes"、"nGenes"
    # 尝试解析常见字段
    # 直接搜索包含 "Total Cells" 或 "Cells"
    tc_line <- grep("Total Cells|Total\\s+Cells|Cells", stat_lines, value = TRUE)
    if (length(tc_line) > 0) {
      # 从行中提取数字
      nums <- as.numeric(gsub("[^0-9]", "", tc_line[1]))
      if (!is.na(nums) && length(nums) > 0) n_cells <- nums
    }
  }

  if (is.na(n_cells) && file.exists(metrics_path)) {
    # 读取 json
    json <- tryCatch(jsonlite::read_json(metrics_path), error = function(e) NULL)
    if (!is.null(json)) {
      if (!is.null(json$summary) && !is.null(json$summary$nCells)) {
        n_cells <- as.numeric(json$summary$nCells)
      } else if (!is.null(json$nCells)) {
        n_cells <- as.numeric(json$nCells)
      }
    }
  }

  # 回退：尝试在 outs 目录下的其他可能文件
  if (is.na(n_cells)) {
    alt <- file.path(sample_dir, "outs", "cell_barcodes.tsv")
    if (file.exists(alt)) {
      n_cells <- length(readLines(alt))
    }
  }

  return(as.integer(n_cells))
}

# ============================================================
# 1. 数据导入 + 空细胞检测 + 环境 RNA 去除 + QC 指标计算
# ============================================================
samples <- list.files(file.path(projectPath, "Data", "QuantifyRawdata"))

seurat_list_raw      <- list()   # 过滤前
seurat_list_filtered <- list()   # 过滤后
sample_ids           <- c()
qc_summary           <- list()

for (sample in samples) {
  tryCatch({
    message(paste("\n========================================"))
    message(paste("Processing sample:", sample))
    message(paste("========================================"))

    # ---------- 路径设置 ----------
    base_path     <- file.path(projectPath, "Data", "QuantifyRawdata", sample)
    filtered_path <- file.path(base_path, "filtered_cell_gene_matrix")
    raw_path      <- file.path(base_path, "raw_cell_gene_matrix")

    if (!dir.exists(filtered_path)) {
      warning(paste("Filtered matrix not found:", filtered_path))
      next
    }

    has_raw <- dir.exists(raw_path)
    message(sprintf("  Raw matrix available: %s", has_raw))

    # ---------- 读取数据 ----------
    filtered_data <- Read10X(data.dir = filtered_path)
    if (has_raw) {
      raw_data <- Read10X(data.dir = raw_path)
    }

    # ==========================================================
    # 步骤 A：空细胞率检测（EmptyDrops）
    # ==========================================================
    if (has_raw) {
      message("  Running emptyDrops on raw matrix...")

      set.seed(42)
      e.out <- emptyDrops(raw_data, lower = EMPTY_DROPS_LOWER)

      is_cell          <- e.out$FDR <= EMPTY_DROPS_FDR & !is.na(e.out$FDR)
      n_empty          <- sum(!is_cell, na.rm = TRUE)
      n_real_cells     <- sum(is_cell, na.rm = TRUE)
      n_total_barcodes <- nrow(e.out)
      empty_rate       <- round((1 - n_real_cells / n_total_barcodes) * 100, 2)

      message(sprintf("  EmptyDrops: %d total barcodes, %d real cells, %.2f%% empty",
                      n_total_barcodes, n_real_cells, empty_rate))

      # 诊断图
      ed_df <- data.frame(
        total_UMI     = e.out$Total,
        neg_log10_FDR = -log10(e.out$FDR),
        is_cell       = is_cell
      )
      ed_df <- ed_df[!is.na(ed_df$neg_log10_FDR), ]

      p_empty <- ggplot(ed_df, aes(x = total_UMI, y = neg_log10_FDR, color = is_cell)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_x_log10() +
        scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray70"),
                           labels = c("TRUE" = "Cell", "FALSE" = "Empty"),
                           name   = "Classification") +
        geom_hline(yintercept = -log10(EMPTY_DROPS_FDR), linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title    = paste0(sample, ": EmptyDrops Detection"),
             subtitle = paste0("Real cells: ", n_real_cells, "  |  Empty rate: ", empty_rate, "%"),
             x = "Total UMI (log10)", y = "-log10(FDR)")
      ggsave(file.path(projectPath, "Output", "EmptyDrops",
                       paste0(sample, "_emptyDrops.png")),
             plot = p_empty, width = 8, height = 6, dpi = 300)

      real_cell_barcodes <- rownames(e.out)[is_cell]
      common_barcodes    <- intersect(colnames(filtered_data), real_cell_barcodes)
      filtered_data      <- filtered_data[, common_barcodes]

      message(sprintf("  After emptyDrops: %d cells retained", ncol(filtered_data)))

    } else {
      message("  No raw matrix. Using UMI threshold for empty cell filtering.")

      col_sums     <- Matrix::colSums(filtered_data)
      n_before     <- ncol(filtered_data)
      keep_cells   <- col_sums >= EMPTY_DROPS_LOWER
      n_empty_like <- sum(!keep_cells)
      empty_rate   <- round(n_empty_like / n_before * 100, 2)

      sorted_umi <- sort(col_sums, decreasing = TRUE)
      knee_df    <- data.frame(rank = seq_along(sorted_umi), UMI = sorted_umi)

      p_knee <- ggplot(knee_df, aes(x = rank, y = UMI)) +
        geom_line(color = "steelblue") +
        scale_x_log10() + scale_y_log10() +
        geom_hline(yintercept = EMPTY_DROPS_LOWER, linetype = "dashed", color = "red") +
        annotate("text", x = max(knee_df$rank) * 0.5, y = EMPTY_DROPS_LOWER * 2,
                 label = paste0("Threshold: ", EMPTY_DROPS_LOWER, " UMI"),
                 color = "red", size = 4) +
        theme_classic() +
        labs(title    = paste0(sample, ": Barcode Rank Plot"),
             subtitle = paste0("Removed: ", n_empty_like, "  |  Empty-like rate: ", empty_rate, "%"),
             x = "Barcode Rank (log10)", y = "Total UMI (log10)")
      ggsave(file.path(projectPath, "Output", "EmptyDrops",
                       paste0(sample, "_kneePlot.png")),
             plot = p_knee, width = 8, height = 6, dpi = 300)

      filtered_data    <- filtered_data[, keep_cells]
      n_real_cells     <- ncol(filtered_data)
      n_total_barcodes <- n_before
      message(sprintf("  After UMI threshold: %d -> %d cells", n_before, n_real_cells))
    }

    # ==========================================================
    # 步骤 B：环境 RNA 污染估算与去除
    # ==========================================================
    if (has_raw) {
      # ---- 方案 A：SoupX ----
      message("  Running SoupX...")

      temp_obj <- CreateSeuratObject(counts = filtered_data, project = sample) %>%
        NormalizeData(verbose = FALSE) %>%
        FindVariableFeatures(verbose = FALSE) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(verbose = FALSE, seed.use = 42) %>%
        FindNeighbors(dims = 1:20, verbose = FALSE) %>%
        FindClusters(resolution = 0.8, verbose = FALSE) %>%
        RunUMAP(dims = 1:20, verbose = FALSE, seed.use = 42)

      sc       <- SoupChannel(tod = raw_data, toc = filtered_data)
      clusters <- setNames(as.character(Idents(temp_obj)), colnames(temp_obj))
      sc       <- setClusters(sc, clusters)
      sc       <- setDR(sc, Embeddings(temp_obj, "umap"), reductName = "UMAP")
      sc       <- autoEstCont(sc, verbose = FALSE)

      contamination_fraction <- sc$fit$rhoEst
      mean_contamination     <- contamination_fraction
      message(sprintf("  SoupX contamination: %.4f (%.2f%%)",
                      contamination_fraction, contamination_fraction * 100))

      corrected_counts <- adjustCounts(sc)

      # SoupX 诊断图
      png(file.path(projectPath, "Output", "AmbientRNA",
                    paste0(sample, "_SoupX_diagnostics.png")),
          width = 1200, height = 800, res = 150)
      top_genes     <- head(sort(Matrix::rowSums(filtered_data), decreasing = TRUE), 6)
      top_gene_name <- names(top_genes)[1]
      if (!is.null(top_gene_name)) {
        tryCatch(plotChangeMap(sc, geneSet = top_gene_name, DR = "UMAP"),
                 error = function(e) message("  SoupX plotChangeMap failed: ", conditionMessage(e)))
      }
      dev.off()

      seurat_obj <- CreateSeuratObject(counts = corrected_counts, project = sample,
                                       min.features = 0, min.cells = 0)
      seurat_obj$ambient_contamination <- contamination_fraction
      seurat_obj$ambient_method        <- "SoupX"

      rm(temp_obj, sc); gc()

    } else {
      # ---- 方案 B：DecontX ----
      message("  Running DecontX...")

      sce <- SingleCellExperiment(assays = list(counts = filtered_data))
      set.seed(42)
      sce <- decontX(sce)

      contamination_scores <- colData(sce)$decontX_contamination
      mean_contamination   <- mean(contamination_scores, na.rm = TRUE)
      median_contamination <- median(contamination_scores, na.rm = TRUE)

      message(sprintf("  DecontX - Mean: %.2f%%, Median: %.2f%%",
                      mean_contamination * 100, median_contamination * 100))

      # 污染分布直方图
      contam_df <- data.frame(contamination = contamination_scores)
      p_contam <- ggplot(contam_df, aes(x = contamination)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
        geom_vline(xintercept = mean_contamination, linetype = "dashed", color = "red", linewidth = 1) +
        geom_vline(xintercept = DECONTX_CONTAMINATION_THRESHOLD,
                   linetype = "dotted", color = "darkred", linewidth = 1) +
        annotate("text", x = mean_contamination + 0.05, y = Inf,
                 label = paste0("Mean: ", round(mean_contamination * 100, 1), "%"),
                 vjust = 2, color = "red", size = 4) +
        annotate("text", x = DECONTX_CONTAMINATION_THRESHOLD + 0.05, y = Inf,
                 label = paste0("Threshold: ", DECONTX_CONTAMINATION_THRESHOLD * 100, "%"),
                 vjust = 4, color = "darkred", size = 3.5) +
        theme_classic() +
        labs(title = paste0(sample, ": DecontX Contamination Distribution"),
             x = "Contamination Fraction", y = "Number of Cells")
      ggsave(file.path(projectPath, "Output", "AmbientRNA",
                       paste0(sample, "_DecontX_contamination_hist.png")),
             plot = p_contam, width = 8, height = 6, dpi = 300)

      # UMAP 上展示污染
      sce <- scater::runPCA(sce, exprs_values = "decontXcounts", seed = 42)
      sce <- scater::runUMAP(sce, dimred = "PCA")

      umap_df <- data.frame(
        UMAP1         = reducedDim(sce, "UMAP")[, 1],
        UMAP2         = reducedDim(sce, "UMAP")[, 2],
        contamination = contamination_scores
      )
      p_umap_contam <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = contamination)) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_viridis_c(option = "inferno", direction = -1, name = "Contamination") +
        theme_classic() +
        labs(title = paste0(sample, ": DecontX Contamination on UMAP"))
      ggsave(file.path(projectPath, "Output", "AmbientRNA",
                       paste0(sample, "_DecontX_UMAP.png")),
             plot = p_umap_contam, width = 9, height = 7, dpi = 300)

      corrected_counts <- decontXcounts(sce)

      seurat_obj <- CreateSeuratObject(counts = corrected_counts, project = sample,
                                       min.features = 0, min.cells = 0)
      seurat_obj$ambient_contamination <- contamination_scores
      seurat_obj$ambient_method        <- "DecontX"
      seurat_obj$high_contamination    <- contamination_scores > DECONTX_CONTAMINATION_THRESHOLD

      n_high_contam <- sum(seurat_obj$high_contamination)
      message(sprintf("  High contamination cells (>%.0f%%): %d (%.2f%%)",
                      DECONTX_CONTAMINATION_THRESHOLD * 100,
                      n_high_contam, n_high_contam / ncol(seurat_obj) * 100))

      rm(sce); gc()
    }

    # ==========================================================
    # 步骤 C：计算 QC 指标（mt + Hb）
    # ==========================================================
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = MT_PATTERN)
    seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = HB_PATTERN)

    # Hb 基因是否匹配
    hb_genes_found <- grep(HB_PATTERN, rownames(seurat_obj), value = TRUE)
    mt_genes_found <- grep(MT_PATTERN, rownames(seurat_obj), value = TRUE)
    message(sprintf("  MT genes matched: %d (%s)", length(mt_genes_found),
                    paste(head(mt_genes_found, 5), collapse = ", ")))
    message(sprintf("  Hb genes matched: %d (%s)", length(hb_genes_found),
                    paste(head(hb_genes_found, 5), collapse = ", ")))
    if (length(hb_genes_found) == 0) {
      warning("  No Hb genes found! Check HB_PATTERN or species. percent.hb will be all 0.")
    }

    # 额外的 ribosome 和 novelty 指标（可选）
    RIBO_PATTERN <- "^RPS|^RPL"
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = RIBO_PATTERN)
    ribo_genes_found <- grep(RIBO_PATTERN, rownames(seurat_obj), value = TRUE)
    message(sprintf("  Ribo genes matched: %d (%s)",
                    length(ribo_genes_found),
                    paste(head(ribo_genes_found, 5), collapse = ", ")))
    seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

    # ---------- 质控前小提琴图 ----------
    vln_nCount   <- VlnPlot(seurat_obj, features = "nCount_RNA",   pt.size = 0) + ggtitle("nCount_RNA")
    vln_nFeature <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + ggtitle("nFeature_RNA")
    vln_mt       <- VlnPlot(seurat_obj, features = "percent.mt",   pt.size = 0) + ggtitle("percent.mt")
    vln_hb       <- VlnPlot(seurat_obj, features = "percent.hb",   pt.size = 0) + ggtitle("percent.hb")
    vln_ribo     <- VlnPlot(seurat_obj, features = "percent.ribo", pt.size = 0) + ggtitle("percent.ribo")
    vln_novelty  <- VlnPlot(seurat_obj, features = "log10GenesPerUMI", pt.size = 0) + ggtitle("Novelty")

    combined_vln_full <- (vln_nCount | vln_nFeature | vln_ribo) /
                         (vln_mt | vln_hb | vln_novelty)
    ggsave(file.path(projectPath, "Output",
                     paste0(sample, ".QC_6Metrics_BeforeFilter.png")),
           plot = combined_vln_full, width = 15, height = 10, dpi = 300)

    # ==========================================================  
    # ★ 先保存原始（未过滤）对象到 seurat_list_raw  
    # ==========================================================  
    seurat_list_raw[[sample]] <- seurat_obj 

    # ==========================================================  
    # 步骤 D：QC 过滤（subset）  
    # ==========================================================  
    cells_before <- ncol(seurat_obj)

    seurat_obj_filtered <- subset(seurat_obj,
                                  subset = nFeature_RNA > MIN_FEATURES &
                                           nFeature_RNA < MAX_FEATURES &
                                           nCount_RNA   > MIN_COUNTS &
                                           percent.mt   < MAX_MT_PERCENT &
                                           percent.hb   < MAX_HB_PERCENT &
                                           log10GenesPerUMI > MIN_NOVELTY)

    cells_after <- ncol(seurat_obj_filtered)

    n_low_feature  <- sum(seurat_obj$nFeature_RNA <= MIN_FEATURES)
    n_high_feature <- sum(seurat_obj$nFeature_RNA >= MAX_FEATURES)
    n_low_count    <- sum(seurat_obj$nCount_RNA   <= MIN_COUNTS)
    n_high_mt      <- sum(seurat_obj$percent.mt   >= MAX_MT_PERCENT)
    n_high_hb      <- sum(seurat_obj$percent.hb   >= MAX_HB_PERCENT)
    n_low_novelty  <- sum(seurat_obj$log10GenesPerUMI <= MIN_NOVELTY, na.rm = TRUE)

    message(sprintf("  [%s] QC filtering: %d -> %d cells (removed %d total)",  
                    sample, cells_before, cells_after, cells_before - cells_after))  
    message(sprintf("    - nFeature < %d:  %d cells", MIN_FEATURES, n_low_feature))  
    message(sprintf("    - nFeature > %d: %d cells", MAX_FEATURES, n_high_feature))  
    message(sprintf("    - nCount   < %d:  %d cells", MIN_COUNTS, n_low_count))  
    message(sprintf("    - %%mt     >= %d%%:  %d cells", MAX_MT_PERCENT, n_high_mt))  
    message(sprintf("    - %%hb     >= %d%%:   %d cells", MAX_HB_PERCENT, n_high_hb))  
    message(sprintf("    - novelty  <= %.1f:     %d cells", MIN_NOVELTY, n_low_novelty))

    # 过滤后统计
    vln_nCount_f   <- VlnPlot(seurat_obj_filtered, features = "nCount_RNA",   pt.size = 0) + ggtitle("nCount_RNA")
    vln_nFeature_f <- VlnPlot(seurat_obj_filtered, features = "nFeature_RNA", pt.size = 0) + ggtitle("nFeature_RNA")
    vln_mt_f       <- VlnPlot(seurat_obj_filtered, features = "percent.mt",   pt.size = 0) + ggtitle("percent.mt")
    vln_hb_f       <- VlnPlot(seurat_obj_filtered, features = "percent.hb",   pt.size = 0) + ggtitle("percent.hb")
    vln_ribo_f     <- VlnPlot(seurat_obj_filtered, features = "percent.ribo", pt.size = 0) + ggtitle("percent.ribo")
    vln_novelty_f  <- VlnPlot(seurat_obj_filtered, features = "log10GenesPerUMI", pt.size = 0) + ggtitle("Novelty")

    combined_vln_f <- (vln_nCount_f | vln_nFeature_f | vln_mt_f) / (vln_hb_f | vln_ribo_f | vln_novelty_f)
    ggsave(file.path(projectPath, "Output",
                     paste0(sample, ".QC_ViolinPlots_AfterFilter.png")),
           plot = combined_vln_f, width = 15, height = 10, dpi = 300)

    # 存入过滤后列表
    seurat_list_filtered[[sample]] <- seurat_obj_filtered
    sample_ids <- c(sample_ids, sample)

    # 记录 QC 摘要
    qc_summary[[sample]] <- data.frame(
      Sample              = sample,
      Raw_Barcodes        = if (!is.na(get_cellranger_cells(base_path)) && get_cellranger_cells(base_path) > 0) get_cellranger_cells(base_path) else NA,
      After_EmptyDrops    = n_real_cells,
      Empty_Rate_Pct      = empty_rate,
      Ambient_Method      = ifelse(has_raw, "SoupX", "DecontX"),
      Ambient_Rate_Pct    = round(mean_contamination * 100, 2),
      Cells_Before_QC     = cells_before,
      Removed_LowFeature  = n_low_feature,
      Removed_HighFeature = n_high_feature,
      Removed_LowCount    = n_low_count,
      Removed_HighMT      = n_high_mt,
      Removed_HighHb      = n_high_hb,
      Cells_After_QC      = cells_after,
      stringsAsFactors    = FALSE
    )

  }, error = function(e) {
    message(paste("Error processing sample:", sample))
    message(conditionMessage(e))
  })
}

# ============================================================
# 输出 QC 汇总表
# ============================================================
# 将 qc_summary 列表转为 DataFrame
qc_summary_df <- do.call(rbind, qc_summary)

# 为了确保 Raw Cells 列存在且与样本对齐，按样本顺序填充 NA
if (!"Raw_Barcodes" %in% names(qc_summary_df)) {
  qc_summary_df$Raw_Barcodes <- NA
}
write.csv(qc_summary_df,
          file.path(projectPath, "Output",
                    paste0(date_tag, "_QC_Summary_Full.csv")),
          row.names = FALSE)
message("\n========== QC Summary ==========")
print(qc_summary_df)

# ============================================================
# ★ 动态样本数
# ============================================================
n_samples_raw      <- length(seurat_list_raw)
n_samples_filtered <- length(seurat_list_filtered)

message(sprintf("Samples in raw list: %d", n_samples_raw))
message(sprintf("Samples in filtered list: %d", n_samples_filtered))

# ============================================================
# ★ 保存：原始（QC 过滤前）对象列表
# ============================================================
save_file_raw <- file.path(projectPath, "Data",
                           make_filename(n_samples_raw,
                                         "AfterEmptyDrops_AmbientRNA_BeforeQC",
                                         prefix = "01a_"))
save(seurat_list_raw, file = save_file_raw)
message(sprintf("Saved raw list -> %s", save_file_raw))

# ============================================================
# ★ 保存：QC 过滤后对象列表
# ============================================================
save_file_filtered <- file.path(projectPath, "Data",
                                make_filename(n_samples_filtered,
                                              "AfterEmptyDrops_AmbientRNA_AfterQC",
                                              prefix = "01b_"))
save(seurat_list_filtered, file = save_file_filtered)
message(sprintf("Saved filtered list -> %s", save_file_filtered))


# ============================================================
# 2. DoubletFinder 双细胞检测
# ============================================================

# 如果需要从中间恢复：
# load('/path/to/your/02_2021216_4Samples_AfterQC_AfterDF.Rdata')

# 确保 seurat_list 已经在全局环境中，且来自上一步 QC 过滤后的结果
# 此处假设 seurat_list 已经是 seurat_list_filtered 的集合
seurat_list <- seurat_list_filtered

umap_pdf <- file.path(projectPath, "Output",
                      paste0(date_tag, "_All_DoubletFinder_UMAP.pdf"))
pdf(umap_pdf, width = 10, height = 8)

doublet_summary <- list()

for (i in seq_along(seurat_list)) {
  obj_name    <- names(seurat_list)[i]
  current_obj <- seurat_list[[obj_name]]
  
  message(sprintf("\nDoubletFinder %d/%d: %s (%d cells)",
                  i, length(seurat_list), obj_name, ncol(current_obj)))
  
  # ---------- 基础预处理 ----------
  current_obj <- current_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE, seed.use = 42) %>%
    RunUMAP(dims = PC_DIMS, verbose = FALSE, seed.use = 42)
  
  # ---------- DoubletFinder 参数优化 ----------
  sweep.res   <- paramSweep(current_obj, PCs = PC_DIMS, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn       <- find.pK(sweep.stats)
  
  # BCmetric vs pK 图
  p_bcmetric <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(x = "pK Value", y = "BCmetric",
         title = paste0(obj_name, ": BCmetric vs. pK"))
  ggsave(file.path(projectPath, "Output", paste0(obj_name, "_BCmetric_vs_pK.png")),
         plot = p_bcmetric, width = 8, height = 6, dpi = 300)
  
  # 选择最佳 pK 值
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message(sprintf("  Best pK = %s", pK))
  
  # 动态估算预期双细胞数（10x 经验值）
  nCells              <- ncol(current_obj)
  doublet_rate_expected <- nCells * 0.8 / 1000 / 100
  nExp_poi            <- max(round(doublet_rate_expected * nCells), 1)
  message(sprintf("  Expected doublets (suggested): %d (rate ~%.2f%%)",
                  nExp_poi, doublet_rate_expected * 100))
  
  # ---------- 运行 DoubletFinder ----------
  current_obj <- doubletFinder(
    current_obj,
    PCs  = PC_DIMS,
    pN   = 0.25,
    pK   = pK,
    nExp = nExp_poi,
    sct  = FALSE
  )
  
  # 提取分类结果
  df.class <- paste("DF.classifications_0.25", pK, nExp_poi, sep = "_")
  current_obj$Doublet_Classification <- current_obj@meta.data[[df.class]]
  
  # 统计
  doublet_count <- sum(current_obj$Doublet_Classification == "Doublet")
  total_cells   <- ncol(current_obj)
  doublet_rate  <- round(doublet_count / total_cells * 100, 2)
  
  # ---------- 可视化 ----------
  doublet_plot <- DimPlot(current_obj,
                          group.by  = "Doublet_Classification",
                          reduction = "umap",
                          cols = c("Doublet" = "red", "Singlet" = "gray90")) +
    labs(title = paste0(obj_name, " Doublet Detection")) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("Doublet Count: ", doublet_count, "  (", doublet_rate, "%)"),
             hjust = 1.1, vjust = 1.1,
             size = 4, color = "darkred", fontface = "bold") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  print(doublet_plot)
  
  ggsave(file.path(projectPath, "Output", "DoubletFinder",
                   paste0(obj_name, "_DoubletFinder_UMAP.png")),
         plot = doublet_plot, width = 10, height = 8, dpi = 300)
  
  # ---------- 过滤双细胞 ----------
  singlet_cells <- colnames(current_obj)[current_obj$Doublet_Classification == "Singlet"]
  current_obj <- subset(current_obj, cells = singlet_cells) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE, seed.use = 42)
  
  message(sprintf("  [%s] After DoubletFinder: %d -> %d cells",
                  obj_name, total_cells, ncol(current_obj)))
  
  seurat_list[[obj_name]] <- current_obj
  
  doublet_summary[[obj_name]] <- data.frame(
    Sample       = obj_name,
    Total_Cells  = total_cells,
    Doublets     = doublet_count,
    Doublet_Rate = doublet_rate,
    Singlets_Kept= ncol(current_obj),
    Best_pK      = pK,
    stringsAsFactors = FALSE
  )
}

dev.off()
message(paste("Doublet UMAP plots saved to:", umap_pdf))

# 将 DoubletFinder 的摘要追加到 qc_summary_df 中
# 假设 qc_summary_df 已经在你的工作区，且包含前面阶段的统计
# 如果没有，可以先创建空的 qc_summary_df，然后用 rbind 添加
if (!exists("qc_summary_df")) {
  qc_summary_df <- data.frame(
    Sample = character(),
    Raw_Barcodes = integer(),
    After_EmptyDrops = integer(),
    Empty_Rate_Pct = numeric(),
    Ambient_Method = character(),
    Ambient_Rate_Pct = numeric(),
    Cells_Before_QC = integer(),
    Removed_LowFeature = integer(),
    Removed_HighFeature = integer(),
    Removed_LowCount = integer(),
    Removed_HighMT = integer(),
    Removed_HighHb = integer(),
    Cells_After_QC = integer(),
    Doublets = integer(),
    Doublet_Rate = numeric(),
    Singlets_Kept = integer(),
    Best_pK = numeric(),
    stringsAsFactors = FALSE
  )
}
# 将当前阶段的 doublet_summary 合并到 qc_summary_df
doublet_summary_df <- do.call(rbind, doublet_summary)
rownames(doublet_summary_df) <- NULL

# 合并 QC 和 Doublet 结果
qc_doublet_merged <- qc_summary_df %>%
  left_join(doublet_summary_df, by = "Sample")
  
# 保存 doublet 汇总
write.csv(
  qc_doublet_merged,
  file.path(projectPath, "Output",
            paste0(date_tag, "_QC_Doublet_Merged.csv")),
  row.names = FALSE
)
message("\n========== Doublet Summary ==========")
print(doublet_summary_df)

# 保存 DoubletFinder 后的对象
n_samples_df <- length(seurat_list)
save_file_df <- file.path(projectPath, "Data",
                          make_filename(n_samples_df,
                                        "AfterQC_AfterDoubletFinder",
                                        prefix = "02_"))
save(seurat_list, file = save_file_df)

message(sprintf("\nProcess Complete! Final list saved to: %s", save_file_df))
message(sprintf("Summary table saved to: %s", final_output_path))


# ============================================================
# 3. 合并样本（未整合）
# ============================================================
if (length(seurat_list) > 1) {
  unintegrated <- merge(seurat_list[[1]],
                        y = seurat_list[2:length(seurat_list)],
                        project = "Merged_Unintegrated")
} else {
  unintegrated <- seurat_list[[1]]
}

unintegrated <- unintegrated %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, seed.use = 42) %>%
  FindNeighbors(dims = PC_DIMS) %>%
  RunUMAP(dims = PC_DIMS, reduction.name = "umap_unintegrated", seed.use = 42)

p_unintegrated <- DimPlot(unintegrated, reduction = "umap_unintegrated",
                          group.by = "orig.ident") +
  ggtitle("Unintegrated UMAP")
ggsave(file.path(projectPath, "Output", "Unintegrated_UMAP_by_sample.png"),
       plot = p_unintegrated, width = 10, height = 8, dpi = 300)

save(unintegrated,
     file = file.path(projectPath, "Data", "03_2021216_Samples_MergedObj_Unintegrated.Rdata"))

message("\n========================================")
message("Pipeline completed successfully!")
message("========================================")

# ============================================================
# 3. 细胞周期评分（Cell Cycle Scoring）
# ============================================================
# 注意：Seurat 自带的是人类基因名
# 犬的基因名通常也是大写（和人一样），所以可以直接用

# Seurat 内置的 S 期和 G2M 期基因
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# 检查基因名匹配情况（犬的基因名可能有部分不匹配）
for (obj_name in names(seurat_list)) {
  current_obj <- seurat_list[[obj_name]]

  # 查看匹配率
  s_matched   <- intersect(s.genes,   rownames(current_obj))
  g2m_matched <- intersect(g2m.genes, rownames(current_obj))
  message(sprintf("  [%s] Cell cycle genes: S-phase %d/%d matched, G2M %d/%d matched",
                  obj_name,
                  length(s_matched),   length(s.genes),
                  length(g2m_matched), length(g2m.genes)))

  # 评分
  current_obj <- NormalizeData(current_obj, verbose = FALSE)
  current_obj <- CellCycleScoring(current_obj,
                                   s.features   = s.genes,
                                   g2m.features = g2m.genes,
                                   set.ident    = FALSE)

  # 查看分布
  cc_table <- table(current_obj$Phase)
  message(sprintf("  [%s] Cell cycle: G1=%d, S=%d, G2M=%d",
                  obj_name,
                  cc_table["G1"], cc_table["S"], cc_table["G2M"]))

  # 可视化
  p_cc <- DimPlot(current_obj %>%
                    RunUMAP(dims = PC_DIMS, verbose = FALSE, seed.use = 42),
                  group.by = "Phase",
                  cols = c("G1" = "#E8E8E8", "S" = "#2196F3", "G2M" = "#FF5722")) +
    ggtitle(paste0(obj_name, ": Cell Cycle Phase"))
  ggsave(file.path(projectPath, "Output","CellCycle",
                   paste0(obj_name, "_CellCycle_UMAP.png")),
         plot = p_cc, width = 9, height = 7, dpi = 300)

  # S.Score 和 G2M.Score 的小提琴图
  p_cc_vln <- VlnPlot(current_obj,
                       features = c("S.Score", "G2M.Score"),
                       pt.size = 0, ncol = 2) +
    plot_annotation(title = paste0(obj_name, ": Cell Cycle Scores"))
  ggsave(file.path(projectPath, "Output", "CellCycle",
                   paste0(obj_name, "_CellCycle_VlnPlot.png")),
         plot = p_cc_vln, width = 10, height = 5, dpi = 300)

  seurat_list[[obj_name]] <- current_obj
}

# 保存包含细胞周期信息的对象
save_file_cc <- file.path(projectPath, "Data",
                          make_filename(length(seurat_list),
                                        "AfterQC_DF_CellCycle",
                                        prefix = "03_"))
save(seurat_list, file = save_file_cc)
message(sprintf("Saved with cell cycle -> %s", save_file_cc))

# ============================================================
# 4. 基因级过滤（可选但推荐）
# ============================================================
# 去除只在极少细胞中表达的基因
MIN_CELLS_PER_GENE <- 3

for (obj_name in names(seurat_list)) {
  current_obj <- seurat_list[[obj_name]]
  genes_before <- nrow(current_obj)

  # 统计每个基因在多少个细胞中表达
  gene_cell_counts <- Matrix::rowSums(GetAssayData(current_obj, slot = "counts") > 0)
  genes_keep <- names(gene_cell_counts[gene_cell_counts >= MIN_CELLS_PER_GENE])

  current_obj <- subset(current_obj, features = genes_keep)
  genes_after <- nrow(current_obj)

  message(sprintf("  [%s] Gene filtering: %d -> %d genes (removed %d)",
                  obj_name, genes_before, genes_after, genes_before - genes_after))

  seurat_list[[obj_name]] <- current_obj
}

# ============================================================
# 5. 合并所有样本（未整合）
# ============================================================
if (length(seurat_list) > 1) {
  merged_obj <- merge(seurat_list[[1]],
                      y = seurat_list[2:length(seurat_list)],
                      add.cell.ids = names(seurat_list),
                      project = "Merged_AllSamples")
} else {
  merged_obj <- seurat_list[[1]]
}

message(sprintf("Merged object: %d cells x %d genes", ncol(merged_obj), nrow(merged_obj)))

# ============================================================
# 6. 合并后整体 QC 对比可视化
# ============================================================
# 按样本对比各 QC 指标
# 假设 merged_obj 已经存在，date_tag、projectPath、Outputs 路径等也已设置

# 要输出的特征列表，逐个输出到同一个 PDF 的单独页面
features_to_plot <- c("nFeature_RNA", "nCount_RNA",
                      "percent.mt", "percent.hb",
                      "percent.ribo", "log10GenesPerUMI")
# 输出 PDF，每页一个特征
pdf_file <- file.path(projectPath, "Output",
                      paste0(date_tag, "_MergedSamples_QC_Comparison_PerFeature.pdf"))
# 打开 PDF 设备
pdf(pdf_file, width = 10, height = 6)
for (feat in features_to_plot) {
  p <- VlnPlot(merged_obj,
               features = feat,
               group.by = "orig.ident",
               pt.size = 0) +
       ggtitle(paste0("QC Metric: ", feat)) +
       NoLegend()
  print(p)  # 每次在同一 PDF 的新页输出
}
dev.off()
message(sprintf("Saved per-feature ViolinPlots to: %s", pdf_file))

# ---------- 未整合 UMAP ----------
merged_obj <- merged_obj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, seed.use = 42) %>%
  FindNeighbors(dims = PC_DIMS, verbose = FALSE) %>%
  RunUMAP(dims = PC_DIMS, reduction.name = "umap_unintegrated", seed.use = 42)

p_unintegrated <- DimPlot(merged_obj, reduction = "umap_unintegrated",
                          group.by = "orig.ident") +
  ggtitle("Unintegrated UMAP (All Samples)")
ggsave(file.path(projectPath, "Output",
                 paste0(date_tag, "_Unintegrated_UMAP.png")),
       plot = p_unintegrated, width = 10, height = 8, dpi = 300)

# 保存合并对象
save_file_merged <- file.path(projectPath, "Data",
                              make_filename(length(seurat_list),
                                            "Merged_Unintegrated",
                                            prefix = "04_"))
save(merged_obj, file = save_file_merged)
message(sprintf("Saved merged object -> %s", save_file_merged))

# ============================================================
# ★ 最终 Pipeline 完成
# ============================================================
message("\n========================================")
message("★ Complete QC Pipeline Finished!")
message("========================================")
message(sprintf("Date: %s", Sys.Date()))
message(sprintf("Samples processed: %d", length(seurat_list)))
message(sprintf("Final merged: %d cells x %d genes", ncol(merged_obj), nrow(merged_obj)))
message("\nSaved files:")
message(sprintf("  01a: %s (before QC)", basename(save_file_raw)))
message(sprintf("  01b: %s (after QC)", basename(save_file_filtered)))
message(sprintf("  02:  %s (after DoubletFinder)", basename(save_file_df)))
message(sprintf("  03:  %s (with CellCycle)", basename(save_file_cc)))
message(sprintf("  04:  %s (merged)", basename(save_file_merged)))
message("========================================")
