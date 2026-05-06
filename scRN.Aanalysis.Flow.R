# ==============================================
# 1. 包加载与全局设置（修正原代码错误+新增依赖包）
# ==============================================
# 基础分析包
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
# 新增：Harmony整合、细胞注释、可视化依赖包
library(harmony)
library(SingleR)
library(celldex)
library(pheatmap)
library(ggrepel)

# 全局随机种子（保证全流程可重复）
set.seed(42)

# 路径设置
projectPath <- "/mnt2/wanggd_group/zjj/ZengMin20260310"
setwd(projectPath)

# 确保输出目录存在
dir.create(file.path(projectPath, "Output"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "EmptyDrops"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "AmbientRNA"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "DoubletFinder"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "CellCycle"), showWarnings = FALSE, recursive = TRUE)
# 新增：整合、聚类、注释相关输出目录
dir.create(file.path(projectPath, "Output", "HarmonyIntegration"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "Clustering"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Output", "CellAnnotation"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(projectPath, "Data"), showWarnings = FALSE, recursive = TRUE)

# 动态日期标签
date_tag <- format(Sys.Date(), "%y%m%d")
message(sprintf("Date tag: %s", date_tag))

# ==============================================
# 2. 核心参数集中管理（修正原代码重复赋值/语法错误）
# ==============================================
# QC阈值
MIN_FEATURES <- 500
MAX_FEATURES <- 5000
MAX_MT_PERCENT <- 10
MAX_HB_PERCENT <- 5
MIN_COUNTS <- 1000
MIN_NOVELTY <- 0.85
# 降维聚类参数
PC_DIMS <- 1:20
CLUSTER_RESOLUTIONS <- seq(0.1, 1, by = 0.05) # 多分辨率梯度
# 物种相关基因Pattern（人源默认，小鼠请切换注释行）
# 人源
MT_PATTERN <- "^MT-"
HB_PATTERN <- "^HBA|^HBB"
RIBO_PATTERN <- "^RPS|^RPL"
# 小鼠（如需使用请注释上面人源行，取消下面注释）
# MT_PATTERN <- "^mt-"
# HB_PATTERN <- "^Hb[ab]-"
# RIBO_PATTERN <- "^Rps|^Rpl"

# EmptyDrops参数
EMPTY_DROPS_FDR <- 0.01
EMPTY_DROPS_LOWER <- 100
# DecontX参数
DECONTX_CONTAMINATION_THRESHOLD <- 0.5
# 基因过滤参数
MIN_CELLS_PER_GENE <- 3
# 标记基因筛选参数
MARKER_LOGFC <- 0.25
MARKER_MIN_PCT <- 0.25

# 动态文件名生成函数
make_filename <- function(n_samples, suffix, prefix = "", ext = "Rdata") {
  fname <- paste0(prefix, date_tag, "_", n_samples, "Samples_", suffix, ".", ext)
  return(fname)
}

# Cell Ranger细胞数读取函数（保留原功能）
get_cellranger_cells <- function(sample_dir) {
  stats_path1 <- file.path(sample_dir, "outs", "stats.txt")
  metrics_path <- file.path(sample_dir, "outs", "metrics.json")
  n_cells <- NA
  if (file.exists(stats_path1)) {
    stat_lines <- readLines(stats_path1, warn = FALSE)
    tc_line <- grep("Total Cells|Total\\s+Cells|Cells", stat_lines, value = TRUE)
    if (length(tc_line) > 0) {
      nums <- as.numeric(gsub("[^0-9]", "", tc_line[1]))
      if (!is.na(nums) && length(nums) > 0) n_cells <- nums
    }
  }
  if (is.na(n_cells) && file.exists(metrics_path)) {
    json <- tryCatch(jsonlite::read_json(metrics_path), error = function(e) NULL)
    if (!is.null(json)) {
      if (!is.null(json$summary) && !is.null(json$summary$nCells)) {
        n_cells <- as.numeric(json$summary$nCells)
      } else if (!is.null(json$nCells)) {
        n_cells <- as.numeric(json$nCells)
      }
    }
  }
  if (is.na(n_cells)) {
    alt <- file.path(sample_dir, "outs", "cell_barcodes.tsv")
    if (file.exists(alt)) n_cells <- length(readLines(alt))
  }
  return(as.integer(n_cells))
}

# ==============================================
# 3. 原QC流程（修正语法/拼写错误，保证流程完整）
# ==============================================
# 样本读取与循环处理
samples <- list.files(file.path(projectPath, "Data", "QuantifyRawdata"))
seurat_list_raw <- list()
seurat_list_filtered <- list()
sample_ids <- c()
qc_summary <- list()

for (sample in samples) {
  tryCatch({
    message(paste("\n========================================"))
    message(paste("Processing sample:", sample))
    message(paste("========================================"))
    
    # 路径设置
    base_path <- file.path(projectPath, "Data", "QuantifyRawdata", sample)
    filtered_path <- file.path(base_path, "filtered_cell_gene_matrix")
    raw_path <- file.path(base_path, "raw_cell_gene_matrix")
    
    if (!dir.exists(filtered_path)) {
      warning(paste("Filtered matrix not found:", filtered_path))
      next
    }
    has_raw <- dir.exists(raw_path)
    message(sprintf("Raw matrix available: %s", has_raw))
    
    # 读取10X数据
    filtered_data <- Read10X(data.dir = filtered_path)
    if (has_raw) raw_data <- Read10X(data.dir = raw_path)

    # ======================
    # 步骤A：空细胞检测EmptyDrops
    # ======================
    if (has_raw) {
      message("Running emptyDrops on raw matrix...")
      set.seed(42)
      e.out <- emptyDrops(raw_data, lower = EMPTY_DROPS_LOWER)
      is_cell <- e.out$FDR <= EMPTY_DROPS_FDR & !is.na(e.out$FDR)
      n_empty <- sum(!is_cell, na.rm = TRUE)
      n_real_cells <- sum(is_cell, na.rm = TRUE)
      n_total_barcodes <- nrow(e.out)
      empty_rate <- round((1 - n_real_cells / n_total_barcodes) * 100, 2)
      message(sprintf("EmptyDrops: %d total barcodes, %d real cells, %.2f%% empty", 
                      n_total_barcodes, n_real_cells, empty_rate))
      
      # 空细胞检测可视化
      ed_df <- data.frame(
        total_UMI = e.out$Total,
        neg_log10_FDR = -log10(e.out$FDR),
        is_cell = is_cell
      )
      ed_df <- ed_df[!is.na(ed_df$neg_log10_FDR), ]
      p_empty <- ggplot(ed_df, aes(x = total_UMI, y = neg_log10_FDR, color = is_cell)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_x_log10() +
        scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray70"),
                           labels = c("TRUE" = "Cell", "FALSE" = "Empty"),
                           name = "Classification") +
        geom_hline(yintercept = -log10(EMPTY_DROPS_FDR), linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title = paste0(sample, ": EmptyDrops Detection"),
             subtitle = paste0("Real cells: ", n_real_cells, " | Empty rate: ", empty_rate, "%"),
             x = "Total UMI (log10)", y = "-log10(FDR)")
      ggsave(file.path(projectPath, "Output", "EmptyDrops", paste0(sample, "_emptyDrops.png")),
             plot = p_empty, width = 8, height = 6, dpi = 300)
      
      # 保留真实细胞barcode
      real_cell_barcodes <- rownames(e.out)[is_cell]
      common_barcodes <- intersect(colnames(filtered_data), real_cell_barcodes)
      filtered_data <- filtered_data[, common_barcodes]
      message(sprintf("After emptyDrops: %d cells retained", ncol(filtered_data)))
    } else {
      message("No raw matrix. Using UMI threshold for empty cell filtering.")
      col_sums <- Matrix::colSums(filtered_data)
      n_before <- ncol(filtered_data)
      keep_cells <- col_sums >= EMPTY_DROPS_LOWER
      n_empty_like <- sum(!keep_cells)
      empty_rate <- round(n_empty_like / n_before * 100, 2)
      sorted_umi <- sort(col_sums, decreasing = TRUE)
      knee_df <- data.frame(rank = seq_along(sorted_umi), UMI = sorted_umi)
      
      # 膝形图可视化
      p_knee <- ggplot(knee_df, aes(x = rank, y = UMI)) +
        geom_line(color = "steelblue") +
        scale_x_log10() + scale_y_log10() +
        geom_hline(yintercept = EMPTY_DROPS_LOWER, linetype = "dashed", color = "red") +
        annotate("text", x = max(knee_df$rank) * 0.5, y = EMPTY_DROPS_LOWER * 2,
                 label = paste0("Threshold: ", EMPTY_DROPS_LOWER, " UMI"),
                 color = "red", size = 4) +
        theme_classic() +
        labs(title = paste0(sample, ": Barcode Rank Plot"),
             subtitle = paste0("Removed: ", n_empty_like, " | Empty-like rate: ", empty_rate, "%"),
             x = "Barcode Rank (log10)", y = "Total UMI (log10)")
      ggsave(file.path(projectPath, "Output", "EmptyDrops", paste0(sample, "_kneePlot.png")),
             plot = p_knee, width = 8, height = 6, dpi = 300)
      
      filtered_data <- filtered_data[, keep_cells]
      n_real_cells <- ncol(filtered_data)
      n_total_barcodes <- n_before
      message(sprintf("After UMI threshold: %d -> %d cells", n_before, n_real_cells))
    }

    # ======================
    # 步骤B：环境RNA去除
    # ======================
    if (has_raw) {
      # 有raw矩阵用SoupX
      message("Running SoupX...")
      temp_obj <- CreateSeuratObject(counts = filtered_data, project = sample) %>%
        NormalizeData(verbose = FALSE) %>%
        FindVariableFeatures(verbose = FALSE) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(verbose = FALSE, seed.use = 42) %>%
        FindNeighbors(dims = PC_DIMS, verbose = FALSE) %>%
        FindClusters(resolution = 0.8, verbose = FALSE) %>%
        RunUMAP(dims = PC_DIMS, verbose = FALSE, seed.use = 42)
      
      sc <- SoupChannel(tod = raw_data, toc = filtered_data)
      clusters <- setNames(as.character(Idents(temp_obj)), colnames(temp_obj))
      sc <- setClusters(sc, clusters)
      sc <- setDR(sc, Embeddings(temp_obj, "umap"), reductName = "UMAP")
      sc <- autoEstCont(sc, verbose = FALSE)
      contamination_fraction <- sc$fit$rhoEst
      mean_contamination <- contamination_fraction
      message(sprintf("SoupX contamination: %.4f (%.2f%%)", contamination_fraction, contamination_fraction * 100))
      
      corrected_counts <- adjustCounts(sc)
      
      # SoupX诊断图
      png(file.path(projectPath, "Output", "AmbientRNA", paste0(sample, "_SoupX_diagnostics.png")),
          width = 1200, height = 800, res = 150)
      top_genes <- head(sort(Matrix::rowSums(filtered_data), decreasing = TRUE), 6)
      top_gene_name <- names(top_genes)[1]
      if (!is.null(top_gene_name)) {
        tryCatch(plotChangeMap(sc, geneSet = top_gene_name, DR = "UMAP"),
                 error = function(e) message("SoupX plotChangeMap failed: ", conditionMessage(e)))
      }
      dev.off()
      
      # 创建Seurat对象
      seurat_obj <- CreateSeuratObject(counts = corrected_counts, project = sample, min.features = 0, min.cells = 0)
      seurat_obj$ambient_contamination <- contamination_fraction
      seurat_obj$ambient_method <- "SoupX"
      rm(temp_obj, sc); gc()
    } else {
      # 无raw矩阵用DecontX
      message("Running DecontX...")
      sce <- SingleCellExperiment(assays = list(counts = filtered_data))
      set.seed(42)
      sce <- decontX(sce)
      contamination_scores <- colData(sce)$decontX_contamination
      mean_contamination <- mean(contamination_scores, na.rm = TRUE)
      median_contamination <- median(contamination_scores, na.rm = TRUE)
      message(sprintf("DecontX - Mean: %.2f%%, Median: %.2f%%", mean_contamination * 100, median_contamination * 100))
      
      # 污染分布可视化
      contam_df <- data.frame(contamination = contamination_scores)
      p_contam <- ggplot(contam_df, aes(x = contamination)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
        geom_vline(xintercept = mean_contamination, linetype = "dashed", color = "red", linewidth = 1) +
        geom_vline(xintercept = DECONTX_CONTAMINATION_THRESHOLD, linetype = "dotted", color = "darkred", linewidth = 1) +
        annotate("text", x = mean_contamination + 0.05, y = Inf,
                 label = paste0("Mean: ", round(mean_contamination * 100, 1), "%"),
                 vjust = 2, color = "red", size = 4) +
        annotate("text", x = DECONTX_CONTAMINATION_THRESHOLD + 0.05, y = Inf,
                 label = paste0("Threshold: ", DECONTX_CONTAMINATION_THRESHOLD * 100, "%"),
                 vjust = 4, color = "darkred", size = 3.5) +
        theme_classic() +
        labs(title = paste0(sample, ": DecontX Contamination Distribution"),
             x = "Contamination Fraction", y = "Number of Cells")
      ggsave(file.path(projectPath, "Output", "AmbientRNA", paste0(sample, "_DecontX_contamination_hist.png")),
             plot = p_contam, width = 8, height = 6, dpi = 300)
      
      # UMAP可视化污染
      sce <- scater::runPCA(sce, exprs_values = "decontXcounts", seed = 42)
      sce <- scater::runUMAP(sce, dimred = "PCA")
      umap_df <- data.frame(
        UMAP1 = reducedDim(sce, "UMAP")[, 1],
        UMAP2 = reducedDim(sce, "UMAP")[, 2],
        contamination = contamination_scores
      )
      p_umap_contam <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = contamination)) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_viridis_c(option = "inferno", direction = -1, name = "Contamination") +
        theme_classic() +
        labs(title = paste0(sample, ": DecontX Contamination on UMAP"))
      ggsave(file.path(projectPath, "Output", "AmbientRNA", paste0(sample, "_DecontX_UMAP.png")),
             plot = p_umap_contam, width = 9, height = 7, dpi = 300)
      
      # 创建Seurat对象
      corrected_counts <- decontXcounts(sce)
      seurat_obj <- CreateSeuratObject(counts = corrected_counts, project = sample, min.features = 0, min.cells = 0)
      seurat_obj$ambient_contamination <- contamination_scores
      seurat_obj$ambient_method <- "DecontX"
      seurat_obj$high_contamination <- contamination_scores > DECONTX_CONTAMINATION_THRESHOLD
      n_high_contam <- sum(seurat_obj$high_contamination)
      message(sprintf("High contamination cells (>%.0f%%): %d (%.2f%%)", 
                      DECONTX_CONTAMINATION_THRESHOLD * 100, n_high_contam, n_high_contam / ncol(seurat_obj) * 100))
      rm(sce); gc()
    }

    # ======================
    # 步骤C：QC指标计算
    # ======================
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = MT_PATTERN)
    seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = HB_PATTERN)
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = RIBO_PATTERN)
    seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    
    # 基因匹配检查
    mt_genes_found <- grep(MT_PATTERN, rownames(seurat_obj), value = TRUE)
    hb_genes_found <- grep(HB_PATTERN, rownames(seurat_obj), value = TRUE)
    ribo_genes_found <- grep(RIBO_PATTERN, rownames(seurat_obj), value = TRUE)
    message(sprintf("MT genes matched: %d", length(mt_genes_found)))
    message(sprintf("Hb genes matched: %d", length(hb_genes_found)))
    message(sprintf("Ribo genes matched: %d", length(ribo_genes_found)))
    if (length(hb_genes_found) == 0) warning("No Hb genes found! Check HB_PATTERN or species.")
    
    # 过滤前QC小提琴图
    vln_nCount <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + ggtitle("nCount_RNA") + NoLegend()
    vln_nFeature <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + ggtitle("nFeature_RNA") + NoLegend()
    vln_mt <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0) + ggtitle("percent.mt") + NoLegend()
    vln_hb <- VlnPlot(seurat_obj, features = "percent.hb", pt.size = 0) + ggtitle("percent.hb") + NoLegend()
    vln_ribo <- VlnPlot(seurat_obj, features = "percent.ribo", pt.size = 0) + ggtitle("percent.ribo") + NoLegend()
    vln_novelty <- VlnPlot(seurat_obj, features = "log10GenesPerUMI", pt.size = 0) + ggtitle("Novelty") + NoLegend()
    combined_vln_full <- (vln_nCount | vln_nFeature | vln_ribo) / (vln_mt | vln_hb | vln_novelty)
    ggsave(file.path(projectPath, "Output", paste0(sample, ".QC_6Metrics_BeforeFilter.png")),
           plot = combined_vln_full, width = 15, height = 10, dpi = 300)
    
    # 保存过滤前原始对象
    seurat_list_raw[[sample]] <- seurat_obj

    # ======================
    # 步骤D：QC过滤
    # ======================
    cells_before <- ncol(seurat_obj)
    seurat_obj_filtered <- subset(seurat_obj,
      subset = nFeature_RNA > MIN_FEATURES &
        nFeature_RNA < MAX_FEATURES &
        nCount_RNA > MIN_COUNTS &
        percent.mt < MAX_MT_PERCENT &
        percent.hb < MAX_HB_PERCENT &
        log10GenesPerUMI > MIN_NOVELTY
    )
    cells_after <- ncol(seurat_obj_filtered)
    
    # 过滤统计
    n_low_feature <- sum(seurat_obj$nFeature_RNA <= MIN_FEATURES)
    n_high_feature <- sum(seurat_obj$nFeature_RNA >= MAX_FEATURES)
    n_low_count <- sum(seurat_obj$nCount_RNA <= MIN_COUNTS)
    n_high_mt <- sum(seurat_obj$percent.mt >= MAX_MT_PERCENT)
    n_high_hb <- sum(seurat_obj$percent.hb >= MAX_HB_PERCENT)
    n_low_novelty <- sum(seurat_obj$log10GenesPerUMI <= MIN_NOVELTY, na.rm = TRUE)
    
    message(sprintf("[%s] QC filtering: %d -> %d cells (removed %d total)", 
                    sample, cells_before, cells_after, cells_before - cells_after))
    message(sprintf("- nFeature < %d: %d cells", MIN_FEATURES, n_low_feature))
    message(sprintf("- nFeature > %d: %d cells", MAX_FEATURES, n_high_feature))
    message(sprintf("- nCount < %d: %d cells", MIN_COUNTS, n_low_count))
    message(sprintf("- %%mt >= %d%%: %d cells", MAX_MT_PERCENT, n_high_mt))
    message(sprintf("- %%hb >= %d%%: %d cells", MAX_HB_PERCENT, n_high_hb))
    message(sprintf("- novelty <= %.1f: %d cells", MIN_NOVELTY, n_low_novelty))
    
    # 过滤后QC小提琴图
    vln_nCount_f <- VlnPlot(seurat_obj_filtered, features = "nCount_RNA", pt.size = 0) + ggtitle("nCount_RNA") + NoLegend()
    vln_nFeature_f <- VlnPlot(seurat_obj_filtered, features = "nFeature_RNA", pt.size = 0) + ggtitle("nFeature_RNA") + NoLegend()
    vln_mt_f <- VlnPlot(seurat_obj_filtered, features = "percent.mt", pt.size = 0) + ggtitle("percent.mt") + NoLegend()
    vln_hb_f <- VlnPlot(seurat_obj_filtered, features = "percent.hb", pt.size = 0) + ggtitle("percent.hb") + NoLegend()
    vln_ribo_f <- VlnPlot(seurat_obj_filtered, features = "percent.ribo", pt.size = 0) + ggtitle("percent.ribo") + NoLegend()
    vln_novelty_f <- VlnPlot(seurat_obj_filtered, features = "log10GenesPerUMI", pt.size = 0) + ggtitle("Novelty") + NoLegend()
    combined_vln_f <- (vln_nCount_f | vln_nFeature_f | vln_mt_f) / (vln_hb_f | vln_ribo_f | vln_novelty_f)
    ggsave(file.path(projectPath, "Output", paste0(sample, ".QC_ViolinPlots_AfterFilter.png")),
           plot = combined_vln_f, width = 15, height = 10, dpi = 300)
    
    # 存入过滤后列表
    seurat_list_filtered[[sample]] <- seurat_obj_filtered
    sample_ids <- c(sample_ids, sample)
    
    # 记录QC摘要
    qc_summary[[sample]] <- data.frame(
      Sample = sample,
      Raw_Barcodes = if (!is.na(get_cellranger_cells(base_path)) && get_cellranger_cells(base_path) > 0) get_cellranger_cells(base_path) else NA,
      After_EmptyDrops = n_real_cells,
      Empty_Rate_Pct = empty_rate,
      Ambient_Method = ifelse(has_raw, "SoupX", "DecontX"),
      Ambient_Rate_Pct = round(mean_contamination * 100, 2),
      Cells_Before_QC = cells_before,
      Removed_LowFeature = n_low_feature,
      Removed_HighFeature = n_high_feature,
      Removed_LowCount = n_low_count,
      Removed_HighMT = n_high_mt,
      Removed_HighHb = n_high_hb,
      Removed_LowNovelty = n_low_novelty,
      Cells_After_QC = cells_after,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(paste("Error processing sample:", sample))
    message(conditionMessage(e))
  })
}

# ======================
# QC结果汇总与保存
# ======================
qc_summary_df <- do.call(rbind, qc_summary)
if (!"Raw_Barcodes" %in% names(qc_summary_df)) qc_summary_df$Raw_Barcodes <- NA
write.csv(qc_summary_df, file.path(projectPath, "Output", paste0(date_tag, "_QC_Summary_Full.csv")), row.names = FALSE)
message("\n========== QC Summary ==========")
print(qc_summary_df)

# 保存QC前后对象
n_samples_raw <- length(seurat_list_raw)
n_samples_filtered <- length(seurat_list_filtered)
save_file_raw <- file.path(projectPath, "Data", make_filename(n_samples_raw, "AfterEmptyDrops_AmbientRNA_BeforeQC", prefix = "01a_"))
save(seurat_list_raw, file = save_file_raw)
save_file_filtered <- file.path(projectPath, "Data", make_filename(n_samples_filtered, "AfterEmptyDrops_AmbientRNA_AfterQC", prefix = "01b_"))
save(seurat_list_filtered, file = save_file_filtered)
message(sprintf("Samples in raw list: %d", n_samples_raw))
message(sprintf("Samples in filtered list: %d", n_samples_filtered))
message(sprintf("Saved raw list -> %s", save_file_raw))
message(sprintf("Saved filtered list -> %s", save_file_filtered))

# ==============================================
# 4. DoubletFinder双细胞去除（修正原代码）
# ==============================================
seurat_list <- seurat_list_filtered
umap_pdf <- file.path(projectPath, "Output/DoubletFinder", paste0(date_tag, "_All_DoubletFinder_UMAP.pdf"))
pdf(umap_pdf, width = 10, height = 8)
doublet_summary <- list()

for (i in seq_along(seurat_list)) {
  obj_name <- names(seurat_list)[i]
  current_obj <- seurat_list[[obj_name]]
  message(sprintf("\nDoubletFinder %d/%d: %s (%d cells)", i, length(seurat_list), obj_name, ncol(current_obj)))
  
  # 基础预处理
  current_obj <- current_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE, seed.use = 42) %>%
    RunUMAP(dims = PC_DIMS, verbose = FALSE, seed.use = 42)
  
  # pK参数优化
  sweep.res <- paramSweep(current_obj, PCs = PC_DIMS, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # pK选择可视化
  p_bcmetric <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_point() + geom_line() +
    theme_classic() +
    labs(x = "pK Value", y = "BCmetric", title = paste0(obj_name, ": BCmetric vs. pK"))
  ggsave(file.path(projectPath, "Output", "DoubletFinder", paste0(obj_name, "_BCmetric_vs_pK.png")),
         plot = p_bcmetric, width = 8, height = 6, dpi = 300)
  
  # 选择最优pK
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message(sprintf("Best pK = %s", pK))
  
  # 预期双细胞率（10X官方经验值：每1000个细胞增加0.8%双细胞率）
  nCells <- ncol(current_obj)
  doublet_rate_expected <- min(nCells * 0.8 / 1000 / 100, 0.08) # 最高不超过8%
  nExp_poi <- max(round(doublet_rate_expected * nCells), 1)
  message(sprintf("Expected doublets: %d (rate ~%.2f%%)", nExp_poi, doublet_rate_expected * 100))
  
  # 运行DoubletFinder
  current_obj <- doubletFinder(
    current_obj,
    PCs = PC_DIMS,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi,
    sct = FALSE
  )
  
  # 提取双细胞分类结果
  df.class <- paste("DF.classifications_0.25", pK, nExp_poi, sep = "_")
  current_obj$Doublet_Classification <- current_obj@meta.data[[df.class]]
  
  # 双细胞统计
  doublet_count <- sum(current_obj$Doublet_Classification == "Doublet")
  total_cells <- ncol(current_obj)
  doublet_rate <- round(doublet_count / total_cells * 100, 2)
  
  # 双细胞可视化
  doublet_plot <- DimPlot(current_obj, group.by = "Doublet_Classification", reduction = "umap",
                           cols = c("Doublet" = "red", "Singlet" = "gray90")) +
    labs(title = paste0(obj_name, " Doublet Detection")) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Doublet Count: ", doublet_count, " (", doublet_rate, "%)"),
             hjust = 1.1, vjust = 1.1, size = 4, color = "darkred", fontface = "bold") +
    theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "right")
  print(doublet_plot)
  ggsave(file.path(projectPath, "Output", "DoubletFinder", paste0(obj_name, "_DoubletFinder_UMAP.png")),
         plot = doublet_plot, width = 10, height = 8, dpi = 300)
  
  # 过滤双细胞，保留单细胞
  singlet_cells <- colnames(current_obj)[current_obj$Doublet_Classification == "Singlet"]
  current_obj <- subset(current_obj, cells = singlet_cells)
  message(sprintf("[%s] After DoubletFinder: %d -> %d cells", obj_name, total_cells, ncol(current_obj)))
  
  seurat_list[[obj_name]] <- current_obj
  doublet_summary[[obj_name]] <- data.frame(
    Sample = obj_name,
    Total_Cells = total_cells,
    Doublets = doublet_count,
    Doublet_Rate_Pct = doublet_rate,
    Singlets_Kept = ncol(current_obj),
    Best_pK = pK,
    stringsAsFactors = FALSE
  )
}
dev.off()
message(paste("Doublet UMAP plots saved to:", umap_pdf))

# 双细胞结果汇总与保存
doublet_summary_df <- do.call(rbind, doublet_summary)
qc_summary_df <- as.data.frame(qc_summary_df)
rownames(doublet_summary_df) <- NULL
qc_doublet_merged <- qc_summary_df %>% left_join(doublet_summary_df, by = "Sample")
write.csv(qc_doublet_merged, file.path(projectPath, "Output", paste0(date_tag, "_QC_Doublet_Merged.csv")), row.names = FALSE)
message("\n========== Doublet Summary ==========")
print(doublet_summary_df)

# 保存去双细胞后的对象
n_samples_df <- length(seurat_list)
save_file_df <- file.path(projectPath, "Data", make_filename(n_samples_df, "AfterQC_AfterDoubletFinder", prefix = "02_"))
save(seurat_list, file = save_file_df)
message(sprintf("Process Complete! Final list saved to: %s", save_file_df))

# ==============================================
# 5. 细胞周期评分（保留原功能，修正代码）
# ==============================================
# 细胞周期基因（人源默认，小鼠请切换注释）
# 人源
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
# 小鼠（如需使用请注释上面两行，取消下面注释）
# library(AnnotationDbi)
# library(org.Mm.eg.db)
# s.genes <- mapIds(org.Mm.eg.db, keys = cc.genes.updated.2019$s.genes, keytype = "SYMBOL", column = "SYMBOL", multiVals = "first") %>% na.omit() %>% as.character()
# g2m.genes <- mapIds(org.Mm.eg.db, keys = cc.genes.updated.2019$g2m.genes, keytype = "SYMBOL", column = "SYMBOL", multiVals = "first") %>% na.omit() %>% as.character()

for (obj_name in names(seurat_list)) {
  current_obj <- seurat_list[[obj_name]]
  # 基因匹配检查
  s_matched <- intersect(s.genes, rownames(current_obj))
  g2m_matched <- intersect(g2m.genes, rownames(current_obj))
  message(sprintf("[%s] Cell cycle genes: S-phase %d/%d matched, G2M %d/%d matched", 
                  obj_name, length(s_matched), length(s.genes), length(g2m_matched), length(g2m.genes)))
  
  # 细胞周期评分
  current_obj <- NormalizeData(current_obj, verbose = FALSE)
  current_obj <- CellCycleScoring(current_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  # 细胞周期分布统计
  cc_table <- table(current_obj$Phase)
  message(sprintf("[%s] Cell cycle: G1=%d, S=%d, G2M=%d", obj_name, cc_table["G1"], cc_table["S"], cc_table["G2M"]))
  
  # 可视化
  current_obj <- RunUMAP(current_obj, dims = PC_DIMS, verbose = FALSE, seed.use = 42)
  p_cc <- DimPlot(current_obj, group.by = "Phase", cols = c("G1" = "#E8E8E8", "S" = "#2196F3", "G2M" = "#FF5722")) +
    ggtitle(paste0(obj_name, ": Cell Cycle Phase"))
  ggsave(file.path(projectPath, "Output", "CellCycle", paste0(obj_name, "_CellCycle_UMAP.png")),
         plot = p_cc, width = 9, height = 7, dpi = 300)
  
  p_cc_vln <- VlnPlot(current_obj, features = c("S.Score", "G2M.Score"), pt.size = 0, ncol = 2) +
    plot_annotation(title = paste0(obj_name, ": Cell Cycle Scores"))
  ggsave(file.path(projectPath, "Output", "CellCycle", paste0(obj_name, "_CellCycle_VlnPlot.png")),
         plot = p_cc_vln, width = 10, height = 5, dpi = 300)
  
  seurat_list[[obj_name]] <- current_obj
}

# 保存细胞周期评分后的对象
save_file_cc <- file.path(projectPath, "Data", make_filename(length(seurat_list), "AfterQC_DF_CellCycle", prefix = "03_"))
save(seurat_list, file = save_file_cc)
message(sprintf("Saved with cell cycle -> %s", save_file_cc))

# ==============================================
# 6. 基因级过滤（保留原功能，修正代码）
# ==============================================
for (obj_name in names(seurat_list)) {
  current_obj <- seurat_list[[obj_name]]
  genes_before <- nrow(current_obj)
  # 保留在至少3个细胞中表达的基因
  gene_cell_counts <- Matrix::rowSums(GetAssayData(current_obj, slot = "counts") > 0)
  genes_keep <- names(gene_cell_counts[gene_cell_counts >= MIN_CELLS_PER_GENE])
  current_obj <- subset(current_obj, features = genes_keep)
  genes_after <- nrow(current_obj)
  message(sprintf("[%s] Gene filtering: %d -> %d genes (removed %d)", 
                  obj_name, genes_before, genes_after, genes_before - genes_after))
  seurat_list[[obj_name]] <- current_obj
}
save_file_cc <- file.path(projectPath, "Data", make_filename(length(seurat_list), "AfterQC_DF_CellCycle_GenesFiltered", prefix = "04_"))
save(seurat_list, file = save_file_cc)
message(sprintf("Saved with cell cycle -> %s", save_file_cc))

# ==============================================
# ========== 新增核心流程开始 ==========
# ==============================================
# ==============================================
# 7. 样本合并与Harmony批次整合（按orig.ident分组）
# ==============================================
message("\n========================================")
message("Starting Harmony Integration by orig.ident")
message("========================================")

# 多样本合并
if (length(seurat_list) > 1) {
  merged_obj <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.ids = names(seurat_list),
    project = "scRNAseq_Harmony"
  )
} else {
  merged_obj <- seurat_list[[1]]
}
message(sprintf("Merged object: %d cells x %d genes", ncol(merged_obj), nrow(merged_obj)))

# 合并后标准化与PCA（Harmony的输入基础）
merged_obj <- merged_obj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE, seed.use = 42)

# PCA碎石图，确认PC维度选择
p_elbow <- ElbowPlot(merged_obj, ndims = 50) +
  geom_vline(xintercept = max(PC_DIMS), linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "PCA Elbow Plot")
ggsave(file.path(projectPath, "Output", "HarmonyIntegration", paste0(date_tag, "_PCA_ElbowPlot.png")),
       plot = p_elbow, width = 8, height = 6, dpi = 300)

# 未整合UMAP（用于批次效应对比）
merged_obj <- RunUMAP(merged_obj, dims = PC_DIMS, reduction.name = "umap_unintegrated", seed.use = 42, verbose = FALSE)
p_unintegrated <- DimPlot(merged_obj, reduction = "umap_unintegrated", group.by = "orig.ident", raster = FALSE) +
  ggtitle("Unintegrated UMAP (by Sample)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(projectPath, "Output", "HarmonyIntegration", paste0(date_tag, "_Unintegrated_UMAP.png")),
       plot = p_unintegrated, width = 10, height = 8, dpi = 300)

# ============================================================
# 合并后整体 QC 对比可视化
# ============================================================
# 按样本对比各 QC 指标
# 假设 merged_obj 已经存在，date_tag、projectPath、Outputs 路径等也已设置
# 要输出的特征列表，逐个输出到同一个 PDF 的单独页面
library(Seurat)
library(ggplot2)
library(RColorBrewer)
# 1. 准备颜色方案 
# 如果样本数超过12个，Set3 或 Paired 会循环使用，颜色更温和
sample_count <- length(unique(merged_obj$orig.ident))
my_colors <- colorRampPalette(brewer.pal(12, "Paired"))(sample_count)
# 设置输出路径
pdf_file <- file.path(projectPath, "Output", 
                      paste0(date_tag, "_Final_QC_Report.pdf"))
pdf(pdf_file, width = 16, height = 8)
# --- 第一部分：优化后的小提琴图 (温和配色) ---
features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                      "percent.hb", "percent.ribo", "log10GenesPerUMI")
for (feat in features_to_plot) {
  p <- VlnPlot(merged_obj, 
               features = feat, 
               group.by = "orig.ident", 
               pt.size = 0, # 必须为0，否则点太多会黑乎乎一片
               cols = my_colors) +
    theme_bw() + # 使用带网格的白色背景，更专业
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      panel.grid.major.x = element_blank(), # 关闭垂直网格线，减少晃眼感
      legend.position = "none"
    ) +
    labs(x = "Sample ID", y = feat) +
    ggtitle(paste0("Refined QC: ", feat))
  print(p)
}
# --- 第二部分：散点图与阈值标记 ---
# 设定建议的阈值（请根据您的具体数据调整这些数值）
thresh_min_nFeature <- 500
thresh_max_nFeature <- 4000
thresh_max_mt <- 15  # 假设 15% 为线粒体上限
# nFeature vs nCount
p1 <- FeatureScatter(merged_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                     group.by = "orig.ident", pt.size = 0.5) +
      scale_color_viridis(discrete = TRUE, alpha = 0.3) +
      geom_hline(yintercept = c(thresh_min_nFeature, thresh_max_nFeature), 
                 linetype = "dashed", color = "red") +
      theme(legend.position = "none") +
      ggtitle("Cell Quality: nFeature vs nCount")

# nFeature vs percent.mt
p2 <- FeatureScatter(merged_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt", 
                     group.by = "orig.ident", pt.size = 0.5) +
      scale_color_viridis(discrete = TRUE, alpha = 0.3) +
      geom_hline(yintercept = thresh_max_mt, linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(thresh_min_nFeature, thresh_max_nFeature), 
                 linetype = "dashed", color = "red") +
      theme(legend.position = "none") +
      ggtitle("Cell Quality: nFeature vs Mitochondrial %")

# 打印散点图
print(p1)
print(p2)
dev.off()
message(sprintf("Refined QC plots and Scatter plots saved to: %s", pdf_file))

    
# ======================
# 核心：Harmony整合（按orig.ident去除批次效应）
# ======================
set.seed(42)
merged_obj <- RunHarmony(
  object = merged_obj,
  group.by.vars = "orig.ident", # 按样本分组去除批次
  reduction = "pca",
  reduction.save = "harmony",
  dims.use = PC_DIMS,
  verbose = TRUE,
  max.iter.harmony = 50
)

# 整合后UMAP/TSNE降维
merged_obj <- merged_obj %>%
  RunUMAP(reduction = "harmony", dims = PC_DIMS, reduction.name = "umap_harmony", seed.use = 42, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = PC_DIMS, reduction.name = "tsne_harmony", seed.use = 42, verbose = FALSE)

# 整合前后对比可视化
p_harmony_sample <- DimPlot(merged_obj, reduction = "umap_harmony", group.by = "orig.ident", raster = FALSE) +
  ggtitle("Harmony Integrated UMAP (by Sample)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p_compare <- p_unintegrated + p_harmony_sample
ggsave(file.path(projectPath, "Output", "HarmonyIntegration", paste0(date_tag, "_Integration_Comparison_UMAP.png")),
       plot = p_compare, width = 20, height = 8, dpi = 300)

# 保存整合后中间对象
save_file_harmony <- file.path(projectPath, "Data", make_filename(length(seurat_list), "Harmony_Integrated", prefix = "04_"))
save(merged_obj, file = save_file_harmony)
message(sprintf("Saved Harmony integrated object -> %s", save_file_harmony))
message("Harmony Integration Completed!")

# ==============================================
# 8. 多分辨率聚类与最优分辨率筛选
# ==============================================
message("\n========================================")
message("Starting Clustering & Resolution Selection")
message("========================================")

# 基于Harmony降维结果构建邻居图
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = PC_DIMS, verbose = FALSE)

# 多分辨率聚类（梯度分辨率，用于筛选最优值）
merged_obj <- FindClusters(
  merged_obj,
  resolution = CLUSTER_RESOLUTIONS,
  verbose = FALSE,
  random.seed = 42
)

# ======================
# 方法1：Clustree树状图评估分群稳定性
# ======================
p_clustree <- clustree(merged_obj, prefix = "RNA_snn_res.") +
  theme_classic() +
  labs(title = "Clustering Tree Across Resolutions")
ggsave(file.path(projectPath, "Output", "Clustering", paste0(date_tag, "_Clustree_Resolution.png")),
       plot = p_clustree, width = 12, height = 10, dpi = 300)

# ======================
# 方法2：轮廓系数评估聚类质量
# ======================
# 提取Harmony降维矩阵
harmony_embeddings <- Embeddings(merged_obj, "harmony")[, PC_DIMS]
dist_matrix <- dist(harmony_embeddings)  
sil_scores <- data.frame(Resolution = numeric(), Mean_Silhouette = numeric(),  
                          N_Clusters = integer())  
for (res in CLUSTER_RESOLUTIONS) {  
  col_name <- paste0("RNA_snn_res.", res)  
  clusters <- as.integer(merged_obj@meta.data[[col_name]])  
  n_clust  <- length(unique(clusters))  
  # Silhouette 要求至少 2 个 cluster  
  if (n_clust >= 2 & n_clust < ncol(merged_obj)) {  
    # 对大数据集采样计算（加速）  
    if (ncol(merged_obj) > 5000) {  
      set.seed(42)  
      sample_idx <- sample(seq_len(ncol(merged_obj)), 5000)  
      sil <- cluster::silhouette(clusters[sample_idx],  
                                  dist(harmony_embeddings[sample_idx, ]))  
    } else {  
      sil <- cluster::silhouette(clusters, dist_matrix)  
    }  
    mean_sil <- mean(sil[, 3])  
  } else {  
    mean_sil <- NA  
  }  
  sil_scores <- rbind(sil_scores,  
                       data.frame(Resolution = res,  
                                  Mean_Silhouette = mean_sil,  
                                  N_Clusters = n_clust))  
}  
# Silhouette 图  
p_sil <- ggplot(sil_scores, aes(x = Resolution, y = Mean_Silhouette)) +  
  geom_line(color = "steelblue", linewidth = 1) +  
  geom_point(color = "steelblue", size = 2.5) +  
  geom_text(aes(label = N_Clusters), vjust = -1, size = 3, color = "gray40") +  
  theme_classic() +  
  labs(title = "Mean Silhouette Score vs Resolution",  
       subtitle = "Numbers indicate cluster count",  
       x = "Resolution", y = "Mean Silhouette Score")  
ggsave(file.path(projectPath, "Output", "Clustering",  
                 paste0(date_tag, "_Silhouette_vs_Resolution.png")),  
       plot = p_sil, width = 10, height = 6, dpi = 300)  
# ======================
# 确定最优分辨率
# 选择规则：轮廓系数最高 + Clustree分群稳定无过度碎片化
# ======================
# 输出建议  
best_res <- sil_scores$Resolution[which.max(sil_scores$Mean_Silhouette)]  
message(sprintf("  Best resolution by silhouette: %.1f (%d clusters, sil=%.3f)",  
                best_res,  
                sil_scores$N_Clusters[which.max(sil_scores$Mean_Silhouette)],  
                max(sil_scores$Mean_Silhouette, na.rm = TRUE)))  

# 保存 Silhouette 结果  
write.csv(sil_scores,  
          file.path(projectPath, "Output", "Clustering",  
                    paste0(date_tag, "_Silhouette_Scores.csv")),  
          row.names = FALSE)  
      
# 将最优分辨率的聚类设为默认ident      
best_res = 0.2
Idents(merged_obj) <- paste0("RNA_snn_res.", best_res)
merged_obj$seurat_clusters <- Idents(merged_obj)

# ★★★ 根据 Clustree + Silhouette + 生物学判断手动选择 ★★★  
# 先用 silhouette 推荐值，运行后根据结果调整  
FINAL_RESOLUTION <- best_res  
message(sprintf("  Using resolution: %.1f", FINAL_RESOLUTION))  

Idents(merged_obj) <- paste0("RNA_snn_res.", FINAL_RESOLUTION)  
merged_obj$seurat_clusters <- Idents(merged_obj)  
# 最终聚类 UMAP  
p_final_cluster <- DimPlot(merged_obj, reduction = "umap_harmony",  
                            label = TRUE, label.size = 5, pt.size = 0.1) +  
  ggtitle(paste0("Final Clustering (Resolution = ", FINAL_RESOLUTION,  
                 ", ", length(unique(Idents(merged_obj))), " clusters)")) +  
  theme(plot.title = element_text(size = 14, hjust = 0.5))  
ggsave(file.path(projectPath, "Output", "Clustering",  
                 paste0(date_tag, "_FinalClusters_UMAP.png")),  
       plot = p_final_cluster, width = 10, height = 8, dpi = 300)  
# tSNE  
p_final_tsne <- DimPlot(merged_obj, reduction = "tsne_harmony",  
                          label = TRUE, label.size = 5, pt.size = 0.1) +  
  ggtitle(paste0("Final Clustering (tSNE, Resolution = ", FINAL_RESOLUTION, ")"))  
ggsave(file.path(projectPath, "Output", "Clustering",  
                 paste0(date_tag, "_FinalClusters_tSNE.png")),  
       plot = p_final_tsne, width = 10, height = 8, dpi = 300)  

# QC 指标在 cluster 上的分布  
p_qc_cluster <- VlnPlot(merged_obj,  
                          features = c("nFeature_RNA", "nCount_RNA",  
                                       "percent.mt", "percent.ribo"),  
                          pt.size = 0, ncol = 2)  
ggsave(file.path(projectPath, "Output", "Clustering",  
                 paste0(date_tag, "_QC_per_Cluster.png")),  
       plot = p_qc_cluster, width = 16, height = 12, dpi = 300)  

# 各样本在 cluster 中的占比  
cluster_sample_table <- table(Idents(merged_obj), merged_obj$orig.ident)  
cluster_prop <- as.data.frame(prop.table(cluster_sample_table, margin = 2))  
colnames(cluster_prop) <- c("Cluster", "Sample", "Proportion")  

p_prop <- ggplot(cluster_prop, aes(x = Sample, y = Proportion, fill = Cluster)) +  
  geom_bar(stat = "identity", position = "fill") +  
  theme_classic() +  
  labs(title = "Cluster Proportion by Sample", y = "Proportion") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave(file.path(projectPath, "Output", "Clustering",  
                 paste0(date_tag, "_ClusterProportion_bySample.png")),  
       plot = p_prop, width = 10, height = 7, dpi = 300)  
      
# 保存聚类后对象
save_file_clustered <- file.path(projectPath, "Data",  
                                  make_filename(length(unique(merged_obj$orig.ident)),  
                                                "Harmony_Clustered", prefix = "07_"))  
save(merged_obj, file = save_file_clustered)  
message(sprintf("Saved clustered object -> %s", save_file_clustered))  


# ==============================================
# 9. 标记基因筛选（FindAllMarkers）
# ==============================================
message("\n========================================")
message("Starting Marker Gene Identification")
message("========================================")


# ==========================================
# 检查当前 Assay 数据结构
message("Current data layers:")
print(Layers(merged_obj, assay = "RNA"))
# 执行 JoinLayers，将 split 的样本层重新合并
# 这样 FindMarkers 才能在所有细胞间进行统计比较
merged_obj <- JoinLayers(merged_obj, assay = "RNA")
message("After JoinLayers:")
print(Layers(merged_obj, assay = "RNA"))      
# 寻找每个cluster的标记基因（仅正相关，过滤低表达基因）
all_markers <- FindAllMarkers(
  merged_obj,
  only.pos = TRUE,
  min.pct = MARKER_MIN_PCT,
  logfc.threshold = MARKER_LOGFC,
  test.use = "wilcox",
  verbose = TRUE
)

# 按cluster和p值排序，保留显著结果
all_markers_filtered <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

# 保存标记基因表
write.csv(all_markers_filtered, file.path(projectPath, "Output", "CellAnnotation", paste0(date_tag, "_AllCluster_MarkerGenes.csv")), row.names = FALSE)

# 每个cluster取top10标记基因做热图
top_markers <- all_markers_filtered %>%
  filter(!grepl("^ENSG\\d", gene)) %>%   # 去除 ENSG 开头的基因
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC)
# 标记基因热图
p_heatmap <- DoHeatmap(merged_obj, features = top_markers$gene,
                        size = 3) +
  theme(axis.text.y = element_text(size = 6)) +
  ggtitle("Top Marker Genes per Cluster")
ggsave(file.path(projectPath, "Output", "CellAnnotation", paste0(date_tag, "_Top_MarkerGenes_Heatmap_.png")),
       plot = p_heatmap, width = 15, height = 12, dpi = 300)
p_marker_heatmap <- DoHeatmap(merged_obj, features = top_markers$gene, group.by = "seurat_clusters", raster = FALSE) +
  scale_fill_viridis_c(option = "inferno") +
  labs(title = "Top Marker Genes per Cluster")
ggsave(file.path(projectPath, "Output", "CellAnnotation", paste0(date_tag, "_Top_MarkerGenes_Heatmap.png")),
       plot = p_marker_heatmap, width = 15, height = 12, dpi = 300)

message(sprintf("Total significant marker genes identified: %d", nrow(all_markers_filtered)))
message("Marker Gene Identification Completed!")

# ==============================================
# 10. 细胞类型注释（SingleR自动注释+手动验证）
# ==============================================
message("\n========================================")
message("Starting Cell Type Annotation")
message("========================================")
dir.create(file.path(projectPath, "Output", "CellAnnotation"),
           showWarnings = FALSE, recursive = FALSE)
# ======================
# 步骤1：加载参考数据集
# ======================
message("  Loading reference datasets...")

# 人源参考（根据物种选择，二选一即可）
ref_hpca      <- HumanPrimaryCellAtlasData()   # 最全面
ref_blueprint <- BlueprintEncodeData()           # 质量最高
ref_dice <- DatabaseImmuneCellExpressionData() # 可选：加入免疫细胞专用数据集

# 小鼠参考（如需使用请注释上面2行，取消下面注释）
# ref_mouse  <- MouseRNAseqData()
# ref_immgen <- ImmGenData()

# ★ 正确合并：取交集基因后再 cbind ★
common_genes_ref <- intersect(rownames(ref_hpca), rownames(ref_blueprint))
ref_combined <- cbind(ref_hpca[common_genes_ref, ], ref_blueprint[common_genes_ref, ])

#
ref_combined = ref_dice
message(sprintf("  Reference: %d genes, %d samples",
                nrow(ref_combined), ncol(ref_combined)))

# ======================
# 步骤2：SingleR 注释
# ======================
# 转为 SCE（SingleR 需要 log-normalized 数据）
sce_query <- as.SingleCellExperiment(merged_obj)

# ---- Cluster 级别注释（更稳健，推荐用于最终展示） ----
message("  Running SingleR at cluster level...")
singler_cluster <- SingleR(test   = sce_query,
                            ref    = ref_combined,
                            labels = ref_combined$label.fine,     #label.main,
                            clusters = merged_obj$seurat_clusters,
                            de.method = "wilcox")

message("  Cluster-level annotation results:")
print(data.frame(Cluster = rownames(singler_cluster),
                  Label   = singler_cluster$pruned.labels))

# ---- Cell 级别注释（精细，但可能有部分 NA） ----
message("  Running SingleR at cell level (may take a while)...")
singler_cell <- SingleR(test   = sce_query,
                          ref    = ref_combined,
                          labels = ref_combined$label.fine,     #label.main,
                          de.method = "wilcox")

# 释放大对象
rm(sce_query); gc()

# ======================
# 步骤3：写入 Seurat 对象（★核心修复区★）
# ======================
message("  Writing annotations to Seurat object...")

# ---- 3a. Cell-level 标签 ----
# pruned.labels 会有 NA（低信度被裁剪掉的），替换为 "Unknown"
cell_labels <- singler_cell$pruned.labels
cell_labels[is.na(cell_labels)] <- "Unknown"

# ★ 用 @meta.data 直接赋值，绕过 Seurat 的 cell name 校验 ★
merged_obj@meta.data$singler_cell_label <- cell_labels
merged_obj@meta.data$singler_cell_score <- apply(singler_cell$scores, 1, max)

message(sprintf("  Cell-level: %d Unknown out of %d cells",
                sum(cell_labels == "Unknown"), length(cell_labels)))

# ---- 3b. Cluster-level 标签映射到每个细胞 ----
cluster_labels <- singler_cluster$pruned.labels
cluster_labels[is.na(cluster_labels)] <- "Unknown"
names(cluster_labels) <- rownames(singler_cluster)
# 获取每个细胞所属的 cluster
cluster_id <- as.character(merged_obj$seurat_clusters)
# 检查名称是否匹配
message(sprintf("  singler_cluster rownames: %s",
                paste(rownames(singler_cluster), collapse = ", ")))
message(sprintf("  seurat_clusters levels:   %s",
                paste(sort(unique(cluster_id)), collapse = ", ")))

# 如果名称不匹配，按排序顺序强制对应
if (!all(unique(cluster_id) %in% names(cluster_labels))) {
  message("  ⚠ Cluster names mismatch! Re-mapping by sorted order...")
  sorted_clusters <- sort(as.numeric(unique(cluster_id)))
  names(cluster_labels) <- as.character(sorted_clusters)
}

# ★ 同样用 @meta.data 直接赋值 ★
merged_obj@meta.data$singler_cluster_label <- cluster_labels[cluster_id]
# 确认没有 NA
message(sprintf("  singler_cell_label NA count:    %d",
                sum(is.na(merged_obj$singler_cell_label))))
message(sprintf("  singler_cluster_label NA count: %d",
                sum(is.na(merged_obj$singler_cluster_label))))

# ======================
# 步骤4：SingleR 诊断图
# ======================
message("  Generating SingleR diagnostic plots...")
# Score Heatmap
png(file.path(projectPath, "Output", "CellAnnotation",
              paste0(date_tag, "_SingleR_ClusterScores_Heatmap.png")),
    width = 1600, height = 1000, res = 150)
plotScoreHeatmap(singler_cluster, show.labels = TRUE)
dev.off()

# Delta 分布图
png(file.path(projectPath, "Output", "CellAnnotation",
              paste0(date_tag, "_SingleR_DeltaDistribution.png")),
    width = 1400, height = 900, res = 150)
plotDeltaDistribution(singler_cell, ncol = 4)
dev.off()

# ======================
# 步骤5：UMAP 可视化注释结果
# ======================
message("  Generating annotation UMAP plots...")
# Cluster 级别
p_anno_cluster <- DimPlot(merged_obj, reduction = "umap_harmony",
                           group.by = "singler_cluster_label",
                           label = TRUE, label.size = 3.5, pt.size = 0.1,
                           repel = TRUE) +
  ggtitle("Cell Type (SingleR Cluster-level)") +
  theme(plot.title = element_text(hjust = 0.5))
# Cell 级别
p_anno_cell <- DimPlot(merged_obj, reduction = "umap_harmony",
                        group.by = "singler_cell_label",
                        label = TRUE, label.size = 3, pt.size = 0.1,
                        repel = TRUE) +
  ggtitle("Cell Type (SingleR Cell-level)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()
# Seurat cluster 编号对照
p_cluster_ref <- DimPlot(merged_obj, reduction = "umap_harmony",
                          group.by = "seurat_clusters",
                          label = TRUE, label.size = 4, pt.size = 0.1) +
  ggtitle("Seurat Clusters") +
  NoLegend()
# 三图拼接
p_annotation_combined <- p_cluster_ref | p_anno_cluster | p_anno_cell
ggsave(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "_Annotation_UMAP_Combined.png")),
       plot = p_annotation_combined, width = 28, height = 8, dpi = 300)
# 单独保存 cluster-level
ggsave(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "_Annotation_ClusterLevel_UMAP.png")),
       plot = p_anno_cluster, width = 12, height = 8, dpi = 300)

# 单独保存 cell-level
ggsave(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "_Annotation_CellLevel_UMAP.png")),
       plot = p_anno_cell, width = 12, height = 8, dpi = 300)

# ======================
# 步骤6：注释汇总表
# ======================
message("  Generating annotation summary table...")

# 用实际 cluster ID 统计
cluster_cell_counts <- table(merged_obj$seurat_clusters)

anno_summary <- data.frame(
  Cluster   = names(cluster_labels),
  Cell_Type = cluster_labels,
  N_Cells   = as.integer(cluster_cell_counts[names(cluster_labels)]),
  row.names = NULL
)

# 添加 delta（如果存在）
if (!is.null(singler_cluster$delta.next)) {
  anno_summary$Delta_Next <- singler_cluster$delta.next
}

# 按 cluster 数字排序
anno_summary <- anno_summary[order(as.numeric(anno_summary$Cluster)), ]

write.csv(anno_summary,
          file.path(projectPath, "Output", "CellAnnotation",
                    paste0(date_tag, "_Annotation_Summary.csv")),
          row.names = FALSE)

message("\n  ===== Annotation Summary =====")
print(anno_summary)

# ======================
# 步骤7：保存最终对象
# ======================
save_file_annotated <- file.path(projectPath, "Data",
                                  make_filename(
                                    length(unique(merged_obj$orig.ident)),
                                    "Harmony_Clustered_Annotated",
                                    prefix = "08_"))
save(merged_obj, file = save_file_annotated)
message(sprintf("  Saved annotated object -> %s", save_file_annotated))

message("\n========================================")
message("Cell Type Annotation Complete!")
message("========================================")

# --- 10.8 已知 Marker 基因验证（根据物种/组织调整） ---  
# ★★★ 以下 marker 基因为人类常见类型示例 ★★★  
# ★★★ 请根据你的实际物种和组织替换 ★★★  
canonical_markers <- list(
  # ==================== 总T细胞 ====================
  "T cells"                = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC1", "TRBC2", "CD28", "TRAT1"),
  # ==================== 幼稚T细胞（共同标记） ====================
  "Naive T cells"    = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),  

  # ==================== CD4+ T细胞亚群 ==================== 
  "CD4+ T naive"          = c("CD4", "CCR7", "SELL", "LEF1", "TCF7", "IL7R", "CD27", "CD28"),                     #CD4+ 幼稚T（核心：CD4 + naive标记）
  "CD4+ T memory"          = c("CD4", "IL7R", "S100A4", "ANXA1", "IL32", "GPR183", "CD69", "AQP3", "ITGB1"),
  "CD4+ T effector"        = c("CD4", "GZMA", "IFNG", "TNF", "CTLA4", "LAG3", "NKG7", "PRF1", "TBX21"),
  #"Treg (调节性T)"         = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF18", "TNFRSF4", "BATF", "LRRC32"),
  #"Th1 cells"              = c("CD4", "TBX21", "IFNG", "CXCR3", "IL12RB2", "TNF", "IL2", "STAT4", "CCR5"),
  #"Th2 cells"              = c("CD4", "GATA3", "IL4", "IL5", "IL13", "CCR4", "CRTH2", "STAT6", "IL17RB"),
  #"Th17 cells"             = c("CD4", "RORC", "IL17A", "IL17F", "CCR6", "IL23R", "STAT3", "RORA", "CD161"),
  # ==================== CD8+ T细胞亚群 ====================
  "CD8+ T naive"           = c("CD8A", "CD8B", "CCR7", "SELL", "LEF1", "TCF7", "IL7R", "S1PR1"),                  # CD8+ 幼稚T（核心：CD8A/B + naive标记）
  "CD8+ T memory"          = c("CD8A", "CD8B", "IL7R", "CD27", "CD28", "GZMK", "ITGA1", "S100A4", "EOMES"),
  "CD8+ T effector"        = c("CD8A", "CD8B", "GZMA", "GZMK", "GZMB", "PRF1", "NKG7", "IFNG", "TBX21"),
  "CD8+ T exhausted"       = c("CD8A", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "ENTPD1", "TOX", "EOMES"),
  # ==================== NK细胞亚群 ====================
  "NK cells"               = c("NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "FCGR3A", "GZMB", "KLRF1", "SPON2"),
  "CD56bright NK"          = c("NCAM1", "XCL1", "XCL2", "GZMK", "IL7R", "CD27", "CD28", "IFNG", "KLRD1"),
  "CD56dim NK"             = c("FCGR3A", "NKG7", "GNLY", "KLRD1", "KLRF1", "GZMB", "PRF1", "KIR2DL1", "KIR3DL1"))

canonical_markers <- list(
  # ==================== B细胞亚群 ====================
  "B cells naive"          = c("CD79A", "CD79B", "MS4A1", "CD19", "IGHD", "IGHM", "TCL1A", "FCER2", "CD24"),
  "B cells memory"         = c("CD79A", "CD79B", "MS4A1", "CD19", "CD27", "TNFRSF13B", "IGHG1", "AIM2", "CR2"),
  "Plasma cells"           = c("JCHAIN", "MZB1", "SDC1", "IGHG1", "IGHG2", "XBP1", "PRDM1", "IRF4", "SLAMF7"))

canonical_markers <- list( 
 # ==================== 单核细胞亚群 ====================
  "Classical Monocytes"    = c("CD14", "LYZ", "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD36", "CSF1R"),
  "Non-classical Monocytes"= c("FCGR3A", "MS4A7", "LST1", "LILRB2", "IFITM3", "AIF1", "CDKN1C", "MTSS1", "ITGAX"),
  "Intermediate Monocytes" = c("CD14", "FCGR3A", "LYZ", "S100A8", "VCAN", "CD36", "ITGAM", "CSF1R", "FCN1"),
 # ==================== 树突状细胞亚群 ====================
  "cDC1 (经典1型DC)"       = c("CLEC9A", "XCR1", "CADM1", "ITGAX", "BATF3", "IRF8", "HLA-DRA", "CD86", "THBD"),
  "cDC2 (经典2型DC)"       = c("CD1C", "FCER1A", "CLEC10A", "ITGAM", "CD1E", "FCER1G", "HLA-DQA1", "CD11c", "IRF4"),
  "pDC (浆细胞样DC)"       = c("LILRA4", "IRF7", "CLEC4C", "IL3RA", "ITM2C", "PLD4", "GZMB", "TCL1A", "SCT"),
  # =================== 血小板/巨核细胞 ====================
  "Platelets"              = c("PPBP", "PF4", "GP9", "ITGA2B", "TUBB1", "TREML1", "CAVIN2", "CMTM5", "ITGB3"),
   # ==================== 造血干/祖细胞 ====================
  #"HSPC"                   = c("CD34", "CRHBP", "SPINK2", "HOPX", "AVP", "CYTL1", "MLLT3", "HLF", "PROM1"),
    # ==================== 粒细胞（PBMC微量） ====================
  #"Neutrophils"            = c("CSF3R", "CEACAM8", "MPO", "ELANE", "CAMP", "PGLYRP1", "ITGAM", "FCGR3B", "CD10"),
    # ==================== 红细胞前体 ====================
  "Erythroid"              = c("HBB", "HBA1", "HBA2", "GYPA", "TFRC", "GATA1", "KLF1", "AHSP", "BPGM")
) 

canonical_markers <- list(
  "Total T Cells" = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC1", "CD28"),
  "T_Naive_Common" = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R", "S1PR1"),
  "T_Memory_Common" = c("IL7R", "CD27", "CD28", "S100A4", "ITGB1"),
  "T_Effector_Common" = c("GZMA", "IFNG", "PRF1", "NKG7", "TBX21"),
  "T_Exhausted_Common" = c("PDCD1", "LAG3", "TIGIT", "CTLA4", "TOX", "EOMES"),
  "Treg_Specific" = c("FOXP3", "IL2RA", "IKZF2", "BATF"),
  "CD4+ T Naive"     = c("CD4", "CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
  "CD4+ T Memory"    = c("CD4", "IL7R", "S100A4", "CD69", "GPR183"),
  "CD4+ T Effector"  = c("CD4", "GZMA", "IFNG", "TNF", "LAG3"),
  "CD4+ T Exhausted" = c("CD4", "PDCD1", "LAG3", "TIGIT"),
  "CD8+ T Naive"     = c("CD8A", "CD8B", "CCR7", "SELL", "LEF1", "TCF7"),
  "CD8+ T Memory"    = c("CD8A", "CD8B", "IL7R", "GZMK", "S100A4"),
  "CD8+ T Effector"  = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7"),
  "CD8+ T Exhausted" = c("CD8A", "CD8B", "PDCD1", "TOX", "EOMES"),
  "NK_Common"        = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "GZMB"),
  "CD56bright NK"    = c("NCAM1", "XCL1", "IL7R", "IFNG"),
  "CD56dim NK"       = c("FCGR3A", "KIR2DL1", "PRF1", "KLRF1")
)

library(Seurat)
library(ggplot2)

# 1. 定义你的原始列表
extended_markers2 <- list(
  "T_cells"       = c("CD3D", "CD3E", "TRAC"),
  "CD4_naive"     = c("CD4", "IL7R", "CCR7", "TCF7", "SELL"),
  "CD4_memory"    = c("IL7R", "S100A4", "ANXA1"),
  "Treg"          = c("FOXP3", "IL2RA", "CTLA4"),
  "CD8_naive"     = c("CD8A", "CCR7", "TCF7"),
  "CD8_cytotoxic" = c("CD8A", "GZMB", "PRF1", "GNLY", "NKG7"),
  "NK_cells"      = c("GNLY", "NKG7", "KLRD1", "NCAM1"),
  "B_naive"       = c("MS4A1", "CD19", "CD79A", "IGHD"),
  "B_memory"      = c("MS4A1", "CD27"),
  "Plasma"        = c("MZB1", "XBP1", "PRDM1", "SDC1"),
  "CD14_Mono"     = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN"),
  "CD16_Mono"     = c("FCGR3A", "MS4A7", "LILRB2"),
  "cDC"           = c("FCER1A", "CD1C", "CLEC10A"),
  "pDC"           = c("LILRA4", "IL3RA", "IRF7"),
  "Platelet"      = c("PPBP", "PF4", "GP9"),
  "HSPC"          = c("CD34", "SPINK2", "HLF")
)

# 2. 【核心修复】自动去重逻辑
# 我们遍历列表，记录已经出现过的基因，只保留第一次出现的基因
extended_markers2=canonical_markers
seen_genes <- c()
unique_markers <- lapply(extended_markers2, function(genes) {
  # 只保留对象中存在的基因
  existing_genes <- genes[genes %in% rownames(merged_obj)]
  # 排除掉之前已经出现过的基因
  new_genes <- setdiff(existing_genes, seen_genes)
  # 更新已出现的基因列表
  seen_genes <<- c(seen_genes, new_genes)
  return(new_genes)
})

# 移除掉因为去重后变为空的分组
unique_markers <- unique_markers[sapply(unique_markers, length) > 0]

# 3. 运行 DotPlot
pdf(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "_CanonicalMarkers_DotPlot_Tcells_.pdf")),
  		 width = 26, height = 10)
p3 <- DotPlot(
  merged_obj,
  features = unique_markers, 
  dot.scale = 8
) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  labs(
    title = "Enhanced Marker Panel (De-duplicated)",
    x = "Unique Canonical Markers",
    y = "Cluster Identity"
  )

# 4. 打印图片
print(p3)     
ggsave(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "_CanonicalMarkers_DotPlot.pdf")),
       plot = p3, width = 26, height = 10, dpi = 300)
# FeaturePlot：每种细胞类型选 1 个代表基因  
representative_genes <- sapply(existing_markers, function(x) x[1])  
representative_genes <- representative_genes[!is.na(representative_genes)]  

p_feature <- FeaturePlot(merged_obj,  
                             features = representative_genes,  
                             reduction = "umap_harmony",  
                             ncol = 4, pt.size = 0.1,  
                             order = TRUE) &  
    scale_color_viridis_c()  
  ggsave(file.path(projectPath, "Output", "CellAnnotation",  
                   paste0(date_tag, "_RepresentativeMarkers_FeaturePlot.pdf")),  
         plot = p_feature,  
         width = 20, height = 5 * ceiling(length(representative_genes) / 4),  
         dpi = 300)  
}  

canonical_markers <- c("TRA", "CD3D", "IL7R", "CCR7", "S100A4", "CD8A", "GZMB", 
                       "NKG7", "GNLY", "MS4A1", "CD19", "CD14", "LYZ", "FCGR3A", 
                       "PPBP", "CD34")

p2 <- DotPlot(merged_obj, features = canonical_markers) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "Verified Canonical Markers by Cell Type")

print(p2)
                   

            
################################################
# ==============================================
三步注释逻辑（先圈 T→再分 CD4/CD8→最后细分亚群）
# ==============================================
################################################
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
seurat_obj <- merged_obj 
Idents(seurat_obj) <- "seurat_clusters"
DefaultAssay(seurat_obj) <- "RNA"
# ==============================================
# 【第一步】区分：T细胞 vs 非T细胞
# ==============================================
# 1.1 定义第一步Marker
step1_markers <- list(
  "T cells" = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1"),
  "B cells" = c("CD79A", "MS4A1"),
  "Monocytes" = c("CD14", "LYZ"),
  "NK cells" = c("NKG7", "GNLY")
)

# 1.2 可视化验证：DotPlot看Marker表达
p1 <- DotPlot(seurat_obj, features = step1_markers, cols = "RdYlBu") + 
  RotatedAxis() + 
  ggtitle("Step 1: Cluster Markers") + 
  theme(plot.title = element_text(hjust = 0.5))

# 1.3 筛选出T细胞（根据DotPlot结果，选择CD3高表达的cluster）
# 假设你看了DotPlot后，确定cluster 0,1,2,3,5是T细胞（替换为你实际的cluster编号）
T_clusters <- c(0, 2, 5, 6) 
T_cells <- subset(seurat_obj, idents = T_clusters)

# 1.4 可视化：UMAP展示T细胞圈选结果
seurat_obj$Step1_Annotation <- ifelse(
  seurat_obj$seurat_clusters %in% T_clusters,
  "T cells",
  "Non-T cells"
)
p2 <- DimPlot(seurat_obj, group.by = "Step1_Annotation", reduction = 'umap_harmony', cols = c("#E63946", "#A8DADC"), label = TRUE) +
  ggtitle("Step 1: T cells vs Non-T cells")

# ==============================================
# 【第二步】在T细胞内区分：CD4+ T vs CD8+ T
# ==============================================
# 2.1 定义第二步Marker
step2_markers <- list(
  "CD4+ T cells" = c("CD4", "IL7R"),
  "CD8+ T cells" = c("CD8A", "CD8B")
)
# 2.2 可视化验证：DotPlot
p3 <- DotPlot(T_cells, features = step2_markers, cols = c("lightgrey", "#457B9D")) +
  RotatedAxis() + ggtitle("Step 2: CD4+ vs CD8+ T cells")
# 2.3 精准区分CD4/CD8（基于基因表达互斥）
gene_expr <- FetchData(T_cells, vars = c("CD4", "CD8A", "CD8B"))
T_cells$CD4_expr <- gene_expr$CD4
T_cells$CD8A_expr <- gene_expr$CD8A
T_cells$CD8B_expr <- gene_expr$CD8B
T_cells$Step2_Annotation <- case_when(
  # CD4+ T：CD4表达 > 0 且 CD8A/CD8B 都不表达
  T_cells$CD4_expr > 0 & T_cells$CD8A_expr == 0 & T_cells$CD8B_expr == 0 ~ "CD4+ T cells", 
  # CD8+ T：CD8A或CD8B表达 > 0 且 CD4 不表达
  (T_cells$CD8A_expr > 0 | T_cells$CD8B_expr > 0) & T_cells$CD4_expr == 0 ~ "CD8+ T cells",
  # 其他情况（双阴/双阳）
  TRUE ~ "Double negative/positive T cells"
)
Idents(T_cells) <- "Step2_Annotation"
table(T_cells$Step2_Annotation)
# 2.4 可视化：UMAP展示CD4/CD8分群
p4 <- DimPlot(T_cells, group.by = "Step2_Annotation", reduction = 'umap_harmony',cols = c("#E63946", "#457B9D", "#999999"), label = TRUE) +
  ggtitle("Step 2: CD4+ vs CD8+ T cells")

########
Idents(T_cells) <- "seurat_clusters"
step02_markers <- list(
  "Naive T cells" = c("CCR7", "SELL", "TCF7", "LEF1"),
  "Regulatory T cells" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Memory T cells" = c("CD27", "CD45RO", "IL7R", "CCR7"),
  "Effector T cells" = c("GZMA", "GZMB", "PRF1", "IFNG", "TNF", "CD4", "CD8A")
)
unique_markers02 <- unique(unlist(step02_markers))
p44 <- DotPlot(
  T_cells,
  features = unique_markers02, 
  dot.scale = 8
) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  labs(
    title = "Enhanced Marker Panel (De-duplicated)",
    x = "Unique Canonical Markers",
    y = "Cluster Identity"
  )

# ==============================================
# 【第三步】在CD4/CD8内细分：功能亚群
# ==============================================
# 3.1 定义第三步Marker
step3_markers <- list(
  "CD4+ Naive" = c("CCR7", "SELL", "LEF1"),
  "CD4+ Memory" = c("IL7R", "S100A4", "CD69"),
  "CD4+ Effector" = c("GZMA", "IFNG", "TNF"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4"),
  "CD8+ Naive" = c("CCR7", "SELL", "LEF1"),
  "CD8+ Memory" = c("IL7R", "GZMK", "S100A4"),
  "CD8+ Effector" = c("GZMB", "PRF1", "NKG7"),
  "CD8+ Exhausted" = c("PDCD1", "LAG3", "TIGIT")
)
unique_markers <- unique(unlist(step3_markers))
Idents(T_cells) <- "seurat_clusters"
# 3.2 可视化验证：DotPlot看所有亚群Marker
p5 <- DotPlot(
  T_cells,
  features = unique_markers, # 替换为去重后的向量
  cols = c("lightgrey", "#1D3557")
) +
  RotatedAxis() +
  ggtitle("Step 3: Fine Annotation Markers") +
  theme(plot.title = element_text(hjust = 0.5))
print(p5)
p55 <- DotPlot(
  T_cells,
  features = unique_markers, 
  dot.scale = 8
) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  labs(
    title = "Enhanced Marker Panel (De-duplicated)",
    x = "Unique Canonical Markers",
    y = "Cluster Identity"
  )


# 3.3 细分亚群（基于Marker表达，这里以CD4为例，CD8同理）
# 先把CD4和CD8分开
CD4_cells <- subset(T_cells, idents = "CD4+ T cells")
CD8_cells <- subset(T_cells, idents = "CD8+ T cells")
# 分别定义 CD4 和 CD8 的 Marker（无重复）
cd4_markers <- list(
  "CD4+ Naive" = c("CCR7", "SELL", "LEF1"),
  "CD4+ Memory" = c("IL7R", "S100A4", "CD69"),
  "CD4+ Effector" = c("GZMA", "IFNG", "TNF"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4")
)
cd8_markers <- list(
  "CD8+ Naive" = c("CCR7", "SELL", "LEF1"),
  "CD8+ Memory" = c("IL7R", "GZMK", "S100A4"),
  "CD8+ Effector" = c("GZMB", "PRF1", "NKG7"),
  "CD8+ Exhausted" = c("PDCD1", "LAG3", "TIGIT")
)
#画 CD4 的 DotPlot
Idents(CD4_cells) = 'seurat_clusters'
p5_cd4 <- DotPlot(
  T_cells,
  features = cd4_markers,
  cols = c("lightgrey", "#E63946")
) +
  RotatedAxis() +
  ggtitle("CD4+ T Cell Subset Markers") +
  theme(plot.title = element_text(hjust = 0.5))
# 画 CD8 的 DotPlot
Idents(CD8_cells) = 'seurat_clusters'
p5_cd8 <- DotPlot(
  T_cells,
  features = cd8_markers,
  cols = c("lightgrey", "#457B9D")
) +
  RotatedAxis() +
  ggtitle("CD8+ T Cell Subset Markers") +
  theme(plot.title = element_text(hjust = 0.5))
# 拼图展示（论文级）
library(patchwork)
#p5_cd4_cd8 = p5_cd4 + p5_cd8

# ==============================================
# 最终汇总与导出
# ==============================================
# 建议：最终的 Fine_Annotation 应该基于 Step2 + seurat_clusters 的综合判断
T_cells$Fine_Annotation <- paste0(T_cells$Step2_Annotation, "_Cluster_", T_cells$seurat_clusters)
# 拼图展示
#final_plot <- (p1 | p2) / (p4) / (p5_cd4 | p5_cd8) + 
#               plot_layout(heights = c(1, 0.8, 1.2))
#print(final_plot)# 查看最终注释统计
table(T_cells$Fine_Annotation)
# 3. 运行 DotPlot
pdf(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "_CanonicalMarkers_DotPlot_Tcells_DotplotDimplotDetailed.pdf")),
  		 width = 26, height = 10)
print(p1)
print(p2)
print(p3)
print(p4)
print(p44)
print(p5)
print(p55)
print(p5_cd4)
print(p5_cd8)
#print(p5_cd4_cd8)
dev.off()










# 3.2 特征表达图 (FeaturePlot) - 查看核心 Marker 在 UMAP 上的分布
DotPlot(
  monocyte_obj,
  features = markers_to_check,
  group.by = "RNA_snn_res.0.05",
  dot.scale = 10,          # 控制点的大小范围
  scale = TRUE,            # 标准化表达量，使颜色对比更明显
  cols = c("lightgrey", "blue"), # 颜色梯度：浅灰（低表达）→ 蓝色（高表达）
  col.min = -1,            # 颜色梯度最小值
  col.max = 1              # 颜色梯度最大值
) +
  scale_y_discrete(labels = c("Classical Monocyte" = "0", "Non-classical Monocyte" = "1")) + # 对应你的y轴0/1
  theme_bw() +
  labs(
    x = "Features",
    y = "Identity",
    fill = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12)
  )







load('/mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia/scRNA_PBMC_DPN/Data/08_260417_6Samples_Harmony_Clustered_Annotated.Rdata')
merged_obj
# 查看对象
merged_obj
table(merged_obj$seurat_clusters)



# ================================================================  
# 11. 手动注释模板（根据 SingleR + Marker 结果调整）  
# ================================================================  

message("\n========================================")  
message("Step 11: Manual Annotation Template")  
message("========================================")  

# ★★★ 查看 SingleR 注释 + marker 后，在此处手动修改 ★★★  
# 示例模板（运行后根据实际结果调整）：  
#  
merged_obj <- RenameIdents(merged_obj,
  "0"  = "CD4+ naive T",
  "1"  = "NK cells",
  "2"  = "CD8+ effector T",
  "3"  = "Classical monocytes",
  "4"  = "NK cells",
  "5"  = "CD8+ naive T",
  "6"  = "CD4+ memory T",
  "7"  = "Plasma cells",
  "8"  = "Non-classical monocytes",
  "9"  = "B cells",
  "10" = "DC cells",
  "11" = "Megakaryocytes"
)

# 应用注释
merged_obj$cell_type_manual <- manual_labels[as.character(merged_obj$seurat_clusters)]

# 查看注释结果
table(merged_obj$cell_type_manual)

# ================================================================
# 绘制 UMAP 图（手动注释）
# ================================================================
p_manual <- DimPlot(
  merged_obj,
  reduction = "umap_harmony",
  group.by = "cell_type_manual",
  label = TRUE,        # 补上了！
  label.size = 4,
  repel = TRUE,
  pt.size = 0.1
) +
  ggtitle("Cell Types (Manual Annotation)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(p_manual)
                              
p_manual <- DimPlot(merged_obj,
                    reduction = "umap_harmony",
                    group.by = "cell_type_manual",
                    label = TRUE,
                    repel = TRUE,
                    label.size = 4,
                    pt.size = 0.1)

print(p_manual)
ggsave(file.path(projectPath, "Output", "CellAnnotation",
                 paste0(date_tag, "Harmony_Clustered_Annotated_CellTypeManual.pdf")),
       plot = p_manual, width = 10, height = 10, dpi = 300)

save_file_annotated <- file.path(projectPath, "Data",
                                  make_filename(
                                    length(unique(merged_obj$orig.ident)),
                                    "Harmony_Clustered_Annotated_CellTypeManual",
                                    prefix = "08_"))
save(merged_obj, file = save_file_annotated)
message("✅ 已保存最终注释对象：", save_file_annotated)
                               
                               
# 流程完成提示
message("\n========================================")
message("★ Full Single-cell RNA-seq Pipeline Finished!")
message("========================================")
message(sprintf("Date: %s", Sys.Date()))
message(sprintf("Samples processed: %d", length(seurat_list)))
message(sprintf("Final cells: %d", ncol(merged_obj)))
message(sprintf("Final genes: %d", nrow(merged_obj)))
message(sprintf("Cell types identified: %d", length(unique(Idents(merged_obj)))))
message("\nFinal saved files:")
message(sprintf("01a: %s (before QC)", basename(save_file_raw)))
message(sprintf("01b: %s (after QC)", basename(save_file_filtered)))
message(sprintf("02: %s (after DoubletFinder)", basename(save_file_df)))
message(sprintf("03: %s (with CellCycle)", basename(save_file_cc)))
message(sprintf("04: %s (Harmony integrated)", basename(save_file_harmony)))
message(sprintf("05: %s (final clustering)", basename(save_file_cluster)))
message(sprintf("06: %s (final annotated)", basename(save_file_final)))
message("========================================")
