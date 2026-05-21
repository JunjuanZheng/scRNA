# ================================================================
# DPN vs Control 单细胞转录组比较分析 完整 Pipeline
# Control: HEA16, HEA17, HEA22
# DPN:     DPN-01, DPN-02, LZF
# ================================================================
# 在脚本第一行加上
# -*- coding: utf-8 -*-
options(encoding = "UTF-8")
Sys.setlocale("LC_ALL", "en_US.UTF-8")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(tidyr)
library(tibble)
library(future)
library(scales)

plan("multicore", workers = 6)
options(future.globals.maxSize = 80 * 1024^3)
set.seed(42)

# ★ 日期标签 ★
date_tag <- format(Sys.Date(), "%y%m%d")

# 动态文件名生成函数
make_filename <- function(n_samples, suffix, prefix = "", ext = "Rdata") {
  fname <- paste0(prefix, date_tag, "_", n_samples, "Samples_", suffix, ".", ext)
  return(fname)
}

# 并行设置
plan("multicore", workers = 6)
options(future.globals.maxSize = 80 * 1024^3)
set.seed(42)

projectPath='/mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia/scRNA_PBMC_DPN'
# 路径
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "01_DataPrep"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "02_CellComp"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "03_DEG"),           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "05_Pathway"),       recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "06_CellChat"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "07_Trajectory"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "08_TF"),            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "09_WGCNA"),         recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "10_CellScore"),     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(projectPath, "Output", "DPN_Analysis", "11_Summary"),       recursive = TRUE, showWarnings = FALSE)


# ================================================================
# 步骤 1：数据准备与分组
# ================================================================
message("\n===== Step 1: Data Preparation =====")

load(file.path(projectPath, "Data",
               "08_260513_6Samples_Harmony_Clustered_Annotated_CellTypeManual.Rdata"))

# ---- 1.1 添加分组标签 ----
group_map <- c(
  "HEA16"  = "Control",
  "HEA17"  = "Control",
  "HEA22"  = "Control",
  "DPN-01" = "DPN",
  "DPN-02" = "DPN",
  "LZF"    = "DPN"
)
merged_obj@meta.data$group <- group_map[merged_obj$orig.ident]
merged_obj$group <- factor(merged_obj$group, levels = c("Control", "DPN"))

# ---- 1.2 设置主要 Ident 为手动注释的细胞类型 ----
Idents(merged_obj) <- 'cell_type_manual'
table(merged_obj$cell_type_manual)

# 定义配色
group_colors     <- c("Control" = "#4878D0", "DPN" = "#EE854A")
cell_type_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(levels(merged_obj$cell_type_manual))),
  levels(merged_obj$cell_type_manual)
)

# ---- 1.3 概览 UMAP（按组）----
p_umap_group <- DimPlot(merged_obj, reduction = "umap_harmony",
                         group.by = "group", pt.size = 0.1,
                         cols = group_colors) +
  ggtitle("UMAP: Control vs DPN") +
  theme(plot.title = element_text(hjust = 0.5))

p_umap_split <- DimPlot(merged_obj, reduction = "umap_harmony",
                          split.by = "group", pt.size = 0.1,
                          cols = cell_type_colors,
                          label = TRUE, label.size = 3) +
  ggtitle("UMAP Split by Group")

p_umap_celltype <- DimPlot(merged_obj, reduction = "umap_harmony",
                             group.by = "cell_type_manual",
                             label = TRUE, label.size = 3.5,
                             pt.size = 0.1, repel = TRUE,
                             cols = cell_type_colors) +
  ggtitle("Cell Type Annotation") +
  NoLegend()

p_step1 <- (p_umap_celltype | p_umap_group) / p_umap_split
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "01_DataPrep",
                 paste0(date_tag, "_Overview_UMAP.png")),
       plot = p_step1, width = 22, height = 16, dpi = 300)

# 样本信息汇总
sample_info <- merged_obj@meta.data %>%
  group_by(orig.ident, group) %>%
  summarise(N_cells = n(), .groups = "drop")
message("  Sample info:")
print(sample_info)
write.csv(sample_info,
          file.path(projectPath, "Output", "DPN_Analysis", "01_DataPrep",
                    paste0(date_tag, "_SampleInfo.csv")),
          row.names = FALSE)


# ================================================================
# 步骤 2：细胞组成差异分析
# ================================================================
message("\n===== Step 2: Cell Composition Analysis =====")

# ---- 2.1 计算各样本各细胞类型比例 ----
comp_df <- merged_obj@meta.data %>%
  group_by(orig.ident, group, cell_type_manual) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()

write.csv(comp_df,
          file.path(projectPath, "Output", "DPN_Analysis", "02_CellComp",
                    paste0(date_tag, "_CellComposition.csv")),
          row.names = FALSE)

# ---- 2.2 堆叠柱状图 ----
p_stack <- ggplot(comp_df, aes(x = orig.ident, y = prop, fill = cell_type_manual)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = cell_type_colors) +
  facet_grid(~ group, scales = "free_x", space = "free") +
  theme_classic() +
  labs(title = "Cell Type Proportion by Sample",
       x = NULL, y = "Proportion (%)", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey90"))
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "02_CellComp",
                 paste0(date_tag, "_Composition_StackBar.png")),
       plot = p_stack, width = 10, height = 7, dpi = 300)

# ---- 2.3 各细胞类型两组比例比较（箱线图 + 统计）----
comp_summary <- comp_df %>%
  group_by(group, cell_type_manual) %>%
  summarise(mean_prop = mean(prop), sd_prop = sd(prop), .groups = "drop")

p_boxcomp <- ggplot(comp_df, aes(x = cell_type_manual, y = prop, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, position = position_dodge(0.75)) +
  geom_jitter(aes(color = group),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  stat_compare_means(aes(group = group), label = "p.signif",
                     method = "wilcox.test", hide.ns = FALSE) +
  theme_classic() +
  labs(title = "Cell Proportion: Control vs DPN",
       x = NULL, y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "02_CellComp",
                 paste0(date_tag, "_Composition_Boxplot.png")),
       plot = p_boxcomp, width = 14, height = 6, dpi = 300)

# ---- 2.4 统计检验（每种细胞类型 Wilcoxon）----
comp_test <- comp_df %>%
  group_by(cell_type_manual) %>%
  summarise(
    mean_Control = mean(prop[group == "Control"]),
    mean_DPN     = mean(prop[group == "DPN"]),
    pval = tryCatch(
      wilcox.test(prop[group == "Control"], prop[group == "DPN"])$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    fold_change = mean_DPN / (mean_Control + 1e-6),
    log2FC      = log2(fold_change),
    padj        = p.adjust(pval, method = "BH"),
    significant = ifelse(!is.na(padj) & padj < 0.05, "Yes", "No")
  )

write.csv(comp_test,
          file.path(projectPath, "Output", "DPN_Analysis", "02_CellComp",
                    paste0(date_tag, "_Composition_Statistics.csv")),
          row.names = FALSE)
message("  Cell composition test results:")
print(comp_test)

# ================================================================
# 全局工具函数：过滤 ENSG 基因（在所有步骤之前定义）
# ================================================================

# ★ 过滤 ENSG 开头的基因的通用函数
filter_ensg <- function(genes) {
  genes[!grepl("^ENSG", genes)]
}


# ================================================================
# 步骤 3：分细胞类型差异表达分析（Pseudo-bulk）
# ================================================================
message("\n===== Step 3: Cell-type Specific DEG Analysis =====")

library(DESeq2)

cell_types    <- levels(merged_obj$cell_type_manual)
all_deg_results <- list()

for (ct in cell_types) {

  message(sprintf("  Processing: %s", ct))

  cells_ct  <- subset(merged_obj, cell_type_manual == ct)
  n_control <- sum(cells_ct$group == "Control")
  n_dpn     <- sum(cells_ct$group == "DPN")

  if (n_control < 30 | n_dpn < 30) {
    message(sprintf("    Skipped %s: too few cells (Control=%d, DPN=%d)",
                    ct, n_control, n_dpn))
    next
  }

  # ---- 3.1 Pseudo-bulk 聚合 ----
  count_mat  <- GetAssayData(cells_ct, slot = "counts", assay = "RNA")
  sample_ids <- cells_ct$orig.ident
  samples_uniq <- unique(sample_ids)

  pb_counts <- sapply(samples_uniq, function(s) {
    idx <- which(sample_ids == s)
    if (length(idx) == 1) as.numeric(count_mat[, idx])
    else Matrix::rowSums(count_mat[, idx])
  })
  rownames(pb_counts) <- rownames(count_mat)
  pb_counts <- as.matrix(pb_counts)

  # 强制整数
  pb_counts <- round(pb_counts)
  storage.mode(pb_counts) <- "integer"
  pb_counts[!is.finite(pb_counts)] <- 0L

  # ★★★ 过滤 ENSG 开头的基因行 ★★★
  ensg_idx  <- grepl("^ENSG", rownames(pb_counts))
  pb_counts <- pb_counts[!ensg_idx, , drop = FALSE]
  message(sprintf("    Removed %d ENSG genes, %d genes remaining",
                  sum(ensg_idx), nrow(pb_counts)))

  # metadata
  pb_meta <- data.frame(
    sample    = samples_uniq,
    group     = group_map[samples_uniq],
    row.names = samples_uniq
  )
  pb_meta$group <- factor(pb_meta$group, levels = c("Control", "DPN"))
  pb_counts     <- pb_counts[, rownames(pb_meta), drop = FALSE]

  # 过滤低表达基因
  keep_genes <- rowSums(pb_counts >= 10) >= 2
  pb_counts  <- pb_counts[keep_genes, , drop = FALSE]

  message(sprintf("    Cells: Control=%d, DPN=%d | Genes after filter: %d",
                  n_control, n_dpn, nrow(pb_counts)))

  # ---- 3.2 DESeq2 ----
  dds <- tryCatch({
    dds_tmp <- DESeqDataSetFromMatrix(
      countData = pb_counts,
      colData   = pb_meta,
      design    = ~ group
    )
    DESeq(dds_tmp, quiet = TRUE)
  }, error = function(e) {
    message(sprintf("    DESeq2 error for %s: %s", ct, e$message))
    return(NULL)
  })
  if (is.null(dds)) next

  res <- results(dds,
                  contrast = c("group", "DPN", "Control"),
                  alpha    = 0.05)

  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    # ★ 双重保险：结果层再次过滤 ENSG ★
    filter(!grepl("^ENSG", gene)) %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(
      cell_type = ct,
      direction = case_when(
        padj < 0.05 & log2FoldChange >  0.5 ~ "Up in DPN",
        padj < 0.05 & log2FoldChange < -0.5 ~ "Down in DPN",
        TRUE ~ "NS"
      )
    )

  all_deg_results[[ct]] <- res_df

  n_up   <- sum(res_df$direction == "Up in DPN")
  n_down <- sum(res_df$direction == "Down in DPN")
  message(sprintf("    DEGs found: Up=%d, Down=%d", n_up, n_down))

  write.csv(res_df,
            file.path(projectPath, "Output", "DPN_Analysis", "03_DEG",
                      paste0(date_tag, "_DEG_", gsub("[/ +]", "_", ct), ".csv")),
            row.names = FALSE)

  # ---- 3.3 Volcano Plot ----
  # 标注基因：优先选有 symbol 的显著基因（ENSG 已被过滤）
  sig_genes <- res_df %>%
    filter(direction != "NS") %>%
    slice_max(abs(log2FoldChange), n = 20)

  p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                               color = direction)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Up in DPN"   = "#E64B35",
                                   "Down in DPN" = "#4DBBD5",
                                   "NS"          = "grey70")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = sig_genes,
                              aes(label = gene), size = 2.5,
                              max.overlaps = 20, show.legend = FALSE) +
    theme_classic() +
    labs(title    = paste0("DEG Volcano: ", ct),
         subtitle = sprintf("Up=%d  Down=%d", n_up, n_down),
         x = "log2 Fold Change (DPN / Control)",
         y = "-log10(padj)")

  ggsave(file.path(projectPath, "Output", "DPN_Analysis", "03_DEG",
                   paste0(date_tag, "_Volcano_", gsub("[/ +]", "_", ct), ".png")),
         plot = p_vol, width = 8, height = 6, dpi = 300)
}

# ================================================================
# 步骤 3：分细胞类型差异表达分析 (Single-cell Level - FindMarkers)
# ================================================================
message("\n===== Step 3: Cell-type Specific DEG Analysis (FindMarkers) =====")
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- 3.0 准备工作：构建联合分组列 ----
# 格式如：DPN_Classical Monocytes, Control_Classical Monocytes
merged_obj$Group_CellTypeMannual <- paste0(merged_obj$group, "_", merged_obj$cell_type_manual)
table(merged_obj$Group_CellTypeMannual)

# 设置当前的 Ident 为新构建的联合分组列
Idents(merged_obj) <- "Group_CellTypeMannual"

cell_types <- levels(as.factor(merged_obj$cell_type_manual))
all_deg_results <- list()

# 确保输出目录存在
deg_out_dir <- file.path(projectPath, "Output", "DPN_Analysis", "03_DEG")
if(!dir.exists(deg_out_dir)) dir.create(deg_out_dir, recursive = TRUE)

groups = levels(as.factor(merged_obj$group))
groups
caseGroup = groups[2]
controlGroup = groups[1]

# 初始化一个用于存放汇总数量的列表
deg_count_summary <- list()

for (ct in cell_types) {
  ident_case <- paste0(caseGroup, "_", ct)
  ident_ctrl <- paste0(controlGroup, "_", ct)
  
  message(sprintf(">>> Analyzing Cell Type: %s", ct))
  
  # 3.1. 细胞数检查
  n_dpn <- sum(merged_obj$Group_CellTypeMannual == ident_case)
  n_ctrl <- sum(merged_obj$Group_CellTypeMannual == ident_ctrl)
  
  if (n_dpn < 10 | n_ctrl < 10) {
    message(sprintf("    [Skipped] %s: Cells too few (DPN=%d, Control=%d)", ct, n_dpn, n_ctrl))
    next
  }
  # 3.2. 差异分析
  res <- FindMarkers(merged_obj, 
                     ident.1 = ident_case, 
                     ident.2 = ident_ctrl,
                     logfc.threshold = 0,
                     min.pct = 0.1, 
                     only.pos = FALSE)
  # 3.3. 结果过滤与清洗
  res_df <- res %>%
    rownames_to_column("gene") %>%
    filter(!grepl("^ENSG", gene)) %>% # 过滤犬类非命名基因
    mutate(
      cell_type = ct,
      direction = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up in DPN",
        p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down in DPN",
        TRUE ~ "NS"
      )
    )
  
  # 3.4. 统计数量供后续画柱状图
  n_up <- sum(res_df$direction == "Up in DPN")
  n_down <- sum(res_df$direction == "Down in DPN")
  deg_count_summary[[ct]] <- data.frame(CellType = ct, Up = n_up, Down = n_down)
  
  # 3.5. 绘制火山图
  file_name_ct <- gsub("[/ +]", "_", ct)
  top_genes <- res_df %>% filter(direction != "NS") %>% slice_max(abs(avg_log2FC), n = 20)
  
  p_vol <- ggplot(res_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c("Up in DPN" = "#E64B35", "Down in DPN" = "#4DBBD5", "NS" = "grey80")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 3, fontface = "italic") +
    theme_minimal() +
    labs(title = paste0("Volcano: ", ct), subtitle = sprintf("Up=%d, Down=%d", n_up, n_down))

  # 强制保存并打印状态
  vol_path <- file.path(deg_out_dir, paste0(date_tag, "_Volcano_", file_name_ct, ".png"))
  ggsave(vol_path, plot = p_vol, width = 7, height = 6, dpi = 300)
  
  if(file.exists(vol_path)) message(sprintf("    [Success] Volcano saved to %s", vol_path))
  
  write.csv(res_df, file.path(deg_out_dir, paste0(date_tag, "_DEG_", file_name_ct, ".csv")), row.names = FALSE)
}

# 3.6. 生成差异基因数量汇总柱状图 (Bar Plot)
message("\n===== Generating High-Quality Summary Bar Plot =====")
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats) # 用于处理因子排序
# 整理数据并计算总数用于排序
summary_df <- do.call(rbind, deg_count_summary) %>%
  mutate(Total = Up + Down) %>%
  # 按照总数排序，如果总数一样，则按 Up 排序
  arrange(desc(Total), desc(Up))
# 转换为长格式并处理绘图逻辑
plot_data <- summary_df %>%
  pivot_longer(cols = c("Up", "Down"), names_to = "Direction", values_to = "Count") %>%
  mutate(
    # 为了镜像效果，将 Down 设置为负值
    PlotCount = ifelse(Direction == "Down", -Count, Count),
    # 强制固定 CellType 的顺序（由 Total 决定）
    CellType = fct_reorder(CellType, Total)
  )
#  绘图
p_bar_polished <- ggplot(plot_data, aes(x = CellType, y = PlotCount, fill = Direction)) +
  # 绘制柱子，设置边框色为白色可以让柱子更有质感
  geom_bar(stat = "identity", width = 0.75, color = "white", size = 0.2) +
  # 添加基准线
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  # 核心：添加数值标签
  # 使用 ifelse 处理标签位置，避免 0 被标出来
  geom_text(aes(label = ifelse(Count > 0, Count, ""), 
                y = PlotCount, 
                hjust = ifelse(Direction == "Up", -0.2, 1.2)), 
            size = 3.5, fontface = "bold") +
  # 配色方案 (选用更具学术感的红蓝)
  scale_fill_manual(values = c("Up" = "#DC0000FF", "Down" = "#3C5488FF"),
                    labels = c("Down" = "Down-regulated", "Up" = "Up-regulated")) +
  # 坐标轴处理：去除负号，并留出足够的空间给标签
  scale_y_continuous(labels = abs, 
                     expand = expansion(mult = c(0.15, 0.15))) + 
  coord_flip() + # 横向排列
  theme_classic() + # 经典学术主题
  labs(title = "Differential Gene Statistics",
       subtitle = "DPN vs Control Group",
       x = NULL, 
       y = "Number of Differentially Expressed Genes",
       fill = "Regulation") +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.y = element_text(face = "bold.italic"), # 细胞类型通常用斜体
    axis.line.y = element_blank(), # 去除纵轴线，通过中间的 hline 区分
    axis.ticks.y = element_blank(),
    legend.position = "top", # 图例放上方更符合你给出的参考图
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30")
  )
final_bar_path <- file.path(deg_out_dir, paste0(date_tag, "_DEG_Summary_BarPlot.png"))
ggsave(final_bar_path, plot = p_bar_polished, width = 9, height = 7, dpi = 300)
message(sprintf(">>> Professional BarPlot saved to: %s", final_bar_path))
# 汇总所有结果到一个大表
#head(all_deg_results)
#all_deg_df <- bind_rows(all_deg_results)
#write.csv(all_deg_df, file.path(deg_out_dir, paste0(date_tag, "_All_CellType_DEGs_Combined.csv")), row.names = FALSE)

# 设置保存路径
save_path <- file.path(projectPath, "Output", "DPN_Analysis", "03_DEG", 
                       paste0(date_tag, "_all_deg_results_list.rds"))

# 存储变量
saveRDS(all_deg_results, file = save_path)

message(sprintf("✅ 变量已成功存储至: %s", save_path))

all_deg_results=readRDS('/mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia/scRNA_PBMC_DPN/Output/DPN_Analysis/03_DEG/260513_all_deg_results_list.rds')
#all_deg_results

# ================================================================
# Step 4: Key Gene Identification
# ================================================================
message("\n===== Step 4: Key Gene Identification =====")

# Check all_deg_results
if (!exists("all_deg_results")) {
  stop("all_deg_results not found! Please run Step 3 first.")
}
if (length(all_deg_results) == 0) {
  stop("all_deg_results is empty! No cell type passed DESeq2 in Step 3.")
}

message(sprintf("  Cell types successfully analyzed in Step 3: %d",
                length(all_deg_results)))
message(sprintf("  Cell type list: %s",
                paste(names(all_deg_results), collapse = ", ")))

# ---- Generate all_deg_df (filter out ENSG genes) ----
all_deg_df <- bind_rows(all_deg_results) %>%
  filter(!grepl("^ENSG", gene))

message(sprintf("  Total DEG entries after merging: %d", nrow(all_deg_df)))
message(sprintf("  Significant DEGs (direction != NS): %d",
                sum(all_deg_df$direction != "NS")))

# Save combined DEG table
write.csv(all_deg_df,
          file.path(projectPath, "Output", "DPN_Analysis", "03_DEG",
                    paste0(date_tag, "_AllDEG_Combined.csv")),
          row.names = FALSE)

# DEG count summary
deg_summary <- all_deg_df %>%
  filter(direction != "NS") %>%
  group_by(cell_type, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)
message("\n  DEG Summary by Cell Type:")
print(deg_summary)

# DEG bar plot
deg_plot_df <- all_deg_df %>%
  filter(direction != "NS") %>%
  group_by(cell_type, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(n_signed = ifelse(direction == "Down in DPN", -n, n))

p_deg_bar <- ggplot(deg_plot_df,
                     aes(x = reorder(cell_type, abs(n_signed)),
                         y = n_signed, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up in DPN"   = "#E64B35",
                                "Down in DPN" = "#4DBBD5")) +
  theme_classic() +
  labs(title = "DEG Count by Cell Type",
       x = NULL, y = "Number of DEGs (up / down)")
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "03_DEG",
                 paste0(date_tag, "_DEG_BarPlot.png")),
       plot = p_deg_bar, width = 10, height = 7, dpi = 300)

# ---- 4.1 Multi-celltype shared DEGs ----
# 过滤显著差异基因（与原代码一致，确保all_deg_df列名与sig_deg相同）
sig_deg <- all_deg_df %>% filter(direction != "NS")

# ==============================================
# 跨细胞类型上调基因统计（DPN中高表达）
# ==============================================
up_genes_per_ct <- sig_deg %>%
  filter(
    direction == "Up in DPN",
    !grepl("^ENSG", gene)  # 过滤ENSG编号基因
  ) %>%
  group_by(gene) %>%
  summarise(
    n_celltypes = n_distinct(cell_type),  # 该基因在多少种细胞类型中上调
    cell_types  = paste(cell_type, collapse = "; "),  # 列出所有细胞类型
    mean_log2FC = mean(avg_log2FC),  # 平均log2倍数变化
    min_padj    = min(p_val_adj),  # 最小校正p值
    .groups = "drop"
  ) %>%
  filter(n_celltypes >= 2) %>%  # 只保留在至少2种细胞类型中上调的基因
  arrange(desc(n_celltypes), desc(mean_log2FC))  # 按细胞类型数和倍数变化排序

# ==============================================
# 跨细胞类型下调基因统计（DPN中低表达）
# ==============================================
dn_genes_per_ct <- sig_deg %>%
  filter(
    direction == "Down in DPN",
    !grepl("^ENSG", gene)
  ) %>%
  group_by(gene) %>%
  summarise(
    n_celltypes = n_distinct(cell_type),
    cell_types  = paste(cell_type, collapse = "; "),
    mean_log2FC = mean(avg_log2FC),
    min_padj    = min(p_val_adj),
    .groups = "drop"
  ) %>%
  filter(n_celltypes >= 2) %>%
  arrange(desc(n_celltypes), desc(abs(mean_log2FC)))  # 下调基因按绝对值排序

message(sprintf("  Multi-celltype UP   genes (>=2 types): %d",
                nrow(up_genes_per_ct)))
message(sprintf("  Multi-celltype DOWN genes (>=2 types): %d",
                nrow(dn_genes_per_ct)))

write.csv(up_genes_per_ct,
          file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes",
                    paste0(date_tag, "_MultiCelltype_UpGenes.csv")),
          row.names = FALSE)
write.csv(dn_genes_per_ct,
          file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes",
                    paste0(date_tag, "_MultiCelltype_DownGenes.csv")),
          row.names = FALSE)

# ---- 4.2 Heatmap ----
top_up        <- head(up_genes_per_ct$gene, 30)
top_down      <- head(dn_genes_per_ct$gene, 20)
top_genes_all <- c(top_up, top_down)

if (length(top_genes_all) > 0) {
  # 构建差异表达矩阵（仅保留显著差异基因）
  deg_matrix <- all_deg_df %>%
    filter(
      gene %in% top_genes_all,
      direction != "NS"  # 双重保险：只保留显著差异基因
    ) %>%
    select(gene, cell_type, avg_log2FC) %>%  # ★ 核心修改：log2FoldChange → avg_log2FC
    pivot_wider(
      names_from  = cell_type,
      values_from = avg_log2FC,  # ★ 核心修改
      values_fill = 0  # 无差异的细胞类型填充为0
    ) %>%
    column_to_rownames("gene")

  # 数据验证：确保矩阵非空
  if (nrow(deg_matrix) == 0 || ncol(deg_matrix) == 0) {
    message("  Warning: No valid DEGs for heatmap, skipping.")
  } else {
    # ★ 新增：数值截断（将超出±3的数值限制在范围内，避免极端值拉偏颜色）
    deg_matrix[deg_matrix > 3] <- 3
    deg_matrix[deg_matrix < -3] <- -3

    # 绘制热图
    pheatmap(
      deg_matrix,
      color         = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
      breaks        = seq(-3, 3, length.out = 101),
      cluster_rows  = TRUE,
      cluster_cols  = TRUE,
      show_rownames = TRUE,
      fontsize_row  = 8,
      fontsize_col  = 10,
      border_color  = NA,  # 去掉单元格边框，更美观
      main          = "Top DEGs Across Cell Types (DPN vs Control)",
      filename      = file.path(projectPath, "Output", "DPN_Analysis",
                                 "04_KeyGenes",
                                 paste0(date_tag, "_KeyGenes_Heatmap.png")),
      width         = 14,
      height        = 12
    )
    message("  Heatmap saved successfully.")
  }
} else {
  message("  Warning: No cross-celltype DEGs found, skipping Heatmap.")
}
  
# ---- 4.3 UpSet Plot ----
# ==============================================
# 绘制UpSet图（优化版）
# ==============================================
if (length(upset_list) >= 3) {
  
  # 关键修正 1：清洗列表名称，防止特殊字符干扰 UpSetR 的列名转换
  # UpSetR 会把名称转为数据框列名，特殊字符会导致匹配失败
  names(upset_list) <- gsub("[/ +\\-]", "_", names(upset_list))
  
  # 调试：检查转换后的矩阵是否为空
  test_df <- fromList(upset_list)
  if (nrow(test_df) == 0 || ncol(test_df) == 0) {
    stop("错误：fromList 转换后的数据框为空，请检查输入基因集是否存在交集。")
  }

  png(
    file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes",
              paste0(date_tag, "_KeyGenes_UpSet.png")),
    width = 1800, 
    height = 1200, 
    res = 150
  )
  
  # 关键修正 2：使用 print() 包裹 upset() 函数
  # 确保绘图对象被正确发送到 png 设备
  p <- upset(
    test_df, # 使用预处理好的数据框
    nsets = length(upset_list),
    nintersects = 40,
    order.by = "freq",
    decreasing = TRUE,
    show.numbers = "yes",
    number.angles = 0,
    point.size = 3,
    line.size = 1,
    text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1),
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Genes per Cell Type",
    main.bar.color = "#E64B35",
    sets.bar.color = "#4DBBD5",
    matrix.color = "#E64B35",
    shade.color = "lightgray"
  )
  
  print(p) # 显式打印
  
  dev.off()
  message("✅ UpSet plot saved successfully.")
} else {
  # ... 原有的警告信息 ...
}
##########
#提取所有细胞类型共同上调的“核心基因”
##########
# 假设你的 upset_list 包含 5 个细胞类型
# 提取所有细胞类型的交集
core_up_genes <- Reduce(intersect, upset_list)
message(sprintf("检测到 %d 个在所有细胞类型中共同上调的核心基因", length(core_up_genes)))
print(core_up_genes)
# 保存到本地
write.table(core_up_genes, 
            file = file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes", "Core_Up_Regulated_Genes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
#############
#生成“基因-细胞类型”对应矩阵（最推荐，最全面）
#############
library(tidyverse)
# 将列表转换为长格式数据框
gene_intersection_df <- stack(upset_list) %>%
  rename(gene = values, cell_type = ind) %>%
  group_by(gene) %>%
  # 将该基因出现的细胞类型合并到一列，用分号隔开
  summarise(
    occurrence_count = n(),
    present_in = paste(sort(cell_type), collapse = "; ")
  ) %>%
  arrange(desc(occurrence_count))

# 保存这个全表，这对应了 UpSet 图里所有的柱子信息
write.csv(gene_intersection_df, 
          file = file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes", "UpSet_Gene_Intersections_Detail.csv"),
          row.names = FALSE)

# 看看那些在 3 个及以上亚群中都上调的基因
high_sharing_genes <- gene_intersection_df %>% filter(occurrence_count >= 3)

#############
#提取特定亚群之间特有的交集
#############
# 定义你感兴趣的组合
target_cts <- c("Classical_Monocytes", "Non_classical_Monocytes")
# 找出这些亚群的并集
all_cts <- names(upset_list)
other_cts <- setdiff(all_cts, target_cts)

# 逻辑：在目标亚群的交集里，且不在其他亚群的并集里
specific_shared_genes <- intersect(upset_list[[target_cts[1]]], upset_list[[target_cts[2]]]) %>%
  setdiff(unlist(upset_list[other_cts]))

message(sprintf("仅在单核细胞亚群间共享的基因数: %d", length(specific_shared_genes)))














# ---- 4.4 VlnPlot ----
key_genes_show <- head(top_up[!grepl("^ENSG", top_up)], 9)
key_genes_show <- key_genes_show[key_genes_show %in% rownames(merged_obj)]

if (length(key_genes_show) > 0) {
  message(sprintf("  VlnPlot genes: %s",
                  paste(key_genes_show, collapse = ", ")))
  p_vln <- VlnPlot(merged_obj,
                    features   = key_genes_show,
                    group.by   = "cell_type_manual",
                    split.by   = "group",
                    split.plot = TRUE,
                    pt.size    = 0,
                    cols       = group_colors,
                    ncol       = 3) &
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
  ggsave(file.path(projectPath, "Output", "DPN_Analysis", "04_KeyGenes",
                   paste0(date_tag, "_KeyGenes_VlnPlot.png")),
         plot = p_vln, width = 18, height = 18, dpi = 300)
  message("  VlnPlot saved.")
} else {
  message("  Warning: No key genes found in Seurat object, skipping VlnPlot.")
}

message("\n  Step 4 complete!")





all_deg_results=readRDS('/mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia/scRNA_PBMC_DPN/Output/DPN_Analysis/03_DEG/260513_all_deg_results_list.rds')
#all_deg_results
# ================================================================
# 步骤 5：功能富集分析（GO / KEGG / GSEA 均过滤 ENSG）
# ================================================================
# ================================================================
# Step 5: GO / KEGG Enrichment Analysis by Cell Type
# ================================================================
message("\n===== Step 5: Pathway Enrichment Analysis =====")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)

# 输出目录
pathway_dir <- file.path(projectPath, "Output", "DPN_Analysis", "05_Pathway")
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

# 检查对象
all_deg_df <- bind_rows(all_deg_results) %>%
  filter(!grepl("^ENSG", gene))

if (!exists("all_deg_df")) {
  stop("all_deg_df not found! Please run DEG analysis first.")
}

sig_deg <- all_deg_df %>%
  filter(direction != "NS",
         !grepl("^ENSG", gene))

cell_types_use <- unique(sig_deg$cell_type)

message(sprintf("  Cell types for enrichment: %d", length(cell_types_use)))
message(paste(cell_types_use, collapse = ", "))


sig_deg <- all_deg_df %>% filter(direction != "NS")


up_genes_per_ct <- sig_deg %>%
  filter(
    direction == "Up in DPN",
    !grepl("^ENSG", gene)  # 过滤ENSG编号基因
  ) %>%
  group_by(gene) %>%
  summarise(
    n_celltypes = n_distinct(cell_type),  # 该基因在多少种细胞类型中上调
    cell_types  = paste(cell_type, collapse = "; "),  # 列出所有细胞类型
    mean_log2FC = mean(avg_log2FC),  # 平均log2倍数变化
    min_padj    = min(p_val_adj),  # 最小校正p值
    .groups = "drop"
  ) %>%
  filter(n_celltypes >= 2) %>%  # 只保留在至少2种细胞类型中上调的基因
  arrange(desc(n_celltypes), desc(mean_log2FC))  # 按细胞类型数和倍数变化排序

# ==============================================
# 跨细胞类型下调基因统计（DPN中低表达）
# ==============================================
dn_genes_per_ct <- sig_deg %>%
  filter(
    direction == "Down in DPN",
    !grepl("^ENSG", gene)
  ) %>%
  group_by(gene) %>%
  summarise(
    n_celltypes = n_distinct(cell_type),
    cell_types  = paste(cell_type, collapse = "; "),
    mean_log2FC = mean(avg_log2FC),
    min_padj    = min(p_val_adj),
    .groups = "drop"
  ) %>%
  filter(n_celltypes >= 2) %>%
  arrange(desc(n_celltypes), desc(abs(mean_log2FC)))  # 下调基因按绝对值排序

message(sprintf("  Multi-celltype UP   genes (>=2 types): %d",
                nrow(up_genes_per_ct)))
message(sprintf("  Multi-celltype DOWN genes (>=2 types): %d",
                nrow(dn_genes_per_ct)))





# 保存所有富集结果
enrich_results <- list()

# ------------------------------------------------
# 富集函数
# ------------------------------------------------
run_enrichment <- function(genes,
                           min_genes = 10,
                           p_cutoff = 0.05,
                           q_cutoff = 0.2) {
  
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  genes <- genes[!grepl("^ENSG", genes)]
  
  if (length(genes) < min_genes) {
    message(sprintf("    Too few genes for enrichment: %d genes", length(genes)))
    return(list(GO_BP = NULL, GO_MF = NULL, GO_CC = NULL, KEGG = NULL))
  }
  
  # SYMBOL -> ENTREZID
  gene_ids <- tryCatch({
    bitr(
      genes,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  }, error = function(e) {
    message("    bitr error: ", e$message)
    return(NULL)
  })
  
  if (is.null(gene_ids) || nrow(gene_ids) == 0) {
    message("    No genes could be converted to ENTREZID.")
    return(list(GO_BP = NULL, GO_MF = NULL, GO_CC = NULL, KEGG = NULL))
  }
  
  gene_ids <- gene_ids %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  entrez_genes <- unique(gene_ids$ENTREZID)
  
  if (length(entrez_genes) < min_genes) {
    message(sprintf("    Too few converted genes: %d ENTREZ IDs", length(entrez_genes)))
    return(list(GO_BP = NULL, GO_MF = NULL, GO_CC = NULL, KEGG = NULL))
  }
  
  # GO BP
  go_bp <- tryCatch({
    enrichGO(
      gene          = entrez_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    )
  }, error = function(e) {
    message("    GO BP error: ", e$message)
    return(NULL)
  })
  
  # GO MF
  go_mf <- tryCatch({
    enrichGO(
      gene          = entrez_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    )
  }, error = function(e) {
    message("    GO MF error: ", e$message)
    return(NULL)
  })
  
  # GO CC
  go_cc <- tryCatch({
    enrichGO(
      gene          = entrez_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    )
  }, error = function(e) {
    message("    GO CC error: ", e$message)
    return(NULL)
  })
  
  # KEGG
  kegg_res <- tryCatch({
    enrichKEGG(
      gene          = entrez_genes,
      organism      = "hsa",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff
    )
  }, error = function(e) {
    message("    KEGG error: ", e$message)
    return(NULL)
  })
  
  # KEGG 转换成 SYMBOL，方便查看
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_res <- tryCatch({
      setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }, error = function(e) {
      kegg_res
    })
  }
  
  return(list(
    GO_BP = go_bp,
    GO_MF = go_mf,
    GO_CC = go_cc,
    KEGG  = kegg_res
  ))
}

# ------------------------------------------------
# 按 cell type 循环富集
# ------------------------------------------------
for (ct in cell_types_use) {
  
  message(sprintf("\n  Processing enrichment for: %s", ct))
  
  ct_safe <- gsub("[/ +]", "_", ct)
  
  up_genes <- sig_deg %>%
    filter(cell_type == ct,
           direction == "Up in DPN") %>%
    pull(gene) %>%
    unique()
  
  dn_genes <- sig_deg %>%
    filter(cell_type == ct,
           direction == "Down in DPN") %>%
    pull(gene) %>%
    unique()
  
  message(sprintf("    Up genes: %d", length(up_genes)))
  message(sprintf("    Down genes: %d", length(dn_genes)))
  
  enrich_up <- run_enrichment(up_genes)
  enrich_dn <- run_enrichment(dn_genes)
  
  enrich_results[[ct]] <- list(
    up   = enrich_up,
    down = enrich_dn
  )
  
  # -------------------------------
  # 内部函数：保存富集结果和图片
  # -------------------------------
  save_enrich_result <- function(enrich_obj, prefix, title_prefix) {
  
    if (is.null(enrich_obj)) return(NULL)
  
  enrich_df <- as.data.frame(enrich_obj)
  
  if (nrow(enrich_df) == 0) {
    message(sprintf("    No enriched terms for %s", prefix))
    return(NULL)
  }
  
  # 保存 CSV
  write.csv(
    enrich_df,
    file.path(pathway_dir, paste0(date_tag, "_", prefix, "_", ct_safe, ".csv")),
    row.names = FALSE
  )
  
  # Dotplot
  p_dot <- dotplot(enrich_obj, showCategory = 20) +
    ggtitle(title_prefix) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(pathway_dir, paste0(date_tag, "_", prefix, "_Dotplot_", ct_safe, ".png")),
    plot = p_dot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Barplot
  p_bar <- barplot(enrich_obj, showCategory = 20) +
    ggtitle(title_prefix) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(pathway_dir, paste0(date_tag, "_", prefix, "_Barplot_", ct_safe, ".png")),
    plot = p_bar,
    width = 10,
    height = 8,
    dpi = 300
  )
  }  
  # Up in DPN
  save_enrich_result(
    enrich_up$GO_BP,
    prefix       = "GO_BP_Up",
    title_prefix = paste0(ct, " - GO BP: Up in DPN")
  )
  
  save_enrich_result(
    enrich_up$GO_MF,
    prefix       = "GO_MF_Up",
    title_prefix = paste0(ct, " - GO MF: Up in DPN")
  )
  
  save_enrich_result(
    enrich_up$GO_CC,
    prefix       = "GO_CC_Up",
    title_prefix = paste0(ct, " - GO CC: Up in DPN")
  )
  
  save_enrich_result(
    enrich_up$KEGG,
    prefix       = "KEGG_Up",
    title_prefix = paste0(ct, " - KEGG: Up in DPN")
  )
  
  # Down in DPN
  save_enrich_result(
    enrich_dn$GO_BP,
    prefix       = "GO_BP_Down",
    title_prefix = paste0(ct, " - GO BP: Down in DPN")
  )
  
  save_enrich_result(
    enrich_dn$GO_MF,
    prefix       = "GO_MF_Down",
    title_prefix = paste0(ct, " - GO MF: Down in DPN")
  )
  
  save_enrich_result(
    enrich_dn$GO_CC,
    prefix       = "GO_CC_Down",
    title_prefix = paste0(ct, " - GO CC: Down in DPN")
  )
  
  save_enrich_result(
    enrich_dn$KEGG,
    prefix       = "KEGG_Down",
    title_prefix = paste0(ct, " - KEGG: Down in DPN")
  )
}

message("\n  Enrichment analysis finished.")


# ================================================================
# Step 5.2: Enrichment for shared multi-celltype DEGs
# ================================================================
message("\n===== Step 5.2: Shared DEG Enrichment =====")

shared_up_genes <- up_genes_per_ct$gene
shared_dn_genes <- dn_genes_per_ct$gene

message(sprintf("  Shared Up genes: %d", length(shared_up_genes)))
message(sprintf("  Shared Down genes: %d", length(shared_dn_genes)))

shared_enrich_up <- run_enrichment(shared_up_genes)
shared_enrich_dn <- run_enrichment(shared_dn_genes)

# 保存函数
save_shared_enrich <- function(enrich_obj, prefix, title_text) {
  
  if (is.null(enrich_obj)) return(NULL)
  
  enrich_df <- as.data.frame(enrich_obj)
  
  if (nrow(enrich_df) == 0) {
    message(sprintf("  No enriched terms for %s", prefix))
    return(NULL)
  }
  
  write.csv(
    enrich_df,
    file.path(pathway_dir, paste0(date_tag, "_Shared_", prefix, ".csv")),
    row.names = FALSE
  )
  
  p_dot <- enrichplot::dotplot(enrich_obj, showCategory = 20) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(pathway_dir, paste0(date_tag, "_Shared_", prefix, "_Dotplot.png")),
    plot = p_dot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Shared Up
save_shared_enrich(
  shared_enrich_up$GO_BP,
  "GO_BP_Up",
  "Shared Up Genes - GO BP"
)

save_shared_enrich(
  shared_enrich_up$KEGG,
  "KEGG_Up",
  "Shared Up Genes - KEGG"
)

# Shared Down
save_shared_enrich(
  shared_enrich_dn$GO_BP,
  "GO_BP_Down",
  "Shared Down Genes - GO BP"
)

save_shared_enrich(
  shared_enrich_dn$KEGG,
  "KEGG_Down",
  "Shared Down Genes - KEGG"
)

# ================================================================
# Step 5.3: GSEA GO-BP by Cell Type
# ================================================================
message("\n===== Step 5.3: GSEA GO-BP Analysis =====")

gsea_dir <- file.path(projectPath, "Output", "DPN_Analysis", "05_Pathway", "GSEA")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

run_gsea_go <- function(deg_df_ct) {
  
  deg_df_ct <- deg_df_ct %>%
    filter(!is.na(avg_log2FC),
           !grepl("^ENSG", gene)) %>%
    arrange(desc(avg_log2FC))
  
  gene_map <- bitr(
    deg_df_ct$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  deg_mapped <- deg_df_ct %>%
    inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    group_by(ENTREZID) %>%
    summarise(avg_log2FC = mean(avg_log2FC), .groups = "drop")
  
  gene_list <- deg_mapped$avg_log2FC
  names(gene_list) <- deg_mapped$ENTREZID
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  if (length(gene_list) < 100) {
    message("    Too few genes for GSEA.")
    return(NULL)
  }
  
  gsea_res <- tryCatch({
    gseGO(
      geneList     = gene_list,
      OrgDb        = org.Hs.eg.db,
      keyType      = "ENTREZID",
      ont          = "BP",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE
    )
  }, error = function(e) {
    message("    GSEA error: ", e$message)
    return(NULL)
  })
  
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_res <- setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  return(gsea_res)
}


# ================================================================
# Step 5.3 (续): 按细胞类型循环运行 GSEA 并保存结果
# ================================================================

gsea_results <- list()

for (ct in cell_types_use) {
  
  message(sprintf("\n  Running GSEA GO-BP for: %s", ct))
  ct_safe <- gsub("[/ +]", "_", ct)
  
  # 提取该细胞类型的全部 DEG（不限方向，GSEA 需要完整排名列表）
  deg_ct <- all_deg_df %>%
    filter(cell_type == ct)
  
  message(sprintf("    Total genes in ranked list: %d", nrow(deg_ct)))
  
  # 运行 GSEA
  gsea_res <- run_gsea_go(deg_ct)
  gsea_results[[ct]] <- gsea_res
  
  if (is.null(gsea_res)) {
    message(sprintf("    ⚠️ No GSEA results for %s, skipping.", ct))
    next
  }
  
  gsea_df <- as.data.frame(gsea_res)
  
  if (nrow(gsea_df) == 0) {
    message(sprintf("    No significant GSEA terms for %s.", ct))
    next
  }
  
  message(sprintf("    ✅ Significant GSEA terms: %d", nrow(gsea_df)))
  
  # ---------------------------
  # 保存 CSV
  # ---------------------------
  write.csv(
    gsea_df,
    file.path(gsea_dir, paste0(date_tag, "_GSEA_GO_BP_", ct_safe, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------
  # 图1: Dotplot（按 NES 着色）
  # ---------------------------
  n_show <- min(20, nrow(gsea_df))
  
  p_dot <- dotplot(gsea_res, showCategory = n_show, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(paste0(ct, " - GSEA GO-BP Dotplot")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(
    file.path(gsea_dir, paste0(date_tag, "_GSEA_Dotplot_", ct_safe, ".png")),
    plot = p_dot,
    width = 12, height = 8, dpi = 300
  )
  
  # ---------------------------
  # 图2: Ridgeplot（密度脊线图）
  # ---------------------------
  p_ridge <- tryCatch({
    ridgeplot(gsea_res, showCategory = min(15, nrow(gsea_df))) +
      ggtitle(paste0(ct, " - GSEA GO-BP Ridgeplot")) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 8)
      )
  }, error = function(e) {
    message("    Ridgeplot error: ", e$message)
    NULL
  })
  
  if (!is.null(p_ridge)) {
    ggsave(
      file.path(gsea_dir, paste0(date_tag, "_GSEA_Ridgeplot_", ct_safe, ".png")),
      plot = p_ridge,
      width = 10, height = 8, dpi = 300
    )
  }
  
  # ---------------------------
  # 图3: 经典 GSEA Running Score Plot（Top 3 通路）
  # ---------------------------
  top_ids <- gsea_df %>%
    arrange(p.adjust) %>%
    head(3) %>%
    pull(ID)
  
  for (i in seq_along(top_ids)) {
    p_gsea <- tryCatch({
      gseaplot2(
        gsea_res,
        geneSetID = top_ids[i],
        title = paste0(ct, "\n", gsea_df$Description[gsea_df$ID == top_ids[i]]),
        pvalue_table = TRUE
      )
    }, error = function(e) {
      message(sprintf("    gseaplot2 error for %s: %s", top_ids[i], e$message))
      NULL
    })
    
    if (!is.null(p_gsea)) {
      ggsave(
        file.path(gsea_dir, paste0(date_tag, "_GSEA_RunningScore_", ct_safe, "_Top", i, ".png")),
        plot = p_gsea,
        width = 8, height = 6, dpi = 300
      )
    }
  }
  
  # ---------------------------
  # 图4: 多通路合并 Running Score（Top 5）
  # ---------------------------
  top5_ids <- gsea_df %>%
    arrange(p.adjust) %>%
    head(min(5, nrow(gsea_df))) %>%
    pull(ID)
  
  if (length(top5_ids) >= 2) {
    p_multi <- tryCatch({
      gseaplot2(
        gsea_res,
        geneSetID = top5_ids,
        pvalue_table = FALSE,
        title = paste0(ct, " - Top GSEA Pathways")
      )
    }, error = function(e) NULL)
    
    if (!is.null(p_multi)) {
      ggsave(
        file.path(gsea_dir, paste0(date_tag, "_GSEA_MultiPath_", ct_safe, ".png")),
        plot = p_multi,
        width = 10, height = 7, dpi = 300
      )
    }
  }
}

# ================================================================
# Step 5.3.2: 汇总所有细胞类型的 GSEA 结果
# ================================================================

# ================================================================
# Step 5.3 (续): 按细胞类型循环运行 GSEA 并保存结果
# ================================================================

gsea_results <- list()

for (ct in cell_types_use) {
  
  message(sprintf("\n  Running GSEA GO-BP for: %s", ct))
  ct_safe <- gsub("[/ +]", "_", ct)
  
  # 提取该细胞类型的全部 DEG（不限方向，GSEA 需要完整排名列表）
  deg_ct <- all_deg_df %>%
    filter(cell_type == ct)
  
  message(sprintf("    Total genes in ranked list: %d", nrow(deg_ct)))
  
  # 运行 GSEA
  gsea_res <- run_gsea_go(deg_ct)
  gsea_results[[ct]] <- gsea_res
  
  if (is.null(gsea_res)) {
    message(sprintf("    ⚠️ No GSEA results for %s, skipping.", ct))
    next
  }
  
  gsea_df <- as.data.frame(gsea_res)
  
  if (nrow(gsea_df) == 0) {
    message(sprintf("    No significant GSEA terms for %s.", ct))
    next
  }
  
  message(sprintf("    ✅ Significant GSEA terms: %d", nrow(gsea_df)))
  
  # ---------------------------
  # 保存 CSV
  # ---------------------------
  write.csv(
    gsea_df,
    file.path(gsea_dir, paste0(date_tag, "_GSEA_GO_BP_", ct_safe, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------
  # 图1: Dotplot（按 NES 着色）
  # ---------------------------
  n_show <- min(20, nrow(gsea_df))
  
  p_dot <- dotplot(gsea_res, showCategory = n_show, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(paste0(ct, " - GSEA GO-BP Dotplot")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(
    file.path(gsea_dir, paste0(date_tag, "_GSEA_Dotplot_", ct_safe, ".png")),
    plot = p_dot,
    width = 12, height = 8, dpi = 300
  )
  
  # ---------------------------
  # 图2: Ridgeplot（密度脊线图）
  # ---------------------------
  p_ridge <- tryCatch({
    ridgeplot(gsea_res, showCategory = min(15, nrow(gsea_df))) +
      ggtitle(paste0(ct, " - GSEA GO-BP Ridgeplot")) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 8)
      )
  }, error = function(e) {
    message("    Ridgeplot error: ", e$message)
    NULL
  })
  
  if (!is.null(p_ridge)) {
    ggsave(
      file.path(gsea_dir, paste0(date_tag, "_GSEA_Ridgeplot_", ct_safe, ".png")),
      plot = p_ridge,
      width = 10, height = 8, dpi = 300
    )
  }
  
  # ---------------------------
  # 图3: 经典 GSEA Running Score Plot（Top 3 通路）
  # ---------------------------
  top_ids <- gsea_df %>%
    arrange(p.adjust) %>%
    head(3) %>%
    pull(ID)
  
  for (i in seq_along(top_ids)) {
    p_gsea <- tryCatch({
      gseaplot2(
        gsea_res,
        geneSetID = top_ids[i],
        title = paste0(ct, "\n", gsea_df$Description[gsea_df$ID == top_ids[i]]),
        pvalue_table = TRUE
      )
    }, error = function(e) {
      message(sprintf("    gseaplot2 error for %s: %s", top_ids[i], e$message))
      NULL
    })
    
    if (!is.null(p_gsea)) {
      ggsave(
        file.path(gsea_dir, paste0(date_tag, "_GSEA_RunningScore_", ct_safe, "_Top", i, ".png")),
        plot = p_gsea,
        width = 8, height = 6, dpi = 300
      )
    }
  }
  # ---------------------------
  # 图4: 多通路合并 Running Score（Top 5）
  # ---------------------------
  top5_ids <- gsea_df %>%
    arrange(p.adjust) %>%
    head(min(5, nrow(gsea_df))) %>%
    pull(ID)
  
  if (length(top5_ids) >= 2) {
    p_multi <- tryCatch({
      gseaplot2(
        gsea_res,
        geneSetID = top5_ids,
        pvalue_table = FALSE,
        title = paste0(ct, " - Top GSEA Pathways")
      )
    }, error = function(e) NULL)
    
    if (!is.null(p_multi)) {
      ggsave(
        file.path(gsea_dir, paste0(date_tag, "_GSEA_MultiPath_", ct_safe, ".png")),
        plot = p_multi,
        width = 10, height = 7, dpi = 300
      )
    }
  }
}
# ================================================================
# Step 5.3.2: 汇总所有细胞类型的 GSEA 结果（稳健版）
# ================================================================
message("\n  Combining all GSEA results...")
all_gsea_df <- data.frame()
for (ct in names(gsea_results)) {
  res <- gsea_results[[ct]]
  if (!is.null(res)) {
    res_df <- as.data.frame(res)
    if (nrow(res_df) > 0) {
      res_df$cell_type <- ct
      all_gsea_df <- rbind(all_gsea_df, res_df)
    }
  }
}
if (nrow(all_gsea_df) > 0) {
  write.csv(
    all_gsea_df,
    file.path(gsea_dir, paste0(date_tag, "_GSEA_GO_BP_AllCellTypes_Combined.csv")),
    row.names = FALSE
  )
  message(sprintf(
    "  ✅ Combined GSEA table saved: %d terms across %d cell types.",
    nrow(all_gsea_df),
    length(unique(all_gsea_df$cell_type))
  ))
  # ------------------------------------------------
  # 选取多个细胞类型中都显著的通路
  # ------------------------------------------------
  pathway_counts <- all_gsea_df %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(
      n_ct = dplyr::n_distinct(cell_type),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_ct >= 2) %>%
    dplyr::arrange(dplyr::desc(n_ct))
  if (nrow(pathway_counts) > 0) {
    top_shared_pathways <- pathway_counts %>%
      head(min(30, nrow(pathway_counts))) %>%
      dplyr::pull(Description)
    # ------------------------------------------------
    # 构建 NES 热图矩阵
    # ------------------------------------------------
  heatmap_data <- all_gsea_df %>%
    dplyr::filter(Description %in% top_shared_pathways) %>%
    dplyr::select(Description, cell_type, NES) %>%
    dplyr::mutate(NES = as.numeric(NES)) %>%
    dplyr::group_by(Description, cell_type) %>%
    dplyr::summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = cell_type,
      values_from = NES
    ) %>%
  tibble::column_to_rownames("Description")
    heatmap_mat <- as.matrix(heatmap_data)
    mode(heatmap_mat) <- "numeric"
    # ------------------------------------------------
    # 清理 NA/Inf，避免 pheatmap 报错
    # 规则：
    # 1) 去掉全 NA 的行
    # 2) 去掉方差为 0 的行
    # 3) 保留至少 2 个非 NA 值的通路
    # ------------------------------------------------
    heatmap_mat <- heatmap_mat[rowSums(!is.na(heatmap_mat)) >= 2, , drop = FALSE]
    if (nrow(heatmap_mat) > 0 && ncol(heatmap_mat) > 1) {
      
      # 对剩余 NA 做行均值填充
      for (i in seq_len(nrow(heatmap_mat))) {
        row_vals <- heatmap_mat[i, ]
        if (any(is.na(row_vals))) {
          row_mean <- mean(row_vals, na.rm = TRUE)
          if (is.finite(row_mean)) {
            row_vals[is.na(row_vals)] <- row_mean
            heatmap_mat[i, ] <- row_vals
          }
        }
      }
      # 再次去掉仍然有问题的行
      heatmap_mat <- heatmap_mat[is.finite(rowMeans(heatmap_mat)), , drop = FALSE]
      heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
      if (nrow(heatmap_mat) > 1 && ncol(heatmap_mat) > 1) {
        if (requireNamespace("pheatmap", quietly = TRUE)) {
          library(pheatmap)
          png(
            file.path(gsea_dir, paste0(date_tag, "_GSEA_NES_Heatmap_SharedPathways.png")),
            width = 12,
            height = 10,
            units = "in",
            res = 300
          )
          pheatmap::pheatmap(
            heatmap_mat,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            color = colorRampPalette(c("blue", "white", "red"))(100),
            na_col = "grey90",
            main = "GSEA NES - Shared Pathways Across Cell Types",
            fontsize_row = 7,
            fontsize_col = 9,
            border_color = NA
          )
          dev.off()
          message("  ✅ NES Heatmap saved.")
        }
      }
    }
  }
}
message("\n===== Step 5.3 GSEA Analysis Complete =====")message("\n===== Step 5.3 GSEA Analysis Complete =====")



                        

# ================================================================
# 步骤 6：细胞间通讯分析（CellChat）
# ================================================================
message("\n===== Step 6: Cell-Cell Communication (CellChat) =====")

library(CellChat)
library(NMF)

run_cellchat <- function(seurat_obj, group_name) {
  message(sprintf("  Creating CellChat object for: %s", group_name))
  data_input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta_input <- data.frame(
    labels   = seurat_obj$cell_type_manual,
    row.names = colnames(seurat_obj)
  )
  cc <- createCellChat(object = data_input, meta = meta_input, group.by = "labels")
  CellChatDB <- CellChatDB.human
  cc@DB <- CellChatDB
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  return(cc)
}

# 分别创建两组的 CellChat 对象
cc_control <- run_cellchat(subset(merged_obj, group == "Control"), "Control")
cc_dpn     <- run_cellchat(subset(merged_obj, group == "DPN"),     "DPN")

# 保存
save(cc_control, file = file.path(projectPath, "Data", "CellChat_Control.Rdata"))
save(cc_dpn,     file = file.path(projectPath, "Data", "CellChat_DPN.Rdata"))

# ---- 6.1 基本通讯强度比较 ----
cc_list <- list(Control = cc_control, DPN = cc_dpn)
cc_merged <- mergeCellChat(cc_list, add.names = names(cc_list))

# 交互数量 & 强度
png(file.path(projectPath, "Output", "DPN_Analysis", "06_CellChat",
              paste0(date_tag, "_CellChat_Count_Strength.png")),
    width = 1400, height = 600, res = 150)
par(mfrow = c(1, 2))
compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2))
compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2),
                    measure = "weight")
dev.off()

# ---- 6.2 差异通讯 Chord Diagram ----
png(file.path(projectPath, "Output", "DPN_Analysis", "06_CellChat",
              paste0(date_tag, "_CellChat_DiffInteraction_Chord.png")),
    width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2))
netVisual_diffInteraction(cc_merged, weight.scale = TRUE)
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight")
dev.off()

# ---- 6.3 气泡图：差异通讯通路 ----
p_bubble <- netVisual_bubble(cc_merged, sources.use = NULL,
                               targets.use = NULL,
                               comparison = c(1, 2), angle.x = 45)
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "06_CellChat",
                 paste0(date_tag, "_CellChat_Bubble.png")),
       plot = p_bubble, width = 16, height = 10, dpi = 300)

# ---- 6.4 通路信号流分析 ----
# 对比两组信号通路强度
rankNet_df <- rankNet(cc_merged, mode = "comparison",
                       stacked = TRUE, do.stat = TRUE)
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "06_CellChat",
                 paste0(date_tag, "_CellChat_PathwayRank.png")),
       plot = rankNet_df, width = 10, height = 8, dpi = 300)


# ================================================================
# 步骤 7：拟时序/轨迹分析（Monocle3）
# ================================================================
message("\n===== Step 7: Trajectory Analysis (Monocle3) =====")

library(monocle3)
library(SeuratWrappers)

# 选择最相关的细胞类型（T 细胞亚群做轨迹）
traj_celltypes <- c("CD4+ T naive", "CD4+ T memory", "CD8+ T naive",
                     "CD8+ T effector", "Treg")
traj_obj <- subset(merged_obj,
                    cell_type_manual %in% traj_celltypes)

# 转换为 CDS 对象
cds <- as.cell_data_set(traj_obj)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

# 以 naive T 细胞作为起点
root_cells <- colnames(traj_obj)[traj_obj$cell_type_manual == "CD4+ T naive"]
cds <- order_cells(cds, root_cells = root_cells)

# ---- 7.1 轨迹图 ----
p_traj_ct <- plot_cells(cds, color_cells_by = "cell_type_manual",
                          label_groups_by_cluster = FALSE,
                          label_leaves = FALSE, label_branch_points = FALSE) +
  scale_color_manual(values = cell_type_colors[traj_celltypes]) +
  ggtitle("T Cell Trajectory by Cell Type")

p_traj_pt <- plot_cells(cds, color_cells_by = "pseudotime",
                          label_cell_groups = FALSE,
                          label_leaves = FALSE) +
  scale_color_viridis_c() +
  ggtitle("T Cell Pseudotime")

p_traj_grp <- plot_cells(cds, color_cells_by = "group",
                           label_cell_groups = FALSE) +
  scale_color_manual(values = group_colors) +
  ggtitle("T Cell Trajectory by Group")

p_traj <- p_traj_ct | p_traj_pt | p_traj_grp
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "07_Trajectory",
                 paste0(date_tag, "_Trajectory_Tcell.png")),
       plot = p_traj, width = 22, height = 7, dpi = 300)

# ---- 7.2 拟时序相关基因检测 ----
traj_obj$pseudotime <- pseudotime(cds)[colnames(traj_obj)]

# 在两组间比较伪时序分布
p_pt_group <- ggplot(traj_obj@meta.data %>% filter(!is.na(pseudotime)),
                      aes(x = pseudotime, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ cell_type_manual, scales = "free") +
  theme_classic() +
  labs(title = "Pseudotime Distribution: Control vs DPN")
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "07_Trajectory",
                 paste0(date_tag, "_Pseudotime_Distribution.png")),
       plot = p_pt_group, width = 14, height = 8, dpi = 300)

# ---- 7.3 Graph-test 找拟时序差异基因 ----
graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
sig_traj_genes <- graph_test_res %>%
  filter(q_value < 0.05, morans_I > 0.1) %>%
  arrange(desc(morans_I))

write.csv(sig_traj_genes,
          file.path(projectPath, "Output", "DPN_Analysis", "07_Trajectory",
                    paste0(date_tag, "_Trajectory_SignificantGenes.csv")),
          row.names = FALSE)


# ================================================================
# 步骤 8：转录因子活性分析（decoupleR）
# ================================================================
message("\n===== Step 8: Transcription Factor Activity (decoupleR) =====")

library(decoupleR)

# 获取 CollecTRI 转录因子调控网络
net <- get_collectri(organism = "human", split_complexes = FALSE)

# 对每种细胞类型计算 TF 活性
for (ct in names(all_deg_results)) {

  deg_ct <- all_deg_results[[ct]] %>%
    filter(!is.na(log2FoldChange)) %>%
    select(gene, log2FoldChange, stat) %>%
    column_to_rownames("gene")

  # 用 deg stat 作为输入矩阵
  stat_mat <- as.matrix(deg_ct[, "stat", drop = FALSE])
  colnames(stat_mat) <- ct

  # 运行 ULM（univariate linear model）
  tf_acts <- tryCatch(
    run_ulm(mat = t(stat_mat), net = net, .source = "source",
             .target = "target", .mor = "mor", minsize = 5),
    error = function(e) NULL
  )
  if (is.null(tf_acts)) next

  # 筛选显著 TF
  sig_tfs <- tf_acts %>%
    filter(statistic == "ulm", p_value < 0.05) %>%
    arrange(desc(abs(score))) %>%
    head(30)

  if (nrow(sig_tfs) == 0) next

  ct_safe <- gsub("[/ +]", "_", ct)
  write.csv(sig_tfs,
            file.path(projectPath, "Output", "DPN_Analysis", "08_TF",
                      paste0(date_tag, "_TF_Activity_", ct_safe, ".csv")),
            row.names = FALSE)

  # 可视化 TF 活性
  p_tf <- ggplot(sig_tfs, aes(x = reorder(source, score), y = score,
                               fill = score > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                       labels = c("Activated in DPN", "Repressed in DPN")) +
    theme_classic() +
    labs(title = paste0("TF Activity: ", ct, " (DPN vs Control)"),
         x = NULL, y = "TF Activity Score", fill = NULL)
  ggsave(file.path(projectPath, "Output", "DPN_Analysis", "08_TF",
                   paste0(date_tag, "_TF_Activity_", ct_safe, ".png")),
         plot = p_tf, width = 9, height = 7, dpi = 300)
}


# ================================================================
# 步骤 9：基因共表达模块分析（hdWGCNA）
# ================================================================
message("\n===== Step 9: Co-expression Module Analysis (hdWGCNA) =====")

library(hdWGCNA)
library(WGCNA)

# 以免疫细胞为例（Monocytes）
wgcna_ct <- "CD14+ Monocytes"
wgcna_obj <- subset(merged_obj, cell_type_manual == wgcna_ct)

# 设置 hdWGCNA
wgcna_obj <- SetupForWGCNA(wgcna_obj,
                             gene_select = "fraction",
                             fraction = 0.05,
                             wgcna_name = wgcna_ct)

# Metacell 聚合
wgcna_obj <- MetacellsByGroups(wgcna_obj,
                                group.by = c("cell_type_manual", "orig.ident"),
                                reduction = "umap_harmony",
                                k = 25, max_shared = 10)
wgcna_obj <- NormalizeMetacells(wgcna_obj)

# 构建表达矩阵
wgcna_obj <- SetDatExpr(wgcna_obj,
                         group_name  = wgcna_ct,
                         group.by    = "cell_type_manual",
                         assay = "RNA", slot = "data")

# 测试软阈值
wgcna_obj <- TestSoftPowers(wgcna_obj, networkType = "signed")
p_sp <- PlotSoftPowers(wgcna_obj)
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "09_WGCNA",
                 paste0(date_tag, "_WGCNA_SoftPower.png")),
       plot = p_sp, width = 12, height = 8, dpi = 300)

# 构建共表达网络
wgcna_obj <- ConstructNetwork(wgcna_obj, soft_power = 9,
                               setDatExpr = FALSE, tom_type = "signed")
wgcna_obj <- ModuleEigengenes(wgcna_obj, group.by.vars = "orig.ident")
wgcna_obj <- ModuleConnectivity(wgcna_obj)

# 关联 DPN 表型
wgcna_obj$DPN_binary <- ifelse(wgcna_obj$group == "DPN", 1, 0)
wgcna_obj <- ModuleTraitCorrelation(wgcna_obj, traits = c("DPN_binary"))

p_mt_cor <- PlotModuleTraitCorrelation(wgcna_obj, label = "fdr",
                                        label_symbol = "stars",
                                        trait_name = "DPN_binary")
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "09_WGCNA",
                 paste0(date_tag, "_WGCNA_ModuleTraitCor.png")),
       plot = p_mt_cor, width = 10, height = 7, dpi = 300)

# 导出 hub genes
hub_genes <- GetHubGenes(wgcna_obj, n_hubs = 10)
write.csv(hub_genes,
          file.path(projectPath, "Output", "DPN_Analysis", "09_WGCNA",
                    paste0(date_tag, "_WGCNA_HubGenes.csv")),
          row.names = FALSE)


# ================================================================
# 步骤 10：细胞状态评分（DPN 疾病特征打分）
# ================================================================
message("\n===== Step 10: Cell State Scoring =====")

# ---- 10.1 用 DEG 构建 DPN 特征基因集评分 ----
# 从步骤4 取跨细胞类型 top 上调基因作为 DPN 特征
dpn_signature <- head(up_genes_per_ct$gene, 50)
dpn_signature <- dpn_signature[dpn_signature %in% rownames(merged_obj)]

ctrl_signature <- head(dn_genes_per_ct$gene, 50)
ctrl_signature <- ctrl_signature[ctrl_signature %in% rownames(merged_obj)]

merged_obj <- AddModuleScore(merged_obj,
                              features = list(dpn_signature),
                              name     = "DPN_Score",
                              seed     = 42)
merged_obj <- AddModuleScore(merged_obj,
                              features = list(ctrl_signature),
                              name     = "Control_Score",
                              seed     = 42)

# ---- 10.2 评分 UMAP 可视化 ----
p_score_umap <- FeaturePlot(merged_obj, reduction = "umap_harmony",
                             features = c("DPN_Score1", "Control_Score1"),
                             split.by = "group",
                             pt.size = 0.1, order = TRUE, ncol = 2) &
  scale_color_gradient2(low = "#4DBBD5", mid = "white",
                         high = "#E64B35", midpoint = 0)
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "10_CellScore",
                 paste0(date_tag, "_DPN_Score_UMAP.png")),
       plot = p_score_umap, width = 18, height = 8, dpi = 300)

# ---- 10.3 评分按细胞类型比较 ----
score_df <- merged_obj@meta.data %>%
  select(cell_type_manual, group, DPN_Score1, Control_Score1)

p_score_vln <- ggplot(score_df, aes(x = cell_type_manual, y = DPN_Score1,
                                     fill = group)) +
  geom_violin(scale = "width", alpha = 0.8, trim = TRUE,
              position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, outlier.size = 0,
               position = position_dodge(0.8)) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(aes(group = group), label = "p.signif",
                     method = "wilcox.test") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DPN Signature Score by Cell Type",
       x = NULL, y = "DPN Score")
ggsave(file.path(projectPath, "Output", "DPN_Analysis", "10_CellScore",
                 paste0(date_tag, "_DPN_Score_Violin.png")),
       plot = p_score_vln, width = 14, height = 6, dpi = 300)


# ================================================================
# 步骤 11：汇总图表
# ================================================================
message("\n===== Step 11: Summary Visualization =====")

# ---- 11.1 整合摘要：关键指标一览 ----
summary_table <- data.frame(
  Analysis         = c("Total cells", "Control cells", "DPN cells",
                        "Cell types identified",
                        "Total DEGs", "Up in DPN (multi-CT)",
                        "Down in DPN (multi-CT)"),
  Value            = c(
    ncol(merged_obj),
    sum(merged_obj$group == "Control"),
    sum(merged_obj$group == "DPN"),
    length(levels(merged_obj$cell_type_manual)),
    nrow(all_deg_df %>% filter(direction != "NS")),
    nrow(up_genes_per_ct),
    nrow(dn_genes_per_ct)
  )
)
write.csv(summary_table,
          file.path(projectPath, "Output", "DPN_Analysis", "11_Summary",
                    paste0(date_tag, "_Analysis_Summary.csv")),
          row.names = FALSE)
message("\n  ===== Analysis Complete =====")
print(summary_table)

# ---- 11.2 保存最终 Seurat 对象 ----
save(merged_obj,
     file = file.path(projectPath, "Data",
                       paste0("09_", date_tag,
                              "_6Samples_Final_WithGroupInfo.Rdata")))

message(sprintf("\n  All outputs saved to: %s",
                file.path(projectPath, "Output", "DPN_Analysis")))
