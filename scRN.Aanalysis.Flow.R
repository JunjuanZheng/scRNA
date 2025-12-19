# -*- coding: UTF-8 -*-  
library(dplyr)  
library(ggplot2)  
library(patchwork)  
library(Seurat)  
library(clustree)  # 用于评估聚类稳定性  
library(DoubletFinder)
library(clustree)  
library(igraph)  
library(cluster)  
library(bluster)
library(SeuratWrappers)
# 设置工作目录  
projectPath <- "/mnt2/wanggd_group/zjj/BGCscRNA/WangHuishan/PJ23062901422"  
setwd(projectPath)  

##### 1.数据导入并创建seurat对象 #####  

# 获取数据文件夹下的所有样本文件列表  
samples <- list.files(paste0(projectPath, "/Data/QuantifyRawdata/"))  
seurat_list <- list()  
sample_ids <- c()  # 改为向量而不是列表  

# 读取每个样本的10x数据并创建Seurat对象  
for (sample in samples) {  
  tryCatch({  
    # 拼接文件路径  
    data.path <- file.path(projectPath, "/Data/QuantifyRawdata", sample,   
                           paste0(sample, "_outs"), "filtered_cell_gene_matrix")  
    
    # 检查路径是否存在  
    if (!dir.exists(data.path)) {  
      warning(paste("Directory not found:", data.path))  
      next  
    }  
    # 读取数据  
    seurat_data<-Read10X(data.dir = data.path)  
    # 提取样本ID  
    sample_id <- sample  
    
    # 创建Seurat对象  
    seurat_obj <- CreateSeuratObject(counts = seurat_data,  
                                    project = sample_id,  
                                    min.features = 0,  
                                    min.cells = 0) 
								
	seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

    vln_plot_nCount <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + 
    ggtitle("nCount")

    vln_plot_nFeature <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + 
    ggtitle("nFeature")

    vln_plot_percent_mt <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0) + 
    ggtitle("percent.mt")

    # 将小提琴图组合起来显示
    combined_vln_plot <- vln_plot_nCount | vln_plot_nFeature | vln_plot_percent_mt
	# 保存小提琴图到文件
    output_file <- file.path(projectPath, "Output", paste0(sample,".nFeature_nCount_percent_mt.ViolinPlots.png"))
    ggsave(output_file, plot = combined_vln_plot, width = 12, height = 12, dpi = 300)
    # 提示保存成功
    message(paste("Violin plots saved to:", output_file))
	
    # 将Seurat对象添加到列表中  
    seurat_list[[sample_id]] <- seurat_obj  
    sample_ids <- c(sample_ids, sample_id)  
    
  }, error = function(e) {  
    message(paste("Error processing sample:", sample))  
    message(e)  
  })  
}  

# 保存初始对象  
save(seurat_list,   
    file = file.path(projectPath, "Data", "01_2021216_4Samples_ObjectOriginal.Rdata"))  

load(file.path(projectPath, "Data", "2021216_4Samples_ObjectOriginal.Rdata"))
# 创建marker_list并设置PDF输出  
marker_list <- list()  
#clustree_pdf <- file.path(projectPath, "Output", "Clustree_Plots.pdf") 
#pdf(clustree_pdf, width = 12, height = 12)
 
umap_pdf <- file.path(projectPath, "Output", "All_doubletFinderPlots_.pdf")  
pdf(umap_pdf, width = 10, height = 8) 
# 创建双细胞检测结果目录
#doublet_dir <- file.path(projectPath, "Output/DoubletFinder")
#dir.create(doublet_dir, showWarnings = FALSE, recursive = TRUE)
#pdf(file.path(doublet_dir, paste0("20250220_45Samples_CanPBMC_Object_doublets.pdf")), width = 10, height = 8)

# 处理每个样本 
for (i in seq_along(seurat_list)) { 
  obj_name <- names(seurat_list)[i]  
  current_obj <- seurat_list[[obj_name]]  
  
  message(sprintf("Processing %d/%d: %s", i, length(seurat_list), obj_name))  
  
  # 基础预处理  
# 基础预处理（到PCA阶段）
  current_obj <- current_obj %>%  
    NormalizeData() %>%  
    FindVariableFeatures() %>%  
    ScaleData() %>%  
    RunPCA(verbose = FALSE)
  
  # 自动参数优化
  sweep.res <- paramSweep(current_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ############## 返回 NULL，意味着它无法确定这个最佳点。 下面画图手动确定拐点
  # 绘制 BCmetric 与 pK 的关系图
  p0 = ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(x = "pK Value", y = "BCmetric", title = "BCmetric vs. pK")
  ggsave(file.path(projectPath, "Output",paste0(obj_name,"_All_doubletFinderPlots.pdf")), plot = p0, width = 12, height = 12, dpi = 300)

  # 展示了 BCmetric（Y轴）随着 pK 值（X轴）变化的趋势。
  # BCmetric 可以理解为模型区分单细胞和双细胞的“信心分数”或“区分度”，
  # 我们的目标就是找到让这个分数最高的 pK 值
   
  # 选择最佳pK值
  
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  #pK <- 0.25
  # 计算预期双细胞数量（按细胞数的2.5%估算，可根据实验调整）
  nExp_poi <- round(0.025 * ncol(current_obj)) 
  
  # 运行DoubletFinder
  current_obj <- doubletFinder(
    current_obj, 
    PCs = 1:10, 
    pN = 0.25, 
    pK = pK, 
    nExp = nExp_poi, 
    #reuse.pANN = FALSE, 
    sct = FALSE
  )
  
  # 提取双细胞分类结果
  df.class <- paste("DF.classifications_0.25", pK, nExp_poi, sep = "_")
  current_obj$Doublet_Classification <- current_obj@meta.data[[df.class]]
  
  # 添加双胞率计算和标注  
  doublet_count <- sum(current_obj$Doublet_Classification == "Doublet")  
  total_cells <- ncol(current_obj)  
  doublet_rate <- round(doublet_count / total_cells * 100, 2)  
  # 生成带统计信息的绘图  
  doublet_plot <- DimPlot(current_obj,   
                         group.by = "Doublet_Classification",  
                         reduction = "pca",  
                         cols = c("gray90", "red")) +  # 设置颜色  
    labs(title = paste0(obj_name, " Doublet Detection")) +  
    annotate("text",  
             x = Inf, y = Inf,  # 右上角定位  
             label = paste0("Doublet Rate: ", doublet_rate, "%\n",  
                          "Predicted: ", nExp_poi, " (", doublet_count, " found)"),  
             hjust = 1.1, vjust = 1.1,  # 微调位置  
             size = 5,  
             color = "darkred",  
             fontface = "bold") +  
    theme(  
      plot.title = element_text(size=14, face="bold"),  
      legend.position = c(0.8, 0.2)  # 调整图例位置  
    )  
  print(doublet_plot)
  # 过滤双细胞
  singlet_cells <- colnames(current_obj)[current_obj$Doublet_Classification == "Singlet"]
  current_obj <- subset(current_obj, cells = singlet_cells) %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)  
  
  seurat_list[[obj_name]] <- current_obj  
}
save(seurat_list,   
    file = file.path(projectPath, "Data", "02_2021216_4Samples_ObjectOriginal_AfterDF.Rdata"))  
  

#unintegrated <- merge(seurat_list[[1]], seurat_list[2:4]) 

unintegrated <- merge(seurat_list[[1]], 
                      y = list(seurat_list[[2]], seurat_list[[4]]))
unintegrated <- NormalizeData(unintegrated)
unintegrated <- FindVariableFeatures(unintegrated)
unintegrated <- ScaleData(unintegrated)
unintegrated <- RunPCA(unintegrated)
unintegrated <- FindNeighbors(unintegrated, dims = 1:30)
unintegrated <- RunUMAP(unintegrated, dims = 1:30,reduction.name = "unintegrated")

save(unintegrated,file = file.path(projectPath, "Data", "03_2021216_3Samples_ObjectOriginal_AfterDF_MergedObj.Rdata"))



#########################################################################
############## ------------------------ 在R里面可以实现的整合方法
integration_info <- list(  
  list(reduction = "integrated.cca", cluster.name = "cca_clusters",method.use ='CCAIntegration',umap.reduction = 'umap.cca'),  
  list(reduction = "integrated.rpca", cluster.name = "rpca_clusters",method.use ='RPCAIntegration',umap.reduction = 'umap.rpca'),  
  list(reduction = "harmony", cluster.name = "harmony_clusters",method.use ='HarmonyIntegration',umap.reduction = 'umap.harmony'),  
  list(reduction = "integrated.jointPCA", cluster.name = "jointPCA_clusters",method.use = 'JointPCAIntegration',umap.reduction = 'umap.jointPCA'),
  list(reduction = "fastMNN", cluster.name = "fastMNN_clusters",method.use ='FastMNNIntegration',umap.reduction = 'umap.fastMNN'),  
  list(reduction = "integrated.scvi", cluster.name = "scVI_clusters",method.use = 'scVIIntegration')    
)

obj=unintegrated
#objIntegrated_List = list() 

i=5
red.use <- integration_info[[i]]$reduction  
clst.name <- integration_info[[i]]$cluster.name  
method.use <- integration_info[[i]]$method.use
umap.reduction <- integration_info[[i]]$umap.reduction

objIntegrated<- IntegrateLayers(
		object = obj,   
		method = method.use,  
		orig.reduction = "pca",       # 原始降维的名称  
		new.reduction = red.use,  
		verbose = FALSE  
    )
# re-join layers after integration
objIntegrated[["RNA"]] <- JoinLayers(objIntegrated[["RNA"]])

objIntegrated <- FindNeighbors(objIntegrated, reduction = red.use, dims = 1:30)
objIntegrated <- RunUMAP(objIntegrated, dims = 1:30, reduction = red.use,reduction.name = umap.reduction)

resolutions <- c(0.025, 0.03,0.04,0.05, 0.1, 0.15, 0.2, 0.25)  
for (j in seq_along(resolutions)) {
  res <- resolutions[j]  
  objIntegrated <- FindClusters(objIntegrated, resolution = 0.04,algorithm = 1)
}
  
clustree_plot <- clustree(objIntegrated, prefix = "RNA_snn_res.") +  
   ggtitle(paste0(obj_name, " Clustering Tree")) +  
   theme(plot.title = element_text(size = 14, face = "bold"))  
print(clustree_plot) 
ggsave(file.path(projectPath, "Output","2021216_3Samples_clustreePlot.pdf"), plot = clustree_plot, width = 12, height = 12, dpi = 300)
 
  


pdf(file.path(projectPath, "Output","2021216_3Samples_DimPlot.pdf"))
Idents(objIntegrated) = objIntegrated$RNA_snn_res.0.1
umap_plot <- DimPlot(  
    objIntegrated,  
    reduction = umap.reduction ,  
    label = TRUE  
  ) +  
    ggtitle(paste0(obj_name, " (Final resolution = ", '0.1', ")"))  
print(umap_plot)  					   

CellDimPlot(
  srt = objIntegrated,
  group.by = c("orig.ident", "RNA_snn_res.0.1"),
  reduction = umap.reduction ,
  #theme_use = "theme_blank"
)

CellDimPlot(
  srt = unintegrated,
  group.by = c("orig.ident"),
  reduction = 'unintegrated',
  #theme_use = "theme_blank"
)
dev.off()

objIntegrated$cell_type = objIntegrated$RNA_snn_res.0.1

#### 注释
library(celldex)
library(SingleR)
library(BiocParallel)
library(ggplot2)

library(DeepCellSeek)
library(Seurat)
markers_df <- FindAllMarkers(object = objIntegrated)
Sys.setenv(USE_DEEPCELLSEEK_API = "TRUE")

annotations <- deepcellseek_celltype(
  input = markers_df,
  tissuename = "Brain",  # Providing tissue context is important
  species = "Mouse",
  model = "gpt-5"  # Choose any supported model
)
objIntegrated@meta.data$LLM_Annotation <- as.factor(annotations[as.character(Idents(objIntegrated))])
DimPlot(objIntegrated, group.by = "LLM_Annotation", label = TRUE, repel = TRUE)



  
pdf(file.path(projectPath, "Output","2021216_3Samples_MarkerGenens.HeatmapPlot.pdf"),height=12) 
objIntegrated.markers<-FindAllMarkers(objIntegrated,only.pos=TRUE,min.pct=0.25)
#top10<-objIntegrated.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)

objIntegrated.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10 
set.seed(123) # 设置随机种子
objIntegrated_sampled <- subset(objIntegrated, downsample = 2000)

DoHeatmap(objIntegrated_sampled, features = top10$gene) + NoLegend()
dev.off()



# 根据 Marker 基因创建细胞类型注释
Idents(objIntegrated) = objIntegrated$RNA_snn_res.0.1

# 定义更新后的 Marker 列表
all_markers <- list(
  # 星形胶质细胞
  Astrocyte = unique(c("Aqp4", "Gja1", "Atp1b2", "Aldoc", "Clu", "Slc1a3", "Mt3", 
                       "Gfap", "S100b", "Slc1a2", "Aldh1l1", "Fabp7", "Glul", "Sox9", "Cd44", "Sparcl1")),
  # 少突胶质细胞
  Oligodendrocytes = unique(c("Mbp", "Plp1", "Mog", "Cldn11", "Mobp", 
                              "Cnp", "Mag", "Rtn4", "Slc44a1", "Opalin", "Qk", "Tf")), 
  # 少突前体细胞
  OPC = unique(c("Pdgfra", "Cspg4", "Olig2", "Sox10")), 
  # 小胶质细胞
  Microglia = unique(c("C1qb", "Ctss", "Tyrobp", "Hexb", "Fcer1g", "Cx3cr1", 
                       "Aif1", "P2ry12", "Tmem119", "Csf1r", "C3", "Itgam", "Trem2", 
                       "Ly86", "Cd68", "Cd14", "Spi1")),
  # 兴奋性神经元
  Excitatory_Neuron = unique(c("Slc17a7", "Slc17a6", "Neurod2", "Tbr1", "Gria1", 
                               "Syt1", "Snap25", "Rbfox3", "Neurod6", "Map2", "Tubb3", 
                               "Grin1", "Camk2a", "Syn1", "Elavl4")),
  # 抑制性神经元（包括 GABA能神经元）
  Inhibitory_Neuron = unique(c("Gad1", "Gad2", "Lhx6", "Nkx2-1", "Npy", "NPY51", "NR2F2")),
  # 神经元前体细胞
  Neuronal_Progenitors = unique(c("Ascl1", "Dcx", "Sox11", "Neurod1", "Prox1", 
                                  "Pax6", "Hes5", "Vim", "Eomes", "Emx2", "Nes")),  
  # 红细胞
  Erythrocyte = unique(c("Hbb-bt", "Hbb-bs", "Hba-a1", "Hba-a2")),
  # 神经干细胞/增殖细胞
  NSC_Proliferating = unique(c("Mki67", "Top2a", "Pclaf", "Rrm2"))
)

# 打印调整后的 Marker 列表，确保没有重复基因
all_markers
p_dot <- DotPlot(objIntegrated, 
                 features = all_markers,
                 group.by = "RNA_snn_res.0.1",  # 或用 "cell_type" 如果已添加注释
                 cols = c("lightgrey", "red"),
                 dot.scale = 6) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.position = "right") +
  ggtitle("Cell Type Markers by Cluster")
pdf(file.path(projectPath, "Output","2021216_3Samples_MarkerGenens.DotPlot.pdf"),height=4,width=25) 
print(p_dot)
dev.off()



table(Idents(objIntegrated))
Idents(objIntegrated) = 'RNA_snn_res.0.1'
cluster_annotation <- c(
  "0" = "C0 unknown",
  "1" = "NSC_Proliferating",     #*
  "2" = "Astrocyte",             #*
  "3" = "Excitatory_Neuron",     #*
  "4" = "Oligodendrocytes",      #* 
  "5" = "Astrocyte",             #*
  "6" = "Microglia",             #*
  "7" = "C7 unknown",
  "8" = "Erythrocyte",           #*
  "9" = "Astrocyte",             #*
  "10" = "Unknown"               #*
)
# 添加到 Seurat 对象
objIntegrated =  RenameIdents(objIntegrated, cluster_annotation)
objIntegrated[['cell_type']] <-Idents(objIntegrated)

# 可视化
p1 <- DimPlot(objIntegrated, 
              group.by = "cell_type", 
              reduction = "umap.fastMNN",
              label = TRUE,
              repel = TRUE,
              label.size = 4) +
  ggtitle("Cell Type Annotation") +
  theme(legend.position = "right")
pdf(file.path(projectPath, "Output","2021216_3Samples_MarkerGenens.Umap.pdf")) 
print(p1)
#用阈值/打分把 astrocyte 先标出来    
astro_markers <- c("Aqp4", "Gja1", "Atp1b2", "Aldoc", "Clu", "Slc1a3", "Mt3", 
                       "Gfap", "S100b", "Slc1a2", "Aldh1l1", "Fabp7", "Glul", "Sox9", "Cd44", "Sparcl1")
objIntegrated <- AddModuleScore(objIntegrated, features = list(astro_markers), name = "AstroScore")
VlnPlot(objIntegrated, features = "AstroScore1", group.by = "RNA_snn_res.0.1")
FeaturePlot(objIntegrated, reduction="umap.fastMNN", features="AstroScore1")
	
############
#在 UMAP 上确认星形 marker 是否集中成簇
FeaturePlot(
  objIntegrated,
  reduction = "umap.fastMNN",
  features = astro_markers,
  ncol = 4
) 
dev.off()     


save(objIntegrated,file = file.path(projectPath, "Data",   
       paste0("04_2021216_3Samples_ObjectOriginal_AfterDF_integratedObjDetails_",method.use,"_Annotationed.Rdata")))
  
objIntegrated

####################
# 加载所需的库
library(Seurat)

# 假设 objIntegrated 是您的 Seurat 对象
# 1. 提取 Astrocyte subset
objSubset <- subset(objIntegrated, idents = "Astrocyte")
objSubset <- NormalizeData(objSubset)
objSubset <- FindVariableFeatures(objSubset)
objSubset <- ScaleData(objSubset)
objSubset <- RunPCA(objSubset)
objSubset <- FindNeighbors(objSubset, dims = 1:10)

objSubset <- FindClusters(objSubset, resolution = 0.05)
objSubset <- RunUMAP(objSubset, dims = 1:10)
DimPlot(objSubset, reduction = "umap", label = TRUE)


astrocyte_subtypes <- list(
  # 胶质纤维形成星形胶质细胞（GFAP+ Astrocytes）
  GFAP_Astro = unique(c("S100b", "Aqp4")),
  # 层状星形胶质细胞
  Laminar_Astro = unique(c("Gs", "Slc1a3")),
  # 少突胶质细胞样星形胶质细胞（OPC-like Astrocytes）
  OPC_like_Astro = unique(c("Pdgfra", "Olig2", "Sox10")),
  # 新生星形胶质细胞（Neonatal Astrocytes）
  Neonatal_Astro= unique(c("Gfp", "Nes")),
  # 反应性星形胶质细胞
  Reactive_Astro = unique(c("Gfap",  "Il6")),
  # 足突星形胶质细胞
  Perivascular_Astro = unique(c("Connexin43", "Gja1"))
)
# 打印 astrocyte_subtypes 列表

p_dot <- DotPlot(objSubset, 
                 features = astrocyte_subtypes,
                 group.by = "RNA_snn_res.0.05",  # 或用 "cell_type" 如果已添加注释
                 cols = c("lightgrey", "red"),
                 dot.scale = 6) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.position = "right") +
  ggtitle("Cell Type Markers of Astrocyte")
    
  
cluster_annotation <- c(
  "0" = "C0 unknown",
  "1" = "OPC_like_Astro",     #*
  "2" = "Astrocyte",             #*
  "3" = "Reactive_Astro",     #*

objSubset.markers<-FindAllMarkers(objSubset,only.pos=TRUE,min.pct=0.25)
#top10<-objIntegrated.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(objSubset.markers,file.path(projectPath, "Data",   
      "2021218_3Samples_subObjectAstrocyte_ReClassed_FindAllMarkers.csv"),quote=F)
objSubset.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup() -> top10 


DoHeatmap(objSubset, features = top10$gene) + NoLegend()


save(objSubset,file = file.path(projectPath, "Data",   
       "04_2021216_3Samples_subObjectAstrocyte_ReClassed_objSubset.Rdata"))
  










ref <- celldex::MouseRNAseqData()

objIntegrated
objIntegrated@meta.data$seurat_clusters = objIntegrated$RNA_snn_res.0.05
query_data <- GetAssayData(objIntegrated, 
                           layer = "data", # 使用 log-normalized counts
                           assay = "RNA") # 确保使用原始 'RNA' Assay
						   
						 
clusters <- objIntegrated@meta.data$seurat_clusters 
pred_results <- SingleR(sc_data = query_data,
                        ref_data = ref,
                        #labels = ref$label.main, # 使用主要标签进行粗略注释
                        clusters = clusters
                      ) 				  				  
objIntegrated$SingleR_label <- pred_results$labels[match(objIntegrated@meta.data$seurat_clusters, rownames(pred_results))]

print(table(objIntegrated$SingleR_label))


avg <- AverageExpression(objIntegrated, assays = "RNA", slot = "data",
                         group.by = "seurat_clusters", verbose = FALSE)$RNA
pred <- SingleR(sc_data = avg,
                        ref_data = ref_data,
                        labels = ref_data$label.main, # 使用主要标签进行粗略注释
                      ) 
cluster_ids <- as.character(objIntegrated$seurat_clusters)
objIntegrated$SingleR_label <- pred$labels[match(cluster_ids, colnames(avg))]

table(objIntegrated$SingleR_label)


  
  # -------------------- 1. 预设参数与对象 --------------------  
  # 尝试不同的resolution值  
  resolutions <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25)  
  metrics_df <- data.frame(  
    resolution = resolutions,  
    modularity = numeric(length(resolutions)),  
    silhouette = numeric(length(resolutions)),  
    n_clusters = integer(length(resolutions))  
  )
  
  # -------------------- 2. 构建SNN图 --------------------  
  # 构建共享最近邻图（复用后续计算）  
  current_obj <- FindNeighbors(
    current_obj, 
	features = VariableFeatures(current_obj),  
    k.param = 20,  # 根据数据量调整（建议范围15-50）  
    prune.SNN = 1/15,  # 保持相同剪枝参数  
	dims = 1:10,
  #  graph.name = "snn_graph"  # 指定SNN图名称  
   )  
  snn_graph <- current_obj@graphs$snn_graph  # 假设使用默认图名称

  # -------------------- 3. 循环评估不同的resolution --------------------   
  # 循环计算不同 resolution 的聚类和指标  
  for (j in seq_along(resolutions)) {  
    # 3.1 获取当前聚类标签  
    current_obj <- FindClusters(  
      current_obj,   
      resolution = res,
      algorithm = 1	  
    )   
    # 3.2 从meta.data中获取聚类标签
    cluster_col <- paste0("RNA_snn_res.", res)
    clusters <- current_obj[[cluster_col, drop = TRUE]]
	
  # 1. 设置需要测试的一系列分辨率  
    resolutions <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25)  

  # 2. 先基于 Seurat 对象构建 SNN 图（供后续循环复用）  
    current_obj <- FindNeighbors(  
      current_obj,  
      features = VariableFeatures(current_obj),  
      dims = 1:10,  
      k.param = 20,       # 根据数据量大小微调  
      prune.SNN = 1/15,   # 与需求保持一致  
   # graph.name = "snn_graph"  # 如需要自定义图名称，可打开此行  
    )  

    # 3. 初始化一个 data.frame 用于记录各分辨率的指标  
    metrics_df <- data.frame(  
      resolution = resolutions,  
      modularity = numeric(length(resolutions)),  
      silhouette = numeric(length(resolutions)),  
      n_clusters = integer(length(resolutions))  
    )  

    # 4. 循环评估不同 resolution  
    for (ii in seq_along(resolutions)) {  
      res <- resolutions[ii]  
  
      # 4.1 运行 FindClusters  
      # 如果上面 graph.name 指定为 "snn_graph" 记得在这里也指定  
      current_obj <- FindClusters(  
        current_obj,  
        resolution = res,  
        algorithm = 1  # Louvain，若想尝试 Leiden 可改为 4  
        # graph.name = "snn_graph"  
      )  
  
      # 4.2 从 meta.data 中获取该分辨率的聚类标签列  
      cluster_col <- paste0("RNA_snn_res.", res)  # 如果实际上是 snn_graph_res.[res]，请改这里  
      clusters <- current_obj[[cluster_col, drop = TRUE]]
  
      # 4.3 计算模块度 Modularity（基于 igraph）  
      #     先获取 SNN 图（默认名称可能是 'RNA_snn'，若自定义则改为 'snn_graph'）  
      snn_graph <- current_obj@graphs$RNA_snn  
      # 转成 igraph 对象  
      library(igraph)  
      ig <- graph_from_adjacency_matrix(  
        snn_graph,  
        mode = "undirected",  
        weighted = TRUE
      ) 
      metrics_df$modularity[ii] <- modularity(ig, membership = clusters)  

      # 4.4 计算轮廓系数 Silhouette（在 PCA 空间）  
      #     抽样部分细胞加快计算  
      library(cluster)  
      set.seed(42)  
      sampled_cells <- sample(colnames(current_obj), size = min(1000, ncol(current_obj)))  
      pca_mat <- Embeddings(current_obj, reduction = "pca")[sampled_cells, 1:30]  
      sil_result <- cluster::silhouette(as.numeric(clusters[sampled_cells]), dist(pca_mat))  
      metrics_df$silhouette[ii] <- mean(sil_result[, "sil_width"])  
  
      # 4.5 记录聚类数  
      metrics_df$n_clusters[ii] <- length(unique(clusters))  
    }  
  # 在循环中添加 clustree 图  
  clustree_plot <- clustree(current_obj, prefix = "RNA_snn_res.") +  
    ggtitle(paste0(obj_name, " Clustering Tree")) +  
    theme(plot.title = element_text(size = 14, face = "bold"))  
  print(clustree_plot)

  # 5. 根据评估指标选最优分辨率  
  #   这里以模块度和轮廓系数的平均排名为例，也可自定义更复杂的方法  
  metrics_df <- metrics_df %>%  
    mutate(  
      # 标准化分数  
      scaled_mod = scale(modularity),  
      scaled_sil = scale(silhouette),  
      # 简单加权综合分数（可根据需要微调权重）  
      combined_score = 0.5 * scaled_mod + 0.5 * scaled_sil  
    )  

  # 6. 找到 combined_score 最高的分辨率  
  best_res_row <- metrics_df[which.max(metrics_df$combined_score), ]  
  best_resolution <- best_res_row$resolution  
  cat("最佳分辨率为:", best_resolution, "\n")  
  print(best_res_row)

  # -------------------- 6. 应用最优resolution聚类并可视化 --------------------  
  current_obj <- FindClusters(current_obj, resolution = best_resolution, algorithm = 1)
  current_obj <- RunUMAP(current_obj, dims = 1:10, verbose = FALSE)
  umap_plot <- DimPlot(  
    current_obj,  
    reduction = "umap",  
    label = TRUE  
  ) +  
    ggtitle(paste0(obj_name, " (Final resolution = ", best_resolution, ")"))  
  print(umap_plot)  

  # -------------------- 7. 找到Marker基因，存储结果 --------------------  
  cluster_markers <- FindAllMarkers(  
    current_obj,  
    only.pos = TRUE,  
    min.pct = 0.25,  
    logfc.threshold = 0.25  
  )  

  # 在结果中记录双细胞信息
  marker_list[[obj_name]] <- list(
    markers = cluster_markers,
    final_resolution  = best_resolution,
    doublet_info = list(
      pK = pK,
      nExp_poi = nExp_poi,
      removed_doublets = nExp_poi,
      remaining_cells = ncol(current_obj)
	  ) 
	)
  # 更新对象  
  seurat_list[[obj_name]] <- current_obj  
  
  # 保存当前进度  
  save(current_obj,  
       file = file.path(projectPath, "Data",   
                       paste0("20250220_45Samples_CanPBMC_Object_", i, ".Rdata")))  
  } 
}  
dev.off()

 
# 保存最终结果  
save(seurat_list, marker_list,  
     file = file.path(projectPath, "Data", "20250227_45Samples_CanPBMC_Object.Rdata"))  

# 创建结果总结报告  
summary_file <- file.path(projectPath, "Output", "clustering_summary.txt")  
sink(summary_file)  
cat("Clustering Analysis Summary\n")  
cat("=========================\n\n")  
for(obj_name in names(seurat_list)) {  
  cat(sprintf("\nSample: %s\n", obj_name)) 
  # 添加双细胞信息
  di <- marker_list[[obj_name]]$doublet_info
  cat(sprintf("Doublet removal:\n"))
  cat(sprintf("  - pK value: %.3f\n", di$pK))
  cat(sprintf("  - Predicted doublets: %d\n", di$nExp_poi))
  cat(sprintf("  - Remaining cells: %d\n", di$remaining_cells))
  cat(sprintf("  - Doublet removal rate: %.1f%%\n", 
              100*di$nExp_poi/(di$nExp_poi + di$remaining_cells)))
			  
  cat(sprintf("Optimal resolution: %.3f\n", marker_list[[obj_name]]$optimal_resolution))  
  cat(sprintf("Number of clusters: %d\n",   
              length(unique(Idents(seurat_list[[obj_name]])))))  
  cat("Top markers per cluster:\n")  
  top_markers <- marker_list[[obj_name]]$markers %>%  
    group_by(cluster) %>%  
    top_n(5, avg_log2FC)  
  print(top_markers)  
  cat("\n-------------------\n")  
}  
sink()


load(file.path(projectPath, "Data", "20250227_45Samples_CanPBMC_Object.Rdata"))
seurat_list
marker_list
markers = marker_list
i=44
names(seurat_list[i]) 
seurat_list[[i]]$seurat_clusters <- seurat_list[[i]]$RNA_snn_res.0.15
save(seurat_list,marker_list
     file = file.path(projectPath, "Data", "20250303_45Samples_CanPBMC_ObjectList.Rdata"))  

A6: 0.15;
C7: 0.15;
D6: 0.1;
K7: 0.15;
M7: 0.15;
F6: 0.15;
G6: 0.15;
H6: 0.25;
H9: 0.15;
I9: 0.1;
J6: 0.2;
J9: 0.2;
K9: 0.1;
M6: 0.2;
E9: 0.2;
M9: 0.1;
D6−2: 0.15

name = 'D6-2'
marker_list[[name]]['final_resolution'] = 0.15

load(file.path(projectPath, "Data", "20250303_45Samples_CanPBMC_ObjectList.Rdata"))
umap_pdf <- file.path(projectPath, "Output", "All_clusterTreeUMAP_Plots.pdf")  
pdf(umap_pdf, width = 10, height = 8) 
for(n in (1:length(seurat_list))) {
  name = names(seurat_list[n])
  current_obj = seurat_list[[name]]
  clustree_plot <- clustree(current_obj, prefix = "RNA_snn_res.") +  
    ggtitle(paste0(name, " Clustering Tree")) +  
    theme(plot.title = element_text(size = 14, face = "bold"))  
  print(clustree_plot)
  best_resolution = marker_list[[name]]['final_resolution']
  umap_plot <- DimPlot(  
    current_obj,  
    reduction = "umap",  
    label = TRUE  
  ) +  
    ggtitle(paste0(name, " (Final resolution = ", best_resolution, ")"))  
  print(umap_plot)  
}
dev.off()

######################
######-----Azimuth注释
load(file.path(projectPath, "Data", "20250303_45Samples_CanPBMC_ObjectList.Rdata"))
library(Seurat)       # >=4.3
library(SeuratDisk)   # 用于加载 .h5seurat 参考文件
library(Azimuth)      # GetAzimuthReference() (若 Seurat <4.3 则不需单独安装)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
#ref_file <- "/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/Data/pbmc3k_Azimuth/AzimuthReference_HumanPBMC.Rds"
#reference <- readRDS(ref_file)   # 加载为 Seurat 对象
seurat_listAzi = list()
for(obj_name in names(seurat_list)) {  
  obj = seurat_list[[obj_name]]
  scRNA1Azi <- RunAzimuth(obj, reference = "pbmcref")
  seurat_listAzi[[obj_name]] = scRNA1Azi
}  
save(seurat_listAzi,
     file = file.path(projectPath, "Data", "20250303_45Samples_CanPBMC_ObjectList_Azimuth.Rdata"))  

load(file = file.path(projectPath, "Data", "20250303_45Samples_CanPBMC_ObjectList_Azimuth.Rdata"))
#################################################################
######-------- 去死样本单独分析
target_samples <- c(
  "C6-1", "C6", "C6-2",
  "D6-1", "D6", "D6-2",
  "H9-1", "H9", "H9-2",
  'B9,C9,D9',   #去死细胞
  'F9','G9','I9'
)
## ------------------------------------------------------------------
## 2. 在 seurat_listAzi 中筛选
## ------------------------------------------------------------------
available <- intersect(target_samples, names(seurat_listAzi))
missing   <- setdiff(target_samples,  names(seurat_listAzi))

# 提取并按 target_samples 的顺序排序
seurat_sub <- seurat_listAzi[match(available, names(seurat_listAzi))]
integration_info <- list(  
  #list(reduction = "integrated.cca", cluster.name = "cca_clusters",method.use ='CCAIntegration',umap.reduction = 'umap.cca'),  
  #list(reduction = "integrated.rpca", cluster.name = "rpca_clusters",method.use ='RPCAIntegration',umap.reduction = 'umap.rpca'),  
  #list(reduction = "harmony", cluster.name = "harmony_clusters",method.use ='HarmonyIntegration',umap.reduction = 'umap.harmony'),  
  list(reduction = "fastMNN", cluster.name = "fastMNN_clusters",method.use ='FastMNNIntegration',umap.reduction = 'umap.fastMNN')  
  #list(reduction = "integrated.jointPCA", cluster.name = "jointPCA_clusters",method.use = 'JointPCAIntegration',umap.reduction = 'umap.jointPCA')
  #list(reduction = "integrated.scvi", cluster.name = "scVI_clusters",method.use = 'scVIIntegration')    
)
seurat_list = seurat_sub
unintegrated <- merge(seurat_list[[1]], seurat_list[2:9]) 
unintegrated <- NormalizeData(unintegrated)
unintegrated <- FindVariableFeatures(unintegrated)
unintegrated <- ScaleData(unintegrated)
unintegrated <- RunPCA(unintegrated)
unintegrated <- FindNeighbors(unintegrated, dims = 1:30)
unintegrated <- RunUMAP(unintegrated, dims = 1:30,reduction.name = "unintegrated")
save(unintegrated,file = file.path('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu',   
                       paste0("20250613_9CanPBMCobject_integratedObjDetails_","unIntegrated.Rdata")))
					   
load(file = file.path('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu',   
                       paste0("20250613_9CanPBMCobject_integratedObjDetails_","unIntegrated.Rdata")))					   
objIntegratedList = list()
obj = unintegrated
library(future)
## Only the *outer* loop is parallel
plan(multisession, workers = 8)          # outer level
options(future.globals.maxSize = 25 * 1024^3)   # safety margin

for (i in seq_along(integration_info)) {
  res <- future({
    red.use  <- integration_info[[i]]$reduction
    clst.name <- integration_info[[i]]$cluster.name
    method.use <- integration_info[[i]]$method.use
    umap.reduction <- integration_info[[i]]$umap.reduction
    print(method.use)

    ## inner part runs sequentially to avoid huge export
    plan(sequential)
    objIntegrated <- IntegrateLayers(
        object = obj,
        method = method.use,
        orig.reduction = "pca",
        new.reduction = red.use,
        verbose = FALSE
    )
    plan(sequential)
    objIntegrated <- JoinLayers(objIntegrated)
    objIntegrated <- FindNeighbors(objIntegrated, reduction = red.use, dims = 1:30)
    objIntegrated <- RunUMAP(objIntegrated, reduction = red.use, dims = 1:30,
                             reduction.name = umap.reduction)
    save(objIntegrated,
         file = file.path("/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu",
                          sprintf("20250613_9CanPBMCobject_integratedObjDetails_%s_objIntegrated.Rdata",
                                  method.use)))
    objIntegrated
  })
  objIntegratedList[[integration_info[[i]]$method.use]] <- value(res)
}
save(objIntegratedList,file = file.path(projectPath, "Data",   
                       paste0("20250613_9CanPBMCobject_integratedObjDetails_objIntegratedList.Rdata")))
					   
			

######整合后的数据用不同Res进行降维聚类，进行ClusterTreePlot
resolutions <-c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25,0.3,0.35,0.4)  
IntegrationObj = objIntegratedList[['HarmonyIntegration']]
for (i in resolutions){
  IntegrationObj = FindClusters(  
    object     = IntegrationObj,  
    resolution = i,  
    algorithm  = 1 
    #cluster.name = clst.name,
  )
}
library(clustree)    # 提供 clustree()
library(ggplot2)  
umap_pdf <- file.path(projectPath, "IntegrationObj_Harmony_AllClusterTreeUMAP_Plots_.pdf")  
pdf(umap_pdf, width = 10, height = 8) 
#for(n in (1:length(seurat_list))) {
  name = 'Harmony' #names(seurat_list[n])
  current_obj =IntegrationObj
  clustree_plot <- clustree(current_obj, prefix = "RNA_snn_res.") +  
    ggtitle(paste0(name, " Clustering Tree")) +  
    theme(plot.title = element_text(size = 14, face = "bold"))  
  print(clustree_plot)
  best_resolution = 'RNA_snn_res.0.2'
  Idents(current_obj) = current_obj$RNA_snn_res.0.2
  umap_plot <- DimPlot(  
    current_obj,  
    reduction = "umap.harmony",  
    label = TRUE  
  ) +  
    ggtitle(paste0(name, " (Final resolution = ", best_resolution, ")"))  
  print(umap_plot) 

  umap_plot01 <- DimPlot(  
    current_obj,  
    reduction = "umap.harmony",
    group.by = 'predicted.celltype.l1',	
    label = TRUE  
  )
  print(umap_plot01) 
  
    umap_plot001 <- DimPlot(  
    current_obj,  
    reduction = "umap.harmony",
    group.by = 'predicted.celltype.l2',	
    label = TRUE  
  )
  print(umap_plot001)
  
    umap_plot0001 <- DimPlot(  
    current_obj,  
    reduction = "umap.harmony",
    group.by = 'orig.ident',	
    label = FALSE  
  )
  print(umap_plot0001)
  
#}
dev.off()
method.use='HarmonyIntegration'
IntegrationObj = current_obj
objIntegratedList[['HarmonyIntegration']] = IntegrationObj
objIntegrated =current_obj
save(objIntegrated,file = file.path('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu',   
                       paste0("20250613_9CanPBMCobject_integratedObjDetails_",method.use,"_objIntegrated_FindClusters.Rdata")))


######-------- 去死样本单独分析
target_samples <- c(
  "C6-2","D6-2","H9-2",         #FlowCytometryTreat
  'B9','C9','D9','C6','D6','H9',   #KitTreat
  'F9','G9','F6' ,"C6-1","D6-1" ,"H9-1"  # WithoutTreat
)
B=load('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu/20250613_9CanPBMCobject_integratedObjDetails_HarmonyIntegration_objIntegrated_FindClusters.Rdata')
Idents(objIntegrated) = 'RNA_snn_res.0.2'
subObjIntegrated = subset(objIntegrated,idents = '7')

A=load('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/Data/20250402_45CanPBMCobject_integratedObjDetails_HarmonyIntegration02.Rdata')
Idents(seurat_obj) = 'orig.ident'
subObject <- subset(seurat_obj, idents = target_samples)

# 创建样本到处理类型的映射
sample_treatment_mapping <- c(
  "C6-2" = "FlowCytometryTreat",
  "D6-2" = "FlowCytometryTreat", 
  "H9-2" = "FlowCytometryTreat",
  "B9" = "KitTreat",
  "C9" = "KitTreat",
  "D9" = "KitTreat",
  "C6" = "KitTreat",
  "D6" = "KitTreat",
  "H9" = "KitTreat",
  "F9" = "WithoutTreat",
  "G9" = "WithoutTreat",
  "F6" = "WithoutTreat",
  "C6-1" = "WithoutTreat",
  "D6-1" = "WithoutTreat",
  "H9-1" = "WithoutTreat"
)

# 添加新的 metadata 列
subObject$Treatment <- data.frame(sample_treatment_mapping[subObject$orig.ident])[,1]















###  Seurat转换格式成Scanpy				   
library(SeuratDisk)  ###SaveH5Seurat-Convert ,but somewhere wrong 
#SaveH5Seurat(unintegrated, filename=file.path(projectPath, "Data", "20250526_3SamplesMSC_ObjectOriginal_AfterDoubletsScaleCellCircle_unintegrated.h5Seurat"))   报错
#Convert(file.path(projectPath, "Data", paste0("20250526_3SamplesMSC_ObjectOriginal_AfterDoubletsScaleCellCircle_unintegrated.h5Seurat")),                               报错
#        dest="h5ad", overwrite=TRUE) 

load(file.path('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu',   
                       paste0("20250613_9CanPBMCobject_ALL.IntegratioMethods_objIntegratedList.Rdata")))

Usedreductions <- data.frame(vapply(objIntegratedList, function(obj) tail(Seurat::Reductions(obj), 1), character(1)))
projectPath = '/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu'
library(Matrix)
library(Matrix.utils)
samples = names(objIntegratedList)
for (ss in (1:length(names(objIntegratedList)))){
  s = names(objIntegratedList)[ss]
  redutionUsed = Usedreductions[ss,1]
  seu = objIntegratedList[[s]]
  seu$barcode <- colnames(seu)
  #seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
  #seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
  write.csv(seu@meta.data, file=file.path(projectPath, paste0("20250613_9CanPBMCobject_integratedObjDetails_",s,"_metadata.csv")), 
      quote=F, row.names=F)
  # write expression counts matrix
  seu = JoinLayers(seu)
  counts_matrix <- GetAssayData(seu[['RNA']], layer='counts')   # slot = counts, data, scale.data
  writeMM(counts_matrix, file=file.path(projectPath, paste0("20250613_9CanPBMCobject_integratedObjDetails_",s,"_counts.mtx")))
  # write gene names
  write.csv(data.frame('gene'=rownames(counts_matrix)),
             file=file.path(projectPath,  paste0("20250613_9CanPBMCobject_integratedObjDetails_",s,"_geneNames.csv")),
	         quote=F,row.names=F,col.names=F)
  # write dimesnionality reduction matrix, in this example case pca matrix
  #write.csv(seu@reductions$pca@cell.embeddings, 
  #          file=file.path(projectPath, paste0("20250613_9CanPBMCobject_integratedObjDetails_",s,"_pca.csv")), 
  #	   	    quote=F, row.names=F)
  write.csv( Embeddings(seu, reduction = redutionUsed), 
            file=file.path(projectPath, paste0( "20250613_9CanPBMCobject_integratedObjDetails_",s,"_UMAP.csv")), 
	    	quote=F, row.names=F)
} 
 
###Python: Scanpy-AnnData, save as h5ad
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
projectPath = '/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu'


for ss in 0:len(samples):
	s = samples[ss]
	redusd = redutionUsed[ss]
	# load sparse matrix:
	X = io.mmread(''.join([pathTemp, s,"_counts.mtx"]))
	# create anndata object
	adata = anndata.AnnData(X=X.transpose().tocsr())
	# load cell metadata:
	cell_meta = pd.read_csv(''.join([pathTemp, s,"_metadata.csv"]))
	# load gene names:
	with open(''.join([pathTemp, s,"._geneNames.csv"]), 'r') as f:
		gene_names = f.read().splitlines()
	# set anndata observations and index obs by barcodes, var by gene names
	adata.obs = cell_meta
	adata.obs.index = adata.obs['barcode']
	adata.var.index = gene_names
	# load dimensional reduction:
	pca = pd.read_csv(''.join([pathTemp, s,"_UMAP.csv"]))
	pca.index = adata.obs.index
	# set pca and umap
	adata.obsm[redusd] = pca.to_numpy()
	#adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
	# plot a UMAP colored by sampleID to test:
	#sc.pl.umap(adata, color=['seurat_clusters'], frameon=False, save=True)
	# save dataset as anndata format
	adata.write(''.join([pathTemp, s,".h5ad"]))
	# reload dataset
	#adata = sc.read_h5ad(''.join([pathTemp, s,".h5ad"]))			   
import os
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
from scipy import io      # io.mmread
samples = ['20250613_9CanPBMCobject_integratedObjDetails_CCAIntegration',
'20250613_9CanPBMCobject_integratedObjDetails_RPCAIntegration', 
'20250613_9CanPBMCobject_integratedObjDetails_HarmonyIntegration', 
'20250613_9CanPBMCobject_integratedObjDetails_FastMNNIntegration', 
'20250613_9CanPBMCobject_integratedObjDetails_JointPCAIntegration',
'20250613_9CanPBMCobject_integratedObjDetails_Unintegrated']
path_temp = '/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu'
redutionUsed = ['CCA','RPCA','Harmony','FastMNN','JointPCA','Unintegrated']

for sample, red_key in zip(samples, redutionUsed):
    # ---------- 1. 读表达矩阵 (cells × genes) ----------
    mtx_path = os.path.join(path_temp, f"{sample}_counts.mtx")
    X = io.mmread(mtx_path).tocsr().T      # 转置后就是 cells × genes
    # ---------- 2. 读基因名 ----------
    gene_path = os.path.join(path_temp, f"{sample}_geneNames.csv")
    gene_names = pd.read_csv(gene_path, header=None,skiprows=1).iloc[:, 0].tolist()
    # ---------- 3. 读细胞元数据 ----------
    meta_path = os.path.join(path_temp, f"{sample}_metadata.csv")
    cell_meta = pd.read_csv(meta_path)
    # 核对维度
    if X.shape[0] != cell_meta.shape[0]:
        raise ValueError(f"{sample}: cell number mismatch "
                         f"({X.shape[0]} in mtx vs {cell_meta.shape[0]} in metadata)")
    if X.shape[1] != len(gene_names):
        raise ValueError(f"{sample}: gene number mismatch "
                         f"({X.shape[1]} in mtx vs {len(gene_names)} in gene list)")
    # ---------- 4. 创建 AnnData ----------
    adata = ad.AnnData(
        X,
        obs = cell_meta.set_index("barcode"),
        var = pd.DataFrame(index = gene_names)
    )
    # ---------- 5. 读 UMAP / PCA 坐标 ----------
    emb_path = os.path.join(path_temp, f"{sample}_UMAP.csv")
    emb = (
        pd.read_csv(emb_path)
          .set_index("barcode")            # 假设第一列就是 barcode
          .loc[adata.obs_names]            # 顺序对齐 AnnData
          .to_numpy()
    )
    adata.obsm[red_key] = emb             # 例如 red_key = "X_umap"
    # ---------- 6. 保存 .h5ad ----------
    out_path = os.path.join(path_temp, f"{sample}.h5ad")
    adata.write(out_path)
    print(f"✓  {sample} saved to {out_path}")					   
	
	

##############################################
###########   占比变化
##############################################

projectPath = '/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu'

seurat_obj = objIntegrated
type = 'RNA_snn_res.0.2'  #'predicted.celltype.l2'#'seurat_clusters'  #SingleR_Labels
Idents(seurat_obj) = seurat_obj@meta.data[[type]]
#  生成对应向量
treat_vec <- ifelse(
  grepl("-1$", seurat_obj$orig.ident),          # 末尾 -1
  "WithoutTreat",
  ifelse(
    grepl("-2$", seurat_obj$orig.ident),        # 末尾 -2
    "FlowCytometryTreat",
    "KitTreat"                                  # 剩余即为 C6 / D6 / H9
  )
)
#  加入 meta.data
seurat_obj <- AddMetaData(
  object   = seurat_obj,
  metadata = treat_vec,
  col.name = "TreatmentDeadCelles"
)
sampleInfo = cbind(seurat_obj$orig.ident,seurat_obj$TreatmentDeadCelles)
colnames(sampleInfo) <- c("Samples","ControlTreatment0102")
sampleInfo <- unique(sampleInfo) 
objIntegrated = seurat_obj
save(objIntegrated,
				file='/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu/20250613_9CanPBMCobject_integratedObjDetails_HarmonyIntegration_objIntegrated_FindClusters.Rdata')
objIntegratedList[['HarmonyIntegration']] = objIntegrated
save(objIntegratedList,file = file.path(projectPath, "Data",   
                       paste0("20250613_9CanPBMCobject_integratedObjDetails_objIntegratedList.Rdata")))
	

#################################################################
# 确保加载所需的包  
library(Seurat)  
library(dplyr)  
library(tidyr) 
library(rstatix)  
library(ggplot2)  
library(tibble)      # 新增这行
library(ggplot2)
library(ggpubr)
# 计算细胞簇在不同样本中的比例 
type = 'seurat_clusters'
celltype_prop <- prop.table(table(seurat_obj@meta.data[[type]], seurat_obj$orig.ident), margin = 2)  
# 将宽格式转换为长格式  
celltype_prop_wide <- as.data.frame.matrix(celltype_prop) %>%   
  rownames_to_column("Cluster")  
celltype_prop_long <- celltype_prop_wide %>%   
  pivot_longer(cols = -Cluster, names_to = "Sample", values_to = "Proportion")  
# 合并样本信息  
merged_data <- merge(celltype_prop_long, sampleInfo, by.x = "Sample", by.y = "Samples")  
# 创建统一分组并转化为因子  
#merged_data$Group <- factor(ifelse(merged_data$ControlTreatment == "Control", "Control", "Treatment"))  
merged_data$Group <- factor(merged_data$ControlTreatment0102)
#merged_data$Group <- factor(merged_data$ControlTreatment)

# 统计检验  
mw_results <- merged_data %>%  
  group_by(Cluster) %>%  
  wilcox_test(Proportion ~ Group) %>%  
  adjust_pvalue(method = "BH") %>%  
  add_significance()  
# 可视化  
print(head(merged_data))  # 打印检查 merged_data 的结构和数据  
print(mw_results)          # 打印检查 mw_results 的结构和数据  
###  Control   Treatment(Treatment01+Treatment02)
# 创建带有显著性标注的箱线图
p = ggplot(merged_data, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Cluster, scales = "free_y", nrow = 2) +
  # 添加显著性标注
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Control", "Treatment")),
    label = "p.signif",
    method.args = list(alternative = "two.sided")
  ) +
  scale_fill_manual(values = c("Control" = "grey50", "Treatment" = "#E69F00")) +
  theme_bw() +
  labs(x = "Treatment Group", y = "Proportion",
       title = "Cluster Proportions Across Treatment Groups")
#pdf(file.path(projectPath, "Output", paste0("45CanPBMCobject_SeuratIntegratedMethods_UMAP_",method.use,"_Proportions__",type,".pdf")),width=15)
pdf(file.path(projectPath, "Output", paste0("45CanPBMCobject_SeuratIntegratedMethods_UMAP_",method.use,"_ReCluster2.Proportions_",type,"ControlTreatment.pdf")),width=15)
print (p)
dev.off()

###  Control   Treatment01    Treatment02
# 可视化  
p0102=ggplot(merged_data, aes(x = Group, y = Proportion, fill = Group)) +  
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 0.5) +  
  facet_wrap(~ Cluster, scales = "free_y", nrow = 2) +  
  stat_compare_means(  
    comparisons = list(  
      c("WithoutTreat", "FlowCytometryTreat"),  #c("Control", "Treatment01"),  
      c("WithoutTreat", "KitTreat"),  #c("Control", "Treatment02"),  
      c("FlowCytometryTreat", "KitTreat")#c("Treatment01", "Treatment02")  
    ),  
    method = "wilcox.test",  
    label = "p.signif"  
  ) +  
  scale_fill_manual(values = c("WithoutTreat" = "grey50",       #"Control" = "grey50",   
                               "FlowCytometryTreat" = "#E69F00",#"Treatment01" = "#E69F00",   
                               "KitTreat" = "#56B4E9"           #"Treatment02" = "#56B4E9"
							   )) +  
  theme_bw() +  
  labs(x = "Treatment Group", y = "Proportion",  
       title = "Cluster Proportions Across Treatment Groups")  
pdf(file.path(projectPath, paste0("20250613_9CanPBMCobject_HarmonyIntegration.Proportions_",type,"ControlTreatment0102.pdf")),width=25,height=10)
print (p0102)
dev.off()


########### 
library(SCP)
load('/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/XiaoLu/20250613_9CanPBMCobject_integratedObjDetails_HarmonyIntegration_objIntegrated_FindClusters.Rdata')

objIntegrated = subObject
objIntegrated$TreatmentDeadCelles = subObject$Treatment

Idents(objIntegrated) = objIntegrated$TreatmentDeadCelles
objIntegrated$harmony_clusters<-factor(x=objIntegrated$harmony_clusters,c("0","1","2","3","4","5","6","7","8"))

sub01_objIntegrated = subset(objIntegrated,ident = c('WithoutTreat','KitTreat'))
sub02_objIntegrated = subset(objIntegrated,ident = c('WithoutTreat','FlowCytometryTreat'))
pdf(file.path(projectPath, "IntegrationObj_Harmony_ZhanBi_Plots_.pdf") )
p1=CellStatPlot(objIntegrated, stat.by = "RNA_snn_res.0.2", group.by = "TreatmentDeadCelles", plot_type = "trend")
print(p1)
Idents(objIntegrated) = objIntegrated$orig.ident
sub001_objIntegrated = subset(objIntegrated,ident = c("C6-1","C6","C6-2"))
sub002_objIntegrated = subset(objIntegrated,ident = c("D6-1","D6","D6-2"))
sub003_objIntegrated = subset(objIntegrated,ident = c("H9-1","H9","H9-2"))
p2=CellStatPlot(sub001_objIntegrated, stat.by = "RNA_snn_res.0.2", group.by = "orig.ident", plot_type = "trend")
print(p2)

p3=CellStatPlot(sub002_objIntegrated, stat.by = "RNA_snn_res.0.2", group.by = "orig.ident", plot_type = "trend")
print(p3)

p4=CellStatPlot(sub003_objIntegrated, stat.by = "RNA_snn_res.0.2", group.by = "orig.ident", plot_type = "trend")
print(p4)
dev.off()


45CanPBMCobject_SeuratIntegratedMethods_UMAP_HarmonyIntegration_Proportions_SingleR_LabelsTreat0102

pdf(file.path(projectPath, paste0("45CanPBMCobject_SeuratIntegratedMethods_UMAP_HarmonyIntegration_Proportions_SingleR_LabelsTreat0102.pdf")),width=15)
print (p0102)
dev.off()
