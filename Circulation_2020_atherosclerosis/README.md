

## single-cell mice Apoe-/- Ldlr-/-

## 1.目录与包导入
```r
setwd("C:/Users/zhusitao/Documents/R_workshop/project/KR/Single-Cell")

library(dplyr) 
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratData)
library(stringr)
####################
````

## 2.文件读取

```r
list.files()
# 使用read.table()函数从txt.gz格式的文件中读取数据，并将第一列作为行名
seurat_data1 <- read.table(gzfile("./GSE155513_RAW/GSM4705592_RPS003_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t") # 0 Week  [ZsGreen1+ Ldlr KO 0 week WD]
seurat_data2 <- read.table(gzfile("./GSE155513_RAW/GSM4705594_RPS011_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t") # 8 Week  [ZsGreen1+ Ldlr KO 8 weeks WD]
seurat_data3 <- read.table(gzfile("./GSE155513_RAW/GSM4705596_RPS007_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t") # 16 Week [ZsGreen1+ Ldlr KO 16 weeks WD]
seurat_data4 <- read.table(gzfile("./GSE155513_RAW/GSM4705598_RPS001_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t") # 26 Week   [ZsGreen1+ Ldlr KO 26 weeks WD]
```

## 3.创建seurat对象
```r
# 读取单个样本
seurat_obj1 <- CreateSeuratObject(counts = seurat_data1,
                                  min.features = 200,
                                  min.cells = 3, 
                                  project = "00w")

seurat_obj2 <- CreateSeuratObject(counts = seurat_data2,
                                  min.features = 200,
                                  min.cells = 3, 
                                  project = "08w")
seurat_obj3 <- CreateSeuratObject(counts = seurat_data3,
                                  min.features = 200,
                                  min.cells = 3, 
                                  project = "16w")
seurat_obj4 <- CreateSeuratObject(counts = seurat_data4,
                                  min.features = 200,
                                  min.cells = 3, 
                                  project = "26w")


data <- list(seurat_obj1, seurat_obj2, seurat_obj3, seurat_obj4)

# 合并数据
seurat_obj <- merge(x=data[[1]], y = data[-1], project = "scBC")
head(seurat_obj)
tail(seurat_obj)
dim(seurat_obj)

seurat_obj <- JoinLayers(seurat_obj)
head(seurat_obj)
tail(seurat_obj)
```






## 4.质量控制
```rrr
# 数据质控
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 数据质控后，使用VlnPlot进行可视化
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


dim(seurat_obj)
# 数据筛选，筛选的是细胞，使得细胞数减少。根据具体情况进行阈值确定
# seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

dim(seurat_obj)

# 标准化，使用的默认参数
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000) 

# 识别高度可变特征
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# Scaling the data标准化数据，使得每个基因的表达量的标准差为1，均值为0
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
```

- 不进行批次矫正
```r
# 不进行批次矫正
# Perform linear dimensional reduction 降维
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), reduction.name = "pca")

# PCA结果可视化
DimPlot(seurat_obj, reduction = "pca") + NoLegend()

# seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, reduction = "pca", reduction.name = "umap")
# seurat_obj <- RunTSNE(seurat_obj, dims = 1:10, reduction = 'tsne', reduction.name = "tsne")
# p1 <- DimPlot(seurat_obj, label.size = 2)

# seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:10)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)

# Look at cluster IDs of the first 5 cells
# head(Idents(seurat_obj), 5)
```

## 5.批次矫正
```r
library(harmony)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", reduction.save = "tsne", dims.use = 1:40)
seurat_obj <- RunTSNE(seurat_obj, reduction = "tsne", dims = 1:10)
DimPlot(seurat_obj, reduction = "tsne", label = T)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)

DimPlot(subset(seurat_obj, subset = orig.ident == "00w"), reduction = 'tsne', group.by = "orig.ident")
DimPlot(subset(seurat_obj, subset = orig.ident == "00w"), reduction = 'tsne', group.by = "seurat_clusters")
DimPlot(subset(seurat_obj, subset = orig.ident == "08w"), reduction = 'tsne', group.by = "seurat_clusters")
DimPlot(subset(seurat_obj, subset = orig.ident == "16w"), reduction = 'tsne', group.by = "seurat_clusters")
DimPlot(subset(seurat_obj, subset = orig.ident == "26w"), reduction = 'tsne', group.by = "seurat_clusters")

```

## 6.结果复现
```r
# 重复B图
FeaturePlot(seurat_obj, features = c("Rnf128", "Lyz2", "Fam107b", "Rnf130"))

# 重复E图
FeaturePlot(subset(seurat_obj, subset = orig.ident == "00w"), features = c("Rnf128"), min.cutoff=0, max.cutoff=0)  
FeaturePlot(subset(seurat_obj, subset = orig.ident == "08w"), features = c("Rnf128"))
FeaturePlot(subset(seurat_obj, subset = orig.ident == "16w"), features = c("Rnf128"))
FeaturePlot(subset(seurat_obj, subset = orig.ident == "26w"), features = c("Rnf128"))

FeaturePlot(subset(seurat_obj, subset = orig.ident == "00w"), features = c("Rnf128", "Lyz2", "Fam107b", "Rnf130"))
FeaturePlot(subset(seurat_obj, subset = orig.ident == "08w"), features = c("Rnf128", "Lyz2", "Fam107b", "Rnf130"))
FeaturePlot(subset(seurat_obj, subset = orig.ident == "16w"), features = c("Rnf128", "Lyz2", "Fam107b", "Rnf130"))
FeaturePlot(subset(seurat_obj, subset = orig.ident == "26w"), features = c("Rnf128", "Lyz2", "Fam107b", "Rnf130"))

# 重复D图
VlnPlot(seurat_obj, features = c("Lyz2", "Spp1"), group.by = "orig.ident")
VlnPlot(seurat_obj, features = c("Lyz2", "Rnf128", "Spp1"), group.by = "orig.ident")
VlnPlot(seurat_obj, features = c("Fam107b", "Rnf130"), group.by = "orig.ident")

# Slfn2 Pqlc1
VlnPlot(seurat_obj, features = c("Lpp","Pde3a","Ppp1r12b","Npnt","Pqlc1","Pdgfrb","Slfn2","Efhd2"), group.by = "orig.ident")
# Trmt1 Slc25a12 Ubr4 Ptcd3 Fam107b
VlnPlot(seurat_obj, features = c("Naif1","Cyp3a41b","Trmt1","Taar7a","Slc25a12","Osbp2","Hrsp12","Ubr4","Ptcd3","Fam107b","Pcdhga9","Mtch2","Smyd1","Olfr1163","Ifnk"), group.by = "orig.ident")

# Lfng
VlnPlot(seurat_obj, features = c("Ctgf","Rgs7bp","Dmpk","Cav1","Lfng"), group.by = "orig.ident")
FeaturePlot(seurat_obj, features = c("Naif1","Cyp3a41b","Trmt1","Taar7a","Slc25a12","Osbp2","Hrsp12","Ubr4","Ptcd3","Fam107b","Pcdhga9","Mtch2","Smyd1","Olfr1163","Ifnk"), label = T)
FeaturePlot(seurat_obj, features = c("Lpp","Pde3a","Ppp1r12b","Npnt","Pqlc1","Pdgfrb","Slfn2","Efhd2"), label.size = 2, label = T)

FeaturePlot(seurat_obj, features = c("Ctgf","Rgs7bp","Dmpk","Cav1","Lfng"), label.size = 2, label = T)
```

## 7.寻找marker基因，便于手动注释
```r
# 找marker基因
seurat_obj.markers <- FindAllMarkers(seurat_obj,  only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 40) %>%
  ungroup() -> top40

DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()


# 手动注释，结合Cell Marker

filter(top40, cluster == 0) # Smooth muscle cell (Cnn1, Myh11)
filter(top40, cluster == 1) # Fibroblast         (Lum)
filter(top40, cluster == 2) # Fibroblast         (Serpinf1, Cygb)
filter(top40, cluster == 3) # Macrophage         (Cd68, C1qb, Lyz2)
filter(top40, cluster == 4) # Monocyte           (C1qb)
filter(top40, cluster == 5) # Fibroblast
filter(top40, cluster == 6) # Fibroblast
filter(top40, cluster == 7) # Fibroblast
filter(top40, cluster == 8) # Endothelial cell  (Egfl7, Gja4)
# 0 ,1,2,3,4,5,6,7,8,9
new.cluster.ids <- c("Smooth muscle cell", "Fibroblast", "Fibroblast", "Macrophage", "Monocyte", "Fibroblast", "Fibroblast", "Fibroblast", "Endothelial cell")
names(new.cluster.ids) <- levels(seurat_obj)

seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave('Harmony_manual.pdf', width = 10, height = 6, units = 'in', dpi=300)

FeaturePlot(seurat_obj, features = c("Abca1", "Slc25a12", "Rnf130", "Apoe", "Myh11", "Fam107a", "Fam107b", "Rnf128", "Lyz2"), label = T)
FeaturePlot(seurat_obj, features = c("Naif1","Cyp3a41b","Trmt1","Taar7a","Slc25a12","Osbp2","Hrsp12","Ubr4","Ptcd3","Fam107b","Pcdhga9","Mtch2","Smyd1","Olfr1163","Ifnk"), label = T)
FeaturePlot(seurat_obj, features = c("Lpp","Pde3a","Ppp1r12b","Npnt","Pqlc1","Pdgfrb","Slfn2","Efhd2"), label.size = 2, label = T)
ggsave('select_DM_26_week.pdf', width = 10, height = 6, units = 'in', dpi=300)

FeaturePlot(seurat_obj, features = c("Ctgf","Rgs7bp","Dmpk","Cav1","Lfng"), label.size = 2, label = T)
ggsave('select_UM_26_week.pdf', width = 10, height = 6, units = 'in', dpi=300)
####################################################################################################################################
```

## 8.自动注释
```r
# 注释细胞类型
library(SingleR)
# hpca.se <- celldex::HumanPrimaryCellAtlasData()  # 加载参考数据集
# bpe.sc <- celldex::BlueprintEncodeData()

hpca.se <- celldex::MouseRNAseqData()
bpe.sc <- celldex::ImmGenData()

# 使用 label.fine 代替 label.main：

#annotations <- SingleR(test = seurat_obj@assays$RNA$data, 
#                       ref = list(BP = bpe.sc, HPCA = hpca.se),
#                       labels =list(bpe.sc$label.fine, hpca.se$label.fine),
#                       clusters = seurat_obj@meta.data$seurat_clusters,
#                       genes = "de",
#                       de.method = "wilcox",
#                       de.n = 1000,
#                       quantile = 0.8,
#                       fine.tune = TRUE,
#                       tune.thresh = 0.05,
#                       )

# 只用这个分析的比较好5个类：de.method = "wilcox", de.n = 1000, quantile = 0.8
annotations <- SingleR(test = seurat_obj@assays$RNA$data, 
                       ref = hpca.se,
                       labels = hpca.se$label.fine,
                       clusters = seurat_obj@meta.data$seurat_clusters,
                       genes = "de",
                       de.method = "wilcox",
                       de.n = 1000,
                       quantile = 0.8
)

unique(hpca.se$label.main)

annotations <- SingleR(test = seurat_obj@assays$RNA$data, 
                       ref = bpe.sc,
                       labels = bpe.sc$label.main,
                       clusters = seurat_obj@meta.data$seurat_clusters,
                       genes = "de",
                       de.method = "classic",
                       de.n = 1000,
                       quantile = 0.8
)


table(annotations$labels)

seurat_obj$cell_type <- annotations$labels

plotScoreHeatmap(annotations, show_colnames = F)
# plotDeltaDistribution(annotations, ncol = 3)


seurat_obj@meta.data$SingleR_label <- unlist(lapply(seurat_obj@meta.data$seurat_clusters, function(x){annotations$labels[x]}))


# DimPlot(seurat_obj, reduction = 'umap', group.by = 'SingleR_label', label = T)
DimPlot(seurat_obj, reduction = 'tsne', group.by = 'SingleR_label', label = T)

FeaturePlot(seurat_obj, features = c("Lyz2"), label = T)
ggsave('GSM4705599_26_week.pdf', width = 10, height = 6, units = 'in', dpi=300)
ggsave('GSM4705598.pdf', width = 10, height = 6, units = 'in', dpi=300)

FeaturePlot(seurat_obj, features = c("Abca1", "Hmgcr", "Ldlr", "Apoe", "Myh11", "Fasn", "Gt(ROSA)26Sor", "Rnf128", "Lyz2"), label = T)

FeaturePlot(seurat_obj, features = c("Abca1", "Hmgcr", "Rnf130", "Apoe", "Myh11", "Fam107a", "Slc30a5", "Rnf128", "Lyz2"), label = T)

```
