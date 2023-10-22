# single-cell RNA-Seq analysis
## 0. 软件安装和数据下载
单细胞RNA-Seq数据分析主要在Rstudio中运行，需要安装的软件主要包括

```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if (!any(rownames(installed.packages()) == "rmarkdown")){
  BiocManager::install("rmarkdown")
}

if (!any(rownames(installed.packages()) == "tinytex")){
  BiocManager::install("tinytex")
}

if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}

if (!any(rownames(installed.packages()) == "hdf5r")){
  BiocManager::install("hdf5r")
}

if (!any(rownames(installed.packages()) == "knitr")){
  BiocManager::install("knitr")
}

if (!any(rownames(installed.packages()) == "kableExtra")){
  BiocManager::install("kableExtra")
}

if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}

if (!any(rownames(installed.packages()) == "dplyr")){
  BiocManager::install("dplyr")
}

if (!any(rownames(installed.packages()) == "reshape2")){
  BiocManager::install("reshape2")
}

if (!any(rownames(installed.packages()) == "biomaRt")){
  BiocManager::install("biomaRt")
}

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  BiocManager::install("org.Hs.eg.db")
}

if (!any(rownames(installed.packages()) == "limma")){
  BiocManager::install("limma")
}

if (!any(rownames(installed.packages()) == "topGO")){
  BiocManager::install("topGO")
}

if (!any(rownames(installed.packages()) == "scran")){
  BiocManager::install("scran")
}

if (!any(rownames(installed.packages()) == "remotes")){
  utils::install.packages("remotes")
}

if (!any(rownames(installed.packages()) == "DoubletFinder")){
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
}

## All of these should now load without error.

library(rmarkdown)
library(tinytex)
library(Seurat)
library(hdf5r)
library(knitr)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(reshape2)
library(biomaRt)
library(limma)
library(topGO)
library(org.Hs.eg.db)
library(scran)

sessionInfo()
```
注意：seurat的安装可以使用上述方式，也可以使用remotes.
```r
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
```
[seurat install](https://satijalab.org/seurat/articles/install)


数据下载
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-July-Single-Cell-RNA-Seq-Analysis/main/data_analysis/scRNA_Workshop-PART1.Rmd", "scRNA_Workshop-PART1.Rmd")
download.file("https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/feb28v7lew62um4/expression_data_cellranger.zip", "expression_data_cellranger.zip")
system("unzip expression_data_cellranger.zip")

```

## 1. Seurat从CellRanger中加载数据到R
Seurat 是一个用于单细胞数据QC、分析和可视化结果的一个R包。Seurat可以鉴定和解释单细胞转录组异质性的来源，并且整合不同类型单细胞数据。
包的开发者在也提供了几个教程[tutorial](https://satijalab.org/seurat/articles/get_started.html)。

```r
# 导入包
library(Seurat)
library(kableExtra)
library(ggplot2)

# 设置实验目录和数据信息
experiment_name = "Colon Cancer"
dataset_loc <- "./expression_data_cellranger"
ids <- c("A001-C-007", "A001-C-104", "B001-A-301")

# 读取cellranger 样本信息矩阵
# cellranger对每个样本信息都具有相应的矩阵，
# 读取, 这里只有一个样本的数据进行展示，ids[1]；后面循环读取采用lapply 
experiment.metrics <- read.csv(file.path(dataset_loc, ids[1], "outs/metrics_summary.csv"))
# 转化成dataframe
sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))
# 设置行名
rownames(sequencing.metrics) <- gsub("\\.", " ", rownames(sequencing.metrics))
# 设置列名
colnames(sequencing.metrics) <- "All samples"
# 使用kableExtra生成表格数据， 在Rstudio中可视化
sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")

# 循环读取多个文件
d10x.metrics <- lapply(ids, function(i){
  metrics <- read.csv(file.path(dataset_loc,paste0(i,"/outs"),"metrics_summary.csv"), colClasses = "character")
})
experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids
sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))
row.names(sequencing.metrics) <- gsub("\\."," ", rownames(sequencing.metrics))
sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")


# 读取cell ranger 矩阵数据，创建基本的Seurat数据对象
# cellranger 提供了一个函数 aggr能够合并多个样本的数据，转化成一个新的数据
# Seurat 提供了一个函数 Read10X 和 Read10X_h5来读取10X数据文件夹，
# 首先我们读取一个文件夹
# 然后，使用没有标准化的数据创建Seurat对象, CreateSeuratObject. 保留至少检测到200个基因的细胞
# 同时给每个细胞提取样本名，计算和添加线粒体数据百分比，最后保存原始的Seurat的对象。

# 使用 Cell Ranger (hdf5)读取，创建一个Seurat对象
d10x.data <- lapply(ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc, i, "/outs","raw_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="_")
  d10x
})
names(d10x.data) <- ids
str(d10x.data)

# 如果没有安装hdf5，则可以直接读取矩阵文件

# d10x.data <- sapply(ids, function(i){
#   d10x <- Read10X(file.path(dataset_loc, i, "/outs","raw_feature_bc_matrix"))
#   colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split="-"), '[[', 1L), i, sep="_")
#   d10x
# })
# names(d10x.data) <- ids

# 创建barcode rank plot from Cell Ranger

# 该代码只画了Features和UMI
plot_cellranger_cells <- function(ind){
  xbreaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6)
  xlabels = c("1","10","100","1000","10k","100K","1M")
  ybreaks = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
  ylabels = c("1","2","5","10","2","5","100","2","5","1000","2","5","10k","2","5","100K","2","5","1M")

  pl1 <- data.frame(index=seq.int(1,ncol(d10x.data[[ind]])),
                    nCount_RNA = sort(Matrix:::colSums(d10x.data[[ind]])+1,decreasing=T),
                    nFeature_RNA = sort(Matrix:::colSums(d10x.data[[ind]]>0)+1,decreasing=T)) %>%
    ggplot() + 
    scale_color_manual(values=c("red2","blue4"), labels=c("Features", "UMI"), name=NULL) +
    ggtitle(paste("CellRanger filltered cells:",ids[ind],sep=" ")) + xlab("Barcodes") + ylab("counts (UMI or Features)") + 
    scale_x_continuous(trans = 'log2', breaks=xbreaks, labels = xlabels) + 
    scale_y_continuous(trans = 'log2', breaks=ybreaks, labels = ylabels) +
    geom_line(aes(x=index, y=nCount_RNA, color = "UMI"), size=1.75) +
    geom_line(aes(x=index, y=nFeature_RNA, color = "Features"), size=1.25)

  return(pl1)
}

plot_cellranger_cells(1)
plot_cellranger_cells(1)
plot_cellranger_cells(1)

# 该代码画Background、Features和UMI

cr_filtered_cells <- as.numeric(gsub(",","",as.character(unlist(sequencing.metrics["Estimated Number of Cells",]))))

plot_cellranger_cells <- function(ind){
  xbreaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6)
  xlabels = c("1","10","100","1000","10k","100K","1M")
  ybreaks = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
  ylabels = c("1","2","5","10","2","5","100","2","5","1000","2","5","10k","2","5","100K","2","5","1M")

  pl1 <- data.frame(index=seq.int(1,ncol(d10x.data[[ind]])),
                    nCount_RNA = sort(Matrix:::colSums(d10x.data[[ind]])+1,decreasing=T),
                    nFeature_RNA = sort(Matrix:::colSums(d10x.data[[ind]]>0)+1,decreasing=T)) %>%
    ggplot() + 
    scale_color_manual(values=c("grey50","red2","blue4"), labels=c("UMI_Background", "Features", "UMI_Cells"), name=NULL) +
    ggtitle(paste("CellRanger filltered cells:",ids[ind],sep=" ")) + xlab("Barcodes") + ylab("counts (UMI or Features)") + 
    scale_x_continuous(trans = 'log2', breaks=xbreaks, labels = xlabels) + 
    scale_y_continuous(trans = 'log2', breaks=ybreaks, labels = ylabels) +
    geom_line(aes(x=index, y=nCount_RNA, color=index<=cr_filtered_cells[ind] , group=1), size=1.75) +
    geom_line(aes(x=index, y=nFeature_RNA, color="Features", group=1), size=1.25)

  return(pl1)
}

plot_cellranger_cells(1)
plot_cellranger_cells(2)
plot_cellranger_cells(3)

# 创建Seurat对象
# 过滤标准：去除少于10个细胞出现的基因，去除少于200个Feature的基因

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = experiment_name,
  min.cells = 10,
  min.features = 300,
  names.field = 2,
  names.delim = "\\_")

experiment.aggregate
## An object of class Seurat 
## 21005 features across 30902 samples within 1 assay 
## Active assay: RNA (21005 features, 0 variable features)

str(experiment.aggregate)

# 比对上线粒体基因组reads的百分比
# 低质量和死亡的细胞通常会出现线粒体污染
# 使用PercentageFeatureSet函数计算线粒体QC矩阵
# 
experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^MT-")
summary(experiment.aggregate$percent.mito)

# Seurat对象简介
# Seurat对象是每个单细胞分析的中心，它存储数据集的全部信息，包括数据，注释，分析等。
# R slotNames函数可以查看这些对象的slot名称
slotNames(experiment.aggregate)
##  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
##  [6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
## [11] "version"      "commands"     "tools"

head(experiment.aggregate[[]])
##                             orig.ident nCount_RNA nFeature_RNA percent.mito
## AAACCCAAGGTCCCTG_A001-C-007 A001-C-007        361          322    2.7700831
## AAACCCAAGTTACGAA_A001-C-007 A001-C-007        459          399    2.6143791
## AAACCCAAGTTATGGA_A001-C-007 A001-C-007       2076         1547    0.5780347
## AAACCCACAACGCCCA_A001-C-007 A001-C-007        854          687    1.5222482
## AAACCCACAAGTAGTA_A001-C-007 A001-C-007        496          410    2.8225806
## AAACCCACAGAAGTTA_A001-C-007 A001-C-007        540          466    1.6666667 

# 最后，保存原始对象并查看
save(experiment.aggregate,file="original_seurat_object.RData")
```


## 2. 质控、过滤和标准化
```r
# 模块加载
library(Seurat)
library(biomaRt)
library(ggplot2)
library(knitr)
library(kableExtra)

# 加载第一步骤生成的seurat数据
load(file="original_seurat_object.RData")
experiment.aggregate
set.seed(12345) # 设置随机种子

# 基本的QA/QC数据分析
## 每个样本中每个细胞的基因数显示5%分位数 以0.05间隔展示数据
## 5% Quantiles of Genes/Cell by Sample
kable(do.call("cbind", tapply(experiment.aggregate$nFeature_RNA, 
                      Idents(experiment.aggregate),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of Genes/Cell by Sample") %>% kable_styling()
## 显示每个样本每个细胞的UMI数量的5%分位数
kable(do.call("cbind", tapply(experiment.aggregate$nCount_RNA, 
                                      Idents(experiment.aggregate),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of UMI/Cell by Sample") %>% kable_styling()
## 样品线粒体百分比的5%分位数
kable(round(do.call("cbind", tapply(experiment.aggregate$percent.mito, Idents(experiment.aggregate),quantile,probs=seq(0,1,0.05))), digits = 3),
      caption = "5% Quantiles of Percent Mitochondria by Sample") %>% kable_styling()

## 画出上面三个统计的小提琴图

VlnPlot(
  experiment.aggregate,
  features = c("nFeature_RNA", "nCount_RNA","percent.mito"),
  ncol = 1, pt.size = 0.3)

## 画出上面三个统计的ridge图，另一种展示方式
RidgePlot(experiment.aggregate, features=c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol = 2)

## 绘制每个基因所代表的细胞数分布。
plot(sort(Matrix::rowSums(GetAssayData(experiment.aggregate) >= 3), decreasing = TRUE) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")

## 基因图，细胞间基因表达的散点图(按样品着色)，在建议的过滤截止点处绘制水平线和垂直线。
FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mito", shuffle = TRUE) + geom_vline(xintercept = c(600,12000)) + geom_hline(yintercept = 8)
FeatureScatter(experiment.aggregate, feature1 = "nFeature_RNA", feature2 = "percent.mito", shuffle = TRUE) + geom_vline(xintercept = c(400,4000)) + geom_hline(yintercept = 8)
FeatureScatter(experiment.aggregate, "nCount_RNA", "nFeature_RNA", pt.size = 0.5, shuffle = TRUE)  + geom_vline(xintercept = c(600,12000)) + geom_hline(yintercept = c(400, 4000))

# 细胞过滤
## 我们使用上面的信息过滤掉单元格。在这里，我们选择那些线粒体基因的百分比(最大为8%)，唯一的UMI计数低于1000或大于12,000，并且其中包含至少700个特征。
table(experiment.aggregate$orig.ident)

experiment.aggregate <- subset(experiment.aggregate, percent.mito <= 8)

experiment.aggregate <- subset(experiment.aggregate, nCount_RNA >= 500 & nCount_RNA <= 12000)

experiment.aggregate <- subset(experiment.aggregate, nFeature_RNA >= 400 & nFeature_RNA < 4000)

experiment.aggregate

table(experiment.aggregate$orig.ident)
## 过滤之后画图
RidgePlot(experiment.aggregate, features=c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol = 2)

## 也可以过滤点额外的基因
## 在创建基础Seurat对象时，我们确实过滤掉了一些基因，记得保留>= 10细胞中表达的所有基因。在过滤细胞之后，你可能想要更积极地使用基因过滤器。Seurat没有提供这样的功能(我可以找到)，所以下面是一个可以这样做的功能，它过滤需要最小值(对数归一化)的基因，在至少400个细胞中，这里表达为1。

experiment.aggregate
FilterGenes <-
 function (object, min.value=1, min.cells = 0, genes = NULL) {
   genes.use <- rownames(object)
   if (!is.null(genes)) {
     genes.use <- intersect(genes.use, genes)
     object@data <- GetAssayData(object)[genes.use, ]
   } else if (min.cells > 0) {
     num.cells <- Matrix::rowSums(GetAssayData(object) > min.value)
     genes.use <- names(num.cells[which(num.cells >= min.cells)])
     object = object[genes.use, ]
   }
  object <- LogSeuratCommand(object = object)
  return(object)
}

experiment.aggregate.genes <- FilterGenes(object = experiment.aggregate, min.value = 1, min.cells = 400)
experiment.aggregate.genes
rm(experiment.aggregate.genes)

# 下一步就是标准化数据
?NormalizeData
experiment.aggregate <- NormalizeData(
  object = experiment.aggregate,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

## 细胞周期基因

# this code is for human samples only!
s.genes <- (cc.genes$s.genes)
g2m.genes <- (cc.genes$g2m.genes)

experiment.aggregate <- CellCycleScoring(experiment.aggregate,
                                         s.features = s.genes,
                                         g2m.features = g2m.genes,
                                         set.ident = TRUE)
table(experiment.aggregate@meta.data$Phase) %>%
  kable(caption = "Number of Cells in each Cell Cycle Stage", col.names = c("Stage", "Count"), align = "c") %>%
  kable_styling()

table(Idents(experiment.aggregate))

# 鉴定可变基因
## 函数FindVariableFeatures通过对log(方差)和log(均值)的关系进行loess平滑拟合来识别最高变异基因（默认为2000个基因），然后使用这些信息对数据进行标准化，并计算标准化数据的方差。这有助于避免选择仅由于表达水平而显得变异的基因。
?FindVariableFeatures

experiment.aggregate <- FindVariableFeatures(
  object = experiment.aggregate,
  selection.method = "vst")

length(VariableFeatures(experiment.aggregate))

top10 <- head(VariableFeatures(experiment.aggregate), 10)

top10

vfp1 <- VariableFeaturePlot(experiment.aggregate)
vfp1 <- LabelPoints(plot = vfp1, points = top10, repel = TRUE)
vfp1

# FindVariableFeatures并不是设置Seurat对象的“可变特性”的唯一方法。另一种合理的方法是选择一组“最低表达”基因。
dim(experiment.aggregate)

min.value = 2
min.cells = 10
num.cells <- Matrix::rowSums(GetAssayData(experiment.aggregate, slot = "count") > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])
length(genes.use)
VariableFeatures(experiment.aggregate) <- genes.use

# 最后保存过滤后和标准化的数据
save(experiment.aggregate, file="pre_sample_corrected.RData")
```


## 3. 整合多个单细胞样本和批次矫正

```r
# 最近发展的大部分将单细胞数据集集成的方法可分为两类。第一类是“锚点”基于的方法。在这种方法中，第一步是选择一个批次作为“锚点”，并将其他批次转换为“锚点”批次。在这种方法中，包括MNN、iMAP、SCALEX和Seurat的整合。这种方法的优点是可以在相同的实验条件下研究不同批次的细胞，缺点是由于每个批次中包含的细胞类型未知，无法完全结合每个批次的特征。第二种方法是将所有批次的数据转换到低维空间以校正批次效应，例如在Scanorama、Harmony、DESC和BBKNN中实现。这种第二种方法的优点是提取具有生物学意义的潜在特征并减少噪声的影响，但不能用于差异基因表达分析。许多现有的方法在批次数据集具有相同细胞类型的情况下效果良好，但在不同数据集中涉及不同细胞类型的情况下会失败。最近提出了一种新方法，利用连接图和生成对抗网络（GAN）来消除数据集批次之间的非生物噪声。这种新方法已被证明在数据集具有相同细胞类型和数据集可能具有不同细胞类型的情况下都能很好地工作。
# 在本研讨会中，我们将介绍Seurat的集成方法。其基本思想是识别处于匹配生物状态(“锚点”)的跨数据集对细胞，并使用它们来纠正数据集之间的技术差异。我们使用的集成方法已经在Seurat中实现，您可以在其出版物中找到该方法的详细信息。

# 导入包
library(Seurat)
# 导入Seurat对象
load(file="pre_sample_corrected.RData")
experiment.aggregate
# 分成单个样本
experiment.split <- SplitObject(experiment.aggregate, split.by = "ident")

# Normalize and find variable features for each individual sample
# 标准化，为每个样本找到变量特征
?NormalizeData
?FindVariableFeatures
#
experiment.split <- lapply(X = experiment.split, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across samples and find integration anchors
features <- SelectIntegrationFeatures(object.list = experiment.split)
anchors <- FindIntegrationAnchors(object.list = experiment.split, anchor.features = features)

# Perform integration
experiment.integrated <- IntegrateData(anchorset = anchors)

# PCA plot before integration
experiment.test <- NormalizeData(object=experiment.integrated, assay="RNA")
experiment.test <- ScaleData(object=experiment.test, assay="RNA")
experiment.test <- FindVariableFeatures(object=experiment.test, assay="RNA")
experiment.test <- RunPCA(object=experiment.test, assay="RNA")
DimPlot(object = experiment.test, group.by="ident", reduction="pca", shuffle=TRUE)

# PCA plot after integration
experiment.test <- ScaleData(object=experiment.integrated, assay="integrated")
experiment.test <- FindVariableFeatures(object=experiment.test, assay="integrated")
experiment.test <- RunPCA(object=experiment.test, assay="integrated")
DimPlot(object = experiment.test, group.by="ident", reduction="pca", shuffle=TRUE)

# saved the integrated data 
save(experiment.integrated, file="sample_integrated.RData")

```


## 4. 主成分分析

## 5. 聚类

## 6. 富集分析、差异表达和细胞类型鉴定

## 7. 添加双态检测
















## 8. Monocle 轨迹分析
monocle 能可视化单细胞测序数据，引入在模拟时间(Pseudotime)内对各个单细胞进行排序的策略，
通过利用单细胞在这写过程中的异步性，将他们匹配到相应的生物过程轨迹上，模拟轨迹分析，解析单细胞数据
在信号通路和发育上的作用。
monocle也可以用于聚类，差异表达分析等。




## 参考
[single-cell RNAseq integration](https://satijalab.org/seurat/articles/integration_introduction)
[Seurate install instructions](https://satijalab.org/seurat/articles/install.html)
[monocle introduction](https://cole-trapnell-lab.github.io/monocle-release/)
