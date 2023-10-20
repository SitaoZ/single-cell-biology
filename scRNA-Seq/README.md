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



## 2. 质控、过滤和标准化

## 3. 整合多个单细胞样本和批次矫正

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
