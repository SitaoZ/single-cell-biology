# single-cell RNA-Seq analysis
## prepare 
使用cellranger 获取fastq文件，
10x genomics 测序下机的文件目录为：
![image](https://github.com/SitaoZ/single-cell-biology/assets/29169319/fecc084d-d44b-46d2-988e-257c2dff971a)

需要使用cellranger mkfastq将单细胞测序的数据转化成fastq文件，mkfastq测序封装了bcl2fastq, 可以将bcl(binary base call)转化成fastq.
使用的脚本如下
```bash
module load bcl2fastq2
cellranger mkfastq --id=test \
                   --run=/path/to/the/run/folder \
                   --csv=test.csv \
                   --jobmode=local \
                   --localmem=40 \
                   --localcores=12
# 生成的test.csv文件是三列用，分隔的文件
Lane,Sample,Index
1,test_sample,SI-GA-A3
```
生成fastq后可以将read比对到参考基因组，使用cellranger count
```bash
cellranger count --id=sample345 \
                   --transcriptome=/opt/refdata-cellranger-GRCh38-3.0.0 \
                   --fastqs=/home/test/outs/fastq_path/HAWT7ADXX/test_sample/ \
                   --sample=mysample \
                   --expect-cells=6000
# 生成的文件在--id指定的sample345目录下，下面有一个outs目录，里面有给Seurat处理的filtered_feature_bc_matrix文件夹，下面有三个文件barcodes.tsv.gz, feature.tsv.gz和matrix.mtx.gz
ls -sh outs/filtered_feature_bc_matrix/
# total 90M
# 60K barcodes.tsv.gz  300K features.tsv.gz   90M matrix.mtx.gz

# 三个文件的作用，
# barcodes.tsv.gz包括cell barcode 用于cellranger filter
zcat barcodes.tsv.gz | head -5
# AAACCCAAGCGCCCAT-1
# AAACCCAAGGTTCCGC-1
# AAACCCACAGAGTTGG-1
# AAACCCACAGGTATGG-1
# AAACCCACATAGTCAC-1
# how many cells (barcodes)?
zcat barcodes.tsv.gz | wc -l
# 11769
# The `features.tsv.gz` contains the ENSEMBLE id and gene symbol
zcat features.tsv.gz | head -5
# ENSG00000243485 MIR1302-2HG     Gene Expression
# ENSG00000237613 FAM138A Gene Expression
# ENSG00000186092 OR4F5   Gene Expression
# ENSG00000238009 AL627309.1      Gene Expression
# ENSG00000239945 AL627309.3      Gene Expression

## how many genes?
zcat features.tsv.gz | wc -l
# 33538

# matrix.mtx.gz is a sparse matrix which contains the non-zero counts
zcat matrix.mtx.gz | head -10
# %%MatrixMarket matrix coordinate integer general
# %metadata_json: {"format_version": 2, "software_version": "3.0.0"}
# 33538 11769 24825783
# 33509 1 1
# 33506 1 4
# 33504 1 2
# 33503 1 10
# 33502 1 5
# 33500 1 20
# 33499 1 9
# 大多数entries是0，因此稀疏矩阵有效地减少了存储空间。
# 矩阵解释：33538 x 11769 表示矩阵的维度，包括的全部的非零元素为24825783个
# 33509 1 是行列索引，行表示基因，列表示细胞；1 表示计数
# 33506 1 表示行列索引，4表示计数
```

cellranger 软件运行非常慢，20个CPU处理小鼠的单细胞RNA-Seq数据需要运行几天，因此需要一些工具来替换，常用的工具主要有，
[Alevin](https://salmon.readthedocs.io/en/latest/alevin.html)
[Kallisto/bustools](https://www.kallistobus.tools/)
[Scumi](https://bitbucket.org/jerry00/scumi-dev/src/master/)

### 公共单细胞数据
[scRNAseq bioc package](https://bioconductor.org/packages/devel/data/experiment/html/scRNAseq.html)
[scRNAseqdb](https://bioinfo.uth.edu/scrnaseqdb/)
[Broad single cell portal](https://singlecell.broadinstitute.org/single_cell)


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

```r
# 导入包
library(Seurat)
library(biomaRt)
library(knitr)
library(ggplot2)

# 导入之前的Seurat数据
load(file="pre_sample_corrected.RData")
experiment.aggregate

# Scale数据
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mito", "nFeature_RNA"))

# 进行PCA降维
?RunPCA
experiment.aggregate <- RunPCA(object = experiment.aggregate, npcs=100)

# 可视化PCA结果
VizDimLoadings(experiment.aggregate, dims = 1, ncol = 1) + theme_minimal(base_size = 8) # 可视化PC_1
VizDimLoadings(experiment.aggregate, dims = 2, ncol = 1) + theme_minimal(base_size = 8) # 可视化PC_2
DimPlot(object = experiment.aggregate, reduction = "pca") # 两个主成分展示
# 绘制聚焦于主成分的热图。细胞和基因都是根据它们的主成分得分来排序的。允许很好地可视化数据集中的异构源。
DimHeatmap(object = experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE) 
DimHeatmap(object = experiment.aggregate, dims = 7:12, cells = 500, balanced = TRUE)

# 选择哪些主成分用于展示
# 为了克服单个基因的技术噪音，Seurat根据其PCA分数将细胞聚类，每个主成分本质上代表了一个元基因，将相关基因集的信息结合起来。因此，在下游确定包含多少个主成分是一项重要的步骤。
# ElbowPlot绘制主成分的标准差（或如果运行PCA快速算法，则为近似奇异值），以便在图中容易识别出拐点。这个拐点通常很好地对应于显著的主成分，并且运行速度更快。这是选择主成分的传统方法。
ElbowPlot(experiment.aggregate, ndims = 100)

# JackStraw函数随机对数据的一个子集进行排列，并计算这些“随机”基因的投影PCA分数，然后将“随机”基因的PCA分数与观察到的PCA分数进行比较，以确定统计学上的显著性。最终结果是每个基因与每个主成分的关联的p值。我们将具有低p值基因的主成分视为显著的主成分。
experiment.aggregate <- JackStraw(object = experiment.aggregate, dims = 100)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:100)
JackStrawPlot(object = experiment.aggregate, dims = 1:100) + theme(legend.position="bottom")

# 保存文件，进行下游分析
save(experiment.aggregate, file="pca_sample_corrected.RData")
```

## 5. 聚类 Clustering
```r
# 导入包
library(Seurat)
library(ggplot2)
library(dplyr)

# 导入Seurat对象
load(file="pca_sample_corrected.RData")
experiment.aggregate

# 根据检验，选择前25
use.pcs = 1:25

# Seurat实现了一种基于图的聚类方法。细胞之间的距离是基于先前确定的主成分进行计算的。
# 在V4版本中，识别k最近邻的默认方法已更改为annoy（“Approximate Nearest Neighbors Oh Yeah!”）。这是一种广泛用于高维分析的近似最近邻方法，包括单细胞分析在内的许多领域都在使用该方法。大量的社区基准测试表# 明，annoy极大地提高了邻居发现的速度和内存需求，对下游结果的影响可以忽略不计。
# Seurat先前的方法受到了近期应用图聚类方法于scRNAseq数据的著作的启发。简单来说，Seurat通过共享最近邻（SNN）模块化优化聚类算法来识别细胞簇。首先计算k最近邻（KNN）并构建SNN图，然后优化模块化函数以确定聚类。关于算法的完整描述，请参阅Waltman和van Eck（2013年）的文章《The European Physical Journal B》。您可以使用nn.method=”rann”将设置切换回之前的默认设置。
# FindClusters函数实现了基于邻居的聚类过程，其中包含一个分辨率参数，用于设置下游聚类的细粒度，增加该值会导致更多的聚类。我倾向于进行一系列的分辨率测试，然后进行调查和选择。

?FindNeighbors
experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction="pca", dims = use.pcs)

experiment.aggregate <- FindClusters(
    object = experiment.aggregate,
    resolution = seq(0.25,4,0.5),
    verbose = FALSE
)

head(experiment.aggregate[[]])

sapply(grep("res",colnames(experiment.aggregate@meta.data),value = TRUE),
       function(x) length(unique(experiment.aggregate@meta.data[,x])))


# Plot TSNE coloring for each resolution
experiment.aggregate <- RunTSNE(
  object = experiment.aggregate,
  reduction.use = "pca",
  dims = use.pcs,
  do.fast = TRUE)

DimPlot(object = experiment.aggregate, group.by=grep("res",colnames(experiment.aggregate@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "tsne", label = T)
DimPlot(object = experiment.aggregate, group.by=grep("res",colnames(experiment.aggregate@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=3.0, reduction = "tsne", label = T)

# Choosing a resolution
# 选择一个分辨率的作图
Idents(experiment.aggregate) <- "RNA_snn_res.1.25"
table(Idents(experiment.aggregate),experiment.aggregate$orig.ident)

DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "tsne", label = T)

# uMAP降维图
experiment.aggregate <- RunUMAP(
  object = experiment.aggregate,
  dims = use.pcs)
DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "umap", label = T)

DimPlot(object = experiment.aggregate, pt.size=0.5, group.by = "Phase", reduction = "umap")


# 是否可以使用特征图来绘制读值元数据，如nUMI、feature count和百分比Mito
FeaturePlot(experiment.aggregate, features = c('nCount_RNA'), pt.size=0.5)
FeaturePlot(experiment.aggregate, features = c('nFeature_RNA'), pt.size=0.5)
FeaturePlot(experiment.aggregate, features = c('percent.mito'), pt.size=0.5)

# 在默认的“Ident”(当前为“RNA_snn_res.1.25”)中建立与每个组的“平均”细胞相关的系统发育树。树的估计基于在基因表达空间或主成分分析空间中构建的距离矩阵。
Idents(experiment.aggregate) <- "RNA_snn_res.1.25"
experiment.aggregate <- BuildClusterTree(
  experiment.aggregate, dims = use.pcs)

PlotClusterTree(experiment.aggregate)

# 创建新的数据
DimPlot(object = experiment.aggregate, pt.size=0.5, label = TRUE, reduction = "umap")
DimPlot(experiment.aggregate, pt.size = 0.5, label = TRUE, reduction = "tsne")
experiment.merged = experiment.aggregate
Idents(experiment.merged) <- "RNA_snn_res.1.25"

experiment.merged <- RenameIdents(
  object = experiment.merged,
  '21' = '0',
  '25' = '3'
)

table(Idents(experiment.merged))
DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "umap")
DimPlot(experiment.merged, pt.size = 0.5, label = TRUE, reduction = "tsne" )
VlnPlot(object = experiment.merged, features = "percent.mito", pt.size = 0.05)

# 记录cluster

experiment.examples <- experiment.merged
levels(experiment.examples@active.ident)
experiment.examples@active.ident <- relevel(experiment.examples@active.ident, "12")
levels(experiment.examples@active.ident)
# now cluster 12 is the "first" factor

DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "umap")

VlnPlot(object = experiment.examples, features = "percent.mito", pt.size = 0.05)

# relevel all the factors to the order I want
neworder <- sample(levels(experiment.examples), replace=FALSE)
Idents(experiment.examples) <- factor(experiment.examples@active.ident, levels=neworder)
levels(experiment.examples@active.ident)

DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "umap")
VlnPlot(object = experiment.examples, features = "percent.mito", pt.size = 0.05)

newIdent = as.character(Idents(experiment.examples))
newIdent[newIdent == '0'] = paste0("R",as.character(experiment.examples$RNA_snn_res.3.75[newIdent == '0']))

Idents(experiment.examples) <- as.factor(newIdent)
table(Idents(experiment.examples))

DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "umap")
DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)

## Pretty umap using alpha
alpha.use <- 2/5
p <- DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
p$layers[[1]]$mapping$alpha <- alpha.use
p + scale_alpha_continuous(range = alpha.use, guide = F)

# create a new tmp object with those removed
experiment.aggregate.tmp <- experiment.aggregate[,-which(Idents(experiment.aggregate) %in% c("23"))]

dim(experiment.aggregate)

dim(experiment.aggregate.tmp)
DimPlot(object = experiment.aggregate.tmp, pt.size=0.5, reduction = "umap", label = T)

# 鉴定marker基因
?FindMarkers
markers = FindMarkers(experiment.aggregate, ident.1=c(4,13), ident.2 = c(6,7))

head(markers)
dim(markers)
table(markers$avg_log2FC > 0)
table(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0)

VlnPlot(object = experiment.aggregate, features = rownames(markers)[1:2], pt.size = 0.05)
FeaturePlot(
    experiment.aggregate,
    features = c("KCNMA1", "LEFTY1"),
    cols = c("lightgrey", "blue"),
    ncol = 2
)

markers_all <- FindAllMarkers(
    object = experiment.merged,
    only.pos = TRUE,
    min.pct = 0.25,
    thresh.use = 0.25
)
dim(markers_all)

head(markers_all)

table(table(markers_all$gene))

markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)

table(table(markers_all_single$gene))

table(markers_all_single$cluster)

head(markers_all_single)
top10 <- markers_all_single %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(
    object = experiment.merged,
    features = top10$gene
)

# Get expression of genes for cells in and out of each cluster
getGeneClusterMeans <- function(gene, cluster){
  x <- GetAssayData(experiment.merged)[gene,]
  m <- tapply(x, ifelse(Idents(experiment.merged) == cluster, 1, 0), mean)
  mean.in.cluster <- m[2]
  mean.out.of.cluster <- m[1]
  return(list(mean.in.cluster = mean.in.cluster, mean.out.of.cluster = mean.out.of.cluster))
}

## for sake of time only using first six (head)
means <- mapply(getGeneClusterMeans, head(markers_all[,"gene"]), head(markers_all[,"cluster"]))
means <- matrix(unlist(means), ncol = 2, byrow = T)

colnames(means) <- c("mean.in.cluster", "mean.out.of.cluster")
rownames(means) <- head(markers_all[,"gene"])
markers_all2 <- cbind(head(markers_all), means)
head(markers_all2)

experiment.clusters <- experiment.aggregate
experiment.clusters <- RenameIdents(
  object = experiment.clusters,
  '0' = 'cell_type_A',
  '1' = 'cell_type_B',
  '2' = 'cell_type_C'
)
# and so on

DimPlot(object = experiment.clusters, pt.size=0.5, label = T, reduction = "tsne")

experiment.merged$finalcluster <- Idents(experiment.merged)
head(experiment.merged[[]])

table(experiment.merged$finalcluster, experiment.merged$orig.ident)
experiment.sample1 <- subset(experiment.merged, orig.ident == "A001-C-007")

DimPlot(object = experiment.sample1, group.by = "RNA_snn_res.0.25", pt.size=0.5, label = TRUE, reduction = "tsne")

experiment.merged$samplecluster = paste(experiment.merged$orig.ident,experiment.merged$finalcluster,sep = '_')

# set the identity to the new variable
Idents(experiment.merged) <- "samplecluster"

markers.comp <- FindMarkers(experiment.merged, ident.1 = c("A001-C-007_12", "A001-C-104_12"), ident.2= "B001-A-301_12")

head(markers.comp)

experiment.subset <- subset(experiment.merged, samplecluster %in%  c( "A001-C-007_12", "A001-C-104_12", "B001-A-301_12" ))
DoHeatmap(object = experiment.subset, features = head(rownames(markers.comp),20))
Idents(experiment.merged) <- "finalcluster"

save(list=grep('experiment', ls(), value = TRUE), file="clusters_seurat_object.RData")
```

## 6. 富集分析、差异表达和细胞类型鉴定

```r
# 导入包
library(Seurat)
library(ggplot2)
library(limma)
library(topGO)

# 导入Seurat对象
load("clusters_seurat_object.RData")
experiment.merged

Idents(experiment.merged) <- "finalcluster"

# 1 在一个cluster中进行GO分析
cluster12 <- subset(experiment.merged, idents = '12')
expr <- as.matrix(GetAssayData(cluster12))

# 在该cluster中去除基因表达量为0基因
bad <- which(rowSums(expr) == 0)
expr <- expr[-bad,]

# 选择在至少一半的细胞中表达量大于0的基因
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# 定义geneList，如果表达了就是1，不表达就是0
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# 创建 topGOdata对象
GOdata <- new("topGOdata",
                ontology = "BP", # use biological process ontology
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

# 使用Fisher测验检验富集分析的结果
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
        GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)

# 2. 在limma包中使用基于建模的差异表达分析

# 选择至少在10%的细胞中表达的基因
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# 设计矩阵用于统计模型分析
cluster12$proper.ident <- make.names(cluster12$orig.ident)
mm <- model.matrix(~0 + proper.ident + S.Score + G2M.Score + percent.mito + nFeature_RNA, data = cluster12[[]])
head(mm)

tail(mm)

# Fit model in limma
fit <- lmFit(expr2, mm)
head(coef(fit))

# 检验'B001-A-301' - 'A001-C-007'两个样本
contr <- makeContrasts(proper.identB001.A.301 - proper.identA001.C.007, levels = colnames(coef(fit)))
contr

fit2 <- contrasts.fit(fit, contrasts = contr)
fit2 <- eBayes(fit2)
out <- topTable(fit2, n = Inf, sort.by = "P")
head(out, 30)

# 使用scMRMA鉴定细胞类型

# 首先安装scMRMA
# install.packages("devtools")
# devtools::install_github("JiaLiVUMC/scMRMA")

suppressPackageStartupMessages(library(scMRMA))
result <- scMRMA(input = experiment.merged,
                 species = "Hs",
                 db = "panglaodb")

table(result$uniformR$annotationResult)

# 将细胞类型注释的结果添加到metadata
experiment.merged <- AddMetaData(experiment.merged, result$uniformR$annotationResult, col.name = "CellType")
table(experiment.merged$CellType, experiment.merged$orig.ident)

table(experiment.merged$CellType, experiment.merged$finalcluster)

DimPlot(experiment.merged, group.by = "CellType", label = TRUE)

```

## 7. Doublet去除
在单细胞RNA测序实验中，双联体(doublets)是由两个细胞生成的人工文库。它们通常是由于细胞分选或捕获中的错误而产生的，特别是在涉及数千个细胞的基于液滴的方案中。
使用DoubletFinder包去除。
```r
# 双联体去除DoubletFinder

# 首先需要安装 DoubletFinder
if (!any(rownames(installed.packages()) == "DoubletFinder")){
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
}

#  导入包
library(DoubletFinder)
library(Seurat)
library(kableExtra)
library(ggplot2)

experiment_name = "Colon Cancer"
dataset_loc <- "./expression_data_cellranger"
ids <- c("A001-C-007", "A001-C-104", "B001-A-301")

d10x.data <- lapply(ids[1], function(i){
  d10x <- Read10X_h5(file.path(dataset_loc, i, "outs","raw_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
names(d10x.data) <- ids[1]

str(d10x.data)

experiment.data <- CreateSeuratObject(
  d10x.data[[1]],
  project = "A001-C-007",
  min.cells = 0,
  min.features = 200,
  names.field = 2,
  names.delim = "\\-")



experiment.data$percent.mito <- PercentageFeatureSet(experiment.data, pattern = "^MT-")
summary(experiment.data$percent.mito)

VlnPlot(
  experiment.data,
  features = c("nFeature_RNA", "nCount_RNA","percent.mito"),
  ncol = 1, pt.size = 0.3)


RidgePlot(experiment.data, features=c("nFeature_RNA","nCount_RNA", "percent.mito"), log=T, ncol = 2)


table(experiment.data$orig.ident)

experiment.data <- subset(experiment.data, percent.mito <= 8)

experiment.data <- subset(experiment.data, nFeature_RNA >= 400 & nFeature_RNA <= 4000)

experiment.data <- subset(experiment.data, nCount_RNA >= 500 & nCount_RNA <= 12000)

experiment.data

table(experiment.data$orig.ident)
RidgePlot(experiment.data, features=c("nFeature_RNA","nCount_RNA", "percent.mito"), log=T, ncol = 2)


experiment.data <- NormalizeData(experiment.data)
experiment.data <- FindVariableFeatures(experiment.data, selection.method = "vst", nfeatures = 2000)
experiment.data <- ScaleData(experiment.data)
experiment.data <- RunPCA(experiment.data)
experiment.data <- FindNeighbors(experiment.data, reduction="pca", dims = 1:20)
experiment.data <- FindClusters(
    object = experiment.data,
    resolution = seq(0.25,4,0.5),
    verbose = FALSE
)
experiment.data <- RunUMAP(experiment.data, dims=1:20)
DimPlot(object = experiment.data, pt.size=0.5, reduction = "umap", label = T)


sweep.res <- paramSweep_v3(experiment.data, PCs = 1:20, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK.set <- unique(sweep.stats$pK)[2]
nExp_poi <- round(0.08*nrow(experiment.data@meta.data))

experiment.data <- doubletFinder_v3(experiment.data, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(pK.set)), nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

annotations <- experiment.data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(experiment.data@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
experiment.data <- doubletFinder_v3(experiment.data, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(pK.set)), nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_142", sct = FALSE)

experiment.data <- subset(experiment.data,  DF.classifications_0.25_0.02_142 == "Singlet")

sessionInfo()

```


## 8. Monocle 轨迹分析
monocle 能可视化单细胞测序数据，引入在模拟时间(Pseudotime)内对各个单细胞进行排序的策略，
通过利用单细胞在这写过程中的异步性，将他们匹配到相应的生物过程轨迹上，模拟轨迹分析，解析单细胞数据
在信号通路和发育上的作用。
monocle也可以用于聚类，差异表达分析等。

[Monocle](https://lishensuo.github.io/posts/bioinfo/008%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--monocle%E8%BD%A8%E8%BF%B9%E5%88%86%E6%9E%90/)
[TSCAN](https://bioconductor.org/books/3.15/OSCA.advanced/trajectory-analysis.html#overview-5)



## 参考
[single-cell RNAseq integration](https://satijalab.org/seurat/articles/integration_introduction)
[Seurate install instructions](https://satijalab.org/seurat/articles/install.html)
[monocle introduction](https://cole-trapnell-lab.github.io/monocle-release/)
