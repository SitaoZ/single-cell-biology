# single-cell RNA-Seq analysis
## 1. 软件安装和数据下载
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

## 2. Seurat从CellRanger中加载数据到R
Seurat 是一个用于单细胞数据QC、分析和可视化结果的一个R包。Seurat可以鉴定和解释单细胞转录组异质性的来源，并且整合不同类型单细胞数据。
包的开发者在也提供了几个教程[tutorial](https://satijalab.org/seurat/articles/get_started.html)。


## 3. 质控、过滤和标准化

## 4. 整合多个单细胞样本和批次矫正

## 5. 主成分分析

## 6. 聚类

## 7. 富集分析、差异表达和细胞类型鉴定

## 8. 添加双态检测


