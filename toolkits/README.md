收录单细胞转录组分析中，常用的一些软件和脚本


## h5ad数据转化成Seurat使用的格式

- step 1
```python
import scanpy as sc
import scipy.io as sio

adata = sc.read_h5ad("input.h5ad")
adata.obs.to_csv("cell_metadata.csv")
adata.var.to_csv("gene_metadata.csv")
sio.mmwrite("matrix.mtx", adata.X.T)
```

- step 2
```r
library(Seurat)
library(Matrix)
counts <- readMM("matrix.mtx")
cells <- read.csv("cell_metadata.csv", row.names = 1)
genes <- read.csv("gene_metadata.csv", row.names = 1)

rownames(counts) <- rownames(genes)
colnames(counts) <- rownames(cells)

seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = cells,
    assay = "RNA"
)
```
