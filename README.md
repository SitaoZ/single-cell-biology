# Single-Cell-Biology
Mass tools in single cell resolution 

- 10X数据命名规则  
[10X说明](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/fastq-input#gex_rightname)  
```bash
  [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
```
```bash
I1: Sample index read (optional) （P7 index，区分样本，8个碱基）
I2: Sample index read (optional)
R1: Read 1 （16个Barcode碱基，然后是12个UMI碱基）
R2: Read 2 （98个碱基，就是转录本的reads）
```
