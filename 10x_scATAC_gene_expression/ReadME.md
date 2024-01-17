## 10x Chromium Single Cell Multiome ATAC + Gene Expression
本文介绍10x单细胞的ATAC结合表达量联合分析


## 一、数据下载
公共数据下载
```bash
# 创建目录
mkdir 10xMultiome
# 下载数据
wget -P 10xMultiome -c https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_fastqs.tar
# 解压数据
tar xf 10xMultiome/e18_mouse_brain_fresh_5k_fastqs.tar -C 10xMultiome/

# 解压缩的数据
tree 10xMultiome/e18_mouse_brain_fresh_5k

10xMultiome/e18_mouse_brain_fresh_5k
├── atac
│   ├── e18_mouse_brain_fresh_5k_S1_L001_I1_001.fastq.gz
│   ├── e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz
│   ├── e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz
│   ├── e18_mouse_brain_fresh_5k_S1_L001_R3_001.fastq.gz
│   ├── e18_mouse_brain_fresh_5k_S1_L002_I1_001.fastq.gz
│   ├── e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz
│   ├── e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz
│   └── e18_mouse_brain_fresh_5k_S1_L002_R3_001.fastq.gz
└── gex
    ├── e18_mouse_brain_fresh_5k_S1_L001_I1_001.fastq.gz
    ├── e18_mouse_brain_fresh_5k_S1_L001_I2_001.fastq.gz
    ├── e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz
    ├── e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz
    ├── e18_mouse_brain_fresh_5k_S1_L002_I1_001.fastq.gz
    ├── e18_mouse_brain_fresh_5k_S1_L002_I2_001.fastq.gz
    ├── e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz
    └── e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz
```
