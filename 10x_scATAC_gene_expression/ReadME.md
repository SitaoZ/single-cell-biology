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

## 二、barcodes 序列下载
barcodes whitelist: 试剂盒中已知的barcodes序列，例如gex_737K-arc-v1.txt含有736,319。
barcodes可以从cellrager下载![10x](https://kb.10xgenomics.com/hc/en-us/articles/4412343032205-Where-can-I-find-the-barcode-whitelist-s-for-Single-Cell-Multiome-ATAC-GEX-product-)，路径为lib/python/cellranger/barcodes和lib/python/atac/barcodes拷贝，也可以从该链接下载。

```bash
# download the GEX expression whitelist
wget -P 10xMultiome/ https://teichlab.github.io/scg_lib_structs/data/gex_737K-arc-v1.txt.gz
gunzip 10xMultiome/gex_737K-arc-v1.txt.gz

# download the atac whitelist
wget -P 10xMultiome/ https://teichlab.github.io/scg_lib_structs/data/atac_737K-arc-v1.txt.gz
gunzip 10xMultiome/atac_737K-arc-v1.txt.gz

# reverse complement the atac whitelist
cat 10xMultiome/atac_737K-arc-v1.txt | \
    rev | tr 'ACGT' 'TGCA' > \
    10xMultiome/atac_737K-arc-v1_rc.txt


```

## 三、表达量数据fastq文件转化成count矩阵
```bash
mkdir -p 10xMultiome/star_outs
# 使用 starsolo 处理 GEX文件
STAR --runThreadN 4 \
     --genomeDir mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix 10xMultiome/star_outs/ \
     --readFilesIn 10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz 10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
     --soloCBwhitelist 10xMultiome/gex_737K-arc-v1.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## 四、ATAC数据处理
```bash
mkdir -p 10xMultiome/chromap_outs
# 使用 chromap软件比对，生成fragment file
chromap -t 4 --preset atac \
        -x mm10/chromap_index/genome.index \
        -r mm10/mm10.fa \
        -1 10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz \
        -2 10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L001_R3_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L002_R3_001.fastq.gz \
        -b 10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz \
        --barcode-whitelist 10xMultiome/atac_737K-arc-v1.txt \
        -o 10xMultiome/chromap_outs/fragments.tsv

# chrom.size
faSize -detailed mm10/mm10.fa | \
    sort -k1,1 > mm10/mm10.chrom.sizes
# 压缩和建立索引
bgzip 10xMultiome/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed 10xMultiome/chromap_outs/fragments.tsv.gz

# 使用 bedClip to remove reads outside the chromosome boundary
# remove reads mapped to the mitochondrial genome (chrM)

zcat 10xMultiome/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin mm10/mm10.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > 10xMultiome/chromap_outs/reads.bed.gz

# /home/zhusitao/anaconda3/envs/py27/bin/macs2
/home/zhusitao/anaconda3/envs/py27/bin/macs2 callpeak -t 10xMultiome/chromap_outs/reads.bed.gz \
               -g mm -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir 10xMultiome/chromap_outs \
               -n aggregate

# format and sort peaks
cut -f 1-4 10xMultiome/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > 10xMultiome/chromap_outs/aggregate_peaks_sorted.bed


# prepare the overlap
bedtools intersect \
    -a 10xMultiome/chromap_outs/aggregate_peaks_sorted.bed \
    -b 10xMultiome/chromap_outs/reads.bed.gz \
    -wo -sorted -g mm10/mm10.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > 10xMultiome/chromap_outs/peak_read_ov.tsv.gz

# create dirctory
mkdir -p 10xMultiome/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 10xMultiome/chromap_outs/aggregate_peaks_sorted.bed > \
    10xMultiome/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat 10xMultiome/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    10xMultiome/chromap_outs/raw_peak_bc_matrix/barcodes.tsv

python generate_csc_mtx.py \
    10xMultiome/chromap_outs/aggregate_peaks_sorted.bed \
    10xMultiome/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    10xMultiome/chromap_outs/peak_read_ov.tsv.gz \
    10xMultiome/chromap_outs/raw_peak_bc_matrix

# trick starsolo to use peaks.bed as features.tsv by creating symlink
ln -s peaks.bed 10xMultiome/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo
STAR --runMode soloCellFiltering \
     10xMultiome/chromap_outs/raw_peak_bc_matrix \
     10xMultiome/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv 10xMultiome/chromap_outs/filtered_peak_bc_matrix/peaks.bed

```
`generate_csc_mtx.py`
```python
# import helper packages
# most entries in the count matrix is 0, so we are going to use a sparse matrix
# since we need to keep updating the sparse matrix, we use lil_matrix from scipy
import sys
import gzip
from scipy.io import mmwrite
from scipy.sparse import lil_matrix

# the unique peak ID is a good renference
# generate a dictionary with peak_id : index_in_the_file
# sys.argv[1] is the 4-column bed file aggregate_peaks_sorted.bed
peak_idx = {}
with open(sys.argv[1]) as fh:
    for i, line in enumerate(fh):
        _, _, _, peak_name = line.strip().split('\t')
        peak_idx[peak_name] = i

# determine and create the dimension of the output matrix
# that is, to calculate the number of peaks and the number of barcodes
# sys.argv[2] is barcodes.tsv
num_peaks = len(peak_idx.keys())
num_cells = len(open(sys.argv[2]).readlines())
mtx = lil_matrix((num_peaks, num_cells), dtype = int)

# read the information from peak_read_ov.tsv.gz
# update the counts into the mtx object
# sys.argv[3] is peak_read_ov.tsv.gz
with gzip.open(sys.argv[3], 'rt') as fh:
    for i, line in enumerate(fh):
        col_idx = i # each column is a cell barcode
        count_info = line.strip().split('\t')[1]
        items = count_info.split(',')
        for pn_count in items:
            pn, count = pn_count.split(':')
            row_idx = peak_idx[pn] # each row is a peak
            mtx[row_idx, col_idx] = int(count)

# get a CSC sparse matrix, which is the same as the 10x Genomics matrix.mtx
mtx = mtx.tocsc()

# sys.argv[4] is the path to the output directory
mmwrite(sys.argv[4] + '/matrix.mtx', mtx, field='integer')
```
