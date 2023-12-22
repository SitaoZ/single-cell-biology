## STAR 
STAR 软件可以用来处理10xgenomics的单细胞测序数据，可以对数据进行barcoder、UMI进行处理，
可以进行定量分析，生成barcodes.tsv、features.tsv、matrix.mtx文件，用于单细胞数据的分析。

## 基因组和注释文件下载

```bash
$ # 下载hg38基因组
$ mkdir hg38
$ wget -P hg38/ -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
$ gunzip hg38/hg38.analysisSet.fa.gz

$ # 下载小鼠基因组
$ mkdir mm10
$ wget -P mm10/ -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
$ gunzip mm10/mm10.fa.gz

$ # 生成混合的基因组文件
$ mkdir mix_hg38_mm10
$ cat <(sed 's/chr/hg38_chr/' hg38/hg38.analysisSet.fa) \
$    <(sed 's/chr/mm10_chr/' mm10/mm10.fa) \
$    > mix_hg38_mm10/genome.fa

$ # hg38注释文件 https://www.gencodegenes.org/human/
$ wget -P hg38/ -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
$ gunzip hg38/gencode.v41.annotation.gtf.gz

$ # 小鼠注释文件 https://www.gencodegenes.org/mouse/
$ wget -P mm10/ -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
$ gunzip mm10/gencode.vM25.annotation.gtf.gz

$ # 创建混合的注释文件
$ cat <(sed 's/^chr/hg38_chr/' hg38/gencode.v41.annotation.gtf) \
$     <(sed 's/^chr/mm10_chr/' mm10/gencode.vM25.annotation.gtf) \
$     > mix_hg38_mm10/annotation.gtf

```

## 建立索引

