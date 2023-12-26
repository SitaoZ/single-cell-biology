## STAR 
STAR 软件可以用来处理10xgenomics的单细胞测序数据，可以对数据进行barcoder、UMI进行处理，
可以进行定量分析，生成barcodes.tsv、features.tsv、matrix.mtx文件，用于单细胞数据的分析。
STRA可以生成和CellRanger类似的结果，可以参考[STAR issue](https://github.com/alexdobin/STAR/issues/2026) [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#how-to-make-starsolo-raw-gene-counts-almost-identical-to-cellrangers)

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

```bash
$ # hg38
$ STAR --runThreadN 4 \
$      --runMode genomeGenerate \
$      --genomeDir hg38/star_index \
$      --genomeFastaFiles hg38/hg38.analysisSet.fa \
$      --sjdbGTFfile hg38/gencode.v41.annotation.gtf

$ # 小鼠
$ STAR --runThreadN 4 \
$      --runMode genomeGenerate \
$      --genomeDir mm10/star_index \
$      --genomeFastaFiles mm10/mm10.fa \
$      --sjdbGTFfile mm10/gencode.vM25.annotation.gtf

$ # 混合
$ STAR --runThreadN 4 \
$      --runMode genomeGenerate \
$      --genomeDir mix_hg38_mm10/star_index \
$      --genomeFastaFiles mix_hg38_mm10/genome.fa \
$      --sjdbGTFfile mix_hg38_mm10/annotation.gtf
```


## 单细胞测序数据建库参考
```bash
$ # 建库参考
$ # https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html
```

## 测序数据下载
```bash
$ 下载SRA文件
$ mkdir -p mereu2020/10xV3
$ wget -P mereu2020/10xV3 -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10587809/SRR10587809_1.fastq.gz
$ wget -P mereu2020/10xV3 -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10587809/SRR10587809_2.fastq.gz

$ # 下载barcodes文件
$ wget -P mereu2020/10xV3 wget https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz
$ gunzip mereu2020/10xV3/3M-february-2018.txt.gz
```



## 比对
```bash 
$ STAR --runThreadN 24 \
$      --genomeDir mix_hg38_mm10/star_index \
$      --readFilesCommand zcat \
$      --outFileNamePrefix mereu2020/star_outs/ \
$      --readFilesIn mereu2020/10xV3/SRR10587809_2.fastq.gz mereu2020/10xV3/SRR10587809_1.fastq.gz \
$      --soloType CB_UMI_Simple \
$      --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
$      --soloCBwhitelist mereu2020/10xV3/3M-february-2018.txt \
$      --soloCellFilter EmptyDrops_CR \
$      --soloStrand Forward \
$      --outSAMattributes CB UB \
$      --outSAMtype BAM SortedByCoordinate

$ # 参数解释
$ # --runThreadN 4 线程数
$ # --genomeDir mix_hg38_mm10/star_index 基因组索引文件
$ # 这里使用混合的star索引，是因为下载的数据来至HCA(human cell atlas)的样本，由人正常外周血单细胞(60%)、小鼠结肠(30%)、
$ # HEK293T细胞系(6%)、NIH3T3细胞系(3%) 和 狗MDCK细胞系(1%)。这里只分析占比高的的人和小鼠。
$ # --readFilesCommand zcat 输入文件是压缩文件
$ # --outFileNamePrefix mereu2020/star_outs/ # 结果输出文件夹
$ # --readFilesIn mereu2020/10xV3/SRR10587809_2.fastq.gz mereu2020/10xV3/SRR10587809_1.fastq.gz # fastq文件
$ # --soloType CB_UMI_Simple # 指定cell barcode和UMI的配置信息
$ # --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 # 指定cell barcode和UMI的信息，CB: cell barcode.
$ # soloCBstart、soloCBlen: cell barcode的起始和cell barcode的长度; soloUMIstart、soloUMIlen: UMI起始和UMI长度
$ # 注意，这些信息的起始都是1-based。
$ # 这四个信息和单细胞文库构建相关，只有清楚了解该文库的构建方式，才知道这些信息。
$ # --soloCBwhitelist mereu2020/10xV3/3M-february-2018.txt # 指定可用的barcode文件，文件的每一行是一个可用的barcode
$ # --soloCellFilter EmptyDrops_CR # 由于实验的效率问题，不是所有的液滴都含有细胞，因此需要指定相应的算法来处理
$ # --soloStrand Forward # 你的文库的cDNA来源的方向，
$ # Forward: cDNA reads are from the same strand as the mRNA (the coding strand),
$ # Reverse: If cDNA reads are from the opposite strand as the mRNA
$ # use Forward for 10x Genomics Single Cell 3’ V3 data
$ # --outSAMattributes CB UB # 指定cell barcode和UMI在SAM文件中的属性，便于后续分析
$ # --outSAMtype BAM SortedByCoordinate # 指定输出文件BAM格式，排序
```
