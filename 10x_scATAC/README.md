## 10x single cell ATAC analysis 
查看每个细胞中的表观组学，判断细胞类型和功能
### 原始数据下载和处理

```bash
mkdir -p 10xscATAC/pbmc500
wget -P 10xscATAC/pbmc500 -c https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fastqs.tar
tar xf 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs.tar -C 10xscATAC/pbmc500/

```

### Prepare Whitelist

```bash
# download the whitelist
wget -P 10xscATAC/pbmc500/ https://teichlab.github.io/scg_lib_structs/data/737K-cratac-v1.txt.gz
gunzip 10xscATAC/pbmc500/737K-cratac-v1.txt.gz

# reverse complement the whitelist
cat 10xscATAC/pbmc500/737K-cratac-v1.txt | \
    rev | tr 'ACGT' 'TGCA' > \
    10xscATAC/pbmc500/737K-cratac-v1_rc.txt
```

### From FastQ To Fragments
chromap 

```bash
mkdir -p 10xscATAC/chromap_outs

map and generate the fragment file

chromap -t 4 --preset atac \
        -x hg38/chromap_index/genome.index \
        -r hg38/hg38.analysisSet.fa \
        -1 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L001_R1_001.fastq.gz,10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L002_R1_001.fastq.gz \
        -2 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L001_R3_001.fastq.gz,10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L002_R3_001.fastq.gz \
        -b 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L001_R2_001.fastq.gz,10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L002_R2_001.fastq.gz \
        --barcode-whitelist 10xscATAC/pbmc500/737K-cratac-v1.txt \
        -o 10xscATAC/chromap_outs/fragments.tsv

compress and index the fragment file

bgzip 10xscATAC/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed 10xscATAC/chromap_outs/fragments.tsv.gz
```

### Fragments To Reads

```bash
we also sort the output by chromosome name
which will be useful later

faSize -detailed hg38/hg38.analysisSet.fa | \
    sort -k1,1 > hg38/hg38.chrom.sizes

# we use bedClip to remove reads outside the chromosome boundary
# we also remove reads mapped to the mitochondrial genome (chrM)

zcat 10xscATAC/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin hg38/hg38.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > 10xscATAC/chromap_outs/reads.bed.gz

```

### Peak Calling By MACS2

```bash
macs2 callpeak -t 10xscATAC/chromap_outs/reads.bed.gz \
               -g hs -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir 10xscATAC/chromap_outs \
               -n aggregate

```


### Find Reads In Peaks Per Cell

```bash
# format and sort peaks

cut -f 1-4 10xscATAC/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > 10xscATAC/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a 10xscATAC/chromap_outs/aggregate_peaks_sorted.bed \
    -b 10xscATAC/chromap_outs/reads.bed.gz \
    -wo -sorted -g hg38/hg38.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > 10xscATAC/chromap_outs/peak_read_ov.tsv.gz
```
### Output The Peak-By-Cell Matrix

```bash
# create dirctory
mkdir -p 10xscATAC/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 10xscATAC/chromap_outs/aggregate_peaks_sorted.bed > \
    10xscATAC/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat 10xscATAC/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    10xscATAC/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
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

```bash
python generate_csc_mtx.py \
    10xscATAC/chromap_outs/aggregate_peaks_sorted.bed \
    10xscATAC/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    10xscATAC/chromap_outs/peak_read_ov.tsv.gz \
    10xscATAC/chromap_outs/raw_peak_bc_matrix
```

### Cell Calling (Filter Cell Barcodes)

```bash
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed 10xscATAC/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     10xscATAC/chromap_outs/raw_peak_bc_matrix \
     10xscATAC/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv 10xscATAC/chromap_outs/filtered_peak_bc_matrix/peaks.bed

```

chromap+macs2 和 cellranger-atac相比，执行更快
| Comparison  | chromap + macs2 | cellranger-atac| Overlap |
| :---        |      :----:     |      :----:    |    ---: |
| Peak number |     71,000      |     65,908     |  61,042 |
| Cell number |      643        |       484      |    482  |
