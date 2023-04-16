# Cancer-Genomics-Data-Analysis
Cancer Genomics Data Analysis exercise for the course Analytical Methods in Cancer Genomics 2023. The goal of the exercise is to generate a read-depth plot for human cancer sample.

## Analysis Pipeline
1. Index the reference genome (reference/hg19) using bwa index.
2. Align the reads to the human reference genome (reference/hg19) using bwa mem.
3. Sort the aligned reads using samtools sort.
4. Index the sorted BAM file using samtools index.
5. Subset the BAM file to the region of interest (chrX:20000000-40000000) using samtools view.
6. Retrieve the read-depth information for the region of interest using samtools depth.
7. Generate a read-depth plot using custom Python script.

## Installation and Dependencies

- bwa (0.7.17-r1188)
- samtools (v1.17)

For Python script:
- pandas (v1.5.3)
- matplotlib (v3.7.1)
- numpy (v1.24.2)


## Usage
The pipeline can be run using the following command:
```
bash scripts/analysis.sh
```

### Notes
- Running genome indexing and alignment is extremely computationally expensive for human genome. The indexing of genome using following command took me ~3 hours (14442.15 sec):

```
bwa index reference/hg19.fa
```

The mapping of reads to the reference genome using following command took another 40 minutes (2362.867 seconds)

```bash
bwa mem -M -t 8 $REF sample_data/tu.r1.fq.gz sample_data/tu.r2.fq.gz
```

Is it how it was supposed to be or I did something wrong? What does it mean to downsample the genome?

- First I tried to use Delly tool but I failed to install it because of the following error:

```
delly: error while loading shared libraries: libboost_iostreams.so.1.67.0: cannot open shared object file: No such file or directory
```
Same issue is discussed [here](https://github.com/jodyphelan/TBProfiler/issues/88).