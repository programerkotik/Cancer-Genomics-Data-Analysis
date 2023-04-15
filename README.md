# Cancer-Genomics-Data-Analysis
Cancer Genomics Data Analysis exercise for the course Analytical Methods in Cancer Genomics 2023. The goal of the exercise is to generate a read-depth plot for human cancer sample.

## Analysis Pipeline
1. Download the paired tumor-normal sample pair data from the given source.
2. Trim the adapter sequences and low-quality reads using fastp.
3. Align the trimmed reads to the human reference genome (GRCh37/hg19) using bwa mem.
4. Sort the aligned reads using samtools sort.
5. Index the sorted BAM file using samtools index.
6. Subset the BAM file to the region of interest (chrX:20000000-40000000) using samtools view.
7. Generate a read-depth plot using samtools depth and gnuplot.

## Repository Structure
```
analysis_pipeline/
├── scripts/
│   ├── align_reads.sh
│   ├── generate_plot.sh
│   └── preprocess_reads.sh
├── output/
│   ├── alignment.bammain
│   ├── alignment.bam.bai
│   └── plot.png
├── README.md
├── requirements.txt
└── sample_data/
    ├── tumor_R1.fastq.gz
    ├── tumor_R2.fastq.gz
    ├── normal_R1.fastq.gz
    └── normal_R2.fastq.gz
```

## Installation and Dependencies

-bwa (v0.7.17)
-samtools (v1.13)
-fastp (v0.20.1)
-gnuplot (v5.2)

All dependencies can be installed using conda by creating a new environment and installing the packages listed in requirements.txt:

```
conda create -n cancer_genomics_env
conda activate cancer_genomics_env
conda install --file requirements.txt
```

## Usage
1. Clone the repository and navigate to the analysis_pipeline directory.
2. Activate the conda environment: conda activate cancer_genomics_env
3. Place the paired tumor-normal sample pair data in the sample_data directory.
4. Run the preprocess_reads.sh script to trim the adapter sequences and low-quality reads:
```
sh scripts/preprocess_reads.sh sample_data/tumor_R1.fastq.gz sample_data/tumor_R2.fastq.gz sample_data/normal_R1.fastq.gz sample_data/normal_R2.fastq.gz
```
5. Run the align_reads.sh script to align the reads to the human reference genome, sort and index the BAM file, and subset the BAM file to the region of interest:
```
sh scripts/align_reads.sh
```
6. Run the generate_plot.sh script to generate the read-depth plot:
```
sh scripts/generate_plot.sh
```
7.The resulting plot can be found in the output directory: output/plot.png.

### Notes
- This pipeline assumes that the data is in paired-end format.
- The align_reads.sh script assumes that the reference genome (GRCh37/hg19) is located in reference/hg19.fa. If using a different reference genome, this script should be modified accordingly.
- The generate_plot.sh script assumes that the BAM file is located in output/alignment.bam. If the BAM file is located elsewhere, this script should be modified accordingly.