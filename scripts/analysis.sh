#!/bin/bash

# Set the path to the reference genome
REF=reference.fasta

# Index the reference genome
bwa index reference.fasta # takes about 3 hours

# Create an output directory for the results
mkdir -p output

# Loop through each plate
for sample in tu wt
do
    # Set the input filenames for this plate
    R1=data/$sample.r1.fq.gz
    R2=data/$sample.r2.fq.gz

    # Align the reads to the reference genome
    bwa mem $REF $R1 $R2 | samtools sort -O BAM -o output/aln_$sample.sorted.bam -

    # Index the alignment file
    samtools index output/aln_$sample.sorted.bam

    # Save only some region
    samtools view -b - chrX:20000000-40000000 > aln_$sample_reg.sorted.bam

    # Index the region file
    samtools index aln_$sample_reg.sorted.bam

    # Get depth of coverage
    samtools depth aln_$sample_reg.sorted.bam > $sample.txt

# Run python script to create plot
python plot.py --tumor $sample.txt --normal $sample.png

done
