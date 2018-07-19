#!/bin/bash

## ************** TEST BEFORE USE ON IMPORTANT FASTQ FILES ******************
# Example usage via condor submission system on hydrogen node7:
# csmit -m 20G -c 1 "bash zcat_fastq_gz_gzip_best.sh HGXF012M_S1 HGXF012M_S4 HGXF012M"

run1=$1
run2=$2
name=$3

if [ ! -f "$name.fastq.gz" ]; then 
    zcat $run1"_L001_R1_001.fastq.gz" $run1"_L002_R1_001.fastq.gz" $run1"_L003_R1_001.fastq.gz" $run1"_L004_R1_001.fastq.gz" \
         $run2"_L001_R1_001.fastq.gz" $run2"_L002_R1_001.fastq.gz" $run2"_L003_R1_001.fastq.gz" $run2"_L004_R1_001.fastq.gz" \
    | gzip -c -k --best > $name.fastq.gz;
else 
    echo "skipping $name"
fi

