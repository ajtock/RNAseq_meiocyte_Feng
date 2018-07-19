#!/bin/bash
# Quantify gene expression levels using salmon 0.9.1

# Example usage via condor submission system on hydrogen node7
# csmit -m 1G -c 24 "bash salmon_quant.sh Col_RNAseq_meiocyte"

prefix=$1

for i in $prefix"_Rep"{1..3};
do
  samp=`basename ${i}`
  echo "Processing sample ${samp}"
  salmon quant -i /projects/ajt200/TAIR10/salmon_transcriptome_index -l A \
               -r ${i}.fastq.gz \
               -p 24 \
               -o ./salmon_quants/${samp}_quant
done

