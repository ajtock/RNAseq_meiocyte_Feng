#!/bin/bash

# Example usage via condor submission system on hydrogen node7
# csmit -m 1G -c 1 "bash fastq_read_counts.sh HGXF012M HGXF012N HGXF012U HGXF012O HGXF012P HGXF012V HGXF012W HGXF012X HGXF012Y"

# Output files with the suffix "fastq.stats" contain read counts,
for i in "$@"
do
( gunzip -k ${i}.fastq.gz
  cat ${i}.fastq | echo $((`wc -l`/4)) > ${i}.fastq.stats ) &
done
wait
