#!/bin/bash
# STAR version 2.5.3a
# samtools version 1.3
# R version 3.3.2

# Example usage via condor submission system on hydrogen node7
# csmit -m 20G -c 1 "bash samtools_PE_RNAseq_mapping_pipeline.sh WT_RNAseq_ATCACG"

i=$1

# Specify directories containing foo.py and multi_unique_extract_pairend.r scripts
pyDir=/projects/ajt200/PyScripts
RDir=/projects/ajt200/mapping_pipeline/multi_unique_extract_Rscripts

# Convert sam to bam
samtools view -bS -o ${i}_STAR.bam ${i}_Aligned.out.sam
# Include header (-h) and retain only properly aligned read pairs (-f 0x02)
samtools view -b -hf 0x02 ${i}_STAR.bam > ${i}_STAR_mapped.bam
# Convert bam to sam and include header
samtools view -h -o ${i}_STAR_mapped.sam ${i}_STAR_mapped.bam
# Write sam header lines to file, to be appended to unique and multiple alignments
samtools view -H ${i}_STAR_mapped.sam > ${i}_STAR_mapped_header.sam

### Identify and extract unique alignments

# Extract alignments for which STAR MAPQ score == 255, which corresponds to uniquely mapped reads \
  # Retain alignments for which the names of both reads in a pair are the same
samtools view -S -q 255 ${i}_STAR_mapped.sam \
  | $pyDir/foo_unixEOL.py > ${i}_STAR_mapped_unique.txt 
# Test whether python script in previous command is required (possibly applicable only to Bowtie2 alignments)
samtools view -S -q 255 ${i}_STAR_mapped.sam > ${i}_STAR_mapped_unique_TEST.txt
diff ${i}_STAR_mapped_unique.txt ${i}_STAR_mapped_unique_TEST.txt > ${i}_STAR_mapped_unique_vs_TEST.diff 2>&1 
# Alternatively, uniquely aligned reads can be obtained with a grep operation on the "NH:i:Nmap" field
  # ([^0-9] matches characters not in the range of 0 to 9; in this has the effect of excluding alignments whose reads map to 1# loci) \
  # Retain alignments for which the names of both reads in a pair are the same \
  # Convert Windows text file to a Unix text file \
samtools view -S ${i}_STAR_mapped.sam \
  | grep -e "NH:i:[1][^0-9]" \
  | $pyDir/foo.py \
  | awk '{ sub("\r$", ""); print }' > ${i}_STAR_mapped_unique_alt.txt
# Check that both methods for extracting unique alignments yield the same result
diff ${i}_STAR_mapped_unique.txt ${i}_STAR_mapped_unique_alt.txt > ${i}_STAR_mapped_unique_vs_alt.diff 2>&1
# Test whether python script in above command is required
samtools view -S ${i}_STAR_mapped.sam \
  | grep -e "NH:i:[1][^0-9]" > ${i}_STAR_mapped_unique_alt_TEST.txt
diff ${i}_STAR_mapped_unique_alt.txt ${i}_STAR_mapped_unique_alt_TEST.txt > ${i}_STAR_mapped_unique_alt_vs_TEST.diff 2>&1

# Concatenate header lines and unique alignments
cat ${i}_STAR_mapped_header.sam ${i}_STAR_mapped_unique.txt > ${i}_STAR_mapped_unique.sam
# Convert sam to bam
samtools view -bS -o ${i}_STAR_mapped_unique.bam ${i}_STAR_mapped_unique.sam
# Sort and index the bam
samtools sort ${i}_STAR_mapped_unique.bam -o ${i}_STAR_mapped_unique_sort.bam
samtools index ${i}_STAR_mapped_unique_sort.bam


### Extract multiply aligned reads and, for each read, extract the alignment with the best score

# Multiply aligned reads can be obtained with a grep -v operation on the "NH:i:Nmap" field
  # ([^0-9] matches characters not in the range of 0 to 9; in this case, this has the effect of excluding alignments whose reads map to 1 loci) \
  # Retain alignments for which the names of both reads in a pair are the same
samtools view -S ${i}_STAR_mapped.sam \
  | grep -v -e "NH:i:[1][^0-9]" \
  | $pyDir/foo_unixEOL.py > ${i}_STAR_mapped_multi.txt
# Test whether python script in previous command is required (possibly applicable only to Bowtie2 alignments)
samtools view -S ${i}_STAR_mapped.sam \
  | grep -v -e "NH:i:[1][^0-9]" >  ${i}_STAR_mapped_multi_TEST.txt
diff ${i}_STAR_mapped_multi.txt ${i}_STAR_mapped_multi_TEST.txt > ${i}_STAR_mapped_multi_vs_TEST.diff 2>&1

# Concatenate header lines and multiple alignments
cat ${i}_STAR_mapped_header.sam ${i}_STAR_mapped_multi.txt > ${i}_STAR_mapped_multi.sam
# Convert sam to bam
samtools view -bS -o ${i}_STAR_mapped_multi.bam ${i}_STAR_mapped_multi.sam

# For each multiply aligned read, extract the primary alignment (that with the best score)
samtools view -S -h -F 0x100 ${i}_STAR_mapped_multi.sam > ${i}_STAR_mapped_multi_primary_alignment.sam
# Convert sam to bam
samtools view -bS -o ${i}_STAR_mapped_multi_primary_alignment.bam ${i}_STAR_mapped_multi_primary_alignment.sam
# Extract primary alignments with a MAPQ score == 3 (one or both reads in alignment map to 2 loci)
  # Retain alignments for which the names of both reads in a pair are the same 
samtools view -S -q 3 ${i}_STAR_mapped_multi_primary_alignment.sam \
  | $pyDir/foo_unixEOL.py > ${i}_STAR_mapped_multi_primary_alignment_fq3.txt 
# Test whether python script in previous command is required (possibly applicable only to Bowtie2 alignments)
samtools view -S -q 3 ${i}_STAR_mapped_multi_primary_alignment.sam > ${i}_STAR_mapped_multi_primary_alignment_fq3_TEST.txt
diff ${i}_STAR_mapped_multi_primary_alignment_fq3.txt ${i}_STAR_mapped_multi_primary_alignment_fq3_TEST.txt > ${i}_STAR_mapped_multi_primary_alignment_fq3_vs_TEST.diff 2>&1

# Concatenate header lines and multiple alignments
cat ${i}_STAR_mapped_header.sam ${i}_STAR_mapped_multi_primary_alignment_fq3.txt > ${i}_STAR_mapped_multi_primary_alignment_fq3.sam
# Convert sam to bam
samtools view -bS -o ${i}_STAR_mapped_multi_primary_alignment_fq3.bam ${i}_STAR_mapped_multi_primary_alignment_fq3.sam
# Sort and index the bam
samtools sort ${i}_STAR_mapped_multi_primary_alignment_fq3.bam -o ${i}_STAR_mapped_multi_primary_alignment_fq3_sort.bam
samtools index ${i}_STAR_mapped_multi_primary_alignment_fq3_sort.bam


### Merge unique and multiple primary alignments

samtools merge ${i}_STAR_mapped_both.bam ${i}_STAR_mapped_unique_sort.bam ${i}_STAR_mapped_multi_primary_alignment_fq3_sort.bam
# Sort and index the bam
samtools sort ${i}_STAR_mapped_both.bam -o ${i}_STAR_mapped_both_sort.bam  
samtools index ${i}_STAR_mapped_both_sort.bam

 
### Calculate alignment summary statistics for fastq files and for each bam file generated by the pipeline

# Output files with the suffix "fastq.stats" or "bam.stats" contain read-pair counts,
# while output files with the suffix "alignments.stats" contain alignment counts
# This is to provide distinct alignment counts where reads map to multiple loci
cat ${i}_L002_R1.fastq | echo $((`wc -l`/4)) > ${i}_L002_R1.fastq.stats
cat ${i}_L002_R3.fastq | echo $((`wc -l`/4)) > ${i}_L002_R3.fastq.stats
samtools view -F 0x4 ${i}_STAR.bam | cut -f 1 | sort | uniq | wc -l > ${i}_STAR.bam.stats
samtools view -F 0x4 ${i}_STAR_mapped.bam | cut -f 1 | sort | uniq | wc -l > ${i}_STAR_mapped.bam.stats
samtools view -F 0x4 ${i}_STAR_mapped_unique_sort.bam | cut -f 1 | sort | uniq | wc -l > ${i}_STAR_mapped_unique_sort.bam.stats
samtools view -c -f 0x40 -F 0x4 ${i}_STAR_mapped_multi.bam > ${i}_STAR_mapped_multi.bam.alignments.stats
samtools view -c -f 0x40 -F 0x4 ${i}_STAR_mapped_multi_primary_alignment.bam > ${i}_STAR_mapped_multi_primary_alignment.bam.alignments.stats
# Sanity check; (read-pair) count generated by the following command should be the same as
# (alignment) count generated by the previous command
samtools view -F 0x4 ${i}_STAR_mapped_multi_primary_alignment.bam | cut -f 1 | sort | uniq | wc -l > ${i}_STAR_mapped_multi_primary_alignment.bam.stats
samtools view -F 0x4 ${i}_STAR_mapped_multi_primary_alignment_fq3_sort.bam | cut -f 1 | sort | uniq | wc -l > ${i}_STAR_mapped_multi_primary_alignment_fq3_sort.bam.stats
samtools view -F 0x4 ${i}_STAR_mapped_both_sort.bam | cut -f 1 | sort | uniq | wc -l > ${i}_STAR_mapped_both_sort.bam.stats

paste -d "\t" ${i}_L002_R1.fastq.stats ${i}_STAR_mapped.bam.stats ${i}_STAR_mapped_unique_sort.bam.stats ${i}_STAR_mapped_multi.bam.alignments.stats ${i}_STAR_mapped_multi_primary_alignment.bam.alignments.stats ${i}_STAR_mapped_multi_primary_alignment.bam.stats ${i}_STAR_mapped_multi_primary_alignment_fq3_sort.bam.stats ${i}_STAR_mapped_both_sort.bam.stats > ${i}_alignment_summary_stats.txt
sed -i '1i Total_sequenced_read_pairs\tAligned,_mismatches<=2\tUniquely_aligned\tMultiply_aligned_read_pair_alignments\tPrimary_multiply_aligned_read_pair_alignments\tPrimary_multiply_aligned\tPrimary_multiply_aligned,_MAPQ==3\tBoth' ${i}_alignment_summary_stats.txt 
[ -d alignment_stats ] || mkdir alignment_stats
mv ${i}*.stats alignment_stats/
mv ${i}_alignment_summary_stats.txt alignment_stats/


### Remove no longer required text and sam files

