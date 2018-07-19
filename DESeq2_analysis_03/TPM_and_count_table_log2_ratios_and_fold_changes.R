#!/applications/R/R-3.3.2/bin/Rscript

inDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_03/"

TPM <- read.table(paste0(inDir, "salmon_TPM_RNAseq_Col_taf4b_Bur0.txt"))
# Remove replicate 2
TPM <- TPM[,-2]
TPM <- cbind(TPM,
             rowMeans(TPM[,1:2]),
             rowMeans(TPM[,3:5]),
             rowMeans(TPM[,6:8]))
colnames(TPM) <- c(colnames(TPM)[1:8], "Col_mean", "taf4b_mean", "Bur0_mean")
TPM <- cbind(TPM,
             log2((TPM$taf4b_mean+1)/(TPM$Col_mean+1)),
             (TPM$taf4b_mean+1)/(TPM$Col_mean+1)-1,
             log2((TPM$Bur0_mean+1)/(TPM$Col_mean+1)),
             (TPM$Bur0_mean+1)/(TPM$Col_mean+1)-1) 
colnames(TPM) <- c(colnames(TPM)[1:11], "taf4b_log2_ratio", "taf4b_fold_change", "Bur0_log2_ratio","Bur0_fold_change")

write.table(as.data.frame(TPM),
            file = paste0(inDir, "salmon_TPM_log2_ratio_and_fold_change_RNAseq_Col_taf4b_Bur0.txt"),
            sep = "\t", quote = F)


counts <- read.table(paste0(inDir, "salmon_counts_RNAseq_Col_taf4b_Bur0.txt"))
# Remove replicate 2
counts <- counts[,-2]
counts <- cbind(counts,
               rowMeans(counts[,1:2]),
               rowMeans(counts[,3:5]),
               rowMeans(counts[,6:8]))
colnames(counts) <- c(colnames(counts)[1:8], "Col_mean", "taf4b_mean", "Bur0_mean")
counts <- cbind(counts,
                log2((counts$taf4b_mean+1)/(counts$Col_mean+1)),
                (counts$taf4b_mean+1)/(counts$Col_mean+1)-1,
                log2((counts$Bur0_mean+1)/(counts$Col_mean+1)),
                (counts$Bur0_mean+1)/(counts$Col_mean+1)-1)
colnames(counts) <- c(colnames(counts)[1:11], "taf4b_log2_ratio", "taf4b_fold_change", "Bur0_log2_ratio","Bur0_fold_change")

write.table(as.data.frame(counts),
            file = paste0(inDir, "salmon_counts_log2_ratio_and_fold_change_RNAseq_Col_taf4b_Bur0.txt"),
            sep = "\t", quote = F)


