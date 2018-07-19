#!/applications/R/R-3.4.0/bin/Rscript

### Create data.frame of differentially expressed genes with TPM values for
### each sample (Col, taf4b, Bur0, leaf, meiocyte)

# R version 3.4.0 (invoke this version by specifying path at command line:
# /applications/R/R-3.4.0/bin/R , or if running as a script:
# /applications/R/R-3.4.0/bin/Rscript)

# Usage: Rscript /projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/taf4b_downReg_TPM_table.R 0.01

DEGdir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/"

args <- commandArgs(trailingOnly = TRUE)
FDR <- as.character(args[1])

# Genes down-regulated in taf4b vs Col
taf4b_downReg <- read.table(paste0(DEGdir,
                                   "res_taf4bVCol_",
                                   FDR, "_lfcShrink_Chr_Sig",
                                   FDR, "_downRegSortedDF.txt"),
                            header = T)
print(dim(taf4b_downReg))
#[1] 1271    5

# Genes up-regulated in meiocytes vs leaf
meioDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/DESeq2_analysis_01/"
meio_upReg <- read.table(paste0(meioDir,
                               "res_meiocyteVleaf_", FDR,
                               "_lfcShrink_Chr_Sig", FDR,
                               "_upRegSortedDF.txt"))
print(dim(meio_upReg))
#[1] 4528    5

Col_taf4b_Bur0_TPM <- read.table("/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_03/salmon_TPM_RNAseq_Col_taf4b_Bur0.txt")
leaf_meio_TPM <- read.table(paste0(meioDir, "salmon_TPM_WT_RNAseq_leaf_and_meiocyte.txt"))

taf4b_downReg_new <- cbind(rownames(taf4b_downReg), taf4b_downReg)
meio_upReg_new <- cbind(rownames(meio_upReg), meio_upReg)
Col_taf4b_Bur0_TPM_new <- cbind(rownames(Col_taf4b_Bur0_TPM), Col_taf4b_Bur0_TPM)
leaf_meio_TPM_new <- cbind(rownames(leaf_meio_TPM), leaf_meio_TPM)

colnames(taf4b_downReg_new) <- c("geneID", paste0("taf4b_", colnames(taf4b_downReg)))
colnames(meio_upReg_new) <- c("geneID", paste0("meio_", colnames(meio_upReg)))
colnames(Col_taf4b_Bur0_TPM_new) <- c("geneID", paste0("TPM_", colnames(Col_taf4b_Bur0_TPM)))
colnames(leaf_meio_TPM_new) <- c("geneID", paste0("TPM_", colnames(leaf_meio_TPM)))

# Obtain intersection of genes down-regulated in taf4b vs Col
# and up-regulated in meiocytes vs leaf tissue
taf4b_downReg_meio_upReg <- merge(x = taf4b_downReg_new,
                                  y = meio_upReg_new,
                                  by.x = "geneID",
                                  by.y = "geneID")
# Append TPM values
taf4b_downReg_meio_upReg_TPM1 <- merge(x = taf4b_downReg_meio_upReg,
                                       y = Col_taf4b_Bur0_TPM_new,
                                       by.x = "geneID",
                                       by.y = "geneID")
taf4b_downReg_meio_upReg_TPM2 <- merge(x = taf4b_downReg_meio_upReg_TPM1,
                                       y = leaf_meio_TPM_new,
                                       by.x = "geneID",
                                       by.y = "geneID")
# Sort by ascending taf4b_log2FoldChange
taf4b_downReg_meio_upReg_TPM <- taf4b_downReg_meio_upReg_TPM2[order(taf4b_downReg_meio_upReg_TPM2$taf4b_log2FoldChange,
                                                                    decreasing = F),]
print(dim(taf4b_downReg_meio_upReg_TPM))
#[1] 646  26
write.table(taf4b_downReg_meio_upReg_TPM,
            file = paste0(DEGdir,
                          "taf4b_downReg_FDR", FDR,
                          "_meio_upReg_FDR", FDR,
                          "_TPM_table_sorted_by_ascending_taf4b_log2FoldChange.txt"),
            quote = F, sep = "\t", row.names = F)

TPM_Col_mean <- rowMeans(cbind(taf4b_downReg_meio_upReg_TPM$TPM_Col_Rep1,
                               taf4b_downReg_meio_upReg_TPM$TPM_Col_Rep3))
TPM_taf4b_mean <- rowMeans(cbind(taf4b_downReg_meio_upReg_TPM$TPM_Col_Rep1,
                                 taf4b_downReg_meio_upReg_TPM$TPM_Col_Rep2,
                                 taf4b_downReg_meio_upReg_TPM$TPM_Col_Rep3))
TPM_leaf_mean <- rowMeans(cbind(taf4b_downReg_meio_upReg_TPM$TPM_leaf_Rep1,
                                taf4b_downReg_meio_upReg_TPM$TPM_leaf_Rep2,
                                taf4b_downReg_meio_upReg_TPM$TPM_leaf_Rep3))
TPM_meiocyte_mean <- rowMeans(cbind(taf4b_downReg_meio_upReg_TPM$TPM_meiocyte_Rep1,
                                    taf4b_downReg_meio_upReg_TPM$TPM_meiocyte_Rep2,
                                    taf4b_downReg_meio_upReg_TPM$TPM_meiocyte_Rep3))
TPM_log2_taf4b_Col <- log2((TPM_taf4b_mean+1)/(TPM_Col_mean+1))
TPM_log2_meiocyte_leaf <- log2((TPM_meiocyte_mean+1)/(TPM_leaf_mean+1))
 
pdf(paste0(DEGdir,
           "taf4b_downReg_FDR", FDR,
           "_meio_upReg_FDR", FDR,
           "_log2FoldChange_scatterplot.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plot(x = taf4b_downReg_meio_upReg_TPM$taf4b_log2FoldChange,
     y = taf4b_downReg_meio_upReg_TPM$meio_log2FoldChange,
     pch = 20, cex = 0.7,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor.test(taf4b_downReg_meio_upReg_TPM$taf4b_log2FoldChange,
                                                             taf4b_downReg_meio_upReg_TPM$meio_log2FoldChange,
                                                             method = "spearman")$estimate[[1]],
                                                             digits = 2))*"," ~ 
                                                             italic("P") ~ " = " ~ .(cor.test(taf4b_downReg_meio_upReg_TPM$taf4b_log2FoldChange,
                                                                                              taf4b_downReg_meio_upReg_TPM$meio_log2FoldChange,
                                                                                              method = "spearman")$p.value)),
     xlab = expression(taf4b~vs~Col~log[2]~fold~change),
     ylab = expression(Meiocyte~vs~leaf~log[2]~fold~change))
dev.off()

pdf(paste0(DEGdir,
           "taf4b_downReg_FDR", FDR,
           "_meio_upReg_FDR", FDR,
           "_adjustedP_scatterplot.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plot(x = taf4b_downReg_meio_upReg_TPM$taf4b_padj,
     y = taf4b_downReg_meio_upReg_TPM$meio_padj,
     pch = 20, cex = 0.7,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor.test(taf4b_downReg_meio_upReg_TPM$taf4b_padj,
                                                             taf4b_downReg_meio_upReg_TPM$meio_padj,
                                                             method = "spearman")$estimate[[1]],
                                                             digits = 2))*"," ~
                                                             italic("P") ~ " = " ~ .(cor.test(taf4b_downReg_meio_upReg_TPM$taf4b_padj,
                                                                                              taf4b_downReg_meio_upReg_TPM$meio_padj,
                                                                                              method = "spearman")$p.value)),
     xlab = expression(taf4b~vs~Col~FDR),
     ylab = expression(Meiocyte~vs~leaf~FDR))
dev.off()

pdf(paste0(DEGdir,
           "taf4b_downReg_FDR", FDR,
           "_meio_upReg_FDR", FDR,
           "_log2TPM_scatterplot.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plot(x = TPM_log2_taf4b_Col,
     y = TPM_log2_meiocyte_leaf,
     pch = 20, cex = 0.7,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor.test(TPM_log2_taf4b_Col,
                                                             TPM_log2_meiocyte_leaf,
                                                             method = "spearman")$estimate[[1]],
                                                             digits = 2))*"," ~
                                                             italic("P") ~ " = " ~ .(cor.test(TPM_log2_taf4b_Col,
                                                                                              TPM_log2_meiocyte_leaf,
                                                                                              method = "spearman")$p.value)),
     xlab = expression(taf4b~vs~Col~log[2]~TPM~fold~change),
     ylab = expression(Meiocyte~vs~leaf~log[2]~TPM~fold~change))
dev.off()


# For all genes down-regulated in taf4b vs Col, append TPM values
taf4b_downReg_TPM1 <- merge(x = taf4b_downReg_new,
                            y = Col_taf4b_Bur0_TPM_new,
                            by.x = "geneID",
                            by.y = "geneID")
taf4b_downReg_TPM2 <- merge(x = taf4b_downReg_TPM1,
                            y = leaf_meio_TPM_new,
                            by.x = "geneID",
                            by.y = "geneID")
# Sort by ascending taf4b_log2FoldChange
taf4b_downReg_TPM <- taf4b_downReg_TPM2[order(taf4b_downReg_TPM2$taf4b_log2FoldChange,
                                              decreasing = F),]
print(dim(taf4b_downReg_TPM))
#[1] 1271   21
write.table(taf4b_downReg_TPM,
            file = paste0(DEGdir,
                          "taf4b_downReg_FDR", FDR,
                          "_TPM_table_sorted_by_ascending_taf4b_log2FoldChange.txt"),
            quote = F, sep = "\t", row.names = F)

# Warning: the following redefines previously defined variables
TPM_Col_mean <- rowMeans(cbind(taf4b_downReg_TPM$TPM_Col_Rep1,
                               taf4b_downReg_TPM$TPM_Col_Rep3))
TPM_taf4b_mean <- rowMeans(cbind(taf4b_downReg_TPM$TPM_Col_Rep1,
                                 taf4b_downReg_TPM$TPM_Col_Rep2,
                                 taf4b_downReg_TPM$TPM_Col_Rep3))
TPM_leaf_mean <- rowMeans(cbind(taf4b_downReg_TPM$TPM_leaf_Rep1,
                                taf4b_downReg_TPM$TPM_leaf_Rep2,
                                taf4b_downReg_TPM$TPM_leaf_Rep3))
TPM_meiocyte_mean <- rowMeans(cbind(taf4b_downReg_TPM$TPM_meiocyte_Rep1,
                                    taf4b_downReg_TPM$TPM_meiocyte_Rep2,
                                    taf4b_downReg_TPM$TPM_meiocyte_Rep3))
TPM_log2_taf4b_Col <- log2((TPM_taf4b_mean+1)/(TPM_Col_mean+1))
TPM_log2_meiocyte_leaf <- log2((TPM_meiocyte_mean+1)/(TPM_leaf_mean+1))

pdf(paste0(DEGdir,
           "taf4b_downReg_FDR", FDR,
           "_log2TPM_scatterplot.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plot(x = TPM_log2_taf4b_Col,
     y = TPM_log2_meiocyte_leaf,
     pch = 20, cex = 0.7,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor.test(TPM_log2_taf4b_Col,
                                                             TPM_log2_meiocyte_leaf,
                                                             method = "spearman")$estimate[[1]],
                                                             digits = 2))*"," ~
                                                             italic("P") ~ " = " ~ .(cor.test(TPM_log2_taf4b_Col,
                                                                                              TPM_log2_meiocyte_leaf,
                                                                                              method = "spearman")$p.value)),
     xlab = expression(taf4b~vs~Col~log[2]~TPM~fold~change),
     ylab = expression(Meiocyte~vs~leaf~log[2]~TPM~fold~change))
dev.off()


pdf(paste0(DEGdir,
           "taf4b_downReg_FDR", FDR,
           "_log2FoldChange_vs_log2TPM_scatterplot.pdf"), height = 4, width = 4)
par(mfrow = c(1, 1), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plot(x = TPM_log2_taf4b_Col,
     y = TPM_log2_meiocyte_leaf,
     pch = 20, cex = 0.7,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor.test(taf4b_downReg_TPM$taf4b_log2FoldChange,
                                                             TPM_log2_meiocyte_leaf,
                                                             method = "spearman")$estimate[[1]],
                                                             digits = 2))*"," ~
                                                             italic("P") ~ " = " ~ .(cor.test(TPM_log2_taf4b_Col,
                                                                                              TPM_log2_meiocyte_leaf,
                                                                                              method = "spearman")$p.value)),
     xlab = expression(taf4b~vs~Col~log[2]~fold~change),
     ylab = expression(Meiocyte~vs~leaf~log[2]~TPM~fold~change))
dev.off()

