#!/applications/R/R-3.4.0/bin
### Perform differential expression analysis using read counts generated by Salmon

# R version 3.4.0 (invoke this version by specifying path at command line: /applications/R/R-3.4.0/bin/R , or if running as a script: /applications/R/R-3.4.0/bin/Rscript)
# DESeq2 version 1.16.1
# Note that this R version or later is required for DESeq2 version 1.16 or later,
# in which "the log2 fold change shrinkage is no longer default for the DESeq and nbinomWaldTest functions".
# DESeq2 version 1.16 introduces "a separate function lfcShrink, which performs log2 fold change shrinkage
# for visualization and ranking of genes." (see https://support.bioconductor.org/p/95695/ and
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#changes)

## Import transcript abundance datasets with tximport package
# Code below based on https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

library(tximport)
print(packageVersion("tximport"))
#[1] ‘1.4.0’

inDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/"
outDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_01/"
plotDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_01/plots/"

# Read in table of sample IDs that will be used to specify paths to count files
samples <- read.table(file.path(inDir, "samples.txt"), header = T)
print(samples)
#  genotype     directory                           sample
#1      Col salmon_quants   Col_RNAseq_meiocyte_Rep1_quant
#2      Col salmon_quants   Col_RNAseq_meiocyte_Rep2_quant
#3      Col salmon_quants   Col_RNAseq_meiocyte_Rep3_quant
#4    taf4b salmon_quants taf4b_RNAseq_meiocyte_Rep1_quant
#5    taf4b salmon_quants taf4b_RNAseq_meiocyte_Rep2_quant
#6    taf4b salmon_quants taf4b_RNAseq_meiocyte_Rep3_quant
#7     Bur0 salmon_quants  Bur0_RNAseq_meiocyte_Rep1_quant
#8     Bur0 salmon_quants  Bur0_RNAseq_meiocyte_Rep2_quant
#9     Bur0 salmon_quants  Bur0_RNAseq_meiocyte_Rep3_quant

# Specify paths to count files
files <- file.path(inDir, samples$genotype, samples$directory, samples$sample, "quant.sf")
# Set 1:#_samples
names(files) <- paste0("sample", 1:9)
all(file.exists(files))
#[1] TRUE

# Create a dataframe of transcript IDs and corresponding gene IDs
transID <- read.table(files[1], colClasses = c(NA, rep("NULL", 4)), header = T)
tx2gene <- data.frame(cbind(as.vector(transID[,1]), substr(transID[,1], 1, 9)))
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Import transcript-level counts, summarised at gene level
# (reads that map to transcript IDs with a common parent gene ID are pooled)
library(readr)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
print(names(txi))
#[1] "abundance"           "counts"              "length"             
#[4] "countsFromAbundance"

# Import transcript-level counts, summarised at transcript level with "txOut = TRUE"
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
# Then summarise to gene level
txi.sum <- summarizeToGene(txi.tx, tx2gene)
# These two approaches should produce identical results
print(all.equal(txi$counts, txi.sum$counts))
#[1] TRUE

print(head(txi$counts))
#          sample1  sample2 sample3 sample4 sample5 sample6 sample7  sample8 sample9
#AT1G01010       2  2.00000       5       0       0       0       6  8.00000       5
#AT1G01020      52 29.12447     171       8       3       6     100 24.00000      59
#AT1G01030       1  0.00000       1       0       0       0       0  0.00000       0
#AT1G01040       2  9.00000      12      21       7      18      21 14.00000      16
#AT1G01050      11 25.00000      22      30       5       7     262 47.00000      87
#AT1G01060      10  7.00000       7       9       8       2      81 25.00434      41


## Use with downstream Bioconductor differential expression package
## DESeq2: http://www.bioconductor.org/help/workflows/rnaseqGene/

library(DESeq2)
print(packageVersion("DESeq2"))
#[1] ‘1.16.1’
sampleTable <- data.frame(sample = c("Col meiocyte RNA-seq Rep1", "Col meiocyte RNA-seq Rep2", "Col meiocyte RNA-seq Rep3",
                                     "taf4b meiocyte RNA-seq Rep1", "taf4b meiocyte RNA-seq Rep2", "taf4b meiocyte RNA-seq Rep3",
                                     "Bur0 meiocyte RNA-seq Rep1", "Bur0 meiocyte RNA-seq Rep2", "Bur0 meiocyte RNA-seq Rep3"),
                          condition = factor(rep(c("Col", "taf4b", "Bur0"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
print(sampleTable)
#                             sample condition
#sample1   Col meiocyte RNA-seq Rep1       Col
#sample2   Col meiocyte RNA-seq Rep2       Col
#sample3   Col meiocyte RNA-seq Rep3       Col
#sample4 taf4b meiocyte RNA-seq Rep1     taf4b
#sample5 taf4b meiocyte RNA-seq Rep2     taf4b
#sample6 taf4b meiocyte RNA-seq Rep3     taf4b
#sample7  Bur0 meiocyte RNA-seq Rep1      Bur0
#sample8  Bur0 meiocyte RNA-seq Rep2      Bur0
#sample9  Bur0 meiocyte RNA-seq Rep3      Bur0

dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, design = ~condition)

# Pre-filter the dataset
print(nrow(dds))
#[1] 27586
# Retain only rows that have more than a single count across all samples
dds <- dds[rowSums(counts(dds)) > 1,]
print(nrow(dds))
#[1] 24117


## The rlog and variance stabilizing transformations
# see http://www.bioconductor.org/help/workflows/rnaseqGene/#the-rlog-and-variance-stabilizing-transformations

rld <- rlog(dds, blind = FALSE)
print(head(assay(rld), 3))

vsd <- vst(dds, blind = FALSE)
print(head(assay(vsd), 3))

# Visualise the effect of transformation
library(dplyr)
library(ggplot2)
library(hexbin)

# For the log2 approach, estimate size factors to account for sequencing depth
# Sequencing-depth correction is done automatically for rlog and vst
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized = T)[,1:2]+1)) %>%
    mutate(transformation = "log2(normalized counts + 1)"),
  as_data_frame(assay(rld)[,1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[,1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("Sample 1 transformed counts", "Sample 2 transformed counts")

plot_transformed_counts <- ggplot(df, aes(x = `Sample 1 transformed counts`, y = `Sample 2 transformed counts`)) +
                             geom_hex(bins = 80) + coord_fixed() + facet_grid(. ~ transformation) +
                             labs(fill = "Occurrences") +
                             theme_classic()
                             #theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
ggsave(plot_transformed_counts,
       file = paste0(plotDir, "Sample1_vs_Sample2_transformed_counts_log2countsPlus1_rlog_vst.pdf"))


## Sample distances

# The plots generated below show that "Col meiocyte RNA-seq Rep2" is closer to Bur0 and taf4b than to other Col replicates
# Analysis repeated from raw data (per-lane fastq files) to confirm this unexpected result

# Sample distances using the rlog-transformed counts
sampleDists <- dist(t(assay(rld)))
print(sampleDists)

library(pheatmap)
library(RColorBrewer)

# Heatmap of sample distances using the rlog-transformed counts
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c("Col meiocyte RNA-seq Rep1", "Col meiocyte RNA-seq Rep2", "Col meiocyte RNA-seq Rep3",
                                "taf4b meiocyte RNA-seq Rep1", "taf4b meiocyte RNA-seq Rep2", "taf4b meiocyte RNA-seq Rep3",
                                "Bur0 meiocyte RNA-seq Rep1", "Bur0 meiocyte RNA-seq Rep2", "Bur0 meiocyte RNA-seq Rep3")
colnames(sampleDistMatrix) <- NULL
mycols <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(paste0(plotDir, "sample_distances_heatmap_rlog.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = mycols)
dev.off()

# Sample distances using the Poisson Distance
library(PoiClaClu)
poisd <- PoissonDistance(t(counts(dds)))
print(poisd)
#Value of alpha used to transform data:  0.515102
#This type of normalization was performed: mle
#Dissimilarity computed for  9  observations.

# Heatmap
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- c("Col meiocyte RNA-seq Rep1", "Col meiocyte RNA-seq Rep2", "Col meiocyte RNA-seq Rep3",
                                    "taf4b meiocyte RNA-seq Rep1", "taf4b meiocyte RNA-seq Rep2", "taf4b meiocyte RNA-seq Rep3",
                                    "Bur0 meiocyte RNA-seq Rep1", "Bur0 meiocyte RNA-seq Rep2", "Bur0 meiocyte RNA-seq Rep3")
colnames(samplePoisDistMatrix) <- NULL
pdf(paste0(plotDir, "sample_distances_heatmap_Poisson.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = mycols)
dev.off()

# PCA plots for visualising sample-to-sample distances
PCAplot_rlog <- DESeq2::plotPCA(rld, intgroup = c("condition", "sample")) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(PCAplot_rlog,
       file = paste0(plotDir, "PCAplot_rlog.pdf"), width = 20, height = 20, units = "cm")

PCAplot_rlog_data <- DESeq2::plotPCA(rld, intgroup = c("condition", "sample"), returnData = T)
PCAplot_rlog_data
# Obtain percentage variance explained by PC1 and PC2 for plotting using ggplot2
percentVar <- round(100 * attr(PCAplot_rlog_data, "percentVar"))

PCAggplot_rlog <- ggplot(PCAplot_rlog_data, aes(x = PC1, y = PC2,
                                                shape = condition,
                                                colour = sample)) +
                    geom_point(size = 3) +
                    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                    theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                    theme_classic() +
                    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                    coord_fixed()
ggsave(PCAggplot_rlog,
       file = paste0(plotDir, "PCAggplot_rlog.pdf"), height = 20, width = 20, units = "cm")

# MDS (multi-dimensional scaling) plots for visualising sample-to-sample distances
# rlog-transformed counts
mds <- as.data.frame(colData(rld)) %>%
         cbind(cmdscale(sampleDistMatrix))
MDSplot_rlog <- ggplot(mds, aes(x = `1`, y = `2`,
                                shape = condition,
                                colour = sample)) +
                  geom_point(size = 3) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(MDSplot_rlog,
       file = paste0(plotDir, "MDSplot_rlog.pdf"), height = 20, width = 20, units = "cm")

# Poisson distance
mdsPois <- as.data.frame(colData(dds)) %>%
             cbind(cmdscale(samplePoisDistMatrix))
MDSplot_Pois <- ggplot(mdsPois, aes(x = `1`, y = `2`,
                                shape = condition,
                                colour = sample)) +
                  geom_point(size = 3) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(MDSplot_Pois,
       file = paste0(plotDir, "MDSplot_Pois.pdf"), height = 20, width = 20, units = "cm")

# Gene clustering
# Clustering of 20 genes with greatest rlog-transformed variance across samples
library(genefilter)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = T), 20)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("sample", "condition")])
pdf(paste0(plotDir, "gene_clustering_rld_topVar20.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(mat, annotation_col = anno)
dev.off()



## Differential expression analysis

# The DESeq developers recommend that the DESeq() function is run to fit a model using samples from all groups,
# followed by use of the "contrast" argument of the results() function to make comparisons between two groups.
# "The model fit by DESeq estimates a single dispersion parameter for each gene, which defines how far we expect
# the observed count for a sample will be from the mean value from the model given its size factor and its condition group."
# (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
# However, given that Col replicates have high within-group variability, as revealed by the above sample-distance analyses,
# greater sensitivity would be achieved by subsetting the DESeqDataSet to only samples from two groups (Col and another)
# before fitting a model by running the DESeq() function.
# This would prevent inflation of the per-gene dispersion estimate.

# Re-import DESeqDataSet to redefine "dds"
dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, design = ~condition)
# Subset DESeqDataSet to only samples from the two conditions to be contrasted
dds_taf4bVCol <- dds[ , dds$condition == "taf4b" | dds$condition == "Col"]
# Remove empty levels
dds_taf4bVCol$condition <- droplevels(dds_taf4bVCol$condition)
dds_taf4bVCol$sample <- droplevels(dds_taf4bVCol$sample)

# Pre-filter the datasets
print(nrow(dds_taf4bVCol))
#[1] 27586
# Retain only rows that have more than a single count across all samples
dds_taf4bVCol <- dds_taf4bVCol[rowSums(counts(dds_taf4bVCol)) > 1,]
print(nrow(dds_taf4bVCol))
#[1] 22179

# Contrast: taf4b vs Col 
dds_taf4bVCol <- DESeq(dds_taf4bVCol)
res_taf4bVCol <- results(dds_taf4bVCol, contrast = c("condition", "taf4b", "Col"))
print(mcols(res_taf4bVCol, use.names = T))
#DataFrame with 6 rows and 2 columns
#                       type                                    description
#                <character>                                    <character>
#baseMean       intermediate      mean of normalized counts for all samples
#log2FoldChange      results log2 fold change (MLE): condition taf4b vs Col
#lfcSE               results         standard error: condition taf4b vs Col
#stat                results         Wald statistic: condition taf4b vs Col
#pvalue              results      Wald test p-value: condition taf4b vs Col
#padj                results                           BH adjusted p-values

print(summary(res_taf4bVCol))
#out of 22179 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 528, 2.4% 
#LFC < 0 (down)   : 1205, 5.4% 
#outliers [1]     : 19, 0.086% 
#low counts [2]   : 8170, 37% 
#(mean count < 3)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

print(table(res_taf4bVCol$padj < 0.1))
#FALSE  TRUE 
#12257  1733 

print(sum(!is.na(res_taf4bVCol$padj)))
#[1] 13990

res_taf4bVCol_.05 <- results(dds_taf4bVCol, contrast = c("condition", "taf4b", "Col"), alpha = 0.05)
print(summary(res_taf4bVCol_.05))
#out of 22179 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 410, 1.8% 
#LFC < 0 (down)   : 974, 4.4% 
#outliers [1]     : 19, 0.086% 
#low counts [2]   : 8600, 39% 
#(mean count < 4)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
print(table(res_taf4bVCol_.05$padj < 0.05))
#FALSE  TRUE 
#12176  1384

res_taf4bVCol_.01 <- results(dds_taf4bVCol, contrast = c("condition", "taf4b", "Col"), alpha = 0.01)
print(summary(res_taf4bVCol_.01))
#out of 22179 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)     : 256, 1.2% 
#LFC < 0 (down)   : 613, 2.8% 
#outliers [1]     : 19, 0.086% 
#low counts [2]   : 9890, 45% 
#(mean count < 6)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
print(table(res_taf4bVCol_.01$padj < 0.01))
#FALSE  TRUE 
#11401   869 

res_taf4bVCol_LFC2 <- results(dds_taf4bVCol, contrast = c("condition", "taf4b", "Col"), lfcThreshold = 2)
print(summary(res_taf4bVCol_LFC2))
#out of 22179 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 2, 0.009% 
#LFC < 0 (down)   : 21, 0.095% 
#outliers [1]     : 19, 0.086% 
#low counts [2]   : 11179, 50% 
#(mean count < 10)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

print(table(res_taf4bVCol_LFC2$padj < 0.1))
#FALSE  TRUE 
#10958    23

res_taf4bVCol_LFC1 <- results(dds_taf4bVCol, contrast = c("condition", "taf4b", "Col"), lfcThreshold = 1)
print(summary(res_taf4bVCol_LFC1))
#out of 22179 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 14, 0.063% 
#LFC < 0 (down)   : 218, 0.98% 
#outliers [1]     : 19, 0.086% 
#low counts [2]   : 9890, 45% 
#(mean count < 6)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
print(table(res_taf4bVCol_LFC1$padj < 0.1))
#FALSE  TRUE 
#12038   232

# Subset results table to significant (FDR < 0.1) genes
res_taf4bVColSig.1 <- subset(res_taf4bVCol, padj < 0.1)
print(dim(res_taf4bVColSig.1))
#[1] 1733    6

# Subset results table to significant (FDR < 0.1) down-regulated genes in taf4b
res_taf4bVColSig.1_downReg <- subset(res_taf4bVCol, padj < 0.1 & log2FoldChange < 0)
print(dim(res_taf4bVColSig.1_downReg))
#[1] 1205    6
# Sort by increasing log2 fold change estimate
# (Genes with strongest down-regulation)
res_taf4bVColSig.1_downRegSorted <- res_taf4bVColSig.1_downReg[order(res_taf4bVColSig.1_downReg$log2FoldChange),]
# Subset results table to significant (FDR < 0.1) up-regulated genes in taf4b
res_taf4bVColSig.1_upReg <- subset(res_taf4bVCol, padj < 0.1 & log2FoldChange > 0)
print(dim(res_taf4bVColSig.1_upReg))
#[1] 528   6
# Sort by decreasing log2 fold change estimate
# (Genes with strongest up-regulation)
res_taf4bVColSig.1_upRegSorted <- res_taf4bVColSig.1_upReg[order(res_taf4bVColSig.1_upReg$log2FoldChange, decreasing = T),]


# MA-plots
# These show the log2(fold changes) in gene expression values for a variable (e.g., treated vs untreated)
# over a the mean of normalized counts for all the samples in the DESeqDataSet

# DESeq2 provides for use of a "Bayesian procedure to moderate (or "shrink") log2 fold changes from genes
# with very low counts and highly variable counts"
# Before generating an MA-plot, the lfcShrink() function should be used to moderate the log2(fold changes) for the contrast
res_taf4bVCol_lfcShrink <- lfcShrink(dds_taf4bVCol, contrast = c("condition", "taf4b", "Col"), res = res_taf4bVCol)
pdf(paste0(plotDir, "MAplot_res_taf4bVCol_lfcShrink.pdf"), height = 10, width = 6)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plotMA(res_taf4bVCol_lfcShrink, ylim = c(-10, 10),
       xlab = "Mean of normalized counts",
       ylab = "log2(fold change)")
plotMA(res_taf4bVCol_lfcShrink, ylim = c(-5, 5),
       xlab = "Mean of normalized counts",
       ylab = "log2(fold change)")
dev.off()

# The following MA-plot is generated without moderating the noisy log2(fold changes)
pdf(paste0(plotDir, "MAplot_res_taf4bVCol.pdf"), height = 10, width = 6)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plotMA(res_taf4bVCol, ylim = c(-10, 10),
       xlab = "Mean of normalized counts",
       ylab = "log2(fold change)")
plotMA(res_taf4bVCol, ylim = c(-5, 5),
       xlab = "Mean of normalized counts",
       ylab = "log2(fold change)")
dev.off()

# Gene selection/subsets
topGene <- rownames(res_taf4bVCol_lfcShrink)[which.min(res_taf4bVCol_lfcShrink$padj)]
topChrGene <- rownames(res_taf4bVCol_lfcShrink)[which.min(res_taf4bVCol_lfcShrink$padj[!grepl("ATM", rownames(res_taf4bVCol_lfcShrink)) &
                                                                                       !grepl("ATC", rownames(res_taf4bVCol_lfcShrink))])]

res_taf4bVCol_lfcShrink_ATM <- res_taf4bVCol_lfcShrink[grep("ATM", rownames(res_taf4bVCol_lfcShrink)), ]
res_taf4bVCol_lfcShrink_ATC <- res_taf4bVCol_lfcShrink[grep("ATC", rownames(res_taf4bVCol_lfcShrink)), ]
res_taf4bVCol_lfcShrink_Chr <- res_taf4bVCol_lfcShrink[!grepl("ATM", rownames(res_taf4bVCol_lfcShrink)) &
                                                       !grepl("ATC", rownames(res_taf4bVCol_lfcShrink)), ]
res_taf4bVCol_lfcShrink_ChrDownReg <- res_taf4bVCol_lfcShrink_Chr[res_taf4bVCol_lfcShrink_Chr$log2FoldChange < 0, ]
res_taf4bVCol_lfcShrink_ChrUpReg <- res_taf4bVCol_lfcShrink_Chr[res_taf4bVCol_lfcShrink_Chr$log2FoldChange > 0, ]

topChrDownRegGene <- rownames(res_taf4bVCol_lfcShrink_ChrDownReg)[which.min(res_taf4bVCol_lfcShrink_ChrDownReg$padj)]
topChrUpRegGene <- rownames(res_taf4bVCol_lfcShrink_ChrUpReg)[which.min(res_taf4bVCol_lfcShrink_ChrUpReg$padj)]

res_taf4bVCol_lfcShrink[res_taf4bVCol_lfcShrink$log2FoldChange > 5,]
#log2 fold change (MAP): condition taf4b vs Col 
#Wald test p-value: condition taf4b vs Col 
#DataFrame with 2 rows and 5 columns
#           baseMean log2FoldChange      stat       pvalue         padj
#          <numeric>      <numeric> <numeric>    <numeric>    <numeric>
#ATMG00370  75.72792       7.436843  7.988849 1.362045e-15 6.351668e-12
#ATMG00910 861.09902       9.980950 11.199823 4.085493e-29 5.715604e-25

res_taf4bVCol_lfcShrink[res_taf4bVCol_lfcShrink$baseMean > 1e+05,]
#log2 fold change (MAP): condition taf4b vs Col 
#Wald test p-value: condition taf4b vs Col 
#DataFrame with 3 rows and 5 columns
#           baseMean log2FoldChange      stat      pvalue       padj
#          <numeric>      <numeric> <numeric>   <numeric>  <numeric>
#AT2G30230  357130.2       1.276435  3.362729 0.000771761 0.01207711
#AT5G10100  526343.7       1.138818  3.055989 0.002243193 0.02755247
#ATMG00030  240353.2       1.130745  2.749023 0.005977323 0.05832868

pdf(paste0(plotDir, "MAplot_res_taf4bVCol_lfcShrink_topChrGenes_UpReg_DownReg.pdf"), height = 10, width = 6)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 2.1, 2.1))

plotMA(res_taf4bVCol_lfcShrink, ylim = c(-10, 10),
       xlab = "Mean of normalized counts",
       ylab = "log2(fold change)")
with(res_taf4bVCol_lfcShrink_ChrDownReg[topChrDownRegGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topChrDownRegGene, pos = 2, col = "dodgerblue")
})
with(res_taf4bVCol_lfcShrink_ChrUpReg[topChrUpRegGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topChrUpRegGene, pos = 2, col = "dodgerblue")
})
with(res_taf4bVCol_lfcShrink[c("ATMG00370", "ATMG00910"), ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, c("ATMG00370", "ATMG00910"), pos = 2, col = "dodgerblue")
})
with(res_taf4bVCol_lfcShrink[c("AT2G30230", "AT5G10100", "ATMG00030"), ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, c("AT2G30230", "AT5G10100", "ATMG00030"), pos = c(1, 3, 2), col = "dodgerblue")
})

plotMA(res_taf4bVCol_lfcShrink, ylim = c(-5, 5),
       xlab = "Mean of normalized counts",
       ylab = "log2(fold change)")
with(res_taf4bVCol_lfcShrink_ChrDownReg[topChrDownRegGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topChrDownRegGene, pos = 2, col = "dodgerblue")
})
with(res_taf4bVCol_lfcShrink_ChrUpReg[topChrUpRegGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topChrUpRegGene, pos = 2, col = "dodgerblue")
})
with(res_taf4bVCol_lfcShrink[c("ATMG00370", "ATMG00910"), ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, c("ATMG00370", "ATMG00910"), pos = 2, col = "dodgerblue")
})
with(res_taf4bVCol_lfcShrink[c("AT2G30230", "AT5G10100", "ATMG00030"), ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, c("AT2G30230", "AT5G10100", "ATMG00030"), pos = c(1, 3, 2), col = "dodgerblue")
})
dev.off()


# Histograms of unadjusted and FDR-adjusted p-values
pdf(paste0(plotDir, "hist_res_taf4bVCol_pvals_unadjusted.pdf"))
hist(res_taf4bVCol$pvalue[res_taf4bVCol$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
dev.off()

pdf(paste0(plotDir, "hist_res_taf4bVCol_pvals_BH_FDR_adjusted.pdf"))
hist(res_taf4bVCol$padj[res_taf4bVCol$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
dev.off()


# Independent filtering
# Create bins
qs <- c(0, quantile(res_taf4bVCol$baseMean[res_taf4bVCol$baseMean > 0], 0:6/6))
# Bin genes by base mean (mean normalized count) using cut
bins <- cut(res_taf4bVCol$baseMean, qs)
# Rename the levels of the bins using the middle point
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(res_taf4bVCol$pvalue, bins, function(p)
                 mean(p < 0.05, na.rm = TRUE))
pdf(paste0(plotDir, "barplot_res_taf4bVCol_fraction_small_pval_genes.pdf"))
barplot(fractionSig, xlab = "Mean normalized count",
                     ylab = "Fraction of small P values")
dev.off()


## Exporting results

# Subset results table to significant (FDR < 0.1) down-regulated genes on chromosomes in taf4b
res_taf4bVCol_lfcShrink_Chr_Sig.1_downReg <- subset(res_taf4bVCol_lfcShrink_Chr, padj < 0.1 & log2FoldChange < 0)
print(dim(res_taf4bVCol_lfcShrink_Chr_Sig.1_downReg))
#[1] 1203    5
# Sort by increasing log2 fold change estimate
# (Genes with strongest down-regulation)
res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSorted <- res_taf4bVCol_lfcShrink_Chr_Sig.1_downReg[order(res_taf4bVCol_lfcShrink_Chr_Sig.1_downReg$log2FoldChange),]
# Subset results table to significant (FDR < 0.1) up-regulated genes on chromosomes in taf4b
res_taf4bVCol_lfcShrink_Chr_Sig.1_upReg <- subset(res_taf4bVCol_lfcShrink_Chr, padj < 0.1 & log2FoldChange > 0)
print(dim(res_taf4bVCol_lfcShrink_Chr_Sig.1_upReg))
#[1] 355   5
# Sort by decreasing log2(fold change) estimate
# (Genes with strongest up-regulation) 
res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSorted <- res_taf4bVCol_lfcShrink_Chr_Sig.1_upReg[order(res_taf4bVCol_lfcShrink_Chr_Sig.1_upReg$log2FoldChange, decreasing = T),]

# Convert DataFrame objects (IRanges package) to data.frame objects that can be processed by write.table
res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSortedDF <- as.data.frame(res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSorted)
res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSortedDF <- as.data.frame(res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSorted)
write.table(res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSortedDF,
            file = paste0(outDir, "res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSortedDF.txt"))
write.table(res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSortedDF,
            file = paste0(outDir, "res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSortedDF.txt"))

res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegGeneIDs <- rownames(res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSortedDF)
res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegGeneIDs <- rownames(res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSortedDF)
write.table(res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegGeneIDs,
            file = paste0(outDir, "res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegGeneIDs.txt"))
write.table(res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegGeneIDs,
            file = paste0(outDir, "res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegGeneIDs.txt"))

# Use ReportingTools bioconductor package to "automatically generate dynamic HTML documents, including links to external databases using gene identifiers and boxplots summarizing the normalized counts across groups. See the ReportingTools vignettes for full details. The simplest version of creating a dynamic ReportingTools report is performed with the following code:"
library(ReportingTools)
setwd("/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_01")
htmlRepDownReg <- HTMLReport(shortName = "taf4bVCol_downReg_report", title = "res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSortedDF report",
                      reportDirectory = "./reports")
publish(res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegSortedDF, htmlRepDownReg)
urlDownReg <- finish(htmlRepDownReg)
browseURL(urlDownReg)

htmlRepUpReg <- HTMLReport(shortName = "taf4bVCol_upReg_report", title = "res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSortedDF report",
                           reportDirectory = "./reports")
publish(res_taf4bVCol_lfcShrink_Chr_Sig.1_upRegSortedDF, htmlRepUpReg)
urlUpReg <- finish(htmlRepUpReg)
browseURL(urlUpReg)

