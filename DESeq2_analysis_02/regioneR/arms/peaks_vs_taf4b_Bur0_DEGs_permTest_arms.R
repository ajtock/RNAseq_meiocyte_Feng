#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if features of
# interest (e.g., peaks) overlap DEGs more than expected by chance

# Usage:
# peaks_vs_taf4b_Bur0_DEGs_permTest_arms.R taf4bVCol_Bur0VCol_.01_lfcShrink_Chr_Sig.01_downRegGeneIDs.txt 

#DEGsFile <- "taf4bVCol_Bur0VCol_.01_lfcShrink_Chr_Sig.01_downRegGeneIDs.txt"

args <- commandArgs(trailingOnly = TRUE)
DEGsFile <- args[1]

DEGsFileBasename <- basename(DEGsFile)
len <- nchar(DEGsFileBasename)
DEGsFileBasename <- substr(DEGsFileBasename, 1, len-11)

inDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/"
outDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/regioneR/arms/"
plotDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/regioneR/arms/plots/FeaturesVsDEGs/"

library(regioneR)

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))

# Load DEG gene IDs as a vector
DEGs <- as.character(read.table(paste0(inDir, DEGsFile))$x)

# Import table of representative genes
genes <- read.table(file = "/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                    header = T)
# Convert to GRanges object
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start, end = genes$end),
                   strand = genes$strand,
                   geneID = substr(genes$gene_model, 1, 9))
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print(length(genesGR))
#[1] 27204

# Retain genes whose IDs match those in DEGs
DEGsGR <- genesGR[genesGR$geneID %in% DEGs]
print(DEGsGR)
# Remove DEGs located within pericentromeric regions
maskDEGsOverlaps <- findOverlaps(mask, DEGsGR,
                                 ignore.strand = TRUE, select = "all") 
DEGsGR <- DEGsGR[-subjectHits(maskDEGsOverlaps)]
print(DEGsGR)

# Import peaks on chromosome arms
load("/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.2_q0.2/idr/sigValRank_TreadsNormCreads/armrangerPeaksGR_REC8_HA_Rep1_REC8_MYC_Rep2_idr0.05.RData")
REC8GR <- armrangerPeaksGR
#seqlevels(REC8GR) <- sub("Chr", "", seqlevels(REC8GR))
maskPeaksOverlaps <- findOverlaps(mask, REC8GR, ignore.strand = TRUE, select = "all")
print(maskPeaksOverlaps)
#REC8GR <- REC8GR[-subjectHits(maskPeaksOverlaps)]
#REC8GR <- REC8GR[width(REC8GR) <= 2000]

load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
nucleRnucsGR <- armPeaksSH99GRmerge
armPeaksSH99GRmerge <- NULL
print(length(nucleRnucsGR))
#[1] 37529

load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/armrangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
rangernucsGR <- armrangerPeaksGRmerge
armrangerPeaksGRmerge <- NULL
rangernucsGR <- rangernucsGR[width(rangernucsGR) <= 500]
print(length(rangernucsGR))
#[1] 42235

load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
SPO11GR <- armrangerPeaksGRmerge
armrangerPeaksGRmerge <- NULL
print(length(SPO11GR))
#[1] 4340

load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
H3K4me3GR <- armrangerPeaksGRmerge
armrangerPeaksGRmerge <- NULL
print(length(H3K4me3GR))
#[1] 12716

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/armrangerPeaksGRmerge_WT_H3K9me2_ChIP_p0.05_q0.05_noMinWidth.RData")
H3K9me2GR <- armrangerPeaksGRmerge
armrangerPeaksGRmerge <- NULL
print(length(H3K9me2GR))
#[1] 12217

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/BCP/p0.05/armbcpPeaksGRmerge_WT_H3K9me2_ChIP_p0.05_noMinWidth.RData")
H3K9me2GRbcp <- armbcpPeaksGRmerge
armbcpPeaksGRmerge <- NULL
print(length(H3K9me2GRbcp))
#[1] 134

load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
COsGR <- COsGRcoords
print(length(COsGR))
#[1] 3320
# Remove COs located within pericentromeric regions
maskCOsOverlaps <- findOverlaps(mask, COsGR, ignore.strand = TRUE, select = "all")
COsGR <- COsGR[-subjectHits(maskCOsOverlaps)]
print(length(COsGR))
#[1] 2452

otherNames <- c("REC8GR", "nucleRnucsGR", "rangernucsGR", "SPO11GR",
                "H3K4me3GR", "H3K9me2GR", "H3K9me2GRbcp", "COsGR")
grl <- GRangesList("REC8GR" = REC8GR, "nucleRnucsGR" = nucleRnucsGR, "rangernucsGR" = rangernucsGR, "SPO11GR" = SPO11GR,
                   "H3K4me3GR" = H3K4me3GR, "H3K9me2GR" = H3K9me2GR, "H3K9me2GRbcp" = H3K9me2GRbcp, "COsGR" = COsGR)

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptDEGsOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = grl[[x]], B = DEGsGR, genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 47)
})
 
for(i in 1:length(ptDEGsOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptDEGsOtherPerChrom[[i]])
}
save(ptDEGsOtherPerChrom,
     file = paste0(outDir,
                   "permTest_REC8_nucleRnucs_rangernucs_SPO11_H3K4me3_H3K9me2_COs_noMinWidth_vs_",
                   DEGsFileBasename,
                   ".RData"))

# summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptDEGsOtherPerChrom)) {
  noOfFeaturesi <- print(length(grl[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptDEGsOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptDEGsOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptDEGsOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptDEGsOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptDEGsOtherPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptDEGsOtherPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_REC8_nucleRnucs_rangernucs_SPO11_H3K4me3_H3K9me2_COs_noMinWidth_vs_",
                          DEGsFileBasename,
                          "_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptDEGsOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_nperm10000_", DEGsFileBasename,"_arms_perChrom.pdf"), width = 10, height = 7) 
  plot(ptDEGsOtherPerChrom[[i]], main = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms"), xlab = "Number of overlaps", ylab = "Relative frequency")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between DEGs and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptDEGsOtherPerChrom[[i]], A = grl[[i]], B = DEGsGR,
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptDEGsOtherPerChrom[[i]], A = grl[[i]], B = DEGsGR,
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptDEGsOtherPerChrom[[i]], A = grl[[i]], B = DEGsGR,
                           window = 10*mean(width(grl[[i]])), step = mean(width(grl[[i]]))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(grl[[i]])))/1000))
  step <- as.character(round(mean(width(grl[[i]]))/2))
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_nperm10000_", DEGsFileBasename, "_arms_w1kb_s50bp_w10kb_s500bp_w", win, "kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms (1-kb shift)"))
  plot(lz_10kb, main = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms (10-kb shift)"))
  plot(lz_custom, main = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(otherNames[i], " vs ", DEGsFileBasename, " on chromosome arms (~", win, "-kb shift)"))
  dev.off()
}



