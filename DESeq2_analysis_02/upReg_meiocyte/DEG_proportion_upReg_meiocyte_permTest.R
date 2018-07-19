#!/applications/R/R-3.4.0/bin/Rscript
### Perform permutation test to determine whether a significant proportion of DEGs
### are up-regulated in meiocytes

# R version 3.4.0 (invoke this version by specifying path at command line:
# /applications/R/R-3.4.0/bin/R , or if running as a script:
# /applications/R/R-3.4.0/bin/Rscript)

# Usage: DEG_proportion_upReg_meiocyte_permTest.R ../taf4bVCol_Bur0VCol_.1_lfcShrink_Chr_Sig.1_downRegGeneIDs.txt upReg_meiocyte 0.01 100000 0.00001

library(methods)
library(parallel)

#geneIDsFile <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/taf4bVCol_Bur0VCol_.1_lfcShrink_Chr_Sig.1_downRegGeneIDs.txt"
#analysisDir <- "upReg_meiocyte"
#meioFDR <- 0.01
#perms <- 10000
#minPval <- 0.0001

args <- commandArgs(trailingOnly = TRUE)
geneIDsFile <- args[1]
analysisDir <- args[2]
meioFDR <- args[3]
# Number of permutations (randomisations) to perform
perms <- as.numeric(args[4])
# Corresponding minimum P-value
minPval <- args[5]

outDir <- paste0(dirname(geneIDsFile), "/", analysisDir, "/")
plotDir <- paste0(dirname(geneIDsFile), "/", analysisDir, "/plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# TAIR10 geneIDs (on chromosomes only)
load(paste0(dirname(geneIDsFile), "/", "TAIR10_geneIDsChr.RData"))
print(length(TAIR10_geneIDsChr))
#[1] 27389

# Query gene IDs (on chromosomes only)
geneIDs <- as.character(read.table(geneIDsFile)$x)

# Genes up-regulated in meiocytes (on chromosomes only)
meioDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/DESeq2_analysis_01/"
upRegMeio_geneIDs <- as.character(read.table(paste0(meioDir,
                                             "res_meiocyteVleaf_", as.character(meioFDR),
                                             "_lfcShrink_Chr_Sig", as.character(meioFDR),
                                             "_upRegGeneIDs.txt"))$x)
print(length(upRegMeio_geneIDs))
#[1] 4528

# Obtain intersection of gene IDs in query gene set and
# gene IDs of genes up-regulated in meiocytes 
geneIDintersection <- intersect(geneIDs, upRegMeio_geneIDs)
# Calculate proportion of query genes up-regulated in meiocytes
proportion_upRegMeio <- length(geneIDintersection)/length(geneIDs)

# Obtain random samples of TAIR10 gene IDs (sample size = length(geneIDs))
TAIR10_geneIDsChr_random <- mclapply(seq(1, perms), function(x) {
  sample(TAIR10_geneIDsChr, size = length(geneIDs), replace = FALSE)
}, mc.cores = 48)

# Obtain intersections of gene IDs in random samples and
# gene IDs of genes up-regulated in meiocytes
randomIDintersections <- mclapply(seq_along(TAIR10_geneIDsChr_random), function(x) {
  intersect(TAIR10_geneIDsChr_random[[x]], upRegMeio_geneIDs)
}, mc.cores = 48)

# Obtain lengths of intersections of gene IDs in random samples and
# gene IDs of genes up-regulated in meiocytes
randomIDintersectionLengths <- mclapply(seq_along(randomIDintersections), function(x) {
  length(randomIDintersections[[x]])
}, mc.cores = 48)

# Calculate proportions of random genes up-regulated in meiocytes
randomIDproportions_upRegMeio <- mclapply(seq_along(randomIDintersections), function(x) {
  length(randomIDintersections[[x]])/length(TAIR10_geneIDsChr_random[[x]])
}, mc.cores = 48)

# Test whether proprtions of random genes up-regulated in meiocytes are
# smaller than the proportion of query genes up-regulated in meiocytes
randomIDproportions_lessThan_proportion_upRegMeio_Bool <- mclapply(seq_along(randomIDproportions_upRegMeio),
  function(x) {
    randomIDproportions_upRegMeio[[x]] < proportion_upRegMeio
}, mc.cores = 48)

# Test whether numbers of random genes up-regulated in meiocytes are
# smaller than the number of query genes up-regulated in meiocytes
randomIDintersectionLengths_lessThan_geneIDintersectionLength_upRegMeio_Bool <- mclapply(seq_along(randomIDintersectionLengths),
  function(x) {
    randomIDintersectionLengths[[x]] < length(geneIDintersection)
}, mc.cores = 48)

# Calculate P-values
PvalProp <- 1-(sum(unlist(randomIDproportions_lessThan_proportion_upRegMeio_Bool))/perms)
PvalNum <- 1-(sum(unlist(randomIDintersectionLengths_lessThan_geneIDintersectionLength_upRegMeio_Bool))/perms)

if(PvalProp == 0) {
  PvalProp <- as.numeric(minPval)
}
if(PvalNum == 0) {
  PvalNum <- as.numeric(minPval)
}

setClass("permTest",
         representation(PvalProp = "numeric",
                        PvalNum = "numeric",
                        proportion_upRegMeio = "numeric",
                        randomIDproportions_upRegMeio = "numeric",
                        geneIDintersectionLength = "numeric",
                        randomIDintersectionLengths = "numeric"))
permTestResults <- new("permTest",
                       PvalProp = PvalProp, 
                       PvalNum = PvalNum,
                       proportion_upRegMeio = proportion_upRegMeio,
                       randomIDproportions_upRegMeio = unlist(randomIDproportions_upRegMeio),
                       geneIDintersectionLength = length(geneIDintersection),
                       randomIDintersectionLengths = unlist(randomIDintersectionLengths))

basenameFile <- basename(geneIDsFile)
len <- nchar(basenameFile)
# WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
basenameFile <- substr(basenameFile, 1, len-4)
save(permTestResults,
     file = paste0(outDir, basenameFile,
                   "_permTestResults_upRegMeioSig", as.character(meioFDR), ".RData"))

library(plotrix)
# Generate histogram
pdf(paste0(plotDir, "hist_", basenameFile,
           "_permTestResults_upRegMeioSig", as.character(meioFDR), ".pdf"),
           height = 4, width = 5)
# Disable scientific notation (e.g., 0.0001 rather than 1e-04)
options(scipen = 100)
# Calculate max density
maxDensityPlus <- max(density(permTestResults@randomIDproportions_upRegMeio)$y)*1.2
alpha0.05 <- quantile(permTestResults@randomIDproportions_upRegMeio, 0.95)[[1]]
#par(mfrow = c(2, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
hist(permTestResults@randomIDproportions_upRegMeio,
     freq = FALSE,
     col = "grey70",
     border = NA,
     lwd = 2,
     xlim = c(pmax(0, min(permTestResults@randomIDproportions_upRegMeio)-.1),
              pmax(permTestResults@proportion_upRegMeio+.1, alpha0.05+.1)),
     ylim = c(0,
              maxDensityPlus),
     xlab = paste0("Proportion of genes \nup-regulated in meiocytes (FDR < ",
                   as.character(meioFDR), ")"),
     ylab = "Density",
     main = bquote(atop(italic("P")~" = "~.(as.character(round(permTestResults@PvalProp,
                                                               digits = 5))),
                   "Permutations = "~.(as.character(perms)))),
     cex.main = 1, cex.lab = 1, cex.axis = 1)
lines(density(permTestResults@randomIDproportions_upRegMeio), lwd = 1.5)
ablineclip(v = mean(permTestResults@randomIDproportions_upRegMeio),
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
ablineclip(v = permTestResults@proportion_upRegMeio,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
ablineclip(v = alpha0.05,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
text(x = c(pmax(0.05, min(permTestResults@randomIDproportions_upRegMeio)-.05),
           mean(permTestResults@randomIDproportions_upRegMeio),
           permTestResults@proportion_upRegMeio,
           alpha0.05),
     y = c(maxDensityPlus*.95,
           maxDensityPlus,
           maxDensityPlus,
           maxDensityPlus*.95),
     labels = c("Permuted genes",
                "Expected",
                "Observed",
                expression(alpha~" = 0.05")),
     col = c("grey70",
             "black",
             "forestgreen",
             "red"),
     cex = 0.7)
box(lwd = 2)
dev.off()

