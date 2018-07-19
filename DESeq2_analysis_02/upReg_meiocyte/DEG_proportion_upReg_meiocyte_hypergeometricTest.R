#!/applications/R/R-3.4.0/bin/Rscript
### Perform hypergeometric test to determine whether a significant proportion of DEGs
### are up-regulated in meiocytes
# P-value is the probability of getting >= length(geneIDintersection) genes
# in a sample size of length(geneIDs) from a total gene set consisting of 
# length(upRegMeio_geneIDs) + (length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs))

# R version 3.4.0 (invoke this version by specifying path at command line:
# /applications/R/R-3.4.0/bin/R , or if running as a script:
# /applications/R/R-3.4.0/bin/Rscript)

# Usage: DEG_proportion_upReg_meiocyte_hypergeometricTest.R ../taf4bVCol_Bur0VCol_.1_lfcShrink_Chr_Sig.1_downRegGeneIDs.txt upReg_meiocyte 0.01 100000 100,000

library(methods)

#geneIDsFile <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/taf4bVCol_Bur0VCol_.1_lfcShrink_Chr_Sig.1_downRegGeneIDs.txt"
#analysisDir <- "upReg_meiocyte"
#meioFDR <- 0.01

#samples <- 100000

args <- commandArgs(trailingOnly = TRUE)
geneIDsFile <- args[1]
analysisDir <- args[2]
meioFDR <- args[3]
# Number of randomisations to perform
samples <- as.numeric(args[4])
samplesChar <- as.character(args[5])

outDir <- paste0(dirname(geneIDsFile), "/", analysisDir, "/")
plotDir <- paste0(dirname(geneIDsFile), "/", analysisDir, "/plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# TAIR10 geneIDs (on chromosomes only)
load(paste0(dirname(geneIDsFile), "/", "TAIR10_geneIDsChr.RData"))
print(length(TAIR10_geneIDsChr))
#[1] 27389

# Genes up-regulated in meiocytes (on chromosomes only)
meioDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/DESeq2_analysis_01/"
upRegMeio_geneIDs <- as.character(read.table(paste0(meioDir,
                                             "res_meiocyteVleaf_", as.character(meioFDR),
                                             "_lfcShrink_Chr_Sig", as.character(meioFDR),
                                             "_upRegGeneIDs.txt"))$x)
print(length(upRegMeio_geneIDs))
#[1] 4528

# Query gene IDs (on chromosomes only)
geneIDs <- as.character(read.table(geneIDsFile)$x)

m <- length(upRegMeio_geneIDs)
print(m)
#[1] 4528
n <- length(TAIR10_geneIDsChr)-m
print(n)
#[1] 22861
k <- length(geneIDs)
print(k)
#[1] 1804

# Obtain intersection of gene IDs in query gene set and
# gene IDs of genes up-regulated in meiocytes 
geneIDintersection <- intersect(geneIDs, upRegMeio_geneIDs)
print(length(geneIDintersection))
#[1] 850
# WARNING: ASSUMES INPUT FILE HAS 3-LETTER EXTENSION
baseName <- basename(geneIDsFile)
baseName <- substr(baseName, 1, nchar(baseName)-4)
write.table(geneIDintersection,
            file = paste0(outDir,
                          baseName,
                          "_res_meiocyteVleaf_", as.character(meioFDR),
                          "_lfcShrink_Chr_Sig", as.character(meioFDR),
                          "_upRegGeneIDs.txt"))

# Calculate proportion of query genes up-regulated in meiocytes
proportion_upRegMeio <- length(geneIDintersection)/length(geneIDs)

# P-value is the probability of drawing >= length(geneIDintersection) [x] genes
# in a sample size of length(geneIDs) [k] from a total gene set consisting of 
# length(upRegMeio_geneIDs) [m] + (length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by sum of the probabilities of drawing
# length(geneIDintersection) to length(geneIDs)
Pval <- sum(dhyper(x = length(geneIDintersection):length(geneIDs),
                   m = length(upRegMeio_geneIDs),
                   n = length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs),
                   k = length(geneIDs)))
print(Pval)
#[1] 2.777458e-219
print(sum(dhyper(x = length(geneIDintersection):k, m = m, n = n, k = k)))
#[1] 2.777458e-219
# Or by 1 minus the sum of probabilities of drawing 0:(length(geneIDintersection)-1)
print(1 - sum(dhyper(x = 0:(length(geneIDintersection)-1),
                     m = length(upRegMeio_geneIDs),
                     n = length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs),
                     k = length(geneIDs))))
#[1] -2.220446e-16
print(1 - sum(dhyper(x = 0:(length(geneIDintersection)-1), m = m, n = n, k = k)))
#[1] -2.220446e-16

# phyper(q, m, n, k) gives the probability of drawing <= q,
# so phyper(q, m, n, k) is the same as sum(dhyper(0:x, m, n, k)).
# phyper(q, m, n, k, lower.tail = FALSE) is the same as 
# 1 - phyper(q, m, n, k), and so is the probablity of drawing >= q+1:
print(phyper(q = length(geneIDintersection)-1,
             m = length(upRegMeio_geneIDs),
             n = length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs),
             k = length(geneIDs),
             lower.tail = FALSE))
#[1] 2.777458e-219
print(1 - phyper(q = length(geneIDintersection)-1,
                 m = length(upRegMeio_geneIDs),
                 n = length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs),
                 k = length(geneIDs)))
#[1] 0

# Sample without replacement
hgDist <- rhyper(nn = samples,
                 m = length(upRegMeio_geneIDs),
                 n = length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs),
                 k = length(geneIDs))
random_proportions_upRegMeio <- hgDist/length(geneIDs)

setClass("hypergeomTest",
         representation(Pval = "numeric",
                        proportion_upRegMeio = "numeric",
                        random_proportions_upRegMeio = "numeric",
                        geneIDintersectionLength = "numeric",
                        hypergeometricDistribution = "numeric"))
hgTestResults <- new("hypergeomTest",
                     Pval = Pval,
                     proportion_upRegMeio = proportion_upRegMeio,
                     random_proportions_upRegMeio = random_proportions_upRegMeio, 
                     geneIDintersectionLength = length(geneIDintersection),
                     hypergeometricDistribution = hgDist)

basenameFile <- basename(geneIDsFile)
len <- nchar(basenameFile)
# WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
basenameFile <- substr(basenameFile, 1, len-4)
save(hgTestResults,
     file = paste0(outDir, basenameFile,
                   "_hypergeometricTestResults_upRegMeioSig", as.character(meioFDR), ".RData"))

library(plotrix)
# Generate histogram
pdf(paste0(plotDir, "hist_", basenameFile,
           "_hypergeometricTestResults_upRegMeioSig", as.character(meioFDR), ".pdf"),
           height = 4, width = 5)
## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
#options(scipen = 100)
# Calculate max density
maxDensityPlus <- max(density(hgTestResults@random_proportions_upRegMeio)$y)*1.2
alpha0.05 <- quantile(hgTestResults@random_proportions_upRegMeio, 0.95)[[1]]
par(mar = c(5.1, 4.1, 4.1, 2.1))
hist(hgTestResults@random_proportions_upRegMeio,
     freq = FALSE,
     col = "grey70",
     border = NA,
     lwd = 2,
     xlim = c(pmax(0, min(hgTestResults@random_proportions_upRegMeio)-.1),
              pmax(hgTestResults@proportion_upRegMeio+.1, alpha0.05+.1)),
     ylim = c(0,
              maxDensityPlus),
     xlab = paste0("Proportion of genes \nup-regulated in meiocytes (FDR < ",
                   as.character(meioFDR), ")"),
     ylab = "Density",
     main = "",
     cex.lab = 1, cex.axis = 1)
titleText <- list(bquote(.(basenameFile)),
                  bquote(italic("P")~" = "~.(hgTestResults@Pval)),
                  bquote("Samples (hypergeometric distribution) = "~.(samplesChar)))
mtext(do.call(expression, titleText), side = 3, line = 2:0, cex = 1)
lines(density(hgTestResults@random_proportions_upRegMeio), lwd = 1.5)
ablineclip(v = mean(hgTestResults@random_proportions_upRegMeio),
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
ablineclip(v = hgTestResults@proportion_upRegMeio,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
ablineclip(v = alpha0.05,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
text(x = c(pmax(0.05, min(hgTestResults@random_proportions_upRegMeio)-.05),
           mean(hgTestResults@random_proportions_upRegMeio),
           hgTestResults@proportion_upRegMeio,
           alpha0.05),
     y = c(maxDensityPlus*.95,
           maxDensityPlus,
           maxDensityPlus,
           maxDensityPlus*.95),
     labels = c("Simulated",
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

