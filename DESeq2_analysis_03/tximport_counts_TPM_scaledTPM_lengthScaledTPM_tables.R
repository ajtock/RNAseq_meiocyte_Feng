#!/applications/R/R-3.4.0/bin/Rscript

# R version 3.4.0 (invoke this version by specifying path at command line:
# /applications/R/R-3.4.0/bin/R , or if running as a script:
# /applications/R/R-3.4.0/bin/Rscript)

# Code below based on:
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# Usage:
# tximport_counts_TPM_scaledTPM_lengthScaledTPM_tables.R ../samples.txt DESeq2_analysis_03 

#sampleFile <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/samples.txt"
#analysisDir <- "DESeq2_analysis_03"

args <- commandArgs(trailingOnly = TRUE)
sampleFile <- args[1]
analysisDir <- args[2]

outDir <- paste0(dirname(sampleFile), "/", analysisDir, "/")

## Import transcript abundance datasets with tximport package

library(tximport)
print(packageVersion("tximport"))
#[1] ‘1.4.0’

# Read in table of sample IDs that will be used to specify paths to count files
samples <- read.table(sampleFile, header = T)
print(samples)

# Specify paths to count files
files <- file.path(dirname(sampleFile),
                   samples$genotype, samples$directory, samples$sample,
                   "quant.sf")
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

txi_scaled <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = "scaledTPM")
txi_lengthScaled <- tximport(files, type = "salmon", tx2gene = tx2gene,
                             countsFromAbundance = "lengthScaledTPM")

txi_counts <- txi$counts
txi_tpm <- txi$abundance
colnames(txi_counts) <- c("Col_Rep1", "Col_Rep2", "Col_Rep3",
                          "taf4b_Rep1", "taf4b_Rep2", "taf4b_Rep3",
                          "Bur0_Rep1", "Bur0_Rep2", "Bur0_Rep3")
colnames(txi_tpm) <- c("Col_Rep1", "Col_Rep2", "Col_Rep3",
                       "taf4b_Rep1", "taf4b_Rep2", "taf4b_Rep3",
                       "Bur0_Rep1", "Bur0_Rep2", "Bur0_Rep3")
write.table(as.data.frame(txi_counts),
            file = paste0(outDir, "salmon_counts_RNAseq_Col_taf4b_Bur0.txt"),
            sep = "\t", quote = F)
write.table(as.data.frame(txi_tpm),
            file = paste0(outDir, "salmon_TPM_RNAseq_Col_taf4b_Bur0.txt"),
            sep = "\t", quote = F)


