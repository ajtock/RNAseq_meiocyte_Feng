#!/applications/R/R-3.4.0/bin/Rscript

# Create heatmap of differentially expressed genes that are up-regulated in meiocytes,
# ordered by decreasing log2(fold change) and coloured according to TPM

# Usage:
# Rscript DEG_heatmap.R taf4b 0.1

library(ComplexHeatmap)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)
genotype1 <- args[1]
FDR <- args[2]

genotypes <- c(genotype1, genotype2)
myColours <- c(colour1, colour2)

inDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/"
res_taf4bVCol_0.1_lfcShrink_Chr_Sig0.1_downRegSortedDF.txt
upRegList <- lapply(seq_along(genotypes), function(x) {
  read.table(paste0(inDir,
                    "res_", genotypes[x], "VCol_",
                                 as.character(FDR), "_lfcShrink_Chr_Sig",
                                 as.character(FDR), "_downRegSortedDF.txt"),
               header = T)$x)
})

upRegMeio <- as.
downRegList <- lapply(seq_along(genotypes), function(x) {
  as.character(read.table(paste0(inDir,
                                 "res_", genotypes[x], "VCol_",
                                 as.character(FDR), "_lfcShrink_Chr_Sig",
                                 as.character(FDR), "_downRegGeneIDs.txt"),
               header = T)$x)
})

df <- data.frame(cbind(length(upRegList[[1]]),
                       length(upRegList[[2]]),
                       length(intersect(upRegList[[1]], upRegList[[2]])),
                       length(downRegList[[1]]),
                       length(downRegList[[2]]),
                       length(intersect(downRegList[[1]], downRegList[[2]]))))
colnames(df) <- c(paste0(genotype1, "_upReg_FDR", as.character(FDR)),
                  paste0(genotype2, "_upReg_FDR", as.character(FDR)),
                  paste0(genotype1, "_", genotype2, "_intersect_upReg_FDR", as.character(FDR)),
                  paste0(genotype1, "_downReg_FDR", as.character(FDR)),
                  paste0(genotype2, "_downReg_FDR", as.character(FDR)),
                  paste0(genotype1, "_", genotype2, "_intersect_downReg_FDR", as.character(FDR)))
write.table(df, paste0("./DEG_numbers_FDR", as.character(FDR), ".txt"),
            row.names = F, sep = "\t", quote = F)

vennPlotFun <- function(geneList, myColours, mainTitle, genotypes) { 
  grid.draw(venn.diagram(x = geneList,
                         filename = NULL,
                         fill = c(myColours),
                         alpha = c(0.5, 0.5),
                         cex = 2,
                         main.cex = 2,
                         cat.cex = 2,
                         fontfamily = "sans",
                         main.fontfamily = "sans",
                         cat.fontfamily = "sans",
                         main.fontface = "bold",
                         cat.fontface = "bold",
                         main = paste0(mainTitle, " (FDR < ", as.character(FDR), ")"),
                         category.names = genotypes,
                         cat.pos = c(-30, 30)))
}

geneList <- list(upRegList, downRegList)
direction <- c("upReg", "downReg")
description <- c("Up-regulated genes", "Down-regulated genes")

lapply(seq_along(direction), function(x) { 
  pdf(paste0("./", "VennDiagram_", direction[x], "_FDR", as.character(FDR), ".pdf"))
  vennPlotFun(geneList = geneList[[x]],
              myColours = myColours,
              mainTitle = description[x],
              genotypes = genotypes)
  dev.off()
})


