#################################################################
# Analyse differentially expressed genes for GO term enrichment #
# compared with all TAIR10 genes                                #
#################################################################

# Example usage:
# Rscript topGO_DEGs.R BP res_taf4bVCol_lfcShrink_Chr_Sig.1_downRegGeneIDs.txt

#options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
ont <- args[1]
target <- args[2]

inDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_01/" 

# load topGO
suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))

degGO <- function(target) {
  # Read in GO annotations for TAIR10 genes to define "gene universe"
  geneID2GO <- readMappings(file = paste0("/projects/ajt200/TAIR10/TAIR10_GO_reshaped_v280317_parentGeneIDs_v310318.txt"))
  geneUniverse <- names(geneID2GO)

  # Define list of genes of interest; file should contain a single column of gene identifiers
  genesOfInterest <- read.table(paste0(inDir, target))
  genesOfInterest <- as.character(genesOfInterest$x)

  # Specify where genes of interest appear in the gene universe vector
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse

  # Build GOdata object in topGO
  capture.output(GOdata <- new("topGOdata", description = "Differentially expressed genes", ontology = ont,
                               allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO),
                 file="/dev/null")

  # Access list of genes of interest
  #sg <- sigGenes(GOdata)
  #print(str(sg))
  #print(numSigGenes(GOdata))

  # Run Fisher's exact tests to determine GO term enrichment
  capture.output(resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultTopGO <- runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
                 file="/dev/null")

  # Count number of results where weight01 gives a P-value <= 0.01
  mySummary <- summary(attributes(resultTopGO)$score <= 0.01)
  numSignif <- as.integer(mySummary[[3]])

  # List significant results and write to file
  capture.output(enrichRes <- GenTable(GOdata, classicFisher = resultClassic,
                                       elimFisher = resultElim,
                                       topGOFisher = resultTopGO,
                                       orderBy = "topGOFisher",
                                       ranksOf = "elimFisher", topNodes = numSignif),
                 file="/dev/null")

  # WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
  basename <- basename(target)
  len <- nchar(basename)
  basename <- substr(basename, 1, len-4)
  
  out_name <- paste(basename, "GO", ont, "enrichment.tsv", sep="_")
  folder <- paste0(dirname(target), "/GO")
  system(paste0("[ -d ", folder, " ] || mkdir ", folder))
  #system(paste0("mkdir ", folder, "/", basename))
  folder2 <- paste0(folder, "/", basename, "_GO_", ont)
  system(paste0("[ -d ", folder2, " ] || mkdir ", folder2))

  capture.output(write.table(enrichRes, file = file.path(folder, out_name), sep = "\t",
                             row.names = FALSE, col.names = TRUE, quote = FALSE),
                 file="/dev/null")
  
  # Visualise the positions of the top 5 statistically significant GO terms in the GO hierarchy
  out_name2 <- paste(basename, "GO", ont, "enrichment", sep="_")
  printGraph(GOdata, resultTopGO, firstSigNodes = 5,
             fn.prefix = file.path(folder, out_name2), useInfo = "all", pdfSW = TRUE)

  # Extract gene IDs annotated with significantly enriched GO terms
  myTerms <- enrichRes$GO.ID
  myGenes <- genesInTerm(GOdata, myTerms)
  for(i in 1:length(myTerms)) {
    myTerm <- myTerms[i]
    myGenesForTerm <- myGenes[myTerm][[1]]
    myFactor <- myGenesForTerm %in% genesOfInterest
    myGenesForTermT <- myGenesForTerm[myFactor == TRUE]
    myGenesForTermT <- paste(myGenesForTermT, collapse = ",")
    myGenesForTermT <- paste(myTerm, myGenesForTermT, sep = "\t")
    out_name3 <- paste0(basename, "_GO_", ont, "_enrichment_", myTerm, ".txt")  
    write(myGenesForTermT, file = file.path(folder2, out_name3))
  }
}

# Apply degGO() function to target (differentially expressed genes file)
degGO(target)


sessionInfo()

