

inDir <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/DESeq2_analysis_02/"

up <- read.csv(paste0(inDir, "up.colvstaf4b.annotate_UNIX2.csv"), header = T)
down <- read.csv(paste0(inDir, "down.colvstaf4b.annotate_UNIX2.csv"), header = T)

colnames(up)[7:8] <- c("leaf", "meiocyte")

up <- up[,1:8]
up <- cbind(up, up$leaf+1, up$meiocyte+1)
colnames(up)[9:10] <- c("leafPlus1", "meiocytePlus1")
up <- cbind(up, up$meiocytePlus1/up$leafPlus1)
colnames(up)[11] <- "meiocytePlus1_leafPlus1_ratio"
print("up mean meiocytePlus1:leafPlus1")
print(mean(up$meiocytePlus1_leafPlus1_ratio))
#[1] 12.96923
print("up mean meiocyte")
print(mean(up$meiocyte))
#[1] 2975.457
print("up mean leaf")
print(mean(up$leaf))
#[1] 1767.028
print("up sum meiocyte")
print(sum(up$meiocyte))
#[1] 2240519
print("up sum leaf")
print(sum(up$leaf))
#[1] 1330572
print("up sum meiocyte:sum leaf")
print((sum(up$meiocyte))/(sum(up$leaf)))
#[1] 1.683877
print("up mean baseMean")
print(mean(up$baseMean))
#[1] 1580.806

down <- down[,1:8]
down <- cbind(down, down$leaf+1, down$meiocyte+1)
colnames(down)[9:10] <- c("leafPlus1", "meiocytePlus1")
down <- cbind(down, down$meiocytePlus1/down$leafPlus1)
colnames(down)[11] <- "meiocytePlus1_leafPlus1_ratio"
print("down mean meiocytePlus1:leafPlus1")
print(mean(down$meiocytePlus1_leafPlus1_ratio))
#[1] 20.93781
print("down mean meiocyte")
print(mean(down$meiocyte))
#[1] 485.4817
print("down mean leaf")
print(mean(down$leaf))
#[1] 235.031
print("down sum meiocyte")
print(sum(down$meiocyte))
#[1] 954457
print("down sum leaf")
print(sum(down$leaf))
#[1] 462071
print("down sum meiocyte:sum leaf")
print((sum(down$meiocyte))/(sum(down$leaf)))
#[1] 2.065607
print("down mean baseMean")
print(mean(down$baseMean))
#[1] 88.08087


