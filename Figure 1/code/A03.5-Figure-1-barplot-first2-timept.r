rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
library(viridis)
library(reshape2)

# ------------------------------------------------------------------------
# set the working directory as the location of the script ---
setwd("/Users/jwhite/Desktop/jhmi-esca-microbiota-2022/code")
scriptPrefix = "A03.5-barplot-first2-timept"
# ------------------------------------------------------------------------
# create output directory
analysisdir = paste("../analysis/", scriptPrefix, sep="")
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)

pdfdir = paste(analysisdir, "/pdf", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ------------------------------------------------------------------------
inDir  = "../analysis/A03.1-add-metadata/";
files  = c("jhmi-esca-microbiota.phylum.txt",
           "jhmi-esca-microbiota.class.txt",
           "jhmi-esca-microbiota.order.txt",
           "jhmi-esca-microbiota.family.txt",
           "jhmi-esca-microbiota.genus.txt",
           "jhmi-esca-microbiota.species-otu.avgGTE0.01Prct.txt")

GroupCols = c("#e4e4e4", "#eb3eab", "#FF9AA2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#8e7cc3", "#f94144","#f8961e","#f9b14a","#f9c74f","#90be6d","#82d3d1","#4d908e","#6ec2e5","#277da1")
numTx     = 19
GroupCols = colorRampPalette(GroupCols)(numTx+1)
for (f in files){
  fpdfdir = paste(pdfdir, "/", gsub("(jhmi-esca-microbiota.|.txt)","",f), sep="")
  dir.create(fpdfdir)
  print(paste(f, "...", sep=""))
  tablefile = paste(inDir,f,sep="")
  A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
  A$pathCR = paste0("PathCR ", as.character(A$pathCR))
  A        = A[A$pathCR != "n/a",]
  A$pathCR = factor(A$pathCR, levels=c("PathCR 0", "PathCR 1"))
  A        = A[A$TimeCat2 != "other",]
  A        = A[A$TimeCat2 != "replicate",]

  colnames(A) = gsub("k_Bacteria.","",colnames(A))

  if (grepl("genus",f)){
    colnames(A) = gsub("p_Firmicutes\\.|p_Actinobacteria\\.|p_Bacteroidetes\\.", "", colnames(A))
  }

  A = A[,c(6,20,26:ncol(A))]

  aggA            = aggregate( . ~ RecordID + pathCR, data=A, FUN=mean)
  cmaggA          = colMeans(aggA[3:ncol(aggA)])
  selectedA       = names(cmaggA[order(cmaggA, decreasing=TRUE)])[1:min(numTx,length(cmaggA))]
  aggA            = aggA[,c("RecordID", "pathCR", selectedA)]
  aggA$Other.taxa = 100-rowSums(aggA[3:ncol(aggA)])

  meltA           = melt(aggA, id.vars = c("RecordID", "pathCR"))
  colnames(meltA) = c("RecordID", "pathCR", "Taxa", "Abundance")

  meltA$pathCR    = factor(meltA$pathCR)

  meltA$Taxa      = factor(meltA$Taxa, levels=rev(levels(meltA$Taxa)))

  customOrder     = aggA[order(rowSums(aggA[,selectedA[1:2]]), decreasing=TRUE),"RecordID"]

  meltA$RecordID  = factor(meltA$RecordID, levels=customOrder)

  p1 <- ggplot(meltA, aes(x=RecordID, y=Abundance, fill=Taxa)) +
  geom_bar(position="stack", stat="identity", colour="white", size=0.25) +
  theme_bw() +
  scale_fill_manual(values=GroupCols) +
  theme(axis.text.x = element_text(size=6, colour="black", angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size=10, colour="black"),
        legend.text      = element_text(size=5),
        legend.title     = element_blank()) +
  xlab(NULL) +
  ylab("% Abundance") +
  facet_grid(~pathCR, space="free_x", scale="free_x") +
  ggtitle("pathCR")
  ggsave(file=paste(fpdfdir, "/",gsub("(jhmi-esca-microbiota.|.txt)","",f), "01.pdf", sep=""), plot=p1, width=10, height=7)

  p1 <- ggplot(meltA, aes(x=RecordID, y=Abundance, fill=Taxa)) +
  geom_bar(position="stack", stat="identity", colour="white", size=0.25) +
  theme_bw() +
  scale_fill_manual(values=GroupCols) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        axis.ticks.x     = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size=10, colour="black"),
        legend.text      = element_text(size=8),
        legend.title     = element_blank()) +
  xlab(NULL) +
  ylab("% Abundance") +
  facet_grid(~pathCR, space="free_x", scale="free_x") +
  scale_y_continuous(limits=c(0,100), expand=c(0,0))
  ggsave(file=paste(fpdfdir, "/",gsub("(jhmi-esca-microbiota.|.txt)","",f), "02.pdf", sep=""), plot=p1, width=12, height=6)


  write.csv(aggA, file=paste(fpdfdir, "/",gsub("(jhmi-esca-microbiota.|.txt)","",f), "01.csv", sep=""), row.names=F)

}
