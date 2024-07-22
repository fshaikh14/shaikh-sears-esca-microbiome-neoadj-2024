rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
library(viridis)
library(reshape2)
library(readr)

# ------------------------------------------------------------------------
# set the working directory as the location of the script ---
# scriptPrefix = per user
# ------------------------------------------------------------------------

# Load CSV
Genus_LDA <- read_csv("Genus_LDA.csv")

GroupCols = c("#e4e4e4", "#eb3eab", "#FF9AA2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#8e7cc3", "#f94144","#f8961e","#f9b14a","#f9c74f","#90be6d","#82d3d1","#4d908e","#6ec2e5","#277da1")
numTx     = 19
GroupCols = colorRampPalette(GroupCols)(numTx+1)

A <- Genus_LDA
A$pathCR = paste0("PathCR ", as.character(A$pathCR))
A$pathCR = factor(A$pathCR, levels=c("PathCR 0", "PathCR 1"))

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
  
  p1
  
  ggsave(file="top.10.genera.barplot.pdf", plot=p1, width=10, height=7)
  
 