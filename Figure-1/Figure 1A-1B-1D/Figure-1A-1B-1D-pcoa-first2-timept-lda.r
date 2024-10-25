rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
library(viridis)
library(vegan)
library(ape)
library(Rtsne)
library(umap)
require(GGally)
require(gtools)
library(ggpubr)

# set the working directory to the exact code base location ------------------

# set the place for the analysis outputs ------------------

# load meta-data ¬
meta <- read.table("../../Data-processing-taxa-assignment/analysis/A03.1-add-metadata/jhmi-esca-microbiota.genus.txt", header=TRUE, sep="\t")
lastmdex = 25
meta     = meta[,1:lastmdex]
meta     = data.frame(meta)

meta = meta[meta$pathCR   != "n/a",]
meta = meta[meta$TimeCat2 != "other",]
meta = meta[meta$TimeCat2 != "replicate",]

meta$RecordID        = as.character(meta$RecordID)
meta$days.from.start = as.numeric(as.character(meta$days.from.start))
meta$Timepoint      = factor(meta$Timepoint,     levels=c("PreTX", "OnTX"))
meta$Arm            = factor(meta$Arm,           levels=c("A", "B", "B2"))
meta$Sex            = factor(meta$Sex,           levels=c("F", "M"))
meta$Age            = as.numeric(as.character(meta$Age))

meta$pathCR90Comp = paste0(meta$pathCR, meta$path90)
meta$pathCR90Comp = gsub("11", "path90+CR+", meta$pathCR90Comp)
meta$pathCR90Comp = gsub("01", "path90+CR-", meta$pathCR90Comp)
meta$pathCR90Comp = gsub("00", "path90-CR-", meta$pathCR90Comp)

meta$pathCR90Comp = factor(meta$pathCR90Comp)
meta$path90       = factor(meta$path90)

meta$pathCR = paste0("PathCR ", as.character(meta$pathCR))

meta$cap.agg = "G0"
meta[meta$CAP.regression.by.central.review %in% c("G1","G2","G3"),"cap.agg"] <- "G1-G3"
meta$cap.agg = factor(meta$cap.agg, levels=c("G0", "G1-G3"))

# BEGIN color scheme for major sample metadata features of interest ---------
mycolors = c()

mycolors$pathCR["PathCR 0"]   = "#ed5565"
mycolors$pathCR["PathCR 1"]   = "#1ab394"

mycolors$path90["0"]   = "#ed5565"
mycolors$path90["1"]   = "#1ab394"
mycolors$path90["n/a"] = "#cccccc"

mycolors$pathCR90Comp["path90-CR-"]   = "#ed5565"
mycolors$pathCR90Comp["path90+CR-"]   = "#95d5c8"
mycolors$pathCR90Comp["path90+CR+"]   = "#1ab394"

mycolors$Sex["F"]   = "#ff9d9d"
mycolors$Sex["M"]   = "#83c5df"

mycolors$Arm["A"]   = "#b17eb8"
mycolors$Arm["B"]   = "#e6cae9"
mycolors$Arm["B2"]  = "#f1e7f2"

mycolors$TimeCat1["first"]     = "#6d8bdb"
mycolors$TimeCat1["other"]     = "#cccccc"
mycolors$TimeCat1["replicate"] = "#efad76"

mycolors$TimeCat2["first2"]     = "#5cb9b6"
mycolors$TimeCat2["other"]      = "#cccccc"
mycolors$TimeCat2["replicate"]  = "#efad76"

mycolors$Timepoint["PreTX"]   = "#5c9fb9"
mycolors$Timepoint["OnTX"]    = "#c7dee7"

mycolors$Smoking["Never"]   = "#cccccc"
mycolors$Smoking["Former"]  = "#ffc36a"

mycolors$CAP.regression.by.central.review["G0"] = "#277da1"
mycolors$CAP.regression.by.central.review["G1"] = "#43aa8b"
mycolors$CAP.regression.by.central.review["G2"] = "#f9c74f"
mycolors$CAP.regression.by.central.review["G3"] = "#f94144"


mycolors$cap.agg["G0"]    = "#277da1"
mycolors$cap.agg["G1-G3"] = "#f9c74f"

# END color scheme for major sample metadata features of interest ---------


## FIGURE 1A, 1B

# begin custom function for PCoA analyses
searshsimmunotxpcoa <- function(file, prefix){
  print(paste("PCoA using distance metric: ", prefix, sep= ""))

  # load and refine distance matrix ------
  bc.D <- read.table(file, header=TRUE, sep="\t", check.names=FALSE)
  rownames(bc.D) = bc.D[,1]
  bc.D           = bc.D[,-c(1)]
  bc.D           = bc.D[rownames(bc.D) %in% meta[,1],  colnames(bc.D) %in% meta[,1] ]
  # execute PCoA ------
  bc.pcoa        = pcoa(bc.D)
  # extract principal coordinates from PCoA results ------
  bc.df          = data.frame(PC1 = bc.pcoa$vectors[,1], PC2 = -1*bc.pcoa$vectors[,2], PC3 = bc.pcoa$vectors[,3])
  # extract variance components from PCoA results ------
  bc.Ax1PrctVar  = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[1])
  bc.Ax2PrctVar  = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[2])
  bc.Ax3PrctVar  = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[3])

  bc.df          = cbind(meta[match(rownames(bc.df), meta[,1]),],bc.df)
  bc.df$SampleID = as.character(bc.df$SampleID)

  bc.df$RecordID        = as.character(bc.df$RecordID)
  bc.df$days.from.start = as.numeric(as.character(bc.df$days.from.start))
  bc.df$Timepoint       = factor(bc.df$Timepoint,      levels=c("PreTX", "OnTX", "PostTX"))
  bc.df$Arm             = factor(bc.df$Arm,            levels=c("A", "B", "B2"))
  bc.df$Sex             = factor(bc.df$Sex,            levels=c("F", "M"))
  bc.df$Age             = as.numeric(as.character(bc.df$Age))

  rownames(bc.df) = bc.df$SampleID

  if (!identical(rownames(bc.df), rownames(bc.D))){
    stop("bc.df and bc.D names don't match")
  }

  # PERMANOVA analysis -----
  capture.output(adonis(bc.D ~ pathCR, data = bc.df, permutations = 1000), file=paste(analysisdir, prefix,".permanova-results.txt",sep=""))

  centroids <- aggregate(cbind(PC1,PC2)~pathCR, data=bc.df, mean)
  bc.df     <- merge(bc.df,centroids,by="pathCR",suffixes=c("",".centroid"))

  #PERMANONA analysis - with correction for Study Arm, remove # to run code
  # capture.output(adonis(bc.D ~ pathCR+Arm, data = bc.df, permutations = 1000), file=paste(analysisdir, prefix,".permanova-results.txt",sep=""))
  # centroids <- aggregate(cbind(PC1,PC2)~pathCR, data=bc.df, mean)
  # bc.df     <- merge(bc.df,centroids,by="pathCR",suffixes=c("",".centroid"))
  
    # begin ggplot code for PCoA
  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=pathCR), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=pathCR), alpha=0.5) +
  geom_point(aes(fill=pathCR, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$pathCR)  +
  scale_color_manual(values=mycolors$pathCR) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.1.1.pdf",sep=""), p1, width=6.5, height=3.5)
  write.csv(bc.df, file=paste(analysisdir, prefix,".pcoa.1.1.csv",sep=""))

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=pathCR), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=pathCR), alpha=0.5) +
  geom_point(aes(fill=pathCR, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$pathCR)  +
  scale_color_manual(values=mycolors$pathCR) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1) +
  facet_wrap(~pathCR, ncol=2) +
  ggtitle(paste0("pathCR")) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.1.2.pdf",sep=""), p1, width=6.5, height=3.5)

  # path 90
  bc.df     = bc.df[,!colnames(bc.df) %in% c("PC1.centroid", "PC2.centroid")]
  centroids <- aggregate(cbind(PC1,PC2)~path90, data=bc.df, mean)
  bc.df     <- merge(bc.df,centroids,by="path90",suffixes=c("",".centroid"))

  # PERMANOVA analysis -----
  capture.output(adonis(bc.D ~ path90, data = bc.df, permutations = 1000), file=paste(analysisdir, prefix,".permanova-results.txt",sep=""), append=TRUE)

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=path90), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=path90), alpha=0.5) +
  geom_point(aes(fill=path90, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$path90)  +
  scale_color_manual(values=mycolors$path90) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1) +
  ggtitle(paste0("path90")) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.2.1.pdf",sep=""), p1, width=6.5, height=3.5)

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=path90), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=path90), alpha=0.5) +
  geom_point(aes(fill=path90, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$path90)  +
  scale_color_manual(values=mycolors$path90) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1) +
  facet_wrap(~path90, ncol=2) +
  ggtitle(paste0("path90")) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.2.2.pdf",sep=""), p1, width=6.5, height=3.5)

  # PERMANOVA analysis -----
  capture.output(adonis(bc.D ~ pathCR90Comp, data = bc.df, permutations = 1000), file=paste(analysisdir, prefix,".permanova-results.txt",sep=""), append=TRUE)

  bc.df     = bc.df[,!colnames(bc.df) %in% c("PC1.centroid", "PC2.centroid")]

  centroids <- aggregate(cbind(PC1,PC2)~pathCR90Comp, data=bc.df, mean)
  bc.df     <- merge(bc.df,centroids,by="pathCR90Comp",suffixes=c("",".centroid"))

  # begin ggplot code for PCoA
  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=pathCR90Comp), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=pathCR90Comp), alpha=0.5) +
  geom_point(aes(fill=pathCR90Comp, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$pathCR90Comp)  +
  scale_color_manual(values=mycolors$pathCR90Comp) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.3.1.pdf",sep=""), p1, width=6.5, height=3.5)

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=pathCR90Comp), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=pathCR90Comp), alpha=0.5) +
  geom_point(aes(fill=pathCR90Comp, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$pathCR90Comp)  +
  scale_color_manual(values=mycolors$pathCR90Comp) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.1) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~pathCR90Comp, ncol=3)
  ggsave(paste(analysisdir, prefix,".pcoa.3.2.pdf",sep=""), p1, width=8, height=5)

  # add by Group G0-G3 status --
  # PERMANOVA analysis -----
  capture.output(adonis(bc.D ~ CAP.regression.by.central.review, data = bc.df, permutations = 1000), file=paste(analysisdir, prefix,".permanova-results.txt",sep=""), append=TRUE)

  bc.df     = bc.df[,!colnames(bc.df) %in% c("PC1.centroid", "PC2.centroid")]

  centroids <- aggregate(cbind(PC1,PC2)~CAP.regression.by.central.review, data=bc.df, mean)
  bc.df     <- merge(bc.df,centroids,by="CAP.regression.by.central.review",suffixes=c("",".centroid"))

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=CAP.regression.by.central.review), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=CAP.regression.by.central.review), alpha=0.5) +
  geom_point(aes(fill=CAP.regression.by.central.review, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$CAP.regression.by.central.review)  +
  scale_color_manual(values=mycolors$CAP.regression.by.central.review) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.4.1.pdf",sep=""), p1, width=8, height=5)

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=CAP.regression.by.central.review), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=CAP.regression.by.central.review), alpha=0.5) +
  geom_point(aes(fill=CAP.regression.by.central.review, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$CAP.regression.by.central.review)  +
  scale_color_manual(values=mycolors$CAP.regression.by.central.review) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.1) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~CAP.regression.by.central.review, ncol=2)
  ggsave(paste(analysisdir, prefix,".pcoa.4.2.pdf",sep=""), p1, width=8, height=5)


  # add by Group G0 vs G1-G3 status --
  # PERMANOVA analysis -----
  capture.output(adonis(bc.D ~ cap.agg, data = bc.df, permutations = 1000), file=paste(analysisdir, prefix,".permanova-results.txt",sep=""), append=TRUE)

  bc.df     = bc.df[,!colnames(bc.df) %in% c("PC1.centroid", "PC2.centroid")]

  centroids <- aggregate(cbind(PC1,PC2)~cap.agg, data=bc.df, mean)
  bc.df     <- merge(bc.df,centroids,by="cap.agg",suffixes=c("",".centroid"))

  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=cap.agg), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=cap.agg), alpha=0.5) +
  geom_point(aes(fill=cap.agg, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$cap.agg)  +
  scale_color_manual(values=mycolors$cap.agg) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.1) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, prefix,".pcoa.5.1.pdf",sep=""), p1, width=8, height=5)

  write.csv(bc.df, file=paste(analysisdir, prefix,".pcoa.5.1.csv",sep=""))


  p1 <- ggplot(bc.df, aes(x=PC1,y=PC2)) +
  geom_point(data=centroids, aes(x=PC1,y=PC2, fill=cap.agg), pch=21, size=1) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=cap.agg), alpha=0.5) +
  geom_point(aes(fill=cap.agg, shape=Arm), color="black", size=4, alpha=0.7) +
  scale_fill_manual(values=mycolors$cap.agg)  +
  scale_color_manual(values=mycolors$cap.agg) +
  scale_shape_manual(values=c(23,21,22)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.1) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~cap.agg, ncol=2)
  ggsave(paste(analysisdir, prefix,".pcoa.5.2.pdf",sep=""), p1, width=8, height=5)

}

# run Bray Curtis file through this new function -----------------------------------------------------------
searshsimmunotxpcoa(file="../analysis/A03.1-add-metadata/jhmi-esca-microbiota.beta-diversity.bray_curtis.txt",        prefix = "bray-curtis")


# FIGURE 1D


# NDMS ---
# load meta-data
meta       <- read.table("../analysis/A03.1-add-metadata/jhmi-esca-microbiota.genus.txt", header=TRUE, sep="\t")
metao       = meta[meta$pathCR != "n/a",]
meta        = meta[meta$TimeCat2 != "n/a",]
meta        = meta[meta$TimeCat2 == "first2",]
meta        = metao[,c(20,26:ncol(metao))]
meta        = data.frame(meta)
meta$pathCR = factor(meta$pathCR)
vals        = meta[,2:ncol(meta)]
input       = cbind(meta[,"pathCR"], vals[,order(colMeans(vals), decreasing=TRUE)[1:10]])
colnames(input)[1] = "pathCR"

# -- lda
ldares <- lda(pathCR ~ ., input)
ldares.values <- predict(ldares, input)

write.csv(input, file=paste0(analysisdir,"lda.input.csv"), row.names=FALSE)

# TSNE plot using all genes --
tinvals   = vals
for (j in 1:ncol(tinvals)){
  #tinvals[,j] = tinvals[,j]+abs(rnorm(nrow(tinvals),0,0.0001))
}
tsneA     = Rtsne(tinvals, perplexity = 5)
umapA     = umap(tinvals)
tsplotDat = cbind(meta,metao$Arm, tsneA$Y[,1:2],umapA$layout[,1:2], ldares.values$x)
colnames(tsplotDat) = c(colnames(meta), "Arm", "tsneX", "tsneY", "umapX", "umapY", "ldaX")

tsplotDat$pathCR = paste0("PathCR ", as.character(tsplotDat$pathCR))

my_comparisons = list(c("PathCR 0","PathCR 1"))

mycolors$pathCR["PathCR 0"]   = "#ed5565"
mycolors$pathCR["PathCR 1"]   = "#1ab394"

p1 <- ggplot(tsplotDat, aes(x=pathCR,y=ldaX)) +
geom_boxplot(aes(color=pathCR), fill=NA, outlier.size=0, outlier.shape=NA, coef=1e100, alpha=0.9) +
geom_point(aes(fill=pathCR, shape=Arm), color="black", size=2.75, alpha=0.45) +
scale_fill_manual(values=mycolors$pathCR)  +
scale_color_manual(values=mycolors$pathCR) +
scale_shape_manual(values=c(23,21,22)) +
xlab(NULL) +
ylab("Linear Discriminant Axis") +
theme_bw() +
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio=1.2) +
guides(fill = guide_legend(override.aes=list(shape=21))) +
scale_y_continuous(expand=c(0.05,0.07)) +
stat_compare_means(aes(label=..p.signif..), comparisons=my_comparisons, label.y.npc='top', method="wilcox.test", size = 3)
ggsave(paste(analysisdir, "lda.4.1.pdf",sep=""), p1, width=4, height=4)
