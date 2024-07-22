# FIGURE 1E - ALL TAXONOMIC LEVLES

rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
require(GGally)
require(gtools)
library(ggpubr)
# ---------------------------------------------------
# configure per user ---
# setwd to file location
analysisdir = "./analysis"
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
taxalevels = c("alpha-diversity",
               "phylum",
               "class",
               "order",
               "family",
               "genus",
               "species-otu.avgGTE0.01Prct",
               "picrust")

for (txlevel in taxalevels){
  # ---------------------------------------------------
  # create output dirs --- ddd
  txdir = paste(analysisdir, "/", txlevel, sep="") #subdirectory in analysisdir
  unlink(txdir, recursive=TRUE)
  dir.create(txdir)
  # ---------------------------------------------------
  pdfdir = paste(txdir, "/pdfs/", sep="") #subdirectory in analysisdir
  unlink(pdfdir, recursive=TRUE)
  dir.create(pdfdir)
  # ---------------------------------------------------
  tablefile = paste("../../Data-processing-taxa-assignment/analysis/A03.1-add-metadata/jhmi-esca-microbiota.", txlevel, ".txt", sep="")
  A         = read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
  Awmeta    = A
  # order specified factors for stats --
  Awmeta              = data.frame(Awmeta)

  Awmeta = Awmeta[Awmeta$pathCR != "n/a",]
  Awmeta = Awmeta[Awmeta$TimeCat2 != "other",]
  Awmeta = Awmeta[Awmeta$TimeCat2 != "replicate",]

  Awmeta$RecordID        = as.character(Awmeta$RecordID)
  Awmeta$days.from.start = as.numeric(as.character(Awmeta$days.from.start))
  Awmeta$Timepoint      = factor(Awmeta$Timepoint,      levels=c("PreTX", "OnTX"))
  Awmeta$Arm            = factor(Awmeta$Arm,            levels=c("A", "B", "B2"))
  Awmeta$Sex            = factor(Awmeta$Sex,            levels=c("F", "M"))
  Awmeta$Age            = as.numeric(as.character(Awmeta$Age))

  Awmeta                = Awmeta[order(Awmeta$pathCR),]
  Awmeta$RecordID       = factor(Awmeta$RecordID, levels=unique(Awmeta$RecordID))

  Awmeta$pathCR = factor(paste0("PathCR ",as.character(Awmeta$pathCR)))

  # --------------------------------------------------------
  # BEGIN color scheme for major sample metadata features of interest ---------
  mycolors = c()

  mycolors$pathCR["PathCR 0"]   = "#ed5565"
  mycolors$pathCR["PathCR 1"]   = "#1ab394"

  mycolors$path90["0"]   = "#ed5565"
  mycolors$path90["1"]   = "#1ab394"
  mycolors$path90["n/a"] = "#cccccc"

  mycolors$Sex["F"]   = "#ff9d9d"
  mycolors$Sex["M"]   = "#83c5df"

  mycolors$Arm["A"]   = "#b17eb8"
  mycolors$Arm["B"]   = "#e6cae9"

  mycolors$Timepoint["PreTX"]   = "#c7dee7"
  mycolors$Timepoint["OnTX"]    = "#5c9fb9"

  mycolors$Smoking["Never"]   = "#cccccc"
  mycolors$Smoking["Former"]  = "#ffc26a"
  # END color scheme for major sample metadata features of interest ---------
  # --------------------------------------------------------

  # create array of features to assess ---
  # colnames of Awmeta are too complicated at this point
  # e.g. "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales;
  # f__Pasteurellaceae; g__Avibacterium; s__Avibacterium_paragallinarum"
  # we will remove semicolons, spaces and simplify the names before visualization
  colnames(Awmeta)[26:ncol(Awmeta)] = gsub("k_Bacteria.", "", colnames(Awmeta)[26:ncol(Awmeta)])
  colnames(Awmeta)[26:ncol(Awmeta)] = gsub(";", ".", colnames(Awmeta)[26:ncol(Awmeta)])
  colnames(Awmeta)[26:ncol(Awmeta)] = gsub(" ", "", colnames(Awmeta)[26:ncol(Awmeta)])
  colnames(Awmeta)[26:ncol(Awmeta)] = gsub("__", ".", colnames(Awmeta)[26:ncol(Awmeta)])
  # completed modifying colnames of taxa ---
  taxaFeatures = colnames(Awmeta)[26:ncol(Awmeta)] # this is the list we will visualize

  # BEGIN STATS CODE -------------------------------------------
  # this will hold all stats results for Q1 / Q2
  txStatsResults = c()

  for (myFeature in taxaFeatures){
    print(paste("Evaluating feature: ", txlevel, " ", myFeature, sep=""))

    if (sum(Awmeta[,myFeature]>0) < 5){
      next
    }

    resultsRow = c(myFeature)

    my_comparisons = list( c("PathCR 0", "PathCR 1"))

    for (c1 in levels(Awmeta$pathCR)){
      nc1        = length(Awmeta[Awmeta$pathCR==c1,myFeature])
      meanc1     = mean(Awmeta[Awmeta$pathCR==c1, myFeature], na.rm=TRUE)
      resultsRow = c(resultsRow, nc1, meanc1)
    }
    for (compi in 1:length(my_comparisons)){
      mwP = wilcox.test(Awmeta[Awmeta$pathCR==my_comparisons[[compi]][1], myFeature],
                        Awmeta[Awmeta$pathCR==my_comparisons[[compi]][2], myFeature])$p.value
      resultsRow = c(resultsRow, mwP)
    }

    txStatsResults  = rbind(txStatsResults, resultsRow)

    if (as.numeric(as.character(resultsRow[6]))> 0.05){
      next
    }

    if (!grepl("alpha", txlevel)){
      Awmeta[,myFeature] = Awmeta[,myFeature] + 0.01
    }

    # END STATS CODE -------------------------------------------

    # BEGIN VIS CODE -------------------------------------------
    # add pseudo count for visualization / stats analysis

    # reorder sample IDs by pathCR and mean value
    aggMeanFeature   = aggregate(as.formula(paste0(myFeature, " ~ pathCR + RecordID")), data = Awmeta, FUN=mean)
    pidCROrder       = aggMeanFeature[order(aggMeanFeature$pathCR, aggMeanFeature[,myFeature], decreasing = TRUE), "RecordID"]

    # reorder sample IDs by pathCR and mean value
    aggMeanFeature   = aggregate(as.formula(paste0(myFeature, " ~ path90 + RecordID")), data = Awmeta, FUN=mean)
    pid90Order       = aggMeanFeature[order(aggMeanFeature$path90, aggMeanFeature[,myFeature], decreasing = TRUE), "RecordID"]

    myTitle = paste(myFeature, sep="")

    p1<-c()
    p2<-c()
    p3<-c()
    p4<-c()
    p5<-c()
    p6<-c()
    if (grepl("(alpha|picrust)", txlevel)){
      give.n <- function(x,min){
        return(c(y=min*0.9, label = length(x)))
      }

      Awmeta$RecordID = factor(Awmeta$RecordID, levels=pidCROrder)

      p1 <- ggplot(Awmeta, aes_string(x="RecordID", y=myFeature)) +
      geom_boxplot(aes(color=pathCR), fill=NA, outlier.size=0, coef=1e100, alpha=0.9) +
      geom_point(aes(fill=pathCR), color="black", pch=21, alpha=0.9, size=2.5, position=position_dodge(width=0.75)) +
      theme_bw() +
      scale_fill_manual(values=mycolors$pathCR) +
      theme(axis.text.x  = element_text(size=6, colour="black", angle=45, vjust=1, hjust=1),
            axis.text.y  = element_text(size=8, colour="black"),
            axis.title.x = element_text(size=12, colour="black"),
            axis.title.y = element_text(size=12, colour="black"),
            plot.title   = element_text(size=5, colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1) +
      xlab(NULL) + # no label for x axis
      ylab(paste(txlevel, " measure", sep="")) +
      ggtitle(myTitle)

      p2 <- ggbarplot(Awmeta, x="RecordID", y=myFeature, add="mean", fill="pathCR") +
        theme_bw() +
        scale_fill_manual(values=mycolors$pathCR) +
        scale_color_manual(values=mycolors$pathCR) +
        theme(axis.text.x      = element_text(size=8, colour="black", angle=45, hjust=1, vjust=1),
              axis.text.y      = element_text(size=10, colour="black"),
              axis.title.y     = element_text(size=9, colour="black"),
              plot.title       = element_text(size=8),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio     = 1.1,
              legend.position  = "none") +
        xlab(NULL)        +
        ylab(paste(txlevel, " measure", sep="")) +
        ggtitle(myTitle)

        p3 <- ggplot(Awmeta, aes_string(x="pathCR", y=myFeature)) +
        geom_boxplot(aes(color=pathCR), fill=NA, outlier.size=0, coef=1e100, alpha=0.9) +
        geom_point(aes(fill=pathCR), color="black", pch=21, alpha=0.9, size=2.5, position=position_dodge(width=0.75)) +
        theme_bw() +
        scale_fill_manual(values=mycolors$pathCR) +
        theme(axis.text.x  = element_text(size=8, colour="black", angle=45, vjust=1, hjust=1),
              axis.text.y  = element_text(size=9, colour="black"),
              axis.title.x = element_text(size=10, colour="black"),
              axis.title.y = element_text(size=10, colour="black"),
              plot.title   = element_text(size=5, colour="black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1.2) +
        xlab(NULL) + # no label for x axis
        ylab(paste(txlevel, " measure", sep="")) +
        ggtitle(myTitle) +
        scale_y_continuous(expand=c(0.05,0.07)) +
        stat_compare_means(aes(label=..p.signif..), comparisons=my_comparisons, label.y.npc='top', method="wilcox.test", size = 3)


    }else{
      give.n <- function(x,min){
        return(c(y=log10(min*0.5), label = length(x)))
      }
      Awmeta$RecordID = factor(Awmeta$RecordID, levels=pidCROrder)

        p3 <- ggplot(Awmeta, aes_string(x="pathCR", y=myFeature)) +
        geom_boxplot(aes(color=pathCR),fill=NA, outlier.size=0, coef=1e100, alpha=0.9) +
        geom_point(aes(fill=pathCR), color="black", pch=21, alpha=0.9, size=2.5, position=position_dodge(width=0.75)) +
        theme_bw() +
        scale_fill_manual(values=mycolors$pathCR) +
        theme(axis.text.x  = element_text(size=8, colour="black"),
              axis.text.y  = element_text(size=9, colour="black"),
              axis.title.x = element_text(size=10, colour="black"),
              axis.title.y = element_text(size=10, colour="black"),
              plot.title   = element_text(size=9, colour="black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1.2) +
        xlab(NULL) + # no label for x axis
        ylab("% Ab.") +
        ggtitle(myTitle) +
        scale_y_log10(breaks=c(0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100), expand=c(0.05,0.07)) +
        ggtitle(myTitle) +
        stat_compare_means(aes(label=..p.signif..), comparisons=my_comparisons, label.y.npc='top', method="wilcox.test", size = 3)

    }
    myFeatureName = substr(myFeature,0,200)
    ggsave(paste(pdfdir, myFeatureName, ".03.pdf", sep=""), plot=p3, width=4, height=3.5)

  } # end of features loop

  # # finalize colnames for txStatsResults ---
  txStatsResultsColNames = c("Feature")
  my_comparisons = list( c("PathCR 0", "PathCR 1"))

  for (c1 in levels(Awmeta$pathCR)){
    txStatsResultsColNames = c(txStatsResultsColNames, paste0(c1,".n"), paste0(c1,".mn"))
  }
  for (compi in 1:length(my_comparisons)){
    txStatsResultsColNames = c(txStatsResultsColNames, paste0(my_comparisons[[compi]][1],".vs.",my_comparisons[[compi]][2]))
  }

  colnames(txStatsResults) = txStatsResultsColNames
  txStatsResults = data.frame(txStatsResults)
  txStatsResults[,6] = as.numeric(as.character(txStatsResults[,6]))
  txStatsResults$p.adj = txStatsResults[,6]
  txStatsResults[txStatsResults$p.adj > 0.97, "p.adj"] <- NA
  txStatsResults$p.adj = p.adjust(txStatsResults$p.adj, method="fdr")

  write.table(x=txStatsResults, file=paste(txdir,"/", txlevel, "-stats-results.csv", sep=""), col.names=TRUE, row.names=FALSE, sep=",")
} # end of files loop
