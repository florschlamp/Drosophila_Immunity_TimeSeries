# Scripts to generate plots in Figure 1
# last updated: March 2nd, 2021
# Florencia Schlamp
# mfs97@cornell.edu


library(ggplot2)
library(colorspace)

# Figure 1A) time points sampling chart

metadata <- read.csv("../data/metadata_table.csv")

time_points <- unique(metadata$hours)
data_base <- data.frame(hours=time_points,y=rep(1,21),timepoint=as.factor(1:21))
divergingx_hcl("Zissou 1",n=(length(unique(data_base$hours))+1))[1:21]

name = "../plots/Figure_1A_timepoints.pdf"
pdf(name,width=14,height=1,paper="special")

ggplot(data_base, aes(hours, y, color=timepoint)) + 
  #geom_line(size=1) +
  geom_point(size=3) +
  scale_color_manual(values = divergingx_hcl("Zissou 1",n=(length(unique(data_base$hours))+1))) +
  #theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position="none") +
  # xlab("Hours after injection",position="bottom") +
  scale_x_continuous(breaks=c(0,24,48,72,96,120),
                     labels = c("LPS injection","Day 1","Day 2","Day 3","Day 4","Day 5"),
                     position="top") +
  scale_y_continuous(breaks=c(1)) +
  geom_text(label=data_base$hours,vjust=2,color="black")

dev.off()


# Figure 1B) PCA plot

library(DESeq2)
library(ggplot2)

raw.counts <- read.delim("../data/raw_counts.csv", header=TRUE, sep=",", row.names=1)
metadata <- read.csv("../data/metadata_table.csv")

ddsMat <- DESeqDataSetFromMatrix(countData = raw.counts,
                                 colData = metadata,
                                 design = ~1) # design doesn't matter for this step
nrow(ddsMat)
# we start wtih 17,736 genes

# remove genes with 0 counts across all samples
ddsMat <- ddsMat[ rowSums(counts(ddsMat)) != 0, ]
nrow(ddsMat)
# 16,813 genes

rld <- rlog(ddsMat, blind=FALSE) #this step can take a long time

data <- plotPCA(rld, intgroup = c("time", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$hours <- factor(metadata$hours,levels=unique(metadata$hours))

name = "../plots/Figure_1B_PCAplot.pdf"
pdf(name,width=6.6,height=5,paper="special")

ggplot(data, aes(PC1, PC2, color=hours, shape=replicate)) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.title.x=element_text(size=12, margin=margin(15,0,0,0)),
        axis.title.y=element_text(size=12, margin=margin(0,15,0,0))) +
  scale_color_manual(values = divergingx_hcl("Zissou 1",n=(length(unique(data$hours))+1)))

dev.off()
