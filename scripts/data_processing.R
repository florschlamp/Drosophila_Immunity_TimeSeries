# Scripts for data processing (getting normalized counts)
# last updated: March 2nd, 2021
# Florencia Schlamp
# mfs97@cornell.edu

library(edgeR)
library(ggplot2)

raw.counts <- read.delim("../data/raw_counts.csv", header=TRUE, 
                         sep=",", row.names=1,check.names = F)
metadata <- read.csv("../data/metadata_table.csv")

dim(raw.counts) # genes vs. samples
# 17736    41

# Remove rows with all zeros
# Genes with 0 counts across all samples are removed (either before
# or after normalization is the same)
keep=rowSums(raw.counts) != 0
raw.counts=raw.counts[keep,]

dim(raw.counts) # genes vs. samples
# 16813    41

# Normalize
libsizes=colSums(raw.counts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(raw.counts)/size.factor)
dim(cnts.norm)
# 16813    41

summary(libsizes/1000000)
# average of 23.5 million mapped reads

data_to_plot <- data.frame(sample_name = metadata$sample_ID,
                           libsize = libsizes)

ggplot(data_to_plot, aes(x = reorder(sample_name, -libsizes/1000000), y=libsizes/1000000)) +
  geom_bar(stat="identity") + coord_flip() +
  labs(y="Library size (in million of reads)", 
       title="Library size per sample (in million of reads)") +
  geom_hline(yintercept=summary(libsizes/1000000)[["Mean"]]) 
# keep an eye on sample 16B, much larger library size than the rest

# save table of normalized counts (before log transformation)
write.csv(cnts.norm,file="../data/normalized_counts.csv",row.names=TRUE)


# Add small yet appropriate positive count number (+1) to all genes 
#before doing log2 transformation, as to avoid -Inf and stabilize 
#the variance at low expression end.
cnts.norm.log=log2(cnts.norm+1)
range(cnts.norm.log)

# remove genes with low counts
hist(cnts.norm.log,breaks=30, main="histogram of normalized counts")

plotMDS(cnts.norm.log, main="MDS plot of normalized counts",
        col=c(rep("red",21),rep("blue",20)))

dim(cnts.norm.log) # genes vs. samples
# 16813    41

# FILTERING STEP
# only keep genes with more than 5 counts in at least 2 samples
keep=rowSums((cnts.norm.log)>5) >= 2
cnts.norm.log.filtered=cnts.norm.log[keep,]
dim(cnts.norm.log.filtered) # genes vs. samples
# 12657    41

hist(cnts.norm.log.filtered,breaks=30, 
     main="histogram of normalized counts after filtering")
# distribution looks more normal now

plotMDS(cnts.norm.log.filtered, main="MDS plot of normalized counts",
        col=c(rep("red",21),rep("blue",20)))

write.csv(cnts.norm.log.filtered,file="../data/normalized_counts_log_filt.csv",row.names=TRUE)

