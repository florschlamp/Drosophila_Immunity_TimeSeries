# Scripts to generate plots in Figure 2
# last updated: March 19th, 2021
# Florencia Schlamp
# mfs97@cornell.edu

library(reshape2)
library(ggplot2)

# Figure 2B and 2C) Plot smoothing vs normalized counts

# read in data:

smooth_short_A <- read.csv("../results/top_genes_limma-voom_short.filtered_A_spline_3_TMM_ALL_SMOOTHED.csv",
                           header=TRUE, row.names=1, as.is=TRUE, check.names=FALSE, sep=',')
smooth_short_B <- read.csv("../results/top_genes_limma-voom_short.filtered_B_spline_3_TMM_ALL_SMOOTHED.csv",
                           header=TRUE, row.names=1, as.is=TRUE, check.names=FALSE, sep=',')
smooth_short_1to7_A <- read.csv("../results/top_genes_limma-voom_short.filtered_A_spline_3_TMM_1to7_ALL_SMOOTHED.csv",
                                header=TRUE, row.names=1, as.is=TRUE, check.names=FALSE, sep=',')
smooth_short_1to7_B <- read.csv("../results/top_genes_limma-voom_short.filtered_B_spline_3_TMM_1to7_ALL_SMOOTHED.csv",
                                header=TRUE, row.names=1, as.is=TRUE, check.names=FALSE, sep=',')

dat_norm <- read.table('../data/normalized_counts_log_filt.csv',
                       header=TRUE, row.names=1, as.is=TRUE, 
                       check.names=FALSE, sep=',')
dat_norm_short <- subset(dat_norm, select = -c(4,19,20,21,39,40,41))

sample_table <- read.csv("../data/metadata_table.csv",header=TRUE)
sample_table_short <- sample_table[-c(4,19,20,21,39,40,41),]

load("../data/translate_geneID_genename.rdata")


# Figure 2B) smoothing of short term patterns
# example with Gale and Galk

genes <- c("Gale","Galk")

list_of_interest <- translation_file_filtered[translation_file_filtered$gene_symbol %in% genes,]$geneID
type = 'name'

genes_from_time_course <- dat_norm_short # choose if dat_norm or dat_norm_short
start_table <- sample_table_short # choose if sample_table or sample_table_short

genes_from_time_course_filt <- genes_from_time_course[list_of_interest,]
smooth_short_A_filt <- smooth_short_A[list_of_interest,]
smooth_short_B_filt <- smooth_short_B[list_of_interest,]
smooth_short_1to7_A_filt <- smooth_short_1to7_A[list_of_interest,]
smooth_short_1to7_B_filt <- smooth_short_1to7_B[list_of_interest,]

rownames(genes_from_time_course_filt) <- as.character(plyr::mapvalues(x=rownames(genes_from_time_course_filt),
                                                                      from=translation_file_filtered$geneID,
                                                                      to=translation_file_filtered$gene_symbol,
                                                                      warn_missing = F))
rownames(smooth_short_A_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_A_filt),
                                                              from=translation_file_filtered$geneID,
                                                              to=translation_file_filtered$gene_symbol,
                                                              warn_missing = F))
rownames(smooth_short_B_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_B_filt),
                                                              from=translation_file_filtered$geneID,
                                                              to=translation_file_filtered$gene_symbol,
                                                              warn_missing = F))
rownames(smooth_short_1to7_A_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_1to7_A_filt),
                                                                   from=translation_file_filtered$geneID,
                                                                   to=translation_file_filtered$gene_symbol,
                                                                   warn_missing = F))
rownames(smooth_short_1to7_B_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_1to7_B_filt),
                                                                   from=translation_file_filtered$geneID,
                                                                   to=translation_file_filtered$gene_symbol,
                                                                   warn_missing = F))


datalist = list()

for(word in genes){
  if(word %in% rownames(genes_from_time_course_filt)){
    gene_counts <- melt(genes_from_time_course_filt[word,])$"value"  -4.5 #value added to match patterns on y-axis
    start_table["counts"] <- gene_counts
    start_table["gene_ID"] <- c(rep(word,length(gene_counts)))
    start_table["type"] <- c(rep("norm_counts",length(gene_counts)))
    start_table["gene"] <- c(rep(word,length(gene_counts)))
    datalist[[word]] <- start_table
    
    start_table <- sample_table_short
    smooth_3_A <- melt(smooth_short_A_filt[word,])$"value"
    smooth_3_B <- melt(smooth_short_B_filt[word,])$"value"
    start_table["counts"] <- c(smooth_3_A,smooth_3_B)
    start_table["gene_ID"] <- c(rep(paste(word,"_smooth",sep=""),length(gene_counts)))
    start_table["type"] <- c(rep("48h_spline",length(gene_counts)))
    start_table["gene"] <- c(rep(word,length(gene_counts)))
    datalist[[paste(word,"_smooth",sep="")]] <- start_table
    
    start_table <- sample_table_short
    smooth_1to7_A <- c(melt(smooth_short_1to7_A_filt[word,])$"value",rep(NA,10))
    smooth_1to7_B <- c(melt(smooth_short_1to7_B_filt[word,])$"value",rep(NA,10))
    start_table["counts"] <- c(smooth_1to7_A,smooth_1to7_B)
    start_table["gene_ID"] <- c(rep(paste(word,"_smooth_1to7",sep=""),length(gene_counts)))
    start_table["type"] <- c(rep("8h_spline",length(gene_counts)))
    start_table["gene"] <- c(rep(word,length(gene_counts)))
    datalist[[paste(word,"_smooth_1to7",sep="")]] <- start_table
  }
}

big_data = do.call(rbind,datalist)

big_data['group'] <- paste(big_data$gene_ID,big_data$replicate,sep="_")

name = "../plots/Figure_2B_smoothing_short_Gale-Galk.pdf"
pdf(name,width=7,height=3,paper="special")

ggplot(big_data, aes(x=hours, y=counts, color=type, group=group)) +
  scale_y_log10() + 
  geom_point(size=2) + 
  geom_line(size=1, aes(linetype=gene)) + 
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.text.y=element_blank(),
        legend.title=element_text(size=12),
        axis.title.x=element_text(size=12, margin=margin(5,0,0,0)),
        axis.title.y = element_blank()) +
  scale_x_continuous(breaks=c(0,6,12,24,36,48)) +
  scale_color_manual(values=c("#5F80F3","#dbc300","#17CE6B")) +
  #ylab("normalized counts") +
  xlab("Hours after injection")

dev.off()



# Figure 2C) smoothing of long term patterns
# example with AttA and DptB

genes <- c("AttA","DptB")

list_of_interest <- translation_file_filtered[translation_file_filtered$gene_symbol %in% genes,]$geneID
type = 'name'

genes_from_time_course <- dat_norm_short # choose if dat_norm or dat_norm_short
start_table <- sample_table_short # choose if sample_table or sample_table_short

genes_from_time_course_filt <- genes_from_time_course[list_of_interest,]
smooth_short_A_filt <- smooth_short_A[list_of_interest,]
smooth_short_B_filt <- smooth_short_B[list_of_interest,]
smooth_short_1to7_A_filt <- smooth_short_1to7_A[list_of_interest,]
smooth_short_1to7_B_filt <- smooth_short_1to7_B[list_of_interest,]

rownames(genes_from_time_course_filt) <- as.character(plyr::mapvalues(x=rownames(genes_from_time_course_filt),
                                                                      from=translation_file_filtered$geneID,
                                                                      to=translation_file_filtered$gene_symbol,
                                                                      warn_missing = F))
rownames(smooth_short_A_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_A_filt),
                                                              from=translation_file_filtered$geneID,
                                                              to=translation_file_filtered$gene_symbol,
                                                              warn_missing = F))
rownames(smooth_short_B_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_B_filt),
                                                              from=translation_file_filtered$geneID,
                                                              to=translation_file_filtered$gene_symbol,
                                                              warn_missing = F))
rownames(smooth_short_1to7_A_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_1to7_A_filt),
                                                                   from=translation_file_filtered$geneID,
                                                                   to=translation_file_filtered$gene_symbol,
                                                                   warn_missing = F))
rownames(smooth_short_1to7_B_filt) <- as.character(plyr::mapvalues(x=rownames(smooth_short_1to7_B_filt),
                                                                   from=translation_file_filtered$geneID,
                                                                   to=translation_file_filtered$gene_symbol,
                                                                   warn_missing = F))


datalist = list()

for(word in genes){
  if(word %in% rownames(genes_from_time_course_filt)){
    gene_counts <- melt(genes_from_time_course_filt[word,])$"value"  -4.5 #value added to match patterns on y-axis
    start_table["counts"] <- gene_counts
    start_table["gene_ID"] <- c(rep(word,length(gene_counts)))
    start_table["type"] <- c(rep("norm_counts",length(gene_counts)))
    start_table["gene"] <- c(rep(word,length(gene_counts)))
    datalist[[word]] <- start_table
    
    start_table <- sample_table_short
    smooth_3_A <- melt(smooth_short_A_filt[word,])$"value"
    smooth_3_B <- melt(smooth_short_B_filt[word,])$"value"
    start_table["counts"] <- c(smooth_3_A,smooth_3_B)
    start_table["gene_ID"] <- c(rep(paste(word,"_smooth",sep=""),length(gene_counts)))
    start_table["type"] <- c(rep("48h_spline",length(gene_counts)))
    start_table["gene"] <- c(rep(word,length(gene_counts)))
    datalist[[paste(word,"_smooth",sep="")]] <- start_table
    
    start_table <- sample_table_short
    smooth_1to7_A <- c(melt(smooth_short_1to7_A_filt[word,])$"value",rep(NA,10))
    smooth_1to7_B <- c(melt(smooth_short_1to7_B_filt[word,])$"value",rep(NA,10))
    start_table["counts"] <- c(smooth_1to7_A,smooth_1to7_B)
    start_table["gene_ID"] <- c(rep(paste(word,"_smooth_1to7",sep=""),length(gene_counts)))
    start_table["type"] <- c(rep("8h_spline",length(gene_counts)))
    start_table["gene"] <- c(rep(word,length(gene_counts)))
    datalist[[paste(word,"_smooth_1to7",sep="")]] <- start_table
  }
}

big_data = do.call(rbind,datalist)

big_data['group'] <- paste(big_data$gene_ID,big_data$replicate,sep="_")

name = "../plots/Figure_2C_smoothing_long_AttA-DptB.pdf"
pdf(name,width=7,height=3,paper="special")

ggplot(big_data, aes(x=hours, y=counts, color=type, group=group)) +
  scale_y_log10() + 
  geom_point(size=2) + 
  geom_line(size=1, aes(linetype=gene)) + 
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.text.y=element_blank(),
        legend.title=element_text(size=12),
        axis.title.x=element_text(size=12, margin=margin(5,0,0,0)),
        axis.title.y = element_blank()) +
  scale_x_continuous(breaks=c(0,6,12,24,36,48)) +
  scale_color_manual(values=c("#5F80F3","#dbc300","#17CE6B")) +
  #ylab("normalized counts") +
  xlab("Hours after injection")

dev.off()

