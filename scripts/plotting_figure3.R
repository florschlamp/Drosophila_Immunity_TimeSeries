# Scripts to generate plots in Figure 3
# last updated: March 16th, 2021
# Florencia Schlamp
# mfs97@cornell.edu

# Figure 3A) heatmap top 91 core DE genes

library(pheatmap)
library(RColorBrewer)
library(dendsort)

all_logFCs <- read.csv("../results/log_fold_change_short.filtered_ALL_spline_3.csv", row.names = 1)
load("../results/gene_lists/pairwiseDE_pval0.01_logFC2_2tp_91genes.rdata")
subset_logFCs <- all_logFCs[core_91,]
mat <- as.matrix(subset_logFCs)
dat <- data.frame(values = as.numeric(mat))
colnames(mat) <- c(1,2,4,5,6,8,10,12,14,16,20,24,30,36,42,48)

load("../data/translate_geneID_genename.rdata")

rownames(mat) <- as.character(plyr::mapvalues(x=rownames(mat),
                                              from=translation_file_filtered$geneID,
                                              to=translation_file_filtered$gene_symbol,
                                              warn_missing = F))

mat_breaks <- seq(-max(mat), max(mat), length.out = 12)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_rows <- sort_hclust(hclust(dist((mat))))

name = "../plots/Figure_3A_core91_heatmap.pdf"
pdf(name,width=6.6,height=10,paper="special")

pheatmap(mat = mat,
         color  = rev(brewer.pal(length(mat_breaks)-1, "PuOr")),
         breaks = mat_breaks,
         border_color = NA,
         show_colnames= T,
         show_rownames= T,
         cluster_rows = mat_cluster_rows,
         cluster_cols = F,
         drop_levels  = TRUE,
         fontsize     = 12,
         fontsize_row = 8,
         fontsize_col = 12,
         angle_col = 0)

dev.off()




# Figure 3B) plot genes up and down and immune bar plot

library(reshape2)
library(ggplot2)

log_fold_change_short_f <- read.csv("../results/log_fold_change_short.filtered_pvalue-cut_0.01_logFC-cut_2_tp-cut_2_spline_3.csv", row.names = 1)
nrow(log_fold_change_short_f)

immune_genes <- as.character(read.csv("../data/List_of_immune_genes_updated.csv")$FlyBaseID)
# note that this file still had old gene names, using flybaseID instead
length(immune_genes)


summary_table <- matrix(ncol=3, nrow=16)
colnames(summary_table) <- c(' up-regulated   ',' down-regulated',' up-regulated immune genes   ')
rownames(summary_table) <- c("1","2","4","5","6","8","10","12","14","16",
                             "20","24","30","36","42","48")

for(timepoint in 1:16){
  summary_table[,' down-regulated'][timepoint] <- -length(log_fold_change_short_f[timepoint][log_fold_change_short_f[timepoint] < 0])
  
  #  summary_table[,'down-regulated immune genes'][timepoint] <- -length(log_fold_change_short_f[timepoint][log_fold_change_short_f[timepoint] < 0 & rownames(log_fold_change_short_f[timepoint]) %in% immune_genes])
  
  summary_table[,' up-regulated immune genes   '][timepoint] <- length(log_fold_change_short_f[timepoint][log_fold_change_short_f[timepoint] > 0 & rownames(log_fold_change_short_f[timepoint]) %in% immune_genes])
  
  summary_table[,' up-regulated   '][timepoint] <- length(log_fold_change_short_f[timepoint][log_fold_change_short_f[timepoint] > 0]) - length(log_fold_change_short_f[timepoint][log_fold_change_short_f[timepoint] > 0 & rownames(log_fold_change_short_f[timepoint]) %in% immune_genes])
}

melted_table <- melt(summary_table)
melted_table$Var1 <- factor(melted_table$Var1, levels = c("1","2","4","5","6","8","10","12","14","16",
                                                          "20","24","30","36","42","48"))
melted_table$labels <- abs(melted_table$value)
colnames(melted_table) <- c("Var1","Genes","value","labels")

melted_table$labels[melted_table$Genes == ' up-regulated   '] <- melted_table$value[melted_table$Genes == ' up-regulated   '] + melted_table$value[melted_table$Genes == ' up-regulated immune genes   ']


melted_table$Genes <- factor(melted_table$Genes, levels = c(" up-regulated   ",
                                                            " up-regulated immune genes   ",
                                                            " down-regulated"))

name = "../plots/Figure_3B_genedistribution.pdf"
pdf(name,width=9,height=4.5,paper="special")

ggplot(melted_table, aes(x=Var1, y=value, fill=Genes)) + 
  geom_bar(position="stack",width=0.75,stat="identity",color='black',size=0.1) +
  ylab("Number of Genes") +
  xlab("Hours after injection") +
  geom_text(data=melted_table[melted_table$Genes == ' up-regulated immune genes   ',],
            aes(label=melted_table[melted_table$Genes == ' up-regulated immune genes   ',]$labels), 
            vjust=-0.3, size=3.5) +
  geom_text(data=melted_table[melted_table$Genes == ' up-regulated   ',],
            aes(y=labels, label=melted_table[melted_table$Genes == ' up-regulated   ',]$labels), 
            vjust=-0.3, size=3.5) +
  geom_text(data=melted_table[melted_table$Genes == ' down-regulated',],
            aes(label=melted_table[melted_table$Genes == ' down-regulated',]$labels), 
            vjust=1.2, size=3.5) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        axis.title.x=element_text(size=12, margin=margin(7,0,0,0)),
        axis.title.y=element_text(size=12, margin=margin(0,5,0,0)),
        legend.position = "bottom") +
  geom_abline(slope=0, color="black") +
  scale_y_continuous(breaks=c(-10,0,10,30,60)) +
  scale_fill_manual(values=c("#F8766D", "#F8766D", "#619CFF"))

dev.off()

