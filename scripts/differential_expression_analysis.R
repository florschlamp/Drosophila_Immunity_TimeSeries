# Scripts for differential expression analyses
# spline fitting and pairwise differential expression
# last updated: March 2nd, 2021
# Florencia Schlamp
# mfs97@cornell.edu

library(limma)
library(edgeR)
library(splines)
library(DESeq2)

# Functions:

LRT_limmavoom <- function(
  dat            # A matrix of (normalized) counts, dimension: no. genes x no. timepts
  , time_ID        # a numeric vector of time pts
  , gene_name     # a character vector of gene names
  , design_type  # type of curve and degree of freedom as a vector, choices: c('spline', *) or c('polynomial', *) 
  , normalize_method = 'TMM'     # set it TMM for now, will add other variations later
  , pval_cutoff = 0.05            
){
  no_gene = nrow(dat); no_timept = ncol(dat)
  if (length(gene_name)!=no_gene){stop("length of gene names do not match data size")}
  if (length(time_ID)!=no_timept){stop("length of timepts does not match data size")}
  
  rownames(dat) = gene_name; colnames(dat) = time_ID
  
  # construct design matrix
  if (design_type[1]=='spline'){
    df = as.numeric(design_type[2])
    design = model.matrix(~ns(time_ID, df))
  }
  if (design_type[1]=='polynomial'){
    df = as.numeric(design_type[2])
    design = model.matrix(~poly(time_ID, df))
  }
  
  # apply normalization
  if (normalize_method == "TMM"){
    dat = DGEList(counts = dat)
    dat = calcNormFactors(dat)
  }
  
  # use voom transform for precision weights
  voomdat = voom(dat, design=design, normalize.method='none', plot=T)
  
  # fit linear models to univariate time course trajectories
  fit = lmFit(voomdat, design = design, weights=NULL)
  dat_smoothed = fit$coefficients %*% t(fit$design)
  fit = eBayes(fit, trend=TRUE)
  top = topTable(fit, coef=seq(2, (df+1)), sort.by="B", number=nrow(dat), p.value=pval_cutoff)
  
  return(list(top=top, dat_smoothed = dat_smoothed))
}


DE_timecourse = function(
  dat1, dat2       # Two matrices of (normalized) counts, each with dimension: no. genes x no. timepts, make sure the gene names are ordered alphabetically (note to sumanta: see if you can change within the function)
  , time_ID          # a numeric vector of time points
  , gene_name        # a character vector of gene names
  , design_type # type of curve and degree of freedom as a vector, choices: c('spline', *) or c('polynomial', *) 
  , normalize_method = 'TMM'     # set it TMM for now, will add other variations later
  , save_output = FALSE
  , out_fname = 'DE_timecourse.RData'
){
  no_gene = nrow(dat1); no_timept = ncol(dat1); T = no_timept-1; p = no_gene
  if (length(gene_name)!=no_gene){stop("length of gene names do not match data size")}
  if (length(time_ID)!=no_timept){stop("length of timepts does not match data size")}
  if (any(dim(dat1)!=dim(dat2))){stop("two replicate datasets have different size")}
  
  # construct design matrix
  if (design_type[1]=='spline'){
    df = as.numeric(design_type[2])
    design = model.matrix(~ns(time_ID, df))
  }
  if (design_type[1]=='polynomial'){
    df = as.numeric(design_type[2])
    design = model.matrix(~poly(time_ID, df))
  }
  
  # apply normalization
  if (normalize_method == "TMM"){
    dat1 = DGEList(counts = dat1); dat1 = calcNormFactors(dat1)
    dat2 = DGEList(counts = dat2); dat2 = calcNormFactors(dat2)
  }
  
  # use voom transform for precision weights
  voomdat1 = voom(dat1, design=design, normalize.method='none', plot=T)
  voomdat2 = voom(dat2, design=design, normalize.method='none', plot=T)
  dat1.voom = voomdat1$E; dat2.voom = voomdat2$E
  
  
  # combine replicates in a single dataset
  dat = cbind(dat1.voom, dat2.voom)
  rownames(dat) = gene_name; colnames(dat) = rep(time_ID, 2)
  
  # fit linear model to aggregate two replicates at each time point
  des_contrast = rbind(diag(T+1), diag(T+1))
  fit = lmFit(dat, des_contrast)
  
  # construct contrast matrix for DE testing (time t - time 0)
  cont = rbind(matrix(rep(-1, T), nrow=1, ncol=T)
               , diag(T))
  ffc = contrasts.fit(fit, cont)
  ffceb = eBayes(ffc)
  
  result = list()
  pval = array(1, c(p, T))
  rownames(pval) = gene_name
  colnames(pval) = paste("t",time_ID[-1], ":t0", sep="")
  logFC = array(0, c(p, T))
  rownames(logFC) = gene_name
  colnames(logFC) = paste("t", time_ID[-1], ":t0", sep="")
  
  
  for (j in 1:T){
    tmp = topTable(ffceb, coef = j, number = p, adjust = "none")
    tmp$gene = rownames(tmp)
    tmp = tmp[order(tmp$gene),]
    #print(head(tmp))
    #print(head(pval))
    
    if (any(rownames(tmp)!=rownames(pval)))
      print(paste('orders of genes do not match for coef ', j))
    result[[j]] = tmp
    pval[,j] <- tmp$adj.P.Val
    logFC[,j] <- tmp$logFC
  }
  
  if (save_output)
    save(pval, logFC, result, file = out_fname)
  
  return(list(pvalue = pval, logFoldChange = logFC))  
  
}

filter_genes <- function(
  pvalue                   # a matrix of size #genes x (#timepts-1), as in the output of DE_timecourse
  , logFC            # a matrix of size #genes x (#timepts-1), as in the output of DE_timecourse
  , p_adjust_method = 'none' # options: 'bonferroni', 'BH', 'BY', see p.adjust() for all available options
  , p_cut
  , logFC_cut
  , rule = 'and'             # options: 'and', 'or'
  , t_cut = 1            # only keep genes which are differential in atleast t_cut timepts
){
  
  if (any(dim(pvalue)!=dim(logFC))){stop("dimensions of pval and logFC do not match!!")}
  out <- logFC
  
  # correct for multiple testing
  if (p_adjust_method!='none'){
    tmp = matrix(p.adjust(as.vector(pvalue), method=p_adjust_method)
                 , nrow(pvalue), ncol(pvalue)
    )
    rownames(tmp) = rownames(pvalue)
    colnames(tmp) = colnames(pvalue)
    pvalue <- tmp
    rm(tmp)			   
  }
  
  # set all logFC with high pval and/or low logFC to zero 
  if (rule == 'and')
    out[!((pvalue <= p_cut) & (abs(logFC) >= logFC_cut))] = 0
  else
    out[!((pvalue <= p_cut) | (abs(logFC) >= logFC_cut))] = 0
  
  # keep only genes differentially expressed in atleast t_cut  points
  keep_genes_idx = (rowSums(!!out) >= t_cut)
  
  return(list(adj.pvalue = pvalue[keep_genes_idx,], logFoldChange = out[keep_genes_idx,]))
}


#### ANALYSIS:
## Step 1: read in data
# A matrix of (normalized) counts, dimension: no. genes x no. timepts
dat_norm <- read.table('../data/normalized_counts.csv',
                       header=TRUE, row.names=1, as.is=TRUE, 
                       check.names=FALSE, sep=',')
dat_norm_log_filt <- read.table('../data/normalized_counts_log_filt.csv',
                                header=TRUE, row.names=1, as.is=TRUE, 
                                check.names=FALSE, sep=',')
dat_short <- subset(dat_norm, select = -c(4,19,20,21,39,40,41))
dat_short_filtered <- dat_short[rownames(dat_short) %in% rownames(dat_norm_log_filt),]

dat_short = data.matrix(dat_short_filtered)
dim(dat_short)
# 12,657 genes in 34 samples (first 48 hours)

sample_ID_short = colnames(dat_short_filtered)
replicate_ID_short = substr(sample_ID_short, start=nchar(sample_ID_short), 
                            stop = nchar(sample_ID_short))
time_ID_short = substr(sample_ID_short, start=1, stop=nchar(sample_ID_short)-1)

hr_short = c(0,1,2,3,4,5,6,8,10,12,14,16,20,24,30,36,42,48)

time_ID_short = hr_short[as.numeric(time_ID_short)]

dat_short_A = dat_short_filtered[,replicate_ID_short == "A"]
dat_short_B = dat_short_filtered[,replicate_ID_short == "B"]

# a numeric vector of time pts
time_ID_short_A = time_ID_short[replicate_ID_short == "A"]
time_ID_short_B = time_ID_short[replicate_ID_short == "B"]

# a character vector of gene names
gene_name <- rownames(dat_short_filtered) 

# type of curve and degree of freedom as a vector
design_type = c('spline', 3) 

# set it TMM for now, will add other variations later
normalize_method = 'TMM'     

# and P-value
pval_cutoff = 0.05


### Step 2: spline analysis (LRT_limmavoom)
# ============== replicate A =============
#window analysis (first 48 hours)
results_short_A <- LRT_limmavoom(dat_short_A,time_ID_short_A,
                                 gene_name,design_type,
                                 normalize_method,pval_cutoff)
top_genes_short_A <- as.data.frame(results_short_A['top'])
nrow(top_genes_short_A)
# 1163 genes at 0.05 
# note: some versions give 1164 genes, but it doesn't change the downstream overlap

# save results
file_name <- "../results/top_genes_limma-voom_short.filtered_A_spline_3_TMM_0.05.csv"
write.csv(top_genes_short_A,file=file_name)

smooth_short_A <- as.data.frame(results_short_A['dat_smoothed'])
nrow(smooth_short_A)
file_name <- "../results/top_genes_limma-voom_short.filtered_A_spline_3_TMM_ALL_SMOOTHED.csv"
write.csv(smooth_short_A,file=file_name)

#window analysis (first 8 hours)
dat_short_A_subset <- dat_short_A[1:7]
time_ID_short_A_subset <- time_ID_short_A[1:7]

results_short_A_window <- LRT_limmavoom(dat_short_A_subset,time_ID_short_A_subset,
                                        gene_name,design_type,
                                        normalize_method,pval_cutoff)
top_genes_short_A_window <- as.data.frame(results_short_A_window['top'])
nrow(top_genes_short_A_window)
# 129 (1to7) at 0.05

# save results
file_name <- "../results/top_genes_limma-voom_short.filtered_A_spline_3_TMM_0.05_1to7.csv"
write.csv(top_genes_short_A_window,file=file_name)

smooth_short_A_window <- as.data.frame(results_short_A_window['dat_smoothed'])
nrow(smooth_short_A_window)
file_name <- "../results/top_genes_limma-voom_short.filtered_A_spline_3_TMM_1to7_ALL_SMOOTHED.csv"
write.csv(smooth_short_A_window,file=file_name)


# ============== replicate B =============
#window analysis (first 48 hours)
results_short_B <- LRT_limmavoom(dat_short_B,time_ID_short_B,
                                 gene_name,design_type,
                                 normalize_method,pval_cutoff)
top_genes_short_B <- as.data.frame(results_short_B['top'])
nrow(top_genes_short_B)
# 845 at 0.05 

# save results
file_name <- "../results/top_genes_limma-voom_short.filtered_B_spline_3_TMM_0.05.csv"
write.csv(top_genes_short_B,file=file_name)

smooth_short_B <- as.data.frame(results_short_B['dat_smoothed'])
nrow(smooth_short_B)
file_name <- "../results/top_genes_limma-voom_short.filtered_B_spline_3_TMM_ALL_SMOOTHED.csv"
write.csv(smooth_short_B,file=file_name)

#windows analysis (first 8 hours)
dat_short_B_subset <- dat_short_B[1:7]
time_ID_short_B_subset <- time_ID_short_B[1:7]

results_short_B_window <- LRT_limmavoom(dat_short_B_subset,time_ID_short_B_subset,
                                        gene_name,design_type,
                                        normalize_method,pval_cutoff)
top_genes_short_B_window <- as.data.frame(results_short_B_window['top'])
nrow(top_genes_short_B_window)
# 145 (1to7) at 0.05

# save results
file_name <- "../results/top_genes_limma-voom_short.filtered_B_spline_3_TMM_0.05_1to7.csv"
write.csv(top_genes_short_B_window,file=file_name)

smooth_short_B_window <- as.data.frame(results_short_B_window['dat_smoothed'])
nrow(smooth_short_B_window)
file_name <- "../results/top_genes_limma-voom_short.filtered_B_spline_3_TMM_1to7_ALL_SMOOTHED.csv"
write.csv(smooth_short_B_window,file=file_name)

#### SUMMARY
##### what is the overlap between spline modeling results for repA and repB #####

## all data (first 48 hours)
union_top_genes_short <- union(rownames(top_genes_short_A),rownames(top_genes_short_B))
length(union_top_genes_short)
# 1628 genes identified by either repA and/or repB
intersect_top_genes_short <- intersect(rownames(top_genes_short_A),rownames(top_genes_short_B))
length(intersect_top_genes_short)
# 380 genes identified by both repA and repB

## smaller window (first 8 hours)
union_top_genes_short_window <- union(rownames(top_genes_short_B_window),rownames(top_genes_short_A_window))
length(union_top_genes_short_window)
# 226 genes identified by either repA and/or repB
intersect_top_genes_short_window <- intersect(rownames(top_genes_short_B_window),rownames(top_genes_short_A_window))
length(intersect_top_genes_short_window)
# 48 genes identified by both repA and rep

### FINAL SUBSET FROM SPLINE MODELING ####
# union of the intersect of replicates for top genes and top genes window
union_of_intersections_limmavoom <- union(intersect_top_genes_short,intersect_top_genes_short_window)
length(union_of_intersections_limmavoom) 
# 411 genes

union_of_unions_limmavoom <- union(union_top_genes_short,union_top_genes_short_window)
length(union_of_unions_limmavoom) 
# 1724 genes


### Step 3: pairwise differential expression (and filtering)
## differential expression analysis and subsetting of genes
## pairwise comparison of each tp against t0
###### ========== DE Analysis =============
dat_short_A = dat_short_A[order(rownames(dat_short_A)),]
dat_short_B = dat_short_B[order(rownames(dat_short_B)),]

out_short = DE_timecourse(dat_short_A, 
                          dat_short_B, 
                          gene_name = rownames(dat_short_A), 
                          time_ID = time_ID_short_B,
                          design_type = design_type)

p_values_short <- as.data.frame(out_short['pvalue'])
log_fold_change_short <- as.data.frame(out_short['logFoldChange'])

## save full tables
file_name = paste("../results/adj.pvalues_short.filtered_ALL_spline_3.csv",sep="")
write.csv(p_values_short,file=file_name)

file_name = paste("../results/log_fold_change_short.filtered_ALL_spline_3.csv",sep="")
write.csv(log_fold_change_short,file=file_name)



# filter results (p-value 0.05, logFC=2, tp=1, total 214 genes)
desired_p_value_cutoff = 0.05 
desired_logFC_cutoff = 2
desired_number_of_timepoints = 1

out_short_f = filter_genes(pvalue = out_short$pvalue, 
                           logFC = out_short$logFoldChange, 
                           p_cut = desired_p_value_cutoff, 
                           logFC_cut = desired_logFC_cutoff, 
                           t_cut = desired_number_of_timepoints,
                           rule = 'and',
                           p_adjust_method = 'BH')

adj.p_values_short_f <- as.data.frame(out_short_f['adj.pvalue'])
log_fold_change_short_f <- as.data.frame(out_short_f['logFoldChange'])
nrow(adj.p_values_short_f)
nrow(log_fold_change_short_f) # 214 genes

DE.214 <- rownames(log_fold_change_short_f)

file_name = paste("../results/adj.pvalues_short.filtered_pvalue-cut_",desired_p_value_cutoff,
                  "_logFC-cut_",desired_logFC_cutoff,
                  "_tp-cut_",desired_number_of_timepoints,
                  "_",design_type[1],"_",design_type[2],".csv",sep="")
write.csv(adj.p_values_short_f,file=file_name)

file_name = paste("../results/log_fold_change_short.filtered_pvalue-cut_",desired_p_value_cutoff,
                  "_logFC-cut_",desired_logFC_cutoff,
                  "_tp-cut_",desired_number_of_timepoints,
                  "_",design_type[1],"_",design_type[2],".csv",sep="")
write.csv(log_fold_change_short_f,file=file_name)


### FINAL SUBSET ######
length(union(union_of_intersections_limmavoom,DE.214))
#551
SUBSET_551 <- union(union_of_intersections_limmavoom,DE.214)
################################





# filter results to core 91 (p-value 0.01, logFC=2, tp=2, total 91 genes)
desired_p_value_cutoff = 0.01
desired_logFC_cutoff = 2
desired_number_of_timepoints = 2

out_short_f = filter_genes(pvalue = out_short$pvalue, 
                           logFC = out_short$logFoldChange, 
                           p_cut = desired_p_value_cutoff, 
                           logFC_cut = desired_logFC_cutoff, 
                           t_cut = desired_number_of_timepoints,
                           rule = 'and',
                           p_adjust_method = 'BH')

adj.p_values_short_f <- as.data.frame(out_short_f['adj.pvalue'])
log_fold_change_short_f <- as.data.frame(out_short_f['logFoldChange'])
nrow(adj.p_values_short_f)
nrow(log_fold_change_short_f) # 91

core_91 <- rownames(log_fold_change_short_f)

file_name = paste("../results/adj.pvalues_short.filtered_pvalue-cut_",desired_p_value_cutoff,
                  "_logFC-cut_",desired_logFC_cutoff,
                  "_tp-cut_",desired_number_of_timepoints,
                  "_",design_type[1],"_",design_type[2],".csv",sep="")
write.csv(adj.p_values_short_f,file=file_name)

file_name = paste("../results/log_fold_change_short.filtered_pvalue-cut_",desired_p_value_cutoff,
                  "_logFC-cut_",desired_logFC_cutoff,
                  "_tp-cut_",desired_number_of_timepoints,
                  "_",design_type[1],"_",design_type[2],".csv",sep="")
write.csv(log_fold_change_short_f,file=file_name)


# filter results to 729 (pval 0.05, logFC 1)
desired_p_value_cutoff = 0.05 
desired_logFC_cutoff = 1
desired_number_of_timepoints = 1

out_short_f = filter_genes(pvalue = out_short$pvalue, 
                           logFC = out_short$logFoldChange, 
                           p_cut = desired_p_value_cutoff, 
                           logFC_cut = desired_logFC_cutoff, 
                           t_cut = desired_number_of_timepoints,
                           rule = 'and',
                           p_adjust_method = 'BH')

adj.p_values_short_f <- as.data.frame(out_short_f['adj.pvalue'])
log_fold_change_short_f <- as.data.frame(out_short_f['logFoldChange'])
nrow(adj.p_values_short_f)
nrow(log_fold_change_short_f) 
# 729

DE.729 <- rownames(log_fold_change_short_f)

file_name = paste("../results/adj.pvalues_short.filtered_pvalue-cut_",desired_p_value_cutoff,
                  "_logFC-cut_",desired_logFC_cutoff,
                  "_tp-cut_",desired_number_of_timepoints,
                  "_",design_type[1],"_",design_type[2],".csv",sep="")
write.csv(adj.p_values_short_f,file=file_name)

file_name = paste("../results/log_fold_change_short.filtered_pvalue-cut_",desired_p_value_cutoff,
                  "_logFC-cut_",desired_logFC_cutoff,
                  "_tp-cut_",desired_number_of_timepoints,
                  "_",design_type[1],"_",design_type[2],".csv",sep="")
write.csv(log_fold_change_short_f,file=file_name)



## save all gene list subsets
# give more descriptive names to some subsets
spline48_380genes <- intersect_top_genes_short
spline8_48genes <- intersect_top_genes_short_window
total_splines_411genes <- union_of_intersections_limmavoom
total_DEGs_951genes <- length(union(DE.729,SUBSET_551))

save(spline48_380genes,file="../results/gene_lists/spline_48h_pval0.05_380genes.rdata")
save(spline8_48genes,file="../results/gene_lists/spline_8h_pval0.05_48genes.rdata")
save(total_splines_411genes,file="../results/gene_lists/both_splines_pval0.05_411genes.rdata")
save(DE.214, file="../results/gene_lists/pairwiseDE_pval0.05_logFC2_1tp_214genes.rdata")
save(DE.729, file="../results/gene_lists/pairwiseDE_pval0.05_logFC1_1tp_729genes.rdata")
save(core_91, file="../results/gene_lists/pairwiseDE_pval0.01_logFC2_2tp_91genes.rdata")
save(total_DEGs_951genes, file="../results/gene_lists/all_DEGs_951genes.rdata")


