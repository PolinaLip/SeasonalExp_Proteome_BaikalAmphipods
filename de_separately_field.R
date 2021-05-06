### Data normalization and DE analysis of field samples (seasonal experiment)
library(tidyverse) 
library(limma) 
library(edgeR)
library(ggfortify) 
library(pheatmap)
library(viridis)

species <- 'Eve'

### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}
path2meta <-
  'labeglo2/MS_results/Field/Eve/Metadata_Proteus.tsv'

meta <- meta_upload(path2meta, species)
meta[which(meta$sex == ''),]$sex <- 'NA'

### 2. Data uploading (proteinGroups file from MaxQuant)
dir <- paste0('labeglo2/MS_results/Field/', species)

proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir, proteinGroups_file), sep = '\t', header = T, 
                     check.names = F) 
pep_annot <- read.csv(file.path(dir, 'annot_proteinGroups_field_eve.csv'), sep = '\t',
                      header = T)

select_data <- function(meta_data, proteinGroups_data){
  dat_ <- proteinGroups_data[,c('Protein IDs', meta_data$measure)]
  rownames(dat_) <- dat_$`Protein IDs`
  dat_ <- dat_[-1]
  colnames(dat_) <- meta_data$sample
  dat_[dat_ == 0] <- NA
  return(dat_)
}
dat <- select_data(meta, dat_init)

dat <- dat[,grepl('BK|pool', colnames(dat))]
dat <- dat[,!grepl('pool_12', colnames(dat))]
meta <- meta[!grepl('PB|pool_12', meta$sample),]

### 3. Get rid of the outliers and samples, which you do not want to analyze
dat <- dat[,!grepl('F/100/BK/3c/2m_2', colnames(dat))]
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')

### 4. Remove rows with all NAs in at least one of the conditions

remove_completeNAs_inCond <- function(data2clean, metafile){
  cond <- as.factor(metafile$condition)
  ixs <- order(cond)
  data2clean <- data2clean[, ixs]
  cond <- cond[ixs]
  
  for (level in levels(cond)) {
    col_ixs <- cond == level
    data2clean <- data2clean[rowSums(is.na(data2clean[,col_ixs])) < sum(col_ixs),]
  }
  return(data2clean)
}
data_raw <- remove_completeNAs_inCond(dat, meta)

### 5. Visualize the raw data
boxplot(log2(data_raw), 
        notch = TRUE, main = 'RAW data',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_raw), 
              main = 'Raw data', legend=F)

### 6. Sample loading normalization

# Eve
exp1_raw <- data_raw[,subset(meta, experiment == 1)$sample]
exp2_raw <- data_raw[,subset(meta, experiment == 2)$sample]
exp3_raw <- data_raw[,subset(meta, experiment == 3)$sample]
exp4_raw <- data_raw[,subset(meta, experiment == 4)$sample]
exp5_raw <- data_raw[,subset(meta, experiment == 5)$sample]
exp6_raw <- data_raw[,subset(meta, experiment == 6)$sample]
exp7_raw <- data_raw[,subset(meta, experiment == 7)$sample]
exp8_raw <- data_raw[,subset(meta, experiment == 8)$sample]
exp9_raw <- data_raw[,subset(meta, experiment == 9)$sample]
exp10_raw <- data_raw[,subset(meta, experiment == 10)$sample]
exp11_raw <- data_raw[,subset(meta, experiment == 11)$sample]
#exp12_raw <- data_raw[,subset(meta, experiment == 12)$sample]

target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
                 colSums(exp3_raw, na.rm = T), colSums(exp4_raw, na.rm = T),
                 colSums(exp5_raw, na.rm = T), colSums(exp6_raw, na.rm = T),
                 colSums(exp7_raw, na.rm = T), colSums(exp8_raw, na.rm = T),
                 colSums(exp9_raw, na.rm = T), colSums(exp10_raw, na.rm = T),
                 colSums(exp11_raw, na.rm = T)),
               na.rm = T)

norm_facs <- target / colSums(exp1_raw, na.rm = T)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw, na.rm = T)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw, na.rm = T)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp4_raw, na.rm = T)
exp4_sl <- sweep(exp4_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp5_raw, na.rm = T)
exp5_sl <- sweep(exp5_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp6_raw, na.rm = T)
exp6_sl <- sweep(exp6_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp7_raw, na.rm = T)
exp7_sl <- sweep(exp7_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp8_raw, na.rm = T)
exp8_sl <- sweep(exp8_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp9_raw, na.rm = T)
exp9_sl <- sweep(exp9_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp10_raw, na.rm = T)
exp10_sl <- sweep(exp10_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp11_raw, na.rm = T)
exp11_sl <- sweep(exp11_raw, 2, norm_facs, FUN = "*")
#norm_facs <- target / colSums(exp12_raw, na.rm = T)
#exp12_sl <- sweep(exp12_raw, 2, norm_facs, FUN = "*")

data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl, exp4_sl, exp5_sl,
                 exp6_sl, exp7_sl, exp8_sl, exp9_sl, exp10_sl,
                 exp11_sl)

boxplot(log2(data_sl), 
        notch = TRUE, main = "Sample Loading (SL) normalized data",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_sl),
              main = "SL normalization", legend = F) # Eve, resolution -> 800/553

### 7. Vizualize normalized data
## 7.1 MDS-plot with groups colored by TMT batch
plotMDS(log2(data_sl), col = c(rep(c('red', 'blue', 'green', 'yellow', 'black',
                                     'pink', 'purple', 'orange', 'grey',
                                     'brown', 'aquamarine', 'cyan'), each = 10)), 
        main = "SL clusters grouped by TMT experiment") # Eve

## 7.2 MDS-plot with groups colored by condition:
col_vec <- ifelse(grepl('pool', colnames(data_sl)), 'darkgreen',
                  ifelse(grepl('BK', colnames(data_sl)), 'orange', 'red'))

plotMDS(log2(data_sl), col = col_vec)

col_vec <- ifelse(grepl('pool', colnames(data_sl)), 'black',
                  ifelse(grepl('/1/', colnames(data_sl)), 'orange', 
                         ifelse(grepl('/2/', colnames(data_sl)), 'blue', 
                                ifelse(grepl('/3c/', colnames(data_sl)), 'red',
                                ifelse(grepl('/4/', colnames(data_sl)), 'darkgreen', 'purple')))))

plotMDS(log2(data_sl), col = col_vec)

## 7.3 PCA-plot
data_sl_t <- na.omit(data_sl) %>% t %>% as.data.frame
pca_res <- prcomp(data_sl_t)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
autoplot(pca_res, data=meta2, colour='condition', size = 3, shape='sex') + theme_light()
# autoplot(pca_res, data=meta2, colour='condition', label = TRUE, label.size = 5) + theme_light()
# meta2[rownames(meta2) == 44,]$sample
# ! outlier found: F/100/BK/3c/2m_2
ggsave(filename = file.path(dir, 'pca_eve_BK_seasonalExp.png'), scale = 2)

### 8. Multiply the normalized data by 10^7 to make it possible to perform DE analysis
data_sl_mult <- data_sl * 10000000 # for PSM-normalized data
data_irs <- data_sl_mult # for PSM-normalized data
data_sl_mult$protein_group <- row.names(data_sl_mult)
data_sl_mult_ <- data.frame(protein_group = data_sl_mult$protein_group, 
                            data_sl_mult[1:length(data_sl_mult)-1])
write.table(data_sl_mult_, file.path(dir, 'intensities_after_slNorm_eve.csv'), 
            sep = '\t', quote = F, row.names = F)

#########
# EdgeR 
#########

### 9. Separate the data on the subset with the same NA pattern and perform EdgeR analysis

# here are two functions to perform de analysis
edger_analysis <- function(datatode, group, sex_group, cond_vector, 
                           protein_annotation){
  
  edge_object <- DGEList(counts = datatode, group = group)
#  if (length(unique(sex_group)) > 1){
#    design <- model.matrix(~0 + group + sex_group)
#    print("Information concerning sex of the animal was used in design")
#  }
#  else {
#    print("Information concerning sex of the animal is not using in design")
#    design <- model.matrix(~0 + group)
#  }
  design <- model.matrix(~0 + group)
  edge_object <- estimateDisp(edge_object, design)
  plotBCV(edge_object, main = "Biological variation SL only (no IRS)")
  edge_res <- exactTest(edge_object, pair = cond_vector)
  
  print(summary(decideTestsDGE(edge_res)))
  
  tt_edge_res <- topTags(edge_res, n = Inf, sort.by = "none")
  tt_edge_res <- tt_edge_res$table
  tt_edge_res_sign <- subset(tt_edge_res, PValue < 0.05)
  
  proteins <- rownames(tt_edge_res)
  samples2annot <- match(proteins, protein_annotation$protein_group)
  sum(is.na(samples2annot))
  
  SignProteinInfo <- data.frame(protein = proteins,
                                geneSymbol = protein_annotation$annotation[samples2annot],
                                logFC = tt_edge_res$logFC,
                                pvalue = tt_edge_res$PValue,
                                FDR = tt_edge_res$FDR)
  
  SignProteinInfo <- SignProteinInfo[order(SignProteinInfo$FDR),]
  return(SignProteinInfo)
}

edger_analysis_wrap <- function(norm_data, protein_annotation, 
                                cond2compare, metafile){
  cond2compare_ <- paste(cond2compare, collapse = '|') # to make it readable for grepl
  samples2take <- metafile[grepl(cond2compare_, metafile$condition),]$sample
  samples2take <- paste(samples2take, collapse = '|')
  data2DE <- norm_data[,grepl(samples2take, colnames(norm_data))] 
  
  sample_order <- match(colnames(data2DE), metafile$sample)
  group <- as.factor(metafile$condition[sample_order])
  
  group_sex <- as.factor(metafile$sex[sample_order])
  
  print("Samples to be analyzed:")
  print(colnames(data2DE))
  
  return(edger_analysis(data2DE, group, group_sex, cond2compare, 
                        protein_annotation))
}

# wo NAs
data_wo_na <- na.omit(data_irs)
SignProteinInfo_0 <- edger_analysis_wrap(norm_data = data_wo_na, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("September_BK", 
                                                          "November_BK"), 
                                         metafile = meta)
### Heatmap for all genes

# Prepare intensities data: scale the data -> transport it back 
intensities <- data_wo_na
scaled_intensities <- apply(intensities, 1, scale)
scaled_intensities <- t(scaled_intensities)
colnames(scaled_intensities) <- meta$sample

# Prepare the annotation table:
annot_col <- data.frame(Months = factor(meta$condition),
                        Sex = factor(meta$sex))
rownames(annot_col) <- colnames(scaled_intensities)
annot_colours <- list(Months = c("December_BK" = "dodgerblue1",
                                 "January_BK" = "cyan1",
                                 "June_BK" = "goldenrod1",
                                 "November_BK" = "brown",
                                 "pool" = "black",
                                 "September_BK" = "darkorange1"),
                      Sex = c("female" = "coral1",
                              "male" = "aquamarine4"))
# Draw the heatmap
pheatmap(scaled_intensities,
         color = plasma(10),
         show_rownames = F,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         annotation_col = annot_col,
         annotation_colors = annot_colours,
         height = 10,
         width = 15,
         filename = file.path(dir, 
                              paste0(species, 
                                     '_heatmap_genes_all', '.png')))
if (dev.cur() != 1) {dev.off()}

# clustering of seasonal profiles and visualizing temporal profiles of proteins



