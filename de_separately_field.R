### Data normalization and DE analysis of field samples (seasonal experiment)
library(tidyverse) 
library(limma) 
library(edgeR)
library(ggfortify) 
library(pheatmap)
library(viridis)

species <- 'Eve'
population <- 'BK'
#dir_to_results <- paste0('labeglo2/MS_results/Field/', species, '/', population)
dir_to_results <- paste0('labeglo2/MS_results/Field/', species, '/', species, 
                         '_refchannels_all')
#dir_to_results <- paste0('labeglo2/MS_results/Field/', species, '/', species, 
#                         '_refchannels_one')
### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}

path2meta <- paste0('labeglo2/MS_results/Field/', species,
                    '/', species, '_refchannels_all', '/Metadata_Proteus.tsv')

meta <- meta_upload(path2meta, species)
meta[which(meta$sex == ''),]$sex <- 'NA'

### 2. Data uploading (proteinGroups file from MaxQuant)
#dir_to_data <- paste0('labeglo2/MS_results/Field/', species)
dir_to_data <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species, '_refchannels_all')
#dir_to_data <- paste0('labeglo2/MS_results/Field/', species, 
#                      '/', species, '_refchannels_one')

proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir_to_data, proteinGroups_file), sep = '\t', 
                     header = T, 
                     check.names = F) 
pep_annot <- read.csv(file.path(dir_to_data, 
                                paste0('annot_proteinGroups_field_', 
                                       tolower(species), '.csv')), 
                      sep = '\t',
                      header = T)

# or you can use updated annotation
pep_annot <- read.csv(file.path(dir_to_data, 
                                  paste0('annot_proteinGroups_field_', tolower(species), 
                                         '_withGOannotation.csv')), 
                        sep = '\t', header = T)

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

# to select PB population
dat <- dat[,grepl('PB|pool', colnames(dat))]
meta <- meta[!grepl('BK', meta$sample),]

### 3. Get rid of the outliers and samples, which you do not want to analyze

dat <- dat[,!grepl('F/100/BK/3c/2m_2', colnames(dat))]
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')

dat <- dat[, !grepl('pool', colnames(dat))] # for Gla
meta <- subset(meta, !grepl('pool', sample)) # for Gla

#dat <- dat[,!grepl('F/104/L/4/5p_1|F/104/L/4/4p_1', colnames(dat))] # it seems that it is not an outlier
#meta <- subset(meta, !grepl('F/104/L/4/5p_1|F/104/L/4/4p_1', sample)) # it seems that it is not an outlier

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

# if use data without remove_completeNAs_inCond
data_raw <- dat
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
exp5_raw <- data_raw[,subset(meta, experiment == 5)$sample] # until here for Gla
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
                 colSums(exp11_raw, na.rm = T) 
                 #,colSums(exp12_raw, na.rm = T)
                 ),
               na.rm = T) # Eve, BK
target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
                   colSums(exp3_raw, na.rm = T), colSums(exp4_raw, na.rm = T),
                   colSums(exp5_raw, na.rm = T), colSums(exp6_raw, na.rm = T),
                   colSums(exp7_raw, na.rm = T)),
                   na.rm = T) # Ecy
target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
                 colSums(exp3_raw, na.rm = T), colSums(exp4_raw, na.rm = T),
                 colSums(exp5_raw, na.rm = T)),
                 na.rm = T) # Gla
target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
                 colSums(exp3_raw, na.rm = T), colSums(exp4_raw, na.rm = T)),
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
                 exp11_sl
                 #, exp12_sl
                 ) # Eve, BK
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl, exp4_sl, exp5_sl,
                 exp6_sl, exp7_sl) # Ecy
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl, exp4_sl, exp5_sl) # Gla
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl, exp4_sl)

boxplot(log2(data_sl), 
        notch = TRUE, main = "Sample Loading (SL) normalized data",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_sl),
              main = "SL normalization", legend = F) # Eve, resolution -> 800/553

####### For Gla: to normalize 5th pool:

med14 <- median((unlist(c(exp1_sl, exp2_sl, exp3_sl, exp4_sl))), na.rm = T)
med5 <- median((unlist(exp5_sl)), na.rm = T)
mult_factor <- med14 / med5
median(unlist(exp5_sl * mult_factor), na.rm = T)
exp5_sl_norm <- sweep(exp5_sl, 2, mult_factor, FUN = "*")

data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl, exp4_sl, exp5_sl_norm)

### 7. Vizualize normalized data
## 7.1 MDS-plot with groups colored by TMT batch
plotMDS(log2(data_sl), col = c(rep(c('red', 'blue', 'green', 'yellow', 'black',
                                     'pink', 'purple', 'orange', 'grey',
                                     'brown', 'aquamarine', 'cyan', 'darkgreen'), 
                                   each = 10)), 
        main = "SL clusters grouped by TMT experiment") # Eve

## 7.2 MDS-plot with groups colored by condition:
col_vec <- ifelse(grepl('pool', colnames(data_sl)), 'darkgreen',
                  ifelse(grepl('BK|PB', colnames(data_sl)), 'orange', 'red'))

plotMDS(log2(data_sl), col = col_vec)

col_vec <- ifelse(grepl('pool', colnames(data_sl)), 'black',
                  ifelse(grepl('/1/', colnames(data_sl)), 'orange', 
                         ifelse(grepl('/2/', colnames(data_sl)), 'blue', 
                         ifelse(grepl('/3c/', colnames(data_sl), ignore.case = T), 'red',
                         ifelse(grepl('/4/', colnames(data_sl)), 'darkgreen', 'purple')))))

plotMDS(log2(data_sl), col = col_vec)

## 7.3 PCA-plot
library(ggthemes)

data_sl_t <- na.omit(data_sl) %>% t %>% as.data.frame
#data_sl_t <- data_sl_t[!grepl('pool', rownames(data_sl_t)),] # to draw wo pools

pca_res <- prcomp(data_sl_t)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
meta2$condition <- sub('_BK', '', meta2$condition)
#meta2 <- subset(meta2, sex != 'NA') # only for the figure with samples with known sex
#data_sl_t <- data_sl_t[meta2$sample,] # only for the figure with samples with known sex 
#meta2 <- subset(meta2, !grepl('pool', sample)) # to draw wo pools
autoplot(pca_res, data=meta2, colour='condition', size = 3, shape='sex') + theme_light()
autoplot(pca_res, data=meta2, colour='condition', size = 3, 
         #shape='sex', 
         frame = T, 
         x = 1,
         y = 2) + 
  theme_light() + 
  scale_color_manual('Condition', 
                     values = palette.colors(n=6, 'Dark2')) +
  scale_fill_manual('Condition', 
                     values = palette.colors(n=6, 'Dark2')) +
  scale_shape_manual('Sex', values = c(16, 17, 10))

# autoplot(pca_res, data=meta2, colour='condition', label = TRUE, label.size = 5) + theme_light()
# meta2[rownames(meta2) == 44,]$sample
# ! outlier found: F/100/BK/3c/2m_2
ggsave(filename = file.path(dir_to_results, paste0('pca_', 
       species,'_seasonalExp_withPools_ONLYSex.png')), 
       scale = 1,
       width = 6, height = 3)

### 8. Multiply the normalized data by 10^7 to make it possible to perform DE analysis
data_sl_mult <- data_sl * 10000000 # for PSM-normalized data
data_irs <- data_sl_mult # for PSM-normalized data
data_sl_mult$protein_group <- row.names(data_sl_mult)
data_sl_mult_ <- data.frame(protein_group = data_sl_mult$protein_group, 
                            data_sl_mult[1:length(data_sl_mult)-1])
write.table(data_sl_mult_, file.path(dir_to_results, 
                                     paste0('intensities_after_slNorm_', 
                                            tolower(species), '.csv')), 
            sep = '\t', quote = F, row.names = F)

data_sl_ <- read.table(file.path(dir_to_results, 
                     paste0('intensities_after_slNorm_', 
                            tolower(species), '.csv')), header = T)

rownames(data_sl_) <- data_sl_$protein_group
data_sl_ <- data_sl_[-1]
data_irs <- data_sl_
colnames(data_irs) <- meta2$sample
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
                                #geneSymbol = protein_annotation$annotation[samples2annot],
                                geneSymbol = protein_annotation$upd_full_annot[samples2annot], # if use updated annotation
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
sign_protein_info_sub_allcomb <- NULL

data_only_dec_and_sep <- data_irs[,subset(meta2, 
                                  grepl('September|June', condition))$sample]
data_wo_na <- na.omit(data_only_dec_and_sep)
#data_wo_na <- na.omit(data_irs)
SignProteinInfo_0 <- edger_analysis_wrap(norm_data = data_wo_na, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("September", 
                                                          "June"), 
                                         metafile = meta2)

sign_protein_info_sub <- subset(SignProteinInfo_0, FDR < 0.05)

sign_protein_info_sub_allcomb <- rbind(sign_protein_info_sub_allcomb,
                                       sign_protein_info_sub)

protein_changed_unique <- unique(sign_protein_info_sub_allcomb$protein)
sign_protein_info_sub_allcomb_ <- unique(sign_protein_info_sub_allcomb[1:2]) 
  
write.table(data_wo_na, file.path(dir_to_results, 
                                  paste0('intensities_after_slNorm_woNA_', 
                                         tolower(species),'.csv')), 
            sep = '\t', quote = F, row.names = T)

write.table(SignProteinInfo_0, 
            file = file.path(paste0(dir_to_data, '/',
                   species, '_edgeResults_JUNvsSEP_proteinGroups.csv')),
            sep = '\t', quote = F, row.names = F)

save(protein_changed_unique, file = paste0(dir_to_results, 
                                           '/sign_changed_proteins_allmonths.rda'))

### Heatmap for all genes ######################################################

# Prepare intensities data: scale the data -> transport it back 
intensities <- data_wo_na
scaled_intensities <- apply(intensities, 1, scale)
scaled_intensities <- t(scaled_intensities)
colnames(scaled_intensities) <- meta$sample

scaled_intensities_only_sign <- 
  scaled_intensities[row.names(scaled_intensities) %in% protein_changed_unique,]

# Prepare the annotation table:
annot_col <- data.frame(Months = factor(meta$condition),
                        Sex = factor(meta$sex))
annot_col <- data.frame(Months = factor(meta2$condition))
rownames(annot_col) <- colnames(scaled_intensities)
# for BK
annot_colours <- list(Months = c("December_BK" = "dodgerblue1",
                                 "January_BK" = "cyan1",
                                 "June_BK" = "goldenrod1",
                                 "November_BK" = "brown",
                                 "pool" = "black",
                                 "September_BK" = "darkorange1"),
                      Sex = c("female" = "coral1",
                              "male" = "aquamarine4",
                              "NA" = "grey"))
annot_colours <- list(Months = c("December" = "dodgerblue1",
                                 "January" = "cyan1",
                                 "June" = "goldenrod1",
                                 "November" = "brown",
                                 "pool" = "black",
                                 "September" = "darkorange1"))
# for PB
annot_colours <- list(Months = c("December_PB" = "dodgerblue1",
                                 "November_PB" = "brown",
                                 "pool" = "black",
                                 "September_PB" = "darkorange1"),
                      Sex = c("female" = "coral1",
                              "male" = "aquamarine4",
                              "NA" = "grey"))
# Draw the heatmap
to_plot <- scaled_intensities
to_plot <- scaled_intensities_only_sign

### !!!! the name !!!!

pheatmap(to_plot,
         color = plasma(10),
         show_rownames = F,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         annotation_col = annot_col,
         annotation_colors = annot_colours,
         height = 10,
         width = 15,
         filename = file.path(dir_to_results, 
                              paste0(species, '_', 
                                     #population,
                                     '_heatmap_genes_sign', '.png')))
if (dev.cur() != 1) {dev.off()}

# clustering of seasonal profiles and visualizing temporal profiles of proteins -> time_series.R

#### Plot proteins (boxplots) for G. lacustris

proteins_to_plot <- subset(sign_protein_info_sub_allcomb_, grepl('hemocyanin', 
                                                geneSymbol, ignore.case = T))

proteins_to_plot <- subset(sign_protein_info_sub_allcomb_, !grepl('hemocyanin', 
                                                geneSymbol, ignore.case = T))

proteins_to_plot <- subset(SignProteinInfo_0, grepl('hemocyanin', 
                          geneSymbol, ignore.case = T)) # to look at all hemocyanins

scaled_intensities_sign_wo_pools <- 
  scaled_intensities_only_sign[,!grepl('pool', colnames(scaled_intensities_only_sign))]
scaled_intensities_sign_wo_pools <- 
  scaled_intensities[,!grepl('pool', colnames(scaled_intensities))] # to look at all hemocyanins

to_plot_1 <- 
  scaled_intensities_sign_wo_pools[match(proteins_to_plot$protein, 
                                         rownames(scaled_intensities_sign_wo_pools)),]
to_plot_1 <- as.data.frame(to_plot_1)
colnames(to_plot_1) <- colnames(scaled_intensities_sign_wo_pools)
to_plot_1$annotation <- proteins_to_plot$geneSymbol
to_plot_1$annotation <- sub('PREDICTED: |-like|, partial|isoform X1', '', to_plot_1$annotation)

to_plot_1_long <- to_plot_1 %>%
  rownames_to_column(var = 'proteins') %>%
  pivot_longer(cols = starts_with('F'), names_to = 'sample',
               values_to = 'intensities')

to_plot_1_long$months <- meta[match(to_plot_1_long$sample, meta$sample),]$condition
to_plot_1_long$time <-
  ifelse(to_plot_1_long$months == 'September', 1,
         ifelse(to_plot_1_long$months == 'November' , 2,
                ifelse(to_plot_1_long$months == 'December', 3,
                       ifelse(to_plot_1_long$months == 'January', 4, 5))))

get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

#to_plot_1_long$go_annotation <- pep_annot_go[match(to_plot_1_long$proteins, 
#                                                   pep_annot_go$protein_group),]$go_annotation
to_plot_1_long$annotation <- str_wrap(to_plot_1_long$annotation,
                                       width = 23)
to_plot_1_long$whole_label <- sprintf('%s|%s', to_plot_1_long$proteins,
                                      to_plot_1_long$annotation)
to_plot_1_long$whole_label2 <- 
  factor(to_plot_1_long$whole_label,
         levels=unique(
         to_plot_1_long$whole_label[order(to_plot_1_long$annotation)]))

months_mean <- aggregate(intensities ~ time, to_plot_1_long, mean)

#to_plot_1_long <- subset(to_plot_1_long, intensities < 6) # if there are samples with intensity value more than 6
to_plot_1_long$Sex <- meta[match(to_plot_1_long$sample, meta$sample),]$sex

ggplot(to_plot_1_long, aes(time, intensities, group = time)) +
  facet_wrap(~ whole_label2, ncol = 5,
             labeller = as_labeller(get_protein_label)
             , scales = 'free_y'
  ) +
  geom_boxplot(outlier.color = NA, 
               #aes(color = time == 4, fill = time == 4), alpha = .3
  ) +
  geom_jitter(
    #aes(color = time == 4), 
    aes(color = Sex),
    size = 0.3) +
  #geom_line(data = months_mean, aes(time, intensities, group = 1), 
  #          color = "firebrick1", alpha = 0.80, 
  #          size = 1.2) +
  #  scale_color_manual(values = c('black', 'dodgerblue4')) +
  scale_color_manual(values = c('deeppink1', 'dodgerblue1', 'grey35') # sex 
                                       ) +
  #  scale_fill_manual(values = c(NA, 'dodgerblue4')) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  guides(fill = F 
         #, color = F
         ) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw()

ggsave(file.path(dir_to_results, paste0('boxplots_', 'all_Hemocyanins_withSex','.png')), 
       scale = 0.7,
       width = 14, height = 7)

scaled_intensities_table <- as.data.frame(scaled_intensities)
hist(scaled_intensities_table$`F/104/L/3C/3_1`)

### Plot chosen proteins:

proteins_to_plot <- 
  subset(SignProteinInfo_0, 
         grepl('GYG1|UGP|GBE|GYS2|sgg|agl1|glycogen|Hex-A|glycogenin|PGM2|Pgm', 
               geneSymbol, ignore.case = T))
proteins_to_plot <- 
  subset(SignProteinInfo_0, 
         grepl('GYG1|UGP|GBE|AGL', 
               geneSymbol, ignore.case = T)) 
scaled_intensities_sign_wo_pools <- 
  scaled_intensities[,!grepl('pool', colnames(scaled_intensities))] # to look at all hemocyanins

to_plot_1 <- 
  scaled_intensities_sign_wo_pools[match(proteins_to_plot$protein, 
                                         rownames(scaled_intensities_sign_wo_pools)),]
to_plot_1 <- as.data.frame(to_plot_1)
colnames(to_plot_1) <- colnames(scaled_intensities_sign_wo_pools)
to_plot_1$annotation <- proteins_to_plot$geneSymbol
to_plot_1$annotation <- sub('PREDICTED: |-like|, partial|isoform X1', '', to_plot_1$annotation)

to_plot_1_long <- to_plot_1 %>%
  rownames_to_column(var = 'proteins') %>%
  pivot_longer(cols = starts_with('F'), names_to = 'sample',
               values_to = 'intensities')

to_plot_1_long$months <- meta2[match(to_plot_1_long$sample, meta2$sample),]$condition
to_plot_1_long$time <-
  ifelse(to_plot_1_long$months == 'September', 1,
         ifelse(to_plot_1_long$months == 'November' , 2,
                ifelse(to_plot_1_long$months == 'December', 3,
                       ifelse(to_plot_1_long$months == 'January', 4, 5))))

get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

#to_plot_1_long$go_annotation <- pep_annot_go[match(to_plot_1_long$proteins, 
#                                                   pep_annot_go$protein_group),]$go_annotation
to_plot_1_long$annotation <- str_wrap(to_plot_1_long$annotation,
                                      width = 23)
to_plot_1_long$whole_label <- sprintf('%s|%s', to_plot_1_long$proteins,
                                      to_plot_1_long$annotation)
to_plot_1_long$whole_label2 <- 
  factor(to_plot_1_long$whole_label,
         levels=unique(
           to_plot_1_long$whole_label[order(to_plot_1_long$annotation)]))

months_mean <- aggregate(intensities ~ time, to_plot_1_long, mean)

#to_plot_1_long <- subset(to_plot_1_long, intensities < 6) # if there are samples with intensity value more than 6
to_plot_1_long$Sex <- meta[match(to_plot_1_long$sample, meta$sample),]$sex

ggplot(to_plot_1_long, aes(time, intensities, group = time)) +
  facet_wrap(~ whole_label2, ncol = 5,
             labeller = as_labeller(get_protein_label)
             , scales = 'free_y'
  ) +
  geom_boxplot(outlier.color = NA, 
               #aes(color = time == 4, fill = time == 4), alpha = .3
  ) +
  geom_jitter(
    #aes(color = time == 4), 
    #aes(color = Sex),
    size = 0.3) +
  geom_line(data = months_mean, aes(time, intensities, group = 1), 
            color = "firebrick1", alpha = 0.80, 
            size = 1.2) +
  scale_color_manual(values = c('black', 'dodgerblue4')) +
  #scale_color_manual(values = c('deeppink1', 'dodgerblue1', 'grey35') # sex 
  #) +
  #  scale_fill_manual(values = c(NA, 'dodgerblue4')) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  guides(fill = F 
         #, color = F
  ) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw()

ggsave(file.path(dir_to_results, paste0('boxplots_', 'all_Hemocyanins_withSex','.png')), 
       scale = 0.7,
       width = 14, height = 7)



### tmp

data_only_dec_and_sep_long <- data_only_dec_and_sep %>% 
  rownames_to_column(var = 'pg') %>%
  pivot_longer(cols = starts_with('F'))

data_only_dec_and_sep_long_sub<- subset(data_only_dec_and_sep_long, 
                                        pg == 'TRINITY_DN10891_c0_g3_i6.p1_Eve')

data_only_dec_and_sep_long_sub$temp <- 
  meta2[match(data_only_dec_and_sep_long_sub$name, meta2$sample),]$condition

ggplot(data_only_dec_and_sep_long_sub, aes(temp, value)) +
  geom_boxplot() +
  geom_jitter() +
  theme_light()



