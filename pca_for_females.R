# to plot PCA for females in December and January and look if they separate in groups (spoiler - they not, we have some outliers but not groups)
library(dplyr)
library(ggfortify) 

species <- 'Eve'
dir_annot <- paste0('labeglo2/MS_results/Field/', species, 
                    '/', species,'_refchannels_all/')
dir_metafile <- paste0('labeglo2/MS_results/Field/', species, 
                       '/', species,'_refchannels_all/')
dir_to_intensities <- paste0('labeglo2/MS_results/Field/', species, 
                             '/', species,'_refchannels_all/')
current_dir <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species,'_refchannels_all/withALLzeros_inAtLeastOneCond/')
pg_annot_complete <- read.csv(file.path(dir_annot, 
                                        paste0('annot_proteinGroups_field_', tolower(species), 
                                               '_withGOannotation_with_organisms.csv')), 
                              sep = '\t', header = T) # eggnog and diamond annotation
metafile <- 'metadata_predictedSex.csv'

###
data_sl <- read.csv(file.path(dir_to_intensities, 
                              paste0('intensities_after_slNorm_withNA', 
                                     tolower(species),'.csv')),
                    header = T, sep = '\t') 
data_sl_ <- data_sl[, !grepl('pool', colnames(data_sl),
                                          ignore.case = T)]

rownames(data_sl_) <- data_sl_[,1]
data_sl_ <- data_sl_[-1]

### Upload metadata
meta_complete <- read.delim(file.path(dir_metafile, metafile), header=TRUE, sep="\t")
meta <- meta_complete
meta <- subset(meta, condition != 'pool')
colnames(data_sl_) <- meta$sample
meta <- subset(meta, sex2 != 'NA')

### choose females and only December:
meta <- subset(meta, sex2 == 'female' & grepl('December|January', condition))
meta <- subset(meta, sex2 == 'female')
###
data_sl_chosen <- 
  data_sl_[,colnames(data_sl_) %in% meta$sample]

data_sl_chosen_cleaned <- data_sl_chosen[rowSums(is.na(data_sl_chosen)) <= 2,]

data_sl_t <- na.omit(data_sl_chosen_cleaned) %>% t %>% as.data.frame
#data_sl_t <- data_sl_t[!grepl('pool', rownames(data_sl_t)),] # to draw wo pools

pca_res <- prcomp(data_sl_t)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
meta2$condition <- sub('_BK', '', meta2$condition)
#meta2 <- subset(meta2, sex != 'NA') # only for the figure with samples with known sex
#data_sl_t <- data_sl_t[meta2$sample,] # only for the figure with samples with known sex 
#meta2 <- subset(meta2, !grepl('pool', sample)) # to draw wo pools
autoplot(pca_res, data=meta2, x = 2, y = 3, colour = 'condition') + theme_light()
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

### Perform clustering for females to look for the groups
sampleTree <- hclust(dist(...), method = "average")




