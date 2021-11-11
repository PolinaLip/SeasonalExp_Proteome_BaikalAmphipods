### To predict sex in June samples
## Conclusion: it seems that it is not possible to do for the samples from June 
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
#meta <- subset(meta, sex2 != 'NA')

### pca
data_sl_t <- na.omit(data_sl_) %>% t %>% as.data.frame
#data_sl_t <- data_sl_t[!grepl('pool', rownames(data_sl_t)),] # to draw wo pools

meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
meta2$condition <- sub('_BK', '', meta2$condition)

pca_res <- prcomp(data_sl_t, rank. = 15)
pca_x <- pca_res$x
rownames(pca_x) <- with(meta2, sprintf('%s_%s', condition, sex2)) %>% make.unique()
tree <- hclust(dist(pca_x))
plot(tree)

library(tidyverse)

meta_pca <- cbind(meta2, pca_x)# %>% select(sample, sex)
meta_train <- filter(meta_pca, !is.na(sex2)) %>%
  mutate(sex3 = ifelse(sex2 == 'male', 1, 0))
meta_test <- filter(meta_pca, is.na(sex2))
model <- glm(sex3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7,
             data = meta_train,
             family = 'binomial')
model$coef
with(summary(model), 1 - deviance/null.deviance)
meta_test$sex3 <- predict(model, meta_test, type = 'response')
meta_test

######
meta_pca <- cbind(meta2, pca_x) %>% mutate(sex = as.factor(sex))
meta_train <- filter(meta_pca, !is.na(sex))
meta_test <- filter(meta_pca, is.na(sex))
model <- glm(sex ~ PC1,# + PC2 + PC3,# + PC4 + PC5,# + PC6 + PC7 + PC8 +
               #PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15,
             data = meta_train,
             family = 'binomial')
model
with(summary(model), 1 - deviance/null.deviance)
meta_test$pred_sex <- levels(meta_pca$sex)[
  1 + round(predict(model, meta_test, type = 'response'))]
meta_test[c('sample', 'sex2', 'pred_sex')]
with(meta_test, sum(!is.na(sex2) & as.character(sex2) == pred_sex) / sum(!is.na(sex2)))

#######

meta_pca <- cbind(meta2, pca_x) %>% mutate(sex2 = as.numeric(sex2 == 'male'))
test_months <- c('November', 'January')
test_months <- 'June'
meta_train <- filter(meta_pca, !is.na(sex2) & !(condition %in% test_months))
meta_test <- filter(meta_pca, condition %in% test_months)
model <- glm(sex2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
             PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15,
             data = meta_train,
             family = 'binomial')
model
with(summary(model), 1 - deviance/null.deviance)
meta_test$pred_sex <- round(predict(model, meta_test, type = 'response'))
meta_test[c('sample', 'sex2', 'pred_sex')]
with(meta_test, sum(!is.na(sex2) & as.character(sex2) == pred_sex) / sum(!is.na(sex2)))

# svm approach
library(e1071)
svm_model <- svm(sex2 ~ PC1 + PC2 + PC3, # + PC4 + PC5 + PC6 + PC7 + PC8, # +
                   #PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15,
                 data = meta_train)

meta_train$pred_sex <- round(predict(svm_model, meta_train[, 9:23], probability = TRUE))
with(meta_train, sum(!is.na(sex2) & as.character(sex2) == pred_sex) / sum(!is.na(sex2)))

meta_test$pred_sex <- round(predict(svm_model, meta_test[, 9:23], probability = TRUE))
meta_test[c('sample', 'sex2', 'pred_sex')]
with(meta_test, sum(!is.na(sex2) & as.character(sex2) == pred_sex) / sum(!is.na(sex2)))
