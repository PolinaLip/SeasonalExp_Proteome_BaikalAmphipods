### To plot PCA for both Field and Gradual Temp Decrease experiments ####

library(stringr)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(ggfortify) 
library(rstatix)

dir_to_save_results <- 'labeglo2/MS_results/Field_Grad_Comparison/'

species <- 'Eve'
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}

path2meta <- paste0('labeglo2/MS_results/Field/', species,
                    '/', species, '_refchannels_all', '/metadata_predictedSex.csv')

meta <- meta_upload(path2meta, species)
meta <- meta[!grepl('PB|pool_12', meta$sample),]
meta <- subset(meta, condition != 'pool')
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')
meta_field <- meta[-6]
meta_field <- meta[-c(6,7)]
meta_field$experiment <- as.numeric(meta_field$experiment)
meta_field <- subset(meta_field, condition != 'pool')

path2meta_grad <-
  paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/',
         species, '/Metadata_Proteus_withPredSex.csv')

meta_g <- meta_upload(path2meta_grad, species)
meta_g <- subset(meta_g, !grepl('PB', sample))
meta_g <- subset(meta_g, !grepl('Pool', sample))

meta_grad <- meta_g

dir_to_field <- paste0('labeglo2/MS_results/Field/', species, '/', species, 
                       '_refchannels_all')
data_field <- read.table(file.path(dir_to_field, 
                                 paste0('intensities_after_slNorm_withNA', 
                                        tolower(species), '.csv')), header = T)

rownames(data_field) <- data_field$protein_group
data_field <- data_field[-1]
data_field <- data_field[,!grepl('pool', colnames(data_field))]
colnames(data_field) <- meta_field$sample

dir_to_grad <-paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/', species)
data_grad <- read.table(file.path(dir_to_grad, 
                                 paste0('intensities_after_slNorm_', 
                                        tolower(species), '.csv')), header = T)
rownames(data_grad) <- data_grad$protein_group
data_grad <- data_grad[-1]

#### upload annotations
prot_annot_field <- read.csv(file.path(dir_to_field, 
                                paste0('annot_proteinGroups_field_', tolower(species), 
                                       '_withGOannotation.csv')), 
                      sep = '\t', header = T)

prot_annot_grad <- 
  read.csv(file.path(dir_to_grad, paste0('annot_proteinGroups_', 
                                       tolower(species), '_withGOannotation.csv')), 
           sep = '\t', header = T) 

#### scale intensities
scaled_intensities_field <- apply(data_field, 1, scale)
scaled_intensities_field <- t(scaled_intensities_field)
colnames(scaled_intensities_field) <- meta_field$sample
scaled_intensities_field <- as.data.frame(scaled_intensities_field)

scaled_intensities_grad <- apply(data_grad, 1, scale)
scaled_intensities_grad <- t(scaled_intensities_grad)
colnames(scaled_intensities_grad) <- meta_grad$sample
scaled_intensities_grad <- as.data.frame(scaled_intensities_grad)

#### 

new_protein_name_field <- 
  lapply(str_split(rownames(scaled_intensities_field), ';'), sort)
scaled_intensities_field$sorted_protein <- 
  sapply(new_protein_name_field, paste0, collapse=';')

new_protein_name_gradual <- 
  lapply(str_split(rownames(scaled_intensities_grad), ';'), sort)
scaled_intensities_grad$sorted_protein <- 
  sapply(new_protein_name_gradual, paste0, collapse=';')

####

protein_name_grid_for_pca <- expand.grid(x=new_protein_name_field, 
                                         y=new_protein_name_gradual)

####
#protein_name_grid_for_pca$match_rate <- apply(protein_name_grid_for_pca, 1,
#  function(row) length(intersect(row$x, row$y)) / min(length(row$x), length(row$y)))
protein_name_grid_for_pca$match_rate <- apply(protein_name_grid_for_pca, 1,
  function(row) length(intersect(row$x, row$y)) / max(length(row$x), length(row$y)))

sum(protein_name_grid_for_pca$match_rate > 0)

###
#protein_name_grid_sub <- subset(protein_name_grid_for_pca, match_rate > 0)
protein_name_grid_sub <- subset(protein_name_grid_for_pca, match_rate == 1) # if take only complete matches in protein group names
protein_name_grid_sub$field_pg <- sapply(protein_name_grid_sub$x, paste0, collapse=';')
protein_name_grid_sub$grad_pg <- sapply(protein_name_grid_sub$y, paste0, collapse=';')

###
prot_single_entry_field <- table(protein_name_grid_sub$field_pg)
prot_single_entry_grad <- table(protein_name_grid_sub$grad_pg)

prot_single_entry_field <- prot_single_entry_field[prot_single_entry_field == 1]
prot_single_entry_grad <- prot_single_entry_grad[prot_single_entry_grad == 1]

protein_name_grid_sub_filtered <- 
  subset(protein_name_grid_sub, field_pg %in% names(prot_single_entry_field))
protein_name_grid_sub_filtered <- 
  subset(protein_name_grid_sub_filtered, grad_pg %in% names(prot_single_entry_grad))

### take only ovellaped proteins from two experiment 

data_field_sub <- 
  subset(scaled_intensities_field, sorted_protein %in% protein_name_grid_sub_filtered$field_pg)
data_grad_sub <- 
  subset(scaled_intensities_grad, sorted_protein %in% protein_name_grid_sub_filtered$grad_pg)

data_field_sub_withAnnot <- data_field_sub
data_grad_sub_withAnnot <- data_grad_sub
data_field_sub_withAnnot$annotation <-
  prot_annot_field[match(rownames(data_field_sub_withAnnot), 
                         prot_annot_field$protein_group),]$upd_full_annot
data_grad_sub_withAnnot$annotation <-
  prot_annot_grad[match(rownames(data_grad_sub_withAnnot), 
                         prot_annot_grad$protein_group),]$upd_full_annot

data_joined <- inner_join(data_field_sub, data_grad_sub, 
                          by = c("sorted_protein" = "sorted_protein"))

rownames(data_joined) <- data_joined$sorted_protein
data_joined <- data_joined[, colnames(data_joined) != 'sorted_protein']
meta_field$sex <- meta_field$sex2 
meta_field <- meta_field[-6]
meta_joined <- rbind(meta_field, meta_grad)

### PCA
rownames(data_joined) <- data_joined$sorted_protein
data_joined <- data_joined[, colnames(data_joined) != 'sorted_protein']
data_joined_t <- na.omit(data_joined) %>% t %>% as.data.frame
#data_sl_t <- data_sl_t[!grepl('pool', rownames(data_sl_t)),] # to draw wo pools

pca_res <- prcomp(data_joined_t)
#meta2 <- meta
#meta2$experiment <- as.factor(meta2$experiment)
meta_joined$condition <- sub('_BK', '', meta_joined$condition)
meta_joined$sex[is.na(meta_joined$sex)] <- 'Unknown'
#meta2 <- subset(meta2, sex != 'NA') # only for the figure with samples with known sex
#data_sl_t <- data_sl_t[meta2$sample,] # only for the figure with samples with known sex 
#meta2 <- subset(meta2, !grepl('pool', sample)) # to draw wo pools
autoplot(pca_res, data=meta_joined, colour='condition', frame = T) + theme_light()
autoplot(pca_res, data=meta_joined, colour='condition', size = 3,
         shape = 'sex',
         frame = F, 
         x = 1,
         y = 2) + 
  theme_light() + 
  scale_color_manual('Condition', 
                     values = palette.colors(n=8, 'Dark2'))

# autoplot(pca_res, data=meta2, colour='condition', label = TRUE, label.size = 5) + theme_light()
# meta2[rownames(meta2) == 44,]$sample
# ! outlier found: F/100/BK/3c/2m_2
dir_to_results <- 'labeglo2/MS_results/Field_Grad_Comparison/'
ggsave(filename = file.path(dir_to_results, paste0('pca_', 
                  species,'_seasonalExp_GradDecrease.png')), 
       scale = 1,
       width = 6, height = 3)

### To plot correlations between different months and conditions from gradual decrease temperature

data_joined_long <- pivot_longer(data_joined %>% rownames_to_column('PG'), 
                            starts_with(c('F','E', 'G')), names_to = 'Samples',
                            values_to = 'intensities')

data_joined_long$condition <- 
  meta_joined[match(data_joined_long$Samples, meta_joined$sample),]$condition
data_joined_long$condition <- sub('_BK|EveBK_|Ecy_|Gla_', '', data_joined_long$condition)

data_joined_mean <- aggregate(intensities ~ PG + condition, 
                              data = data_joined_long, median, na.rm = T)

data_joined_mean_wide <- pivot_wider(data_joined_mean, names_from = condition,
                                     values_from = intensities)
colnames(data_joined_mean_wide)

cond1 <- "10.5"
cond2 <- "June"

cor_test_res <- 
  cor.test(c(data_joined_mean_wide[,cond1][[1]]), 
           c(data_joined_mean_wide[,cond2][[1]]))
  
ggplot(data_joined_mean_wide, aes(`10.5`, `June`)) +
  geom_point(alpha=.7, color='gray70') +
  geom_smooth(method='lm') +
  annotate(geom='text', size = 3,
           x = 1.1, y = 2.1, hjust = 0,
           label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                          'p-value = ', p_format(cor_test_res$p.value,
                                               accuracy = 0.005
                                               #accuracy = cor_test_res$p.value * 10
                                               ))) +
  scale_x_continuous(limits = c(-2.3, 2)) +
  scale_y_continuous(limits = c(-1.5, 2.3)) +
  theme_light()

ggsave(file.path(dir_to_save_results, 
                 paste0(species, '_corPlot_', cond1, '_', cond2, '', '.png')),
       width = 6, height = 4, scale = .7)


