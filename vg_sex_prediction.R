library(metR)
library(e1071)
library(tidyr)
library(tibble)
library(ggplot2)
### upload intensities
species <- 'Eve'
dir_to_intensities <- paste0('labeglo2/MS_results/Field/', species, '/', species, 
                             '_refchannels_all')
intensities <- read.csv(file.path(dir_to_intensities, 
                                     paste0('intensities_after_slNorm_withNA', 
                                            tolower(species),'.csv')),
                           header = T, sep = '\t') 

intensities <- intensities[, !grepl('pool', colnames(intensities),
                                          ignore.case = T)]

rownames(intensities) <- intensities[,1]
intensities <- intensities[-1]

pep_annot <- read.csv(file.path(dir_to_intensities, 
                                paste0('annot_proteinGroups_field_', tolower(species), 
                                       '_withGOannotation.csv')), 
                      sep = '\t', header = T)

### upload metafile
### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}

path2meta <- paste0('labeglo2/MS_results/Field/', species,
                    '/', species, '_refchannels_all', '/metadata_predictedSex.csv') # with predicted sex

meta <- meta_upload(path2meta, species)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
meta2$condition <- sub('_BK', '', meta2$condition)
meta2 <- subset(meta2, !grepl('pool', sample))
meta2$sex2[is.na(meta2$sex2)] <- 'NA'
meta2$condition <- factor(meta2$condition, levels = unique(meta2$condition),
                          labels = unique(meta2$condition))

### scale intensities:
scaled_intensities <- apply(intensities, 1, scale)
scaled_intensities <- t(scaled_intensities)
colnames(scaled_intensities) <- meta2$sample

scaled_intensities <- 
  scaled_intensities[,!grepl('pool', colnames(scaled_intensities))]
scaled_intensities <- as.data.frame(scaled_intensities)
scaled_intensities_long <- scaled_intensities %>%
  rownames_to_column(var = 'proteins') %>%
  pivot_longer(cols = starts_with('F'), 
               names_to = 'sample', 
               values_to = 'intensities')

scaled_intensities_long$annotation <- 
  pep_annot[match(scaled_intensities_long$proteins, 
                     pep_annot$protein_group),]$upd_full_annot

### Take only Vg
vg_only <- subset(scaled_intensities_long, annotation == 'Vg')

vg_only$Month <- meta2[match(vg_only$sample, meta2$sample),]$condition
#vg_only$Month <- sub('_BK', '', vg_only$Month)
vg_only$time <-
  ifelse(vg_only$Month == 'September', 1,
         ifelse(vg_only$Month == 'November' , 2,
                ifelse(vg_only$Month == 'December', 3,
                       ifelse(vg_only$Month == 'January', 4, 5))))
get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

vg_only$whole_label <- sprintf('%s|%s', vg_only$proteins,
                                        vg_only$annotation)
vg_only$whole_label <- sprintf('%s|%s', vg_only$proteins,
                               vg_only$annotation)
vg_only$whole_label2 <- 
  factor(vg_only$whole_label,
         levels=unique(
           vg_only$whole_label[order(vg_only$annotation)]))

vg_only$Sex <- meta[match(vg_only$sample, meta$sample),]$sex

ggplot(vg_only, aes(time, intensities, group = time)) +
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

unique(vg_only$proteins)[c(4, 7, 9)]

vg_to_plot <- subset(vg_only, 
                     proteins == unique(vg_only$proteins)[4] | 
                       proteins == unique(vg_only$proteins)[7] | 
                       proteins == unique(vg_only$proteins)[9])

ggplot(vg_to_plot, aes(time, intensities, group = time)) +
  facet_wrap(~ whole_label2, ncol = 4
             #,
             #labeller = as_labeller(get_protein_label)
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

#vg3 <- subset(vg_to_plot, proteins)

### pca using 3 VGs
chosen_month <- c('December', 'January', 'November', 'September')

data_sl_t_vg <- data_sl_t[,colnames(data_sl_t) == unique(vg_only$proteins)[3] | 
                            colnames(data_sl_t) == unique(vg_only$proteins)[5] | 
                            colnames(data_sl_t) == unique(vg_only$proteins)[6]]
data_sl_t_vg <- data_sl_t[,colnames(data_sl_t) == unique(vg_only$proteins)[4] | 
                            colnames(data_sl_t) == unique(vg_only$proteins)[7] | 
                            colnames(data_sl_t) == unique(vg_only$proteins)[9]]
data_sl_t_vg <- data_sl_t_vg[!grepl('pool', rownames(data_sl_t_vg)),]
data_sl_t_vg <- data_sl_t_vg[match(meta2[meta2$condition %in% chosen_month,]$sample,
                                   rownames(data_sl_t_vg)),]

pca_res <- prcomp(data_sl_t_vg)
meta_vg <- subset(meta2, condition != 'pool')
meta_vg <- subset(meta2, condition %in% chosen_month)
meta_vg$sex[is.na(meta_vg$sex)] <- 'NA'
autoplot(pca_res, data=meta_vg, 
         colour='condition', 
         alpha = 'sex',
         size = 2.5, 
         color = 0,
         shape='sex'
         ) +
  scale_color_manual('Months', 
                     values = palette.colors(n=6, 'Dark2')) +
  scale_shape_manual('Sex', values = c(16, 17, 4)) +
  scale_alpha_manual('Sex', values = c(1, 1, 0.4)) + 
  theme_light()

ggsave(filename = file.path(dir_to_results, paste0('pca_', 
       species,'Vg_SepNovDecJan.png')), 
       scale = 1,
       width = 6, height = 3)

res <- cbind(meta_vg, as.data.frame(pca_res$x))
tmp <- data_sl_t_vg
colnames(tmp) <- c('vg1', 'vg2', 'vg3')
res <- cbind(meta_vg, tmp)
ggplot(res) +
  geom_histogram(aes(PC1, fill=sex), binwidth=0.02) +
  theme_bw()

#### svm

res$sex_fact <- with(res, factor(ifelse(sex == 'NA', NA, sex)))
res_part <- res[c('sex_fact', 'PC1', 'PC2')]
res_part <- res[c('sex_fact', 'vg1', 'vg2')]
svm_model <- svm(sex_fact ~ PC1, data = res_part, cost = 1, probability = T)
svm_model <- svm(sex_fact ~ vg1 + vg2, data = res_part, cost = 1, probability = T)

grange = apply(dplyr::select(res_part, starts_with('PC')), 2,
               function(x) { r <- range(x); c(r, r[2] - r[1]) }) %>%
  as.data.frame

grange = apply(dplyr::select(res_part, starts_with('vg')), 2,
               function(x) { r <- range(x); c(r, r[2] - r[1]) }) %>%
  as.data.frame

grid <- expand.grid(
  PC1 = seq(grange$PC1[1] - grange$PC1[3] * 0.05,
            grange$PC1[2] + grange$PC1[3] * 0.05, length.out = 1000),
  PC2 = seq(grange$PC2[1] - grange$PC2[3] * 0.05,
            grange$PC2[2] + grange$PC2[3] * 0.05, length.out = 100))

cor(res$vg2, res$vg3, method='pearson')
summary(pca_res)


grid <- expand.grid(
  vg1 = seq(grange$vg1[1] - grange$vg1[3] * 0.05,
            grange$vg1[2] + grange$vg1[3] * 0.05, length.out = 100),
  vg2 = seq(grange$vg2[1] - grange$vg2[3] * 0.05,
            grange$vg2[2] + grange$vg2[3] * 0.05, length.out = 100)
)

grid$sex_fact <- predict(svm_model, grid)
grid$sex_prob <- (predict(svm_model, grid, probability = T) %>%
                    attr('probabilities'))[,1]

# ggplot(grid, aes(PC1, PC2)) +
#   geom_tile(aes(fill = sex_prob)) +
#   #geom_point(data = res, size=2) +
#   theme_bw()

ggplot(grid, aes(vg1, vg2)) +
  geom_tile(aes(fill = sex_prob)) +
  geom_contour(aes(z = sex_prob), color='gray30') +
  geom_text_contour(aes(z = sex_prob), stroke = 0.1,
                    label.placer = label_placer_fraction(frac = 0.3)) +
  geom_point(aes(shape = sex, alpha = sex, color = sex), data = res, size = 6) +
  scale_fill_gradientn('Probability\nto be male', colours = c('pink', 'white', 'dodgerblue'),
                       values = c(0, 0.5, 1)) +
  scale_shape_manual('Sex', values = c("\u2640", "\u2642", '*'), 
                     labels = c('Female', 'Male', 'Unknown')) +
  scale_alpha_manual('Sex', values = c(1, 1, 0.7), 
                     labels = c('Female', 'Male', 'Unknown')) +
  scale_color_manual('Sex', values = c('black', 'black', 'coral3'), 
                     labels = c('Female', 'Male', 'Unknown')) +
  #scale_fill_manual(values = colorRampPalette(c("pink", "white", "blue"))(10)) +
  #coord_cartesian(xlim = c(-0.185, 0.35)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw()

ggsave(filename = file.path(dir_to_results, paste0('pca_', 
       species,'_afterSvm_to_predict_sex.png')), 
       scale = 0.65,
       width = 12, height = 6)

ggplot(grid, aes(PC1, PC2)) +
  geom_contour_filled(aes(z = sex_prob), binwidth = 0.1) +
  geom_point(aes(shape = sex2), data = res, size = 2) +
  #scale_fill_gradientn(colours = c('red', 'white', 'blue')) +
  scale_fill_manual(values = colorRampPalette(c("pink", "white", "blue"))(10)) +
  theme_bw()

####

meta_dec <- subset(meta, grepl('December', condition))
meta_dec$sex_prediction <- ifelse(res[meta_dec$sample,]$PC1 < 0.1, 
                                  'female', 'male')
dec_sexes <- meta_dec$sex_prediction
names(dec_sexes) <- meta_dec$sample

meta_tmp <- left_join(meta, meta_dec, by=colnames(meta))

meta_tmp$sex2 <- 
  ifelse(meta_tmp$sex == 'female', 'female', 
         ifelse(meta_tmp$sex == 'male', 'male',
                ifelse(meta_tmp$sex_prediction == 'female', 'female',
                       ifelse(meta_tmp$sex_prediction == 'male', 'male', 'NA'))))

#### if you decided to assign sex to January samples
# I assigned sex to the January samples based on the Vg boxplots
meta_tmp[meta_tmp$sample == 'F/100/BK/4/1p',]$sex2 <- 'male'
meta_tmp[meta_tmp$sample == 'F/100/BK/4/2p',]$sex2 <- 'female'
meta_tmp[meta_tmp$sample == 'F/100/BK/4/3p',]$sex2 <- 'female'
meta_tmp[meta_tmp$sample == 'F/100/BK/4/4p',]$sex2 <- 'female'
meta_tmp[meta_tmp$sample == 'F/100/BK/4/5p',]$sex2 <- 'female'
####
data_sl_t <- na.omit(data_sl) %>% t %>% as.data.frame
data_sl_t <- data_sl_t[!grepl('pool', rownames(data_sl_t)),] # to draw wo pools
pca_res <- prcomp(data_sl_t)
meta2 <- meta_tmp
meta2$experiment <- as.factor(meta2$experiment)
meta2$condition <- sub('_BK', '', meta2$condition)
meta2 <- subset(meta2, condition != 'pool')
meta2$sex2[is.na(meta2$sex2)] <- 'NA'
#meta2 <- subset(meta2, sex2 != 'NA') # only for the figure with samples with known sex
#data_sl_t <- data_sl_t[meta2$sample,] # only for the figure with samples with known sex 

autoplot(pca_res, data=meta2, colour='condition', size = 3, 
         shape='sex2', 
         frame = F, 
         x = 1,
         y = 2) + 
  theme_light() + 
  scale_color_manual('Condition', 
                     values = palette.colors(n=6, 'Dark2')) +
  scale_fill_manual('Condition', 
                    values = palette.colors(n=6, 'Dark2')) +
  scale_shape_manual('Sex', values = c(16, 17, 4))

ggsave(filename = file.path(dir_to_results, paste0('pca_', 
            species,'_seasonalExp_withPredictedSexDecJan_ONLYsex.png')), 
       scale = 1.4,
       width = 6, height = 3)


write.table(meta_tmp, 
          '~/labeglo2/MS_results/Field/Eve/Eve_refchannels_all/metadata_predictedSex.csv', 
          sep = '\t', row.names = F)


