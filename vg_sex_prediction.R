scaled_intensities <- 
  scaled_intensities[,!grepl('pool', colnames(scaled_intensities))]
scaled_intensities <- as.data.frame(scaled_intensities)
scaled_intensities_long <- scaled_intensities %>%
  rownames_to_column(var = 'proteins') %>%
  pivot_longer(cols = starts_with('F'), 
               names_to = 'sample', 
               values_to = 'intensities')

scaled_intensities_long$annotation <- 
  pep_annot_go[match(scaled_intensities_long$proteins, 
                     pep_annot_go$protein_group),]$annotation

vg_only <- subset(scaled_intensities_long, annotation == 'Vg')

vg_only$Month <- meta[match(vg_only$sample, meta$sample),]$condition
vg_only$Month <- sub('_BK', '', vg_only$Month)
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
vg_only$whole_label2 <- 
  factor(vg_only$whole_label,
         levels=unique(
           vg_only$whole_label[order(vg_only$annotation)]))

vg_only$Sex <- meta[match(vg_only$sample, meta$sample),]$sex

ggplot(vg_only, aes(time, intensities, group = time)) +
  facet_wrap(~ whole_label2, ncol = 4,
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

unique(vg_only$proteins)[c(3, 5, 6)]

vg_to_plot <- subset(vg_only, 
                     proteins == unique(vg_only$proteins)[3] | 
                       proteins == unique(vg_only$proteins)[5] | 
                       proteins == unique(vg_only$proteins)[6])

ggplot(vg_to_plot, aes(time, intensities, group = time)) +
  facet_wrap(~ whole_label2, ncol = 4,
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

#vg3 <- subset(vg_to_plot, proteins)

### pca using 3 VGs
choosen_month <- 'December'
data_sl_t_vg <- data_sl_t[,colnames(data_sl_t) == unique(vg_only$proteins)[3] | 
                            colnames(data_sl_t) == unique(vg_only$proteins)[5] | 
                            colnames(data_sl_t) == unique(vg_only$proteins)[6]]
data_sl_t_vg <- data_sl_t_vg[!grepl('pool', rownames(data_sl_t_vg)),]
data_sl_t_vg <- data_sl_t_vg[match(meta[grepl(choosen_month, meta$condition),]$sample,
                                   rownames(data_sl_t_vg)),]

pca_res <- prcomp(data_sl_t_vg)
meta_vg <- subset(meta2, condition != 'pool')
meta_vg <- subset(meta2, grepl(choosen_month, condition))
autoplot(pca_res, data=meta_vg, colour='condition', size = 3, shape='sex') + theme_light()

res <- pca_res$x
#res <- as_tibble(res)
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


