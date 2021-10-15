#### Check for parasites ####

## upload annotation
dir_annot <- paste0('labeglo2/MS_results/Field/', species, 
                    '/', species,'_refchannels_all/')
pg_annot_complete <- read.csv(file.path(dir_annot, 
                              paste0('annot_proteinGroups_field_', tolower(species), 
                                               '_withGOannotation_with_organisms.csv')), 
                              sep = '\t', header = T) 

## upload and work on metafile
dir_metafile <- paste0('labeglo2/MS_results/Field/', species, 
                       '/', species,'_refchannels_all/')
metafile <- 'Metadata_Proteus.tsv'
meta_complete <- read.delim(file.path(dir_metafile, metafile), header=TRUE, sep="\t")
meta <- meta_complete
population <- 'BK'
meta <- subset(meta, condition != 'pool')
meta[which(meta$sex == ''),]$sex <- NA
meta <- meta[grepl(population, meta$sample),]
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')


## upload and clean intensities 
dir_to_intensities <- paste0('labeglo2/MS_results/Field/', species, 
                             '/', species,'_refchannels_all/')
current_dir <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species,'_refchannels_all/withALLzeros_inAtLeastOneCond/')

data_with_NAs <- read.csv(file.path(dir_to_intensities, 
                                     paste0('intensities_after_slNorm_withNA', 
                                            tolower(species),'.csv')),
                           header = T, sep = '\t') 
data_with_NAs <- data_with_NAs[, !grepl('pool', colnames(data_with_NAs),
                                          ignore.case = T)]
rownames(data_with_NAs) <- data_with_NAs[,1]
data_with_NAs <- data_with_NAs[-1]
colnames(data_with_NAs) <- meta$sample

##### 

data_parasites <- subset(pg_annot_complete, grepl('Annelida|Nemato|Tremato|Cesto', organism_diamond))

intensity_scaled <- apply(data_with_NAs, 1, scale)
intensity_scaled <- t(intensity_scaled)
intensity_scaled <- data.frame(intensity_scaled)
colnames(intensity_scaled) <- rownames(metadata)

intensity_long <- intensity_scaled %>% 
  rownames_to_column(var = 'protein') %>%
  pivot_longer(-protein, values_to = 'intensity', names_to = 'sample')

intensity_long$annotation <- 
  pg_annot_complete[match(intensity_long$protein, 
                          pg_annot_complete$protein_group),]$upd_full_annot
intensity_long$organism <- 
  pg_annot_complete[match(intensity_long$protein, 
                          pg_annot_complete$protein_group),]$organism_diamond
intensity_long$exp <- meta[match(intensity_long$sample, meta$sample),]$experiment

intensity_long_subset <- subset(intensity_long, grepl('Annelida|Nema', organism))

intensity_long_subset$missed_value <- 
  ifelse(is.na(intensity_long_subset$intensity), 1, 0)

dat_to_plot <- aggregate(missed_value ~ sample, intensity_long_subset, sum)
dat_to_plot$exp <- meta[match(dat_to_plot$sample, meta$sample),]$experiment
dat_to_plot <- dat_to_plot[order(dat_to_plot$exp),]
dat_to_plot$exp <- as.factor(dat_to_plot$exp)
dat_to_plot$sample <- factor(dat_to_plot$sample, levels = dat_to_plot$sample)

ggplot(dat_to_plot, aes(sample, missed_value)) +
  geom_point(aes(color = exp)) +
  geom_line(aes(group = exp, color = exp)) +
  scale_y_continuous(limits = c(0, 24)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(current_dir, 'annelide_NAs_distribution.png'), scale = 1.9)

# Eve: it seems that the number of NAs only depends on TMT batch. It seems that we do not have Annelides in our samples



