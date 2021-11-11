#### Look at the chosen proteins and top proteins from a module ####
# can be used after wgcna script
library(ggplot2)
get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

toplot <- intensity_long[grepl('HSPD1|cct', intensity_long$upd_full_annotation, 
                               ignore.case = T),]
toplot <- scaled_intensities_long[grepl('Vg', scaled_intensities_long$annotation, 
                               ignore.case = T),]
toplot <- intensity_long[grepl('Vg', intensity_long$upd_full_annotation, 
                               ignore.case = T),]
toplot <- intensity_long[grepl('pheromone|ATIC|Vg|VMO1', intensity_long$upd_full_annotation, 
                               ignore.case = T),]
toplot <- intensity_long[grepl('ABC|pix', intensity_long$upd_full_annotation, 
                               ignore.case = T),]

# if want to look at top proteins of a module #
chosen_module <- module
chosen_trait <- trait
#toplot <- subset(intensity_long, module == chosen_module)
module_threshold <- 0.7
trait_threshold <- 0.4
mm_pvalue_threshold <- 0.01
gs_pvalue_threshold <- 0.01
toplot <- subset(intensity_long, MM > module_threshold & pMM < mm_pvalue_threshold)
toplot <- subset(intensity_long, GS > trait_threshold & pGS < gs_pvalue_threshold)
max(toplot$pMM)
max(toplot$pGS)
#toplot <- subset(toplot, MM > module_threshold & pMM < mm_pvalue_threshold)
###

toplot$month <- meta[match(toplot$sample, meta$sample),]$condition

toplot$month <- sub('_BK', '', toplot$month)

toplot$time <- ifelse(toplot$month == 'September', 1,
                      ifelse(toplot$month == 'November' , 2,
                             ifelse(toplot$month == 'December', 3,
                                    ifelse(toplot$month == 'January', 4, 5))))
toplot$upd_full_annotation <- sub('PREDICTED: |-like| isoform X\\d+', '', toplot$upd_full_annotation)
toplot$upd_full_annotation <- sapply(
  toplot$upd_full_annotation, function(x)
    paste(strwrap(x, width=26), collapse='\n'))

toplot$whole_label <- sprintf('%s|%s', toplot$protein,
                              toplot$upd_full_annotation)
#toplot$whole_label <- sprintf('%s|%s', toplot$proteins,
#                              toplot$annotation)
toplot$whole_label2 <- factor(toplot$whole_label,
                              levels=unique(
                              toplot$whole_label[order(toplot$upd_full_annotation)]))
#toplot$whole_label2 <- factor(toplot$whole_label,
#                              levels=unique(
#                                toplot$whole_label[order(toplot$annotation)]))


## add mean value for all proteins
## remove hypothetical proteins (with hypothetical to suppl)
toplot <- subset(toplot, !grepl('hypothetical|uncharacterized', upd_full_annotation, ignore.case = T))

median_per_month <- aggregate(intensity ~ month, toplot, median) 
mean_per_month <- aggregate(intensity ~ month, toplot, mean) 
average_per_month <- median_per_month
average_per_month$time <- ifelse(average_per_month$month == 'September', 1,
                           ifelse(average_per_month$month == 'November', 2,
                             ifelse(average_per_month$month == 'December', 3,
                              ifelse(average_per_month$month == 'January', 4, 5))))

ggplot(toplot, aes(time, intensity, group = time)) +
  facet_wrap(~ whole_label2, ncol = 4,
             labeller = as_labeller(get_protein_label)
  ) +
  #, scales = 'free_y') +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(size = 0.3
              #, aes(color = sex) # sex
  ) +
  geom_line(data = average_per_month, aes(time, intensity, group = 1), 
            color = "firebrick1", alpha = 0.80, 
            size = 1.2) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  #  scale_color_manual('Sex', values = c('deeppink1', 'dodgerblue1'), # sex 
  #                     na.value = 'grey35') + # sex
  guides(fill = F
         #, color = F) + # sex
  ) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

ggsave(file.path(current_dir, paste0('boxplots_', chosen_module, 
                                     '_more', module_threshold, 'pvalue', 
                                     mm_pvalue_threshold,
                                     '.png')),
       width = 8, height = 6, scale = 0.8)
ggsave(file.path(current_dir, paste0('boxplots_', chosen_module, 
                                     '_more', module_threshold, '_woHypothetical.png')),
       width = 6, height = 4.5, scale = 0.8)

#### plot proteins with sex info ####

toplot$sex <- meta[match(toplot$sample, meta$sample),]$sex
#toplot <- subset(toplot, intensity < 5)
#toplot <- subset(toplot, intensities < 5)
toplot$sex[toplot$sex == 'NA'] <- NA

#toplot <- 
#  subset(toplot, 
#         protein %in% unique(toplot$protein)[c(4,7,9)])

ggplot(toplot, aes(time, intensity, group = time)) +
  facet_wrap(~ whole_label2, ncol = 5,
             labeller = as_labeller(get_protein_label)
  ) +
  #, scales = 'free_y') +
  geom_boxplot(outlier.color = NA, size = 0.3) +
  geom_jitter(size = 0.8, alpha = 0.7
              , aes(color = sex) # sex
  ) +
  #geom_line(data = average_per_month, aes(time, intensity, group = 1), 
  #          color = "firebrick1", alpha = 0.80, 
  #          size = 1.2) +
  scale_x_continuous('Months', breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  scale_color_manual('Sex: ', values = c('deeppink1', 'dodgerblue1'), 
                     labels = c('Female', 'Male', 'Unknown'), 
                       na.value = 'grey35') + # sex
  guides(fill = F
         #, color = F) + # sex
  ) +
  ylab('Scaled intensities') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = 'top',
        legend.margin = margin(0,0,-10,0),
        legend.key.width = unit(5, 'pt'),
        legend.text = element_text(margin = margin(0, 5, 0, 0)))

ggsave(file.path(current_dir, paste0('boxplots_', chosen_module, 
                                     '_more', module_threshold, 'pvalue', 
                                     mm_pvalue_threshold,
                                     '_withSexAssignment.png')),
       width = 10, height = 4.8, scale = 0.8)
ggsave(file.path(current_dir, paste0('boxplots_Vg_withSexAssignment_onlyUsedForPrediction',
                                     '.png')),
       width = 6, height = 2.7, scale = 0.8)

#### Plot proteins with predicted sex ####

toplot$sex <- meta[match(toplot$sample, 
                   meta$sample),]$sex2
#toplot$sex[is.na(toplot$sex)] <- 'male'

toplot_f <- subset(toplot, sex == 'female')
toplot_m <- subset(toplot, sex == 'male')
months_mean_female <- aggregate(intensity ~ month, toplot_f, median)
months_mean_male <- aggregate(intensity ~ month, toplot_m, median)
months_mean_female$time <- ifelse(months_mean_female$month == 'September', 1,
                                  ifelse(months_mean_female$month == 'November', 2,
                                         ifelse(months_mean_female$month == 'December', 3,
                                                ifelse(months_mean_female$month == 'January', 4, 5))))
months_mean_male$time <- ifelse(months_mean_male$month == 'September', 1,
                                ifelse(months_mean_male$month == 'November', 2,
                                       ifelse(months_mean_male$month == 'December', 3,
                                              ifelse(months_mean_male$month == 'January', 4, 5))))

toplot <- subset(toplot, intensity > -2 & intensity < 4) 

#intensity_long_toplot$sex <- meta[match(intensity_long_toplot$sample,
#                                        meta$sample),]$sex2

toplot <- subset(toplot, upd_full_annotation != '*')

(ggplot(toplot, aes(factor(time), intensity)) +
    facet_wrap(~ whole_label2, ncol = 3,
               labeller = as_labeller(get_protein_label), 
               #scales = 'free_y') +
    )+
    geom_boxplot(outlier.color = NA, aes(color = sex)) +
    geom_jitter(size = 0.3,
                aes(color = sex), # sex
                position=position_jitterdodge(jitter.width = 0.3)) +
    geom_line(data = months_mean_female, aes(time, intensity, group = 1), 
              color = "deeppink1", alpha = 0.40, 
              size = 1.2) +
    geom_line(data = months_mean_male, aes(time, intensity, group = 1), 
              color = "dodgerblue1", alpha = 0.40, 
              size = 1.2) +
    scale_x_discrete(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
    scale_color_manual('Sex', values = c('deeppink1', 'dodgerblue1'), # sex 
                       na.value = 'grey35') + # sex
    #  guides(fill = F
    #         , color = F) + # sex
    #) +
    xlab('Months') +
    ylab('Scaled intensities') +
    theme_bw() +
    theme(axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          legend.position = 'bottom',
          legend.margin = margin(-10, 0, 0, 0)))

ggsave(file.path(current_dir, paste0('boxplots_', chosen_module, 
                                     '_more', module_threshold, 'pvalue', 
                                     mm_pvalue_threshold,
                                     '_withPredSex.png')),
       width = 6, height = 3.2, scale = 0.9)

ggsave(file.path(current_dir, paste0('boxplots_', chosen_module, 
                                     '_more', module_threshold, 'pvalue', 
                                     mm_pvalue_threshold,
                                     '_withPredSex_withJUNE.png')),
       width = 6, height = 4, scale = 0.9)

ggsave(file.path(current_dir, paste0('boxplots_', chosen_trait, 
                                     '_more', trait_threshold, 'pvalue', 
                                     gs_pvalue_threshold,
                                     '_withPredSex.png')),
       width = 8, height = 5, scale = 0.9)

ggsave(file.path(current_dir, paste0('boxplots_', 'ABCpix_',
                                     '_withPredSex_withJune.png')),
       width = 3, height = 2.2, scale = 0.9)

#write.table(x = toplot, file.path(current_dir, 'females_specific_proteins_pinkmodule_intensities.csv'))


