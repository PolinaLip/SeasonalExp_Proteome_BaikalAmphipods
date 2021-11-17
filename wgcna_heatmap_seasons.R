### Plot heatmap for wgcna correlation plot
# for the plotting, wgcna script should be running until the moment with correlation plots
library(ggthemes)

dir_to_save <- 'labeglo2/MS_results/Field/Eve/Eve_refchannels_all/withALLzeros_inAtLeastOneCond/'

moduleTraitCor_df <- data.frame(moduleTraitCor)
moduleTraitPvalue_df <- data.frame(moduleTraitPvalue) 

mtc_long <- moduleTraitCor_df %>%
  rownames_to_column(var = 'Module') %>%
  pivot_longer(-Module, values_to = 'Correlation', names_to = 'Condition')

mtp_long <- moduleTraitPvalue_df %>%
  rownames_to_column(var = 'Module') %>%
  pivot_longer(-Module, values_to = 'p-value', names_to = 'Condition')

mtc_long$pvalue <- mtp_long$`p-value`

mtc_long$Condition <- factor(mtc_long$Condition, levels = unique(mtc_long$Condition))

mtc_long$Module <- sub('ME', '', mtc_long$Module)
mtc_long$Module <- factor(mtc_long$Module, levels = rev(unique(mtc_long$Module)))

mtc_long$pvalue_ast <- ifelse(mtc_long$pvalue <= 0.1 & mtc_long$pvalue > 0.05, '*', 
                       ifelse(mtc_long$pvalue <= 0.05 & mtc_long$pvalue > 0.01, '**', 
                       ifelse(mtc_long$pvalue <= 0.01 & mtc_long$pvalue > 0.001, '***', 
                              ifelse(mtc_long$pvalue <= 0.001, '****', ''))))

mtc_long$text_label <- ifelse(mtc_long$pvalue_ast != '',
                              paste0(round(mtc_long$Correlation, digits = 2), 
                                     '\n', 
                                     mtc_long$pvalue_ast), 
                              #paste0(round(mtc_long$Correlation, digits = 2), '\n'))
                              '')

mtc_long <- subset(mtc_long, Module != 'grey')
mtc_long$Condition <- factor(mtc_long$Condition, levels = levels(mtc_long$Condition),
                             #labels = c(levels(mtc_long$Condition)[1:6], # wo June
                             labels = c(levels(mtc_long$Condition)[1:7], 
                                        intToUtf8(9792),
                                        intToUtf8(9794), 
                             paste0('Sep, ', intToUtf8(9792)),
                             paste0('Sep, ', intToUtf8(9794)), 
                             paste0('Nov, ', intToUtf8(9792)),
                             paste0('Nov, ', intToUtf8(9794)), 
                             paste0('Dec, ', intToUtf8(9792)),
                             paste0('Dec, ', intToUtf8(9794)), 
                             paste0('Jan, ', intToUtf8(9792)),
                             paste0('Jan, ', intToUtf8(9794))
                             #,
                             #paste0('Jun, ', intToUtf8(9794))
                             #,
                             #paste0('Amplexus, ', intToUtf8(9792)),
                             #paste0('Amplexus, ', intToUtf8(9794)))
                             ))

### plot all modules
ggplot(mtc_long, aes(Condition, Module, fill = Correlation)) +
  geom_tile() +
  #geom_tile(color = 'grey10') +
  geom_text(aes(label = text_label), size = 3.7) +
  ylab('Modules') +
  xlab('Traits') +
  scale_fill_gradientn('Correlation:  ', colors = blueWhiteRed(50), limits = c(-1, 1)) +
  # scale_fill_gradientn(colors = blueWhiteRed(4)[c(1, 1, 2, 4, 4)], limits = c(-1, 1)) +
  scale_x_discrete(expand = expansion(add=0.51)) +
  theme_tufte(base_family="Helvetica") + 
  theme(axis.ticks=element_blank(),
        axis.text=element_text(size=15),
        axis.title = element_text(size=16),
        legend.title = element_text(size=13),
        legend.text = element_text(size=10),
        legend.position = 'top',
        legend.justification = 'left',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_colourbar(barwidth = 7,
                                barheight = 0.8))

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', MCH, '_MMS', MMS,'_DS', DS,
                         '_', species, 'withATLEAST5notNAsInCond_withJune.png'), 
       scale = 0.9, width = 11, height = 8)
ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', MCH, '_MMS', MMS,'_DS', DS,
                         '_', species, 'withATLEAST5notNAsInCond_withPredSex_AMPL.png'), 
       scale = 0.9, width = 12, height = 9)
ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', MCH, '_MMS', MMS,'_DS', DS,
                         '_', species, 'withATLEAST5notNAsInCond_withPredSex_AMPL_withJUNE.png'), 
       scale = 0.9, width = 13, height = 8)

### plot only sign modules

sign_mtc_logical <- mtc_long %>%
  group_by(Module) %>%
  summarise(has_sign = any(grepl('\\*{3}', pvalue_ast))) %>%
  filter(has_sign)

sign_mtc <- mtc_long[mtc_long$Module %in% sign_mtc_logical$Module,]
#sign_mtc <- subset(mtc_long, Module == 'pink' | Module == 'cyan' |
#                     Module == 'lightcyan' | Module == 'darkgreen' |
#                     Module == 'brown' | Module == 'blue')

ggplot(sign_mtc, aes(Condition, Module, fill = Correlation)) +
  geom_tile() +
  #geom_tile(color = 'grey10') +
  geom_text(aes(label = text_label), size = 3.7) +
  geom_vline(xintercept = 2.5) +
  geom_vline(xintercept = 7.5) +
  geom_vline(xintercept = 9.5) +
  geom_vline(xintercept = 18.5) + 
  ylab('Modules') +
  xlab('Traits') +
  scale_fill_gradientn('Correlation:  ', colors = blueWhiteRed(50), limits = c(-1, 1)) +
  # scale_fill_gradientn(colors = blueWhiteRed(4)[c(1, 1, 2, 4, 4)], limits = c(-1, 1)) +
  scale_x_discrete(expand = expansion(add=0.51)) +
  theme_tufte(base_family="Helvetica") + 
  theme(axis.ticks=element_blank(),
        axis.text=element_text(size=15),
        axis.title = element_text(size=16),
        legend.title = element_text(size=13),
        legend.text = element_text(size=10),
        legend.position = 'top',
        legend.justification = 'left',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_colourbar(barwidth = 7,
                                barheight = 0.8))

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', 0.15, '_MMS', 25,'_DS', 4,
                         '_', species, 'onlySignModules_less001_4.png'), 
       scale = 0.9, width = 5.8, height = 5.5)

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', MCH, '_MMS', MMS,'_DS', DS,
                         '_', species, 'onlySignModules_less001_withPredSex_AMPL2.png'), 
       scale = 0.9, width = 12, height = 9)

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', MCH, '_MMS', MMS,'_DS', DS,
                         '_', species, 'onlySignModules_less001_withPredSex_withJUNE.png'), 
       scale = 0.9, width = 10.5, height = 8)

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', 0.15, '_MMS', 25,'_DS', 4,
                         '_', species, 'choosenModules.png'), 
       scale = 0.9, width = 3.7, height = 4.5)


