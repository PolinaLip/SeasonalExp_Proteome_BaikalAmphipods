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
mtc_long$pvalue_ast <- ifelse(mtc_long$pvalue <= 0.2 & mtc_long$pvalue > 0.05, '*', 
                       ifelse(mtc_long$pvalue <= 0.05 & mtc_long$pvalue > 0.01, '**', 
                       ifelse(mtc_long$pvalue <= 0.01, '***', '')))

mtc_long$text_label <- ifelse(mtc_long$pvalue_ast != '',
                              paste0(round(mtc_long$Correlation, digits = 2), 
                                     '\n', 
                                     mtc_long$pvalue_ast), 
                              paste0(round(mtc_long$Correlation, digits = 2), '\n'))
mtc_long <- subset(mtc_long, Module != 'grey')
### plot all modules
ggplot(mtc_long, aes(Condition, Module, fill = Correlation)) +
  #geom_tile() +
  geom_tile(color = 'grey10') +
  geom_text(aes(label = text_label), size = 3.7) +
  scale_fill_gradientn(colors = blueWhiteRed(50), limits = c(-1, 1)) +
  # scale_fill_gradientn(colors = blueWhiteRed(4)[c(1, 1, 2, 4, 4)], limits = c(-1, 1)) +
  scale_x_discrete(expand = expansion(add=0.51)) +
  theme_tufte(base_family="Helvetica") + 
  theme(axis.ticks=element_blank(),
        axis.text=element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size=14))

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', 0.15, '_MMS', 25,'_DS', 4,
                         '_', species, '.png'), 
       scale = 0.9, width = 6.5, height = 8)

### plot only sign modules

sign_mtc <- subset(mtc_long, pvalue_ast != '')
sign_mtc <- subset(mtc_long, Module == 'pink' | Module == 'cyan' |
                     Module == 'lightcyan' | Module == 'darkgreen' |
                     Module == 'brown' | Module == 'blue')

ggplot(sign_mtc, aes(Condition, Module, fill = Correlation)) +
  geom_tile(color = 'grey10') +
  geom_text(aes(label = text_label), size = 3.7) +
  scale_fill_gradientn(colors = blueWhiteRed(50), limits = c(-1, 1)) +
  # scale_fill_gradientn(colors = blueWhiteRed(4)[c(1, 1, 2, 4, 4)], limits = c(-1, 1)) +
  scale_x_discrete(expand = expansion(add=0.51)) +
  theme_tufte(base_family="Helvetica") + 
  theme(axis.ticks=element_blank(),
        axis.text=element_text(size=12),
        axis.title = element_text(size=15),
        legend.title = element_text(size=10)) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 4))

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', 0.15, '_MMS', 25,'_DS', 4,
                         '_', species, 'onlySignModules.png'), 
       scale = 0.9, width = 3.7, height = 2.4)

ggsave(filename = paste0(dir_to_save, 'cor_heatmap_power', 
                         sth_power,'_MCH', 0.15, '_MMS', 25,'_DS', 4,
                         '_', species, 'choosenModules.png'), 
       scale = 0.9, width = 3.7, height = 4.5)


