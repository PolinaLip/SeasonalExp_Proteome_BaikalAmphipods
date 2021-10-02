### find same proteins from field sampling and gradual temperature decrease experiment
### Then, correlate them by lfc (Dec/Sep - 1.5/10.5)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(rstatix)

species <- 'Eve'
month1 <- 'JUN'
month2 <- 'SEP'

dir_to_save_res <- paste0('labeglo2/MS_results/Field_Grad_Comparison/', species)

# 1. Upload dataframes after EdgeR (de_separately...R)
dir_field <- paste0('labeglo2/MS_results/Field/', species,'/', species, '_refchannels_all')
proteins_field <- 
  read.csv(paste0(dir_field, '/', species,
                  '_edgeResults_', month1, 'vs', month2, '_proteinGroups.csv'),
           sep = '\t', header = T)

dir_grad <- paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/', species)

## Eve
proteins_grad <- 
  read.csv(paste0(dir_grad, '/', species,
                  'BK_edgeResults_10.5vs1.5_proteinGroups.csv'),
           sep = '\t', header = T) # Eve

## Ecy and Gla

cond1 <- '1.5'
cond1 <- '6'
cond2 <- '10.5ALL'
cond2 <- '10.5parall'
cond2 <- '10.5FP'

proteins_grad <- 
  read.csv(paste0(dir_grad, '/', species,
                  '_edgeResults_', cond1, 'vs', cond2, '_proteinGroups.csv'),
           sep = '\t', header = T) # Ecy, Gla

# 2. "Open" protein group names and sort contig names inside them, them "close" them back
new_protein_name <- lapply(str_split(proteins_field$protein, ';'), sort)
proteins_field$sorted_protein <- sapply(new_protein_name, paste0, collapse=';')

new_protein_name_grad <- lapply(str_split(proteins_grad$protein, ';'), sort)
proteins_grad$sorted_protein <- sapply(new_protein_name_grad, paste0, collapse=';')

# 3. Make matrix of all possible pairs of protein groups (for two dataframes)
protein_name_grid <- expand.grid(x=new_protein_name, y=new_protein_name_grad)

# 4. Intersect the sets of contig names (composing protein group) to find most similar protein groups: 
#protein_name_grid$match_rate <- apply(protein_name_grid, 1,
#      function(row) length(intersect(row$x, row$y)) / min(length(row$x), length(row$y)))

protein_name_grid$match_rate <- apply(protein_name_grid, 1,
      function(row) length(intersect(row$x, row$y)) / max(length(row$x), length(row$y)))

sum(protein_name_grid$match_rate > 0) # so many pairs are similar of even identical

# 5. Take only pairs with similarity > 0
#protein_name_grid_sub <- subset(protein_name_grid, match_rate > 0)
protein_name_grid_sub <- subset(protein_name_grid, match_rate == 1) # if take only complete matches in protein group names
protein_name_grid_sub$field_pg <- sapply(protein_name_grid_sub$x, paste0, collapse=';')
protein_name_grid_sub$grad_pg <- sapply(protein_name_grid_sub$y, paste0, collapse=';')

# 6. Assign annotations to protein groups
protein_name_grid_sub$protein_name_field <- 
  proteins_field[match(protein_name_grid_sub$field_pg, 
                       proteins_field$sorted_protein),]$geneSymbol

protein_name_grid_sub$protein_name_grad <- 
  proteins_grad[match(protein_name_grid_sub$grad_pg, 
                       proteins_grad$sorted_protein),]$geneSymbol

# 7. Some of the protein groups are duplicated but it is not possible to say 
# what is the true match. Thus, take only those pairs protein groups of which 
# meet only once through dataframe.
prot_single_entry_field <- table(protein_name_grid_sub$field_pg)
prot_single_entry_grad <- table(protein_name_grid_sub$grad_pg)

prot_single_entry_field <- prot_single_entry_field[prot_single_entry_field == 1]
prot_single_entry_grad <- prot_single_entry_grad[prot_single_entry_grad == 1]

protein_name_grid_sub_filtered <- 
  subset(protein_name_grid_sub, field_pg %in% names(prot_single_entry_field))
protein_name_grid_sub_filtered <- 
  subset(protein_name_grid_sub_filtered, grad_pg %in% names(prot_single_entry_grad))

# 8. Add LFC and pvalue to the dataframe

protein_name_grid_sub_filtered$LFC_field <- 
  proteins_field[match(protein_name_grid_sub_filtered$field_pg, 
                       proteins_field$sorted_protein),]$logFC
protein_name_grid_sub_filtered$pv_field <- 
  proteins_field[match(protein_name_grid_sub_filtered$field_pg, 
                       proteins_field$sorted_protein),]$pvalue

protein_name_grid_sub_filtered$LFC_grad <- 
  proteins_grad[match(protein_name_grid_sub_filtered$grad_pg, 
                       proteins_grad$sorted_protein),]$logFC
protein_name_grid_sub_filtered$pv_grad <- 
  proteins_grad[match(protein_name_grid_sub_filtered$grad_pg, 
                      proteins_grad$sorted_protein),]$pvalue

# 9. Clean dataframe a bit
#protein_name_grid_sub_filtered <- 
#  protein_name_grid_sub_filtered[,!colnames(protein_name_grid_sub_filtered) %in% c('x','y')]

# 10. Merge protein names (from two datasets (field and grad))  
protein_name_grid_sub_filtered$common_name <- 
  ifelse(protein_name_grid_sub_filtered$protein_name_field == protein_name_grid_sub_filtered$protein_name_grad,
         protein_name_grid_sub_filtered$protein_name_field,
         sprintf('%s/%s', protein_name_grid_sub_filtered$protein_name_field, 
                 protein_name_grid_sub_filtered$protein_name_grad))

protein_name_grid_sub_filtered$common_name <- 
  sub('PREDICTED: |-like|, partial|LOW QUALITY PROTEIN: ', 
      '', protein_name_grid_sub_filtered$common_name)

# 11. 

protein_name_grid_sub_filtered$sign <- ifelse(
  protein_name_grid_sub_filtered$pv_field < 0.05 & 
    protein_name_grid_sub_filtered$pv_grad < 0.05, 
  "< 0.05 (both)",
  ifelse(protein_name_grid_sub_filtered$pv_field < 0.05, "< 0.05 (Field)", 
         ifelse(protein_name_grid_sub_filtered$pv_grad < 0.05, 
                "< 0.05 (Lab)", "> 0.05 (both)")))

sign_changed <- subset(protein_name_grid_sub_filtered, sign == "< 0.05 (both)")
sign_changed <- subset(protein_name_grid_sub_filtered, 
                abs(LFC_field) > 0.4 & pv_field < 0.05 | 
                  abs(LFC_grad) > 0.4 & pv_grad < 0.05)

sign_changed <- subset(protein_name_grid_sub_filtered, 
  abs(LFC_field) > 0.5 & pv_field < 0.05 | abs(LFC_grad) > 0.5 & pv_grad < 0.05 |
    sign == "< 0.05 (both)")

cor_test_res <- cor.test(protein_name_grid_sub_filtered$LFC_field, 
                         protein_name_grid_sub_filtered$LFC_grad)

my_geom_text_repel <- 
  geom_text_repel(aes(label = 
                 ifelse(grepl('uncharacterized|predicted protein|hypothetical|\\*',
                 common_name), '', common_name)),
                 data = sign_changed, 
                 max.overlaps = Inf,
                 segment.colour = 'grey50', box.padding = .6, size = 2)

ggplot(protein_name_grid_sub_filtered, aes(LFC_field, LFC_grad)) +
  annotate(geom='text', x = -1.9, y = 1.3, hjust = 0, size = 5,
           label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                          'p-value ', p_format(cor_test_res$p.value,
                                               #accuracy = 0.05
                                               accuracy = cor_test_res$p.value * 10
                                               ))) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_point(aes(shape = sign, color = sign, alpha = sign)) +
  scale_shape_manual('p-value:', values = c(16, 0, 2, 3)) + # Eve, Ecy
  #scale_shape_manual('p-value:', values = c(0, 2, 3)) + # Gla
  scale_color_manual('p-value:', 
                     values = c('firebrick1', 'olivedrab','olivedrab','grey50')) + # Eve, Ecy
  #scale_color_manual('p-value:', 
  #                   values = c('olivedrab','olivedrab','grey50')) + # Gla
  scale_alpha_manual('p-value:', values = c(1, .8, .8, 0.3)) + # Eve, Ecy
  geom_smooth(method = 'lm', color = 'cornflowerblue', fill = 'grey85', size = .6) +
 
  #my_geom_text_repel +
  
  #xlab('log2FC (December, 2.3 °C/September, 10.5 °C; 71 days) ') +
  #xlab('log2FC (January, 0.3 °C/September, 10.5 °C; 118 days) ') +
  #xlab('log2FC (November, 6.4 °C/September, 10.5 °C; 38 days) ') +
  xlab('log2FC (June, ~4-6 °C/September, 10.5 °C; 254 days) ') +
  ylab('log2FC (Lab: 1.5 °C/10.5 °C; 24 days) ') +
  
  #ylab('log2FC (Lab: 1.5 °C/10.5 °C, all controls; 24 days) ') +
  #ylab('log2FC (Lab: 6 °C/10.5 °C, all controls; 12 days) ') +
  #ylab('log2FC (Lab: 1.5 °C/10.5 °C, parallel control; 24 days) ') +
  #ylab('log2FC (Lab: 6 °C/10.5 °C, parallel control; 12 days) ') +
  #ylab('log2FC (Lab: 1.5 °C/10.5 °C, first point; 12 days) ') +
  #ylab('log2FC (Lab: 6 °C/10.5 °C, first point; 12 days) ') +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  theme_light() +
  theme(axis.title = element_text(size = 8.8))

ggsave(file.path(dir, paste0('FieldGradComp_', species, '_', month1, '_', month2, 
                             '_onlyCompletePGmatches.png')), 
       scale = 1.2, width = 7, height = 5.5)

ggsave(file.path(dir_to_save_res, paste0('FieldGradComp_', species, '_', month1, '_', 
                                         month2, 
                             '_1.5vs10.5parall',
                             '_onlyCompletePGmatches.png')), 
       scale = 1.2, width = 8, height = 5.5)

ggsave(file.path(dir_to_save_res, paste0('FieldGradComp_', species, '_', month1, '_', 
                                         month2, 
                                         '_1.5vs10.5parall',
                                         '_onlyCompletePGmatches_woNames.png')), 
       scale = 0.9, width = 7, height = 5.5)

### if adjust all plots to one size
ggsave(file.path(dir_to_save_res, paste0('FieldGradComp_', species, '_', month1, '_', 
                                         month2, 
                                         '_', cond1, 'vs', cond2,
                                         '_onlyCompletePGmatches.png')), 
       scale = 1.2, width = 8, height = 3.5)

ggsave(file.path(dir_to_save_res, paste0('FieldGradComp_', species, '_', month1, '_', 
                                         month2, 
                                         '_', cond1, 'vs', cond2,
                                         '_onlyCompletePGmatches_woNames.png')), 
       scale = 0.91, width = 8.5, height = 3.5)

