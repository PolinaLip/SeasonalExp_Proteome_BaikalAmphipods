library(stringr)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(ggfortify) 
library(rstatix)
get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}
dir_to_save_results <- 'labeglo2/MS_results/Field_Grad_Comparison/'

species <- 'Eve'
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}

### Upload metafile for field
path2meta <- paste0('labeglo2/MS_results/Field/', species,
                    '/', species, '_refchannels_all', '/metadata_predictedSex.csv')

meta <- meta_upload(path2meta, species)
meta <- meta[!grepl('PB|pool_12', meta$sample),]
meta <- subset(meta, condition != 'pool')
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')
#meta_field <- meta[-6]
meta_field <- meta[-c(6,7)]
meta_field$experiment <- as.numeric(meta_field$experiment)
meta_field <- subset(meta_field, condition != 'pool')

### Upload metafile for gradual decrease experiment
path2meta_grad <-
  paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/',
         species, '/Metadata_Proteus_withPredSex.csv')

meta_g <- meta_upload(path2meta_grad, species)
meta_g <- subset(meta_g, !grepl('PB', sample))
meta_g <- subset(meta_g, !grepl('Pool', sample))

meta_grad <- meta_g

### Combine metafiles
meta_field$sex <- meta_field$sex2 
meta_field <- meta_field[-6]
meta_joined <- rbind(meta_field, meta_grad)
meta_joined$condition <- sub('_BK|EveBK_', '', meta_joined$condition)

### Upload intensities for field sampling
dir_to_field <- paste0('labeglo2/MS_results/Field/', species, '/', species, 
                       '_refchannels_all')
data_field <- read.table(file.path(dir_to_field, 
                                   paste0('intensities_after_slNorm_withNA', 
                                          tolower(species), '.csv')), header = T)

rownames(data_field) <- data_field$protein_group
data_field <- data_field[-1]
data_field <- data_field[,!grepl('pool', colnames(data_field))]
colnames(data_field) <- meta_field$sample

### Upload intensities for gradual decrease experiment
dir_to_grad <-paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/', species)
data_grad <- read.table(file.path(dir_to_grad, 
                                  paste0('intensities_after_slNorm_woNA_', # i take data wo NAs in the case of grad decrease exp
                                         tolower(species), '.csv')), header = T)
#rownames(data_grad) <- data_grad$protein_group
#data_grad <- data_grad[-1]

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

#### Sort each protein group name (inside itself to get as much equal protein groups as possible)  

new_protein_name_field <- 
  lapply(str_split(rownames(scaled_intensities_field), ';'), sort)
scaled_intensities_field$sorted_protein <- 
  sapply(new_protein_name_field, paste0, collapse=';')

new_protein_name_gradual <- 
  lapply(str_split(rownames(scaled_intensities_grad), ';'), sort)
scaled_intensities_grad$sorted_protein <- 
  sapply(new_protein_name_gradual, paste0, collapse=';')

#### Try oversect the pg by the first name in it
splitted_name_field <- lapply(rownames(scaled_intensities_field), str_split, ';')
splitted_name_field <- sapply(splitted_name_field, "[[", 1)
scaled_intensities_field$first_name <- sapply(splitted_name_field, "[[", 1)
scaled_intensities_field$annotation <- 
  prot_annot_field[match(rownames(scaled_intensities_field), 
                   prot_annot_field$protein_group),]$upd_full_annot
  
splitted_name_grad <- lapply(rownames(scaled_intensities_grad), str_split, ';')
splitted_name_grad <- sapply(splitted_name_grad, "[[", 1)
scaled_intensities_grad$first_name <- sapply(splitted_name_grad, "[[", 1)

data_joined <- inner_join(scaled_intensities_field, scaled_intensities_grad, 
                          by = c("first_name" = "first_name"))
rownames(data_joined) <- data_joined$first_name

data_joined <- data_joined[,colnames(data_joined) != c('sorted_protein.x','sorted_protein.y')]
data_joined_long <- data_joined %>%
  pivot_longer(-c(first_name, annotation), values_to = 'intensity', names_to = 'sample')

data_joined_long$Condition <- 
  meta_joined[match(data_joined_long$sample, meta_joined$sample),]$condition

data_joined_long$Sex <- 
  meta_joined[match(data_joined_long$sample, meta_joined$sample),]$sex
### To plot chosen proteins

toplot <- data_joined_long[grepl('ATIC|Cpq|pheromone|PRSS21|lipid-transport|Vg|ALDH4A1|F12|Gstm2|TKT', 
                               data_joined_long$annotation, 
                               ignore.case = T),]
toplot <- data_joined_long[grepl('Amyrel|14-3-3|CAPN1|DLST1|NDRG|SDCBP|smc|SYNE1|UTP18|H1F0|OGDHL|UBC|Rab6', 
                                 data_joined_long$annotation, 
                                 ignore.case = T),]

toplot$Condition <- factor(toplot$Condition, levels = unique(toplot$Condition),
                           labels = c('Sep', 'Nov', 'Dec', 'Jan', 'June', 
                                      '10.5 °C', '1.5 °C'))
toplot$annotation <- sub('PREDICTED: |-like| isoform X\\d+|, partial', '', 
                                  toplot$annotation)

toplot$annotation <- sapply(
  toplot$annotation, function(x)
    paste(strwrap(x, width=24), collapse='\n'))

toplot$whole_label <- sprintf('%s|%s', toplot$first_name,
                              toplot$annotation)
#toplot$whole_label <- sprintf('%s|%s', toplot$proteins,
#                              toplot$annotation)
toplot$whole_label2 <- factor(toplot$whole_label,
                              levels=unique(
                                toplot$whole_label[order(toplot$annotation)]))
#toplot <- subset(toplot, !grepl('uncharacterized', upd_full_annotation))

ggplot(toplot, aes(Condition, intensity, group = Condition)) +
  facet_wrap(~ whole_label2, ncol = 6, 
             labeller = as_labeller(get_protein_label)) +
  geom_boxplot(outlier.color = NA 
  ) +
               #,aes(color = temp)) +
               #aes(fill = temp), alpha = 0.3, color = 'grey70') + # Vg
  geom_jitter(
    #aes(color = temp)
    ) +
  ylab('Scaled intensities') +
  #scale_color_manual('Temperature',
  #                   values = as.character(palette.colors(n=9)[c(7, 3)])) +
  #scale_fill_manual('Temperature',
  #                  values = as.character(palette.colors(n=9)[c(7, 3)])) + # Vg
  theme_bw() +
  guides(color = F, fill = F) +
  theme(axis.title.x = element_blank())

toplot <- subset(toplot, intensity < 4) 
(ggplot(toplot, aes(factor(Condition), intensity)) +
    facet_wrap(~ whole_label2, ncol = 4,
               labeller = as_labeller(get_protein_label), 
               #scales = 'free_y') +
    )+
    geom_boxplot(outlier.color = NA, aes(fill = Sex), color = 'grey57', alpha = 0.3) +
    geom_jitter(size = 0.9,
                aes(color = Sex), # sex
                position=position_jitterdodge(jitter.width = 0.3)) +
    geom_vline(xintercept = 5.5, linetype="dotted") +
    #geom_line(data = temp_mean_female, aes(time, intensity, group = 1), 
    #          color = "deeppink1", alpha = 0.40, 
    #          size = 1.2) +
    #geom_line(data = temp_mean_male, aes(time, intensity, group = 1), 
    #          color = "dodgerblue1", alpha = 0.40, 
    #          size = 1.2) +
    scale_fill_manual('Sex', values = c('deeppink1', 'dodgerblue1'), # sex 
                      na.value = 'grey35') + # sex
    scale_color_manual('Sex', values = c('deeppink1', 'dodgerblue1'), # sex 
                       na.value = 'grey35') +
    #  guides(fill = F
    #         , color = F) + # sex
    #) +
    ylab('Scaled intensities') +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 13),
          legend.position = 'bottom',
          legend.margin = margin(-10, 0, 0, 0)))

ggsave(file.path(dir_to_save_results,
                 paste0('boxplots_topGreenModuleGrad', 
                        #chosen_module, 
                        #'_more', module_threshold, 'pvalue', 
                        #mm_pvalue_threshold,
                        '_GradwoNAs_FieldwithNAs.png')),
       width = 13.5, height = 8, scale = 0.9)

#### Make a grid from all field pg names against all pg names from grad decrease experiment  

protein_name_grid_for_pca <- expand.grid(x=new_protein_name_field, 
                                         y=new_protein_name_gradual)

####
#protein_name_grid_for_pca$match_rate <- apply(protein_name_grid_for_pca, 1,
#  function(row) length(intersect(row$x, row$y)) / min(length(row$x), length(row$y)))
protein_name_grid_for_pca$match_rate <- 
  apply(protein_name_grid_for_pca, 1,
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
