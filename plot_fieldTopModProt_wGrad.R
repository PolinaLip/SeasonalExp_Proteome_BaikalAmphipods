#### to plot top proteins from the field modules with intensities from gradual decrease experiment
library(stringr)

order_recode <- function(x, ...) {
  args <- unlist(rlang::list2(...))
  old_levels = as.vector(args)
  new_levels = names(args)
  new_levels[new_levels == ""] <- old_levels[new_levels == ""]
  
  res <- factor(x, levels = old_levels, labels = new_levels)
  stopifnot('Not all factor levels are supplied' = all(is.na(x) == is.na(res)))
  res
}

dir_to_save_results <- 'labeglo2/MS_results/Field_Grad_Comparison/'

### Upload metafile for gradual decrease experiment
path2meta_grad <-
  paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/',
         species, '/Metadata_Proteus_withPredSex.csv')

meta_g <- meta_upload(path2meta_grad, species)
meta_g <- subset(meta_g, !grepl('PB', sample))
meta_g <- subset(meta_g, !grepl('Pool', sample))

meta_grad <- meta_g

### Upload intensities for gradual decrease experiment
dir_to_grad <-paste0('labeglo2/MS_results/GradDecrease_2/separately_assembled/', species)
data_grad <- read.table(file.path(dir_to_grad, 
                                  paste0('intensities_after_slNorm_', 
                                         tolower(species), '.csv')), header = T)
rownames(data_grad) <- data_grad[,1] # if take data with NAs
data_grad <- data_grad[,-1] # if take data with NAs

data_grad <- read.table(file.path(dir_to_grad, 
                                  paste0('intensities_after_slNorm_woNA_', # i take data wo NAs in the case of grad decrease exp
                                         tolower(species), '.csv')), header = T)

### scale intensities from gradual experiment

scaled_intensities_grad <- apply(data_grad, 1, scale)
scaled_intensities_grad <- t(scaled_intensities_grad)
colnames(scaled_intensities_grad) <- meta_grad$sample
scaled_intensities_grad <- as.data.frame(scaled_intensities_grad)

###########################
toplot_ <- 
  toplot[,colnames(toplot) %in% c('protein', 'intensity', 'month', 'sex')]
colnames(toplot_) <- c('protein', 'intensity', 'condition', 'sex')

### choose only proteins to plot
sorted_protein_list_field <- str_split(unique(toplot_$protein), ';')
  # lapply(str_split(unique(toplot_$protein), ';'), sort)

sorted_protein_list_grad <- str_split(rownames(scaled_intensities_grad), ';')
  # lapply(str_split(rownames(scaled_intensities_grad), ';'), sort)
#scaled_intensities_grad$sorted_protein <- 
#  sapply(sorted_protein_list_grad, paste0, collapse=';')

protein_name_grid <- expand.grid(x=sorted_protein_list_field, 
                                 y=sorted_protein_list_grad)

protein_name_grid$match_rate <- 
  apply(protein_name_grid, 1,
        function(row) length(intersect(row$x, row$y)))

protein_name_grid_match <- filter(protein_name_grid, match_rate > 0) %>%
  mutate(x_name = sapply(x, paste0, collapse = ';'),
         y_name = sapply(y, paste0, collapse = ';'))

sum(protein_name_grid$match_rate > 0)
prot_to_plot <- sapply(protein_name_grid[protein_name_grid$match_rate > 0,]$y,
                       paste0, collapse=';')

scaled_intensities_grad_toplot <- 
  subset(scaled_intensities_grad, rownames(scaled_intensities_grad) %in% prot_to_plot)

grad_toplot_long <- scaled_intensities_grad_toplot %>%
  rownames_to_column(var = 'init_protein') %>%
  pivot_longer(cols = -init_protein, values_to = 'intensity', names_to = 'condition')

grad_toplot_long$protein <- 
  protein_name_grid_match[match(grad_toplot_long$init_protein, 
                                protein_name_grid_match$y_name),]$x_name

grad_toplot_long_ <- grad_toplot_long[,c('protein', 'intensity', 'condition')]


#### 
grad_toplot_long_ <- aggregate(intensity ~ protein + condition, mean, 
                                data = grad_toplot_long_)
grad_toplot_long_ <- grad_toplot_long_[,c(1,3,2)]
#### add info about sex

grad_toplot_long_$sex <- 
  meta_grad[match(grad_toplot_long_$condition, meta_grad$sample),]$sex

grad_toplot_long_$condition <- 
  meta_grad[match(grad_toplot_long_$condition, meta_grad$sample),]$condition
grad_toplot_long_$condition <- sub('EveBK_', '', grad_toplot_long_$condition)

grad_toplot_long_$condition <- sub('5', '5 °C', grad_toplot_long_$condition)

### 

all_toplot <- rbind(toplot_, grad_toplot_long_)
all_toplot$annotation <- 
  pg_annot_complete[match(all_toplot$protein, 
                          pg_annot_complete$protein_group),]$upd_full_annot
###

all_toplot$condition <- order_recode(all_toplot$condition,
    'Sep' = 'September', 'Nov' = 'November', 'Dec' = 'December',
    'Jan' = 'January', 'Jun' = 'June', '10.5 °C', '1.5 °C')
all_toplot$annotation <- sub('PREDICTED: |-like| isoform X\\d+|, partial', '', 
                             all_toplot$annotation)

all_toplot$annotation <- sapply(
  all_toplot$annotation, function(x)
    paste(strwrap(x, width=24), collapse='\n'))

all_toplot$whole_label <- sprintf('%s|%s', all_toplot$protein,
                                  all_toplot$annotation)
#toplot$whole_label <- sprintf('%s|%s', toplot$proteins,
#                              toplot$annotation)
all_toplot$whole_label2 <- factor(all_toplot$whole_label,
                              levels=unique(
                                all_toplot$whole_label[order(all_toplot$annotation)]))

(ggplot(all_toplot, aes(condition, intensity)) +
    facet_wrap(~ whole_label2, ncol = 4,
               labeller = as_labeller(get_protein_label), 
               #scales = 'free_y') +
    )+
    geom_boxplot(outlier.color = NA, aes(fill = sex), color = 'grey45', alpha = 0.3) +
    geom_jitter(size = 0.3,
                aes(color = sex), # sex
                position=position_jitterdodge(jitter.width = 0.3)) +
    geom_vline(xintercept = 5.5, linetype="dotted") +
    geom_line(data = months_mean_female, aes(time, intensity, group = 1), 
              color = "deeppink1", alpha = 0.40, 
              size = 1.2) +
    geom_line(data = months_mean_male, aes(time, intensity, group = 1), 
              color = "dodgerblue1", alpha = 0.40, 
              size = 1.2) +
    #geom_line(data = temp_mean_female, aes(time, intensity, group = 1), 
    #          color = "deeppink1", alpha = 0.40, 
    #          size = 1.2) +
    #geom_line(data = temp_mean_male, aes(time, intensity, group = 1), 
    #          color = "dodgerblue1", alpha = 0.40, 
    #          size = 1.2) +
    scale_fill_manual('Sex: ', values = c('deeppink1', 'dodgerblue1'), # sex 
                      labels = c('Female', 'Male', 'Unknown'),
                      na.value = 'grey35') + # sex
    scale_color_manual('Sex: ', values = c('deeppink1', 'dodgerblue1'), # sex 
                       labels = c('Female', 'Male', 'Unknown'),
                       na.value = 'grey35') +
    #  guides(fill = F
    #         , color = F) + # sex
    #) +
    ylab('Scaled intensities') +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 13),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          legend.position = 'top',
          legend.margin = margin(0,0,-10,0),
          legend.key.width = unit(5, 'pt'),
          legend.text = element_text(margin = margin(0, 5, 0, 0))))

ggsave(file.path(dir_to_save_results,
                 paste0('boxplots_', chosen_module, 
                        '_more', module_threshold, '_MAXpvalue', 
                        max(toplot$pMM),
                        '_fieldwNA_GradWoNAs',
                        '_withPredSex_withJUNE.png')),
       width = 9, height = 7, scale = 0.87)

# 4x2: width = 9, height = 4.7, scale = 0.87

