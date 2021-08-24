library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
species <- 'Gla'
current_dir <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species,'_refchannels_all/')
current_dir <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species,'_refchannels_all/BK/') # Eve

data_for_analysis <- read.csv(file.path(current_dir, 
                                     paste0('intensities_after_slNorm_woNA_', 
                                            tolower(species),'.csv')),
                           header = T, sep = '\t')

meta_complete <- read.delim(file.path(current_dir, 'Metadata_Proteus.tsv'), 
                            header=TRUE, sep="\t")

pg_annot_complete <- read.csv(file.path(current_dir, 
                              paste0('annot_proteinGroups_field_', 
                                     tolower(species), 
                                     '_withGOannotation.csv')), 
                              sep = '\t', header = T)

data_for_analysis <- data_for_analysis[, !grepl('pool', colnames(data_for_analysis),
                                          ignore.case = T)]
meta_complete <- subset(meta_complete, condition != 'pool')
meta_complete <- meta_complete[grepl('BK', meta_complete$sample),] 
meta_complete <- subset(meta_complete, sample != 'F/100/BK/3c/2m_2')

colnames(data_for_analysis) <- meta_complete$sample

intensity_scaled <- apply(data_for_analysis, 1, scale)
intensity_scaled <- t(intensity_scaled)
intensity_scaled <- data.frame(intensity_scaled)
colnames(intensity_scaled) <- meta_complete$sample

intensity_long <- intensity_scaled %>% 
  rownames_to_column(var = 'protein') %>%
  pivot_longer(-protein, values_to = 'intensity', names_to = 'sample')

intensity_long$annotation <- pg_annot_complete[match(intensity_long$protein,
                                                     pg_annot_complete$protein_group),]$annotation


intensity_long$month <- meta_complete[match(intensity_long$sample, 
                                            meta_complete$sample),]$condition
intensity_long_sd <- aggregate(intensity ~ protein + month, intensity_long,
                               function(x) { r <- range(x); r[2] - r[1]} )

ggplot(intensity_long_sd) +
  geom_boxplot(aes(factor(month), intensity))

ggplot(intensity_long_sd, aes(intensity)) +
  facet_wrap(~ month, ncol = 5) +
  geom_histogram()

### to draw all three species
intensity_long_sd$species <- species
#all_species <- NULL
all_species <- rbind(all_species, intensity_long_sd)
all_species$month <- sub('_BK', '', all_species$month)

ggplot(all_species, aes(factor(month), intensity)) +
  facet_wrap(~ species, ncol = 1) +
  geom_boxplot()

ggplot(all_species, aes(intensity)) +
  facet_wrap(~ species + month, ncol = 5) +
  geom_histogram(binwidth=.2)

### Brown–Forsythe test to test equality of group variances
library(onewaytests)

outp_df <- NULL
for (pr in unique(intensity_long$protein)) {
  current_data <- subset(intensity_long, protein == pr)
  current_out <- bf.test(intensity ~ month, current_data, verbose = F)
  curr_row <- data.frame(protein = pr,
                         pvalue = current_out$p.value)
  outp_df <- rbind(outp_df, curr_row)
}

outp_df$padj <- p.adjust(outp_df$pvalue, method = 'fdr')
outp_sign <- subset(outp_df, padj < 0.05)
outp_sign$annotation <- 
  pg_annot_complete[match(outp_sign$protein, 
                          pg_annot_complete$protein_group),]$annotation
mean(outp_sign$padj)

to_plot <- subset(intensity_long, protein == outp_sign$protein[1])

ggplot(to_plot, aes(month, intensity)) +
  geom_boxplot()
