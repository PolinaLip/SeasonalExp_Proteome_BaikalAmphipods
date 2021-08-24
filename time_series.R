### Performing kmeans clustering and building plots based on the obtained clusters ###

library(TSrepr)
library(ggplot2)
library(data.table)
library(dplyr)
library(topGO)
library(readr)

species <- "Eve"
population <- "BK"
dir <- paste0('labeglo2/MS_results/Field/', species, '/', 
              species, '_refchannels_all/', population) # Eve
dir <- paste0('labeglo2/MS_results/Field/', species, '/', species, '_refchannels_all')
data_wo_na <- read.csv(file.path(dir, 
                        paste0("intensities_after_slNorm_woNA_", tolower(species), ".csv")),
                        sep = '\t')
### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}
path2meta <-
  paste0('labeglo2/MS_results/Field/', species, '/', 
         species,'_refchannels_all/Metadata_Proteus.tsv')

meta <- meta_upload(path2meta, species)
meta[which(meta$sex == ''),]$sex <- 'NA'
meta <- meta[!grepl('PB|pool_12', meta$sample),]
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')

# Prepare data 
intensities <- data_wo_na # from de_separately_field.R
scaled_intensities <- apply(intensities, 1, scale)
scaled_intensities <- t(scaled_intensities)
colnames(scaled_intensities) <- meta$sample
data_time_series <- scaled_intensities[,!grepl('pool', colnames(scaled_intensities))]

# if you want to use only proteins significantly changed at least once:
load('~/labeglo2/MS_results/Field/Eve/BK/sign_changed_proteins_allmonths.rda')
data_time_series <- 
  data_time_series[row.names(data_time_series) %in% protein_changed_unique,]

# Obtain medians for every month
ts_1 <- data_time_series[,subset(meta, 
                                 condition == paste0('September_', population))$sample]
ts_2 <- data_time_series[,subset(meta, 
                                 condition == paste0('November_', population))$sample]
ts_3 <- data_time_series[,subset(meta, 
                                 condition == paste0('December_', population))$sample]
ts_4 <- data_time_series[,subset(meta, condition == 'January_BK')$sample]
ts_5 <- data_time_series[,subset(meta, condition == 'June_BK')$sample]

ts_1 <- data_time_series[,subset(meta, 
                                 condition == 'September')$sample]
ts_2 <- data_time_series[,subset(meta, 
                                 condition == 'November')$sample]
ts_3 <- data_time_series[,subset(meta, 
                                 condition == 'December')$sample]
ts_4 <- data_time_series[,subset(meta, condition == 'January')$sample]
ts_5 <- data_time_series[,subset(meta, condition == 'June')$sample]

ts_medians <- data.frame("M1_Sep" = apply(ts_1, 1, median), 
                         "M2_Nov" = apply(ts_2, 1, median),
                         "M3_Dec" = apply(ts_3, 1, median),
                         "M4_Jan" = apply(ts_4, 1, median),
                         "M5_Jun" = apply(ts_5, 1, median))

ts_medians <- data.frame("M1_Sep" = apply(ts_1, 1, median), 
                         "M2_Nov" = apply(ts_2, 1, median),
                         "M3_Dec" = apply(ts_3, 1, median))

# Extract patterns of proteins using k-means clustering method
k <- 10 # number of wanted clusters
set.seed(999)
res_km <- kmeans(ts_medians, k, nstart = 10)

# Prepare data for plotting
data_plot <- melt(data.table(ID = 1:nrow(ts_medians),
                             cluster = res_km$cluster,
                             ts_medians),
                  id.vars = c("ID", "cluster"),
                  variable.name = "Time",
                  variable.factor = FALSE
)

#data_plot[, Time := as.integer(gsub("V", "", Time))]

# Prepare centroids
centers <- melt(data.table(ID = 1:nrow(res_km$centers),
                           cluster = 1:nrow(res_km$centers),
                           res_km$centers),
                id.vars = c("ID", "cluster"),
                variable.name = "Time",
                variable.factor = FALSE
)

#centers[, Time := as.integer(gsub("V", "", Time))]

# Plot the results
ggplot(data_plot, aes(Time, value, group = ID)) +
  facet_wrap(~cluster, ncol = 5) +
  geom_line(color = "grey10", alpha = 0.4) +
  geom_line(data = centers, aes(Time, value), color = "firebrick1", alpha = 0.80, 
            size = 1.2) +
  labs(x = "Months", y = "Scaled intensities") +
  scale_x_discrete(labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
#  scale_x_discrete(labels = c('Sep', 'Nov', 'Dec')) +
  theme_bw()

ggsave(file.path(dir, paste0('kmeans', k,'_', species, '.png')), 
       scale = 0.8, width = 13, height = 5)
# k = 12: width = 12, height = 6, ncol = 4
# k = 5: width = 13, height = 2.7, ncol = 5
# k = 10: width = 13, height = 5, ncol = 5
# k = 18, 20:  width = 13, height = 10, ncol = 5

#############################################################################
####### GO ENRICHMENT ANALYSIS ##############################################
#############################################################################

### Upload annotation of the protein groups to prepare data for GO annotation
dir_with_data <- paste0('labeglo2/MS_results/Field/', species, '/', 
                        species,'_refchannels_all/')
protein_annotation <- read.csv(file.path(dir_with_data, 
                               paste0('annot_proteinGroups_field_', tolower(species),
                                      '_withGOannotation.csv')), 
                               sep = '\t',
                               header = T)

### Select the cluster of interest for GO enrichment analysis
cluster_number <- 9
data_plot_cluster <- subset(data_plot, cluster == cluster_number)
protein_ids <- unique(data_plot_cluster$ID) # to get ids of the protein from the cluster
pr_groups <- rownames(ts_medians[protein_ids,]) # to get their pr_group names

# annotate all proteins
protein_names_all <- protein_annotation[match(row.names(ts_medians), 
                              protein_annotation$protein_group),]$go_annotation
ts_medians_annot <- cbind(ts_medians, protein_names_all)
ts_medians_annot$id <- seq(1:length(ts_medians_annot$M1_Sep)) # assign ids to all proteins
# remove proteins with no annotation (suitable for the go enrichment analysis)
ts_medians_annot_sub <- subset(ts_medians_annot, protein_names_all != '*')

# get the true-false vector with true for the proteins from the chosen cluster
#true_false_vector <- match(protein_ids, ts_medians_annot_sub$id)
true_false_vector <- ts_medians_annot_sub$id %in% protein_ids  

################### Calculate distances to the centroids ######################

curr_center <- centers[centers$cluster == cluster_number,]
data_plot$cluster_value <- curr_center[match(data_plot$Time, 
                                             curr_center$Time),]$value
data_plot$diff <- with(data_plot, value - cluster_value)
sq_diffs <- aggregate(diff ~ ID, data_plot, function(x) sqrt(sum(x ^ 2)))
sq_diffs <- sq_diffs[sq_diffs$ID %in% ts_medians_annot_sub$id,]
sq_diffs$in_cluster <- true_false_vector

ggplot(sq_diffs) +
  geom_histogram(aes(diff, fill=in_cluster), binwidth=.1) +
  theme_bw()

###############################################################################

proteins <- ts_medians_annot_sub$protein_names_all
proteins <- toupper(proteins)

# prepare list with true-false vector and names of the proteins for it 
proteins_list <- true_false_vector
names(proteins_list) <- proteins
proteins_list <- factor(proteins_list)

# OR to take distance to the centroids of the cluster as a measure:
proteins_list <- sq_diffs$diff
names(proteins_list) <- proteins
sq_diffs$protein_names <- proteins
sq_diffs$prot_groups <- row.names(ts_medians_annot_sub)

# upload go annotation file
gene_to_go <- read_delim('~/labeglo2/MS_results/annotation/go_gaf/combined.gaf',
                      '\t', comment = '!', col_names = F, col_types = cols(.default = "c"))
# create named list with genes and corresponded to them go terms:
gene_go <- aggregate(X5 ~ X3, gene_to_go, function(x) list(unique(x)))
gene_go$X4 <- toupper(gene_go$X3)
gene_go_ <- gene_go$X5
names(gene_go_) <- gene_go$X4

proteins_list2 <- proteins_list[names(proteins_list) %in% gene_go$X4]
length(proteins_list2)

# run topgo
t <- 0.7
ontology <- 'BP'
GOdata <- new("topGOdata", description = "Simple session", 
              ontology = ontology,
              allGenes = proteins_list2, 
              #geneSel = topDiffGenes,
              geneSel = function(x) x <= t,
              nodeSize = 8,
              annot = annFUN.gene2GO, gene2GO = gene_go_)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, 
                   elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", 
                   topNodes = 50,
                   numChar = 100)

allRes_ordered <- allRes[order(allRes$Annotated, decreasing = T),]
Res_sign <- subset(allRes_ordered, classicKS < 0.05 | elimKS < 0.05)

#### Extract gene names of significant observations from go analysis results
allGO <- genesInTerm(GOdata)
#allGO['GO:0006090']

go_table_w_genes <- NULL

for (row in 1:nrow(allRes)) {
 current_row <- allRes[row,]
 current_go <- as.character(current_row[1])
 wanted_ <- proteins_list2[names(proteins_list2) %in% allGO[[current_go]]] 
 genes_names <- paste0(names(wanted_[wanted_ < 0.7]), collapse = '', sep = ';')
 names(genes_names) <- 'significant_genes'
 go_table_w_genes <- rbind(go_table_w_genes, 
                           as.data.frame(c(current_row, genes_names)))
}

write.table(go_table_w_genes, 
            file = file.path(dir, 
                             paste0('go_enrich_cluster', cluster_number, 
                                    '_outof', k, '_thres', t, 
                                    '_ontology', ontology, 
                                    '.csv')), row.names = F)

### Plot significant proteins from chosen GO terms (boxplots) 
data_time_series_ <- data.frame(data_time_series)
colnames(data_time_series_) <- colnames(data_time_series)

intensity_long <- data_time_series_ %>% 
  rownames_to_column(var = 'protein') %>%
  pivot_longer(-protein, values_to = 'intensity', names_to = 'sample')

intensity_long_sub <- subset(intensity_long, protein %in% sq_diffs$prot_groups)
intensity_long_sub$go_annotation <- 
  sq_diffs[match(intensity_long_sub$protein, 
                 sq_diffs$prot_groups),]$protein_names

get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

intensity_long_sub$whole_label <- sprintf('%s|%s', intensity_long_sub$protein,
                                             intensity_long_sub$go_annotation)
intensity_long_sub$whole_label2 <- factor(intensity_long_sub$whole_label,
    levels=unique(
    intensity_long_sub$whole_label[order(intensity_long_sub$go_annotation)]))
intensity_long_sub$diff <- sq_diffs[match(intensity_long_sub$protein, 
                                          sq_diffs$prot_groups),]$diff
intensity_long_sub$annot_long <- 
  protein_annotation[match(intensity_long_sub$protein,
                           protein_annotation$protein_group),]$diamond_annot

GO2draw <- 'GO:0006888'
wanted <- proteins_list2[names(proteins_list2) %in% allGO[[GO2draw]]] 
proteints2plot <- wanted[wanted <= t]
intensity_long_sub_toplot <- subset(intensity_long_sub, diff %in% proteints2plot)
intensity_long_sub_toplot$month <- meta[match(intensity_long_sub_toplot$sample,
                                              meta$sample),]$condition

intensity_long_sub_toplot$time <- 
  ifelse(grepl('September', intensity_long_sub_toplot$month), 1,
  ifelse(grepl('November', intensity_long_sub_toplot$month), 2,
  ifelse(grepl('December', intensity_long_sub_toplot$month), 3,
  ifelse(grepl('January', intensity_long_sub_toplot$month), 4, 5))))

centers$time <- ifelse(centers$Time == 'M1_Sep', 1,
                ifelse(centers$Time == 'M2_Nov', 2,
                ifelse(centers$Time == 'M3_Dec', 3,
                ifelse(centers$Time == 'M4_Jan', 4, 5))))
centers_to_draw <- subset(centers, cluster == cluster_number)
intensity_long_sub_toplot <- subset(intensity_long_sub_toplot, intensity < 4)
ggplot(intensity_long_sub_toplot, aes(time, intensity, group = time)) +
  facet_wrap(~ whole_label2, ncol = 3,
             labeller = as_labeller(get_protein_label)
             #, scales = 'free_y'
             ) +
  geom_boxplot(outlier.color = NA, 
               #aes(color = time == 4, fill = time == 4), alpha = .3
  ) +
  geom_jitter(
    #aes(color = time == 4), 
    size = 0.3) +
  geom_line(data = centers_to_draw, aes(time, value, group = 1), 
            color = "firebrick1", alpha = 0.80, 
            size = 1.2) +
  #  scale_color_manual(values = c('black', 'dodgerblue4')) +
  #  scale_fill_manual(values = c(NA, 'dodgerblue4')) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  guides(fill = F, color = F) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw()

ggsave(file.path(dir, paste0('boxplots_timeseries_cluster', cluster_number,
                            'outof', k, '_', GO2draw, 
                            '_ontology', ontology,'.png')), scale = 0.8,
       width = 6, height = 4.5)

intensity_long_sub_toplot <- subset(intensity_long_sub, 
                                    grepl('COP', go_annotation, ignore.case = T))
intensity_long_sub_toplot <- subset(intensity_long_sub, 
                                    grepl('coatomer', annot_long, ignore.case = T))

ggsave(file.path(dir, paste0('boxplots_', 
                              'COP', '.png')), 
       scale = 0.8,
       width = 9.5, height = 3)
