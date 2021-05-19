### Performing kmeans clustering and building plots based on the obtained clusters ###

library(TSrepr)
library(ggplot2)
library(data.table)
library(dplyr)

species <- "Eve"
population <- "BK"
dir <- paste0('labeglo2/MS_results/Field/', species, '/', population)
data_wo_na <- read.csv(file.path(dir, "intensities_after_slNorm_woNA_eve.csv"),
                        sep = '\t')
### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  return(meta)
}
path2meta <-
  'labeglo2/MS_results/Field/Eve/Metadata_Proteus.tsv'

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

# Obtain medians for every month
ts_1 <- data_time_series[,subset(meta, 
                                 condition == paste0('September_', population))$sample]
ts_2 <- data_time_series[,subset(meta, 
                                 condition == paste0('November_', population))$sample]
ts_3 <- data_time_series[,subset(meta, 
                                 condition == paste0('December_', population))$sample]
ts_4 <- data_time_series[,subset(meta, condition == 'January_BK')$sample]
ts_5 <- data_time_series[,subset(meta, condition == 'June_BK')$sample]
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

ggsave(file.path(dir, paste0('kmeans', k,'_eve.png')), width = 13, height = 5)
# k = 12: width = 12, height = 6, ncol = 4
# k = 5: width = 13, height = 2.7, ncol = 5
# k = 10: width = 13, height = 5, ncol = 5
# k = 18, 20:  width = 13, height = 10, ncol = 5

### Upload annotation of the protein groups to prepare data for GO annotation
dir_with_data <- 'labeglo2/MS_results/Field/Eve/'
protein_annotation <- read.csv(file.path(dir_with_data, 
                                         'annot_proteinGroups_field_eve_withGOannotation.csv'), 
                               sep = '\t',
                               header = T)

### Select the cluster of interest
cluster_number <- 9
data_plot_cluster <- subset(data_plot, cluster == cluster_number)
protein_ids <- unique(data_plot_cluster$ID)

pr_groups <- rownames(ts_medians[protein_ids,])

protein_names <- protein_annotation[protein_annotation$protein_group %in% 
                                      pr_groups,]$go_annotation
protein_names_ <- subset(protein_names, protein_names != '*')
true_false_vector <- protein_annotation$protein_group %in% pr_groups





