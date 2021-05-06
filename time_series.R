### Performing kmeans clustering and building plots based on the obtained clusters ###

library(TSrepr)
library(ggplot2)
library(data.table)
library(dplyr)

# Prepare data 
intensities <- data_wo_na # from de_separately_field.R
scaled_intensities <- apply(intensities, 1, scale)
scaled_intensities <- t(scaled_intensities)
colnames(scaled_intensities) <- meta$sample
data_time_series <- scaled_intensities[,!grepl('pool', colnames(scaled_intensities))]

# Obtain medians for every month
ts_1 <- data_time_series[,subset(meta, condition == 'September_BK')$sample]
ts_2 <- data_time_series[,subset(meta, condition == 'November_BK')$sample]
ts_3 <- data_time_series[,subset(meta, condition == 'December_BK')$sample]
ts_4 <- data_time_series[,subset(meta, condition == 'January_BK')$sample]
ts_5 <- data_time_series[,subset(meta, condition == 'June_BK')$sample]
ts_medians <- data.frame("M1_Sep" = apply(ts_1, 1, median), 
                         "M2_Nov" = apply(ts_2, 1, median),
                         "M3_Dec" = apply(ts_3, 1, median),
                         "M4_Jan" = apply(ts_4, 1, median),
                         "M5_Jun" = apply(ts_5, 1, median))

# Extract patterns of proteins using k-means clustering method
set.seed(999)
k <- 12
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
  facet_wrap(~cluster, ncol = 4) +
  geom_line(color = "grey10", alpha = 0.4) +
  geom_line(data = centers, aes(Time, value), color = "firebrick1", alpha = 0.80, 
            size = 1.2) +
  labs(x = "Months", y = "Scaled intensities") +
  scale_x_discrete(labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  theme_bw()

ggsave(file.path(dir, paste0('kmeans', k,'_eve.png')), width = 12, height = 7)
# k = 12: width = 12, height = 6, col = 4
# k = 5: width = 13, height = 2.7, col = 5
# k = 10: width = 13, height = 5, col = 5

### Upload annotation of the protein groups to prepare data for GO annotation

protein_annotation <- read.csv(file.path(dir, 'annot_proteinGroups_field_eve.csv'), 
                               sep = '\t',
                               header = T)






