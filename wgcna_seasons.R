library(WGCNA)
library(ggplot2)
library(pheatmap)
library(tempR)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tibble)
options(stringsAsFactors = FALSE)
brighten <- function(..., alpha) {
  inner <- function(col) {
    rgb(red=col[1] + (255 - col[1]) * alpha,
        green=col[2] + (255 - col[2]) * alpha,
        blue=col[3] + (255 - col[3]) * alpha,
        maxColorValue=255)
  }
  new_colors <- col2rgb(col=unlist(list(...)), alpha=F)
  apply(new_colors, 2, inner)
}

species <- 'Eve'
# 1. Set the paths and upload annotation file
dir_annot <- paste0('labeglo2/MS_results/Field/', species, 
                    '/', species,'_refchannels_all/')
dir_metafile <- paste0('labeglo2/MS_results/Field/', species, 
                       '/', species,'_refchannels_all/')
current_dir <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species,'_refchannels_all/')
pg_annot_complete <- read.csv(file.path(dir_annot, 
                              paste0('annot_proteinGroups_field_', tolower(species), 
                                     '_withGOannotation.csv')), 
                              sep = '\t', header = T) # eggnog and diamond annotation
metafile <- 'Metadata_Proteus.tsv'

# 2. Upload normalized intensities
data_for_wgcna <- read.csv(file.path(current_dir, 
                                     paste0('intensities_after_slNorm_woNA_', 
                                            tolower(species),'.csv')),
                           header = T, sep = '\t')
data_for_wgcna <- data_for_wgcna[, !grepl('pool', colnames(data_for_wgcna),
                                          ignore.case = T)]

# 3. Data preprocessing - get rid of missing data 
data_transp <- as.data.frame(t(data_for_wgcna))
gsp <- goodSamplesGenes(data_transp, verbose = 1)
gsp$allOK

if (!gsp$allOK) {
  dat_wo_missing <- data_transp[gsp$goodSamples, gsp$goodGenes]
} else {
  dat_wo_missing <- data_transp
}

sampleTree <- hclust(dist(dat_wo_missing), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# 4. Choosing the soft-thresholding power: analysis of network topology
soft_threshold <- function(data_wo_miss){
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(data_wo_miss, powerVector = powers, verbose = 5)
  #save(sft, file = 'results_soft_thresh_WGCNA.RData')
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 1.3
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"), cex.axis=1.5, cex.lab=1.3, cex.main=1.5)
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.90, col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"), cex.axis=1.5, cex.lab=1.3, cex.main=1.5)
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
}
soft_threshold(dat_wo_missing)

# Choose the power, which is the lowest power, for which the scale-free topology fit
# From the wgcna manual: index curve flattens out upon reaching a high value 

# 5. Upload metadata
meta_complete <- read.delim(file.path(dir_metafile, metafile), header=TRUE, sep="\t")
meta <- meta_complete
population <- 'BK'
meta <- subset(meta, condition != 'pool')
meta[which(meta$sex == ''),]$sex <- NA
meta <- meta[grepl(population, meta$sample),]
meta <- subset(meta, sample != 'F/100/BK/3c/2m_2')

colnames(data_for_wgcna) <- meta$sample

# 6. Prepare table with traits
metadata <- meta
row.names(metadata) <- metadata$sample
metadata$temperature <- ifelse(grepl('September', metadata$condition), 10.5,
                        ifelse(grepl('November', metadata$condition), 6.4,
                        ifelse(grepl('December', metadata$condition), 2.3,
                        ifelse(grepl('January', metadata$condition), 0.3, 7))))

metadata$temperature <- ifelse(grepl('September', metadata$condition), 10.2,
                        ifelse(grepl('November', metadata$condition), 2.1,
                        ifelse(grepl('December', metadata$condition), 0.3,
                        ifelse(grepl('January', metadata$condition), 0.7, 8)))) # Gla

sep_t_20days <- mean(c(9.6, 9.7, rep(10.5, 3), 10.4, 10.1, 10.8, 11.1, 10.7, 10.7, 10.1,
                       12.1, 12.4, 12.5, 12.3, 12.5, 12.7, 12.8, 13.4))
sep_t_20days_sd <- sd(c(9.6, 9.7, rep(10.5, 3), 10.4, 10.1, 10.8, 11.1, 10.7, 10.7, 10.1,
                        12.1, 12.4, 12.5, 12.3, 12.5, 12.7, 12.8, 13.4))
nov_t_20days <- mean(c(4.7, 5.2, 5.2, 5, 5.8, 5.7, 5.7, 5.5, 5.2, 5.9,
                       5.7, 5.5, 5.5, 5.6, 5.8, 6.7, 6.7, 6.2, 6.9, 6.9))
nov_t_20days_sd <- sd(c(4.7, 5.2, 5.2, 5, 5.8, 5.7, 5.7, 5.5, 5.2, 5.9,
                       5.7, 5.5, 5.5, 5.6, 5.8, 6.7, 6.7, 6.2, 6.9, 6.9))
dec_t_20days <- mean(c(rep(2.8, 3), 3.2, 3.2, 3.3, 3.1, rep(3.4, 3), 
                       3.5, 3.3, 3.1, 3.4, 3.5, 2.8, 3.4, 3.3, 3.1, 3.2))
dec_t_20days_sd <- sd(c(rep(2.8, 3), 3.2, 3.2, 3.3, 3.1, rep(3.4, 3), 
                       3.5, 3.3, 3.1, 3.4, 3.5, 2.8, 3.4, 3.3, 3.1, 3.2))
jan_t_20days <- mean(c(1.3, 1.4, 1.4, 1.5, 1.5, 1.6, 1.6, rep(1.7, 4),
                       rep(1.8, 4), 1.9, 2, 2, 2.1, 2.1))
jan_t_20days_sd <- sd(c(1.3, 1.4, 1.4, 1.5, 1.5, 1.6, 1.6, rep(1.7, 4),
                          rep(1.8, 4), 1.9, 2, 2, 2.1, 2.1))
jun_t_20days <- mean(c(3.3, 3.5, 3.7, 3.8, 3.8, rep(3.9, 4), 3.8,
                       4, 3.9, 4.1, 3.9, 3.6, 3.7, 3.7, 3, 3.4, 3))
jun_t_20days_sd <- sd(c(3.3, 3.5, 3.7, 3.8, 3.8, rep(3.9, 4), 3.8,
                       4, 3.9, 4.1, 3.9, 3.6, 3.7, 3.7, 3, 3.4, 3))

metadata$temp_20days <- ifelse(grepl('September', metadata$condition), sep_t_20days,
                        ifelse(grepl('November', metadata$condition), nov_t_20days,
                        ifelse(grepl('December', metadata$condition), dec_t_20days,
                        ifelse(grepl('January', metadata$condition), jan_t_20days, 
                                                                    jun_t_20days))))

metadata$sd_20days <- ifelse(grepl('September', metadata$condition), sep_t_20days_sd,
                        ifelse(grepl('November', metadata$condition), nov_t_20days_sd,
                        ifelse(grepl('December', metadata$condition), dec_t_20days_sd,
                        ifelse(grepl('January', metadata$condition), jan_t_20days_sd, 
                                                                     jun_t_20days_sd))))

metadata$cc_perc <- ifelse(grepl('September', metadata$condition), 100,
                    ifelse(grepl('November', metadata$condition), 0,
                    ifelse(grepl('December', metadata$condition), 100,
                    ifelse(grepl('January', metadata$condition), 0, 
                                                25))))

metadata$cc_points <- ifelse(grepl('September', metadata$condition), 9,
                           ifelse(grepl('November', metadata$condition), 6,
                                  ifelse(grepl('December', metadata$condition), 9,
                                         ifelse(grepl('January', metadata$condition), 9, 
                                                4))))

metadata$wind <- ifelse(grepl('September', metadata$condition), 4,
                           ifelse(grepl('November', metadata$condition), 1,
                                  ifelse(grepl('December', metadata$condition), 3,
                                         ifelse(grepl('January', metadata$condition), 0, 
                                                2))))
metadata$pa <- ifelse(grepl('September', metadata$condition), 1009.2,
                           ifelse(grepl('November', metadata$condition), 1026.9,
                                  ifelse(grepl('December', metadata$condition), 1022.4,
                                         ifelse(grepl('January', metadata$condition), 1025.2, 
                                                1005.5))))

#metadata$ice <- ifelse(grepl('January', metadata$condition), 1, 0)
metadata$`fouling/amplexus` <- ifelse(grepl('September|November', 
                                            metadata$condition), 1, 0)

metadata$daylight <- ifelse(grepl('September', metadata$condition), 11.78,
                        ifelse(grepl('November', metadata$condition), 9.37,
                        ifelse(grepl('December', metadata$condition), 7.92,
                        ifelse(grepl('January', metadata$condition), 8.73, 16.6))))

metadata$blooming <- ifelse(grepl('September', metadata$condition), 1.5,
                            ifelse(grepl('June', metadata$condition), 0.3, 0.2))

#metadata$amplexus <- ifelse(grepl('September|November', metadata$condition), 1, 0)
metadata$time <- ifelse(grepl('September', metadata$condition), 1,
                 ifelse(grepl('November', metadata$condition), 2,
                 ifelse(grepl('December', metadata$condition), 3,
                 ifelse(grepl('January', metadata$condition), 4, 5))))
metadata$Sep <- ifelse(grepl('September', metadata$condition), 1, 0)
metadata$Nov <- ifelse(grepl('November', metadata$condition), 1, 0)
metadata$Dec <- ifelse(grepl('December', metadata$condition), 1, 0)
metadata$Jan <- ifelse(grepl('January', metadata$condition), 1, 0)
metadata$Jun <- ifelse(grepl('June', metadata$condition), 1, 0)

#metadata$male <- ifelse(metadata$sex == 'male', 1, 0)
#metadata[is.na(metadata$sex),]$male <- 0
#metadata$female <- ifelse(metadata$sex == 'female', 1, 0)
#metadata[is.na(metadata$sex),]$female <- 0
metadata <- metadata[-c(1, 2, 3, 4, 5, 6)] # Eve, Gla
metadata <- metadata[-c(1, 2, 3, 4, 5)]

# 7. Constructing the gene network and identifying modules: 
sth_power <- 4 # soft-threshold power
sth_power <- 7
net <- blockwiseModules(dat_wo_missing, power = sth_power,
                        TOMType = "signed", minModuleSize = 25,
                        reassignThreshold = 0, mergeCutHeight = 0.15,
                        maxBlockSize = 12000,
                        deepSplit = 4, 
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        #saveTOMFileBase = "EveTOM_peptides",
                        verbose = 3)

table(net$colors) # number of modules

# 8. Plot the dendrogram
dendrogram <- function(net_to_draw){
  sizeGrWindow(12, 9)
  mergedColors <- labels2colors(net_to_draw$colors)
  plotDendroAndColors(net_to_draw$dendrograms[[1]], 
                      mergedColors[net_to_draw$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}
dendrogram(net)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 9. ##### Relating modules to external trait
### Quantifying module–trait associations 
nGenes <- ncol(dat_wo_missing)
nSamples <- nrow(dat_wo_missing)
MEs0 <- moduleEigengenes(dat_wo_missing, moduleColors)$eigengenes # Recalculate MEs with color labels
MEs <- orderMEs(MEs0)

# 10. ### Correlate eigengenes with external traits and look for 
# the most significant associations
dim(MEs)
moduleTraitCor <- cor(MEs, metadata, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

correlation_heatmap <- function(dir=current_dir, plot_width=900, plot_height=900,
                                ModuleTraitCorMatrix=moduleTraitCor,
                                ModuleTraitPvalueMatrix=moduleTraitPvalue,
                                meta=metadata, ME=MEs, power, MCH, MMS, ds, species, nameOfexp){
  sizeGrWindow(10,10)
  png(file=file.path(dir,
                     paste0('cor_heatmap_power', power,'_MCH', MCH, '_MMS', MMS,'_DS', 
                            ds,'_', species, nameOfexp, '_Proteins_IRSnorm.png')),
      width=plot_width, height=plot_height, res=200, pointsize=10)
  textMatrix <- paste(signif(ModuleTraitCorMatrix, 2), " (",
                      signif(ModuleTraitPvalueMatrix, 1), ")", sep = "");
  dim(textMatrix) <- dim(ModuleTraitCorMatrix)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = ModuleTraitCorMatrix,
                 xLabels = names(meta),
                 yLabels = names(ME),
                 ySymbols = names(ME),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
}

correlation_heatmap(power=sth_power, MCH=0.15, MMS=25, ds=4, species = species,
                    plot_width=1150, plot_height=900,
                    nameOfexp = 'Seasons_short')

### Module Membership (MM) and Gene Significance (GS) calculation
trait <- 'Jun'
trait_data <- as.data.frame(metadata[trait]) # take the trait of interest
names(trait_data) <- trait
modNames <- substring(names(MEs), 3) # take module names without 'ME'-preposition
# calculate correlations between the protein intensities and values of module eigengenes
geneModuleMembership <- as.data.frame(cor(dat_wo_missing, MEs, use = "p"))
# calculate asymptotic p-value for the given correlations of protein intensities and
# values of module eigengenes (after cor function above)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) 
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

# calculate correlation between the protein intensities and the trait
geneTraitSignificance <- as.data.frame(cor(dat_wo_missing, trait_data, use = "p"))
# calculate asymptotic p-values for the given correlation of protein intensities and traits
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(trait_data), sep="")
names(GSPvalue) <-paste("p.GS.", names(trait_data), sep="")

### Plot the genes from one module

module <- "turquoise"
column <- match(module, modNames)
moduleGenes <- moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for 1.5C",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, pch = 21,
                   bg = module)

### Summary output of network analysis results

proteins <- names(dat_wo_missing)
samples2annot <- match(proteins, pg_annot_complete$protein_group)
sum(is.na(samples2annot))

proteinInfo <- data.frame(protein = proteins,
                          geneSymbol = pg_annot_complete$annotation[samples2annot],
                          moduleColor = moduleColors,
                          geneTraitSignificance,
                          GSPvalue)

modOrder <- order(-abs(cor(MEs, trait_data, use = "p"))) # order modules by their significance for temperature

for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(proteinInfo)
  proteinInfo = data.frame(proteinInfo, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
  names(proteinInfo) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

proteinOrder <- order(proteinInfo$moduleColor, 
                      -abs(proteinInfo[sprintf('GS.%s', trait)]))
proteinInfo_final <- proteinInfo[proteinOrder,]

### Draw all hub genes on the plot

GS_treshold <- 0.6
MM_treshold <- 0.6
proteinInfo_final_hub <- proteinInfo_final[abs(proteinInfo_final[sprintf('GS.%s', 
                                                                         trait)]) 
                                           > GS_treshold,]

hubsMM <- NULL
for (protein in 1:nrow(proteinInfo_final_hub)){
  temp <- proteinInfo_final_hub[sprintf(c('MM.%s', 'p.MM.%s'), 
                  as.character(proteinInfo_final_hub[protein,][3]))][protein,]
  colnames(temp) <- c('MM', 'pMM')
  hubsMM <- rbind(hubsMM, temp)
}
proteinInfo_hub_reduced <- data.frame(proteinInfo_final_hub[c(1:5)], hubsMM)
proteinInfo_hub_reduced <- subset(proteinInfo_hub_reduced, abs(MM) > MM_treshold)

rgb.in <- col2rgb(proteinInfo_hub_reduced$moduleColor)
hsv_colors <- rgb2hsv(rgb.in)
#hsv_colors[2:3, ] <- 0.5
hub_colors_bright <- apply(hsv_colors, 2, function(row) hsv(row[1], row[2], row[3]))

proteinInfo_hub_reduced$geneSymbol <- 
  sub('PREDICTED: |LOW QUALITY PROTEIN: |-like|, partial', '', 
  proteinInfo_hub_reduced$geneSymbol)

proteinInfo_hub_reduced$geneSymbolWrap <- sapply(
  proteinInfo_hub_reduced$geneSymbol, function(x)
    paste(strwrap(x, width=20), collapse='\n'))

ggplot(proteinInfo_hub_reduced,
       aes_string(x = 'abs(MM)', y = sprintf('GS.%s', trait))) +
  geom_point(fill = hub_colors_bright,
             size = 3, shape = 21, stroke = 0.2) +
  xlab('Module Membership (abs. value)') +
  ylab('Gene Trait Significance') +
  geom_text_repel(aes(label = geneSymbolWrap), 
                  #ifelse(nchar(as.character(geneSymbol)) < 10, 
                  #as.character(geneSymbol), '')), 
                  min.segment.length = 0,
                  seed = 555, box.padding = 0.5, point.padding = 0.4, size = 3,
                  segment.alpha = 0.5, max.overlaps = 33) +
  theme_light() + 
  scale_y_continuous(breaks = seq(-1.00, 1.10, 0.25), limits = c(-1,1.1)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(), 
        panel.grid = element_line(color = 'white'), 
        panel.border = element_rect(color = 'white'))

ggsave(file.path(current_dir, paste0('hub_proteins_', trait, '_GS',
                                     GS_treshold, '.png')), scale = 1.7)

##############################################################
intensity_scaled <- apply(data_for_wgcna, 1, scale)
intensity_scaled <- t(intensity_scaled)
intensity_scaled <- data.frame(intensity_scaled)
colnames(intensity_scaled) <- rownames(metadata)

intensity_long <- intensity_scaled %>% 
  rownames_to_column(var = 'protein') %>%
  pivot_longer(-protein, values_to = 'intensity', names_to = 'sample')

pr_id <- match(intensity_long$protein, proteinInfo_final$protein)
intensity_long$module <- proteinInfo_final[pr_id,]$moduleColor
intensity_long$annotation <- proteinInfo_final[pr_id,]$geneSymbol
intensity_long$GS <- proteinInfo_final[pr_id,]$GS.temp_20days
intensity_long$GS <- proteinInfo_final[pr_id,]$GS.daylight
intensity_long$GS <- proteinInfo_final[pr_id,]$GS.time
intensity_long$GS <- proteinInfo_final[pr_id,]$GS.Jun

#intensity_long <- subset(intensity_long, abs(GS)>0.5)

intensity_long$temp_20days <- metadata[match(intensity_long$sample, 
                                          rownames(metadata)),]$temp_20days
intensity_long$daylight <- metadata[match(intensity_long$sample, 
                                             rownames(metadata)),]$daylight
intensity_long$time <- metadata[match(intensity_long$sample, 
                                          rownames(metadata)),]$time
#intensity_long$daylight <- factor(intensity_long$daylight,
#                                  levels = c('11.78', '9.37', '7.92', '8.73', '16.6'))
intensity_long_sub <- subset(intensity_long, module == 'yellow' | module == 'green' |
                               module == 'brown' | module == 'turquoise' |
                               module == 'purple' | module == 'tan' |
                               module == 'blue')
## OR
intensity_long_sub <- subset(intensity_long, module != 'grey')

## OR 
intensity_long_sub <- intensity_long

ggplot(intensity_long_sub, aes(temp_20days, intensity, group = protein)) +
  geom_line() +
  theme_light()

ggplot(intensity_long, aes(temp_20days, intensity, group = protein)) +
  geom_line(color = 'grey80', alpha = 0.6) +
  geom_line(data = intensity_long_sub, aes(temp_20days, intensity), 
            color = 'forestgreen') +
  theme_light()

ggplot(intensity_long_sub, aes(temp_20days, intensity)) +
  geom_jitter()

intensity_long_sub$condition <- 
  meta[match(intensity_long_sub$sample, meta$sample),]$condition

intensity_long_sub_med <- aggregate(intensity ~ protein + condition, intensity_long_sub, 
                                    median)
intensity_long_sub_med$daylight <- 
  ifelse(grepl('September', intensity_long_sub_med$condition), 11.78,
      ifelse(grepl('November', intensity_long_sub_med$condition), 9.37,
          ifelse(grepl('December', intensity_long_sub_med$condition), 7.92,
              ifelse(grepl('January', intensity_long_sub_med$condition), 8.73, 
                        16.6))))

intensity_long_sub_med$temperature <- 
  ifelse(grepl('September', intensity_long_sub_med$condition), 10.5,
    ifelse(grepl('November', intensity_long_sub_med$condition), 6.4,
      ifelse(grepl('December', intensity_long_sub_med$condition), 2.3,
         ifelse(grepl('January', intensity_long_sub_med$condition), 0.3, 7))))

intensity_long_sub_med$temp_20days <- 
  ifelse(grepl('September', intensity_long_sub_med$condition), sep_t_20days,
                        ifelse(grepl('November', intensity_long_sub_med$condition), nov_t_20days,
                        ifelse(grepl('December', intensity_long_sub_med$condition), dec_t_20days,
                        ifelse(grepl('January', intensity_long_sub_med$condition), jan_t_20days, 
                                                    jun_t_20days))))

intensity_long_sub_med$time <- 
  ifelse(grepl('September', intensity_long_sub_med$condition), 'Sep',
         ifelse(grepl('November', intensity_long_sub_med$condition), 'Nov',
                ifelse(grepl('December', intensity_long_sub_med$condition), 'Dec',
                       ifelse(grepl('January', intensity_long_sub_med$condition), 'Jan', 
                              'Jun'))))
intensity_long_sub_med$time <- factor(intensity_long_sub_med$time,
                                      levels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun'))

intensity_long_sub_med$annotation <- 
  pg_annot_complete[match(intensity_long_sub_med$protein, 
                          pg_annot_complete$protein_group),]$annotation
intensity_long_sub_med$go_annotation <- 
  pg_annot_complete[match(intensity_long_sub_med$protein, 
                          pg_annot_complete$protein_group),]$go_annotation
intensity_long_sub_med$module <- proteinInfo_final[match(intensity_long_sub_med$protein,
                                       proteinInfo_final$protein),]$moduleColor
intensity_long_sub_med$GS <- proteinInfo_final[match(intensity_long_sub_med$protein,
                             proteinInfo_final$protein),]$GS.temp_20days
intensity_long_sub_med$GS <- proteinInfo_final[match(intensity_long_sub_med$protein,
                             proteinInfo_final$protein),]$GS.daylight
intensity_long_sub_med$GS <- proteinInfo_final[match(intensity_long_sub_med$protein,
                             proteinInfo_final$protein),]$GS.time
intensity_long_sub_med$GS <- proteinInfo_final[match(intensity_long_sub_med$protein,
                                                     proteinInfo_final$protein),]$GS.Jun

little <- data.frame(temp = c(11.27, 5.77, 3.2, 1.72, 3.69),
                     month = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun'))
ggplot(intensity_long_sub_med, aes(temp_20days, intensity, group = protein)) +
  geom_line(aes(color = module)) +
  scale_color_manual('Module', 
                     values = sort(unique(intensity_long_sub_med$module))
                       #c('black', 'blue', 'brown', 'green', 'magenta', 
                      #          'deeppink', 'orange'), 
                     #labels=c('Pink', 'Yellow')
                     ) +
  geom_vline(xintercept = little$temp, color = 'slategray4') +
  xlab('Mean water temperature for 20 days before sampling') +
  ylab('Scaled intensities') +
  geom_text(data = little, aes(x = temp - 0.2, y = -1.5, label = month), angle = 90) +
  theme_light()

ggsave(file.path(current_dir, 'mean_temp_20days_pinkandyellowModules.png'))

## To plot 
little <- data.frame(temp = c(sep_t_20days, nov_t_20days, dec_t_20days, 
                              jan_t_20days, jun_t_20days),
                     month = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) # 20 days temp
little <- data.frame(daylight = c(11.78, 9.37, 7.92, 8.73, 16.6),
                     month = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) # daylight

gs_treshold <- 0.7
inten_sub_sub <- subset(intensity_long_sub_med, abs(GS) > gs_treshold)
protein_names <- subset(inten_sub_sub, temp_20days == max(temp_20days))
protein_names <- subset(inten_sub_sub, daylight == max(daylight))
protein_names <- subset(inten_sub_sub, time == 'Jun')

protein_names$annotation <- 
  sub('PREDICTED: |LOW QUALITY PROTEIN: |LOW QUALITY PROTEIN: |-like|, partial|isoform X1', 
                                '',
                                protein_names$annotation)

protein_names$ix <- 0
protein_names$ix[order(protein_names$intensity)] <- 1:nrow(protein_names)
protein_names$ix01 <- (protein_names$ix - 1) / (nrow(protein_names) - 1)
min_intens <- min(intensity_long_sub_med$intensity)
max_intens <- max(intensity_long_sub_med$intensity)
protein_names$y <- min_intens + protein_names$ix01 * (max_intens - min_intens)

### FOR temp_20days:
ggplot(intensity_long_sub_med, aes(temp_20days, intensity, group = protein)) +
  geom_line(color = 'grey80', alpha = .7) +
  geom_line(data = inten_sub_sub, 
            aes(temp_20days, intensity, group = protein, color = module)) +
  geom_vline(xintercept = little$temp, color = 'slategray4') +
  xlab('Mean water temperature for 20 days before sampling') +
  ylab('Scaled intensities') +
  geom_text(data = little, 
            aes(x = temp - 0.1, y = -1.9, label = month), angle = 90) +
  geom_label(data = protein_names,
             aes(x = 12, y = y, 
                 label = go_annotation, 
                 #label = annotation,
                 fill=module, hjust = 0),
             size = 3) +
  geom_segment(data = protein_names, 
             aes(x = 11.3, xend = 11.95, y = intensity, yend = y, color = module)) +
  scale_fill_manual(values = brighten('dodgerblue3', 'chartreuse3', 'darkorchid1', 
                                      'turquoise2', 'gold1', alpha = .6)) +
  scale_color_manual('Module', 
                     values = c('dodgerblue3','chartreuse3', 'darkorchid1', 
                                'turquoise2', 'gold1')) +
  scale_x_continuous(limits = c(NA, 12.5)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'meanTemp20days_6Modules_goAnnot.png'), width = 13,
       height = 7)
ggsave(file.path(current_dir, 'meanTemp20days_6Modules2.png'), width = 13,
       height = 7)

### FOR daylight:
ggplot(intensity_long_sub_med, aes(daylight, intensity, group = protein)) +
  geom_line(color = 'grey80', alpha = .7) +
  geom_line(data = inten_sub_sub, 
            aes(daylight, intensity, group = protein, color = module)) +
  geom_vline(xintercept = little$daylight, color = 'slategray4') +
  xlab('Daylight length, h') +
  ylab('Scaled intensities') +
  geom_text(data = little, 
            aes(x = daylight - 0.1, y = -1.9, label = month), angle = 90) +
  geom_label(data = protein_names,
             aes(x = 17.5, y = y, 
                 label = go_annotation, 
                 #label = annotation,
                 fill=module, hjust = 0),
             size = 2) +
  geom_segment(data = protein_names, 
               aes(x = 16.7, xend = 17.4, y = intensity, yend = y, color = module)) +
  scale_fill_manual(values = brighten('sienna4', 'firebrick2',  
                                      'turquoise2', 'gold1', alpha = .6)) +
  scale_color_manual('Module', 
                     values = c('sienna4', 'firebrick2',  
                                'turquoise2', 'gold1')) +
  #scale_fill_manual(values = brighten('black', 'dodgerblue3', 'sienna4', 
  #                                    'chartreuse3', 'tan1', 
  #                                    'turquoise2', alpha = .6)) +
  #scale_color_manual('Module', 
  #                   values = c('black', 'dodgerblue3', 'sienna4',
  #                              'chartreuse3', 'tan1', 
  #                              'turquoise2')) +
  scale_x_continuous(limits = c(NA, 17.7)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'daylight_7Modules_goAnnot.png'), width = 13,
       height = 8)
ggsave(file.path(current_dir, 'daylight_7Modules.png'), width = 15,
       height = 8)

my_colors <- c('brown' = 'saddlebrown', 'red' = 'firebrick2', 
               'turquoise' = 'turquoise2', 'yellow' = 'goldenrod1',
               'grey' = 'gray26', 'blue' = 'dodgerblue3', 
               'pink' = 'hotpink', 'magenta' = 'magenta4', 'purple' = 'purple4',
               'tan' = 'tan')
### FOR time:
#intensity_long_sub_med$as.numeric(intensity_long_sub_med$time)
ggplot(intensity_long_sub_med, aes(as.numeric(time), intensity, group = protein)) +
  geom_line(color = 'grey80', alpha = .7) +
  geom_line(data = inten_sub_sub, 
            aes(as.numeric(time), intensity, group = protein, color = module), 
            alpha = .4) +
  geom_vline(xintercept = c(1, 2, 3, 4, 5), color = 'slategray4') +
  geom_point(data = inten_sub_sub, 
             aes(as.numeric(time), intensity, group = protein, color = module)) +
  xlab('Months') +
  ylab('Scaled intensities') +
  geom_label(data = protein_names,
             aes(x = 5.4, y = y, 
                 label = go_annotation, 
                 #label = annotation,
                 fill=module, hjust = 0),
             size = 3) +
  geom_segment(data = protein_names, 
               aes(x = 5.05, xend = 5.35, y = intensity, yend = y, color = module)) +
  scale_fill_manual(values = brighten(my_colors, alpha = .6)) +
  scale_color_manual('Module', 
                     values = c(my_colors)) +
  scale_x_continuous(breaks = 1:length(levels(intensity_long_sub_med$time)),
                     labels = levels(intensity_long_sub_med$time),
                     minor_breaks = NULL, 
                     limits = c(1, 5.7)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))

ggsave(file.path(current_dir, 
                 paste0('Trait', trait, gs_treshold, 
                        'time_AllModules.png')),
       width = 14,
       height = 7,
       scale = 0.9)
ggsave(file.path(current_dir, 'time_7Modules.png'), width = 18,
       height = 8)

### FOR time, only one module (if take absGS > 0.5)
module2draw <- 'red'
inten_sub_sub <- subset(inten_sub_sub, module == module2draw)
protein_names <- subset(inten_sub_sub, time == 'Jun')
protein_names$annotation <- sub('PREDICTED: |LOW QUALITY PROTEIN: |LOW QUALITY PROTEIN: |-like', 
                                '',
                                protein_names$annotation)
protein_names$ix <- 0
protein_names$ix[order(protein_names$intensity)] <- 1:nrow(protein_names)
protein_names$ix01 <- (protein_names$ix - 1) / (nrow(protein_names) - 1)
min_intens <- min(intensity_long_sub_med$intensity)
max_intens <- max(intensity_long_sub_med$intensity)
protein_names$y <- min_intens + protein_names$ix01 * (max_intens - min_intens)

ggplot(intensity_long_sub_med, aes(as.numeric(time), intensity, group = protein)) +
  geom_line(color = 'grey80', alpha = .7) +
  geom_line(data = inten_sub_sub, 
            aes(as.numeric(time), intensity, group = protein, color = module)) +
  geom_vline(xintercept = c(1, 2, 3, 4, 5), color = 'slategray4') +
  xlab('Months') +
  ylab('Scaled intensities') +
  geom_label(data = protein_names,
             aes(x = 5.4, y = y, 
                 #label = go_annotation, 
                 label = annotation,
                 fill=module, hjust = 0),
             size = 3) +
  geom_segment(data = protein_names, 
               aes(x = 5.05, xend = 5.35, y = intensity, yend = y, color = module)) +
  scale_fill_manual(values = brighten(my_colors, 
                                      alpha = .6)) +
  scale_color_manual('Module', 
                     values = my_colors) +
  scale_x_continuous(breaks = 1:length(levels(intensity_long_sub_med$time)),
                     labels = levels(intensity_long_sub_med$time),
                     minor_breaks = NULL, 
                     limits = c(1, 5.7)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'time_7Modules_yellow_goAnnot_GS05.png'), width = 13,
       height = 8)

##### GO enrichment
library(topGO)
library(readr)
# upload go annotation file
gene_to_go <- read_delim('~/labeglo2/MS_results/annotation/go_gaf/combined.gaf',
                         '\t', comment = '!', col_names = F, 
                         col_types = cols(.default = "c"))
gene_to_go_bonus <- read_delim('~/labeglo2/MS_results/annotation/go_gaf/combined_bonus.csv',
                               '\t') # output from goterms_enrichment.py
# create named list with genes and corresponded go terms:
gene_go <- aggregate(X5 ~ X3, gene_to_go, function(x) list(unique(x)))
gene_go$X4 <- toupper(gene_go$X3)
gene_go_ <- gene_go$X5
names(gene_go_) <- gene_go$X4

# create named the same list but for the proteins with alternative names (bonus from goterms_enrichment.py) 
gene_go_bonus <- aggregate(go_term ~ gene_name, gene_to_go_bonus, 
                           function(x) list(unique(x)))
gene_go_bonus$up_gene_name <- toupper(gene_go_bonus$gene_name)
gene_go_bonus_ <- gene_go_bonus$go_term
names(gene_go_bonus_) <- gene_go_bonus$up_gene_name

#####
module <- 'magenta'
intensity_long$MM <- proteinInfo_final[pr_id, paste0('MM.', module)]

intensity_long$pMM <- proteinInfo_final[pr_id,]$p.MM.green
intensity_long$MM <- proteinInfo_final[pr_id,]$MM.green
intensity_long$pGS <- proteinInfo_final[pr_id,]$p.GS.daylight

intensity_long$go_annotation <- 
  pg_annot_complete[match(intensity_long$protein, 
                          pg_annot_complete$protein_group),]$go_annotation

# create the dataframe with proteins and information to them:
unique_proteins <- group_by(intensity_long, protein) %>% 
  slice_head() %>%
  ungroup()

# take wanted statistics 
proteins_list <- unique_proteins$MM

proteins_list <- unique_proteins$pMM
proteins_list <- unique_proteins$pGS
proteins_list <- unique_proteins$GS

# name vector with statistics by the corresponed protein:
names(proteins_list) <- toupper(unique_proteins$go_annotation)
proteins_list2 <- proteins_list[names(proteins_list) != '*'] # avoid not annotated proteins
proteins_list3 <- proteins_list2[names(proteins_list2) %in% gene_go$X4] # take only proteins with known go terms

# enrich the analysis by the observations with alternative names:
proteins_not_in_main_list <- proteins_list2[!names(proteins_list2) %in% gene_go$X4]
gene_go_from_bonus <- gene_go_bonus_[names(gene_go_bonus_) %in%
                                     names(proteins_not_in_main_list)]

# update go term list
gene_go_all_ <- c(gene_go_, gene_go_from_bonus)

# update the list with proteins in analysis
proteins_list4 <- proteins_list2[names(proteins_list2) %in% names(gene_go_all_)]

# run topgo
threshold <- 0.6
ontology <- 'BP'
GOdata <- new("topGOdata", description = "Simple session", 
              ontology = ontology,
              allGenes = proteins_list4, 
              #geneSel = function(x) x < 0.01, # pMM, pGS
              geneSel = function(x) abs(x) > threshold, # MM, GS
              nodeSize = 8,
              annot = annFUN.gene2GO, 
              gene2GO = gene_go_all_)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", 
                   topNodes = 70,
                   numChar = 100)

allGO <- genesInTerm(GOdata)
#allGO['GO:0006090']

go_table_w_genes <- NULL

for (row in 1:nrow(allRes)) {
  current_row <- allRes[row,]
  current_go <- as.character(current_row[1])
  wanted_ <- proteins_list4[names(proteins_list4) %in% allGO[[current_go]]] 
  genes_names <- paste0(names(wanted_[wanted_ > threshold]), collapse = '', sep = ';')
  names(genes_names) <- 'significant_genes'
  go_table_w_genes <- rbind(go_table_w_genes, 
                            as.data.frame(c(current_row, genes_names)))
}

write.table(go_table_w_genes, 
            file = file.path(current_dir, 
                             paste0('topgoResults_', module, 
                                    'Module_', threshold, 
                                    '_ontology', ontology,'.csv')), row.names = F)

### Plot significant proteins from chosen GO terms (boxplots) 

intensity_long_sub$go_annotation <- 
  intensity_long[match(intensity_long_sub$protein, 
                 intensity_long$protein),]$go_annotation
intensity_long_sub$annot_long <- 
  pg_annot_complete[match(intensity_long_sub$protein,
                          pg_annot_complete$protein_group),]$diamond_annot

get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

intensity_long_sub$whole_label <- sprintf('%s|%s', intensity_long_sub$protein,
                                          intensity_long_sub$go_annotation)
intensity_long_sub$whole_label2 <- factor(intensity_long_sub$whole_label,
                                          levels=unique(
                                            intensity_long_sub$whole_label[order(intensity_long_sub$go_annotation)]))
intensity_long_sub$MM <- intensity_long[match(intensity_long_sub$protein, 
                                          intensity_long$protein),]$MM

GO2draw <- 'GO:0006096'
wanted <- proteins_list4[names(proteins_list4) %in% allGO[[GO2draw]]] 
proteints2plot <- wanted[abs(wanted) > threshold]
intensity_long_sub_toplot <- subset(intensity_long_sub, MM %in% proteints2plot)

intensity_long_sub_toplot$month <- meta[match(intensity_long_sub_toplot$sample,
                                              meta$sample),]$condition

intensity_long_sub_toplot$month <- sub('_BK', '', intensity_long_sub_toplot$month)

intensity_long_sub_toplot$time <- ifelse(intensity_long_sub_toplot$month == 'September', 1,
                                         ifelse(intensity_long_sub_toplot$month == 'November' , 2,
                                                ifelse(intensity_long_sub_toplot$month == 'December', 3,
                                                       ifelse(intensity_long_sub_toplot$month == 'January', 4, 5))))

months_mean <- aggregate(intensity ~ condition, intensity_long_sub_toplot, mean)
months_mean$condition <- sub('_BK', '', months_mean$condition)
months_mean$time <- ifelse(months_mean$condition == 'September', 1,
                           ifelse(months_mean$condition == 'November', 2,
                                  ifelse(months_mean$condition == 'December', 3,
                                         ifelse(months_mean$condition == 'January', 4, 5))))

intensity_long_sub_toplot <- subset(intensity_long_sub_toplot, intensity < 4) 

intensity_long_sub_toplot$sex <- meta[match(intensity_long_sub_toplot$sample,
                                            meta$sample),]$sex

ggplot(intensity_long_sub_toplot, aes(time, intensity, group = time)) +
  facet_wrap(~ whole_label2, ncol = 5,
             labeller = as_labeller(get_protein_label)
             ) +
             #, scales = 'free_y') +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(size = 0.3
              #, aes(color = sex) # sex
              ) +
  geom_line(data = months_mean, aes(time, intensity, group = 1), 
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

ggsave(file.path(current_dir, paste0('boxplots_', module, 
                             'Module_MM', threshold,
                             '_ontology', ontology, '_',  GO2draw,'.png')), 
       scale = 0.8,
       width = 10, height = 6)

##### To look at certain proteins:
intensity_long_sub_toplot <- subset(intensity_long_sub, 
                                    grepl('myh|myl', annotation, ignore.case = T))
intensity_long_sub_toplot <- subset(intensity_long_sub, 
                                    grepl('hemocyanin', annot_long, ignore.case = T))
# then -> go to the steps above the plot

ggsave(file.path(current_dir, paste0('boxplots_', 
                                     'bt', '.png')), 
       scale = 0.8,
       width = 8.7, height = 2.5)

#### facet with boxplots with proteins with high GS for certain trait

intensity_long_gs_sub <- subset(intensity_long, GS > 0.5)
intensity_long_gs_sub <- 
  intensity_long_gs_sub[toupper(intensity_long_gs_sub$go_annotation) %in% names(gene_go_all_),]

get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

intensity_long_gs_sub$whole_label <- sprintf('%s|%s', intensity_long_gs_sub$protein,
                                             intensity_long_gs_sub$go_annotation)
intensity_long_gs_sub$whole_label2 <- factor(intensity_long_gs_sub$whole_label,
    levels=unique(
          intensity_long_gs_sub$whole_label[order(intensity_long_gs_sub$go_annotation)]))
lm_res <- lm(intensity ~ time, intensity_long_gs_sub)
lm_res <- lm(intensity ~ daylight, intensity_long_gs_sub)
# time
ggplot(intensity_long_gs_sub, aes(time, intensity, group = time)) +
  facet_wrap(~ whole_label2, ncol = 6,
             labeller = as_labeller(get_protein_label),
             scales = 'free_y') +
  geom_boxplot(outlier.color = NA, 
               #aes(color = time == 4, fill = time == 4), alpha = .3
               ) +
  geom_jitter(
    #aes(color = time == 4), 
    size = 0.3) +
  geom_abline(intercept = lm_res$coefficients[1], slope = lm_res$coefficients[2],
              color = 'firebrick') +
#  scale_color_manual(values = c('black', 'dodgerblue4')) +
#  scale_fill_manual(values = c(NA, 'dodgerblue4')) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  guides(fill = F, color = F) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw()

ggsave(file.path(current_dir, 'gsMinus50time_proteins.png'), scale = 2.1)

# daylight
little <- data.frame(daylight = c(11.78, 9.37, 7.92, 8.73, 16.6),
                     month = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) # daylight
intensity_long_gs_sub$month <- little[match(intensity_long_gs_sub$daylight, 
                                              little$daylight),]$month
intensity_long_gs_sub$month <- factor(intensity_long_gs_sub$month,
                                      levels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun'))
ggplot(intensity_long_gs_sub, aes(daylight, intensity, group = daylight)) +
  facet_wrap(~ whole_label2, ncol = 6,
             labeller = as_labeller(get_protein_label),
             scales = 'free_y') +
  geom_boxplot(outlier.color = NA, 
               aes(color = month)) +
  geom_jitter(aes(color = month),
    size = 0.2) +
  geom_abline(intercept = lm_res$coefficients[1], slope = lm_res$coefficients[2],
              color = 'firebrick') +
  scale_color_manual('Month',
                     values = c('#ff6361', '#bc5090', '#58508d', '#003f5c', '#ffa600')) +
  xlab('Daylight, h') +
  ylab('Scaled intensities') +
  theme_bw()

ggsave(file.path(current_dir, 'gs50daylight_proteins.png'), scale = 5)

### 
large_ribosome <- subset(intensity_long, 
                         grepl('Rpl', go_annotation, ignore.case = T)) # 60S
large_ribosome <- subset(intensity_long, 
                         grepl('^Rps|sup|sop', go_annotation, ignore.case = T)) # 40S
large_ribosome$month <- 
  ifelse(large_ribosome$daylight == 11.78, 'Sep',
         ifelse(large_ribosome$daylight == 9.37, 'Nov',
                ifelse(large_ribosome$daylight == 7.92, 'Dec',
                       ifelse(large_ribosome$daylight == 8.73, 'Jan', 
                              'Jun'))))
large_ribosome_med <- aggregate(intensity ~ protein + month, large_ribosome, 
                                median)
large_ribosome_med$daylight <- 
  ifelse(grepl('Sep', large_ribosome_med$month), 11.78,
         ifelse(grepl('Nov', large_ribosome_med$month), 9.37,
                ifelse(grepl('Dec', large_ribosome_med$month), 7.92,
                         ifelse(grepl('Jan', large_ribosome_med$month), 8.73, 
                              16.6))))
large_ribosome_med$annotation <- 
  pg_annot_complete[match(large_ribosome_med$protein, 
                          pg_annot_complete$protein_group),]$annotation
large_ribosome_med$go_annotation <- 
  pg_annot_complete[match(large_ribosome_med$protein, 
                          pg_annot_complete$protein_group),]$go_annotation
large_ribosome_med$module <- proteinInfo_final[match(large_ribosome_med$protein,
                          proteinInfo_final$protein),]$moduleColor

protein_names <- subset(large_ribosome_med, daylight == max(daylight))

protein_names$annotation <- sub('PREDICTED: |LOW QUALITY PROTEIN: |LOW QUALITY PROTEIN: |-like', 
                                '',
                                protein_names$annotation)

protein_names$ix <- 0
protein_names$ix[order(protein_names$intensity)] <- 1:nrow(protein_names)
protein_names$ix01 <- (protein_names$ix - 1) / (nrow(protein_names) - 1)
min_intens <- min(intensity_long_sub_med$intensity)
max_intens <- max(intensity_long_sub_med$intensity)
protein_names$y <- min_intens + protein_names$ix01 * (max_intens - min_intens)

ggplot(large_ribosome_med, aes(daylight, intensity, group = protein)) +
  geom_line(aes(color = module)) +
  geom_vline(xintercept = little$daylight, color = 'slategray4') +
  xlab('Daylight length, h') +
  ylab('Scaled intensities') +
  geom_text(data = little, 
            aes(x = daylight - 0.1, y = -1.9, label = month), angle = 90) +
  geom_label(data = protein_names,
             aes(x = 17.5, y = y, 
                 label = go_annotation, 
                 #label = annotation,
                 fill=module, hjust = 0),
             size = 2) +
  geom_segment(data = protein_names, 
               aes(x = 16.7, xend = 17.4, y = intensity, yend = y, color = module)) +
  scale_fill_manual(values = brighten('sienna4', 'chartreuse3', 
                                      'turquoise2', alpha = .6)) +
  scale_color_manual('Module', 
                     values = c('sienna4','chartreuse3',
                                'turquoise2')) +
  scale_x_continuous(limits = c(NA, 17.7)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

large_ribosome_med$time <- ifelse(large_ribosome_med$month == 'Sep', 1,
                           ifelse(large_ribosome_med$month == 'Nov', 2, 
                           ifelse(large_ribosome_med$month == 'Dec', 3,
                           ifelse(large_ribosome_med$month == 'Jan', 4, 5))))
large_ribosome_med$month <- factor(large_ribosome_med$month, levels = c('Sep', 'Nov', 'Dec', 'Jan',
                                                                  'Jun'))
ggplot(large_ribosome_med, aes(time, intensity, group = protein)) +
  geom_line(aes(color = module)) +
#  geom_vline(xintercept = c(1, 2, 3, 4, 5), color = 'slategray4') +
  xlab('Months') +
  ylab('Scaled intensities') +
  geom_label(data = protein_names,
             aes(x = 5.4, y = y, 
                 label = go_annotation, 
                 #label = annotation,
                 fill=module, hjust = 0),
             size = 3) +
  geom_segment(data = protein_names, 
               aes(x = 5.05, xend = 5.35, y = intensity, yend = y, color = module)) +
  scale_fill_manual(values = brighten('sienna4', 'chartreuse4','turquoise3', 
                                      'gold1', alpha = .6)) +
  scale_color_manual('Module', 
                     values = c('sienna4','chartreuse4','turquoise3', 
                                'gold1')) +
  scale_x_continuous(breaks = 1:length(levels(large_ribosome_med$month)),
                     labels = levels(large_ribosome_med$month),
                     minor_breaks = NULL, 
                     limits = c(1, 5.7)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'RPSs_months_allRpss.png'), scale = 2)

#### TRiC complex

tric <- subset(intensity_long, 
               grepl('cct|tcp|HSPD1', go_annotation, ignore.case = T))
tric$month <- metadata[match(tric$sample, rownames(metadata)),]$time
tric_med <- aggregate(intensity ~ protein + month, tric, 
                      median)

tric_med$annotation <- 
  pg_annot_complete[match(tric_med$protein, 
                          pg_annot_complete$protein_group),]$annotation
tric_med$go_annotation <- 
  pg_annot_complete[match(tric_med$protein, 
                          pg_annot_complete$protein_group),]$go_annotation
tric_med$module <- proteinInfo_final[match(tric_med$protein,
                                           proteinInfo_final$protein),]$moduleColor

tric_med$time <- ifelse(tric_med$month == 1, 'Sep',
                        ifelse(tric_med$month == 2, 'Nov', 
                               ifelse(tric_med$month == 3, 'Dec',
                                      ifelse(tric_med$month == 4, 'Jan', 'Jun'))))
tric_med$time <- factor(tric_med$time, levels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun'))
tric_med$mitoch <- ifelse(tric_med$go_annotation == 'HSPD1', 'Mitochondria', 'Cytoplasm')

protein_names <- subset(tric_med, month == max(month))

protein_names$annotation <- 
  sub('PREDICTED: |LOW QUALITY PROTEIN: |LOW QUALITY PROTEIN: |-like', '',
                                protein_names$annotation)

protein_names$ix <- 0
protein_names$ix[order(protein_names$intensity)] <- 1:nrow(protein_names)
protein_names$ix01 <- (protein_names$ix - 1) / (nrow(protein_names) - 1)
min_intens <- min(intensity_long_sub_med$intensity)
max_intens <- max(intensity_long_sub_med$intensity)
protein_names$y <- min_intens + protein_names$ix01 * (max_intens - min_intens)

ggplot(tric_med, aes(month, intensity, group = protein)) +
  geom_line(aes(color = mitoch)) +
  #  geom_vline(xintercept = c(1, 2, 3, 4, 5), color = 'slategray4') +
  xlab('Months') +
  ylab('Scaled intensities') +
  geom_label(data = protein_names,
             aes(x = 5.4, y = y, 
                 label = go_annotation, hjust = 0, 
                 fill = mitoch),
             size = 3) +
  geom_segment(data = protein_names, 
               aes(x = 5.05, xend = 5.35, y = intensity, yend = y, color = mitoch),
               size = 0.4) +
  scale_color_manual('Location', 
                     values = c('#58508d', '#bc5090')) +
  scale_fill_manual('Location', 
                    values = c(brighten('#58508d', '#bc5090', alpha = .7))) +
  scale_x_continuous(breaks = 1:length(levels(tric_med$time)),
                     labels = levels(tric_med$time),
                     minor_breaks = NULL, 
                     limits = c(1, 5.7)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'TriC_months.png'), scale = 1.6)

### hemocyanins
hemo <- subset(intensity_long, 
               grepl('hemocyanin', annotation, ignore.case = T))
hemo$month <- metadata[match(hemo$sample, rownames(metadata)),]$time
hemo_med <- aggregate(intensity ~ protein + month, hemo, 
                      median)

hemo_med$annotation <- 
  pg_annot_complete[match(hemo_med$protein, 
                          pg_annot_complete$protein_group),]$annotation
hemo_med$go_annotation <- 
  pg_annot_complete[match(hemo_med$protein, 
                          pg_annot_complete$protein_group),]$go_annotation
hemo_med$module <- proteinInfo_final[match(hemo_med$protein,
                                           proteinInfo_final$protein),]$moduleColor

hemo_med$time <- ifelse(hemo_med$month == 1, 'Sep',
                        ifelse(hemo_med$month == 2, 'Nov', 
                               ifelse(hemo_med$month == 3, 'Dec',
                                      ifelse(hemo_med$month == 4, 'Jan', 'Jun'))))
hemo_med$time <- factor(hemo_med$time, levels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun'))
hemo_med$temp20 <- metadata[match(hemo_med$month, metadata$time),]$temp_20days

protein_names <- subset(hemo_med, month == max(month)) # if x - month
protein_names <- subset(hemo_med, temp20 == max(temp20)) # if x - temp

protein_names$annotation <- sub('PREDICTED: |LOW QUALITY PROTEIN: |LOW QUALITY PROTEIN: |-like', 
                                '',
                                protein_names$annotation)

protein_names$ix <- 0
protein_names$ix[order(protein_names$intensity)] <- 1:nrow(protein_names)
protein_names$ix01 <- (protein_names$ix - 1) / (nrow(protein_names) - 1)
min_intens <- min(intensity_long_sub_med$intensity)
max_intens <- max(intensity_long_sub_med$intensity)
protein_names$y <- min_intens + protein_names$ix01 * (max_intens - min_intens)

# x - month
ggplot(hemo_med, aes(month, intensity, group = protein)) +
  geom_line(aes(color = module)) +
  #  geom_vline(xintercept = c(1, 2, 3, 4, 5), color = 'slategray4') +
  xlab('Months') +
  ylab('Scaled intensities') +
  geom_label(data = protein_names,
             aes(x = 5.4, y = y, 
                 label = annotation, hjust = 0, 
                 fill = module),
             size = 3) +
  geom_segment(data = protein_names, 
               aes(x = 5.05, xend = 5.35, y = intensity, yend = y, color = module),
               size = 0.4) +
  scale_color_manual('Module', 
                     #values = c('grey40', '#bc5090')) +
                     values = c('brown', 'green', 'grey40', '#bc5090')) +
  scale_fill_manual('Module', 
                    #values = c(brighten('grey40', '#bc5090', alpha = .7))) +
                    values = c(brighten('brown', 'green', 'grey40', '#bc5090', alpha = .7))) +
  scale_x_continuous(breaks = 1:length(levels(hemo_med$time)),
                     labels = levels(hemo_med$time),
                     minor_breaks = NULL, 
                     limits = c(1, 6.5)) +
  guides(fill = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'Hemocyanin_months.png'), scale = 1.7)

# x - temp
little <- data.frame(temp = c(sep_t_20days, nov_t_20days, dec_t_20days, 
                              jan_t_20days, jun_t_20days),
                     month = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) # 20 days temp
hemo_med$corr_w_temp <- 
  ifelse(grepl('NODE_10114_length_2944_cov_11444.185838_g4667_i1.p1_Eve;TRINITY_DN25900_c0_g1_i4.p1_Eve', 
               hemo_med$protein), T,
  ifelse(grepl('NODE_84781_length_880_cov_3856.855596_g66860_i0.p1_Gla;NODE_15570_length_2415_cov',
               hemo_med$protein), T, F))
ggplot(hemo_med, aes(temp20, intensity, group = protein)) +
  geom_line(aes(color = module, linetype = corr_w_temp)) +
  geom_vline(xintercept = little$temp, color = 'slategray4') +
  xlab('Mean water temperature for 20 days before sampling') +
  ylab('Scaled intensities') +
  geom_text(data = little, 
            aes(x = temp - 0.2, y = -1.9, label = month), angle = 90) +
  geom_label(data = protein_names,
             aes(x = 12, y = y, 
                 #label = go_annotation, 
                 label = annotation,
                 fill=module, hjust = 0),
             size = 3) +
  geom_segment(data = protein_names, 
               aes(x = 11.3, xend = 11.95, y = intensity, yend = y, color = module),
               size = 0.4) +
  scale_fill_manual(values = brighten('grey40', '#bc5090', alpha = .6)) +
  scale_color_manual('Module', 
                     values = c('grey40', '#bc5090')) +
  scale_linetype_manual(values = c('dotted', 'solid')) +
  scale_x_continuous(limits = c(NA, 14.5)) +
  guides(fill = F, linetype = F) +
  theme_light() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(file.path(current_dir, 'Hemocyanin_temp20_withCorrWithTemp20.png'), scale = 1.7)
  
