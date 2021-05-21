library(WGCNA)
library(ggplot2)
library(pheatmap)
library(tempR)
library(ggrepel)
library(dplyr)
options(stringsAsFactors = FALSE)

# 1. Set the paths and upload annotation file
dir_annot <- 'labeglo2/MS_results/Field/Eve/'
dir_metafile <- 'labeglo2/MS_results/Field/Eve/'
current_dir <- 'labeglo2/MS_results/Field/Eve/BK/'
pg_annot_complete <- read.csv(file.path(dir_annot, 
                              'annot_proteinGroups_field_eve_withGOannotation.csv'), 
                              sep = '\t', header = T) # eggnog and diamond annotation
metafile <- 'Metadata_Proteus.tsv'

# 2. Upload normalized intensities
data_for_wgcna <- read.csv(file.path(current_dir, 'intensities_after_slNorm_woNA_eve.csv'),
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

metadata$ice <- ifelse(grepl('January', metadata$condition), 1, 0)
metadata$fouling <- ifelse(grepl('September|November', metadata$condition), 1, 0)

metadata$daylight <- ifelse(grepl('September', metadata$condition), 11.78,
                        ifelse(grepl('November', metadata$condition), 9.37,
                        ifelse(grepl('December', metadata$condition), 7.92,
                        ifelse(grepl('January', metadata$condition), 8.73, 16.6))))

metadata$blooming <- ifelse(grepl('September', metadata$condition), 1.5,
                            ifelse(grepl('June', metadata$condition), 0.3, 0.2))

metadata$amplexus <- ifelse(grepl('September|November', metadata$condition), 1, 0)
metadata$time <- ifelse(grepl('September', metadata$condition), 1,
                 ifelse(grepl('November', metadata$condition), 2,
                 ifelse(grepl('December', metadata$condition), 3,
                 ifelse(grepl('January', metadata$condition), 4, 5))))
metadata$Sep <- ifelse(grepl('September', metadata$condition), 1, 0)
metadata$Nov <- ifelse(grepl('November', metadata$condition), 1, 0)
metadata$Dec <- ifelse(grepl('December', metadata$condition), 1, 0)
metadata$Jan <- ifelse(grepl('January', metadata$condition), 1, 0)
metadata$Jun <- ifelse(grepl('June', metadata$condition), 1, 0)
metadata$male <- ifelse(metadata$sex == 'male', 1, 0)
metadata[is.na(metadata$sex),]$male <- 0
metadata$female <- ifelse(metadata$sex == 'female', 1, 0)
metadata[is.na(metadata$sex),]$female <- 0
metadata <- metadata[-c(1, 2, 3, 4, 5, 6)]

# 7. Constructing the gene network and identifying modules: 
sth_power <- 4 # soft-threshold power 
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
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
}

correlation_heatmap(power=sth_power, MCH=0.15, MMS=25, ds=4, species = 'Eve',
                    plot_width=1800, plot_height=1000,
                    nameOfexp = 'Seasons_BK_june7C')

### Gene Significance (GS) and Module Membership (MM)
trait <- 'daylight'
temperature <- as.data.frame(metadata[trait]) # take the trait of interest
names(temperature) <- trait
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(dat_wo_missing, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
geneTraitSignificance <-  as.data.frame(cor(dat_wo_missing, temperature, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(temperature), sep="")
names(GSPvalue) <-paste("p.GS.", names(temperature), sep="")

### Plot the genes from one module

module <- "green"
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

modOrder <- order(-abs(cor(MEs, temperature, use = "p"))) # order modules by their significance for temperature

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

proteinInfo_final_hub <- proteinInfo_final[abs(proteinInfo_final[sprintf('GS.%s', 
                                                                         trait)]) 
                                           > 0.60,]

hubsMM <- NULL
for (protein in 1:nrow(proteinInfo_final_hub)){
  temp <- proteinInfo_final_hub[sprintf(c('MM.%s', 'p.MM.%s'), 
                  as.character(proteinInfo_final_hub[protein,][3]))][protein,]
  colnames(temp) <- c('MM', 'pMM')
  hubsMM <- rbind(hubsMM, temp)
}
proteinInfo_hub_reduced <- data.frame(proteinInfo_final_hub[c(1:5)], hubsMM)
proteinInfo_hub_reduced <- subset(proteinInfo_hub_reduced, abs(MM) > 0.60)

rgb.in <- col2rgb(proteinInfo_hub_reduced$moduleColor)
hsv_colors <- rgb2hsv(rgb.in)
#hsv_colors[2:3, ] <- 0.5
hub_colors_bright <- apply(hsv_colors, 2, function(row) hsv(row[1], row[2], row[3]))

proteinInfo_hub_reduced$geneSymbol <- sub('PREDICTED: |LOW QUALITY PROTEIN: ', '', 
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
                  seed = 555, box.padding = 0.5, point.padding = 0.4, size = 3.5,
                  segment.alpha = 0.5, max.overlaps = 15) +
  theme_light() + 
  scale_y_continuous(breaks = seq(-1.00, 1.00, 0.25), limits = c(-1,1)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(), 
        panel.grid = element_line(color = 'white'), 
        panel.border = element_rect(color = 'white'))

##############################################################
intensity_scaled <- apply(data_for_wgcna, 1, scale)
intensity_scaled <- t(intensity_scaled)
intensity_scaled <- data.frame(intensity_scaled)
colnames(intensity_scaled) <- rownames(metadata)

intensity_long <- intensity_scaled %>% 
  rownames_to_column(var = 'protein') %>%
  pivot_longer(-protein, values_to = 'intensity')

pr_id <- match(intensity_long$protein, proteinInfo_final$protein)
intensity_long$module <- proteinInfo_final[pr_id,]$moduleColor
intensity_long$annotation <- proteinInfo_final[pr_id,]$geneSymbol
intensity_long$GS <- proteinInfo_final[pr_id,]$GS.daylight
intensity_long <- subset(intensity_long, abs(GS)>0.6)

intensity_long$daylight <- metadata[match(intensity_long$name, 
                                          rownames(metadata)),]$daylight
#intensity_long$daylight <- factor(intensity_long$daylight,
#                                  levels = c('11.78', '9.37', '7.92', '8.73', '16.6'))
intensity_long_sub <- subset(intensity_long, module == 'green')

ggplot(intensity_long_sub, aes(daylight, intensity, group = protein)) +
  geom_line() +
  theme_light()

ggplot(intensity_long, aes(daylight, intensity, group = protein)) +
  geom_line(color = 'grey80', alpha = 0.6) +
  geom_line(data = intensity_long_sub, aes(daylight, intensity), 
            color = 'forestgreen') +
  theme_light()

ggplot(intensity_long_sub, aes(daylight, intensity)) +
  geom_jitter()

aggregate()

