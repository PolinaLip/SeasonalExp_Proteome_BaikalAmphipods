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
dir_to_intensities <- paste0('labeglo2/MS_results/Field/', species, 
                             '/', species,'_refchannels_all/')
current_dir <- paste0('labeglo2/MS_results/Field/', species, 
                      '/', species,'_refchannels_all/withALLzeros_inAtLeastOneCond/')
pg_annot_complete <- read.csv(file.path(dir_annot, 
                                        paste0('annot_proteinGroups_field_', tolower(species), 
                                               '_withGOannotation_with_organisms.csv')), 
                              sep = '\t', header = T) # eggnog and diamond annotation
metafile <- 'metadata_predictedSex.csv'

# 2. Upload normalized intensities
data_for_wgcna <- read.csv(file.path(current_dir, 
                                     paste0('intensities_after_slNorm_woNA_', 
                                            tolower(species),'.csv')),
                           header = T, sep = '\t')
# if to take into analysis data with missing values (zeros)
data_for_wgcna <- read.csv(file.path(dir_to_intensities, 
                                     paste0('intensities_after_slNorm_withNA', 
                                            tolower(species),'.csv')),
                           header = T, sep = '\t') 

data_for_wgcna <- data_for_wgcna[, !grepl('pool', colnames(data_for_wgcna),
                                          ignore.case = T)]

### If you work with data with missing values: ####

rownames(data_for_wgcna) <- data_for_wgcna[,1]
data_for_wgcna <- data_for_wgcna[-1]

# 3. Upload metadata
meta_complete <- read.delim(file.path(dir_metafile, metafile), header=TRUE, sep="\t")
meta <- meta_complete
meta <- subset(meta, condition != 'pool')
colnames(data_for_wgcna) <- meta$sample
#meta <- subset(meta, sex2 != 'NA')

### choose sex if you want to analyse only one sex
#meta <- subset(meta, sex2 == 'male')
#meta <- subset(meta, sex2 == 'female')
###
data_for_wgcna <- 
  data_for_wgcna[,colnames(data_for_wgcna) %in% meta$sample]

# 4. Choose only proteins with at least 5 not NA samples in each condition
remove_prot_with_someNA_inCond <- function(data2clean, metafile, notNAnumber = 5){
  cond <- as.factor(metafile$condition)
  ixs <- order(cond)
  data2clean <- data2clean[, ixs]
  cond <- cond[ixs]
  
  for (level in levels(cond)) {
    col_ixs <- cond == level
    data2clean <- data2clean[rowSums(!is.na(data2clean[,col_ixs])) >= notNAnumber,]
  }
  return(data2clean)
}

data_for_wgcna_ <- remove_prot_with_someNA_inCond(data_for_wgcna, meta)
col_indx <- match(colnames(data_for_wgcna), colnames(data_for_wgcna_))
data_for_wgcna <- data_for_wgcna_[col_indx] # to return the order of columns

# 4. Data preprocessing - get rid of missing data 
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

# 5. Choosing the soft-thresholding power: analysis of network topology
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

# 6. Prepare table with traits
metadata <- meta
row.names(metadata) <- metadata$sample
metadata$temperature <- ifelse(grepl('September', metadata$condition), 10.5,
                               ifelse(grepl('November', metadata$condition), 6.4,
                                      ifelse(grepl('December', metadata$condition), 2.3,
                                             ifelse(grepl('January', metadata$condition), 0.3, 7))))
metadata$time <- ifelse(grepl('September', metadata$condition), 1,
                        ifelse(grepl('November', metadata$condition), 2,
                               ifelse(grepl('December', metadata$condition), 3,
                                      ifelse(grepl('January', metadata$condition), 4, 5))))
metadata$Sep <- ifelse(grepl('September', metadata$condition), 1, 0)
metadata$Nov <- ifelse(grepl('November', metadata$condition), 1, 0)
metadata$Dec <- ifelse(grepl('December', metadata$condition), 1, 0)
metadata$Jan <- ifelse(grepl('January', metadata$condition), 1, 0)
#metadata$Jun <- ifelse(grepl('June', metadata$condition), 1, 0)

metadata$female <- ifelse(metadata$sex2 == 'female', 1, 0)
metadata$male <- ifelse(metadata$sex2 == 'male', 1, 0)
#metadata$female <- ifelse(metadata$sex2 != 'male' & !is.na(metadata$sex2), 1, 0)
#metadata$male <- ifelse(metadata$sex2 == 'male' | metadata$condition == 'June_BK', 1, 0)
metadata$Sep_f <- 
  ifelse(grepl('September', metadata$condition) & metadata$sex2 == 'female', 1, 0)
metadata$Sep_m <- 
  ifelse(grepl('September', metadata$condition) & metadata$sex2 == 'male', 1, 0)
metadata$Nov_f <- 
  ifelse(grepl('November', metadata$condition) & metadata$sex2 == 'female', 1, 0)
metadata$Nov_m <- 
  ifelse(grepl('November', metadata$condition) & metadata$sex2 == 'male', 1, 0)
metadata$Dec_f <- 
  ifelse(grepl('December', metadata$condition) & metadata$sex2 == 'female', 1, 0)
metadata$Dec_m <- 
  ifelse(grepl('December', metadata$condition) & metadata$sex2 == 'male', 1, 0)
metadata$Jan_f <- 
  ifelse(grepl('January', metadata$condition) & metadata$sex2 == 'female', 1, 0)
metadata$Jan_m <- 
  ifelse(grepl('January', metadata$condition) & metadata$sex2 == 'male', 1, 0)
#metadata$Jun_m <- 
#  ifelse(grepl('June', metadata$condition), 1, 0)

metadata$ampl_f <- ifelse(metadata$Sep_f == 1 | metadata$Nov_f == 1, 1, 0)
metadata$ampl_m <- ifelse(metadata$Sep_m == 1 | metadata$Nov_m == 1, 1, 0)

metadata <- metadata[-c(1, 2, 3, 4, 5, 6, 7, 8)]
# 7. Constructing the gene network and identifying modules: 
sth_power <- 4 # soft-threshold power
sth_power <- 6 # both sexes
sth_power <- 5
sth_power <- 5 # male
sth_power <- 7 # female
MMS <- 25
MCH <- 0.15
DS <- 4
net <- blockwiseModules(dat_wo_missing, power = sth_power,
                        TOMType = "signed", minModuleSize = MMS,
                        reassignThreshold = 0, mergeCutHeight = MCH,
                        maxBlockSize = 12000,
                        deepSplit = DS, 
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

correlation_heatmap(power=sth_power, MCH=MCH, MMS=MMS, ds=DS, species = species,
                    plot_width=2200, plot_height=800,
                    nameOfexp = 'SeasonsWithKnownSex_MonthSex_June2')

### Module Membership (MM) and Gene Significance (GS) calculation
trait <- 'female'
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

### Summary output of network analysis results

proteins <- names(dat_wo_missing)
samples2annot <- match(proteins, pg_annot_complete$protein_group)
sum(is.na(samples2annot))

proteinInfo <- data.frame(protein = proteins,
                          geneSymbol = pg_annot_complete$upd_full_annot[samples2annot],
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
                 GS_treshold, '_MM', MM_treshold, '.png')), scale = 1.7)

##################### Prepare long table #########################################
intensity_scaled <- apply(data_for_wgcna, 1, scale)
#intensity_scaled <- apply(data_for_wgcna, 1, scale)
intensity_scaled <- t(intensity_scaled)
intensity_scaled <- data.frame(intensity_scaled)
colnames(intensity_scaled) <- rownames(metadata)

intensity_long <- intensity_scaled %>% 
  rownames_to_column(var = 'protein') %>%
  pivot_longer(-protein, values_to = 'intensity', names_to = 'sample')

pr_id <- match(intensity_long$protein, proteinInfo_final$protein)
intensity_long$module <- proteinInfo_final[pr_id,]$moduleColor
intensity_long$annotation <- proteinInfo_final[pr_id,]$geneSymbol

##### GO enrichment #############################################################
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
module <- 'cyan'
intensity_long$MM <- proteinInfo_final[pr_id, paste0('MM.', module)]

intensity_long$pMM <- proteinInfo_final[pr_id, paste0('p.MM.', module)]

intensity_long$pGS <- proteinInfo_final[pr_id, paste0('p.GS.', trait)]
intensity_long$GS <- proteinInfo_final[pr_id, paste0('GS.', trait)]

intensity_long$go_annotation <- 
  pg_annot_complete[match(intensity_long$protein, 
                          pg_annot_complete$protein_group),]$go_annotation
intensity_long$upd_full_annotation <- 
  pg_annot_complete[match(intensity_long$protein, 
                          pg_annot_complete$protein_group),]$upd_full_annot
# create the dataframe with proteins and information to them:
unique_proteins <- group_by(intensity_long, protein) %>% 
  slice_head() %>%
  ungroup()
unique_proteins <- subset(unique_proteins, !is.na(module))
# take wanted statistics 
proteins_list <- unique_proteins$MM
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
threshold <- 0.5
ontology <- 'BP'
GOdata <- new("topGOdata", description = "Simple session", 
              ontology = ontology,
              allGenes = proteins_list4, 
              #geneSel = function(x) x < 0.01, # pMM, pGS
              #geneSel = function(x) abs(x) > threshold, # MM, GS
              geneSel = function(x) x > threshold, # MM
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
  genes_names <- paste0(names(wanted_[wanted_ < -threshold]), collapse = '', sep = ';')
  names(genes_names) <- 'significant_genes'
  go_table_w_genes <- rbind(go_table_w_genes, 
                            as.data.frame(c(current_row, genes_names)))
}

write.table(go_table_w_genes, 
            file = file.path(current_dir, 
                             paste0('topgo_ONLYpredSEX_', module, 
                                    'Module_', threshold, 
                                    '_ontology', ontology,'.csv')), row.names = F)

more_or_less <- 'lessMinus'
write.table(go_table_w_genes, 
            file = file.path(current_dir, 
                             paste0('topgo_ONLYpredSEX_', trait, 
                                    '_Trait_', more_or_less, threshold, 
                                    '_ontology', ontology,'.csv')), row.names = F)

### Plot significant proteins from chosen GO terms (boxplots) 

intensity_long$annot_long <- 
  pg_annot_complete[match(intensity_long$protein,
                          pg_annot_complete$protein_group),]$diamond_annot

get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}

intensity_long$whole_label <- sprintf('%s|%s', intensity_long$protein,
                                          intensity_long$go_annotation)
intensity_long$whole_label2 <- factor(intensity_long$whole_label,
                                          levels=unique(
                                            intensity_long$whole_label[order(intensity_long$go_annotation)]))

GO2draw <- 'GO:0006165'
wanted <- proteins_list4[names(proteins_list4) %in% allGO[[GO2draw]]] 
proteints2plot <- wanted[abs(wanted) > threshold]
intensity_long_toplot <- subset(intensity_long, MM %in% proteints2plot)

intensity_long_toplot$month <- meta[match(intensity_long_toplot$sample,
                                              meta$sample),]$condition

intensity_long_toplot$month <- sub('_BK', '', intensity_long_toplot$month)

intensity_long_toplot$time <- ifelse(intensity_long_toplot$month == 'September', 1,
                                         ifelse(intensity_long_toplot$month == 'November' , 2,
                                                ifelse(intensity_long_toplot$month == 'December', 3,
                                                       ifelse(intensity_long_toplot$month == 'January', 4, 5))))
intensity_long_toplot$sex <- meta[match(intensity_long_toplot$sample, 
                                        meta$sample),]$sex2

intensity_long_toplot_f <- subset(intensity_long_toplot, sex == 'female')
intensity_long_toplot_m <- subset(intensity_long_toplot, sex == 'male')
months_mean_female <- aggregate(intensity ~ month, intensity_long_toplot_f, mean)
months_mean_male <- aggregate(intensity ~ month, intensity_long_toplot_m, mean)
months_mean_female$time <- ifelse(months_mean_female$month == 'September', 1,
                           ifelse(months_mean_female$month == 'November', 2,
                           ifelse(months_mean_female$month == 'December', 3,
                           ifelse(months_mean_female$month == 'January', 4, 5))))
months_mean_male$time <- ifelse(months_mean_male$month == 'September', 1,
                           ifelse(months_mean_male$month == 'November', 2,
                           ifelse(months_mean_male$month == 'December', 3,
                           ifelse(months_mean_male$month == 'January', 4, 5))))

#intensity_long_toplot <- subset(intensity_long_toplot, intensity < 4) 

intensity_long_toplot$sex <- meta[match(intensity_long_toplot$sample,
                                            meta$sample),]$sex2

(ggplot(intensity_long_toplot, aes(factor(time), intensity)) +
  facet_wrap(~ whole_label2, ncol = 4,
             labeller = as_labeller(get_protein_label), scales = 'free_y') +
  geom_boxplot(outlier.color = NA, aes(color = sex)) +
  geom_jitter(size = 0.3,
              aes(color = sex), # sex
              position=position_jitterdodge(jitter.width = 0.3)) +
  geom_line(data = months_mean_female, aes(time, intensity, group = 1), 
            color = "deeppink1", alpha = 0.40, 
            size = 1.2) +
  geom_line(data = months_mean_male, aes(time, intensity, group = 1), 
            color = "dodgerblue1", alpha = 0.40, 
            size = 1.2) +
  scale_x_discrete(breaks = c(1, 2, 3, 4, 5),
                   labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
  scale_color_manual('Sex', values = c('deeppink1', 'dodgerblue1'), # sex 
                       na.value = 'grey35') + # sex
#  guides(fill = F
#         , color = F) + # sex
  #) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = 'bottom',
        legend.margin = margin(-10, 0, 0, 0)))

ggsave(file.path(current_dir, paste0('boxplots_ONLYpredSEX', module, 
                                     'Module_MM', threshold,
                                     '_ontology', ontology, '_',  GO2draw,'.png')), 
       scale = 0.8,
       width = 9, height = 5)

#############################################################################
intensity_long$GS <- proteinInfo_final[pr_id,]$GS.Sep_m
intensity_long$time <- metadata[match(intensity_long$sample, 
                                      rownames(metadata)),]$time
intensity_long$Sex <- meta[match(intensity_long$sample, 
                                 meta$sample),]$sex2

th <- 0.5
intensity_long_gs_sub <- subset(intensity_long, GS > th)
#intensity_long_gs_sub <- 
#  intensity_long_gs_sub[toupper(intensity_long_gs_sub$go_annotation) %in% names(gene_go_all_),]
intensity_long_gs_sub <- to_plot
get_protein_label <- function(x) {
  sapply(strsplit(as.character(x), '|', fixed=T), `[`, 2)
}
intensity_long_gs_sub$annotation <- 
  sub('PREDICTED: |LOW QUALITY PROTEIN: |-like|, partial| isoform X\\d+', '', 
      intensity_long_gs_sub$annotation)
intensity_long_gs_sub$annotationWrap <- sapply(
  intensity_long_gs_sub$annotation, function(x)
    paste(strwrap(x, width=26), collapse='\n'))

intensity_long_gs_sub$whole_label <- sprintf('%s|%s', intensity_long_gs_sub$protein,
                                             intensity_long_gs_sub$annotationWrap)
intensity_long_gs_sub$whole_label2 <- 
  factor(intensity_long_gs_sub$whole_label,
    levels=unique(
    intensity_long_gs_sub$whole_label[order(intensity_long_gs_sub$annotationWrap)]))
#lm_res <- lm(intensity ~ time, intensity_long_gs_sub)
#lm_res <- lm(intensity ~ daylight, intensity_long_gs_sub)
# time
ggplot(intensity_long_gs_sub, aes(factor(time), intensity)) +
  facet_wrap(~ whole_label2, ncol = 4,
             labeller = as_labeller(get_protein_label),
             scales = 'free_y') +
  geom_boxplot(outlier.color = NA, 
               aes(color = Sex)
               #aes(color = time == 4, fill = time == 4), alpha = .3
  ) +
  geom_jitter(
    aes(color = Sex), position = position_jitterdodge(jitter.width = 0.3),
    size = 0.3) +
#  geom_abline(intercept = lm_res$coefficients[1], slope = lm_res$coefficients[2],
#              color = 'firebrick') +
  scale_color_manual('Sex', values = c('deeppink1', 'dodgerblue1'), # sex 
                     na.value = 'grey35') +
  #  scale_fill_manual(values = c(NA, 'dodgerblue4')) + 
  scale_x_discrete(breaks = c(1, 2, 3, 4, 5),
                     labels = c('Sep', 'Nov', 'Dec', 'Jan', 'Jun')) +
#  guides(fill = F, color = F) +
  xlab('Months') +
  ylab('Scaled intensities') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.margin = margin(-10, 0, 0, 0))

ggsave(file.path(current_dir, paste0('boxplots_ONLYpredSEX_', trait,'Trait_', 
                                     th, 'GS',
                                     '.png')), 
       scale = 0.8,
       width = 10, height = 8)


to_plot <- subset(intensity_long, annotation == 'MDH2')



