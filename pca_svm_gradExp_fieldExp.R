# for sex prediction in grad exp. 
# take only Vgs, then leave only the first transcript names from the protein group -> overlap two tables, 
# then make pca, svm as well.
library(dplyr)
library(ggfortify) 
library(e1071)
library(ggplot2)
library(metR)
vg1 <- 'NODE_1944_length_5232_cov_28.357901_g1272_i0.p1_Gla'
vg2 <- 'TRINITY_DN3055_c0_g1_i31.p1_Eve'
vg3 <- 'TRINITY_DN5221_c0_g1_i3.p1_Eve'
vgs <- c(vg1, vg2, vg3)
scaled_intensities_field_vg <- 
  scaled_intensities_field[grepl(paste(vgs, collapse = '|'), 
                                 rownames(scaled_intensities_field)),]
scaled_intensities_field_vg <- 
  scaled_intensities_field_vg[,colnames(scaled_intensities_field_vg) != 'sorted_protein']
scaled_intensities_grad_vg <- 
  scaled_intensities_grad[grepl(paste(vgs, collapse = '|'), 
                                rownames(scaled_intensities_grad)),]
scaled_intensities_grad_vg <- 
  scaled_intensities_grad_vg[,colnames(scaled_intensities_grad_vg) != 'sorted_protein']


sc_int_all_samples <- cbind(scaled_intensities_field_vg, scaled_intensities_grad_vg)

### PCA
sc_int_all_samples <- sc_int_all_samples[,!grepl('June', meta_joined$condition)]
sc_int_all_samples_t <- na.omit(sc_int_all_samples) %>% t %>% as.data.frame
pca_res <- prcomp(sc_int_all_samples_t)
meta_joined_ <- subset(meta_joined, condition != 'June')
meta_joined_$condition <- sub('EveBK_1.5', 'Lab: 1.5 °C', meta_joined_$condition)
meta_joined_$condition <- sub('EveBK_10.5', 'Lab: 10.5 °C', meta_joined_$condition)
autoplot(pca_res, data=meta_joined_, colour='condition', size = 3,
         shape = 'sex',
         frame = F, 
         x = 1,
         y = 2) + 
  theme_light() + 
  scale_color_manual('Condition', 
                     values = palette.colors(n=8, 'Dark2'))

dir_to_results <- 'labeglo2/MS_results/Field_Grad_Comparison/'
ggsave(filename = file.path(dir_to_results, paste0('pca_', 
                  species,'_seasonalExp_GradDecreaseANDfield.png')), 
       scale = 1.2,
       width = 6, height = 3)

###
meta_joined_$sex_before_pred <- 
  ifelse(meta_joined_$condition %in% c('September', 'November'), meta_joined_$sex, 'NA')
res <- cbind(meta_joined_, as.data.frame(pca_res$x))
tmp <- sc_int_all_samples_t
colnames(tmp) <- c('vg1', 'vg2', 'vg3')
res <- cbind(meta_, tmp)
ggplot(res) +
  geom_histogram(aes(PC1, fill=sex_before_pred), binwidth=0.05) +
  theme_bw()

#### svm

res$sex_fact <- with(res, factor(ifelse(sex_before_pred == 'NA', NA, sex)))
res_part <- res[c('sex_fact', 'PC1', 'PC2')]
#res_part <- res[c('sex_fact', 'vg1', 'vg2')]
svm_model <- svm(sex_fact ~ PC1, data = res_part, cost = 1, probability = T)
#svm_model <- svm(sex_fact ~ vg1 + vg2, data = res_part, cost = 1, probability = T)

grange = apply(dplyr::select(res_part, starts_with('PC')), 2,
               function(x) { r <- range(x); c(r, r[2] - r[1]) }) %>%
  as.data.frame

grid <- expand.grid(
  PC1 = seq(grange$PC1[1] - grange$PC1[3] * 0.05,
            grange$PC1[2] + grange$PC1[3] * 0.05, length.out = 1000),
  PC2 = seq(grange$PC2[1] - grange$PC2[3] * 0.05,
            grange$PC2[2] + grange$PC2[3] * 0.05, length.out = 100))

cor(res$vg2, res$vg3, method='pearson')
summary(pca_res)

grid$sex_fact <- predict(svm_model, grid)
grid$sex_prob <- (predict(svm_model, grid, probability = T) %>%
                    attr('probabilities'))[,1]

# ggplot(grid, aes(PC1, PC2)) +
#   geom_tile(aes(fill = sex_prob)) +
#   #geom_point(data = res, size=2) +
#   theme_bw()

ggplot(grid, aes(PC1, PC2)) +
  geom_tile(aes(fill = sex_prob)) +
  geom_contour(aes(z = sex_prob), color='gray30') +
  geom_text_contour(aes(z = sex_prob), stroke = 0.1,
                    label.placer = label_placer_fraction(frac = 0.3)) +
  geom_point(aes(shape = sex_before_pred, 
                 alpha = sex_before_pred, 
                 color = sex_before_pred), data = res, size = 6) +
  scale_fill_gradientn('Probability\nto be male', colours = c('pink', 'white', 'dodgerblue'),
                       values = c(0, 0.5, 1)) +
  scale_shape_manual('Sex', values = c("\u2640", "\u2642", '*'), 
                     labels = c('Female', 'Male', 'Unknown')) +
  scale_alpha_manual('Sex', values = c(1, 1, 0.7), 
                     labels = c('Female', 'Male', 'Unknown')) +
  scale_color_manual('Sex', values = c('black', 'black', 'coral3'), 
                     labels = c('Female', 'Male', 'Unknown')) +
  #scale_fill_manual(values = colorRampPalette(c("pink", "white", "blue"))(10)) +
  #coord_cartesian(xlim = c(-0.185, 0.35)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw()

ggsave(filename = file.path(dir_to_results, paste0('pca_', 
       species,'_afterSvm_to_predict_sex_FieldPlusGrad.png')), 
       scale = 0.65,
       width = 13, height = 6)

ggplot(grid, aes(PC1, PC2)) +
  geom_contour_filled(aes(z = sex_prob), binwidth = 0.1) +
  geom_point(aes(shape = sex2), data = res, size = 2) +
  #scale_fill_gradientn(colours = c('red', 'white', 'blue')) +
  scale_fill_manual(values = colorRampPalette(c("pink", "white", "blue"))(10)) +
  theme_bw()



