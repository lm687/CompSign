### test/debugging scripts for the package
### once done, move to vignette knitr

rm(list = ls())

## .rs.restartR()

library(CompSign)
library(compositions)
data("three_comp_linear_relation")
data("two_normal_pops")
data("two_normal_pops_extended")
plot(three_comp_linear_relation)

merged_compositional_to_sign(two_normal_pops)


plot(acomp(two_normal_pops@count_matrix),
     col=as.factor(two_normal_pops@df[,1]),
     pch=4)

plot(hclust(dist(acomp(two_normal_pops@count_matrix))))

two_normal_pops_downsampling <- two_normal_pops
downsampling_idx <- sample(nrow(two_normal_pops_downsampling@count_matrix),
                           size = 0.1*nrow(two_normal_pops_downsampling@count_matrix))
#two_normal_pops_downsampling@count_matrix <- two_normal_pops_downsampling@count_matrix[downsampling_idx,]
#two_normal_pops_downsampling@df <- data.frame(col=two_normal_pops_downsampling@df[downsampling_idx,])

# png('../CDA_in_Cancer/text/other_presentations/20181015/ex_cluster1.png')
# plot(acomp(two_normal_pops_downsampling@count_matrix),
#      col=as.factor(two_normal_pops_downsampling@df[,1]),
#      pch=4)
# dev.off()
#
# # dendrogram
# library(dendextend)
# tmp_toplot_dendro <- as.dendrogram(hclust(dist(acomp(two_normal_pops_downsampling@count_matrix))))
# tmp_toplot_dendro_nonComp <- as.dendrogram(hclust(dist(two_normal_pops_downsampling@count_matrix)))
# labels_colors(tmp_toplot_dendro) <- c('black', 'red')[as.factor(two_normal_pops_downsampling@df[,1][labels(tmp_toplot_dendro)])]
# labels_colors(tmp_toplot_dendro_nonComp) <- c('black', 'red')[as.factor(two_normal_pops_downsampling@df[,1][labels(tmp_toplot_dendro_nonComp)])]
# labels(tmp_toplot_dendro) <- labels(tmp_toplot_dendro_nonComp) <-  rep('.', nrow(two_normal_pops_downsampling@df))
# png('../CDA_in_Cancer/text/other_presentations/20181015/ex_cluster2.png')
# plot(tmp_toplot_dendro, main='Aitchison distance')
# dev.off()
# png('../CDA_in_Cancer/text/other_presentations/20181015/ex_cluster3.png')
# plot(tmp_toplot_dendro_nonComp, main='Euclidean distance')
# dev.off()

createDendrogram(two_normal_pops, bool_comparison = TRUE)
createDendrogram(two_normal_pops_downsampling, bool_comparison = TRUE)

var.acomp(acomp(two_normal_pops@count_matrix))
var.acomp(acomp(rdirichlet(1e4, c(1,1,1))))

## test for normality
qqnorm.acomp(acomp(two_normal_pops@count_matrix), pch=19, cex=0.2)
qqnorm.acomp(acomp(two_normal_pops@count_matrix[1:1000,]), pch=19, cex=0.2)

hist(alr(acomp(two_normal_pops@count_matrix[1:1000,]))[,1], breaks=100)
hist(alr(acomp(two_normal_pops@count_matrix[1:1000,]))[,2], breaks=100)

## example where the presence of two mutations is equivalent
## TODO

###############################
## cat variable as predictor ##
contrasts(as.factor(unlist(two_normal_pops_extended@df)))
boxplot(acomp(two_normal_pops_extended@count_matrix), unlist(two_normal_pops_extended@df))

outer(X = 1:8, Y = unlist(two_normal_pops_extended@df),
      FUN = function(X,Y) boxplot(acomp(two_normal_pops_extended@count_matrix)[,X],
                    Y))

model_discrete <- lm(ilr(acomp(two_normal_pops_extended@count_matrix))~unlist(two_normal_pops_extended@df))
model_discrete_with_alr <- lm(alr(acomp(two_normal_pops_extended@count_matrix))~unlist(two_normal_pops_extended@df))
anova(model_discrete)
anova(model_discrete_with_alr)
coef(model_discrete)
coef(model_discrete_with_alr)
plot(coef(model_discrete_with_alr)[2,], type='l')
acomp(two_normal_pops@count_matrix)

model_discrete_scrambled <- lm(ilr(acomp(two_normal_pops@count_matrix))~sample(c(0, 1), length(unlist(two_normal_pops@df)), TRUE))
anova(model_discrete_scrambled)

## with LDA
library(MASS)
library(scales)
png('../CDA_in_Cancer/text/other_presentations/20181015/ex_lda1.png')
plot(acomp(two_normal_pops_extended@count_matrix),
     col=alpha(c('black', 'red')[as.factor(unlist(two_normal_pops_extended@df))],
               0.01), cex=0.2, pch=19)
dev.off()

lda_res <- lda(x=data.frame(ilr(acomp(two_normal_pops_extended@count_matrix))),
               grouping=unlist(two_normal_pops_extended@df))
ilrInv(lda_res$means, orig=acomp(two_normal_pops_extended@count_matrix))
V <- ilrBase(acomp(two_normal_pops_extended@count_matrix))
rownames(V) <- colnames(acomp(two_normal_pops_extended@count_matrix))
t(ilr2clr(t(lda_res$scaling), V=V))
plot3D(acomp(two_normal_pops_extended@count_matrix))

lda_res_alr <- lda(x=data.frame(alr(acomp(two_normal_pops_extended@count_matrix))),
               grouping=unlist(two_normal_pops_extended@df))

png('../CDA_in_Cancer/text/other_presentations/20181015/ex_lda2.png')
biplot(princomp(acomp(two_normal_pops_extended@count_matrix)))
dev.off()
biplot(princomp(alr(acomp(two_normal_pops_extended@count_matrix))))
biplot(princomp(ilr(acomp(two_normal_pops_extended@count_matrix))))

image(var.acomp(acomp(two_normal_pops_extended@count_matrix)))


tmp_pca <- princomp(acomp(two_normal_pops_extended@count_matrix))
tmp_pca$sdev^2/sum(tmp_pca$sdev^2)


