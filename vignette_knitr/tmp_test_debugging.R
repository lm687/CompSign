### test/debugging scripts for the package
### once done, move to vignette knitr

rm(list = ls())

library(CompSign)
library(compositions)
data("three_comp_linear_relation")
data("two_normal_pops")
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

png('../CDA_in_Cancer/text/other_presentations/20181015/ex_cluster1.png')
plot(acomp(two_normal_pops_downsampling@count_matrix),
     col=as.factor(two_normal_pops_downsampling@df[,1]),
     pch=4)
dev.off()

## dendrogram
library(dendextend)
tmp_toplot_dendro <- as.dendrogram(hclust(dist(acomp(two_normal_pops_downsampling@count_matrix))))
tmp_toplot_dendro_nonComp <- as.dendrogram(hclust(dist(two_normal_pops_downsampling@count_matrix)))
labels_colors(tmp_toplot_dendro) <- c('black', 'red')[as.factor(two_normal_pops_downsampling@df[,1][labels(tmp_toplot_dendro)])]
labels_colors(tmp_toplot_dendro_nonComp) <- c('black', 'red')[as.factor(two_normal_pops_downsampling@df[,1][labels(tmp_toplot_dendro_nonComp)])]
labels(tmp_toplot_dendro) <- labels(tmp_toplot_dendro_nonComp) <-  rep('.', nrow(two_normal_pops_downsampling@df))
png('../CDA_in_Cancer/text/other_presentations/20181015/ex_cluster2.png')
plot(tmp_toplot_dendro, main='Aitchison distance')
dev.off()
png('../CDA_in_Cancer/text/other_presentations/20181015/ex_cluster3.png')
plot(tmp_toplot_dendro_nonComp, main='Euclidean distance')
dev.off()

var.acomp(acomp(two_normal_pops@count_matrix))
