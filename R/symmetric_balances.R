##' Implementing symmetric balances from
##' Correlation Between Compositional Parts Based on Symmetric Balances
##' Petra Kyncˇlová · Karel Hron · Peter Filzmoser
##'
##' which is a compositional data analysis-like correlation coefficient

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(MCMCpack)
library(circlize)
library(compositions)
library(ComplexHeatmap)


##
pdf(file = "~/Desktop/symmetric_balances.pdf")
for(i in list.files("../../CompSign/data", pattern = "*rda", full.names = TRUE)){load(i)}
dataset <- list(CNA_12K_TCGA=CNA_12K_TCGA,
                CNA_12K_TCGA_v2_QP=CNA_12K_TCGA_v2_QP,
                CNA_12K_TCGA_v2_YAPSA=CNA_12K_TCGA_v2_YAPSA)
version_dataset <- 'CNA_12K_TCGA_v2_QP'
exposures_tmp <- count_matrix(dataset[[version_dataset]])
print(plotcomputeRho(x = exposures_tmp, pseudocount =  1e-7,
               names_sigs = paste0('s', 1:7),
               column_title='CNA TCGA QP dataset - rho')) ## why is it all so extremely high?

print(plotcomputeclrcor(x = exposures_tmp, pseudocount =  1e-7,
               names_sigs = paste0('s', 1:7),
               column_title='CNA TCGA QP dataset - cor clr'))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

Ncomp <- 5
Nsamps <- 100
X <- MCMCpack::rdirichlet(Nsamps, rep(1,Ncomp))
print(plotcomputeRho(X, names_sigs = paste0('S', 1:Ncomp),
      column_title='CNA TCGA QP dataset - rho'))
print(plotcomputeclrcor(X, names_sigs = paste0('S', 1:Ncomp),
      column_title='CNA TCGA QP dataset - cor clr'))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

print(plotcomputeRho(count_matrix(Breast560), pseudocount = 1e-7,
               names_sigs = colnames(count_matrix(Breast560)),
               column_title = '560 Breast cancers (Nik-Zainal) - rho'))

print(plotcomputeclrcor(count_matrix(Breast560), pseudocount = 1e-7,
                  names_sigs = colnames(count_matrix(Breast560)),
                  column_title = '560 Breast cancers (Nik-Zainal) - cor clr'))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dev.off()

