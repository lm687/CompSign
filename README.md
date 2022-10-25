![logo simplex](compsign.png ""))

CompSign: An R package for differential abundace of compositional mutational signatures

# Installation

    library(devtools)
    devtools::install_github("lm687/CompSign")



setwd("/Users/morril01/Documents/PhD/CompSign")
setwd("/Users/morril01/Documents/PhD/")
system("R CMD build CompSign/")
system("R CMD check CompSign_0.1.0.tar.gz")
system("R CMD INSTALL CompSign_0.1.0.tar.gz")
