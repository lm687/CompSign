<!-- ![logo simplex](compsign2.png "") -->
<img src="compsign2.png" alt="logo" width="400"/>

CompSign: An R package for differential abundace of compositional mutational signatures

# Installation

    library(devtools)
    devtools::install_github("lm687/CompSign")



setwd("/Users/morril01/Documents/PhD/CompSign")
setwd("/Users/morril01/Documents/PhD/")
system("R CMD build CompSign/")
system("R CMD check CompSign_0.1.0.tar.gz")
system("R CMD INSTALL CompSign_0.1.0.tar.gz")

# Display
## Display A: files that contain exposures
X, Z matrices



