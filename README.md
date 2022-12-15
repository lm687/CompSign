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

# Variations of the models

| name of model  | description  |cpp file (no need to use)   |   |   |
|---|---|---|---|---|
|   | DM with independent RE and one lambda  | diagRE_dirichletmultinomial_single_lambda  |   |   |
|   | DM with independent RE and two lambdas  | diagRE_ME_dirichletmultinomial  |   |   |
|   | M with independent RE  | diagRE_ME_multinomial  |   |   |
|   | DM with no RE and one lambda  | FE_dirichletmultinomial_single_lambda  |   |   |
|   | DM with no RE and two lambdas  | FE_dirichletmultinomial  |   |   |
|   | DM with independent RE and two lambdas  | fullRE_dirichletmultinomial_single_lambda  |   |   |
|   | DM assuming that there is no overdispersion in the first group (fixed lambda=1)  | fullRE_ME_dirichletmultinomial_onefixedlambda  |   |   |
|   | DM with correlated RE and two lambdas  | fullRE_ME_dirichletmultinomial  |   |   |
|   | M with correlated RE  | fullRE_ME_multinomial  |   |   |
|   | DM with a single RE intercept and two lambdas  | singleRE_dirichlet_multinomial  |   |   |




