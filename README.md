<!-- ![logo simplex](compsign2.png "") -->
<img src="compsign3.png" alt="logo" width="400"/>

*CompSign*: An R package for differential abundance of compositional mutational signatures

Lena Morrill 2022

# Installation

    library(devtools)
    devtools::install_github("lm687/CompSign", build_vignettes = TRUE)
    
# Vignette

`vignette('CompSign')`

# Variations of the model

The first row is the `<model>` argument in the function `wrapper_run_TMB()`.
| name of model (for user) | description  |cpp file (no need to use)   |
|---|---|---|
| diagREDMsinglelambda  | DM with independent RE and one lambda  | diagRE_dirichletmultinomial_single_lambda  |
| diagRE_DM  | DM with independent RE and two lambdas  | diagRE_ME_dirichletmultinomial  |
| diagRE_M  | M with independent RE  | diagRE_ME_multinomial  |
| FEDMsinglelambda  | DM with no RE and one lambda  | FE_dirichletmultinomial_single_lambda  |
| FE_DM  | DM with no RE and two lambdas  | FE_dirichletmultinomial  |
| fullREDMsinglelambda  | DM with independent RE and two lambdas  | fullRE_dirichletmultinomial_single_lambda  |
| fullRE_DMonefixedlambda  | DM assuming that there is no overdispersion in the first group (fixed lambda=1)  | fullRE_ME_dirichletmultinomial_onefixedlambda  |
| fullRE_DM  | DM with correlated RE and two lambdas  | fullRE_ME_dirichletmultinomial  |
| fullRE_M  | M with correlated RE  | fullRE_ME_multinomial  |
| singleRE_DM  | DM with a single RE intercept and two lambdas  | singleRE_dirichlet_multinomial  |




