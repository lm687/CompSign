<!-- ![logo simplex](compsign2.png "") -->
<img src="compsign3.png" alt="logo" width="400"/>

**CompSign**: An *R* package for differential abundance of compositional mutational signatures

Lena Morrill Gavarró 2024

## Installation

[Pandoc](https://pandoc.org/installing.html) is required to build R Markdown vignettes - please install it before installing the package as follows:

    install.packages('devtools') ## if needed
    library(devtools)
    devtools::install_github("lm687/CompSign", build_vignettes = TRUE)

## Vignettes
Several examples of input data, how to run the models, and how to interpret the results are found in the vignette:

    browseVignettes("CompSign")

This package as, applied to study the differences in mutational signatures between clonal and subclonal mutations in the PCAWG dataset. The github repository reproducing these results can be [found here](https://github.com/lm687/CompSign-results).

## Brief summary of the package
We refer to the vignettes for a more in-depth explanation of models and the functioning of the package, but a minimal example is found here:


## How to run the model

The function `wrapper_run_TMB()` is used to run all variations of the model.

### Input dataset

The input dataset is the argument `object`, which is a list with the following structure:
- `x`: covariate matrix (`p x n`)
- `z`: matrix of random effects indicating which patient-specific subsample belongs to which patient (`n x N`)
- `Y` (`n x d`)
- `d`


### Estimating the parameters

A minimal example is:

```
diagDM_no_small_sigs <- wrapper_run_TMB(object = simplified_object,
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
```

in which `simplified_object` is a list containing `simplified_object$x`, `simplified_object$z`, `simplified_object$Y`. In particular, `simplified_object` corresponds to the nucleotide changes in the Lung-AdenoCA cohort. You can load `simplified_object`, as well as several other example datasets:
```
data(package='CompSign')
data(simplified_object)
```

To test for differential abundance, a generalised Wald test can be used with the function `wald_TMB_wrapper`, which gives a p-value as output:

```
wald_TMB_wrapper(diagDM_no_small_sigs)
```

The estimated coefficients can be visualised as follows:
```
library(ggplot2)
plot_betas(diagDM_no_small_sigs, return_ggplot = T)
```

If signatures have not yet been extracted, this can be done with the wrapper function `extract_sigs_TMB_obj` and using as input a mutation matrix (commonly, a matrix with trinucleotide mutations; see more details in the Vignette).

### Variations of the model


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
| diagDMpatientlambda  | DM with independent RE and one lambda for each patient  | diagREpatientlambda_ME_dirichletmultinomial  |
| fullDMpatientlambda  | DM with correlated RE and one lambda for each patient  | fullREpatientlambda_ME_dirichletmultinomial  |

