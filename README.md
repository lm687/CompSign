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
We recommend browing the vignettes for a more in-depth explanation of models and the functioning of the package, several examples of input data, how to run several models, and how to interpret the results. You can find the vignette here:

    browseVignettes("CompSign")

## Brief summary of the package
This is a minimal example for running the model and extracting the estimated coefficients.

### Input dataset

The input dataset is the argument `object`, which is a list with the following structure:
- `x`: covariate matrix (`p x n`)
- `z`: matrix of random effects indicating which patient-specific subsample belongs to which patient (`n x N`)
- `Y` (`n x d`)

### Fitting the model

The function `wrapper_run_TMB()` is used to run all variations of the model.

A minimal example is:

```
diagDM_no_small_sigs <- wrapper_run_TMB(object = simplified_object,
                                        model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
```

in which `simplified_object` is a list containing `simplified_object$x`, `simplified_object$z`, `simplified_object$Y`. The `object` input for `wrapper_run_TMB` is always a list containing a matrix of covariates `x`, a matrix of observations `Y` and a matrix matching observations patients and the labels corresponding to the random effects. All three matrices should have the same number of rows - the number of observations. In this case, `simplified_object` corresponds to the nucleotides abundances in clonal and subclonal mutations in the Lung-AdenoCA cohort. `simplified_object$x`, contains two columns (a column of 1 for the intercept and a column of 0 and 1 indicating whether the observation corresponds to clonal (0) or subclonal (1) mutations. `simplified_object$z` contains a matrix where each observation is paired with its corresponding patient, `simplified_object$Y` is a matrix of counts for the nucleotide abundances for each patient and (clonal/subclonal) group.

Here we are using the model `diagRE_DM`, which is a mixed-effect Dirichlet-multinomial model with non-correlated random effects and two dispersion parameters (one per group). See below for the list of alternative models available.

### Data included in CompSign
You can load `simplified_object`, as well as several other example datasets, using the `data` funtion:

```
data(package='CompSign')
data(simplified_object)
```

### Testing for differential abundance

Depending on the model, the output of `wrapper_run_TMB` contains the estimated values of different coefficients. The parameters of interest are `beta_0`, `beta_1` and `lambda`. Both `beta` are values in log-ratio space in with respect to the last category included in the input `object` for `wrapper_run_TMB` (i.e. the last colun of `object$Y`). `beta_0` indicates the coefficient for the baseline - in the two-group example of `simplified_object`, this baseline reflects the abundance in the first (clonal) group. `beta_1` indicates the cofficient for the second column - here representing the difference between clonal and subclonal mutations. If `beta_1` is a vector of 0, this indicates that none of the log-ratios of abundances have changed between the two groups, and hence that there is no differential abundance.

To test for differential abundance, a generalised Wald test can be used with the function `wald_TMB_wrapper`, which gives a p-value as output:

```
wald_TMB_wrapper(diagDM_no_small_sigs)
```

### Visualising the output

The estimated `beta` coefficients can be visualised as follows:

```
library(ggplot2)
plot_betas(diagDM_no_small_sigs, return_ggplot = T)
```

The interest of the user is generally on which categories increase or decrease in absolute terms. This is not possible in compositional data, but, assuming that the abundance of most categories does not change, we can visualise which categories deviate highly from the median coefficient. This visualisation can be generated with `plot_betas_minimal_perturbation`. The red line is the line of “zero perturbation” and for each estimate its 95% confidence interval is shown in blue. We can consider a signature to be differentially abundant if the line of “zero perturbation” is not included in the confidence interval of its coefficient for differential abundance. Note that we draw conclusions for the signatures included in the numerator of the log-ratios, but we are unable to get equivalent results for the baseline signature (although it is used to compute the line of “zero perturbation”).

```
plot_betas_minimal_perturbation(diagDM_no_small_sigs)
```

### Usage for mutational signatures

In the example above, nucleotide abundances are used as input. Alternatively, mutational signatures (as counts) and be used instead. If signatures have not yet been extracted, this can be done with the wrapper function `extract_sigs_TMB_obj` and using as input a mutation matrix (commonly, a matrix with trinucleotide mutations; see more details in the Vignette).

`CompSign` has been applied to study the differences in mutational signatures between clonal and subclonal mutations in the PCAWG dataset. The github repository reproducing these results can be [found here](https://github.com/lm687/CompSign-results).


### Variations of the model

There are several variations on the model, each with a particular combination of the following:
- Base model: Dirichlet-multinomial in most cases, but two simpler models -- multinomial models `diagRE_M` and `fullRE_M` are implemented for comparison.
- Random effects: correlated or not (e.g. uncorrelated in `diagRE_DM`, correlated in `fullRE_DM`).
- Number of dispersion parameters: generally one per group (e.g. in `diagRE_DM`), but possibly fewer or more if specified in the name (e.g. one single dipersion parameter in `diagRE_DM_singlelambda`, one dispersion parameter per patient in `diagRE_DM_patientlambda`).

The names of the model correspond to `{random effects}_{base model}_{particularities of the dispersion parameter, if any}`.

The first column is the `<model>` argument in the function `wrapper_run_TMB()`.
| Name of model (for user) | Description and usage  |
|---|---|
| `FE_DM_singlelambda`  | Dirichlet-multinomial  with no RE and a single `lambda` which can be used in a two-group comparison (if we consider the dipersion to be the same in both groups), a multi-group comparison, or any regression setting.  | 
| `FE_DM`  | Dirichlet-multinomial  with no random effects and two `lambda`, used for the comparison between two groups in which we want to account for group-specific dispersion. Model slightly more complex than `FE_DM_singlelambda` | 
| `diagRE_M`  |  Multinomial with non-correlated multivariate effects. Simpler model than `diagRE_DM`, as no dispersion is included | 
| `diagRE_DM_singlelambda`  |  Dirichlet-multinomial with non-correlated multivariate RE and one: equivalent of diagRE_DM but with a shared . It can be used in a two-group comparison (if we consider the dipersion to be the same in both groups), a multi-group comparison, or any regression setting.  | 
| `singleRE_DM`  | Dirichlet-multinomial with a single RE intercept and two `lambda`: simple model in which random effects are not multivariate, but where group-specific dispersion is needed, in a two-group comparison. | 
| `diagRE_DM`  |  Dirichlet-multinomial with independent RE and two `lambda`: used most commonly throughout the paper. The data are matched to warrant multivariate RE, but correlations between categories are not explicitly modelled. Faster than `fullRE_DM` with often comparable results for `beta_1`.   | 
| `fullRE_DM_singlelambda`  | Dirichlet-multinomial with correlated RE and two lambdab: equivalent of `diagR_DM_singlelambda` that can be used if categories have strong correlations. | 
| `fullRE_M`  | Multinomial with correlated RE: assuming no overdispersion, it can be used in any regression setting  | 
| `fullRE_DM`  | Dirichlet-multinomial with correlated RE and two `lambda`: model with multivariate RE with correlations modelled explicitly. As a `(K-1)*(K-1)` covariance matrix is estimated, this model is not recommended where the ratio of signatures to samples is high.  |
| `diagRE_DM_patientlambda`  | Dirichlet-multinomial with non-correlated RE and one lambda$ per patient: model more complex than `diag_RE_DM` that can be used when there are several observations per patient and we wish to include a patient-specific dispersion parameter  | 
| `fullRE_DM_patientlambda`  |  Dirichlet-multinomial with correlated RE and one  per patient: equivalent of `diagRE_DM_patientlambda`, but modelling correlations between categories | 


<!---
[comment]: <>  |---|---|---|
[comment]: <>  diagREDMsinglelambda  | DM with independent RE and one lambda  | diagRE_dirichletmultinomial_single_lambda  |
[comment]: <> | diagRE_DM  | DM with independent RE and two lambdas  | diagRE_ME_dirichletmultinomial  |
[comment]: <> | diagRE_M  | M with independent RE  | diagRE_ME_multinomial  |
[comment]: <> | FEDMsinglelambda  | DM with no RE and one lambda  | FE_dirichletmultinomial_single_lambda  |
[comment]: <> | FE_DM  | DM with no RE and two lambdas  | FE_dirichletmultinomial  |
[comment]: <> | fullREDMsinglelambda  | DM with independent RE and two lambdas  | fullRE_dirichletmultinomial_single_lambda  |
[comment]: <> | fullRE_DMonefixedlambda  | DM assuming that there is no overdispersion in the first group (fixed lambda=1)  | fullRE_ME_dirichletmultinomial_onefixedlambda  |
[comment]: <> | fullRE_DM  | DM with correlated RE and two lambdas  | fullRE_ME_dirichletmultinomial  |
[comment]: <> | fullRE_M  | M with correlated RE  | fullRE_ME_multinomial  |
[comment]: <> | singleRE_DM  | DM with a single RE intercept and two lambdas  | singleRE_dirichlet_multinomial  |
[comment]: <> | diagDMpatientlambda  | DM with independent RE and one lambda for each patient  | diagREpatientlambda_ME_dirichletmultinomial  |
[comment]: <> | fullDMpatientlambda  | DM with correlated RE and one lambda for each patient  | fullREpatientlambda_ME_dirichletmultinomial  |
--->

