##https://github.com/lgatto/TeachingMaterial/blob/master/_roo/roo.pdf

require(compositions)
#require(MCMCpack)
#require(NMF)

## before calling nmf
setClass("count_data",
         representation = representation(
           id="character",
           id_samples = "character",
           id_categories = "character",
           count_matrix = "matrix",
           modified = "logical"
         ))

## extracted signatures
setClass("sign",
         representation = representation(
           id="character",
           id_samples = "character",
           id_signatures = "character",
           count_matrix = "matrix",
           modified = "logical"#,
           #ilr = "rmult" TODO: add ilr
         ))

setClass("metadata",
         representation = representation(
           id="character",
           id_samples="character",
           df= "data.frame"
         ))

## both metadata and the mutational signatures
setClass("merged_compositional",
         representation = representation(
           id="character",
           id_samples="character",
           id_signatures = "character",
           count_matrix = "matrix",
           df= "data.frame",
           pseudocounts="logical"
         ))

## alternative:
## foo <- function(x) {
##   if (!is.numeric(x)) stop("X must be numeric")
##   structure(list(x), class = "foo")
## }
to_sign <- function(x){
  x_sign_obj <- new("sign",
      id=deparse(substitute(x)), ## get name
      id_samples=rownames(x),
      id_signatures= colnames(x),
      count_matrix=x,
      modified=FALSE)
  x_sign_obj
}

##########################
comp_lda <- function(x, indices_response){
  if(class(x)[1] != 'merged_compositional') stop('Input must be of class <merged_compositional>')

}


###################################
############# CURRENT #############

comp_lm <- function(x, indices_predictor){
  ## the composition is the response

  if(class(x)[1] != 'merged_compositional') stop('Input must be of class <merged_compositional>')

  #require(compositions)
  apply(compositions::ilr(compositions::acomp(x@count_matrix)), 2, function(it_ilr){
   model_compReg <- compositions::lm(unlist(it_ilr~(x@df)[,indices_predictor]))
   list(summary(model_compReg), anova(model_compReg))
  })
}
comp_lm <- function(x, indices_predictor){
  ## the composition is the response

  if(class(x)[1] != 'merged_compositional') stop('Input must be of class <merged_compositional>')

  #require(compositions)
  model_compReg <- compositions::lm(compositions::ilr(compositions::acomp(x@count_matrix))~(x@df)[,indices_predictor])
  list(summary(model_compReg)#,
       #anova(model_compReg)
       )
}

############# CURRENT #############
###################################


##########################
####### Dummy data #######
##########################
# dum_sign <- new("sign",
#                 id="dum1",
#                 id_samples=c("sam1", "sam2", "sam3"),
#                 id_signatures= c('s1', 's2', 's3', 's4'), ## signature names
#                 count_matrix=MCMCpack::rdirichlet(3, c(1,1,1,1)),
#                 modified=FALSE)
# d <- readRDS("~/Desktop/tcga_ov_signature_exposures_20180926.rds")
# many_tcga <- new("sign",
#                  id='many_tcga',
#                  id_samples=rownames(d),
#                  id_signatures= colnames(d), ## signature names
#                  count_matrix=d,
#                  modified=FALSE)
#
# ## TODO: add types of response (categorical, numerical, ...)
# dum_metadata <- new("metadata",
#                     id="dum1",
#                     id_samples=c("sam1", "sam2", "sam3"),
#                     df=data.frame(age=sample(1L:100L,3),
#                                   some_cat=as.character(letters[sample(1L:length(letters), 3)]),
#                                   stringsAsFactors = FALSE)
# )

### Example of matrix transformed into sign object
# aaa <- matrix(runif(100), 4)
# colnames(aaa) <- paste0('s', 1:25); rownames(aaa) <- paste0('sam', 1:4)
# to_sign(aaa)


## put together the matrix and its row/col names
setGeneric("add_together_matrix", function(obj, ...) standardGeneric("add_together_matrix"))
setMethod("add_together_matrix",
          "sign",
          function(obj){
            colnames(obj@count_matrix) <- obj@id_signatures
            rownames(obj@count_matrix) <- obj@id_samples
            obj@modified <- TRUE

            ## TODO: add pseudocounts by default and throw warning
            return(obj)
          })

##########################
## Function 1 summarise ##
##########################

setGeneric("summarise", function(obj, ...) standardGeneric("summarise"))
setMethod("summarise", "sign", function(obj, ...){
  list(General=paste0("Object of class ", class(obj)),
  `Number of signatures` = ncol(obj@count_matrix),
  `Number of samples` = nrow(obj@count_matrix),
  `Geometric means of signatures`= sort(apply(obj@count_matrix, 2, function(x) exp(mean(log(x)))), decreasing=TRUE),
  `Covariance` = compositions::var.acomp(obj@count_matrix))
})

# summarise(add_together_matrix(dum_sign))
# summarise(many_tcga)
#
# many_tcga

################################
## Function 2 signature-co-oc ##
################################
setGeneric("within_signature_analysis", function(obj, ...) standardGeneric("summarise"))
#setMethod("summarise", "sign", function(obj, ...){
  ### calculate the co-ocurrence of mutationl signatures in samples
#})

############################################
## Function 3 comparison with other signs ##
############################################

# ?data

## NMF
## NMF::nmf(dum_sign@count_matrix, rank = 3)
## below: TODO: install sigfit
## devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)


##############################################
## Function 4 comparison with clinical data ##
##############################################
## boxplots for alr and each predictor:
## in exploring_correlations.R

## lm
## discriminant analysis
# dum_sign@id_samples
# dum_metadata@df
#
# tmp_merged_compositional <- new("merged_compositional",
#                                 id='adas',
#                                 id_samples=c("sam1", "sam2", "sam3"),
#                                 id_signatures= c('s1', 's2', 's3', 's4'), ## signature names
#                                 count_matrix=MCMCpack::rdirichlet(3, c(1,1,1,1)),
#                                 df=data.frame(a=c(3,4,1), b=c(10, 10, 10)))

link_to_clinical_data <- function(predictors, response){
  ## the predictors are the signatures (compositional)
  ## the response is a non-compositional dataframe
  if( !all((predictors@id_samples == response@id_samples)) ) stop('Data and metadata contain different samples')
  if( !all((predictors@id == response@id)) ) stop('Data and metadata seem to be different datasets (discordant id)')

  sap_types_col <- sapply(response@df, typeof)
  idx_types_col <- lapply(unique(sap_types_col),
               function(x) which (sap_types_col == x))
  names(idx_types_col) <- unique(sap_types_col)
  idx_types_col

  merged_object <-  new("merged_compositional",
                        id=predictors@id,
                        id_samples=predictors@id_samples,
                        id_signatures= predictors@id_signatures,
                        count_matrix=predictors@count_matrix,
                        df=response@df,
                        pseudocounts=FALSE
                        )

  print(typeof(merged_object))

  ## character response ##
  if(!is.null(idx_types_col$character)){
    for(idx in idx_types_col$character){
      comp_lda(merged_object)
    }
  }

  ## integer response ##
  if(!is.null(idx_types_col$integer)){
    #for(idx in idx_types_col$integer){
     comp_lm(merged_object)
    #}
  }


}

# link_to_clinical_data(dum_sign, dum_metadata)
#
# comp_lm(tmp_merged_compositional, 1)
#

dum_sign <- new("sign",
                id="dum1",
                id_samples=c("sam1", "sam2", "sam3"),
                id_signatures= c('s1', 's2', 's3', 's4'), ## signature names
                count_matrix=MCMCpack::rdirichlet(3, c(1,1,1,1)),
                modified=FALSE)
d <- readRDS("~/Desktop/tcga_ov_signature_exposures_20180926.rds")
many_tcga <- new("sign",
                 id='many_tcga',
                 id_samples=rownames(d),
                 id_signatures= colnames(d), ## signature names
                 count_matrix=d,
                 modified=FALSE)

## TODO: add types of response (categorical, numerical, ...)
dum_metadata <- new("metadata",
                    id="dum1",
                    id_samples=c("sam1", "sam2", "sam3"),
                    df=data.frame(age=sample(1L:100L,3),
                                  some_cat=as.character(letters[sample(1L:length(letters), 3)]),
                                  stringsAsFactors = FALSE)
)

### Example of matrix transformed into sign object
aaa <- matrix(runif(100), 4)
colnames(aaa) <- paste0('s', 1:25); rownames(aaa) <- paste0('sam', 1:4)
to_sign(aaa)

input_dummy <- matrix(runif(100), 4)
colnames(input_dummy) <- paste0('s', 1:25); rownames(input_dummy) <- paste0('sam', 1:4)
sign_dummy <- to_sign(input_dummy)
add_together_matrix(sign_dummy)
summarise(add_together_matrix(sign_dummy))

obj <- add_together_matrix(sign_dummy)
list(General=paste0("Object of class ", class(obj)),
     `Number of signatures` = ncol(obj@count_matrix),
     `Number of samples` = nrow(obj@count_matrix),
     `Geometric means of signatures`= sort(apply(obj@count_matrix, 2, function(x) exp(mean(log(x)))), decreasing=TRUE),
     `Covariance` = compositions::var.acomp(obj@count_matrix))
summarise(obj)


summarise(to_sign(aaa))
