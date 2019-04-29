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
         ))

setClass("metadata",
         representation = representation(
           id = "character",
           id_samples = "character",
           df = "data.frame"
         ))

#' Class whiuch inludes features of both metadata and
#'  mutational signatures.
#'  types_metadata is a vector indicating which type of data each
#'  column in df (metadata) is.
setClass("merged_compositional",
         representation = representation(
           id="character",
           id_samples="character",
           id_signatures = "character",
           count_matrix = "matrix",
           df= "data.frame",
           pseudocounts="logical",
           types_metadata = "character",
           bool_empty_metadata = "logical"
         ))

## alternative:
## foo <- function(x) {
##   if (!is.numeric(x)) stop("X must be numeric")
##   structure(list(x), class = "foo")
## }

#' Converts a matrix to a sign object
to_sign <- function(x){
  x_sign_obj <- new("sign",
      id=deparse(substitute(x)), ## get name
      id_samples=rownames(x),
      id_signatures= colnames(x),
      count_matrix=x,
      modified=FALSE)
  x_sign_obj
}

#' Convert a 'sign' object to a 'merged_compositional_to_sign'
#' (note loss of information)
merged_compositional_to_sign <- function(x){
  x_sign_obj <- new("sign",
                    id=x@id, ## get name
                    id_samples=x@id_samples,
                    id_signatures= x@id_signatures,
                    count_matrix=x@count_matrix)
  x_sign_obj
}

###################################
############# CURRENT #############
#' It computes the LDA for categorical responses (??? need to change? they are actually the predictors)
comp_lda <- function(x, indices_response){
  if(class(x)[1] != 'merged_compositional') stop('Input must be of class <merged_compositional>')
  model_compReg <- lm(compositions::ilr(compositions::acomp(x@count_matrix))~as.matrix((x@df)[,indices_predictor]))
  list(summary(model_compReg))
}
############# CURRENT #############
###################################
#' (ongoing) Linear regression with compositional data as the response
#' It computes a linear regression with some numerical value as the predictor and the compositions as response
comp_lm <- function(x, indices_predictor){
  ## the composition is the response

  if(class(x)[1] != 'merged_compositional') stop('Input must be of class <merged_compositional>')

  #require(compositions)
  model_compReg <- lm(compositions::ilr(compositions::acomp(x@count_matrix))~as.matrix((x@df)[,indices_predictor]))
  list(summary(model_compReg)#,
       #anova(model_compReg)
       )
}

#' Perform logistic regression
#'
#' @param merged_obj a merged_compositional object
#' @param colname either the name or the index of the column in metadata(merged_compositional) with which to perform logistic regression
#' @return The summary of a logistic regression
deprecated_comp_logistic <- function(merged_obj, colname){
  lab <- (as.numeric(as.factor(metadata(merged_obj)[,colname]))-1)
  if(all(sort(lab) != c(0,1))){stop('Check labels of metadata!')}
  summary(glm(formula = lab ~ ilr(acomp(count_matrix(merged_obj))),
              family = binomial(link = "logit")))
}

#' Performs logistic regression, as explained in the reference
#' @param m_obj merged_compositional object
#' @param col_idx name or index of column of interest in metadata(m_obj)
#' @references Van den Boogaart, K. Gerald, and Raimon Tolosana-Delgado. Analyzing compositional data with R. Vol. 122. Berlin: Springer, 2013.
comp_logistic <- function(compData, binaryLabels, relax_binary_assumption=FALSE){
  require(nnet)
  ## potential next: arguments are a merged_compositional object and the column

  ## transform if necessary
  if((length(unique(binaryLabels))>2) & !relax_binary_assumption){
    stop('There must be only two labels. Use cleanObject() if necessarys')
  }
  dat <- data.frame(binaryLabels, compData)
  ncol(compData)
  if(is.null(colnames(compData))){
    frm <- paste0(deparse(substitute(binaryLabels)), '~', paste(paste0('X',1:ncol(compData)), collapse = '+'))
  }else{
    frm <- paste0(deparse(substitute(binaryLabels)), '~', paste(colnames(compData), collapse = '+'))
  }
  res <- multinom(formula = as.formula(frm), data = dat)
  list(coefTransformed=coef(res),
       ## coefsSimplex=ilr2clr(coef(res)[-1], x=compData), ## need to debugged
       summary=summary(res),
       res=res,
       table(as.factor(binaryLabels), FP_FN_table=as.factor(round(fitted(res)[,1]))))
}


#' To be added in the package
#' @param merged_obj a merged object
#' @param colname column of interest
#' @return a merged object of smaller size, with only the rows in which the column of interest is part of the binary classification
cleanObject <- function(merged_obj, colname, expected_labels='', verbose=FALSE){
  tmp <- merged_obj@df[,colname]
  table(tmp)
  if(expected_labels==''){
    if(verbose) cat("No expected labels specified.\nLabels used: ")
    used <- unique(unique(tmp))
    not_used <- c('--', 'no_data_supplied', 'not reported')
    used <- used[! (used %in% not_used)]
    if(verbose) cat(paste0(used, collapse=', ' ))
    if(verbose) cat("\nLabels not used: ")
    if(verbose) cat(paste0(not_used, collapse=', ' ))
    if(verbose) cat("\n")
  }
  merged_obj@df <- merged_obj@df[tmp %in% used,]
  merged_obj@count_matrix <- merged_obj@count_matrix[tmp %in% used,]
  merged_obj
}

###########################################
##### Functions to manipulate objects #####
###########################################

#' Retrieve count matrix from object of class sign.
#' @param X object of class sign
count_matrix <- function(X){
  return(X@count_matrix)
}

#' Assign count matrix from object of class sign
setGeneric("count_matrix<-",
           function(object, value) standardGeneric("count_matrix<-"))
setReplaceMethod("count_matrix",
                 signature(object="sign",
                           value="matrix"),
                 function(object, value){
                   object@count_matrix <- value
                   object@id_samples <- rownames(value)
                   object@id_signatures <- colnames(value)
                   return(object)
                 })
setReplaceMethod("count_matrix",
                 signature(object="merged_compositional",
                           value="matrix"),
                 function(object, value){
                   object@count_matrix <- value
                   object@id_samples <- rownames(value)
                   object@id_signatures <- colnames(value)
                   return(object)
                 })

#' Add name to object of class sign or merged_compositional
setGeneric("id<-",
           function(object, value) standardGeneric("id<-"))
setReplaceMethod("id",
                 signature(object="sign",
                           value="character"),
                 function(object, value){
                   object@id <- value
                   return(object)
                 })
setReplaceMethod("id",
                 signature(object="merged_compositional",
                           value="character"),
                 function(object, value){
                   object@id <- value
                   return(object)
                 })

#' Retrieve metadata from object of class sign.
#' @param X object of class sign
metadata <- function(X){
  return(X@df)
}

#' Retrieve names of metadata from object of class sign.
#' Used to be called metadata_types()
#' @param X object of class sign
names_metadata <- function(X){
  return(colnames(X@df))
}

#' Assign metadata to object of class sign
setGeneric("metadata<-",
           function(object, value) standardGeneric("metadata<-"))

setReplaceMethod("metadata",
                 signature(object="merged_compositional",
                           value="data.frame"),
                 function(object, value){
                   object@df <- value
                   return(object)
                 })


#' Get types of metadata, as a char vector
types_metadata <- function(mergedComp_obj){
  mergedComp_obj@types_metadata
}

#' Add pseudocounts; renormalise
addPseudoCounts <- function(count_matrix, pseudocount=1e-5){
  count_matrix <- count_matrix + pseudocount
  count_matrix <- sweep(count_matrix, 1,
                        rowSums(count_matrix), '/')
  cat('Pseudocount of', pseudocount, 'added\n')
  count_matrix
}

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
setGeneric("within_signature_analysis", function(obj, ...) standardGeneric("within_signature_analysis"))
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

#' Compare two groups
compare_populations <- function(predictors, response, ...){
  if(length(unique(response))>2){stop('Only two categories in response')}
  require(Compositional)
  tmp_response <- response
  tmp_response <- factor(tmp_response)
  levels(tmp_response) <- c(1,2)
  Compositional::comp.test(x = predictors[,-1],
                           ina = tmp_response, ...)
}

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
  ## or mauybe it should be logistic regression?
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


###################################################################
## Function 5 comparison with clinical data: dendrogram of samples
###################################################################

#' Create a dendrogram, using Aithchison distance, of the samples in a merged
#' object. The labels are coloured according to the one of the columns of its
#' metadata dataframe (to be specified in name_clinical). There is also an
#' option (bool_comparison) to add a second dendrogram using Euclidean distance,
#' for comparison. Arguments for plot() are inherited (not tested).
#' WARNING! colours for the second dendrogram need fixing
createDendrogram <- function(merged_object, name_clinical, bool_comparison, ...){
  require(dendextend)
  require(RColorBrewer)
  # to add: possibility to extend colour vector
  col_vec <- brewer.pal(n = 8, name = "Set2")
  tmp_toplot_dendro <- as.dendrogram(hclust(dist(acomp(merged_object@count_matrix))))
  labels_colors(tmp_toplot_dendro) <- col_vec[as.factor(merged_object@df[,name_clinical][labels(tmp_toplot_dendro)])]
  labels(tmp_toplot_dendro) <- rep('.', nrow(merged_object@df))
  if(bool_comparison) par(mfrow=c(1,2))
  plot(tmp_toplot_dendro, main='Aitchison distance')
  if(bool_comparison){
    tmp_toplot_dendro_nonComp <- as.dendrogram(hclust(dist(merged_object@count_matrix)))
    labels_colors(tmp_toplot_dendro_nonComp) <- col_vec[as.factor(merged_object@df[,name_clinical][labels(tmp_toplot_dendro_nonComp)])]
    labels(tmp_toplot_dendro_nonComp) <-  rep('.', nrow(merged_object@df))
    plot(tmp_toplot_dendro_nonComp, main='Euclidean distance')
  }
}

#########################################
############### DEBUGGING ###############
#########################################
#
# link_to_clinical_data(dum_sign, dum_metadata)
#
# comp_lm(tmp_merged_compositional, 1)
#
#
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
#
# ## Example of matrix transformed into sign object
# aaa <- matrix(runif(100), 4)
# colnames(aaa) <- paste0('s', 1:25); rownames(aaa) <- paste0('sam', 1:4)
# to_sign(aaa)
#
# input_dummy <- matrix(runif(100), 4)
# colnames(input_dummy) <- paste0('s', 1:25); rownames(input_dummy) <- paste0('sam', 1:4)
# sign_dummy <- to_sign(input_dummy)
# add_together_matrix(sign_dummy)
# summarise(add_together_matrix(sign_dummy))
#
# summarise(to_sign(aaa))
#
#
# ## example of comp_lm
# tmp_merged_compositional <- new("merged_compositional",
#                                 id='adas',
#                                 id_samples=c("sam1", "sam2", "sam3"),
#                                 id_signatures= c('s1', 's2', 's3', 's4'), ## signature names
#                                 count_matrix=MCMCpack::rdirichlet(3, c(1,1,1,1)),
#                                 df=data.frame(a=c(3,4,1), b=c(10, 10, 10)))
# comp_lm(tmp_merged_compositional)
