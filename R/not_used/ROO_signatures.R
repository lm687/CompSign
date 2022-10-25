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

##' Plot ggtern (ternary plot)
plot_ggtern <- function(exposures, colours, title='CNA_12K_TCGA'){
  require(ggtern)
  stopifnot(nrow(exposures) == nrow(colours))
  stopifnot(ncol(exposures) == 3)
  tmp_df <- cbind.data.frame(exposures,
                             col=colours)
  if(is.null(colnames(exposures))){
    lab1 <- paste0('S', 1); lab2 <- paste0('S', 2); lab3 <- paste0('S', 3)
  }else{
    lab1 <- colnames(exposures)[1]; lab2 <- colnames(exposures)[2]; lab3 <- colnames(exposures)[3]
  }
  names(tmp_df) <- c('s1', 's2', 's3', 'grp')
  ggtern(data=tmp_df,
         aes(x=s1, y=s2, z=s3, grp=as.numeric(grp)), aes(x,y,z)) +
    geom_point(aes(color=grp),shape=20, alpha=1)+
    labs(x=lab1, y=lab2, z=lab3)#+
  # scale_fill_gradient(low = "blue", high="gold",
  #                     space = "Lab", na.value = "white", guide = "colourbar",
  #                     aesthetics = "colour")
  #ggtitle(paste("Cancer type: ",cancer_type, "\nSignatures: ", paste0(c(sig1, sig2, sig3), collapse=', '),
  #              "\nColour: immune score for ", gsub('_', ' ', col_immune), "\n"))
}

#' exposures: matrix of exposures, with exposures in the columns and samples in rows
#' keep_cols: columns to keep
close_data <- function(exposures, keep_cols){
  .exposures <- exposures[,keep_cols]
  sweep(.exposures, 1, rowSums(.exposures), '/')
}

#' Transform an acomp object to a matrix
acomp_to_matrix <- function(acomp_object){
  .ncol <- ncol(acomp_object)
  .res <- matrix(acomp_object, ncol=.ncol)
  if(!is.null(colnames(acomp_object))){
    colnames(.res) <- colnames(acomp_object)
  }
  if(!is.null(rownames(acomp_object))){
    rownames(.res) <- rownames(acomp_object)
  }
  .res
}

###################################################################
## Functions 6 Equivalent of cor coef with symmetric balances    ##
###################################################################

#' Compute a coefficient based on D which is used for
#' both compute_z1s() and compute_z2s()
coef_D <- function(D){
  sqrt( ((D-1)+sqrt(D*(D-2)))/(2*D) )
}

#' x: compositional vector
#' idx1: index of the first component
#' idx2: index of the second component
compute_zs <- function(x, idx1, idx2){
  ## equation 13 of Kynclova et al.
  .D <- ncol(x)
  .notidx <- (1:.D)[!((1:.D) %in% c(idx1, idx2))]
  coef <- coef_D(.D)

  ## exponents
  .expnt_a <- 1/((.D-1)+sqrt(.D*(.D-2)))
  .expnt_b <- ( sqrt(.D-2)+sqrt(.D) ) / ( sqrt(.D-2)*(.D - 1 + sqrt(.D*(.D-2))) )
  ## denominator
  denom <- ( (x[,idx2]**.expnt_a) * (apply(x[,.notidx], 1, prod))**.expnt_b )

  return(coef*log(x[idx1]/denom))
}

computeRho <- function(z1s, z2s){
  return(cov(z1s, z2s)/(sqrt(var(z1s)*var(z2s))))
}

#' x: compositional vector
#' computes the matrix of rhos for all pairwise
#' combinations
computeRhoWrapper <- function(x){
  .D <- ncol(x)
  outer(1:.D, 1:.D, Vectorize(function(i,j){
    z1s <- compute_zs(x = x, idx1 = i, idx2 = j)
    z2s <- compute_zs(x, j, i)
    computeRho(z1s, z2s)
  }))
}

#' Compute rho and plot the heatmap

plotcomputeRho <- function(x, pseudocount = 0, names_sigs, column_title='',  return_mat=FALSE){
  .mat <- computeRhoWrapper(addPseudoCounts(x, pseudocount = pseudocount))
  colnames(.mat) <- rownames(.mat) <- names_sigs
  if(return_mat){
    .mat
  }else{
    ComplexHeatmap::Heatmap(.mat,  col = circlize::colorRamp2(c(min(.mat), median(.mat), max(.mat)), c("#e6cb1f", "white", "#921fe6")),
                            column_title = column_title)
  }
}

#' Compute clr correlation and plot the heatmap
plotcomputeclrcor <- function(x, pseudocount = 0, names_sigs, column_title='', return_mat=FALSE){
  .mat <- compositions::clr(addPseudoCounts(x, pseudocount = pseudocount))
  .mat <- matrix(.mat, ncol = ncol(x))
  .mat <- cor(.mat)
  colnames(.mat) <- rownames(.mat) <- names_sigs
  if(return_mat){
    .mat
  }else{
    ComplexHeatmap::Heatmap(.mat,  col = circlize::colorRamp2(c(min(.mat), median(.mat), max(.mat)), c("#e6cb1f", "white", "#921fe6")),
                            column_title = column_title)
  }
}


#' Compute matrix of total variation
total_variation <- function(x, pseudocount = 0, remove_zeroes=FALSE){
  x <- addPseudoCounts(x, pseudocount)
  r <- outer(1:ncol(x), 1:ncol(x), Vectorize(function(i,j){
    if(remove_zeroes){
      which_keep <- which(apply(x[,c(i,j)], 1, function(k) all(k > 0)))
      if(length(which_keep) < 2){
        NA
      }else{
        subsetx <- x[which_keep,]
        var(log(subsetx[,i]/subsetx[,j]))
      }
    }else{
      var(log(x[,i]/x[,j]))
    }
  }))
  colnames(r) <- rownames(r) <- colnames(x)
  r
}

#' Remove columns and rows which contain all NA in a matrix
remove_na_columns <- function(x){
  which_rm_col <- apply(x, 2, function(i) all(is.na(i)))
  x <- x[,!which_rm_col]
  which_rm_row <- apply(x, 1, function(i) all(is.na(i)))
  x <- x[!which_rm_col,]
  x
}

give_all_combinat <- function(Nsig, exclude_complement){
  require(combinat)
  it_partitions <- c()
  if(exclude_complement){
    for(k in 1:floor(Nsig/2)){
      it_partitions <- c(it_partitions, lapply(1:ncol(combinat::combn(1:Nsig, k)), function(x) combinat::combn(1:Nsig, k)[,x]))
    }
  }else{
#    stop('Not implemented yet')
    for(k in 1:floor(Nsig/2)){
      it_partitions <- c(it_partitions, lapply(1:ncol(combinat::combn(1:Nsig, k)), function(x) combinat::combn(1:Nsig, k)[,x]))
    }
    it_partitions <- c(it_partitions, sapply(it_partitions, function(i) list(c(1:Nsig)[!(1:Nsig %in% i)])))
  }
  c(list(1:Nsig), unique(it_partitions)) ## sloppy. there are duplicates if Nsig is even (not if odd)
}

remove_some_signature <- function(mat, which_sig=1){
  if(which_sig != 'None'){
    mat <- mat[,-which_sig]
  }
  sweep(mat, 1, rowSums(mat), '/')
}

# aitch_distance <- function(x, y){
#   d <- length(x) # same as y
#   tmp <- outer(1:d, 1:d, function(i, j) (log(x[i]/x[j]) - log(y[i]/y[j]) )**2 )
#   (sum(tmp[upper.tri(tmp)]))**(1/2)
# }

aitch_distance <- function(x, y, verbose=TRUE){
  if(verbose) 'Previous results were off by a constant'
  d <- length(x) # same as y
  dist(rbind(compositions::clr(x), compositions::clr(y)), method = 'euclidean')
}

dist_Aitch <- function(x){
  require(usedist)
  dist_make(x, aitch_distance, "Aitchison Distance")
}

AitchDistphyloseq = function(physeq, ...){
  require(phyloseq)
  # PoiClaClu expects samples as rows and taxa/variables as columns
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  x = as(phyloseq::otu_table(physeq), "matrix")
  dd = usedist::dist_make(x = x, distance_fcn = CompSign::aitch_distance, method = "euclidean")
  attributes(dd)$Labels <- sample_names(physeq)
  return(dd)
}

AitchDist = function(physeq, ...){
  physeq <- t(physeq)
  x = physeq
  dd = dist_Aitch(x, ...)
  return(dd)
}

tomatrix <- function(mat){
  .ncol <- ncol(mat)
  .bool_row <- !is.null(rownames(mat))
  .bool_col <- !is.null(colnames(mat))
  if(.bool_row)  .tmp_rownames <- rownames(mat)
  if(.bool_col)  .tmp_colnames <- colnames(mat)
  mat <- matrix(sapply(mat, as.numeric), ncol=.ncol)
  if(.bool_row)  rownames(mat) <- .tmp_rownames
  if(.bool_col)  colnames(mat) <- .tmp_colnames
  mat
}

amalgamate <- function(x, which_to_amalgamate, shorten_name=FALSE, length_shorten=4){
  stopifnot(sapply(rowSums(x), all.equal, target = 1))
  require(compositions) ## tomatrix
  .cnames <- colnames(x)
  .xtmp <- cbind(x[,-which_to_amalgamate], rowSums(x[,which_to_amalgamate]))
  name_amalg <- paste0(.cnames[which_to_amalgamate], collapse = '+')
  if(shorten_name){
    name_amalg <- substr(name_amalg, 1, length_shorten)
  }
  colnames(.xtmp) <- c(.cnames[-which_to_amalgamate],
                       name_amalg)
  .xtmp
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

normalise_rw <- function(x){
  ## normalise row-wise
  sweep(x, 1, rowSums(x), '/')
}

normalise_cl <- function(x){
  ## normalise col-wise
  t(sweep(x, 2, colSums(x), '/'))
}

giveTree <- function(nleaves){
  current_tree <- lapply(1:nleaves, function(x) x)
  current_row_signary <- lapply(1:nleaves, function(x) rep(0, nleaves))
  while(length(current_tree) > 1){
    idx <- sample(1:length(current_tree), 2)
    current_tree <- c(current_tree,
                      paste0('(', current_tree[idx[1]],
                             ',', current_tree[idx[2]], ')'))
    if(is.integer(current_tree[[idx[1]]]) & is.integer(current_tree[[idx[2]]])){
      for(j in 1:2){
        current_row_signary[[idx[j]]][idx[1]] <- 1
        current_row_signary[[idx[j]]][idx[2]] <- -1
      }
    }else{
      
    }
    
    current_tree[idx] <- NULL
  }
  current_tree[[1]]
}

createSignary <- function(nleaves, tree=NULL){
  
  if(is.null(tree)){
    tmp_tree <- giveTree(nleaves)
  }else{
    tmp_tree <- tree
  }
  signary <- matrix(NA, ncol = nleaves, nrow = (nleaves-1))
  ct <- 1
  
  tmp_split <- strsplit(tmp_tree, "(?=[/(|/)/,])", perl = TRUE)[[1]]
  idx_open <- which(tmp_split == "(")
  idx_close <- which(tmp_split == ")")
  if(length(idx_close) != length(idx_open)){stop('This is not a tree')}
  
  ## find those four indices apart
  idx_immeditely <- which(outer(idx_close, idx_open, "-") == 4, arr.ind = TRUE)[1,]    
  
  while(is.na(sum(rowSums(signary)))){ ## while the signary is not completed
    ## create signary for all of these cases
    ## (whenever there is a (x,y) tree, just put 1 -1 and 0 elsewhere)
    
    X <- idx_immeditely[1];  Y <- idx_immeditely[2]
    
    tmp_group <- c(tmp_split[idx_open[Y]+1], tmp_split[idx_close[X]-1])
    signary[ct,] <- 0
    ifill <- -1
    for(i in tmp_group){
      if(nchar(i)==(1)){
        signary[ct,as.numeric(i)] <- ifill
      }else{ ## merged nodes
        signary[ct,as.numeric(strsplit(i, '&')[[1]])] <- ifill
      }
      ifill <- 1
    }
    
    if(is.na(sum(rowSums(signary)))){
      ## merge together and repeat
      if((idx_close[X]+1) > length(tmp_split)){
        ## last round
        tmp_split <- c(tmp_split[1:(idx_open[Y]-1)], paste(tmp_group, collapse = "&"), tmp_split[length(tmp_split)])
      }else{
        tmp_split <- c(tmp_split[1:(idx_open[Y]-1)], paste(tmp_group, collapse = "&"), tmp_split[(idx_close[X]+1):length(tmp_split)])
      }
      
      idx_open <- which(tmp_split == "(")
      idx_close <- which(tmp_split == ")")
      
      idx_immeditely <- which(outer(idx_close, idx_open, "-") == 4, arr.ind = TRUE)[1,]
      inverse = FALSE
      
      ct <- ct+1
    }
    
  }
  signary[nrow(signary):1,]
}

give_labels_tree= function(tree){
 .x =  strsplit(tree, "(?=[/(|/)/,])", perl = TRUE)[[1]]
 .x
}

give_idx_binary_split = function(split_labels_tree){
  ## to find where it ends, count the number of parentehesis that open and stop where the number of parenthesis that close matches it
  num_close = 0; num_open = 1; idx = 2 ## init
  while(num_close != num_open & idx <= length(split_labels_tree) ){
    .label_it = split_labels_tree[idx]
    if(.label_it == ")"){
      num_close = num_close + 1
    }else if(.label_it == "("){
      num_open = num_open + 1
    }
    idx = idx + 1
  }
  idx
}

## as_numeric_bool: are the labels numeric?
standardise_tree = function(tree, as_numeric_bool){
  stop('Deprecated. See <standardise_list_tree>')
#   ## put trees in a common shape, i.e. sorted, so that two trees can be easily compared
#   split_labels_tree = give_labels_tree(tree)
#   num_labels_tree = sum(!(split_labels_tree %in% c(")", "(", ",")))
#   num_labels_tree
#   
#   if(num_labels_tree > 2){
#     ## continue splitting further
#     
#     ## remove outer parentheses
#     split_labels_tree = split_labels_tree[-c(1, length(split_labels_tree))]
#     
#     ## open the first parenthesis, and split into two
#     idx_split = give_idx_binary_split(split_labels_tree)
#     
#     if(idx_split == (length(split_labels_tree)+1) | idx_split == (length(split_labels_tree))){
#       ## it's 1 vs all others
#       split_labels_tree_naked = split_labels_tree
#       idx_comma = which(c(split_labels_tree_naked[2], split_labels_tree_naked[length(split_labels_tree_naked)-1]) == ",")
#       part_a_full = c("(", split_labels_tree[1:(idx_comma)], ")") #paste0(c(paste0(split_labels_tree[1:(idx_comma+1)], collapse = ""),  ")"), collapse = "")
#       part_b_full = c("(", split_labels_tree[(idx_comma+3):length(split_labels_tree)], ")")
#     }else{
#       part_a_full = split_labels_tree[1:(idx_split-1)]
#       part_b_full = split_labels_tree[(idx_split+1):length(split_labels_tree)]
#     }
#     
#     ## change the order depending on where is the smallest
#     part_a = part_a_full[! (part_a_full %in% c(")", "(", ","))]
#     part_b = part_b_full[! (part_b_full %in% c(")", "(", ","))]
#     if(as_numeric_bool){
#       bool_min_in_a = min(as.numeric(part_a)) < min(as.numeric(part_b))
#     }else{
#       bool_min_in_a = min((part_a)) < min((part_b))
#     }
#     
#     cat(part_a_full, '\t\t', part_b_full, '\n')
#     
#     ## we'll make part a be the one with smallest labels
#     if(!bool_min_in_a){
#       ## need to change it
#       a = standardise_tree(paste0(part_b_full, collapse = ""), as_numeric_bool)
#       b = standardise_tree(paste0(part_a_full, collapse = ""), as_numeric_bool)
#     }else{
#       ## good as it is
#       a = standardise_tree(paste0(part_a_full, collapse = ""), as_numeric_bool)
#       b = standardise_tree(paste0(part_b_full, collapse = ""), as_numeric_bool)
#     }
#     
#     return(paste0(a,b, collapse=""))
# 
#   }else{
#     ## end
#     if(num_labels_tree == 1){
#       return(tree)
#     }else if(num_labels_tree == 2){
#       ## sort it. Here we know what the indices are as it's always in the format "(X,Y)"
#       part_a = split_labels_tree[2]
#       part_b = split_labels_tree[4]
#       if(as_numeric_bool){
#         bool_min_in_a = min(as.numeric(part_a)) < min(as.numeric(part_b))
#       }else{
#         bool_min_in_a = min((part_a)) < min((part_b))
#       }
#       
#       if(!bool_min_in_a){
#         ## change it
#         return(paste0("(", split_labels_tree[4], ",", split_labels_tree[2], ")"))
#       }else{
#         ## leave it as it is
#         return(tree)
#       }
#       
#     }
#   }
}

standardise_list_tree = function(list_tree){
  if(length(list_tree[[1]]) == 2 | length(list_tree[[2]]) == 2){ ## tree continues
    null_bool = sapply(list_tree, is.null)
    if(sum(!null_bool) == 2){
      .idx_change = (which.min(sapply(list_tree, function(i) min(unlist(i)))))-1
      return(list(standardise_list_tree(list_tree[[.idx_change+1]]), standardise_list_tree(list_tree[[(1-.idx_change)+1]])))
    }else{
      return(list_tree)
    }
  }else{
    ## it's not a list, just a value
    null_bool = sapply(list_tree, is.null)
    if(sum(!null_bool) == 2){
      .idx_change = (which.min(sapply(list_tree, function(i) min(unlist(i)))))-1
      return(list(list_tree[[.idx_change+1]], list_tree[[(1-.idx_change)+1]]))
    }else{
      return(list_tree)
    }
  }
}

tomatrix <- function(mat){
  .ncol <- ncol(mat)
  .bool_row <- !is.null(rownames(mat))
  .bool_col <- !is.null(colnames(mat))
  if(.bool_row)  .tmp_rownames <- rownames(mat)
  if(.bool_col)  .tmp_colnames <- colnames(mat)
  mat <- matrix(sapply(mat, as.numeric), ncol=.ncol)
  if(.bool_row)  rownames(mat) <- .tmp_rownames
  if(.bool_col)  colnames(mat) <- .tmp_colnames
  mat
}

give_orthonormal_basis_rw <- function(signary_rw){
  r = sum(signary_rw == 1)
  s = sum(signary_rw == (-1))
  x = sqrt((r*s)/(r+s))
  signary_rw[signary_rw == 1] = 1/r*x
  signary_rw[signary_rw == -1] = -1/s*x
  signary_rw
}
