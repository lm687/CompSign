#' pheatmap without clustering or names
pheatmap0 <- function(x, noclust=TRUE, nonames=TRUE, ...){
  if(noclust){
    if(nonames){
      pheatmap(x, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, ...)
    }else{
      pheatmap(x, cluster_rows = FALSE, cluster_cols = FALSE, ...)
    }
  }else{
    if(noclust){
      pheatmap(x, show_rownames = FALSE, show_colnames = FALSE, ...)
    }else{
      pheatmap(x, ...)
    }
  }
}

#' from a matrix of cosine similarities (both dimensions have to coincide), pair them.
#' This has to be done jointly with the function below sequential_to_original()
pair_signatures <- function(xxxxx, nsignatures){
  working_xxxxx <- xxxxx
  list_pairs <- c()
  k <- 1
  while(k < nsignatures){
    wm <- which.max(working_xxxxx) ## bycol

    ## true signatures are in the y axis; NMF-derived are in the x axis
    idx1 <- (wm %% nrow(working_xxxxx))    ## true signatures
    if(idx1 == 0){idx1 <- nrow(working_xxxxx)}
    idx2 <- (ceiling(wm/(nrow(working_xxxxx)))) ## NMF-derived

    list_pairs <- rbind(list_pairs, c(idx1, idx2))
    if(!all(dim(working_xxxxx) == c(2,2))){
      working_xxxxx <- working_xxxxx[-idx1,]
      working_xxxxx <- working_xxxxx[,-idx2]
    }
    k <- k + 1
  }
  list_pairs <- rbind(list_pairs, c(1, 1))
  list_pairs
}

#' from the sequential pairing of signatures, get the indices of the original pairs
sequential_to_original <- function(list_pairs){
  nsignatures <- nrow(list_pairs)
  list_pairs_new <- matrix(NA, nrow = nsignatures, ncol = 2)
  list_pairs_new[1,] <- list_pairs[1,]
  remaining <- list()
  for(rw in 1:nrow(list_pairs)){
    remaining[[rw]] <- list()
    if(rw == 1){
      for(cl in 1:ncol(list_pairs)){
        remaining[[rw]][[cl]] <- 1:nsignatures
        remaining[[rw]][[cl]] <- remaining[[rw]][[cl]][-list_pairs[rw,cl]]
      }
    }else{
      for(cl in 1:ncol(list_pairs)){
        #print('here')
        remaining[[rw]][[cl]] <- remaining[[rw-1]][[cl]]
        cur <- remaining[[rw]][[cl]][list_pairs[rw,cl]] ## current
        remaining[[rw]][[cl]] <- remaining[[rw]][[cl]][-list_pairs[rw,cl]]
        if(is.na(cur)){stop()}
        list_pairs_new[rw,cl] <- cur
      }
    }
  }
  list_pairs_new
}

cosineSimilarity <- (function(A, B, verbatim=TRUE){
  A <- as.vector(A); B <- as.vector(B)
  if(verbatim) cat('Note that previous versions of the code might have been erroneous for input vectors of class <acomp>\n')
  ## A, B vectors
  sum(A*B)/sqrt(sum(A^2)*sum(B^2))
})

wrapperCosineSimilarity <- Vectorize(function(X, Y, matrixA, matrixB){
  ## X, Y indices
  cosineSimilarity(matrixA[X,], matrixB[Y,])
}, vectorize.args = c('X', 'Y'))

outerCosineSimilarity <- function(matrixA, matrixB){
  ## given two matrices, where the rows are the signatures and the columns the categories of mutations
  ## with which we define the signatures, output a matrix which has the cosine similarity for any pair
  ## of mutational signatures
  print(dim(matrixA))
  print(dim(matrixB))
  if(ncol(matrixA) != ncol(matrixB)){stop('Matrices need to have the same number of columns')}
  return(outer(1:nrow(matrixA),
               1:nrow(matrixB),
               (wrapperCosineSimilarity), matrixA=matrixA, matrixB=matrixB ))
}


outerNNLS <- function(matrixA, matrixB){
  ## given two matrices, where the rows are the signatures and the columns the categories of mutations
  ## with which we define the signatures, output a matrix which has the NNLS for any pair
  ## of mutational signatures. Based on outerCosineSimilarity
  require(nnls)
  if(ncol(matrixA) != ncol(matrixB)){stop('Matrices need to have the same number of columns')}
  return(do.call('rbind', lapply(1:nrow(matrixB),
               function(Y) coef(nnls::nnls(t(matrixA), matrixB[Y,]) ))))
}
