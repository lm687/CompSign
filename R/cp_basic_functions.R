## Lena Morrill
## Copied on 12 April 2019

#' vector_dim: which dimensions to plot
# plotPCA <- function(df, groups=NA, labels=FALSE, vector_dim=c(1,2), ...){
#   require(stats); require(ggplot2)
#   if(labels) require(ggrepel)
#
#   zeroes <- which(rowSums(df) == 0)
#   nan <- apply(df, 1, function(a) any(is.nan(a)))
#   if(length(zeroes)>0 | sum(nan) > 0){
#     cat('Removing rows ', zeroes, 'which only contain zeroes\n')
#     df <- df[-zeroes,]
#     cat('Removing rows ', which(nan), 'which only contain NaNs\n')
#     df <- df[!nan,]
#   }
#
#   pca <- prcomp(df,
#                 center = TRUE,
#                 scale. = TRUE)
#   xlabel <- paste0('PC', vector_dim[1], ' (', round(pca$sdev[vector_dim[1]]**2/sum(pca$sdev**2), 3)*100, '%)')
#   ylabel <- paste0('PC', vector_dim[2],' (', round(pca$sdev[vector_dim[2]]**2/sum(pca$sdev**2), 3)*100, '%)')
#   if(!all(is.na(groups))){
#
#     if(length(unique(groups))> 8){
#       library(RColorBrewer)
#       n <- 60
#       qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#       col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#       col_vector = unique(col_vector)
#       set_swatch(col_vector)
#     }
#
#     df_plot <- data.frame(PC1=pca$x[,vector_dim[1]], PC2=pca$x[,vector_dim[2]], groups=groups)
#     out <- ggplot(data = df_plot, aes(x=PC1, y=PC2, col=groups, label=groups))+
#       geom_point()+
#       labs(x=xlabel,
#            y=ylabel)
#     if(labels) out + geom_text_repel(aes(label=groups))
#     else out
#   }else{
#     plot(pca$x[,vector_dim],
#          xlab=xlabel,
#          ylab=ylabel, ...)
#   }
# }

#' Select at random
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


#' create random signary
#' example:
#' createSignary(5, "((((1,2),3),4),5)")
createSignary <- function(nleaves, tree=NULL){
  if(is.null(tree)){
    tmp_tree <- giveTree(nleaves)
  }else{
    tmp_tree <- tree
  }
  signary <- matrix(NA, ncol = nleaves, nrow = (nleaves-1))
  ct <- 1

  tmp_split <- strsplit(tmp_tree, '')[[1]]
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
      if(nchar(i)==1){
        signary[ct,as.numeric(i)] <- ifill
      }else{ ## merged nodes
        signary[ct,as.numeric(strsplit(i, '&')[[1]])] <- ifill
      }
      ifill <- 1
    }

    if(is.na(sum(rowSums(signary)))){
      ## merge together and repeat
      tmp_split <- c(tmp_split[1:(idx_open[Y]-1)], paste(tmp_group, collapse = "&"), tmp_split[(idx_close[X]+1):length(tmp_split)])

      idx_open <- which(tmp_split == "(")
      idx_close <- which(tmp_split == ")")
      idx_immeditely <- which(outer(idx_close, idx_open, "-") == 4, arr.ind = TRUE)[1,]
      ct <- ct+1
    }

  }
  signary
}

#' clean dataset based on some metadata missing etc
clean_dataset <- function(metadata_name, expected_categories, exposures_matrix, metadata_matrix){
  group_vector <- metadata_matrix[,metadata_name]
  group_vector_bool <- group_vector %in% expected_categories
  group_vector_bool_clean <- group_vector[group_vector_bool]
  exposures_clean <- exposures_matrix[group_vector_bool,]
  list(exposures_clean=exposures_clean, group_vector_bool_clean=group_vector_bool_clean)
}

remove_some_signature <- function(mat, which_sig=1){
  if(which_sig != 'None'){
    mat <- mat[,-which_sig]
  }
  sweep(mat, 1, rowSums(mat), '/')
}
